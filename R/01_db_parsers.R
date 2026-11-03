#' Import compound records from an NIST/MoNA-style `.msp` text file
#'
#' Reads the entire file with [readLines()], splits records on blank lines or
#' lines starting with `Name:`, extracts key metadata fields, parses `Num
#' Peaks:` regions into paired m/z and intensity values, and returns tables
#' ready for [build_nts_database()].
#'
#' When `Ion_mode:` / `Ionization_mode:` / `Polarity:` is present, records whose
#' mode does not match `polarity` are dropped. If absent, `polarity` is
#' assigned to every kept record.
#'
#' Neutral `mass` is derived from `PrecursorMZ` (adduct m/z) using `polarity`,
#' or taken from `Precursor_mass` / `ExactMass` when present. When
#' `PrecursorMZ` is available, `target_mz` is also filled for downstream use.
#'
#' Ion mode lines may appear as \code{Ion_mode:}, \code{IONMODE:}, \code{IonMode:},
#' etc. Precursor adduct lines may appear as \code{Precursor_type:},
#' \code{PRECURSORTYPE:}, or \code{PrecursorType:} (case-insensitive).
#'
#' `SMILES` is read from dedicated `SMILES:` / `Smiles:` lines when present;
#' otherwise MoNA-style quoted tokens in `Comments:` are parsed (`"SMILES=..."`
#' preferred over `"computed SMILES=..."`).
#'
#' Additional columns when present: `Formula`, `Precursor_type`,
#' `Instrument_type`, `Collision_energy` (from key lines or quoted
#' `"collision energy=..."` in `Comments:`), `CAS`, `PubChem_CID` (from quoted
#' tokens in the block text).
#'
#' @section mzCloud support:
#' When importing mzCloud offline library MSP exports, the following
#' additional fields are recognized:
#' \itemize{
#'   \item \code{COMPOUNDCLASS:} -> compound class annotation
#'   \item \code{SMILES:} -> SMILES structure
#'   \item \code{INCHIKEY:} -> InChIKey identifier
#'   \item \code{INSTRUMENT:} -> instrument name (e.g., "Q Exactive Plus Orbitrap")
#'   \item \code{MASSANALYZER:} -> mass analyzer type (e.g., "FT")
#'   \item \code{FRAGMENTATIONMODE:} -> fragmentation mode (e.g., "HCD", "CID")
#'   \item \code{COLLISIONENERGY:} / \code{Collision_energy:} -> collision energy (eV)
#'   \item \code{NEUTRALMASS:} -> neutral exact mass (overrides derived mass)
#'   \item \code{SYNONYM:} -> compound synonyms
#' }
#'
#' @param msp_file Path to the `.msp` file.
#' @param polarity `"NEG"`, `"POS"`, or `"EI"` - filter when ion mode is stated;
#'   default assigned when not. `"EI"` skips ion-mode filtering (for GC-EI).
#'
#' @return A named `list` with `metadata` and `spectra_long` (`data.table`s).
#'   `metadata` includes `Lab_ID`, `name`, `mass`, `Polarity`, `SMILES`,
#'   `InChIKey`, `Formula`, `Precursor_type`, `Collision_energy`,
#'   `Instrument_type`, `CAS`, `PubChem_CID` (missing values are `NA`), plus
#'   `target_mz` when precursor m/z was parsed. For mzCloud imports additional
#'   columns include `Compound_class`, `Fragmentation_mode`, `Instrument`,
#'   `Mass_analyzer`, `Synonyms`, `Spectrum_type`. `spectra_long` has `Lab_ID`,
#'   `mz`, `intensity`.
#'
#' @examples
#' \dontrun{
#' # Import mzCloud offline library MSP export (LC-MS/MS, ESI positive)
#' mzcloud_imp <- import_msp("path/to/mzCloud_Reference_2023.msp",
#'                            polarity = "POS")
#' # Build reference database
#' mzcloud_lib <- build_nts_database(mzcloud_imp$metadata,
#'                                    mzcloud_imp$spectra_long,
#'                                    version_str = "2023.1")
#' }
#' @export
import_msp <- function(msp_file, polarity = "NEG") {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  pol_user <- toupper(as.character(polarity)[1L])
  if (!pol_user %in% c("NEG", "POS", "EI")) {
    stop("`polarity` must be \"NEG\", \"POS\", or \"EI\".", call. = FALSE)
  }

  msp_file <- path.expand(as.character(msp_file)[1L])
  .nts_check(msp_file, "file_exists", .arg = "msp_file")

  message("import_msp: reading ", msp_file)
  
  dt <- data.table::fread(
    msp_file, sep = "", header = FALSE, col.names = "line",
    blank.lines.skip = TRUE, quote = "", strip.white = FALSE,
    encoding = "UTF-8", stringsAsFactors = FALSE
  )
  
  if (!nrow(dt)) {
    message("import_msp: empty file.")
    return(list(
      metadata = data.table::data.table(),
      spectra_long = data.table::data.table(Lab_ID = character(), mz = numeric(), intensity = numeric())
    ))
  }
  
  message("import_msp: identifying blocks...")
  is_name <- grepl("^Name\\s*:", dt$line, ignore.case = TRUE)
  dt[, block_id := cumsum(is_name)]
  dt <- dt[block_id > 0L]
  
  if (!nrow(dt)) {
    message("import_msp: no records kept.")
    return(list(
      metadata = data.table::data.table(),
      spectra_long = data.table::data.table(Lab_ID = character(), mz = numeric(), intensity = numeric())
    ))
  }
  
  message("import_msp: extracting metadata...")
  # Filter peak lines to massively reduce regex workload (from 50M rows down to 1M rows)
  is_peak <- grepl("^[ \t]*[0-9]", dt$line)
  peaks <- dt[is_peak]
  dt_meta <- dt[!is_peak]
  
  # Free up memory to prevent swap hanging
  rm(dt, is_peak, is_name)
  gc(FALSE)
  
  extract_kv <- function(pattern) {
    idx <- grep(pattern, dt_meta$line, ignore.case = TRUE)
    if (!length(idx)) return(data.table::data.table(block_id = integer(), val = character()))
    val <- sub(paste0("(?i)", pattern, "\\s*"), "", dt_meta$line[idx], perl = TRUE)
    res <- data.table::data.table(block_id = dt_meta$block_id[idx], val = trimws(val))
    res <- res[!duplicated(block_id)]
    res
  }
  
  names_dt <- extract_kv("^Name\\s*:")
  ion_dt <- extract_kv(
    "^(Ionization_mode|Ionization|Ion[_ ]?mode|Polarity)\\s*:"
  )
  prec_mz_dt <- extract_kv("^(PrecursorMZ|Precursor_m/z|precursor_mz|PRECURSORMZ)\\s*:")
  prec_mass_dt <- extract_kv("^(Precursor_mass|PrecursorMass|ExactMass|EXACTMASS|MW|Monoisotopic_mass)\\s*:")
  lab_dt <- extract_kv("^(DB#|DBNO|DBno|SpectrumID|SPECTRUMID|Spectrum\\s+ID|MSP_FILE|FILENAME|ID)\\s*:")
  smi_dt <- extract_kv("^(Smiles|SMILES)\\s*:")
  iki_dt <- extract_kv("^(InChIKey|INCHIKEY)\\s*:")
  form_dt <- extract_kv("^Formula\\s*:")
  ## MoNA: "Precursor_type:", "Precursor type:"; many exports: "PRECURSORTYPE:" / "PrecursorType:"
  pt_dt <- extract_kv("^(?i)precursor[_\\s-]*type\\s*:")
  it_dt <- extract_kv("^Instrument[-_ ]type\\s*:")
  ce_dt <- extract_kv("^(Collision[-_ ]?energy|COLLISIONENERGY|CE)\\s*:")
  cas_dt <- extract_kv("^CAS#?\\s*:")
  rt_dt <- extract_kv("^(RetentionTime|Retention[_ ]?time|RT)\\s*:")
  ## mzCloud-specific fields
  neutmass_dt <- extract_kv("^(NEUTRALMASS|NEUTRAL[_ ]?MASS)\\s*:")
  cc_dt <- extract_kv("^COMPOUNDCLASS\\s*:")
  fm_dt <- extract_kv("^FRAGMENTATIONMODE\\s*:")
  instr_dt <- extract_kv("^INSTRUMENT\\s*:")
  ma_dt <- extract_kv("^MASSANALYZER\\s*:")
  syn_dt <- extract_kv("^SYNONYM\\s*:")
  st_dt <- extract_kv("^SPECTRUMTYPE\\s*:")
  pcid_dt <- extract_kv("^PUBCHEM\\s*:")

  cm_idx <- grep("^Comments\\s*:", dt_meta$line, ignore.case = TRUE)
  cm_dt <- data.table::data.table(block_id = integer(), smi = character(), ce = character(), cas = character(), pcid = character())
  if (length(cm_idx)) {
    cm_lines <- dt_meta$line[cm_idx]
    cm_b <- dt_meta$block_id[cm_idx]
    
    get_quoted <- function(lines, regex) {
      m <- regexec(regex, lines, perl = TRUE)
      g <- regmatches(lines, m)
      vapply(g, function(x) if (length(x) >= 2L) trimws(x[2L]) else NA_character_, "")
    }
    
    smi_c1 <- get_quoted(cm_lines, "(?i)\"SMILES=([^\"]+)\"")
    smi_c2 <- get_quoted(cm_lines, "(?i)\"computed\\s+SMILES=([^\"]+)\"")
    smi_c <- ifelse(!is.na(smi_c1), smi_c1, smi_c2)
    
    ce_c <- get_quoted(cm_lines, "(?i)\"collision\\s+energy=([^\"]+)\"")
    cas_c <- get_quoted(cm_lines, "(?i)\"cas=([^\"]+)\"")
    pcid_c <- get_quoted(cm_lines, "(?i)\"pubchem\\s+cid=([^\"]+)\"")
    
    cm_dt <- data.table::data.table(block_id = cm_b, smi = smi_c, ce = ce_c, cas = cas_c, pcid = pcid_c)
  }
  
  meta <- data.table::data.table(block_id = unique(dt_meta$block_id))
  
  meta <- merge(meta, names_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "name")
  
  meta <- merge(meta, ion_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "ion_mode")
  
  meta[, Polarity := pol_user]
  if (pol_user != "EI") {
  meta[!is.na(ion_mode) & grepl("neg|negative|n\\-|anion|^n$", ion_mode, ignore.case=TRUE), Polarity := "NEG"]
      meta[!is.na(ion_mode) & grepl("pos|positive|p\\+|cation|^p$", ion_mode, ignore.case=TRUE), Polarity := "POS"]
      
      meta <- meta[Polarity == pol_user]
  }
  if (!nrow(meta)) {
    message("import_msp: no records kept.")
    return(list(
      metadata = data.table::data.table(),
      spectra_long = data.table::data.table(Lab_ID = character(), mz = numeric(), intensity = numeric())
    ))
  }
  
  meta <- merge(meta, prec_mz_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "prec_mz")
  meta[, prec_mz := suppressWarnings(as.numeric(prec_mz))]
  
  meta <- merge(meta, prec_mass_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "prec_mass")
  meta[, prec_mass := suppressWarnings(as.numeric(prec_mass))]
  
  proton <- 1.007276
  meta[, mass := NA_real_]
  meta[is.finite(prec_mass), mass := prec_mass]
  meta[is.na(mass) & is.finite(prec_mz) & Polarity == "NEG", mass := prec_mz + proton]
  meta[is.na(mass) & is.finite(prec_mz) & Polarity == "POS", mass := prec_mz - proton]
  
  meta <- meta[is.finite(mass)]
  if (!nrow(meta)) {
    message("import_msp: no records kept.")
    return(list(
      metadata = data.table::data.table(),
      spectra_long = data.table::data.table(Lab_ID = character(), mz = numeric(), intensity = numeric())
    ))
  }
  
  meta <- merge(meta, lab_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Lab_ID")
  meta[is.na(Lab_ID) | !nzchar(Lab_ID), Lab_ID := sprintf("MoNA_%05d", block_id)]
  
  meta <- merge(meta, smi_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "SMILES")
  meta <- merge(meta, iki_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "InChIKey")
  meta <- merge(meta, form_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Formula")
  meta <- merge(meta, pt_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Precursor_type")
  meta <- merge(meta, it_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Instrument_type")
  meta <- merge(meta, ce_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Collision_energy")
  meta <- merge(meta, cas_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "CAS_standalone")
  meta <- merge(meta, rt_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "RT")

  ## mzCloud: neutral mass (used to correct/override derived mass)
  meta <- merge(meta, neutmass_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "neutral_mass_raw")
  meta[, neutral_mass_raw := suppressWarnings(as.numeric(neutral_mass_raw))]
  meta[is.finite(neutral_mass_raw) & neutral_mass_raw > 0, mass := neutral_mass_raw]

  ## mzCloud: compound metadata
  meta <- merge(meta, cc_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Compound_class")
  meta <- merge(meta, fm_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Fragmentation_mode")
  meta <- merge(meta, instr_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Instrument")
  meta <- merge(meta, ma_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Mass_analyzer")
  meta <- merge(meta, syn_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Synonyms")
  meta <- merge(meta, st_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "Spectrum_type")
  meta <- merge(meta, pcid_dt, by="block_id", all.x=TRUE)
  data.table::setnames(meta, "val", "PubChem_line")
  
  if (nrow(cm_dt)) {
    meta <- merge(meta, cm_dt, by="block_id", all.x=TRUE)
    meta[is.na(SMILES) & !is.na(smi), SMILES := smi]
    meta[is.na(Collision_energy) & !is.na(ce), Collision_energy := ce]
    meta[, CAS_from_cm := cas]
    meta[, PubChem_CID := pcid]
    meta[, c("smi", "ce", "cas", "pcid") := NULL]
  } else {
    meta[, CAS_from_cm := NA_character_]
    meta[, PubChem_CID := NA_character_]
  }
  # Prefer standalone CAS# field over Comments-embedded CAS
  meta[, CAS := NA_character_]
  meta[!is.na(CAS_standalone) & nzchar(CAS_standalone), CAS := CAS_standalone]
  meta[is.na(CAS) & !is.na(CAS_from_cm) & nzchar(CAS_from_cm), CAS := CAS_from_cm]
  meta[, CAS_standalone := NULL]
  meta[, CAS_from_cm := NULL]

  meta[is.na(PubChem_CID) & !is.na(PubChem_line) & nzchar(PubChem_line), PubChem_CID := PubChem_line]
  meta[, PubChem_line := NULL]
  
  meta[, target_mz := prec_mz]
  
  keep_cols <- c("block_id", "Lab_ID", "name", "mass", "Polarity", "SMILES", "InChIKey",
                 "Formula", "Precursor_type", "Collision_energy", "Instrument_type",
                 "CAS", "PubChem_CID", "target_mz", "RT",
                 "Compound_class", "Fragmentation_mode", "Instrument",
                 "Mass_analyzer", "Synonyms", "Spectrum_type")
  meta <- meta[, ..keep_cols]
  
  message("import_msp: parsing peaks...")
  
  if (nrow(peaks)) {
    parsed_peaks <- data.table::fread(
      text = peaks$line, select = 1:2, col.names = c("mz", "intensity"),
      fill = TRUE, na.strings = c("NA", ""), quote = ""
    )
    parsed_peaks[, block_id := peaks$block_id]
    
    spec <- merge(parsed_peaks, meta[, .(block_id, Lab_ID)], by="block_id", all.x=TRUE)
    spec <- spec[!is.na(mz) & !is.na(intensity), .(Lab_ID, mz, intensity = as.numeric(intensity))]
    
    # Free peak memory
    rm(peaks, parsed_peaks)
    gc(FALSE)
  } else {
    spec <- data.table::data.table(Lab_ID=character(), mz=numeric(), intensity=numeric())
  }
  
  meta[, block_id := NULL]
  
  message(
    "import_msp: metadata rows = ", nrow(meta),
    "; spectra_long rows = ", nrow(spec)
  )
  
  list(metadata = meta, spectra_long = spec)
}

#' @rdname import_msp
#' @export
import_mona_msp <- import_msp

#' Import peak lists from MGF (`BEGIN IONS` / `END IONS`)
#'
#' Reads the file once, locates ion blocks with vectorized `grepl`, and parses
#' peak numbers per block with a single [scan()] over concatenated peak lines
#' (fast for ~1e5+ blocks). `TITLE=` (case-insensitive) supplies `Lab_ID`;
#' missing titles become `MGF_00000001`, etc.
#'
#' @param mgf_file Path to the `.mgf` file.
#'
#' @return A `data.table` with columns `Lab_ID`, `mz`, `intensity`.
#'
#' @export
import_mgf <- function(mgf_file) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  mgf_file <- path.expand(as.character(mgf_file)[1L])
  .nts_check(mgf_file, "file_exists", .arg = "mgf_file")

  message("import_mgf: reading ", mgf_file)
  
  # Read entire file as single column using C-level fread
  dt <- data.table::fread(
    mgf_file, sep = "", header = FALSE, col.names = "line",
    blank.lines.skip = TRUE, quote = "", strip.white = FALSE,
    encoding = "UTF-8", stringsAsFactors = FALSE
  )
  if (!nrow(dt)) {
    return(data.table::data.table(Lab_ID = character(), mz = numeric(), intensity = numeric()))
  }

  # Fast vectorized pattern matching
  is_begin <- grepl("^BEGIN IONS", dt$line, ignore.case = TRUE)
  dt[, block_id := cumsum(is_begin)]
  
  # Remove preamble before first block
  dt <- dt[block_id > 0L]
  if (!nrow(dt)) {
    return(data.table::data.table(Lab_ID = character(), mz = numeric(), intensity = numeric()))
  }
  
  # Extract TITLE lines for Lab_ID mapping
  is_title <- grepl("^TITLE\\s*=", dt$line, ignore.case = TRUE)
  titles <- dt[is_title, .(block_id, Lab_ID = trimws(sub("^TITLE\\s*=\\s*", "", line, ignore.case = TRUE)))]
  
  # Filter peak lines: drop BEGIN, END, and generic key=value lines
  is_kv <- grepl("^[A-Za-z][A-Za-z0-9_]*[[:space:]]*=", dt$line)
  is_end <- grepl("^END IONS", dt$line, ignore.case = TRUE)
  
  peaks <- dt[!(is_begin | is_end | is_kv)]
  if (!nrow(peaks)) {
    message("import_mgf: no spectra extracted.")
    return(data.table::data.table(Lab_ID = character(), mz = numeric(), intensity = numeric()))
  }
  
  # Parse mz and intensity instantly using fread on the character vector
  parsed_peaks <- data.table::fread(
    text = peaks$line, select = 1:2, col.names = c("mz", "intensity"),
    fill = TRUE, na.strings = c("NA", ""), quote = ""
  )
  parsed_peaks[, block_id := peaks$block_id]
  
  # Assemble the final mapping
  blocks_all <- data.table::data.table(block_id = unique(dt$block_id))
  titles_all <- merge(blocks_all, titles, by = "block_id", all.x = TRUE)
  titles_all[is.na(Lab_ID) | !nzchar(Lab_ID), Lab_ID := sprintf("MGF_%08d", block_id)]
  
  out <- merge(parsed_peaks, titles_all, by = "block_id", all.x = TRUE)
  out <- out[!is.na(mz) & !is.na(intensity), .(Lab_ID, mz, intensity = as.numeric(intensity))]
  
  message(
    "import_mgf: blocks parsed = ", nrow(blocks_all),
    "; peak rows = ", nrow(out)
  )
  
  out[]
}
