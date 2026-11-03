#' @keywords internal
#' @noRd
.nts_resolve_csv_col <- function(nm, candidates) {
  if (!length(nm)) {
    return(NULL)
  }
  ln <- tolower(nm)
  for (cand in candidates) {
    w <- match(tolower(cand), ln)
    if (!is.na(w)) {
      return(nm[w])
    }
  }
  NULL
}

#' @keywords internal
#' @noRd
.nts_spectraverse_normalize_polarity_vec <- function(pol) {
  pa <- as.character(pol)
  empty <- is.na(pa) | !nzchar(trimws(pa))
  x <- toupper(pa)
  x[grepl("NEG", x, ignore.case = TRUE)] <- "NEG"
  x[grepl("POS", x, ignore.case = TRUE)] <- "POS"
  x[!x %in% c("POS", "NEG")] <- "NEG"
  x[empty] <- "NEG"
  x
}

#' Normalise an ion-mode vector to `"POS"` / `"NEG"`, or fail loudly.
#'
#' Accepts the spellings that appear in vendor/metadata exports
#' (`POS`, `POSITIVE`, `+`, `POSITIVE ION`, ...). Anything non-empty that cannot
#' be classified is an **error**: silently defaulting an unrecognised mode to
#' `"NEG"` turns a positive-mode library into a negative-mode one without any
#' signal.
#'
#' @param pol Vector to normalise.
#' @param arg Name used in error messages.
#' @param default Value for `NA`/empty entries, or `NA_character_` to leave them
#'   `NA` for the caller to handle.
#'
#' @keywords internal
#' @noRd
.nts_normalize_polarity <- function(pol, arg = "Polarity",
                                    default = NA_character_) {
  pa <- as.character(pol)
  raw <- trimws(pa)
  empty <- is.na(pa) | !nzchar(raw)
  x <- toupper(raw)
  word_pos <- grepl("POS", x, fixed = TRUE)
  word_neg <- grepl("NEG", x, fixed = TRUE)
  neither <- !word_pos & !word_neg
  out <- rep(NA_character_, length(x))
  out[word_pos & !word_neg] <- "POS"
  out[word_neg & !word_pos] <- "NEG"
  # Bare signs count only when no word matched, so "POS-ESI" stays POS.
  out[neither & grepl("+", x, fixed = TRUE) & !grepl("-", x, fixed = TRUE)] <- "POS"
  out[neither & grepl("-", x, fixed = TRUE) & !grepl("+", x, fixed = TRUE)] <- "NEG"
  bad <- is.na(out) & !empty
  if (any(bad)) {
    stop(
      "`", arg, "` must describe positive or negative ion mode; ",
      "unrecognised value(s): ",
      paste(utils::head(unique(raw[bad]), 5L), collapse = ", "),
      call. = FALSE
    )
  }
  if (any(empty) && !is.na(default)) {
    out[empty] <- as.character(default)[1L]
  }
  out
}

#' @keywords internal
#' @noRd
.nts_normalize_adduct_token <- function(x) {
  s <- trimws(tolower(as.character(x)))
  s <- gsub("\\s+", "", s, perl = TRUE)
  s <- gsub("\\[|\\]", "", s, perl = TRUE)
  s
}

#' Neutral monoisotopic mass from observed precursor m/z: `mass = prec_mz + delta`.
#'
#' @keywords internal
#' @noRd
.nts_precursor_mz_to_neutral_delta <- function(adduct_token, polarity_fallback) {
  # Deltas are expressed against the *adduct ion* m/z. A singly charged ion
  # carries one electron fewer than the neutral atoms it is built from (one
  # more, for an anion), so the tabulated ATOM masses below need the electron
  # term:
  #   [M+Na]+  : m/z = M + m(Na)  - m(e)
  #   [M+K]+   : m/z = M + m(K)   - m(e)
  #   [M+NH4]+ : m/z = M + m(NH4) - m(e)
  #   [M+Cl]-  : m/z = M + m(Cl)  + m(e)
  # The proton is the exception: the *proton mass* already equals m(1H) - m(e),
  # so [M+H]+ / [M-H]- use it as-is and must NOT subtract a second electron.
  # Omitting the electron term for the atom adducts biases the inferred neutral
  # mass by ~0.55 mDa.
  e_m <- 0.00054858
  proton <- 1.00727646                # m(1H) - m(e)
  na_m <- 22.98976928 - e_m
  k_m <- 38.96370668 - e_m
  cl_m <- 34.96885268 + e_m
  nh4_m <- 18.034374 - e_m            # 14.003074 + 4 * 1.007825 = 18.034374(neutral)
  h2o_m <- 18.01056468346
  a <- .nts_normalize_adduct_token(adduct_token)
  if (is.na(a) || !nzchar(a)) {
    pol <- toupper(as.character(polarity_fallback)[1L])
    if (!pol %in% c("POS", "NEG")) {
      pol <- "NEG"
    }
    return(if (pol == "NEG") proton else -proton)
  }
  if (a %in% c("m-h", "m-h-", "[m-h]-", "negative", "neg")) {
    return(proton)
  }
  if (a %in% c("m+h", "m+h+", "[m+h]+", "positive", "pos")) {
    return(-proton)
  }
  if (a %in% c("m+na", "m+na+", "[m+na]+", "na+")) {
    return(-na_m)
  }
  if (a %in% c("m+k", "m+k+", "[m+k]+", "k+")) {
    return(-k_m)
  }
  if (a %in% c("m+nh4", "m+nh4+", "[m+nh4]+", "nh4+")) {
    return(-nh4_m)
  }
  if (a %in% c("m+cl", "m+cl-", "[m+cl]-", "cl-")) {
    return(-cl_m)
  }
  if (a %in% c("m-h-h2o", "m-h-h2o-", "[m-h-h2o]-", "m-h2o-h")) {
    return(proton + h2o_m)
  }
  if (a %in% c("m+h-h2o", "m+h-h2o+", "[m+h-h2o]+")) {
    return(-proton + h2o_m)
  }
  pol <- toupper(as.character(polarity_fallback)[1L])
  if (!pol %in% c("POS", "NEG")) {
    pol <- "NEG"
  }
  if (pol == "NEG") proton else -proton
}

#' @keywords internal
#' @noRd
.nts_default_features_long_mz_digits <- function() {
  c(DF = 5L, NL = 5L, MD = 3L)
}

#' @keywords internal
#' @noRd
.nts_resolve_features_long_mz_digits <- function(x) {
  def <- .nts_default_features_long_mz_digits()
  if (is.null(x)) {
    return(def)
  }
  if (length(x) == 1L && is.null(names(x))) {
    d <- suppressWarnings(as.integer(x)[1L])
    if (!is.finite(d) || d < 0L) {
      d <- 5L
    }
    return(c(DF = d, NL = d, MD = d))
  }
  xu <- unlist(x, recursive = FALSE, use.names = TRUE)
  nm <- names(xu)
  if (is.null(nm) || !any(nzchar(nm))) {
    stop(
      "`features_long_mz_digits` must be a single integer or a named vector ",
      "with keys among DF, NL, MD.",
      call. = FALSE
    )
  }
  out <- def
  for (i in seq_along(nm)) {
    key <- toupper(trimws(nm[i]))
    if (!nzchar(key) || !key %in% names(out)) {
      next
    }
    di <- suppressWarnings(as.integer(xu[[i]])[1L])
    if (is.finite(di) && di >= 0L) {
      out[key] <- di
    }
  }
  out
}

#' @keywords internal
#' @noRd
.nts_apply_features_long_mz_digits <- function(features_tb, digit_map) {
  if (!nrow(features_tb) || !"Feature_Type" %in% names(features_tb)) {
    return(features_tb)
  }
  ftu <- toupper(as.character(features_tb$Feature_Type))
  nm <- toupper(names(digit_map))
  d <- digit_map[match(ftu, nm)]
  d[is.na(d)] <- digit_map[["DF"]]
  features_tb[, mz_value := round(as.numeric(mz_value), d)]
  features_tb
}

#' @keywords internal
#' @noRd
.nts_mz_assign_bin <- function(mz, bin_width = NULL, digits = NULL) {
  mz <- as.numeric(mz)
  bw <- suppressWarnings(as.numeric(bin_width)[1L])
  if (is.finite(bw) && bw > 0) {
    return(round(mz / bw) * bw)
  }
  d <- suppressWarnings(as.integer(digits)[1L])
  if (!is.finite(d) || d < 0L) {
    d <- 4L
  }
  round(mz, d)
}

#' Build [NTSDatabase] from Spectraverse-style CSV metadata + MGF spectra
#'
#' Reads compound metadata with [data.table::fread()], peak lists with
#' [import_mgf()], harmonises `Lab_ID` / `mass` / `Polarity`, then calls
#' [build_nts_database()]. CSV identifiers are matched to MGF `TITLE=` values.
#'
#' @param csv_file Path to metadata `.csv` (or other `fread`-readable table).
#' @param mgf_file Path to `.mgf` from the same export.
#' @param id_col Optional column name for spectrum ID (matched to MGF `TITLE`).
#'   If `NULL`, the first match among common names is used.
#' @param mass_col Optional column for neutral mass. If `NULL`, auto-detect; if
#'   only a precursor m/z column exists, neutral mass and `target_mz` are
#'   derived using per-row polarity.
#' @param polarity_col Optional polarity column; if `NULL`, auto-detect; if
#'   absent, `default_polarity` is used for all rows.
#' @param default_polarity `"NEG"` or `"POS"` when the table has no polarity
#'   column.
#' @param adduct_col Optional column name for adduct labels when inferring
#'   neutral `mass` from precursor m/z (e.g. `"[M-H]-"`, `"M+Na"`). Unrecognised
#'   labels fall back to `[M+H]+` / `[M-H]-` from the row `Polarity`.
#' @param noise_threshold,top_k,version_str Passed to [build_nts_database()].
#' @param spectra_mz_bin_digits,features_long_mz_digits Passed to
#'   [build_nts_database()] (consensus bin width and per-feature-type rounding
#'   for `features_long$mz_value`).
#'
#' @return An [NTSDatabase] S4 object.
#'
#' @export
build_spectraverse_db <- function(
    csv_file,
    mgf_file,
    id_col = NULL,
    mass_col = NULL,
    polarity_col = NULL,
    default_polarity = "NEG",
    adduct_col = NULL,
    noise_threshold = 0.01,
    top_k = 30L,
    version_str = "1.0.0",
    spectra_mz_bin_digits = 4L,
    features_long_mz_digits = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }

  dp <- toupper(as.character(default_polarity)[1L])
  if (!dp %in% c("NEG", "POS")) {
    stop("`default_polarity` must be \"NEG\" or \"POS\".", call. = FALSE)
  }

  csv_file <- path.expand(as.character(csv_file)[1L])
  .nts_check(csv_file, "file_exists", .arg = "csv_file")

  md <- data.table::fread(
    csv_file,
    data.table = TRUE,
    na.strings = c("", "NA", "N/A", "null")
  )
  if (!nrow(md)) {
    stop("`csv_file` has no rows.", call. = FALSE)
  }

  nm <- names(md)
  id_src <- id_col
  if (is.null(id_src) || !nzchar(as.character(id_src)[1L])) {
    id_src <- .nts_resolve_csv_col(
      nm,
      c(
        "Lab_ID",
        "lab_id",
        "spectrum_id",
        "Spectrum_ID",
        "SPECTRUM_ID",
        "spectrum_uuid",
        "ID",
        "id",
        "TITLE",
        "title",
        "ACCESSION",
        "accession",
        "name",
        "Name"
      )
    )
  } else {
    id_src <- as.character(id_col)[1L]
    if (!id_src %in% nm) {
      stop("`id_col` not found in CSV: ", id_src, call. = FALSE)
    }
  }
  if (is.null(id_src)) {
    stop(
      "Could not find an ID column; set `id_col` (e.g. spectrum id / TITLE).",
      call. = FALSE
    )
  }
  md[, Lab_ID := as.character(trimws(as.character(.SD[[1L]]))), .SDcols = id_src]

  pol_src <- polarity_col
  if (is.null(pol_src) || !nzchar(as.character(pol_src)[1L])) {
    pol_src <- .nts_resolve_csv_col(
      nm,
      c(
        "Polarity",
        "polarity",
        "ion_mode",
        "Ion_mode",
        "ION_MODE",
        "Ionization",
        "ionization",
        "Mode",
        "mode"
      )
    )
  } else {
    pol_src <- as.character(polarity_col)[1L]
    if (!pol_src %in% nm) {
      stop("`polarity_col` not found in CSV: ", pol_src, call. = FALSE)
    }
  }
  if (is.null(pol_src)) {
    md[, Polarity := dp]
  } else {
    md[, Polarity := .nts_spectraverse_normalize_polarity_vec(as.character(.SD[[1L]])), .SDcols = pol_src]
  }

  mass_src <- mass_col
  if (is.null(mass_src) || !nzchar(as.character(mass_src)[1L])) {
    mass_src <- .nts_resolve_csv_col(
      nm,
      c(
        "mass",
        "neutral_mass",
        "NeutralMass",
        "MonoisotopicMass",
        "exact_mass",
        "ExactMass",
        "EXACT_MASS",
        "neutral_mass_gnps"
      )
    )
  } else {
    mass_src <- as.character(mass_col)[1L]
    if (!mass_src %in% nm) {
      stop("`mass_col` not found in CSV: ", mass_src, call. = FALSE)
    }
  }

  prec_src <- .nts_resolve_csv_col(
    nm,
    c(
      "PrecursorMZ",
      "precursor_mz",
      "PEPMASS",
      "pepmass",
      "Precursor_mz",
      "precursorMz",
      "PRECURSORMZ"
    )
  )

  proton <- 1.007276
  ad_src <- adduct_col
  if (!is.null(ad_src) && nzchar(as.character(ad_src)[1L])) {
    ad_src <- as.character(ad_src)[1L]
    if (!ad_src %in% names(md)) {
      stop("`adduct_col` not found in CSV: ", ad_src, call. = FALSE)
    }
  } else {
    ad_src <- NULL
  }

  if (!is.null(mass_src)) {
    md[, mass := suppressWarnings(as.numeric(.SD[[1L]])), .SDcols = mass_src]
    if (!is.null(prec_src)) {
      md[, target_mz := suppressWarnings(as.numeric(.SD[[1L]])), .SDcols = prec_src]
    }
  } else if (!is.null(prec_src)) {
    pmz <- suppressWarnings(as.numeric(md[[prec_src]]))
    polv <- md[["Polarity"]]
    if (!is.null(ad_src)) {
      ad_vec <- md[[ad_src]]
      delta <- mapply(
        function(a, pol) .nts_precursor_mz_to_neutral_delta(a, pol),
        ad_vec,
        polv,
        SIMPLIFY = TRUE,
        USE.NAMES = FALSE
      )
      md[, mass := pmz + as.numeric(delta)]
    } else {
      md[, mass := ifelse(
        polv == "NEG",
        pmz + proton,
        pmz - proton
      )]
    }
    md[, target_mz := pmz]
  } else {
    stop(
      "CSV must contain a neutral mass column or a precursor m/z column ",
      "(e.g. mass, PrecursorMZ, PEPMASS); set `mass_col` if non-standard.",
      call. = FALSE
    )
  }

  sl <- import_mgf(mgf_file)
  if (!nrow(sl)) {
    warning("import_mgf returned no peaks; database may be empty.", call. = FALSE)
  } else {
    sl[, intensity := (as.numeric(intensity) / max(as.numeric(intensity), na.rm = TRUE)) * 1000, by = .(Lab_ID)]
    sl <- sl[intensity >= 10]
  }

  n_md <- nrow(md)
  n_sl_lab <- if (nrow(sl)) data.table::uniqueN(sl$Lab_ID) else 0L
  n_both <- if (nrow(sl)) {
    data.table::uniqueN(sl[Lab_ID %in% md$Lab_ID, Lab_ID])
  } else {
    0L
  }
  message(
    "build_spectraverse_db: CSV rows = ",
    n_md,
    "; MGF spectra (unique Lab_ID) = ",
    n_sl_lab,
    "; overlapping Lab_ID = ",
    n_both
  )

  db <- build_nts_database(
    md,
    sl,
    noise_threshold = noise_threshold,
    top_k = top_k,
    version_str = version_str,
    spectra_mz_bin_digits = spectra_mz_bin_digits,
    features_long_mz_digits = features_long_mz_digits,
    db_type = "online"
  )
  
  meta_dt <- data.table::as.data.table(methods::slot(db, "metadata"))
  spec_dt <- data.table::as.data.table(methods::slot(db, "spectra_long"))
  feat_dt <- data.table::as.data.table(methods::slot(db, "features_long"))
  
  has_inchikey <- "INCHIKEY" %in% names(meta_dt)
  has_smiles <- "SMILES" %in% names(meta_dt)
  
  if (has_inchikey || has_smiles) {
    message("build_spectraverse_db: Applying deduplication and polarity mapping...")
    
    # 1. Sync polarity from metadata to spectra and features
    polarity_map <- meta_dt[, .(Lab_ID, correct_Polarity = Polarity)]
    
    spec_dt <- merge(spec_dt, polarity_map, by = "Lab_ID", all.x = TRUE)
    spec_dt[, Polarity := correct_Polarity]
    spec_dt[, correct_Polarity := NULL]
    
    feat_dt <- merge(feat_dt, polarity_map, by = "Lab_ID", all.x = TRUE)
    feat_dt[, Polarity := correct_Polarity]
    feat_dt[, correct_Polarity := NULL]
    
    # 2. Group by compound
    compound_keys <- character(nrow(meta_dt))
    if (has_inchikey && has_smiles) {
      compound_keys <- ifelse(!is.na(meta_dt$INCHIKEY) & nzchar(as.character(meta_dt$INCHIKEY)), 
                              as.character(meta_dt$INCHIKEY), 
                              as.character(meta_dt$SMILES))
    } else if (has_inchikey) {
      compound_keys <- as.character(meta_dt$INCHIKEY)
    } else {
      compound_keys <- as.character(meta_dt$SMILES)
    }
    
    meta_dt[, compound_key := compound_keys]
    
    id_mapping <- meta_dt[, .(
      new_Lab_ID = paste0("NTS_C", sprintf("%06d", .GRP))
    ), by = .(compound_key)]
    
    meta_dt <- merge(meta_dt, id_mapping, by = "compound_key", all.x = TRUE)
    id_lookup <- unique(meta_dt[, .(Lab_ID, new_Lab_ID)])
    
    spec_dt <- merge(spec_dt, id_lookup, by = "Lab_ID", all.x = TRUE)
    spec_dt[, Lab_ID := new_Lab_ID]
    spec_dt[, new_Lab_ID := NULL]
    
    feat_dt <- merge(feat_dt, id_lookup, by = "Lab_ID", all.x = TRUE)
    feat_dt[, Lab_ID := new_Lab_ID]
    feat_dt[, new_Lab_ID := NULL]
    
    # 3. Aggregate metadata
    get_first_valid <- function(x) {
      valid_x <- stats::na.omit(x[nzchar(as.character(x))])
      if (length(valid_x) > 0) return(as.character(valid_x)[1L]) else return(NA_character_)
    }
    get_all_unique <- function(x) {
      valid_x <- unique(stats::na.omit(x[nzchar(as.character(x))]))
      if (length(valid_x) > 0) return(paste(valid_x, collapse = " | ")) else return(NA_character_)
    }
    
    for (cname in c("CHARGE", "ADDUCT", "INSTRUMENT_TYPE", "FORMULA", "INCHIKEY", "SMILES")) {
      if (!cname %in% names(meta_dt)) meta_dt[, (cname) := NA_character_]
    }
    
    clean_meta <- meta_dt[, .(
      mass = round(mean(as.numeric(mass), na.rm = TRUE), 4),
      SMILES = get_first_valid(SMILES),
      InChIKey = get_first_valid(INCHIKEY),
      Formula = get_first_valid(FORMULA),
      Charge = get_first_valid(CHARGE),
      Polarity = paste(sort(unique(as.character(Polarity))), collapse = "/"),
      Adduct = get_first_valid(ADDUCT),
      Instrument_type = get_all_unique(INSTRUMENT_TYPE)
    ), by = .(Lab_ID = new_Lab_ID)]
    
    # Resolve conflicting polarities using Charge
    clean_meta[Charge == 1 | Charge == "1", Polarity := "POS"]
    clean_meta[Charge == -1 | Charge == "-1", Polarity := "NEG"]
    
    final_polarity_map <- clean_meta[, .(Lab_ID, final_Pol = Polarity)]
    
    spec_dt <- merge(spec_dt, final_polarity_map, by = "Lab_ID", all.x = TRUE)
    spec_dt[, Polarity := final_Pol]
    spec_dt[, final_Pol := NULL]
    
    feat_dt <- merge(feat_dt, final_polarity_map, by = "Lab_ID", all.x = TRUE)
    feat_dt[, Polarity := final_Pol]
    feat_dt[, final_Pol := NULL]
    
    # No source RT is available for this ETL path: store a true missing value
    # (`NA_real_`) rather than a constant `1`, so downstream `is.na()` checks on
    # the metadata RT behave correctly (a constant would sit outside every RT
    # window and silently disable RT-based filtering/sorting).
    clean_meta[, RT := NA_real_]

    std_cols <- intersect(c("Lab_ID", "mass", "RT"), names(clean_meta))
    other_cols <- setdiff(names(clean_meta), std_cols)
    data.table::setcolorder(clean_meta, c(std_cols, other_cols))
    
    methods::slot(db, "metadata") <- clean_meta
    methods::slot(db, "spectra_long") <- spec_dt
    methods::slot(db, "features_long") <- feat_dt
  }
  
  db
}

#' @keywords internal
#' @noRd
.nts_build_nts_database_precursor_mz <- function(md) {
  proton <- 1.007276
  pol <- toupper(as.character(md[["Polarity"]]))
  pol[is.na(pol) | pol == ""] <- "NEG"
  pol[!pol %in% c("POS", "NEG")] <- "NEG"
  mass <- as.numeric(md[["mass"]])
  ifelse(pol == "NEG", mass - proton, mass + proton)
}

#' Screening precursor m/z from metadata (`target_mz` > polarity mz > mass).
#'
#' Aligns with [extract_mix_standards()] and in-house build scripts: when
#' `target_mz` is set (ISF / adduct / override ions), it wins; otherwise use
#' `mz_neg` / `mz_pos` (or `mz.neg` / `mz.pos`) for the active ion mode.
#'
#' @param metadata_tb `data.frame` or `data.table` with `Lab_ID` optional.
#' @param polarity `"neg"` or `"pos"` for mode-specific fallbacks.
#' @return Numeric vector length `nrow(metadata_tb)` (may contain `NA`).
#' @keywords internal
#' @noRd
.nts_metadata_screening_mz <- function(metadata_tb, polarity = c("neg", "pos")) {
  polarity <- match.arg(polarity)
  meta <- data.table::as.data.table(metadata_tb)
  n <- nrow(meta)
  if (!n) {
    return(numeric(0))
  }
  proton <- 1.007276
  pol_neg <- identical(polarity, "neg")
  out <- rep(NA_real_, n)

  if ("target_mz" %in% names(meta)) {
    tm <- suppressWarnings(as.numeric(meta$target_mz))
    ok <- is.finite(tm) & tm > 0
    out[ok] <- tm[ok]
  }

  miss <- !is.finite(out) | out <= 0
  if (any(miss)) {
    pol_cols <- if (pol_neg) {
      c("mz_neg", "mz.neg")
    } else {
      c("mz_pos", "mz.pos")
    }
    for (col in pol_cols) {
      if (!col %in% names(meta)) {
        next
      }
      v <- suppressWarnings(as.numeric(meta[[col]]))
      ok <- miss & is.finite(v) & v > 0
      out[ok] <- v[ok]
      miss <- !is.finite(out) | out <= 0
    }
    if (any(miss) && all(c("mz_pos", "mz_neg") %in% names(meta))) {
      v <- if (pol_neg) {
        suppressWarnings(as.numeric(meta$mz_neg))
      } else {
        suppressWarnings(as.numeric(meta$mz_pos))
      }
      ok <- miss & is.finite(v) & v > 0
      out[ok] <- v[ok]
      miss <- !is.finite(out) | out <= 0
    }
  }

  if (any(miss) && "Precursor_mz" %in% names(meta)) {
    v <- suppressWarnings(as.numeric(meta$Precursor_mz))
    ok <- miss & is.finite(v) & v > 0
    out[ok] <- v[ok]
    miss <- !is.finite(out) | out <= 0
  }
  if (any(miss) && "mz" %in% names(meta)) {
    v <- suppressWarnings(as.numeric(meta$mz))
    ok <- miss & is.finite(v) & v > 0
    out[ok] <- v[ok]
    miss <- !is.finite(out) | out <= 0
  }
  if (any(miss) && "mass" %in% names(meta)) {
    m <- suppressWarnings(as.numeric(meta$mass))
    v <- if (pol_neg) m - proton else m + proton
    ok <- miss & is.finite(v) & v > 0
    out[ok] <- v[ok]
  }
  out
}

#' Build a unified [NTSDatabase] from metadata and long MS2 peaks
#'
#' Single entry point for in-house and external-derived tables: purifies
#' fragment lists (base-peak normalization, noise floor, top-\code{k}),
#' derives simple neutral-loss and diagnostic-fragment features, appends an
#' `NCE == \"Merged\"` consensus layer compatible with
#' [compile_reference_database()], and returns an S4 object.
#'
#' @param metadata `data.table` or `data.frame` with at least `Lab_ID` and
#'   `mass`. `Polarity` is normalised to `"POS"` / `"NEG"` (accepting
#'   `POSITIVE` / `NEGATIVE` / `+` / `-`); empty values default to `"NEG"` and
#'   an unrecognised non-empty value is an **error**. If the column is missing
#'   entirely, a warning is emitted and every row is assumed `"NEG"` -- pass a
#'   `Polarity` column (e.g. from [extract_mix_standards()]) when building a
#'   positive-mode library. Optional `target_mz` (e.g. from
#'   [extract_mix_standards()]) overrides theoretical precursor m/z for
#'   neutral-loss features.
#' @param spectra_long `data.table` with `Lab_ID`, `mz`, `intensity`.
#' @param noise_threshold Relative intensity threshold after per-`Lab_ID`
#'   base-peak normalization (default `0.01` = 1%).
#' @param top_k Maximum number of peaks to retain per `Lab_ID` after noise
#'   filtering (default `30`).
#' @param version_str Schema version stored in the S4 object (default
#'   `"1.0.0"`).
#' @param spectra_mz_bin_digits Decimal places for m/z binning (`round(m/z, d)`)
#'   when `spectra_mz_bin_width` is not set (default `4`).
#' @param spectra_mz_bin_width Optional fixed m/z bin width in Da (e.g. `0.005`).
#'   When set, takes precedence over `spectra_mz_bin_digits`.
#' @param add_merged_layer If `TRUE` (default), append an `NCE == \"Merged\"`
#'   consensus layer across collision energies. If `FALSE`, only per-NCE spectra
#'   are kept.
#' @param features_mz_bin_width Optional fixed bin width (Da) for
#'   `features_long$mz_value` after digit rounding; e.g. `0.005`.
#' @param features_long_mz_digits Rounding for `features_long$mz_value` by
#'   `Feature_Type`: default `DF = 5`, `NL = 5`, `MD = 3`. Pass a single
#'   non-negative integer to use the same width for all types, or a named
#'   vector/list (e.g. `c(DF = 5, NL = 4)`).
#' @param predict_rt If `TRUE`, predict a retention time per `Lab_ID` with
#'   `rt_model` (see [build_rt_model]) and attach it to the result. Default
#'   `FALSE`.
#' @param rt_model An `nts_rt_model`; required when `predict_rt = TRUE`.
#' @param db_type Database provenance stamped onto every compound as a
#'   `db_type` column: one of `"inhouse"`, `"online"` or `"prediction"`. An
#'   explicit `metadata$db_type` column takes precedence. This value is what
#'   [assign_confidence_levels] uses to select the confidence tiers
#'   (`inhouse` -> L1, `online` -> L2a, `prediction` -> L2b). Leave
#'   `NA` (the default) when the provenance is unknown; the column is then
#'   omitted and [assign_confidence_levels] falls back to `Lab_ID` pattern
#'   detection.
#'
#' @return An [NTSDatabase] instance.
#'
#' @export
build_nts_database <- function(
    metadata,
    spectra_long,
    noise_threshold = 0.01,
    top_k = 30L,
    version_str = "1.0.0",
    spectra_mz_bin_digits = 4L,
    spectra_mz_bin_width = NULL,
    add_merged_layer = TRUE,
    features_long_mz_digits = NULL,
    features_mz_bin_width = NULL,
    predict_rt = FALSE,
    rt_model = NULL,
    db_type = NA_character_) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }

  md <- data.table::as.data.table(metadata)
  md <- .nts_dot_to_underscore_names(md)
  req_m <- c("Lab_ID", "mass")
  miss_m <- setdiff(req_m, names(md))
  if (length(miss_m)) {
    stop(
      "`metadata` must contain columns: ",
      paste(req_m, collapse = ", "),
      call. = FALSE
    )
  }

  if (!"SMILES" %in% names(md)) {
    warning(
      "\u26a0\ufe0f Metadata \u7f3a\u5c11 SMILES \u5217\uff0c\u5c06\u65e0\u6cd5\u53c2\u4e0e\u53cc\u7f51\u62d3\u6251\u63a8\u6f14\uff01",
      call. = FALSE
    )
    md[, SMILES := NA_character_]
  }

  if (!"Polarity" %in% names(md)) {
    warning(
      "`metadata` has no `Polarity` column; assuming \"NEG\" for every entry. ",
      "Add a `Polarity` column (e.g. from extract_mix_standards()) so that a ",
      "positive-mode library is not mislabelled as negative.",
      call. = FALSE
    )
    md[, Polarity := "NEG"]
  } else {
    md[, Polarity := .nts_normalize_polarity(Polarity,
      arg = "metadata$Polarity", default = "NEG")]
  }

  # Database provenance ("inhouse" / "online" / "prediction") drives the
  # confidence-level tiers (L1 = in-house, L2a = online, L2b =
  # prediction). An explicit `metadata$db_type` column wins; otherwise the
  # caller's `db_type` argument is stamped onto every compound. When neither is
  # available the column is left absent, so that downstream confidence
  # assignment can fall back to ID-pattern detection and warn (rather than
  # silently treating an in-house library as "online").
  dbt <- tolower(as.character(db_type)[1L])
  if (!dbt %in% c("inhouse", "online", "prediction")) dbt <- NA_character_
  if ("db_type" %in% names(md)) {
    md[, db_type := tolower(as.character(db_type))]
    md[!db_type %in% c("inhouse", "online", "prediction"), db_type := NA_character_]
  } else if (!is.na(dbt)) {
    md[, db_type := dbt]
  }

  md[, Lab_ID := as.character(Lab_ID)]
  md[, mass := as.numeric(mass)]
  if (any(!is.finite(md$mass))) {
    stop("`metadata$mass` must be finite numeric.", call. = FALSE)
  }

  if ("target_mz" %in% names(md)) {
    md[, prec_mz_work := as.numeric(target_mz)]
  } else {
    md[, prec_mz_work := .nts_build_nts_database_precursor_mz(md)]
  }
  md[is.na(prec_mz_work) | prec_mz_work <= 0, prec_mz_work := NA_real_]

  sl <- data.table::as.data.table(spectra_long)
  need_s <- c("Lab_ID", "mz", "intensity")
  miss_s <- setdiff(need_s, names(sl))
  if (length(miss_s)) {
    stop(
      "`spectra_long` must contain columns: ",
      paste(need_s, collapse = ", "),
      call. = FALSE
    )
  }
  sl[, Lab_ID := as.character(Lab_ID)]
  sl[, mz := as.numeric(mz)]
  sl[, intensity := as.numeric(intensity)]
  sl <- sl[is.finite(mz) & is.finite(intensity)]

  nth <- suppressWarnings(as.integer(top_k)[1L])
  if (!is.finite(nth) || nth < 1L) {
    nth <- 30L
  }
  nt <- suppressWarnings(as.numeric(noise_threshold)[1L])
  if (!is.finite(nt) || nt < 0) {
    nt <- 0.01
  }

  bin_d <- suppressWarnings(as.integer(spectra_mz_bin_digits)[1L])
  if (!is.finite(bin_d) || bin_d < 0L) {
    bin_d <- 4L
  }
  spec_bin_w <- suppressWarnings(as.numeric(spectra_mz_bin_width)[1L])
  if (!is.finite(spec_bin_w) || spec_bin_w <= 0) {
    spec_bin_w <- NA_real_
  }
  feat_bin_w <- suppressWarnings(as.numeric(features_mz_bin_width)[1L])
  if (!is.finite(feat_bin_w) || feat_bin_w <= 0) {
    feat_bin_w <- NA_real_
  }
  add_merged <- isTRUE(add_merged_layer)
  feat_dmap <- .nts_resolve_features_long_mz_digits(features_long_mz_digits)

  pol_join <- unique(md[, .(Lab_ID, Polarity)])
  sl <- sl[Lab_ID %in% pol_join$Lab_ID]
  sl <- pol_join[sl, on = "Lab_ID"]

  # Resolve the collision-energy label *before* grouping. Baseline
  # normalization and the top-k peak selection below must be per
  # (Lab_ID, Polarity, NCE) -- the same groups optimize_spectra() uses -- so
  # that peaks measured at one CE are not renormalized/truncated by the
  # intensities of another CE.
  if (!"NCE" %in% names(sl)) {
    sl[, NCE := "InHouse"]
  } else {
    sl[, NCE := as.character(NCE)]
    sl[is.na(NCE) | !nzchar(NCE), NCE := "InHouse"]
  }

  grp_cols <- c("Lab_ID", "Polarity", "NCE")
  sl[, imax := max(intensity, na.rm = TRUE), by = grp_cols]
  sl <- sl[imax > 0]
  sl[, int_rel := intensity / imax]
  sl <- sl[int_rel >= nt]
  sl[, rk := data.table::frank(-intensity, ties.method = "first"), by = grp_cols]
  sl <- sl[rk <= nth]
  sl[, c("imax", "int_rel", "rk") := NULL]

  if (!add_merged) {
    sl <- sl[NCE != "Merged"]
  }

  spectra_tb <- sl[, .(
    Lab_ID = Lab_ID,
    Polarity = Polarity,
    NCE = NCE,
    frag_mz = mz,
    intensity = intensity
  )]

  if (nrow(spectra_tb)) {
    if (add_merged) {
      merged_src <- spectra_tb[NCE != "Merged"]
      merged_dt <- merged_src[
        ,
        .(
          NCE = "Merged",
          frag_mz = if (all(is.na(frag_mz))) {
            NA_real_
          } else {
            mean(frag_mz, na.rm = TRUE)
          },
          intensity = if (all(is.na(intensity))) {
            NA_real_
          } else {
            max(intensity, na.rm = TRUE)
          }
        ),
        by = .(
          Lab_ID,
          Polarity,
          mz_bin = .nts_mz_assign_bin(frag_mz, bin_width = spec_bin_w, digits = bin_d)
        )
      ]
      merged_dt[, mz_bin := NULL]
      data.table::setcolorder(merged_dt, names(spectra_tb))
      spectra_tb <- data.table::rbindlist(
        list(spectra_tb, merged_dt),
        use.names = TRUE
      )
    }
    data.table::setkeyv(spectra_tb, c("Lab_ID", "Polarity", "NCE"))
  } else {
    spectra_tb <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      NCE = character(),
      frag_mz = numeric(),
      intensity = numeric()
    )
  }

  prec_map <- unique(md[, .(Lab_ID, prec_mz_work)])
  sp_inhouse <- spectra_tb[NCE != "Merged"]
  sp_inhouse <- prec_map[sp_inhouse, on = "Lab_ID"]

  nl_ok <- sp_inhouse[
    is.finite(prec_mz_work) & is.finite(frag_mz) & frag_mz < prec_mz_work
  ]
  nl_ok[, int_max := max(intensity, na.rm = TRUE), by = Lab_ID]
  nl_ok[int_max > 0, Frequency := intensity / int_max]
  nl_ok[int_max <= 0 | is.na(int_max), Frequency := NA_real_]
  nl_dt <- nl_ok[
    ,
    .(
      Lab_ID = Lab_ID,
      Polarity = Polarity,
      Feature_Type = "NL",
      mz_value = prec_mz_work - frag_mz,
      Frequency = Frequency
    )
  ]

  df_ok <- copy(sp_inhouse)
  df_ok[, int_max := max(intensity, na.rm = TRUE), by = Lab_ID]
  df_ok[int_max > 0, Frequency := intensity / int_max]
  df_ok[int_max <= 0 | is.na(int_max), Frequency := NA_real_]
  df_dt <- df_ok[
    ,
    .(
      Lab_ID = Lab_ID,
      Polarity = Polarity,
      Feature_Type = "DF",
      mz_value = frag_mz,
      Frequency = Frequency
    )
  ]

  md_ok <- copy(sp_inhouse)
  md_dt <- md_ok[, {
    n <- .N
    if (n > 1) {
      grid <- expand.grid(i = seq_len(n), j = seq_len(n))
      grid <- grid[grid$i < grid$j, ]
      diffs <- abs(frag_mz[grid$i] - frag_mz[grid$j])
      freqs <- sqrt(intensity[grid$i] * intensity[grid$j])
      list(mz_value = diffs, raw_freq = freqs)
    } else {
      list(mz_value = numeric(), raw_freq = numeric())
    }
  }, by = .(Lab_ID, Polarity, NCE)]

  if (nrow(md_dt)) {
    md_dt[, int_max := max(raw_freq, na.rm = TRUE), by = Lab_ID]
    md_dt[int_max > 0, Frequency := raw_freq / int_max]
    md_dt[int_max <= 0 | is.na(int_max), Frequency := NA_real_]
    md_dt <- md_dt[!is.na(Frequency)]
    
    # Keep top 50% per Lab_ID
    md_dt[, rk := data.table::frank(-Frequency, ties.method = "first"), by = Lab_ID]
    md_dt[, max_rk := .N, by = Lab_ID]
    md_dt <- md_dt[rk <= ceiling(max_rk * 0.5)]
    
    md_dt <- md_dt[
      ,
      .(
        Lab_ID = Lab_ID,
        Polarity = Polarity,
        Feature_Type = "MD",
        mz_value = mz_value,
        Frequency = Frequency
      )
    ]
  } else {
    md_dt <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      Feature_Type = character(),
      mz_value = numeric(),
      Frequency = numeric()
    )
  }

  features_tb <- data.table::rbindlist(list(nl_dt, df_dt, md_dt), use.names = TRUE)
  if (nrow(features_tb)) {
    features_tb[, Lab_ID := as.character(Lab_ID)]
    features_tb[, Polarity := as.character(Polarity)]
    features_tb[, Feature_Type := as.character(Feature_Type)]
    features_tb[, mz_value := as.numeric(mz_value)]
    features_tb[, Frequency := as.numeric(Frequency)]
    data.table::setkeyv(
      features_tb,
      c("Lab_ID", "Polarity", "Feature_Type")
    )
    features_tb <- .nts_apply_features_long_mz_digits(features_tb, feat_dmap)
    if (is.finite(feat_bin_w)) {
      features_tb[, mz_value := .nts_mz_assign_bin(mz_value, bin_width = feat_bin_w)]
    }
    features_tb <- features_tb[
      ,
      .(Frequency = max(Frequency, na.rm = TRUE)),
      by = .(Lab_ID, Polarity, Feature_Type, mz_value)
    ]
  } else {
    features_tb <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      Feature_Type = character(),
      mz_value = numeric(),
      Frequency = numeric()
    )
  }

  md[, prec_mz_work := NULL]

  ver <- as.character(version_str)[1L]
  if (!nzchar(ver)) {
    ver <- "1.0.0"
  }

  # Optional RT prediction
  if (isTRUE(predict_rt)) {
    if (!inherits(rt_model, "nts_rt_model")) {
      stop("predict_rt = TRUE requires a valid rt_model (nts_rt_model object).",
           call. = FALSE)
    }
    md <- .nts_add_predicted_rt(md, rt_model)
  }

  message(
    "build_nts_database: spectra_long rows = ",
    nrow(spectra_tb),
    "; features_long rows = ",
    nrow(features_tb)
  )

  methods::new(
    "NTSDatabase",
    metadata = md,
    spectra_long = spectra_tb,
    features_long = features_tb,
    version = ver
  )
}

#' Set the `db_type` provenance of an [NTSDatabase]
#'
#' Stamps a `db_type` value onto every row of `@metadata`
#' (`"inhouse"` / `"online"` / `"prediction"`). This drives the confidence
#' tiers assigned by [assign_confidence_levels]: `inhouse` -> L1,
#' `online` -> L2a, `prediction` -> L2b. Mainly useful for retro-fitting
#' provenance onto a library that was built before `db_type` existed (in which
#' case [assign_confidence_levels] would otherwise fall back to `Lab_ID`
#' pattern detection and could silently demote an in-house library to the
#' online tiers).
#'
#' @param db An [NTSDatabase].
#' @param db_type One of `"inhouse"`, `"online"` or `"prediction"`.
#'
#' @return `db` with an updated `db_type` column in `@metadata`.
#'
#' @seealso [build_nts_database], [assign_confidence_levels]
#' @export
set_db_type <- function(db, db_type) {
  if (!.is_reference_spectral_db(db)) {
    stop("`db` must be an NTSDatabase object", call. = FALSE)
  }
  dbt <- tolower(as.character(db_type)[1L])
  if (!dbt %in% c("inhouse", "online", "prediction")) {
    stop("`db_type` must be one of \"inhouse\", \"online\", \"prediction\".",
         call. = FALSE)
  }
  md <- data.table::copy(data.table::as.data.table(methods::slot(db, "metadata")))
  if (!nrow(md)) {
    stop("`db` has an empty metadata table; nothing to stamp.", call. = FALSE)
  }
  md[, db_type := dbt]
  methods::slot(db, "metadata") <- md
  db
}

#' Subset an [NTSDatabase] by compound identifiers
#'
#' Keeps rows in `metadata` whose `Lab_ID` is in `lab_ids`, and restricts
#' `spectra_long` / `features_long` to the same identifiers. Useful when
#' carving a target list out of a large Spectraverse or merged master library.
#'
#' @param db An [NTSDatabase].
#' @param lab_ids Character vector of `Lab_ID` values to retain.
#'
#' @return A new [NTSDatabase] (same `version` as `db`).
#'
#' @export
subset_nts_database <- function(db, lab_ids) {
  if (!.is_reference_spectral_db(db)) {
    stop("`db` must be an NTSDatabase object", call. = FALSE)
  }
  keep <- unique(as.character(lab_ids))
  if (!length(keep)) {
    stop("`lab_ids` must be a non-empty character vector", call. = FALSE)
  }
  md <- data.table::as.data.table(methods::slot(db, "metadata"))
  if (!nrow(md) || !"Lab_ID" %in% names(md)) {
    stop("db@metadata is empty or has no Lab_ID column", call. = FALSE)
  }
  md[, Lab_ID := as.character(Lab_ID)]
  md <- md[Lab_ID %in% keep]
  if (!nrow(md)) {
    warning("No metadata rows matched `lab_ids`; returning empty tables.", call. = FALSE)
  }
  lids <- unique(md$Lab_ID)
  sp <- data.table::as.data.table(methods::slot(db, "spectra_long"))
  if (nrow(sp) && "Lab_ID" %in% names(sp)) {
    sp <- sp[as.character(Lab_ID) %in% lids]
  }
  fl <- data.table::as.data.table(methods::slot(db, "features_long"))
  if (nrow(fl) && "Lab_ID" %in% names(fl)) {
    fl <- fl[as.character(Lab_ID) %in% lids]
  }
  methods::new(
    "NTSDatabase",
    metadata = md,
    spectra_long = sp,
    features_long = fl,
    version = methods::slot(db, "version")
  )
}

#' ---- merged from 01_pfas_database_compile.R ----
#' @keywords internal
#' @noRd
.nts_dot_to_underscore_names <- function(x) {
  nm <- names(x)
  if (length(nm)) {
    names(x) <- gsub(".", "_", nm, fixed = TRUE)
  }
  x
}

#' @keywords internal
#' @noRd
.nts_spectra_polarity_from_key <- function(pol_key) {
  k <- as.character(pol_key)[1L]
  if (!nzchar(k)) {
    return(NA_character_)
  }
  if (identical(k, "Spectra.positive")) {
    return("POS")
  }
  if (identical(k, "Spectra.negative")) {
    return("NEG")
  }
  if (toupper(k) %in% c("POS", "NEG")) {
    return(toupper(k))
  }
  NA_character_
}

#' @keywords internal
#' @noRd
.nts_nce_key_to_char <- function(nce_key) {
  x <- as.character(nce_key)[1L]
  if (!nzchar(x)) {
    return(x)
  }
  if (grepl("^NCE", x, ignore.case = FALSE)) {
    x <- sub("^NCE", "", x)
  }
  x
}

#' @keywords internal
#' @noRd
.nts_feature_type_map <- function(ft_key) {
  m <- c(
    neutral_loss = "NL",
    mass_difference = "MD",
    diagnostic_fragments = "DF"
  )
  k <- as.character(ft_key)[1L]
  hit <- unname(m[k])
  if (length(hit) == 1L && !is.na(hit)) {
    return(hit)
  }
  k
}

#' @keywords internal
#' @noRd
.nts_rename_first_present <- function(dt, from, to) {
  for (f in from) {
    if (f %in% names(dt)) {
      if (f != to) {
        data.table::setnames(dt, f, to, skip_absent = TRUE)
      }
      return(invisible(TRUE))
    }
  }
  invisible(FALSE)
}

#' Unify heterogeneous fragment feature column names (NL/MD/DF/mz -> mz_value).
#' @keywords internal
#' @noRd
.nts_fragment_feat_normalize_cols <- function(dt_feat) {
  dt_feat <- data.table::as.data.table(dt_feat)
  nm <- names(dt_feat)
  mz_col <- grep("^(NL|MD|DF|mz)$", nm, ignore.case = TRUE, value = TRUE)
  if (length(mz_col) > 0L) {
    data.table::setnames(dt_feat, mz_col[1L], "mz_value", skip_absent = TRUE)
  } else if (ncol(dt_feat) > 0L) {
    cn1 <- names(dt_feat)[1L]
    if (nzchar(cn1)) {
      data.table::setnames(dt_feat, cn1, "mz_value", skip_absent = TRUE)
    }
  }
  nm2 <- names(dt_feat)
  freq_col <- grep("^(Frequency|Freq|int)$", nm2, ignore.case = TRUE, value = TRUE)
  if (length(freq_col) > 0L) {
    data.table::setnames(dt_feat, freq_col[1L], "Frequency", skip_absent = TRUE)
  } else if (ncol(dt_feat) >= 2L) {
    cn2 <- names(dt_feat)[2L]
    if (nzchar(cn2)) {
      data.table::setnames(dt_feat, cn2, "Frequency", skip_absent = TRUE)
    }
  }
  dt_feat
}

#' Compile nested list reference data into a flat [NTSDatabase]
#'
#' ETL entry point: collects rows in lists and binds once with
#' [data.table::rbindlist] to avoid repeated allocation and copies in tight
#' loops. After flattening MS2 peaks, builds an NCE = `"Merged"` consensus layer
#' per `Lab_ID`, `Polarity`, and m/z bin (`round(frag_mz, spectra_mz_bin_digits)`).
#'
#' @param info_df Compound metadata (`data.frame` or `data.table`), e.g. legacy
#'   `spectra.info`. Column names get `.` replaced by `_`.
#' @param spectra_list Nested list: `Spectra.positive` / `Spectra.negative` (or
#'   `POS` / `NEG`) -> `Lab_ID` -> NCE name -> `data.frame` with m/z and
#'   intensity columns.
#' @param fragment_list Nested list: `POS` / `NEG` -> `Lab_ID` -> feature key
#'   (`neutral_loss`, `mass_difference`, `diagnostic_fragments`, ...) ->
#'   `data.frame` with m/z and frequency-like columns.
#' @param version_str Schema version string stored in the S4 object (length 1).
#' @param spectra_mz_bin_digits Passed through to consensus binning (same meaning
#'   as [build_nts_database()]).
#' @param features_long_mz_digits Rounding map for `features_long$mz_value`
#'   (same meaning as [build_nts_database()]).
#' @param predict_rt If `TRUE`, predict a retention time per `Lab_ID` with
#'   `rt_model` (see [build_rt_model]) and attach it to the result. Default
#'   `FALSE`.
#' @param rt_model An `nts_rt_model`; required when `predict_rt = TRUE`.
#'
#' @return An [NTSDatabase] instance.
#' @export
compile_reference_database <- function(info_df,
                                  spectra_list,
                                  fragment_list,
                                  version_str = "1.0.0",
                                  spectra_mz_bin_digits = 4L,
                                  features_long_mz_digits = NULL,
    predict_rt = FALSE,
    rt_model = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  if (!is.list(spectra_list)) {
    stop("spectra_list must be a list", call. = FALSE)
  }
  if (!is.list(fragment_list)) {
    stop("fragment_list must be a list", call. = FALSE)
  }
  ver <- as.character(version_str)[1L]
  if (!nzchar(ver)) {
    ver <- "1.0.0"
  }
  bin_d <- suppressWarnings(as.integer(spectra_mz_bin_digits)[1L])
  if (!is.finite(bin_d) || bin_d < 0L) {
    bin_d <- 4L
  }
  feat_dmap <- .nts_resolve_features_long_mz_digits(features_long_mz_digits)

  if (is.null(info_df)) {
    info_dt <- data.table::data.table()
  } else {
    info_dt <- data.table::as.data.table(info_df)
    info_dt <- .nts_dot_to_underscore_names(info_dt)
  }

  dt_list <- list()
  li <- 0L
  for (pol_key in names(spectra_list)) {
    pol_short <- .nts_spectra_polarity_from_key(pol_key)
    if (is.na(pol_short)) {
      next
    }
    by_lab <- spectra_list[[pol_key]]
    if (!is.list(by_lab)) {
      next
    }
    for (lab_id in names(by_lab)) {
      by_nce <- by_lab[[lab_id]]
      if (!is.list(by_nce)) {
        next
      }
      for (nce_name in names(by_nce)) {
        chunk <- by_nce[[nce_name]]
        if (is.null(chunk) || !is.data.frame(chunk)) {
          next
        }
        if (!nrow(chunk)) {
          next
        }
        tb <- data.table::as.data.table(chunk)
        tb[, Lab_ID := as.character(lab_id)]
        tb[, Polarity := pol_short]
        tb[, NCE := .nts_nce_key_to_char(nce_name)]
        li <- li + 1L
        dt_list[[li]] <- tb
      }
    }
  }

  if (!li) {
    spectra_tb <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      NCE = character(),
      frag_mz = numeric(),
      intensity = numeric()
    )
  } else {
    spectra_tb <- data.table::rbindlist(dt_list, fill = TRUE)
    .nts_rename_first_present(spectra_tb, c("mz", "MZ", "m.z", "mass"), "frag_mz")
    .nts_rename_first_present(
      spectra_tb,
      c("intensity", "Intensity", "int", "rel_intensity", "relative_intensity"),
      "intensity"
    )
    need <- c("Lab_ID", "Polarity", "NCE", "frag_mz", "intensity")
    miss <- setdiff(need, names(spectra_tb))
    if (length(miss)) {
      stop(
        "spectra_long is missing required column(s) after normalisation: ",
        paste(miss, collapse = ", "),
        call. = FALSE
      )
    }
    spectra_tb <- spectra_tb[, need, with = FALSE]
    spectra_tb[, Lab_ID := as.character(Lab_ID)]
    spectra_tb[, Polarity := as.character(Polarity)]
    spectra_tb[, NCE := as.character(NCE)]
    spectra_tb[, frag_mz := as.numeric(frag_mz)]
    spectra_tb[, intensity := as.numeric(intensity)]
  }

  data.table::setkeyv(spectra_tb, c("Lab_ID", "Polarity", "NCE"))

  if (nrow(spectra_tb)) {
    merged_src <- spectra_tb[NCE != "Merged"]
    merged_dt <- merged_src[
      ,
      .(
        NCE = "Merged",
        frag_mz = if (all(is.na(frag_mz))) {
          NA_real_
        } else {
          mean(frag_mz, na.rm = TRUE)
        },
        intensity = if (all(is.na(intensity))) {
          NA_real_
        } else {
          max(intensity, na.rm = TRUE)
        }
      ),
      by = .(Lab_ID, Polarity, mz_bin = round(frag_mz, bin_d))
    ]
    merged_dt[, mz_bin := NULL]
    data.table::setcolorder(merged_dt, names(spectra_tb))
    spectra_tb <- data.table::rbindlist(
      list(spectra_tb, merged_dt),
      use.names = TRUE
    )
    data.table::setkeyv(spectra_tb, c("Lab_ID", "Polarity", "NCE"))
  }

  f_list <- list()
  fi <- 0L
  for (pol_key in names(fragment_list)) {
    pol_short <- toupper(as.character(pol_key)[1L])
    if (!pol_short %in% c("POS", "NEG")) {
      next
    }
    by_lab <- fragment_list[[pol_key]]
    if (!is.list(by_lab)) {
      next
    }
    for (lab_id in names(by_lab)) {
      by_ft <- by_lab[[lab_id]]
      if (!is.list(by_ft)) {
        next
      }
      for (ft_key in names(by_ft)) {
        chunk <- by_ft[[ft_key]]
        if (is.null(chunk) || !is.data.frame(chunk)) {
          next
        }
        if (!nrow(chunk)) {
          next
        }
        tb <- .nts_fragment_feat_normalize_cols(chunk)
        tb[, Lab_ID := as.character(lab_id)]
        tb[, Polarity := pol_short]
        tb[, Feature_Type := .nts_feature_type_map(ft_key)]
        fi <- fi + 1L
        f_list[[fi]] <- tb
      }
    }
  }

  if (!fi) {
    features_tb <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      Feature_Type = character(),
      mz_value = numeric(),
      Frequency = numeric()
    )
  } else {
    features_tb <- data.table::rbindlist(f_list, fill = TRUE)
    .nts_rename_first_present(
      features_tb,
      c("mz", "MZ", "m.z", "mass", "NL", "MD", "DF", "nl", "md", "df"),
      "mz_value"
    )
    .nts_rename_first_present(
      features_tb,
      c("frequency", "freq", "Frequency", "rel_frequency", "int", "Int"),
      "Frequency"
    )
    fneed <- c("Lab_ID", "Polarity", "Feature_Type", "mz_value", "Frequency")
    fmiss <- setdiff(fneed, names(features_tb))
    if (length(fmiss)) {
      stop(
        "features_long is missing required column(s) after normalisation: ",
        paste(fmiss, collapse = ", "),
        call. = FALSE
      )
    }
    features_tb <- features_tb[, fneed, with = FALSE]
    features_tb[, Lab_ID := as.character(Lab_ID)]
    features_tb[, Polarity := as.character(Polarity)]
    features_tb[, Feature_Type := as.character(Feature_Type)]
    features_tb[, mz_value := as.numeric(mz_value)]
    features_tb[, Frequency := as.numeric(Frequency)]
    features_tb <- .nts_apply_features_long_mz_digits(features_tb, feat_dmap)
  }

  data.table::setkeyv(features_tb, c("Lab_ID", "Polarity", "Feature_Type"))

  # Optional RT prediction
  if (isTRUE(predict_rt)) {
    if (!inherits(rt_model, "nts_rt_model")) {
      stop("predict_rt = TRUE requires a valid rt_model (nts_rt_model object).",
           call. = FALSE)
    }
    info_dt <- .nts_add_predicted_rt(info_dt, rt_model)
  }

  methods::new(
    "NTSDatabase",
    metadata = info_dt,
    spectra_long = spectra_tb,
    features_long = features_tb,
    version = ver
  )
}

#' ---- merged from 02_pfas_database_optimize.R ----
#' @keywords internal
#' @noRd
.nts_spectra_require_dt <- function(spectra_tb) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  if (!data.table::is.data.table(spectra_tb)) {
    spectra_tb <- data.table::as.data.table(spectra_tb)
  }
  req <- c("Lab_ID", "Polarity", "NCE", "frag_mz", "intensity")
  miss <- setdiff(req, names(spectra_tb))
  if (length(miss)) {
    stop(
      "spectra_tb is missing column(s): ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  spectra_tb
}

#' @keywords internal
#' @noRd
.nts_metadata_join_precursor <- function(dt_ce, metadata_tb) {
  meta <- data.table::as.data.table(metadata_tb)
  if (!nrow(meta)) {
    dt_ce[, prec_mz := NA_real_]
    return(dt_ce)
  }
  if ("Lab.ID" %in% names(meta) && !"Lab_ID" %in% names(meta)) {
    data.table::setnames(meta, "Lab.ID", "Lab_ID", skip_absent = TRUE)
  }
  if (!"Lab_ID" %in% names(meta)) {
    stop(
      "metadata_tb must contain Lab_ID (or Lab.ID) for precursor join",
      call. = FALSE
    )
  }
  keep <- c(
    "Lab_ID",
    intersect(
      names(meta),
      c("Precursor_mz", "target_mz", "mz_pos", "mz_neg", "mz.neg", "mz.pos", "mz", "mass")
    )
  )
  keep <- unique(keep)
  meta_sub <- unique(meta[, ..keep])
  out <- merge(dt_ce, meta_sub, by = "Lab_ID", all.x = TRUE, sort = FALSE)
  out[, prec_mz := NA_real_]
  if ("target_mz" %in% names(out)) {
    tm <- suppressWarnings(as.numeric(out$target_mz))
    ok <- is.finite(tm) & tm > 0
    out[ok, prec_mz := tm[ok]]
  }
  miss <- !is.finite(out$prec_mz) | out$prec_mz <= 0
  if (any(miss) && "Precursor_mz" %in% names(out)) {
    v <- suppressWarnings(as.numeric(out$Precursor_mz))
    ok <- miss & is.finite(v) & v > 0
    out[ok, prec_mz := v[ok]]
    miss <- !is.finite(out$prec_mz) | out$prec_mz <= 0
  }
  if (any(miss) && all(c("mz_pos", "mz_neg") %in% names(out))) {
    out[miss, prec_mz := ifelse(
      Polarity == "POS",
      as.numeric(mz_pos),
      as.numeric(mz_neg)
    )]
    miss <- !is.finite(out$prec_mz) | out$prec_mz <= 0
  } else if (any(miss)) {
    for (pol in c("NEG", "POS")) {
      col <- if (pol == "NEG") {
        intersect(c("mz_neg", "mz.neg"), names(out))[1L]
      } else {
        intersect(c("mz_pos", "mz.pos"), names(out))[1L]
      }
      if (is.na(col)) {
        next
      }
      v <- suppressWarnings(as.numeric(out[[col]]))
      ok <- miss & out$Polarity == pol & is.finite(v) & v > 0
      out[ok, prec_mz := v[ok]]
      miss <- !is.finite(out$prec_mz) | out$prec_mz <= 0
    }
  }
  if (any(miss) && "mz" %in% names(out)) {
    v <- suppressWarnings(as.numeric(out$mz))
    ok <- miss & is.finite(v) & v > 0
    out[ok, prec_mz := v[ok]]
    miss <- !is.finite(out$prec_mz) | out$prec_mz <= 0
  }
  if (any(miss) && "mass" %in% names(out)) {
    proton <- 1.007276
    m <- suppressWarnings(as.numeric(out$mass))
    out[miss & out$Polarity == "POS", prec_mz := m + proton]
    out[miss & out$Polarity != "POS", prec_mz := m - proton]
  }
  drop <- intersect(
    c("Precursor_mz", "target_mz", "mz_pos", "mz_neg", "mz.neg", "mz.pos", "mz", "mass"),
    names(out)
  )
  if (length(drop)) {
    out[, (drop) := NULL]
  }
  out
}

#' Denoise MS2 long table and rebuild consensus spectrum
#'
#' Drops prior `NCE == "Merged"` layers, micro-centroids duplicate m/z bins
#' (same `Lab_ID` / `Polarity` / `NCE`: mean `frag_mz`, max `intensity`), then
#' filters by relative abundance (>= 1%), retains top 15 peaks per group, and
#' recomputes a merged spectrum per compound and polarity.
#'
#' @param spectra_long A `data.table` with `Lab_ID`, `Polarity`, `NCE`,
#'   `frag_mz`, `intensity`.
#' @param mz_bin_digits Decimal places for `round(m/z, d)` when `mz_bin_width`
#'   is not set (default `4`).
#' @param mz_bin_width Optional fixed m/z bin width in Da (e.g. `0.005`).
#' @param add_merged_layer If `TRUE` (default), append `NCE == \"Merged\"`.
#'
#' @return A keyed `data.table` (same schema as input).
#' @export
optimize_spectra <- function(
    spectra_long,
    mz_bin_digits = 4L,
    mz_bin_width = NULL,
    add_merged_layer = TRUE) {
  dt <- .nts_spectra_require_dt(spectra_long)
  bin_d <- suppressWarnings(as.integer(mz_bin_digits)[1L])
  if (!is.finite(bin_d) || bin_d < 0L) {
    bin_d <- 4L
  }
  spec_bin_w <- suppressWarnings(as.numeric(mz_bin_width)[1L])
  if (!is.finite(spec_bin_w) || spec_bin_w <= 0) {
    spec_bin_w <- NA_real_
  }
  add_merged <- isTRUE(add_merged_layer)
  dt <- dt[NCE != "Merged"]
  if (!nrow(dt)) {
    empty <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      NCE = character(),
      frag_mz = numeric(),
      intensity = numeric()
    )
    data.table::setkeyv(empty, c("Lab_ID", "Polarity", "NCE"))
    return(empty)
  }

  dt[, mz_bin := .nts_mz_assign_bin(frag_mz, bin_width = spec_bin_w, digits = bin_d)]
  dt <- dt[
    ,
    .(
      frag_mz = if (all(is.na(frag_mz))) {
        NA_real_
      } else {
        mean(frag_mz, na.rm = TRUE)
      },
      intensity = if (all(is.na(intensity))) {
        NA_real_
      } else {
        max(intensity, na.rm = TRUE)
      }
    ),
    by = .(Lab_ID, Polarity, NCE, mz_bin)
  ]
  dt[, mz_bin := NULL]

  dt[, max_int := max(intensity, na.rm = TRUE), by = .(Lab_ID, Polarity, NCE)]
  dt <- dt[
    is.finite(max_int) & max_int > 0 & is.finite(intensity) &
      ((intensity / max_int) >= 0.01)
  ]
  if (!nrow(dt)) {
    empty <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      NCE = character(),
      frag_mz = numeric(),
      intensity = numeric()
    )
    data.table::setkeyv(empty, c("Lab_ID", "Polarity", "NCE"))
    return(empty)
  }

  dt <- dt[
    order(Lab_ID, Polarity, NCE, -intensity),
    head(.SD, 15L),
    by = .(Lab_ID, Polarity, NCE)
  ]
  dt[, max_int := NULL]

  if (add_merged) {
    merged_dt <- dt[
      ,
      .(
        NCE = "Merged",
        frag_mz = if (all(is.na(frag_mz))) {
          NA_real_
        } else {
          mean(frag_mz, na.rm = TRUE)
        },
        intensity = if (all(is.na(intensity))) {
          NA_real_
        } else {
          max(intensity, na.rm = TRUE)
        }
      ),
      by = .(
        Lab_ID,
        Polarity,
        mz_bin = .nts_mz_assign_bin(frag_mz, bin_width = spec_bin_w, digits = bin_d)
      )
    ]
    merged_dt[, mz_bin := NULL]
    data.table::setcolorder(merged_dt, names(dt))
    out <- data.table::rbindlist(list(dt, merged_dt), use.names = TRUE)
  } else {
    out <- dt
  }
  data.table::setkeyv(out, c("Lab_ID", "Polarity", "NCE"))
  out
}

#' Mine DF / NL / MD features with CE-frequency via `uniqueN(NCE)`
#'
#' Uses only rows with `NCE != "Merged"`. For each final feature key
#' (`Lab_ID`, `Polarity`, `Feature_Type`, `mz_value`), **`NCE_Count`** is
#' `uniqueN(NCE)` over contributing rows (i.e. the number of distinct collision
#' energies in which that feature appeared, matching the original nested-loop
#' counter without per-CE iteration in R).
#'
#' @details The count column is named `NCE_Count` (not `Frequency`) on purpose.
#'   `[build_nts_database()]` and [compile_reference_database()] emit a
#'   `features_long$Frequency` column that is a **relative intensity** in
#'   `[0, 1]`, and this function's raw output is a **count** of collision
#'   energies. They are different quantities, so they must not share a name.
#'
#' @param spectra_long Cleaned `data.table` from [optimize_spectra].
#' @param metadata_tb Compound table with `Lab_ID` and precursor columns:
#'   `target_mz` when present, else `Precursor_mz`, or `mz_pos` / `mz_neg` (by
#'   `Polarity`), or `mz` / `mass`.
#' @param mz_value_digits Same convention as [build_nts_database()]'s
#'   `features_long_mz_digits` (default DF/NL = 5, MD = 3).
#'
#' @return Keyed `data.table` with `Lab_ID`, `Polarity`, `Feature_Type`,
#'   `mz_value`, `NCE_Count`.
#' @export
mine_features <- function(spectra_long, metadata_tb, mz_value_digits = NULL) {
  dt_ce <- .nts_spectra_require_dt(spectra_long)
  dig <- .nts_resolve_features_long_mz_digits(mz_value_digits)
  dt_ce <- dt_ce[NCE != "Merged"]
  if (!nrow(dt_ce)) {
    feat <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      Feature_Type = character(),
      mz_value = numeric(),
      NCE_Count = numeric()
    )
    data.table::setkeyv(feat, c("Lab_ID", "Polarity", "Feature_Type"))
    return(feat)
  }

  dt_ce <- .nts_metadata_join_precursor(dt_ce, metadata_tb)

  df_dt <- dt_ce[
    is.finite(frag_mz) & is.finite(prec_mz) & frag_mz < prec_mz,
    .(
      Lab_ID,
      Polarity,
      NCE,
      Feature_Type = "DF",
      mz_value = round(frag_mz, dig[["DF"]])
    )
  ]

  nl_dt <- dt_ce[
    is.finite(frag_mz) & is.finite(prec_mz),
    .(
      Lab_ID,
      Polarity,
      NCE,
      Feature_Type = "NL",
      mz_value = round(prec_mz - frag_mz, dig[["NL"]])
    )
  ]
  nl_dt <- nl_dt[is.finite(mz_value) & mz_value > 0]

  dt_md_base <- unique(dt_ce[
    is.finite(frag_mz),
    .(Lab_ID, Polarity, NCE, frag_mz)
  ])
  md_dt <- merge(
    dt_md_base,
    dt_md_base,
    by = c("Lab_ID", "Polarity", "NCE"),
    allow.cartesian = TRUE,
    suffixes = c("_1", "_2")
  )
  md_dt <- md_dt[frag_mz_1 > frag_mz_2]
  md_dt[, mz_value := round(frag_mz_1 - frag_mz_2, dig[["MD"]])]
  md_dt[, Feature_Type := "MD"]
  md_dt <- md_dt[, .(Lab_ID, Polarity, NCE, Feature_Type, mz_value)]

  all_features <- data.table::rbindlist(
    list(df_dt, nl_dt, md_dt),
    use.names = TRUE
  )
  if (!nrow(all_features)) {
    feat <- data.table::data.table(
      Lab_ID = character(),
      Polarity = character(),
      Feature_Type = character(),
      mz_value = numeric(),
      NCE_Count = numeric()
    )
    data.table::setkeyv(feat, c("Lab_ID", "Polarity", "Feature_Type"))
    return(feat)
  }

  feat_tb <- all_features[
    ,
    .(NCE_Count = uniqueN(NCE)),
    by = .(Lab_ID, Polarity, Feature_Type, mz_value)
  ]
  feat_tb <- feat_tb[
    order(Lab_ID, Polarity, Feature_Type, -NCE_Count)
  ]
  data.table::setkeyv(feat_tb, c("Lab_ID", "Polarity", "Feature_Type"))
  feat_tb
}

#' Optimize a [NTSDatabase] in memory (spectra + mined features)
#'
#' Runs [optimize_spectra] on `db@spectra_long`, then [mine_features] using
#' `db@metadata`, and returns a new S4 object with updated `spectra_long` and
#' `features_long` (metadata and `version` are preserved). Note that
#' `features_long` from [mine_features] carries an `NCE_Count` column (number of
#' collision energies a feature was seen in), **not** the relative-intensity
#' `Frequency` produced by [build_nts_database].
#'
#' @param db An [NTSDatabase] instance.
#' @param spectra_mz_bin_digits,features_long_mz_digits Passed to
#'   [optimize_spectra()] and [mine_features()] respectively.
#' @param spectra_mz_bin_width,add_merged_layer Passed to [optimize_spectra()].
#' @param features_mz_bin_width If set (e.g. `0.005`), re-bins `features_long$mz_value`
#'   after [mine_features()] the same way as [build_nts_database()].
#'
#' @return A new [NTSDatabase].
#' @export
optimize_reference_database <- function(db,
                                          spectra_mz_bin_digits = 4L,
                                          spectra_mz_bin_width = NULL,
                                          add_merged_layer = TRUE,
                                          features_long_mz_digits = NULL,
                                          features_mz_bin_width = NULL) {
  if (!.is_reference_spectral_db(db)) {
    stop("`db` must be an NTSDatabase (or legacy PFASDatabase) object", call. = FALSE)
  }
  clean_spectra <- optimize_spectra(
    methods::slot(db, "spectra_long"),
    mz_bin_digits = spectra_mz_bin_digits,
    mz_bin_width = spectra_mz_bin_width,
    add_merged_layer = add_merged_layer
  )
  new_features <- mine_features(
    clean_spectra,
    methods::slot(db, "metadata"),
    mz_value_digits = features_long_mz_digits
  )
  feat_bin_w <- suppressWarnings(as.numeric(features_mz_bin_width)[1L])
  if (is.finite(feat_bin_w) && feat_bin_w > 0 && nrow(new_features)) {
    new_features <- data.table::as.data.table(new_features)
    new_features[, mz_value := .nts_mz_assign_bin(mz_value, bin_width = feat_bin_w)]
    new_features <- new_features[
      ,
      .(NCE_Count = max(NCE_Count, na.rm = TRUE)),
      by = .(Lab_ID, Polarity, Feature_Type, mz_value)
    ]
    data.table::setkeyv(new_features, c("Lab_ID", "Polarity", "Feature_Type"))
  }
  methods::new(
    "NTSDatabase",
    metadata = methods::slot(db, "metadata"),
    spectra_long = clean_spectra,
    features_long = new_features,
    version = methods::slot(db, "version")
  )
}

#' ---- merged from 13_structure_annotation.R ----
#' @keywords internal
#' @noRd
.nts_rcdk_stop <- function() {
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop(
      "rcdk is required: install.packages(\"rcdk\") (see package README for CDK/JVM setup).",
      call. = FALSE
    )
  }
  if (!requireNamespace("fingerprint", quietly = TRUE)) {
    stop("fingerprint package is required (rcdk fingerprint similarity).", call. = FALSE)
  }
}

#' Reconstruct [fingerprint::fingerprint-class] from stored on-bit indices + length.
#' @keywords internal
#' @noRd
.nts_fp_from_onbits <- function(onbits, nbit) {
  nb <- as.integer(nbit)[1L]
  if (!is.finite(nb) || nb < 1L) {
    nb <- 1024L
  }
  ob <- as.numeric(onbits)
  if (!length(ob)) {
    ob <- numeric(0)
  }
  methods::new("fingerprint", nbit = nb, bits = ob, provider = "CDK", name = "")
}

#' Similarity between two on-bit representations (same `nbit` as library build).
#' @keywords internal
#' @noRd
.nts_fingerprint_similarity <- function(onbits_a, onbits_b, nbit, method = "dice") {
  fa <- .nts_fp_from_onbits(onbits_a, nbit)
  fb <- .nts_fp_from_onbits(onbits_b, nbit)
  tryCatch(
    as.numeric(fingerprint::distance(fa, fb, method = method)),
    error = function(e) NA_real_
  )
}

#' Try full SMILES then dot-separated fragments (salts).
#' @keywords internal
#' @noRd
.nts_mol_first_from_smiles_rcdk <- function(smiles) {
  parts <- trimws(strsplit(smiles, ".", fixed = TRUE)[[1L]])
  parts <- parts[nzchar(parts)]
  cand <- unique(c(smiles, parts))
  cand <- cand[nzchar(cand)]
  for (s in cand) {
    mol_list <- tryCatch(rcdk::parse.smiles(s), error = function(e) NULL)
    if (is.null(mol_list) || length(mol_list) == 0L) {
      next
    }
    mol <- mol_list[[1L]]
    if (!is.null(mol)) {
      return(mol)
    }
  }
  NULL
}

# ============================================================================
# Unified Fingerprint / Similarity Interface (Plan A)
# All fingerprint & similarity operations go through these two entry points.
# ============================================================================

#' Compute fingerprint for one SMILES with optional timeout protection.
#'
#' Returns `list(exact_mass, onbits, fp_nbits)` or `NULL` on failure/timeout.
#' Uses `.nts_mol_first_from_smiles_rcdk` for salt-tolerant SMILES parsing.
#'
#' @param smiles Character scalar.
#' @param fp_type rcdk fingerprint type (e.g. "circular", "maccs", "substructure").
#' @param timeout_sec Seconds before aborting. `0` or `NULL` = no timeout.
#' @keywords internal
#' @noRd
.nts_fp_compute_one <- function(smiles, fp_type = "circular", timeout_sec = 30L) {
  .do <- function() {
    mol <- .nts_mol_first_from_smiles_rcdk(smiles)
    if (is.null(mol)) return(NULL)
    mass <- tryCatch(rcdk::get.exact.mass(mol), error = function(e) NA_real_)
    if (!is.finite(mass)) return(NULL)
    fp <- tryCatch(
      rcdk::get.fingerprint(mol, type = fp_type),
      error = function(e) NULL
    )
    if (is.null(fp)) return(NULL)
    nbit <- methods::slot(fp, "nbit")
    bits <- methods::slot(fp, "bits")
    on_bits <- as.integer(bits)
    list(exact_mass = mass, onbits = on_bits, fp_nbits = as.integer(nbit))
  }

  ts <- suppressWarnings(as.numeric(timeout_sec)[1L])
  if (is.finite(ts) && ts > 0 && requireNamespace("R.utils", quietly = TRUE)) {
    tryCatch(
      R.utils::withTimeout(.do(), timeout = ts, onTimeout = "error"),
      error = function(e) NULL
    )
  } else {
    tryCatch(.do(), error = function(e) NULL)
  }
}

#' Compute fingerprint object (fingerprint-class) for one SMILES with timeout.
#'
#' Returns a `fingerprint::fingerprint` object or `NULL`.
#' Used by pairwise comparison functions (e.g. `infer_by_structure`).
#'
#' @param smiles Character scalar.
#' @param fp_type rcdk fingerprint type.
#' @param timeout_sec Seconds before aborting.
#' @keywords internal
#' @noRd
.nts_fp_object_one <- function(smiles, fp_type = "circular", timeout_sec = 30L) {
  .do <- function() {
    if (!nzchar(trimws(as.character(smiles)))) return(NULL)
    .nts_rcdk_stop()
    mol <- rcdk::parse.smiles(trimws(as.character(smiles)))
    if (length(mol) > 1L) mol <- mol[[1L]]
    rcdk::get.fingerprint(mol, type = as.character(fp_type)[1L])
  }

  ts <- suppressWarnings(as.numeric(timeout_sec)[1L])
  if (is.finite(ts) && ts > 0 && requireNamespace("R.utils", quietly = TRUE)) {
    tryCatch(
      R.utils::withTimeout(.do(), timeout = ts, onTimeout = "error"),
      error = function(e) NULL
    )
  } else {
    tryCatch(.do(), error = function(e) NULL)
  }
}

#' Compute pairwise similarity between two fingerprint objects.
#'
#' Wraps `fingerprint::distance` with a configurable method.
#' Falls back to NA on error.
#'
#' @param fp1,fp2 `fingerprint::fingerprint` objects.
#' @param method "tanimoto", "dice", or "cosine".
#' @keywords internal
#' @noRd
.nts_fp_similarity_pair <- function(fp1, fp2, method = "tanimoto") {
  if (is.null(fp1) || is.null(fp2)) return(NA_real_)
  tryCatch(
    as.numeric(fingerprint::distance(fp1, fp2, method = method)),
    error = function(e) NA_real_
  )
}

#' Compute similarity from dense bit vectors (matrix multiplication path).
#'
#' Used by `search_asymmetric_smiles` and `build_structural_library`
#' (compute_similarity_matrix mode). Supports dice, tanimoto, cosine.
#'
#' @param Q Query matrix (M x L, 0/1).
#' @param T_mat Target matrix (N x L, 0/1).
#' @param method "dice", "tanimoto", or "cosine".
#' @return M x N similarity matrix.
#' @keywords internal
#' @noRd
.nts_sim_matrix_from_bits <- function(Q, T_mat, method = "dice") {
  storage.mode(Q) <- "double"
  storage.mode(T_mat) <- "double"
  inter_mat <- Q %*% t(T_mat)
  query_sums <- rowSums(Q)
  target_sums <- rowSums(T_mat)

  sim_meth <- tolower(trimws(as.character(method)[1L]))
  if (sim_meth %in% c("tanimoto", "jaccard")) {
    denom <- outer(query_sums, target_sums, "+") - inter_mat
    sim_mat <- ifelse(denom > 0, inter_mat / denom, 0)
  } else if (sim_meth == "cosine") {
    denom <- outer(sqrt(query_sums), sqrt(target_sums), "*")
    sim_mat <- ifelse(denom > 0, inter_mat / denom, 0)
  } else {
    # Default: Dice
    denom_base <- outer(query_sums, target_sums, "+")
    sim_mat <- ifelse(denom_base > 0, 2 * inter_mat / denom_base, 0)
  }
  sim_mat[!is.finite(sim_mat)] <- 0
  sim_mat
}

#' Legacy wrapper - delegates to `.nts_fp_compute_one`.
#' @keywords internal
#' @noRd
.nts_fp_info_from_smiles_rcdk <- function(smiles, fingerprint_type) {
  mol <- .nts_mol_first_from_smiles_rcdk(smiles)
  if (is.null(mol)) return(NULL)
  mass <- tryCatch(rcdk::get.exact.mass(mol), error = function(e) NA_real_)
  if (!is.finite(mass)) return(NULL)
  fp <- tryCatch(
    rcdk::get.fingerprint(mol, type = fingerprint_type),
    error = function(e) NULL
  )
  if (is.null(fp)) return(NULL)
  nbit <- methods::slot(fp, "nbit")
  bits <- methods::slot(fp, "bits")
  on_bits <- as.integer(bits)
  list(exact_mass = mass, onbits = on_bits, fp_nbits = as.integer(nbit))
}

#' Normalize one SMILES cell to a single string or NA (placeholders removed).
#' @keywords internal
#' @noRd
.nts_normalize_smiles_string <- function(x) {
  if (length(x) != 1L) {
    return(NA_character_)
  }
  x <- as.character(x)[1L]
  if (is.na(x)) {
    return(NA_character_)
  }
  x <- trimws(x)
  if (!nzchar(x)) {
    return(NA_character_)
  }
  if (x %in% c("-", ".", "--", "\u2014", "\u2013")) {
    return(NA_character_)
  }
  ul <- toupper(x)
  if (ul %in% c("NA", "N/A", "NONE", "NULL", "UNKNOWN")) {
    return(NA_character_)
  }
  x
}

#' Build a static structure library with CDK fingerprints (rcdk)
#'
#' Reads a CSV of SMILES, parses each structure with **rcdk** (CDK), computes
#' exact mass and a user-selected binary fingerprint (`fingerprint_type`), and
#' saves a `data.table` to disk. Fingerprint is stored as **on-bit indices**
#' (`fp_onbits` list column) for portable RDS and later
#' [annotate_by_dual_network]. `fp_nbits` matches the fingerprinter size for
#' the chosen type (e.g. substructure pattern count from CDK).
#'
#' Attribute `fingerprint_type` is set on the returned table (and preserved
#' in `output_rds`) so propagation uses the same fingerprinter as the library.
#'
#' Parsing tries the full SMILES string, then dot-separated fragments (salts).
#'
#' Requires CRAN packages **rcdk** and **fingerprint**.
#'
#' @param csv_path Path to CSV.
#' @param smiles_col Column name holding the SMILES string.
#' @param id_col Column name holding the compound identifier.
#' @param fallback_smiles_col Optional second column (e.g. `"SMILES"`) used when
#'   `smiles_col` is empty, common for OECD lists where `QSAR_READY_SMILES` is
#'   missing for some rows.
#' @param fingerprint_type Passed to [rcdk::get.fingerprint] `type` (default
#'   `"substructure"`). Examples: `"standard"`, `"maccs"`, `"circular"`, etc.
#' @param output_rds Path for `saveRDS` output (parent dirs are created if needed).
#' @param max_smi_len Maximum SMILES length (characters); longer rows are dropped
#'   before parsing.
#' @param timeout_sec Per-molecule timeout (seconds) for fingerprint computation.
#' @param compute_similarity_matrix If `TRUE`, also precompute an all-vs-all
#'   similarity matrix for the library. Default `FALSE`.
#' @param similarity_method Similarity metric used when
#'   `compute_similarity_matrix = TRUE` (e.g. `"dice"`).
#'
#' @details
#' Columns are read with [data.table::fread], then copied to a plain
#' `data.frame` and extracted as **atomic character vectors** (`id_vec`,
#' `smi_vec`) before the rcdk loop, avoiding `data.table` row-subset edge cases.
#' Cells equal to `NA`, `""`, `"-"`, `"."`, etc. are dropped before parsing.
#' After load, prints a short message. In interactive sessions, a
#' [utils::txtProgressBar] (style 3) tracks rows with usable SMILES; it is
#' suppressed when `interactive()` is false (e.g. `devtools::test()`, CRAN
#' checks) to avoid log noise.
#'
#' @return A `data.table` with columns `ID`, `SMILES`, `ExactMass`, `fp_onbits`
#'   (list of integer vectors), `fp_nbits` (integer). Attribute
#'   `fingerprint_type` records the fingerprinter used.
#' @export
build_structural_library <- function(csv_path,
                                          smiles_col = "SMILES",
                                          id_col = "ID",
                                          fallback_smiles_col = NULL,
                                          fingerprint_type = "circular",
                                          output_rds = "OECD_Structural_DB.rds",
                                          max_smi_len = 500L,
                                          timeout_sec = 30L,
                                          compute_similarity_matrix = FALSE,
                                          similarity_method = "dice") {
  csv_path <- path.expand(as.character(csv_path)[1L])
  .nts_check(csv_path, "file_exists", .arg = "csv_path")
  .nts_rcdk_stop()
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  fp_type <- as.character(fingerprint_type)[1L]
  if (!nzchar(fp_type)) {
    stop("fingerprint_type must be a non-empty string", call. = FALSE)
  }

  max_len <- suppressWarnings(as.integer(max_smi_len)[1L])
  if (!is.finite(max_len) || max_len < 1L) max_len <- 500L
  ts <- suppressWarnings(as.numeric(timeout_sec)[1L])
  if (!is.finite(ts) || ts < 0) ts <- 30

  dt <- data.table::fread(csv_path, na.strings = c("", "NA", "N/A"))
  if (!nrow(dt)) {
    stop("CSV has no rows", call. = FALSE)
  }
  if (!smiles_col %in% names(dt)) {
    stop("Column smiles_col not found: ", smiles_col, call. = FALSE)
  }
  if (!id_col %in% names(dt)) {
    stop("Column id_col not found: ", id_col, call. = FALSE)
  }
  if (!is.null(fallback_smiles_col)) {
    fb <- as.character(fallback_smiles_col)[1L]
    if (!nzchar(fb) || !fb %in% names(dt)) {
      stop("fallback_smiles_col not found in CSV: ", fallback_smiles_col, call. = FALSE)
    }
    if (identical(fb, smiles_col)) {
      fallback_smiles_col <- NULL
    }
  }

  df <- as.data.frame(dt, stringsAsFactors = FALSE)
  id_vec <- as.character(df[[id_col]])
  smi_vec <- as.character(df[[smiles_col]])
  name_col <- grep(
    "^Compound\\.name$|^compound_name$|^Name$|^Chemical_Name$|^Compound_Name$",
    names(df),
    ignore.case = TRUE,
    value = TRUE
  )[1L]
  cname_vec <- rep(NA_character_, length(id_vec))
  if (!is.na(name_col) && nzchar(name_col) && name_col %in% names(df)) {
    cname_vec <- trimws(as.character(df[[name_col]]))
  }
  n <- length(id_vec)
  if (length(smi_vec) != n) {
    stop("Internal error: id_col and smiles_col length mismatch.", call. = FALSE)
  }
  fb_vec <- NULL
  if (!is.null(fallback_smiles_col)) {
    fb_vec <- as.character(df[[fallback_smiles_col]])
    if (length(fb_vec) != n) {
      stop("Internal error: fallback_smiles_col length mismatch.", call. = FALSE)
    }
  }

  n_primary_ok <- 0L
  n_fallback_only <- 0L
  valid_idx <- integer(n)
  n_valid <- 0L
  for (j in seq_len(n)) {
    pr <- .nts_normalize_smiles_string(smi_vec[j])
    if (!is.na(pr)) {
      n_primary_ok <- n_primary_ok + 1L
      n_valid <- n_valid + 1L
      valid_idx[n_valid] <- j
    } else if (!is.null(fb_vec)) {
      fb <- .nts_normalize_smiles_string(fb_vec[j])
      if (!is.na(fb)) {
        n_fallback_only <- n_fallback_only + 1L
        n_valid <- n_valid + 1L
        valid_idx[n_valid] <- j
      }
    }
  }
  valid_idx <- valid_idx[seq_len(n_valid)]

  message(sprintf("Structural library data loaded. (%d rows, fingerprint='%s')",
                  n, fp_type))
  use_pb <- interactive()
  n_skip <- 0L

  # -- Pre-filter: build vectors of valid (index, id, smiles) --
  kept_idx <- integer(n)
  kept_smi <- character(n)
  kept_id  <- character(n)
  kept_cname <- character(n)
  nk <- 0L
  for (i in seq_len(n)) {
    smi <- .nts_normalize_smiles_string(smi_vec[i])
    if (is.na(smi) && !is.null(fb_vec)) {
      smi <- .nts_normalize_smiles_string(fb_vec[i])
    }
    if (is.na(smi)) {
      n_skip <- n_skip + 1L
      next
    }
    nch <- suppressWarnings(nchar(smi, allowNA = TRUE))
    if (is.na(nch) || nch > max_len) {
      n_skip <- n_skip + 1L
      next
    }
    nk <- nk + 1L
    kept_idx[nk] <- i
    kept_smi[nk] <- smi
    kept_id[nk]  <- id_vec[i]
    cn <- cname_vec[i]
    if (is.na(cn) || !nzchar(cn)) {
      cn <- id_vec[i]
    }
    kept_cname[nk] <- cn
  }
  kept_idx <- kept_idx[seq_len(nk)]
  kept_smi <- kept_smi[seq_len(nk)]
  kept_id  <- kept_id[seq_len(nk)]
  kept_cname <- kept_cname[seq_len(nk)]

  if (use_pb) {
    pb <- utils::txtProgressBar(min = 0L, max = nk, style = 3L)
    on.exit({ close(pb); cat("\n") }, add = TRUE)
  }

  # -- Chunk-based batch processing --
  chunk_sz <- 100L
  out_rows <- vector("list", nk)
  chunk_starts <- seq.int(1L, nk, by = chunk_sz)

  for (cs in chunk_starts) {
    ce <- min(nk, cs + chunk_sz - 1L)
    idx_chunk <- cs:ce
    smi_chunk <- kept_smi[idx_chunk]

    # Batch parse SMILES (single JNI call per chunk)
    mols <- tryCatch(
      rcdk::parse.smiles(smi_chunk),
      error = function(e) {
        # Fallback: parse one by one
        lapply(smi_chunk, function(s) {
          tryCatch({
            m <- rcdk::parse.smiles(s)
            if (is.list(m)) m[[1L]] else m
          }, error = function(e2) NULL)
        })
      }
    )
    if (!is.list(mols)) mols <- list(mols)

    for (j in seq_along(idx_chunk)) {
      k <- idx_chunk[j]
      mol <- mols[[j]]
      if (is.null(mol)) {
        n_skip <- n_skip + 1L
        if (use_pb) utils::setTxtProgressBar(pb, k)
        next
      }

      info <- tryCatch({
        mass <- rcdk::get.exact.mass(mol)
        if (!is.finite(mass)) stop("no mass")
        fp <- rcdk::get.fingerprint(mol, type = fp_type)
        nbit <- methods::slot(fp, "nbit")
        bits <- methods::slot(fp, "bits")
        list(exact_mass = mass, onbits = as.integer(bits), fp_nbits = as.integer(nbit))
      }, error = function(e) NULL)

      if (is.null(info)) {
        n_skip <- n_skip + 1L
      } else {
        out_rows[[k]] <- data.table::data.table(
          ID = kept_id[k],
          Compound.name = kept_cname[k],
          SMILES = kept_smi[k],
          ExactMass = info$exact_mass,
          fp_onbits = list(info$onbits),
          fp_nbits = info$fp_nbits
        )
      }
      if (use_pb) utils::setTxtProgressBar(pb, k)
    }

    # JVM + R garbage collection between chunks
    rm(mols)
    tryCatch({
      sys <- rJava::J("java/lang/System")
      rJava::.jcall(sys, "V", "gc")
    }, error = function(e) invisible(NULL))
    gc(verbose = FALSE)
  }
  bridge_note <- ""
  if (n_skip > 0L) {
    warning(
      paste0(
        sprintf(
          paste0(
            "build_structural_library: skipped %d of %d row(s) ",
            "(empty/placeholder SMILES, too long, timeout, or rcdk error). ",
            "Pre-scan: %d usable primary column cell(s), %d only via fallback."
          ),
          n_skip,
          n,
          n_primary_ok,
          n_fallback_only
        ),
        bridge_note
      ),
      call. = FALSE
    )
  }
  out <- data.table::rbindlist(out_rows[!vapply(out_rows, is.null, logical(1L))])
  if (!nrow(out)) {
    stop(
      paste0(
        "No valid SMILES rows in CSV (all ",
        n_skip,
        " row(s) failed). Rows=",
        n,
        "; usable in ",
        smiles_col,
        "=",
        n_primary_ok,
        "; usable only via fallback=",
        n_fallback_only,
        ". Check column names, placeholders (-, .), and rcdk/Java.",
        bridge_note
      ),
      call. = FALSE
    )
  }
  attr(out, "fingerprint_type") <- fp_type

  # -- Optional NxN similarity matrix --
  if (isTRUE(compute_similarity_matrix) && nrow(out) > 0L) {
    message("Computing N x N similarity matrix...")
    max_nbit <- max(out$fp_nbits, na.rm = TRUE)
    L <- max(max_nbit, 1024L)
    mat <- matrix(0L, nrow = nrow(out), ncol = L)
    for (k in seq_len(nrow(out))) {
      bits <- out$fp_onbits[[k]]
      bits <- bits[is.finite(bits) & bits >= 1L & bits <= L]
      if (length(bits)) mat[k, bits] <- 1L
    }
    sim_mat_out <- .nts_sim_matrix_from_bits(mat, mat, method = similarity_method)
    rownames(sim_mat_out) <- out$ID
    colnames(sim_mat_out) <- out$ID
    attr(out, "similarity_matrix") <- sim_mat_out
    attr(out, "similarity_method") <- similarity_method
    message(sprintf("Similarity matrix: %d x %d (%s)", nrow(sim_mat_out), ncol(sim_mat_out), similarity_method))
  }

  od <- dirname(output_rds)
  if (nzchar(od) && od != ".") {
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
  }
  saveRDS(out, output_rds)
  out
}

#' @keywords internal
#' @noRd
.nts_anchor_fp_onbits <- function(node_id, nodes_dt, db_dt, fingerprint_type) {
  nr <- nodes_dt[as.character(id) == as.character(node_id), ]
  cand_ids <- character(0)
  if (nrow(nr)) {
    if ("propagated_ID" %in% names(nr)) {
      x <- as.character(nr$propagated_ID[1L])
      if (nzchar(x) && !is.na(x)) {
        cand_ids <- c(cand_ids, x)
      }
    }
    if ("Library_ID" %in% names(nr)) {
      x <- as.character(nr$Library_ID[1L])
      if (nzchar(x) && !is.na(x)) {
        cand_ids <- c(cand_ids, x)
      }
    }
  }
  x <- as.character(node_id)
  if (nzchar(x) && !is.na(x)) {
    cand_ids <- c(cand_ids, x)
  }
  cand_ids <- unique(cand_ids[nzchar(cand_ids) & !is.na(cand_ids)])
  for (pid in cand_ids) {
    hit <- db_dt[as.character(ID) == as.character(pid), ]
    if (nrow(hit) == 1L && is.list(hit$fp_onbits) && length(hit$fp_onbits[[1L]])) {
      return(hit$fp_onbits[[1L]])
    }
  }
  if (!nrow(nr)) {
    return(NULL)
  }
  smi <- NA_character_
  if ("SMILES" %in% names(nr)) {
    smi <- as.character(nr$SMILES[1L])
  } else if ("smiles" %in% names(nr)) {
    smi <- as.character(nr$smiles[1L])
  }
  if (is.na(smi) || !nzchar(trimws(smi))) {
    return(NULL)
  }
  info <- tryCatch(
    .nts_fp_info_from_smiles_rcdk(trimws(smi), fingerprint_type),
    error = function(e) NULL
  )
  if (is.null(info)) {
    return(NULL)
  }
  info$onbits
}

#' Compound display name from a structural-library row (falls back to `ID`).
#' @keywords internal
#' @noRd
.nts_struct_lib_row_compound_name <- function(dt, k) {
  if (!is.data.frame(dt) && !data.table::is.data.table(dt)) {
    return(NA_character_)
  }
  for (cn in c(
    "Compound.name", "compound_name", "Name", "name",
    "Chemical_Name", "chemical_name", "Compound_Name"
  )) {
    if (cn %in% names(dt)) {
      v <- trimws(as.character(dt[[cn]][k]))
      if (!is.na(v) && nzchar(v)) {
        return(v)
      }
    }
  }
  if ("ID" %in% names(dt)) {
    return(trimws(as.character(dt[["ID"]][k])))
  }
  NA_character_
}

#' Ensure network node identity columns exist as character vectors.
#' @keywords internal
#' @noRd
.nts_network_nodes_ensure_identity_cols <- function(nodes) {
  for (col in c("ref_id", "Compound.name", "SMILES", "Library_ID", "Annotation")) {
    if (!col %in% names(nodes)) {
      nodes[, (col) := NA_character_]
    } else {
      nodes[, (col) := as.character(get(col))]
    }
  }
  nodes
}

#' @keywords internal
#' @noRd
.nts_edge_ms2_lookup_build <- function(edges_dt) {
  E <- data.table::as.data.table(edges_dt)
  if (!nrow(E)) {
    out <- data.table::data.table(
      ka = character(0),
      kb = character(0),
      sc = numeric(0)
    )
    data.table::setkey(out, ka, kb)
    return(out)
  }
  fa <- as.character(E$from)
  tb <- as.character(E$to)
  col_nm <- NA_character_
  if ("ms2_similarity" %in% names(E)) {
    col_nm <- "ms2_similarity"
  } else if ("MS2_Similarity" %in% names(E)) {
    col_nm <- "MS2_Similarity"
  } else if ("edge_score" %in% names(E)) {
    col_nm <- "edge_score"
  } else if ("weight" %in% names(E)) {
    col_nm <- "weight"
  }
  if (is.na(col_nm)) {
    sc <- rep(-Inf, nrow(E))
  } else {
    sc <- suppressWarnings(as.numeric(E[[col_nm]]))
    sc[!is.finite(sc)] <- -Inf
  }
  ka <- pmin(fa, tb)
  kb <- pmax(fa, tb)
  out <- data.table::data.table(ka = ka, kb = kb, sc = sc)
  out <- out[, .(sc = max(sc)), by = .(ka, kb)]
  out[!is.finite(sc), sc := -Inf]
  data.table::setkey(out, ka, kb)
  out
}

#' @keywords internal
#' @noRd
.nts_edge_ms2_lookup_score <- function(lk, node_a, node_b) {
  if (!nrow(lk)) {
    return(-Inf)
  }
  a <- as.character(node_a)[1L]
  b <- as.character(node_b)[1L]
  if (!nzchar(a) || !nzchar(b)) {
    return(-Inf)
  }
  x <- pmin(a, b)
  y <- pmax(a, b)
  row <- lk[list(x, y)]
  if (nrow(row) && is.finite(row$sc[1L])) {
    return(as.numeric(row$sc[1L]))
  }
  -Inf
}

#' @keywords internal
#' @noRd
.nts_lineage_node_label <- function(nodes_dt, seed_id_chr) {
  sid <- as.character(seed_id_chr)[1L]
  if ("name" %in% names(nodes_dt)) {
    v <- nodes_dt[id_chr == sid, as.character(name)][1L]
    if (length(v) && !is.na(v[1L]) && nzchar(v[1L])) {
      return(v[1L])
    }
  }
  sid
}

#' @keywords internal
#' @noRd
.nts_propagation_pick_better <- function(ms2_a, fp_a, ms2_b, fp_b) {
  if (!is.finite(ms2_a)) {
    ms2_a <- -Inf
  }
  if (!is.finite(ms2_b)) {
    ms2_b <- -Inf
  }
  if (!is.finite(fp_a)) {
    fp_a <- -Inf
  }
  if (!is.finite(fp_b)) {
    fp_b <- -Inf
  }
  if (ms2_a > ms2_b + 1e-12) {
    return(TRUE)
  }
  if (ms2_b > ms2_a + 1e-12) {
    return(FALSE)
  }
  fp_a > fp_b + 1e-12
}

#' @title Propagate structural library hits on a dual (MS1+MS2) molecular network
#'
#' @description
#' Performs breadth-first (layered) propagation: seed nodes on the graph are
#' matched against rows of a structural fingerprint RDS ([build_structural_library])
#' for each frontier-connected unknown, using neutral-mass filtering and
#' fingerprint similarity. Lineage from the winning seed is recorded in
#' `propagated_from_node` when multiple seeds compete (edge MS2 score, then
#' fingerprint score).
#'
#' @details
#' **Seed (layer 0) rule:** a node is a seed if any of: `annotation_level` matches
#' `"Database Confirmed"` (substring); or non-empty `Annotation`; or non-empty
#' `Library_ID`. The `SMILES` column on nodes is **not** used for seed detection.
#'
#' **Neutral mass:** `target_neutral_mass = mz - mass_offset` with
#' `mass_offset = -1.007276` for `"NEG"` (`[M-H]-`) and `+1.007276` for
#' positive mode, matching `ExactMass` in the RDS within `mz_tol_ppm`.
#'
#' **Output columns:** `Propagation_Level`, `network_propagation`,
#' `propagated_ID`, `propagated_SMILES`, `propagated_similarity`,
#' `propagated_from_node`. For each successfully propagated node, the same hit
#' is also copied into **`ref_id`**, **`Compound.name`** (from the structural
#' RDS name column when present, else `ID`), and **`SMILES`**, plus **`Library_ID`**
#' / **`Annotation`** when those columns exist on `nodes`. Existing `propagated_*` columns are
#' coerced to character / numeric before writing to avoid logical coercion bugs.
#'
#' **NSE / imports:** Uses `data.table` assignment by reference (`:=`) and joins
#' on `nodes` and `edges`; structural RDS must contain `ID`, `SMILES`,
#' `ExactMass`, `fp_onbits`, `fp_nbits`.
#'
#' @param net_data `list(nodes, edges)` from [generate_network_data].
#' @param structural_db_path Path to RDS from [build_structural_library].
#' @param mz_tol_ppm Ppm tolerance: `abs(ExactMass - neutral) / neutral * 1e6`.
#' @param min_similarity Minimum similarity (threshold for `similarity_method`).
#' @param similarity_method Passed to [fingerprint::distance] (e.g. `"dice"`,
#'   `"tanimoto"`). Default `"dice"`.
#' @param fingerprint_type Fingerprinter `type` for on-the-fly SMILES parsing
#'   when the node is not keyed in the library. Default: `attr(db, "fingerprint_type")`
#'   or `"substructure"`.
#' @param polarity `"auto"` (default: sniff from `net_data$polarity`,
#'   `net_data$info$polarity`, `net_data@polarity` when `net_data` is an
#'   [NTSProject] with that slot, or `project@polarity` when `project` is set),
#'   `"NEG"`, or `"POS"` for MS1 adduct to neutral mass. Manual values are
#'   normalised to a 3-letter prefix.
#' @param max_iterations Maximum BFS depth (layers `1..max_iterations`).
#' @param project Optional [NTSProject]; used for polarity sniffing only when
#'   `polarity = "auto"` and no polarity was found on `net_data`.
#'
#' @return The input `net_data` list: `nodes` is replaced with an updated
#'   `data.frame` (propagation columns in **Details**); `edges` is unchanged.
#'
#' @importFrom data.table := as.data.table is.data.table rbindlist
#' @export
annotate_by_dual_network <- function(net_data,
                                     structural_db_path,
                                     mz_tol_ppm = 5,
                                     min_similarity = 0.5,
                                     similarity_method = "dice",
                                     fingerprint_type = NULL,
                                     polarity = "auto",
                                     max_iterations = 3L,
                                     project = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }

  pol_in <- as.character(polarity)[1L]
  if (identical(tolower(trimws(pol_in)), "auto")) {
    inferred_pol <- NULL
    if (is.list(net_data) && !is.null(net_data$polarity)) {
      inferred_pol <- net_data$polarity
    } else if (is.list(net_data) && is.list(net_data$info) && !is.null(net_data$info$polarity)) {
      inferred_pol <- net_data$info$polarity
    } else if (methods::is(net_data, "NTSProject") && "polarity" %in% methods::slotNames(net_data)) {
      inferred_pol <- methods::slot(net_data, "polarity")
    } else if (!is.null(project) && methods::is(project, "NTSProject") && "polarity" %in% methods::slotNames(project)) {
      inferred_pol <- methods::slot(project, "polarity")
    }
    if (!is.null(inferred_pol)) {
      inferred_pol <- as.character(inferred_pol)[1L]
      if (nzchar(inferred_pol) && !is.na(inferred_pol)) {
        polarity <- toupper(substr(inferred_pol, 1L, 3L))
        message(sprintf("Auto-detected polarity: %s", polarity))
      } else {
        warning("Could not infer polarity from net_data; falling back to default: NEG", call. = FALSE)
        polarity <- "NEG"
      }
    } else {
      warning("Could not infer polarity from net_data; falling back to default: NEG", call. = FALSE)
      polarity <- "NEG"
    }
  } else {
    polarity <- toupper(substr(pol_in, 1L, 3L))
  }

  if (!is.list(net_data) || !all(c("nodes", "edges") %in% names(net_data))) {
    stop("net_data must be a list with nodes and edges", call. = FALSE)
  }

  structural_db_path <- path.expand(as.character(structural_db_path)[1L])
  .nts_check(structural_db_path, "file_exists", .arg = "structural_db_path")
  .nts_rcdk_stop()

  nodes <- data.table::as.data.table(net_data$nodes)
  edges <- data.table::as.data.table(net_data$edges)
  edge_ms2_lk <- .nts_edge_ms2_lookup_build(edges)
  if (!nrow(nodes) || !"id" %in% names(nodes)) {
    stop("net_data$nodes must be non-empty with column id", call. = FALSE)
  }
  if (!nrow(edges) || !all(c("from", "to") %in% names(edges))) {
    stop("net_data$edges must contain from and to", call. = FALSE)
  }
  if (!"annotation_level" %in% names(nodes)) {
    nodes[, annotation_level := NA_character_]
  }
  if (!"mz" %in% names(nodes)) {
    stop("net_data$nodes must contain mz", call. = FALSE)
  }

  db_dt <- readRDS(structural_db_path)
  if (!data.table::is.data.table(db_dt)) {
    db_dt <- data.table::as.data.table(db_dt)
  }
  req <- c("ID", "SMILES", "ExactMass", "fp_onbits", "fp_nbits")
  if (!all(req %in% names(db_dt))) {
    stop(
      "Structural DB must contain columns: ",
      paste(req, collapse = ", "),
      call. = FALSE
    )
  }
  if (!"Compound.name" %in% names(db_dt)) {
    db_dt[, Compound.name := trimws(as.character(ID))]
  }

  # Provenance-safe display name for a propagated hit. A structural library may
  # carry only IDs (or an empty name cell), in which case the propagation would
  # otherwise stamp `Compound.name` (and `Annotation`) with a bare ref_id such as
  # "PFAS_L53" -- which then masquerades as a database match in node labels /
  # tooltips. Resolve the human name from the project's `@db_annotations` by
  # `ref_id` when available; return NA when nothing is known so the caller leaves
  # the name blank instead of writing an ID into a name field.
  .name_from_ref <- local({
    nm <- NULL
    if (!is.null(project) && methods::is(project, "NTSProject")) {
      da <- tryCatch(
        as.data.frame(methods::slot(project, "db_annotations"), stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      if (is.data.frame(da) && nrow(da) && "ref_id" %in% names(da)) {
        rid <- trimws(as.character(da$ref_id))
        cn <- if ("Compound.name" %in% names(da)) {
          trimws(as.character(da$Compound.name))
        } else {
          rep(NA_character_, nrow(da))
        }
        ok <- !is.na(rid) & nzchar(rid) & !is.na(cn) &
          nzchar(cn) & toupper(cn) != "NA"
        if (any(ok)) {
          nm <- stats::setNames(cn[ok], rid[ok])
          nm <- nm[!duplicated(names(nm))]
        }
      }
    }
    function(rid) {
      rid <- trimws(as.character(rid)[1L])
      if (!is.null(nm) && length(nm) && rid %in% names(nm)) {
        as.character(nm[[rid]])
      } else {
        NA_character_
      }
    }
  })

  fp_type <- fingerprint_type
  if (is.null(fp_type) || !nzchar(as.character(fp_type)[1L])) {
    fp_type <- attr(db_dt, "fingerprint_type")
  }
  if (is.null(fp_type) || !nzchar(as.character(fp_type)[1L])) {
    fp_type <- "substructure"
  }
  fp_type <- as.character(fp_type)[1L]

  sim_method <- as.character(similarity_method)[1L]
  if (!nzchar(sim_method)) {
    stop("similarity_method must be non-empty", call. = FALSE)
  }

  .nts_check(mz_tol_ppm, "positive_finite", .arg = "mz_tol_ppm")
  ppm <- as.numeric(mz_tol_ppm)[1L]
  mt <- as.numeric(min_similarity)[1L]
  if (!is.finite(mt) || mt < 0 || mt > 1) {
    stop("min_similarity must be in [0, 1]", call. = FALSE)
  }

  pol <- gsub("[^A-Za-z]", "", as.character(polarity)[1L])
  pol <- toupper(pol)
  if (!nzchar(pol)) {
    pol <- "NEG"
  }
  mass_offset <- if (startsWith(pol, "NEG")) -1.007276 else 1.007276

  mi <- suppressWarnings(as.integer(max_iterations)[1L])
  if (!is.finite(mi) || mi < 1L) {
    stop("max_iterations must be a positive integer", call. = FALSE)
  }

  nodes[, id_chr := as.character(id)]
  nodes <- .nts_network_nodes_ensure_identity_cols(nodes)

  al <- as.character(nodes$annotation_level)
  al[is.na(al)] <- ""
  seed_mask <- grepl("Database Confirmed", al, fixed = FALSE)
  if ("Annotation" %in% names(nodes)) {
    a <- trimws(as.character(nodes$Annotation))
    seed_mask <- seed_mask | (!is.na(a) & nzchar(a))
  }
  if ("Library_ID" %in% names(nodes)) {
    lib <- trimws(as.character(nodes$Library_ID))
    seed_mask <- seed_mask | (!is.na(lib) & nzchar(lib))
  }
  if (!any(seed_mask)) {
    warning(
      "No seed nodes (Annotation, Library_ID, or Database Confirmed in annotation_level).",
      call. = FALSE
    )
    if (!"Propagation_Level" %in% names(nodes)) {
      nodes[, Propagation_Level := NA_integer_]
    }
    if (!"propagated_from_node" %in% names(nodes)) {
      nodes[, propagated_from_node := NA_character_]
    }
    nodes[, id_chr := NULL]
    net_data$nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
    return(net_data)
  }

  nodes[, Propagation_Level := NA_integer_]
  nodes[seed_mask, Propagation_Level := 0L]

  e1 <- edges[, .(a = as.character(from), b = as.character(to))]
  e2 <- edges[, .(a = as.character(to), b = as.character(from))]
  eb <- unique(data.table::rbindlist(list(e1, e2)))

  if (!"network_propagation" %in% names(nodes)) {
    nodes[, network_propagation := NA_character_]
  } else {
    nodes[, network_propagation := as.character(network_propagation)]
  }
  if (!"propagated_ID" %in% names(nodes)) {
    nodes[, propagated_ID := NA_character_]
  } else {
    nodes[, propagated_ID := as.character(propagated_ID)]
  }
  if (!"propagated_SMILES" %in% names(nodes)) {
    nodes[, propagated_SMILES := NA_character_]
  } else {
    nodes[, propagated_SMILES := as.character(propagated_SMILES)]
  }
  if (!"propagated_similarity" %in% names(nodes)) {
    nodes[, propagated_similarity := NA_real_]
  } else {
    nodes[, propagated_similarity := suppressWarnings(as.numeric(propagated_similarity))]
  }
  if (!"propagated_from_node" %in% names(nodes)) {
    nodes[, propagated_from_node := NA_character_]
  } else {
    nodes[, propagated_from_node := as.character(propagated_from_node)]
  }

  sim_lab <- paste0(toupper(substring(sim_method, 1L, 1L)), substring(sim_method, 2L))

  for (it in seq_len(mi)) {
    frontier <- nodes[Propagation_Level == it - 1L, unique(id_chr)]
    frontier <- frontier[!is.na(frontier) & nzchar(frontier)]
    if (!length(frontier)) {
      break
    }
    cand_targets <- unique(
      c(eb[a %in% frontier, b], eb[b %in% frontier, a])
    )
    cand_targets <- cand_targets[
      !is.na(cand_targets) & nzchar(cand_targets)
    ]
    open <- nodes[is.na(Propagation_Level), unique(id_chr)]
    targets <- intersect(cand_targets, open)
    if (!length(targets)) {
      break
    }

    updates <- vector("list", length(targets))
    names(updates) <- targets

    for (uid in targets) {
      nei <- unique(c(eb[a == uid, b], eb[b == uid, a]))
      seeds_conn <- intersect(nei, frontier)
      if (!length(seeds_conn)) {
        next
      }
      u_mz <- nodes[id_chr == uid, as.numeric(mz)][1L]
      if (!is.finite(u_mz) || u_mz <= 0) {
        next
      }
      target_neutral <- u_mz - mass_offset
      if (!is.finite(target_neutral) || target_neutral <= 0) {
        next
      }
      cand <- db_dt[
        is.finite(ExactMass) &
          abs(ExactMass - target_neutral) / target_neutral * 1e6 <= ppm
      ]
      if (!nrow(cand)) {
        next
      }

      best <- NULL
      for (sid in seeds_conn) {
        ms2_sid <- .nts_edge_ms2_lookup_score(edge_ms2_lk, uid, sid)
        ob_a <- .nts_anchor_fp_onbits(sid, nodes, db_dt, fp_type)
        if (is.null(ob_a) || !length(ob_a)) {
          next
        }
        for (k in seq_len(nrow(cand))) {
          ob_c <- cand$fp_onbits[[k]]
          if (!length(ob_c)) {
            next
          }
          nb <- as.integer(cand$fp_nbits[k])
          if (!is.finite(nb) || nb < 1L) {
            nb <- 1024L
          }
          sim <- tryCatch(
            .nts_fingerprint_similarity(ob_a, ob_c, nb, method = sim_method),
            error = function(e) NA_real_
          )
          sim <- suppressWarnings(as.numeric(sim)[1L])
          if (!is.finite(sim) || sim < mt) {
            next
          }
          if (
            is.null(best) ||
              .nts_propagation_pick_better(ms2_sid, sim, best$ms2_edge, best$similarity)
          ) {
            rid <- trimws(as.character(cand[["ID"]][k]))
            rsmi <- trimws(as.character(cand[["SMILES"]][k]))
            rcname <- .nts_struct_lib_row_compound_name(cand, k)
            if (!nzchar(rcname)) {
              rcname <- rid
            }
            best <- list(
              similarity = sim,
              pid = rid,
              psmi = rsmi,
              compound_name = rcname,
              seed_id = sid,
              ms2_edge = ms2_sid,
              seed_label = .nts_lineage_node_label(nodes, sid)
            )
          }
        }
      }

      if (!is.null(best)) {
        updates[[uid]] <- best
      }
    }

    hit_names <- names(updates)[!vapply(updates, is.null, logical(1L))]
    if (!length(hit_names)) {
      break
    }
    for (uid in hit_names) {
      hit <- updates[[uid]]
      pid_hit <- trimws(as.character(hit[["pid"]][1L]))
      smi_hit <- trimws(as.character(hit[["psmi"]][1L]))
      cname_hit <- trimws(as.character(hit[["compound_name"]][1L]))
      if (is.na(cname_hit) || !nzchar(cname_hit) || identical(cname_hit, pid_hit)) {
        alt <- .name_from_ref(pid_hit)
        cname_hit <- if (!is.na(alt) && nzchar(alt)) alt else NA_character_
      }
      sim_hit <- suppressWarnings(as.numeric(hit[["similarity"]][1L]))
      if (!is.finite(sim_hit)) {
        sim_hit <- NA_real_
      }
      seed_hit <- as.character(hit[["seed_label"]][1L])
      lab <- sprintf(
        "Network Propagated: %s (%s %.3f)",
        pid_hit,
        sim_lab,
        sim_hit
      )
      nodes[id_chr == uid, `:=`(
        propagated_ID = pid_hit,
        propagated_SMILES = smi_hit,
        propagated_similarity = sim_hit,
        propagated_from_node = seed_hit,
        network_propagation = lab,
        annotation_level = lab,
        Propagation_Level = as.integer(it),
        ref_id = pid_hit,
        Compound.name = cname_hit,
        SMILES = smi_hit
      )]
      if ("Library_ID" %in% names(nodes)) {
        nodes[id_chr == uid, Library_ID := pid_hit]
      }
      if ("Annotation" %in% names(nodes)) {
        nodes[id_chr == uid, Annotation := cname_hit]
      }
    }
  }

  net_data$nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  net_data
}

#' Pre-compile a CSV Compound Library into an Indexed RDS
#'
#' Reads a large CSV file with compound information, validates required columns,
#' sorts by ExactMass for binary-search indexing, and saves as a compact RDS.
#' This RDS can be passed directly to \code{annotate_neighbors_from_csv()} for
#' dramatically faster mass-guided lookups.
#'
#' @param csv_path Path to the CSV file (must contain ExactMass, SMILES, ID columns).
#' @param output_rds Path for the output RDS file. If NULL, uses \code{csv_path} with
#'        \code{.rds} extension.
#' @param smiles_col Column name for SMILES (auto-detected if NULL).
#' @param id_col Column name for compound ID (auto-detected if NULL).
#'
#' @return Invisible path to the created RDS file.
#' @export
precompile_csv_library <- function(csv_path,
                                   output_rds = NULL,
                                   smiles_col = NULL,
                                   id_col = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")

  csv_path <- path.expand(as.character(csv_path)[1L])
  if (!file.exists(csv_path)) stop("CSV file not found: ", csv_path)

  message(sprintf("Reading CSV: %s", basename(csv_path)))
  dt <- data.table::fread(csv_path, na.strings = c("", "NA", "N/A"))
  message(sprintf("  Loaded %d rows x %d columns.", nrow(dt), ncol(dt)))

  # Validate ExactMass
  if (!"ExactMass" %in% names(dt)) stop("CSV must contain an ExactMass column.")
  dt[, ExactMass := as.numeric(ExactMass)]
  dt <- dt[is.finite(ExactMass) & ExactMass > 0]

  # Validate SMILES
  if (is.null(smiles_col)) {
    smiles_col <- grep("^SMILES$|^QSAR_READY_SMILES$", names(dt),
                       ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(smiles_col)) stop("Could not find SMILES column.")
  }
  dt <- dt[!is.na(dt[[smiles_col]]) & nzchar(trimws(as.character(dt[[smiles_col]])))]

  # Validate ID
  if (is.null(id_col)) {
    id_col <- grep("^ID$|^DTXSID$|^Lab_ID$", names(dt),
                   ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(id_col)) id_col <- names(dt)[1L]
  }

  # Sort and key by ExactMass for binary search
  data.table::setorder(dt, ExactMass)
  data.table::setkey(dt, ExactMass)

  # Store column name metadata as attributes
  attr(dt, "smiles_col") <- smiles_col
  attr(dt, "id_col") <- id_col

  if (is.null(output_rds)) {
    output_rds <- sub("\\.[^.]+$", ".rds", csv_path)
  }
  saveRDS(dt, output_rds)
  message(sprintf("Indexed RDS saved: %s (%d compounds, %.1f MB)",
                  basename(output_rds), nrow(dt),
                  file.size(output_rds) / 1024^2))
  invisible(output_rds)
}

#' Targeted Structural Annotation via Mass-Guided BFS Propagation
#'
#' Performs BFS-layered structure propagation on a molecular network by:
#' (1) filtering an external compound library by MS1 exact mass (binary search),
#' (2) computing fingerprint similarity between seed SMILES and candidates,
#' (3) optionally computing Maximum Common Substructure (MCS) score with
#'     atom-count pre-screening,
#' (4) using MS2 edge weight + combined structural score for competition.
#'
#' @param net_data A list with `nodes` and `edges` (from `generate_network_data`).
#' @param library_path Path to CSV or pre-compiled RDS from `precompile_csv_library()`.
#' @param mz_tol_ppm Mass tolerance in ppm (default 10).
#' @param min_similarity Minimum fingerprint similarity threshold (default 0.5).
#' @param fingerprint_type rcdk fingerprint type (default "circular").
#' @param similarity_method "dice", "tanimoto", or "cosine" (default "dice").
#' @param use_mcs If TRUE, compute MCS score for candidates passing FP threshold.
#' @param mcs_fp_gate Minimum FP similarity to trigger MCS calculation (default 0.6).
#' @param polarity "neg" or "pos" (default "neg").
#' @param max_iterations Maximum BFS depth (default 3).
#' @param smiles_col Column name for SMILES (auto-detected if NULL).
#' @param id_col Column name for ID (auto-detected if NULL).
#'
#' @return Updated `net_data` with propagation columns on nodes.
#' @export
annotate_neighbors_from_csv <- function(net_data,
                                        library_path,
                                        mz_tol_ppm = 10,
                                        min_similarity = 0.5,
                                        fingerprint_type = "circular",
                                        similarity_method = "dice",
                                        use_mcs = TRUE,
                                        mcs_fp_gate = 0.4,
                                        polarity = "neg",
                                        max_iterations = 3L,
                                        smiles_col = NULL,
                                        id_col = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  .nts_rcdk_stop()

  # --- 1. Load compound library (CSV or pre-compiled RDS) ---
  library_path <- path.expand(as.character(library_path)[1L])
  if (!file.exists(library_path)) stop("Library file not found: ", library_path)

  is_rds <- grepl("\\.rds$", library_path, ignore.case = TRUE)

  if (is_rds) {
    message("Loading pre-compiled RDS library...")
    dt_list <- readRDS(library_path)
    if (!data.table::is.data.table(dt_list)) {
      dt_list <- data.table::as.data.table(dt_list)
    }
    # Recover column metadata from RDS attributes
    if (is.null(smiles_col)) smiles_col <- attr(dt_list, "smiles_col")
    if (is.null(id_col)) id_col <- attr(dt_list, "id_col")
  } else {
    message("Loading CSV library...")
    dt_list <- data.table::fread(library_path, na.strings = c("", "NA", "N/A"))
  }

  # Auto-detect columns if still NULL
  if (is.null(smiles_col) || is.na(smiles_col)) {
    smiles_col <- grep("^SMILES$|^QSAR_READY_SMILES$", names(dt_list),
                       ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(smiles_col)) stop("Could not find SMILES column in library.")
  }
  if (is.null(id_col) || is.na(id_col)) {
    id_col <- grep("^ID$|^DTXSID$|^Lab_ID$", names(dt_list),
                   ignore.case = TRUE, value = TRUE)[1L]
    if (is.na(id_col)) id_col <- names(dt_list)[1L]
  }
  if (!"ExactMass" %in% names(dt_list)) {
    stop("Library must contain an ExactMass column.")
  }

  # [P0] Ensure ExactMass is numeric and indexed
  dt_list[, ExactMass := as.numeric(ExactMass)]
  dt_list <- dt_list[is.finite(ExactMass) & ExactMass > 0]
  if (!identical(data.table::key(dt_list), "ExactMass")) {
    data.table::setorder(dt_list, ExactMass)
    data.table::setkey(dt_list, ExactMass)
  }

  # [P0] Mass-range pre-filter: only keep compounds within network mz range
  nodes_mz <- suppressWarnings(as.numeric(net_data$nodes$mz))
  nodes_mz <- nodes_mz[is.finite(nodes_mz) & nodes_mz > 0]
  if (length(nodes_mz)) {
    mass_lo <- min(nodes_mz) - 2
    mass_hi <- max(nodes_mz) + 2
    n_before <- nrow(dt_list)
    dt_list <- dt_list[ExactMass >= mass_lo & ExactMass <= mass_hi]
    message(sprintf("  Library: %d -> %d compounds (mass range %.0f-%.0f Da)",
                    n_before, nrow(dt_list), mass_lo, mass_hi))
  }

  message(sprintf("  Library ready: %d compounds with ExactMass index.", nrow(dt_list)))

  # --- 2. Prepare network tables ---
  nodes <- data.table::as.data.table(net_data$nodes)
  edges <- data.table::as.data.table(net_data$edges)
  nodes[, id_chr := as.character(id)]
  nodes <- .nts_network_nodes_ensure_identity_cols(nodes)

  # [P3] Set key on nodes for O(log n) lookups
  data.table::setkey(nodes, id_chr)

  # Build MS2 edge-weight lookup
  edge_ms2_lk <- .nts_edge_ms2_lookup_build(edges)

  # [P2] Build hash adjacency list for O(1) neighbor queries
  e1 <- edges[, .(a = as.character(from), b = as.character(to))]
  e2 <- edges[, .(a = as.character(to),   b = as.character(from))]
  eb <- unique(data.table::rbindlist(list(e1, e2)))

  adj_env <- new.env(hash = TRUE, parent = emptyenv())
  for (r in seq_len(nrow(eb))) {
    a_node <- eb$a[r]; b_node <- eb$b[r]
    if (exists(a_node, envir = adj_env, inherits = FALSE)) {
      assign(a_node, c(get(a_node, envir = adj_env, inherits = FALSE), b_node),
             envir = adj_env)
    } else {
      assign(a_node, b_node, envir = adj_env)
    }
  }

  pol <- tolower(trimws(as.character(polarity)[1L]))
  mass_offset <- if (startsWith(pol, "neg")) -1.007276 else 1.007276
  mt <- as.numeric(min_similarity)[1L]
  mcs_gate <- as.numeric(mcs_fp_gate)[1L]
  ppm_val <- as.numeric(mz_tol_ppm)[1L]
  mi <- suppressWarnings(as.integer(max_iterations)[1L])
  if (!is.finite(mi) || mi < 1L) mi <- 3L

  sim_lab <- paste0(toupper(substring(similarity_method, 1L, 1L)),
                    substring(similarity_method, 2L))

  # --- 3. Initialize propagation columns ---
  # Initialize with proper types to avoid data.table warnings
  if (!"propagated_ID" %in% names(nodes)) nodes[, propagated_ID := NA_character_]
  if (!"propagated_SMILES" %in% names(nodes)) nodes[, propagated_SMILES := NA_character_]
  if (!"propagated_similarity" %in% names(nodes)) nodes[, propagated_similarity := NA_real_]
  if (!"propagated_mcs_score" %in% names(nodes)) nodes[, propagated_mcs_score := NA_real_]
  if (!"propagated_from_node" %in% names(nodes)) nodes[, propagated_from_node := NA_character_]
  if (!"network_propagation" %in% names(nodes)) nodes[, network_propagation := NA_character_]
  if (!"Propagation_Level" %in% names(nodes)) nodes[, Propagation_Level := NA_integer_]

  # Detect seed nodes (layer 0)
  al <- as.character(nodes$annotation_level)
  al[is.na(al)] <- ""
  seed_mask <- grepl("Database Confirmed", al, fixed = FALSE)
  if ("Annotation" %in% names(nodes)) {
    a <- trimws(as.character(nodes$Annotation))
    seed_mask <- seed_mask | (!is.na(a) & nzchar(a))
  }
  if ("Library_ID" %in% names(nodes)) {
    lib <- trimws(as.character(nodes$Library_ID))
    seed_mask <- seed_mask | (!is.na(lib) & nzchar(lib))
  }
  if (!any(seed_mask)) {
    seed_mask <- (!is.na(nodes$SMILES) & nzchar(as.character(nodes$SMILES)))
  }
  if (!any(seed_mask)) {
    message("No seed nodes found. Returning net_data unchanged.")
    nodes[, id_chr := NULL]
    net_data$nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
    return(net_data)
  }

  nodes[, Propagation_Level := NA_integer_]
  nodes[seed_mask, Propagation_Level := 0L]

  # --- 4. Fingerprint & molecule cache ---
  fp_cache  <- new.env(hash = TRUE, parent = emptyenv())
  mol_cache <- new.env(hash = TRUE, parent = emptyenv())

  .get_cached_fp <- function(smi) {
    if (exists(smi, envir = fp_cache, inherits = FALSE)) {
      return(get(smi, envir = fp_cache, inherits = FALSE))
    }
    fp <- .nts_fp_info_from_smiles_rcdk(smi, fingerprint_type)
    assign(smi, fp, envir = fp_cache)
    fp
  }

  .get_cached_mol <- function(smi) {
    if (exists(smi, envir = mol_cache, inherits = FALSE)) {
      return(get(smi, envir = mol_cache, inherits = FALSE))
    }
    mol <- .nts_mol_first_from_smiles_rcdk(smi)
    assign(smi, mol, envir = mol_cache)
    mol
  }

  # --- 5. BFS propagation (optimized) ---
  for (it in seq_len(mi)) {
    frontier <- nodes[Propagation_Level == it - 1L, unique(id_chr)]
    frontier <- frontier[!is.na(frontier) & nzchar(frontier)]
    if (!length(frontier)) break

    # [P2] Use hash adj_env for fast target discovery
    cand_targets <- unique(unlist(lapply(frontier, function(f) {
      if (exists(f, envir = adj_env, inherits = FALSE)) {
        get(f, envir = adj_env, inherits = FALSE)
      }
    })))
    cand_targets <- cand_targets[!is.na(cand_targets) & nzchar(cand_targets)]
    open <- nodes[is.na(Propagation_Level), unique(id_chr)]
    targets <- intersect(cand_targets, open)
    if (!length(targets)) break

    message(sprintf("\n[BFS Layer %d] %d frontier seeds -> %d candidate targets",
                    it, length(frontier), length(targets)))
    pb <- utils::txtProgressBar(min = 0, max = length(targets), style = 3)

    updates <- vector("list", length(targets))
    names(updates) <- targets

    for (i in seq_along(targets)) {
      utils::setTxtProgressBar(pb, i)
      uid <- targets[i]

      # [P3] Keyed lookup on nodes
      u_mz <- nodes[uid, as.numeric(mz)][1L]
      if (!is.finite(u_mz) || u_mz <= 0) next
      target_neutral <- u_mz - mass_offset
      if (!is.finite(target_neutral) || target_neutral <= 0) next

      # [P0] Binary-search mass filter
      target_lower <- target_neutral * (1 - ppm_val * 1e-6)
      target_upper <- target_neutral * (1 + ppm_val * 1e-6)
      cand <- dt_list[ExactMass >= target_lower & ExactMass <= target_upper]
      if (!nrow(cand)) next

      # [P2] Hash-based neighbor lookup
      nei <- if (exists(uid, envir = adj_env, inherits = FALSE)) {
        get(uid, envir = adj_env, inherits = FALSE)
      } else character(0)
      seeds_conn <- intersect(nei, frontier)
      if (!length(seeds_conn)) next

      best <- NULL

      for (sid in seeds_conn) {
        # [P3] Keyed lookup
        s_smi <- nodes[sid, as.character(SMILES)][1L]
        if (is.na(s_smi) || !nzchar(s_smi)) next

        ms2_edge <- .nts_edge_ms2_lookup_score(edge_ms2_lk, uid, sid)

        s_fp <- .get_cached_fp(s_smi)
        if (is.null(s_fp)) next
        s_mol <- if (isTRUE(use_mcs)) .get_cached_mol(s_smi) else NULL

        for (k in seq_len(nrow(cand))) {
          c_smi <- cand[[smiles_col]][k]
          if (is.na(c_smi) || !nzchar(c_smi)) next

          c_fp <- .get_cached_fp(c_smi)
          if (is.null(c_fp)) next

          # Fingerprint similarity
          sim <- .nts_fingerprint_similarity(
            s_fp$onbits, c_fp$onbits, s_fp$fp_nbits, method = similarity_method
          )
          sim <- suppressWarnings(as.numeric(sim)[1L])
          if (!is.finite(sim) || sim < mt) next

          # [P1] MCS scoring with atom-count pre-screening
          mcs_score <- 0.0
          if (isTRUE(use_mcs) && !is.null(s_mol) && sim >= mcs_gate) {
            c_mol <- .get_cached_mol(c_smi)
            if (!is.null(c_mol)) {
              # Pre-screen: skip MCS if atom counts differ > 2x
              n_s <- tryCatch(rcdk::get.atom.count(s_mol), error = function(e) NA)
              n_c <- tryCatch(rcdk::get.atom.count(c_mol), error = function(e) NA)
              if (is.finite(n_s) && is.finite(n_c) &&
                  min(n_s, n_c) / max(n_s, n_c) >= 0.4) {
                mcs_score <- tryCatch({
                  mcs <- rcdk::get.mcs(s_mol, c_mol)
                  if (is.null(mcs)) 0 else {
                    n_mcs <- rcdk::get.atom.count(mcs)
                    n_mcs / max(n_s, n_c)
                  }
                }, error = function(e) 0)
              } else {
                mcs_score <- 0
              }
            }
          }

          # Combined structural score
          combined <- if (isTRUE(use_mcs) && is.finite(mcs_score)) {
            0.6 * sim + 0.4 * mcs_score
          } else {
            sim
          }

          # Competition: MS2 edge weight first, then combined structural score
          if (is.null(best) ||
              .nts_propagation_pick_better(ms2_edge, combined,
                                          best$ms2_edge, best$combined)) {
            rcname <- .nts_struct_lib_row_compound_name(cand, k)
            if (!nzchar(rcname)) {
              rcname <- trimws(as.character(cand[[id_col]][k]))
            }
            best <- list(
              pid = cand[[id_col]][k],
              psmi = c_smi,
              compound_name = rcname,
              similarity = sim,
              mcs = mcs_score,
              combined = combined,
              seed_id = sid,
              ms2_edge = ms2_edge,
              seed_label = .nts_lineage_node_label(nodes, sid)
            )
          }
        }
      }

      if (!is.null(best)) {
        updates[[uid]] <- best
      }
    }
    close(pb)

    # [P4] Batch apply updates
    hit_names <- names(updates)[!vapply(updates, is.null, logical(1L))]
    if (!length(hit_names)) break

    # Build update data.table for batch join
    upd_dt <- data.table::rbindlist(lapply(hit_names, function(uid) {
      hit <- updates[[uid]]
      pid_hit  <- trimws(as.character(hit$pid[1L]))
      smi_hit  <- trimws(as.character(hit$psmi[1L]))
      sim_hit  <- suppressWarnings(as.numeric(hit$similarity[1L]))
      mcs_hit  <- suppressWarnings(as.numeric(hit$mcs[1L]))
      seed_hit <- as.character(hit$seed_label[1L])

      mcs_text <- if (isTRUE(use_mcs) && is.finite(mcs_hit)) {
        sprintf("; MCS %.2f", mcs_hit)
      } else ""
      lab <- sprintf("Network Propagated: %s (%s %.3f%s)",
                     pid_hit, sim_lab, sim_hit, mcs_text)

      cname_hit <- trimws(as.character(hit$compound_name[1L]))
      if (is.na(cname_hit) || !nzchar(cname_hit) || identical(cname_hit, pid_hit)) {
        ## Unlike annotate_by_dual_network(), this CSV-neighbour path has no
        ## NTSProject whose @db_annotations could resolve a human-readable name
        ## (and the library row already had no name, which is why
        ## compound_name fell back to the bare id). Leave the field blank rather
        ## than stamping a ref_id into a name column; the caller treats NA as
        ## "no name known".
        cname_hit <- NA_character_
      }
      out_row <- data.table::data.table(
        id_chr = uid,
        propagated_ID = pid_hit,
        propagated_SMILES = smi_hit,
        propagated_similarity = sim_hit,
        propagated_mcs_score = mcs_hit,
        propagated_from_node = seed_hit,
        network_propagation = lab,
        annotation_level = lab,
        Propagation_Level = as.integer(it),
        ref_id = pid_hit,
        Compound.name = cname_hit,
        SMILES = smi_hit
      )
      if ("Library_ID" %in% names(nodes)) {
        out_row[, Library_ID := pid_hit]
      }
      if ("Annotation" %in% names(nodes)) {
        out_row[, Annotation := cname_hit]
      }
      out_row
    }))

    # Batch update via keyed join
    for (col in setdiff(names(upd_dt), "id_chr")) {
      nodes[upd_dt$id_chr, (col) := upd_dt[[col]]]
    }

    message(sprintf("[BFS Layer %d] Propagated %d nodes.", it, length(hit_names)))
  }

  # --- 6. Cleanup and return ---
  nodes[, id_chr := NULL]
  net_data$nodes <- as.data.frame(nodes, stringsAsFactors = FALSE)
  return(net_data)
}

#' Build a unified NTSDatabase from a MoNA export MSP file
#'
#' This function extracts both negative and positive mode spectra from a single MoNA
#' MSP file, cleans and filters the metadata (keeping only HRMS instruments and 
#' valid M-H/M+H precursors), normalizes and prunes noise from MS2 spectra, and 
#' deduplicates records using InChIKey or SMILES.
#'
#' @param msp_file Path to the MoNA MSP file.
#' @param noise_threshold Threshold for MS2 noise filtering.
#' @param top_k Max peaks to keep per spectrum.
#' @param ... Additional arguments passed to \code{build_nts_database}.
#'
#' @return An \code{NTSDatabase} object.
#' @export
build_mona_db <- function(msp_file, noise_threshold = 0.01, top_k = 30L, ...) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  
  msp_file <- path.expand(as.character(msp_file)[1L])
  .nts_check(msp_file, "file_exists")
  
  message("build_mona_db: Extracting Negative mode data...")
  res_neg <- suppressWarnings(import_msp(msp_file, polarity = "NEG"))
  
  message("build_mona_db: Extracting Positive mode data...")
  res_pos <- suppressWarnings(import_msp(msp_file, polarity = "POS"))
  
  all_metadata <- data.table::rbindlist(
    list(res_neg$metadata, res_pos$metadata),
    fill = TRUE,
    use.names = TRUE
  )
  
  dt_spec_neg <- data.table::as.data.table(res_neg$spectra_long)
  if (nrow(dt_spec_neg) > 0) dt_spec_neg[, Polarity := "NEG"]
  
  dt_spec_pos <- data.table::as.data.table(res_pos$spectra_long)
  if (nrow(dt_spec_pos) > 0) dt_spec_pos[, Polarity := "POS"]
  
  all_spectra <- data.table::rbindlist(list(dt_spec_neg, dt_spec_pos), fill = TRUE, use.names = TRUE)
  
  if (!nrow(all_metadata) || !nrow(all_spectra)) {
    stop("build_mona_db: No valid spectra extracted from MSP file.", call. = FALSE)
  }
  
  # Filter SMILES
  if (!"SMILES" %in% names(all_metadata)) {
      all_metadata[, SMILES := NA_character_]
  }
  all_metadata <- all_metadata[!is.na(SMILES) & nzchar(trimws(as.character(SMILES)))]
  
  if (!nrow(all_metadata)) stop("build_mona_db: No records with valid SMILES remaining.", call. = FALSE)
  
  # Filter Precursor Type
  m_minus_h_variants <- c("[M-H]-", "M-H", "[M-H]", "M-1", "[M-H]1-", "deprotonated molecule", "Deprotonated molecule")
  m_plus_h_variants  <- c("[M+H]+", "M+H", "[M+H]", "M+1", "(M+H)+", "Protonated")
  
  if ("Precursor_type" %in% names(all_metadata)) {
    all_metadata <- all_metadata[Precursor_type %in% c(m_minus_h_variants, m_plus_h_variants)]
    all_metadata[Precursor_type %in% m_minus_h_variants, Precursor_type := "[M-H]-"]
    all_metadata[Precursor_type %in% m_plus_h_variants, Precursor_type := "[M+H]+"]
  }
  
  if (!nrow(all_metadata)) stop("build_mona_db: No records matching standard precursor types.", call. = FALSE)
  
  # Filter HRMS (only when instrument metadata is present on at least one row;
  # curated MSPs may omit Instrument_type entirely - do not drop those records.)
  if ("Instrument_type" %in% names(all_metadata)) {
    it_chr <- trimws(as.character(all_metadata$Instrument_type))
    has_instr <- !is.na(it_chr) & nzchar(it_chr)
    if (any(has_instr)) {
      hrms_pattern <- "(?i)TOF|Orbitrap|FT|Exactive|Tribrid|Maxis|Impact"
      keep_rows <- !has_instr | grepl(hrms_pattern, it_chr, perl = TRUE)
      all_metadata <- all_metadata[keep_rows]
    }
  }
  
  if (!nrow(all_metadata)) stop("build_mona_db: No records matching HRMS instrument patterns.", call. = FALSE)
  
  # Sync spectra
  all_spectra <- all_spectra[Lab_ID %in% all_metadata$Lab_ID]
  
  # Normalize and denoise spectra
  all_spectra[, intensity := (as.numeric(intensity) / max(as.numeric(intensity), na.rm = TRUE)) * 1000, by = .(Lab_ID)]
  all_spectra <- all_spectra[intensity >= 10]
  
  valid_labs <- unique(all_spectra$Lab_ID)
  all_metadata <- all_metadata[Lab_ID %in% valid_labs]
  
  if (!nrow(all_metadata)) stop("build_mona_db: No valid spectra remaining after denoising.", call. = FALSE)
  
  # Deduplicate
  if (!"InChIKey" %in% names(all_metadata)) all_metadata[, InChIKey := NA_character_]
  
  all_metadata[, compound_key := ifelse(!is.na(InChIKey) & nzchar(trimws(as.character(InChIKey))), as.character(InChIKey), as.character(SMILES))]
  
  unique_compounds <- all_metadata[, .(
    new_Lab_ID = paste0("MONA_C", sprintf("%06d", .GRP))
  ), by = .(compound_key)]
  
  all_metadata <- merge(all_metadata, unique_compounds, by = "compound_key", all.x = TRUE)
  
  get_first_valid <- function(x) {
    valid_x <- stats::na.omit(x[nzchar(trimws(as.character(x)))])
    if (length(valid_x) > 0) return(as.character(valid_x)[1L]) else return(NA_character_)
  }
  get_all_unique <- function(x) {
    valid_x <- unique(stats::na.omit(x[nzchar(trimws(as.character(x)))]))
    if (length(valid_x) > 0) return(paste(valid_x, collapse = " | ")) else return(NA_character_)
  }
  
  for (cname in c("CAS", "PubChem_CID", "Instrument_type", "name", "Formula")) {
    if (!cname %in% names(all_metadata)) all_metadata[, (cname) := NA_character_]
  }
  
  clean_metadata <- all_metadata[, .(
    name = get_all_unique(name),
    mass = round(mean(as.numeric(mass), na.rm = TRUE), 4),
    SMILES = get_first_valid(SMILES),
    InChIKey = get_first_valid(InChIKey),
    Formula = get_first_valid(Formula),
    Polarity = paste(sort(unique(as.character(Polarity))), collapse = "/"),
    CAS = get_all_unique(CAS),
    PubChem_CID = get_all_unique(PubChem_CID),
    Instrument_type = get_all_unique(Instrument_type),
    target_mz = as.numeric(stats::na.omit(as.numeric(target_mz))[1L])
  ), by = .(new_Lab_ID)]
  
  data.table::setnames(clean_metadata, "new_Lab_ID", "Lab_ID")
  # MoNA is a public/online spectral library -> provenance tag "online".
  clean_metadata[, db_type := "online"]

  id_map <- unique(all_metadata[, .(Lab_ID, new_Lab_ID)])
  all_spectra <- merge(all_spectra, id_map, by = "Lab_ID", all.x = TRUE)
  all_spectra <- all_spectra[!is.na(new_Lab_ID)]
  all_spectra[, Lab_ID := new_Lab_ID]
  all_spectra[, new_Lab_ID := NULL]
  
  build_nts_database(
    metadata = clean_metadata,
    spectra_long = all_spectra,
    noise_threshold = noise_threshold,
    top_k = top_k,
    ...
  )
}

