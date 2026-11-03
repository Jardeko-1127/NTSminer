# GC-QTOF data I/O - EI spectrum extraction, TIC reading, PCDL import ---------

#' Extract GC-EI full-scan spectra at feature apex RT
#'
#' **Deprecated**: use [detect_gc_compounds] instead for erah-style
#' full-chromatogram deconvolution that directly outputs compound spectra
#' without requiring prior CentWave peak picking.
#'
#' For each feature in an NTSproject that has been processed with
#' \code{extract_features()}, opens the raw mzML files and extracts the
#' full-scan EI spectrum at the feature apex retention time.
#'
#' @param project An [NTSProject] containing \code{@xdata} (XCMSnExp)
#'   and \code{@feature_table}.
#' @param method Character, method to select scans. Currently only
#'   \code{"apex"} is supported: uses the apex RT from chromPeaks.
#' @param merge_scans Integer, number of scans on each side of the apex to
#'   merge (takes the maximum intensity per m/z across scans). Default 1
#'   (apex scan only).
#' @param bg_subtract Logical, whether to perform background subtraction
#'   using scans at the feature boundaries (rtmin/rtmax). Default TRUE.
#' @param BPPARAM A \code{BiocParallelParam} object for parallel processing.
#'   Default \code{BiocParallel::SerialParam()}.
#'
#' @return The project with \code{@gc_spectra} populated: a named list
#'   where names correspond to \code{feature_table$name} and each element
#'   is a \code{data.frame} with columns \code{mz} and \code{intensity}.
#'
#' @seealso [extract_features], [extract_ms2]
#' @importFrom mzR openMSfile header spectra close
#' @export
extract_gc_spectra <- function(project,
                               method = c("apex"),
                               merge_scans = 1L,
                               bg_subtract = TRUE,
                               BPPARAM = BiocParallel::SerialParam()) {
  .nts_check(project, "project", .arg = "project")
  if (!requireNamespace("mzR", quietly = TRUE)) stop("mzR is required", call. = FALSE)
  if (!requireNamespace("xcms", quietly = TRUE)) stop("xcms is required", call. = FALSE)

  method <- match.arg(method)
  merge_scans <- max(1L, as.integer(merge_scans)[1L])

  ft <- get_feature_table(project)
  if (!nrow(ft)) {
    warning("feature_table is empty; nothing to extract", call. = FALSE)
    project@gc_spectra <- list()
    return(project)
  }
  if (!all(c("name", "mz", "rt") %in% names(ft))) {
    stop("feature_table must contain columns: name, mz, rt", call. = FALSE)
  }

  xd <- project@xdata
  if (!inherits(xd, "XCMSnExp")) {
    stop("@xdata must be an XCMSnExp from extract_features()", call. = FALSE)
  }

  cpeaks <- xcms::chromPeaks(xd)
  fnames  <- xd@processingData@files

  if (!nrow(cpeaks)) {
    warning("No chromPeaks found in @xdata; using feature_table mz/rt directly", call. = FALSE)
    cpeaks <- NULL
  }

  ft_names <- as.character(ft$name)
  n <- length(ft_names)

  .extract_one <- function(i) {
    nm <- ft_names[i]
    ft_mz <- as.numeric(ft$mz[i])
    ft_rt <- as.numeric(ft$rt[i])

    # Determine target RT and file from chromPeaks (best intensity row)
    target_rt  <- ft_rt
    target_rtmin <- target_rt - 5
    target_rtmax <- target_rt + 5
    target_file <- fnames[1L]

    if (!is.null(cpeaks) && nrow(cpeaks)) {
      cp_idx <- integer(0L)
      # Match by m/z proximity first
      if ("mz" %in% colnames(cpeaks)) {
        mz_diff <- abs(cpeaks[, "mz"] - ft_mz)
        cp_idx  <- which(mz_diff < 0.02 & cpeaks[, "sample"] %in% seq_along(fnames))
      } else {
        cp_idx <- which(cpeaks[, "sample"] %in% seq_along(fnames))
      }
      if (length(cp_idx)) {
        cp_sub <- cpeaks[cp_idx, , drop = FALSE]
        if ("into" %in% colnames(cp_sub)) {
          best <- which.max(cp_sub[, "into"])
        } else {
          best <- 1L
        }
        if (length(best) && is.finite(best)) {
          target_rt  <- as.numeric(cp_sub[best, "rt"])
          target_rtmin <- as.numeric(cp_sub[best, "rtmin"])
          target_rtmax <- as.numeric(cp_sub[best, "rtmax"])
          si <- as.integer(cp_sub[best, "sample"])
          if (si >= 1L && si <= length(fnames)) target_file <- fnames[si]
        }
      }
    }

    if (!nzchar(target_file) || !file.exists(target_file)) {
      warning("File not found for feature ", nm, call. = FALSE)
      return(data.frame(mz = numeric(0L), intensity = numeric(0L)))
    }

    ms <- try(mzR::openMSfile(target_file), silent = TRUE)
    if (inherits(ms, "try-error")) {
      warning("Cannot open mzML for feature ", nm, call. = FALSE)
      return(data.frame(mz = numeric(0L), intensity = numeric(0L)))
    }
    on.exit(try(mzR::close(ms), silent = TRUE), add = TRUE)

    hdr <- try(mzR::header(ms), silent = TRUE)
    if (inherits(hdr, "try-error") || !is.data.frame(hdr) || !nrow(hdr)) {
      warning("Cannot read header for feature ", nm, call. = FALSE)
      return(data.frame(mz = numeric(0L), intensity = numeric(0L)))
    }

    rt_hdr <- as.numeric(hdr$retentionTime)
    n_scans <- length(rt_hdr)

    # Find apex scan
    apex_idx <- which.min(abs(rt_hdr - target_rt))
    if (!length(apex_idx)) {
      return(data.frame(mz = numeric(0L), intensity = numeric(0L)))
    }

    # Collect scans to merge
    scan_start <- max(1L, apex_idx - merge_scans)
    scan_end   <- min(n_scans, apex_idx + merge_scans)
    scan_seq   <- scan_start:scan_end

    spec_list <- lapply(scan_seq, function(si) {
      sp <- try(mzR::spectra(ms, si)[[1L]], silent = TRUE)
      if (inherits(sp, "try-error") || !is.matrix(sp)) return(NULL)
      data.frame(mz = as.numeric(sp[, "mz"]), intensity = as.numeric(sp[, "intensity"]))
    })
    spec_list <- spec_list[!vapply(spec_list, is.null, logical(1L))]
    if (!length(spec_list)) {
      return(data.frame(mz = numeric(0L), intensity = numeric(0L)))
    }

    # Merge: combine all peaks, bin by m/z (0.01 Da), take max intensity
    spec_all <- do.call(rbind, spec_list)
    spec_all$mz_bin <- round(spec_all$mz, 2L)
    spec_merged <- stats::aggregate(
      intensity ~ mz_bin, data = spec_all, FUN = max
    )
    names(spec_merged)[names(spec_merged) == "mz_bin"] <- "mz"
    spec_merged <- spec_merged[order(spec_merged$mz), ]

    # Background subtraction
    if (bg_subtract && nrow(spec_merged)) {
      bg_spectra <- list()
      # Left background: scan just before rtmin
      if (is.finite(target_rtmin)) {
        left_idx <- which(rt_hdr < target_rtmin)
        if (length(left_idx)) {
          left_scan <- left_idx[length(left_idx)]
          sp <- try(mzR::spectra(ms, left_scan)[[1L]], silent = TRUE)
          if (!inherits(sp, "try-error") && is.matrix(sp) && nrow(sp)) {
            bg_spectra[[length(bg_spectra) + 1L]] <- sp
          }
        }
      }
      # Right background: scan just after rtmax
      if (is.finite(target_rtmax)) {
        right_idx <- which(rt_hdr > target_rtmax)
        if (length(right_idx)) {
          right_scan <- right_idx[1L]
          sp <- try(mzR::spectra(ms, right_scan)[[1L]], silent = TRUE)
          if (!inherits(sp, "try-error") && is.matrix(sp) && nrow(sp)) {
            bg_spectra[[length(bg_spectra) + 1L]] <- sp
          }
        }
      }
      if (length(bg_spectra)) {
        bg_all <- do.call(rbind, lapply(bg_spectra, function(m) {
          data.frame(mz = as.numeric(m[, "mz"]), intensity = as.numeric(m[, "intensity"]))
        }))
        bg_all$mz_bin <- round(bg_all$mz, 2L)
        bg_avg <- stats::aggregate(intensity ~ mz_bin, data = bg_all, FUN = mean)
        names(bg_avg)[names(bg_avg) == "mz_bin"] <- "mz"
        mz_match <- match(round(spec_merged$mz, 2L), round(bg_avg$mz, 2L))
        spec_merged$bg_intensity <- bg_avg$intensity[mz_match]
        spec_merged$bg_intensity[is.na(spec_merged$bg_intensity)] <- 0
        spec_merged$intensity <- spec_merged$intensity - spec_merged$bg_intensity
        spec_merged$intensity[spec_merged$intensity < 0] <- 0
        spec_merged$bg_intensity <- NULL
      }
    }

    spec_merged[, c("mz", "intensity")]
  }

  .nts_assert_parallel_ntsminer()
  slist <- BiocParallel::bplapply(
    seq_len(n), .extract_one, BPPARAM = BPPARAM
  )
  names(slist) <- ft_names

  project@gc_spectra <- slist
  validObject(project)
  project
}


#' Read total ion chromatogram from GC-QTOF mzML files
#'
#' Extracts the total ion current (TIC) from mzML file headers
#' without reading full spectral data.
#'
#' @param mzml_files Character vector of mzML file paths.
#'
#' @return A `data.table` with columns: `filename`, `rt` (seconds),
#'   `intensity` (total ion current).
#'
#' @importFrom data.table data.table rbindlist
#' @export
read_gc_tic <- function(mzml_files) {
  if (!requireNamespace("mzR", quietly = TRUE)) stop("mzR is required", call. = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)

  mzml_files <- as.character(mzml_files)
  mzml_files <- mzml_files[file.exists(mzml_files)]
  if (!length(mzml_files)) stop("No valid mzML files provided", call. = FALSE)

  tic_list <- lapply(mzml_files, function(f) {
    ms <- mzR::openMSfile(f)
    on.exit(try(mzR::close(ms), silent = TRUE))
    hdr <- mzR::header(ms)
    if (!is.data.frame(hdr) || !nrow(hdr)) {
      return(data.table::data.table(
        filename = basename(f), rt = numeric(0L), intensity = numeric(0L)
      ))
    }
    data.table::data.table(
      filename  = basename(f),
      rt        = as.numeric(hdr$retentionTime),
      intensity = as.numeric(hdr$totIonCurrent)
    )
  })

  data.table::rbindlist(tic_list)
}


# PCDL import and GC library building ------------------------------------------

#' Import an Agilent PCDL database for GC-QTOF spectral matching
#'
#' PCDL (Personal Compound Database and Library) files are Agilent's
#' proprietary format (SQL Server Compact 3.5, encrypted with password
#' \code{Stevens_Creek_95051}). This function attempts to convert them
#' to MSP format usable by [import_msp] and [build_gc_library].
#'
#' @param pcdl_path Character, path to a \code{.cdb} PCDL file.
#' @param export_msp Logical, if TRUE attempt to use Agilent tools to
#'   export to MSP format.
#' @param msp_dir Character, directory to save exported MSP files.
#'   If NULL, uses a temporary directory.
#'
#' @return A character vector of paths to the exported MSP files,
#'   or an error if conversion is not possible.
#'
#' @section PCDL Conversion Strategies:
#' \enumerate{
#'   \item Agilent MassHunter PCDL Manager CLI (if \code{PCDLManager.exe} is
#'     found on \code{PATH}, in the usual install locations, or via the
#'     \code{NTS_PCDLMANAGER} environment variable): runs \code{--export-msp}
#'     against the cdb file.
#'   \item SQLite direct read (PCDL v2+): attempts \code{DBI/RSQLite}
#'     connection. Not yet fully implemented.
#'   \item External extraction (recommended for Windows): use PowerShell
#'     with \code{System.Data.SqlServerCe} (\code{ai_tasks/extract_pcdl.ps1})
#'     or Python with \code{pythonnet/clr} (\code{ai_tasks/extract_pcdl.py}).
#'     This reads the \code{Compounds} and \code{Spectra} tables directly
#'     and writes CSV + MSP output.
#'   \item Agilent PCDL Manager GUI: \code{File > Export > To CEF} then
#'     parse with [import_cef]; or \code{File > Export > Compounds/Spectra}.
#' }
#'
#' @seealso [import_cef], [build_gc_library], [import_msp]
#' @export
import_pcdl <- function(pcdl_path, export_msp = TRUE, msp_dir = NULL) {
  pcdl_path <- path.expand(as.character(pcdl_path)[1L])

  # Strategy 1: Agilent PCDL Manager CLI export
  if (export_msp) {
    exe <- .nts_resolve_executable(
      "PCDLManager",
      candidates = c(
        "C:/Program Files/Agilent/MassHunter/PCDLManager/PCDLManager.exe",
        "C:/Program Files (x86)/Agilent/MassHunter/PCDLManager/PCDLManager.exe",
        "C:/Agilent/MassHunter/PCDLManager/PCDLManager.exe"
      ),
      envvar = "NTS_PCDLMANAGER"
    )
    if (!is.null(exe)) {
      msp_dir <- if (is.null(msp_dir)) tempfile("pcdl_msp_") else msp_dir
      dir.create(msp_dir, showWarnings = FALSE, recursive = TRUE)
      result <- try(system2(
        exe,
        args = c("--export-msp", shQuote(pcdl_path), shQuote(msp_dir)),
        stdout = TRUE, stderr = TRUE
      ), silent = TRUE)
      msp_files <- list.files(msp_dir, pattern = "\\.msp$", full.names = TRUE, ignore.case = TRUE)
      if (length(msp_files)) {
        message("PCDL exported ", length(msp_files), " MSP file(s) to ", msp_dir)
        return(msp_files)
      }
    }
  }

  # Strategy 2: SQLite direct read
  if (requireNamespace("DBI", quietly = TRUE) && requireNamespace("RSQLite", quietly = TRUE)) {
    con <- try(DBI::dbConnect(RSQLite::SQLite(), pcdl_path), silent = TRUE)
    if (!inherits(con, "try-error")) {
      on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)
      tables <- try(DBI::dbListTables(con), silent = TRUE)
      if (!inherits(tables, "try-error") && length(tables)) {
        message("PCDL SQLite tables found: ", paste(tables, collapse = ", "))
        message("SQLite direct read not yet implemented. Please export manually.")
      }
    }
  }

  stop(
    "Cannot parse PCDL directly. ",
    "Please export to MSP using Agilent MassHunter Qualitative:\n",
    "  File > Export > MSP format\n",
    "Then use build_gc_library() on the exported .msp file.",
    call. = FALSE
  )
}


#' Build an NTSDatabase-compatible GC-EI spectral library from MSP files
#'
#' Convenience wrapper that reads one or more MSP files and assembles
#' them into an [NTSDatabase] object for use with [score_ei_spectra].
#'
#' @param msp_files Character vector of paths to MSP-format spectral
#'   library files.
#' @param polarity Character, "EI" by default (GC-EI has no polarity).
#' @param ... Additional arguments passed to [import_msp].
#'
#' @return An [NTSDatabase] object.
#'
#' @seealso [import_msp], [NTSDatabase], [score_ei_spectra]
#' @export
build_gc_library <- function(msp_files, polarity = "EI", ...) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)

  msp_files <- as.character(msp_files)
  msp_files <- msp_files[file.exists(msp_files)]
  if (!length(msp_files)) stop("No valid MSP files provided", call. = FALSE)

  pol_msp <- toupper(polarity[1L])

  meta_list <- list()
  spec_list <- list()

  for (f in msp_files) {
    res <- import_msp(f, polarity = pol_msp, ...)
    if (nrow(res$metadata)) {
      meta_list[[length(meta_list) + 1L]] <- res$metadata
      spec_list[[length(spec_list) + 1L]] <- res$spectra_long
    }
  }

  if (!length(meta_list)) {
    stop("No valid spectra imported from MSP files", call. = FALSE)
  }

  meta_all <- data.table::rbindlist(meta_list, fill = TRUE)
  spec_all <- data.table::rbindlist(spec_list, fill = TRUE)

  # Ensure required columns
  if (!"Lab_ID" %in% names(meta_all)) {
    meta_all[, Lab_ID := as.character(seq_len(.N))]
  }
  if (!"Lab_ID" %in% names(spec_all) && nrow(spec_all)) {
    stop("spectra_long must contain Lab_ID column", call. = FALSE)
  }
  # Map metadata Lab_ID to spectra_long lab_id
  data.table::setnames(meta_all, "Lab_ID", "lab_id", skip_absent = TRUE)
  data.table::setnames(spec_all, "Lab_ID", "lab_id", skip_absent = TRUE)

  methods::new(
    "NTSDatabase",
    metadata      = meta_all,
    spectra_long  = spec_all,
    features_long = data.table::data.table(),
    version       = "0.1.0"
  )
}


#' Import an Agilent CEF (Compound Exchange Format) XML file
#'
#' Parses Agilent's CEF XML export format (from PCDL Manager GUI:
#' \code{File > Export > To CEF}) and returns a list with
#' \code{metadata} and \code{spectra_long} \code{data.table}s
#' suitable for constructing an [NTSDatabase] via
#' \code{methods::new("NTSDatabase", ...)}.
#'
#' CEF is Agilent's XML-based compound exchange format containing
#' compound names, CAS numbers, formulas, exact masses, retention
#' times, and full spectral peak lists.
#'
#' @param cef_path Character, path to a \code{.cef} XML file.
#'
#' @return A list with two \code{data.table}s: \code{metadata}
#'   (columns: \code{lab_id, name, cas, formula, exactmass, rt})
#'   and \code{spectra_long} (columns: \code{lab_id, mz, intensity}).
#'
#' @section CEF Format:
#' The CEF XML structure expected:
#' \preformatted{
#'   <CEF>
#'     <CompoundList>
#'       <Compound>
#'         <Name>...</Name>
#'         <CAS>...</CAS>
#'         <Formula>...</Formula>
#'         <ExactMass>...</ExactMass>
#'         <RT>...</RT>
#'         <Spectrum>
#'           <MSPeaks>
#'             <p x="mz" y="intensity"/>
#'           </MSPeaks>
#'         </Spectrum>
#'       </Compound>
#'     </CompoundList>
#'   </CEF>
#' }
#'
#' @seealso [import_pcdl], [import_msp], [build_gc_library]
#' @importFrom xml2 read_xml xml_find_all xml_find_first xml_text xml_attr xml_ns
#' @export
import_cef <- function(cef_path) {
  if (!requireNamespace("xml2", quietly = TRUE))
    stop("xml2 is required. install.packages('xml2')", call. = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table is required", call. = FALSE)

  cef_path <- path.expand(as.character(cef_path)[1L])
  .nts_check(cef_path, "file_exists", .arg = "cef_path")

  message("import_cef: parsing ", cef_path)
  doc <- xml2::read_xml(cef_path)
  ns <- xml2::xml_ns(doc)

  compounds <- xml2::xml_find_all(doc, "//d1:Compound", ns)
  if (!length(compounds))
    compounds <- xml2::xml_find_all(doc, "//Compound")
  message("import_cef: found ", length(compounds), " compounds")

  .txt <- function(node, xpath) {
    child <- xml2::xml_find_first(node, xpath, ns)
    if (length(child) && !is.na(child)) {
      val <- xml2::xml_text(child)
      if (nchar(val) > 0) val else NA_character_
    } else NA_character_
  }

  meta_list <- vector("list", length(compounds))
  spec_list <- vector("list", length(compounds))

  for (i in seq_along(compounds)) {
    cmp <- compounds[[i]]

    name_val    <- .txt(cmp, "d1:Name")
    cas_val     <- .txt(cmp, "d1:CAS")
    formula_val <- .txt(cmp, "d1:Formula")
    emass_val   <- .txt(cmp, "d1:ExactMass")
    rt_val      <- .txt(cmp, "d1:RT")
    mw_val      <- .txt(cmp, "d1:MolecularWeight")

    lab_id <- as.character(i)
    rt_num <- suppressWarnings(as.numeric(rt_val))
    if (is.na(rt_num)) rt_num <- NA_real_
    emass_num <- suppressWarnings(as.numeric(emass_val))
    if (is.na(emass_num)) {
      mw_num <- suppressWarnings(as.numeric(mw_val))
      if (!is.na(mw_num)) emass_num <- mw_num
    }

    meta_list[[i]] <- data.table::data.table(
      lab_id    = lab_id,
      name      = if (is.na(name_val)) paste0("Compound_", i) else name_val,
      cas       = if (is.na(cas_val)) NA_character_ else cas_val,
      formula   = if (is.na(formula_val)) NA_character_ else formula_val,
      exactmass = emass_num,
      rt        = rt_num
    )

    peaks <- xml2::xml_find_all(cmp, ".//d1:p", ns)
    if (!length(peaks))
      peaks <- xml2::xml_find_all(cmp, ".//p")

    if (length(peaks)) {
      mz_vec  <- numeric(length(peaks))
      int_vec <- numeric(length(peaks))
      for (j in seq_along(peaks)) {
        mz_vec[j]  <- suppressWarnings(as.numeric(xml2::xml_attr(peaks[[j]], "x")))
        int_vec[j] <- suppressWarnings(as.numeric(xml2::xml_attr(peaks[[j]], "y")))
      }
      keep <- !is.na(mz_vec) & !is.na(int_vec)
      if (any(keep)) {
        spec_list[[i]] <- data.table::data.table(
          lab_id    = lab_id,
          mz        = mz_vec[keep],
          intensity = int_vec[keep]
        )
      }
    }
  }

  meta <- data.table::rbindlist(meta_list, fill = TRUE)
  spec <- data.table::rbindlist(
    spec_list[!vapply(spec_list, is.null, logical(1L))], fill = TRUE
  )
  if (!nrow(spec)) {
    spec <- data.table::data.table(
      lab_id = character(), mz = numeric(), intensity = numeric()
    )
  }

  message("import_cef: ", nrow(meta), " compounds, ",
          length(unique(spec$lab_id)), " spectra, ",
          nrow(spec), " peaks")
  list(metadata = meta, spectra_long = spec)
}


# GC-EI deconvolution - Biller-Biemann algorithm ---------------------------------

#' Deconvolve GC-EI co-eluting spectra using Biller-Biemann algorithm
#'
#' For each feature in @gc_spectra, reads the raw scan data in the
#' chromatographic window and separates co-eluting components by
#' clustering m/z channels with shared elution apexes.
#'
#' @param project An [NTSProject] with populated \code{@xdata} and
#'   \code{@gc_spectra}.
#' @param min_snr Numeric, minimum signal-to-noise ratio for an EIC
#'   apex to be considered a valid ion. Default 3.
#' @param apex_tol Numeric, RT tolerance in seconds for grouping m/z
#'   channels into the same component. Default 1.5.
#' @param min_ions Integer, minimum number of m/z channels required
#'   for a valid deconvolved component. Default 5.
#' @param edge_trim Numeric, RT margin in seconds to trim from the
#'   window edges (avoids boundary artifacts). Default 2.
#' @param BPPARAM A \code{BiocParallelParam} object. Default serial.
#'
#' @return The project with \code{@gc_deconvolved} populated: a named
#'   list where each element is a list of component lists. Each
#'   component has \code{$spectrum} (data.frame: mz, intensity),
#'   \code{$apex_rt} (numeric), and \code{$n_ions} (integer).
#'
#' @references Biller, J. E., & Biemann, K. (1974). Reconstructed
#'   mass spectra, a novel approach for the utilization of gas
#'   chromatograph-mass spectrometer data. Analytical Letters, 7(7),
#'   515-528.
#'
#' @seealso [extract_gc_spectra], [score_ei_spectra]
#' @importFrom mzR openMSfile header peaks close
#' @export
deconvolve_gc_spectra <- function(project,
                                  min_snr    = 3,
                                  apex_tol   = 1.5,
                                  min_ions   = 5L,
                                  edge_trim  = 2,
                                  BPPARAM    = BiocParallel::SerialParam()) {
  .nts_check(project, "project", .arg = "project")
  if (!length(project@gc_spectra)) stop("@gc_spectra is empty; run extract_gc_spectra() first", call. = FALSE)
  if (!requireNamespace("mzR", quietly = TRUE)) stop("mzR is required", call. = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)

  xd <- project@xdata
  if (!inherits(xd, "XCMSnExp")) stop("@xdata must be an XCMSnExp from extract_features()", call. = FALSE)

  cpeaks <- xcms::chromPeaks(xd)
  fnames <- xd@processingData@files
  ft <- get_feature_table(project)
  ft_names <- as.character(ft$name)
  min_ions <- max(1L, as.integer(min_ions)[1L])
  n <- length(ft_names)

  .deconv_one <- function(i) {
    nm <- ft_names[i]
    ft_mz <- as.numeric(ft$mz[i])
    ft_rt <- as.numeric(ft$rt[i])

    rtmin <- ft_rt - 5
    rtmax <- ft_rt + 5
    target_file <- fnames[1L]

    if (!is.null(cpeaks) && nrow(cpeaks)) {
      mz_diff <- abs(cpeaks[, "mz"] - ft_mz)
      cp_idx <- which(mz_diff < 0.02 & cpeaks[, "sample"] %in% seq_along(fnames))
      if (length(cp_idx)) {
        cp_sub <- cpeaks[cp_idx, , drop = FALSE]
        if ("into" %in% colnames(cp_sub)) {
          best <- which.max(cp_sub[, "into"])
        } else {
          best <- 1L
        }
        if (length(best) && is.finite(best)) {
          rtmin <- as.numeric(cp_sub[best, "rtmin"])
          rtmax <- as.numeric(cp_sub[best, "rtmax"])
          si <- as.integer(cp_sub[best, "sample"])
          if (si >= 1L && si <= length(fnames)) target_file <- fnames[si]
        }
      }
    }

    rtmin <- rtmin + edge_trim
    rtmax <- rtmax - edge_trim
    if (rtmin >= rtmax) return(list())

    ms <- try(mzR::openMSfile(target_file), silent = TRUE)
    if (inherits(ms, "try-error")) return(list())
    on.exit(try(mzR::close(ms), silent = TRUE), add = TRUE)

    hdr <- try(mzR::header(ms), silent = TRUE)
    if (inherits(hdr, "try-error")) return(list())

    scan_rt <- as.numeric(hdr$retentionTime)
    win_idx <- which(scan_rt >= rtmin & scan_rt <= rtmax)
    if (length(win_idx) < 5) return(list())

    peak_list <- lapply(win_idx, function(si) {
      sp <- try(mzR::peaks(ms, si), silent = TRUE)
      if (inherits(sp, "try-error") || !is.matrix(sp) || !nrow(sp)) return(NULL)
      data.table::data.table(
        scan_idx = si,
        mz_bin   = round(as.numeric(sp[, "mz"]), 1L),
        intensity = as.numeric(sp[, "intensity"])
      )
    })
    peak_list <- peak_list[!vapply(peak_list, is.null, logical(1L))]
    if (!length(peak_list)) return(list())

    dt <- data.table::rbindlist(peak_list)
    dt <- dt[, .(intensity = max(intensity)), by = .(scan_idx, mz_bin)]

    scan_order <- match(dt$scan_idx, win_idx)
    dt[, scan_order := scan_order]

    n_scans <- length(win_idx)
    mz_bins <- sort(unique(dt$mz_bin))

    eic_list <- vector("list", length(mz_bins))
    names(eic_list) <- as.character(mz_bins)
    for (j in seq_along(mz_bins)) {
      mb <- mz_bins[j]
      sub <- dt[mz_bin == mb]
      eic <- numeric(n_scans)
      eic[sub$scan_order] <- sub$intensity
      eic_list[[j]] <- eic
    }

    apex_rows <- vector("list", length(mz_bins))
    ari <- 0L
    for (j in seq_along(mz_bins)) {
      eic <- eic_list[[j]]
      apexes <- .nts_find_eic_apexes(eic, min_snr)
      if (!length(apexes)) next
      for (apex in apexes) {
        ari <- ari + 1L
        noise_sd <- .nts_eic_noise_sd(eic)
        if (noise_sd < 1e-6) noise_sd <- 1e-6
        apex_rows[[ari]] <- data.table::data.table(
          mz_bin      = mz_bins[j],
          apex_order  = apex,
          apex_rt     = scan_rt[win_idx[apex]],
          intensity   = eic[apex],
          snr         = eic[apex] / noise_sd
        )
      }
    }
    if (!ari) return(list())

    apex_dt <- data.table::rbindlist(apex_rows[1L:ari])
    data.table::setorder(apex_dt, apex_rt)

    apex_dt[, cluster := .nts_cluster_apexes(apex_rt, apex_tol)]

    clusters <- unique(apex_dt$cluster)
    comps <- vector("list", length(clusters))
    ci <- 0L
    for (cl in clusters) {
      crows <- apex_dt[cluster == cl]
      if (nrow(crows) < min_ions) next

      best <- which.max(crows$snr)
      apex_order <- crows$apex_order[best]
      apex_rt_val <- crows$apex_rt[best]

      spec_mz <- numeric(nrow(crows))
      spec_int <- numeric(nrow(crows))
      for (k in seq_len(nrow(crows))) {
        mb <- crows$mz_bin[k]
        eic <- eic_list[[as.character(mb)]]
        spec_mz[k] <- mb
        spec_int[k] <- eic[apex_order]
      }

      keep <- spec_int > 0
      if (sum(keep) < min_ions) next
      spec_mz <- spec_mz[keep]
      spec_int <- spec_int[keep]

      spec_int <- spec_int / max(spec_int) * 100

      ci <- ci + 1L
      comps[[ci]] <- list(
        spectrum = data.frame(mz = spec_mz, intensity = spec_int, stringsAsFactors = FALSE),
        apex_rt  = apex_rt_val,
        n_ions   = length(spec_mz)
      )
    }

    comps[seq_len(ci)]
  }

  .nts_assert_parallel_ntsminer()
  results <- BiocParallel::bplapply(seq_len(n), .deconv_one, BPPARAM = BPPARAM)
  names(results) <- ft_names

  n_components <- sum(vapply(results, length, integer(1L)))
  n_deconv <- sum(vapply(results, function(x) length(x) > 0L, logical(1L)))
  message("Deconvolution: ", n_components, " components from ", n_deconv, " features")

  project@gc_deconvolved <- results
  validObject(project)
  project
}


# Internal helpers for Biller-Biemann deconvolution ----------------------------

#' Estimate noise standard deviation from an EIC vector
#'
#' Uses the lower portion of non-zero intensity values (below median + 2*MAD)
#' to estimate background noise SD.
#'
#' @param eic Numeric vector of EIC intensities.
#' @return A non-negative noise SD estimate.
#' @noRd
.nts_eic_noise_sd <- function(eic) {
  low <- eic[eic > 0]
  if (length(low) < 5) return(1)
  threshold <- stats::median(low) + 2 * stats::mad(low)
  if (!is.finite(threshold)) return(stats::sd(low))
  noise_vals <- low[low <= threshold]
  if (length(noise_vals) < 3) return(stats::sd(low))
  stats::sd(noise_vals)
}


#' Find local maxima in an EIC vector
#'
#' @param eic Numeric vector of EIC intensities.
#' @param min_snr Minimum signal-to-noise ratio for a valid apex.
#' @return Integer vector of 1-based indices of valid apexes.
#' @noRd
.nts_find_eic_apexes <- function(eic, min_snr = 3) {
  eic <- pmax(eic, 0)
  n <- length(eic)
  if (n < 3) return(integer(0L))

  d <- diff(eic)
  pos <- which(d[-length(d)] > 0 & d[-1L] < 0) + 1L

  if (!length(pos)) return(integer(0L))

  noise_sd <- .nts_eic_noise_sd(eic)
  if (noise_sd < 1e-6) noise_sd <- 1e-6

  snr <- eic[pos] / noise_sd
  pos <- pos[snr >= min_snr]

  if (length(pos) > 1L) {
    gaps <- diff(pos)
    keep <- rep(TRUE, length(pos))
    for (k in seq_len(length(pos) - 1L)) {
      if (gaps[k] <= 2L) {
        if (eic[pos[k]] >= eic[pos[k + 1L]]) {
          keep[k + 1L] <- FALSE
        } else {
          keep[k] <- FALSE
        }
      }
    }
    pos <- pos[keep]
  }

  pos
}


#' Greedy RT clustering of EIC apexes
#'
#' @param apex_rt Numeric vector of apex retention times, sorted ascending.
#' @param apex_tol Maximum RT difference (seconds) for grouping into one cluster.
#' @return Integer vector of cluster IDs (same length as apex_rt).
#' @noRd
.nts_cluster_apexes <- function(apex_rt, apex_tol) {
  n <- length(apex_rt)
  if (n == 0L) return(integer(0L))
  if (n == 1L) return(rep(1L, 1L))

  cluster <- integer(n)
  cluster[1L] <- 1L
  cid <- 1L
  for (i in seq(2L, n)) {
    if (apex_rt[i] - apex_rt[i - 1L] > apex_tol) {
      cid <- cid + 1L
    }
    cluster[i] <- cid
  }
  cluster
}


#' Cluster m/z channels by EIC shape similarity (not just RT proximity)
#'
#' In GC-EI, all fragments from the same compound share an identical
#' elution profile. Co-eluting compounds have fragments with different
#' EIC shapes even when RT difference is < 1s. This function uses
#' greedy seed-based Pearson correlation clustering to group ions
#' by elution profile similarity.
#'
#' Algorithm:
#' 1. Pick the m/z bin with highest total intensity as seed
#' 2. Compute Pearson r between seed EIC and all remaining EICs
#' 3. Group all bins with r >= min_cor into a cluster
#' 4. Remove clustered bins, repeat with next strongest remaining
#'
#' @param eic_matrix Numeric matrix, rows = scans, cols = m/z bins.
#' @param min_cor Numeric, minimum Pearson r for same-compound grouping.
#' @param min_ions Integer, minimum ions per compound.
#' @param intensity Optional numeric vector of total intensities per bin.
#' @return List of integer vectors, each = column indices of one compound.
#' @noRd
.nts_eic_shape_cluster <- function(eic_matrix, min_cor = 0.85,
                                   min_ions = 5L, intensity = NULL) {
  n_bins <- ncol(eic_matrix)
  if (n_bins < min_ions) return(list())

  if (is.null(intensity)) {
    intensity <- colSums(eic_matrix)
  }

  remaining <- seq_len(n_bins)
  clusters  <- list()
  ci <- 0L

  while (length(remaining) >= min_ions) {
    idx <- which.max(intensity[remaining])
    seed_col <- remaining[idx]
    seed_eic <- eic_matrix[, seed_col]

    cors <- apply(eic_matrix[, remaining, drop = FALSE], 2,
                  function(x) {
                    r <- suppressWarnings(stats::cor(seed_eic, x))
                    if (is.na(r)) 0 else r
                  })

    cluster_cols <- remaining[cors >= min_cor]
    if (length(cluster_cols) < min_ions) break

    ci <- ci + 1L
    clusters[[ci]] <- cluster_cols
    remaining <- setdiff(remaining, cluster_cols)
  }

  clusters
}


#' Compute FWHM (full width at half maximum) of an EIC peak
#'
#' Searches left and right from the apex until the signal falls
#' below half of the apex intensity.
#'
#' @param eic Numeric vector, the extracted ion chromatogram.
#' @param apex_idx Integer, index of the peak apex in the EIC.
#' @return Named numeric: c(fwhm = scans, left = idx, right = idx).
#' @noRd
.nts_peak_fwhm <- function(eic, apex_idx) {
  n <- length(eic)
  if (apex_idx < 1L || apex_idx > n) return(c(fwhm = NA_real_, left = NA_real_, right = NA_real_))

  apex_val <- eic[apex_idx]
  if (apex_val <= 0) return(c(fwhm = NA_real_, left = NA_real_, right = NA_real_))
  half_max <- apex_val * 0.5

  left <- apex_idx
  while (left > 1L && eic[left] > half_max) left <- left - 1L
  right <- apex_idx
  while (right < n && eic[right] > half_max) right <- right + 1L

  fwhm <- max(1L, right - left)
  c(fwhm = as.numeric(fwhm), left = as.numeric(left), right = as.numeric(right))
}


#' Compute peak symmetry (ratio of left to right half-width)
#'
#' Perfect symmetry = 1.0. Values < 1 indicate fronting, > 1 tailing.
#' Returns symmetry in [0, 1] where 1 = perfect and 0 = completely asymmetric.
#'
#' @param eic Numeric vector, the extracted ion chromatogram.
#' @param apex_idx Integer, index of the peak apex.
#' @return Numeric in [0, 1].
#' @noRd
.nts_peak_symmetry <- function(eic, apex_idx) {
  fwhm_info <- .nts_peak_fwhm(eic, apex_idx)
  if (any(is.na(fwhm_info))) return(0)

  left_hw  <- apex_idx - fwhm_info["left"]
  right_hw <- fwhm_info["right"] - apex_idx
  if (left_hw <= 0 || right_hw <= 0) return(0)

  min(left_hw, right_hw) / max(left_hw, right_hw)
}


# erah-style full-chromatogram deconvolution -----------------------------------

#' Detect compounds in GC-EI data by TIC-guided multi-component deconvolution
#'
#' Uses TIC peak detection to identify chromatographic regions of interest,
#' then applies Biller-Biemann deconvolution within each region to separate
#' co-eluting compounds into pure component spectra. A per-m/z-bin pre-filter
#' (max/median intensity ratio) excludes noise-only channels for efficiency.
#' This is the primary feature detection step for GC-EI data, replacing
#' CentWave-based peak picking.
#'
#' @param mzml_files Character vector of mzML file paths.
#' @param min_snr Numeric, minimum signal-to-noise ratio for an EIC apex
#'   to be considered valid. Agilent default 2. Default 2.
#' @param apex_tol Numeric, RT tolerance in seconds for merging adjacent
#'   EIC apexes within the same m/z channel. No longer used for cross-m/z
#'   clustering (which now uses EIC shape correlation). Default 1.0.
#' @param min_ions Integer, minimum number of m/z channels per compound.
#'   Agilent default 3. Default 3.
#' @param avoid_mz Numeric vector of m/z values to exclude (column bleed,
#'   e.g. 73, 74, 75, 147, 148, 149 for derivatized samples, or
#'   c(207, 281) for siloxanes). Default c(73, 74, 75, 147, 148, 149).
#' @param mz_bin_width Numeric, m/z bin width in Da for spectrum binning.
#'   Default 0.1 (unit-mass). Use 0.01 for QTOF high-resolution data.
#' @param shape_min_cor Numeric, minimum Pearson correlation for grouping
#'   m/z channels into the same compound by EIC shape similarity.
#'   Default 0.85. Lower to 0.7 for looser grouping, raise to 0.95 for
#'   stricter separation.
#' @param peak_ratio_min Numeric, minimum max/median intensity ratio for
#'   an m/z channel to pass the pre-filter (chromatographic peak shape).
#'   Default 2. Lower values allow more ions into deconvolution.
#' @param peak_sym_min Numeric in `[0,1]`, minimum symmetry of the base peak
#'   EIC. Compounds whose base peak symmetry < this are rejected.
#'   Default 0.6.
#' @param delta_apex_max Numeric, maximum RT spread (seconds) of individual
#'   ion apexes within a compound cluster. Larger spreads suggest mixed
#'   compounds. Default 2.0.
#' @param delta_fwhm_max Numeric, maximum ratio of widest to median FWHM
#'   within a cluster. Ions with FWHM > median * this are excluded.
#'   Default 5.0.
#' @param min_peak_width Numeric, minimum chromatographic peak width in
#'   seconds for TIC region detection. Default 1.
#' @param baseline_correct Logical, apply running-median baseline
#'   correction. Default TRUE.
#' @param smooth_window Integer, moving-average smoothing half-window
#'   size in scans. 0 disables smoothing. Default 3.
#' @param baseline_window Integer, running-median baseline correction
#'   window width in scans. Passed to internal \code{.nts_baseline_correct}.
#'   Default 31.
#' @param BPPARAM A BiocParallelParam object. Default SerialParam.
#'
#' @return A data.table with columns:
#'   \item{compound_id}{Unique compound identifier (e.g. "C001")}
#'   \item{filename}{Source mzML filename}
#'   \item{rt}{Retention time in seconds (component apex)}
#'   \item{basepeak_mz}{m/z of the most intense ion}
#'   \item{basepeak_int}{Intensity of the base peak}
#'   \item{spectrum}{List-column: data.frame(mz, intensity), the pure
#'     EI spectrum normalized to 0-100}
#'   \item{n_ions}{Number of m/z channels in the deconvolved spectrum}
#'   \item{tic_apex}{TIC intensity at apex (proxy for abundance)}
#'
#' @section erah integration:
#' For improved deconvolution accuracy, consider using [run_erah_detection]
#' which wraps erah's PCA-OSD algorithm. In benchmarks on GC-QTOF pesticide
#' standards, erah achieved 52% identification rate (cos>=0.7) vs 13% for
#' this self-implemented pipeline.
#'
#' @references
#' Biller, J. E. & Biemann, K. (1974). Analytical Letters, 7(7), 515-528.
#' Domingo-Almenara, X. et al. (2016). Analytical Chemistry, 88(19), 9821-9829.
#'
#' @importFrom mzR openMSfile header peaks close
#' @importFrom data.table data.table rbindlist setorder setnames
#' @export
detect_gc_compounds <- function(mzml_files,
                                min_snr           = 2,
                                apex_tol          = 1.0,
                                min_ions          = 3L,
                                avoid_mz          = c(73, 74, 75, 147, 148, 149),
                                mz_bin_width      = 0.1,
                                shape_min_cor     = 0.85,
                                peak_ratio_min    = 2,
                                peak_sym_min      = 0.6,
                                delta_apex_max    = 5.0,
                                delta_fwhm_max    = 2.0,
                                min_peak_width    = 1,
                                baseline_correct  = TRUE,
                                smooth_window     = 3L,
                                baseline_window   = 31L,
                                BPPARAM           = BiocParallel::SerialParam()) {
  if (!requireNamespace("mzR", quietly = TRUE)) stop("mzR is required", call. = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)

  mzml_files <- as.character(mzml_files)
  mzml_files <- mzml_files[file.exists(mzml_files)]
  if (!length(mzml_files)) stop("No valid mzML files provided", call. = FALSE)

  min_ions   <- max(1L, as.integer(min_ions)[1L])
  avoid_mz   <- round(as.numeric(avoid_mz), 1L)
  avoid_mz   <- avoid_mz[is.finite(avoid_mz)]
  baseline_window <- max(5L, as.integer(baseline_window)[1L])
  if (baseline_window %% 2L == 0L) baseline_window <- baseline_window + 1L
  mz_digits  <- max(0L, as.integer(round(-log10(mz_bin_width[1L]))))
  shape_min_cor <- max(0.5, min(0.99, as.numeric(shape_min_cor[1L])))
  peak_ratio_min <- max(1, as.numeric(peak_ratio_min[1L]))
  peak_sym_min   <- max(0, min(1, as.numeric(peak_sym_min[1L])))
  delta_apex_max <- max(0.05, as.numeric(delta_apex_max[1L]))
  delta_fwhm_max <- max(1.1, as.numeric(delta_fwhm_max[1L]))

  .process_one_file <- function(f) {
    ms <- mzR::openMSfile(f)
    on.exit(try(mzR::close(ms), silent = TRUE), add = TRUE)

    hdr <- mzR::header(ms)
    if (!is.data.frame(hdr) || !nrow(hdr)) {
      return(data.table::data.table(
        compound_id = character(0L), filename = character(0L),
        rt = numeric(0L), basepeak_mz = numeric(0L),
        basepeak_int = numeric(0L), spectrum = list(),
        n_ions = integer(0L), tic_apex = numeric(0L)
      ))
    }

    scan_rt   <- as.numeric(hdr$retentionTime)
    scan_tic  <- as.numeric(hdr$totIonCurrent)
    n_scans   <- length(scan_rt)

    # --- Step 1: Load all spectra into sparse data.table ---
    peak_list <- vector("list", n_scans)
    for (si in seq_len(n_scans)) {
      sp <- try(mzR::peaks(ms, si), silent = TRUE)
      if (inherits(sp, "try-error") || !is.matrix(sp) || !nrow(sp)) next
      sp_mz <- as.numeric(sp[, "mz"])
      sp_int <- as.numeric(sp[, "intensity"])
      keep <- sp_mz >= 35 & sp_mz <= 800 & sp_int > 0
      if (!any(keep)) next
      sp_mz <- round(sp_mz[keep], mz_digits)
      sp_int <- sp_int[keep]
      if (length(avoid_mz)) {
        ok <- !sp_mz %in% avoid_mz
        sp_mz <- sp_mz[ok]
        sp_int <- sp_int[ok]
      }
      if (!length(sp_mz)) next
      peak_list[[si]] <- data.table::data.table(
        scan_order = si,
        mz_bin     = sp_mz,
        intensity  = sp_int
      )
    }

    dt <- data.table::rbindlist(peak_list)
    rm(peak_list)
    if (!nrow(dt)) return(data.table::data.table(
      compound_id = character(0L), filename = character(0L),
      rt = numeric(0L), basepeak_mz = numeric(0L),
      basepeak_int = numeric(0L), spectrum = list(),
      n_ions = integer(0L), tic_apex = numeric(0L)
    ))

    # Collapse duplicate mz_bin per scan
    dt <- dt[, .(intensity = max(intensity)), by = .(scan_order, mz_bin)]

    # --- Step 2: TIC-based region detection (coarse pre-filter) ---
    # Use TIC peaks to divide the chromatogram into regions of interest.
    # Unlike the old "TIC first" approach, TIC is a pre-filter, not a gate:
    # gaps between TIC peaks are also scanned for low-abundance compounds.
    tic_full <- scan_tic
    tic_full[!is.finite(tic_full)] <- 0

    tic_smoothed <- .nts_smooth_eic(tic_full, 2L)
    tic_apexes <- .nts_find_eic_apexes(tic_smoothed, min_snr = 3)

    # Build TIC peak regions (FWHM boundaries)
    tic_regions <- list()
    if (length(tic_apexes)) {
      tic_regions <- lapply(tic_apexes, function(apex) {
        peak_height <- tic_smoothed[apex]
        half_h <- peak_height * 0.5
        left <- apex
        while (left > 1L && tic_smoothed[left] > half_h) left <- left - 1L
        right <- apex
        while (right < n_scans && tic_smoothed[right] > half_h) right <- right + 1L
        width <- right - left + 1L
        if (width < 6L) {
          extra <- 6L - width
          left  <- max(1L, left - extra %/% 2L)
          right <- min(n_scans, right + extra %/% 2L)
        }
        max_width <- 85L
        if (right - left > max_width) {
          mid <- (left + right) %/% 2L
          left  <- max(1L, mid - max_width %/% 2L)
          right <- min(n_scans, mid + max_width %/% 2L)
        }
        c(left, right)
      })

      # Merge overlapping TIC regions
      tic_regions <- tic_regions[order(vapply(tic_regions, `[`, numeric(1L), 1L))]
      if (length(tic_regions) > 1L) {
        dedup <- list(tic_regions[[1L]])
        for (k in seq(2L, length(tic_regions))) {
          prev <- dedup[[length(dedup)]]
          curr <- tic_regions[[k]]
          if (curr[1L] <= prev[2L]) {
            dedup[[length(dedup)]][2L] <- max(prev[2L], curr[2L])
          } else {
            dedup[[length(dedup) + 1L]] <- curr
          }
        }
        tic_regions <- dedup
      }
    }

    # --- Step 3: Per-TIC-peak multi-component deconvolution ---
    # Within each TIC peak region, find EIC apexes per m/z channel and
    # cluster by RT -> one cluster = one compound (Biller-Biemann).
    # TIC is used as a pre-filter to divide the chromatogram into regions
    # with meaningful signal; each region is deconvolved independently.
    all_compounds <- vector("list", length(tic_regions) * 4L)
    comp_count <- 0L

    for (region in tic_regions) {
      r_start <- region[1L]
      r_end   <- region[2L]
      if (r_end - r_start < 4L) next

      dt_sub <- dt[scan_order >= r_start & scan_order <= r_end]
      if (!nrow(dt_sub)) next

      # Pre-filter m/z bins: only keep those with chromatographic peak shape
      # (max intensity at least 5x the median, indicating a real peak vs noise)
      mz_stats <- dt_sub[, .(
        mx = max(intensity),
        md = stats::median(intensity)
      ), by = mz_bin]
      mz_stats[, peak_ratio := mx / pmax(md, 0.01)]
      mz_list <- mz_stats[peak_ratio >= peak_ratio_min, mz_bin]
      if (!length(mz_list)) next
      mz_list <- sort(mz_list)
      n_scans_sub <- r_end - r_start + 1L

      apex_rows <- vector("list", length(mz_list))
      eic_cache <- list()
      ari <- 0L

      # Build EIC for all pre-filtered m/z bins (always cache for shape clustering)
      # and detect apexes per bin for spectrum building.
      for (j in seq_along(mz_list)) {
        mb <- mz_list[j]
        mz_bin_key <- as.character(mb)
        sub <- dt_sub[mz_bin == mb]
        eic <- numeric(n_scans_sub)
        eic[sub$scan_order - r_start + 1L] <- sub$intensity

        if (baseline_correct) eic <- .nts_baseline_correct(eic, baseline_window)

        # Always cache EIC (for shape correlation)
        eic_cache[[mz_bin_key]] <- eic

        # Detect apexes
        apexes <- .nts_find_eic_apexes(eic, min_snr)
        if (!length(apexes)) next

        noise_sd <- .nts_eic_noise_sd(eic)
        if (noise_sd < 1e-6) noise_sd <- 1e-6

        for (apex in apexes) {
          ari <- ari + 1L
          apex_scan <- r_start + apex - 1L
          apex_rows[[ari]] <- data.table::data.table(
            mz_bin      = mb,
            apex_scan   = apex_scan,
            apex_order  = apex,
            apex_rt     = scan_rt[apex_scan],
            intensity   = eic[apex],
            snr         = eic[apex] / noise_sd
          )
        }
      }

      if (length(eic_cache) < min_ions) next

      # Build EIC matrix: rows = scans, cols = cached m/z bins
      cached_names <- names(eic_cache)
      n_cached <- length(cached_names)
      eic_mat <- matrix(0, nrow = n_scans_sub, ncol = n_cached)
      colnames(eic_mat) <- cached_names
      intensities <- numeric(n_cached)
      for (j in seq_len(n_cached)) {
        eic <- eic_cache[[cached_names[j]]]
        eic_mat[, j] <- eic
        intensities[j] <- sum(eic)
      }

      # --- EIC shape clustering (replaces RT-only clustering) ---
      shape_clusters <- .nts_eic_shape_cluster(
        eic_mat, min_cor = shape_min_cor, min_ions = min_ions,
        intensity = intensities
      )

      if (!length(shape_clusters)) next

      # Build apex lookup: mz_bin -> apex_row index (fast lookup)
      apex_dt <- if (ari > 0L) data.table::rbindlist(apex_rows[1L:ari]) else NULL

      for (sc in shape_clusters) {
        if (length(sc) < min_ions) next

        cl_names <- cached_names[sc]
        cl_mz <- as.numeric(cl_names)
        n_cl <- length(cl_names)

        cl_apex <- if (!is.null(apex_dt)) apex_dt[mz_bin %in% cl_mz] else data.table::data.table()
        if (nrow(cl_apex) == 0L) next

        # --- Post-cluster quality validation ---
        # (A) Compute FWHM, symmetry, and apex RT for every ion in the cluster
        cl_fwhm    <- numeric(n_cl)
        cl_sym     <- numeric(n_cl)
        cl_apex_rt <- numeric(n_cl)
        cl_intensity <- numeric(n_cl)
        cl_apex_order <- integer(n_cl)

        for (k in seq_len(n_cl)) {
          eic <- eic_cache[[cl_names[k]]]
          this_apex <- cl_apex[mz_bin == cl_mz[k]]
          if (nrow(this_apex) == 0L) {
            cl_fwhm[k] <- NA_real_
            cl_sym[k]  <- 0
            cl_apex_rt[k] <- NA_real_
            cl_intensity[k] <- 0
            next
          }
          ao <- this_apex$apex_order[1L]
          cl_apex_order[k] <- ao
          fwhm_info <- .nts_peak_fwhm(eic, ao)
          cl_fwhm[k]    <- fwhm_info["fwhm"]
          cl_sym[k]     <- .nts_peak_symmetry(eic, ao)
          cl_apex_rt[k] <- this_apex$apex_rt[1L]
          cl_intensity[k] <- this_apex$intensity[1L]
        }

        # (B) Base peak symmetry check
        bp_local <- which.max(cl_intensity)
        if (cl_sym[bp_local] < peak_sym_min) next

        # (C) Cluster-wide apex RT alignment
        valid_apex <- is.finite(cl_apex_rt)
        if (sum(valid_apex) < min_ions) next
        apex_range <- diff(range(cl_apex_rt[valid_apex]))
        if (apex_range > delta_apex_max) next

        # (D) FWHM consistency filter
        fwhm_pos <- cl_fwhm > 0 & is.finite(cl_fwhm)
        med_fwhm <- if (sum(fwhm_pos) > 0L) stats::median(cl_fwhm[fwhm_pos]) else NA_real_
        if (is.na(med_fwhm) || med_fwhm <= 0) next
        fwhm_ok <- fwhm_pos & cl_fwhm <= med_fwhm * delta_fwhm_max

        # Keep only validated ions
        keep_ions <- valid_apex & fwhm_ok & cl_intensity > 0
        if (sum(keep_ions) < min_ions) next

        cl_names_clean <- cl_names[keep_ions]
        cl_mz_clean    <- cl_mz[keep_ions]
        n_cl_clean <- length(cl_names_clean)

        # (E) Reference scan: base peak apex
        ref_scan_local <- cl_apex_order[bp_local]
        apex_rt_val    <- cl_apex_rt[bp_local]

        # --- Biller-Biemann spectrum extraction ---
        spec_mz <- numeric(n_cl_clean)
        spec_int <- numeric(n_cl_clean)
        for (k in seq_len(n_cl_clean)) {
          eic <- eic_cache[[cl_names_clean[k]]]
          spec_mz[k] <- cl_mz_clean[k]
          this_apex <- cl_apex[mz_bin == cl_mz_clean[k]]
          if (nrow(this_apex) > 0L) {
            spec_int[k] <- this_apex$intensity[1L]
          } else {
            spec_int[k] <- if (ref_scan_local <= length(eic)) eic[ref_scan_local] else 0
          }
        }

        # 1% base peak filter + normalize 0-100
        keep <- spec_int > 0
        if (sum(keep) < min_ions) next
        spec_mz <- spec_mz[keep]
        spec_int <- spec_int[keep]

        bp <- max(spec_int)
        keep <- spec_int >= bp * 0.01
        if (sum(keep) < min_ions) next
        spec_mz <- spec_mz[keep]
        spec_int <- spec_int[keep]

        spec_int <- spec_int / bp * 100

        bp_idx <- which.max(spec_int)
        apex_scan_ref <- r_start + ref_scan_local - 1L
        comp_count <- comp_count + 1L
        all_compounds[[comp_count]] <- list(
          spectrum     = data.frame(
            mz = spec_mz, intensity = spec_int, stringsAsFactors = FALSE
          ),
          rt           = apex_rt_val,
          basepeak_mz  = spec_mz[bp_idx],
          basepeak_int = spec_int[bp_idx] * max(spec_int) / 100 *
            if (apex_scan_ref <= length(scan_tic)) scan_tic[apex_scan_ref] / 100 else 1,
          n_ions       = length(spec_mz),
          tic_apex     = if (apex_scan_ref <= length(scan_tic)) scan_tic[apex_scan_ref] else NA_real_
        )
      }
    }

    if (!comp_count) return(data.table::data.table(
      compound_id = character(0L), filename = character(0L),
      rt = numeric(0L), basepeak_mz = numeric(0L),
      basepeak_int = numeric(0L), spectrum = list(),
      n_ions = integer(0L), tic_apex = numeric(0L)
    ))

    # Deduplicate: remove near-duplicates (same RT + same base peak m/z)
    comp_dt <- data.table::rbindlist(lapply(all_compounds[1L:comp_count], function(x) {
      data.table::data.table(
        rt = x$rt, basepeak_mz = x$basepeak_mz,
        basepeak_int = x$basepeak_int,
        spectrum = list(x$spectrum),
        n_ions = x$n_ions, tic_apex = x$tic_apex
      )
    }))

    data.table::setorder(comp_dt, -basepeak_int)
    comp_dt[, rt_bin := round(rt, 0L)]
    comp_dt[, bp_bin := round(basepeak_mz, 0L)]
    comp_dt <- comp_dt[!duplicated(comp_dt[, .(rt_bin, bp_bin)])]
    comp_dt[, `:=`(rt_bin = NULL, bp_bin = NULL)]
    data.table::setorder(comp_dt, rt)

    comp_dt[, compound_id := sprintf("C%04d", seq_len(.N))]
    comp_dt[, filename := basename(f)]

    comp_dt[, .(compound_id, filename, rt, basepeak_mz, basepeak_int,
                spectrum, n_ions, tic_apex)]
  }

  .nts_assert_parallel_ntsminer()
  results <- BiocParallel::bplapply(mzml_files, .process_one_file, BPPARAM = BPPARAM)
  out <- data.table::rbindlist(results)
  ## compound_id is assigned per-file inside .process_one_file, so it restarts at
  ## C0001 for every file and collides across files. Re-index the combined table
  ## so ids are globally unique; otherwise name-keyed lookups downstream
  ## (gc_spectra list names, and the feature_table intensity fill in
  ## run_gc_detection via match(compound_id, name)) silently hit the wrong
  ## file's rows. Single-file input is unchanged (ids 1..N either way).
  out[, compound_id := sprintf("C%04d", seq_len(.N))]
  data.table::setorder(out, filename, rt)
  out
}


# Internal helpers for erah-style pipeline -------------------------------------

#' Running-median baseline correction for an EIC vector
#' @param eic Numeric vector of intensities.
#' @param window Integer, running median window width in scans. Default 31.
#' @return Baseline-corrected (and zero-clamped) EIC.
#' @noRd
.nts_baseline_correct <- function(eic, window = 31L) {
  n <- length(eic)
  if (n < window) return(eic)
  baseline <- stats::runmed(eic, k = min(window, n %/% 2L * 2L + 1L))
  corrected <- eic - baseline
  pmax(corrected, 0)
}


#' Simple moving-average smoothing for an EIC vector
#' @param eic Numeric vector of intensities.
#' @param half_window Integer, number of points each side. Default 3.
#' @return Smoothed EIC.
#' @noRd
.nts_smooth_eic <- function(eic, half_window = 3L) {
  hw <- as.integer(half_window[1L])
  if (hw < 1L) return(eic)
  n <- length(eic)
  out <- numeric(n)
  for (i in seq_len(n)) {
    lo <- max(1L, i - hw)
    hi <- min(n, i + hw)
    out[i] <- mean(eic[lo:hi])
  }
  out
}


#' Safe quantile (returns fallback on short input)
#' @param x Numeric vector.
#' @param prob Quantile probability.
#' @return Quantile value or 0.
#' @noRd
.nts_try_quantile <- function(x, prob = 0.25) {
  if (length(x) < 3L) return(0)
  q <- tryCatch(stats::quantile(x, prob, na.rm = TRUE), error = function(e) 0)
  if (!is.finite(q)) 0 else q
}


# NTSProject integration -------------------------------------------------------

#' Run GC compound detection on project mzML files
#'
#' Wrapper around [detect_gc_compounds] that operates on an NTSProject,
#' filling the \code{gc_spectra} slot with deconvolved compound spectra and
#' the \code{feature_table} slot with compound-level features.
#'
#' @param project An [NTSProject] with \code{raw_data_path} pointing to mzML files.
#' @param ... Arguments passed to [detect_gc_compounds].
#'
#' @return The project with \code{gc_spectra}, \code{feature_table},
#'   and \code{ms1_annotations}
#'   populated with GC compound data.
#'
#' @export
run_gc_detection <- function(project, ...) {
  .nts_check(project, "project", .arg = "project")

  mzml_files <- list.files(
    project@raw_data_path,
    pattern = "\\.mzML$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(mzml_files)) stop("No mzML files found in raw_data_path", call. = FALSE)

  compounds <- detect_gc_compounds(mzml_files, ...)

  if (!nrow(compounds)) {
    warning("No compounds detected in ", length(mzml_files), " mzML file(s)", call. = FALSE)
    return(project)
  }

  # Fill @gc_spectra
  spec_list <- lapply(seq_len(nrow(compounds)), function(i) compounds$spectrum[[i]])
  names(spec_list) <- compounds$compound_id
  project@gc_spectra <- spec_list

  # Build feature_table (compound-level)
  ft <- data.table::data.table(
    name = compounds$compound_id,
    mz   = compounds$basepeak_mz,
    rt   = compounds$rt
  )
  for (fnm in unique(compounds$filename)) {
    col_name <- tools::file_path_sans_ext(basename(fnm))
    ft_sub <- compounds[filename == fnm]
    ft[[col_name]] <- NA_real_
    ft[match(ft_sub$compound_id, ft$name), (col_name) := ft_sub$basepeak_int]
  }
  project@feature_table <- .nts_round_feature_intensities(as.data.frame(ft, stringsAsFactors = FALSE))

  # Build ms1_annotations (minimal columns)
  ms1 <- data.frame(
    name = compounds$compound_id,
    mz = compounds$basepeak_mz,
    rt = compounds$rt,
    stringsAsFactors = FALSE
  )
  project@ms1_annotations <- .nts_ms1_partition(.nts_ms1_sanitize_and_order(ms1))

  message(sprintf(
    "GC detection: %d compounds from %d file(s)",
    nrow(compounds), length(mzml_files)
  ))
  validObject(project)
  project
}


# ==== erah-based GC-EI deconvolution ====

#' Run GC-EI compound detection using erah PCA-OSD deconvolution + alignment
#'
#' Wraps erah's \code{deconvolveComp()} and \code{alignComp()} to perform
#' spectral deconvolution and cross-sample alignment on GC-QTOF mzML files,
#' then converts the output to NTSminer's standard data structures.
#'
#' After alignment, each compound in \code{@feature_table} has a \code{found_in}
#' column (integer count of samples where detected) and per-sample columns
#' containing chromatographic peak areas (from erah's \code{Area} field).
#' The \code{@gc_spectra} slot contains consensus spectra (averaged across
#' samples) for each aligned compound.
#'
#' @param project An NTSProject with valid \code{@raw_data_path}.
#' @param min.peak.width Numeric, minimum peak width in seconds (default 1.0).
#' @param min.peak.height Numeric, minimum peak height (default 1000).
#' @param noise.threshold Numeric, noise threshold (default 500).
#' @param avoid.processing.mz Numeric vector, m/z to exclude (default column bleed).
#' @param analysis.time Numeric vector of length 2, RT window in minutes (default c(5,40)).
#' @param align.min_cor Numeric, minimum spectral correlation for
#'   cross-sample alignment (0-1). Default 0.70.
#' @param align.max_rt_diff Numeric, maximum RT difference in seconds
#'   for cross-sample alignment. Default 0.50.
#' @param align.mz_range Integer vector, m/z range for consensus
#'   spectrum binning. Default 40:600.
#' @param cache_dir Character, optional directory for erah caching.
#'   Default NULL.
#' @param compression.coef Numeric, erah compression coefficient. Default 1.
#'
#' @return The NTSProject with \code{@gc_spectra} (consensus spectra),
#'   \code{@feature_table} (aligned features with \code{found_in} and
#'   per-sample peak areas), and \code{@ms1_annotations} populated.
#'
#' @seealso \code{\link[erah]{deconvolveComp}}, \code{\link[erah]{alignComp}}
#' @export
run_erah_detection <- function(project,
                                min.peak.width      = 1.0,
                                min.peak.height     = 1000,
                                noise.threshold     = 500,
                                avoid.processing.mz = c(35:69, 73:75, 147:149),
                                analysis.time       = c(5, 40),
                                align.min_cor       = 0.70,
                                align.max_rt_diff   = 0.50,
                                align.mz_range      = 40:600,
                                cache_dir           = NULL,
                                compression.coef    = 1) {

  .nts_check(project, "project", .arg = "project")

  if (!requireNamespace("erah", quietly = TRUE))
    stop("erah is required. Install with: remotes::install_github('xdomingoal/erah-devel')", call. = FALSE)

  mzml_files <- list.files(
    project@raw_data_path,
    pattern = "\\.mzML$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (!length(mzml_files)) stop("No mzML files found in raw_data_path", call. = FALSE)

  # ---- Create erah experiment ----
  tmp_dir <- file.path(tempdir(), "erah_ntsminer")
  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  inst_csv <- file.path(tmp_dir, "instrumental.csv")
  pheno_csv <- file.path(tmp_dir, "phenotype.csv")

  inst_df <- data.frame(
    sampleID  = tools::file_path_sans_ext(basename(mzml_files)),
    filename  = normalizePath(mzml_files, winslash = "/"),
    date      = rep("2024-01-01", length(mzml_files)),
    time      = rep("12:00:00", length(mzml_files)),
    class     = rep("sample", length(mzml_files)),
    stringsAsFactors = FALSE
  )
  write.table(inst_df, inst_csv, sep = ";", row.names = FALSE, quote = FALSE)

  pheno_df <- data.frame(
    sampleID = inst_df$sampleID,
    class    = "sample",
    stringsAsFactors = FALSE
  )
  write.table(pheno_df, pheno_csv, sep = ";", row.names = FALSE, quote = FALSE)

  ex <- erah::newExp(
    instrumental = inst_csv,
    phenotype    = pheno_csv,
    info         = c(project = "NTSminer")
  )

  dec_par <- erah::setDecPar(
    min.peak.width       = min.peak.width,
    min.peak.height      = min.peak.height,
    noise.threshold      = noise.threshold,
    avoid.processing.mz  = avoid.processing.mz,
    compression.coef     = compression.coef,
    analysis.time        = analysis.time
  )

  # ---- Run erah deconvolution ----
  message("Running erah PCA-OSD deconvolution on ", length(mzml_files), " file(s)...")
  t0 <- Sys.time()
  ex <- erah::deconvolveComp(ex, decParameters = dec_par)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("erah deconvolution: %.0f sec", elapsed))

  # ---- Run erah cross-sample alignment ----
  message("Running erah cross-sample alignment...")
  t0_align <- Sys.time()
  ex <- erah::alignComp(ex, alParameters = erah::setAlPar(
    min.spectra.cor = align.min_cor,
    max.time.dist   = align.max_rt_diff,
    mz.range        = align.mz_range
  ))
  elapsed_align <- as.numeric(difftime(Sys.time(), t0_align, units = "secs"))
  message(sprintf("erah alignment: %.0f sec", elapsed_align))

  # ---- Convert aligned erah output to NTSminer format ----
  # After alignment:
  #   ex@Results@Alignment: data.frame with AlignID, Factor, tmean, FoundIn,
  #     Spectra (consensus MSP string), and per-sample columns (Peak Height).
  #   ex@Data@FactorList: per-sample data.frames with AlignID column added;
  #     columns: ID, RT, Area, Peak.Height, Spectra, Profile, AlignID
  align_res <- ex@Results@Alignment
  if (!is.data.frame(align_res) || !nrow(align_res)) {
    warning("No aligned compounds found", call. = FALSE)
    return(project)
  }

  factor_list <- ex@Data@FactorList
  sample_names <- names(factor_list)

  spec_list <- list()
  ft_list <- list()

  for (i in seq_len(nrow(align_res))) {
    aid <- align_res$AlignID[i]
    compound_id <- sprintf("EA%04d", aid)

    sp <- .nts_parse_erah_spectrum(align_res$Spectra[i])
    if (is.null(sp) || nrow(sp) < 3) next

    bp_idx <- which.max(sp$intensity)
    rt_sec <- align_res$tmean[i] * 60
    found_in <- as.integer(align_res$FoundIn[i])

    spec_list[[compound_id]] <- sp

    ft_row <- list(
      name     = compound_id,
      mz       = sp$mz[bp_idx],
      rt       = rt_sec,
      found_in = found_in
    )

    # Per-sample chromatographic peak Area from FactorList
    for (sname in sample_names) {
      fl <- factor_list[[sname]]
      area_val <- NA_real_
      if (is.data.frame(fl) && nrow(fl) && "AlignID" %in% colnames(fl)) {
        match_rows <- which(fl[["AlignID"]] == aid)
        if (length(match_rows)) {
          a <- fl[["Area"]][match_rows[1L]]
          if (!is.na(a) && a > 0) area_val <- as.numeric(a)
        }
      }
      ft_row[[sname]] <- area_val
    }

    ft_list[[length(ft_list) + 1L]] <- data.table::as.data.table(ft_row)
  }

  if (!length(ft_list)) {
    warning("No valid aligned compounds after parsing", call. = FALSE)
    return(project)
  }

  project@gc_spectra <- spec_list

  ft <- data.table::rbindlist(ft_list, fill = TRUE)
  project@feature_table <- .nts_round_feature_intensities(as.data.frame(ft))

  ms1 <- data.frame(
    name = ft$name,
    mz   = ft$mz,
    rt   = ft$rt,
    stringsAsFactors = FALSE
  )
  project@ms1_annotations <- .nts_ms1_partition(.nts_ms1_sanitize_and_order(ms1))

  message(sprintf(
    "erah GC detection: %d aligned compounds from %d file(s) (deconv: %.0f sec, align: %.0f sec)",
    nrow(ft), length(mzml_files), elapsed, elapsed_align
  ))

  validObject(project)
  project
}


# Internal: parse erah spectrum string "mz,int mz,int ..."
#' @noRd
#' @keywords internal
.nts_parse_erah_spectrum <- function(s) {
  if (is.na(s) || is.null(s) || s == "" || nchar(s) < 5) return(NULL)
  tokens <- tryCatch(strsplit(s, "\\s+")[[1]], error = function(e) NULL)
  if (is.null(tokens) || length(tokens) < 3) return(NULL)

  mz <- numeric(length(tokens))
  int <- integer(length(tokens))
  k <- 0L

  for (i in seq_along(tokens)) {
    p <- tryCatch(strsplit(tokens[i], ",")[[1]], error = function(e) character(0))
    if (length(p) != 2) next
    mi <- suppressWarnings(as.numeric(p[1]))  # mz (preserve high-res)
    ii <- suppressWarnings(as.integer(p[2]))  # int (second)
    if (is.na(mi) || is.na(ii)) next
    if (mi > 0 && ii > 0 && mi < 1000) {
      k <- k + 1L
      mz[k] <- mi
      int[k] <- ii
    }
  }
  if (k < 3) return(NULL)
  if (k < length(tokens) * 0.3) {
    warning("Low parse yield (", k, "/", length(tokens), " tokens); format may have changed", call. = FALSE)
  }
  mx <- max(int[1:k])
  data.frame(
    mz        = mz[1:k],
    intensity = round(int[1:k] / mx * 100, 1),
    stringsAsFactors = FALSE
  )
}


#'' Replace erah integer-mass spectra with high-res spectra from mzML
#'
#' erah deconvolution bins GC-EI spectra to nominal (integer) masses.
#' Re-extracts high-resolution spectra from original mzML at apex RT.
#'
#' @param project An [NTSProject] with @gc_spectra and @raw_data_path.
#' @param rt_window Numeric, RT half-window in seconds. Default 1.0.
#' @param top_n_scans Integer, scans to co-add. Default 1.
#' @return The project with updated high-res spectra.
#' @importFrom mzR openMSfile header peaks close
#' @export
refine_gc_spectra <- function(project, rt_window = 1.0, top_n_scans = 1L) {
  if (!requireNamespace("mzR", quietly = TRUE)) stop("mzR is required", call. = FALSE)
  .nts_check(project, "project", .arg = "project")
  gc <- project@gc_spectra
  if (!length(gc)) stop("@gc_spectra is empty", call. = FALSE)
  ft <- data.table::as.data.table(project@feature_table)
  sample_cols <- setdiff(names(ft), c("name", "mz", "rt", "found_in"))
  n_total <- length(gc); n_replaced <- 0L; n_skipped <- 0L
  mzml_cache <- list()
  on.exit({ for (ms in mzml_cache) try(mzR::close(ms), silent = TRUE) }, add = TRUE)
  message(sprintf("Refining %d spectra...", n_total))
  for (nm in names(gc)) {
    sp_old <- gc[[nm]]
    if (!is.data.frame(sp_old) || nrow(sp_old) < 3) { n_skipped <- n_skipped + 1L; next }
    ft_row <- ft[name == nm]
    if (!nrow(ft_row)) { n_skipped <- n_skipped + 1L; next }
    rt_sec <- ft_row$rt[1L]
    if (is.na(rt_sec)) { n_skipped <- n_skipped + 1L; next }
    target_sample <- NULL
    for (sc in sample_cols) { if (!is.na(ft_row[[sc]]) && ft_row[[sc]] > 0) { target_sample <- sc; break } }
    if (is.null(target_sample)) { n_skipped <- n_skipped + 1L; next }
    mzml_file <- file.path(project@raw_data_path, paste0(target_sample, ".mzML"))
    if (!file.exists(mzml_file)) { n_skipped <- n_skipped + 1L; next }
    if (is.null(mzml_cache[[target_sample]])) mzml_cache[[target_sample]] <- mzR::openMSfile(mzml_file)
    ms <- mzml_cache[[target_sample]]
    hdr <- tryCatch(data.table::as.data.table(mzR::header(ms)), error = function(e) NULL)
    if (is.null(hdr)) { n_skipped <- n_skipped + 1L; next }
    rt_col <- if ("retentionTime" %in% names(hdr)) "retentionTime" else "rt"
    hdr[, delta := abs(get(rt_col) - rt_sec)]
    hdr <- hdr[delta <= rt_window]
    if (!nrow(hdr)) { n_skipped <- n_skipped + 1L; next }
    data.table::setorder(hdr, delta)
    best_scans <- head(hdr, top_n_scans)
    all_peaks <- lapply(best_scans$seqNum, function(sn) {
      p <- tryCatch(mzR::peaks(ms, sn), error = function(e) NULL)
      if (is.null(p) || !is.matrix(p)) return(NULL)
      data.table::data.table(mz = p[, "mz"], intensity = p[, "intensity"])
    })
    all_peaks <- all_peaks[!vapply(all_peaks, is.null, logical(1L))]
    if (!length(all_peaks)) { n_skipped <- n_skipped + 1L; next }
    merged <- data.table::rbindlist(all_peaks)
    merged[, mz_bin := round(mz, 6)]
    merged <- merged[, .(intensity = max(intensity)), by = mz_bin]
    data.table::setnames(merged, "mz_bin", "mz")
    data.table::setorder(merged, mz)
    if (nrow(merged) < 3) { n_skipped <- n_skipped + 1L; next }
    merged[, intensity := round(intensity / max(intensity) * 100, 1)]
    merged <- merged[intensity >= 0.5]
    if (nrow(merged) < 3) { n_skipped <- n_skipped + 1L; next }
    gc[[nm]] <- as.data.frame(merged)
    bp_idx <- which.max(merged$intensity)
    ft[name == nm, mz := merged$mz[bp_idx]]
    n_replaced <- n_replaced + 1L
  }
  project@feature_table <- .nts_round_feature_intensities(as.data.frame(ft))
  project@gc_spectra <- gc
  message(sprintf("Refined: %d / %d", n_replaced, n_total))
  project
}

#'' Build NTSDatabase from NIST MSP export files
#'
#' Reads NIST MSP files with NIST# and InChIKey extraction.
#' @param msp_dir Directory containing NIST MSP files.
#' @param pattern File pattern. Default \code{"\\.msp$"}.
#' @param polarity Ionization polarity. Default \code{"EI"}.
#' @param compress Compression for saveRDS. Default \code{"xz"}.
#' @param output_rds Optional RDS output path.
#' @return An NTSDatabase.
#' @export
build_nist_library <- function(msp_dir, pattern = "\\.msp$", polarity = "EI",
                                compress = "xz", output_rds = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)
  msp_files <- list.files(msp_dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (!length(msp_files)) stop("No MSP files found in: ", msp_dir, call. = FALSE)
  message(sprintf("Found %d MSP files", length(msp_files)))
  meta_list <- list(); spec_list <- list(); counter <- 0L
  for (i in seq_along(msp_files)) {
    f <- msp_files[i]
    message(sprintf("[%d/%d] %s", i, length(msp_files), basename(f)))
    res <- tryCatch(import_msp(f, polarity = polarity), error = function(e) { message("  ERROR"); return(NULL) })
    if (is.null(res) || !nrow(res$metadata)) next
    n <- nrow(res$metadata)
    msp_lines <- readLines(f, warn = FALSE)
    nist_lines <- grep("^NIST#\\s*:", msp_lines, value = TRUE, ignore.case = TRUE)
    nist_nums <- suppressWarnings(as.integer(sub("^NIST#\\s*:\\s*", "", nist_lines, ignore.case = TRUE)))
    syn_lines <- grep("^Synonym\\s*:", msp_lines, value = TRUE, ignore.case = TRUE)
    inchikeys <- rep(NA_character_, n)
    for (j in seq_len(min(n, length(syn_lines)))) {
      ik <- regmatches(syn_lines[j], regexpr("\\$:[0-9]*[A-Z]{14}-[A-Z]{10}-[A-Z]", syn_lines[j]))
      if (length(ik) > 0 && nchar(ik) > 0) inchikeys[j] <- sub("^\\$:[0-9]*", "", ik)
    }
    new_ids <- sprintf("NIST23_%06d", counter + seq_len(n))
    id_map <- setNames(new_ids, res$metadata$Lab_ID)
    res$metadata[, Lab_ID := id_map[Lab_ID]]
    res$spectra_long[, Lab_ID := id_map[Lab_ID]]
    if (length(nist_nums) == n) res$metadata[, NIST_number := nist_nums] else res$metadata[, NIST_number := NA_integer_]
    res$metadata[, InChIKey := inchikeys]
    counter <- counter + n
    meta_list[[length(meta_list) + 1L]] <- res$metadata
    spec_list[[length(spec_list) + 1L]] <- res$spectra_long
    message(sprintf("  %d compounds (total %d)", n, counter))
  }
  if (!length(meta_list)) stop("No valid compounds", call. = FALSE)
  meta_all <- data.table::rbindlist(meta_list, fill = TRUE)
  spec_all <- data.table::rbindlist(spec_list, fill = TRUE)
  data.table::setnames(meta_all, "Lab_ID", "lab_id", skip_absent = TRUE)
  data.table::setnames(spec_all, "Lab_ID", "lab_id", skip_absent = TRUE)
  nb_m <- nrow(meta_all); nb_s <- nrow(spec_all)
  if ("CAS" %in% names(meta_all)) meta_all <- unique(meta_all, by = c("name", "CAS"))
  spec_all <- unique(spec_all, by = c("lab_id", "mz", "intensity"))
  if (nrow(meta_all) < nb_m || nrow(spec_all) < nb_s)
    message(sprintf("Dedup: meta %d->%d, spec %d->%d", nb_m, nrow(meta_all), nb_s, nrow(spec_all)))
  n_tot <- nrow(meta_all)
  n_nist <- sum(!is.na(meta_all$NIST_number))
  n_iki <- sum(!is.na(meta_all$InChIKey) & nzchar(meta_all$InChIKey))
  med <- as.integer(median(spec_all[, .N, by = lab_id]$N, na.rm = TRUE))
  message(sprintf("NIST Library: %d cpds, %d peaks, median %d ions", n_tot, nrow(spec_all), med))
  message(sprintf("  NIST# %d%%, InChIKey %d%%", 100*n_nist/n_tot, 100*n_iki/n_tot))
  lib <- methods::new("NTSDatabase", metadata = meta_all, spectra_long = spec_all,
    features_long = data.table::data.table(), version = "1.0.0")
  if (!is.null(output_rds)) { saveRDS(lib, output_rds, compress = compress); message("Saved: ", output_rds) }
  invisible(lib)
}

#'' Build GC-QTOF RT-indexed library from NIST + measured data
#'
#' Combines NIST metadata with QTOF RT and spectra.
#' @param project NTSProject with refined @gc_spectra.
#' @param nist_lib NTSDatabase from build_nist_library.
#' @param cas_list CAS numbers to include.
#' @param min_cosine Minimum cosine. Default 0.7.
#' @param mz_tol m/z tolerance. Default 0.5.
#' @return NTSDatabase with NIST metadata + RT + QTOF spectra.
#' @export
build_gc_rt_library <- function(project, nist_lib, cas_list, min_cosine = 0.7, mz_tol = 0.5) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)
  .nts_check(project, "project", .arg = "project")
  gc <- project@gc_spectra
  if (!length(gc)) stop("@gc_spectra is empty", call. = FALSE)
  ft <- data.table::as.data.table(project@feature_table)
  cas_vec <- unique(as.character(cas_list))
  cas_vec <- cas_vec[!is.na(cas_vec) & nzchar(cas_vec) & cas_vec != "-"]
  nist_meta <- data.table::as.data.table(nist_lib@metadata)
  if (!"CAS" %in% names(nist_meta)) stop("nist_lib@metadata must have CAS column", call. = FALSE)
  nist_sub_meta <- nist_meta[CAS %in% cas_vec]
  if (!nrow(nist_sub_meta)) stop("No CAS matches in NIST", call. = FALSE)
  nist_sub_spec <- nist_lib@spectra_long[nist_lib@spectra_long$lab_id %in% nist_sub_meta$lab_id, ]
  nist_sub <- methods::new("NTSDatabase", metadata = nist_sub_meta, spectra_long = nist_sub_spec,
    features_long = data.table::data.table(), version = "1.0.0")
  message(sprintf("Sub-library: %d cpds from %d CAS", nrow(nist_sub_meta), length(cas_vec)))
  scores <- score_ei_spectra(project, nist_sub, mz_tol = mz_tol, min_score = min_cosine,
    method = "cosine", prefilter_threshold = 1, best_only = FALSE)
  if (!nrow(scores)) { warning("No matches", call. = FALSE); return(methods::new("NTSDatabase",
    metadata = data.table::data.table(), spectra_long = data.table::data.table(),
    features_long = data.table::data.table(), version = "1.0.0")) }
  ft_rt <- setNames(ft$rt / 60, ft$name)
  lib_meta <- list(); lib_spec <- list(); n_inc <- 0L
  scores[, gid := .GRP, by = ref_id]
  ## NB: the loop variable must NOT be named `gid` -- inside `scores[i]` the
  ## data.table column `gid` would shadow it, making `gid == gid` always TRUE
  ## and returning every row on every iteration (all entries then take the
  ## globally-best reference). Use a distinct name and an explicit join key.
  for (grp in unique(scores$gid)) {
    group <- scores[gid == grp][order(-cosine)][1]
    feat <- group$feature_id; ref_id <- group$ref_id
    sp <- gc[[feat]]
    if (!is.data.frame(sp) || nrow(sp) < 3) next
    meta_row <- nist_sub_meta[lab_id == ref_id]
    if (!nrow(meta_row)) next
    n_inc <- n_inc + 1L; lid <- sprintf("NRT_%04d", n_inc)
    lib_meta[[n_inc]] <- data.table::data.table(lab_id = lid, name = meta_row$name[1L],
      Formula = meta_row$Formula[1L] %||% NA_character_, CAS = meta_row$CAS[1L] %||% NA_character_,
      NIST_number = meta_row$NIST_number[1L] %||% NA_integer_,
      mass = meta_row$mass[1L] %||% NA_real_, RT = ft_rt[feat],
      Polarity = "EI", match_cosine = group$cosine)
    lib_spec[[n_inc]] <- data.table::data.table(lab_id = lid, mz = sp$mz, intensity = sp$intensity)
  }
  if (!n_inc) { warning("No valid spectra", call. = FALSE); return(methods::new("NTSDatabase",
    metadata = data.table::data.table(), spectra_long = data.table::data.table(),
    features_long = data.table::data.table(), version = "1.0.0")) }
  meta_all <- data.table::rbindlist(lib_meta)
  spec_all <- data.table::rbindlist(lib_spec)
  message(sprintf("GC-RT Library: %d cpds, %d peaks, RT %d/%d",
    nrow(meta_all), nrow(spec_all), sum(!is.na(meta_all$RT)), nrow(meta_all)))
  methods::new("NTSDatabase", metadata = meta_all, spectra_long = spec_all,
    features_long = data.table::data.table(), version = "1.0.0")
}

#'' Extract fragment-ion EIC traces for a GC-EI compound
#'
#' @param project NTSProject.
#' @param compound_id Compound identifier.
#' @param sample Sample name (auto-select if NULL).
#' @param mz_tol m/z tolerance. Default 0.02.
#' @param rt_window RT half-window in minutes. Default 0.5.
#' @param top_n Top-N ions. Default 10.
#' @return data.table with columns mz, rt, intensity, relative, ion_rank.
#' @importFrom mzR openMSfile header peaks close
#' @export
extract_gc_fragment_eic <- function(project, compound_id, sample = NULL,
                                     mz_tol = 0.02, rt_window = 0.5, top_n = 10L) {
  if (!requireNamespace("mzR", quietly = TRUE)) stop("mzR is required", call. = FALSE)
  .nts_check(project, "project", .arg = "project")
  gc <- project@gc_spectra
  if (!length(gc)) stop("@gc_spectra is empty", call. = FALSE)
  sp <- gc[[compound_id]]
  if (is.null(sp) || !is.data.frame(sp) || !nrow(sp)) stop("Compound not found", call. = FALSE)
  ft <- data.table::as.data.table(project@feature_table)
  ft_row <- ft[name == compound_id]
  if (!nrow(ft_row)) stop("Not in feature_table", call. = FALSE)
  rt_sec <- ft_row$rt[1L]
  sample_cols <- setdiff(names(ft), c("name", "mz", "rt", "found_in"))
  if (is.null(sample)) {
    for (sc in sample_cols) { if (!is.na(ft_row[[sc]]) && ft_row[[sc]] > 0) { sample <- sc; break } }
    if (is.null(sample)) stop("Not detected in any sample", call. = FALSE)
  }
  mzml_file <- file.path(project@raw_data_path, paste0(sample, ".mzML"))
  if (!file.exists(mzml_file)) stop("mzML not found", call. = FALSE)
  sp_sorted <- sp[order(-sp$intensity), , drop = FALSE]
  if (is.finite(top_n) && nrow(sp_sorted) > top_n) sp_sorted <- sp_sorted[1L:top_n, ]
  ms <- mzR::openMSfile(mzml_file); on.exit(mzR::close(ms), add = TRUE)
  hdr <- data.table::as.data.table(mzR::header(ms))
  hdr[, rt_min := retentionTime / 60]
  rt_lo <- max(0, rt_sec/60 - rt_window); rt_hi <- rt_sec/60 + rt_window
  all_eics <- list(); ei <- 0L
  for (i in seq_len(nrow(sp_sorted))) {
    tgt <- sp_sorted$mz[i]
    sw <- hdr[rt_min >= rt_lo & rt_min <= rt_hi]
    if (!nrow(sw)) next
    rows <- lapply(sw$seqNum, function(sn) {
      p <- tryCatch(mzR::peaks(ms, sn), error = function(e) NULL)
      if (is.null(p) || !is.matrix(p)) return(NULL)
      ir <- p[, "mz"] >= (tgt - mz_tol) & p[, "mz"] <= (tgt + mz_tol)
      if (!any(ir)) return(NULL)
      data.table::data.table(mz = tgt, rt = hdr[seqNum == sn, rt_min], intensity = sum(p[ir, "intensity"]))
    })
    dt <- data.table::rbindlist(rows)
    if (!nrow(dt)) next
    dt[, relative := round(intensity / max(intensity) * 100, 1)]
    dt[, ion_rank := i]
    ei <- ei + 1L; all_eics[[ei]] <- dt
  }
  if (!ei) stop("No EIC data", call. = FALSE)
  out <- data.table::rbindlist(all_eics)
  data.table::setorder(out, ion_rank, rt); out[]
}


#'' Plot fragment-ion EICs for a GC-EI compound
#'
#' @param project NTSProject.
#' @param compound_id Compound identifier.
#' @param ... Passed to extract_gc_fragment_eic.
#' @param overlay Logical. Default FALSE.
#' @return Invisibly, EIC data.table.
#' @export
plot_gc_fragment_eic <- function(project, compound_id, ..., overlay = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
  eic <- extract_gc_fragment_eic(project, compound_id, ...)
  eic[, mz_label := sprintf("m/z %.4f (r%d)", mz, ion_rank)]
  if (overlay) {
    p <- ggplot2::ggplot(eic, ggplot2::aes(x = rt, y = relative, color = mz_label, group = mz_label)) +
      ggplot2::geom_line(linewidth = 0.6) + ggplot2::scale_color_viridis_d() +
      ggplot2::labs(title = compound_id, x = "RT (min)", y = "%", color = "Ion") + ggplot2::theme_minimal()
  } else {
    p <- ggplot2::ggplot(eic, ggplot2::aes(x = rt, y = relative)) +
      ggplot2::geom_line(linewidth = 0.5, color = "#2563eb") + ggplot2::facet_wrap(~ mz_label, ncol = 2, scales = "free_y") +
      ggplot2::labs(title = compound_id, x = "RT (min)", y = "%") + ggplot2::theme_minimal()
  }
  print(p); invisible(eic)
}