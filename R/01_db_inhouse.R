#' @keywords internal
#' @noRd
.nts_normalize_mzr_header_names <- function(h) {
  dt <- data.table::as.data.table(h)
  nm <- names(dt)
  map <- c(
    precursorMZ = "precursorMZ",
    precursorMz = "precursorMZ",
    precursor_mz = "precursorMZ",
    precursorIntensity = "precursorIntensity",
    precursorintensity = "precursorIntensity",
    retentionTime = "retentionTime",
    retentiontime = "retentionTime",
    acquisitionNum = "acquisitionNum",
    acquisitionnum = "acquisitionNum",
    msLevel = "msLevel",
    mslevel = "msLevel"
  )
  hit <- nm %in% names(map)
  if (any(hit)) {
    newn <- nm
    newn[hit] <- unname(map[nm[hit]])
    data.table::setnames(dt, newn)
  }
  dt
}

#' Extract MS2 scans for mixed standards from mzML
#'
#' Reads MS2-level entries from an mzML file via **mzR**, matches each compound's
#' theoretical precursor m/z (from neutral mass and polarity), clusters matched
#' scans in retention time to separate peak groups (e.g. isomers), and records
#' the highest-`precursorIntensity` scan per group. The global best scan (max
#' group apex intensity) is written to `best_*` columns. If a second peak group
#' has apex intensity greater than 10% of the strongest group, `isomer_warning`
#' is set to `TRUE` and representative RTs of all groups are listed in
#' `isomer_rts`.
#'
#' Intended as a front-end step before curating in-house reference MS2 for
#' [compile_reference_database()] or related workflows.
#'
#' Requires the **mzR** package ([Bioconductor](https://bioconductor.org/packages/mzR/));
#' install with `BiocManager::install("mzR")` if missing.
#'
#' @param cmpd_list `data.table` or `data.frame` with columns `Lab_ID`, `name`,
#'   and `mass` (exact neutral mass).
#' @param mzml_file Character path(s) to mzML/mzXML file(s) (all must exist).
#' @param polarity `"NEG"` for `[M-H]-` (`mass - 1.007276`) or `"POS"` for
#'   `[M+H]+` (`mass + 1.007276`).
#' @param mz_tol_ppm Precursor m/z tolerance in ppm (default `10`).
#' @param mz_tol_da Optional absolute precursor m/z tolerance in Da. When set,
#'   matching uses the wider of ppm-derived and absolute tolerances.
#' @param rt_gap_tol Minimum retention-time gap (same units as in the file,
#'   typically seconds for mzML) to start a new peak group (default `15`).
#'
#' @return A `data.table` equal to `cmpd_list` with added columns:
#'   `Polarity` (the normalised `"POS"` / `"NEG"` ion mode, so the result can be
#'   passed straight to [build_nts_database()] without the polarity being
#'   defaulted), `target_mz`, `best_rt`, `best_scan_num` (`mzR::header()`
#'   `acquisitionNum`), `best_source_file`, `best_intensity`, `isomer_warning`,
#'   `isomer_rts`.
#'
#' @export
extract_mix_standards <- function(
    cmpd_list,
    mzml_file,
    polarity = "NEG",
    mz_tol_ppm = 10,
    mz_tol_da = NULL,
    rt_gap_tol = 15) {
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("Package \"mzR\" is required. Install with BiocManager::install(\"mzR\").", call. = FALSE)
  }

  req <- c("Lab_ID", "name", "mass")
  cn <- names(cmpd_list)
  miss <- setdiff(req, cn)
  if (length(miss)) {
    stop("`cmpd_list` must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
  }

  pol <- toupper(as.character(polarity)[1L])
  if (!pol %in% c("NEG", "POS")) {
    stop("`polarity` must be \"NEG\" or \"POS\".", call. = FALSE)
  }

  mzml_files <- unique(path.expand(as.character(mzml_file)))
  mzml_files <- mzml_files[nzchar(trimws(mzml_files))]
  if (!length(mzml_files)) {
    stop("`mzml_file` must contain at least one file path.", call. = FALSE)
  }
  miss_f <- mzml_files[!file.exists(mzml_files)]
  if (length(miss_f)) {
    stop("Missing mzML/mzXML file(s): ", paste(miss_f, collapse = ", "), call. = FALSE)
  }

  proton <- 1.007276
  out <- data.table::copy(data.table::as.data.table(cmpd_list))
  out[["mass"]] <- as.numeric(out[["mass"]])
  if (any(!is.finite(out[["mass"]]))) {
    stop("`cmpd_list$mass` must be finite numeric.", call. = FALSE)
  }
  dmz <- if (pol == "NEG") -proton else proton
  nm_low <- tolower(names(out))
  pick_col <- function(cands) {
    hit <- match(tolower(cands), nm_low)
    idx <- hit[!is.na(hit)][1L]
    if (is.na(idx)) return(NULL)
    names(out)[idx]
  }
  
  target_col <- pick_col(c("target_mz", "target.mz"))
  if (is.null(target_col)) {
    target_col <- if (identical(pol, "NEG")) {
      pick_col(c("mz.neg", "mz_neg", "mzneg", "mz_negative", "mznegative"))
    } else {
      pick_col(c("mz.pos", "mz_pos", "mzpos", "mz_positive", "mzpositive"))
    }
  }
  if (!is.null(target_col)) {
    out[["target_mz"]] <- suppressWarnings(as.numeric(out[[target_col]]))
    bad <- !is.finite(out[["target_mz"]]) | out[["target_mz"]] <= 0
    if (any(bad)) {
      out[["target_mz"]][bad] <- out[["mass"]][bad] + dmz
    }
  } else {
    out[["target_mz"]] <- out[["mass"]] + dmz
  }

  rt_col <- pick_col(c("rt", "retentiontime", "retention_time", "tr"))
  if (!is.null(rt_col)) {
    out[["target_rt"]] <- suppressWarnings(as.numeric(out[[rt_col]]))
  }

  need <- c("msLevel", "precursorMZ", "retentionTime", "acquisitionNum", "precursorIntensity", "seqNum")
  ms2_chunks <- vector("list", length(mzml_files))
  n_ms2_total <- 0L
  for (i in seq_along(mzml_files)) {
    f <- mzml_files[i]
    message(sprintf("Opening file [%d/%d]: %s", i, length(mzml_files), basename(f)))
    ms <- mzR::openMSfile(f)
    # Using tryCatch for mzR::close inside loop to avoid on.exit pileup
    h <- mzR::header(ms)
    
    if (is.null(h) || !nrow(h)) {
      try(mzR::close(ms), silent = TRUE)
      next
    }
    hdt <- .nts_normalize_mzr_header_names(h)
    miss_h <- setdiff(need, names(hdt))
    if (length(miss_h)) {
      try(mzR::close(ms), silent = TRUE)
      stop("mzML header is missing required columns after normalization in file ", basename(f), ": ", paste(miss_h, collapse = ", "), call. = FALSE)
    }
    ms2 <- hdt[as.integer(hdt[["msLevel"]]) == 2L]
    n_ms2 <- nrow(ms2)
    n_zero_pmz <- sum(ms2[["precursorMZ"]] == 0, na.rm = TRUE)
    message("  [DIAG] Total MS2 spectra: ", n_ms2)
    message("  [DIAG] precursorMZ == 0 count: ", n_zero_pmz)
    if (n_ms2 > n_zero_pmz) {
      non_zero_pmz <- ms2[["precursorMZ"]][ms2[["precursorMZ"]] != 0 & is.finite(ms2[["precursorMZ"]])]
      if (length(non_zero_pmz) > 0) {
        # Deterministic diagnostics: take the first 3 values instead of
        # `sample()`. This is only a `[DIAG]` message, but `sample()` would
        # perturb the global RNG stream, silently breaking a user's `set.seed()`
        # reproducibility for downstream modelling.
        smpl <- utils::head(non_zero_pmz, 3L)
        message("  [DIAG] Sample non-zero precursorMZ: ", paste(smpl, collapse = ", "))
      }
    }
    if (!n_ms2) {
      try(mzR::close(ms), silent = TRUE)
      next
    }
    ms2[["source_file"]] <- f
    if ("collisionEnergy" %in% names(ms2)) {
      ce_val <- suppressWarnings(as.numeric(ms2[["collisionEnergy"]]))
      ms2[["nce_level"]] <- ifelse(is.finite(ce_val) & ce_val > 0, paste0("CE", round(ce_val)), NA_character_)
    } else {
      ms2[["nce_level"]] <- NA_character_
    }

    bad_nce <- is.na(ms2[["nce_level"]])
    if (any(bad_nce)) {
      nce_match <- regexpr("(?i)(C|NC)E\\d+", basename(f), perl = TRUE)
      fallback_val <- if (nce_match > 0) toupper(regmatches(basename(f), nce_match)) else "CE_Unknown"
      ms2[["nce_level"]][bad_nce] <- fallback_val
    }
    ms2[["orig_pmz"]] <- ms2[["precursorMZ"]]
    ms2[["orig_rt"]] <- ms2[["retentionTime"]]
    ms2_chunks[[i]] <- ms2
    n_ms2_total <- n_ms2_total + nrow(ms2)
    try(mzR::close(ms), silent = TRUE)
  }
  if (!n_ms2_total) {
    warning("No MS2 header rows found across input files.", call. = FALSE)
  }
  hd_ms2 <- data.table::rbindlist(ms2_chunks, use.names = TRUE, fill = TRUE)
  message("MS2 spectra in headers: ", nrow(hd_ms2), " across ", length(mzml_files), " files.")

  ncmp <- nrow(out)
  message("Compounds in list: ", ncmp)

  out_ids <- out[, .(Lab_ID, target_mz)]
  if ("target_rt" %in% names(out)) {
    out_ids[, target_rt := out[["target_rt"]]]
  } else {
    out_ids[, target_rt := NA_real_]
  }

  tol_da_val <- if (is.null(mz_tol_da) || !is.finite(mz_tol_da)) NA_real_ else as.numeric(mz_tol_da)
  out_ids[, mz_tol_abs := max(target_mz * mz_tol_ppm / 1e6, tol_da_val, na.rm = TRUE), by = seq_len(nrow(out_ids))]
  out_ids[, mz_min := target_mz - mz_tol_abs]
  out_ids[, mz_max := target_mz + mz_tol_abs]

  has_rt <- "target_rt" %in% names(out_ids) && any(is.finite(out_ids$target_rt))
  if (has_rt) {
    out_ids[, rt_min := target_rt - rt_gap_tol]
    out_ids[, rt_max := target_rt + rt_gap_tol]
    message(sprintf("[DIAG] Target RT Range: %.2f - %.2f", min(out_ids$target_rt, na.rm=TRUE), max(out_ids$target_rt, na.rm=TRUE)))
    message(sprintf("[DIAG] mzXML RT Range: %.2f - %.2f", min(hd_ms2$retentionTime, na.rm=TRUE), max(hd_ms2$retentionTime, na.rm=TRUE)))
    message("Performing 2D (MZ + RT) Non-equi Join...")
    matched_hits <- hd_ms2[out_ids, 
                           on = .(precursorMZ >= mz_min, precursorMZ <= mz_max, 
                                  retentionTime >= rt_min, retentionTime <= rt_max),
                           nomatch = NULL, allow.cartesian = TRUE]
  } else {
    message("No valid RT column found. Performing 1D (MZ) Non-equi Join...")
    matched_hits <- hd_ms2[out_ids, 
                           on = .(precursorMZ >= mz_min, precursorMZ <= mz_max),
                           nomatch = NULL, allow.cartesian = TRUE]
  }

  if (nrow(matched_hits) > 0) {
    data.table::setorderv(matched_hits, "precursorIntensity", -1L, na.last = TRUE)
    best_hits <- matched_hits[, .SD[1L], by = .(Lab_ID, nce_level)]
  } else {
    best_hits <- data.table::data.table(Lab_ID = character(), nce_level = character(), source_file = character(), acquisitionNum = integer())
  }
  
  if (!"diagnostic_log" %in% names(best_hits)) {
     best_hits[, diagnostic_log := ""]
  }

  # Fallback for 0.5 Da
  missing_ids <- setdiff(out_ids$Lab_ID, best_hits$Lab_ID)
  if (length(missing_ids) > 0) {
    message("Fallback (0.5 Da window) for ", length(missing_ids), " unassigned targets...")
    fallback_ids <- out_ids[Lab_ID %in% missing_ids]
    fallback_ids[, mz_min := target_mz - 0.5]
    fallback_ids[, mz_max := target_mz + 0.5]
    
    if (has_rt) {
      fallback_hits <- hd_ms2[fallback_ids, 
                              on = .(precursorMZ >= mz_min, precursorMZ <= mz_max, 
                                     retentionTime >= rt_min, retentionTime <= rt_max),
                              nomatch = NULL, allow.cartesian = TRUE]
    } else {
      fallback_hits <- hd_ms2[fallback_ids, 
                              on = .(precursorMZ >= mz_min, precursorMZ <= mz_max),
                              nomatch = NULL, allow.cartesian = TRUE]
    }
    
    if (nrow(fallback_hits) > 0) {
      data.table::setorderv(fallback_hits, "precursorIntensity", -1L, na.last = TRUE)
      best_fallback <- fallback_hits[, .SD[1L], by = .(Lab_ID, nce_level)]
      best_fallback[, diagnostic_log := "Dependent on wide window matching (0.5 Da)"]
      best_hits <- rbind(best_hits, best_fallback, fill = TRUE)
    }
  }

  # Extract Peaks & Hardcore Background Purification
  best_hits[, peak_data := list()]
  if (nrow(best_hits) > 0) {
    unique_files <- unique(best_hits$source_file)
    for (f in unique_files) {
      f_hits_idx <- which(best_hits$source_file == f)
      if (length(f_hits_idx) == 0) next
      
      ms <- mzR::openMSfile(f)
      pk_list <- mzR::peaks(ms)
      mzR::close(ms)
      
      seqs <- best_hits$seqNum[f_hits_idx]
      
      for (j in seq_along(seqs)) {
        seq_idx <- seqs[j]
        # Ultimate out-of-bounds protection
        if (is.na(seq_idx) || seq_idx < 1 || seq_idx > length(pk_list)) next
        
        mat <- pk_list[[seq_idx]]
        # Core null-pointer prevention
        if (is.matrix(mat) && nrow(mat) > 0) {
          mz_vec <- as.numeric(mat[, 1L])
          int_vec <- as.numeric(mat[, 2L])
          max_int <- max(int_vec, na.rm = TRUE)
          if (is.finite(max_int) && max_int > 0) {
            rel_int <- int_vec / max_int * 1000
            keep <- rel_int >= 5
            if (any(keep)) {
              df <- data.table::data.table(
                mz = mz_vec[keep], 
                intensity = int_vec[keep], 
                rel_intensity = rel_int[keep]
              )
              best_hits$peak_data[[f_hits_idx[j]]] <- df
            }
          }
        }
      }
    }
  }

  # Structure into Nested List-column
  data.table::setDT(out)
  out[, Polarity := pol]
  out$ms2_spectra <- vector("list", nrow(out))
  out[, diagnostic_log := ""]
  out[, matched_rt := NA_real_]
  out[, matched_intensity := NA_real_]
  out[, rt_drift := NA_real_]
  
  for (k in seq_len(nrow(out))) {
    lab <- out$Lab_ID[k]
    hits <- best_hits[Lab_ID == lab]
    if (nrow(hits) > 0) {
      res_list <- list()
      diag_logs <- c()
      
      max_int_val <- max(hits$precursorIntensity, na.rm = TRUE)
      if (is.finite(max_int_val)) {
        out$matched_intensity[k] <- max_int_val
        max_idx <- which.max(hits$precursorIntensity)
        out$matched_rt[k] <- hits$orig_rt[max_idx]
        if ("target_rt" %in% names(out)) {
          out$rt_drift[k] <- out$matched_rt[k] - out$target_rt[k]
        }
      }
      
      for (h in seq_len(nrow(hits))) {
        nce <- hits$nce_level[h]
        pd <- hits$peak_data[[h]]
        if (!is.null(pd) && is.data.frame(pd) && nrow(pd) > 0) {
          res_list[[nce]] <- pd
        }
        if (hits$diagnostic_log[h] != "") {
          diag_logs <- c(diag_logs, hits$diagnostic_log[h])
        }
      }
      if (length(res_list) > 0) {
        out$ms2_spectra[[k]] <- res_list
      }
      if (length(diag_logs) > 0) {
        out$diagnostic_log[k] <- paste(unique(diag_logs), collapse = "; ")
      }
    }
  }

  message("Done. Extracted MS2 for ", sum(lengths(out$ms2_spectra) > 0), " / ", ncmp, " compounds.")
  out[]
}

#' Extract MS2 peak lists for selected scans from mzML
#'
#' Uses **mzR** to read centroid (or profile) peak matrices for spectra whose
#' `acquisitionNum` matches `best_scan_num` in `metadata`. Row order follows
#' `mzR::header()` / `mzR::peaks()` (1-based spectrum index), not necessarily
#' the numeric order of `acquisitionNum`.
#'
#' @param metadata `data.table` or `data.frame` with `Lab_ID` and `best_scan_num`
#'   (as returned by [extract_mix_standards()]).
#' @param mzml_file Path(s) to mzML/mzXML file(s). When `metadata` contains
#'   `best_source_file`, extraction is restricted to that file per row.
#'
#' @return A `data.table` with columns `Lab_ID`, `source_file`, `mz`,
#'   `intensity` (one row per peak). Compounds with missing/invalid scans
#'   contribute no rows.
#'
#' @export
extract_spectra_from_scans <- function(metadata, mzml_file) {
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop(
      "Package \"mzR\" is required for extract_spectra_from_scans(). ",
      "Install with BiocManager::install(\"mzR\").",
      call. = FALSE
    )
  }
  md <- data.table::as.data.table(metadata)
  req <- c("Lab_ID", "best_scan_num")
  miss <- setdiff(req, names(md))
  if (length(miss)) {
    stop(
      "`metadata` must contain columns: ",
      paste(req, collapse = ", "),
      call. = FALSE
    )
  }

  mzml_files <- unique(path.expand(as.character(mzml_file)))
  mzml_files <- mzml_files[nzchar(trimws(mzml_files))]
  if (!length(mzml_files)) {
    stop("`mzml_file` must contain at least one file path.", call. = FALSE)
  }
  miss_f <- mzml_files[!file.exists(mzml_files)]
  if (length(miss_f)) {
    stop("Missing mzML/mzXML file(s): ", paste(miss_f, collapse = ", "), call. = FALSE)
  }

  file_maps <- vector("list", length(mzml_files))
  names(file_maps) <- mzml_files
  # Single function-level safety net: closes whatever handle is currently open
  # if the function exits via an error. Per-iteration handles are closed
  # explicitly below. (The older code used `on.exit(..., add = TRUE)` *inside*
  # the loop, which piled up one deferred close per file and could conflict with
  # the explicit closes, leaking mzR handles on large batches.)
  cur_ms <- NULL
  on.exit(
    if (!is.null(cur_ms)) try(mzR::close(cur_ms), silent = TRUE),
    add = TRUE
  )
  for (i in seq_along(mzml_files)) {
    f <- mzml_files[i]
    message(sprintf("extract_spectra_from_scans: opening [%d/%d] %s", i, length(mzml_files), basename(f)))
    ms <- mzR::openMSfile(f)
    cur_ms <- ms
    hd <- mzR::header(ms)
    if (is.null(hd) || !nrow(hd)) {
      try(mzR::close(ms), silent = TRUE)
      cur_ms <- NULL
      next
    }
    hdt <- .nts_normalize_mzr_header_names(hd)
    if (!"acquisitionNum" %in% names(hdt)) {
      stop("mzML header has no acquisitionNum column in file: ", basename(f), call. = FALSE)
    }
    acq <- suppressWarnings(as.integer(hdt[["acquisitionNum"]]))
    pk <- mzR::peaks(ms)
    # Peaks are now in memory; release the file handle immediately.
    try(mzR::close(ms), silent = TRUE)
    cur_ms <- NULL
    if (!is.list(pk) || !length(pk)) {
      next
    }
    file_maps[[i]] <- list(
      file = f,
      file_base = basename(f),
      acq = acq,
      pk = pk,
      n_spec = length(pk)
    )
    if (length(pk) != nrow(hdt)) {
      message(
        "Note: length(peaks) (",
        length(pk),
        ") differs from nrow(header) (",
        nrow(hdt),
        ") in file ",
        basename(f),
        "; indexing uses header row order."
      )
    }
  }

  n_md <- nrow(md)
  chunks <- vector("list", n_md)
  n_ok <- 0L
  for (k in seq_len(n_md)) {
    lab <- as.character(md[["Lab_ID"]][k])
    sn <- suppressWarnings(as.integer(md[["best_scan_num"]][k]))
    if (!is.finite(sn)) {
      next
    }
    src_pref <- if ("best_source_file" %in% names(md)) as.character(md[["best_source_file"]][k]) else NA_character_
    src_pref <- trimws(src_pref)
    cands <- file_maps
    if (nzchar(src_pref)) {
      cands <- Filter(
        function(x) !is.null(x) && (identical(x$file, src_pref) || identical(x$file_base, basename(src_pref))),
        file_maps
      )
      if (!length(cands)) {
        cands <- file_maps
      }
    }
    idx <- NA_integer_
    pk_hit <- NULL
    src_hit <- NA_character_
    for (m in cands) {
      if (is.null(m)) {
        next
      }
      idx0 <- match(sn, m$acq)
      if (is.na(idx0) || idx0 < 1L || idx0 > m$n_spec) {
        next
      }
      idx <- idx0
      pk_hit <- m$pk
      src_hit <- m$file
      break
    }
    if (!is.finite(idx) || is.null(pk_hit)) {
      warning(
        "No spectrum index for Lab_ID=",
        lab,
        " best_scan_num=",
        sn,
        call. = FALSE,
        immediate. = TRUE
      )
      next
    }
    mat <- pk_hit[[idx]]
    if (is.null(mat)) {
      next
    }
    if (!is.matrix(mat)) {
      mat <- tryCatch(
        as.matrix(mat),
        error = function(e) matrix(numeric(), nrow = 0L, ncol = 2L)
      )
    }
    if (!nrow(mat)) {
      next
    }
    nc <- ncol(mat)
    if (nc < 2L) {
      next
    }
    mz <- as.numeric(mat[, 1L])
    int <- as.numeric(mat[, 2L])
    ok <- is.finite(mz) & is.finite(int)
    if (!any(ok)) {
      next
    }
    n_ok <- n_ok + 1L
    chunks[[n_ok]] <- data.table::data.table(
      Lab_ID = lab,
      source_file = src_hit,
      mz = mz[ok],
      intensity = int[ok]
    )
  }

  if (!n_ok) {
    message("extract_spectra_from_scans: no peaks extracted.")
    return(data.table::data.table(
      Lab_ID = character(),
      source_file = character(),
      mz = numeric(),
      intensity = numeric()
    ))
  }
  out <- data.table::rbindlist(chunks[seq_len(n_ok)], use.names = TRUE)
  message(
    "extract_spectra_from_scans: ",
    nrow(out),
    " peak rows from ",
    n_ok,
    " compounds."
  )
  out[]
}
