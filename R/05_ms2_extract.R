#' Slice an MS2-only experiment for one file, or read it from disk
#'
#' Reuses an already-loaded MSnbase experiment (e.g. `project@xdata`) when it
#' covers `f`, so callers that already read the same files (MS1 / MS2 in one
#' file) do not pay for a second full `readMSData()`. Any failure falls back to
#' the original on-disk read, so the fast path can never change results.
#' @keywords internal
#' @noRd
.nts_ms2_raw_for_file <- function(f, xdata = NULL) {
  if (!is.null(xdata) && methods::is(xdata, "MSnExp")) {
    got <- tryCatch({
      fns <- basename(MSnbase::fileNames(xdata))
      w <- which(fns == basename(f))
      if (length(w)) {
        sub <- MSnbase::filterFile(xdata, w[1L])
        sub <- MSnbase::filterMsLevel(sub, msLevel = 2L)
        if (length(sub)) sub else NULL
      } else {
        NULL
      }
    }, error = function(e) NULL)
    if (!is.null(got)) return(got)
  }
  MSnbase::readMSData(files = f, mode = "onDisk", msLevel = 2L)
}

#' Per-scan precursor isolation purity from the file's own MS1 scans
#'
#' For each MS2 scan, locates the closest preceding MS1 scan (in RT) in the same
#' file and returns the fraction of summed MS1 intensity inside the scan's
#' isolation window that is attributable to the target precursor m/z:
#' `I(precursor +/- mz_tol) / I(isolation window)`. A value near 1 means the
#' precursor was isolated cleanly; low values flag co-isolation (isobaric /
#' isomeric interferences), the main driver of chimeric MS2 in DDA.
#'
#' Returns `NA` for a scan when no MS1 scan is available in the file, the window
#' is unknown, or the in-window intensity is zero. Purely diagnostic: never
#' changes which scans are matched.
#' @keywords internal
#' @noRd
.nts_iso_purity_for_file <- function(f,
                                     scan_rt,
                                     iso_lo,
                                     iso_hi,
                                     prec_mz,
                                     mz_tol_ppm = 15,
                                     xdata = NULL) {
  n <- length(scan_rt)
  out <- rep(NA_real_, n)
  if (!n) return(out)

  ms1 <- NULL
  if (!is.null(xdata) && methods::is(xdata, "MSnExp")) {
    ms1 <- tryCatch({
      fns <- basename(MSnbase::fileNames(xdata))
      w <- which(fns == basename(f))
      if (length(w)) {
        sub <- MSnbase::filterFile(xdata, w[1L])
        sub <- MSnbase::filterMsLevel(sub, msLevel = 1L)
        if (length(sub)) sub else NULL
      } else {
        NULL
      }
    }, error = function(e) NULL)
  }
  if (is.null(ms1)) {
    ms1 <- tryCatch(
      MSnbase::readMSData(files = f, mode = "onDisk", msLevel = 1L),
      error = function(e) NULL
    )
  }
  if (is.null(ms1) || !length(ms1)) return(out)

  ms1_rt <- suppressWarnings(as.numeric(MSnbase::rtime(ms1)))
  ok <- is.finite(ms1_rt)
  if (!any(ok)) return(out)
  ms1 <- ms1[ok]
  ms1_rt <- ms1_rt[ok]
  ord <- order(ms1_rt)
  ms1 <- ms1[ord]
  ms1_rt <- ms1_rt[ord]

  tol <- suppressWarnings(as.numeric(mz_tol_ppm)[1L]) * 1e-6
  if (!is.finite(tol) || tol <= 0) tol <- 15e-6

  for (i in seq_len(n)) {
    rt_i <- scan_rt[i]
    pm_i <- prec_mz[i]
    if (!is.finite(rt_i) || !is.finite(pm_i)) next
    j <- findInterval(rt_i, ms1_rt)
    if (j < 1L) j <- 1L

    sp <- ms1[[j]]
    mzv <- suppressWarnings(as.numeric(MSnbase::mz(sp)))
    iv <- suppressWarnings(as.numeric(MSnbase::intensity(sp)))
    k <- is.finite(mzv) & is.finite(iv) & iv > 0
    if (!any(k)) next
    mzv <- mzv[k]
    iv <- iv[k]

    lo <- iso_lo[i]
    hi <- iso_hi[i]
    win <- if (is.finite(lo) && is.finite(hi) && hi > lo) {
      mzv >= lo & mzv <= hi
    } else {
      rep(TRUE, length(mzv))
    }
    denom <- sum(iv[win])
    if (!is.finite(denom) || denom <= 0) next
    tgt <- abs(mzv - pm_i) <= pm_i * tol
    num <- sum(iv[tgt & win])
    out[i] <- max(0, min(1, num / denom))
  }
  out
}


#' Attach consensus MS2 spectra to an [NTSProject]
#'
#' Reads DDA mzML/mzXML files one-by-one (with progress messages), performs vectorized
#' m/z-range overlap matching (`data.table::foverlaps`) plus RT bounds filter,
#' lazily loads only matched spectra indices, and stores **all** cleaned MS2
#' spectra per feature as a nested list. Peaks below 2.5% relative intensity are
#' removed; retained intensities are scaled to 0-1.
#'
#' Per-scan acquisition metadata is captured into each `@ms2data[[feature]][[i]]`
#' element regardless of `scoring`: `ce` (collision energy), `polarity`
#' (`0 = neg`, `1 = pos`, `NA = unknown`), `isolation_lo` / `isolation_hi`
#' (precursor isolation window edges in Da; `NA` if not stored by the vendor).
#'
#' @param project An [NTSProject] with `ms1_annotations` containing `name`, `mz`,
#'   and `rt`.
#' @param ms2_dir Directory of centroid MS2 files: **mzML** and **mzXML** are
#'   scanned recursively (case-insensitive extensions).
#' @param mz_tol_ppm,rt_tol Matching tolerances (ppm and RT units consistent
#'   with `rt` in `ms1_annotations`).
#' @param scoring One of `"none"` (default, fully backward compatible),
#'   `"soft"` (compute multi-evidence match score per scan, keep all hits),
#'   or `"strict"` (compute score and drop hits below `min_match_score`).
#'   When enabled, each scan element gains `match_score`, `s_mz`, `s_rt`,
#'   `s_iso`, `s_purity`, `s_tic`, `s_npeaks` fields. `s_iso` measures how
#'   centrally the feature m/z sits inside the scan's isolation window;
#'   `s_purity` is the precursor isolation purity read from the file's own MS1
#'   scans (fraction of in-window intensity belonging to the target precursor).
#'   `s_purity` is `NA` when the file exposes no MS1 scans.
#' @param weights Named numeric vector of component weights summed in the
#'   composite `match_score`. NA components are skipped and remaining weights
#'   are renormalised on the fly.
#' @param min_match_score Threshold used when `scoring = "strict"`. Default
#'   `0.4` (weighted score is in `[0, 1]`).
#' @param isolation_width_da Optional fallback isolation window full width
#'   (Da) used only when the vendor file omits per-scan isolation offsets;
#'   `NULL` leaves them as `NA`.
#' @param xdata Optional pre-read MSnbase experiment (e.g. `project@xdata`)
#'   covering the same files. When supplied, its MS2 spectra are sliced per
#'   file instead of re-reading each file with [MSnbase::readMSData]. This
#'   avoids a duplicate disk read when the MS1 and MS2 data live in the same
#'   files. Any mismatch silently falls back to the on-disk read.
#'
#' @return Updated `project`: `@ms2data` keeps **only** features with at least
#'   one non-empty MS2 scan; `@ms2_annotations` holds one row per surviving
#'   feature with `name`, `mz`, `rt`, and `n_fragments` (raw peak count),
#'   joined from `@ms1_annotations` with `feature_table` mz/rt fallback when
#'   needed.
#' @export
extract_ms2 <- function(project,
                        ms2_dir,
                        mz_tol_ppm = 15,
                        rt_tol = 30,
                        scoring = c("none", "soft", "strict"),
                        weights = c(mz = 0.30, rt = 0.20, iso = 0.20,
                                    purity = 0.15, tic = 0.05, npeaks = 0.10),
                        min_match_score = 0.4,
                        isolation_width_da = NULL,
                        xdata = NULL) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  scoring <- match.arg(scoring)
  ms1 <- get_ms1_annotations(project)
  if (!is.data.frame(ms1) || !nrow(ms1)) {
    stop("ms1_annotations is empty", call. = FALSE)
  }
  if (!all(c("name", "mz", "rt") %in% names(ms1))) {
    stop("ms1_annotations must contain columns name, mz, rt", call. = FALSE)
  }

  ft_dt <- data.table::as.data.table(ms1)
  ft_dt <- ft_dt[, .(feature_name = as.character(name), mz = as.numeric(mz), rt = as.numeric(rt))]
  ft_dt <- ft_dt[is.finite(mz) & is.finite(rt)]
  if (!nrow(ft_dt)) {
    stop("ms1_annotations has no finite mz/rt rows", call. = FALSE)
  }
  ppm <- mz_tol_ppm / 1e6
  ft_dt[, `:=`(
    mz_lo = mz * (1 - ppm),
    mz_hi = mz * (1 + ppm),
    rt_lo = rt - rt_tol,
    rt_hi = rt + rt_tol
  )]

  ms2_dir <- normalizePath(ms2_dir, mustWork = TRUE)
  files <- list.files(
    ms2_dir,
    pattern = "\\.(mzml|mzxml)$",
    ignore.case = TRUE,
    full.names = TRUE,
    recursive = TRUE
  )
  if (!length(files)) {
    stop("No mzML/mzXML files found in ms2_dir", call. = FALSE)
  }

  ms2_map <- data.table::data.table(
    feature_name = character(0),
    scan_idx = integer(0),
    source_file = character(0),
    tic = numeric(0),
    precursorMz = numeric(0),
    rt = numeric(0),
    ce = numeric(0),
    polarity = integer(0),
    isolation_lo = numeric(0),
    isolation_hi = numeric(0),
    purity = numeric(0),
    mz = I(list()),
    intensity = I(list())
  )

  n_files <- length(files)
  n_dropped_no_prec <- 0L
  for (i in seq_len(n_files)) {
    f <- files[i]
    message(sprintf("[%d/%d] Parsing and matching MS2: %s", i, n_files, basename(f)))

    raw <- .nts_ms2_raw_for_file(f, xdata)
    n_sp <- length(raw)
    if (!n_sp) {
      next
    }

    pm <- as.numeric(MSnbase::precursorMz(raw))
    rt <- as.numeric(MSnbase::rtime(raw))
    meta <- .nts_extract_scan_meta(
      raw,
      n_sp,
      default_isolation_width_da = isolation_width_da
    )
    idx <- seq_len(n_sp)
    sp_dt <- data.table::data.table(
      scan_idx = idx,
      precursorMz = pm,
      rt = rt
    )
    n_before_pmz <- nrow(sp_dt)
    sp_dt <- sp_dt[is.finite(precursorMz) & is.finite(rt)]
    if (nrow(sp_dt) < n_before_pmz) {
      n_dropped_no_prec <- n_dropped_no_prec + (n_before_pmz - nrow(sp_dt))
    }
    if (!nrow(sp_dt)) {
      next
    }
    sp_dt[, `:=`(mz_lo = precursorMz, mz_hi = precursorMz)]

    sp_ov <- sp_dt[, .(mz_lo, mz_hi, scan_idx, precursorMz, rt)]
    ft_ov <- ft_dt[, .(mz_lo, mz_hi, feature_name, rt_lo, rt_hi)]
    data.table::setkey(ft_ov, mz_lo, mz_hi)

    hits <- data.table::foverlaps(sp_ov, ft_ov, nomatch = 0L, type = "within")
    if (!nrow(hits)) {
      next
    }
    hits <- hits[rt >= rt_lo & rt <= rt_hi]
    if (!nrow(hits)) {
      next
    }

    matched_indices <- sort(unique(hits$scan_idx))
    sp_list <- MSnbase::spectra(raw[matched_indices])
    if (!length(sp_list)) {
      next
    }

    proc_list <- lapply(sp_list, function(sp) {
      mz <- as.numeric(MSnbase::mz(sp))
      intensity <- as.numeric(MSnbase::intensity(sp))
      proc <- .nts_process_ms2_consensus(mz, intensity)
      if (!nrow(proc)) {
        return(list(mz = numeric(0), intensity = numeric(0), tic = 0))
      }
      list(
        mz = proc[, 1L],
        intensity = proc[, 2L],
        tic = sum(proc[, 2L], na.rm = TRUE)
      )
    })
    names(proc_list) <- as.character(matched_indices)

    pm_match <- pm[matched_indices]
    iso_lo_off <- meta$iso_lo[matched_indices]
    iso_hi_off <- meta$iso_hi[matched_indices]
    iso_lo_edge <- pm_match - iso_lo_off
    iso_hi_edge <- pm_match + iso_hi_off

    # Precursor isolation purity (diagnostic; only needed for scoring modes).
    # Requires the file's own MS1 scans, so skip the extra read when scoring is
    # off and the value would never be used.
    purity <- rep(NA_real_, length(matched_indices))
    if (!identical(scoring, "none")) {
      purity <- .nts_iso_purity_for_file(
        f,
        scan_rt = rt[matched_indices],
        iso_lo = iso_lo_edge,
        iso_hi = iso_hi_edge,
        prec_mz = pm_match,
        mz_tol_ppm = mz_tol_ppm,
        xdata = xdata
      )
    }

    scan_info <- data.table::data.table(
      scan_idx = matched_indices,
      tic = vapply(proc_list, function(z) z$tic, numeric(1L)),
      ce = meta$ce[matched_indices],
      polarity = meta$polarity[matched_indices],
      isolation_lo = iso_lo_edge,
      isolation_hi = iso_hi_edge,
      purity = purity,
      mz = I(lapply(proc_list, function(z) z$mz)),
      intensity = I(lapply(proc_list, function(z) z$intensity))
    )
    scan_info <- scan_info[tic > 0]
    if (!nrow(scan_info)) {
      next
    }

    cand <- merge(
      hits[, .(feature_name, scan_idx, precursorMz, rt)],
      scan_info,
      by = "scan_idx",
      all.x = FALSE,
      all.y = FALSE
    )
    if (!nrow(cand)) {
      next
    }
    cand[, source_file := basename(f)]
    ms2_map <- data.table::rbindlist(
      list(
        ms2_map,
        cand[, .(feature_name, scan_idx, source_file, tic, precursorMz, rt,
                 ce, polarity, isolation_lo, isolation_hi, purity, mz, intensity)]
      ),
      use.names = TRUE
    )
  }

  if (n_dropped_no_prec > 0L) {
    message(sprintf(
      paste0(
        "extract_ms2: skipped %d MS2 scan(s) without a finite precursor m/z ",
        "or RT (e.g. DIA / missing precursor); they cannot be matched to ",
        "features. If you expected hits, check the acquisition mode."
      ),
      n_dropped_no_prec
    ))
  }

  if (scoring != "none" && nrow(ms2_map)) {
    ms2_map <- .nts_score_dda_match(
      ms2_map,
      ft_dt,
      weights = weights,
      mz_tol_ppm = mz_tol_ppm,
      rt_tol = rt_tol
    )
    if (identical(scoring, "strict")) {
      keep <- is.na(ms2_map$match_score) | ms2_map$match_score >= min_match_score
      ms2_map <- ms2_map[keep]
    }
  }

  has_score <- isTRUE(scoring != "none") && "match_score" %in% names(ms2_map)

  out <- list()
  if (nrow(ms2_map)) {
    grp_idx <- split(seq_len(nrow(ms2_map)), ms2_map$feature_name)
    out <- lapply(names(grp_idx), function(fid) {
      idx <- grp_idx[[fid]]
      sub <- ms2_map[idx]
      lapply(seq_len(nrow(sub)), function(i) {
        el <- list(
          feature_id = fid,
          id = fid,
          scan_idx = sub$scan_idx[i],
          source_file = sub$source_file[i],
          precursorMz = sub$precursorMz[i],
          rt = sub$rt[i],
          ce = sub$ce[i],
          polarity = sub$polarity[i],
          isolation_lo = sub$isolation_lo[i],
          isolation_hi = sub$isolation_hi[i],
          mz = sub$mz[[i]],
          intensity = sub$intensity[[i]]
        )
        if (has_score) {
          el$match_score <- sub$match_score[i]
          el$s_mz <- sub$s_mz[i]
          el$s_rt <- sub$s_rt[i]
          el$s_iso <- sub$s_iso[i]
          el$s_purity <- sub$s_purity[i]
          el$s_tic <- sub$s_tic[i]
          el$s_npeaks <- sub$s_npeaks[i]
        }
        el
      })
    })
    names(out) <- names(grp_idx)
  }

  out <- .nts_ms2data_drop_empty(out)
  methods::slot(project, "ms2data") <- out
  methods::slot(project, "clean_ms2") <- list()
  project <- .nts_build_ms2_annotations_from_ms2data(project, out)
  if (length(methods::slot(project, "ms2data")) > nrow(get_ms1_annotations(project))) {
    warning(
      "extract_ms2: more MS2 feature groups in @ms2data (",
      length(methods::slot(project, "ms2data")), ") than rows in ms1_annotations (",
      nrow(get_ms1_annotations(project)),
      "). Check name alignment across ms1 / feature_table / MS2.",
      call. = FALSE
    )
  }
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "extract_ms2",
    list(
      ms2_dir = ms2_dir,
      mz_tol_ppm = mz_tol_ppm,
      rt_tol = rt_tol,
      scoring = scoring,
      weights = as.list(weights),
      min_match_score = if (identical(scoring, "strict")) min_match_score else NA_real_,
      isolation_width_da = if (is.null(isolation_width_da)) NA_real_ else isolation_width_da,
      n_files = length(files),
      n_ms2 = sum(vapply(out, length, integer(1L))),
      n_features_with_ms2 = length(out)
    )
  )
  .nts_ms2_strip_mz_bounds(project)
}

# Extract per-scan MS2 acquisition metadata from an MSnbase onDisk handle.
# Returns a named list of length-n_sp atomic vectors. All fields fall back to
# NA when the vendor / parser does not expose them.
#' @noRd
.nts_extract_scan_meta <- function(raw, n_sp, default_isolation_width_da = NULL) {
  n_sp <- as.integer(n_sp)[1L]
  if (!is.finite(n_sp) || n_sp < 1L) {
    return(list(
      ce = numeric(0),
      polarity = integer(0),
      iso_lo = numeric(0),
      iso_hi = numeric(0)
    ))
  }

  ce <- tryCatch(
    as.numeric(MSnbase::collisionEnergy(raw)),
    error = function(e) rep(NA_real_, n_sp)
  )
  if (!length(ce) || length(ce) != n_sp) {
    ce <- rep(NA_real_, n_sp)
  }

  pol <- tryCatch(
    as.integer(MSnbase::polarity(raw)),
    error = function(e) rep(NA_integer_, n_sp)
  )
  if (!length(pol) || length(pol) != n_sp) {
    pol <- rep(NA_integer_, n_sp)
  }

  iso_lo <- rep(NA_real_, n_sp)
  iso_hi <- rep(NA_real_, n_sp)
  hdr <- tryCatch(MSnbase::header(raw), error = function(e) NULL)
  if (is.data.frame(hdr) && nrow(hdr) == n_sp) {
    if ("isolationWindowLowerOffset" %in% names(hdr)) {
      iso_lo <- suppressWarnings(as.numeric(hdr$isolationWindowLowerOffset))
    }
    if ("isolationWindowUpperOffset" %in% names(hdr)) {
      iso_hi <- suppressWarnings(as.numeric(hdr$isolationWindowUpperOffset))
    }
  }

  if (!is.null(default_isolation_width_da)) {
    w <- suppressWarnings(as.numeric(default_isolation_width_da)[1L])
    if (is.finite(w) && w > 0) {
      half <- w / 2
      iso_lo[!is.finite(iso_lo)] <- half
      iso_hi[!is.finite(iso_hi)] <- half
    }
  }

  list(ce = ce, polarity = pol, iso_lo = iso_lo, iso_hi = iso_hi)
}

# Compute DDA scan-to-feature multi-evidence match scores.
# Inputs:
#   ms2_map: data.table with columns feature_name, source_file, precursorMz,
#            rt, tic, mz (list), intensity (list).
#   ft_dt  : data.table with feature_name, mz (feature precursor), rt.
# Component sub-scores (each in [0, 1], NA when not computable):
#   s_mz, s_rt, s_iso, s_purity, s_tic, s_npeaks.
# Composite match_score is a weight-normalised mean over non-NA components.
# s_iso is derived from the scan's isolation window; s_purity is carried in from
# extract_ms2 (per-scan MS1 isolation purity) and is NA when unavailable.
#' @noRd
.nts_score_dda_match <- function(ms2_map,
                                 ft_dt,
                                 weights = c(mz = 0.30, rt = 0.20, iso = 0.20,
                                             purity = 0.15, tic = 0.05, npeaks = 0.10),
                                 mz_tol_ppm = 15,
                                 rt_tol = 30) {
  dt <- data.table::as.data.table(ms2_map)
  if (!nrow(dt)) {
    return(dt)
  }
  if (!all(c("feature_name", "precursorMz", "rt", "tic", "mz", "source_file") %in% names(dt))) {
    return(dt)
  }

  ft <- data.table::as.data.table(ft_dt)
  if (!all(c("feature_name", "mz", "rt") %in% names(ft))) {
    return(dt)
  }
  ft_meta <- unique(
    ft[, .(feature_name = as.character(feature_name),
           feat_mz = as.numeric(mz),
           feat_rt = as.numeric(rt))],
    by = "feature_name"
  )
  dt <- ft_meta[dt, on = "feature_name"]

  mz_tol <- suppressWarnings(as.numeric(mz_tol_ppm)[1L])
  if (!is.finite(mz_tol) || mz_tol <= 0) {
    mz_tol <- 15
  }
  rt_tol_v <- suppressWarnings(as.numeric(rt_tol)[1L])
  if (!is.finite(rt_tol_v) || rt_tol_v <= 0) {
    rt_tol_v <- 30
  }

  dt[, ppm_err := ifelse(
    is.finite(feat_mz) & feat_mz > 0,
    abs(precursorMz - feat_mz) / feat_mz * 1e6,
    NA_real_
  )]
  dt[, rt_err := abs(rt - feat_rt)]
  dt[, s_mz := pmax(0, pmin(1, 1 - ppm_err / mz_tol))]
  dt[, s_rt := pmax(0, pmin(1, 1 - rt_err / rt_tol_v))]

  dt[, n_peaks_proc := vapply(mz, function(x) sum(is.finite(x)), integer(1L))]
  dt[, s_npeaks := 1 - exp(-as.numeric(n_peaks_proc) / 5)]

  dt[, file_tic_med := stats::median(tic, na.rm = TRUE), by = source_file]
  dt[, s_tic := pmin(1, ifelse(
    is.finite(file_tic_med) & file_tic_med > 0,
    tic / file_tic_med,
    NA_real_
  ))]

  # s_iso: how centrally the feature m/z sits inside the scan's isolation
  # window (1 = dead centre, 0 = at the edge). A feature hugging the window edge
  # signals likely co-isolation.
  if (all(c("isolation_lo", "isolation_hi") %in% names(dt))) {
    half <- (as.numeric(dt$isolation_hi) - as.numeric(dt$isolation_lo)) / 2
    dt[, s_iso := ifelse(
      is.finite(half) & half > 0 & is.finite(feat_mz),
      pmax(0, pmin(1, 1 - abs(precursorMz - feat_mz) / half)),
      NA_real_
    )]
  } else {
    dt[, s_iso := NA_real_]
  }

  # s_purity: precursor isolation purity carried in from extract_ms2 (fraction
  # of in-window MS1 intensity attributable to the target precursor). NA when the
  # file exposes no MS1 scans or purity was not computed.
  if ("purity" %in% names(dt)) {
    dt[, s_purity := pmax(0, pmin(1, as.numeric(purity)))]
  } else {
    dt[, s_purity := NA_real_]
  }

  comp_cols <- c("s_mz", "s_rt", "s_iso", "s_purity", "s_tic", "s_npeaks")
  comp_names <- c("mz", "rt", "iso", "purity", "tic", "npeaks")
  w_vec <- suppressWarnings(as.numeric(weights[comp_names]))
  w_vec[!is.finite(w_vec) | w_vec < 0] <- 0
  if (sum(w_vec) <= 0) {
    w_vec <- c(0.30, 0.20, 0.20, 0.15, 0.05, 0.10)
  }

  mat <- as.matrix(dt[, ..comp_cols])
  num <- rowSums(t(t(mat) * w_vec), na.rm = TRUE)
  den <- as.numeric((!is.na(mat)) %*% w_vec)
  dt[, match_score := ifelse(den > 0, num / den, NA_real_)]

  dt[, c("feat_mz", "feat_rt", "ppm_err", "rt_err",
         "n_peaks_proc", "file_tic_med") := NULL]
  dt
}

#' @noRd
.nts_process_ms2_consensus <- function(mz, intensity) {
  ok <- is.finite(mz) & is.finite(intensity) & intensity > 0
  mz <- mz[ok]
  intensity <- intensity[ok]
  if (!length(mz)) {
    return(matrix(numeric(0), ncol = 2L))
  }
  mx <- max(intensity, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) {
    return(matrix(numeric(0), ncol = 2L))
  }
  rel <- intensity / mx
  keep <- rel >= 0.0025
  mz <- mz[keep]
  intensity <- intensity[keep]
  if (!length(mz)) {
    return(matrix(numeric(0), ncol = 2L))
  }
  intensity <- (intensity / max(intensity, na.rm = TRUE)) * 1000
  cbind(mz, intensity)
}

#' Purify MS2 spectra and store per-feature spectrum lists on the project
#'
#' Runs per-scan relative-intensity filtering (default 5% floor), m/z binning,
#' cross-scan consensus (absolute `min_occurrence` and/or fractional
#' `min_freq`), and precursor-region removal (`neutral_loss < 1` Da),
#' then fills `@clean_ms2` with a **named list** of **only** features that have
#' at least one raw scan and a non-empty consensus `spectrum`.
#'
#' When `by_ce = TRUE` and at least one scan element carries a finite `ce`
#' (collision energy, written by [extract_ms2]), spectra are stratified by
#' CE bin first; the **total** `spectrum` is then assembled per
#' `merge_strategy`:
#'
#' * `"ce_top_intensity"` (default): each CE bin is purified independently,
#'   then for every aligned `mz` bin the **most intense** CE sub-consensus
#'   wins. Fragment visibility increases when one CE is far better suited.
#' * `"union_then_filter"`: pool all scans first (legacy behaviour, exactly
#'   equivalent to versions < 0.1.3 of this function).
#'
#' Each `@clean_ms2[[k]]` element now exposes:
#' `feature_id`, `precursorMZ`, `rt`, `n_raw_scans`, `spectrum`
#' (`data.frame` with `mz`, `intensity`, `frequency`), and -- only when CE
#' information is available -- `ce_levels` (numeric vector),
#' `ce_subspectra` (named list of per-CE `mz/intensity/frequency` frames),
#' `ce_summary` (`data.frame` with `ce`, `n_scans`, `n_peaks`).
#'
#' @param project An [NTSProject] with non-empty `@ms2data` (unless `ms2_list`
#'   is passed) and `ms1_annotations` with `name` and `mz` (precursor m/z).
#' @param mz_tol_da Bin width (Da) for alignment: `round(mz / mz_tol_da) * mz_tol_da`.
#' @param min_occurrence Minimum scan-count for a binned ion to be kept (OR gate).
#' @param min_freq Minimum fraction of scans (`occurrence >= n_scans * min_freq`, OR gate).
#'   Default `0` keeps single-scan ions when `min_occurrence` is not met.
#' @param rel_int_floor Minimum relative intensity within each scan (default `0.05`).
#' @param ms2_list Optional named list of MS2 scan lists; `NULL` uses `@ms2data`.
#' @param by_ce Logical. When `TRUE` (default) and CE metadata is available,
#'   build per-CE sub-consensus spectra and use them when assembling the total
#'   `spectrum`. Falls back to the legacy single-group pipeline when no scan
#'   carries finite CE.
#' @param ce_bin Numeric. CE bucket size used to group scans, in CE units
#'   (default `5`). NA CE values are gathered into a single `"CE_NA"` group.
#' @param merge_strategy One of `"ce_top_intensity"` (default) or
#'   `"union_then_filter"`; see Details.
#' @param keep_ce_subspectra Logical. Store per-CE sub-consensus spectra
#'   on each `@clean_ms2[[k]]` element (default `TRUE`). Set to `FALSE` to
#'   reduce memory footprint while keeping the total `spectrum`.
#' @param precursor_filter Either `NULL` (legacy behaviour: drop
#'   `neutral_loss < 1 Da` only) or a named list of rules consumed by the
#'   internal cleaner. Recognised keys: `drop_neutral_loss_below` (Da),
#'   `drop_isotope_envelope` (`TRUE`/`FALSE`), `drop_window_da` (length-2
#'   numeric, both sides of `prec_mz`), `drop_low_mass_below` (Da), and
#'   `mz_tol_da` (isotope-envelope tolerance, default `0.01`).
#'   Use [precursor_filter_template()] for ready-made presets.
#'
#' @return Updated `project` with `@clean_ms2` populated (named `list`) and
#'   `@ms2_annotations` decorated with `n_fragments`, `n_ce_levels`,
#'   `ce_min`, `ce_max` for surviving features.
#' @export
purify_ms2_spectra <- function(project,
                               mz_tol_da = 0.01,
                               min_occurrence = 2L,
                               min_freq = 0.0,
                               rel_int_floor = 0.05,
                               ms2_list = NULL,
                               by_ce = TRUE,
                               ce_bin = 5,
                               merge_strategy = c("ce_top_intensity",
                                                  "union_then_filter"),
                               keep_ce_subspectra = TRUE,
                               precursor_filter = NULL) {
  merge_strategy <- match.arg(merge_strategy)
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  if (is.null(ms2_list)) {
    ms2_list <- methods::slot(project, "ms2data")
  }
  ms1 <- get_ms1_annotations(project)
  if (!is.data.frame(ms1) || !nrow(ms1) || !all(c("name", "mz") %in% names(ms1))) {
    stop(
      "ms1_annotations must be non-empty with columns name and mz (precursor m/z).",
      call. = FALSE
    )
  }

  if (!length(ms2_list)) {
    warning("ms2data / ms2_list is empty; @clean_ms2 cleared.", call. = FALSE)
    methods::slot(project, "clean_ms2") <- list()
    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "purify_ms2_spectra",
      list(
        note = "empty_ms2",
        mz_tol_da = mz_tol_da,
        min_occurrence = min_occurrence,
        min_freq = min_freq
      )
    )
    return(project)
  }

  m2_seed <- get_ms2_annotations(project)
  if (!is.data.frame(m2_seed) || !nrow(m2_seed) || !"name" %in% names(m2_seed)) {
    project <- .nts_build_ms2_annotations_from_ms2data(project, ms2_list)
    ms2_list <- methods::slot(project, "ms2data")
  }

  ms2_names <- names(ms2_list)
  if (is.null(ms2_names)) {
    ms2_names <- as.character(seq_along(ms2_list))
    names(ms2_list) <- ms2_names
  }
  ms1_names <- unique(as.character(ms1$name))
  keep_idx <- ms2_names %in% ms1_names
  ms2_list <- ms2_list[keep_idx]
  feat_names <- unique(names(ms2_list))
  if (!length(feat_names)) {
    warning("ms2_list has no names matching ms1_annotations; @clean_ms2 cleared.", call. = FALSE)
    methods::slot(project, "clean_ms2") <- list()
    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "purify_ms2_spectra",
      list(
        note = "ms2_names_not_in_ms1",
        mz_tol_da = mz_tol_da,
        min_occurrence = min_occurrence,
        min_freq = min_freq
      )
    )
    return(project)
  }

  ms2_raw <- .nts_ms2_list_to_long_scans(ms2_list, int_threshold = 0)
  raw_scan_n <- if (nrow(ms2_raw)) {
    ms2_raw[, .(n_raw_scans = data.table::uniqueN(scan_id)), by = name]
  } else {
    data.table::data.table(name = character(0), n_raw_scans = integer(0))
  }

  use_ce <- isTRUE(by_ce) &&
    identical(merge_strategy, "ce_top_intensity") &&
    "ce" %in% names(ms2_raw) &&
    nrow(ms2_raw) > 0L &&
    any(is.finite(ms2_raw$ce))

  prec_dt <- unique(
    data.table::as.data.table(ms1)[, .(name = as.character(name), prec_mz = as.numeric(mz))],
    by = "name"
  )

  ce_sub_dt <- NULL
  ce_nscan <- NULL
  if (use_ce) {
    ce_step <- as.numeric(ce_bin)[1L]
    if (!is.finite(ce_step) || ce_step <= 0) {
      ce_step <- 1
    }
    ms2_raw[, ce_grp := ifelse(
      is.finite(ce),
      sprintf("CE%g", round(ce / ce_step) * ce_step),
      "CE_NA"
    )]
    ce_nscan <- ms2_raw[, .(n_scans_ce = data.table::uniqueN(scan_id)),
                       by = .(name, ce_grp)]
    raw_sub <- ms2_raw[, .(
      name = paste(name, ce_grp, sep = "|||"),
      scan_id = scan_id,
      mz = mz,
      intensity = intensity
    )]
    clean_core_sub <- .nts_purify_raw_ms2(
      raw_sub,
      mz_tol_da = mz_tol_da,
      min_occurrence = min_occurrence,
      min_freq = min_freq,
      rel_int_min = rel_int_floor
    )
    if (nrow(clean_core_sub)) {
      sp <- data.table::tstrsplit(clean_core_sub$name, "|||", fixed = TRUE)
      clean_core_sub[, name := sp[[1L]]]
      clean_core_sub[, ce_grp := sp[[2L]]]
      ce_sub_dt <- .nts_ms2_clean_with_precursor(
        clean_core_sub,
        prec_dt,
        rules = precursor_filter
      )
    } else {
      ce_sub_dt <- clean_core_sub[0L]
      ce_sub_dt[, ce_grp := character(0)]
    }
    if (nrow(ce_sub_dt)) {
      ce_sub_dt <- ce_nscan[ce_sub_dt, on = c("name", "ce_grp")]
      ce_sub_dt[is.na(n_scans_ce), n_scans_ce := 0L]
      ce_sub_dt[, frequency := ifelse(
        n_scans_ce > 0L,
        as.numeric(occurrence) / as.numeric(n_scans_ce),
        NA_real_
      )]
      ce_sub_dt[, mz_bin_key := round(as.numeric(mz) / mz_tol_da) * mz_tol_da]
      data.table::setorderv(
        ce_sub_dt,
        c("name", "mz_bin_key", "intensity"),
        c(1L, 1L, -1L)
      )
      clean_dt <- unique(ce_sub_dt, by = c("name", "mz_bin_key"))[
        ,
        .(name, mz, intensity, frequency)
      ]
    } else {
      clean_dt <- data.table::data.table(
        name = character(0),
        mz = numeric(0),
        intensity = numeric(0),
        frequency = numeric(0)
      )
    }
  } else {
    clean_core <- .nts_purify_raw_ms2(
      ms2_raw[, .(name, scan_id, mz, intensity)],
      mz_tol_da = mz_tol_da,
      min_occurrence = min_occurrence,
      min_freq = min_freq,
      rel_int_min = rel_int_floor
    )
    clean_dt <- .nts_ms2_clean_with_precursor(
      clean_core,
      prec_dt,
      rules = precursor_filter
    )
    clean_dt <- raw_scan_n[clean_dt, on = "name"]
    clean_dt[is.na(n_raw_scans), n_raw_scans := 0L]
    clean_dt[, frequency := ifelse(
      n_raw_scans > 0L,
      as.numeric(occurrence) / as.numeric(n_raw_scans),
      NA_real_
    )]
  }

  ms1_dt <- data.table::as.data.table(ms1)
  rt_lookup <- unique(
    ms1_dt[, .(name = as.character(name), precursorMZ = as.numeric(mz), rt = as.numeric(rt))],
    by = "name"
  )

  ns_tbl <- raw_scan_n[data.table::data.table(name = feat_names), on = "name"]
  ns_tbl[is.na(n_raw_scans), n_raw_scans := 0L]

  ce_split <- if (use_ce && !is.null(ce_sub_dt) && nrow(ce_sub_dt)) {
    split(ce_sub_dt, by = "name", keep.by = FALSE)
  } else {
    list()
  }

  out_list <- vector("list", length(feat_names))
  names(out_list) <- feat_names
  n_peaks <- 0L
  clean_by_name <- if (nrow(clean_dt)) split(clean_dt, clean_dt$name) else list()
  for (k in seq_along(feat_names)) {
    fn <- feat_names[k]
    sub <- clean_by_name[[fn]]
    if (is.null(sub)) {
      sub <- clean_dt[0L]
    }
    row <- rt_lookup[name == fn]
    prec <- if (nrow(row)) row$precursorMZ[1L] else NA_real_
    rt_val <- if (nrow(row) && "rt" %in% names(row)) row$rt[1L] else NA_real_
    n_raw <- as.integer(ns_tbl[k, n_raw_scans])
    if (!length(n_raw) || is.na(n_raw)) {
      n_raw <- 0L
    }
    if (nrow(sub)) {
      sp <- data.frame(
        mz = sub$mz,
        intensity = sub$intensity,
        frequency = sub$frequency,
        stringsAsFactors = FALSE
      )
      n_peaks <- n_peaks + nrow(sp)
    } else {
      sp <- data.frame(
        mz = numeric(0),
        intensity = numeric(0),
        frequency = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    el <- list(
      feature_id = fn,
      precursorMZ = prec,
      rt = rt_val,
      n_raw_scans = n_raw,
      spectrum = sp
    )

    if (use_ce) {
      ce_info <- .nts_assemble_ce_outputs(
        feature_name = fn,
        ce_feat_dt = ce_split[[fn]],
        ce_nscan = ce_nscan,
        keep_ce_subspectra = keep_ce_subspectra
      )
      el$ce_levels <- ce_info$ce_levels
      el$ce_subspectra <- ce_info$ce_subspectra
      el$ce_summary <- ce_info$ce_summary
    }
    out_list[[k]] <- el
  }

  keep_el <- vapply(out_list, function(z) {
    nr <- suppressWarnings(as.integer(z$n_raw_scans)[1L])
    spf <- z$spectrum
    (is.finite(nr) && nr > 0L) && is.data.frame(spf) && nrow(spf) > 0L
  }, logical(1L))
  out_list <- out_list[keep_el]

  methods::slot(project, "clean_ms2") <- out_list
  project <- .nts_update_ms2_annotations_n_fragments(project)
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "purify_ms2_spectra",
    list(
      mz_tol_da = mz_tol_da,
      min_occurrence = min_occurrence,
      min_freq = min_freq,
      rel_int_floor = rel_int_floor,
      by_ce = isTRUE(by_ce),
      ce_bin = ce_bin,
      merge_strategy = merge_strategy,
      keep_ce_subspectra = isTRUE(keep_ce_subspectra),
      ce_used = use_ce,
      precursor_filter = if (is.null(precursor_filter)) "legacy_nl_ge_1" else precursor_filter,
      n_peaks = n_peaks,
      n_features = length(out_list)
    )
  )
  project
}

# Build per-feature CE outputs: ce_levels (numeric, sorted ascending),
# ce_subspectra (named list keyed by "CExx"/"CE_NA"), ce_summary
# (`data.frame(ce, n_scans, n_peaks)`). Empty atoms when no peaks survive.
#' @noRd
.nts_assemble_ce_outputs <- function(feature_name,
                                     ce_feat_dt,
                                     ce_nscan,
                                     keep_ce_subspectra = TRUE) {
  empty <- list(
    ce_levels = numeric(0),
    ce_subspectra = list(),
    ce_summary = data.frame(
      ce = numeric(0),
      n_scans = integer(0),
      n_peaks = integer(0),
      stringsAsFactors = FALSE
    )
  )
  if (is.null(ce_feat_dt) || !nrow(ce_feat_dt)) {
    return(empty)
  }

  data.table::setorderv(ce_feat_dt, c("ce_grp", "mz"))
  grp_split <- split(ce_feat_dt, by = "ce_grp", keep.by = FALSE)
  ce_strs <- names(grp_split)

  ce_nums <- suppressWarnings(as.numeric(sub("^CE", "", ce_strs)))
  ce_nums[ce_strs == "CE_NA"] <- NA_real_

  ord <- order(ce_nums, na.last = TRUE)
  ce_strs <- ce_strs[ord]
  ce_nums <- ce_nums[ord]
  grp_split <- grp_split[ord]

  sub_list <- lapply(grp_split, function(d) {
    data.frame(
      mz = as.numeric(d$mz),
      intensity = as.numeric(d$intensity),
      frequency = as.numeric(d$frequency),
      stringsAsFactors = FALSE
    )
  })

  nscan_vec <- integer(length(ce_strs))
  if (!is.null(ce_nscan) && nrow(ce_nscan)) {
    sub_nscan <- ce_nscan[name == feature_name]
    idx <- match(ce_strs, sub_nscan$ce_grp)
    nscan_vec <- as.integer(sub_nscan$n_scans_ce[idx])
    nscan_vec[is.na(nscan_vec)] <- 0L
  }

  summary_df <- data.frame(
    ce = ce_nums,
    n_scans = nscan_vec,
    n_peaks = vapply(sub_list, nrow, integer(1L)),
    stringsAsFactors = FALSE
  )

  if (!isTRUE(keep_ce_subspectra)) {
    sub_list <- list()
  }

  list(
    ce_levels = ce_nums,
    ce_subspectra = sub_list,
    ce_summary = summary_df
  )
}

# ---------------------------------------------------------------------------
# Internal: MS2 list flattening and consensus (used by purify_ms2_spectra &
# network MS2 edge mining in 07_network_builder.R)
# ---------------------------------------------------------------------------
# MS2 list -> long table with scan_id (internal)
# When scan elements carry `ce` (added by extract_ms2), it is captured into a
# `ce` column (NA when absent). Downstream callers that don't need CE simply
# ignore the extra column.
#' @noRd
.nts_ms2_list_to_long_scans <- function(ms2_list, int_threshold = 0) {
  empty_dt <- data.table::data.table(
    name = character(0),
    scan_id = integer(0),
    mz = numeric(0),
    intensity = numeric(0),
    ce = numeric(0)
  )
  if (!length(ms2_list)) {
    return(empty_dt)
  }
  nms <- names(ms2_list)
  if (is.null(nms)) {
    nms <- as.character(seq_along(ms2_list))
  }
  thr <- as.numeric(int_threshold)[1L]
  if (!is.finite(thr)) {
    thr <- 0
  }
  parts <- vector("list", 0L)
  for (i in seq_along(ms2_list)) {
    nm <- nms[i]
    scans <- ms2_list[[i]]
    if (!length(scans)) {
      next
    }
    for (j in seq_along(scans)) {
      sp <- scans[[j]]
      mz <- as.numeric(sp$mz)
      int <- as.numeric(sp$intensity)
      ok <- is.finite(mz) & is.finite(int) & int > thr
      mz <- mz[ok]
      int <- int[ok]
      if (!length(mz)) {
        next
      }
      ce_val <- suppressWarnings(as.numeric(sp$ce)[1L])
      if (!length(ce_val)) {
        ce_val <- NA_real_
      }
      parts[[length(parts) + 1L]] <- data.table::data.table(
        name = nm,
        scan_id = j,
        mz = mz,
        intensity = int,
        ce = ce_val
      )
    }
  }
  if (!length(parts)) {
    return(empty_dt)
  }
  data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
}

# Per-scan rel-int floor, binning, cross-scan consensus (internal)
#
# Retention (per name, mz_bin): keep if occurrence >= min_occurrence OR
# occurrence >= n_scans * min_freq. Supports multi-CE / pooled data (absolute
# replicate count) while min_freq = 0 preserves single-scan ions.
#' @noRd
.nts_purify_raw_ms2 <- function(raw_ms2_long,
                                mz_tol_da = 0.01,
                                min_occurrence = 2L,
                                min_freq = 0.0,
                                rel_int_min = 0.05) {
  dt <- data.table::as.data.table(raw_ms2_long)
  empty_purify <- function() {
    data.table::data.table(
      name = character(0),
      mz = numeric(0),
      intensity = numeric(0),
      occurrence = integer(0),
      n_scans = integer(0)
    )
  }
  if (!nrow(dt)) {
    return(empty_purify())
  }
  req <- c("name", "mz", "intensity")
  if (!all(req %in% names(dt))) {
    stop(".nts_purify_raw_ms2 needs columns: name, mz, intensity", call. = FALSE)
  }
  if (!"scan_id" %in% names(dt)) {
    dt[, scan_id := 1L]
  } else {
    dt[, scan_id := as.integer(scan_id)]
  }
  dt <- dt[is.finite(mz) & is.finite(intensity) & intensity > 0]
  if (!nrow(dt)) {
    return(empty_purify())
  }
  dt[, rel_int := intensity / max(intensity, na.rm = TRUE), by = .(name, scan_id)]
  dt <- dt[is.finite(rel_int) & rel_int >= rel_int_min]
  if (!nrow(dt)) {
    return(empty_purify())
  }
  tol <- as.numeric(mz_tol_da)[1L]
  .nts_check(tol, "positive_finite", .arg = "mz_tol_da")
  dt[, mz_bin := round(mz / tol) * tol]
  nscan <- dt[, .(n_scans = data.table::uniqueN(scan_id)), by = name]
  agg <- dt[, .(
    occurrence = .N,
    intensity = mean(intensity),
    mz = mean(mz)
  ), by = .(name, mz_bin)]
  agg <- nscan[agg, on = "name"]
  mo <- suppressWarnings(as.integer(min_occurrence)[1L])
  if (!is.finite(mo)) {
    mo <- 2L
  }
  mo <- max(0L, mo)
  mf <- as.numeric(min_freq)[1L]
  if (!is.finite(mf)) {
    mf <- 0
  }
  mf <- max(0, min(1, mf))
  thr <- as.numeric(agg$n_scans) * mf
  agg <- agg[
    is.finite(occurrence) & is.finite(n_scans) &
      (occurrence >= mo | occurrence >= thr)
  ]
  agg[, .(name, mz, intensity, occurrence, n_scans)]
}

# Precursor join + configurable precursor-region cleaning (internal)
#
# When `rules = NULL` the legacy behaviour is preserved exactly:
# drop rows with `neutral_loss < 1.0 Da`. When `rules` is a named list, the
# following keys are honoured (any subset; NULL/FALSE entries are ignored):
#   * drop_neutral_loss_below : numeric Da; drop rows where prec_mz - mz < x
#   * drop_isotope_envelope   : TRUE/FALSE; drop peaks within `mz_tol_da` of
#                               prec_mz - n * 1.003355 for n = 1, 2
#   * drop_window_da          : length-2 numeric c(lo, hi) absolute m/z window
#                               (peaks with mz in [prec_mz + lo, prec_mz + hi]
#                               are removed). Use e.g. c(-0.5, 5) to remove
#                               isolation neighbours.
#   * drop_low_mass_below     : numeric Da; drop peaks with mz < x
#   * mz_tol_da               : tolerance used by isotope-envelope rule
#                               (default 0.01)
#' @noRd
.nts_ms2_clean_with_precursor <- function(clean_long, prec_dt, rules = NULL) {
  if (!nrow(clean_long) || !nrow(prec_dt)) {
    return(clean_long[0L])
  }
  if (!all(c("name", "mz") %in% names(clean_long))) {
    return(clean_long[0L])
  }
  if (!all(c("name", "prec_mz") %in% names(prec_dt))) {
    stop("prec_dt needs name, prec_mz", call. = FALSE)
  }
  x <- prec_dt[clean_long, on = "name", nomatch = 0L]
  x <- x[is.finite(mz) & is.finite(prec_mz)]
  x[, neutral_loss := prec_mz - mz]
  if (is.null(rules)) {
    return(x[neutral_loss >= 1.0])
  }
  .nts_apply_precursor_rules(x, rules)
}

# Apply a precursor-region rule list to a data.table that already carries
# `mz`, `prec_mz` (and `neutral_loss = prec_mz - mz`). Returns the filtered
# table.
#' @noRd
.nts_apply_precursor_rules <- function(x, rules) {
  if (!nrow(x)) {
    return(x)
  }
  if (!is.list(rules)) {
    stop("`precursor_filter` must be a named list of rules.", call. = FALSE)
  }
  keep <- rep(TRUE, nrow(x))
  nl_below <- suppressWarnings(as.numeric(rules$drop_neutral_loss_below)[1L])
  if (length(nl_below) && is.finite(nl_below)) {
    keep <- keep & (x$neutral_loss >= nl_below)
  }
  low_mass <- suppressWarnings(as.numeric(rules$drop_low_mass_below)[1L])
  if (length(low_mass) && is.finite(low_mass)) {
    keep <- keep & (x$mz >= low_mass)
  }
  mz_tol_da <- suppressWarnings(as.numeric(rules$mz_tol_da)[1L])
  if (!length(mz_tol_da) || !is.finite(mz_tol_da) || mz_tol_da <= 0) {
    mz_tol_da <- 0.01
  }
  if (isTRUE(rules$drop_isotope_envelope)) {
    iso_step <- 1.003355
    for (n in 1:2) {
      target <- x$prec_mz - n * iso_step
      keep <- keep & !(abs(x$mz - target) <= mz_tol_da)
    }
  }
  win <- suppressWarnings(as.numeric(rules$drop_window_da))
  if (length(win) == 2L && all(is.finite(win))) {
    lo <- min(win)
    hi <- max(win)
    keep <- keep & !(x$mz >= (x$prec_mz + lo) & x$mz <= (x$prec_mz + hi))
  }
  x[keep]
}

#' Pre-built precursor-region cleaning rule templates
#'
#' Helper for [purify_ms2_spectra]'s `precursor_filter` argument. Returns a
#' named list compatible with the internal rule engine; pass it directly to
#' [purify_ms2_spectra]. Custom users can build their own list with the same
#' key names instead.
#'
#' @param type One of `"default"` (legacy, drop `neutral_loss < 1 Da`),
#'   `"qtof"` (default + isotope envelope), `"orbitrap"` (qtof + 0.5-5 Da
#'   precursor window), `"permissive"` (drop only `neutral_loss < 0.2 Da`),
#'   or `"strict"` (low-mass cut + isotope envelope + +/-5 Da window).
#'
#' @return Named list. Recognised keys: `drop_neutral_loss_below`,
#'   `drop_isotope_envelope`, `drop_window_da`, `drop_low_mass_below`,
#'   `mz_tol_da`.
#'
#' @examples
#' precursor_filter_template("qtof")
#' precursor_filter_template("strict")
#' @export
precursor_filter_template <- function(type = c("default", "qtof", "orbitrap",
                                               "permissive", "strict")) {
  type <- match.arg(type)
  switch(
    type,
    default = list(
      drop_neutral_loss_below = 1.0,
      mz_tol_da = 0.01
    ),
    qtof = list(
      drop_neutral_loss_below = 1.0,
      drop_isotope_envelope = TRUE,
      mz_tol_da = 0.01
    ),
    orbitrap = list(
      drop_neutral_loss_below = 0.5,
      drop_isotope_envelope = TRUE,
      drop_window_da = c(0.5, 5),
      mz_tol_da = 0.005
    ),
    permissive = list(
      drop_neutral_loss_below = 0.2,
      mz_tol_da = 0.01
    ),
    strict = list(
      drop_neutral_loss_below = 1.0,
      drop_isotope_envelope = TRUE,
      drop_window_da = c(-5, 5),
      drop_low_mass_below = 70,
      mz_tol_da = 0.005
    )
  )
}

# Flatten @clean_ms2 in one pass for MS2 network joins: peaks + precursor table.
# Optional `ms2_anchor_names` restricts to features listed in `@ms2_annotations`.
#' @noRd
.nts_clean_ms2_network_tables <- function(clean_ms2, ms2_anchor_names = NULL) {
  empty_clean <- data.table::data.table(
    name = character(0),
    mz = numeric(0),
    intensity = numeric(0),
    frequency = numeric(0)
  )
  empty_prec <- data.table::data.table(name = character(0), prec_mz = numeric(0))
  if (!is.list(clean_ms2) || !length(clean_ms2)) {
    return(list(clean_dt = empty_clean, prec_dt = empty_prec))
  }
  kn <- if (length(ms2_anchor_names)) {
    unique(as.character(ms2_anchor_names))
  } else {
    character(0)
  }
  kn <- kn[!is.na(kn) & nzchar(kn)]
  nms <- names(clean_ms2)
  n_el <- length(clean_ms2)
  peak_parts <- vector("list", n_el)
  prec_parts <- vector("list", n_el)
  for (i in seq_len(n_el)) {
    el <- clean_ms2[[i]]
    if (!is.list(el)) {
      next
    }
    nm <- if (!is.null(nms) && length(nms) >= i && nzchar(nms[i])) {
      nms[i]
    } else {
      ""
    }
    fid <- el$feature_id
    if (length(fid) && nzchar(as.character(fid)[1L])) {
      nm <- as.character(fid)[1L]
    }
    if (!nzchar(nm)) {
      next
    }
    if (length(kn) && !(nm %in% kn)) {
      next
    }
    pm <- suppressWarnings(as.numeric(el$precursorMZ)[1L])
    prec_parts[[i]] <- data.table::data.table(name = nm, prec_mz = pm)
    sp <- el$spectrum
    if (!is.data.frame(sp) || !nrow(sp)) {
      next
    }
    if (!all(c("mz", "intensity") %in% names(sp))) {
      next
    }
    mz_c <- as.numeric(sp$mz)
    int_c <- as.numeric(sp$intensity)
    if ("frequency" %in% names(sp)) {
      fr_c <- suppressWarnings(as.numeric(sp$frequency))
    } else {
      fr_c <- rep(NA_real_, length(mz_c))
    }
    ok <- is.finite(mz_c) & is.finite(int_c)
    if (!any(ok)) {
      next
    }
    if (length(fr_c) != length(mz_c)) {
      fr_c <- rep(NA_real_, length(mz_c))
    }
    peak_parts[[i]] <- data.table::data.table(
      name = nm,
      mz = mz_c[ok],
      intensity = int_c[ok],
      frequency = fr_c[ok]
    )
  }
  clean_dt <- data.table::rbindlist(peak_parts[!vapply(peak_parts, is.null, logical(1L))])
  prec_dt <- data.table::rbindlist(prec_parts[!vapply(prec_parts, is.null, logical(1L))])
  if (nrow(prec_dt)) {
    prec_dt <- prec_dt[is.finite(prec_mz)]
    prec_dt <- unique(prec_dt, by = "name")
  } else {
    prec_dt <- empty_prec
  }
  if (!nrow(clean_dt)) {
    clean_dt <- empty_clean
  }
  list(clean_dt = clean_dt, prec_dt = prec_dt)
}

# Flatten @clean_ms2 spectrum list to long peaks for MS2 network joins

# Drop m/z bins present in > max_bg_freq fraction of features (background ions)
#' @noRd
.nts_clean_ms2_filter_bg_mz <- function(clean_dt, max_bg_freq) {
  dt <- data.table::as.data.table(clean_dt)
  if (!nrow(dt)) {
    return(dt)
  }
  mbf <- max_bg_freq[1L]
  if (!is.numeric(mbf) || !is.finite(mbf) || mbf <= 0 || mbf >= 1) {
    return(dt)
  }
  n_feat <- dt[, data.table::uniqueN(name)]
  if (n_feat < 3L) {
    return(dt)
  }
  dt[, mz_k := round(as.numeric(mz), 2L)]
  ub <- dt[, .(nf = data.table::uniqueN(name)), by = mz_k]
  bad <- ub[nf / n_feat > mbf, mz_k]
  if (length(bad)) {
    dt <- dt[!mz_k %in% bad]
  }
  dt[, mz_k := NULL]
  dt
}

