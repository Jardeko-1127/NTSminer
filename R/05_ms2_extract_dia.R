# Waters MSE / DIA MS2 extraction (co-elution assignment) ----------------------------------
#
# DIA entry point that is mutually exclusive with [extract_ms2()]; currently only
# Waters MSE is implemented. A Waters MSE file usually has **the whole file at
# msLevel == 1**, with low/high energy written alternately and distinguished only
# by the function= field in `spectrumId` (1 = low-energy MS1, 2 = high-energy
# MSE/MS2). This function splits on that field and co-elution-assigns each feature
# by the **Pearson correlation between the MS1 EIC and candidate fragment EICs**,
# then writes back to `@ms2data`; downstream `purify_ms2_spectra()` /
# `compute_ms2_qc()` need no change.


# Parse the function=N field in spectrumId ------------------------------------
#' @noRd
.nts_dia_parse_function <- function(spectrum_id) {
  s <- as.character(spectrum_id)
  out <- suppressWarnings(as.integer(sub(".*function=([0-9]+).*", "\\1", s)))
  out[!is.finite(out)] <- NA_integer_
  out
}


# Read the mzR header and split scans into MS1 / high-energy ------------------------------
# Returns list(ms1 = data.table, ms2 = data.table, ce_lo = scalar, ce_hi = scalar)
#' @noRd
.nts_dia_split_scans <- function(file, mode = "waters_mse") {
  if (!requireNamespace("mzR", quietly = TRUE)) {
    stop("mzR is required for extract_ms2_dia(); install Bioconductor mzR.",
         call. = FALSE)
  }
  m <- mzR::openMSfile(file)
  on.exit(try(mzR::close(m), silent = TRUE), add = TRUE)
  hdr <- mzR::header(m)
  if (!is.data.frame(hdr) || !nrow(hdr)) {
    stop("Empty header for ", file, call. = FALSE)
  }
  fn <- .nts_dia_parse_function(hdr$spectrumId)
  if (all(is.na(fn))) {
    stop(
      "Cannot detect Waters MSE function= field in spectrumId; ",
      "is this really a Waters MSE mzML? File: ", file,
      call. = FALSE
    )
  }
  is_ms1 <- fn == 1L & !is.na(fn)
  is_ms2 <- fn == 2L & !is.na(fn)
  if (!any(is_ms1) || !any(is_ms2)) {
    stop(
      "MSE split failed: function=1 rows=", sum(is_ms1, na.rm = TRUE),
      " function=2 rows=", sum(is_ms2, na.rm = TRUE), call. = FALSE
    )
  }
  ms1 <- data.table::data.table(
    seqNum = as.integer(hdr$seqNum[is_ms1]),
    rt = as.numeric(hdr$retentionTime[is_ms1])
  )
  ms2 <- data.table::data.table(
    seqNum = as.integer(hdr$seqNum[is_ms2]),
    rt = as.numeric(hdr$retentionTime[is_ms2])
  )
  ms1 <- ms1[is.finite(rt)]
  ms2 <- ms2[is.finite(rt)]
  data.table::setorder(ms1, rt)
  data.table::setorder(ms2, rt)

  data.table::set(ms1, j = "peaks", value = mzR::peaks(m, ms1$seqNum))
  data.table::set(ms2, j = "peaks", value = mzR::peaks(m, ms2$seqNum))

  list(ms1 = ms1, ms2 = ms2, mode = mode, file = file)
}


# Flatten the scan list into a long data.table(scan_idx, rt, mz, intensity) -------
#' @noRd
.nts_dia_flatten <- function(scans) {
  if (!nrow(scans)) {
    return(data.table::data.table(
      scan_idx = integer(0), rt = numeric(0),
      mz = numeric(0), intensity = numeric(0)
    ))
  }
  npeaks <- vapply(scans$peaks, function(p) if (is.null(p)) 0L else nrow(p), integer(1L))
  if (sum(npeaks) == 0L) {
    return(data.table::data.table(
      scan_idx = integer(0), rt = numeric(0),
      mz = numeric(0), intensity = numeric(0)
    ))
  }
  big <- do.call(rbind, scans$peaks)
  data.table::data.table(
    scan_idx = rep.int(seq_len(nrow(scans)), npeaks),
    rt = rep.int(scans$rt, npeaks),
    mz = as.numeric(big[, 1L]),
    intensity = as.numeric(big[, 2L])
  )
}


# Single-feature co-elution extraction (matrix version) ----------------------------------------
# Within the RT window the MS1 EIC is one vector; candidate MS2 m/z bins form an
# (n_scan_ms2 x n_bin) matrix X in one shot. The Pearson correlation of each column
# of X against y (the MS1 EIC interpolated onto ms2_rt) is done in one step via a
# col-centered crossprod. In Waters MSE the MS1 and MS2 RTs interleave strictly
# 1:1, so when length(ms1) == length(ms2) they pair directly with no approx().
#' @noRd
.nts_dia_extract_one_feature <- function(feat,
                                          ms1_long, ms2_long,
                                          ms1_grid_dt,
                                          ms2_grid_dt,
                                          eic_tol_ppm,
                                          coelution_window_s,
                                          min_coelution,
                                          min_frag_intensity,
                                          frag_count_max,
                                          source_file) {
  rt_lo <- feat$rt - coelution_window_s
  rt_hi <- feat$rt + coelution_window_s
  ms1_grid <- ms1_grid_dt[rt >= rt_lo & rt <= rt_hi]
  ms2_grid <- ms2_grid_dt[rt >= rt_lo & rt <= rt_hi]
  if (nrow(ms1_grid) < 3L || nrow(ms2_grid) < 3L) {
    return(NULL)
  }

  tol_pre <- feat$mz * eic_tol_ppm / 1e6
  ms1_in <- ms1_long[
    scan_idx %in% ms1_grid$scan_idx &
      mz >= feat$mz - tol_pre & mz <= feat$mz + tol_pre
  ]
  if (!nrow(ms1_in)) {
    return(NULL)
  }
  pre_eic <- ms1_in[, .(intensity = max(intensity, na.rm = TRUE)), by = scan_idx]
  ms1_idx <- ms1_grid$scan_idx
  ms1_rt_seq <- ms1_grid$rt
  pre_vec_full <- numeric(length(ms1_idx))
  m_pre <- match(pre_eic$scan_idx, ms1_idx)
  pre_vec_full[m_pre] <- pre_eic$intensity
  # min_frag_intensity is a *fragment* apex threshold and lives on a different
  # intensity scale than the MS1 precursor EIC, so it must not gate the
  # precursor here. Require only that the precursor EIC has signal in at least
  # three scans; the co-elution test below does the real discrimination.
  if (sum(pre_vec_full > 0) < 3L) {
    return(NULL)
  }

  ms2_idx <- ms2_grid$scan_idx
  ms2_rt_seq <- ms2_grid$rt
  if (length(ms2_idx) == length(ms1_idx)) {
    pre_on_ms2 <- pre_vec_full
  } else {
    pre_on_ms2 <- stats::approx(
      x = ms1_rt_seq, y = pre_vec_full,
      xout = ms2_rt_seq, method = "linear",
      rule = 2, ties = mean
    )$y
    pre_on_ms2[!is.finite(pre_on_ms2)] <- 0
  }
  sd_pre <- stats::sd(pre_on_ms2)
  if (!is.finite(sd_pre) || sd_pre == 0) {
    return(NULL)
  }

  ms2_in <- ms2_long[
    scan_idx %in% ms2_idx &
      intensity >= min_frag_intensity &
      mz <= feat$mz - 1
  ]
  if (!nrow(ms2_in)) {
    return(NULL)
  }
  bin_da <- max(0.01, eic_tol_ppm * feat$mz / 1e6)
  ms2_in[, mz_bin := round(mz / bin_da) * bin_da]
  per_scan <- ms2_in[, .(intensity = max(intensity, na.rm = TRUE),
                         mz_apex = mz[which.max(intensity)]),
                     by = .(mz_bin, scan_idx)]
  bin_summary <- per_scan[, .(
    n_pos = .N,
    apex = max(intensity, na.rm = TRUE),
    rep_mz = mz_apex[which.max(intensity)]
  ), by = mz_bin]
  bin_summary <- bin_summary[n_pos >= 3L & apex >= min_frag_intensity]
  if (!nrow(bin_summary)) {
    return(NULL)
  }
  data.table::setorder(bin_summary, -apex)
  cand_cap <- max(frag_count_max * 4L, 32L)
  if (nrow(bin_summary) > cand_cap) {
    bin_summary <- bin_summary[seq_len(cand_cap)]
  }
  per_scan <- per_scan[mz_bin %in% bin_summary$mz_bin]

  bin_levels <- bin_summary$mz_bin
  k <- length(bin_levels)
  n <- length(ms2_idx)
  bin_pos <- match(per_scan$mz_bin, bin_levels)
  scan_pos <- match(per_scan$scan_idx, ms2_idx)
  ok <- !is.na(bin_pos) & !is.na(scan_pos)
  if (!any(ok)) {
    return(NULL)
  }
  X <- matrix(0, nrow = n, ncol = k)
  X[cbind(scan_pos[ok], bin_pos[ok])] <- per_scan$intensity[ok]

  apex_use <- matrixStats_colMax(X)
  col_means <- .colMeans(X, n, k)
  Xc <- X - rep(col_means, each = n)
  yc <- pre_on_ms2 - mean(pre_on_ms2)
  num <- as.numeric(crossprod(Xc, yc))
  den_X <- sqrt(.colSums(Xc * Xc, n, k))
  den_y <- sqrt(sum(yc * yc))
  bad <- !is.finite(den_X) | den_X == 0
  scores <- ifelse(bad, NA_real_, num / (den_X * den_y))

  rep_mz <- bin_summary$rep_mz
  keep <- is.finite(scores) & scores >= min_coelution & apex_use >= min_frag_intensity
  if (!any(keep)) {
    return(NULL)
  }
  rep_mz <- rep_mz[keep]
  scores <- scores[keep]
  apex_use <- apex_use[keep]
  if (length(rep_mz) > frag_count_max) {
    ord_top <- order(-apex_use)[seq_len(frag_count_max)]
    rep_mz <- rep_mz[ord_top]
    scores <- scores[ord_top]
    apex_use <- apex_use[ord_top]
  }
  ord <- order(rep_mz)
  rep_mz <- rep_mz[ord]
  scores <- scores[ord]
  apex_use <- apex_use[ord]

  list(list(
    feature_id = as.character(feat$name),
    id = as.character(feat$name),
    scan_idx = NA_integer_,
    source_file = source_file,
    precursorMz = as.numeric(feat$mz),
    rt = as.numeric(feat$rt),
    ce = NA_real_,
    polarity = NA_integer_,
    isolation_lo = NA_real_,
    isolation_hi = NA_real_,
    mz = as.numeric(rep_mz),
    intensity = as.numeric(apex_use),
    is_dia = TRUE,
    dia_mode = "waters_mse",
    coelution_score = as.numeric(scores)
  ))
}


# Tiny colMax without matrixStats dependency -----------------------------
#' @noRd
matrixStats_colMax <- function(X) {
  if (!ncol(X)) {
    return(numeric(0))
  }
  apply(X, 2L, max)
}


#' Extract Waters MSE (DIA) MS2 spectra by EIC co-elution
#'
#' Independent DIA entry to be used **instead of** [extract_ms2()] when the data
#' was acquired in Waters MSE mode. For each feature in `@ms1_annotations`, the
#' MS1 EIC at `feat$mz` is correlated (Pearson) with the EIC of every candidate
#' high-energy fragment ion in the same RT window; fragments above
#' `min_coelution` are kept and written to `@ms2data` in the **same shape** as
#' [extract_ms2()] so downstream [purify_ms2_spectra()] / [compute_ms2_qc()]
#' / network steps need no changes.
#'
#' @details
#' Waters MSE mzML files (msconvert default) report every spectrum as
#' \code{msLevel = 1}; low-energy and high-energy scans are alternated and
#' distinguished only by the \code{function=N} field embedded in
#' \code{spectrumId}. This function parses that field; it errors out if the
#' field is absent (e.g. for true DDA / non-Waters files).
#'
#' Co-elution metric: per-feature Pearson correlation between the MS1 EIC and
#' the high-energy EIC, both sampled at the MS1 retention-time grid via
#' linear interpolation. Fragment intensity reported in `@ms2data` is the
#' interpolated apex intensity inside the RT window.
#'
#' @param project An [NTSProject] after [extract_features()] and an MS1
#'   annotation (e.g. from `flag_isotopes` / `flag_isf`).
#' @param ms2_files Character vector of Waters MSE mzML files. **MS1 and MS2
#'   live in the same file** (Waters MSE convention).
#' @param mode Currently only `"waters_mse"`.
#' @param mz_tol_ppm Feature precursor m/z tolerance (currently informative).
#' @param rt_tol Reserved for future use; currently `coelution_window_s`
#'   controls the RT window.
#' @param eic_tol_ppm m/z tolerance (ppm) used when extracting MS1 / MS2 EICs.
#' @param coelution_window_s Half-width (seconds in RT units of `feature_table`)
#'   of the RT window used for the EIC pair.
#' @param min_coelution Minimum Pearson correlation between the MS1 EIC and a
#'   candidate fragment EIC for that fragment to be kept (default `0.7`).
#' @param min_frag_intensity Minimum apex intensity for fragment retention.
#' @param frag_count_max Maximum number of fragments per feature (top by apex
#'   intensity).
#' @param verbose Progress messages.
#'
#' @return Updated `project` with `@ms2data` populated. Each element has the
#'   fields produced by [extract_ms2()] plus `is_dia = TRUE`, `dia_mode`, and
#'   `coelution_score` (numeric vector aligned with `mz`). `@clean_ms2` is
#'   reset; run [purify_ms2_spectra()] next.
#'
#' @seealso [extract_ms2()], [purify_ms2_spectra()]
#' @export
extract_ms2_dia <- function(project,
                            ms2_files,
                            mode = c("waters_mse"),
                            mz_tol_ppm = 15,
                            rt_tol = 30,
                            eic_tol_ppm = 10,
                            coelution_window_s = 6,
                            min_coelution = 0.7,
                            min_frag_intensity = 100,
                            frag_count_max = 200,
                            verbose = TRUE) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  mode <- match.arg(mode)
  ms1 <- get_ms1_annotations(project)
  if (!is.data.frame(ms1) || !nrow(ms1)) {
    stop("ms1_annotations is empty.", call. = FALSE)
  }
  if (!all(c("name", "mz", "rt") %in% names(ms1))) {
    stop("ms1_annotations must contain columns name, mz, rt.", call. = FALSE)
  }
  ft_dt <- data.table::as.data.table(ms1)[, .(
    name = as.character(name),
    mz = as.numeric(mz),
    rt = as.numeric(rt)
  )]
  ft_dt <- ft_dt[is.finite(mz) & is.finite(rt)]
  if (!nrow(ft_dt)) {
    stop("No finite mz/rt rows in ms1_annotations.", call. = FALSE)
  }

  ms2_files <- as.character(ms2_files)
  if (!length(ms2_files)) {
    stop("ms2_files is empty.", call. = FALSE)
  }
  for (f in ms2_files) {
    if (!file.exists(f)) {
      stop("MSE mzML not found: ", f, call. = FALSE)
    }
  }

  acc <- list()
  for (fi in seq_along(ms2_files)) {
    f <- ms2_files[fi]
    if (isTRUE(verbose)) {
      message(sprintf("[%d/%d] Reading Waters MSE: %s",
                      fi, length(ms2_files), basename(f)))
    }
    sp <- .nts_dia_split_scans(f, mode = mode)

    if (isTRUE(verbose)) {
      message(sprintf("  flattening peaks: MS1 scans=%d  MS2 scans=%d",
                      nrow(sp$ms1), nrow(sp$ms2)))
    }
    ms1_long <- .nts_dia_flatten(sp$ms1)
    ms2_long <- .nts_dia_flatten(sp$ms2)
    data.table::setkey(ms1_long, scan_idx, mz)
    data.table::setkey(ms2_long, scan_idx, mz)
    ms1_grid_dt <- data.table::data.table(
      scan_idx = seq_len(nrow(sp$ms1)),
      rt = as.numeric(sp$ms1$rt)
    )
    ms2_grid_dt <- data.table::data.table(
      scan_idx = seq_len(nrow(sp$ms2)),
      rt = as.numeric(sp$ms2$rt)
    )

    rt_min <- min(ms1_grid_dt$rt, ms2_grid_dt$rt, na.rm = TRUE)
    rt_max <- max(ms1_grid_dt$rt, ms2_grid_dt$rt, na.rm = TRUE)
    feats <- ft_dt[rt >= rt_min - coelution_window_s &
                     rt <= rt_max + coelution_window_s]
    if (!nrow(feats)) {
      next
    }
    if (isTRUE(verbose)) {
      message(sprintf("  features in RT range=%d   ms1_long=%d rows  ms2_long=%d rows",
                      nrow(feats), nrow(ms1_long), nrow(ms2_long)))
    }
    n_feats <- nrow(feats)
    next_msg <- max(1L, as.integer(n_feats / 10L))
    for (i in seq_len(n_feats)) {
      r <- feats[i]
      one <- .nts_dia_extract_one_feature(
        feat = list(name = r$name, mz = r$mz, rt = r$rt),
        ms1_long = ms1_long, ms2_long = ms2_long,
        ms1_grid_dt = ms1_grid_dt, ms2_grid_dt = ms2_grid_dt,
        eic_tol_ppm = eic_tol_ppm,
        coelution_window_s = coelution_window_s,
        min_coelution = min_coelution,
        min_frag_intensity = min_frag_intensity,
        frag_count_max = frag_count_max,
        source_file = basename(f)
      )
      if (!is.null(one)) {
        key <- as.character(r$name)
        acc[[key]] <- c(acc[[key]], one)
      }
      if (isTRUE(verbose) && (i %% next_msg == 0L)) {
        message(sprintf("    progress %d/%d  features with MS2 so far=%d",
                        i, n_feats, length(acc)))
      }
    }
    if (isTRUE(verbose)) {
      message(sprintf("  cumulative features with MS2: %d", length(acc)))
    }
  }

  acc <- .nts_ms2data_drop_empty(acc)
  methods::slot(project, "ms2data") <- acc
  methods::slot(project, "clean_ms2") <- list()
  project <- .nts_build_ms2_annotations_from_ms2data(project, acc)
  if (length(methods::slot(project, "ms2data")) > nrow(get_ms1_annotations(project))) {
    stop("Fatal Error: Extracted more DIA MS2 spectra than valid MS1 annotations!",
         call. = FALSE)
  }
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "extract_ms2_dia",
    list(
      mode = mode,
      n_files = length(ms2_files),
      mz_tol_ppm = mz_tol_ppm,
      eic_tol_ppm = eic_tol_ppm,
      coelution_window_s = coelution_window_s,
      min_coelution = min_coelution,
      min_frag_intensity = min_frag_intensity,
      n_features_with_ms2 = length(acc),
      n_frag_total = sum(vapply(acc, function(z) {
        sum(vapply(z, function(s) length(s$mz), integer(1L)))
      }, integer(1L)))
    )
  )
  .nts_ms2_strip_mz_bounds(project)
}
