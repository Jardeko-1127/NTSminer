# Wide-window EIC trapezoid areas for feature_table (MassHunter-aligned optional path).

#' @keywords internal
#' @noRd
.nts_trap_integrate <- function(rt_s, y) {
  if (length(rt_s) < 2L || length(y) < 2L) {
    return(NA_real_)
  }
  sum(diff(rt_s) * (head(y, -1L) + tail(y, -1L)) / 2)
}

#' @keywords internal
#' @noRd
.nts_validate_feature_area_method <- function(method) {
  method <- match.arg(
    method,
    c("xcms_default", "eic_rt30", "apex_into")
  )
  method
}

#' @keywords internal
#' @noRd
.nts_validate_feature_area_rt_center <- function(rt_center) {
  match.arg(rt_center, c("apex", "rtmed"))
}

#' Per-feature per-sample apex chromPeak stats from `peakidx`.
#' @keywords internal
#' @noRd
.nts_apex_peak_stats <- function(xdata, sample_names = NULL) {
  fd <- xcms::featureDefinitions(xdata)
  n_f <- nrow(fd)
  if (!n_f) {
    return(list(into = matrix(), rt = matrix(), maxo = matrix()))
  }
  feat_names <- rownames(fd)
  if (is.null(feat_names) || !length(feat_names)) {
    feat_names <- sprintf("FT%06d", seq_len(n_f))
  }
  xdata_sn <- MSnbase::sampleNames(xdata)
  if (is.null(sample_names)) {
    sample_names <- xdata_sn
  } else {
    sample_names <- intersect(as.character(sample_names), xdata_sn)
  }
  n_s <- length(sample_names)
  into_mat <- matrix(NA_real_, n_f, n_s)
  rt_mat <- matrix(NA_real_, n_f, n_s)
  maxo_mat <- matrix(NA_real_, n_f, n_s)
  colnames(into_mat) <- sample_names
  rownames(into_mat) <- feat_names
  colnames(rt_mat) <- sample_names
  rownames(rt_mat) <- feat_names
  colnames(maxo_mat) <- sample_names
  rownames(maxo_mat) <- feat_names

  pkix <- fd[["peakidx"]]
  cp <- xcms::chromPeaks(xdata)
  cp_n <- nrow(cp)
  if (!cp_n || is.null(pkix)) {
    return(list(into = into_mat, rt = rt_mat, maxo = maxo_mat))
  }

  for (fi in seq_len(n_f)) {
    idx <- suppressWarnings(as.integer(pkix[[fi]]))
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= cp_n]
    if (!length(idx)) {
      next
    }
    sub <- cp[idx, , drop = FALSE]
    for (j in seq_len(n_s)) {
      si <- match(sample_names[j], xdata_sn)
      if (is.na(si)) {
        next
      }
      rows <- which(sub[, "sample"] == si)
      if (!length(rows)) {
        next
      }
      ss <- sub[rows, , drop = FALSE]
      into_vals <- ss[, "into"]
      best <- which.max(into_vals)
      into_mat[fi, j] <- into_vals[best]
      rt_mat[fi, j] <- ss[best, "rt"]
      maxo_mat[fi, j] <- ss[best, "maxo"]
    }
  }
  list(into = into_mat, rt = rt_mat, maxo = maxo_mat)
}

#' RT center matrix for EIC integration (`rtmed` or apex RT per sample).
#' @keywords internal
#' @noRd
.nts_feature_rt_center_matrix <- function(xdata,
                                         rt_center = c("apex", "rtmed"),
                                         sample_names = NULL) {
  rt_center <- .nts_validate_feature_area_rt_center(rt_center)
  fd <- xcms::featureDefinitions(xdata)
  n_f <- nrow(fd)
  if (!n_f) {
    return(matrix())
  }
  feat_names <- rownames(fd)
  if (is.null(feat_names) || !length(feat_names)) {
    feat_names <- sprintf("FT%06d", seq_len(n_f))
  }
  xdata_sn <- MSnbase::sampleNames(xdata)
  if (is.null(sample_names)) {
    sample_names <- xdata_sn
  }
  n_s <- length(sample_names)
  mat <- matrix(NA_real_, n_f, n_s)
  colnames(mat) <- sample_names
  rownames(mat) <- feat_names
  rtmed <- as.numeric(fd[["rtmed"]])
  if (rt_center == "rtmed") {
    for (j in seq_len(n_s)) {
      mat[, j] <- rtmed
    }
    return(mat)
  }
  apex <- .nts_apex_peak_stats(xdata, sample_names)
  mat <- apex$rt
  miss <- !is.finite(mat)
  if (any(miss)) {
    mat[miss] <- rtmed[row(mat)[miss]]
  }
  mat
}

#' Batch EIC trap integration for one sample (one scan loop, many features).
#' @keywords internal
#' @noRd
.nts_eic_trap_one_sample <- function(xsub,
                                    mz_vec,
                                    rt_center_vec,
                                    ppm = 5,
                                    rt_half = 30) {
  n_f <- length(mz_vec)
  out <- rep(NA_real_, n_f)
  if (!n_f) {
    return(out)
  }
  valid <- is.finite(mz_vec) & is.finite(rt_center_vec)
  if (!any(valid)) {
    return(out)
  }
  vi <- which(valid)
  tol <- mz_vec[vi] * ppm * 1e-6
  mz_range <- cbind(mz_vec[vi] - tol, mz_vec[vi] + tol)
  rt_range <- cbind(rt_center_vec[vi] - rt_half, rt_center_vec[vi] + rt_half)

  chrs <- tryCatch(
    MSnbase::chromatogram(xsub, rt = rt_range, mz = mz_range,
                          aggregationFun = "max", missing = 0),
    error = function(e) NULL
  )
  if (is.null(chrs)) {
    return(out)
  }

  for (k in seq_along(vi)) {
    fi <- vi[k]
    ints <- MSnbase::intensity(chrs[k, 1L])
    rts  <- MSnbase::rtime(chrs[k, 1L])
    out[fi] <- .nts_trap_integrate(rts, ints)
  }
  out
}

#' Feature x sample area matrix for `feature_table` sample columns.
#' @keywords internal
#' @noRd
.nts_feature_area_matrix <- function(xdata,
                                    method = c("xcms_default", "eic_rt30", "apex_into"),
                                    ppm = 5,
                                    rt_half = 30,
                                    rt_center = c("apex", "rtmed"),
                                    samples = NULL,
                                    BPPARAM = NULL) {
  method <- .nts_validate_feature_area_method(method)
  rt_center <- .nts_validate_feature_area_rt_center(rt_center)
  ppm <- suppressWarnings(as.numeric(ppm)[1L])
  rt_half <- suppressWarnings(as.numeric(rt_half)[1L])
  if (!is.finite(ppm) || ppm <= 0) {
    stop("`ppm` must be a positive finite number.", call. = FALSE)
  }
  if (!is.finite(rt_half) || rt_half <= 0) {
    stop("`rt_half` must be a positive finite number.", call. = FALSE)
  }

  fd <- xcms::featureDefinitions(xdata)
  n_f <- nrow(fd)
  xdata_sn <- MSnbase::sampleNames(xdata)
  sample_names <- xdata_sn
  if (!is.null(samples)) {
    sample_names <- intersect(as.character(samples), xdata_sn)
    if (!length(sample_names)) {
      stop("No `samples` match sampleNames(xdata).", call. = FALSE)
    }
  }
  n_s <- length(sample_names)
  feat_names <- rownames(fd)
  if (is.null(feat_names) || !length(feat_names)) {
    feat_names <- sprintf("FT%06d", seq_len(n_f))
  }

  fv_xcms <- as.matrix(xcms::featureValues(xdata, value = "into"))
  cn <- colnames(fv_xcms)
  if (is.null(cn) || anyNA(cn) || any(cn == "")) {
    cn <- xdata_sn[seq_len(ncol(fv_xcms))]
    colnames(fv_xcms) <- cn
  }
  rownames(fv_xcms) <- feat_names

  if (method == "xcms_default") {
    return(fv_xcms[, sample_names, drop = FALSE])
  }

  mat <- matrix(NA_real_, n_f, n_s)
  colnames(mat) <- sample_names
  rownames(mat) <- feat_names

  if (method == "apex_into") {
    apex <- .nts_apex_peak_stats(xdata, sample_names)
    return(apex$into)
  }

  mz_vec <- as.numeric(fd[["mzmed"]])
  valid_mz <- is.finite(mz_vec)
  if (!any(valid_mz)) {
    return(mat)
  }

  rt_centers <- .nts_feature_rt_center_matrix(
    xdata,
    rt_center = rt_center,
    sample_names = sample_names
  )
  use_parallel <- !is.null(BPPARAM) &&
    BiocParallel::bpnworkers(BPPARAM) > 1L &&
    n_s > 1L

  raw_files <- MSnbase::fileNames(xdata)

  if (use_parallel) {
    results <- BiocParallel::bplapply(seq_len(n_s), function(j) {
      si <- match(sample_names[j], xdata_sn)
      if (is.na(si)) {
        return(rep(NA_real_, n_f))
      }
      xraw <- MSnbase::readMSData(raw_files[si], msLevel. = 1L,
                                  mode = "onDisk")
      .nts_eic_trap_one_sample(
        xraw,
        mz_vec = mz_vec,
        rt_center_vec = rt_centers[, j],
        ppm = ppm,
        rt_half = rt_half
      )
    }, BPPARAM = BPPARAM)
    for (j in seq_len(n_s)) {
      mat[, j] <- results[[j]]
    }
  } else {
    for (j in seq_len(n_s)) {
      si <- match(sample_names[j], xdata_sn)
      if (is.na(si)) {
        next
      }
      xraw <- MSnbase::readMSData(raw_files[si], msLevel. = 1L,
                                  mode = "onDisk")
      mat[, j] <- .nts_eic_trap_one_sample(
        xraw,
        mz_vec = mz_vec,
        rt_center_vec = rt_centers[, j],
        ppm = ppm,
        rt_half = rt_half
      )
    }
  }
  mat
}

#' Resolve primary + xcms-default area matrices from `XcmsParam` settings.
#' @keywords internal
#' @noRd
.nts_resolve_feature_areas <- function(xdata, param) {
  method <- param$feature_area_method
  if (is.null(method) || !nzchar(as.character(method)[1L])) {
    method <- "xcms_default"
  }
  method <- .nts_validate_feature_area_method(method)
  ppm <- if (is.null(param$feature_area_ppm)) 5 else param$feature_area_ppm
  rt_half <- if (is.null(param$feature_area_rt_half)) 30 else param$feature_area_rt_half
  rt_center <- if (is.null(param$feature_area_rt_center)) {
    "apex"
  } else {
    param$feature_area_rt_center
  }
  samples <- param$feature_area_samples
  bpp <- param$feature_area_BPPARAM

  fv_xcms <- as.matrix(xcms::featureValues(xdata, value = "into"))
  cn <- colnames(fv_xcms)
  if (is.null(cn) || anyNA(cn) || any(cn == "")) {
    cn <- MSnbase::sampleNames(xdata)[seq_len(ncol(fv_xcms))]
    colnames(fv_xcms) <- cn
  }

  primary <- if (method == "xcms_default") {
    fv_xcms
  } else if (is.null(samples)) {
    .nts_feature_area_matrix(
      xdata,
      method = method,
      ppm = ppm,
      rt_half = rt_half,
      rt_center = rt_center,
      samples = NULL,
      BPPARAM = bpp
    )
  } else {
    out <- fv_xcms
    sub <- .nts_feature_area_matrix(
      xdata,
      method = method,
      ppm = ppm,
      rt_half = rt_half,
      rt_center = rt_center,
      samples = samples,
      BPPARAM = bpp
    )
    out[, colnames(sub)] <- sub
    out
  }

  list(
    matrix = primary,
    xcms_default = fv_xcms,
    method = method,
    ppm = as.numeric(ppm)[1L],
    rt_half = as.numeric(rt_half)[1L],
    rt_center = .nts_validate_feature_area_rt_center(rt_center)
  )
}

#' Attach feature-area provenance attributes on a feature bounds table.
#' @keywords internal
#' @noRd
.nts_set_feature_area_attrs <- function(ft, area_bundle) {
  if (is.null(area_bundle) || !is.list(area_bundle)) {
    return(ft)
  }
  attr(ft, "feature_area_method") <- area_bundle$method
  attr(ft, "feature_area_ppm") <- area_bundle$ppm
  attr(ft, "feature_area_rt_half") <- area_bundle$rt_half
  attr(ft, "feature_area_rt_center") <- area_bundle$rt_center
  if (!identical(area_bundle$method, "xcms_default")) {
    attr(ft, "xcms_default_intensity") <- area_bundle$xcms_default
  }
  ft
}

#' Refresh sample-column peak areas on `@feature_table`
#'
#' Recomputes per-feature per-sample intensities from `project@xdata` using
#' an alternate area definition aligned with wide-window EIC integration
#' (see [XcmsParam] fields `feature_area_*`). Default xcms values are kept
#' on attribute `xcms_default_intensity` when `method != "xcms_default"`.
#'
#' @param project An [NTSProject] with `@xdata` from [extract_features].
#' @param xdata Optional `XCMSnExp`; defaults to `project@xdata`.
#' @param method `"eic_rt30"` (default), `"apex_into"`, or `"xcms_default"`.
#' @param ppm m/z tolerance (ppm) for `eic_rt30` (default `5`).
#' @param rt_half Half-width of RT integration window in seconds (default `30`).
#' @param rt_center RT center for `eic_rt30`: `"apex"` (default) or `"rtmed"`.
#' @param samples Optional sample column names to update; default all samples
#'   present in both `feature_table` and `xdata`. For large studies, update QC
#'   or a subset first (`samples = c("QC-20-1.mzML", ...)`).
#' @param threads Parallel worker count when `BPPARAM` is `NULL` (default `1`,
#'   serial). Uses [get_robust_bpparam]; on Windows/RStudio this is
#'   [BiocParallel::SnowParam]. Parallelism is by sample; each worker holds a
#'   copy of `xdata`.
#' @param BPPARAM Optional [BiocParallel::BiocParallelParam]. When non-`NULL`,
#'   overrides `threads`. **Requires NTSminer installed** (`remotes::install_local()`);
#'   `load_all()` alone does not work on SnowParam workers.
#' @param keep_xcms_default If `TRUE` and `method != "xcms_default"`, store the
#'   previous sample intensities on `attr(feature_table, "xcms_default_intensity")`.
#' @details
#' `method = "eic_rt30"` extracts per-feature EICs via batched [MSnbase::chromatogram]
#' calls (500 features per batch) for a balance of speed and memory. Use
#' `threads > 1` to parallelise across samples.
#' @return Updated `project`.
#' @export
#' @seealso [XcmsParam], [extract_features]
refresh_feature_areas <- function(project,
                                  xdata = NULL,
                                  method = c("eic_rt30", "apex_into", "xcms_default"),
                                  ppm = 5,
                                  rt_half = 30,
                                  rt_center = c("apex", "rtmed"),
                                  samples = NULL,
                                  threads = 1L,
                                  BPPARAM = NULL,
                                  keep_xcms_default = TRUE) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  method <- .nts_validate_feature_area_method(match.arg(method))
  rt_center <- .nts_validate_feature_area_rt_center(match.arg(rt_center))
  threads <- suppressWarnings(as.integer(threads)[1L])
  if (!is.finite(threads) || threads < 1L) {
    threads <- 1L
  }
  bpp <- BPPARAM
  if (is.null(bpp)) {
    bpp <- get_robust_bpparam(threads)
  }

  xd <- xdata
  if (is.null(xd)) {
    xd <- methods::slot(project, "xdata")
  }
  if (is.null(xd)) {
    stop(
      "No xcms experiment: run extract_features() first or pass `xdata`.",
      call. = FALSE
    )
  }

  ft <- get_feature_table(project)
  if (!nrow(ft)) {
    stop("feature_table is empty.", call. = FALSE)
  }
  int_cols <- .nts_feature_intensity_colnames(ft)
  xdata_sn <- MSnbase::sampleNames(xd)
  if (is.null(samples)) {
    samples <- intersect(int_cols, xdata_sn)
  } else {
    samples <- intersect(as.character(samples), int_cols)
    samples <- intersect(samples, xdata_sn)
  }
  if (!length(samples)) {
    stop("No sample columns to update.", call. = FALSE)
  }

  nw <- BiocParallel::bpnworkers(bpp)
  message(
    "refresh_feature_areas: method=", method,
    if (method == "eic_rt30") {
      paste0(
        " (ppm=", ppm, ", RT+/-", rt_half,
        " s, center=", rt_center, ", n_samples=", length(samples),
        ", workers=", nw, ")"
      )
    } else {
      paste0(" (n_samples=", length(samples), ", workers=", nw, ")")
    }
  )

  xcms_prev <- NULL
  if (isTRUE(keep_xcms_default) && !identical(method, "xcms_default")) {
    xcms_prev <- attr(ft, "xcms_default_intensity")
    if (is.null(xcms_prev)) {
      xcms_prev <- as.matrix(ft[, samples, drop = FALSE])
      rownames(xcms_prev) <- ft$name
    }
  }

  param <- list(
    feature_area_method = method,
    feature_area_ppm = ppm,
    feature_area_rt_half = rt_half,
    feature_area_rt_center = rt_center,
    feature_area_samples = samples,
    feature_area_BPPARAM = bpp
  )
  area_bundle <- .nts_resolve_feature_areas(xd, param)
  new_mat <- area_bundle$matrix

  for (samp in samples) {
    if (!samp %in% colnames(new_mat)) {
      next
    }
    ft[[samp]] <- new_mat[match(ft$name, rownames(new_mat)), samp]
  }

  if (isTRUE(keep_xcms_default) && !identical(method, "xcms_default")) {
    area_bundle$xcms_default <- xcms_prev
  }
  ft <- .nts_set_feature_area_attrs(ft, area_bundle)
  set_feature_table(project, ft)
}

#' Apex chromPeak stats for a subset of feature names (fast path).
#' @keywords internal
#' @noRd
.nts_apex_peak_stats_for_names <- function(xdata,
                                          feature_names,
                                          sample_names = NULL) {
  fd <- xcms::featureDefinitions(xdata)
  feat_names <- rownames(fd)
  if (is.null(feat_names) || !length(feat_names)) {
    feat_names <- sprintf("FT%06d", seq_len(nrow(fd)))
  }
  want <- unique(as.character(feature_names))
  fi <- match(want, feat_names)
  fi <- fi[is.finite(fi) & fi >= 1L]
  if (!length(fi)) {
    return(list(into = matrix(), rt = matrix(), maxo = matrix()))
  }
  xdata_sn <- MSnbase::sampleNames(xdata)
  if (is.null(sample_names)) {
    sample_names <- xdata_sn
  }
  n_s <- length(sample_names)
  into_mat <- matrix(NA_real_, length(fi), n_s)
  rt_mat <- matrix(NA_real_, length(fi), n_s)
  maxo_mat <- matrix(NA_real_, length(fi), n_s)
  rownames(into_mat) <- feat_names[fi]
  colnames(into_mat) <- sample_names
  rownames(rt_mat) <- rownames(into_mat)
  colnames(rt_mat) <- sample_names
  rownames(maxo_mat) <- rownames(into_mat)
  colnames(maxo_mat) <- sample_names
  pkix <- fd[["peakidx"]]
  cp <- xcms::chromPeaks(xdata)
  cp_n <- nrow(cp)
  for (k in seq_along(fi)) {
    fidx <- fi[k]
    idx <- suppressWarnings(as.integer(pkix[[fidx]]))
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= cp_n]
    if (!length(idx)) {
      next
    }
    sub <- cp[idx, , drop = FALSE]
    for (j in seq_len(n_s)) {
      si <- match(sample_names[j], xdata_sn)
      if (is.na(si)) {
        next
      }
      rows <- which(sub[, "sample"] == si)
      if (!length(rows)) {
        next
      }
      ss <- sub[rows, , drop = FALSE]
      into_vals <- ss[, "into"]
      best <- which.max(into_vals)
      into_mat[k, j] <- into_vals[best]
      rt_mat[k, j] <- ss[best, "rt"]
      maxo_mat[k, j] <- ss[best, "maxo"]
    }
  }
  list(into = into_mat, rt = rt_mat, maxo = maxo_mat)
}

#' Resolve EIC m/z and RT center for db_annotations rows.
#' @keywords internal
#' @noRd
.nts_prepare_db_annotation_quant_tasks <- function(ann,
                                                  ft,
                                                  rt_center = c("exp_rt", "feature_rt", "lib_rt", "apex")) {
  rt_center <- match.arg(rt_center)
  if (!is.data.frame(ann) || !nrow(ann)) {
    return(data.frame())
  }
  if (!"name" %in% names(ann)) {
    stop("annotations must contain column `name` (feature id).", call. = FALSE)
  }
  ft_mz <- setNames(as.numeric(ft$mz), ft$name)
  ft_rt <- setNames(as.numeric(ft$rt), ft$name)
  lib_mz <- if ("lib_mz" %in% names(ann)) {
    suppressWarnings(as.numeric(ann[["lib_mz"]]))
  } else {
    rep(NA_real_, nrow(ann))
  }
  lib_rt <- if ("lib_rt" %in% names(ann)) {
    suppressWarnings(as.numeric(ann[["lib_rt"]]))
  } else {
    rep(NA_real_, nrow(ann))
  }
  exp_rt <- if ("exp_rt" %in% names(ann)) {
    suppressWarnings(as.numeric(ann[["exp_rt"]]))
  } else {
    rep(NA_real_, nrow(ann))
  }
  feat <- as.character(ann[["name"]])
  eic_mz <- lib_mz
  miss_mz <- !is.finite(eic_mz)
  if (any(miss_mz)) {
    eic_mz[miss_mz] <- ft_mz[feat[miss_mz]]
  }
  rt_feat <- ft_rt[feat]
  rt_center_val <- switch(
    rt_center,
    exp_rt = ifelse(is.finite(exp_rt), exp_rt, rt_feat),
    feature_rt = rt_feat,
    lib_rt = ifelse(is.finite(lib_rt), lib_rt, rt_feat),
    apex = rt_feat
  )
  data.frame(
    ann_row = seq_len(nrow(ann)),
    name = feat,
    ref_id = if ("ref_id" %in% names(ann)) {
      as.character(ann[["ref_id"]])
    } else if ("Lab_ID" %in% names(ann)) {
      as.character(ann[["Lab_ID"]])
    } else {
      NA_character_
    },
    compound = if ("Compound.name" %in% names(ann)) {
      as.character(ann[["Compound.name"]])
    } else {
      NA_character_
    },
    lib_mz = lib_mz,
    lib_rt = lib_rt,
    exp_rt = exp_rt,
    feature_mz = as.numeric(ft_mz[feat]),
    feature_rt = as.numeric(rt_feat),
    eic_mz = as.numeric(eic_mz),
    rt_center = as.numeric(rt_center_val),
    rt_center_source = rt_center,
    stringsAsFactors = FALSE
  )
}

#' EIC quantification for one sample (db_annotations subset).
#' @keywords internal
#' @noRd
.nts_quantify_db_annotations_one_sample <- function(j,
                                                   samp,
                                                   si,
                                                   raw_files,
                                                   tasks,
                                                   apex,
                                                   ft_into,
                                                   ppm,
                                                   rt_half,
                                                   rt_center_mode) {
  if (is.na(si)) {
    return(NULL)
  }
  rt_vec <- tasks$rt_center
  if (rt_center_mode == "apex") {
    ap_rt <- apex$rt[tasks$name, samp, drop = TRUE]
    rt_vec <- ap_rt
    miss <- !is.finite(rt_vec)
    if (any(miss)) {
      rt_vec[miss] <- tasks$feature_rt[miss]
    }
  }
  valid <- is.finite(tasks$eic_mz) & is.finite(rt_vec)
  if (!any(valid)) {
    return(NULL)
  }
  vi <- which(valid)
  tol <- tasks$eic_mz[vi] * ppm * 1e-6
  mz_range <- cbind(tasks$eic_mz[vi] - tol, tasks$eic_mz[vi] + tol)
  rt_range <- cbind(rt_vec[vi] - rt_half, rt_vec[vi] + rt_half)

  xraw <- MSnbase::readMSData(raw_files[si], msLevel. = 1L, mode = "onDisk")
  chrs <- tryCatch(
    MSnbase::chromatogram(xraw, rt = rt_range, mz = mz_range,
                          aggregationFun = "max", missing = 0),
    error = function(e) NULL
  )

  areas <- rep(NA_real_, nrow(tasks))
  if (!is.null(chrs)) {
    for (k in seq_along(vi)) {
      ti <- vi[k]
      ints <- MSnbase::intensity(chrs[k, 1L])
      rts  <- MSnbase::rtime(chrs[k, 1L])
      areas[ti] <- .nts_trap_integrate(rts, ints)
    }
  }

  ap_into <- apex$into[tasks$name, samp, drop = TRUE]
  ft_val <- ft_into[cbind(match(tasks$name, rownames(ft_into)), j)]
  data.frame(
    name = tasks$name,
    ref_id = tasks$ref_id,
    compound = tasks$compound,
    lib_mz = tasks$lib_mz,
    lib_rt = tasks$lib_rt,
    exp_rt = tasks$exp_rt,
    feature_mz = tasks$feature_mz,
    feature_rt = tasks$feature_rt,
    eic_mz = tasks$eic_mz,
    rt_center = rt_vec,
    rt_center_source = tasks$rt_center_source,
    sample = samp,
    eic_rt30_area = areas,
    feature_table_into = ft_val,
    apex_max_into = ap_into,
    stringsAsFactors = FALSE
  )
}

#' bplapply with txtProgressBar updates on the master (batched by worker count).
#' @keywords internal
#' @noRd
.nts_bplapply_with_progress <- function(X, FUN, BPPARAM, label = "progress") {
  n <- length(X)
  if (!n) {
    return(list())
  }
  bpp <- BPPARAM
  if (methods::is(bpp, "SnowParam") || methods::is(bpp, "MulticoreParam")) {
    BiocParallel::bpprogressbar(bpp) <- TRUE
  }
  workers <- max(1L, BiocParallel::bpnworkers(bpp))
  batch_n <- max(1L, min(n, workers))
  pb <- .nts_progress_start(n, label)
  on.exit(.nts_progress_done(pb), add = TRUE)
  if (n <= batch_n) {
    res <- BiocParallel::bplapply(X, FUN, BPPARAM = bpp)
    .nts_progress_update(pb, n)
    return(res)
  }
  res <- vector("list", n)
  pos <- 1L
  while (pos <= n) {
    end <- min(pos + batch_n - 1L, n)
    batch <- X[pos:end]
    batch_res <- BiocParallel::bplapply(batch, FUN, BPPARAM = bpp)
    res[pos:end] <- batch_res
    .nts_progress_update(pb, end)
    pos <- end + 1L
  }
  res
}

#' Wide-window EIC areas for [batch_annotate] hits (MassHunter-aligned)
#'
#' Computes **RT+/-30 s** (default) EIC trapezoid areas for rows in `@db_annotations`
#' without modifying the [NTSProject]. Intended use: one line after
#' [batch_annotate], e.g.
#' `quant_areas <- quantify_db_annotations(project)`.
#'
#' Only **annotated hits** are processed (typically tens of rows), not the full
#' feature matrix, so runtime is comparable to the standalone `mh_peak_check`
#' export scripts.
#'
#' @param project An [NTSProject] with `@xdata` and non-empty `@db_annotations`.
#' @param annotations Optional `data.frame`; defaults to [get_db_annotations].
#' @param xdata Optional `XCMSnExp`; defaults to `project@xdata`.
#' @param samples Sample columns to quantify; default all intensity columns in
#'   [get_feature_table] that match `sampleNames(xdata)`.
#' @param ppm m/z tolerance for EIC extraction (default `5`).
#' @param rt_half RT half-window in seconds (default `30`).
#' @param rt_center RT center: `"exp_rt"` (default, from annotation then feature),
#'   `"feature_rt"`, `"lib_rt"`, or `"apex"` (per-sample apex chromPeak RT).
#' @param threads Parallel worker count when `BPPARAM` is `NULL` (default `1`,
#'   serial). Uses [get_robust_bpparam]; on Windows/RStudio this is
#'   [BiocParallel::SnowParam].
#' @param BPPARAM Optional [BiocParallel::BiocParallelParam]. When non-`NULL`,
#'   overrides `threads`. Parallelism is **by sample** (`filterFile` per worker);
#'   each worker holds a copy of `xdata`, so keep `threads` modest on Windows.
#'   **Requires NTSminer installed** (`remotes::install_local()`); `load_all()`
#'   alone does not work on SnowParam workers. Progress uses `txtProgressBar`
#'   (updated after each worker batch).
#' @return Long `data.frame`: one row per annotation x sample with
#'   `eic_rt30_area`, `feature_table_into`, `apex_max_into`, and annotation ids.
#' @export
#' @seealso [batch_annotate], [refresh_feature_areas], [XcmsParam]
quantify_db_annotations <- function(project,
                                    annotations = NULL,
                                    xdata = NULL,
                                    samples = NULL,
                                    ppm = 5,
                                    rt_half = 30,
                                    rt_center = c("exp_rt", "feature_rt", "lib_rt", "apex"),
                                    threads = 1L,
                                    BPPARAM = NULL) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  rt_center <- match.arg(rt_center)
  ppm <- suppressWarnings(as.numeric(ppm)[1L])
  rt_half <- suppressWarnings(as.numeric(rt_half)[1L])
  threads <- suppressWarnings(as.integer(threads)[1L])
  if (!is.finite(ppm) || ppm <= 0) {
    stop("`ppm` must be a positive finite number.", call. = FALSE)
  }
  if (!is.finite(rt_half) || rt_half <= 0) {
    stop("`rt_half` must be a positive finite number.", call. = FALSE)
  }
  if (!is.finite(threads) || threads < 1L) {
    stop("`threads` must be a positive integer.", call. = FALSE)
  }
  bpp <- BPPARAM
  if (is.null(bpp)) {
    bpp <- get_robust_bpparam(threads)
  }

  ann <- annotations
  if (is.null(ann)) {
    ann <- get_db_annotations(project)
  }
  if (!is.data.frame(ann) || !nrow(ann)) {
    return(data.frame(
      name = character(0),
      ref_id = character(0),
      compound = character(0),
      sample = character(0),
      eic_rt30_area = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  xd <- xdata
  if (is.null(xd)) {
    xd <- methods::slot(project, "xdata")
  }
  if (is.null(xd)) {
    stop("No xcms experiment: run extract_features() first or pass `xdata`.", call. = FALSE)
  }

  ft <- get_feature_table(project)
  int_cols <- .nts_feature_intensity_colnames(ft)
  xdata_sn <- MSnbase::sampleNames(xd)
  if (is.null(samples)) {
    samples <- intersect(int_cols, xdata_sn)
  } else {
    samples <- intersect(as.character(samples), int_cols)
    samples <- intersect(samples, xdata_sn)
  }
  if (!length(samples)) {
    stop("No sample columns to quantify.", call. = FALSE)
  }
  use_parallel <- BiocParallel::bpnworkers(bpp) > 1L && length(samples) > 1L

  tasks <- .nts_prepare_db_annotation_quant_tasks(ann, ft, rt_center = rt_center)
  if (!nrow(tasks)) {
    return(data.frame())
  }

  apex <- .nts_apex_peak_stats_for_names(xd, unique(tasks$name), samples)
  ft_into <- as.matrix(ft[, samples, drop = FALSE])
  rownames(ft_into) <- ft$name

  if (rt_center == "apex") {
    tasks$rt_center <- NA_real_
  }

  raw_files <- MSnbase::fileNames(xd)
  sample_si <- match(samples, xdata_sn)
  quant_one <- function(j) {
    if (use_parallel) {
      .nts_parallel_call_internal(
        ".nts_quantify_db_annotations_one_sample",
        j = j,
        samp = samples[j],
        si = sample_si[j],
        raw_files = raw_files,
        tasks = tasks,
        apex = apex,
        ft_into = ft_into,
        ppm = ppm,
        rt_half = rt_half,
        rt_center_mode = rt_center
      )
    } else {
      .nts_quantify_db_annotations_one_sample(
        j = j,
        samp = samples[j],
        si = sample_si[j],
        raw_files = raw_files,
        tasks = tasks,
        apex = apex,
        ft_into = ft_into,
        ppm = ppm,
        rt_half = rt_half,
        rt_center_mode = rt_center
      )
    }
  }

  if (use_parallel) {
    .nts_assert_parallel_ntsminer(c(
      ".nts_quantify_db_annotations_one_sample",
      ".nts_trap_integrate"
    ))
    res <- .nts_bplapply_with_progress(
      seq_along(samples),
      quant_one,
      BPPARAM = bpp,
      label = "quantify_db_annotations"
    )
    out_list <- res[!vapply(res, is.null, logical(1L))]
  } else {
    out_list <- vector("list", length(samples))
    ok <- 0L
    pb <- .nts_progress_start(length(samples), "quantify_db_annotations")
    on.exit(.nts_progress_done(pb), add = TRUE)
    for (j in seq_along(samples)) {
      chunk <- quant_one(j)
      if (!is.null(chunk)) {
        ok <- ok + 1L
        out_list[[ok]] <- chunk
      }
      .nts_progress_update(pb, j)
    }
    if (!ok) {
      out_list <- list()
    } else {
      out_list <- out_list[seq_len(ok)]
    }
  }

  if (!length(out_list)) {
    return(data.frame())
  }
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

#' Pivot EIC areas to a wide matrix for comparison
#'
#' Reshapes the long-format output of [quantify_db_annotations] into a wide
#' matrix suitable for direct comparison with manual integration results
#' (e.g. MassHunter). Two formats are available:
#' - `"compound_sample"` (default): rows = compounds, columns = samples.
#'   Optionally joins feature metadata (`name`, `mz`, `rt`) from a
#'   `feature_table`-like `data.frame`.
#' - `"feature_compound"`: rows = features, columns = compounds.
#'
#' @param quant_areas `data.frame` from [quantify_db_annotations].
#' @param format `"compound_sample"` or `"feature_compound"`.
#' @param value_col Column for cell values (default `"eic_rt30_area"`).
#' @param feature_table Optional `data.frame` with columns `name`, `mz`, `rt`
#'   (e.g. [get_feature_table]). When non-`NULL` and `format =
#'   `"compound_sample"`, left-joins these columns onto each compound row.
#' @return `data.frame` in wide format.
#' @export
#' @seealso [quantify_db_annotations]
#' @examples
#' \dontrun{
#' quant_areas <- quantify_db_annotations(project)
#' # Compound x sample (like feature_table)
#' wide <- pivot_quant_areas(quant_areas)
#' # Compound x sample + metadata
#' wide <- pivot_quant_areas(quant_areas, feature_table = get_feature_table(project))
#' # Feature x compound
#' wide <- pivot_quant_areas(quant_areas, format = "feature_compound")
#' }
pivot_quant_areas <- function(quant_areas,
                              format = c("compound_sample", "feature_compound"),
                              value_col = "eic_rt30_area",
                              feature_table = NULL) {
  if (!is.data.frame(quant_areas) || !nrow(quant_areas)) {
    return(data.frame(name = character(0L)))
  }
  format <- match.arg(format)
  needed <- c("name", "compound", "sample", value_col)
  miss <- setdiff(needed, names(quant_areas))
  if (length(miss)) {
    stop("quant_areas missing columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  .agg <- function(x) {
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else mean(x)
  }

  dt <- data.table::as.data.table(quant_areas)

  if (format == "feature_compound") {
    frm <- stats::as.formula("name ~ compound")
    out <- data.table::dcast(dt, frm,
      value.var = value_col,
      fun.aggregate = .agg,
      fill = NA_real_
    )
  } else {
    frm <- stats::as.formula("compound ~ sample")
    out <- data.table::dcast(dt, frm,
      value.var = value_col,
      fun.aggregate = .agg,
      fill = NA_real_
    )

    if (!is.null(feature_table) && is.data.frame(feature_table)) {
      ft <- data.table::as.data.table(feature_table)
      meta <- c("name", "mz", "rt")
      have <- intersect(meta, names(ft))
      if (length(have)) {
        fmap <- unique(dt[, list(name, compound)])
        fmap <- fmap[!is.na(name) & !is.na(compound)]
        fmap <- fmap[fmap[, .I[1L], by = compound]$V1]
        ft_sub <- ft[, c("name", setdiff(have, "name")), with = FALSE]
        fmap <- merge(fmap, ft_sub, by = "name", all.x = TRUE, sort = FALSE)
        out <- merge(fmap, out, by = "compound", all.y = TRUE, sort = FALSE)
        col_order <- c("compound", have, setdiff(names(out), c("compound", have)))
        out <- out[, col_order, with = FALSE]
      }
    }
  }

  as.data.frame(out, stringsAsFactors = FALSE)
}
