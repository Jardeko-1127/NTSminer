#' Extract XCMSnExp from project (with fallback to on-disk reload)
#' @noRd
.nts_get_xdata <- function(project, files = NULL) {
  xd <- methods::slot(project, "xdata")
  if (!is.null(xd)) {
    return(xd)
  }
  if (is.null(files) || !length(files)) {
    files <- .nts_list_ms_files(methods::slot(project, "raw_data_path"))
  }
  if (!length(files)) {
    stop("No mzML/mzXML files found and @xdata is empty", call. = FALSE)
  }
  message("> @xdata is empty; reloading from mzML (no RT correction in this pass)...")
  MSnbase::readMSData(files = files, mode = "onDisk", msLevel = 1L)
}

#' Chromatogram object -> tidy data.frame
#' @noRd
.chrom_to_df <- function(chrom_obj, files = NULL) {
  nc <- ncol(chrom_obj)
  frames <- vector("list", nc)
  for (j in seq_len(nc)) {
    chr <- chrom_obj[1, j]
    rt_v <- MSnbase::rtime(chr)
    int_v <- MSnbase::intensity(chr)
    ok <- is.finite(rt_v) & is.finite(int_v)
    if (!any(ok)) next
    sn <- if (!is.null(files) && length(files) >= j) {
      basename(files[j])
    } else {
      paste0("sample_", j)
    }
    sn <- sub("\\.[Mm][Zz][MmXx][Ll]$", "", sn)
    frames[[j]] <- data.frame(
      rt = rt_v[ok],
      intensity = int_v[ok],
      sample = sn,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, frames[!vapply(frames, is.null, logical(1L))])
}

#' Plot Base Peak Chromatogram (BPC)
#'
#' Uses `project@xdata` (inherits RT alignment from [extract_features]). Falls
#' back to on-disk reload if `@xdata` is empty.
#'
#' @param project An [NTSProject].
#' @param files Optional character vector (only used for fallback / labels).
#'
#' @return A `ggplot2` object.
#' @export
plot_bpc <- function(project, files = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_bpc()", call. = FALSE)
  }
  stopifnot(inherits(project, "NTSProject"))

  xd <- .nts_get_xdata(project, files)
  src_files <- tryCatch(
    MSnbase::fileNames(xd),
    error = function(e) files
  )

  message("> Extracting base peak chromatogram (BPC)...")
  bpc <- xcms::chromatogram(xd, aggregationFun = "max")

  df <- .chrom_to_df(bpc, src_files)
  if (is.null(df) || !nrow(df)) {
    stop("BPC extraction returned no data", call. = FALSE)
  }

  n_samples <- length(unique(df$sample))
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["rt"]],
    y = .data[["intensity"]],
    colour = .data[["sample"]]
  )) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.85) +
    ggplot2::labs(
      title = "Base Peak Chromatogram (BPC)",
      x = "Retention Time (s)",
      y = "Intensity"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = if (n_samples > 8) "none" else "right",
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )
  p
}

#' Build a single EIC trace spec from an `ms1_annotations` / `feature_table` row
#' @noRd
.nts_eic_trace_from_row <- function(row,
                                    mz_tol,
                                    rt_tol,
                                    label,
                                    use_bounds,
                                    is_target) {
  pick <- function(nm) {
    if (!nm %in% names(row)) return(NA_real_)
    suppressWarnings(as.numeric(row[[nm]])[1L])
  }
  mz <- pick("mz")
  rt <- pick("rt")
  mzmin <- pick("mzmin")
  mzmax <- pick("mzmax")
  rtmin <- pick("rtmin")
  rtmax <- pick("rtmax")
  if (!isTRUE(use_bounds) ||
        !is.finite(mzmin) || !is.finite(mzmax) || mzmax <= mzmin) {
    mzmin <- mz - mz_tol
    mzmax <- mz + mz_tol
  }
  peak_rtmin <- if (isTRUE(use_bounds) && is.finite(rtmin) && is.finite(rtmax) &&
                     rtmax > rtmin) {
    rtmin
  } else {
    NA_real_
  }
  peak_rtmax <- if (isTRUE(use_bounds) && is.finite(rtmin) && is.finite(rtmax) &&
                     rtmax > rtmin) {
    rtmax
  } else {
    NA_real_
  }
  if (is.finite(peak_rtmin) && is.finite(peak_rtmax)) {
    rtmin_eff <- peak_rtmin - rt_tol
    rtmax_eff <- peak_rtmax + rt_tol
  } else {
    rtmin_eff <- rt - rt_tol
    rtmax_eff <- rt + rt_tol
  }
  list(
    name = as.character(row[["name"]])[1L],
    label = label,
    mz = mz,
    rt = rt,
    mz_min = mzmin,
    mz_max = mzmax,
    rt_min = rtmin_eff,
    rt_max = rtmax_eff,
    peak_rt_min = peak_rtmin,
    peak_rt_max = peak_rtmax,
    is_target = isTRUE(is_target)
  )
}

#' Collect EIC traces for `feature_name` + (optional) isotope siblings
#' @noRd
.nts_eic_collect_traces <- function(project,
                                    feature_name,
                                    mz_tol,
                                    rt_tol,
                                    show_isotopes,
                                    include_parent_if_child,
                                    use_bounds) {
  feature_name <- as.character(feature_name)[1L]
  ft <- get_feature_table(project)
  if (!nrow(ft) || !"name" %in% names(ft)) {
    stop("feature_table is empty or lacks 'name' column", call. = FALSE)
  }
  ann <- get_ms1_annotations(project)
  has_ann <- is.data.frame(ann) && nrow(ann) && "name" %in% names(ann)

  ft_idx <- match(feature_name, as.character(ft$name))
  if (is.na(ft_idx)) {
    stop("feature_name '", feature_name,
         "' not found in feature_table", call. = FALSE)
  }

  if (has_ann) {
    aidx <- match(feature_name, as.character(ann$name))
    if (is.na(aidx)) {
      row_target <- ft[ft_idx, , drop = FALSE]
      target_iso <- FALSE
      parent_name <- feature_name
      target_label <- "M (target)"
    } else {
      row_target <- ann[aidx, , drop = FALSE]
      target_iso <- isTRUE(as.logical(row_target[["is_isotope"]])[1L] %in% TRUE)
      raw_lab <- as.character(row_target[["isotope_label"]])[1L]
      if (target_iso && !is.na(raw_lab) && nzchar(raw_lab)) {
        target_label <- paste0(raw_lab, " (target)")
        pn <- as.character(row_target[["parent_name"]])[1L]
        parent_name <- if (!is.na(pn) && nzchar(pn)) pn else feature_name
      } else {
        target_label <- "M (target)"
        parent_name <- feature_name
      }
    }
  } else {
    row_target <- ft[ft_idx, , drop = FALSE]
    target_iso <- FALSE
    parent_name <- feature_name
    target_label <- "M (target)"
  }

  traces <- list(
    .nts_eic_trace_from_row(
      row_target, mz_tol, rt_tol,
      label = target_label,
      use_bounds = use_bounds,
      is_target = TRUE
    )
  )

  if (!isTRUE(show_isotopes) || !has_ann) {
    return(traces)
  }

  iso_cols_present <- all(c("is_isotope", "parent_name") %in% names(ann))
  if (!iso_cols_present) {
    return(traces)
  }

  if (target_iso && isTRUE(include_parent_if_child) &&
        parent_name != feature_name) {
    pidx <- match(parent_name, as.character(ann$name))
    if (!is.na(pidx)) {
      traces[[length(traces) + 1L]] <- .nts_eic_trace_from_row(
        ann[pidx, , drop = FALSE], mz_tol, rt_tol,
        label = "M (parent)",
        use_bounds = use_bounds,
        is_target = FALSE
      )
    }
  }

  is_iso <- as.logical(ann[["is_isotope"]]) %in% TRUE
  same_parent <- as.character(ann[["parent_name"]]) == parent_name
  not_self <- as.character(ann[["name"]]) != feature_name
  child_idx <- which(is_iso & same_parent & not_self &
                       !is.na(ann[["parent_name"]]))
  if (length(child_idx)) {
    shifts <- suppressWarnings(as.numeric(ann[["isotope_shift_da"]][child_idx]))
    ord <- order(shifts, na.last = TRUE)
    child_idx <- child_idx[ord]
    for (k in child_idx) {
      lab <- as.character(ann[["isotope_label"]][k])
      if (is.na(lab) || !nzchar(lab)) {
        sh <- suppressWarnings(as.numeric(ann[["isotope_shift_da"]][k]))
        kind <- as.character(ann[["isotope_kind"]][k])
        lab <- if (!is.na(sh) && is.finite(sh)) {
          sprintf("\u0394%.4f Da%s", sh,
                  if (!is.na(kind) && nzchar(kind)) paste0(" (", kind, ")") else "")
        } else {
          paste0("isotope:", as.character(ann[["name"]][k]))
        }
      }
      traces[[length(traces) + 1L]] <- .nts_eic_trace_from_row(
        ann[k, , drop = FALSE], mz_tol, rt_tol,
        label = lab,
        use_bounds = use_bounds,
        is_target = FALSE
      )
    }
  }
  traces
}

#' Centered moving-average / Savitzky-Golay smoother for an EIC trace
#'
#' No new package dependency: moving-average is implemented with `cumsum`
#' (linear time); Savitzky-Golay uses the optional `signal` package when
#' available, otherwise falls back to moving-average. Reflect-padding at
#' the boundaries keeps the smoothed series the same length as `y`.
#'
#' @noRd
.nts_smooth_intensity <- function(y, method = "moving", n = 5L) {
  method <- match.arg(method, c("none", "moving", "savgol"))
  if (identical(method, "none")) return(as.numeric(y))
  y <- as.numeric(y)
  if (length(y) < 3L) return(y)
  n <- as.integer(n)[1L]
  if (!is.finite(n) || n < 3L) n <- 3L
  if (n %% 2L == 0L) n <- n + 1L
  if (length(y) < n) return(y)
  if (identical(method, "savgol") && requireNamespace("signal", quietly = TRUE)) {
    out <- tryCatch(
      suppressWarnings(signal::sgolayfilt(y, p = 3L, n = n)),
      error = function(e) NULL
    )
    if (!is.null(out)) {
      out[!is.finite(out)] <- y[!is.finite(out)]
      return(pmax(as.numeric(out), 0))
    }
  }
  k <- (n - 1L) %/% 2L
  pad <- c(rev(y[seq_len(k) + 1L]), y, rev(y[length(y) - seq_len(k)]))
  cs <- cumsum(c(0, pad))
  out <- (cs[(n + 1L):length(cs)] - cs[seq_len(length(cs) - n)]) / n
  out <- out[seq_along(y)]
  pmax(out, 0)
}

#' Extract feature or ROI chromatograms for EIC plotting.
#'
#' @return `list(cm = MChromatograms or XChromatograms or NULL, engine_used = character)`.
#' @noRd
.nts_eic_extract_chromatograms <- function(xd,
                                            traces,
                                            rt_tol = 0,
                                            engine = "auto") {
  engine <- match.arg(engine, c("auto", "feature", "chromatogram"))
  if (!length(traces)) {
    return(list(cm = NULL, engine_used = "none"))
  }
  ft_names <- vapply(traces, `[[`, character(1L), "name")

  use_feature <- engine %in% c("auto", "feature")
  if (use_feature) {
    fd <- tryCatch(xcms::featureDefinitions(xd), error = function(e) NULL)
    fd_names <- if (!is.null(fd)) BiocGenerics::rownames(fd) else character(0)
    if (!length(fd_names) || !all(ft_names %in% fd_names)) {
      if (identical(engine, "feature")) {
        stop(
          "Some traces are not features in featureDefinitions(xdata).",
          call. = FALSE
        )
      }
      use_feature <- FALSE
    }
  }

  bpp <- BiocParallel::SerialParam()
  cm <- NULL
  engine_used <- "chromatogram"
  if (use_feature) {
    feat_idx <- match(ft_names, BiocGenerics::rownames(
      xcms::featureDefinitions(xd)
    ))
    cm <- tryCatch(
      xcms::featureChromatograms(
        xd,
        features = feat_idx,
        expandRt = as.numeric(rt_tol)[1L],
        filled = TRUE,
        BPPARAM = bpp
      ),
      error = function(e) {
        message(
          "  featureChromatograms failed (", conditionMessage(e),
          "); falling back to chromatogram()."
        )
        NULL
      }
    )
    if (!is.null(cm)) engine_used <- "featureChromatograms"
  }
  if (is.null(cm)) {
    mz_mat <- do.call(rbind, lapply(traces, function(t) c(t$mz_min, t$mz_max)))
    rt_lo <- min(vapply(traces, `[[`, numeric(1L), "rt_min"), na.rm = TRUE)
    rt_hi <- max(vapply(traces, `[[`, numeric(1L), "rt_max"), na.rm = TRUE)
    cm <- tryCatch(
      xcms::chromatogram(
        xd,
        mz = mz_mat,
        rt = c(rt_lo, rt_hi),
        BPPARAM = bpp
      ),
      error = function(e) NULL
    )
    if (!is.null(cm)) engine_used <- "chromatogram"
  }
  list(cm = cm, engine_used = engine_used)
}

#' Per-sample colours for EIC: named by basename(file) or length `n`.
#' @noRd
.nts_resolve_eic_sample_colors <- function(sample_col_map, file_names, n_col) {
  bn <- basename(as.character(file_names)[seq_len(n_col)])
  if (!is.null(sample_col_map) && length(sample_col_map)) {
    nms <- names(sample_col_map)
    if (!is.null(nms) && any(nzchar(nms))) {
      out <- character(n_col)
      for (j in seq_len(n_col)) {
        hit <- match(bn[j], nms)
        out[j] <- if (!is.na(hit)) as.character(sample_col_map[hit])[1L] else NA_character_
      }
      miss <- is.na(out) | !nzchar(out)
      if (any(miss)) {
        nmiss <- sum(miss)
        fill <- tryCatch(
          grDevices::hcl.colors(nmiss, palette = "Dark 2"),
          error = function(e) grDevices::rainbow(nmiss, alpha = 0.85)
        )
        out[miss] <- fill
      }
      return(out)
    }
    if (length(sample_col_map) >= n_col) {
      return(as.character(sample_col_map)[seq_len(n_col)])
    }
  }
  if (n_col <= 1L) {
    return("#1F77B4")
  }
  tryCatch(
    grDevices::hcl.colors(n_col, palette = "Dark 3"),
    error = function(e) grDevices::rainbow(n_col, alpha = 0.85)
  )
}

#' Pull EICs for a batch of trace specs in **one** xcms call.
#'
#' Picks the engine in this order:
#' 1. [xcms::featureChromatograms] when every trace `name` is present in
#'    `featureDefinitions(xdata)` (this is by far the fastest path: xcms
#'    aggregates per-feature chromPeak m/z windows + `expandRt` so we read
#'    only the minimal scan slice).
#' 2. Otherwise [xcms::chromatogram] with an m/z matrix and the union RT
#'    window across all traces (still one call, no per-trace IO).
#'
#' [BiocParallel::SerialParam] is forced for the call to avoid worker
#' fork overhead that dominates Windows `SnowParam` chromatogram timings.
#'
#' @noRd
.nts_eic_traces_to_df <- function(xd,
                                  traces,
                                  src_files = NULL,
                                  rt_tol = 0,
                                  engine = "auto",
                                  smooth = "none",
                                  smooth_window = 5L) {
  ext <- .nts_eic_extract_chromatograms(xd, traces, rt_tol = rt_tol, engine = engine)
  cm <- ext$cm
  engine_used <- ext$engine_used
  if (is.null(cm)) return(NULL)

  n_trace <- nrow(cm)
  n_samp  <- ncol(cm)
  out <- vector("list", n_trace * n_samp)
  pos <- 1L
  for (i in seq_len(n_trace)) {
    tr <- traces[[i]]
    for (j in seq_len(n_samp)) {
      chr <- cm[i, j]
      rt_v <- MSnbase::rtime(chr)
      int_v <- MSnbase::intensity(chr)
      ok <- is.finite(rt_v) & is.finite(int_v)
      if (!any(ok)) {
        pos <- pos + 1L
        next
      }
      rt_v <- rt_v[ok]
      int_v <- int_v[ok]
      int_sm <- .nts_smooth_intensity(int_v, method = smooth, n = smooth_window)
      sn <- if (!is.null(src_files) && length(src_files) >= j) {
        basename(src_files[j])
      } else {
        paste0("sample_", j)
      }
      sn <- sub("\\.[Mm][Zz][MmXx][Ll]$", "", sn)
      out[[pos]] <- data.frame(
        rt = rt_v,
        intensity = int_sm,
        intensity_raw = int_v,
        sample = sn,
        series = sprintf("%s  m/z %.4f", tr$label, tr$mz),
        series_order = i,
        is_target = tr$is_target,
        peak_rt_min = tr$peak_rt_min,
        peak_rt_max = tr$peak_rt_max,
        feature = tr$name,
        stringsAsFactors = FALSE
      )
      pos <- pos + 1L
    }
  }
  out <- out[!vapply(out, is.null, logical(1L))]
  if (!length(out)) return(NULL)
  df <- do.call(rbind, out)
  attr(df, "engine") <- engine_used
  df
}

#' Plot Extracted Ion Chromatogram (EIC) for a feature with its isotope peaks
#'
#' For the requested `feature_name` (and, by default, every isotope peak
#' flagged by [flag_isotopes] that points back to the same monoisotopic
#' parent), pulls an EIC from `project@xdata` (the `XCMSnExp` produced by
#' [extract_features] after `findChromPeaks -> adjustRtime ->
#' groupChromPeaks -> fillChromPeaks`, so traces reflect the aligned
#' retention-time axis and gap-filled peaks). When every requested feature
#' is present in `featureDefinitions(xdata)`, [xcms::featureChromatograms]
#' is used (one call across all traces and samples, m/z window taken from
#' the chromPeaks of that feature, padded by `rt_tol` on the RT axis);
#' otherwise the function falls back to [xcms::chromatogram] with an m/z
#' matrix and the union RT window. The xcms call is forced onto
#' [BiocParallel::SerialParam] to avoid the multi-worker fork overhead
#' that dominates `SnowParam` chromatogram timings on Windows.
#'
#' If the requested feature itself is flagged as an isotope child
#' (`is_isotope == TRUE`), the monoisotopic parent is also drawn when
#' `include_parent_if_child = TRUE` (default), making it convenient to QC
#' a single isotope satellite.
#'
#' Samples are overlaid on the same panel and coloured by `sample`;
#' the monoisotopic / isotope traces are distinguished by **`linetype`**
#' (target = solid, isotopes = dashed). When `show_integration = TRUE`
#' (default), the per-trace integration range (`rtmin/rtmax` from
#' [extract_features], when available) is filled as a translucent ribbon
#' under each curve so the integration window is visually obvious.
#'
#' With `graphics = "xcms"`, draws the native **xcms** `plot()` method for
#' `XChromatograms` (base graphics) so **chromatographic peak integration**
#' is shown like the xcms vignette:
#' `peakBg = colours[chromPeaks(object)[, "sample"]]`.
#' This requires an `XChromatograms` result from [xcms::featureChromatograms]
#' (all isotope traces must exist as grouped features in `xdata`). If the
#' extractor fell back to [xcms::chromatogram] (`MChromatograms` without
#' attached `chromPeaks`), switch to `graphics = "ggplot"` or ensure isotope
#' `name`s are present in `featureDefinitions(xdata)`.
#'
#' @param project An [NTSProject] with non-empty `@feature_table` and (for
#'   isotope overlay) `@ms1_annotations` populated by [flag_isotopes].
#' @param feature_name Character scalar matching a `name` in
#'   `@feature_table` / `@ms1_annotations`.
#' @param mz_tol Symmetric m/z window in **Da** (default `0.01`), used when
#'   `mzmin/mzmax` are missing or `use_bounds = FALSE`.
#' @param rt_tol Symmetric RT padding in **seconds** (default `15`)
#'   applied around `rtmin/rtmax` (or the feature centroid). Smaller =
#'   faster.
#' @param show_isotopes If `TRUE` (default), overlay every isotope child
#'   pointing to the same monoisotopic parent.
#' @param include_parent_if_child If `TRUE` (default) and `feature_name` is
#'   itself an isotope child, also draw the monoisotopic parent trace.
#' @param use_bounds If `TRUE` (default), prefer the CentWave `mzmin/mzmax`
#'   and `rtmin/rtmax` on `@ms1_annotations`; set `FALSE` to force symmetric
#'   `mz_tol` / `rt_tol` windows around the feature centroid.
#' @param show_integration If `TRUE` (default), shade the integration
#'   window (per-trace `rtmin/rtmax`) under each curve with a translucent
#'   ribbon coloured by sample (**ggplot** only).
#' @param graphics `"ggplot"` (default) or `"xcms"`. The latter calls
#'   `plot(XChromatograms, col = ..., peakBg = ...)` and returns the object
#'   invisibly.
#' @param sample_col_map Optional **named** character vector mapping
#'   `basename(fileNames(xdata))` to colours (one name per file). If `NULL`,
#'   colours are chosen automatically (`ggplot` and `xcms`).
#' @param peak_type Passed as `peakType` to `plot(object)` when
#'   `graphics = "xcms"` (see **xcms** `plot,XChromatograms` in \code{?plot}).
#' @param engine xcms extraction engine: `"auto"` (default) prefers
#'   [xcms::featureChromatograms] when all traces are in
#'   `featureDefinitions(xdata)`, falling back to [xcms::chromatogram];
#'   `"feature"` forces the former and errors if any trace is not a
#'   feature; `"chromatogram"` forces the latter.
#' @param smooth EIC smoothing (**ggplot only**): `"moving"` (default, centered
#'   moving average), `"savgol"` (Savitzky-Golay; uses optional `signal`
#'   package; falls back to moving average), or `"none"`.
#' @param smooth_window Smoothing window (odd integer; default `5`).
#' @param files Optional; only used for fallback / label derivation when
#'   `@xdata` is empty.
#' @param ... Passed to `plot(x, ...)` when `graphics = "xcms"`.
#'
#' @return A **ggplot2** object when `graphics = "ggplot"`; **invisibly** the
#'   `XChromatograms` when `graphics = "xcms"`.
#' @export
plot_feature_eic <- function(project,
                             feature_name,
                             mz_tol = 0.01,
                             rt_tol = 15,
                             show_isotopes = TRUE,
                             include_parent_if_child = TRUE,
                             use_bounds = TRUE,
                             show_integration = TRUE,
                             graphics = c("ggplot", "xcms"),
                             sample_col_map = NULL,
                             peak_type = c("polygon", "point", "rectangle", "none"),
                             engine = c("auto", "feature", "chromatogram"),
                             smooth = c("moving", "savgol", "none"),
                             smooth_window = 5L,
                             files = NULL,
                             ...) {
  stopifnot(inherits(project, "NTSProject"))
  graphics <- match.arg(graphics)
  peak_type <- match.arg(peak_type)
  engine <- match.arg(engine)
  smooth <- match.arg(smooth)

  traces <- .nts_eic_collect_traces(
    project,
    feature_name = feature_name,
    mz_tol = mz_tol,
    rt_tol = rt_tol,
    show_isotopes = show_isotopes,
    include_parent_if_child = include_parent_if_child,
    use_bounds = use_bounds
  )
  if (!length(traces)) {
    stop("No EIC trace could be assembled for ", feature_name, call. = FALSE)
  }
  target <- traces[[1L]]
  feat_rt <- target$rt
  feat_mz <- target$mz

  xd <- .nts_get_xdata(project, files)

  t0 <- proc.time()[[3L]]
  message(
    "\u25b6 Extracting EIC: ", feature_name,
    " (m/z=", round(feat_mz, 4), ", RT=", round(feat_rt, 1),
    ") + ", length(traces) - 1L, " isotope trace(s); engine='", engine,
    "', graphics='", graphics,
    if (identical(graphics, "ggplot")) paste0("', smooth='", smooth, "'") else "'",
    "..."
  )

  if (identical(graphics, "xcms")) {
    ext <- .nts_eic_extract_chromatograms(
      xd, traces, rt_tol = rt_tol, engine = engine
    )
    cm <- ext$cm
    eng_used <- ext$engine_used
    if (is.null(cm)) {
      stop("EIC extraction returned no data for ", feature_name, call. = FALSE)
    }
    if (!methods::is(cm, "XChromatograms")) {
      stop(
        "graphics = \"xcms\" needs XChromatograms from featureChromatograms(); ",
        "got ", class(cm)[[1L]], " (", eng_used, "). ",
        "Isotope traces may be missing from featureDefinitions(xdata), or ",
        "use graphics = \"ggplot\", or engine = \"feature\" after ensuring ",
        "all traces are grouped features.",
        call. = FALSE
      )
    }
    fn <- MSnbase::fileNames(xd)
    n_col <- ncol(cm)
    n_tr <- nrow(cm)
    title <- if (n_tr > 1L) {
      sprintf(
        "EIC (xcms): %s + %d trace(s) | m/z %.4f RT %.1f",
        feature_name, n_tr - 1L, feat_mz, feat_rt
      )
    } else {
      sprintf("EIC (xcms): %s | m/z %.4f RT %.1f", feature_name, feat_mz, feat_rt)
    }
    cols <- .nts_resolve_eic_sample_colors(sample_col_map, fn, n_col)
    pks <- tryCatch(xcms::chromPeaks(cm), error = function(e) NULL)
    plot_args <- list(
      col = cols,
      peakType = peak_type,
      main = title,
      xlab = "Retention time (s)",
      ylab = "Intensity"
    )
    if (is.matrix(pks) && nrow(pks) && "sample" %in% colnames(pks)) {
      si <- suppressWarnings(as.integer(pks[, "sample"]))
      peak_bg <- cols[pmax(1L, pmin(n_col, si))]
      bad <- !is.finite(si) | si < 1L | si > n_col
      if (any(bad)) {
        peak_bg[bad] <- grDevices::adjustcolor("#888888", alpha.f = 0.35)
      }
      plot_args$peakBg <- peak_bg
    }
    udots <- list(...)
    if (length(udots)) {
      for (nm in names(udots)) plot_args[[nm]] <- udots[[nm]]
    }
    do.call(plot, c(list(cm), plot_args))
    message(sprintf(
      "  done in %.2fs via %s (xcms plot, %d traces \u00d7 %d samples)",
      proc.time()[[3L]] - t0,
      eng_used,
      n_tr,
      n_col
    ))
    return(invisible(cm))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_feature_eic() when graphics = \"ggplot\"",
         call. = FALSE)
  }

  src_files <- tryCatch(
    MSnbase::fileNames(xd),
    error = function(e) files
  )

  df <- .nts_eic_traces_to_df(
    xd, traces,
    src_files = src_files,
    rt_tol = rt_tol,
    engine = engine,
    smooth = smooth,
    smooth_window = smooth_window
  )
  if (is.null(df) || !nrow(df)) {
    stop("EIC extraction returned no data for ", feature_name, call. = FALSE)
  }
  eng_used <- attr(df, "engine")
  if (is.null(eng_used)) eng_used <- "chromatogram"
  message(sprintf(
    "  done in %.2fs via %s (%d points, %d traces, %d samples)",
    proc.time()[[3L]] - t0,
    eng_used,
    nrow(df),
    length(unique(df$series)),
    length(unique(df$sample))
  ))

  df$series <- factor(
    df$series,
    levels = unique(df$series[order(df$series_order)])
  )
  series_levels <- levels(df$series)
  lt_vec <- rep("dashed", length(series_levels))
  lt_vec[1L] <- "solid"
  names(lt_vec) <- series_levels

  n_samples <- length(unique(df$sample))
  n_traces  <- nlevels(df$series)
  title <- if (n_traces > 1L) {
    sprintf(
      "EIC: %s + %d isotope trace(s)  (m/z %.4f, RT %.1f)",
      feature_name, n_traces - 1L, feat_mz, feat_rt
    )
  } else {
    sprintf("EIC: %s  (m/z %.4f, RT %.1f)", feature_name, feat_mz, feat_rt)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["rt"]],
    y = .data[["intensity"]],
    colour = .data[["sample"]],
    linetype = .data[["series"]],
    group = interaction(.data[["sample"]], .data[["series"]])
  ))

  if (isTRUE(show_integration)) {
    in_peak <- is.finite(df$peak_rt_min) & is.finite(df$peak_rt_max) &
      df$rt >= df$peak_rt_min & df$rt <= df$peak_rt_max
    df_fill <- df[in_peak, , drop = FALSE]
    if (nrow(df_fill)) {
      p <- p + ggplot2::geom_ribbon(
        data = df_fill,
        ggplot2::aes(
          x = .data[["rt"]],
          ymin = 0,
          ymax = .data[["intensity"]],
          fill = .data[["sample"]],
          group = interaction(.data[["sample"]], .data[["series"]])
        ),
        colour = NA,
        alpha = 0.15,
        inherit.aes = FALSE
      )
    }
  }

  p <- p +
    ggplot2::geom_line(linewidth = 0.6, alpha = 0.95) +
    ggplot2::geom_vline(
      xintercept = feat_rt,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.3
    ) +
    ggplot2::scale_linetype_manual(values = lt_vec, drop = FALSE) +
    ggplot2::labs(
      title = title,
      x = "Retention Time (s)",
      y = "Intensity",
      colour = "Sample",
      fill = "Sample",
      linetype = "Trace"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      legend.position = "right"
    )
  if (n_samples > 12) {
    p <- p + ggplot2::guides(colour = "none", fill = "none")
  }
  p
}

#' Annotation row lookup for a feature
#' @noRd
.nts_get_annotation_row <- function(project, feature_name) {
  ann <- methods::slot(project, "db_annotations")
  if (!nrow(ann) || !"name" %in% names(ann)) {
    stop("db_annotations is empty or lacks 'name' column", call. = FALSE)
  }
  sub <- ann[as.character(ann$name) == feature_name, , drop = FALSE]
  if (!nrow(sub)) {
    stop("feature_name '", feature_name, "' not found in db_annotations", call. = FALSE)
  }
  if ("score_total" %in% names(sub)) {
    sub <- sub[order(sub$score_total, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  }
  sub[1L, , drop = FALSE]
}

#' Normalize a peak matrix for mirror plotting
#' @noRd
.nts_scale_mirror_matrix <- function(mat, sign = 1) {
  if (is.null(mat) || !is.matrix(mat) || !nrow(mat)) {
    return(data.frame(mz = numeric(0), intensity = numeric(0), source = character(0)))
  }
  mz <- as.numeric(mat[, 1L])
  intensity <- as.numeric(mat[, 2L])
  ok <- is.finite(mz) & is.finite(intensity) & intensity > 0
  if (!any(ok)) {
    return(data.frame(mz = numeric(0), intensity = numeric(0), source = character(0)))
  }
  mz <- mz[ok]
  intensity <- intensity[ok]
  intensity <- 100 * intensity / max(intensity, na.rm = TRUE)
  data.frame(
    mz = mz,
    intensity = sign * intensity,
    stringsAsFactors = FALSE
  )
}

#' Get the best experimental MS2 spectrum for a feature
#'
#' Selects the exact experimental spectrum the scoring engine used. Preference
#' order: (1) an explicit `best_ms2_index` written by
#' [annotate_by_spectral_match] (0.1.30+); (2) a label match against
#' `exp_label` - the `db_annotations$exp_CE` value, e.g.
#' `"CE10#1835@90ppb2"` - which is what older projects carry; (3) subset #1.
#' Matching by label recovers the winner for projects mined before
#' `best_ms2_index` existed, where a feature routinely offered several MS2
#' scans and the highest-scoring one was often not the first.
#' @noRd
.nts_get_best_exp_ms2 <- function(project, feature_name, best_ms2_index = NULL,
                                  exp_label = NULL) {
  # Build matrices through the same helper the scoring engine uses, so the
  # subset order (and its labels) is identical to the one that produced the
  # recorded `best_ms2_index` / `exp_CE`.
  pl <- .nts_feature_ms2_peak_list(project, feature_name, "ms2data",
                                   with_labels = TRUE)
  subset <- pl$mats
  labs <- pl$labels
  if (!length(subset)) {
    # Fall back to a direct read for the (rare) case the helper cannot serve.
    sp <- methods::slot(project, "ms2data")[[feature_name]]
    if (is.null(sp)) {
      stop("No MS2 spectra found for feature: ", feature_name, call. = FALSE)
    }
    subset <- .nts_exp_ms2_to_peak_matrices(sp)
    labs <- names(subset)
  }
  if (!length(subset)) {
    stop("No valid MS2 subset found for feature: ", feature_name, call. = FALSE)
  }

  idx <- NA_integer_
  # (2) label match first when an explicit index is unavailable - the label is
  # what older projects recorded, and it is unambiguous across source files.
  if (is.null(best_ms2_index) || !length(best_ms2_index) ||
      !is.finite(suppressWarnings(as.integer(best_ms2_index)[1L]))) {
    lab <- tryCatch(as.character(exp_label)[1L], error = function(e) NA_character_)
    if (!is.na(lab) && nzchar(lab) && !is.null(labs)) {
      hit <- which(as.character(labs) == lab)
      if (length(hit)) idx <- hit[[1L]]
    }
  }
  # (1) explicit index.
  if (!is.finite(idx)) {
    idx <- suppressWarnings(as.integer(best_ms2_index)[1L])
  }
  if (!is.finite(idx) || idx < 1L || idx > length(subset)) {
    idx <- 1L
  }
  subset[[idx]]
}

#' First present column among a list of candidate names
#' @noRd
.nts_pick_col <- function(df, candidates) {
  nm <- names(df)
  if (is.null(nm)) return(NULL)
  hit <- candidates[candidates %in% nm]
  if (length(hit)) hit[[1L]] else NULL
}

#' Reference spectra (per NCE) from a flat [NTSDatabase] library
#'
#' `NTSDatabase` - the flat, row-oriented library NTSminer itself builds - keeps
#' every MS2 peak in `@spectra_long` (`Lab_ID`, `Polarity`, `NCE`, `frag_mz`,
#' `intensity`). This reshapes one compound's rows back into the nested
#' `list(NCE -> data.frame(mz, intensity))` shape that
#' [.nts_lib_consensus_matrix] consumes, so mirror plotting works with the flat
#' libraries as well as with legacy [databaseClass] ones. Column names are
#' resolved leniently (`Lab_ID`/`lab_id`, `frag_mz`/`mz`, `Polarity`/`polarity`,
#' `NCE`/`nce`/`CE`), and a missing `NCE` column collapses to a single layer.
#'
#' @param database An [NTSDatabase].
#' @param ref_id `Lab_ID` of the reference compound.
#' @param ion_mode One of `"negative"` / `"positive"`.
#' @param nce Optional single NCE level to keep (e.g. `"CE10"`). When supplied,
#'   only rows whose `NCE` equals it are returned, so the mirror plot draws the
#'   one collision-energy level that was actually scored rather than a
#'   cross-CE consensus. `NULL` (default) keeps every level.
#' @return Named `list` of `data.frame(mz, intensity)` per NCE level, or `NULL`
#'   when the compound carries no spectrum of that polarity (or none at `nce`).
#' @noRd
.nts_ref_levels_from_flat <- function(database, ref_id, ion_mode, nce = NULL) {
  sl <- as.data.frame(methods::slot(database, "spectra_long"))
  if (!nrow(sl)) return(NULL)
  id_col  <- .nts_pick_col(sl, c("Lab_ID", "lab_id", "Lab.ID", "lab.id", "id"))
  mz_col  <- .nts_pick_col(sl, c("frag_mz", "mz", "mz_value", "Mz", "MZ"))
  int_col <- .nts_pick_col(sl, c("intensity", "Intensity", "int", "Int"))
  if (is.null(id_col) || is.null(mz_col) || is.null(int_col)) return(NULL)
  keep <- as.character(sl[[id_col]]) == as.character(ref_id)
  pol_col <- .nts_pick_col(sl, c("Polarity", "polarity", "Ion_mode", "ion_mode", "Ion_Mode"))
  if (!is.null(pol_col)) {
    want <- if (identical(ion_mode, "positive")) "POS" else "NEG"
    kk <- toupper(substr(as.character(sl[[pol_col]]), 1L, 3L))
    keep <- keep & !is.na(kk) & (kk == want)
  }
  nce_col <- .nts_pick_col(sl, c("NCE", "nce", "CE", "collision_energy"))
  # Restrict to one collision-energy level when the caller names it. Done on the
  # full frame (before the `sub` slice) so the NCE column is still visible even
  # if a lenient name was resolved above.
  if (!is.null(nce) && length(nce) && !is.na(nce[1L]) &&
      nzchar(as.character(nce[1L])) && !is.null(nce_col)) {
    keep <- keep & as.character(sl[[nce_col]]) == as.character(nce[1L])
  }
  sub <- sl[keep, , drop = FALSE]
  if (!nrow(sub)) return(NULL)
  if (!is.null(nce_col) && nce_col %in% names(sub)) {
    grp <- as.character(sub[[nce_col]])
  } else {
    grp <- rep("Merged", nrow(sub))
  }
  grp[is.na(grp) | !nzchar(grp)] <- "Merged"
  lv <- split(sub, grp)
  lapply(lv, function(d) {
    data.frame(mz = as.numeric(d[[mz_col]]), intensity = as.numeric(d[[int_col]]))
  })
}

#' Resolve reference MS2 matrix for a given annotation
#'
#' Works with both reference-library formats NTSminer handles: the flat
#' [NTSDatabase] (peaks in `@spectra_long`) that all shipped `lc_libraries`
#' files use, and the legacy [databaseClass] (nested `@spectra.data`).
#' `ion_mode = "auto"` tries negative first, then positive.
#'
#' `lib_ce` narrows the reference to a single collision-energy level - the
#' `db_annotations$lib_CE` value the scoring engine picked - so the mirror plot's
#' reference half is the same spectrum the score came from, not a consensus of
#' every CE. It falls back to the full (all-CE) consensus whenever that level is
#' missing, keeping older projects and single-level libraries working.
#'
#' @noRd
.nts_get_ref_ms2 <- function(database, ref_id, ion_mode = c("auto", "negative", "positive"),
                             lib_ce = NULL) {
  ion_mode <- match.arg(ion_mode)
  modes <- if (identical(ion_mode, "auto")) c("negative", "positive") else ion_mode
  ce <- tryCatch(as.character(lib_ce)[1L], error = function(e) NA_character_)
  if (is.na(ce) || !nzchar(ce)) ce <- NULL
  is_classic <- methods::is(database, "databaseClass")
  is_flat <- !is_classic && methods::is(database, "NTSDatabase")
  if (!is_classic && !is_flat) {
    stop("`database` must be a databaseClass or NTSDatabase object", call. = FALSE)
  }

  ## ---- flat NTSDatabase: reshape @spectra_long ----
  if (is_flat) {
    for (im in modes) {
      # Prefer the exact CE level; fall back to the all-level consensus.
      mat <- if (!is.null(ce)) {
        .nts_lib_consensus_matrix(.nts_ref_levels_from_flat(database, ref_id, im, nce = ce))
      } else NULL
      if (is.null(mat) || !nrow(mat)) {
        mat <- .nts_lib_consensus_matrix(.nts_ref_levels_from_flat(database, ref_id, im))
      }
      if (!is.null(mat) && nrow(mat)) return(mat)
    }
    stop("Reference spectrum not found for ref_id: ", ref_id,
         " (NTSDatabase, ion_mode=", paste(modes, collapse = "/"), ")", call. = FALSE)
  }

  ## ---- legacy databaseClass: nested @spectra.data ----
  spec_data <- methods::slot(database, "spectra.data")
  slots <- paste0("Spectra.", modes)
  ref_levels <- NULL
  for (sn in slots) {
    root <- spec_data[[sn]]
    if (is.null(root)) next
    ref_levels <- root[[ref_id]]
    if (is.null(ref_levels)) {
      ref_levels <- root[[as.character(ref_id)]]
    }
    if (!is.null(ref_levels)) break
  }
  if (is.null(ref_levels)) {
    stop("Reference spectrum not found for ref_id: ", ref_id, call. = FALSE)
  }
  # Narrow to the scored CE when the nested list is keyed by CE.
  if (!is.null(ce) && is.list(ref_levels) && !is.null(names(ref_levels)) &&
      ce %in% names(ref_levels)) {
    ref_levels <- ref_levels[ce]
  }
  ref_mat <- .nts_lib_consensus_matrix(ref_levels)
  if (is.null(ref_mat) || !nrow(ref_mat)) {
    stop("Reference spectrum is empty for ref_id: ", ref_id, call. = FALSE)
  }
  ref_mat
}

#' Shared KMD plot filtering logic
#' @noRd
.nts_filter_kmd_plot_df <- function(df,
                                    filter_isf = TRUE,
                                    only_show_series = FALSE,
                                    anchor_annotated = FALSE,
                                    target_kmds = NULL,
                                    kmd_tol = 0.01) {
  if (isTRUE(filter_isf) && "is_isf" %in% names(df)) {
    df <- df[is.na(df$is_isf) | !df$is_isf, , drop = FALSE]
  }
  if (isTRUE(only_show_series) && "kmd_series_id" %in% names(df)) {
    df <- df[!is.na(df$kmd_series_id) & nzchar(as.character(df$kmd_series_id)), , drop = FALSE]
  }
  if (isTRUE(anchor_annotated)) {
    valid_series <- unique(
      df$kmd_series_id[
        !is.na(df$kmd_series_id) &
          !is.na(df$Compound.name) &
          is.finite(df$score_total) &
          df$score_total > 0
      ]
    )
    valid_series <- valid_series[!is.na(valid_series) & nzchar(as.character(valid_series))]
    df <- df[df$kmd_series_id %in% valid_series, , drop = FALSE]
  }
  if (!is.null(target_kmds)) {
    target_kmds <- as.numeric(target_kmds)
    target_kmds <- target_kmds[is.finite(target_kmds)]
    if (length(target_kmds)) {
      kmd_tol <- as.numeric(kmd_tol)[1L]
      if (!is.finite(kmd_tol) || kmd_tol < 0) {
        stop("kmd_tol must be a non-negative finite number", call. = FALSE)
      }
      kmd_dist <- abs(outer(df$KMD, target_kmds, "-"))
      keep <- rowSums(kmd_dist <= kmd_tol, na.rm = TRUE) > 0
      df <- df[keep, , drop = FALSE]
    }
  }
  df
}

#' Drop mirror peaks below a relative-intensity floor
#'
#' Intensities are already normalized to percent of the base peak, so
#' `min_rel_intensity` compares to them directly. `0` (or a non-finite value)
#' keeps every peak.
#' @noRd
.nts_drop_low_rel_peaks <- function(df, min_rel_intensity) {
  if (!nrow(df)) return(df)
  thr <- suppressWarnings(as.numeric(min_rel_intensity)[1L])
  if (!is.finite(thr) || thr <= 0) return(df)
  df[is.finite(df$intensity) & abs(df$intensity) >= thr, , drop = FALSE]
}

#' Mark mirror-plot peaks that have a counterpart on the other half
#'
#' A "matched fragment ion" on a mirror plot is an experimental peak whose m/z is
#' within `tol_da` of some reference peak (and vice versa). Matching is done
#' directly on the peaks being drawn - after the noise floor - so the labels can
#' never disagree with the figure. Each returned data frame gains a `matched_mz`
#' column (the peak's own m/z when matched, `NA` otherwise).
#' @noRd
.nts_mark_matched_peaks <- function(exp_df, ref_df, tol_da) {
  exp_df$matched_mz <- rep(NA_real_, nrow(exp_df))
  ref_df$matched_mz <- rep(NA_real_, nrow(ref_df))
  tol_da <- suppressWarnings(as.numeric(tol_da)[1L])
  if (!is.finite(tol_da) || tol_da < 0) tol_da <- 0.02
  if (!nrow(exp_df) || !nrow(ref_df)) {
    return(list(exp = exp_df, ref = ref_df))
  }
  em <- exp_df$mz
  rm <- ref_df$mz
  ok_e <- is.finite(em)
  ok_r <- is.finite(rm)
  # Nearest partner on the opposite half, greedily and independently per peak.
  # Two peaks that round to the same reported m/z are rare within one spectrum,
  # so a one-sided nearest match is enough and keeps exp/lib m/z both visible.
  de <- rep(Inf, length(em))
  dr <- rep(Inf, length(rm))
  if (any(ok_e) && any(ok_r)) {
    de[ok_e] <- abs(rm[vapply(em[ok_e], function(x) which.min(abs(rm - x))[1L], integer(1L))] - em[ok_e])
    dr[ok_r] <- abs(em[vapply(rm[ok_r], function(x) which.min(abs(em - x))[1L], integer(1L))] - rm[ok_r])
  }
  exp_df$matched_mz[de <= tol_da] <- em[de <= tol_da]
  ref_df$matched_mz[dr <= tol_da] <- rm[dr <= tol_da]
  list(exp = exp_df, ref = ref_df)
}

#' Spread rotated m/z labels that would otherwise overlap
#'
#' On a mirror plot two matched fragments 1 Da apart sit at almost the same x, so
#' their rotated m/z numbers collide. Peaks closer than `sep` (m/z units) are
#' pushed to alternating vertical levels `gap` (intensity units) apart. Adds a
#' `y_lab` column - the anchor for the rotated text, `dir * level * gap` outward
#' from the tip - and a `lvl` column so the caller can reserve headroom.
#' @noRd
.nts_stagger_labels <- function(lab, dir = 1, sep, gap = 34) {
  lab$y_lab <- lab$intensity
  lab$lvl   <- 0L
  if (nrow(lab) < 2L) return(lab)
  lab <- lab[order(lab$mz), , drop = FALSE]
  # Last m/z placed at each of the two levels, to test a real overlap (not just
  # against the immediately preceding peak, which would mis-level a run of 3+).
  last <- c(NA_real_, NA_real_)
  for (i in seq_len(nrow(lab))) {
    m  <- lab$mz[i]
    lv <- if (is.finite(last[1L]) && m - last[1L] < sep) 1L else 0L
    if (lv == 1L && is.finite(last[2L]) && m - last[2L] < sep) lv <- 0L
    last[lv + 1L] <- m
    lab$lvl[i] <- lv
  }
  lab$y_lab <- lab$intensity + dir * lab$lvl * gap
  lab
}

#' Compose score-decomposition subtitle lines for a mirror plot
#'
#' Reads the per-component scores stored in `db_annotations` and renders two
#' text lines: how `score_total` splits into the MS1 m/z, RT and MS2 terms, and
#' how the MS2 term itself splits into ndp / nl / md / diagnostic. Missing
#' columns render as `NA`, so the helper stays safe for older annotations.
#'
#' `score_total` is `w_mz*score_mz + w_rt*score_rt + w_ms2*score_ms2`. The flat
#' library path ([annotate_by_spectral_match]) stores raw components with
#' weights `0.25 / 0.25 / 0.50`, whereas the legacy [databaseClass] path stores
#' m/z and RT already weighted (`1 / 1 / 1`); the scheme that reproduces
#' `score_total` is chosen automatically.
#' @noRd
.nts_score_breakdown_lines <- function(ann_row) {
  comp <- function(col) {
    if (col %in% names(ann_row)) {
      suppressWarnings(as.numeric(ann_row[[col]][[1L]]))
    } else {
      NA_real_
    }
  }
  fmt <- function(x) if (is.finite(x)) sprintf("%.3f", x) else "NA"

  tot <- comp("score_total")
  mz  <- comp("score_mz")
  rt  <- comp("score_rt")
  ms2 <- comp("score_ms2")

  weights <- c(0.25, 0.25, 0.5)
  if (all(is.finite(c(tot, mz, rt, ms2)))) {
    candidates <- list(c(0.25, 0.25, 0.5), c(1, 1, 1))
    err <- vapply(candidates, function(w) abs(sum(w * c(mz, rt, ms2)) - tot), numeric(1))
    weights <- candidates[[which.min(err)]]
  }
  line1 <- sprintf(
    "score_total %s  =  m/z %s + RT %s + MS2 %s",
    fmt(tot), fmt(weights[1L] * mz), fmt(weights[2L] * rt), fmt(weights[3L] * ms2)
  )

  ndp <- comp("score_ndp")
  nl  <- comp("score_nl")
  md  <- comp("score_md")
  # The stored `score_diagnostic` column is ambiguous (the flat path may keep the
  # raw [0,1] coverage while the composite uses it weighted by 0.15), so the
  # diagnostic term is recovered as the residual that makes the four terms sum to
  # the MS2 composite exactly.
  dg <- if (all(is.finite(c(ms2, ndp, nl, md)))) max(0, ms2 - (ndp + nl + md)) else NA_real_
  extra <- character(0)
  ndp_raw <- comp("ndp_raw")
  npk <- comp("n_matched_peaks")
  # `score_diagnostic` is the *raw* diagnostic-fragment coverage kept by the flat
  # path; it is reported here for transparency because it does not itself enter the
  # MS2 composite (the composite's own diag term is the residual above).
  diag_cov <- comp("score_diagnostic")
  if (is.finite(ndp_raw)) extra <- c(extra, sprintf("raw ndp %s", fmt(ndp_raw)))
  if (is.finite(npk)) extra <- c(extra, sprintf("matched peaks %d", as.integer(npk)))
  if (is.finite(diag_cov)) extra <- c(extra, sprintf("diag coverage %s", fmt(diag_cov)))
  line2 <- sprintf(
    "MS2 %s  =  ndp %s + nl %s + md %s + diag %s%s",
    fmt(ms2), fmt(ndp), fmt(nl), fmt(md), fmt(dg),
    if (length(extra)) paste0("   [", paste(extra, collapse = ", "), "]") else ""
  )
  c(line1, line2)
}

#' Plot an experimental vs reference MS2 mirror spectrum
#'
#' Uses the top-ranked annotation for `feature_name` and draws the two halves of
#' the same m/z axis: experimental up, reference down. Both halves are the exact
#' spectra that produced the reported score, so the figure and the numbers agree:
#' the experimental half is the winning MS2 subset (the one named by
#' `@db_annotations$exp_CE`, e.g. `"CE10#1835@90ppb2"` - or, for projects mined
#' from 0.1.30 on, the explicit `best_ms2_index`), and the reference half is the
#' library collision-energy level `@db_annotations$lib_CE` rather than a
#' consensus of every CE on file. When those columns are absent (older projects)
#' or a named level cannot be found, it falls back to subset #1 and the all-CE
#' consensus, so nothing regresses.
#'
#' @param project An [NTSProject] with non-empty `@db_annotations`.
#' @param database Reference library: a flat [NTSDatabase] (peaks in
#'   `@spectra_long`, the format all `lc_libraries` files use) or a legacy
#'   [databaseClass] (nested `@spectra.data`).
#' @param feature_name Feature identifier present in `project@db_annotations$name`.
#' @param ion_mode Which library slot to use. `"auto"` tries negative first,
#'   then positive.
#' @param min_rel_intensity Noise floor as a **percentage of the base peak**
#'   (default `2.5`, matching the `int_thresh = 0.025` denoising used by the
#'   scoring engine, so the mirror shows the peaks that were actually compared).
#'   Peaks whose relative intensity is below this are dropped from both traces
#'   before plotting; pass `0` to keep every peak.
#' @param label_matched_mz Logical; if `TRUE` (default) annotate every *matched*
#'   fragment ion - an experimental peak within `label_mz_tol` of a reference
#'   peak, and vice versa - with its own m/z, printed beside the peak tip on both
#'   halves. Matching is done on the drawn peaks (after `min_rel_intensity`), so a
#'   label can never point at a peak that is not in the figure. Set `FALSE` for
#'   the plain unlabelled mirror.
#' @param label_mz_tol m/z tolerance (Da) for calling two peaks a matched pair
#'   (default `0.02`, the `ms2_tol_da` the scoring engine used).
#' @param label_digits Decimal places used when printing the labelled m/z
#'   (default `2`).
#'
#' @return A `ggplot2` object.
#' @export
plot_ms2_mirror <- function(project,
                            database,
                            feature_name,
                            ion_mode = c("auto", "negative", "positive"),
                            min_rel_intensity = 2.5,
                            label_matched_mz = TRUE,
                            label_mz_tol = 0.02,
                            label_digits = 2) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_ms2_mirror()", call. = FALSE)
  }
  stopifnot(inherits(project, "NTSProject"))
  ion_mode <- match.arg(ion_mode)
  min_rel_intensity <- suppressWarnings(as.numeric(min_rel_intensity)[1L])
  if (is.na(min_rel_intensity) || min_rel_intensity < 0) {
    stop("`min_rel_intensity` must be a non-negative finite percentage", call. = FALSE)
  }
  label_digits <- suppressWarnings(as.integer(label_digits)[1L])
  if (is.na(label_digits) || label_digits < 0) label_digits <- 2L
  do_label <- isTRUE(label_matched_mz)

  ann_row <- .nts_get_annotation_row(project, feature_name)
  ref_id <- if ("ref_id" %in% names(ann_row)) {
    as.character(ann_row$ref_id[[1L]])
  } else {
    stop("db_annotations lacks 'ref_id' column", call. = FALSE)
  }
  best_idx <- if ("best_ms2_index" %in% names(ann_row)) ann_row$best_ms2_index[[1L]] else NULL
  exp_lab  <- if ("exp_CE" %in% names(ann_row)) as.character(ann_row$exp_CE[[1L]]) else NULL
  lib_ce   <- if ("lib_CE" %in% names(ann_row)) as.character(ann_row$lib_CE[[1L]]) else NULL
  # Both halves must be the exact spectra that produced the reported score:
  # the experimental subset `exp_CE` names (fall back to `best_ms2_index`) and
  # the library CE level `lib_CE`. Otherwise a feature with several MS2 scans -
  # or a library stored at several collision energies - is drawn against the
  # wrong peaks and the mirror looks far worse than the match really is.
  exp_mat <- .nts_get_best_exp_ms2(project, feature_name, best_idx, exp_label = exp_lab)
  ref_mat <- .nts_get_ref_ms2(database, ref_id, ion_mode = ion_mode, lib_ce = lib_ce)

  exp_df <- .nts_scale_mirror_matrix(exp_mat, sign = 1)
  ref_df <- .nts_scale_mirror_matrix(ref_mat, sign = -1)
  exp_df <- .nts_drop_low_rel_peaks(exp_df, min_rel_intensity)
  ref_df <- .nts_drop_low_rel_peaks(ref_df, min_rel_intensity)
  # Mark the fragments that pair across the two halves so they (and only they)
  # can be labelled. Done on the noise-floored peaks actually drawn, so a label
  # can never name a peak that is missing from the figure.
  if (do_label) {
    mk <- .nts_mark_matched_peaks(exp_df, ref_df, label_mz_tol)
    exp_df <- mk$exp
    ref_df <- mk$ref
  } else {
    exp_df$matched_mz <- rep(NA_real_, nrow(exp_df))
    ref_df$matched_mz <- rep(NA_real_, nrow(ref_df))
  }
  exp_df$source <- "Experimental"
  ref_df$source <- "Reference"
  df <- rbind(exp_df, ref_df)
  if (!nrow(df)) {
    stop("No valid peaks available for mirror plot", call. = FALSE)
  }

  hit_name <- if ("Compound.name" %in% names(ann_row)) as.character(ann_row$Compound.name[[1L]]) else ref_id

  # Header line: compound + (when present) the experimental / library CE, since a
  # CE mismatch is the usual reason a good match looks visually poor. The full
  # experimental label is kept (e.g. "CE10#1835@90ppb2") because it names the one
  # scan now drawn on the top half, not just its collision energy.
  ce_bits <- character(0)
  if ("exp_CE" %in% names(ann_row)) {
    ce_bits <- c(ce_bits, paste0("exp ", as.character(ann_row$exp_CE[[1L]])))
  }
  if ("lib_CE" %in% names(ann_row)) {
    ce_bits <- c(ce_bits, paste0("lib ", as.character(ann_row$lib_CE[[1L]])))
  }
  header <- paste0(
    "Match: ", hit_name,
    if (length(ce_bits)) paste0("  |  ", paste(ce_bits, collapse = " vs ")) else "",
    "  |  noise<", format(min_rel_intensity, trim = TRUE), "%"
  )
  subtitle <- paste(c(header, .nts_score_breakdown_lines(ann_row)), collapse = "\n")

  # Matched peaks whose own m/z gets printed next to the tip. Labels on the
  # experimental half read upward from the tip; reference-half labels read
  # downward, so the two halves never collide.
  lab_exp <- exp_df[!is.na(exp_df$matched_mz), , drop = FALSE]
  lab_ref <- ref_df[!is.na(ref_df$matched_mz), , drop = FALSE]
  if (nrow(lab_exp)) lab_exp$lab <- formatC(lab_exp$matched_mz, format = "f", digits = label_digits)
  if (nrow(lab_ref)) lab_ref$lab <- formatC(lab_ref$matched_mz, format = "f", digits = label_digits)
  # Peaks closer than a glyph's worth of m/z get moved to a second vertical level
  # so their rotated numbers do not sit on top of each other. The m/z separation
  # is judged against the panel width, so it adapts to the feature's mass range.
  lbl_gap <- 34
  x_span <- diff(range(df$mz, na.rm = TRUE))
  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  sep_mz <- max(0.3, 3 * x_span / 150)
  if (do_label && nrow(lab_exp)) lab_exp <- .nts_stagger_labels(lab_exp, dir = 1, sep = sep_mz, gap = lbl_gap)
  if (do_label && nrow(lab_ref)) lab_ref <- .nts_stagger_labels(lab_ref, dir = -1, sep = sep_mz, gap = lbl_gap)
  # Headroom for the printed m/z: a rotated label is about as long as its text
  # width (~7 characters for "1000.00"), so the tip needs a good margin past 100,
  # plus one more level wherever the stagger pushed a label outward.
  lvl_max <- max(c(lab_exp$lvl, lab_ref$lvl, 0L))
  y_pad <- if (do_label && (nrow(lab_exp) || nrow(lab_ref))) 30 + lbl_gap * lvl_max else 0
  # Label layers, kept in a list so the main plot expression stays a plain `+`
  # chain (adding a list of layers is equivalent to adding them one by one).
  lab_layers <- list()
  if (do_label && nrow(lab_exp)) {
    lab_layers$exp <- ggplot2::geom_text(
      data = lab_exp,
      ggplot2::aes(x = .data[["mz"]], y = .data[["y_lab"]], label = .data[["lab"]]),
      angle = 90, hjust = 0, vjust = 0.5, size = 2.6, colour = "#B71C1C"
    )
  }
  if (do_label && nrow(lab_ref)) {
    lab_layers$ref <- ggplot2::geom_text(
      data = lab_ref,
      ggplot2::aes(x = .data[["mz"]], y = .data[["y_lab"]], label = .data[["lab"]]),
      angle = 90, hjust = 1, vjust = 0.5, size = 2.6, colour = "#0D47A1"
    )
  }

  ggplot2::ggplot(df) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
    ggplot2::geom_segment(
      data = df[df$source == "Experimental", , drop = FALSE],
      ggplot2::aes(
        x = .data[["mz"]], xend = .data[["mz"]],
        y = 0, yend = .data[["intensity"]]
      ),
      colour = "#D62728",
      linewidth = 0.45
    ) +
    ggplot2::geom_segment(
      data = df[df$source == "Reference", , drop = FALSE],
      ggplot2::aes(
        x = .data[["mz"]], xend = .data[["mz"]],
        y = 0, yend = .data[["intensity"]]
      ),
      colour = "#1F77B4",
      linewidth = 0.45
    ) +
    lab_layers +
    ggplot2::scale_y_continuous(
      limits = c(-105 - y_pad, 105 + y_pad),
      breaks = c(-100, -50, 0, 50, 100),
      labels = function(x) abs(x)
    ) +
    ggplot2::labs(
      title = paste0("MS2 Mirror Plot: ", feature_name),
      subtitle = subtitle,
      x = "m/z",
      y = "Relative Intensity"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 9, colour = "grey35"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

#' Plot a Kendrick Mass Defect (KMD) scatter plot
#'
#' Draws `mz` vs `KMD` after optional ISF filtering. Merges `ms1_annotations`
#' (KMD, ISF flags) and best-per-feature `db_annotations` scores for sizing.
#' Unknown features whose KMD aligns with library-annotated series are highlighted.
#'
#' @param project An [NTSProject] with `feature_table` and populated
#'   `ms1_annotations` after [calculate_kmd].
#' @param filter_isf Logical; if `TRUE`, removes rows with `is_isf == TRUE`.
#' @param only_show_series Logical; if `TRUE`, keeps only rows with
#'   non-missing `kmd_series_id`.
#' @param anchor_annotated Logical; if `TRUE`, keeps only homologous series that
#'   contain at least one annotated feature with `score_total > 0`.
#' @param target_kmds Optional numeric vector of target KMD values to retain.
#' @param kmd_tol Absolute tolerance around `target_kmds` (default `0.01`).
#' @param kmd_highlight_tol Absolute KMD tolerance used to highlight unknown
#'   points near library-annotated features (default `0.002`).
#'
#' @return A `ggplot2` object.
#' @export
plot_kmd <- function(project,
                     filter_isf = TRUE,
                     only_show_series = FALSE,
                     anchor_annotated = FALSE,
                     target_kmds = NULL,
                     kmd_tol = 0.01,
                     kmd_highlight_tol = 0.002) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_kmd()", call. = FALSE)
  }
  stopifnot(inherits(project, "NTSProject"))

  ft <- get_feature_table(project)
  if (!nrow(ft) || !"mz" %in% names(ft)) {
    stop("feature_table is empty or lacks mz; run calculate_kmd() first.", call. = FALSE)
  }
  ms1 <- get_ms1_annotations(project)
  if (!nrow(ms1) || !"KMD" %in% names(ms1)) {
    stop("ms1_annotations is empty or lacks KMD; run calculate_kmd() first.", call. = FALSE)
  }

  df <- .nts_ft_ms1_view(project)
  if (!"name" %in% names(df)) {
    df$name <- sprintf("FT%05d", seq_len(nrow(df)))
  }

  ann <- methods::slot(project, "db_annotations")
  if (nrow(ann) && "name" %in% names(ann)) {
    keep_cols <- intersect(c("name", "score_total", "Compound.name"), names(ann))
    ann_sub <- ann[, keep_cols, drop = FALSE]
    ann_sub <- ann_sub[!duplicated(as.character(ann_sub$name)), , drop = FALSE]
    df <- merge(df, ann_sub, by = "name", all.x = TRUE, sort = FALSE)
  }
  df <- .nts_filter_kmd_plot_df(
    df,
    filter_isf = filter_isf,
    only_show_series = only_show_series,
    anchor_annotated = anchor_annotated,
    target_kmds = target_kmds,
    kmd_tol = kmd_tol
  )
  if (!nrow(df)) {
    message("No points left after filtering; skipping plot.")
    return(invisible(NULL))
  }

  num_cols <- names(df)[vapply(df, is.numeric, logical(1L))]
  abundance_cols <- setdiff(
    num_cols,
    c("mz", "rt", "KM", "NKM", "KMD", "parent_mz", "score_total", "standard_md")
  )
  if ("score_total" %in% names(df) && any(is.finite(df$score_total))) {
    df$point_size <- ifelse(is.finite(df$score_total), df$score_total, 0)
  } else if (length(abundance_cols)) {
    df$point_size <- rowMeans(as.matrix(df[, abundance_cols, drop = FALSE]), na.rm = TRUE)
  } else {
    df$point_size <- 1
  }

  identified <- rep(FALSE, nrow(df))
  if ("score_total" %in% names(df) && "Compound.name" %in% names(df)) {
    identified <- is.finite(df$score_total) &
      df$score_total > 0.5 &
      !is.na(df$Compound.name) &
      nzchar(as.character(df$Compound.name))
  }

  aligned_unknown <- rep(FALSE, nrow(df))
  if (any(identified)) {
    ref_kmd <- df$KMD[identified]
    ref_kmd <- ref_kmd[is.finite(ref_kmd)]
    if (length(ref_kmd)) {
      kmd_dist <- abs(outer(df$KMD, ref_kmd, "-"))
      aligned_unknown <- rowSums(kmd_dist <= kmd_highlight_tol, na.rm = TRUE) > 0
      aligned_unknown <- aligned_unknown & !identified
    }
  }

  if (!"standard_md" %in% names(df)) {
    df$standard_md <- df$mz - round(df$mz)
  }
  md_ok <- is.finite(df$standard_md) & df$standard_md >= -0.15 & df$standard_md <= 0.15
  in_series <- "kmd_series_id" %in% names(df) &
    !is.na(df$kmd_series_id) &
    nzchar(as.character(df$kmd_series_id))

  df$point_type <- "Isolated Features"
  df$point_type[!md_ok] <- "Filtered: High Mass Defect"
  df$point_type[md_ok & in_series] <- "Series-resonant"
  df$point_type[identified] <- "Library annotated"
  df$point_type <- factor(
    df$point_type,
    levels = c(
      "Library annotated",
      "Series-resonant",
      "Isolated Features",
      "Filtered: High Mass Defect"
    )
  )

  total_n <- nrow(df)
  point_counts <- table(df$point_type)
  annotated_n <- unname(point_counts["Library annotated"])
  series_resonant_n <- unname(point_counts["Series-resonant"])
  isolated_n <- unname(point_counts["Isolated Features"])
  filtered_n <- unname(point_counts["Filtered: High Mass Defect"])
  annotated_n <- if (length(annotated_n)) annotated_n else 0L
  series_resonant_n <- if (length(series_resonant_n)) series_resonant_n else 0L
  isolated_n <- if (length(isolated_n)) isolated_n else 0L
  filtered_n <- if (length(filtered_n)) filtered_n else 0L
  filtered_pct <- if (total_n > 0) 100 * filtered_n / total_n else 0
  view_tag <- paste0(
    if (isTRUE(filter_isf)) "ISF-filtered view" else "All-features view",
    if (isTRUE(only_show_series)) "; series-only" else "",
    if (isTRUE(anchor_annotated)) "; anchored" else "",
    if (!is.null(target_kmds) && length(target_kmds)) paste0("; target KMD +/- ", kmd_tol) else ""
  )
  subtitle_text <- sprintf(
    "%s | Total: %d | Library annotated: %d | Series-resonant: %d | Isolated: %d | Filtered background: %d (%.1f%%)",
    view_tag,
    total_n,
    annotated_n,
    series_resonant_n,
    isolated_n,
    filtered_n,
    filtered_pct
  )

  color_map <- c(
    "Library annotated" = "#D62728",
    "Series-resonant" = "#111111",
    "Isolated Features" = "#BDBDBD",
    "Filtered: High Mass Defect" = "#E6E6E6"
  )
  alpha_map <- c(
    "Library annotated" = 1.0,
    "Series-resonant" = 0.8,
    "Isolated Features" = 0.4,
    "Filtered: High Mass Defect" = 0.1
  )
  size_map <- c(
    "Library annotated" = 3.8,
    "Series-resonant" = 2.7,
    "Isolated Features" = 2.0,
    "Filtered: High Mass Defect" = 1.6
  )

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["mz"]],
    y = .data[["KMD"]],
    colour = .data[["point_type"]],
    alpha = .data[["point_type"]],
    size = .data[["point_type"]]
  )) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = color_map, drop = FALSE) +
    ggplot2::scale_alpha_manual(values = alpha_map, drop = FALSE) +
    ggplot2::scale_size_manual(values = size_map, drop = FALSE) +
    ggplot2::labs(
      title = "Kendrick mass defect (KMD) plot",
      subtitle = subtitle_text,
      x = "m/z",
      y = "KMD",
      colour = "Feature Class",
      alpha = "Feature Class",
      size = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey35"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top"
    ) +
    ggplot2::guides(
      alpha = "none",
      size = "none"
    )
}

#' Plot homologous-series m/z-RT trend after KMD grouping
#'
#' Keeps only features assigned to a `kmd_series_id`, optionally removes ISF,
#' and visualizes the retention trend of each homologous series in `mz-rt`
#' space.
#'
#' @param project An [NTSProject] with `feature_table` containing `mz`, `rt`,
#'   and `kmd_series_id`.
#' @param filter_isf Logical; if `TRUE`, removes rows with `is_isf == TRUE`.
#' @param anchor_annotated Logical; if `TRUE`, keeps only homologous series that
#'   contain at least one annotated feature with `score_total > 0`.
#' @param target_kmds Optional numeric vector of target KMD values to retain.
#' @param kmd_tol Absolute tolerance around `target_kmds` (default `0.01`).
#'
#' @return A `ggplot2` object.
#' @export
plot_kmd_rt_series <- function(project,
                               filter_isf = TRUE,
                               anchor_annotated = FALSE,
                               target_kmds = NULL,
                               kmd_tol = 0.01) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_kmd_rt_series()", call. = FALSE)
  }
  stopifnot(inherits(project, "NTSProject"))

  ft <- get_feature_table(project)
  ms1 <- get_ms1_annotations(project)
  if (!nrow(ft) || !nrow(ms1) || !all(c("mz", "rt", "kmd_series_id") %in% names(ms1))) {
    stop("feature_table / ms1_annotations missing mz/rt/kmd_series_id; run calculate_kmd() first.", call. = FALSE)
  }

  df <- .nts_ft_ms1_view(project)
  if (!"name" %in% names(df)) {
    df$name <- sprintf("FT%05d", seq_len(nrow(df)))
  }
  ann <- methods::slot(project, "db_annotations")
  if (nrow(ann) && "name" %in% names(ann)) {
    keep_cols <- intersect(c("name", "score_total", "Compound.name"), names(ann))
    ann_sub <- ann[, keep_cols, drop = FALSE]
    ann_sub <- ann_sub[!duplicated(as.character(ann_sub$name)), , drop = FALSE]
    df <- merge(df, ann_sub, by = "name", all.x = TRUE, sort = FALSE)
  }
  df <- .nts_filter_kmd_plot_df(
    df,
    filter_isf = filter_isf,
    only_show_series = TRUE,
    anchor_annotated = anchor_annotated,
    target_kmds = target_kmds,
    kmd_tol = kmd_tol
  )
  if (!nrow(df)) {
    message("No points left after filtering; skipping plot.")
    return(invisible(NULL))
  }

  df$kmd_series_id <- factor(df$kmd_series_id)
  df <- df[order(df$kmd_series_id, df$mz), , drop = FALSE]
  grp_n <- table(df$kmd_series_id)
  line_df <- df[df$kmd_series_id %in% names(grp_n)[grp_n >= 2L], , drop = FALSE]
  if (!nrow(line_df)) {
    message("No KMD series with >=2 points for trend lines; skipping plot.")
    return(invisible(NULL))
  }
  line_df <- line_df[order(line_df$kmd_series_id, line_df$mz), , drop = FALSE]

  ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["mz"]],
    y = .data[["rt"]],
    colour = .data[["kmd_series_id"]],
    group = .data[["kmd_series_id"]]
  )) +
    ggplot2::geom_line(
      data = line_df,
      alpha = 0.55,
      linewidth = 0.45
    ) +
    ggplot2::geom_point(size = 2.2, alpha = 0.9) +
    ggplot2::geom_smooth(
      data = line_df,
      method = "lm",
      se = FALSE,
      linetype = "dashed",
      linewidth = 0.5,
      alpha = 0.7
    ) +
    ggplot2::labs(
      title = "KMD Homologous Series m/z-RT Trend",
      subtitle = if (isTRUE(filter_isf)) {
        "Only non-ISF features assigned to homologous series are shown"
      } else {
        "Only features assigned to homologous series are shown"
      },
      x = "m/z",
      y = "Retention Time (s)",
      colour = "Feature Class"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10, colour = "grey35"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
}

#' MS2 extraction efficiency diagnostic (funnel chart)
#'
#' Reads MS2 mzML headers to count raw MS1/MS2 scans, compares against the
#' number of consensus MS2 spectra actually attached to the project, and plots
#' a descending bar chart showing the data loss funnel.
#'
#' @param project An [NTSProject] with `@ms2data` populated.
#' @param ms2_dir Directory containing MS2 mzML files. Defaults to a `MS2`
#'   sibling of `@raw_data_path`.
#'
#' @return A `ggplot2` object (invisibly also prints a console summary).
#' @export
plot_ms2_efficiency <- function(project, ms2_dir = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_ms2_efficiency()", call. = FALSE)
  }
  stopifnot(inherits(project, "NTSProject"))

  if (is.null(ms2_dir)) {
    base <- methods::slot(project, "raw_data_path")
    candidate <- file.path(dirname(base), "MS2")
    if (dir.exists(candidate)) {
      ms2_dir <- candidate
    } else {
      stop("ms2_dir not specified and no MS2 sibling found next to @raw_data_path",
           call. = FALSE)
    }
  }

  ms2_files <- .nts_list_ms_files(ms2_dir)
  if (!length(ms2_files)) {
    stop("No mzML/mzXML files found in: ", ms2_dir, call. = FALSE)
  }

  message("> Scanning headers of ", length(ms2_files), " MS2 file(s)...")
  total_ms1 <- 0L
  total_ms2 <- 0L
  for (i in seq_along(ms2_files)) {
    hdr <- tryCatch({
      raw <- MSnbase::readMSData(ms2_files[i], mode = "onDisk")
      MSnbase::fData(raw)
    }, error = function(e) NULL)
    if (is.null(hdr) || !"msLevel" %in% names(hdr)) next
    lvl <- hdr$msLevel
    total_ms1 <- total_ms1 + sum(lvl == 1L)
    total_ms2 <- total_ms2 + sum(lvl == 2L)
  }

  n_features <- nrow(get_feature_table(project))
  n_extracted <- length(methods::slot(project, "ms2data"))
  rate <- if (total_ms2 > 0) round(n_extracted / total_ms2 * 100, 2) else NA_real_

  cat("\n", strrep("=", 50), "\n")
  cat("  MS2 Extraction Efficiency Report\n")
  cat(strrep("=", 50), "\n")
  cat("  MS2 files scanned     :", length(ms2_files), "\n")
  cat("  Raw MS1 scans (in MS2):", total_ms1, "\n")
  cat("  Raw MS2 scans         :", total_ms2, "\n")
  cat("  MS1 features          :", n_features, "\n")
  cat("  Extracted MS2 (unique):", n_extracted, "\n")
  cat("  Extraction rate       :", if (is.na(rate)) "N/A" else paste0(rate, "%"), "\n")
  cat(strrep("=", 50), "\n\n")

  df <- data.frame(
    stage = factor(
      c("Raw MS1 Scans", "Raw MS2 Scans", "Extracted MS2"),
      levels = c("Raw MS1 Scans", "Raw MS2 Scans", "Extracted MS2")
    ),
    count = c(total_ms1, total_ms2, n_extracted),
    stringsAsFactors = FALSE
  )

  fills <- c("#4E79A7", "#F28E2B", "#59A14F")

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x = .data[["stage"]],
    y = .data[["count"]],
    fill = .data[["stage"]]
  )) +
    ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
    ggplot2::geom_text(
      ggplot2::aes(label = format(.data[["count"]], big.mark = ",")),
      vjust = -0.5,
      size = 4.5,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = fills) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.15)),
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    ggplot2::labs(
      title = "DDA MS2 Extraction Funnel",
      subtitle = sprintf("Extraction rate: %s", if (is.na(rate)) "N/A" else paste0(rate, "%%")),
      x = NULL,
      y = "Count"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(colour = "grey40", size = 11),
      axis.text.x = ggplot2::element_text(face = "bold", size = 11),
      panel.grid.major.x = ggplot2::element_blank()
    )
  p
}

#' Evidence-based plot groups (see `.nts_network_nodes_evidence_group`)
#' @noRd
.nts_network_nodes_plot_group <- function(project, nodes_dt) {
  .nts_network_nodes_evidence_group(project, nodes_dt)
}

#' Strict node evidence tiers (data.table). Seed (DB Matched) when
#' \code{Library_ID} or \code{Compound.name} is non-empty after trim.
#'
#' Priority: Seed (DB Matched) > Propagated > MS2-Rich (Diagnostic) >
#' Homologue Series > ISF > Unknown.
#' @noRd
.nts_network_nodes_evidence_group <- function(project, nodes_dt) {
  force(project)
  if (!nrow(nodes_dt)) {
    return(character(0))
  }
  .nts_check(NULL, "data.table_ns")
  DT <- data.table::copy(data.table::as.data.table(nodes_dt))
  if (!"Library_ID" %in% names(DT)) {
    DT[, Library_ID := NA_character_]
  }
  if (!"Compound.name" %in% names(DT)) {
    DT[, Compound.name := NA_character_]
  }
  if (!"Propagation_Level" %in% names(DT)) {
    DT[, Propagation_Level := NA_integer_]
  }
  if (!"diagnostic_fragments" %in% names(DT)) {
    DT[, diagnostic_fragments := NA_character_]
  }
  if (!"n_fragments" %in% names(DT)) {
    DT[, n_fragments := 0]
  }
  if (!"mz_remainder" %in% names(DT)) {
    DT[, mz_remainder := NA_character_]
  }
  if (!"is_isf" %in% names(DT)) {
    DT[, is_isf := FALSE]
  }
  DT[, pl_i := suppressWarnings(as.integer(Propagation_Level))]
  DT[
    ,
    group := data.table::fcase(
      (!is.na(Library_ID) & nzchar(trimws(as.character(Library_ID)))) |
        (!is.na(Compound.name) & nzchar(trimws(as.character(Compound.name)))) |
        grepl("Database Confirmed", trimws(as.character(annotation_level)), fixed = FALSE),
      "Seed (DB Matched)",
      !is.na(pl_i) & pl_i > 0L,
      "Propagated",
      (!is.na(diagnostic_fragments) &
        nzchar(trimws(as.character(diagnostic_fragments)))) |
        (is.finite(suppressWarnings(as.numeric(n_fragments))) &
          suppressWarnings(as.numeric(n_fragments)) > 1),
      "MS2-Rich (Diagnostic)",
      !is.na(mz_remainder) & nzchar(trimws(as.character(mz_remainder))),
      "Homologue Series",
      is_isf %in% TRUE,
      "ISF",
      default = "Unknown"
    )
  ]
  out <- as.character(DT$group)
  DT[, pl_i := NULL]
  out
}

#' Fixed semantic groups and colors; Seed / legacy long labels use anchor red.
#' @noRd
.nts_network_group_color_palette <- function() {
  c(
    "Seed (DB Matched)" = "#B2182B",
    "Propagated" = "#FF7F0E",
    "MS2-Rich (Diagnostic)" = "#9467BD",
    "MS2-Rich" = "#9467BD",
    "Homologue Series" = "#1F78B4",
    "Homologue" = "#1F78B4",
    "ISF" = "#33A02C",
    "Unknown" = "#D9D9D9",
    "RootNode" = "#FF8C00",
    "Subnet_Root" = "#FF8C00"
  )
}

#' @noRd
.nts_network_nodes_apply_group_colors <- function(nodes_df) {
  if (!nrow(nodes_df) || !"group" %in% names(nodes_df)) {
    return(nodes_df)
  }
  pal <- .nts_network_group_color_palette()
  bg <- pal[as.character(nodes_df$group)]
  bg[is.na(bg)] <- "#B2182B"
  nodes_df$color.background <- unname(bg)
  nodes_df$color.border <- "#404040"
  nodes_df
}

#' `+`-joined edge / interaction token match (merged-layer collapse)
#' @noRd
.nts_edge_has_token <- function(vec, tok) {
  vapply(vec, function(s) {
    if (length(s) != 1L) {
      return(FALSE)
    }
    s <- as.character(s)
    if (is.na(s) || !nzchar(s)) {
      return(FALSE)
    }
    tok %in% trimws(strsplit(s, " + ", fixed = TRUE)[[1L]])
  }, logical(1L))
}

#' @noRd
.nts_edge_vis_dash_one <- function(s) {
  if (length(s) != 1L || is.na(s) || !nzchar(as.character(s))) {
    return(FALSE)
  }
  toks <- trimws(strsplit(as.character(s), " + ", fixed = TRUE)[[1L]])
  pick <- function(nm) isTRUE(nm %in% toks)
  if (pick("MS2_Shared_Fragment")) {
    return(c(2, 4))
  }
  if (pick("MS2_Shared_Neutral_Loss")) {
    return(c(3, 3))
  }
  if (pick("MS2_Shared_Mass_Diff")) {
    return(c(4, 2))
  }
  if (pick("MS2_Similarity")) {
    return(c(2, 4))
  }
  if (pick("MS2_Semantic")) {
    return(c(5, 5))
  }
  if (pick("KMD_Homologue")) {
    return(FALSE)
  }
  if (pick("Transformation")) {
    return(TRUE)
  }
  if (pick("ISF_Link")) {
    return(c(5, 3))
  }
  FALSE
}

#' @noRd
.nts_edge_vis_arrow_one <- function(s) {
  if (length(s) != 1L || is.na(s) || !nzchar(as.character(s))) {
    return("to")
  }
  toks <- trimws(strsplit(as.character(s), " + ", fixed = TRUE)[[1L]])
  undir <- c(
    "MS2_Similarity",
    "MS2_Shared_Fragment",
    "MS2_Shared_Neutral_Loss",
    "MS2_Shared_Mass_Diff",
    "MS2_Semantic",
    "ISF_Link"
  )
  if (any(undir %in% toks)) {
    # NOTE: vis-network >= 2.x throws
    #   "The arrow newOptions can only be an object or a string"
    # on a per-edge `arrows = FALSE` (it validates against object|string),
    # which aborts the *entire* network render. Use the empty string, which is
    # a valid string and parses to all arrowheads disabled (undirected look).
    return("")
  }
  "to"
}

#' @noRd
.nts_edge_network_vis <- function(edges_df, color_by_edge = TRUE) {
  if (!nrow(edges_df)) {
    return(edges_df)
  }
  # Tolerate arbitrary edge tables (e.g. a bare list(nodes, edges) without the
  # canonical layer columns): the token checks below must never see a
  # zero-length or all-NA `Interaction` (that would yield a length-0 / NA
  # `width` column and break visNetwork).
  if (!"Interaction" %in% names(edges_df)) {
    edges_df$Interaction <- ""
  }
  if (!"label" %in% names(edges_df)) {
    edges_df$label <- ""
  }
  intr <- as.character(edges_df$Interaction)
  intr[is.na(intr)] <- ""
  if (!isTRUE(color_by_edge)) {
    edges_df$color <- "#AAAAAA"
  } else {
    et <- if ("edge_type" %in% names(edges_df)) {
      as.character(edges_df$edge_type)
    } else {
      rep(NA_character_, nrow(edges_df))
    }
    mt <- if ("match_type" %in% names(edges_df)) {
      as.character(edges_df$match_type)
    } else {
      rep(NA_character_, nrow(edges_df))
    }
    colv <- rep("#999999", nrow(edges_df))
    # MS1 edge_type colours (supports Merged-layer `a + b` strings)
    w <- !is.na(et) & .nts_edge_has_token(et, "Homologue")
    colv[w] <- "#2CA02C"
    w <- !is.na(et) & .nts_edge_has_token(et, "Transformation")
    colv[w] <- "#D62728"
    w <- !is.na(et) & .nts_edge_has_token(et, "ISF")
    colv[w] <- "#FF9896"
    # MS2 match_type colours
    w <- !is.na(mt) & .nts_edge_has_token(mt, "fragment")
    colv[w] <- "#7F7F7F"
    w <- !is.na(mt) & .nts_edge_has_token(mt, "neutral_loss")
    colv[w] <- "#9467BD"
    w <- !is.na(mt) & .nts_edge_has_token(mt, "mass_difference")
    colv[w] <- "#17BECF"
    w <- !is.na(mt) & .nts_edge_has_token(mt, "semantic")
    colv[w] <- "#9370DB"
    # Fallback: Interaction
    is_ms2_intr <- .nts_edge_has_token(intr, "MS2_Similarity") |
      .nts_edge_has_token(intr, "MS2_Shared_Fragment") |
      .nts_edge_has_token(intr, "MS2_Shared_Neutral_Loss") |
      .nts_edge_has_token(intr, "MS2_Shared_Mass_Diff") |
      .nts_edge_has_token(intr, "MS2_Semantic")
    fb <- colv == "#999999"
    colv[fb & .nts_edge_has_token(intr, "KMD_Homologue")] <- "#2CA02C"
    colv[fb & .nts_edge_has_token(intr, "Transformation")] <- "#D62728"
    colv[fb & .nts_edge_has_token(intr, "ISF_Link")] <- "#FF9896"
    colv[fb & .nts_edge_has_token(intr, "MS2_Shared_Fragment")] <- "#7F7F7F"
    colv[fb & .nts_edge_has_token(intr, "MS2_Shared_Neutral_Loss")] <- "#9467BD"
    colv[fb & .nts_edge_has_token(intr, "MS2_Shared_Mass_Diff")] <- "#17BECF"
    colv[fb & .nts_edge_has_token(intr, "MS2_Semantic")] <- "#9370DB"
    colv[fb & is_ms2_intr] <- "#7F7F7F"
    edges_df$color <- colv
  }
  intr <- as.character(edges_df$Interaction)
  intr[is.na(intr)] <- ""
  mt_w <- if ("match_type" %in% names(edges_df)) {
    as.character(edges_df$match_type)
  } else {
    rep(NA_character_, nrow(edges_df))
  }
  is_ms2 <- .nts_edge_has_token(mt_w, "fragment") |
    .nts_edge_has_token(mt_w, "neutral_loss") |
    .nts_edge_has_token(mt_w, "mass_difference") |
    .nts_edge_has_token(mt_w, "semantic") |
    grepl("MS2_", intr, fixed = TRUE)
  edges_df$width <- ifelse(is_ms2, 1, 2)
  edges_df$dashes <- I(lapply(intr, .nts_edge_vis_dash_one))
  edges_df$arrows <- I(lapply(intr, .nts_edge_vis_arrow_one))
  ttl <- rep("", nrow(edges_df))
  if ("shared_fragments_mz" %in% names(edges_df)) {
    sf <- as.character(edges_df$shared_fragments_mz)
    ok <- !is.na(sf) & nzchar(sf)
    ttl[ok] <- paste0(
      ttl[ok],
      "<b>Shared fragments (m/z)</b><br>",
      sf[ok],
      "<br>"
    )
  }
  if ("shared_nls_mz" %in% names(edges_df)) {
    sn <- as.character(edges_df$shared_nls_mz)
    ok <- !is.na(sn) & nzchar(sn)
    ttl[ok] <- paste0(
      ttl[ok],
      "<b>Shared neutral losses (m/z)</b><br>",
      sn[ok],
      "<br>"
    )
  }
  if ("shared_mds_mz" %in% names(edges_df)) {
    mm <- as.character(edges_df$shared_mds_mz)
    intr_md <- as.character(edges_df$Interaction)
    mt_md <- if ("match_type" %in% names(edges_df)) {
      as.character(edges_df$match_type)
    } else {
      rep(NA_character_, nrow(edges_df))
    }
    okm <- !is.na(mm) & nzchar(mm) & (
      .nts_edge_has_token(intr_md, "MS2_Shared_Mass_Diff") |
        .nts_edge_has_token(mt_md, "mass_difference")
    )
    ttl[okm] <- paste0(
      ttl[okm],
      "<b>Shared mass differences (Da)</b><br>",
      mm[okm],
      "<br>"
    )
  }
  if ("semantic_similarity" %in% names(edges_df)) {
    ss <- suppressWarnings(as.numeric(edges_df$semantic_similarity))
    ok_ss <- is.finite(ss) & ss > 0
    ttl[ok_ss] <- paste0(
      ttl[ok_ss],
      "<b>Semantic similarity</b><br>",
      sprintf("%.3f", ss[ok_ss]),
      "<br>"
    )
  }
  lab <- as.character(edges_df$label)
  okl <- !is.na(lab) & nzchar(lab)
  ttl[okl] <- paste0(ttl[okl], "<b>Label</b><br>", lab[okl])
  edges_df$title <- ttl
  edges_df
}
#' @noRd
.nts_network_apply_show_iteration <- function(net_data, show_iteration) {
  si <- trimws(as.character(show_iteration)[1L])
  if (!nzchar(si) || tolower(si) == "all") {
    return(net_data)
  }
  nodes_df <- net_data$nodes
  edges_df <- net_data$edges
  if (!"Propagation_Level" %in% names(nodes_df)) {
    warning(
      "plot_molecular_network: show_iteration is not \"all\" but nodes lack Propagation_Level; showing full network.",
      call. = FALSE,
      immediate. = TRUE
    )
    return(net_data)
  }
  pl <- suppressWarnings(as.integer(nodes_df$Propagation_Level))
  nid <- as.character(nodes_df$id)
  keep <- rep(FALSE, nrow(nodes_df))
  sil <- tolower(si)
  if (sil %in% c("seeds_only", "iteration_0")) {
    keep <- !is.na(pl) & pl == 0L
  } else if (grepl("^iteration_[0-9]+$", sil)) {
    k <- as.integer(sub("^iteration_([0-9]+)$", "\\1", sil))
    keep <- !is.na(pl) & pl >= 0L & pl <= k
  } else {
    warning(
      sprintf("plot_molecular_network: unknown show_iteration=%s; using all.", si),
      call. = FALSE,
      immediate. = TRUE
    )
    return(net_data)
  }
  if (!any(keep)) {
    net_data$nodes <- nodes_df[FALSE, , drop = FALSE]
    net_data$edges <- edges_df[FALSE, , drop = FALSE]
    return(net_data)
  }
  vid <- nid[keep]
  kn <- nodes_df[keep, , drop = FALSE]
  if (!nrow(edges_df) || !all(c("from", "to") %in% names(edges_df))) {
    net_data$nodes <- kn
    net_data$edges <- edges_df
    return(net_data)
  }
  ff <- as.character(edges_df$from)
  tt <- as.character(edges_df$to)
  ke <- ff %in% vid & tt %in% vid
  net_data$nodes <- kn
  net_data$edges <- edges_df[ke, , drop = FALSE]
  net_data
}

#' Add a second "Iteration" dropdown to a visNetwork widget
#'
#' [visNetwork::visOptions] `selectedBy` renders a single selection list, so a
#' **second** list is injected with [htmlwidgets::onRender]. It enumerates the
#' distinct non-`NA` integer values of `iteration_col` on the nodes and, on
#' change, hides every node whose level exceeds the chosen value (vis.js
#' `hidden`, which also hides incident edges). Selecting the header (empty value)
#' restores the full graph.
#'
#' No-op — the widget is returned unchanged — when `iteration_col` is `NULL`,
#' absent from the node table, non-integer, or carries fewer than two distinct
#' levels. Graphs whose `Propagation_Level` is entirely `NA` (the historical
#' default) therefore render exactly as before.
#'
#' @param p A **visNetwork** htmlwidget.
#' @param nodes_dt Node `data.frame` used to build `p` (level column read here).
#' @param iteration_col Column name holding the propagation level, or `NULL`.
#' @return `p`, optionally with an `onRender` hook appended.
#' @noRd
.nts_network_add_iteration_selector <- function(p, nodes_dt, iteration_col) {
  if (is.null(p) || is.null(iteration_col) || !length(iteration_col)) {
    return(p)
  }
  iteration_col <- as.character(iteration_col)[[1L]]
  # visNetworkdataframeToD3 splits column names on ".", so only flat names are
  # safe as a JS-readable node field.
  if (!nzchar(iteration_col) ||
      !grepl("^[A-Za-z_][A-Za-z0-9_]*$", iteration_col) ||
      !iteration_col %in% names(nodes_dt)) {
    return(p)
  }
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
    return(p)
  }
  lv <- suppressWarnings(as.integer(nodes_dt[[iteration_col]]))
  if (length(unique(lv[!is.na(lv)])) < 2L) {
    return(p)
  }
  js <- sprintf(
    "function(el, x) {
      var root = document.getElementById(el.id) || el;
      var gdiv = document.getElementById('graph' + el.id);
      var chart = (gdiv && gdiv.chart) || el.chart;
      if (!chart || !chart.body || !chart.body.data) { return; }
      var nodesDS = chart.body.data.nodes;
      var all = nodesDS.get();
      var VAR = '%s';
      var lvMap = {};
      var levels = [];
      for (var i = 0; i < all.length; i++) {
        var v = all[i][VAR];
        if (v === null || v === undefined || v === '') { continue; }
        var n = parseInt(v, 10);
        if (isNaN(n)) { continue; }
        lvMap[all[i].id] = n;
        if (levels.indexOf(n) === -1) { levels.push(n); }
      }
      if (levels.length < 2) { return; }
      levels.sort(function(a, b) { return a - b; });
      var sel = document.createElement('select');
      sel.className = 'dropdown';
      sel.id = 'iterationSelect' + el.id;
      sel.style.marginLeft = '6px';
      sel.style.display = 'inline';
      var o0 = document.createElement('option');
      o0.value = ''; o0.text = 'Iteration (propagation)';
      sel.appendChild(o0);
      for (var j = 0; j < levels.length; j++) {
        var o = document.createElement('option');
        o.value = String(levels[j]);
        o.text = levels[j] === 0 ? '0 (seeds)' : ('\\u2264 ' + levels[j]);
        sel.appendChild(o);
      }
      var anchor = document.getElementById('maindiv' + el.id);
      if (anchor && anchor.parentNode === root) { root.insertBefore(sel, anchor); }
      else { root.appendChild(sel); }
      sel.onchange = function() {
        var lim = sel.value === '' ? null : parseInt(sel.value, 10);
        var upd = [];
        for (var k = 0; k < all.length; k++) {
          var id = all[k].id;
          var lv = lvMap[id];
          var hidden = (lim === null) ? false : !(lv !== undefined && lv <= lim);
          upd.push({ id: id, hidden: hidden });
        }
        nodesDS.update(upd);
        if (chart.fit) { chart.fit(); }
      };
    }", iteration_col)
  htmlwidgets::onRender(p, htmlwidgets::JS(js))
}

#' @title Plot a molecular network (visNetwork)
#'
#' @description
#' Builds an interactive **visNetwork** graph from precomputed `nodes` / `edges`
#' tables: node tooltips (`title`), evidence-based colours when `group` is
#' present, undirected-style edges, and `forceAtlas2Based` physics with
#' stabilization then physics disabled via [visNetwork::visEvents] for a stable
#' layout drag experience. Edges are styled per layer via
#' `.nts_edge_network_vis` (Homologue / Transformation / ISF / MS2 colour,
#' dashes, arrow direction and hover titles) unless `color_by_edge = FALSE`.
#' Node tooltips also show a **PFAS ladder** line when a node carries an MS2 CF2
#' fragment ladder (`pfas_ladder_len >= 3`, from
#' [detect_pfas_fragment_series]), including the `pfas_subclass_hint` when set,
#' and a **KMD series** line when the node carries a non-empty `kmd_series_id`.
#' Node **labels** fall back from the homologue inference label (`homologue_label`)
#' to the database-annotation name (`Compound.name`, then `Annotation`), so
#' DB-matched seeds are labelled in the MS2 feature and merged graphs too; the
#' tooltip gains a **DB match** line (`name (Library_ID)`) when those columns are
#' populated.
#' When `iteration_col` carries two or more distinct levels, a second dropdown
#' steps through them by hiding higher-level nodes (see
#' `.nts_network_add_iteration_selector`).
#' With `remove_isolated = TRUE`, node ids listed in `keep_nodes` are retained
#' even without an incident edge, so annotated seeds are never hidden merely
#' because they form no edge in this layer.
#'
#' @param x Network `list` from [get_network_data], [generate_network_data], or
#'   [extract_subnetwork], with elements `nodes` and `edges` (coerced with
#'   [as.data.frame]). An [NTSProject] is **not** accepted here; call
#'   [generate_network_data] or [plot_nts_network] first.
#' @param remove_isolated Logical; when `TRUE` (default), drops nodes with no
#'   incident edge after both endpoints are matched on `id`.
#' @param keep_nodes Optional character vector of node `id`s that are **never**
#'   dropped by `remove_isolated`, even when they have no incident edge. Use it
#'   to keep seed / annotated nodes visible in a sparse graph (e.g. pass the
#'   database-annotation feature ids so every matched seed is shown, rather than
#'   only the subset that happens to form an edge). Ids not present in `nodes`
#'   are ignored.
#' @param search_by_name Logical; when `TRUE`, enables [visNetwork::visOptions]
#'   `nodesIdSelection` so users can jump to a vertex by `id` (visNetwork
#'   `label` defaults to `id`). When `FALSE`, the ID selector is hidden.
#' @param color_by_edge Logical; when `TRUE` (default), styles edges by layer via
#'   `.nts_edge_network_vis` (colour / dashes / arrows / hover `title`). When
#'   `FALSE`, draws all edges in a flat grey. [plot_nts_network] forwards this
#'   argument to [plot_molecular_network].
#' @param height,width CSS size passed to [visNetwork::visNetwork]. Defaults match
#'   the historical widget size (`"900px"`); callers embedding the widget in a
#'   fixed-height container (e.g. Shiny) should pass the container height.
#' @param iteration_col Name of a node column holding a propagation level
#'   (default `"Propagation_Level"`). When it is present and carries at least two
#'   distinct non-`NA` integer values, a **second** dropdown ("Iteration
#'   (propagation)") is added next to the group selector: choosing level `k` hides
#'   every node whose level exceeds `k` (the header restores the full graph). Set
#'   to `NULL` to suppress. Because the historical `Propagation_Level` column is
#'   `NA` throughout, graphs built straight from [generate_network_data] are
#'   unaffected; populate the column (e.g. with a per-node hop count) to enable it.
#' @param ... Absorbs legacy unused arguments (e.g. `ms2_edge_types`) for
#'   backward-compatible scripts.
#'
#' @return A **visNetwork** htmlwidget. When no rows remain in `nodes` after
#'   filtering, returns `NULL` `invisibly`.
#'
#' @importFrom visNetwork visEdges visEvents visNetwork visNodes visOptions visPhysics
#' @export
plot_molecular_network <- function(x, remove_isolated = TRUE, keep_nodes = NULL,
                                   search_by_name = TRUE,
                                   color_by_edge = TRUE, height = "900px",
                                   width = "100%", iteration_col = "Propagation_Level",
                                   ...) {
  if (!requireNamespace("visNetwork", quietly = TRUE)) {
    stop("install.packages(\"visNetwork\") is required for plot_molecular_network()", call. = FALSE)
  }
  if (!is.list(x) || !all(c("nodes", "edges") %in% names(x))) {
    stop(
      "x must be a list with elements `nodes` and `edges` (e.g. from get_network_data() or extract_subnetwork()).",
      call. = FALSE
    )
  }

  nodes_dt <- as.data.frame(x$nodes, stringsAsFactors = FALSE)
  edges_dt <- as.data.frame(x$edges, stringsAsFactors = FALSE)

  if (!nrow(nodes_dt)) {
    return(invisible(NULL))
  }

  # 1. Ensure core ID columns exist
  if (!"id" %in% colnames(nodes_dt) && "name" %in% colnames(nodes_dt)) {
    nodes_dt$id <- nodes_dt$name
  }
  if (!"id" %in% colnames(nodes_dt)) {
    stop("nodes must contain column `id` or `name`.", call. = FALSE)
  }

  # 2. Drop isolated nodes (minimal safe filter). `keep_nodes` (e.g. database
  # matched seed ids) are exempt: an annotated seed stays visible even when it
  # forms no edge in this particular layer, so the displayed seed count matches
  # the annotation count instead of the (sparser) edge-connected subset.
  if (isTRUE(remove_isolated) && nrow(edges_dt) > 0 && all(c("from", "to") %in% colnames(edges_dt))) {
    connected_ids <- unique(c(edges_dt$from, edges_dt$to))
    if (length(keep_nodes)) {
      kp <- as.character(keep_nodes)
      connected_ids <- union(connected_ids, kp)
    }
    nodes_dt <- nodes_dt[nodes_dt$id %in% connected_ids, , drop = FALSE]
  }

  if (!nrow(nodes_dt)) {
    return(invisible(NULL))
  }

  # 3. Coerce types for JSON (avoid Factor columns breaking visNetwork/JS)
  nodes_dt$id <- as.character(nodes_dt$id)
  if (nrow(edges_dt) && all(c("from", "to") %in% colnames(edges_dt))) {
    edges_dt$from <- as.character(edges_dt$from)
    edges_dt$to <- as.character(edges_dt$to)
  }

  # 3b. Style edges by layer (colour / dashes / arrows / hover title). This is
  # the same enrichment plot_nts_network relies on; running it here means every
  # caller (e.g. the Shiny Network tab, which calls plot_molecular_network()
  # directly) renders Homologue / Transformation / ISF / MS2 layers distinctly
  # instead of as a uniform grey web (P1-1). `.nts_edge_network_vis` is a no-op
  # on an empty edge table.
  if (nrow(edges_dt) && all(c("from", "to") %in% colnames(edges_dt))) {
    edges_dt <- .nts_edge_network_vis(edges_dt, color_by_edge = color_by_edge)
  }

  if (!"group" %in% colnames(nodes_dt)) {
    nodes_dt$group <- "Isolated Unknown"
  }
  nodes_dt <- .nts_network_nodes_apply_group_colors(nodes_dt)

  # 4. Trim columns and attach tooltip fields
  # Node label: prefer the homologue inference label ("X - n CF2"), then the
  # database-annotation name (Compound.name, then Annotation), then the bare id.
  nodes_dt$label <- nodes_dt$id
  if ("homologue_label" %in% colnames(nodes_dt)) {
    hl <- trimws(as.character(nodes_dt$homologue_label))
    ok_lbl <- !is.na(hl) & nzchar(hl) & toupper(hl) != "NA"
    nodes_dt$label[ok_lbl] <- hl[ok_lbl]
  }
  # Seed / DB-matched nodes: name them from the database annotation so the MS2
  # feature networks and the merged graph label their seeds too (the homologue
  # graph already gets this through homologue_label). Only nodes still showing the
  # bare id are overwritten; Compound.name wins over Annotation.
  bare <- nodes_dt$label == nodes_dt$id
  for (.cn_col in c("Compound.name", "Annotation")) {
    if (!.cn_col %in% colnames(nodes_dt)) next
    .v <- trimws(as.character(nodes_dt[[.cn_col]]))
    .ok <- bare & !is.na(.v) & nzchar(.v) & toupper(.v) != "NA"
    if (any(.ok)) {
      nodes_dt$label[.ok] <- .v[.ok]
      bare[.ok] <- FALSE
    }
  }
  nodes_dt$title <- paste0(
    "<div style='font-family:Arial; font-size:12px; padding:5px;'>",
    "<b>ID:</b> ", nodes_dt$id, "<br>",
    "<b>Group:</b> ", ifelse(is.na(nodes_dt$group), "Isolated Unknown", nodes_dt$group)
  )
  if ("mz" %in% colnames(nodes_dt)) {
    nodes_dt$title <- paste0(nodes_dt$title, "<br><b>m/z:</b> ", round(as.numeric(nodes_dt$mz), 4L))
  }
  # Database annotation (seed identity): show the matched compound name and its
  # library accession, so every graph names its DB-matched seeds, not just the
  # homologue graph. Compound.name wins over Annotation; Library_ID is optional.
  if (any(c("Compound.name", "Annotation") %in% colnames(nodes_dt))) {
    cn <- if ("Compound.name" %in% colnames(nodes_dt)) {
      trimws(as.character(nodes_dt$Compound.name))
    } else {
      rep(NA_character_, nrow(nodes_dt))
    }
    if ("Annotation" %in% colnames(nodes_dt)) {
      an <- trimws(as.character(nodes_dt$Annotation))
      na_cn <- is.na(cn) | !nzchar(cn) | toupper(cn) == "NA"
      use_an <- na_cn & !is.na(an) & nzchar(an) & toupper(an) != "NA"
      cn[use_an] <- an[use_an]
    }
    ok_db <- !is.na(cn) & nzchar(cn) & toupper(cn) != "NA"
    if (any(ok_db)) {
      lib <- if ("Library_ID" %in% colnames(nodes_dt)) {
        trimws(as.character(nodes_dt$Library_ID))
      } else {
        rep(NA_character_, nrow(nodes_dt))
      }
      ok_lib <- !is.na(lib) & nzchar(lib) & toupper(lib) != "NA"
      db_tag <- cn
      db_tag[ok_db & ok_lib] <- paste0(cn[ok_db & ok_lib], " (", lib[ok_db & ok_lib], ")")
      nodes_dt$title[ok_db] <- paste0(
        nodes_dt$title[ok_db], "<br><b>DB match:</b> ", db_tag[ok_db])
    }
  }
  if ("homologue_label" %in% colnames(nodes_dt)) {
    hl2 <- trimws(as.character(nodes_dt$homologue_label))
    ok2 <- !is.na(hl2) & nzchar(hl2)
    if (any(ok2)) {
      nodes_dt$title[ok2] <- paste0(
        nodes_dt$title[ok2],
        "<br><b>Homologue:</b> ", hl2[ok2]
      )
    }
  }
  if ("evidence_source" %in% colnames(nodes_dt)) {
    es <- trimws(as.character(nodes_dt$evidence_source))
    ok3 <- !is.na(es) & nzchar(es)
    if (any(ok3)) {
      nodes_dt$title[ok3] <- paste0(
        nodes_dt$title[ok3],
        "<br><b>Evidence:</b> ", es[ok3]
      )
    }
  }
  # (P2-4) Surface KMD homologue-series membership. Path-A nodes carry
  # `kmd_series_id` but no `homologue_label` (that column is only produced by
  # the infer_homologue_structures path), so without this line the App shows the
  # "Homologue Series" colour but never the series id.
  if ("kmd_series_id" %in% colnames(nodes_dt)) {
    ks <- trimws(as.character(nodes_dt$kmd_series_id))
    ok5 <- !is.na(ks) & nzchar(ks) & toupper(ks) != "NA"
    if (any(ok5)) {
      nodes_dt$title[ok5] <- paste0(
        nodes_dt$title[ok5],
        "<br><b>KMD series:</b> ", ks[ok5]
      )
    }
  }
  # (P2-3) Surface the library-free MS2 CF2 ladder as node-tooltip context.
  if ("pfas_ladder_len" %in% colnames(nodes_dt)) {
    ll <- suppressWarnings(as.numeric(nodes_dt$pfas_ladder_len))
    ok4 <- is.finite(ll) & ll >= 3
    if (any(ok4)) {
      lab <- paste0(as.integer(ll[ok4]), " CF2 ladder")
      if ("pfas_subclass_hint" %in% colnames(nodes_dt)) {
        sh <- trimws(as.character(nodes_dt$pfas_subclass_hint[ok4]))
        has <- !is.na(sh) & nzchar(sh) & toupper(sh) != "NA"
        lab[has] <- paste0(lab[has], " (", sh[has], ")")
      }
      nodes_dt$title[ok4] <- paste0(
        nodes_dt$title[ok4],
        "<br><b>PFAS ladder:</b> ", lab
      )
    }
  }
  nodes_dt$title <- paste0(nodes_dt$title, "</div>")

  # 5. Physics layout (forceAtlas2Based; not visIgraphLayout)
  p <- visNetwork::visNetwork(nodes_dt, edges_dt, width = width, height = height) |>
    visNetwork::visNodes(
      shape = "dot",
      size = 20,
      borderWidth = 2,
      shadow = list(enabled = TRUE, size = 5)
    ) |>
    visNetwork::visEdges(
      smooth = FALSE,
      color = list(color = "rgba(150,150,150,0.6)", highlight = "#ff4b4b"),
      width = 2
    ) |>
    visNetwork::visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      selectedBy = list(variable = "group", main = "Select Group"),
      nodesIdSelection = if (isTRUE(search_by_name)) {
        list(enabled = TRUE, main = "Select Node")
      } else {
        FALSE
      }
    ) |>
    visNetwork::visPhysics(
      solver = "forceAtlas2Based",
      forceAtlas2Based = list(gravitationalConstant = -50, centralGravity = 0.01),
      stabilization = list(enabled = TRUE, iterations = 1000)
    ) |>
    visNetwork::visEvents(
      stabilizationIterationsDone = "function () {this.setOptions( { physics: false } );}"
    )

  # Optional second dropdown: step through propagation levels (see
  # .nts_network_add_iteration_selector); a no-op unless the column carries >= 2
  # distinct non-NA levels.
  p <- .nts_network_add_iteration_selector(p, nodes_dt, iteration_col)

  p
}

#' @rdname plot_molecular_network
#' @aliases plot_nts_network
#' @param project An [NTSProject] or `list(nodes, edges)` from
#'   [generate_network_data] / [annotate_by_dual_network]; when an
#'   [NTSProject], network tables are built with [generate_network_data] using
#'   `ms2_edge_types` and tolerances below. Recommended entry for dual-network
#'   workflows: pass an annotated `list(nodes, edges)` to preserve
#'   `Propagation_Level` and related columns.
#' @inheritParams generate_network_data
#' @param focus_nodes,max_steps Passed to [extract_subnetwork] when set.
#' @param show_iteration Passed to `.nts_network_apply_show_iteration`.
#' @export
plot_nts_network <- function(project,
                              ms2_edge_types = c("fragment", "neutral_loss", "mass_difference"),
                              ms1_edge_types = c("homologue", "transformation", "isf"),
                              mz_tol_ppm = 10,
                              min_shared_ms2 = 2,
                              max_bg_freq = 0.20,
                              min_frag_freq = 0.5,
                              max_frag_freq = NULL,
                              isf_mz_tol_ppm = 20,
                              md_min_da = 1.0,
                              min_ms2_score = NULL,
                              top_k = Inf,
                              max_edges = Inf,
                              min_degree = 0L,
                              cf2_ladder_boost = TRUE,
                              include_semantic = FALSE,
                              min_semantic_score = 0.70,
                              semantic_top_k = 3L,
                              embeddings = NULL,
                              python_path = NULL,
                              focus_nodes = NULL,
                              max_steps = 3L,
                              remove_isolated = TRUE,
                              color_by_edge = TRUE,
                              show_iteration = "all",
                              search_by_name = TRUE,
                              ...) {
  use_focus <- !is.null(focus_nodes)
  if (use_focus) {
    .nts_check(NULL, "igraph_ns", .arg = "igraph")
  }
  if (!is.null(max_frag_freq)) {
    warning("`max_frag_freq` has been renamed to `min_frag_freq` ",
            "(it keeps consensus peaks with frequency >= the value).", call. = FALSE)
    min_frag_freq <- max_frag_freq
  }
  if (inherits(project, "NTSProject")) {
    net_data <- generate_network_data(
      project,
      ms2_edge_types = ms2_edge_types,
      ms1_edge_types = ms1_edge_types,
      mz_tol_ppm = mz_tol_ppm,
      min_shared_ms2 = min_shared_ms2,
      max_bg_freq = max_bg_freq,
      min_frag_freq = min_frag_freq,
      isf_mz_tol_ppm = isf_mz_tol_ppm,
      md_min_da = md_min_da,
      min_ms2_score = min_ms2_score,
      top_k = top_k,
      max_edges = max_edges,
      min_degree = min_degree,
      cf2_ladder_boost = cf2_ladder_boost,
      include_semantic = include_semantic,
      min_semantic_score = min_semantic_score,
      semantic_top_k = semantic_top_k,
      embeddings = embeddings,
      python_path = python_path
    )
    if (use_focus) {
      sn <- extract_subnetwork(net_data, focus_nodes, max_steps = max_steps)
      net_data <- list(nodes = as.data.frame(sn$nodes), edges = as.data.frame(sn$edges))
    }
  } else if (is.list(project) && all(c("nodes", "edges") %in% names(project))) {
    net_data <- list(
      nodes = as.data.frame(project$nodes, stringsAsFactors = FALSE),
      edges = as.data.frame(project$edges, stringsAsFactors = FALSE)
    )
    if (use_focus) {
      sn <- extract_subnetwork(net_data, focus_nodes, max_steps = max_steps)
      net_data <- list(nodes = as.data.frame(sn$nodes), edges = as.data.frame(sn$edges))
    }
  } else {
    stop(
      "plot_nts_network: `project` must be an NTSProject or list(nodes, edges).",
      call. = FALSE
    )
  }
  net_data <- .nts_network_apply_show_iteration(net_data, show_iteration)
  if (!nrow(net_data$nodes)) {
    if (use_focus) {
      stop("Target node(s) not found or subnetwork is empty.", call. = FALSE)
    }
    message("> plot_molecular_network: no nodes after filtering - returning NULL.")
    return(invisible(NULL))
  }
  # Edge styling now lives in plot_molecular_network (P1-1); do not style here to
  # avoid applying .nts_edge_network_vis twice.
  plot_molecular_network(
    net_data,
    color_by_edge = color_by_edge,
    remove_isolated = remove_isolated,
    search_by_name = search_by_name,
    ...
  )
}

#' Export batch MS2 mirror plots to PDF
#'
#' Iterates through scored annotations and writes one mirror plot per feature to
#' a multi-page PDF.
#'
#' @param project An [NTSProject] with populated `@db_annotations`.
#' @param database Reference library: a flat [NTSDatabase] (peaks in
#'   `@spectra_long`) or a legacy [databaseClass] (nested `@spectra.data`).
#' @param output_pdf Output PDF path.
#' @param top_n Optional integer limiting export to the top-N scored features.
#' @param ion_mode Which library slot to use. `"auto"` tries negative first.
#' @param min_rel_intensity Noise floor as a percentage of the base peak
#'   (default `2.5`, matching the scoring engine's denoising); passed through to
#'   [plot_ms2_mirror]. Pass `0` to keep every peak.
#' @param label_matched_mz,label_mz_tol,label_digits Passed through to
#'   [plot_ms2_mirror]: label each matched fragment ion with its m/z
#'   (`label_matched_mz`, default `TRUE`), the matching tolerance in Da
#'   (`label_mz_tol`, default `0.02`) and the printed decimal places
#'   (`label_digits`, default `2`).
#'
#' @return Invisibly returns the normalized PDF path.
#' @export
export_mirror_report <- function(project,
                                 database,
                                 output_pdf = "mirror_plots.pdf",
                                 top_n = NULL,
                                 ion_mode = c("auto", "negative", "positive"),
                                 min_rel_intensity = 2.5,
                                 label_matched_mz = TRUE,
                                 label_mz_tol = 0.02,
                                 label_digits = 2) {
  stopifnot(inherits(project, "NTSProject"))
  ion_mode <- match.arg(ion_mode)
  ann <- methods::slot(project, "db_annotations")
  if (!nrow(ann) || !"name" %in% names(ann) || !"score_total" %in% names(ann)) {
    stop("db_annotations is empty or lacks required columns", call. = FALSE)
  }
  ann <- ann[is.finite(ann$score_total) & ann$score_total > 0, , drop = FALSE]
  if (!nrow(ann)) {
    stop("No positive-scoring annotations to export", call. = FALSE)
  }
  ann <- ann[order(ann$score_total, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  ann <- ann[!duplicated(as.character(ann$name)), , drop = FALSE]
  if (!is.null(top_n)) {
    top_n <- max(1L, as.integer(top_n)[1L])
    ann <- utils::head(ann, top_n)
  }

  grDevices::pdf(output_pdf, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  n_done <- 0L
  for (i in seq_len(nrow(ann))) {
    fname <- as.character(ann$name[[i]])
    p <- tryCatch(
      plot_ms2_mirror(project, database, fname, ion_mode = ion_mode,
                      min_rel_intensity = min_rel_intensity,
                      label_matched_mz = label_matched_mz,
                      label_mz_tol = label_mz_tol,
                      label_digits = label_digits),
      error = function(e) {
        warning("Skipping mirror plot for ", fname, ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(p)) next
    print(p)
    n_done <- n_done + 1L
  }

  message("> Exported ", n_done, " mirror plot(s) to: ", normalizePath(output_pdf, winslash = "/", mustWork = FALSE))
  invisible(normalizePath(output_pdf, winslash = "/", mustWork = FALSE))
}
