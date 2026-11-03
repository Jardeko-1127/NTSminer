# Source order matches AGENTS typical pipeline: [XcmsParam] to [extract_features];
# vendor [convert_mzml] sits last as optional preprocessing to produce mzML.

#' Parameters for xcms CentWave feature extraction
#'
#' `threads` is consumed by [get_robust_bpparam] (not passed to
#' [xcms::CentWaveParam]). Grouping parameters `binSize`, `bw`, and
#' `min_fraction` control [xcms::PeakDensityParam]; `obiwarp_binSize` sets the
#' **m/z** bin size for [xcms::ObiwarpParam] (see [extract_features]).
#'
#' @param threads Parallel worker count for [get_robust_bpparam].
#' @param ppm,peakwidth,noise,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,scanrange
#'   Passed to [xcms::CentWaveParam] when present. `mzdiff` controls
#'   minimum m/z difference for co-eluting peaks (`-0.001` forces separation).
#' @param snthresh,sn_thresh Signal-to-noise threshold for CentWave peak
#'   detection (default `10`). `sn_thresh` is an alias; `snthresh` takes
#'   precedence if both are supplied.
#' @param binSize m/z bin width for [xcms::PeakDensityParam] (default `0.008`,
#'   tuned for HRMS; too-large values cause over-merging of close m/z features).
#' @param bw RT bandwidth (seconds) for [xcms::PeakDensityParam] (default `20`).
#' @param min_fraction Minimum fraction of samples in which a grouped feature
#'   must be present (passed to [xcms::PeakDensityParam]; default `0.001`,
#'   matching NTS single-injection workflows). The same default is used by
#'   [suggest_xcms_params] and the shipped workflow templates.
#' @param obiwarp_binSize **m/z** bin size (Da) for [xcms::ObiwarpParam] in
#'   [xcms::adjustRtime] -- matching `xcms::ObiwarpParam`, this is in the m/z
#'   dimension, **not** seconds. `NULL` (default) uses `0.6` Da.
#'   [suggest_xcms_params] may return `obiwarp_binSize` (also in Da) for the
#'   same batch.
#' @param eic_metrics Logical, default `TRUE`. When `TRUE`, [extract_features]
#'   for an [NTSProject] derives a small set of **per-feature EIC-shape
#'   proxies** from [xcms::chromPeaks] (no extra raw-data reading) and writes
#'   them to `@ms1_annotations`. Columns added:
#'   \itemize{
#'     \item `apex_pos_ratio_median` - median of `(rt - rtmin) / (rtmax - rtmin)`
#'       across sample-level peaks (0.5 = symmetric; <0.5 leading apex; >0.5 tail).
#'     \item `into_maxo_ratio_median` - median of `into / maxo` (area/height).
#'     \item `peakwidth_median_s` - median per-sample peak width in seconds.
#'     \item `sn_median` - median CentWave signal-to-noise.
#'   }
#'   Set to `FALSE` to skip and keep `@ms1_annotations` lean.
#' @param peak_pick_serial If `TRUE`, run [xcms::findChromPeaks] with
#'   [BiocParallel::SerialParam] (lower **peak-pick** memory); if `FALSE` (default),
#'   use [get_robust_bpparam] with `threads` for that step.
#'   alignment and grouping still use the **registered** backend from
#'   `get_robust_bpparam` (in xcms **4.2.x**, [xcms::adjustRtime] and
#'   [xcms::groupChromPeaks] do not accept `BPPARAM`; verified error:
#'   `unused argument (BPPARAM)`). Per-step backends as in older tutorials thus
#'   cannot be applied after peak picking except via [xcms::findChromPeaks].
#' @param read_ms_cache Optional value for [MSnbase::readMSData] argument `cache.`
#'   (only applied when not overridden in `...`). On very large studies, consult
#'   `?MSnbase::readMSData` for cache levels. **Note:** tutorial snippets with
#'   `backend = MsBackendHdf5()` and `BPPARAM =` on `readMSData` depend on MSnbase/Spectra versions; pass them
#'   through `...` when your installation supports them.
#' @param feature_area_method How sample-column intensities are computed for
#'   `feature_table`: `"xcms_default"` ( [xcms::featureValues] `into`, default),
#'   `"eic_rt30"` (ppm-window max-per-scan EIC + RT trapezoid, MH-aligned), or
#'   `"apex_into"` (max `into` chromPeak per feature and sample). Non-default
#'   methods store the xcms matrix on attribute `xcms_default_intensity`.
#' @param feature_area_ppm m/z tolerance (ppm) for `feature_area_method = "eic_rt30"`
#'   (default `5`).
#' @param feature_area_rt_half RT half-window in seconds for `"eic_rt30"` (default `30`).
#' @param feature_area_rt_center RT center for `"eic_rt30"`: `"apex"` (default, per
#'   sample apex chromPeak) or `"rtmed"`.
#' @param feature_area_samples Optional sample names to recompute; `NULL` (default)
#'   means all samples in the xcms experiment.
#'
#' @details
#' **Construction checks**: `threads`, `ppm`, `peakwidth`, `snthresh`, `noise`,
#' `prefilter`, `integrate`, `mzdiff`, `scanrange`, `binSize`, `bw`,
#' `min_fraction`, and `obiwarp_binSize` (if non-`NULL`) are validated so invalid
#' values fail here with [stop] instead of inside [xcms::CentWaveParam].
#'
#' CentWave tuning: some examples use `CentWaveParam(maxScans = ...)`; **xcms
#' 4.2.x** `CentWaveParam` has no `maxScans` slot-raise `prefilter` second element
#' or tighten `prefilter` scans to reduce ROI work instead.
#'
#' @return An object of class `XcmsParam` (a tagged list).
#' @export
XcmsParam <- function(threads = 1L,
                      ppm = 10,
                      peakwidth = c(20, 50),
                      snthresh = NULL,
                      sn_thresh = 10,
                      noise = 0,
                      prefilter = c(3, 80),
                      mzCenterFun = "wMean",
                      integrate = 1L,
                      mzdiff = -0.001,
                      fitgauss = FALSE,
                      scanrange = integer(),
                      binSize = 0.008,
                      bw = 20,
                      min_fraction = 0.001,
                      obiwarp_binSize = NULL,
                      eic_metrics = TRUE,
                      peak_pick_serial = FALSE,
                      read_ms_cache = NULL,
                      feature_area_method = c("xcms_default", "eic_rt30", "apex_into"),
                      feature_area_ppm = 5,
                      feature_area_rt_half = 30,
                      feature_area_rt_center = c("apex", "rtmed"),
                      feature_area_samples = NULL) {
  snt <- if (!is.null(snthresh)) snthresh else sn_thresh
  threads <- suppressWarnings(as.integer(threads)[1L])
  if (!is.finite(threads) || threads < 1L) {
    stop("`threads` must be a positive integer.", call. = FALSE)
  }
  ppm <- suppressWarnings(as.numeric(ppm)[1L])
  if (!is.finite(ppm) || ppm <= 0) {
    stop("`ppm` must be a positive finite number.", call. = FALSE)
  }
  peakwidth <- suppressWarnings(as.numeric(peakwidth))
  if (length(peakwidth) != 2L || !all(is.finite(peakwidth)) ||
        peakwidth[1L] <= 0 || peakwidth[2L] < peakwidth[1L]) {
    stop(
      "`peakwidth` must be c(lower, upper) with positive finite bounds and ",
      "upper >= lower.",
      call. = FALSE
    )
  }
  snt <- suppressWarnings(as.numeric(snt)[1L])
  if (!is.finite(snt) || snt < 0) {
    stop("`snthresh` / `sn_thresh` must be a non-negative finite number.", call. = FALSE)
  }
  noise <- suppressWarnings(as.numeric(noise)[1L])
  if (!is.finite(noise) || noise < 0) {
    stop("`noise` must be a non-negative finite number.", call. = FALSE)
  }
  prefilter <- suppressWarnings(as.numeric(prefilter))
  if (length(prefilter) != 2L || !all(is.finite(prefilter)) ||
        any(prefilter < 0)) {
    stop("`prefilter` must be two non-negative finite numbers.", call. = FALSE)
  }
  mzCenterFun <- as.character(mzCenterFun)[1L]
  if (is.na(mzCenterFun) || !nzchar(mzCenterFun)) {
    stop("`mzCenterFun` must be a non-empty string.", call. = FALSE)
  }
  integrate <- suppressWarnings(as.integer(integrate)[1L])
  if (!is.finite(integrate) || integrate < 1L) {
    stop("`integrate` must be a positive integer (see ?xcms::CentWaveParam).", call. = FALSE)
  }
  mzdiff <- suppressWarnings(as.numeric(mzdiff)[1L])
  if (!is.finite(mzdiff)) {
    stop("`mzdiff` must be a finite number.", call. = FALSE)
  }
  fitgauss <- isTRUE(fitgauss)
  scanrange <- as.integer(scanrange)
  if (!length(scanrange)) {
    scanrange <- integer()
  } else if (length(scanrange) != 2L) {
    stop("`scanrange` must be integer() or a vector of two integers.", call. = FALSE)
  }
  binSize <- suppressWarnings(as.numeric(binSize)[1L])
  if (!is.finite(binSize) || binSize <= 0) {
    stop("`binSize` must be a positive finite number.", call. = FALSE)
  }
  bw <- suppressWarnings(as.numeric(bw)[1L])
  if (!is.finite(bw) || bw <= 0) {
    stop("`bw` must be a positive finite number.", call. = FALSE)
  }
  min_fraction <- suppressWarnings(as.numeric(min_fraction)[1L])
  if (!is.finite(min_fraction) || min_fraction <= 0 || min_fraction > 1) {
    stop("`min_fraction` must be in (0, 1].", call. = FALSE)
  }
  if (!is.null(obiwarp_binSize)) {
    obi <- suppressWarnings(as.numeric(obiwarp_binSize)[1L])
    if (!is.finite(obi) || obi <= 0) {
      stop("`obiwarp_binSize` must be NULL or a positive finite number.", call. = FALSE)
    }
    obiwarp_binSize <- obi
  }

  feature_area_method <- match.arg(feature_area_method)
  feature_area_ppm <- suppressWarnings(as.numeric(feature_area_ppm)[1L])
  if (!is.finite(feature_area_ppm) || feature_area_ppm <= 0) {
    stop("`feature_area_ppm` must be a positive finite number.", call. = FALSE)
  }
  feature_area_rt_half <- suppressWarnings(as.numeric(feature_area_rt_half)[1L])
  if (!is.finite(feature_area_rt_half) || feature_area_rt_half <= 0) {
    stop("`feature_area_rt_half` must be a positive finite number.", call. = FALSE)
  }
  feature_area_rt_center <- match.arg(feature_area_rt_center)
  if (!is.null(feature_area_samples)) {
    feature_area_samples <- as.character(feature_area_samples)
    if (!length(feature_area_samples)) {
      stop("`feature_area_samples` must be NULL or a non-empty character vector.", call. = FALSE)
    }
  }

  structure(
    list(
      threads = threads,
      ppm = ppm,
      peakwidth = peakwidth,
      snthresh = snt,
      noise = noise,
      prefilter = prefilter,
      mzCenterFun = mzCenterFun,
      integrate = integrate,
      mzdiff = mzdiff,
      fitgauss = fitgauss,
      scanrange = scanrange,
      binSize = binSize,
      bw = bw,
      min_fraction = min_fraction,
      obiwarp_binSize = obiwarp_binSize,
      eic_metrics = isTRUE(eic_metrics),
      peak_pick_serial = isTRUE(peak_pick_serial),
      read_ms_cache = read_ms_cache,
      feature_area_method = feature_area_method,
      feature_area_ppm = feature_area_ppm,
      feature_area_rt_half = feature_area_rt_half,
      feature_area_rt_center = feature_area_rt_center,
      feature_area_samples = feature_area_samples
    ),
    class = c("XcmsParam", "list")
  )
}


.centwave_from_xcmsparam <- function(object) {
  f <- names(formals(xcms::CentWaveParam))
  f <- f[f != "..."]
  pl <- as.list(object)
  pl <- pl[names(pl) %in% f]
  keep <- vapply(
    names(pl),
    function(nm) {
      z <- pl[[nm]]
      if (is.null(z)) {
        return(FALSE)
      }
      if (length(z) == 0L && (is.numeric(z) || is.integer(z))) {
        return(FALSE)
      }
      TRUE
    },
    logical(1L)
  )
  pl <- pl[keep]
  do.call(xcms::CentWaveParam, pl)
}


#' Extract chromatographic features (S3 generic)
#'
#' @param object [XcmsParam] **or** an [NTSProject]. For [NTSProject], pass
#'   `param = XcmsParam(...)` and optionally `files` (defaults to all mzML/mzXML
#'   under `@raw_data_path`).
#' @param files Character vector of MS data files; ignored for [NTSProject] if
#'   `NULL` (auto-discover).
#' @param msLevel MS level passed to [MSnbase::readMSData].
#' @param param,method `param` is an [XcmsParam] instance when `object` is
#'   [NTSProject]. `method` is accepted as an alias of `param` (PRD-style API).
#' @param ... Passed to [MSnbase::readMSData] (e.g. centroiding options, or
#'   version-specific `backend` / `cache.` where supported). If you omit `cache.`
#'   here, [XcmsParam] field `read_ms_cache` is forwarded as `cache.`.
#'
#' @return `data.frame` of peak metrics for [XcmsParam]; updated [NTSProject]
#'   with filled `@feature_table` and `@log` for project input.
#' @export
extract_features <- function(object,
                             files = NULL,
                             msLevel = 1L,
                             param = NULL,
                             method = NULL,
                             ...) {
  if (!is.null(method) && is.null(param)) {
    param <- method
  }
  if (inherits(object, "NTSProject")) {
    if (is.null(param) || !inherits(param, "XcmsParam")) {
      stop(
        "For NTSProject, supply `param = XcmsParam(...)` (or `method =`), e.g. ",
        "extract_features(project, param = XcmsParam(threads = 1, ppm = 10))",
        call. = FALSE
      )
    }
    if (is.null(files) || !length(files)) {
      files <- .nts_list_ms_files(methods::slot(object, "raw_data_path"))
    }
    if (!length(files)) {
      stop(
        "No mzML/mzXML under raw_data_path; pass `files` explicitly.",
        call. = FALSE
      )
    }
    peaks <- extract_features.XcmsParam(param, files = files, msLevel = msLevel, ...)
    xd <- attr(peaks, ".xdata")
    attr(peaks, ".xdata") <- NULL
    ft_core <- .nts_feature_table_core(peaks)
    object <- set_feature_table(object, ft_core)
    ms1_df <- .nts_init_ms1_from_ft(ft_core, bounds_df = peaks)
    if (isTRUE(param$eic_metrics)) {
      if (is.null(xd)) {
        message(
          "extract_features: no xcms result attached (.xdata); EIC-shape metrics skipped. ",
          "Re-run with the package extract_features implementation, or use refresh_eic_shape_metrics() ",
          "once `project@xdata` is available."
        )
      } else {
        ms1_df <- .nts_apply_eic_shape_to_ms1_df(ms1_df, xd)
      }
    }
    object <- set_ms1_annotations(object, ms1_df)
    methods::slot(object, "ms2data") <- list()
    object <- set_ms2_annotations(object, list())
    methods::slot(object, "clean_ms2") <- list()
    methods::slot(object, "db_annotations") <- data.frame()
    methods::slot(object, "physchem_proxy") <- data.frame()
    if (!is.null(xd)) {
      methods::slot(object, "xdata") <- xd
    }
    methods::slot(object, "log") <- .nts_log_append(
      methods::slot(object, "log"),
      "extract_features",
      list(param = param, n_files = length(files), msLevel = msLevel)
    )
    return(object)
  }
  UseMethod("extract_features", object)
}


#' Refresh EIC-shape columns on MS1 annotations
#'
#' Writes or updates the four CentWave-derived shape proxies (see
#' [XcmsParam] field `eic_metrics`) in `@ms1_annotations` by re-reading
#' [xcms::chromPeaks] from the `XCMSnExp` in `@xdata`.
#' Use when those columns are missing but `extract_features` already filled
#' `@xdata` (e.g. session resumed, or an older package version skipped the merge).
#'
#' @param project An [NTSProject].
#' @param xdata Optional xcms experiment object; defaults to `project@xdata`.
#' @return The updated `project`.
#' @export
refresh_eic_shape_metrics <- function(project, xdata = NULL) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  xd <- xdata
  if (is.null(xd)) xd <- methods::slot(project, "xdata")
  if (is.null(xd)) {
    stop(
      "No xcms experiment object: pass `xdata` or run extract_features(project, param = XcmsParam(...)) to fill `@xdata`.",
      call. = FALSE
    )
  }
  ms1 <- get_ms1_annotations(project)
  if (!nrow(ms1) || !"name" %in% names(ms1)) {
    stop("ms1_annotations is empty or has no `name` column.", call. = FALSE)
  }
  ms1 <- .nts_apply_eic_shape_to_ms1_df(ms1, xd)
  set_ms1_annotations(project, ms1)
}


#' @export
extract_features.default <- function(object,
                                     files,
                                     msLevel = 1L,
                                     param = NULL,
                                     method = NULL,
                                     ...) {
  stop(
    "No extract_features method for class: ",
    paste(class(object), collapse = "/"),
    call. = FALSE
  )
}


#' @rdname extract_features
#' @export
extract_features.XcmsParam <- function(object,
                                       files,
                                       msLevel = 1L,
                                       param = NULL,
                                       method = NULL,
                                       ...) {
  if (!length(files)) {
    stop("files must be a non-empty character vector", call. = FALSE)
  }

  bpp <- get_robust_bpparam(object$threads)
  old_bpp <- BiocParallel::bpparam()
  BiocParallel::register(bpp)
  on.exit(BiocParallel::register(old_bpp), add = TRUE)

  bpp_peak <- if (isTRUE(object$peak_pick_serial)) {
    BiocParallel::SerialParam()
  } else {
    bpp
  }

  cwp <- .centwave_from_xcmsparam(object)

  rd_dots <- list(...)
  if (!is.null(object$read_ms_cache) && !"cache." %in% names(rd_dots)) {
    rd_dots$cache. <- object$read_ms_cache
  }
  message("> Reading raw data (", length(files), " file(s))...")
  raw <- do.call(
    MSnbase::readMSData,
    c(list(files = files, mode = "onDisk", msLevel = msLevel), rd_dots)
  )

  message(
    "> Peak detection (CentWave, ppm=", object$ppm,
    ", snthresh=", object$snthresh, ", mzdiff=", object$mzdiff,
    if (isTRUE(object$peak_pick_serial)) "; BPPARAM=SerialParam" else "",
    ")..."
  )
  xdata <- xcms::findChromPeaks(raw, param = cwp, BPPARAM = bpp_peak)

  if (length(files) >= 2L) {
    obi_bs <- object$obiwarp_binSize
    if (is.null(obi_bs) || !is.finite(obi_bs) || obi_bs <= 0) obi_bs <- 0.6
    message(
      "> Retention-time alignment (Obiwarp, binSize=", obi_bs, ")..."
    )
    xdata <- xcms::adjustRtime(
      xdata,
      param = xcms::ObiwarpParam(binSize = obi_bs)
    )
  }

  min_frac <- object$min_fraction
  if (is.null(min_frac) || !is.finite(min_frac)) min_frac <- 0.001
  bin_sz <- object$binSize
  if (is.null(bin_sz) || !is.finite(bin_sz)) bin_sz <- 0.008
  bw_val <- object$bw
  if (is.null(bw_val) || !is.finite(bw_val)) bw_val <- 20
  message(
    "> Peak grouping (PeakDensity, binSize=", bin_sz,
    ", bw=", bw_val, ", minFraction=", min_frac, ")..."
  )
  xdata <- xcms::groupChromPeaks(
    xdata,
    param = xcms::PeakDensityParam(
      sampleGroups = rep(1L, length(files)),
      minFraction = min_frac,
      binSize = bin_sz,
      bw = bw_val
    )
  )

  message("> Gap filling (missing peaks)...")
  xdata <- xcms::fillChromPeaks(xdata)

  message("> Building feature matrix...")
  fd_obj <- xcms::featureDefinitions(xdata)
  rid <- BiocGenerics::rownames(fd_obj)
  if (is.null(rid) || !length(rid)) {
    rid <- sprintf("FT%06d", seq_len(nrow(fd_obj)))
  }
  # xcms 4.2.x: `featureDefinitions()` is often a **DFrame**; `as.data.frame(fd_obj)`
  # can corrupt mzmed/rtmed into nested data.frames. Index columns with [[.
  fd_cols <- colnames(fd_obj)
  stopifnot(all(c("mzmed", "rtmed") %in% fd_cols))
  fd <- data.frame(
    name = rid,
    mzmed = as.numeric(fd_obj[["mzmed"]]),
    mzmin = as.numeric(fd_obj[["mzmin"]]),
    mzmax = as.numeric(fd_obj[["mzmax"]]),
    rtmed = as.numeric(fd_obj[["rtmed"]]),
    rtmin = as.numeric(fd_obj[["rtmin"]]),
    rtmax = as.numeric(fd_obj[["rtmax"]]),
    npeaks = as.integer(fd_obj[["npeaks"]]),
    stringsAsFactors = FALSE
  )
  extra_fd <- setdiff(fd_cols, colnames(fd))
  for (cl in extra_fd) {
    fd[[cl]] <- fd_obj[[cl]]
  }
  fd$mz <- fd$mzmed
  ## xcms reports rt in **seconds** (same order of magnitude as MS file `rtime()`).
  fd$rt <- fd$rtmed

  area_bundle <- .nts_resolve_feature_areas(xdata, object)
  fv_m <- area_bundle$matrix
  cn <- colnames(fv_m)
  if (is.null(cn) || anyNA(cn) || any(cn == "")) {
    cn <- basename(files)[seq_len(ncol(fv_m))]
    colnames(fv_m) <- cn
  }
  fv_df <- as.data.frame(fv_m, stringsAsFactors = FALSE, check.names = FALSE)

  out <- cbind(fd, fv_df)
  rownames(out) <- NULL

  if (!identical(area_bundle$method, "xcms_default")) {
    message(
      "> Feature areas: method=", area_bundle$method,
      " (ppm=", area_bundle$ppm,
      ", RT+/-", area_bundle$rt_half, " s, center=", area_bundle$rt_center, ")"
    )
  }
  out <- .nts_set_feature_area_attrs(out, area_bundle)

  message("> Feature extraction done: ", nrow(out), " features x ", ncol(fv_m), " samples")
  attr(out, ".xdata") <- xdata
  out
}


#' Batch convert vendor raw files to mzML via ProteoWizard `msconvert`
#'
#' Scans `input_dir` recursively for Thermo `.raw`, Waters `.raw`, Agilent `.d`
#' folders, and `.wiff` files, then runs centroiding peak picking per file.
#'
#' @param input_dir Root directory to search.
#' @param output_dir Directory for mzML outputs (created if needed).
#' @param mode Conversion mode: `"centroid"` adds default peak-picking
#'   (`peakPicking true 1-`) when `filters = NULL`; `"profile"` keeps profile
#'   spectra by default (no auto filters).
#' @param exe `msconvert` executable path. `NULL` (default) auto-detects it:
#'   the `NTS_MSCONVERT` environment variable, then `msconvert` on `PATH`, then
#'   a `PROTEOWIZARD_HOME` install and the usual ProteoWizard locations.
#' @param intensity_filter Optional list to build an msconvert threshold filter.
#'   Example: `list(type = "count", value = 200)` becomes
#'   `"threshold count 200 most-intense"`.
#' @param filters Character vector of ProteoWizard filter expressions passed as
#'   repeated `--filter` arguments. Set `NULL` to use defaults from `mode`; pass
#'   explicit empty `character(0)` to disable all filters.
#' @param stdout,stderr Passed to [base::system2].
#' @param project Optional [NTSProject]. When supplied, the conversion is
#'   recorded in its `@log` (one entry with the input/output dirs, mode,
#'   resolved `exe`, filters and the converted file list) and the updated
#'   project is returned invisibly instead of the file vector.
#'
#' @return Invisibly, a character vector of the input paths that were
#'   processed. When `project` is supplied, invisibly the updated
#'   [NTSProject] (with a new `convert_mzml` log entry) instead.
#' @export
convert_mzml <- function(input_dir,
                         output_dir,
                         mode = c("centroid", "profile"),
                         exe = NULL,
                         intensity_filter = NULL,
                         filters = NULL,
                         stdout = FALSE,
                         stderr = "",
                         project = NULL) {
  mode <- match.arg(mode)
  input_dir <- normalizePath(input_dir, mustWork = TRUE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(exe) || !nzchar(as.character(exe)[1L])) {
    pw_home <- Sys.getenv("PROTEOWIZARD_HOME", unset = "")
    exe <- .nts_resolve_executable(
      "msconvert",
      candidates = c(
        if (nzchar(pw_home)) file.path(pw_home, "msconvert.exe") else character(0),
        "C:/Program Files/ProteoWizard/msconvert.exe",
        "C:/Program Files (x86)/ProteoWizard/msconvert.exe"
      ),
      envvar = "NTS_MSCONVERT"
    )
    if (is.null(exe)) {
      stop(
        "Could not locate `msconvert`. Pass `exe =` explicitly, put ",
        "`msconvert` on your PATH, or set the `NTS_MSCONVERT` environment ",
        "variable to its full path.",
        call. = FALSE
      )
    }
  }
  exe <- as.character(exe)[1L]
  if (!nzchar(exe)) {
    stop("`exe` must be a non-empty msconvert executable path.", call. = FALSE)
  }
  if (grepl("[/\\\\]", exe) && !file.exists(exe)) {
    stop("`exe` does not exist: ", exe, call. = FALSE)
  }
  if (is.null(filters)) {
    filters <- if (identical(mode, "centroid")) {
      "peakPicking true 1-"
    } else {
      character(0)
    }
  }
  if (!is.null(intensity_filter)) {
    filters <- c(filters, .nts_msconvert_intensity_filter(intensity_filter))
  }
  filters <- as.character(filters)
  filters <- filters[nzchar(trimws(filters))]

  # Shared exit: optionally record the conversion in a project's @log, then
  # return either the file vector (default) or the updated project.
  finish <- function(files) {
    if (is.null(project)) {
      return(invisible(files))
    }
    .nts_check(project, "project", .arg = "project")
    methods::slot(project, "log") <- .nts_log_append(
      get_project_log(project), "convert_mzml",
      list(
        input_dir = input_dir, output_dir = output_dir, mode = mode,
        exe = exe, filters = filters, intensity_filter = intensity_filter,
        n_files = length(files), files = basename(files)
      )
    )
    invisible(project)
  }

  f_raw <- list.files(
    input_dir,
    pattern = "\\.(raw|wiff)$",
    ignore.case = TRUE,
    full.names = TRUE,
    recursive = TRUE
  )
  f_d <- list.dirs(input_dir, recursive = TRUE, full.names = TRUE)
  f_d <- f_d[grepl("\\.d$", f_d, ignore.case = TRUE)]
  files <- unique(c(sort(f_raw), sort(f_d)))
  n <- length(files)
  if (!n) {
    warning("No .raw, .wiff, or .d inputs found under input_dir.", call. = FALSE)
    return(finish(character(0)))
  }

  for (i in seq_len(n)) {
    f <- files[i]
    message(sprintf("[%d/%d] Converting %s", i, n, basename(f)))
    args <- c(
      shQuote(f),
      "-o", shQuote(output_dir),
      "--mzML"
    )
    if (length(filters)) {
      filter_args <- unlist(
        lapply(filters, function(fl) c("--filter", shQuote(fl))),
        use.names = FALSE
      )
      args <- c(args, filter_args)
    }
    exe_cmd <- if (grepl(" ", exe, fixed = TRUE)) shQuote(exe) else exe
    status <- system2(exe_cmd, args = args, stdout = stdout, stderr = stderr)
    if (!identical(as.integer(status), 0L)) {
      stop(
        "msconvert failed on ", f, " (exit code ", status, ").",
        call. = FALSE
      )
    }
  }
  finish(files)
}


#' @keywords internal
#' @noRd
.nts_msconvert_intensity_filter <- function(intensity_filter) {
  if (!is.list(intensity_filter)) {
    stop("`intensity_filter` must be a list, e.g. list(type='count', value=200).", call. = FALSE)
  }
  type <- if (!is.null(intensity_filter$type)) {
    tolower(trimws(as.character(intensity_filter$type)[1L]))
  } else {
    "count"
  }
  if (identical(type, "custom")) {
    expr <- if (!is.null(intensity_filter$expr)) {
      as.character(intensity_filter$expr)[1L]
    } else {
      ""
    }
    expr <- trimws(expr)
    if (!nzchar(expr)) {
      stop("`intensity_filter$type='custom'` requires non-empty `expr`.", call. = FALSE)
    }
    return(expr)
  }

  val <- suppressWarnings(as.numeric(intensity_filter$value)[1L])
  if (!is.finite(val) || val <= 0) {
    stop("`intensity_filter$value` must be a positive number.", call. = FALSE)
  }
  orientation <- if (!is.null(intensity_filter$orientation)) {
    tolower(trimws(as.character(intensity_filter$orientation)[1L]))
  } else {
    "most-intense"
  }
  if (!orientation %in% c("most-intense", "least-intense")) {
    stop("`intensity_filter$orientation` must be 'most-intense' or 'least-intense'.", call. = FALSE)
  }

  if (identical(type, "count")) {
    return(sprintf("threshold count %d %s", as.integer(round(val)), orientation))
  }
  if (identical(type, "absolute")) {
    v <- format(val, scientific = FALSE, trim = TRUE)
    return(sprintf("threshold absolute %s %s", v, orientation))
  }

  stop(
    "`intensity_filter$type` must be one of: 'count', 'absolute', 'custom'.",
    call. = FALSE
  )
}


#' @keywords internal
#' @noRd
.nts_apply_eic_shape_to_ms1_df <- function(ms1_df, xd) {
  if (!is.data.frame(ms1_df) || !nrow(ms1_df) || !"name" %in% names(ms1_df)) {
    return(ms1_df)
  }
  if (is.null(xd)) {
    return(.nts_ms1_sanitize_and_order(ms1_df))
  }
  shape_df <- tryCatch(
    .nts_compute_eic_shape_metrics(xd),
    error = function(e) {
      warning(
        "EIC-shape metrics failed; leaving columns unchanged (",
        conditionMessage(e), ")",
        call. = FALSE
      )
      data.frame()
    }
  )
  if (!nrow(shape_df) || !"name" %in% names(shape_df)) {
    return(.nts_ms1_sanitize_and_order(ms1_df))
  }
  idx <- match(as.character(ms1_df$name), as.character(shape_df$name))
  for (cl in setdiff(names(shape_df), "name")) {
    ms1_df[[cl]] <- shape_df[[cl]][idx]
  }
  .nts_ms1_sanitize_and_order(ms1_df)
}


# Per-feature EIC-shape proxies derived from [xcms::chromPeaks] (CentWave-side
# scalars; no extra raw-data IO). Aggregated across all sample-level peaks
# belonging to a feature via `featureDefinitions(xdata)$peakidx`.
#' @keywords internal
#' @noRd
.nts_compute_eic_shape_metrics <- function(xdata) {
  if (is.null(xdata)) return(data.frame())
  if (!requireNamespace("xcms", quietly = TRUE)) return(data.frame())

  fd <- tryCatch(xcms::featureDefinitions(xdata), error = function(e) NULL)
  if (is.null(fd) || !nrow(fd)) return(data.frame())

  pks <- tryCatch(xcms::chromPeaks(xdata), error = function(e) NULL)
  rid <- BiocGenerics::rownames(fd)
  if (is.null(rid) || !length(rid)) {
    rid <- sprintf("FT%06d", seq_len(nrow(fd)))
  }
  if (is.null(pks) || !nrow(pks)) {
    out <- data.frame(name = rid, stringsAsFactors = FALSE)
    for (nm in .nts_ms1_eic_shape_colnames()) {
      out[[nm]] <- rep(NA_real_, length(rid))
    }
    return(out)
  }

  pcols <- colnames(pks)
  has_rt    <- all(c("rt", "rtmin", "rtmax") %in% pcols)
  has_ratio <- all(c("into", "maxo") %in% pcols)
  has_sn    <- "sn" %in% pcols

  n_feat <- nrow(fd)
  ridx_list <- tryCatch(fd[["peakidx"]], error = function(e) NULL)
  if (is.null(ridx_list) || length(ridx_list) != n_feat) {
    if (is.null(ridx_list)) {
      warning(
        "featureDefinitions(xdata)$peakidx is NULL; ",
        "EIC-shape proxies will be NA. (xcms version / pipeline state.)",
        call. = FALSE
      )
    } else {
      warning(
        "featureDefinitions(xdata)$peakidx length (", length(ridx_list),
        ") != nrow(features) (", n_feat, "); EIC-shape proxies will be NA.",
        call. = FALSE
      )
    }
    ridx_list <- NULL
  }

  apex <- rep(NA_real_, n_feat)
  ar   <- rep(NA_real_, n_feat)
  pw   <- rep(NA_real_, n_feat)
  snm  <- rep(NA_real_, n_feat)

  for (i in seq_len(n_feat)) {
    idx <- integer(0)
    if (!is.null(ridx_list)) {
      idx <- suppressWarnings(as.integer(ridx_list[[i]]))
    }
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= nrow(pks)]
    if (!length(idx)) next
    sub <- pks[idx, , drop = FALSE]

    if (has_rt) {
      span <- as.numeric(sub[, "rtmax"]) - as.numeric(sub[, "rtmin"])
      span[!is.finite(span) | span <= 0] <- NA_real_
      ratio <- (as.numeric(sub[, "rt"]) - as.numeric(sub[, "rtmin"])) / span
      apex[i] <- stats::median(ratio, na.rm = TRUE)
      pw[i]   <- stats::median(span, na.rm = TRUE)
    }
    if (has_ratio) {
      maxo <- as.numeric(sub[, "maxo"])
      maxo[!is.finite(maxo) | maxo <= 0] <- NA_real_
      ar[i] <- stats::median(as.numeric(sub[, "into"]) / maxo, na.rm = TRUE)
    }
    if (has_sn) {
      snm[i] <- stats::median(
        as.numeric(sub[, "sn"])[is.finite(as.numeric(sub[, "sn"]))],
        na.rm = TRUE
      )
    }
  }

  apex[!is.finite(apex)] <- NA_real_
  ar[!is.finite(ar)]     <- NA_real_
  pw[!is.finite(pw)]     <- NA_real_
  snm[!is.finite(snm)]   <- NA_real_

  eic_c <- .nts_ms1_eic_shape_colnames()
  cbind(
    data.frame(name = rid, stringsAsFactors = FALSE),
    stats::setNames(
      data.frame(apex, ar, pw, snm, stringsAsFactors = FALSE),
      eic_c
    )
  )
}
