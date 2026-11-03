# Parse `function=N` from the mzML spectrumId via mzR.
# The result is restricted to MS1 spectra so it aligns element-for-element with
# `MSnbase::readMSData(msLevel = 1L)` (used by both Sonar and the BO refiner).
# Returns an integer vector (NA where absent), or NULL on any failure (so
# callers can gracefully fall back to DDA behaviour).
#' @noRd
.sonar_try_parse_function <- function(file) {
  if (!requireNamespace("mzR", quietly = TRUE)) return(NULL)
  m <- try(mzR::openMSfile(file), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  hdr <- try(mzR::header(m), silent = TRUE)
  try(mzR::close(m), silent = TRUE)
  if (inherits(hdr, "try-error") || !is.data.frame(hdr) || !nrow(hdr) ||
      is.null(hdr$spectrumId)) {
    return(NULL)
  }
  lvl <- if (is.null(hdr$msLevel)) rep(1L, nrow(hdr)) else hdr$msLevel
  sel <- which(lvl == 1L)
  if (!length(sel)) return(NULL)
  fn <- suppressWarnings(as.integer(
    sub(".*function=([0-9]+).*", "\\1", as.character(hdr$spectrumId[sel]))
  ))
  fn[!is.finite(fn)] <- NA_integer_
  fn
}


# Paramounter-style weighted moving-average smoothing of an EIC segment.
# `level` = 0 disables smoothing. Fallbacks keep short / constant inputs
# unchanged rather than erroring.
#' @noRd
.sonar_peak_smooth <- function(x, level = 0) {
  n <- as.integer(level)[1L]
  if (is.na(n) || n <= 0L || length(x) < 2L * n) return(x)
  if (length(unique(x)) == 1L) return(x)

  weights <- c(seq_len(n + 1L), n:1L)
  y <- numeric(length(x))
  for (i in seq_along(x)) {
    start <- max(1L, i - n)
    end <- min(length(x), i + n)
    w_start <- n + 1L - (i - start)
    w_end <- n + 1L + (end - i)
    curr_w <- weights[w_start:w_end]
    y[i] <- sum(x[start:end] * curr_w) / sum(curr_w)
  }
  y
}


# Paramounter-style knee search for the background noise floor.
# Ascends the non-zero intensities and locates the first point where the
# running mean + 3*sd of the lower bulk is exceeded -- the "knee" that separates
# background from real signal. If no knee is found (uniform / tightly-clustered /
# near-constant / very short vectors), fall back to the median of the non-zero
# values instead of `max(nz)`: returning the maximum would drive the global
# pre-filter `dt[int > global_noise]` to retain 0% of points and Sonar would
# later abort with a misleading "no peaks" error. The median fallback guarantees
# >= 50% retention.
#' @noRd
.sonar_estimate_noise <- function(int_vec) {
  nz <- sort(int_vec[int_vec > 0])
  n <- length(nz)
  if (n == 0L) return(0)
  if (n < 10L) return(stats::median(nz))

  cutoff <- NA_real_
  if (n > 10L) {
    # `x` ascends in steps of 10; stop before `x == n` so `nz[x + 1L]` is valid.
    for (x in seq.int(10L, n - 1L, by = 10L)) {
      s_dev <- stats::sd(nz[seq_len(x)])
      thres <- mean(nz[seq_len(x)]) + 3 * s_dev
      if (nz[x + 1L] >= thres) { cutoff <- nz[x]; break }
    }
  }
  if (!is.finite(cutoff)) cutoff <- stats::median(nz)
  cutoff
}


# Turn the per-ROI statistics table (`res`, one row per detected chromatographic
# peak: columns mz, ppm, mz_span, pw, scans, height, noise, sn) into the
# suggested XCMS parameter list. Extracted from `suggest_xcms_params` so the
# quantile / clamping logic can be unit-tested without real mzML files.
#
# Robustness notes:
#   * `ppm` uses a Tukey IQR upper fence before the quantile. Without it, a
#     handful of pathological ROIs (overlapped bins, overloaded standard peaks)
#     inflate the upper quantile -- e.g. a clean PFAS standard shows a raw 95th
#     percentile of ~150 ppm against a median of ~5 ppm, so the value was always
#     pinned to the `max_suggested_ppm` ceiling. The fence restores a realistic
#     value while keeping the high-quantile design intent.
#   * every quantile is guarded for an empty / all-NA input so it cannot raise
#     "missing value where TRUE/FALSE needed" on degenerate ROIs.
#   * `prefilter[2]` is anchored to the detection floor (`noise`) rather than to
#     peak heights: a height quantile is a peak-height statistic and
#     systematically over-estimates the per-scan intensity floor.
#' @noRd
.sonar_derive_params <- function(res, cutoff_quantile = 0.95,
                                 max_suggested_ppm = 80L,
                                 mzdiff_range = c(-0.001, 0.01),
                                 prefilter_scans = c(2L, 4L),
                                 prefilter_intensity = c(50, 1e5),
                                 min_fraction_default = 0.001) {
  psr <- as.integer(prefilter_scans)
  pir <- as.numeric(prefilter_intensity)
  roi_ppm_cap <- 500

  # --- ppm ---------------------------------------------------------------
  res_clean <- res[is.finite(res$ppm) & res$ppm > 0]
  ppm_pool <- res_clean$ppm[res_clean$ppm <= roi_ppm_cap]
  if (length(ppm_pool) < 20L) ppm_pool <- res_clean$ppm
  if (length(ppm_pool) >= 5L) {
    pq <- stats::quantile(ppm_pool, c(0.25, 0.75), na.rm = TRUE)
    fence <- as.numeric(pq[2L]) + 1.5 * (as.numeric(pq[2L]) - as.numeric(pq[1L]))
    trimmed <- ppm_pool[ppm_pool <= fence]
    if (length(trimmed) >= 5L) ppm_pool <- trimmed
  }
  suggested_ppm <- if (length(ppm_pool)) {
    ceiling(stats::quantile(ppm_pool, cutoff_quantile, na.rm = TRUE))
  } else NA_real_
  if (!is.finite(suggested_ppm) || suggested_ppm < 5) suggested_ppm <- 5
  if (!is.null(max_suggested_ppm)) {
    cap <- suppressWarnings(as.integer(max_suggested_ppm)[1L])
    if (is.finite(cap) && cap >= 5L && suggested_ppm > cap) suggested_ppm <- cap
  }

  # --- peakwidth (FWHM) --------------------------------------------------
  pw_clean <- res$pw[is.finite(res$pw) & res$pw > 0 & res$pw < 500]
  if (length(pw_clean) >= 5L) {
    suggested_min_pw <- floor(stats::quantile(pw_clean, 0.05, na.rm = TRUE))
    suggested_max_pw <- ceiling(stats::quantile(pw_clean, cutoff_quantile, na.rm = TRUE))
  } else {
    suggested_min_pw <- NA_real_
    suggested_max_pw <- NA_real_
  }
  if (!is.finite(suggested_min_pw) || suggested_min_pw < 3) suggested_min_pw <- 3
  if (!is.finite(suggested_max_pw) || suggested_max_pw <= suggested_min_pw) {
    suggested_max_pw <- 30
  }

  # --- noise & sn --------------------------------------------------------
  noise_ok <- res$noise[is.finite(res$noise)]
  suggested_noise <- if (length(noise_ok)) {
    floor(stats::quantile(noise_ok, 1 - cutoff_quantile, na.rm = TRUE))
  } else NA_real_
  if (!is.finite(suggested_noise) || suggested_noise < 100) suggested_noise <- 100

  sn_ok <- res$sn[is.finite(res$sn) & res$sn > 0]
  suggested_sn <- if (length(sn_ok)) {
    floor(stats::quantile(sn_ok, 1 - cutoff_quantile, na.rm = TRUE))
  } else NA_real_
  if (!is.finite(suggested_sn) || suggested_sn < 3) suggested_sn <- 3

  # --- prefilter ---------------------------------------------------------
  suggested_min_scan <- floor(stats::quantile(res$scans, 0.05, na.rm = TRUE))
  if (!is.finite(suggested_min_scan)) suggested_min_scan <- psr[1L]
  suggested_min_scan <- as.integer(max(psr[1L], min(psr[2L], suggested_min_scan)))
  int_floor <- stats::quantile(res$noise, 0.05, na.rm = TRUE)
  suggested_min_int <- if (is.finite(int_floor)) floor(int_floor) else pir[1L]
  suggested_min_int <- as.integer(max(pir[1L], min(pir[2L], suggested_min_int)))

  # --- mzdiff ------------------------------------------------------------
  mz_span_q <- res$mz_span[is.finite(res$mz_span) & res$mz_span >= 0]
  if (length(mz_span_q) >= 5L) {
    q50 <- stats::quantile(mz_span_q, 0.5, na.rm = TRUE)
    suggested_mzdiff <- 0.5 * as.numeric(q50)
  } else {
    suggested_mzdiff <- mzdiff_range[1L]
  }
  suggested_mzdiff <- max(mzdiff_range[1L], min(mzdiff_range[2L], suggested_mzdiff))
  if (!is.finite(suggested_mzdiff) || abs(suggested_mzdiff) < 1e-6) {
    suggested_mzdiff <- mzdiff_range[1L]
  }
  suggested_mzdiff <- round(suggested_mzdiff, 4L)

  # --- peakdensity / obiwarp --------------------------------------------
  mz_med <- stats::median(res$mz[is.finite(res$mz)], na.rm = TRUE)
  if (!is.finite(mz_med) || mz_med <= 0) mz_med <- 300
  suggested_bin <- 0.5 * mz_med * (suggested_ppm / 1e6)
  suggested_bin <- max(0.003, min(0.05, suggested_bin))

  med_pw <- stats::median(pw_clean, na.rm = TRUE)
  if (!is.finite(med_pw)) med_pw <- 10
  suggested_bw <- max(10, min(30, med_pw * 3))
  # ObiwarpParam@binSize is in the m/z dimension (Da); it is the profile-matrix
  # bin and cannot be inferred from MS1 EICs, so return the shipped safe default.
  suggested_obiwarp_bin <- 0.6

  list(
    ppm = suggested_ppm,
    peakwidth = c(suggested_min_pw, suggested_max_pw),
    snthresh = suggested_sn,
    noise = suggested_noise,
    prefilter = c(suggested_min_scan, suggested_min_int),
    mzdiff = suggested_mzdiff,
    min_fraction = min_fraction_default,
    binSize = suggested_bin,
    bw = suggested_bw,
    obiwarp_binSize = suggested_obiwarp_bin
  )
}


#' Data-Driven XCMS Parameter Optimization Engine (Sonar)
#'
#' This function implements a modernized, high-performance version of the 
#' Paramounter algorithm (Guo & Huan, 2021). It performs an automated 
#' "sonar scan" of raw MS1 data to detect optimal parameters for XCMS 
#' CentWave peak picking (ppm, peakwidth, snthresh, noise, prefilter, mzdiff) 
#' by analyzing the statistical properties of Ion Chromatograms (EICs).
#' It also returns **light-touch** suggestions for \code{PeakDensityParam}
#' (\code{min_fraction}, \code{binSize}, \code{bw}) and \code{ObiwarpParam}
#' (\code{obiwarp_binSize}) so per-batch QC files can drive a single Sonar pass
#' without full Bayesian optimization.
#'
#' @param ms1_files Character vector of paths to raw MS1 files (mzXML, mzML).
#'   Typical use: 3--4 QC files for the current batch.
#' @param mode One of \code{"auto"} / \code{"dda"} / \code{"dia"} / \code{"gc"}.
#'   \code{"dda"} (current behaviour) feeds all MS1 scans to Sonar.
#'   \code{"dia"} targets Waters MSE mzML where low/high collision energies are
#'   both written under \code{msLevel = 1} and distinguished only by
#'   \code{function=N} in \code{spectrumId}; in this mode Sonar uses only
#'   \code{function=1} (low-energy) scans, avoiding fragment-induced inflation
#'   of \code{ppm}/\code{peakwidth}. \code{"auto"} (default) probes
#'   \code{spectrumId} via \pkg{mzR}: dia if the field is present, dda otherwise.
#'   \code{"gc"} is **deprecated**: use [detect_gc_compounds] instead for
#'   erah-style full-chromatogram deconvolution that replaces CentWave
#'   peak picking for GC-EI data.
#' @param sample_n Number of files to sample for detection. \code{NULL} (the
#'   default) processes all files. Supplied values are clamped to
#'   \code{[1, length(ms1_files)]} so that \code{sample_n = NULL} or an
#'   oversized value can never truncate the run.
#' @param cutoff_quantile Statistical quantile for parameter detection (default 0.95).
#' @param mass_sd_range Number of standard deviations for mass stability (default 2).
#' @param smooth Smooth level for EIC (0 = no smoothing).
#' @param roi_gap_da Maximum m/z gap (Da) used to chain points into one ROI
#'   before profiling; default \code{0.02} (see S-5 note in the evaluation
#'   report). Lower it if distinct species are being merged.
#' @param max_suggested_ppm Hard ceiling on suggested CentWave \code{ppm}. Sonar
#'   derives ppm from per-ROI m/z scatter; outliers (common with Waters MSE
#'   mzML in \code{"dda"} fallback mode) can drive quantiles to absurd values
#'   (e.g. 10^4 ppm) without this cap. Default \code{80} is reasonable for HRMS
#'   under \code{"dda"}; under \code{"dia"} you can safely lower it to ~15--25.
#'   Use \code{NULL} for no ceiling (not recommended).
#' @param mzdiff_range Numeric vector of length 2, \code{c(lower, upper)} for
#'   suggested CentWave \code{mzdiff} (Da). Default \code{c(-0.001, 0.01)} matches
#'   typical bounds used in grid/BO tuning scripts.
#' @param prefilter_scans_range Integerish length-2 vector: \code{c(min, max)} scan
#'   count for \code{prefilter[1]}. Defaults \code{c(2L, 4L)} (same window as typical
#'   CentWave BO scripts); pass \code{c(2L, 8L)} for noisier environmental matrices.
#' @param prefilter_intensity_range Numeric length-2 vector: \code{c(min, max)} for
#'   \code{prefilter[2]} (minimum intensity). Default \code{c(50, 1e5)}; floor
#'   follows common practice (see \code{XcmsParam} default 100, BO examples 50+).
#' @param min_fraction_default Numeric in \code{[0, 1]} returned as the
#'   \code{min_fraction} suggestion. Default \code{0.001} matches
#'   \code{\link{XcmsParam}} and the shipped workflow templates (NTS workflows
#'   with few or no replicates). **Sonar cannot infer this from QC files**---if
#'   you have replicate groups, override it in \code{\link{XcmsParam}}.
#' @param seed Integer seed used to draw the \code{sample_n} files, so the
#'   Sonar suggestion is reproducible (default \code{123L}).
#'
#' @return A named list: CentWave fields \code{ppm}, \code{peakwidth},
#'   \code{snthresh}, \code{noise}, \code{prefilter}, \code{mzdiff}; plus
#'   \code{min_fraction}, \code{binSize}, \code{bw} for \code{PeakDensityParam},
#'   and \code{obiwarp_binSize} for \code{ObiwarpParam}. Note that
#'   \code{obiwarp_binSize} is expressed in the \strong{m/z dimension (Da)},
#'   matching \code{xcms::ObiwarpParam}, \strong{not} in seconds; it is a safe
#'   default (Sonar cannot infer the profile-matrix bin from MS1 EICs), much
#'   like \code{min_fraction}. All fields can be passed into
#'   \code{\link{XcmsParam}}.
#' @export
suggest_xcms_params <- function(ms1_files, 
                                mode = c("auto", "dda", "dia", "gc"),
                                sample_n = NULL, 
                                cutoff_quantile = 0.95, 
                                mass_sd_range = 2,
                                smooth = 0,
                                roi_gap_da = 0.02,
                                max_suggested_ppm = 80L,
                                mzdiff_range = c(-0.001, 0.01),
                                prefilter_scans_range = c(2L, 4L),
                                prefilter_intensity_range = c(50, 1e5),
                                min_fraction_default = 0.001,
                                seed = 123L) {
  
  if (!requireNamespace("xcms", quietly = TRUE)) stop("xcms is required")
  if (!requireNamespace("MSnbase", quietly = TRUE)) stop("MSnbase is required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required")
  
  mode <- match.arg(mode)

  if (mode == "gc") {
    centwave_params <- list(
      ppm          = 20,
      peakwidth    = c(1, 8),
      snthresh     = 5,
      prefilter    = c(3, 100),
      mzCenterFun  = "wMean",
      integrate    = 1L,
      mzdiff       = -0.001,
      noise        = 100,
      fitgauss     = FALSE
    )
    out <- list(
      algorithm   = "centWave",
      centwave    = centwave_params,
      peakdensity = list(
        minFraction = 0.5,
        bw          = 2,
        binSize     = 0.01
      ),
      obiwarp = list(
        binSize = 0.6
      ),
      mode = "gc",
      notes = "GC-QTOF EI parameters: narrow peaks (1-5s), TOF mass accuracy (~20 ppm), no polarity filtering"
    )
    return(out)
  }

  min_fraction_default <- as.numeric(min_fraction_default[1L])
  if (!is.finite(min_fraction_default) || min_fraction_default < 0 ||
      min_fraction_default > 1) {
    stop("`min_fraction_default` must be a finite numeric in [0, 1]",
         call. = FALSE)
  }
  
  mzdiff_range <- sort(as.numeric(mzdiff_range)[1:2])
  if (length(mzdiff_range) != 2L || !all(is.finite(mzdiff_range))) {
    stop("`mzdiff_range` must be two finite numerics", call. = FALSE)
  }
  psr <- suppressWarnings(as.integer(round(prefilter_scans_range[1:2])))
  if (anyNA(psr) || min(psr) < 1L) {
    stop("`prefilter_scans_range` must be two positive integers", call. = FALSE)
  }
  if (psr[1] > psr[2]) psr <- sort(psr)
  pir <- sort(as.numeric(prefilter_intensity_range)[1:2])
  if (length(pir) != 2L || pir[1] < 0 || pir[1] > pir[2]) {
    stop("`prefilter_intensity_range` must be sorted non-negative numerics", call. = FALSE)
  }
  
  # --- 1. Preparation ---
  ms1_files <- path.expand(as.character(ms1_files))
  ms1_files <- ms1_files[file.exists(ms1_files)]
  if (length(ms1_files) == 0) stop("No valid MS1 files provided.")
  
  # `sample_n = NULL` is the documented default and means "process all files".
  # Guard against `length(x) > NULL` -> logical(0) -> "argument is of length zero".
  seed <- suppressWarnings(as.integer(seed)[1L])
  if (is.na(seed)) seed <- 123L
  roi_gap_da <- suppressWarnings(as.numeric(roi_gap_da)[1L])
  if (!is.finite(roi_gap_da) || roi_gap_da <= 0) roi_gap_da <- 0.02
  # Seed the whole run: the global-noise estimate below subsamples the intensity
  # vector with `sample()`, so without this the suggestion would differ between
  # identical calls.
  set.seed(seed)
  n_files <- length(ms1_files)
  if (is.null(sample_n)) {
    target_files <- ms1_files
  } else {
    sample_n <- suppressWarnings(as.integer(sample_n)[1L])
    if (is.na(sample_n) || sample_n < 1L) sample_n <- n_files
    if (sample_n >= n_files) {
      target_files <- ms1_files
    } else {
      target_files <- sample(ms1_files, sample_n)
    }
  }
  
  message(sprintf("Starting Sonar Scan on %d file(s)... (mode=%s)",
                  length(target_files), mode))
  
  all_stats <- list()
  
  for (i in seq_along(target_files)) {
    f <- target_files[i]
    
    # --- DIA branch: decide per-file whether to subset to function=1 ---
    use_dia_subset <- FALSE
    fn_vec <- NULL
    if (mode != "dda") {
      fn_vec <- .sonar_try_parse_function(f)
      has_fn <- !is.null(fn_vec) && any(!is.na(fn_vec))
      if (mode == "dia") {
        if (!has_fn) {
          stop(
            sprintf("mode='dia' but `function=N` not found in spectrumId of %s",
                    basename(f)),
            call. = FALSE
          )
        }
        use_dia_subset <- TRUE
      } else {
        use_dia_subset <- has_fn
      }
    }
    message(sprintf("[%d/%d] Analyzing: %s%s",
                    i, length(target_files), basename(f),
                    if (use_dia_subset) "  [DIA: function=1 only]" else ""))
    
    # 2. Fast Data Loading
    raw_data <- MSnbase::readMSData(f, msLevel = 1, mode = "onDisk")
    rts <- MSnbase::rtime(raw_data)
    mzs <- MSnbase::mz(raw_data)
    ints <- MSnbase::intensity(raw_data)
    
    if (use_dia_subset) {
      n_total <- length(rts)
      # MSnbase honours mzML spectrum order; align by length with mzR header.
      if (length(fn_vec) == n_total) {
        keep <- which(fn_vec == 1L & !is.na(fn_vec))
      } else {
        warning(sprintf(
          "%s: mzR header length (%d) != MSnbase scan count (%d); skip DIA subset",
          basename(f), length(fn_vec), n_total
        ), call. = FALSE)
        keep <- seq_len(n_total)
      }
      if (length(keep) < 5L) {
        warning(sprintf("%s: too few function=1 scans (%d); skip DIA subset",
                        basename(f), length(keep)), call. = FALSE)
      } else {
        rts <- rts[keep]
        mzs <- mzs[keep]
        ints <- ints[keep]
      }
    }
    
    # Flatten to long data.table for vector processing
    dt <- data.table::data.table(
      mz = unlist(mzs),
      int = unlist(ints),
      scan = rep(seq_along(mzs), lengths(mzs))
    )
    dt[, rt := rts[scan]]
    
    # Preliminary noise floor filter (speed up)
    global_noise <- .sonar_estimate_noise(dt$int[sample(nrow(dt), min(nrow(dt), 50000))])
    dt <- dt[int > global_noise]
    
    # 3. Dynamic ROI Clustering (Smart binning)
    data.table::setorder(dt, mz)
    # ROIs are single-linkage runs of points whose m/z gaps stay below
    # `roi_gap_da` (configurable; 0.02 Da is a safe separator for nominal-
    # resolution-ish grouping). NOTE: single-linkage can chain-merge distinct
    # species when small gaps line up -- narrow `roi_gap_da` if that matters.
    dt[, mz_group := cumsum(c(0, diff(mz) > roi_gap_da))]

    groups <- unique(dt$mz_group)
    pb <- txtProgressBar(min = 0, max = length(groups), style = 3)

    # Process each m/z ROI
    roi_results_list <- dt[, {
      if (.GRP %% 100 == 0) setTxtProgressBar(pb, .GRP)

      res_dt <- data.table::data.table()
      if (.N >= 5) {
        local_cutoff <- .sonar_estimate_noise(int)

        # Identify peaks in time axis. Index relative to `scan_range` rather
        # than allocating a full length-`max(scan)` vector per ROI: with many
        # ROIs the full-length version is O(n_roi x n_scan) memory.
        scan_range <- min(scan):max(scan)
        eic_segment <- numeric(length(scan_range))
        eic_segment[scan - scan_range[1L] + 1L] <- int
        
        # Smoothing
        eic_s <- .sonar_peak_smooth(eic_segment, smooth)
        
        # Peak segmentation
        above_idx <- which(eic_s > local_cutoff)
        if (length(above_idx) >= 5) {
          # Split into contiguous segments
          segs <- split(above_idx, cumsum(c(1, diff(above_idx) != 1)))
          
          # Calculate stats for each segment
          seg_stats <- lapply(segs, function(idx) {
            if (length(idx) >= 5) {
              actual_idx <- idx + min(scan_range) - 1
              # Find points belonging to this segment
              p_mask <- scan %in% actual_idx
              p_mz <- mz[p_mask]
              p_int <- int[p_mask]
              p_rt <- rt[p_mask]
              
              if (length(p_mz) < 5) return(NULL)
              
              ref_mz <- p_mz[which.max(p_int)]
              bg_mask <- p_int < (local_cutoff * 1.5)
              bg_sd <- if (sum(bg_mask) > 3) sd(p_int[bg_mask]) else (local_cutoff / 3)
              bg_mean <- if (sum(bg_mask) > 3) mean(p_int[bg_mask]) else 0
              
              sn_val <- (max(p_int) - bg_mean) / max(bg_sd, 1e-6)

              # Peak width at half height (FWHM), not the full baseline-to-
              # baseline span: CentWave's `peakwidth` is a half-height width, so
              # `max(rt) - min(rt)` over the whole above-cutoff segment
              # systematically over-estimates it. When the half-height band
              # spans fewer than two scans we fall back to the baseline span
              # scaled by 0.588 (the Gaussian FWHM / base-width ratio).
              pw_base <- max(p_rt) - min(p_rt)
              half_level <- local_cutoff + 0.5 * (max(p_int) - local_cutoff)
              half_rt <- p_rt[p_int >= half_level]
              pw_fwhm <- if (length(half_rt) >= 2L) {
                max(half_rt) - min(half_rt)
              } else {
                0.588 * pw_base
              }

              list(
                mz = ref_mz,
                ppm = (sd(p_mz) / ref_mz) * 1e6 * mass_sd_range,
                mz_span = max(p_mz) - min(p_mz),
                pw = pw_fwhm,
                pw_base = pw_base,
                scans = length(idx),
                height = max(p_int),
                noise = local_cutoff,
                sn = as.numeric(sn_val)
              )
            } else NULL
          })
          res_dt <- data.table::rbindlist(seg_stats)
        }
      }
      # Return a list containing the data.table to ensure consistent structure
      list(res = list(res_dt))
    }, by = .(mz_group)]
    
    setTxtProgressBar(pb, length(groups))
    close(pb)
    
    # 3b. Combine all ROI data for this file
    all_stats[[f]] <- data.table::rbindlist(roi_results_list$res, fill = TRUE)
  }
  
  res <- data.table::rbindlist(all_stats, fill = TRUE)
  
  if (is.null(res) || nrow(res) == 0) {
    stop(sprintf(
      paste0("Sonar failed: no valid chromatographic peaks detected in %s. ",
             "This usually means the data-driven noise floor removed every ",
             "point; check the file contents, or pass mode = 'dia' for Waters ",
             "MSE data (function=1 low-energy scans only)."),
      basename(target_files[1L])
    ), call. = FALSE)
  }
  
  # --- 4. Statistical Derivation ---
  params <- .sonar_derive_params(
    res,
    cutoff_quantile     = cutoff_quantile,
    max_suggested_ppm   = max_suggested_ppm,
    mzdiff_range        = mzdiff_range,
    prefilter_scans     = psr,
    prefilter_intensity = pir,
    min_fraction_default = min_fraction_default
  )

  message("DONE: Sonar Scan Complete.")
  message(sprintf("  > Detected optimal ppm: %s", params$ppm))
  message(sprintf("  > PeakWidth range: [%s, %s] s", params$peakwidth[1L], params$peakwidth[2L]))
  message(sprintf("  > Baseline noise: %s", params$noise))
  message(sprintf(
    "  > min_fraction = %.2f (placeholder; override in XcmsParam() if replicates exist)",
    params$min_fraction
  ))

  params
}
