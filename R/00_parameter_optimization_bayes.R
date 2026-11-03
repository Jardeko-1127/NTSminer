# Bayesian-optimisation-based CentWave parameter refinement.
#
# Wraps the ParBayesianOptimization workflow originally drafted in
# `tests/20260405_01_ParBayes_Pick_optimize.R`. Intended as a heavier-weight
# alternative to [suggest_xcms_params]: when the Sonar initial guess is not
# good enough, users can call [refine_xcms_params_bayes] on a small set of
# QC files (1--3 typical) to search the CentWave space directly against a
# peak-shape + simple M+1 isotope objective.

# Peak-shape score over a `chromPeaks(xpd)` matrix (one row per detected peak).
#
# Every peak is scored against absolute, non-self-referential criteria, so a set
# of "uniform but wrong" peaks can no longer score highly:
#   * symmetry     -- the apex should sit near the middle of the elution window;
#   * Gaussian law -- for an ideal Gaussian, area/height = sigma*sqrt(2*pi) and
#                     FWHM = 2.3548*sigma, hence (area/height)/FWHM = 1.0645.
#                     This is checked per peak against that constant, NOT
#                     against the run's own median. The previous version graded
#                     "how alike the peaks are" (1 - |x - median(x)|/median(x)),
#                     which is invariant to any global rescaling and awarded a
#                     perfect score to degenerate, uniform peak sets.
#   * plausibility -- an LC peak should occupy only a small fraction of the
#                     gradient. `rt_span` (seconds, the full acquisition range)
#                     anchors the score to an absolute scale, so globally
#                     inflated widths are penalised. This is the
#                     scale-awareness the previous objective lacked (it returned
#                     an identical score when all widths were multiplied by 10).
#' @noRd
.bo_peak_shape_score <- function(pks, rt_span = NULL) {
  if (!is.matrix(pks) || nrow(pks) < 1L) return(0)
  rt    <- pks[, "rt"]
  rtmin <- pks[, "rtmin"]
  rtmax <- pks[, "rtmax"]
  into  <- pks[, "into"]
  maxo  <- pks[, "maxo"]

  width <- pmax(rtmax - rtmin, 1e-9)

  sym <- (rt - rtmin) / width
  sym_score <- pmax(0, 1 - 2 * abs(sym - 0.5))

  ar_ratio <- into / pmax(maxo, 1e-9)
  gauss_ratio <- (ar_ratio / width) / 1.0645
  gauss_score <- pmax(0, 1 - abs(gauss_ratio - 1))

  if (is.finite(rt_span) && rt_span > 0) {
    width_score <- pmax(0, 1 - (width / rt_span) / 0.25)
  } else {
    width_score <- rep(1, length(width))
  }

  shape <- 0.5 * sym_score + 0.3 * gauss_score + 0.2 * width_score
  res <- mean(shape, na.rm = TRUE)
  if (is.nan(res) || is.na(res)) 0 else res
}

# M+1 isotope-pairing score over a `chromPeaks(xpd)` matrix.
#
# A peak is "corroborated" when a co-eluting peak (within `rt_win` seconds)
# exists at mz + mass_diff within `ppm_tol` and its intensity ratio is
# physically plausible for a 13C M+1 isotopologue (roughly 1.1% per carbon,
# i.e. ~0.005 at C2 up to ~0.5 at C45). We deliberately do NOT derive a
# "theoretical" ratio from m/z via round(mz/14) as the previous version did:
# that overestimates the carbon count for heteroatom-rich analytes (e.g. PFAS,
# where every F adds 19 Da), and its min(0.1, ...) cap collapsed every analyte
# above ~133 Da to an identical target ratio, erasing the information.
#
# Uses a sorted index + findInterval, so the pairing search is O(n log n)
# rather than O(n^2).
#' @noRd
.bo_isotope_score <- function(pks,
                              mass_diff = 1.00335,
                              rt_win = 1.5,
                              ppm_tol = 15,
                              ratio_range = c(0.005, 0.5)) {
  if (!is.matrix(pks) || nrow(pks) < 1L) return(0)
  mz  <- pks[, "mz"]
  int <- pks[, "into"]
  rt  <- pks[, "rt"]

  o     <- order(mz)
  mz_s  <- mz[o]
  int_s <- int[o]
  rt_s  <- rt[o]

  iso_scores <- numeric(0)
  for (i in seq_along(mz)) {
    m1  <- mz[i] + mass_diff
    tol <- m1 * ppm_tol * 1e-6
    lo  <- findInterval(m1 - tol, mz_s) + 1L
    hi  <- findInterval(m1 + tol, mz_s)
    if (hi < lo) next
    cand <- lo:hi
    cand <- cand[abs(rt_s[cand] - rt[i]) < rt_win]
    if (!length(cand)) next
    err    <- abs(mz_s[cand] - m1) / m1 * 1e6
    r_obs  <- int_s[cand] / int[i]
    in_band <- r_obs >= ratio_range[1] & r_obs <= ratio_range[2]
    s <- pmax(0, 1 - err / (2 * ppm_tol)) * in_band
    iso_scores <- c(iso_scores, s)
  }
  if (!length(iso_scores)) return(0)
  res <- mean(iso_scores, na.rm = TRUE)
  if (is.nan(res) || is.na(res)) 0 else res
}


#' Bayesian-Optimisation Refinement of CentWave Parameters
#'
#' Runs `ParBayesianOptimization::bayesOpt()` over the CentWave parameter space
#' against a composite objective (peak-shape + simple M+1 isotope score). This
#' is a **heavier, optional** alternative to [suggest_xcms_params]: each
#' iteration triggers a full [xcms::findChromPeaks] call, so wall-clock cost
#' scales with `(init_points + n_iters) x len(ms1_files) x CentWave_time`.
#' Typical use: 1--3 QC files, 20+20 iterations, total 20--90 minutes.
#'
#' Bounds default to those used in `tests/20260405_01_ParBayes_Pick_optimize.R`
#' and can be narrowed around a [suggest_xcms_params] result to speed
#' convergence (see Examples).
#'
#' Requires the **Suggests** package `ParBayesianOptimization`. Install via
#' `install.packages("ParBayesianOptimization")`.
#'
#' @param ms1_files Character vector of MS1 mzML / mzXML files (typically QC).
#'   Loaded once with [MSnbase::readMSData] and reused for every BO iteration.
#' @param mode One of `"auto"` / `"dda"` / `"dia"`. `"dia"` restricts the search
#'   to `function=1` (low-energy) MS1 scans, mirroring [suggest_xcms_params]
#'   and avoiding fragment contamination of Waters MSE data; `"dda"` uses every
#'   MS1 scan; `"auto"` (default) uses the DIA subset only when a `function=N`
#'   field is actually present in the `spectrumId`.
#' @param bounds Named list of length-2 numeric vectors. Required fields:
#'   `min_peakw`, `max_peakw`, `ppm`, `mzdiff`, `snthresh`, `noise`, `prefilter`,
#'   `val_prefilter`. Pass `NULL` for the reference defaults.
#' @param init_points Number of random initial evaluations (`bayesOpt(initPoints)`).
#' @param n_iters Number of subsequent GP-guided iterations (`bayesOpt(iters.n)`).
#' @param acq The **Bayesian-optimisation acquisition function** passed to
#'   `ParBayesianOptimization::bayesOpt()` -- one of `"ei"`, `"ucb"`, `"poi"`.
#'   Default `"ei"` (expected improvement). This is *not* an MS acquisition
#'   mode; see the separate `mode` argument for DDA/DIA selection.
#' @param eps Early-stop tolerance passed to `bayesOpt(eps = ...)`.
#' @param weights Named numeric, `c(shape = ..., isotope = ...)`. Normalised
#'   to sum 1 internally. Default `c(0.5, 0.5)` (equal mix).
#' @param min_peaks Penalty threshold: if `nrow(chromPeaks(xpd)) < min_peaks`
#'   the iteration returns score `0.001` (effectively a dead trial). If *every*
#'   iteration returns the floor, a warning is emitted at the end because the
#'   search then carries no information (lower `min_peaks` for small QC files).
#' @param iso_mass_diff,iso_rt_window,iso_ppm_tol Tuning of the M+1 isotope
#'   pairing rule (Da / s / ppm).
#' @param iso_ratio_range Length-2 numeric, plausible `M+1/M` intensity-ratio
#'   band (default `c(0.005, 0.5)`) used by the isotope corroboration score.
#' @param bpparam Optional `BiocParallelParam` passed to [xcms::findChromPeaks]
#'   for each iteration. Default `NULL` -> [BiocParallel::SerialParam] (safer
#'   for repeated calls on small file sets).
#' @param log_file Optional path to a text file for per-iteration logs.
#' @param seed Integer RNG seed for reproducibility (passed to [set.seed]).
#'   Reproducibility holds for `parallel = FALSE`; with `parallel = TRUE` the
#'   evaluation order is non-deterministic and results are not reproducible.
#' @param parallel Forwarded to `bayesOpt(parallel = ...)`. Most reliable as
#'   `FALSE`; CentWave parallelism is already controlled by `bpparam`.
#' @param verbose `bayesOpt` verbosity *and* whether to print per-iteration
#'   diagnostic `message()`s. Use `0L` for silent runs.
#' @param keep_opt_res If `TRUE` (default) include the full `bayesOpt` object
#'   in the return value (useful for `plot()` / diagnostics). Set `FALSE` for
#'   compact downstream serialisation.
#'
#' @return A named list. The CentWave fields (`ppm`, `peakwidth`, `snthresh`,
#'   `noise`, `prefilter`, `mzdiff`) are available both at the top level and
#'   inside `params`, so the result is drop-in compatible with the historical
#'   top-level access (`refined$ppm`) *and* with the documented
#'   `refined$params` idiom. Also contains a `bo` sub-list (best score,
#'   iteration count, `score_summary`, and optionally the full `opt_res`).
#'   Note: this routine does **not** suggest PeakDensity / Obiwarp parameters
#'   (use [suggest_xcms_params] or set them manually in [XcmsParam]).
#'
#' @seealso [suggest_xcms_params] for the lightweight Sonar initial guess;
#'   [XcmsParam] for plugging the result into [extract_features].
#'
#' @examples
#' \dontrun{
#' qc_files <- list.files("data/QC", pattern = "\\.mzML$", full.names = TRUE)
#'
#' # (Optional) start from a Sonar guess and narrow the bounds around it
#' sonar <- suggest_xcms_params(qc_files, sample_n = length(qc_files))
#' refined <- refine_xcms_params_bayes(
#'   qc_files,
#'   bounds = list(
#'     min_peakw     = sonar$peakwidth[1] + c(-3, 5),
#'     max_peakw     = sonar$peakwidth[2] + c(-10, 10),
#'     ppm           = c(max(3, sonar$ppm * 0.5), sonar$ppm * 1.5),
#'     mzdiff        = c(-0.001, 0.01),
#'     snthresh      = c(3, 10),
#'     noise         = c(0, max(1000, sonar$noise * 2)),
#'     prefilter     = c(2, 4),
#'     val_prefilter = c(50, max(100, sonar$prefilter[2]))
#'   ),
#'   init_points = 15, n_iters = 15
#' )
#'
#' param <- XcmsParam(
#'   ppm = refined$params$ppm, peakwidth = refined$params$peakwidth,
#'   snthresh = refined$params$snthresh, noise = refined$params$noise,
#'   prefilter = refined$params$prefilter, mzdiff = refined$params$mzdiff,
#'   binSize = sonar$binSize, bw = sonar$bw,
#'   min_fraction = sonar$min_fraction,
#'   obiwarp_binSize = sonar$obiwarp_binSize
#' )
#' }
#' @export
refine_xcms_params_bayes <- function(ms1_files,
                                     mode = c("auto", "dda", "dia"),
                                     bounds = NULL,
                                     init_points = 20L,
                                     n_iters = 20L,
                                     acq = "ei",
                                     eps = 0.01,
                                     weights = c(shape = 0.5, isotope = 0.5),
                                     min_peaks = 200L,
                                     iso_mass_diff = 1.00335,
                                     iso_rt_window = 1.5,
                                     iso_ppm_tol = 15,
                                     iso_ratio_range = c(0.005, 0.5),
                                     bpparam = NULL,
                                     log_file = NULL,
                                     seed = 123L,
                                     parallel = FALSE,
                                     verbose = 1L,
                                     keep_opt_res = TRUE) {
  if (!requireNamespace("ParBayesianOptimization", quietly = TRUE)) {
    stop(
      "`refine_xcms_params_bayes()` requires the suggested package ",
      "'ParBayesianOptimization'. Install it via ",
      "`install.packages('ParBayesianOptimization')`.",
      call. = FALSE
    )
  }
  if (!requireNamespace("xcms", quietly = TRUE)) stop("xcms is required")
  if (!requireNamespace("MSnbase", quietly = TRUE)) stop("MSnbase is required")

  mode <- match.arg(mode)
  acq  <- match.arg(acq, c("ei", "ucb", "poi"))

  ms1_files <- path.expand(as.character(ms1_files))
  ms1_files <- ms1_files[file.exists(ms1_files)]
  if (!length(ms1_files)) {
    stop("No valid MS1 files provided.", call. = FALSE)
  }

  w_names <- c("shape", "isotope")
  if (is.null(names(weights))) {
    if (length(weights) != 2L) {
      stop("`weights` must be length 2 (shape, isotope)", call. = FALSE)
    }
    names(weights) <- w_names
  }
  if (!all(w_names %in% names(weights))) {
    stop("`weights` must contain `shape` and `isotope`", call. = FALSE)
  }
  w <- as.numeric(weights[w_names])
  if (anyNA(w) || any(w < 0) || sum(w) <= 0) {
    stop("`weights` must be non-negative with positive sum", call. = FALSE)
  }
  w <- w / sum(w)

  iso_ratio_range <- as.numeric(iso_ratio_range)[1:2]
  if (length(iso_ratio_range) != 2L || anyNA(iso_ratio_range) ||
      iso_ratio_range[1] <= 0 || iso_ratio_range[1] >= iso_ratio_range[2]) {
    stop("`iso_ratio_range` must be two positive numbers with lower < upper",
         call. = FALSE)
  }

  default_bounds <- list(
    min_peakw     = c(5, 25),
    max_peakw     = c(30, 60),
    ppm           = c(5, 25),
    mzdiff        = c(-0.001, 0.01),
    snthresh      = c(3, 10),
    noise         = c(0, 1000),
    prefilter     = c(2, 4),
    val_prefilter = c(50, 100)
  )
  if (is.null(bounds)) {
    bounds <- default_bounds
  } else {
    miss <- setdiff(names(default_bounds), names(bounds))
    if (length(miss)) {
      stop("`bounds` missing fields: ", paste(miss, collapse = ", "), call. = FALSE)
    }
    bounds <- bounds[names(default_bounds)]
    for (nm in names(bounds)) {
      v <- as.numeric(bounds[[nm]])
      if (length(v) != 2L || anyNA(v) || v[1] >= v[2]) {
        stop(sprintf("`bounds$%s` must be a length-2 numeric with lower < upper", nm),
             call. = FALSE)
      }
      bounds[[nm]] <- v
    }
  }

  bpp_used <- bpparam
  if (is.null(bpp_used)) bpp_used <- BiocParallel::SerialParam()

  if (verbose >= 1L) {
    message(sprintf("Loading %d file(s) for Bayesian refinement...",
                    length(ms1_files)))
  }
  raw_data <- MSnbase::readMSData(ms1_files, msLevel = 1L, mode = "onDisk")

  # --- DIA branch: restrict to function=1 (low-energy) MS1 scans ---------------
  if (mode != "dda") {
    fnv_list <- lapply(ms1_files, .sonar_try_parse_function)
    has_fn <- vapply(fnv_list,
                     function(v) !is.null(v) && any(!is.na(v)), logical(1))
    if (mode == "dia" && !all(has_fn)) {
      stop("mode = 'dia' but `function=N` was not found in the spectrumId of: ",
           paste(basename(ms1_files[!has_fn]), collapse = ", "), call. = FALSE)
    }
    if (any(has_fn)) {
      fidx <- MSnbase::fromFile(raw_data)
      per_file <- split(seq_along(fidx), fidx)
      keep <- integer(0)
      for (fi in seq_along(ms1_files)) {
        pos <- per_file[[as.character(fi)]]
        if (is.null(pos) || !length(pos)) next
        fnv <- fnv_list[[fi]]
        if (is.null(fnv) || length(fnv) != length(pos)) {
          keep <- c(keep, pos)
          next
        }
        k <- pos[!is.na(fnv) & fnv == 1L]
        keep <- c(keep, if (length(k) >= 5L) k else pos)
      }
      if (length(keep) >= 5L && length(keep) < length(fidx)) {
        if (verbose >= 1L) {
          message(sprintf("  > DIA mode: kept %d / %d MS1 spectra (function=1)",
                          length(keep), length(fidx)))
        }
        raw_data <- raw_data[keep]
      }
    }
  }

  # Full acquisition time range (seconds): anchors the absolute-scale width term
  # of the shape objective.
  rts_all <- MSnbase::rtime(raw_data)
  rt_span <- if (length(rts_all) > 1L) diff(range(rts_all)) else NA_real_

  iter_count <- 0L
  if (!is.null(log_file)) {
    cat(sprintf("=== BO start: %s ===\n", Sys.time()), file = log_file)
  }

  fn <- function(min_peakw, max_peakw, ppm, mzdiff,
                 snthresh, noise, prefilter, val_prefilter) {
    iter_count <<- iter_count + 1L

    min_peakw     <- as.integer(round(min_peakw))
    max_peakw     <- as.integer(round(max_peakw))
    snthresh      <- as.integer(round(snthresh))
    noise         <- as.integer(round(noise))
    prefilter     <- as.integer(round(prefilter))
    val_prefilter <- as.integer(round(val_prefilter))
    ppm           <- round(ppm, 1)
    mzdiff        <- round(mzdiff, 4)

    if (max_peakw <= min_peakw) return(list(Score = 0.001))

    cwp <- xcms::CentWaveParam(
      peakwidth = c(min_peakw, max_peakw),
      ppm = ppm, mzdiff = mzdiff,
      snthresh = snthresh, noise = noise,
      prefilter = c(prefilter, val_prefilter),
      mzCenterFun = "wMean", integrate = 1L, fitgauss = FALSE
    )

    score <- tryCatch({
      xpd <- xcms::findChromPeaks(raw_data, param = cwp, BPPARAM = bpp_used)
      pks <- xcms::chromPeaks(xpd)
      if (NROW(pks) < min_peaks) {
        if (verbose >= 1L) {
          message(sprintf("[BO %d] only %d peaks (<%d); Score=0.001",
                          iter_count, NROW(pks), min_peaks))
        }
        0.001
      } else {
        ss <- .bo_peak_shape_score(pks, rt_span = rt_span)
        is_ <- .bo_isotope_score(pks, iso_mass_diff, iso_rt_window, iso_ppm_tol,
                                 iso_ratio_range)
        s <- w[1] * ss + w[2] * is_
        s <- pmax(0.001, pmin(0.999, s + stats::runif(1L, -1e-5, 1e-5)))
        if (verbose >= 1L) {
          message(sprintf(
            "[BO %d] shape=%.3f iso=%.3f total=%.3f | pw=[%d,%d] ppm=%.1f mzdiff=%.4f sn=%d noise=%d pre=c(%d,%d)",
            iter_count, ss, is_, s, min_peakw, max_peakw, ppm, mzdiff,
            snthresh, noise, prefilter, val_prefilter
          ))
        }
        if (!is.null(log_file)) {
          cat(sprintf(
            "[%d] shape=%.3f iso=%.3f total=%.3f | pw=[%d,%d] ppm=%.1f mzdiff=%.4f sn=%d noise=%d pre=c(%d,%d)\n",
            iter_count, ss, is_, s, min_peakw, max_peakw, ppm, mzdiff,
            snthresh, noise, prefilter, val_prefilter
          ), file = log_file, append = TRUE)
        }
        s
      }
    }, error = function(e) {
      if (verbose >= 1L) {
        message(sprintf("[BO %d] CentWave error: %s",
                        iter_count, conditionMessage(e)))
      }
      if (!is.null(log_file)) {
        cat(sprintf("[%d] ERROR: %s\n", iter_count, conditionMessage(e)),
            file = log_file, append = TRUE)
      }
      0.001
    })

    list(Score = score)
  }

  set.seed(seed)
  if (verbose >= 1L) {
    message(sprintf(
      "Running ParBayesianOptimization (init=%d, iter=%d, acq=%s)...",
      init_points, n_iters, acq
    ))
  }
  opt_res <- ParBayesianOptimization::bayesOpt(
    FUN        = fn,
    bounds     = bounds,
    initPoints = init_points,
    iters.n    = n_iters,
    acq        = acq,
    eps        = eps,
    parallel   = parallel,
    verbose    = verbose
  )

  best <- ParBayesianOptimization::getBestPars(opt_res)
  best_score <- suppressWarnings(max(opt_res$scoreSummary$Score, na.rm = TRUE))
  if (!is.finite(best_score)) best_score <- NA_real_

  all_scores <- opt_res$scoreSummary$Score
  if (!any(is.finite(all_scores) & all_scores > 0.0011)) {
    warning(
      "Every BO iteration returned the score floor (0.001): no parameter set ",
      "in the search space produced at least `min_peaks` peaks. The result is ",
      "not meaningful -- lower `min_peaks`, widen `bounds`, or check the data.",
      call. = FALSE
    )
  }

  cw_params <- list(
    ppm       = round(best$ppm, 1),
    peakwidth = c(as.integer(round(best$min_peakw)),
                  as.integer(round(best$max_peakw))),
    snthresh  = as.integer(round(best$snthresh)),
    noise     = as.integer(round(best$noise)),
    prefilter = c(as.integer(round(best$prefilter)),
                  as.integer(round(best$val_prefilter))),
    mzdiff    = round(best$mzdiff, 4)
  )
  # CentWave fields are exposed both at the top level (historical access such as
  # `refined$ppm`) and inside `params` (the documented idiom, `refined$params`).
  out <- c(cw_params, list(
    params = cw_params,
    bo = list(
      method        = "ParBayesianOptimization",
      best_score    = best_score,
      n_iters_total = init_points + n_iters,
      acq           = acq,
      mode          = mode,
      weights       = w,
      bounds        = bounds,
      score_summary = opt_res$scoreSummary
    )
  ))
  if (isTRUE(keep_opt_res)) out$bo$opt_res <- opt_res

  if (!is.null(log_file)) {
    cat(sprintf("=== BO end: %s | best=%.4f ===\n",
                Sys.time(), best_score),
        file = log_file, append = TRUE)
  }

  if (verbose >= 1L) {
    message("===== BO done =====")
    message(sprintf("Best score: %.4f", best_score))
    message(sprintf(
      "ppm=%s peakwidth=[%d,%d] snthresh=%d noise=%d prefilter=c(%d,%d) mzdiff=%g",
      out$ppm, out$peakwidth[1], out$peakwidth[2], out$snthresh, out$noise,
      out$prefilter[1], out$prefilter[2], out$mzdiff
    ))
  }

  out
}
