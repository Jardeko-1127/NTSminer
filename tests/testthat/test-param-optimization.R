# Unit tests for the data-driven parameter-optimisation functions
#   suggest_xcms_params()          (Sonar)  -- R/00_parameter_optimization.R
#   refine_xcms_params_bayes()     (BO)     -- R/00_parameter_optimization_bayes.R
#
# Focus: the pure helpers behind both functions, so the quantile / clamping /
# scoring logic is covered without requiring mzML files or the optional
# ParBayesianOptimization package. Covers the regressions fixed for 0.1.17:
#   * empty / all-NA quantile pools no longer raise (E-2)
#   * the noise estimator falls back to the median, not max (S-4)
#   * the ppm estimator is not pinned to the cap by a heavy tail (S-2/S-6)
#   * the BO shape score is no longer invariant to peak width (S-9)

# --------------------------------------------------------------------------
# .sonar_estimate_noise
# --------------------------------------------------------------------------
test_that(".sonar_estimate_noise returns 0 for an all-zero / empty vector", {
  expect_identical(NTSminer:::.sonar_estimate_noise(numeric(0)), 0)
  expect_identical(NTSminer:::.sonar_estimate_noise(c(0, 0, 0)), 0)
})

test_that(".sonar_estimate_noise returns the median for short vectors", {
  v <- c(1, 2, 3, 4, 5)
  expect_equal(NTSminer:::.sonar_estimate_noise(v), stats::median(v))
})

test_that(".sonar_estimate_noise does not error at n == 10 (regression)", {
  # seq.int(10, n - 1, 10) with n == 10 used to raise "wrong sign in 'by'".
  v <- c(seq_len(9) * 10, 1e6)
  expect_error(NTSminer:::.sonar_estimate_noise(v), NA)
  expect_true(is.finite(NTSminer:::.sonar_estimate_noise(v)))
})

test_that(".sonar_estimate_noise finds the knee between noise and signal", {
  set.seed(1)
  v <- c(stats::runif(200, 1, 50), stats::runif(200, 5000, 9000))
  est <- NTSminer:::.sonar_estimate_noise(v)
  expect_true(est > 40 && est < 2000)
})

test_that(".sonar_estimate_noise falls back to the median, never max (S-4)", {
  set.seed(2)
  v <- stats::runif(5000, 1, 100)          # no knee: running mean + 3sd never crossed
  est <- NTSminer:::.sonar_estimate_noise(v)
  expect_true(est < max(v))                # would equal max under the old fallback
  expect_equal(est, stats::median(v))
})

# --------------------------------------------------------------------------
# .sonar_peak_smooth
# --------------------------------------------------------------------------
test_that(".sonar_peak_smooth is a no-op for level <= 0 or short input", {
  x <- c(1, 5, 2, 8, 3)
  expect_identical(NTSminer:::.sonar_peak_smooth(x, 0), x)
  expect_identical(NTSminer:::.sonar_peak_smooth(x, -1), x)
  expect_identical(NTSminer:::.sonar_peak_smooth(x, 5), x)   # length(x) < 2n
})

test_that(".sonar_peak_smooth preserves length and leaves constants alone", {
  x <- c(1, 100, 2, 90, 3, 80, 4)
  y <- NTSminer:::.sonar_peak_smooth(x, 1)
  expect_length(y, length(x))
  expect_identical(NTSminer:::.sonar_peak_smooth(rep(4, 10), 2), rep(4, 10))
})

# --------------------------------------------------------------------------
# .sonar_derive_params  (pure table -> params)
# --------------------------------------------------------------------------
make_roi <- function(n = 300, ppm = NULL, pw = NULL, noise = 1000) {
  data.table::data.table(
    mz      = stats::runif(n, 200, 600),
    ppm     = if (is.null(ppm)) stats::runif(n, 2, 6) else ppm,
    mz_span = stats::runif(n, 0.001, 0.01),
    pw      = if (is.null(pw)) stats::runif(n, 4, 12) else pw,
    scans   = sample(5:20, n, replace = TRUE),
    height  = stats::runif(n, 5e3, 5e5),
    noise   = rep(noise, n),
    sn      = stats::runif(n, 3, 50)
  )
}

test_that(".sonar_derive_params returns the full named contract", {
  p <- NTSminer:::.sonar_derive_params(make_roi())
  expect_named(p, c("ppm", "peakwidth", "snthresh", "noise", "prefilter",
                    "mzdiff", "min_fraction", "binSize", "bw", "obiwarp_binSize"))
  expect_length(p$peakwidth, 2L)
  expect_length(p$prefilter, 2L)
  expect_true(p$peakwidth[1L] >= 3)
  expect_true(p$peakwidth[2L] > p$peakwidth[1L])
  expect_true(p$noise >= 100)
  expect_identical(p$obiwarp_binSize, 0.6)
  expect_identical(p$min_fraction, 0.001)
})

test_that(".sonar_derive_params survives an empty / all-NA peakwidth pool (E-2)", {
  p <- NTSminer:::.sonar_derive_params(make_roi(n = 50, pw = rep(NA_real_, 50)))
  expect_identical(p$peakwidth, c(3, 30))
  p2 <- NTSminer:::.sonar_derive_params(make_roi(n = 50, pw = rep(-1, 50)))
  expect_identical(p2$peakwidth, c(3, 30))
})

test_that(".sonar_derive_params floors a degenerate ppm pool at 5", {
  p <- NTSminer:::.sonar_derive_params(make_roi(n = 50, ppm = rep(NA_real_, 50)))
  expect_identical(p$ppm, 5)
})

test_that(".sonar_derive_params is not pinned to the cap by a heavy tail (S-2/S-6)", {
  # 200 clean ROIs at ~5 ppm + 20 pathological ones at ~150 ppm. The raw 95th
  # percentile sits in the tail and used to force the value to `max_suggested_ppm`.
  set.seed(3)
  ppm <- c(stats::rnorm(200, 5, 1), rep(150, 20))
  ppt <- abs(ppm[ppm > 0])
  p <- NTSminer:::.sonar_derive_params(make_roi(n = 220, ppm = ppt))
  expect_true(p$ppm >= 5)
  expect_true(p$ppm < 80)          # not pinned to the default ceiling
  expect_true(p$ppm <= 20)         # close to the clean core of the distribution
})

test_that(".sonar_derive_params honours max_suggested_ppm when genuinely exceeded", {
  set.seed(4)
  p <- NTSminer:::.sonar_derive_params(make_roi(n = 100, ppm = rep(30, 100)),
                                       max_suggested_ppm = 25L)
  expect_identical(p$ppm, 25L)
})

test_that(".sonar_derive_params anchors prefilter intensity to the noise floor", {
  p <- NTSminer:::.sonar_derive_params(make_roi(n = 100, noise = 800))
  expect_true(p$prefilter[2L] >= 800 && p$prefilter[2L] <= 1000)
  expect_true(p$prefilter[1L] >= 2L && p$prefilter[1L] <= 4L)
})

# --------------------------------------------------------------------------
# .bo_peak_shape_score  (S-9: must NOT be invariant to peak width)
# --------------------------------------------------------------------------
make_peaks <- function(fwhm = 5, maxo = 1000, center = 100, n = 5) {
  sigma <- fwhm / 2.3548
  into  <- rep(maxo * sigma * sqrt(2 * pi), n)   # ideal Gaussian area
  rt    <- rep(center, n)
  matrix(c(rt, rt - fwhm / 2, rt + fwhm / 2, into, rep(maxo, n)),
         ncol = 5L, byrow = FALSE,
         dimnames = list(NULL, c("rt", "rtmin", "rtmax", "into", "maxo")))
}

test_that(".bo_peak_shape_score rewards clean symmetric Gaussian peaks", {
  s <- NTSminer:::.bo_peak_shape_score(make_peaks(), rt_span = 1000)
  expect_true(s > 0.9 && s < 1)          # high, but not a degenerate perfect 1.0
})

test_that(".bo_peak_shape_score is sensitive to peak width (S-9 regression)", {
  s_narrow <- NTSminer:::.bo_peak_shape_score(make_peaks(fwhm = 5), rt_span = 1000)
  s_wide   <- NTSminer:::.bo_peak_shape_score(make_peaks(fwhm = 50), rt_span = 1000)
  expect_false(isTRUE(all.equal(s_narrow, s_wide)))
  expect_true(s_wide < s_narrow)         # 10x over-wide peaks must score worse
})

test_that(".bo_peak_shape_score does not return a degenerate 1.0 for identical peaks", {
  s <- NTSminer:::.bo_peak_shape_score(make_peaks(n = 300), rt_span = 1000)
  expect_true(s < 1)
})

test_that(".bo_peak_shape_score returns 0 for an empty / non-matrix input", {
  expect_identical(NTSminer:::.bo_peak_shape_score(NULL), 0)
  expect_identical(NTSminer:::.bo_peak_shape_score(matrix(numeric(0), 0, 5)), 0)
})

# --------------------------------------------------------------------------
# .bo_isotope_score
# --------------------------------------------------------------------------
make_iso <- function(mz1 = 300, ratio = 0.1, delta = 1.00335, rt = 100) {
  matrix(c(mz1, mz1 + delta, 1000, 1000 * ratio, rt, rt),
         ncol = 3L, dimnames = list(NULL, c("mz", "into", "rt")))
}

test_that(".bo_isotope_score detects a correctly-spaced M+1 partner", {
  s <- NTSminer:::.bo_isotope_score(make_iso(mz1 = 300, ratio = 0.1))
  expect_true(s > 0)
})

test_that(".bo_isotope_score works at high m/z (no constant-ratio degeneracy, S-10)", {
  expect_true(NTSminer:::.bo_isotope_score(make_iso(mz1 = 800, ratio = 0.08)) > 0)
  expect_true(NTSminer:::.bo_isotope_score(make_iso(mz1 = 133, ratio = 0.06)) > 0)
})

test_that(".bo_isotope_score is 0 without a co-eluting partner", {
  expect_identical(NTSminer:::.bo_isotope_score(make_iso(delta = 5)), 0)
  # partner present but far in RT -> not co-eluting
  m <- make_iso()
  m[2, "rt"] <- m[2, "rt"] + 100
  expect_identical(NTSminer:::.bo_isotope_score(m), 0)
})

test_that(".bo_isotope_score rejects out-of-band intensity ratios", {
  expect_identical(NTSminer:::.bo_isotope_score(make_iso(ratio = 0.9)), 0)  # > 0.5
  expect_identical(NTSminer:::.bo_isotope_score(make_iso(ratio = 0.001)), 0) # < 0.005
})

# --------------------------------------------------------------------------
# refine_xcms_params_bayes argument contract
# (skipped when the optional ParBayesianOptimization is absent)
# --------------------------------------------------------------------------
test_that("refine_xcms_params_bayes validates its arguments", {
  skip_if_not_installed("ParBayesianOptimization")
  f <- tempfile(fileext = ".mzML"); file.create(f)
  on.exit(unlink(f), add = TRUE)
  # NB: match.arg() reports the error as "'arg' should be one of ...", so match
  # on that (the argument name itself does not appear in the message).
  expect_error(refine_xcms_params_bayes(f, acq = "nope"), "should be one of")
  expect_error(refine_xcms_params_bayes(f, weights = c(1, 2, 3)), "weights")
  expect_error(refine_xcms_params_bayes(f, iso_ratio_range = c(0.5, 0.05)),
               "iso_ratio_range")
})
