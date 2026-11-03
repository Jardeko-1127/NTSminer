# Tests for Kendrick mass-defect (KMD) series detection ------------------------

.cf2 <- 49.9968064  # CF2 (CF2 = C F2) exact mass

.kmd_series_ft <- function() {
  # A synthetic 5-member PFAS-like CF2 homologous series, with rt and abundance
  # increasing monotonically along the series so the RT/Spearman filters pass.
  mz <- 400.0 + .cf2 * (0:4)
  data.frame(
    name = paste0("S", 1:5),
    mz = mz,
    rt = 100 + 10 * (0:4),
    s1 = c(100, 200, 300, 400, 500),
    s2 = c(200, 400, 600, 800, 1000),
    s3 = c(400, 800, 1200, 1600, 2000),
    stringsAsFactors = FALSE
  )
}


test_that("calculate_kmd assigns one CF2 series to a homologous ladder", {
  res <- calculate_kmd(.kmd_series_ft(), base_formulas = "CF2")
  expect_true("kmd_series_id" %in% names(res))
  expect_equal(sum(!is.na(res$kmd_series_id)), 5L)
  # All members share a single series id.
  expect_length(unique(stats::na.omit(res$kmd_series_id)), 1L)
})


test_that("calculate_kmd KMD columns are internally consistent", {
  res <- calculate_kmd(.kmd_series_ft(), base_formulas = "CF2")
  expect_true(all(c("KM", "NKM", "KMD", "mz_remainder") %in% names(res)))
  # KMD = NKM - KM by construction.
  expect_equal(res$KMD, res$NKM - res$KM, tolerance = 1e-9)
})


test_that("calculate_kmd excludes isotope-labelled rows from series", {
  ft <- .kmd_series_ft()
  ft$is_isotope <- c(FALSE, FALSE, TRUE, FALSE, FALSE)  # S3 is an isotope
  res <- calculate_kmd(ft, base_formulas = "CF2")

  expect_true(is.na(res$kmd_series_id[res$name == "S3"]))
  # A gap of one isotope row is tolerated (max_gap = 2), so the remaining four
  # members still form one series.
  expect_equal(sum(!is.na(res$kmd_series_id)), 4L)
})


test_that("calculate_kmd honours min_series_size", {
  ft <- .kmd_series_ft()
  res <- calculate_kmd(ft, base_formulas = "CF2", min_series_size = 6)
  expect_equal(sum(!is.na(res$kmd_series_id)), 0L)
})


test_that("calculate_kmd rejects an out-of-range mass-defect series", {
  ft <- .kmd_series_ft()
  # Shrink md_range so the (negative) CF2 mass defects fall outside it.
  res <- calculate_kmd(ft, base_formulas = "CF2", md_range = c(-0.001, 0.001))
  expect_equal(sum(!is.na(res$kmd_series_id)), 0L)
})


test_that("calculate_kmd requires mz and rt columns", {
  ft <- data.frame(name = "x", mz = 400)
  expect_error(calculate_kmd(ft), "rt")
})


test_that("calculate_kmd validates kmd_tol", {
  expect_error(
    calculate_kmd(.kmd_series_ft(), base_formulas = "CF2", kmd_tol = -1),
    "kmd_tol"
  )
})
