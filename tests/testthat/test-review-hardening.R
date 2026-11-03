# Regression / hardening from 0.1.4 review checklist (no network / no large raw data).

test_that("XcmsParam validates numeric CentWave-related fields", {
  expect_error(XcmsParam(peakwidth = c("a", 50)), "peakwidth")
  expect_error(XcmsParam(ppm = -1), "ppm")
  expect_error(XcmsParam(threads = 0L), "threads")
  expect_s3_class(XcmsParam(), "XcmsParam")
})

test_that("calculate_exact_mass does not mis-check length(formula) after coercion", {
  skip_if_not_installed("Rdisop")
  m <- calculate_exact_mass(c("H2O", "CO2"))
  expect_true(is.finite(m))
  expect_true(m > 0)
})

test_that("extended NTSProject accessors round-trip", {
  d <- tempfile("nts_test_proj_")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  p <- create_nts_project(d)
  expect_true(is.list(get_project_log(p)))

  p2 <- set_db_annotations(
    p,
    data.frame(name = "a", stringsAsFactors = FALSE)
  )
  expect_equal(nrow(get_db_annotations(p2)), 1L)

  p3 <- set_ms2data(p2, list(FT1 = list()))
  expect_equal(length(get_ms2data(p3)), 1L)
  expect_named(get_ms2data(p3), "FT1")
})
