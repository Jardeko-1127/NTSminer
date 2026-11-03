test_that(".nts_trap_integrate computes trapezoid area", {
  rt <- c(0, 1, 2)
  y <- c(0, 10, 0)
  area <- NTSminer:::.nts_trap_integrate(rt, y)
  expect_equal(area, 10)
})

test_that("refresh_feature_areas validates project input", {
  expect_error(
    refresh_feature_areas(data.frame()),
    "project"
  )
})

test_that("XcmsParam accepts feature_area_method", {
  p <- XcmsParam(feature_area_method = "eic_rt30", feature_area_rt_center = "rtmed")
  expect_equal(p$feature_area_method, "eic_rt30")
  expect_equal(p$feature_area_rt_center, "rtmed")
  expect_equal(p$feature_area_ppm, 5)
  expect_error(XcmsParam(feature_area_ppm = 0), "feature_area_ppm")
})

test_that(".nts_validate_feature_area_method rejects unknown values", {
  expect_error(
    NTSminer:::.nts_validate_feature_area_method("invalid"),
    "should be one of"
  )
})

test_that("quantify_db_annotations returns empty frame without annotations", {
  skip_if_not(requireNamespace("methods", quietly = TRUE))
  proj <- methods::new("NTSProject")
  out <- quantify_db_annotations(proj, annotations = data.frame())
  expect_equal(nrow(out), 0L)
})

test_that("quantify_db_annotations validates threads", {
  skip_if_not(requireNamespace("methods", quietly = TRUE))
  proj <- methods::new("NTSProject")
  expect_error(quantify_db_annotations(proj, threads = 0), "threads")
})
