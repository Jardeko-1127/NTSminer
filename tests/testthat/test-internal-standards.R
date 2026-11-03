make_is_ft <- function() {
  data.frame(
    name = c("IS_early", "IS_late", "FT_mid", "FT_late"),
    mz = c(100, 300, 200, 300.0001),
    rt = c(100, 500, 250, 480),
    sample_A = c(1000, 2000, 500, 800),
    sample_B = c(1000, 4000, 600, 1200),
    stringsAsFactors = FALSE
  )
}

test_that("read_internal_standards parses CSV", {
  p <- testthat::test_path("fixtures", "internal_standards.csv")
  is_tbl <- read_internal_standards(p)
  expect_equal(nrow(is_tbl), 2L)
  expect_true(all(c("is_name", "mz", "rt") %in% names(is_tbl)))
})

test_that("match_internal_standards finds features by mz/rt", {
  ft <- make_is_ft()
  is_tbl <- data.frame(
    is_name = c("IS_early", "IS_late"),
    mz = c(100, 300),
    rt = c(100, 500),
    mz_tol_ppm = 15,
    rt_tol_sec = 60,
    amount = c(50, 50),
    stringsAsFactors = FALSE
  )
  hit <- match_internal_standards(ft, is_tbl, pick = "closest_rt")
  expect_equal(nrow(hit), 2L)
  expect_equal(hit$feature_name, c("IS_early", "IS_late"))
  expect_true(all(hit$match_status == "ok"))
})

test_that("normalize_by_IS nearest_rt uses RT-local IS", {
  ft <- make_is_ft()
  out <- normalize_by_IS(
    ft,
    is_names = c("IS_early", "IS_late"),
    method = "nearest_rt",
    scale_factor = 1
  )
  expect_true(all(out$is_normalise_method == "nearest_rt"))
  expect_equal(out$is_anchor[ft$name == "FT_mid"], "IS_early")
  expect_equal(out$is_anchor[ft$name == "FT_late"], "IS_late")
  expect_equal(out$sample_A[ft$name == "FT_mid"], 500 / 1000)
  expect_equal(out$sample_A[ft$name == "FT_late"], 800 / 2000)
})

test_that("normalize_by_IS geom_mean matches legacy behaviour", {
  ft <- data.frame(
    name = c("IS1", "IS2", "FT01"),
    mz = c(100, 150, 200),
    rt = c(10, 15, 20),
    s1 = c(1000, 4000, 500),
    stringsAsFactors = FALSE
  )
  out <- normalize_by_IS(ft, is_names = c("IS1", "IS2"), method = "geom_mean", scale_factor = 1)
  geo_mean <- exp(mean(log(c(1000, 4000))))
  expect_equal(out$s1[3], 500 / geo_mean, tolerance = 1e-6)
})

test_that("run_semi_quantification returns estimated_conc on feature table", {
  ft <- make_is_ft()
  is_tbl <- data.frame(
    is_name = c("IS_early", "IS_late"),
    mz = c(100, 300),
    rt = c(100, 500),
    mz_tol_ppm = 20,
    rt_tol_sec = 60,
    amount = c(50, 50),
    stringsAsFactors = FALSE
  )
  res <- run_semi_quantification(
    ft,
    is_spec = is_tbl,
    normalize_method = "nearest_rt",
    scale_factor = 1,
    default_rie = 1
  )
  expect_true("estimated_conc" %in% names(res$feature_table))
  expect_true("rie" %in% names(res$feature_table))
  expect_equal(nrow(res$is_match), 2L)
})
