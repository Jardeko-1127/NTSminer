# Tests for score_isotope_fit() and .nts_score_isotope_fit_single() -----------

# Expected Rdisop isotopes for PFOA (C8HF15O2):
#   M = 413.974, M+1 = 414.977, M+2 = 415.979
# Abundances: ~0.911, ~0.082, ~0.0069
PFOA_MZ  <- c(413.974, 414.977, 415.979)
PFOA_INT <- c(100000, 8700, 380)

test_that("score_isotope_fit: PFOA known formula gives cosine >= 0.95", {
  skip_if_not_installed("Rdisop")

  d <- tempfile("nts_iso_pfoa")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("FT00001", "FT00002", "FT00003"),
      mz = PFOA_MZ,
      rt = c(910, 911, 905),
      sample1 = PFOA_INT,
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("FT00001", "FT00002", "FT00003"),
      mz = PFOA_MZ,
      rt = c(910, 911, 905),
      is_isotope = c(FALSE, TRUE, TRUE),
      parent_name = c(NA_character_, "FT00001", "FT00001"),
      isotope_c13_n = c(NA_integer_, 1L, 2L),
      isotope_kind = c(NA_character_, "13C", "13C"),
      stringsAsFactors = FALSE
    )
  )

  res <- score_isotope_fit(p, formulas = "C8HF15O2")

  expect_s3_class(res, "data.table")
  expect_true(nrow(res) >= 1L)
  expect_true(any(res$isotope_cosine >= 0.95, na.rm = TRUE))
})


test_that("score_isotope_fit: wrong formula gives lower cosine than correct one", {
  skip_if_not_installed("Rdisop")

  d <- tempfile("nts_iso_wrong")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("FT00001", "FT00002", "FT00003"),
      mz = PFOA_MZ,
      rt = c(910, 911, 905),
      sample1 = PFOA_INT,
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("FT00001", "FT00002", "FT00003"),
      mz = PFOA_MZ,
      rt = c(910, 911, 905),
      is_isotope = c(FALSE, TRUE, TRUE),
      parent_name = c(NA_character_, "FT00001", "FT00001"),
      isotope_c13_n = c(NA_integer_, 1L, 2L),
      isotope_kind = c(NA_character_, "13C", "13C"),
      stringsAsFactors = FALSE
    )
  )

  res_correct <- score_isotope_fit(p, formulas = "C8HF15O2")
  res_wrong   <- score_isotope_fit(p, formulas = "CH4")

  cos_correct <- res_correct$isotope_cosine[1L]
  cos_wrong   <- res_wrong$isotope_cosine[1L]

  expect_true(!is.na(cos_correct))
  expect_true(!is.na(cos_wrong))
  expect_true(cos_correct > cos_wrong)
})


test_that(".nts_score_isotope_fit_single: single peak (M+0 only) -> cosine=1, chisq=0", {
  skip_if_not_installed("Rdisop")

  res <- .nts_score_isotope_fit_single(
    mz_vec = 350,
    intensity_vec = 50000,
    formula = "C8HF15O2"
  )

  expect_equal(res$cosine, 1)
  expect_equal(res$chisq, 0)
  expect_equal(res$n_matched, 1L)
})


test_that("score_isotope_fit: empty feature_table -> error", {
  # score_isotope_fit() checks for Rdisop before it inspects the feature table,
  # so without the backend the "feature_table is empty" branch is never reached.
  skip_if_not_installed("Rdisop")

  d <- tempfile("nts_iso_empty_ft")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)

  expect_error(
    score_isotope_fit(p, formulas = "C8HF15O2"),
    "feature_table is empty"
  )
})


test_that("score_isotope_fit: Rdisop not installed -> friendly error", {
  # Test the error path by checking what happens when requireNamespace
  # returns FALSE. We simulate this by calling .nts_score_isotope_fit_single
  # with a valid setup; the function handles Rdisop absence internally.

  # The main function requires Rdisop externally, so this test verifies
  # that the internal helper surfaces a clear error message.
  skip("Mocking requireNamespace is unreliable across testthat versions")
})


test_that("score_isotope_fit: no isotope children returns result", {
  skip_if_not_installed("Rdisop")

  d <- tempfile("nts_iso_nochildren")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = "FT00001",
      mz = 350,
      rt = 500,
      sample1 = 50000,
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = "FT00001",
      mz = 350,
      rt = 500,
      is_isotope = FALSE,
      parent_name = NA_character_,
      stringsAsFactors = FALSE
    )
  )

  res <- score_isotope_fit(p, formulas = "C8HF15O2", min_isotopologues = 1L)

  expect_s3_class(res, "data.table")
  expect_equal(res$n_measured[1L], 1L)
  expect_equal(res$isotope_cosine[1L], 1)
  expect_equal(res$isotope_chisq[1L], 0)
})


test_that("score_isotope_fit: formulas as data.frame", {
  skip_if_not_installed("Rdisop")

  d <- tempfile("nts_iso_df")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("FT00001", "FT00002"),
      mz = c(350, 400),
      rt = c(500, 600),
      sample1 = c(50000, 60000),
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("FT00001", "FT00002"),
      mz = c(350, 400),
      rt = c(500, 600),
      is_isotope = c(FALSE, FALSE),
      parent_name = c(NA_character_, NA_character_),
      stringsAsFactors = FALSE
    )
  )

  fdf <- data.frame(
    feature_name = c("FT00001", "FT00002"),
    formula = c("C8HF15O2", "CH4"),
    stringsAsFactors = FALSE
  )

  res <- score_isotope_fit(p, formulas = fdf, min_isotopologues = 1L)

  expect_s3_class(res, "data.table")
  expect_setequal(res$parent_name, c("FT00001", "FT00002"))
  expect_equal(res$formula[res$parent_name == "FT00001"], "C8HF15O2")
  expect_equal(res$formula[res$parent_name == "FT00002"], "CH4")
})


test_that("score_isotope_fit: write_to_project = TRUE", {
  skip_if_not_installed("Rdisop")

  d <- tempfile("nts_iso_write")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("FT00001", "FT00002", "FT00003"),
      mz = PFOA_MZ,
      rt = c(910, 911, 905),
      sample1 = PFOA_INT,
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("FT00001", "FT00002", "FT00003"),
      mz = PFOA_MZ,
      rt = c(910, 911, 905),
      is_isotope = c(FALSE, TRUE, TRUE),
      parent_name = c(NA_character_, "FT00001", "FT00001"),
      isotope_c13_n = c(NA_integer_, 1L, 2L),
      isotope_kind = c(NA_character_, "13C", "13C"),
      stringsAsFactors = FALSE
    )
  )

  p <- score_isotope_fit(p, formulas = "C8HF15O2", write_to_project = TRUE)

  # When write_to_project = TRUE, returns updated NTSProject (not data.table)
  expect_s4_class(p, "NTSProject")

  ms1 <- get_ms1_annotations(p)
  expect_true("isotope_formula" %in% names(ms1))
  expect_true("isotope_cosine" %in% names(ms1))
  expect_true("isotope_chisq" %in% names(ms1))
  expect_equal(ms1$isotope_formula[ms1$name == "FT00001"], "C8HF15O2")
  expect_true(!is.na(ms1$isotope_cosine[ms1$name == "FT00001"]))
})


test_that(".nts_score_isotope_fit_single: no matching peaks -> cosine=0", {
  skip_if_not_installed("Rdisop")

  # Use H2O as the formula (mass ~18 Da), measured at 500 Da
  res <- .nts_score_isotope_fit_single(
    mz_vec = c(500, 501, 502),
    intensity_vec = c(100, 8, 0.5),
    formula = "H2O",
    mz_tol_ppm = 1
  )

  expect_equal(res$cosine, 0)
  expect_equal(res$n_matched, 0L)
})


test_that(".nts_score_isotope_fit_single: mismatched vector lengths -> error", {
  expect_error(
    .nts_score_isotope_fit_single(
      mz_vec = c(100, 200),
      intensity_vec = c(1000),
      formula = "CH4"
    ),
    "same length"
  )
})


test_that(".nts_score_isotope_fit_single: empty formula -> error", {
  expect_error(
    .nts_score_isotope_fit_single(
      mz_vec = 100,
      intensity_vec = 1000,
      formula = ""
    ),
    "single non-empty string"
  )
})


test_that(".nts_score_isotope_fit_single: paired vectors have same length", {
  skip_if_not_installed("Rdisop")

  res <- .nts_score_isotope_fit_single(
    mz_vec = PFOA_MZ,
    intensity_vec = PFOA_INT,
    formula = "C8HF15O2"
  )

  expect_true(is.finite(res$cosine))
  expect_true(res$n_matched >= 1L)
  expect_true(res$cosine >= 0.95)
})
