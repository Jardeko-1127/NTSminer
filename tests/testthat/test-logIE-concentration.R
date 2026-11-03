test_that("predict_rie with logIE_model sets rie column", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")

  m <- load_logIE_model(pos_path)
  ft <- data.frame(
    SMILES = c("c1ccccc1", "CC(=O)OC1=CC=CC=C1C(=O)O", "CCCC"),
    mz = c(78, 180, 58), rt = c(600, 450, 200),
    intensity = c(1000, 2000, 500),
    stringsAsFactors = FALSE
  )
  out <- predict_rie(ft, logIE_model = m)
  expect_true("rie" %in% names(out))
  expect_equal(nrow(out), 3L)
  expect_true(all(out$rie > 0))
})

test_that("predict_rie without model falls back to default_rie", {
  ft <- data.frame(
    SMILES = c("c1ccccc1", "CCCC"),
    mz = c(78, 58), rt = c(600, 200),
    intensity = c(1000, 500),
    stringsAsFactors = FALSE
  )
  out <- predict_rie(ft, default_rie = 2.5)
  expect_true("rie" %in% names(out))
  # Without LogP column and no model, defaults to default_rie
  expect_equal(unique(out$rie), 2.5)
})

test_that("predict_rie auto-loads model by polarity", {
  skip_if_not_installed("xgboost")  # in-house models are xgboost-based
  ft <- data.frame(
    SMILES = c("c1ccccc1", "CCCC"),
    mz = c(78, 58), rt = c(600, 200),
    intensity = c(1000, 500),
    stringsAsFactors = FALSE
  )
  out_pos <- predict_rie(ft, polarity = "positive")
  out_neg <- predict_rie(ft, polarity = "negative")
  expect_true("rie" %in% names(out_pos))
  expect_true("rie" %in% names(out_neg))
  # POS and NEG give different RIE
  expect_false(isTRUE(all.equal(out_pos$rie, out_neg$rie)))
})

test_that("predict_rie handles NA / invalid SMILES", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")

  m <- load_logIE_model(pos_path)
  ft <- data.frame(
    SMILES = c("c1ccccc1", NA_character_, "CCCC"),
    mz = c(78, 100, 58), rt = c(600, 400, 200),
    intensity = c(1000, 1500, 500),
    stringsAsFactors = FALSE
  )
  out <- predict_rie(ft, logIE_model = m, default_rie = 1.0)
  expect_equal(nrow(out), 3L)
  expect_true(all(is.finite(out$rie)))
  # NA SMILES gets default_rie
  expect_equal(out$rie[2], 1.0)
})

test_that("predict_rie detects SMILES column automatically", {
  ft <- data.frame(
    CanonicalSMILES = c("c1ccccc1", "CCCC"),
    mz = c(78, 58), rt = c(600, 200),
    intensity = c(1000, 500),
    stringsAsFactors = FALSE
  )
  out <- suppressWarnings(predict_rie(ft, polarity = "positive"))
  expect_true("rie" %in% names(out))
  expect_true(all(out$rie > 0))
})

test_that("predict_rie falls back to LogP when no SMILES", {
  ft <- data.frame(
    LogP = c(1.5, 3.2, -0.5),
    mz = c(78, 180, 58),
    rt = c(600, 450, 200),
    intensity = c(1000, 2000, 500),
    stringsAsFactors = FALSE
  )
  out <- predict_rie(ft, default_rie = 1.0)
  expect_true("rie" %in% names(out))
  # LogP-based: values should vary (not all default)
  expect_false(isTRUE(all.equal(out$rie, rep(1.0, 3))))
})

test_that("calculate_concentration errors without rie column", {
  ft <- data.frame(mz = 78, rt = 600, intensity = 1000)
  expect_error(calculate_concentration(ft), "rie")
})

test_that("NTSProject path works with predict_rie", {
  ft <- data.frame(
    name = c("A", "B"),
    SMILES = c("c1ccccc1", "CCCC"),
    mz = c(78, 58), rt = c(600, 200),
    intensity = c(1000, 500),
    stringsAsFactors = FALSE
  )
  out <- predict_rie(ft, polarity = "positive")
  expect_true("rie" %in% names(out))
})

# ---- regressions: model auto-load path and metadata column naming -----------

test_that(".nts_resolve_logIE_model_path finds the in-house models", {
  # Models live in data/logIE_models/, NOT data/ directly.
  pos <- NTSminer:::.nts_resolve_logIE_model_path("pos")
  neg <- NTSminer:::.nts_resolve_logIE_model_path("neg")
  if (is.null(pos)) skip("in-house POS model not available on this machine")
  expect_true(file.exists(pos))
  expect_match(basename(pos), "^logIE_model_inhouse_pos\\.rds$")
  expect_equal(load_logIE_model(pos)$polarity, "positive")
  if (!is.null(neg)) {
    expect_true(file.exists(neg))
    expect_equal(load_logIE_model(neg)$polarity, "negative")
  }
})

test_that(".nts_logIE_meta_matrix names the RT column as the model expects", {
  # Regression: the RT column was hard-coded to "RT_avg", but the POS model's
  # feature_cols uses "RT" -> it was renamed and the model's "RT" feature was
  # silently filled with 0.
  mm <- NTSminer:::.nts_logIE_meta_matrix(c("RT", "MW"), rt = 123, mw = 58)
  expect_identical(colnames(mm), c("RT", "MW"))
  expect_equal(unname(mm[1L, "RT"]), 123)
  expect_equal(unname(mm[1L, "MW"]), 58)

  mm2 <- NTSminer:::.nts_logIE_meta_matrix(c("RT_avg", "Organic_Percentage"),
                                           rt = 77, organic_pct = 0.4)
  expect_identical(colnames(mm2), c("RT_avg", "Organic_Percentage"))
  expect_equal(unname(mm2[1L, "RT_avg"]), 77)
  expect_false("RT" %in% colnames(mm2))

  # Columns the model does not use are omitted; NULL inputs are skipped.
  expect_null(NTSminer:::.nts_logIE_meta_matrix(c("MW"), rt = 123))
  expect_null(NTSminer:::.nts_logIE_meta_matrix(c("RT", "MW")))
})
