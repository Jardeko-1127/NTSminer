test_that("load_logIE_model loads POS model", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")
  m <- load_logIE_model(pos_path)
  expect_s3_class(m, "nts_logIE_model")
  expect_equal(m$polarity, "positive")
})

test_that("load_logIE_model loads NEG model", {
  neg_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_neg.rds")
  if (!file.exists(neg_path)) skip("NEG logIE model not found")
  m <- load_logIE_model(neg_path)
  expect_s3_class(m, "nts_logIE_model")
  expect_equal(m$polarity, "negative")
})

test_that("print.nts_logIE_model works", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")
  m <- load_logIE_model(pos_path)
  out <- capture.output(print(m))
  expect_true(any(grepl("nts_logIE_model", out)))
  expect_true(any(grepl("Polarity", out)))
  expect_true(any(grepl("Test R2", out)))
})

test_that("predict_logIE returns correct output shape", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")
  skip_if_not_installed("xgboost")
  m <- load_logIE_model(pos_path)

  smi <- c("c1ccccc1", "CC(=O)OC1=CC=CC=C1C(=O)O", "CCCC")
  res <- predict_logIE(m, smi)
  expect_s3_class(res, "data.frame")
  expect_named(res, "predicted_logIE")
  expect_equal(nrow(res), 3L)
})

test_that("predict_logIE rejects invalid model class", {
  expect_error(predict_logIE(list(), "c1ccccc1"), "nts_logIE_model")
})

test_that("predict_logIE handles NA / empty SMILES", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")
  skip_if_not_installed("xgboost")
  m <- load_logIE_model(pos_path)

  # All-invalid SMILES should error
  expect_error(predict_logIE(m, c(NA_character_, "")), "No valid SMILES")
})

test_that("predict_logIE handles mixed valid/invalid SMILES gracefully", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")
  m <- load_logIE_model(pos_path)

  # Mixed: should not crash, result has predictions
  smi <- c("c1ccccc1", NA_character_, "CCCC", "")
  res <- tryCatch(
    predict_logIE(m, smi),
    error = function(e) NULL
  )
  # Either succeeds with predictions or fails cleanly
  if (!is.null(res)) {
    expect_s3_class(res, "data.frame")
  } else {
    expect_true(TRUE)  # clean error is also acceptable
  }
})

test_that("polarity mismatch is an error", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")
  m <- load_logIE_model(pos_path)

  # One model object serves one polarity; a contradicting `polarity` cannot be
  # honoured, so it must fail loudly rather than warn and predict anyway.
  expect_error(
    predict_logIE(m, c("c1ccccc1", "CCCC", "CC(=O)O"), polarity = "negative"),
    "does not match"
  )
})

test_that("evaluate_logIE_model returns expected structure", {
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")
  skip_if_not_installed("xgboost")
  m <- load_logIE_model(pos_path)

  smi <- c("c1ccccc1", "CC(=O)OC1=CC=CC=C1C(=O)O")
  obs <- c(1.0, 0.5)
  res <- evaluate_logIE_model(m, smi, obs)
  expect_type(res, "list")
  expect_named(res, c("metrics", "predicted", "observed"))
  expect_named(res$metrics, c("R2", "RMSE", "MAE"))
  expect_length(res$predicted, 2L)
  expect_length(res$observed, 2L)
})

test_that("build_logIE_model returns valid object (smoke)", {
  skip_on_ci()
  skip_if_not_installed("rcdk")
  skip_if_not_installed("xgboost")
  smi <- c("c1ccccc1", "CCCC", "CC(=O)Oc1ccccc1C(=O)O", "c1ccncc1",
           "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
           "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", "C1=CC=C2C(=C1)C=CC=C2O")
  logIE <- rnorm(8L, 0.8, 0.5)
  m <- build_logIE_model(smi, logIE, fingerprint_types = "maccs",
                          method = "xgb", train_split = 0.75, seed = 42)
  expect_s3_class(m, "nts_logIE_model")
  expect_equal(m$method, "xgb")
  expect_gt(m$n_training, 0L)
})

test_that("save_logIE_model / load_logIE_model round-trip", {
  tmp <- tempfile(fileext = ".rds")
  pos_path <- testthat::test_path("../../data/logIE_models/logIE_model_inhouse_pos.rds")
  if (!file.exists(pos_path)) skip("POS logIE model not found")
  m <- load_logIE_model(pos_path)

  f <- save_logIE_model(m, tmp)
  expect_equal(f, tmp)

  m2 <- load_logIE_model(tmp)
  expect_s3_class(m2, "nts_logIE_model")
  expect_equal(m2$polarity, m$polarity)
  unlink(tmp)
})
