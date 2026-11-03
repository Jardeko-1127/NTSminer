# Regression tests for the 0.1.17 P1 round (see the whole-pipeline
# evaluation report under docs/analysis/, sections 14-15).
#
# Each test pins one of the eight fixes and is written to fail on the pre-fix
# code. Where a fix lives in a code path that needs a modelling backend, the
# test uses a stub S3 method so the *contract* is checked without training
# anything; the xgboost early-stopping control flow (7-2) can only be observed
# with an instrumented backend and is covered by the standalone reproduction in
# tests/rawdata_pfas/verify_p1_fixes_0.1.17.R.

# ---- 7-1: predict_rt must not re-apply original-space pruning indices --------

test_that("predict_rt keeps the post-preprocessing width (7-1)", {
  # A stub backend that reports how many columns it was handed. The model below
  # carries nzv_idx / cor_idx that index the *original* matrix; before the fix
  # predict_rt applied them to the already-pruned matrix, so the stub would have
  # seen 3 columns instead of 6 (and real xgboost would have errored outright).
  predict.width_probe <- function(object, newdata, ...) {
    structure(list(predictions = rep(ncol(newdata), nrow(newdata))))
  }
  registerS3method("predict", "width_probe", predict.width_probe,
                   envir = asNamespace("stats"))

  X6 <- matrix(seq_len(60 * 6) / 60, ncol = 6,
               dimnames = list(NULL, paste0("f", 1:6)))
  model <- structure(list(
    model = structure(list(), class = "width_probe"),
    method = "rf",
    fp_types = character(0), fp_nbits = integer(0),
    nzv_removed = 2L, cor_removed = 1L,
    nzv_idx = c(1L, 2L), cor_idx = 3L,
    preproc_center = NULL, preproc_scale = NULL,
    kept_colnames = colnames(X6),
    descriptor_cols = character(0),
    test_metrics = list(RMSE = 1, R2 = 0.5, MAE = 1),
    p_features = 6L, p_original = 6L,
    calib_shift = NULL, calib_scale = NULL, model_version = "0.2.0"
  ), class = "nts_rt_model")

  pred <- predict_rt(model, smiles = rep("C", 3),
                     fingerprint_matrix = X6[1:3, , drop = FALSE])
  expect_equal(pred$predicted_rt, rep(6, 3))

  # The bug itself, inlined: the stale indices do collapse 6 columns to 3.
  rm_idx <- unique(c(model$nzv_idx, model$cor_idx))
  expect_equal(ncol(X6[, -rm_idx, drop = FALSE]), 3L)

  expect_false(grepl("nzv_idx", paste(deparse(body(predict_rt)), collapse = "\n")))
  expect_false(grepl("cor_idx", paste(deparse(body(predict_rt)), collapse = "\n")))
})

# ---- 7-2: early stopping must not watch the held-out test set ----------------

test_that("build_logIE_model early-stops on an inner validation split (7-2)", {
  src <- paste(deparse(body(build_logIE_model)), collapse = "\n")
  expect_true(grepl("valid = dva", src, fixed = TRUE))
  expect_false(grepl("test = dte", src, fixed = TRUE))

  skip_if_not_installed("xgboost")
  set.seed(7)
  X <- matrix(rnorm(80 * 5), ncol = 5, dimnames = list(NULL, paste0("f", 1:5)))
  y <- X[, 1] * 1.5 + rnorm(80) * 0.3

  m <- build_logIE_model(logIE = y, method = "xgb", train_split = 0.8,
                         seed = 11, fingerprint_matrix = X)

  expect_s3_class(m, "nts_logIE_model")
  expect_equal(m$n_training + m$n_test, 80L)
  expect_true(is.finite(m$best_nrounds))
  expect_gte(m$best_nrounds, 1L)
  expect_lte(m$best_nrounds, 500L)
  expect_true(is.finite(m$test_metrics$RMSE))
  expect_true(is.finite(m$test_metrics$R2))
})

# ---- 7-3: missing feature columns are imputed with training means -----------

test_that("predict_logIE imputes missing columns with training means, not 0 (7-3)", {
  skip_if_not_installed("ranger")

  # Stub backend: predictions depend linearly on every column, so the exact
  # imputed value is observable from the output.
  predict.stub_sums <- function(object, newdata, ...) {
    structure(list(predictions = as.numeric(rowSums(newdata))))
  }
  registerS3method("predict", "stub_sums", predict.stub_sums,
                   envir = asNamespace("stats"))

  model <- structure(list(
    model = structure(list(), class = "stub_sums"),
    method = "rf",
    feature_set = "user",
    feature_cols = c("a", "b", "c"),
    n_features = 3L,
    hp = "stub",
    best_nrounds = NULL,
    train_metrics = list(R2 = 1, RMSE = 0, MAE = 0),
    test_metrics = list(R2 = 1, RMSE = 0, MAE = 0),
    cv_metrics = list(R2 = NA_real_, RMSE = NA_real_),
    polarity = "positive",
    n_training = 10L, n_test = 2L,
    feature_centers = c(a = 1, b = 2, c = 3),
    training_date = Sys.time(),
    version = "0.1.0"
  ), class = "nts_logIE_model")

  Xab <- matrix(c(10, 20, 30, 40), ncol = 2,
                dimnames = list(NULL, c("a", "b")))   # column `c` is missing

  expect_warning(
    pred <- predict_logIE(model, smiles = c("x", "y"), fingerprint_matrix = Xab),
    "imputed with the training-set mean"
  )
  expect_equal(pred$predicted_logIE, rowSums(cbind(Xab, c = c(3, 3))))
  expect_false(isTRUE(all.equal(pred$predicted_logIE, rowSums(cbind(Xab, c = c(0, 0))))))

  # A model built before the field existed must still work, and say why.
  legacy <- model
  legacy$feature_centers <- NULL
  expect_warning(
    pred_legacy <- predict_logIE(legacy, smiles = c("x", "y"),
                                 fingerprint_matrix = Xab),
    "no training means"
  )
  expect_equal(pred_legacy$predicted_logIE, rowSums(cbind(Xab, c = c(0, 0))))
})

# ---- 8-1: abundance-column pickers must exclude shape/coordinate columns -----

test_that("abundance-column pickers exclude feature shape columns (8-1)", {
  ft <- data.frame(
    name = c("F1", "F2"),
    mz = c(100.1, 200.2), mzmed = c(100.1, 200.2),
    mzmin = c(100.0, 200.1), mzmax = c(100.2, 200.3),
    rt = c(60, 120), rtmed = c(60, 120),
    rtmin = c(55, 115), rtmax = c(65, 125),
    npeaks = c(12L, 9L),
    kmd = c(0.11, 0.22),
    S1 = c(1000, 2000), S2 = c(1500, 2500),
    stringsAsFactors = FALSE
  )
  sample_cols <- c("S1", "S2")

  expect_setequal(NTSminer:::.nts_abundance_colnames(ft), sample_cols)
  expect_setequal(NTSminer:::.smq_abundance_cols(ft), sample_cols)
  expect_setequal(NTSminer:::.nts_isotope_fit_sample_cols(ft), sample_cols)
  expect_equal(unname(NTSminer:::.nts_isotope_row_intensity(ft)),
               unname(rowMeans(ft[, sample_cols])))

  proj <- methods::new("NTSProject")
  proj@feature_table <- ft
  expect_equal(unname(NTSminer:::.nts_isotope_mean_intensity(proj)),
               unname(rowMeans(ft[, sample_cols])))

  shape <- NTSminer:::.nts_feature_shape_cols
  expect_true(all(c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                    "npeaks") %in% shape))

  # The pickers used to carry their own copies of the reserved set and drifted
  # apart; they must all defer to the shared constant.
  ff <- paste(readLines(file.path("..", "..", "R", "04_feature_filter.R")),
              collapse = "\n")
  expect_gte(length(gregexpr(".nts_feature_shape_cols", ff, fixed = TRUE)[[1]]), 3L)
})

# ---- 8-2: GNPS export column names / units / FEATURE_ID ----------------------

test_that("export_to_gnps writes FBMN-conformant column names (8-2)", {
  proj <- methods::new("NTSProject")
  proj@feature_table <- data.frame(
    name = c("F1", "F2"),
    mz = c(100.1, 200.2),
    rt = c(120, 300),                       # seconds (xcms)
    mzmin = c(100.0, 200.1), mzmax = c(100.2, 200.3),
    rtmin = c(115, 295), rtmax = c(125, 305), npeaks = c(5L, 6L),
    S1 = c(1000, 2000), S2 = c(1500, 2500),
    stringsAsFactors = FALSE
  )
  proj@ms2data <- list(
    list(feature_id = "F1", precursorMz = 100.1, rt = 120,
         mz = c(50.1, 60.2), intensity = c(500, 300)),
    list(feature_id = "F2", precursorMz = 200.2, rt = 300,
         mz = c(80.3, 90.4), intensity = c(700, 400))
  )
  out_dir <- file.path(tempdir(), "nts_gnps_p1")
  unlink(out_dir, recursive = TRUE)

  suppressWarnings(export_to_gnps(proj, out_dir))

  qt <- utils::read.csv(file.path(out_dir, "quantification_table.csv"),
                        check.names = FALSE)
  expect_equal(names(qt)[1:3],
               c("row ID", "row m/z", "row retention time"))
  expect_equal(names(qt)[-(1:3)], c("S1", "S2"))
  expect_equal(qt[["row retention time"]], c(2, 5))   # 120 s / 300 s -> minutes

  mgf <- readLines(file.path(out_dir, "spectra.mgf"))
  expect_true(any(grepl("^FEATURE_ID=", mgf)))
  expect_false(any(grepl("^TITLE=", mgf)))
  expect_equal(sub("^FEATURE_ID=", "", grep("^FEATURE_ID=", mgf, value = TRUE)),
               qt[["row ID"]])
})

# ---- 9-1: baseline normalisation / top-k must group by NCE -------------------

test_that("build_nts_database normalises and top-k's per (Lab_ID, Polarity, NCE) (9-1)", {
  md <- data.frame(Lab_ID = "L1", mass = 300, Polarity = "NEG",
                   SMILES = "CCCC", stringsAsFactors = FALSE)
  sl <- data.frame(
    Lab_ID = "L1",
    NCE = c(rep("CE10", 4), rep("CE20", 4)),
    mz = c(100, 150, 200, 250, 100, 150, 200, 250),
    intensity = c(1000, 500, 250, 100, 10, 5, 2.5, 1)
  )
  db <- suppressMessages(
    build_nts_database(md, sl, noise_threshold = 0.01, top_k = 30L,
                       add_merged_layer = FALSE)
  )
  expect_equal(nrow(db@spectra_long[NCE == "CE10"]), 4L)
  expect_equal(nrow(db@spectra_long[NCE == "CE20"]), 4L)

  # The old single-group rule starves the quieter collision energy.
  old <- data.table::as.data.table(sl)
  old[, imax := max(intensity, na.rm = TRUE), by = "Lab_ID"]
  old[, int_rel := intensity / imax]
  expect_lt(nrow(old[int_rel >= 0.01 & NCE == "CE20"]), 4L)
})

# ---- 3-1: in-house polarity must be carried, never silently NEG --------------

test_that("polarity is normalised tolerantly and never silently defaulted (3-1)", {
  pol <- NTSminer:::.nts_normalize_polarity
  expect_equal(pol(c("POS", "POSITIVE", "pos", "+", "Pos-ESI")), rep("POS", 5))
  expect_equal(pol(c("NEG", "NEGATIVE", "neg", "-")), rep("NEG", 4))
  expect_equal(pol(c("", NA), default = "NEG"), c("NEG", "NEG"))
  expect_error(pol("banana"), "positive or negative")

  sl <- data.frame(
    Lab_ID = "L1", NCE = "CE10", mz = c(100, 150), intensity = c(1000, 500)
  )
  md_pos <- data.frame(Lab_ID = "L1", mass = 300, Polarity = "POSITIVE",
                       SMILES = "CCCC", stringsAsFactors = FALSE)
  db_pos <- suppressMessages(
    build_nts_database(md_pos, sl, add_merged_layer = FALSE)
  )
  expect_true(all(db_pos@metadata$Polarity == "POS"))
  expect_true(all(db_pos@spectra_long$Polarity == "POS"))

  expect_error(
    suppressMessages(suppressWarnings(build_nts_database(
      data.frame(Lab_ID = "L1", mass = 300, Polarity = "weird", SMILES = "CCCC"),
      sl, add_merged_layer = FALSE
    ))),
    "positive or negative"
  )
  expect_warning(
    suppressMessages(build_nts_database(
      data.frame(Lab_ID = "L1", mass = 300, SMILES = "CCCC"), sl,
      add_merged_layer = FALSE
    )),
    "no `Polarity` column"
  )

  # extract_mix_standards must emit the column so the warning is unreachable
  # for the in-house workflow. deparse() rewrites `out[, Polarity := pol]` as
  # `` out[, `:=`(Polarity, pol)] ``.
  src <- paste(deparse(body(extract_mix_standards)), collapse = "\n")
  expect_true(grepl("`:=`(Polarity, pol)", src, fixed = TRUE))
})

# ---- 9-3: adduct constants must carry the charge-carrier electron ------------

test_that("adduct deltas account for the electron (9-3)", {
  d <- NTSminer:::.nts_precursor_mz_to_neutral_delta
  e_m <- 0.00054858

  # The proton *mass* already equals m(1H) - m(e): no second subtraction.
  expect_equal(d("M+H", "POS"), -1.00727646, tolerance = 1e-9)
  expect_equal(d("M-H", "NEG"), 1.00727646, tolerance = 1e-9)

  # Atom-derived adducts: cations lose an electron, chloride gains one.
  expect_equal(d("M+Na", "POS"), -(22.98976928 - e_m), tolerance = 1e-9)
  expect_equal(d("M+K", "POS"), -(38.96370668 - e_m), tolerance = 1e-9)
  expect_equal(d("M+NH4", "POS"), -(18.034374 - e_m), tolerance = 1e-9)
  expect_equal(d("M+Cl", "NEG"), -(34.96885268 + e_m), tolerance = 1e-9)

  # Pre-fix values were the neutral-atom masses: ~0.55 mDa too heavy.
  expect_gt(abs(d("M+NH4", "POS") + 18.034374), 5e-4)
  expect_gt(abs(d("M+Na", "POS") + 22.98976928), 1e-6)
})
