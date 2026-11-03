# Regression tests for the 0.1.33 RT / logIE prediction review
# (see docs/analysis/NTSminer_RT与logIE预测审阅_0.1.32.md).
#
# Each test pins one fix and is written to fail on the pre-fix code. Where a fix
# lives behind a modelling backend, the test registers a stub S3 `predict`
# method so the *contract* is observable without training anything (neither
# xgboost nor ranger is required).

# ---- P1-1: the RT feature transform is fitted on the training rows only ------

test_that("RT preprocessing is fitted on the training rows only (P1-1)", {
  set.seed(11)
  n <- 40L; p <- 12L
  X <- matrix(rnorm(n * p), nrow = n, ncol = p,
              dimnames = list(NULL, paste0("f", seq_len(p))))
  tr <- 1:20; te <- 21:40
  X[tr, 1L] <- 5                    # constant within train, varied overall
  X[te, 1L] <- rnorm(20L)

  prep <- NTSminer:::.nts_rt_preprocess(X, fit_idx = tr)
  old  <- NTSminer:::.nts_rt_preprocess(X)      # fit on everything (pre-fix)

  # Column 1 is near-zero-variance *within the training rows*, so it must be
  # screened out; a full-matrix fit cannot see that (the test rows vary).
  expect_true(1L %in% prep$nzv_idx)
  expect_false(1L %in% old$nzv_idx)

  kept <- setdiff(seq_len(p), prep$removed_idx)
  expect_equal(ncol(prep$matrix), length(kept))
  expect_equal(colnames(prep$matrix), colnames(X)[kept])

  # Centre / scale are the *training* constants, and the returned matrix is every
  # row pushed through that transform.
  expect_equal(prep$center, unname(apply(X[tr, kept, drop = FALSE], 2L, mean)),
               tolerance = 1e-10)
  manual <- scale(X[, kept, drop = FALSE], center = prep$center, scale = prep$scale)
  expect_equal(unname(prep$matrix), unname(manual), tolerance = 1e-10)

  # The held-out rows must not influence the transform at all.
  X2 <- X
  X2[te, ] <- X2[te, ] * 3 + 10
  prep2 <- NTSminer:::.nts_rt_preprocess(X2, fit_idx = tr)
  expect_equal(prep2$center, prep$center)
  expect_equal(prep2$scale, prep$scale)
  expect_equal(prep2$removed_idx, prep$removed_idx)
})


test_that(".nts_rt_apply_preproc reuses the learned transform (P1-1)", {
  set.seed(2)
  X <- matrix(rnorm(30 * 8), nrow = 30, ncol = 8,
              dimnames = list(NULL, paste0("f", 1:8)))
  prep <- NTSminer:::.nts_rt_preprocess(X, fit_idx = 1:20)

  expect_equal(unname(NTSminer:::.nts_rt_apply_preproc(X, prep)),
               unname(prep$matrix), tolerance = 1e-12)

  Xn  <- X[25:30, , drop = FALSE]
  got <- NTSminer:::.nts_rt_apply_preproc(Xn, prep)
  kept <- setdiff(seq_len(8), prep$removed_idx)
  expect_equal(colnames(got), colnames(X)[kept])
  expect_equal(unname(got),
               unname(scale(Xn[, kept, drop = FALSE],
                            center = prep$center, scale = prep$scale)),
               tolerance = 1e-12)
})


test_that("build_rt_model / build_logIE_model split before they preprocess (P1-1, P1-3)", {
  rt_src <- paste(deparse(body(build_rt_model)), collapse = "\n")
  expect_true(grepl("preproc <- .nts_rt_preprocess(X_raw, fit_idx = train_idx)",
                    rt_src, fixed = TRUE))
  spl <- regexpr("split_obj <- .nts_rt_group_split", rt_src, fixed = TRUE)
  pre <- regexpr("preproc <- .nts_rt_preprocess", rt_src, fixed = TRUE)
  expect_true(spl > 0L && pre > 0L && spl < pre)

  li_src <- paste(deparse(body(build_logIE_model)), collapse = "\n")
  expect_true(grepl("prep <- .nts_rt_preprocess(X, fit_idx = tr_idx)",
                    li_src, fixed = TRUE))
})


# ---- P2-2: missing RT model columns are imputed with the training centre -----

test_that("predict_rt imputes missing columns with the training centre (P2-2)", {
  predict.width_probe2 <- function(object, newdata, ...) {
    structure(list(predictions = as.numeric(rowSums(newdata))))
  }
  registerS3method("predict", "width_probe2", predict.width_probe2,
                   envir = asNamespace("stats"))

  model <- structure(list(
    model = structure(list(), class = "width_probe2"),
    method = "rf",
    fp_types = character(0), fp_nbits = integer(0),
    nzv_removed = 0L, cor_removed = 0L,
    nzv_idx = integer(0), cor_idx = integer(0),
    preproc_center = c(g1 = 0, g2 = 0, k1 = 10),
    preproc_scale  = c(g1 = 1, g2 = 1, k1 = 2),
    kept_colnames = c("g1", "g2", "k1"),
    descriptor_cols = character(0),
    test_metrics = list(RMSE = 1, R2 = 0.5, MAE = 1),
    p_features = 3L, p_original = 3L,
    calib_shift = NULL, calib_scale = NULL, model_version = "0.2.0"
  ), class = "nts_rt_model")

  X <- matrix(c(1, 2), ncol = 2, dimnames = list(NULL, c("g1", "g2")))
  expect_warning(
    pred <- predict_rt(model, smiles = "C", fingerprint_matrix = X),
    "imputed with the training mean"
  )
  # `k1` is filled with its training centre (10), so after standardisation it
  # contributes exactly 0: (1-0)/1 + (2-0)/1 + (10-10)/2 = 3.
  expect_equal(pred$predicted_rt, 3)
  # The pre-fix behaviour filled raw 0 and standardised afterwards, which maps to
  # z = (0-10)/2 = -5 and gives 1 + 2 - 5 = -2.
  expect_false(isTRUE(all.equal(pred$predicted_rt, -2)))
})


# ---- P2-5: intervals are withheld outside the applicability domain -----------

test_that("predict_rt withholds intervals for out-of-domain compounds (P2-5)", {
  predict.const_probe <- function(object, newdata, ...) {
    structure(list(predictions = rep(100, nrow(newdata))))
  }
  registerS3method("predict", "const_probe", predict.const_probe,
                   envir = asNamespace("stats"))

  cols <- c("maccs_1", "maccs_2", "maccs_3", "maccs_4")
  model <- structure(list(
    model = structure(list(), class = "const_probe"),
    method = "rf", fp_types = "maccs", fp_nbits = 4L,
    nzv_removed = 0L, cor_removed = 0L,
    nzv_idx = integer(0), cor_idx = integer(0),
    preproc_center = NULL, preproc_scale = NULL,
    kept_colnames = cols, descriptor_cols = character(0),
    test_metrics = list(RMSE = 10, R2 = 0.5, MAE = 5),
    p_features = 4L, p_original = 4L,
    calib_shift = NULL, calib_scale = NULL, model_version = "0.2.0",
    app_domain = list(fp_bin = rbind(c(1, 1, 0, 0), c(1, 1, 0, 0)),
                      sim_threshold = 0.6)
  ), class = "nts_rt_model")

  Q <- rbind(c(0, 0, 1, 1),      # disjoint from every training fingerprint
             c(1, 1, 0, 0))      # identical to one of them
  colnames(Q) <- cols
  pred <- predict_rt(model, smiles = c("a", "b"), fingerprint_matrix = Q,
                     return_interval = TRUE)

  expect_equal(pred$applicability, c("out_of_domain", "in_domain"))
  expect_true(is.na(pred$lower[1L]) && is.na(pred$upper[1L]))
  expect_equal(pred$upper[2L] - pred$lower[2L], 40)   # 2 * RMSE on each side

  # A model without an applicability domain keeps the interval everywhere.
  no_dom <- model
  no_dom$app_domain <- NULL
  pred2 <- predict_rt(no_dom, smiles = c("a", "b"), fingerprint_matrix = Q,
                      return_interval = TRUE)
  expect_true(all(is.na(pred2$applicability)))
  expect_false(anyNA(pred2$lower))
})


# ---- P1-2: model provenance is reported, not silently absent -----------------

test_that("print.nts_rt_model reports the applicability domain (P1-2)", {
  base <- structure(list(
    model = list(), method = "xgb", fp_types = "maccs", fp_nbits = 166L,
    nzv_removed = 0L, cor_removed = 0L, nzv_idx = integer(0), cor_idx = integer(0),
    preproc_center = NULL, preproc_scale = NULL, kept_colnames = "maccs_1",
    descriptor_cols = character(0),
    test_metrics = list(RMSE = 10, R2 = 0.8, MAE = 5),
    p_features = 1L, p_original = 1L,
    calib_shift = NULL, calib_scale = NULL, model_version = "0.2.0"
  ), class = "nts_rt_model")

  out <- capture.output(print(base))
  expect_true(any(grepl("App. domain: none", out)))

  dom <- base
  dom$app_domain <- list(fp_bin = matrix(0, nrow = 7, ncol = 4), sim_threshold = 0.6)
  out2 <- capture.output(print(dom))
  expect_true(any(grepl("App. domain: Tanimoto >= 0.6", out2)))
  expect_true(any(grepl("7 training fingerprints", out2)))
})


test_that("load_rt_model accepts the current model versions only (P1-2)", {
  m <- structure(list(
    model = list(), method = "xgb", fp_types = "maccs", fp_nbits = 166L,
    test_metrics = list(RMSE = 1, R2 = 0.5, MAE = 1),
    model_version = "0.2.0"
  ), class = "nts_rt_model")

  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)

  saveRDS(m, f)
  expect_silent(load_rt_model(f))

  m$model_version <- "0.1.0"          # legacy format, still readable
  saveRDS(m, f)
  expect_silent(load_rt_model(f))

  m$model_version <- "9.9.9"
  saveRDS(m, f)
  expect_warning(load_rt_model(f), "may not be compatible")

  m$model_version <- NULL
  saveRDS(m, f)
  expect_warning(load_rt_model(f), "may not be compatible")
})


# ---- P1-3: logIE models carry split / standardisation provenance -------------

test_that("logIE S3 models record split and standardisation provenance (P1-3)", {
  m <- NTSminer:::new_nts_logIE_model(
    model = list(), method = "xgb", feature_set = "user",
    feature_cols = c("a", "b"), n_features = 2L, hp = "",
    train_metrics = list(R2 = 0, RMSE = 0, MAE = 0),
    test_metrics  = list(R2 = 0, RMSE = 0, MAE = 0),
    cv_metrics    = list(R2 = NA_real_, RMSE = NA_real_),
    polarity = "positive", n_training = 10L, n_test = 2L,
    .preproc_center = c(a = 0, b = 1), .preproc_scale = c(a = 1, b = 2),
    .split_strategy = "grouped", .n_groups = 3L
  )
  expect_equal(m$preproc_center, c(a = 0, b = 1))
  expect_equal(m$preproc_scale, c(a = 1, b = 2))
  expect_equal(m$split_strategy, "grouped")
  expect_equal(m$n_groups, 3L)

  out <- capture.output(print(m))
  expect_true(any(grepl("grouped", out)))
  expect_true(any(grepl("3 groups", out)))
  expect_true(any(grepl("Scaled", out)))

  # Externally trained models carry no such fields; print stays quiet about them.
  raw <- structure(list(model_type = "xgboost", model = list(),
                        feature_cols = c("a"), n_features = 1L,
                        polarity = "negative"),
                   class = "list")
  up <- NTSminer:::.upgrade_raw_logIE_model(raw)
  expect_null(up$preproc_center)
  expect_null(up$split_strategy)
})


test_that("predict_logIE applies the training standardisation (P1-3)", {
  # predict_logIE gates on the backend *package* being installed; the actual
  # prediction goes through a stub below, so nothing from ranger is used.
  skip_if_not_installed("ranger")

  predict.stub_sums3 <- function(object, newdata, ...) {
    structure(list(predictions = as.numeric(rowSums(newdata))))
  }
  registerS3method("predict", "stub_sums3", predict.stub_sums3,
                   envir = asNamespace("stats"))

  model <- structure(list(
    model = structure(list(), class = "stub_sums3"),
    method = "rf", feature_set = "user", feature_cols = c("a", "b"),
    n_features = 2L, hp = "stub", best_nrounds = NULL,
    train_metrics = list(R2 = 1, RMSE = 0, MAE = 0),
    test_metrics  = list(R2 = 1, RMSE = 0, MAE = 0),
    cv_metrics    = list(R2 = NA_real_, RMSE = NA_real_),
    polarity = "positive", n_training = 10L, n_test = 2L,
    feature_centers = c(a = 0, b = 0),
    preproc_center  = c(a = 1, b = 2),
    preproc_scale   = c(a = 2, b = 4),
    split_strategy  = "grouped", n_groups = 4L,
    training_date = Sys.time(), version = "0.1.0"
  ), class = "nts_logIE_model")

  X <- matrix(c(5, 3, 9, 7), ncol = 2, dimnames = list(NULL, c("a", "b")))
  res <- predict_logIE(model, smiles = c("x", "y"), fingerprint_matrix = X)

  expect_equal(res$predicted_logIE,
               as.numeric(rowSums(scale(X, center = c(1, 2), scale = c(2, 4)))))
  expect_false(isTRUE(all.equal(res$predicted_logIE,
                                as.numeric(rowSums(X)))))
})


# ---- P0: the bundled logIE models' metadata features are addressable ---------

test_that("bundled logIE models declare the metadata features predict_logIE fills (P0)", {
  pos_path <- NTSminer:::.nts_resolve_logIE_model_path("pos")
  neg_path <- NTSminer:::.nts_resolve_logIE_model_path("neg")
  skip_if(is.null(pos_path) || is.null(neg_path), "bundled logIE models absent")

  mp <- load_logIE_model(pos_path)
  mn <- load_logIE_model(neg_path)

  # POS: retention time column is literally "RT" (not "RT_avg"), plus MW and the
  # run-level organic fraction.
  expect_true(all(c("RT", "MW", "Organic_Percentage") %in% mp$feature_cols))
  # NEG: "RT_avg" + MW, no organic fraction (that is a POS-gradient parameter).
  expect_true(all(c("RT_avg", "MW") %in% mn$feature_cols))
  expect_false("Organic_Percentage" %in% mn$feature_cols)
})


test_that(".nts_logIE_meta_matrix names the model's own RT column (P0)", {
  fcols <- c("maccs_1", "RT", "MW", "Organic_Percentage")
  mm <- NTSminer:::.nts_logIE_meta_matrix(fcols, rt = c(1, 2), mw = c(100, 200),
                                          organic_pct = c(50, 60))
  expect_setequal(colnames(mm), c("RT", "MW", "Organic_Percentage"))
  expect_equal(unname(mm[, "RT"]), c(1, 2))

  mm2 <- NTSminer:::.nts_logIE_meta_matrix(c("RT_avg", "MW"), rt = c(3, 4),
                                           mw = c(1, 1))
  expect_true("RT_avg" %in% colnames(mm2))

  # A model with no metadata features gets no metadata block.
  expect_null(NTSminer:::.nts_logIE_meta_matrix("maccs_1", rt = c(1, 2)))

  # Missing inputs omit the column (predict_logIE then warns) rather than
  # inventing a value.
  mm3 <- NTSminer:::.nts_logIE_meta_matrix(fcols, rt = c(1, 2))
  expect_false("MW" %in% colnames(mm3))
})


# ---- P3-1 / P2-7: bundled model paths resolve and configs agree --------------

test_that("bundled RT / logIE model paths resolve (P3-1)", {
  rt_pos <- NTSminer:::.nts_resolve_rt_model_path("pos")
  rt_neg <- NTSminer:::.nts_resolve_rt_model_path("neg")
  li_pos <- NTSminer:::.nts_resolve_logIE_model_path("pos")
  li_neg <- NTSminer:::.nts_resolve_logIE_model_path("neg")

  for (p in list(rt_pos, rt_neg, li_pos, li_neg)) {
    skip_if(is.null(p), "bundled models absent")
    expect_true(file.exists(p))
    expect_match(p, "\\.rds$")
  }
  expect_null(NTSminer:::.nts_resolve_rt_model_path("does_not_exist"))

  # The RT models on disk are loadable and report the format version the
  # resolver's siblings expect.
  m_rt <- load_rt_model(rt_pos)
  expect_s3_class(m_rt, "nts_rt_model")
  expect_true(m_rt$model_version %in% c("0.1.0", "0.2.0"))
})


test_that("each bundled model config.json describes the model beside it (P2-7)", {
  skip_if_not_installed("jsonlite")
  paths <- c(NTSminer:::.nts_resolve_rt_model_path("pos"),
             NTSminer:::.nts_resolve_rt_model_path("neg"),
             NTSminer:::.nts_resolve_logIE_model_path("pos"),
             NTSminer:::.nts_resolve_logIE_model_path("neg"))
  skip_if(any(vapply(paths, is.null, logical(1L))), "bundled models absent")

  for (rds in paths) {
    cfg <- sub("\\.rds$", "_config.json", rds)
    expect_true(file.exists(cfg), info = cfg)
    js <- jsonlite::fromJSON(cfg)
    expect_equal(basename(rds), js$file, info = cfg)
  }

  # The polarity recorded next to a model must match the model itself: the POS
  # config used to name the pubchem-only variant and the NEG config used to
  # describe a different (881-feature) model than the one shipped.
  for (abbr in c("pos", "neg")) {
    p <- NTSminer:::.nts_resolve_logIE_model_path(abbr)
    m <- load_logIE_model(p)
    js <- jsonlite::fromJSON(sub("\\.rds$", "_config.json", p))
    expect_equal(js$polarity, m$polarity)
    expect_equal(as.numeric(js$feature_count), as.numeric(m$n_features))
    expect_equal(as.numeric(js$n_training), as.numeric(m$n_training))
    expect_equal(as.numeric(js$n_test), as.numeric(m$n_test))
    expect_equal(as.numeric(js$test_r2), round(m$test_metrics$R2, 4),
                 tolerance = 5e-4)
  }
})
