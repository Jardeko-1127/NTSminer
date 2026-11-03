# RT model grouped split + applicability domain (finding B5) ------------------
# A naive random split leaks structural analogues across train/test; and a model
# applied outside its training chemistry predicts confidently but meaninglessly.
# These tests cover the grouping, the group-aware split and the domain scoring.

test_that(".nts_rt_tanimoto computes Jaccard similarity and aligns columns", {
  f <- NTSminer:::.nts_rt_tanimoto
  A <- matrix(c(1, 0, 0, 1,
                1, 0, 1, 0,
                0, 1, 0, 1), nrow = 3, byrow = TRUE)
  colnames(A) <- paste0("f", 1:4)

  S <- f(A)
  expect_equal(diag(S), c(1, 1, 1))
  # {f1,f4} vs {f1,f3}: 1 shared / 3 union
  expect_equal(S[1, 2], 1 / 3)
  # {f1,f4} vs {f2,f4}: 1 shared / 3 union
  expect_equal(S[1, 3], 1 / 3)
  # {f1,f3} vs {f2,f4}: no overlap
  expect_equal(S[2, 3], 0)

  # Column order must not matter (aligned by name).
  B <- A[, c("f4", "f1", "f3", "f2"), drop = FALSE]
  expect_equal(f(A, B), S)
})

test_that(".nts_rt_resolve_groups honours a user vector and fills gaps", {
  f <- NTSminer:::.nts_rt_resolve_groups
  g <- f(group = c("a", "a", NA, "b"), n = 4)
  expect_length(g, 4)
  expect_equal(g[1], g[2])
  expect_false(any(is.na(g)))
  expect_equal(length(unique(g)), 3)
})

test_that(".nts_rt_resolve_groups clusters identical fingerprints together", {
  f <- NTSminer:::.nts_rt_resolve_groups
  fp <- rbind(
    c(1, 1, 1, 0, 0, 0),
    c(1, 1, 1, 0, 0, 0),
    c(0, 0, 0, 1, 1, 1),
    c(0, 0, 0, 1, 1, 1)
  )
  colnames(fp) <- paste0("f", 1:6)
  g <- f(n = 4, fp = fp, cluster_similarity = 0.8)
  # Two tight clusters -> two groups.
  expect_equal(length(unique(g)), 2L)
  expect_equal(g[1], g[2])
  expect_equal(g[3], g[4])
  expect_false(g[1] == g[3])
})

test_that(".nts_rt_group_split keeps whole groups on one side", {
  f <- NTSminer:::.nts_rt_group_split
  groups <- c("a", "a", "b", "b", "c", "c")
  sp <- f(groups, n = 6, train_split = 0.5, seed = 1)
  expect_equal(sp$strategy, "grouped")
  expect_true(length(sp$train_idx) > 0)
  expect_true(length(sp$test_idx) > 0)
  expect_setequal(c(sp$train_idx, sp$test_idx), seq_len(6))
  # No group may straddle the split.
  tr_g <- unique(groups[sp$train_idx])
  te_g <- unique(groups[sp$test_idx])
  expect_length(intersect(tr_g, te_g), 0L)
})

test_that(".nts_rt_group_split falls back to a random split without groups", {
  f <- NTSminer:::.nts_rt_group_split
  sp <- f(as.character(seq_len(10)), n = 10, train_split = 0.8, seed = 2)
  expect_true(sp$strategy %in% c("random", "random_single_group"))
  expect_equal(length(sp$train_idx), 8L)
  expect_length(sp$test_idx, 2L)
})

test_that(".nts_rt_check_group rejects a mis-sized group vector", {
  expect_error(
    NTSminer:::.nts_rt_check_group(c("a", "b"), n = 3),
    "one value per input compound"
  )
})

test_that(".nts_rt_build_app_domain captures training fingerprints", {
  f <- NTSminer:::.nts_rt_build_app_domain
  expect_null(f(NULL, c(1, 2)))
  fp <- matrix(c(1, 0, 1, 0, 1, 1), nrow = 3, byrow = TRUE)
  colnames(fp) <- c("a", "b")
  ad <- f(fp, train_idx = c(1, 3), app_min_similarity = 0.7)
  expect_equal(nrow(ad$fp_bin), 2L)
  expect_equal(ad$sim_threshold, 0.7)
  # Binary-ised.
  expect_true(all(ad$fp_bin %in% c(0, 1)))
})

test_that("predict_rt exposes applicability columns (NA without a domain)", {
  # Minimal stand-in model: a registered predict method lets us exercise the
  # applicability branch without training a real forest.
  fake_fit <- structure(list(), class = "nts_rt_fake_fit")
  registerS3method(
    "predict", "nts_rt_fake_fit",
    function(object, newdata, ...) list(predictions = rep(0, nrow(newdata))),
    envir = asNamespace("NTSminer")
  )

  fp_train <- matrix(c(1, 1, 0, 0,
                       0, 0, 1, 1), nrow = 2, byrow = TRUE)
  colnames(fp_train) <- paste0("f", 1:4)
  model <- structure(list(
    model = fake_fit, method = "rf",
    fp_types = "base", descriptor_cols = character(0),
    kept_colnames = colnames(fp_train),
    preproc_center = NULL, preproc_scale = NULL,
    p_original = 4L, p_features = 4L, p_features_final = 4L,
    calib_shift = NULL, calib_scale = NULL,
    test_metrics = list(RMSE = 5, R2 = 0.9, MAE = 3),
    app_domain = NTSminer:::.nts_rt_build_app_domain(fp_train, c(1, 2), 0.6)
  ), class = "nts_rt_model")

  # Query close to training row 1 (Tanimoto 1) and far from both (0).
  fp_query <- matrix(c(1, 1, 0, 0,
                       0, 0, 0, 0), nrow = 2, byrow = TRUE)
  colnames(fp_query) <- paste0("f", 1:4)
  pred <- predict_rt(model, smiles = c("a", "b"),
                     fingerprint_matrix = fp_query)
  expect_true(all(c("max_train_similarity", "applicability") %in% names(pred)))
  expect_equal(pred$max_train_similarity[1], 1)
  expect_equal(pred$applicability[1], "in_domain")
  expect_equal(pred$applicability[2], "out_of_domain")

  # A model carrying no domain still returns the columns as NA.
  model$app_domain <- NULL
  pred2 <- predict_rt(model, smiles = c("a", "b"),
                      fingerprint_matrix = fp_query)
  expect_true(all(is.na(pred2$max_train_similarity)))
  expect_true(all(is.na(pred2$applicability)))
})
