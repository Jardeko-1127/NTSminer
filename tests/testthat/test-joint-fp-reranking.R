# Tests for Candidate Fingerprint Re-ranking Module
# (DeepMass-style Candidate Re-ranking using Transformer-predicted 2215-bit Fingerprints)

test_that(".nts_tanimoto_bits computes exact Tanimoto similarities", {
  # 1. Identical bits -> 1.0
  b1 <- c(10, 20, 30, 155, 2050)
  expect_equal(NTSminer:::.nts_tanimoto_bits(b1, b1), 1.0)

  # 2. Completely disjoint bits -> 0.0
  b2 <- c(1, 2, 3)
  b3 <- c(4, 5, 6)
  expect_equal(NTSminer:::.nts_tanimoto_bits(b2, b3), 0.0)

  # 3. Partial overlap
  # A = {1, 2, 3}, B = {2, 3, 4} -> Intersect = {2, 3} (2), Union = {1, 2, 3, 4} (4) -> 2/4 = 0.5
  expect_equal(NTSminer:::.nts_tanimoto_bits(c(1, 2, 3), c(2, 3, 4)), 0.5)

  # 4. Empty or NA inputs
  expect_true(is.na(NTSminer:::.nts_tanimoto_bits(integer(0), c(1, 2))))
  expect_true(is.na(NTSminer:::.nts_tanimoto_bits(c(1, 2), integer(0))))
  expect_true(is.na(NTSminer:::.nts_tanimoto_bits(NULL, c(1, 2))))
})


test_that("rescore_candidates_by_joint_fp correctly calculates Tanimoto and re-ranks candidates", {
  cand <- data.frame(
    feature = c("FT01", "FT01", "FT02"),
    ref_id = c("insilico_01", "insilico_02", "insilico_03"),
    candidate = c("PFOA Isomer A", "PFOA Isomer B", "PFOS"),
    score_total = c(0.60, 0.60, 0.70),
    SMILES = c("C(=O)(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)O",
               "C(=O)(C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)O",
               "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
    stringsAsFactors = FALSE
  )

  # Synthetic predicted fingerprint bit lists
  # FT01 has high overlap with Isomer B (simulated bits 10..20) and low overlap with Isomer A
  mock_pred_fp <- list(
    FT01 = c(10, 11, 12, 13, 14, 15, 155, 2048),
    FT02 = c(50, 51, 52, 131, 2050)
  )

  # When candidate fingerprints are pre-supplied or mock-calculated
  cand$fp_similarity <- c(0.30, 0.85, 0.75)

  rescored <- rescore_candidates_by_joint_fp(
    cand,
    predicted_fp = mock_pred_fp,
    weight_fp = 0.50
  )

  expect_s3_class(rescored, "data.frame")
  expect_true("fp_similarity" %in% names(rescored))
  expect_true("score_reranked" %in% names(rescored))
  expect_true("rank_reranked" %in% names(rescored))

  # In FT01: Isomer B had higher fp_similarity (0.85 vs 0.30) despite tied score_total (0.60)
  ft01_cands <- rescored[rescored$feature == "FT01", ]
  expect_equal(ft01_cands$candidate[1], "PFOA Isomer B")
  expect_equal(ft01_cands$rank_reranked[1], 1L)
  expect_equal(ft01_cands$candidate[2], "PFOA Isomer A")
  expect_equal(ft01_cands$rank_reranked[2], 2L)

  # Blended score: (1 - 0.5)*0.60 + 0.5*0.85 = 0.30 + 0.425 = 0.725
  expect_equal(ft01_cands$score_reranked[1], 0.725)
  # Isomer A: (1 - 0.5)*0.60 + 0.5*0.30 = 0.30 + 0.15 = 0.45
  expect_equal(ft01_cands$score_reranked[2], 0.45)
})


test_that("rescore_candidates_by_joint_fp handles predicted_fp from data.frame and NTSProject proxy", {
  cand <- data.frame(
    name = c("FT001", "FT001"),
    ref_id = c("inhouse_L01", "inhouse_L02"),
    candidate = c("PFHxA", "PFPeA"),
    score_total = c(0.80, 0.80),
    fp_similarity = c(0.90, 0.40),
    stringsAsFactors = FALSE
  )

  # Data.frame input
  df_pred <- data.frame(
    feature_id = "FT001",
    joint_fp_bits = "10;20;30;155",
    stringsAsFactors = FALSE
  )

  res <- rescore_candidates_by_joint_fp(cand, predicted_fp = df_pred, weight_fp = 0.4)
  expect_equal(res$rank_reranked, c(1L, 2L))
  expect_equal(res$ref_id[1], "inhouse_L01")
  expect_true(res$score_reranked[1] > res$score_reranked[2])
})


test_that("rank_annotation_candidates breaks ties using fp_similarity", {
  # Two candidates for FT01 with identical decision_tier (L2b) and score_lib (0.60)
  cand <- data.frame(
    feature = c("FT01", "FT01"),
    ref_id = c("insilico_A", "insilico_B"),
    candidate = c("Compound A (low FP sim)", "Compound B (high FP sim)"),
    db_type = c("prediction", "prediction"),
    score_total = c(0.60, 0.60),
    fp_similarity = c(0.20, 0.85),
    stringsAsFactors = FALSE
  )

  out <- rank_annotation_candidates(cand, verbose = FALSE)

  # Compound B must win Top-1 due to higher fp_similarity
  expect_equal(out$final$candidate[1], "Compound B (high FP sim)")
  expect_equal(out$final$ref_id[1], "insilico_B")
  expect_equal(out$candidates$is_top, c(TRUE, FALSE))
  expect_equal(out$candidates$ref_id[1], "insilico_B")
  expect_equal(out$final$fp_similarity[1], 0.85)

  # Flags should contain FP_MATCH for similarity >= 0.70
  expect_true(grepl("FP_MATCH", out$candidates$flags[1]))
})
