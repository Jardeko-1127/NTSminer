# tests/testthat/test-forward-spectrum-prediction.R
# Unit tests for FraG-FIORA Forward MS/MS Prediction, Rescoring, and In-Silico Annotation

test_that(".nts_match_spectral_peaks calculates exact cosine and explained intensity", {
  exp_df <- data.frame(
    mz = c(68.9952, 112.9856, 168.9880),
    intensity = c(100.0, 50.0, 25.0)
  )

  # Exact match
  pred_exact <- data.frame(
    mz = c(68.9950, 112.9860, 168.9875),
    intensity = c(100.0, 50.0, 25.0)
  )

  res_exact <- .nts_match_spectral_peaks(exp_df, pred_exact, mz_tol = 0.02)
  expect_equal(res_exact$matched_peaks, 3L)
  expect_gt(res_exact$cosine_sqrt, 0.99)
  expect_gt(res_exact$cosine_linear, 0.99)
  expect_equal(res_exact$explained_intensity, 1.0)

  # Partial match
  pred_partial <- data.frame(
    mz = c(68.9950, 200.0000),
    intensity = c(100.0, 80.0)
  )
  res_partial <- .nts_match_spectral_peaks(exp_df, pred_partial, mz_tol = 0.02)
  expect_equal(res_exact$matched_peaks, 3L)
  expect_equal(res_partial$matched_peaks, 1L)
  expect_gt(res_partial$cosine_sqrt, 0.40)
  expect_lt(res_partial$explained_intensity, 1.0)

  # Completely disjoint
  pred_none <- data.frame(
    mz = c(300.0, 400.0),
    intensity = c(100.0, 50.0)
  )
  res_none <- .nts_match_spectral_peaks(exp_df, pred_none, mz_tol = 0.02)
  expect_equal(res_none$matched_peaks, 0L)
  expect_equal(res_none$cosine_sqrt, 0.0)
  expect_equal(res_none$explained_intensity, 0.0)
})

test_that("predict_ms2_spectrum degrades gracefully when Python is unavailable", {
  # Fake python path to verify graceful degradation
  expect_warning(
    res_fake <- predict_ms2_spectrum("C(=O)(C(F)(F)F)O", python_path = "non_existent_py_binary_12345.exe"),
    "Python with RDKit is required"
  )
  expect_s3_class(res_fake, "nts_predicted_spectrum")
  expect_s3_class(res_fake, "data.frame")
  expect_equal(nrow(res_fake), 0L)
  expect_true(all(c("mz", "intensity", "formula", "annotation", "type") %in% names(res_fake)))
})

test_that("predict_ms2_spectrum runs on real molecules when Python is available", {
  py_exe <- .nts_resolve_forward_python()
  skip_if(is.null(py_exe), "Python with RDKit is not available on this environment.")

  # 1. Single PFAS molecule (TFA, automatic routing to pfas_inhouse)
  tfa_smi <- "C(=O)(C(F)(F)F)O"
  res_tfa <- predict_ms2_spectrum(tfa_smi, nce = "20", polarity = "negative", model = "auto")

  expect_s3_class(res_tfa, "nts_predicted_spectrum")
  expect_s3_class(res_tfa, "data.frame")
  expect_gt(nrow(res_tfa), 0L)
  expect_true(all(c("mz", "intensity", "formula", "annotation", "type") %in% names(res_tfa)))
  expect_equal(attr(res_tfa, "model"), "pfas_inhouse")
  expect_equal(attr(res_tfa, "polarity"), "negative")
  expect_equal(attr(res_tfa, "nce"), "20")

  # Contains characteristic TFA / CF3 fragments (e.g. 68.995 Da CF3- or 112.985 Da M-H)
  expect_true(any(abs(res_tfa$mz - 68.995) < 0.05))
  expect_true(any(abs(res_tfa$mz - 112.985) < 0.05))

  # Test print method executes without error
  expect_output(print(res_tfa), "FraG-FIORA In-Silico Predicted MS/MS Spectrum")

  # 2. General organic molecule (Sulfadimethoxine, no fluorine -> routes to general model)
  sulfa_smi <- "COc1cc(nc(n1)NS(=O)(=O)c2ccc(cc2)N)OC"
  res_sulfa <- predict_ms2_spectrum(sulfa_smi, nce = "20", polarity = "negative", model = "auto")
  expect_s3_class(res_sulfa, "nts_predicted_spectrum")
  expect_gt(nrow(res_sulfa), 0L)
  expect_equal(attr(res_sulfa, "model"), "general")

  # 3. Vector of multiple SMILES
  smis <- c(TFA = tfa_smi, Sulfa = sulfa_smi)
  res_batch <- predict_ms2_spectrum(smis, nce = "Merged", polarity = "negative")
  expect_type(res_batch, "list")
  expect_equal(length(res_batch), 2L)
  expect_equal(names(res_batch), c("TFA", "Sulfa"))
  expect_s3_class(res_batch$TFA, "nts_predicted_spectrum")
  expect_s3_class(res_batch$Sulfa, "nts_predicted_spectrum")

  # 4. Top-k filtering and relative intensity threshold
  res_top3 <- predict_ms2_spectrum(tfa_smi, nce = "20", polarity = "negative", top_k = 3L)
  expect_lte(nrow(res_top3), 3L)
})

test_that("rescore_candidates_by_predicted_spectrum correctly prioritizes true structure", {
  py_exe <- .nts_resolve_forward_python()
  skip_if(is.null(py_exe), "Python with RDKit is not available on this environment.")

  # Experimental spectrum resembling TFA (CF3- and [M-H]-)
  exp_tfa <- data.frame(
    mz = c(68.9952, 112.9855),
    intensity = c(100.0, 65.0)
  )

  cands <- c(
    TFA = "C(=O)(C(F)(F)F)O",
    Acetic_Acid = "CC(=O)O",
    Formic_Acid = "C(=O)O",
    PFOA = "C(=O)(C(C(C(C(C(C(C(F)(F)F)(F)(F))(F)(F))(F)(F))(F)(F))(F)(F))(F)(F))O"
  )

  ranked <- rescore_candidates_by_predicted_spectrum(
    exp_ms2 = exp_tfa,
    candidate_smiles = cands,
    polarity = "negative",
    nce = "20"
  )

  expect_s3_class(ranked, "data.frame")
  expect_equal(nrow(ranked), 4L)
  expect_true(all(c("candidate_id", "smiles", "rank", "score_pred", "cosine_sqrt", "explained_intensity") %in% names(ranked)))

  # TFA must rank #1 among these candidates
  expect_equal(ranked$candidate_id[1L], "TFA")
  expect_equal(ranked$rank[1L], 1L)
  expect_gt(ranked$score_pred[1L], ranked$score_pred[2L])
  expect_gt(ranked$cosine_sqrt[1L], 0.35)
  expect_gt(ranked$explained_intensity[1L], 0.80)
})

test_that("annotate_by_insilico_prediction works on NTSProject and feature tables", {
  py_exe <- .nts_resolve_forward_python()
  skip_if(is.null(py_exe), "Python with RDKit is not available on this environment.")

  # Build mock NTSProject
  proj <- new("NTSProject")

  ft <- data.frame(
    name = c("FT001", "FT002"),
    mz = c(112.9856, 300.1234), # FT001 = TFA [M-H]-
    rt = c(120.0, 350.0),
    stringsAsFactors = FALSE
  )
  proj <- set_feature_table(proj, ft)

  # Mock experimental MS2 for FT001
  clean_ms2_mock <- list(
    FT001 = list(
      spectrum = data.frame(
        mz = c(68.9952, 112.9855),
        intensity = c(100.0, 65.0)
      )
    )
  )
  proj <- set_clean_ms2(proj, clean_ms2_mock)

  # Mock suspect table containing TFA
  suspect_df <- data.frame(
    ID = c("SUS_TFA", "SUS_OTHER"),
    Name = c("Trifluoroacetic acid", "Hypothetical compound"),
    SMILES = c("C(=O)(C(F)(F)F)O", "CCCCCCCC(=O)O"),
    Formula = c("C2HF3O2", "C8H16O2"),
    ExactMass = c(113.9929, 144.1150),
    stringsAsFactors = FALSE
  )

  # Run in-silico annotation
  proj_ann <- annotate_by_insilico_prediction(
    object = proj,
    suspect_table = suspect_df,
    polarity = "negative",
    nce = "20",
    mz_tol_ms1 = 15,
    min_cosine = 0.30,
    verbose = FALSE
  )

  expect_s4_class(proj_ann, "NTSProject")
  db_ann <- proj_ann@db_annotations
  expect_s3_class(db_ann, "data.frame")
  expect_gt(nrow(db_ann), 0L)

  hit_tfa <- db_ann[db_ann$name == "FT001", , drop = FALSE]
  expect_equal(nrow(hit_tfa), 1L)
  expect_equal(hit_tfa$ref_id, "SUS_TFA")
  expect_equal(hit_tfa$db_type, "prediction")
  expect_gt(hit_tfa$score_ms2, 0.35)

  # Check project log was recorded
  expect_true(any(vapply(proj_ann@log, function(x) x$step == "annotate_by_insilico_prediction", logical(1L))))

  # Existing L1/L2a annotation protection test
  # If FT001 already has an L1 row, skip_annotated = TRUE must protect it
  l1_mock <- data.frame(
    name = "FT001",
    ref_id = "INHOUSE_001",
    db_type = "inhouse",
    Compound.name = "TFA Confirmed Standard",
    SMILES = "C(=O)(C(F)(F)F)O",
    score_total = 0.99,
    score_ms2 = 0.98,
    stringsAsFactors = FALSE
  )
  proj_l1 <- proj
  methods::slot(proj_l1, "db_annotations") <- l1_mock

  proj_l1_ann <- annotate_by_insilico_prediction(
    object = proj_l1,
    suspect_table = suspect_df,
    skip_annotated = TRUE,
    verbose = FALSE
  )
  db_ann_prot <- proj_l1_ann@db_annotations
  # FT001 remains solely the in-house hit
  expect_equal(nrow(db_ann_prot), 1L)
  expect_equal(db_ann_prot$db_type, "inhouse")
})

test_that("plot_predicted_mirror produces valid ggplot object", {
  py_exe <- .nts_resolve_forward_python()
  skip_if(is.null(py_exe), "Python with RDKit is not available on this environment.")

  exp_tfa <- data.frame(
    mz = c(68.9952, 112.9855),
    intensity = c(100.0, 65.0)
  )

  p <- plot_predicted_mirror(
    exp_spectrum = exp_tfa,
    smiles = "C(=O)(C(F)(F)F)O",
    nce = "20",
    polarity = "negative"
  )

  expect_s3_class(p, "ggplot")
  expect_true(length(p$layers) >= 3L) # hline, exp segments, pred segments
})
