# Tests for Joint 2215-bit Fingerprint & Ecotoxicity Prediction Module
# (MS2 -> Joint 2215-bit Fingerprint -> Downstream QSAR Ecotoxicity)

test_that("load_joint_fp_model_config loads and validates bundled configuration", {
  cfg <- load_joint_fp_model_config()
  expect_s3_class(cfg, "nts_joint_fp_config")
  expect_type(cfg, "list")

  # Architecture
  expect_equal(cfg$architecture$d_model, 256L)
  expect_equal(cfg$architecture$num_heads, 8L)
  expect_equal(cfg$architecture$num_encoder_layers, 4L)
  expect_equal(cfg$architecture$num_decoder_layers, 4L)
  expect_equal(cfg$architecture$vocab_size, 2218L)
  expect_equal(cfg$architecture$num_mz_vocabs, 148003L)

  # Training summary
  expect_equal(cfg$training_summary$training_spectra, 124908L)
  expect_equal(cfg$training_summary$training_scaffolds, 60502L)
  expect_equal(cfg$training_summary$validation_spectra, 20980L)
  expect_true(cfg$training_summary$train_loss < 1.5)
  expect_true(cfg$training_summary$val_tanimoto > 0.65)

  # PFAS benchmark performance
  expect_equal(cfg$pfas_benchmark_performance$isolated_scaffolds, 174L)
  expect_equal(cfg$pfas_benchmark_performance$evaluated_spectra, 3858L)
  expect_true(cfg$pfas_benchmark_performance$joint_2215_tanimoto_mean > 0.60)
  expect_true(cfg$pfas_benchmark_performance$maccs167_tanimoto_mean > 0.80)

  # Downstream toxicity fidelity
  tox <- cfg$pfas_benchmark_performance$downstream_toxicity_fidelity
  expect_true(tox$spearman_rho > 0.70)
  expect_true(tox$pearson_r > 0.75)
  expect_true(tox$mae < 0.20)

  # Print method does not crash
  out_txt <- capture.output(print(cfg))
  expect_true(any(grepl("IDSL_MINT", out_txt) | grepl("Joint 2215-bit", out_txt)))
  expect_true(any(grepl("174 Scaffolds", out_txt)))
})


test_that(".nts_resolve_joint_fp_model_paths resolves bundled assets correctly", {
  paths <- NTSminer:::.nts_resolve_joint_fp_model_paths()
  expect_type(paths, "list")
  expect_true(file.exists(paths$weights))
  expect_true(file.exists(paths$config))
  expect_true(file.exists(paths$tox))
  expect_true(grepl("joint_fp_model_weights", paths$weights))
  expect_true(grepl("joint_fp_model_config.json", paths$config))
  expect_true(grepl("tox_model.joblib", paths$tox))

  script <- NTSminer:::.nts_resolve_joint_fp_engine_script()
  expect_true(file.exists(script))
  expect_true(grepl("joint_fp_engine.py", script))
})


test_that("predict_joint_fingerprint validates inputs and handles missing python gracefully", {
  # 1. Non-existent Python path throws informative error
  fake_py <- file.path(tempdir(), "non_existent_python_binary_xyz.exe")
  expect_error(
    predict_joint_fingerprint(
      matrix(c(100, 200, 50, 100), ncol = 2),
      precursor_mz = 412.9664,
      python_path = fake_py
    ),
    "Python executable not found"
  )

  # 2. Missing precursor_mz for raw matrix throws informative error
  expect_error(
    predict_joint_fingerprint(
      matrix(c(100, 200, 50, 100), ncol = 2),
      precursor_mz = NULL,
      python_path = "some_python"
    )
  )

  # 3. Invalid object type throws error
  expect_error(
    predict_joint_fingerprint(
      object = 12345,
      precursor_mz = 100,
      python_path = "some_python"
    )
  )
})


test_that("predict_joint_fingerprint data exchange and NTSProject update work with mock python", {
  # Create a mock Python script that generates a valid prediction response
  mock_py_script <- tempfile(fileext = ".bat")
  mock_r_worker <- tempfile(fileext = ".R")

  # Write an R worker script to simulate Python CLI output
  cat('
    args <- commandArgs(trailingOnly = TRUE)
    in_idx <- which(args == "--input_json") + 1L
    out_idx <- which(args == "--output_json") + 1L
    in_file <- args[in_idx]
    out_file <- args[out_idx]
    
    in_data <- jsonlite::fromJSON(in_file, simplifyVector = FALSE)
    out_data <- lapply(in_data, function(x) {
      list(
        id = x$id,
        precursor_mz = x$precursor_mz,
        n_total_bits = 32L,
        n_ecfp_bits = 11L,
        n_maccs_bits = 21L,
        joint_fp_on_bits = c(12L, 45L, 2048L, 2055L),
        ecfp_on_bits = c(12L, 45L),
        maccs_on_bits = c(0L, 7L),
        toxicity_score = 0.1726,
        toxicity_risk = "Low",
        detected_functional_groups = list("Fluorine atom (F)", "Carboxylic acid (-COOH) / Carbonyl")
      )
    })
    jsonlite::write_json(out_data, out_file, auto_unbox = TRUE)
  ', file = mock_r_worker)

  # On Windows, wrap with a batch file; on Unix, a shell script
  if (.Platform$OS.type == "windows") {
    r_exe <- file.path(R.home("bin"), "Rscript.exe")
    cat(sprintf('@echo off\n"%s" "%s" %%*\n', r_exe, mock_r_worker), file = mock_py_script)
  } else {
    r_exe <- file.path(R.home("bin"), "Rscript")
    cat(sprintf('#!/bin/sh\n"%s" "%s" "$@"\n', r_exe, mock_r_worker), file = mock_py_script)
    Sys.chmod(mock_py_script, "755")
  }

  pks <- matrix(c(68.995, 118.992, 168.989, 25.0, 100.0, 85.0), ncol = 2)
  colnames(pks) <- c("mz", "intensity")

  # 1. Test data.frame output for single spectrum
  res_df <- predict_joint_fingerprint(
    object = pks,
    precursor_mz = 412.9664,
    python_path = mock_py_script,
    verbose = FALSE
  )

  expect_s3_class(res_df, "data.frame")
  expect_equal(nrow(res_df), 1L)
  expect_equal(res_df$precursor_mz, 412.9664)
  expect_equal(res_df$n_total_bits, 32L)
  expect_equal(res_df$toxicity_score, 0.1726)
  expect_equal(res_df$toxicity_risk, "Low")
  expect_true(grepl("Fluorine atom", res_df$detected_functional_groups))

  # 2. Test predict_toxicity_from_ms2 wrapper
  tox_df <- predict_toxicity_from_ms2(
    object = pks,
    precursor_mz = 412.9664,
    python_path = mock_py_script,
    verbose = FALSE
  )
  expect_equal(tox_df$toxicity_score, 0.1726)

  # 3. Test NTSProject update
  proj <- new("NTSProject")
  proj@feature_table <- data.frame(
    name = c("FT0001", "FT0002"),
    mz = c(412.9664, 227.1078),
    rt = c(300, 150),
    stringsAsFactors = FALSE
  )
  proj@clean_ms2 <- list(
    FT0001 = pks,
    FT0002 = matrix(c(133.065, 212.084, 100.0, 80.0), ncol = 2)
  )

  proj_updated <- predict_joint_fingerprint(
    object = proj,
    python_path = mock_py_script,
    write_to_project = TRUE,
    verbose = FALSE
  )

  expect_s4_class(proj_updated, "NTSProject")
  proxy <- get_physchem_proxy(proj_updated)
  expect_true("MS2_Toxicity_Score" %in% names(proxy))
  expect_true("MS2_Toxicity_Risk" %in% names(proxy))
  expect_true("MS2_Functional_Groups" %in% names(proxy))
  expect_equal(nrow(proxy), 2L)
  expect_equal(proxy$MS2_Toxicity_Score[1], 0.1726)
  expect_equal(proxy$MS2_Toxicity_Risk[1], "Low")

  # Verify project log was appended
  log_entries <- proj_updated@log
  expect_true(any(sapply(log_entries, function(e) e$step == "predict_joint_fingerprint")))

  # Cleanup temporary mock scripts
  unlink(c(mock_py_script, mock_r_worker))
})
