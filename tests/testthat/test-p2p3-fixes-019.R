# Regression tests for the 0.1.19 P2/P3 whole-pipeline review fixes.
# Each block pins one item from docs/analysis/NTSminer_全流程函数评估_0.1.17.md (section 14).

ns <- asNamespace("NTSminer")

# ---- 9-2: mine_features reports an integer NCE_Count, not a relative Frequency -----

test_that("9-2 mine_features exposes NCE_Count (not the relative-intensity Frequency)", {
  skip_if_not_installed("data.table")
  sl <- data.table::data.table(
    Lab_ID = character(), Polarity = character(), NCE = character(),
    frag_mz = numeric(), intensity = numeric(), prec_mz = numeric()
  )
  meta <- data.table::data.table(Lab_ID = character(), Polarity = character())
  out <- mine_features(sl, meta)
  expect_true("NCE_Count" %in% names(out))
  expect_false("Frequency" %in% names(out))
})

# ---- 5-4: has_rt requires agreement within rt_tol_sec ----------------------------
# (0.1.20: has_rt was tightened from "lib_rt present" to "RT agrees within tol".)

test_that("5-4 .nts_confidence_has_rt_vec needs RT agreement within rt_tol_sec", {
  f <- get(".nts_confidence_has_rt_vec", envir = ns)
  # Numeric-looking character lib_rt still parses; agreement gates the flag.
  db <- data.frame(
    lib_rt     = c("120.5", "120.5", "abc", NA, ""),
    rt_error_s = c(5, 60, 5, 5, 5)
  )
  expect_equal(f(db), c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(f(data.frame(lib_rt = c(10, NA), rt_error_s = c(1, 1))),
               c(TRUE, FALSE))
  # A column without lib_rt is never "has RT".
  expect_equal(f(data.frame(other = 1:2)), c(FALSE, FALSE))
  # Without a measured RT error we cannot confirm agreement.
  expect_equal(f(data.frame(lib_rt = c(10, 20))), c(FALSE, FALSE))
})

# ---- 5-3: min_score default is a documented signature value of 0.5 ----------------

test_that("5-3 batch_annotate defaults min_score to 0.5", {
  expect_equal(formals(NTSminer::batch_annotate)$min_score, 0.5)
})

# ---- 6-1: homologue expansion follows edges in both directions by default ----------

test_that("6-1 build_homologue_network defaults to bidirectional expansion", {
  expect_true(isTRUE(formals(NTSminer::build_homologue_network)$bidirectional))
})

# ---- 4-1 / 4-2: extract_ms2 can reuse a pre-read experiment and reports drops ----

test_that("4-1 extract_ms2 gained an optional pre-read xdata argument", {
  expect_true("xdata" %in% names(formals(NTSminer::extract_ms2)))
})

test_that("4-2 extract_ms2 counts MS2 scans dropped for a missing precursor", {
  src <- paste(deparse(body(NTSminer::extract_ms2)), collapse = "\n")
  expect_true(grepl("n_dropped_no_prec", src, fixed = TRUE))
})

# ---- 3-2 / 3-3: scan diagnostics leave the RNG and mzR handles alone --------------

test_that("3-2 scan diagnostics do not perturb the global RNG", {
  src <- paste(deparse(body(NTSminer:::extract_spectra_from_scans)), collapse = "\n")
  expect_false(grepl("sample(", src, fixed = TRUE))
})

test_that("3-3 per-file mzR handles are closed without leaking on.exit handlers", {
  src <- paste(deparse(body(NTSminer:::extract_spectra_from_scans)), collapse = "\n")
  expect_false(grepl(
    "on.exit(try(mzR::close(ms), silent = TRUE), add = TRUE)", src, fixed = TRUE
  ))
})

# ---- 0-1 / 0-2 / 0-3: uniform object validation -----------------------------------

test_that("0-1 NTSProject has a validity method and bare new() is valid", {
  d <- tempfile("p2p3"); dir.create(d)
  expect_true(isTRUE(methods::validObject(create_nts_project(d))))
  expect_true(isTRUE(methods::validObject(methods::new("NTSProject"))))
})

test_that("0-2 NTSProject setters validate the result", {
  d <- tempfile("p2p3"); dir.create(d)
  p <- create_nts_project(d)
  expect_true(isTRUE(methods::validObject(set_project_log(p, list()))))
  b <- methods::new("NTSProject")
  b <- set_feature_table(b, data.frame(name = "F1", mz = 1, rt = 1, s = 1))
  expect_true(isTRUE(methods::validObject(b)))
})

test_that("0-3 get_ms2_score_grid returns an empty grid on a fresh project", {
  d <- tempfile("p2p3"); dir.create(d)
  p <- create_nts_project(d)
  expect_equal(nrow(NTSminer::get_ms2_score_grid(p)), 0L)
})

# ---- E-2: external executables resolve without machine-specific paths -------------

test_that("E-2 executables resolve via env/PATH/candidates and convert_mzml auto-detects", {
  f <- get(".nts_resolve_executable", envir = ns)
  expect_null(f("definitely_not_a_real_tool_xyz", candidates = character(0)))
  # env var takes precedence and is used verbatim
  old <- Sys.getenv("NTS_TEST_TOOL", unset = NA)
  Sys.setenv(NTS_TEST_TOOL = "sentinel-path")
  on.exit(
    if (is.na(old)) Sys.unsetenv("NTS_TEST_TOOL") else Sys.setenv(NTS_TEST_TOOL = old),
    add = TRUE
  )
  expect_equal(f("definitely_not_a_real_tool_xyz", envvar = "NTS_TEST_TOOL"), "sentinel-path")
  # an existing candidate path is returned when neither env nor PATH match
  cand <- tempfile(); file.create(cand)
  expect_equal(f("definitely_not_a_real_tool_xyz", candidates = cand), cand)
  expect_true(is.null(formals(NTSminer::convert_mzml)$exe))
})

# ---- E-4: the abandoned MSP/MGF parser island is gone -----------------------------

test_that("E-4 dead MSP/MGF/standards helpers were removed", {
  for (nm in c(
    ".nts_msp_trim", ".nts_msp_parse_peak_lines", ".nts_mgf_parse_block_lines",
    ".nts_mix_standards_cluster_ms2", ".nts_is_unknown_node",
    ".nts_is_anchor_node", ".nts_build_eic_traces"
  )) {
    expect_false(exists(nm, envir = ns, inherits = FALSE), info = nm)
  }
})

# ---- E-5: no R/*.bak build artefacts remain ---------------------------------------

test_that("E-5 no R/*.bak build artefacts remain in the source tree", {
  expect_length(list.files(file.path("..", "..", "R"), pattern = "\\.bak$"), 0L)
})

# ---- E-7: fingerprint bit-widths match the documented totals ----------------------

test_that("E-7 fingerprint base/full bit-widths match the roxygen documentation", {
  bits <- get(".FP_BITS", envir = ns)
  base <- get(".FP_BASE_TYPES", envir = ns)
  full <- get(".FP_FULL_TYPES", envir = ns)
  expect_equal(sum(bits[base]), 4119)
  expect_equal(sum(bits[full]), 10310)
})

# ---- E-8: convert_mzml folds its outcome into the project log ---------------------

test_that("E-8 convert_mzml records a log entry and returns the project when given one", {
  d <- tempfile("p2p3"); dir.create(d)
  p <- create_nts_project(d)
  out <- suppressWarnings(
    convert_mzml(d, tempfile("o"), exe = "echo", project = p)
  )
  expect_true(methods::is(out, "NTSProject"))
  lg <- NTSminer::get_project_log(out)
  last <- lg[[length(lg)]]
  expect_equal(last$step, "convert_mzml")
  expect_equal(last$params$n_files, 0L)
})

# ---- 7-4: print.nts_rt_model no longer prints a never-populated train slot --------

test_that("7-4 print.nts_rt_model renders without a train_metrics element", {
  m <- structure(list(
    model_version = "0.2.0", method = "rf",
    fp_types = character(0), descriptor_cols = character(0),
    p_features = 6L, p_original = 6L, n_training = 10L,
    test_metrics = list(R2 = 0.5, RMSE = 1, MAE = 1),
    training_date = Sys.time(), calib_shift = NULL, calib_scale = NULL
  ), class = "nts_rt_model")
  txt <- paste(capture.output(print(m)), collapse = "\n")
  expect_match(txt, "Test R2")
  expect_false(grepl("Train R2", txt, fixed = TRUE))
})
