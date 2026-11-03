# Default-on diagnostic-fragment validation (packaged PFAS expert table) -----

test_that(".nts_apply_default_diagnostic_fragments tags matched PFAS fragments", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_diag_default")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  p <- create_nts_project(td)
  p <- set_feature_table(
    p,
    data.frame(name = c("F1", "F2"), mz = c(413, 500), rt = c(100, 200),
               stringsAsFactors = FALSE)
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("F1", "F2"),
      mz = c(413, 500),
      rt = c(100, 200),
      annotation_level = c("Unknown candidate", "Unknown candidate"),
      stringsAsFactors = FALSE
    )
  )
  # F1 carries CF2- (49.9973564) and SO3- (79.957367); F2 carries neither.
  p@ms2data <- list(
    F1 = list(list(mz = c(49.9974, 79.9574, 118.9), intensity = c(1000, 800, 50))),
    F2 = list(list(mz = c(200.1, 300.2), intensity = c(500, 400)))
  )

  p <- .nts_apply_default_diagnostic_fragments(p, mz_tol_ppm = 10)

  m2 <- get_ms2_annotations(p)
  expect_true("diagnostic_fragments" %in% names(m2))
  expect_match(m2$diagnostic_fragments[m2$name == "F1"], "CF2-")
  expect_match(m2$diagnostic_fragments[m2$name == "F1"], "SO3-")
  expect_true(is.na(m2$diagnostic_fragments[m2$name == "F2"]))

  m1 <- get_ms1_annotations(p)
  expect_match(m1$annotation_level[m1$name == "F1"], "MS2 Verified", fixed = TRUE)
  expect_false(grepl("MS2 Verified", m1$annotation_level[m1$name == "F2"], fixed = TRUE))
})


test_that(".nts_apply_default_diagnostic_fragments no-ops without MS2 data", {
  td <- tempfile("nts_diag_noms2")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  p <- create_nts_project(td)
  p <- set_feature_table(p, data.frame(name = "F1", mz = 300, rt = 1,
                                       stringsAsFactors = FALSE))

  out <- .nts_apply_default_diagnostic_fragments(p)
  # No ms2data -> unchanged, no error
  expect_s4_class(out, "NTSProject")
})


test_that("annotate_by_spectral_match runs diagnostic validation by default (toggleable)", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_diag_ann")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)

  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = "L1", target_mz = 300, RT = 1, name = "TestPFAS", SMILES = "CC"
    ),
    spectra_long = data.table::data.table(
      Lab_ID = "L1", Polarity = "NEG", NCE = "CE10",
      frag_mz = c(100, 150), intensity = c(50, 30)
    ),
    features_long = data.table::data.table(),
    version = "t"
  )
  p <- create_nts_project(td)
  p <- set_feature_table(p, data.frame(name = "F1", mz = 300, rt = 1,
                                       stringsAsFactors = FALSE))
  p@ms2data <- list(
    F1 = list(list(mz = c(100, 150, 49.9974), intensity = c(50, 30, 20)))
  )

  p_on <- annotate_by_spectral_match(
    p, db, mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.2, polarity = "neg"
  )
  expect_true("diagnostic_fragments" %in% names(get_ms2_annotations(p_on)))

  p_off <- annotate_by_spectral_match(
    p, db, mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.2, polarity = "neg",
    validate_diagnostic_fragments = FALSE
  )
  expect_false("diagnostic_fragments" %in% names(get_ms2_annotations(p_off)))
})
