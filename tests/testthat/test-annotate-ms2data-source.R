test_that("annotate_by_spectral_match uses ms2data when clean_ms2 is present", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_ann_ms2")
  dir.create(td)
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = "L1",
      target_mz = 300,
      RT = 1,
      name = "TestPFAS",
      SMILES = "CC"
    ),
    spectra_long = data.table::data.table(
      Lab_ID = "L1",
      Polarity = "NEG",
      NCE = "CE10",
      frag_mz = c(100, 150),
      intensity = c(50, 30)
    ),
    features_long = data.table::data.table(),
    version = "t"
  )
  p <- create_nts_project(td)
  p <- set_feature_table(
    p,
    data.frame(name = "F1", mz = 300, rt = 1, stringsAsFactors = FALSE)
  )
  p@ms2data <- list(
    F1 = list(
      list(mz = c(100, 150), intensity = c(50, 30))
    )
  )
  p@clean_ms2 <- list(
    list(
      feature_id = "F1",
      spectrum = data.frame(mz = 999, intensity = 1)
    )
  )

  p2 <- annotate_by_spectral_match(
    p,
    db,
    mz_tol_ppm = 50,
    ms2_tol_da = 0.05,
    min_score = 0.2,
    polarity = "neg"
  )
  expect_true(nrow(p2@db_annotations) >= 1L)
  expect_equal(as.character(p2@db_annotations$ref_id[1]), "L1")
})
