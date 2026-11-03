test_that("annotate_by_spectral_match keeps one hit per feature per Lab_ID", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_ann_multi")
  dir.create(td)
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = c("L1", "L2"),
      target_mz = c(300, 300),
      RT = c(1, 1),
      name = c("CmpA", "CmpB"),
      SMILES = c("CC", "CN")
    ),
    spectra_long = data.table::data.table(
      Lab_ID = rep(c("L1", "L2"), each = 2L),
      Polarity = "NEG",
      NCE = "CE10",
      frag_mz = c(100, 150, 100, 180),
      intensity = c(50, 30, 50, 40)
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
      list(mz = c(100, 150), intensity = c(50, 30)),
      list(mz = c(100, 180), intensity = c(50, 40))
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
  ann <- p2@db_annotations
  expect_equal(nrow(ann), 2L)
  expect_equal(sort(as.character(ann$ref_id)), c("L1", "L2"))
  expect_true(all(ann$name == "F1"))
})
