test_that("per-CE reference rescues a spectrum lost by a lossy Merged layer", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_ann_perce")
  dir.create(td)
  # Merged layer deliberately drops the CE10 fragments (simulating a lossy
  # consensus); the per-CE strategy must still match the CE10 experimental scan.
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
      Lab_ID = c("L1", "L1", "L1", "L1", "L1", "L1"),
      Polarity = "NEG",
      NCE = c("CE10", "CE10", "CE20", "CE20", "Merged", "Merged"),
      frag_mz = c(100, 150, 200, 250, 200, 250),
      intensity = c(100, 50, 100, 50, 100, 50)
    ),
    features_long = data.table::data.table(),
    version = "t"
  )
  p <- create_nts_project(td)
  p <- set_feature_table(
    p,
    data.frame(name = "F1", mz = 300, rt = 1, stringsAsFactors = FALSE)
  )
  # A single experimental scan carrying the CE10 fragments.
  p@ms2data <- list(
    F1 = list(list(mz = c(100, 150), intensity = c(100, 50)))
  )

  # per_ce (default): matches CE10 and reports which CE won.
  p_ce <- annotate_by_spectral_match(
    p, db,
    mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.5, polarity = "neg"
  )
  expect_equal(nrow(p_ce@db_annotations), 1L)
  expect_equal(as.character(p_ce@db_annotations$ref_id[1]), "L1")
  expect_equal(as.character(p_ce@db_annotations$lib_CE[1]), "CE10")
  expect_gt(p_ce@db_annotations$score_cosine_ms2[1], 0.99)

  # Legacy merged strategy only sees the two CE20 fragments -> no match.
  p_mg <- annotate_by_spectral_match(
    p, db,
    mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.5, polarity = "neg",
    ref_spectra = "merged"
  )
  expect_equal(nrow(p_mg@db_annotations), 0L)
})


test_that("per-CE grid scores each experimental scan against each library CE", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_ann_grid")
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
      Lab_ID = c("L1", "L1", "L1", "L1"),
      Polarity = "NEG",
      NCE = c("CE10", "CE10", "CE20", "CE20"),
      frag_mz = c(100, 150, 200, 250),
      intensity = c(100, 50, 100, 50)
    ),
    features_long = data.table::data.table(),
    version = "t"
  )
  p <- create_nts_project(td)
  p <- set_feature_table(
    p,
    data.frame(name = "F1", mz = 300, rt = 1, stringsAsFactors = FALSE)
  )
  # Two scans: the first matches CE10, the second matches CE20.
  p@ms2data <- list(
    F1 = list(
      list(mz = c(100, 150), intensity = c(100, 50)),
      list(mz = c(200, 250), intensity = c(100, 50))
    )
  )

  p2 <- annotate_by_spectral_match(
    p, db,
    mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.5, polarity = "neg"
  )
  expect_equal(nrow(p2@db_annotations), 1L)
  # Both scans score 1.0 against their own CE; the grid must surface a CE hit.
  expect_true(as.character(p2@db_annotations$lib_CE[1]) %in% c("CE10", "CE20"))
  expect_gt(p2@db_annotations$score_cosine_ms2[1], 0.99)
})


test_that("db_annotations keeps every passing compound but nodes stay one per feature", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_ann_multi")
  dir.create(td)
  # Two compounds share the precursor m/z. L1 reproduces the experimental
  # spectrum exactly; L2 shares the same two fragment m/z but with a very
  # different intensity ratio, so it is a genuine (weaker) candidate.
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = c("L1", "L2"),
      target_mz = c(300, 300),
      RT = c(1, 1),
      name = c("Win", "Runner"),
      SMILES = c("CC", "CCC")
    ),
    spectra_long = data.table::data.table(
      Lab_ID = c("L1", "L1", "L2", "L2"),
      Polarity = "NEG",
      NCE = c("CE10", "CE10", "CE10", "CE10"),
      frag_mz = c(100, 150, 100, 150),
      intensity = c(100, 50, 100, 1000)
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
    F1 = list(list(mz = c(100, 150), intensity = c(100, 50)))
  )

  p2 <- annotate_by_spectral_match(
    p, db,
    mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.2, polarity = "neg"
  )
  # Both compounds pass the (low) cosine gate and are retained.
  expect_equal(nrow(p2@db_annotations), 2L)
  expect_setequal(as.character(p2@db_annotations$ref_id), c("L1", "L2"))
  # The per-CE grid discriminates on fragment intensity: the exact match wins.
  ann <- p2@db_annotations
  cos_l1 <- ann$score_cosine_ms2[ann$ref_id == "L1"]
  cos_l2 <- ann$score_cosine_ms2[ann$ref_id == "L2"]
  expect_gt(cos_l1, cos_l2)
  # Rows are ordered by score_total descending within a feature.
  expect_true(all(diff(ann$score_total) <= 1e-9))
  # Annotating no longer stashes a legacy `@network_results` stub (0.1.28):
  # scoring writes `db_annotations` / `ms2_annotations` only, and the network
  # slot stays empty so a mined project carries no network content.
  expect_length(p2@network_results, 0L)
  # The top-scoring compound for the feature is the first `db_annotations` row.
  expect_equal(as.character(ann$Compound.name[1]), "Win")
})


test_that("batch_annotate inherits the per_ce default (no silent merged fallback)", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_ann_batch")
  dir.create(td)
  # Same lossy Merged layer as the first test: only the per-CE strategy can
  # match the CE10 experimental scan. If batch_annotate fell back to "merged"
  # (its old hard-coded default) this would return zero annotations.
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
      Lab_ID = rep("L1", 6L),
      Polarity = "NEG",
      NCE = c("CE10", "CE10", "CE20", "CE20", "Merged", "Merged"),
      frag_mz = c(100, 150, 200, 250, 200, 250),
      intensity = c(100, 50, 100, 50, 100, 50)
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
    F1 = list(list(mz = c(100, 150), intensity = c(100, 50)))
  )

  p2 <- batch_annotate(
    p, db,
    mz_tol = 50, rt_tol = 150, feature_tol = 0.05,
    min_score = 0.5, polarity = "neg"
  )
  expect_equal(nrow(p2@db_annotations), 1L)
  expect_equal(as.character(p2@db_annotations$ref_id[1]), "L1")
  expect_equal(as.character(p2@db_annotations$lib_CE[1]), "CE10")
})
