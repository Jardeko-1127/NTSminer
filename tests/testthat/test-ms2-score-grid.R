test_that("ms2_score_grid starts empty and follows the per-cell contract", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_grid_ce")
  dir.create(td)
  # One compound, two library CEs; the two experimental scans each reproduce
  # exactly one of the CE spectra, so the grid is a full CE x spectrum cross.
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
  p@ms2data <- list(
    F1 = list(
      list(mz = c(100, 150), intensity = c(100, 50), ce = 10, scan_idx = 1L),
      list(mz = c(200, 250), intensity = c(100, 50), ce = 20, scan_idx = 2L)
    )
  )

  # Empty before any annotation.
  expect_true(is.data.frame(get_ms2_score_grid(p)))
  expect_equal(nrow(get_ms2_score_grid(p)), 0L)

  p2 <- annotate_by_spectral_match(
    p, db,
    mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.9, polarity = "neg"
  )
  g <- get_ms2_score_grid(p2)

  # One row per (feature x ref_id x library CE x experimental spectrum).
  expect_equal(nrow(g), 4L)
  expect_setequal(unique(g$name), "F1")
  expect_setequal(unique(g$ref_id), "L1")
  expect_setequal(unique(g$lib_CE), c("CE10", "CE20"))
  expect_equal(nrow(unique(g[, c("lib_CE", "exp_spectrum")])), 4L)

  # Columns promised by the docs are all present.
  expect_true(all(c(
    "name", "ref_id", "lib_CE", "exp_spectrum", "score_cosine",
    "n_matched_peaks", "passes_min_score", "is_winner", "score_ms2",
    "score_ndp", "ndp_raw", "score_nl", "score_md", "score_diagnostic"
  ) %in% names(g)))

  # Cells below the gate are kept, not dropped.
  expect_true(is.logical(g$passes_min_score))
  expect_equal(sum(g$passes_min_score), 2L)
  expect_equal(sum(!g$passes_min_score), 2L)
  # A non-passing cell still carries its (raw) cosine but no composite.
  failed <- g[!g$passes_min_score, , drop = FALSE]
  expect_true(all(failed$score_cosine < 0.9))
  expect_true(all(is.na(failed$score_ms2)))
  # Passing cells carry a finite composite.
  expect_true(all(is.finite(g$score_ms2[g$passes_min_score])))

  # Exactly one winner per scored (feature, compound) pair.
  expect_equal(sum(g$is_winner), 1L)
  expect_equal(nrow(p2@db_annotations), 1L)

  # The winner row mirrors the merged db_annotations result exactly.
  w <- g[g$is_winner, , drop = FALSE]
  expect_equal(as.character(w$lib_CE), as.character(p2@db_annotations$lib_CE))
  expect_equal(as.character(w$exp_spectrum), as.character(p2@db_annotations$exp_CE))
  expect_equal(w$score_cosine, p2@db_annotations$score_cosine_ms2)
})


test_that("is_winner has one row per db_annotations row across features/compounds", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_grid_multi")
  dir.create(td)
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
    data.frame(name = c("F1", "F2"), mz = c(300, 300), rt = c(1, 1))
  )
  p@ms2data <- list(
    F1 = list(list(mz = c(100, 150), intensity = c(100, 50))),
    F2 = list(list(mz = c(100, 150), intensity = c(100, 50)))
  )

  p2 <- annotate_by_spectral_match(
    p, db,
    mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.2, polarity = "neg"
  )
  g <- get_ms2_score_grid(p2)
  # 2 features x 2 compounds x 1 CE x 1 spectrum = 4 cells.
  expect_equal(nrow(g), 4L)
  # Each feature x compound contributes exactly one winner.
  expect_equal(sum(g$is_winner), nrow(p2@db_annotations))
  wins <- g[g$is_winner, , drop = FALSE]
  expect_equal(
    nrow(unique(wins[, c("name", "ref_id")])),
    sum(g$is_winner)
  )
  # Winners align cell-for-cell with db_annotations.
  key_g <- paste(wins$name, wins$ref_id, sep = "|")
  key_a <- paste(p2@db_annotations$name, p2@db_annotations$ref_id, sep = "|")
  expect_setequal(key_g, key_a)
})


test_that("return_score_grid = FALSE leaves the grid empty but still annotates", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_grid_off")
  dir.create(td)
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = "L1", target_mz = 300, RT = 1, name = "TestPFAS", SMILES = "CC"
    ),
    spectra_long = data.table::data.table(
      Lab_ID = c("L1", "L1"), Polarity = "NEG", NCE = c("CE10", "CE10"),
      frag_mz = c(100, 150), intensity = c(100, 50)
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
    mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.5, polarity = "neg",
    return_score_grid = FALSE
  )
  expect_equal(nrow(get_ms2_score_grid(p2)), 0L)
  expect_equal(nrow(p2@db_annotations), 1L)
})


test_that("ms2_score_grid survives a copy_project deep copy", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_grid_copy")
  dir.create(td)
  p <- create_nts_project(td)
  p <- set_ms2_score_grid(
    p,
    data.frame(
      name = "F1", ref_id = "L1", lib_CE = "CE10", exp_spectrum = "CE10",
      score_cosine = 0.9, passes_min_score = TRUE, is_winner = TRUE,
      stringsAsFactors = FALSE
    )
  )
  pc <- copy_project(p)
  expect_equal(nrow(get_ms2_score_grid(pc)), 1L)
  # Mutating the copy must not touch the original.
  pc <- set_ms2_score_grid(pc, data.frame())
  expect_equal(nrow(get_ms2_score_grid(p)), 1L)
  expect_equal(nrow(get_ms2_score_grid(pc)), 0L)
})


test_that("is_winner stays unique when two experimental spectra share a label", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_grid_duplab")
  dir.create(td)
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = "L1", target_mz = 300, RT = 1, name = "TestPFAS", SMILES = "CC"
    ),
    spectra_long = data.table::data.table(
      Lab_ID = c("L1", "L1"), Polarity = "NEG", NCE = c("CE10", "CE10"),
      frag_mz = c(100, 150), intensity = c(100, 50)
    ),
    features_long = data.table::data.table(),
    version = "t"
  )
  p <- create_nts_project(td)
  p <- set_feature_table(
    p,
    data.frame(name = "F1", mz = 300, rt = 1, stringsAsFactors = FALSE)
  )
  # `scan_idx` is only unique within one file: two scans from two files share it.
  p@ms2data <- list(
    F1 = list(
      list(mz = c(100, 150), intensity = c(100, 50), ce = 10, scan_idx = 5L,
           source_file = "/data/runA.mzML"),
      list(mz = c(100, 150), intensity = c(100, 50), ce = 10, scan_idx = 5L,
           source_file = "/data/runB.mzML")
    )
  )

  p2 <- annotate_by_spectral_match(
    p, db,
    mz_tol_ppm = 50, ms2_tol_da = 0.05, min_score = 0.5, polarity = "neg"
  )
  g <- get_ms2_score_grid(p2)
  expect_equal(nrow(g), 2L)             # 1 ref x 1 CE x 2 spectra
  expect_equal(sum(g$is_winner), 1L)    # exactly one winner, not one per label
  # Labels are disambiguated by the source file so cells stay distinguishable.
  expect_equal(length(unique(g$exp_spectrum)), 2L)
  expect_true(all(grepl("@run[AB]$", g$exp_spectrum)))
  # The winner's label is the one recorded in db_annotations.
  expect_equal(
    as.character(g$exp_spectrum[g$is_winner]),
    as.character(p2@db_annotations$exp_CE)
  )
})

test_that("get_ms2_score_grid folds a pre-refactor project into the list", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_grid_legacy")
  dir.create(td)
  p <- create_nts_project(td)
  # A fresh project reads as an empty grid, not an error.
  expect_true(is.data.frame(get_ms2_score_grid(p)))
  expect_equal(nrow(get_ms2_score_grid(p)), 0L)
  # Pre-refactor objects kept the grid in a standalone `@ms2_score_grid` slot
  # and `@ms2_annotations` as a flat data.frame. Simulate that stored shape.
  methods::slot(p, "ms2_annotations") <-
    data.frame(name = "F1", n_fragments = 3L, stringsAsFactors = FALSE)
  attr(p, "ms2_score_grid") <- data.frame(
    name = "F1", ref_id = "L1", lib_CE = "CE10", exp_spectrum = "CE10",
    score_cosine = 0.9, passes_min_score = TRUE, is_winner = TRUE,
    stringsAsFactors = FALSE
  )
  # The getter folds both legacy shapes into the list.
  expect_equal(nrow(get_ms2_score_grid(p)), 1L)
  expect_equal(nrow(get_ms2_annotations(p)), 1L)
  expect_equal(nrow(get_ms2_annotations(p, "score_grid")), 1L)
  # The setter migrates the object onto the new list layout (legacy attr gone).
  p2 <- set_ms2_score_grid(p, get_ms2_score_grid(p))
  expect_equal(nrow(get_ms2_score_grid(p2)), 1L)
  expect_null(attr(p2, "ms2_score_grid"))
})

test_that("batch_annotate keeps grid winners in step with the deduped db_annotations", {
  skip_if_not_installed("data.table")

  td <- tempfile("nts_grid_batch"); dir.create(td)
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = "L1", target_mz = 300, RT = 1, name = "TestPFAS", SMILES = "CC"
    ),
    spectra_long = data.table::data.table(
      Lab_ID = c("L1", "L1"), Polarity = "NEG", NCE = c("CE10", "CE10"),
      frag_mz = c(100, 150), intensity = c(100, 50)
    ),
    features_long = data.table::data.table(),
    version = "t"
  )
  p <- create_nts_project(td)
  p <- set_feature_table(
    p,
    data.frame(name = c("F1", "F2"), mz = c(300, 300), rt = c(1, 1))
  )
  # Two features compete for the SAME library compound: annotate writes 2 rows
  # (and 2 grid winners), then batch_annotate dedups db_annotations to 1.
  p@ms2data <- list(
    F1 = list(list(mz = c(100, 150), intensity = c(100, 50))),
    F2 = list(list(mz = c(100, 150), intensity = c(100, 50)))
  )

  p2 <- batch_annotate(
    p, db,
    mz_tol = 50, feature_tol = 0.05, min_score = 0.2, polarity = "neg"
  )
  ann <- get_db_annotations(p2)
  g <- get_ms2_score_grid(p2)

  expect_equal(nrow(ann), 1L)                 # deduped to one row per ref_id
  expect_equal(nrow(g), 2L)                   # both cells kept in the grid
  expect_equal(sum(g$is_winner), nrow(ann))   # exactly as many winners as rows
  expect_equal(sum(g$is_winner), 1L)
  # The single winner is the surviving feature"s cell.
  expect_equal(as.character(g$name[g$is_winner]), as.character(ann$name))
  expect_equal(as.character(g$ref_id[g$is_winner]), as.character(ann$ref_id))
})
