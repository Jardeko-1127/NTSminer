## Unit tests for infer_isf_parents: ISF fragment-to-database candidate lookup.

test_that("infer_isf_parents fills candidate columns on ISF fragment rows", {
  skip_if_not_installed("data.table")

  # Build mock NTSDatabase
  meta <- data.table::data.table(
    Lab_ID = c("C001", "C002", "C003"),
    mass = c(500.0, 500.1, 600.0),
    Polarity = "NEG",
    SMILES = c("C1", "C2", "C3"),
    Name = c("Compound1", "Compound2", "Compound3"),
    stringsAsFactors = FALSE
  )
  feat <- data.table::data.table(
    Lab_ID = c("C001", "C001", "C002", "C003"),
    Polarity = "NEG",
    Feature_Type = "DF",
    mz_value = c(218.985, 280.995, 300.0, 400.0),
    Frequency = c(0.9, 0.5, 0.8, 0.7)
  )
  db <- methods::new("NTSDatabase",
    metadata = meta,
    spectra_long = data.table::data.table(),
    features_long = feat,
    version = "1.0.0"
  )

  # Build input table with ISF fragments pointing to a parent at m/z 500
  dt <- data.table::data.table(
    name = c("FT00001", "FT00002", "FT00003"),
    mz = c(500.001, 218.985, 280.995),
    rt = c(300, 301, 302),
    is_isf = c(FALSE, TRUE, TRUE),
    parent_mz = c(NA_real_, 500.001, 500.001),
    isf_annotation = c(NA_character_, "ISF of FT00001 (-HF)", "ISF of FT00001 (-CO2)")
  )

  out <- infer_isf_parents(dt, db, mz_tol_ppm = 20)

  # ISF fragment rows should have candidates
  expect_false(is.na(out$isf_candidate_ids[2]))
  expect_true(grepl("C001", out$isf_candidate_ids[2]))
  expect_true(grepl("Compound1", out$isf_candidate_names[2]))

  # Parent row (FT00001) should also get candidates via m/z matching
  expect_false(is.na(out$isf_candidate_ids[1]))
})

test_that("infer_isf_parents returns NA candidates when no ISF features exist", {
  skip_if_not_installed("data.table")

  meta <- data.table::data.table(
    Lab_ID = "C001", mass = 500.0, Polarity = "NEG",
    SMILES = "C1", Name = "Compound1"
  )
  feat <- data.table::data.table(
    Lab_ID = "C001", Polarity = "NEG", Feature_Type = "DF",
    mz_value = 200.0, Frequency = 0.9
  )
  db <- methods::new("NTSDatabase",
    metadata = meta, spectra_long = data.table::data.table(),
    features_long = feat, version = "1.0.0"
  )

  dt <- data.table::data.table(
    name = "FT00001", mz = 500.0, rt = 300,
    is_isf = FALSE, parent_mz = NA_real_
  )

  out <- infer_isf_parents(dt, db)
  expect_true(all(is.na(out$isf_candidate_ids)))
  expect_true(all(is.na(out$isf_candidate_names)))
})

test_that("infer_isf_parents respects polarity filter", {
  skip_if_not_installed("data.table")

  meta <- data.table::data.table(
    Lab_ID = c("C001", "C002"),
    mass = c(500.0, 500.1),
    Polarity = c("NEG", "POS"),
    SMILES = c("C1", "C2"),
    Name = c("Compound1", "Compound2")
  )
  feat <- data.table::data.table(
    Lab_ID = c("C001", "C002"),
    Polarity = c("NEG", "POS"),
    Feature_Type = "DF",
    mz_value = c(218.985, 218.985),
    Frequency = c(0.9, 0.9)
  )
  db <- methods::new("NTSDatabase",
    metadata = meta, spectra_long = data.table::data.table(),
    features_long = feat, version = "1.0.0"
  )

  dt <- data.table::data.table(
    name = c("FT00001", "FT00002"),
    mz = c(500.0, 218.985),
    rt = c(300, 301),
    is_isf = c(FALSE, TRUE),
    parent_mz = c(NA_real_, 500.0)
  )

  out_neg <- infer_isf_parents(dt, db, polarity = "NEG")
  expect_true(grepl("C001", out_neg$isf_candidate_ids[2]))
  expect_false(grepl("C002", out_neg$isf_candidate_ids[2]))

  out_pos <- infer_isf_parents(dt, db, polarity = "POS")
  expect_false(grepl("C001", out_pos$isf_candidate_ids[2]))
  expect_true(grepl("C002", out_pos$isf_candidate_ids[2]))
})

test_that("infer_isf_parents ranks candidates by matching fragment count", {
  skip_if_not_installed("data.table")

  meta <- data.table::data.table(
    Lab_ID = c("C001", "C002"),
    mass = c(500.0, 500.1),
    Polarity = "NEG",
    SMILES = c("C1", "C2"),
    Name = c("Compound1", "Compound2")
  )
  # C001 matches BOTH ISF fragments; C002 matches only one
  feat <- data.table::data.table(
    Lab_ID = c("C001", "C001", "C002"),
    Polarity = "NEG",
    Feature_Type = "DF",
    mz_value = c(218.985, 280.995, 218.985),
    Frequency = c(0.9, 0.5, 0.8)
  )
  db <- methods::new("NTSDatabase",
    metadata = meta, spectra_long = data.table::data.table(),
    features_long = feat, version = "1.0.0"
  )

  dt <- data.table::data.table(
    name = c("FT00001", "FT00002", "FT00003"),
    mz = c(500.0, 218.985, 280.995),
    rt = c(300, 301, 302),
    is_isf = c(FALSE, TRUE, TRUE),
    parent_mz = c(NA_real_, 500.0, 500.0)
  )

  out <- infer_isf_parents(dt, db)

  # C001 should come first (2 fragments matched) before C002 (1 fragment)
  ids <- strsplit(out$isf_candidate_ids[2], "; ")[[1]]
  expect_equal(ids[1], "C001")
  expect_equal(ids[2], "C002")

  counts <- strsplit(out$isf_candidate_n_fragments[2], "; ")[[1]]
  expect_equal(as.integer(counts), c(2L, 1L))
})

test_that("infer_isf_parents writes candidates to parent feature via mz match", {
  skip_if_not_installed("data.table")

  meta <- data.table::data.table(
    Lab_ID = "C001", mass = 500.0, Polarity = "NEG",
    SMILES = "C1", Name = "Compound1"
  )
  feat <- data.table::data.table(
    Lab_ID = "C001", Polarity = "NEG", Feature_Type = "DF",
    mz_value = 218.985, Frequency = 0.9
  )
  db <- methods::new("NTSDatabase",
    metadata = meta, spectra_long = data.table::data.table(),
    features_long = feat, version = "1.0.0"
  )

  dt <- data.table::data.table(
    name = c("FT00001", "FT00002"),
    mz = c(500.0, 218.985),
    rt = c(300, 301),
    is_isf = c(FALSE, TRUE),
    parent_mz = c(NA_real_, 500.0)
  )

  out <- infer_isf_parents(dt, db)

  # Parent row (non-ISF at m/z matching parent_mz) gets candidates
  expect_false(is.na(out$isf_candidate_ids[1]))
  expect_true(grepl("C001", out$isf_candidate_ids[1]))
  expect_true(grepl("Compound1", out$isf_candidate_names[1]))
})

test_that("infer_isf_parents handles empty ISF gracefully", {
  skip_if_not_installed("data.table")

  meta <- data.table::data.table(
    Lab_ID = "C001", mass = 500.0, Polarity = "NEG",
    SMILES = "C1", Name = "Compound1"
  )
  feat <- data.table::data.table(
    Lab_ID = "C001", Polarity = "NEG", Feature_Type = "DF",
    mz_value = 200.0, Frequency = 0.9
  )
  db <- methods::new("NTSDatabase",
    metadata = meta, spectra_long = data.table::data.table(),
    features_long = feat, version = "1.0.0"
  )

  dt <- data.table::data.table(
    name = c("FT00001", "FT00002"),
    mz = c(500.0, 218.985),
    is_isf = c(FALSE, FALSE),  # no ISF
    parent_mz = c(NA_real_, NA_real_)
  )

  out <- infer_isf_parents(dt, db)
  expect_true(all(is.na(out$isf_candidate_ids)))
  expect_true("isf_candidate_ids" %in% names(out))
})
