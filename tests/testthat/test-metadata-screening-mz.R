test_that(".nts_metadata_screening_mz prefers target_mz over mz_neg", {
  meta <- data.frame(
    Lab_ID = c("A", "B"),
    target_mz = c(292.983, NA_real_),
    mz_neg = c(400, 500.1),
    mass = c(300, 400),
    stringsAsFactors = FALSE
  )
  mz <- NTSminer:::.nts_metadata_screening_mz(meta, polarity = "neg")
  expect_equal(mz[1L], 292.983, tolerance = 1e-6)
  expect_equal(mz[2L], 500.1, tolerance = 1e-6)
})

test_that(".nts_metadata_screening_mz uses mz_pos in positive mode", {
  meta <- data.frame(
    Lab_ID = "C",
    mz_pos = 501.2,
    mz_neg = 499.0,
    mass = 500,
    stringsAsFactors = FALSE
  )
  mz <- NTSminer:::.nts_metadata_screening_mz(meta, polarity = "pos")
  expect_equal(mz[1L], 501.2, tolerance = 1e-6)
})

test_that(".nts_db_ms1_table resolves target_mz before polarity/mass fallbacks", {
  skip_if_not_installed("data.table")
  info <- data.frame(
    Lab.ID = c("L1", "L2", "L3"),
    rt = c(1, 2, 3),
    target_mz = c(118.9926, NA_real_, NA_real_),
    mz.neg = c(200, 399, NA_real_),
    mz = c(299, 399, 499),
    stringsAsFactors = FALSE
  )
  db <- methods::new(
    "databaseClass",
    spectra.info = info,
    spectra.data = list(),
    fragments.data = list()
  )
  out <- NTSminer:::.nts_db_ms1_table(db, polarity = "neg")
  # L1: target_mz wins over mz.neg; L2: mz.neg; L3: mz - proton (last resort).
  expect_equal(out$mz[1L], 118.9926, tolerance = 1e-6)
  expect_equal(out$mz[2L], 399, tolerance = 1e-6)
  expect_equal(out$mz[3L], 499 - 1.007276, tolerance = 1e-6)
  expect_equal(out$exact_mass, c(299, 399, 499))
})
