# MS2 DDA match-score components: s_iso and s_purity (finding C4) -------------
# `s_iso` and `s_purity` were reserved-but-always-NA. They are now computed:
# s_iso from the scan's isolation window, s_purity from the file's own MS1 scan.

test_that(".nts_score_dda_match computes s_iso and consumes s_purity", {
  f <- NTSminer:::.nts_score_dda_match

  ms2_map <- data.table::data.table(
    feature_name = c("F1", "F2"),
    source_file = "a.mzML",
    precursorMz = c(100.0, 200.0),
    rt = c(10, 20),
    tic = c(1000, 1000),
    isolation_lo = c(99.5, 199.0),   # half-width 0.5 (F1), 1.0 (F2)
    isolation_hi = c(100.5, 201.0),
    purity = c(0.9, 0.3)
  )
  ms2_map[, mz := list(c(50, 60), c(70, 80))]
  ms2_map[, intensity := list(c(100, 50), c(100, 50))]

  ft <- data.table::data.table(
    feature_name = c("F1", "F2"),
    mz = c(100.0, 200.5),           # F2 sits 0.5 Da off the window centre
    rt = c(10, 20)
  )

  out <- f(ms2_map, ft, mz_tol_ppm = 15, rt_tol = 30)
  out <- out[order(out$feature_name)]

  # F1 precursor is dead-centre in its window -> s_iso == 1
  expect_equal(out$s_iso[out$feature_name == "F1"], 1)
  # F2 is 0.5 Da off a 1.0 Da half-width -> s_iso == 0.5
  expect_equal(out$s_iso[out$feature_name == "F2"], 0.5)
  # purity is carried through, clamped to [0, 1]
  expect_equal(out$s_purity[out$feature_name == "F1"], 0.9)
  expect_equal(out$s_purity[out$feature_name == "F2"], 0.3)
})

test_that(".nts_score_dda_match leaves components NA when inputs are absent", {
  f <- NTSminer:::.nts_score_dda_match

  ms2_map <- data.table::data.table(
    feature_name = "F1",
    source_file = "a.mzML",
    precursorMz = 100.0,
    rt = 10,
    tic = 1000
  )
  ms2_map[, mz := list(c(50, 60))]
  ms2_map[, intensity := list(c(100, 50))]

  ft <- data.table::data.table(feature_name = "F1", mz = 100.0, rt = 10)

  out <- f(ms2_map, ft, mz_tol_ppm = 15, rt_tol = 30)
  expect_true(is.na(out$s_iso))
  expect_true(is.na(out$s_purity))
})

test_that(".nts_iso_purity_for_file returns NA without accessible MS1", {
  g <- NTSminer:::.nts_iso_purity_for_file
  out <- suppressWarnings(
    g("/no/such/file.mzML",
      scan_rt = c(1, 2), iso_lo = c(0.9, 1.9), iso_hi = c(1.1, 2.1),
      prec_mz = c(1.0, 2.0), mz_tol_ppm = 10)
  )
  expect_equal(out, c(NA_real_, NA_real_))
})
