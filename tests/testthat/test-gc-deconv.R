test_that("@gc_deconvolved slot exists and is accessible", {
  p <- create_nts_project(tempdir())
  expect_true("gc_deconvolved" %in% methods::slotNames(p))
  expect_type(get_gc_deconvolved(p), "list")
  expect_length(get_gc_deconvolved(p), 0L)
})

test_that("set_gc_deconvolved round-trips", {
  p <- create_nts_project(tempdir())
  dummy <- list(
    "FT001" = list(
      list(spectrum = data.frame(mz = 1:3, intensity = 10:12),
           apex_rt = 10.5, n_ions = 3L)
    )
  )
  p <- set_gc_deconvolved(p, dummy)
  expect_equal(names(get_gc_deconvolved(p)), "FT001")
  expect_equal(get_gc_deconvolved(p)[["FT001"]][[1]]$apex_rt, 10.5)
  expect_equal(nrow(get_gc_deconvolved(p)[["FT001"]][[1]]$spectrum), 3L)
})

test_that("deconvolve_gc_spectra handles empty gc_spectra", {
  p <- create_nts_project(tempdir())
  expect_error(deconvolve_gc_spectra(p), "empty|extract_gc_spectra")
})

test_that(".nts_eic_noise_sd estimates noise from quiet regions", {
  # Signal with known noise floor
  eic <- c(2, 1, 3, 2, 100, 3, 1, 2, 1, 0, 2, 1, 3)
  noise_sd <- NTSminer:::.nts_eic_noise_sd(eic)
  expect_true(noise_sd > 0)
  expect_true(noise_sd < 5)  # noise region has low SD
})

test_that(".nts_eic_noise_sd handles all-zero input", {
  nsd <- NTSminer:::.nts_eic_noise_sd(numeric(10))
  expect_true(is.finite(nsd))
})

test_that(".nts_find_eic_apexes finds correct local maxima", {
  eic <- c(0, 0, 5, 20, 100, 50, 10, 0, 3, 30, 60, 20, 0)
  apexes <- NTSminer:::.nts_find_eic_apexes(eic, min_snr = 2)
  expect_true(5L %in% apexes)   # value 100
  expect_true(11L %in% apexes)  # value 60
})

test_that(".nts_find_eic_apexes filters by min_snr", {
  eic <- c(0, 0, 2, 5, 2, 0, 0, 0, 3, 8, 3, 0)
  apexes_low <- NTSminer:::.nts_find_eic_apexes(eic, min_snr = 1)
  apexes_high <- NTSminer:::.nts_find_eic_apexes(eic, min_snr = 100)
  expect_true(length(apexes_low) >= 1L)
  expect_equal(length(apexes_high), 0L)
})

test_that(".nts_find_eic_apexes merges adjacent apexes", {
  # Two close peaks: at indices 4 and 5, should merge to keep the higher (index 4)
  eic <- c(0, 0, 3, 100, 80, 2, 0, 0)
  apexes <- NTSminer:::.nts_find_eic_apexes(eic, min_snr = 2)
  expect_equal(length(apexes), 1L)
  expect_equal(apexes[1L], 4L)  # higher of the two adjacent
})

test_that(".nts_cluster_apexes groups by RT tolerance", {
  rt <- c(10.0, 10.8, 11.2, 14.0, 14.5)
  clusters <- NTSminer:::.nts_cluster_apexes(rt, apex_tol = 1.5)
  expect_equal(clusters, c(1L, 1L, 1L, 2L, 2L))
})

test_that(".nts_cluster_apexes handles single element", {
  expect_equal(NTSminer:::.nts_cluster_apexes(15, 1.5), 1L)
})

test_that(".nts_cluster_apexes handles empty input", {
  expect_equal(length(NTSminer:::.nts_cluster_apexes(numeric(0), 1.5)), 0L)
})

test_that("score_ei_spectra use_deconvolved switches source correctly", {
  skip_if_not_installed("OrgMassSpecR")

  p <- create_nts_project(tempdir())
  p <- set_gc_spectra(p, list(
    "FT001" = data.frame(mz = 50:55, intensity = 10:15)
  ))
  p <- set_gc_deconvolved(p, list(
    "FT001" = list(
      list(spectrum = data.frame(mz = 77:80, intensity = 20:23),
           apex_rt = 10.5, n_ions = 4L)
    )
  ))

  lib <- methods::new("NTSDatabase",
    metadata = data.table::data.table(lab_id = "R1", name = "Test"),
    spectra_long = data.table::data.table(
      lab_id = "R1", mz = 77:80, intensity = 20:23
    ),
    features_long = data.table::data.table(),
    version = "0.1.0"
  )

  r1 <- score_ei_spectra(p, lib, min_score = 0, use_deconvolved = FALSE, prefilter_threshold = 0)
  r2 <- score_ei_spectra(p, lib, min_score = 0, use_deconvolved = TRUE, prefilter_threshold = 0)

  expect_true(nrow(r2) > 0)
  # Deconvolved spectrum (77-80) should match reference (77-80) better
  # than the mixed spectrum (50-55) vs reference (77-80)
  expect_true(r2$cosine[1L] > r1$cosine[1L])
})

test_that("score_ei_spectra use_deconvolved=TRUE with multiple components", {
  skip_if_not_installed("OrgMassSpecR")

  p <- create_nts_project(tempdir())
  p <- set_gc_spectra(p, list(
    "FT001" = data.frame(mz = c(50, 55), intensity = c(10, 20))
  ))
  p <- set_gc_deconvolved(p, list(
    "FT001" = list(
      list(spectrum = data.frame(mz = c(77, 78), intensity = c(30, 40)),
           apex_rt = 10.0, n_ions = 2L),
      list(spectrum = data.frame(mz = c(105, 120), intensity = c(100, 50)),
           apex_rt = 12.0, n_ions = 2L)
    )
  ))

  lib <- methods::new("NTSDatabase",
    metadata = data.table::data.table(lab_id = c("R1", "R2"), name = c("A", "B")),
    spectra_long = data.table::data.table(
      lab_id = c("R1", "R1", "R2", "R2"),
      mz = c(77, 78, 105, 120),
      intensity = c(30, 40, 100, 50)
    ),
    features_long = data.table::data.table(),
    version = "0.1.0"
  )

  r <- score_ei_spectra(p, lib, min_score = 0, use_deconvolved = TRUE, prefilter_threshold = 0)
  expect_s3_class(r, "data.table")
  # Should have matches from both components
  expect_true(any(grepl("FT001_c1", r$feature_id)))
  expect_true(any(grepl("FT001_c2", r$feature_id)))
})

test_that("deconvolve_gc_spectra returns empty list without valid xdata", {
  p <- create_nts_project(tempdir())
  p <- set_gc_spectra(p, list("FT001" = data.frame(mz = 1:3, intensity = 10:12)))
  expect_error(deconvolve_gc_spectra(p), "XCMSnExp")
})

test_that("detect_gc_compounds runs on single mzML", {
  a500 <- file.path(getwd(), "data", "DLL", "mzml", "A500.mzML")
  skip_if_not(file.exists(a500), "A500.mzML not available")

  result <- detect_gc_compounds(a500, min_snr = 5, min_ions = 10)
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_true(all(c(
    "compound_id", "filename", "rt", "basepeak_mz",
    "basepeak_int", "spectrum", "n_ions"
  ) %in% names(result)))

  spec1 <- result$spectrum[[1]]
  expect_s3_class(spec1, "data.frame")
  expect_true(all(c("mz", "intensity") %in% names(spec1)))
  expect_true(max(spec1$intensity) <= 100)

  cat(sprintf("Detected %d compounds from A500\n", nrow(result)))
})

test_that("detect_gc_compounds respects avoid_mz", {
  a500 <- file.path(getwd(), "data", "DLL", "mzml", "A500.mzML")
  skip_if_not(file.exists(a500), "A500.mzML not available")

  result <- detect_gc_compounds(a500,
    min_snr = 5, min_ions = 5,
    avoid_mz = c(50:500)
  )
  expect_true(nrow(result) == 0L || all(!result$basepeak_mz %in% 50:500))
})

test_that("detect_gc_compounds handles non-existent file gracefully", {
  expect_error(
    detect_gc_compounds("nonexistent.mzML"),
    "No valid"
  )
})

test_that("run_gc_detection populates project slots", {
  skip_on_ci()
  a500 <- file.path(getwd(), "data", "DLL", "mzml", "A500.mzML")
  skip_if_not(file.exists(a500), "A500.mzML not available")

  prj <- create_nts_project(dirname(a500))
  prj <- run_gc_detection(prj, min_snr = 5, min_ions = 10)

  expect_true(length(prj@gc_spectra) > 0)
  expect_true(nrow(prj@feature_table) > 0)
  expect_equal(nrow(prj@feature_table), length(prj@gc_spectra))
})

test_that(".nts_baseline_correct removes low-frequency baseline", {
  eic <- c(10, 11, 9, 12, 100, 50, 8, 10, 11, 9)
  corrected <- NTSminer:::.nts_baseline_correct(eic, window = 5L)
  expect_true(all(corrected >= 0))
  expect_true(max(corrected) > 50)
})

test_that(".nts_smooth_eic reduces noise", {
  eic <- c(0, 0, 10, 5, 100, 8, 0, 0)
  smoothed <- NTSminer:::.nts_smooth_eic(eic, half_window = 2L)
  expect_equal(length(smoothed), length(eic))
  expect_true(max(smoothed) < max(eic))
})

# Phase 16: Full-chromatogram deconvolution tests ------------------------------

test_that("per-m/z EIC apex detection finds peaks invisible in summed TIC", {
  # Compound B has a weak but clear single-m/z peak that would be
  # invisible in a summed TIC dominated by compound A's many fragments
  n_scans <- 100
  eic_B <- numeric(n_scans)
  eic_B[65:75] <- c(0, 2, 5, 15, 50, 80, 50, 15, 5, 2, 0)

  apexes <- NTSminer:::.nts_find_eic_apexes(eic_B, min_snr = 3)
  expect_true(length(apexes) > 0)
  expect_equal(apexes[1], 70)
})

test_that(".nts_cluster_apexes separates two compounds at different RT", {
  # Compound A fragments cluster around RT 10.0, compound B around 15.0
  rt <- c(9.8, 10.0, 10.2, 14.8, 15.0, 15.3)
  clusters <- NTSminer:::.nts_cluster_apexes(rt, apex_tol = 1.0)
  expect_equal(clusters, c(1L, 1L, 1L, 2L, 2L, 2L))
})

test_that(".nts_cluster_apexes with tighter apex_tol separates close peaks", {
  # Two compounds 1.5s apart should be separate clusters with tol=1.0
  rt <- c(10.0, 10.3, 11.6, 11.9)
  clusters <- NTSminer:::.nts_cluster_apexes(rt, apex_tol = 1.0)
  expect_equal(clusters, c(1L, 1L, 2L, 2L))
})

test_that("detect_gc_compounds validates baseline_window parameter", {
  a500 <- file.path(getwd(), "data", "DLL", "mzml", "A500.mzML")
  skip_if_not(file.exists(a500), "A500.mzML not available")

  result <- detect_gc_compounds(a500,
    min_snr = 5, min_ions = 10, baseline_window = 51L)
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)

  spec1 <- result$spectrum[[1]]
  expect_true(max(spec1$intensity) <= 100)
})

test_that("detect_gc_compounds with apex_tol=1.0 runs on single file", {
  a500 <- file.path(getwd(), "data", "DLL", "mzml", "A500.mzML")
  skip_if_not(file.exists(a500), "A500.mzML not available")

  result <- detect_gc_compounds(a500,
    min_snr = 3, apex_tol = 1.0, min_ions = 5)
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_true(all(c(
    "compound_id", "filename", "rt", "basepeak_mz",
    "basepeak_int", "spectrum", "n_ions"
  ) %in% names(result)))
})

# Phase 17: EIC shape clustering tests -----------------------------------------

test_that(".nts_eic_shape_cluster separates compounds by EIC shape", {
  # Two compounds with different elution profiles
  # Compound 1: peaks at scan 20, Compound 2: peaks at scan 50
  n_scans <- 80
  eic_mat <- matrix(0, nrow = n_scans, ncol = 6)
  # Compound 1: 3 fragments peaking at scan 20
  pk1 <- c(1, 2, 5, 10, 20, 50, 100, 50, 20, 10, 5, 2, 1)
  for (c in 1:3) {
    eic_mat[14:26, c] <- pk1 * runif(1, 0.5, 1.5)
  }
  # Compound 2: 3 fragments peaking at scan 50 (different elution time)
  pk2 <- c(1, 3, 8, 15, 30, 60, 90, 60, 30, 15, 8, 3, 1)
  for (c in 4:6) {
    eic_mat[44:56, c] <- pk2 * runif(1, 0.5, 1.5)
  }

  clusters <- NTSminer:::.nts_eic_shape_cluster(
    eic_mat, min_cor = 0.85, min_ions = 2L
  )

  expect_true(length(clusters) >= 2L)
  # Verify cluster sizes
  sizes <- vapply(clusters, length, integer(1L))
  expect_true(all(sizes >= 2L))
  # Clusters should be disjoint
  all_cols <- unlist(clusters)
  expect_equal(length(all_cols), length(unique(all_cols)))
})

test_that(".nts_eic_shape_cluster groups identical EIC shapes together", {
  # Two perfectly correlated fragments should be grouped
  n_scans <- 50
  eic_mat <- matrix(0, nrow = n_scans, ncol = 4)
  peak <- c(0, 2, 5, 10, 20, 50, 100, 50, 20, 10, 5, 2, 0)
  start <- 20
  for (c in 1:2) {
    eic_mat[start:(start + length(peak) - 1), c] <- peak * runif(1, 0.5, 1.5)
  }
  # Add uncorrelated noise bins
  eic_mat[30:35, 3] <- c(1, 2, 3, 2, 1, 0)
  eic_mat[10:15, 4] <- c(2, 4, 6, 4, 2, 0)

  clusters <- NTSminer:::.nts_eic_shape_cluster(
    eic_mat, min_cor = 0.85, min_ions = 2L
  )

  # The two correlated fragments should be in the same cluster
  found_together <- any(vapply(clusters, function(cl) all(c(1L, 2L) %in% cl), logical(1L)))
  expect_true(found_together)
})

test_that(".nts_eic_shape_cluster handles edge cases", {
  # Empty matrix (< min_ions columns)
  empty_mat <- matrix(0, nrow = 10, ncol = 2)
  result <- NTSminer:::.nts_eic_shape_cluster(
    empty_mat, min_cor = 0.85, min_ions = 5L
  )
  expect_equal(length(result), 0L)

  # Single-column matrix (should return empty if min_ions > 1)
  single_mat <- matrix(c(1:10), nrow = 10, ncol = 1)
  result <- NTSminer:::.nts_eic_shape_cluster(
    single_mat, min_cor = 0.85, min_ions = 5L
  )
  expect_equal(length(result), 0L)

  # All-identical EICs → one cluster
  n_scans <- 20
  all_same <- matrix(0, nrow = n_scans, ncol = 5)
  pk <- c(1, 2, 5, 10, 25, 50, 80, 50, 25, 10, 5, 2, 1)
  for (c in 1:5) {
    all_same[4:16, c] <- pk * runif(1, 0.5, 1.5)
  }
  clusters <- NTSminer:::.nts_eic_shape_cluster(
    all_same, min_cor = 0.85, min_ions = 2L
  )
  expect_equal(length(clusters), 1L)
  expect_equal(length(clusters[[1L]]), 5L)
})

test_that("detect_gc_compounds accepts mz_bin_width and shape_min_cor", {
  a500 <- file.path(getwd(), "data", "DLL", "mzml", "A500.mzML")
  skip_if_not(file.exists(a500), "A500.mzML not available")

  result <- detect_gc_compounds(a500,
    min_snr = 5, min_ions = 10,
    mz_bin_width = 0.1, shape_min_cor = 0.80
  )
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_true(all(c(
    "compound_id", "filename", "rt", "basepeak_mz",
    "basepeak_int", "spectrum", "n_ions"
  ) %in% names(result)))
})

# Phase 18: Peak quality tests ------------------------------------------------

test_that(".nts_peak_fwhm computes correct width for Gaussian peak", {
  # Gaussian-like EIC centered at scan 25
  n <- 50
  x <- seq_len(n)
  center <- 25; sigma <- 3
  eic <- exp(-(x - center)^2 / (2 * sigma^2)) * 100
  apex_idx <- which.max(eic)

  fwhm_info <- NTSminer:::.nts_peak_fwhm(eic, apex_idx)
  expect_true(is.finite(fwhm_info["fwhm"]))
  # FWHM should be roughly 2.355 * sigma ≈ 7 scans
  expect_true(fwhm_info["fwhm"] >= 4 && fwhm_info["fwhm"] <= 12)
  expect_true(fwhm_info["left"] < apex_idx)
  expect_true(fwhm_info["right"] > apex_idx)
})

test_that(".nts_peak_fwhm handles edge cases", {
  # Peak at edge of EIC
  eic <- c(5, 20, 100, 60, 20, 5, 0, 0, 0, 0)
  fwhm_info <- NTSminer:::.nts_peak_fwhm(eic, 1L)
  # Apex at index 1, left boundary can't extend
  expect_true(is.na(fwhm_info["fwhm"]) || is.finite(fwhm_info["fwhm"]))

  # Zero apex value
  fwhm_zero <- NTSminer:::.nts_peak_fwhm(c(0, 0, 0), 2L)
  expect_true(all(is.na(fwhm_zero)))

  # Out-of-bounds apex
  fwhm_oob <- NTSminer:::.nts_peak_fwhm(1:5, 10L)
  expect_true(all(is.na(fwhm_oob)))
})

test_that(".nts_peak_symmetry gives ~1 for symmetric peak", {
  n <- 50
  x <- seq_len(n)
  center <- 25; sigma <- 3
  eic <- exp(-(x - center)^2 / (2 * sigma^2)) * 100
  apex_idx <- which.max(eic)

  sym <- NTSminer:::.nts_peak_symmetry(eic, apex_idx)
  expect_true(sym >= 0.7)  # near-symmetric
})

test_that(".nts_peak_symmetry detects tailing peak", {
  # Right-tailing peak: sharp rise, slow fall
  n <- 30
  eic <- c(0, 1, 3, 8, 20, 50, 100, 80, 60, 45, 35, 28, 22, 18, 14, 12, 10,
           8, 7, 6, 5, 4, 3, 3, 2, 2, 1, 1, 0, 0)
  apex_idx <- which.max(eic)  # index 7, value 100

  sym <- NTSminer:::.nts_peak_symmetry(eic, apex_idx)
  # Left half (7→4) ≈ 3, Right half (7→30) ≈ 23
  # symmetry ≈ 3/23 ≈ 0.13, which is < 0.6
  expect_true(sym < 0.6)
  expect_true(sym > 0)
})

test_that("detect_gc_compounds runs with new Phase 18 defaults", {
  a500 <- file.path(getwd(), "data", "DLL", "mzml", "A500.mzML")
  skip_if_not(file.exists(a500), "A500.mzML not available")

  result <- detect_gc_compounds(a500,
    min_snr = 2, min_ions = 3, peak_ratio_min = 2,
    mz_bin_width = 0.1
  )
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  # With relaxed params, should detect more compounds
  expect_true(all(c(
    "compound_id", "filename", "rt", "basepeak_mz",
    "basepeak_int", "spectrum", "n_ions"
  ) %in% names(result)))
})

# erah integration tests -------------------------------------------------------

test_that("run_erah_detection populates project slots", {
  skip_if_not_installed("erah")
  skip_on_ci()
  mzml_dir <- file.path(getwd(), "data", "DLL", "mzml")
  skip_if_not(dir.exists(mzml_dir), "mzML data directory not available")

  prj <- create_nts_project(mzml_dir)
  prj <- run_erah_detection(prj)

  expect_true(length(prj@gc_spectra) > 0)
  expect_true(nrow(prj@feature_table) > 0)
  expect_equal(nrow(prj@feature_table), length(prj@gc_spectra))
})

test_that(".nts_parse_erah_spectrum parses correctly", {
  s <- "91,1000 92,394 105,385 106,138 175,7"
  sp <- NTSminer:::.nts_parse_erah_spectrum(s)
  expect_s3_class(sp, "data.frame")
  expect_equal(names(sp), c("mz", "intensity"))
  expect_equal(sp$mz[1], 91)
  expect_equal(sp$intensity[1], 100)  # 1000/10 rounded
})

test_that(".nts_parse_erah_spectrum handles edge cases", {
  expect_null(NTSminer:::.nts_parse_erah_spectrum(""))
  expect_null(NTSminer:::.nts_parse_erah_spectrum(NA_character_))
  expect_null(NTSminer:::.nts_parse_erah_spectrum("1,2"))
})
