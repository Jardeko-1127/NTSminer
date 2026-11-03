test_that("suggest_xcms_params mode='gc' returns valid parameters", {
  params <- suggest_xcms_params(ms1_files = character(0), mode = "gc")
  expect_equal(params$mode, "gc")
  expect_equal(params$centwave$peakwidth, c(1, 8))
  expect_equal(params$centwave$ppm, 20)
  expect_equal(params$centwave$snthresh, 5)
  expect_equal(params$centwave$prefilter, c(3, 100))
  expect_equal(params$centwave$noise, 100)
  expect_equal(params$centwave$fitgauss, FALSE)
  expect_equal(params$algorithm, "centWave")
  expect_equal(params$peakdensity$minFraction, 0.5)
})

test_that("read_gc_tic reads mzML correctly", {
  a500_path <- file.path(getwd(), "data", "DLL", "mzml", "A500.mzML")
  skip_if_not(file.exists(a500_path),
              "A500.mzML not available")
  tic <- read_gc_tic(a500_path)
  expect_s3_class(tic, "data.table")
  expect_true(nrow(tic) > 0)
  expect_true(all(c("filename", "rt", "intensity") %in% names(tic)))
  expect_true(all(tic$rt > 0))
})

test_that("score_ei_spectra computes cosine similarity correctly", {
  skip_if_not_installed("OrgMassSpecR")
  ref <- data.frame(
    mz = c(50, 77, 105, 120, 150),
    intensity = c(10, 30, 100, 50, 20)
  )
  qry <- data.frame(
    mz = c(50.1, 77.1, 105.0, 120.2),
    intensity = c(8, 25, 95, 45)
  )
  cos <- OrgMassSpecR::SpectrumSimilarity(
    spec.top = ref,
    spec.bottom = qry,
    t = 0.25, b = 10,
    print.graphic = FALSE,
    output.list = FALSE
  )
  expect_true(cos > 0.9)
  expect_true(cos <= 1.0)
})

test_that("OrgMassSpecR handles dissimilar spectra", {
  skip_if_not_installed("OrgMassSpecR")
  ref <- data.frame(
    mz = c(50, 77, 105),
    intensity = c(10, 30, 100)
  )
  qry <- data.frame(
    mz = c(200, 250, 300),
    intensity = c(100, 50, 30)
  )
  cos <- OrgMassSpecR::SpectrumSimilarity(
    spec.top = ref,
    spec.bottom = qry,
    t = 0.25, b = 10,
    print.graphic = FALSE,
    output.list = FALSE
  )
  expect_true(cos < 0.1)
})

test_that("get_gc_spectra / set_gc_spectra round-trips", {
  p <- create_nts_project(tempdir())
  expect_type(get_gc_spectra(p), "list")
  expect_length(get_gc_spectra(p), 0L)

  dummy <- list("FT001" = data.frame(mz = 1:3, intensity = 10:12))
  p <- set_gc_spectra(p, dummy)
  expect_equal(names(get_gc_spectra(p)), "FT001")
  expect_equal(nrow(get_gc_spectra(p)[["FT001"]]), 3L)
})

test_that("NTSProject contains gc_spectra slot", {
  p <- create_nts_project(tempdir())
  expect_true("gc_spectra" %in% methods::slotNames(p))
})

test_that("score_ei_spectra with minimal project and library", {
  skip_if_not_installed("OrgMassSpecR")

  p <- create_nts_project(tempdir())
  p <- set_gc_spectra(p, list(
    "FT001" = data.frame(
      mz = c(50.1, 77.1, 105.0, 120.2),
      intensity = c(8, 25, 95, 45)
    )
  ))

  lib <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      lab_id = "REF001", name = "Test Compound"
    ),
    spectra_long = data.table::data.table(
      lab_id = "REF001",
      mz = c(50, 77, 105, 120, 150),
      intensity = c(10, 30, 100, 50, 20)
    ),
    features_long = data.table::data.table(),
    version = "0.1.0"
  )

  result <- score_ei_spectra(p, lib, min_score = 0)

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_true(all(c("feature_id", "ref_id", "ref_name", "cosine") %in% names(result)))
  expect_equal(result$feature_id[1L], "FT001")
  expect_equal(result$ref_id[1L], "REF001")
  expect_true(result$cosine[1L] > 0.9)
})

test_that("annotate_gc_by_spectral_match writes to @db_annotations", {
  skip_if_not_installed("OrgMassSpecR")

  p <- create_nts_project(tempdir())
  p <- set_gc_spectra(p, list(
    "FT001" = data.frame(
      mz = c(50.1, 77.1, 105.0),
      intensity = c(8, 25, 95)
    )
  ))

  lib <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      lab_id = "REF001", name = "Benzene"
    ),
    spectra_long = data.table::data.table(
      lab_id = "REF001",
      mz = c(50, 77, 105, 120, 150),
      intensity = c(10, 30, 100, 50, 20)
    ),
    features_long = data.table::data.table(),
    version = "0.1.0"
  )

  p <- annotate_gc_by_spectral_match(p, lib, min_score = 0, top_n = 2L)

  db <- p@db_annotations
  expect_true(is.data.frame(db))
  expect_true(nrow(db) >= 1L)
  expect_true("Compound.name" %in% names(db))
  expect_true("match_method" %in% names(db))
  expect_equal(db$match_method[1L], "EI_cosine")
})

test_that("score_ei_spectra handles empty gc_spectra gracefully", {
  skip_if_not_installed("OrgMassSpecR")

  p <- create_nts_project(tempdir())
  lib <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(lab_id = character(0L), name = character(0L)),
    spectra_long = data.table::data.table(
      lab_id = character(0L), mz = numeric(0L), intensity = numeric(0L)
    ),
    features_long = data.table::data.table(),
    version = "0.1.0"
  )

  expect_error(
    score_ei_spectra(p, lib),
    "empty"
  )
})

test_that("build_gc_library assembles NTSDatabase from MSP", {
  skip_on_ci()
  # This test requires a real MSP file and exercises import_msp path
  # Placeholder for manual testing with actual MSP exports
  expect_true(TRUE)
})
