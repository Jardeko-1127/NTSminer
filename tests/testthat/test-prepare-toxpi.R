test_that("validate_toxpi_input passes fixture tables", {
  expect_true(
    NTSminer::validate_toxpi_input(.toxpi_fixture_epi(), .toxpi_fixture_admet())
  )
})

test_that("validate_toxpi_input stops on missing EPI columns", {
  epi <- .toxpi_fixture_epi()
  epi$Fish_LC50 <- NULL
  expect_error(
    NTSminer::validate_toxpi_input(epi, .toxpi_fixture_admet()),
    "Missing required columns"
  )
})

test_that("Name coalesces from lowercase name when episuite_Name is all NA", {
  epi <- .toxpi_fixture_epi()
  epi$episuite_Name <- NA_character_
  epi$name <- c("TFA", "PFPrA", "PFOA")
  res <- NTSminer::prepare_toxpi(
    epi,
    .toxpi_fixture_admet(),
    verbose = FALSE
  )
  expect_equal(res$data$Name, c("TFA", "PFPrA", "PFOA"))
  expect_equal(res$data$SMILES, .toxpi_fixture_epi()$SMILES)
})

test_that("compound_ids overrides missing episuite_Name", {
  epi <- .toxpi_fixture_epi()
  epi$episuite_Name <- NA_character_
  res <- NTSminer::prepare_toxpi(
    epi,
    .toxpi_fixture_admet(),
    compound_ids = c("A", "B", "C"),
    verbose = FALSE
  )
  expect_equal(res$data$Name, c("A", "B", "C"))
})

test_that("prepare_toxpi returns six slices for environmental scenario", {
  res <- NTSminer::prepare_toxpi(
    .toxpi_fixture_epi(),
    .toxpi_fixture_admet(),
    scenario = "environmental",
    tox21_agg = "max",
    verbose = FALSE
  )
  expect_equal(nrow(res$data), 3L)
  expect_true(all(c(
    "Persistence", "BioDegradation", "AquaticTox",
    "SystemicTox", "Genotox", "Tox21"
  ) %in% names(res$data)))
  expect_equal(sum(res$weights), 1, tolerance = 1e-8)
  expect_equal(as.numeric(res$weights["Persistence"]), 0.25, tolerance = 1e-8)
  expect_true(all(res$data$Persistence >= 0 & res$data$Persistence <= 1))
})

test_that("prepare_toxpi split tox21 adds NR and SR slices", {
  res <- NTSminer::prepare_toxpi(
    .toxpi_fixture_epi(),
    .toxpi_fixture_admet(),
    tox21_agg = "split",
    verbose = FALSE
  )
  expect_true(all(c("Tox21_NR", "Tox21_SR") %in% names(res$data)))
  expect_equal(length(res$slices$Tox21_NR), 7L)
  expect_equal(length(res$slices$Tox21_SR), 5L)
})

test_that("prepare_toxpi writes csv format", {
  out <- tempfile(fileext = ".csv")
  on.exit(unlink(out), add = TRUE)
  res <- NTSminer::prepare_toxpi(
    .toxpi_fixture_epi(),
    .toxpi_fixture_admet(),
    output_format = "csv",
    outfile = out,
    verbose = FALSE
  )
  expect_true(file.exists(out))
  lines <- readLines(out, warn = FALSE)
  expect_true(grepl("^Name,", lines[1L]))
  expect_true(grepl("^Weights,", lines[2L]))
  expect_equal(length(lines), 2L + nrow(res$data))
})

test_that("aquatic tox ranks higher toxicity for lower LC50", {
  res <- NTSminer::prepare_toxpi(
    .toxpi_fixture_epi(),
    .toxpi_fixture_admet(),
    verbose = FALSE
  )
  aq <- stats::setNames(res$data$AquaticTox, res$data$SMILES)
  expect_gt(aq["CCC"], aq["CCO"])
  expect_gt(aq["CC"], aq["CCO"])
})

test_that("order merge keeps nrow when EPI rewrites SMILES", {
  epi <- .toxpi_fixture_epi()
  admet <- .toxpi_fixture_admet()
  epi <- rbind(epi, epi[3L, , drop = FALSE])
  admet <- rbind(admet, admet[3L, , drop = FALSE])
  epi$SMILES[1L] <- "C(C)O"
  epi$SMILES[3L] <- "C(=O)(O)CCC"
  query <- c("CCO", "CC", "CCC", "CCC")

  res <- NTSminer::prepare_toxpi(
    epi,
    admet,
    query_smiles = query,
    align_by = "order",
    verbose = FALSE
  )
  expect_equal(nrow(res$data), 4L)
  expect_false(any(is.na(res$data$SystemicTox)))
  expect_false(any(is.na(res$data$AquaticTox)))
})

test_that("Biowin3 timeframe uses discrete EPI risk bins", {
  expect_equal(NTSminer:::.nts_biowin3_to_risk(5), 0)
  expect_equal(NTSminer:::.nts_biowin3_to_risk(4), 0.25, tolerance = 1e-8)
  expect_equal(NTSminer:::.nts_biowin3_to_risk(3), 0.5, tolerance = 1e-8)
  expect_equal(NTSminer:::.nts_biowin3_to_risk(1), 1, tolerance = 1e-8)
  out1 <- NTSminer:::.nts_toxpi_norm_biowin(0.8, "Biowin1", verbose = FALSE)
  expect_equal(out1, 0.2, tolerance = 1e-8)
})

test_that("Ready biodegradability risk uses raw Biowin3 and Biowin5 thresholds", {
  expect_equal(NTSminer:::.nts_biowin_ready_biodeg_risk(4.0, 0.85), 0)
  expect_equal(NTSminer:::.nts_biowin_ready_biodeg_risk(2.0, 0.85), 1)
  expect_equal(NTSminer:::.nts_biowin_ready_biodeg_risk(4.0, 0.4), 1)
})

test_that("BioDegradation slice has seven Biowins; Ready is separate diagnostic", {
  res <- NTSminer::prepare_toxpi(
    .toxpi_fixture_epi(),
    .toxpi_fixture_admet(),
    verbose = FALSE
  )
  expect_equal(length(res$slices$BioDegradation), 7L)
  expect_false("Ready_biodeg_risk" %in% res$slices$BioDegradation)
  expect_true("Ready_biodeg_risk" %in% names(res$data))
})

test_that("Genotox slice uses _norm column names", {
  res <- NTSminer::prepare_toxpi(
    .toxpi_fixture_epi(),
    .toxpi_fixture_admet(),
    verbose = FALSE
  )
  expect_equal(
    res$slices$Genotox,
    c("AMES_norm", "Carcinogens_Lagunin_norm")
  )
})

test_that("normalize minmax handles flat vectors", {
  x <- rep(5, 4)
  out <- NTSminer:::.nts_toxpi_normalize_minmax(x, "positive")
  expect_equal(out, rep(0.5, 4))
})
