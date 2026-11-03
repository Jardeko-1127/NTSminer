.nts_mock_episuite_response <- function() {
  jsonlite::fromJSON(
    testthat::test_path("fixtures", "episuite_cco.json"),
    simplifyVector = FALSE
  )
}

.nts_with_mock_episuite <- function(expr) {
  fx_body <- .nts_mock_episuite_response()
  testthat::with_mocked_bindings(
    .nts_episuite_call = function(smiles, ...) {
      fx_body
    },
    .package = "NTSminer",
    expr
  )
}

test_that("EPI Suite response parser produces expected columns", {
  fx <- .nts_mock_episuite_response()
  row <- NTSminer:::.nts_episuite_row_from_response(fx, "CCO")
  expect_true("LogKow_est" %in% names(row))
  expect_true("Biowin1" %in% names(row))
  expect_true("Fish_LC50" %in% names(row))
  expect_equal(row$LogP_pred, row$LogKow_est)
  expect_equal(row$SMILES, "CCO")
})

test_that("analyze_episuite_biodegradation classifies fixture BIOWIN models", {
  fx <- .nts_mock_episuite_response()
  bio <- NTSminer::analyze_episuite_biodegradation(fx)
  expect_equal(nrow(bio$bio_results), 7L)
  expect_equal(bio$bio_results$Classification[1L], "Biodegrades fast")
  expect_equal(bio$bio_results$Classification[3L], "About Days")
  expect_true(is.data.frame(bio$prediction_result))
  expect_match(bio$prediction_result$Results[1L], "YES")
})

test_that("predict_episuite works with mocked EPI call", {
  out <- .nts_with_mock_episuite({
    NTSminer::predict_episuite("CCO", progress = FALSE, delay_sec = c(0, 0))
  })

  expect_equal(nrow(out), 1L)
  expect_true(is.finite(out$LogKow_est[1L]))
})

test_that("predict_properties accepts SMILES character vector", {
  out <- .nts_with_mock_episuite({
    NTSminer::predict_properties("CCO", progress = FALSE, delay_sec = c(0, 0))
  })
  expect_true(is.data.frame(out))
  expect_true("LogKoc_est" %in% names(out))
})

test_that("predict_admet_risk merges backends structure", {
  merged <- NTSminer:::.nts_merge_admet_backends(
    data.frame(SMILES = "CCO", AMES = 0.1, admetlab_backend = "ADMETlab3", stringsAsFactors = FALSE),
    data.frame(SMILES = "CCO", hERG = 0.2, admet_ai_backend = "ADMET-AI", stringsAsFactors = FALSE)
  )
  expect_equal(nrow(merged), 1L)
  expect_true(all(c("AMES", "hERG") %in% names(merged)))
})

test_that("admet_ai row parser handles native R list like toxPi cbind", {
  mock <- list(AMES = 0.1, hERG = 0.2, DILI = 0.05)
  out <- NTSminer:::.nts_admet_ai_row_from_predict(mock, "CCO")
  expect_equal(out$SMILES, "CCO")
  expect_equal(out$AMES, 0.1)
  expect_equal(out$hERG, 0.2)
})

test_that("NTSProject predict_properties writes physchem_proxy when SMILES present", {
  d <- tempfile("nts_epi_proj")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- NTSminer::create_nts_project(d)
  p <- NTSminer::set_feature_table(
    p,
    data.frame(name = "F1", mz = 46, rt = 1, stringsAsFactors = FALSE)
  )
  ms1 <- data.frame(
    name = "F1", mz = 46, rt = 1,
    master_smiles = "CCO",
    stringsAsFactors = FALSE
  )
  p <- NTSminer::set_ms1_annotations(p, ms1)

  p2 <- .nts_with_mock_episuite({
    NTSminer::predict_properties(p, progress = FALSE, delay_sec = c(0, 0))
  })

  proxy <- NTSminer::get_physchem_proxy(p2)
  expect_true("LogKow_est" %in% names(proxy))
  expect_true(is.finite(proxy$LogKow_est[proxy$name == "F1"]))
})
