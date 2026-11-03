.nts_mock_classyfire_query <- function() {
  jsonlite::fromJSON(
    testthat::test_path("fixtures", "classyfire_cco.json"),
    simplifyVector = FALSE
  )
}

.nts_with_mock_classyfire <- function(expr) {
  fx <- .nts_mock_classyfire_query()
  testthat::with_mocked_bindings(
    .nts_classyfire_submit = function(smiles, ids = NULL, ...) {
      list(id = fx$id, label = fx$label)
    },
    .nts_classyfire_wait_query = function(query_id, ...) {
      fx
    },
    .package = "NTSminer",
    expr
  )
}

test_that("classyfire entity parser extracts ChemOnt levels", {
  fx <- .nts_mock_classyfire_query()
  out <- NTSminer:::.nts_classyfire_parse_query_json(fx, input_ids = "MOL1")
  expect_equal(nrow(out), 1L)
  expect_equal(out$identifier, "MOL1")
  expect_equal(out$Kingdom, "Organic compounds")
  expect_equal(out$Subclass, "Primary alcohols")
  expect_equal(out$classyfire_status, "classified")
})

test_that("classify_classyfire works with mocked API", {
  out <- .nts_with_mock_classyfire({
    NTSminer::classify_classyfire(
      c(MOL1 = "CCO"),
      canonicalize = FALSE,
      progress = FALSE,
      prefer_inchikey_lookup = FALSE
    )
  })
  expect_equal(nrow(out), 1L)
  expect_equal(out$Superclass, "Organooxygen compounds")
})

test_that("unclassified inputs are retained in output", {
  fx <- .nts_mock_classyfire_query()
  out <- NTSminer:::.nts_classyfire_parse_query_json(fx, input_ids = c("MOL1", "BAD1"))
  expect_equal(nrow(out), 2L)
  expect_equal(out$classyfire_status[out$identifier == "BAD1"], "unclassified")
  expect_equal(out$cf_status[out$identifier == "BAD1"], "not_found")
})

test_that("InChIKey lookup skips POST when all compounds found", {
  fx <- .nts_mock_classyfire_query()
  ent <- fx$entities[[1L]]
  submit_called <- FALSE
  out <- testthat::with_mocked_bindings(
    .nts_toxpi_inchikeys = function(smiles) {
      rep("LFQSCWFLJHTTHZ-UHFFFAOYSA-N", length(smiles))
    },
    .nts_classyfire_entity_by_inchikey = function(inchi_key, ...) ent,
    .nts_classyfire_submit = function(...) {
      submit_called <<- TRUE
      stop("POST should not run")
    },
    .package = "NTSminer",
    {
      NTSminer::classify_classyfire(
        c(MOL1 = "CCO"),
        canonicalize = FALSE,
        progress = FALSE,
        prefer_inchikey_lookup = TRUE,
        fallback_single = FALSE
      )
    }
  )
  expect_false(submit_called)
  expect_equal(out$cf_status, "found")
  expect_equal(out$Kingdom, "Organic compounds")
})

test_that("POST failure with fallback_single FALSE returns lookup hits", {
  fx <- .nts_mock_classyfire_query()
  ent <- fx$entities[[1L]]
  run_cf <- function() {
    NTSminer::classify_classyfire(
      c(A = "CCO", B = "invalid"),
      compound_ids = c("A", "B"),
      canonicalize = FALSE,
      progress = FALSE,
      prefer_inchikey_lookup = TRUE,
      fallback_single = FALSE
    )
  }
  with_cf_mocks <- function(expr) {
    testthat::with_mocked_bindings(
      .nts_toxpi_inchikeys = function(smiles) {
        c("KEY_FOUND", rep(NA_character_, max(0L, length(smiles) - 1L)))
      },
      .nts_classyfire_entity_by_inchikey = function(inchi_key, ...) {
        if (inchi_key == "KEY_FOUND") ent else NULL
      },
      .nts_classyfire_submit = function(...) {
        structure(list(ok = FALSE, http_status = 500L), class = "nts_classyfire_submit_failed")
      },
      .package = "NTSminer",
      expr
    )
  }
  out <- with_cf_mocks(run_cf())
  testthat::expect_warning(
    with_cf_mocks(run_cf()),
    "fallback_single = FALSE"
  )
  expect_equal(out$cf_status[out$identifier == "A"], "found")
  expect_equal(out$cf_status[out$identifier == "B"], "not_found")
})
