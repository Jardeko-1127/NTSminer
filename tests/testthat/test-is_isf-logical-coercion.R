## Regression tests: logical / NA handling for `is_isf` and related graph logic.
## Motivation: `is_isf == TRUE` vs `is_isf %in% TRUE` and CSV-derived types differ;
## `NA` must not create spurious ISF links; evidence tiers must stay stable.

test_that("%in% TRUE never yields NA (safe data.table `i` rows)", {
  x <- c(TRUE, NA, FALSE, 1L, 0L)
  r <- x %in% TRUE
  expect_false(anyNA(r))
  expect_identical(r, c(TRUE, FALSE, FALSE, TRUE, FALSE))
})

test_that(".nts_build_isf_edges_dt links only is_isf %in% TRUE with finite parent_mz", {
  skip_if_not_installed("data.table")
  prec <- 300.1234
  dt <- data.table::data.table(
    name = c("parent", "child_na", "child_ok", "child_false"),
    mz = c(prec, 280.1, prec - 0.0005, 270.0),
    is_isf = c(FALSE, NA, TRUE, FALSE),
    parent_mz = c(NA_real_, prec, prec, prec)
  )
  ed <- NTSminer:::.nts_build_isf_edges_dt(dt, mz_tol_ppm = 20)
  expect_s3_class(ed, "data.table")
  expect_true(nrow(ed) >= 1L)
  expect_true(any(ed$to == "child_ok"))
  expect_false(any(ed$to == "child_na"))
  expect_false(any(ed$to == "child_false"))
})

test_that(".nts_build_isf_edges_dt: character 'TRUE' matches via %in% (coercion); prefer as.logical for CSV", {
  skip_if_not_installed("data.table")
  prec <- 400.0
  dt <- data.table::data.table(
    name = c("parent", "child_str"),
    mz = c(prec, 380.0),
    is_isf = c("FALSE", "TRUE"),
    parent_mz = c(NA_real_, prec)
  )
  ed <- NTSminer:::.nts_build_isf_edges_dt(dt, mz_tol_ppm = 20)
  ## `%in%` uses match(); character "TRUE" matches logical TRUE after coercion.
  expect_true(nrow(ed) >= 1L)
  expect_true(any(ed$to == "child_str"))
})

test_that(".nts_network_nodes_evidence_group uses strict ISF tier (is_isf %in% TRUE)", {
  skip_if_not_installed("data.table")
  nodes_dt <- data.table::data.table(
    id = paste0("n", 1:4),
    name = paste0("n", 1:4),
    Library_ID = NA_character_,
    Compound.name = NA_character_,
    annotation_level = NA_character_,
    Propagation_Level = NA_integer_,
    diagnostic_fragments = NA_character_,
    n_fragments = 0,
    mz_remainder = NA_character_,
    is_isf = c(FALSE, NA, TRUE, FALSE)
  )
  grp <- NTSminer:::.nts_network_nodes_evidence_group(NULL, nodes_dt)
  expect_identical(length(grp), nrow(nodes_dt))
  expect_identical(grp[[3L]], "ISF")
  expect_false(grp[[2L]] == "ISF")
  expect_false(grp[[1L]] == "ISF")
  expect_false(grp[[4L]] == "ISF")
})
