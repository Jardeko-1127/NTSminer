test_that("run_toxpi works with toxpiR when installed", {
  skip_if_not_installed("toxpiR")
  prep <- NTSminer::prepare_toxpi(
    .toxpi_fixture_epi(),
    .toxpi_fixture_admet(),
    scenario = "environmental",
    verbose = FALSE
  )
  res <- NTSminer::run_toxpi(prep, weight_multiplier = 1)
  sc <- NTSminer::toxpi_scores(res)
  expect_equal(length(sc), nrow(prep$data))
  expect_equal(names(sc), prep$data$Name)
  mat <- NTSminer::toxpi_slice_scores(res, adjusted = FALSE)
  expect_equal(nrow(mat), nrow(prep$data))
  expect_equal(ncol(mat), length(names(prep$weights)))
})
