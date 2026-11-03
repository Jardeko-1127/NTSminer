test_that("summarize_subclass_fragment_features prevalence and plot", {
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = c("A", "B", "C"),
      Second_Class = c("X", "X", "Y")
    ),
    spectra_long = data.table::data.table(),
    features_long = data.table::data.table(
      Lab_ID = c("A", "A", "B", "C"),
      Polarity = "NEG",
      Feature_Type = c("DF", "DF", "DF", "NL"),
      mz_value = c(100, 100.0001, 100, 50),
      Frequency = c(0.9, 0.8, 0.7, 0.6)
    ),
    version = "test"
  )
  s <- summarize_subclass_fragment_features(
    db,
    class_col = "Second_Class",
    min_prevalence = 0.5,
    mz_tol_ppm = NULL,
    mz_digits = 3L
  )
  expect_true(nrow(s$by_feature) >= 1L)
  expect_true("subclass_coverage" %in% names(s))
  expect_no_error(plot_subclass_fragment_features(s, print_plot = FALSE))
})

test_that("summarize_subclass_fragment_features: lowercase column names normalization", {
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      lab_id = c("A", "B", "C"),
      second_class = c("X", "X", "Y")
    ),
    spectra_long = data.table::data.table(),
    features_long = data.table::data.table(
      lab_id = c("A", "A", "B", "C"),
      polarity = "NEG",
      feature_type = c("df", "df", "df", "nl"),
      mz = c(100, 100.0001, 100, 50),
      frequency = as.factor(c("0.9", "0.8", "0.7", "0.6"))
    ),
    version = "test"
  )
  s <- summarize_subclass_fragment_features(
    db,
    class_col = "second_class",
    lab_id_col = "lab_id",
    min_prevalence = 0.5,
    mz_tol_ppm = NULL,
    mz_digits = 3L
  )
  expect_true(nrow(s$by_feature) >= 1L)
  expect_true(all(s$by_feature$prevalence >= 0.5))
})

test_that("summarize_subclass_fragment_features: factor Frequency handled correctly", {
  db <- methods::new(
    "NTSDatabase",
    metadata = data.table::data.table(
      Lab_ID = c("A", "B"),
      Second_Class = c("X", "X")
    ),
    spectra_long = data.table::data.table(),
    features_long = data.table::data.table(
      Lab_ID = c("A", "B"),
      Polarity = "NEG",
      Feature_Type = c("DF", "DF"),
      mz_value = c(68.995, 68.995),
      Frequency = factor(c("0.95", "0.85"))
    ),
    version = "test"
  )
  s <- summarize_subclass_fragment_features(
    db,
    min_prevalence = 0.5,
    mz_tol_ppm = 10
  )
  expect_true(nrow(s$by_feature) >= 1L)
  expect_equal(s$by_feature$prevalence[1L], 1)
})
