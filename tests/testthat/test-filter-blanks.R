test_that("filter_blanks keeps rows when all blank columns are NA", {
  ft <- data.frame(
    name = c("FT1", "FT2"),
    mz = c(300, 400),
    rt = c(100, 200),
    Sample_A = c(1e5, NA_real_),
    ME1 = c(NA_real_, 1e4),
    ME2 = c(NA_real_, 1e4),
    stringsAsFactors = FALSE
  )
  out <- filter_blanks(
    ft,
    sample_cols = "Sample_A",
    blank_cols = c("ME1", "ME2"),
    threshold = 5,
    mode = "any"
  )
  expect_equal(out$name, "FT1")
})

test_that("filter_blanks still drops sample below blank fold-change", {
  ft <- data.frame(
    name = "FT1",
    mz = 300,
    rt = 100,
    Sample_A = 1e4,
    ME1 = 1e4,
    stringsAsFactors = FALSE
  )
  out <- filter_blanks(
    ft,
    sample_cols = "Sample_A",
    blank_cols = "ME1",
    threshold = 5,
    mode = "any"
  )
  expect_equal(nrow(out), 0L)
})
