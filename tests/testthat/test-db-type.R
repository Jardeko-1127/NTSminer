# db_type provenance: builder stamping, set_db_type(), confidence detection ---
# Regression coverage for finding A1 of the 0.1.19 accuracy review: an in-house
# library with no `db_type` column and non-matching Lab_IDs was silently judged
# "online", making the in-house tiers L1/L2 unreachable.

test_that("build_nts_database stamps db_type only when provenance is known", {
  md <- data.frame(Lab_ID = "L1", mass = 300, Polarity = "NEG",
                   SMILES = "CCCC", stringsAsFactors = FALSE)
  sl <- data.frame(Lab_ID = "L1", NCE = "CE10",
                   mz = c(100, 150), intensity = c(1000, 500))

  db_i <- suppressMessages(
    build_nts_database(md, sl, add_merged_layer = FALSE, db_type = "inhouse")
  )
  expect_true("db_type" %in% names(db_i@metadata))
  expect_true(all(db_i@metadata$db_type == "inhouse"))

  # Unknown provenance -> column omitted (NOT silently defaulted to inhouse).
  db_u <- suppressMessages(
    build_nts_database(md, sl, add_merged_layer = FALSE)
  )
  expect_false("db_type" %in% names(db_u@metadata))

  # An explicit metadata column takes precedence over the argument.
  md2 <- md
  md2$db_type <- "online"
  db_c <- suppressMessages(
    build_nts_database(md2, sl, add_merged_layer = FALSE, db_type = "inhouse")
  )
  expect_true(all(db_c@metadata$db_type == "online"))

  # Unrecognised argument values are ignored -> column omitted.
  db_bad <- suppressMessages(
    build_nts_database(md, sl, add_merged_layer = FALSE, db_type = "bogus")
  )
  expect_false("db_type" %in% names(db_bad@metadata))
})

test_that("set_db_type validates and stamps provenance", {
  md <- data.frame(Lab_ID = "PFAS_L01", mass = 300, Polarity = "NEG",
                   SMILES = "CCCC", stringsAsFactors = FALSE)
  sl <- data.frame(Lab_ID = "PFAS_L01", NCE = "CE10",
                   mz = c(100, 150), intensity = c(1000, 500))
  db <- suppressMessages(build_nts_database(md, sl, add_merged_layer = FALSE))

  db2 <- set_db_type(db, "inhouse")
  expect_true("db_type" %in% names(db2@metadata))
  expect_true(all(db2@metadata$db_type == "inhouse"))
  # Original is untouched (no by-reference mutation).
  expect_false("db_type" %in% names(db@metadata))

  expect_error(set_db_type(db, "bogus"), "must be one of")
  expect_error(set_db_type("not a db", "inhouse"), "NTSDatabase")
})

test_that(".nts_confidence_db_type prefers the column, falls back to ID pattern", {
  dbt <- NTSminer:::.nts_confidence_db_type

  # Explicit column wins even when the ID looks like an online entry.
  expect_equal(
    dbt(data.frame(Lab_ID = "MoNA_1", db_type = "inhouse")),
    "inhouse"
  )
  # Column present but NA/unknown -> fall back to the ID pattern.
  expect_equal(
    dbt(data.frame(Lab_ID = c("NTS_InHouse_PFAS_01", "pred_x"),
                   db_type = c(NA_character_, NA_character_))),
    c("inhouse", "prediction")
  )

  # No db_type column at all -> warn once, then infer from Lab_ID.
  old <- getOption("NTSminer.db_type_warned")
  on.exit(options(NTSminer.db_type_warned = old), add = TRUE)
  options(NTSminer.db_type_warned = FALSE)

  expect_warning(
    out <- dbt(data.frame(Lab_ID = c("PFAS_L01", "NTS_InHouse_X"))),
    "no `db_type` column"
  )
  expect_equal(out, c("online", "inhouse"))

  # The warning is one-shot: a second call in the same session is silent.
  expect_silent(dbt(data.frame(Lab_ID = c("PFAS_L01"))))
})

test_that("in-house db_type reaches L1 via propagation path", {
  # Mirrors the real db_annotations shape: a db_type column carried through from
  # the library metadata. A non-matching Lab_ID (PFAS_L*) must still hit L1,
  # regardless of RT agreement (library tiers are provenance-only).
  d <- tempfile("nts_conf_dbt")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(name = c("A", "B"), mz = c(100, 200), rt = c(1, 2),
               stringsAsFactors = FALSE)
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      Lab_ID = c("PFAS_L01", "PFAS_L02"),
      db_type = c("inhouse", "inhouse"),
      lib_rt = c(1.2, 2.1),
      rt_error_s = c(0.5, 400),
      score_total = c(0.95, 0.95),
      stringsAsFactors = FALSE
    )
  )
  p <- assign_confidence_levels(p, rt_tol_sec = 30)
  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L1")  # RT agrees
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L1")  # RT off -> still L1
})
