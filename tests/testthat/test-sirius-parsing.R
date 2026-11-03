# Tests for SIRIUS parsing / CLI hardening (0.1.34).
#
# Scope: CANOPUS NPC column resolution + feature-keyed merge, the
# `explained_peaks` / `confidence` column mapping, `batch_size` validation, the
# `output_dir` unlink protection and the CLI failure path.
#
# NONE of these tests start SIRIUS, start a REST service, or touch the network:
# the CLI is stubbed at `.nts_system2_capture()`. Fixtures are the real
# SIRIUS 6.3.4 outputs copied into fixtures/sirius/ (read via test_path()).

sirius_fixture <- function(...) testthat::test_path("fixtures", "sirius", ...)

read_sirius_tsv <- function(...) {
  utils::read.table(sirius_fixture(...), header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE, quote = "\"",
                    comment.char = "", check.names = FALSE)
}

# Read a fixture as a data.table with the lower-cased column names that the
# parsing code works with.
read_sirius_dt <- function(...) {
  dt <- data.table::as.data.table(read_sirius_tsv(...))
  names(dt) <- tolower(names(dt))
  dt
}

# Same, plus the key-column renaming that compute_sirius_fingerprints() /
# parse_canopus_results() apply before merging.
canopus_fixture <- function(...) {
  dt <- read_sirius_dt(...)
  if ("formularank" %in% names(dt)) data.table::setnames(dt, "formularank", "rank")
  if ("molecularformula" %in% names(dt)) {
    data.table::setnames(dt, "molecularformula", "molecular_formula")
  }
  dt
}

# ---------------------------------------------------------------------------
# NPC delimiter regex
# ---------------------------------------------------------------------------

test_that("NPC delimiter regex splits on ':' ';' '/' with optional whitespace", {
  x <- strsplit("Alkaloids: Polyketides; Terpenoids/Lipids", "[\\:/;]\\s*")[[1]]
  expect_length(x, 4L)
  expect_identical(x, c("Alkaloids", "Polyketides", "Terpenoids", "Lipids"))
  # trimming the pieces yields the three usable levels
  expect_identical(trimws(x[1:3]), c("Alkaloids", "Polyketides", "Terpenoids"))
})

# ---------------------------------------------------------------------------
# NPC column resolution
# ---------------------------------------------------------------------------

test_that(".nts_canopus_npc_columns resolves multi / single / none column shapes", {
  hdr <- names(read_sirius_tsv("canopus_formula_summary.tsv"))

  r6 <- NTSminer:::.nts_canopus_npc_columns(hdr)
  expect_equal(r6$mode, "multi")
  expect_equal(r6$pathway, "NPC#pathway")
  expect_equal(r6$superclass, "NPC#superclass")
  expect_equal(r6$class, "NPC#class")
  # SIRIUS 6 has no NPC#subclass -> NA (expected, not a failure)
  expect_true(is.na(r6$subclass))

  r1 <- NTSminer:::.nts_canopus_npc_columns(c("name", "NPC", "ionMass"))
  expect_equal(r1$mode, "single")
  expect_equal(r1$across, "NPC")

  r0 <- NTSminer:::.nts_canopus_npc_columns(c("foo", "bar"))
  expect_equal(r0$mode, "none")

  # exact matching: a bare `npc` is checked before the `npc#...` columns and
  # `npc#pathway probability` must never be picked as a level
  rmix <- NTSminer:::.nts_canopus_npc_columns(
    c("NPC#pathway", "NPC#pathway Probability", "NPC#superclass", "NPC#class", "npc"))
  expect_equal(rmix$mode, "single")
  expect_equal(rmix$across, "npc")
  rprob <- NTSminer:::.nts_canopus_npc_columns(c("NPC#pathway Probability"))
  expect_true(is.na(rprob$pathway))
})

# ---------------------------------------------------------------------------
# NPC values (multi-column shape)
# ---------------------------------------------------------------------------

test_that("multi-column CANOPUS fills npc/superclass/class and leaves subclass NA", {
  dt <- read_sirius_dt("canopus_formula_summary_with_rows.tsv")
  cols <- NTSminer:::.nts_canopus_npc_columns(names(dt))
  expect_equal(cols$mode, "multi")

  dt <- NTSminer:::.nts_sirius_apply_canopus_npc(dt, cols)

  expect_identical(dt$canopus_npc[1],
                   "Alkaloids; Alkaloids; Alkaloids and derivatives")
  expect_identical(dt$canopus_superclass[1], "Alkaloids")
  expect_identical(dt$canopus_class[1], "Alkaloids and derivatives")

  expect_identical(dt$canopus_npc[2],
                   "Lipids and lipid-like molecules; Fatty Acyls; Fatty acids and conjugates")
  expect_identical(dt$canopus_superclass[2], "Fatty Acyls")
  expect_identical(dt$canopus_class[2], "Fatty acids and conjugates")

  # row 3 has empty NPC levels -> the joined string is NA
  expect_true(is.na(dt$canopus_npc[3]))

  # SIRIUS 6 emits no NPC#subclass -> the whole column is NA
  expect_length(dt$canopus_subclass, 3L)
  expect_true(all(is.na(dt$canopus_subclass)))
})

# ---------------------------------------------------------------------------
# Empty CANOPUS table
# ---------------------------------------------------------------------------

test_that("a header-only CANOPUS table is reported, not silently dropped", {
  tabs <- NULL
  expect_warning(
    tabs <- NTSminer:::.nts_sirius_read_canopus(
      sirius_fixture("canopus_formula_summary.tsv")),
    "no data rows")
  expect_length(tabs, 0L)

  # an empty CANOPUS table must not invent classifications
  empty_dt <- canopus_fixture("canopus_formula_summary.tsv")
  # a header-only TSV reads as all-logical columns; give the join keys the types
  # they would have if the table had data rows
  empty_dt[, molecular_formula := as.character(molecular_formula)]
  empty_dt[, adduct := as.character(adduct)]
  empty_dt[, rank := as.integer(rank)]
  out <- data.table::data.table(
    name = "FT0001", adduct = "[M - H]-",
    molecular_formula = "C7H7NO3", rank = 1L)
  merged <- suppressWarnings(
    NTSminer:::.nts_sirius_merge_canopus(out, empty_dt, "FT0001"))
  expect_equal(nrow(merged), 1L)
  for (cc in c("canopus_npc", "canopus_superclass", "canopus_class", "canopus_subclass")) {
    expect_true(cc %in% names(merged))
    expect_true(is.na(merged[[cc]][1]))
  }
})

# ---------------------------------------------------------------------------
# Score column mapping
# ---------------------------------------------------------------------------

test_that("explained_peaks maps to numExplainedPeaks, never to explainedIntensity", {
  hdr <- names(read_sirius_tsv("formula_identifications.tsv"))
  cm <- suppressWarnings(NTSminer:::.nts_sirius_col_map(hdr))
  expect_identical(cm$explained_peaks, "numExplainedPeaks")
  expect_identical(cm$sirius_score, "SiriusScore")
  expect_identical(cm$rank, "formulaRank")

  # the integer count (not the 0-1 fraction) is what ends up in the slot
  meta <- data.frame(
    numExplainedPeaks  = c(0, 3, 21),
    explainedIntensity = c(0.87, 1, 0.4),
    check.names = FALSE)
  got <- as.integer(meta[[cm$explained_peaks[1L]]])
  expect_identical(got, c(0L, 3L, 21L))
})

test_that("confidence is absent for SIRIUS 6 formula output and warns once", {
  hdr <- names(read_sirius_tsv("formula_identifications.tsv"))
  cm <- NULL
  expect_warning(cm <- NTSminer:::.nts_sirius_col_map(hdr), "confidence")
  expect_length(cm$confidence, 0L)

  # if a confidence score column does appear, prefer ConfidenceScoreApproximate
  cm2 <- NTSminer:::.nts_sirius_col_map(
    c("SiriusScore", "ConfidenceScoreExact", "ConfidenceScoreApproximate"))
  expect_identical(cm2$confidence, "ConfidenceScoreApproximate")
})

# ---------------------------------------------------------------------------
# CANOPUS merge: feature-keyed, no cross-talk, no row growth
# ---------------------------------------------------------------------------

test_that("CANOPUS merge keys on the feature, so same-formula features do not cross-talk", {
  dt <- canopus_fixture("canopus_formula_summary_with_rows.tsv")
  out <- data.table::data.table(
    name = c("FT0001", "FT0002"),
    adduct = c("[M - H]-", "[M - H]-"),
    molecular_formula = c("C7H7NO3", "C7H7NO3"),
    rank = c(1L, 1L))

  merged <- suppressWarnings(
    NTSminer:::.nts_sirius_merge_canopus(out, dt, c("FT0001", "FT0002")))

  expect_equal(nrow(merged), 2L)   # no row growth
  expect_identical(merged$canopus_superclass[merged$name == "FT0001"], "Alkaloids")
  expect_identical(merged$canopus_superclass[merged$name == "FT0002"], "Fatty Acyls")
  expect_identical(merged$canopus_class[merged$name == "FT0002"],
                   "Fatty acids and conjugates")
})

# ---------------------------------------------------------------------------
# set_sirius_fingerprints() accepts a plain data.frame (0-row path)
#
# Regression: the slot is declared `data.table` while the validator accepts a
# `data.frame`. The 0-row path of compute_sirius_fingerprints() used to pass a
# data.frame straight to the slot and die inside checkSlotAssignment() with
# "is(value, \"data.table\") is not TRUE" -- which masked the real warning
# ("SIRIUS may have produced no results") behind a confusing error.
# ---------------------------------------------------------------------------

test_that("set_sirius_fingerprints coerces a 0-row data.frame to data.table", {
  d <- tempfile("nts_sirius_empty")
  dir.create(d)
  p <- create_nts_project(d)

  empty_df <- data.frame(
    name = character(0), adduct = character(0),
    molecular_formula = character(0), rank = integer(0),
    sirius_score = numeric(0), zodiac_score = numeric(0),
    confidence = numeric(0), explained_peaks = integer(0),
    fingerprint = I(list()), fp_version = character(0),
    canopus_npc = character(0), canopus_superclass = character(0),
    canopus_class = character(0), canopus_subclass = character(0),
    stringsAsFactors = FALSE
  )

  expect_silent(p2 <- set_sirius_fingerprints(p, empty_df))
  sf <- methods::slot(p2, "sirius_fingerprints")
  expect_true(data.table::is.data.table(sf))
  expect_equal(nrow(sf), 0L)
  expect_setequal(names(sf), names(empty_df))

  # a non-empty data.frame is coerced too, and the rows survive
  df1 <- data.frame(
    name = "FT0001", adduct = "[M-H]-", molecular_formula = "C2F3O2",
    rank = 1L, sirius_score = 0.8, zodiac_score = 0.7, confidence = 0.5,
    explained_peaks = 12L, fp_version = "1",
    canopus_npc = NA_character_, canopus_superclass = NA_character_,
    canopus_class = NA_character_, canopus_subclass = NA_character_,
    stringsAsFactors = FALSE
  )
  df1$fingerprint <- I(list(c(0.1, 0.2)))
  p3 <- set_sirius_fingerprints(p, df1)
  expect_equal(nrow(methods::slot(p3, "sirius_fingerprints")), 1L)
  expect_identical(get_sirius_fingerprints(p3)$name, "FT0001")
})

test_that("a CANOPUS table without mappingFeatureId falls back with a warning", {
  dt <- canopus_fixture("canopus_formula_summary_with_rows.tsv")
  dt[, mappingfeatureid := NULL]
  out <- data.table::data.table(
    name = "FT0001", adduct = "[M - H]-",
    molecular_formula = "C7H7NO3", rank = 1L)

  expect_warning(
    NTSminer:::.nts_sirius_merge_canopus(out, dt, "FT0001"),
    "mappingFeatureId")
})

# ---------------------------------------------------------------------------
# batch_size validation
# ---------------------------------------------------------------------------

test_that(".nts_check_batch_size rejects invalid input and coerces valid input", {
  f <- NTSminer:::.nts_check_batch_size

  expect_error(f(0, 100))
  expect_error(f(-1, 100))
  expect_error(f(NA, 100))
  expect_error(f(Inf, 100))
  expect_error(f("20", 100))
  expect_error(f(c(1, 2), 100))

  expect_identical(f(NULL, 100), 100L)
  expect_identical(f(20.0, 100), 20L)
  expect_identical(f(20.9, 100), 20L)
})

# ---------------------------------------------------------------------------
# output_dir unlink protection
# ---------------------------------------------------------------------------

test_that(".nts_sirius_prepare_output_dir only deletes directories it created", {
  prep <- NTSminer:::.nts_sirius_prepare_output_dir

  # 1. does not exist -> created, marker written
  d1 <- tempfile("nts_sirius_")
  on.exit(unlink(d1, recursive = TRUE), add = TRUE)
  prep(d1)
  expect_true(dir.exists(d1))
  expect_true(file.exists(file.path(d1, ".nts_sirius_output")))

  # 2. exists WITH marker -> wiped and re-created (stale files gone)
  writeLines("stale", file.path(d1, "stale.txt"))
  suppressMessages(prep(d1))
  expect_true(dir.exists(d1))
  expect_false(file.exists(file.path(d1, "stale.txt")))
  expect_true(file.exists(file.path(d1, ".nts_sirius_output")))

  # 3. exists WITHOUT marker, overwrite = FALSE -> error, directory untouched
  d2 <- tempfile("nts_sirius_user_")
  on.exit(unlink(d2, recursive = TRUE), add = TRUE)
  dir.create(d2)
  writeLines("my data", file.path(d2, "keep.txt"))
  expect_error(prep(d2), "not created by NTSminer")
  expect_true(file.exists(file.path(d2, "keep.txt")))

  # 4. same directory, overwrite = TRUE -> warns and wipes
  expect_warning(
    suppressMessages(prep(d2, overwrite = TRUE)),
    "not created by NTSminer")
  expect_false(file.exists(file.path(d2, "keep.txt")))
  expect_true(file.exists(file.path(d2, ".nts_sirius_output")))
})

# ---------------------------------------------------------------------------
# CLI failure path (mocked - no process is spawned)
# ---------------------------------------------------------------------------

test_that(".nts_sirius_run_cli surfaces the exit status and stops/warns on failure", {
  run_cli <- NTSminer:::.nts_sirius_run_cli

  testthat::local_mocked_bindings(
    .nts_system2_capture = function(exe, args) list(status = 0L, output = "ok"),
    .package = "NTSminer"
  )
  r <- run_cli("exe", c("a"), quiet = TRUE)
  expect_equal(r$status, 0L)

  testthat::local_mocked_bindings(
    .nts_system2_capture = function(exe, args) list(status = 7L, output = "boom"),
    .package = "NTSminer"
  )
  # run_sirius() uses on_fail = "stop"; compute_sirius_fingerprints() uses "warn"
  expect_error(run_cli("exe", c("a"), quiet = TRUE, on_fail = "stop"), "exit status 7")
  expect_warning(run_cli("exe", c("a"), quiet = TRUE, on_fail = "warn"), "exit status 7")
  # the log file is written even on failure
  lf <- tempfile(fileext = ".log")
  on.exit(unlink(lf), add = TRUE)
  suppressWarnings(run_cli("exe", c("a"), log_file = lf, quiet = TRUE, on_fail = "none"))
  expect_true(file.exists(lf))
  expect_identical(readLines(lf), "boom")
})

# ---------------------------------------------------------------------------
# API fingerprint back-fill row selection
# ---------------------------------------------------------------------------

test_that(".nts_api_pick_fp_row maps the returned candidate to its exact row", {
  f <- NTSminer:::.nts_api_pick_fp_row
  fp <- data.table::data.table(
    name = c("FT0001", "FT0001", "FT0001"),
    rank = c(1L, 2L, 3L),
    molecular_formula = c("C7H7NO3", "C7H7NO3", "C7H7NO3"),
    adduct = c("[M - H]-", "[M - H]-", "[M - H]-"))
  fp[, fingerprint := vector("list", 3L)]   # all-empty fingerprint vectors

  # API returns rank 2 -> row 2 (not row 1) is selected
  cand <- f(fp, "FT0001",
            list(rank = 2L, molecularFormula = "C7H7NO3", adduct = "[M - H]-"))
  expect_identical(cand, 2L)

  # an already-filled row is never overwritten
  data.table::set(fp, 2L, "fingerprint", list(list(c(1, 2, 3))))
  expect_length(f(fp, "FT0001", list(rank = 2L)), 0L)
  expect_identical(f(fp, "FT0001", list(rank = 1L)), 1L)

  # an unknown feature name yields no candidate
  expect_length(f(fp, "FT9999", list(rank = 1L)), 0L)
})

# ---------------------------------------------------------------------------
# rbindlist(fill = TRUE) across differing batch column sets
# ---------------------------------------------------------------------------

test_that("rbindlist(fill = TRUE) tolerates differing column sets across batches", {
  a <- data.frame(name = "FT0001", SiriusScore = 1, stringsAsFactors = FALSE)
  b <- data.frame(name = "FT0002", SiriusScore = 2, extra = 9, stringsAsFactors = FALSE)
  res <- data.table::rbindlist(list(a, b), fill = TRUE)
  expect_equal(nrow(res), 2L)
  expect_true("extra" %in% names(res))
  expect_true(is.na(res$extra[1]))
  expect_equal(res$extra[2], 9)
})

# ---------------------------------------------------------------------------
# elements flag (item #1, sandbox-checkable parts only)
# ---------------------------------------------------------------------------

test_that(".nts_sirius_elements_flag keeps symbol case and warns when S is absent", {
  f <- NTSminer:::.nts_sirius_elements_flag

  # default value passes through verbatim (case-sensitive two-letter symbols)
  expect_identical(f("SClBrSeI"), c("--elements-considered", "SClBrSeI"))
  # an explicit token list is joined back into the same symbol string
  expect_identical(f("S,Cl,Br,Se,I"), c("--elements-considered", "SClBrSeI"))
  # C/H/N/O/P alone reduce to "no flag"
  expect_identical(f(c("C", "H", "N", "O", "P")), character(0))
  expect_identical(f(NULL), character(0))

  # missing S is warned about (PFAS sulfonates/sulfonamides)
  expect_warning(f("FClBr"), "does not include S")
  expect_identical(suppressWarnings(f("FClBr")),
                   c("--elements-considered", "FClBr"))
  expect_identical(suppressWarnings(f("CHNOPFClBr")),
                   c("--elements-considered", "CHNOPFClBr"))
})
