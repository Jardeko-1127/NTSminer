# Tests for the library-free MS2 CF2 fragment-ladder / series detector.
# Design: docs/analysis/NTSminer_MS2_CF2碎片阶梯检测设计方案_v0.1.22.md

CF3 <- 68.9957596
C2F5 <- 118.992566
C3F7 <- 168.9893724
C4F9 <- 218.9861788
C5F11 <- 268.9829852
CF2 <- 49.9968064

# ---------------------------------------------------------------------------
# S1: pure function
# ---------------------------------------------------------------------------

test_that("find_pfas_fragment_series: clean 3-member CnF(2n+1)- ladder", {
  sp <- data.frame(mz = c(CF3, C2F5, C3F7), intensity = c(100, 80, 60))
  r <- find_pfas_fragment_series(sp)

  expect_type(r, "list")
  expect_equal(r$ladder$len, 3L)
  expect_equal(r$ladder$span, 2L)
  expect_equal(r$ladder$gaps, c(1L, 1L))
  expect_equal(r$series$max_n, 3L)
  expect_equal(r$series$hits, c(1L, 2L, 3L))
  expect_true(r$is_hit)
  expect_equal(unname(r$chain[["from_series"]]), 3L)
  expect_equal(unname(r$chain[["from_ladder"]]), 3L)
  expect_gt(r$score, 0)
})

test_that("find_pfas_fragment_series: missing rungs (decision #4, loose vs strict)", {
  sp <- data.frame(mz = c(CF3, C4F9), intensity = c(100, 90))

  rl <- find_pfas_fragment_series(sp, mode = "loose")
  expect_equal(rl$ladder$len, 2L)
  expect_equal(rl$ladder$span, 3L)
  expect_equal(rl$ladder$gaps, 3L)
  expect_lt(rl$ladder$len, 3L) # below min_len -> not a qualifying ladder

  rs <- find_pfas_fragment_series(sp, mode = "strict")
  expect_equal(rs$ladder$len, 0L)
  expect_equal(rs$ladder$span, 0L)
})

test_that("find_pfas_fragment_series: boundaries never error", {
  expect_equal(find_pfas_fragment_series(data.frame(mz = 100, intensity = 1))$ladder$len, 0L)
  expect_equal(
    find_pfas_fragment_series(data.frame(mz = numeric(0), intensity = numeric(0)))$ladder$len,
    0L
  )
  expect_equal(find_pfas_fragment_series(NULL)$ladder$len, 0L)
  # tolerance outside: no edge
  sp <- data.frame(mz = c(CF3, CF3 + CF2 + 0.5), intensity = c(1, 1))
  expect_equal(find_pfas_fragment_series(sp)$ladder$len, 0L)
  # input shapes
  expect_equal(
    find_pfas_fragment_series(list(mz = c(CF3, C2F5, C3F7), intensity = c(1, 1, 1)))$ladder$len,
    3L
  )
  expect_equal(
    find_pfas_fragment_series(cbind(c(CF3, C2F5, C3F7), c(1, 1, 1)))$ladder$len,
    3L
  )
  # unsorted input is sorted internally
  expect_equal(
    find_pfas_fragment_series(data.frame(mz = c(C3F7, CF3, C2F5), intensity = c(60, 100, 80)))$ladder$len,
    3L
  )
})

test_that("find_pfas_fragment_series: relative-intensity floor shrinks the ladder", {
  sp <- data.frame(mz = c(CF3, C2F5, C3F7), intensity = c(1000, 1e-6, 1000))
  r <- find_pfas_fragment_series(sp, min_rel_int = 0.01)
  expect_equal(r$ladder$len, 2L) # the weak C2F5- is dropped
})

test_that("find_pfas_fragment_series: coincidental ~50 Da spacing is not a PFAS hit", {
  sp <- data.frame(mz = c(100, 150, 200), intensity = c(1, 1, 1))
  r <- find_pfas_fragment_series(sp)
  expect_equal(r$ladder$len, 3L)      # naive ladder does form
  expect_equal(r$series$n_hits, 0L)   # but no CnF(2n+1)- series
  expect_lt(r$score, 0.5)
  expect_false(r$is_hit)              # joint-evidence gate (design §8.7)
})

test_that("find_pfas_fragment_series: polarity gates the anion series", {
  sp <- data.frame(mz = c(CF3, C2F5, C3F7), intensity = c(1, 1, 1))
  expect_equal(find_pfas_fragment_series(sp, polarity = "pos")$series$n_hits, 0L)
  expect_equal(find_pfas_fragment_series(sp, polarity = "neg")$series$n_hits, 3L)
})

test_that("find_pfas_fragment_series: longer ladder does not score lower", {
  short <- find_pfas_fragment_series(
    data.frame(mz = c(CF3, C2F5, C3F7), intensity = c(1, 1, 1))
  )$score
  long <- find_pfas_fragment_series(
    data.frame(mz = c(CF3, C2F5, C3F7, C4F9, C5F11, C5F11 + CF2), intensity = rep(1, 6))
  )$score
  expect_gte(long, short)
})

test_that("pfas_diagnostic_rules: schema and semantics", {
  rr <- pfas_diagnostic_rules()
  expect_s3_class(rr, "data.frame")
  expect_true(all(
    c("subclass", "ion_formula", "ion_mz", "polarity", "type", "required", "label") %in%
      names(rr)
  ))
  expect_true(all(c("PFCA", "PFSA", "FTS") %in% rr$subclass))
  # exactly one ladder-requirement row, belonging to FTS
  expect_equal(sum(rr$type == "series"), 1L)
  expect_equal(rr$subclass[rr$type == "series"], "FTS")
  # sulfonate head-group mass present
  expect_true(any(abs(rr$ion_mz - 79.957367) < 1e-9, na.rm = TRUE))
  # no positive-mode preset rules (PFAS sub-classes are negative mode)
  expect_equal(nrow(pfas_diagnostic_rules("pos")), 0L)
})

test_that("sub-class hints follow the rule semantics", {
  # PFCA: any CnF(2n+1)- member
  expect_true("PFCA" %in% find_pfas_fragment_series(
    data.frame(mz = c(CF3, C2F5, C3F7), intensity = c(1, 1, 1))
  )$diagnostic$subclass)

  # PFSA: sulfonate head group, no ladder required
  sa <- data.frame(mz = c(79.957367, 98.9557702), intensity = c(1, 1))
  expect_true("PFSA" %in% find_pfas_fragment_series(sa)$diagnostic$subclass)

  # FTS needs head group AND a CF2 ladder; head group alone must not fire
  expect_false("FTS" %in% find_pfas_fragment_series(sa)$diagnostic$subclass)
  fts <- data.frame(mz = c(CF3, C2F5, C3F7, 79.957367), intensity = rep(1, 4))
  expect_true("FTS" %in% find_pfas_fragment_series(fts)$diagnostic$subclass)
})

test_that("diagnostic expert table is matched when supplied explicitly", {
  diag_tbl <- data.frame(
    fragment_formula = c("CF3- ", "SO3- "),
    fragment_mass = c(CF3, 79.957367),
    polarity = c("neg", "neg"),
    stringsAsFactors = FALSE
  )
  r <- find_pfas_fragment_series(
    data.frame(mz = c(CF3, C2F5, C3F7), intensity = c(1, 1, 1)),
    diag_csv = diag_tbl
  )
  expect_true(any(trimws(r$diagnostic$hits) == "CF3-"))
})

test_that("the packaged diagnostic CSV resolves and matches", {
  csv <- default_diagnostic_fragments_csv()
  skip_if(is.na(csv) || !nzchar(csv), "packaged diagnostic table unavailable")
  r <- find_pfas_fragment_series(
    data.frame(mz = c(CF3, C2F5, C3F7), intensity = c(1, 1, 1)),
    diag_csv = csv
  )
  expect_true(any(trimws(r$diagnostic$hits) == "CF3-"))
  expect_true(any(trimws(r$diagnostic$hits) == "C3F7-"))
  expect_true("PFCA" %in% r$diagnostic$subclass)
})

# ---------------------------------------------------------------------------
# S2: project-level integration
# ---------------------------------------------------------------------------

.nts_make_series_project <- function() {
  td <- tempfile("nts_pfas_series")
  dir.create(td)
  p <- create_nts_project(td)
  p <- set_feature_table(
    p,
    data.frame(name = c("F1", "F2"), mz = c(300, 400), rt = c(1, 2),
               stringsAsFactors = FALSE)
  )
  p <- set_ms1_annotations(
    p,
    data.frame(name = c("F1", "F2"), mz = c(300, 400), rt = c(1, 2),
               annotation_level = c("Isolated Unknown", "Isolated Unknown"),
               stringsAsFactors = FALSE)
  )
  p <- set_ms2_annotations(
    p,
    data.frame(name = c("F1", "F2"),
               diagnostic_fragments = c("SO3-", NA_character_),
               stringsAsFactors = FALSE)
  )
  p@clean_ms2 <- list(
    F1 = list(
      feature_id = "F1", precursorMZ = 300, rt = 1, n_raw_scans = 3L,
      spectrum = data.frame(
        mz = c(CF3, C2F5, C3F7), intensity = c(100, 80, 60)
      )
    ),
    F2 = list(
      feature_id = "F2", precursorMZ = 400, rt = 2, n_raw_scans = 3L,
      spectrum = data.frame(mz = c(123.4, 187.7, 250.2), intensity = c(10, 9, 8))
    )
  )
  p
}

test_that("detect_pfas_fragment_series writes the new columns and tags hits", {
  skip_if_not_installed("data.table")
  p <- .nts_make_series_project()

  p2 <- detect_pfas_fragment_series(p, quiet = TRUE)
  m2 <- get_ms2_annotations(p2)

  expect_true(all(.nts_pfas_series_columns() %in% names(m2)))

  i1 <- match("F1", m2$name)
  i2 <- match("F2", m2$name)
  expect_equal(m2$pfas_ladder_len[i1], 3L)
  expect_equal(m2$pfas_series_max_n[i1], 3L)
  expect_true(grepl("PFCA", m2$pfas_subclass_hint[i1], fixed = TRUE))
  expect_gt(m2$pfas_series_score[i1], 0)

  expect_equal(m2$pfas_ladder_len[i2], 0L)
  expect_true(is.na(m2$pfas_subclass_hint[i2]) || !nzchar(m2$pfas_subclass_hint[i2]))

  # existing columns are preserved untouched
  expect_equal(m2$diagnostic_fragments[i1], "SO3-")
})

test_that("detect_pfas_fragment_series tags annotation_level only for hits", {
  skip_if_not_installed("data.table")
  p <- .nts_make_series_project()
  p2 <- detect_pfas_fragment_series(p, quiet = TRUE)

  ms1 <- get_ms1_annotations(p2)
  lvl <- ms1$annotation_level[match(c("F1", "F2"), ms1$name)]
  expect_true(grepl("[MS2 Verified]", lvl[1L], fixed = TRUE))
  expect_false(grepl("[MS2 Verified]", lvl[2L], fixed = TRUE))
})

test_that("detect_pfas_fragment_series appends a log step", {
  skip_if_not_installed("data.table")
  p <- .nts_make_series_project()
  p2 <- detect_pfas_fragment_series(p, quiet = TRUE)
  steps <- vapply(p2@log, function(z) as.character(z$step)[1L], character(1L))
  expect_true("detect_pfas_fragment_series" %in% steps)
})

test_that("detect_pfas_fragment_series is a safe no-op without clean_ms2", {
  td <- tempfile("nts_pfas_empty")
  dir.create(td)
  p <- create_nts_project(td)
  p2 <- detect_pfas_fragment_series(p, quiet = TRUE)
  expect_s4_class(p2, "NTSProject")
})

# ---------------------------------------------------------------------------
# S4: structural confidence guardrail keeps PFAS_Series rows structural
# ---------------------------------------------------------------------------

test_that("PFAS_Series is treated as a structural annotation_method token", {
  skip_if_not_installed("data.table")
  # A db row carrying the PFAS_Series token must not receive a library tier
  # even when it also carries a db_type: PFAS_Series is a structural marker.
  db <- data.frame(
    name = c("A", "B"),
    annotation_method = c("PFAS_Series", "Spectral_DB"),
    db_type = c("inhouse", "inhouse"),
    stringsAsFactors = FALSE
  )
  ft <- data.frame(
    name = c("A", "B"), mz = c(1, 2), rt = c(1, 2),
    annotation_level = c("x", "y"),
    stringsAsFactors = FALSE
  )
  res <- .nts_confidence_build_feature_table(
    ft, db, ms2_cosine_threshold = 0.5, propagation_max_level = 3L
  )
  # A (structural) keeps the default L5 -> priority 10; B (spectral) gets L1.
  pri <- res$priority[match(c("A", "B"), res$name)]
  expect_equal(pri[2L], 1L)
  expect_gt(pri[1L], 1L)
})
