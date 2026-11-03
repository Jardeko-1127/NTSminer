## @ms1_annotations is a named list of per-step sub-tables (see NEWS 0.1.15).
## These tests pin the slot shape, the flat-view / per-step accessors, the
## lossless round-trip, and backward compatibility with the legacy flat slot.

make_proj <- function() {
  d <- tempfile("nts_ms1list")
  dir.create(d)
  p <- create_nts_project(d)
  ft <- data.frame(
    name = c("F1", "F2", "F3"),
    mz = c(400.1, 414.1, 428.1),
    rt = c(60, 70, 80),
    s1 = c(1000.4, NA, 2000.6),
    s2 = c(1500.5, 2000.5, NA),
    stringsAsFactors = FALSE
  )
  set_feature_table(p, ft)
}

make_flat <- function() {
  data.frame(
    name = c("F1", "F2", "F3"),
    mz = c(400.1, 414.1, 428.1),
    rt = c(60, 70, 80),
    is_isf = c(TRUE, FALSE, NA),
    parent_mz = c(NA_real_, 380.0, NA_real_),
    kmd_series_id = c("S1", "S1", NA_character_),
    KMD = c(-0.1, -0.2, NA_real_),
    my_custom_col = c("a", "b", "c"),
    stringsAsFactors = FALSE
  )
}

test_that("the slot is a named list keyed by processing step", {
  p <- set_ms1_annotations(make_proj(), make_flat())
  sl <- methods::slot(p, "ms1_annotations")
  expect_true(is.list(sl))
  expect_false(is.data.frame(sl))
  expect_true(all(names(sl) %in% names(NTSminer:::.nts_ms1_groups())))
  ## only steps that actually carry data are stored
  expect_setequal(names(sl), c("base", "isotopes", "isf", "kmd", "extra"))
  expect_false(any(c("eic", "library", "adduct", "master", "network",
                     "annotation") %in% names(sl)))
})

test_that("steps whose columns are absent or all-NA are pruned from the slot", {
  ## `adduct` / `network` / `annotation` columns exist (created all-NA by the
  ## canonical schema) but carry no data here, so they must not be stored.
  p <- set_ms1_annotations(make_proj(), make_flat())
  sl <- methods::slot(p, "ms1_annotations")
  expect_true("isf" %in% names(sl))
  expect_false(any(c("adduct", "master", "network", "annotation") %in% names(sl)))
  ## a pruned step is still a valid `step`, just empty
  expect_equal(nrow(get_ms1_annotations(p, "adduct")), 0L)
  ## populating the step makes it reappear in the slot
  p2 <- set_ms1_annotations(
    p,
    data.frame(name = c("F1", "F2", "F3"), adduct_type = c("[M-H]-", NA, NA)),
    step = "adduct"
  )
  expect_true("adduct" %in% names(methods::slot(p2, "ms1_annotations")))
  expect_identical(get_ms1_annotations(p2)$adduct_type, c("[M-H]-", NA, NA))
  ## the merged flat view still exposes the full canonical column set
  fv <- get_ms1_annotations(p2)
  expect_true(all(c("adduct_type", "neutral_mass", "inferred_formula",
                    "inferred_name", "annotation_level") %in% names(fv)))
})

test_that("each step holds only its own columns (plus the name key)", {
  p <- set_ms1_annotations(make_proj(), make_flat())
  flat <- get_ms1_annotations(p)

  isf <- get_ms1_annotations(p, "isf")
  expect_true(is.data.frame(isf))
  expect_setequal(names(isf), c("name", intersect(NTSminer:::.nts_ms1_groups()$isf, names(flat))))
  expect_false("KMD" %in% names(isf))
  expect_false("is_isotope" %in% names(isf))

  kmd <- get_ms1_annotations(p, "kmd")
  expect_setequal(names(kmd), c("name", intersect(NTSminer:::.nts_ms1_groups()$kmd, names(flat))))

  ## a column that belongs to no step is parked in `extra`, never dropped
  expect_setequal(names(get_ms1_annotations(p, "extra")), c("name", "my_custom_col"))
  ## a step whose columns were never created is an empty table
  expect_equal(nrow(get_ms1_annotations(p, "library")), 0L)
  expect_error(get_ms1_annotations(p, "nope"), "Unknown `ms1_annotations` step")
  expect_error(set_ms1_annotations(p, data.frame(name = "F1"), step = "nope"),
               "Unknown `ms1_annotations` step")
})

test_that("the flat view is lossless and keeps the legacy column set", {
  p <- set_ms1_annotations(make_proj(), make_flat())
  fv <- get_ms1_annotations(p)
  expect_true(is.data.frame(fv))
  expect_identical(as.character(fv$name), c("F1", "F2", "F3"))
  expect_identical(fv$my_custom_col, c("a", "b", "c"))
  expect_identical(fv$KMD, c(-0.1, -0.2, NA_real_))
  expect_true(is.na(fv$is_isf[3]))
  expect_identical(names(fv)[1:7],
                   c("name", "mzmin", "mzmax", "rtmin", "rtmax", "mz", "rt"))
  ## the eager head/tail columns of the legacy flat schema are still created
  expect_true(all(c("is_isotope", "isotope_shift_da", "isf_annotation",
                    "standard_md", "KM", "NKM", "mz_remainder", "adduct_type",
                    "neutral_mass", "inferred_formula", "annotation_level",
                    "my_custom_col") %in% names(fv)))

  ## flat -> list -> flat is value-lossless
  p2 <- set_ms1_annotations(p, fv)
  fv2 <- get_ms1_annotations(p2)
  expect_identical(dim(fv2), dim(fv))
  expect_identical(as.character(fv2$name), as.character(fv$name))
  for (cl in names(fv)) expect_equal(fv2[[cl]], fv[[cl]])
})

test_that("set(step = ) replaces only that step", {
  p <- set_ms1_annotations(make_proj(), make_flat())
  p2 <- set_ms1_annotations(
    p, data.frame(name = c("F1", "F2", "F3"), KMD = c(9, 8, 7)), step = "kmd"
  )
  f2 <- get_ms1_annotations(p2)
  expect_identical(f2$KMD, c(9, 8, 7))
  expect_identical(f2$kmd_series_id, c("S1", "S1", NA_character_))
  expect_identical(f2$my_custom_col, c("a", "b", "c"))
  expect_identical(get_ms1_annotations(p2, "isf")$is_isf, c(TRUE, FALSE, NA))
})

test_that("forbidden physchem / MS2 columns are still stripped on flat writes", {
  p <- set_ms1_annotations(make_proj(), make_flat())
  p3 <- set_ms1_annotations(
    p, transform(get_ms1_annotations(p), SMILES = "CCO", formula = "C2H6")
  )
  expect_false(any(c("SMILES", "formula") %in% names(get_ms1_annotations(p3))))
})

test_that("an empty project keeps an empty list slot", {
  d <- tempfile("nts_ms1empty")
  dir.create(d)
  pe <- create_nts_project(d)
  expect_identical(methods::slot(pe, "ms1_annotations"), list())
  expect_true(is.data.frame(get_ms1_annotations(pe)))
  expect_equal(nrow(get_ms1_annotations(pe)), 0L)
  expect_identical(methods::slot(set_ms1_annotations(pe, data.frame()), "ms1_annotations"),
                   list())
  expect_null(get_ms1_annotations(pe, "kmd"))
})

test_that(".nts_init_ms1_from_ft still returns the flat view, stored as a list", {
  ft <- data.frame(name = "F1", mz = 400, rt = 60, s1 = 10, stringsAsFactors = FALSE)
  init <- NTSminer:::.nts_init_ms1_from_ft(ft)
  expect_true(is.data.frame(init))
  expect_true(all(c("name", "is_isotope", "is_isf", "KMD") %in% names(init)))
  p <- set_ms1_annotations(make_proj(), init)
  sl <- methods::slot(p, "ms1_annotations")
  expect_true(is.list(sl) && !is.data.frame(sl))
})

test_that("copy_project deep-copies the sub-tables", {
  p <- set_ms1_annotations(make_proj(), make_flat())
  cp <- copy_project(p)
  expect_true(is.list(methods::slot(cp, "ms1_annotations")))
  methods::slot(cp, "ms1_annotations")[["isf"]][["is_isf"]][1] <- FALSE
  expect_true(get_ms1_annotations(p, "isf")$is_isf[1])
})

test_that("a legacy flat data.frame pushed straight into the slot still reads", {
  p <- set_ms1_annotations(make_proj(), make_flat())
  flat <- make_flat()
  ## S4 coerces a data.frame to a bare list of same-length columns
  methods::slot(p, "ms1_annotations") <- flat
  expect_true(is.data.frame(get_ms1_annotations(p)))
  expect_equal(nrow(get_ms1_annotations(p)), 3L)
  expect_identical(get_ms1_annotations(p)$my_custom_col, c("a", "b", "c"))
  expect_setequal(names(get_ms1_annotations(p, "isf")), c("name", "is_isf", "parent_mz"))
})

test_that("feature_table sample columns are stored as integers (NA preserved)", {
  p <- make_proj()
  ft <- get_feature_table(p)
  expect_true(is.integer(ft$s1))
  expect_true(is.integer(ft$s2))
  expect_true(is.na(ft$s1[2]))
  expect_identical(ft$s1[1], 1000L)
  expect_false(is.integer(ft$mz))
})

test_that("get_sirius_fingerprints(as_matrix = TRUE) resolves mz via the accessor", {
  p <- set_ms1_annotations(make_proj(), make_flat())
  fp <- data.table::data.table(
    name = "F1", adduct = "[M-H]-", molecular_formula = "C2H6O", rank = 1L,
    sirius_score = 0.9, zodiac_score = 0.9, confidence = "high",
    explained_peaks = 10L, fingerprint = list(c(0.1, 0.2)), fp_version = "1",
    canopus_npc = NA_character_, canopus_superclass = NA_character_,
    canopus_class = NA_character_, canopus_subclass = NA_character_
  )
  m <- get_sirius_fingerprints(set_sirius_fingerprints(p, fp), as_matrix = TRUE)
  expect_true(is.matrix(m))
  expect_match(rownames(m)[1], "mz=400.1")
})
