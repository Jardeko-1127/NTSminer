# Tests for assign_confidence_levels() (v3 — 10-level provenance scheme) -----
# The three library tiers are keyed on db_type ONLY:
#   in-house -> L1, online -> L2a, prediction -> L2b
# RT agreement splits only the network-propagation tier (L3a vs L3b).

test_that("empty db_annotations -> all L5", {
  d <- tempfile("nts_conf")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B", "C"),
      mz = c(100, 200, 300),
      rt = c(1, 2, 3),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_true("confidence_level" %in% names(ms1))
  expect_true("rt_confidence" %in% names(ms1))
  expect_equal(ms1$confidence_level, c("L5", "L5", "L5"))
  expect_equal(ms1$rt_confidence, rep("unavailable", 3L))
})


test_that("in-house library always -> L1 (regardless of RT)", {
  d <- tempfile("nts_conf_l1")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      ref_id = c("STD_001", "STD_002"),
      Library_ID = c("NTS_InHouse_PFAS", "NTS_InHouse_PFAS"),
      db_type = c("inhouse", "inhouse"),
      lib_rt = c(1.2, NA_real_),
      rt_error_s = c(0.5, NA_real_),
      score_total = c(0.95, 0.95),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  # Both A (RT present) and B (no RT) are in-house -> L1
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L1")
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L1")
})


test_that("online library always -> L2a (regardless of RT)", {
  d <- tempfile("nts_conf_l2a")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      ref_id = c("MoNA_001", "MoNA_002"),
      Library_ID = c("MoNA_PFAS", "MoNA_PFAS"),
      db_type = c("online", "online"),
      lib_rt = c(1.2, NA_real_),
      rt_error_s = c(0.5, NA_real_),
      score_total = c(0.91, 0.91),
      score_cosine_ms2 = c(0.90, 0.90),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  # Both A (RT present) and B (no RT) are online -> L2a
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L2a")
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L2a")
})


test_that("prediction library always -> L2b (regardless of RT)", {
  d <- tempfile("nts_conf_l2b")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      ref_id = c("Pred_001", "Pred_002"),
      Library_ID = c("Pred_CFMID", "Pred_EPISUITE"),
      db_type = c("prediction", "prediction"),
      lib_rt = c(1.2, 2.1),
      rt_error_s = c(0.5, 400),   # A agrees; B is far beyond rt_tol_sec
      score_total = c(0.8, 0.8),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p, rt_tol_sec = 30)

  ms1 <- get_ms1_annotations(p)
  # Both A and B are prediction matches -> L2b (RT no longer splits this tier)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L2b")
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L2b")
})


test_that("L3a/L3b: network propagation +/- RT", {
  d <- tempfile("nts_conf_l3ab")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B", "C"),
      mz = c(100, 200, 300),
      rt = c(1, 2, 3),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B", "C"),
      ref_id = c("NET_A", "NET_B", "NET_C"),
      annotation_method = c("Dual_Network", "Propagation", "Dual_Network"),
      propagation_level = c(1L, 3L, 5L),
      db_type = c("online", "online", "online"),
      lib_rt = c(1.1, NA_real_, 3.1),
      rt_error_s = c(0.3, NA_real_, 0.4),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p, propagation_max_level = 3L)

  ms1 <- get_ms1_annotations(p)
  # A: propagation + RT -> L3a
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L3a")
  # B: propagation - RT -> L3b
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L3b")
  # C: propagation_level > max -> L5
  expect_equal(ms1$confidence_level[ms1$name == "C"], "L5")
})


test_that("L3c: homologue series evidence", {
  d <- tempfile("nts_conf_l3c")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      homologous_series_id = c("S1", NA_character_),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L3c")
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L5")
})


test_that("L3c from db annotation_method", {
  d <- tempfile("nts_conf_l3c_db")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A"),
      mz = c(100),
      rt = c(1),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = "HM_A",
      annotation_method = "Homologue",
      db_type = "online",
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L3c")
})


test_that("L3d: transformation evidence", {
  d <- tempfile("nts_conf_l3d")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      transformed_formula = c("C8HF15O2", NA_character_),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L3d")
})


test_that("L3e: multi-method consensus", {
  d <- tempfile("nts_conf_l3e")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A"),
      mz = c(100),
      rt = c(1),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = "MM_A",
      # Two structural methods agreeing (no spectral component) -> L3e. A row
      # that also carried a Spectral_DB hit would instead be a library match.
      annotation_method = "Network_DF; Network_MD",
      db_type = "online",
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L3e")
})


test_that("L4: isotope_cosine >= 0.9 in db_annotations", {
  d <- tempfile("nts_conf_l4_iso")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      ref_id = c("DB_A", "DB_B"),
      score_total = c(0.6, 0.6),
      isotope_cosine = c(0.95, 0.5),
      # Structural (non-spectral) rows, so the isotope evidence -- not a library
      # tier -- decides the level.
      annotation_method = c("Network_MD", "Network_MD"),
      db_type = c("online", "online"),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L4")
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L5")
})


test_that("L4 legacy fallback: isotope_kind non-NA", {
  d <- tempfile("nts_conf_l4_legacy")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      isotope_kind = c("13C", NA_character_),
      isotope_c13_n = c(2L, NA_integer_),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L4")
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L5")
})


test_that("priority: L1 beats L2a", {
  d <- tempfile("nts_conf_prio_l1")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A"),
      mz = c(100),
      rt = c(1),
      stringsAsFactors = FALSE
    )
  )
  # Both inhouse and online matches exist for A
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "A"),
      ref_id = c("STD_001", "MoNA_001"),
      Library_ID = c("NTS_InHouse_PFAS", "MoNA_PFAS"),
      db_type = c("inhouse", "online"),
      lib_rt = c(1.2, 1.3),
      rt_error_s = c(0.5, 0.5),
      score_total = c(0.95, 0.91),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  # L1 (in-house) beats L2a (online)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L1")
})


test_that("priority: L2a beats L3c", {
  d <- tempfile("nts_conf_prio_l2a_l3c")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A"),
      mz = c(100),
      rt = c(1),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = "MoNA_001",
      Library_ID = "MoNA_PFAS",
      db_type = "online",
      lib_rt = 1.2,
      rt_error_s = 0.5,
      score_total = 0.91,
      annotation_method = "Spectral_DB; Homologue",
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = "A",
      mz = 100,
      rt = 1,
      homologous_series_id = "S1",
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  # L2a beats L3c and L3e
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L2a")
})


test_that("priority: L3c beats L4 (homologue > isotope)", {
  d <- tempfile("nts_conf_prio_l3c_l4")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A"),
      mz = c(100),
      rt = c(1),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = "HM_A",
      annotation_method = "Homologue",
      isotope_cosine = 0.95,
      db_type = "online",
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  # L3c (priority 6) beats L4 (priority 9)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L3c")
})


test_that("priority: L3c > L3d > L3e", {
  d <- tempfile("nts_conf_prio_l3_order")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A"),
      mz = c(100),
      rt = c(1),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = "MULTI_A",
      annotation_method = "Homologue; Transformation",
      db_type = "online",
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  # L3c (Homologue) beats L3d (Transformation) and L3e (multi-method)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L3c")
})


test_that("auto-detect db_type from Library_ID pattern", {
  d <- tempfile("nts_conf_autodetect")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B", "C"),
      mz = c(100, 200, 300),
      rt = c(1, 2, 3),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B", "C"),
      ref_id = c("IH_001", "MoNA_001", "Pred_001"),
      Library_ID = c("NTS_InHouse_v2", "MoNA_PFAS_v1", "Pred_EPISUITE"),
      lib_rt = c(1.2, 2.1, 3.0),
      rt_error_s = c(0.5, 0.3, 0.4),
      score_total = c(0.95, 0.85, 0.75),
      stringsAsFactors = FALSE
    )
  )

  # No db_type column here on purpose: this test exercises the ID-pattern
  # fallback. That path emits a one-shot warning (covered in test-db-type.R),
  # silenced here to keep the suite output clean.
  p <- suppressWarnings(assign_confidence_levels(p))

  ms1 <- get_ms1_annotations(p)
  # A: Library_ID contains "InHouse" -> in-house -> L1
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L1")
  # B: Library_ID contains "MoNA" -> online -> L2a
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L2a")
  # C: Library_ID contains "Pred" -> prediction -> L2b
  expect_equal(ms1$confidence_level[ms1$name == "C"], "L2b")
})


test_that("db_type column takes priority over Library_ID pattern", {
  d <- tempfile("nts_conf_dbtype")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A"),
      mz = c(100),
      rt = c(1),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = "ANY_001",
      Library_ID = "MoNA_PFAS",
      db_type = "inhouse",
      lib_rt = 1.2,
      rt_error_s = 0.5,
      score_total = 0.95,
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  # db_type="inhouse" -> L1, even though Library_ID looks like MoNA
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L1")
})


test_that("library tiers ignore RT agreement (provenance-only)", {
  # RT agreement no longer gates the library tiers: an in-house match with a
  # far-off RT is still L1, and an online match with a far-off RT is still L2a.
  d <- tempfile("nts_conf_rt_tol")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B", "C", "D"),
      mz = c(100, 200, 300, 400),
      rt = c(1, 2, 3, 4),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B", "C", "D"),
      ref_id = c("STD_001", "STD_002", "MoNA_001", "MoNA_002"),
      Library_ID = c("NTS_InHouse", "NTS_InHouse", "MoNA_PFAS", "MoNA_PFAS"),
      db_type = c("inhouse", "inhouse", "online", "online"),
      lib_rt = c(1.2, 2.1, 3.2, 4.1),
      rt_error_s = c(5, 400, 5, 400),   # B and D are 400 s off
      score_total = c(0.95, 0.95, 0.90, 0.90),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p, rt_tol_sec = 30)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L1")  # in-house, RT agrees
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L1")  # in-house, RT far off -> still L1
  expect_equal(ms1$confidence_level[ms1$name == "C"], "L2a") # online, RT agrees
  expect_equal(ms1$confidence_level[ms1$name == "D"], "L2a") # online, RT far off -> still L2a
})


test_that("db row with NA fields still gets its provenance tier (online -> L2a)", {
  d <- tempfile("nts_conf_na")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = "A",
      mz = 100,
      rt = 1,
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = NA_character_,
      score_total = NA_real_,
      annotation_method = NA_character_,
      db_type = "online",
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  # The presence of a db row with online provenance is enough for L2a; only a
  # feature with NO db row (and no other evidence) falls to L5.
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L2a")
})


test_that("empty ms1_annotations -> warning", {
  d <- tempfile("nts_conf_empty")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)

  expect_warning(
    p <- assign_confidence_levels(p),
    "Empty ms1_annotations"
  )
})


test_that("non-NTSProject input -> error", {
  expect_error(
    assign_confidence_levels(list()),
    "must be an NTSProject"
  )
})


test_that("rt_confidence always unavailable", {
  d <- tempfile("nts_conf_rt")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = "STD_A",
      Library_ID = "NTS_InHouse",
      db_type = "inhouse",
      lib_rt = 1.2,
      rt_error_s = 0.5,
      score_total = 0.95,
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$rt_confidence, rep("unavailable", 2L))
})


test_that("feature_id column in db_annotations works", {
  d <- tempfile("nts_conf_featid")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A"),
      mz = c(100),
      rt = c(1),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      feature_id = "A",
      ref_id = "STD_A",
      db_type = "inhouse",
      lib_rt = 1.2,
      rt_error_s = 0.5,
      score_total = 0.95,
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L1")
})


test_that("kmd_series_id -> L3c", {
  d <- tempfile("nts_conf_kmd")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_ms1_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      kmd_series_id = c("CF2_001", NA_character_),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p)

  ms1 <- get_ms1_annotations(p)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L3c")
})


test_that("online tier is provenance-only: low MS2 cosine still L2a", {
  # The old cosine gate (>= ms2_cosine_threshold) was removed with the switch
  # to provenance-only library tiers: an online hit is L2a regardless of how
  # low its MS2 cosine is.
  d <- tempfile("nts_conf_cosine_gate")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_feature_table(
    p,
    data.frame(
      name = c("A", "B"),
      mz = c(100, 200),
      rt = c(1, 2),
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      ref_id = c("MoNA_001", "MoNA_002"),
      db_type = c("online", "online"),
      score_total = c(0.95, 0.95),
      score_cosine_ms2 = c(0.40, 0.90),
      stringsAsFactors = FALSE
    )
  )

  p <- assign_confidence_levels(p, ms2_cosine_threshold = 0.7)

  ms1 <- get_ms1_annotations(p)
  # A: low cosine 0.40 still L2a (gate removed)
  expect_equal(ms1$confidence_level[ms1$name == "A"], "L2a")
  # B: cosine 0.90 also L2a
  expect_equal(ms1$confidence_level[ms1$name == "B"], "L2a")
})
