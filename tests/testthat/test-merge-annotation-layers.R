test_that("merge_annotation_layers combines sources with semicolon collapse", {
  d <- tempfile("nts_merge_ann")
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
      name = "A",
      ref_id = "DB_A",
      Compound.name = "Compound A",
      score_total = 0.91,
      stringsAsFactors = FALSE
    )
  )

  net_frag <- list(
    nodes = data.frame(
      id = c("A", "B"),
      ref_id = c(NA_character_, "NET_B"),
      propagated_ID = c(NA_character_, "NET_B"),
      Compound.name = c(NA_character_, "Compound B"),
      Propagation_Level = c(NA_integer_, 1L),
      propagated_similarity = c(NA_real_, 0.55),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = "A", to = "B", stringsAsFactors = FALSE)
  )

  out <- merge_annotation_layers(
    p,
    net_frag_ann = net_frag,
    include_all_features = TRUE
  )

  expect_true(all(c("feature", "annotation_method", "ref_id") %in% names(out)))
  expect_equal(nrow(out), 3L)

  row_a <- out[out$feature == "A", , drop = FALSE]
  expect_true(grepl("Spectral_DB", row_a$annotation_method, fixed = TRUE))
  expect_equal(row_a$ref_id, "DB_A")

  row_b <- out[out$feature == "B", , drop = FALSE]
  expect_true(grepl("Network_DF", row_b$annotation_method, fixed = TRUE))
  expect_equal(row_b$ref_id, "NET_B")

  expect_true(is.na(out$ref_id[out$feature == "C"]))
})

test_that("apply_annotation_master_to_ms1 writes merged columns by name", {
  d <- tempfile("nts_apply_master")
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
      stringsAsFactors = FALSE
    )
  )
  p <- set_db_annotations(
    p,
    data.frame(
      name = "A",
      ref_id = "DB_A",
      Compound.name = "Compound A",
      SMILES = "C(F)(F)F",
      score_total = 0.9,
      stringsAsFactors = FALSE
    )
  )

  master <- merge_annotation_layers(p, include_all_features = TRUE, write_to_ms1 = TRUE)
  p2 <- attr(master, "project")
  expect_s4_class(p2, "NTSProject")
  ms1 <- get_ms1_annotations(p2)

  expect_true("annotation_method" %in% names(ms1))
  expect_true("ref_id" %in% names(ms1))
  expect_true("master_smiles" %in% names(ms1))
  expect_equal(ms1$ref_id[ms1$name == "A"], "DB_A")
  expect_equal(ms1$master_smiles[ms1$name == "A"], "C(F)(F)F")
  expect_true(is.na(ms1$ref_id[ms1$name == "B"]))

  p3 <- apply_annotation_master_to_ms1(p, annotation_master = master)
  expect_true("ref_id" %in% names(get_ms1_annotations(p3)))
  expect_false(identical(p, p3))
})

test_that("count_annotation_layers returns per-method counts", {
  d <- tempfile("nts_count_ann")
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  p <- create_nts_project(d)
  p <- set_db_annotations(
    p,
    data.frame(
      name = c("A", "B"),
      ref_id = c("DB_A", "DB_B"),
      stringsAsFactors = FALSE
    )
  )
  net <- list(
    nodes = data.frame(
      id = "C",
      ref_id = "NET_C",
      Propagation_Level = 2L,
      stringsAsFactors = FALSE
    ),
    edges = data.frame(from = "A", to = "C", stringsAsFactors = FALSE)
  )

  cnt <- count_annotation_layers(p, net_md_ann = net)
  expect_equal(nrow(cnt), 5L)
  expect_equal(cnt$n_features[cnt$method == "Spectral_DB"], 2L)
  expect_equal(cnt$n_features[cnt$method == "Network_DF"], 0L)
  expect_equal(cnt$n_features[cnt$method == "Network_MD"], 1L)
  expect_equal(cnt$n_propagated[cnt$method == "Network_MD"], 1L)
  expect_true(is.na(cnt$n_level0[cnt$method == "Spectral_DB"]))
})
