# tests/testthat/test-semantic-scaffold-inference.R
# Unit tests for DreaMS-paradigm Analog and Scaffold Retrieval

test_that("default_transformation_rules returns expected schema and rules", {
  rules <- default_transformation_rules()
  expect_s3_class(rules, "data.table")
  expect_true(nrow(rules) >= 25L)
  expect_true(all(c("rule_name", "mass_shift", "delta_formula", "category") %in% names(rules)))

  # Check key environmental rules
  rule_names <- rules$rule_name
  expect_true(any(grepl("Isomer", rule_names)))
  expect_true(any(grepl("\\+CF2", rule_names)))
  expect_true(any(grepl("-CF2", rule_names)))
  expect_true(any(grepl("\\+CH2", rule_names)))
  expect_true(any(grepl("Hydroxylation", rule_names)))
  expect_true(any(grepl("Decarboxylation", rule_names)))
  expect_true(any(grepl("Sulfonation", rule_names)))
})

test_that(".nts_match_transformation correctly identifies mass shifts", {
  rules <- .nts_default_environmental_transformations()

  # 1. Isomer / Identical mass (delta = 0)
  m_iso <- .nts_match_transformation(0.0001, rules, mz_tol_ppm = 10)
  expect_true(m_iso$matched)
  expect_equal(m_iso$rule_name, "Isomer / Identical mass")
  expect_equal(m_iso$category, "Isomer")

  # 2. +CF2 (49.9968 Da)
  m_cf2 <- .nts_match_transformation(49.9968, rules, mz_tol_ppm = 10)
  expect_true(m_cf2$matched)
  expect_true(grepl("\\+CF2", m_cf2$rule_name))
  expect_equal(m_cf2$delta_formula, "+CF2")

  # 3. -CF2 (-49.9968 Da)
  m_mcf2 <- .nts_match_transformation(-49.9968, rules, mz_tol_ppm = 10)
  expect_true(m_mcf2$matched)
  expect_true(grepl("-CF2", m_mcf2$rule_name))
  expect_equal(m_mcf2$delta_formula, "-CF2")

  # 4. +O (15.9949 Da)
  m_o <- .nts_match_transformation(15.9949, rules, mz_tol_ppm = 10)
  expect_true(m_o$matched)
  expect_true(grepl("Hydroxylation", m_o$rule_name))
  expect_equal(m_o$delta_formula, "+O")

  # 5. -CO2 (-43.9898 Da)
  m_co2 <- .nts_match_transformation(-43.9898, rules, mz_tol_ppm = 10)
  expect_true(m_co2$matched)
  expect_true(grepl("Decarboxylation", m_co2$rule_name))
  expect_equal(m_co2$delta_formula, "-CO2")

  # 6. Unmatched shift
  m_un <- .nts_match_transformation(123.4567, rules, mz_tol_ppm = 10)
  expect_false(m_un$matched)
  expect_equal(m_un$category, "Unclassified")
  expect_true(grepl("123.4567", m_un$rule_name))
})

test_that(".nts_extract_murcko_scaffold_batch handles cyclic and acyclic molecules", {
  smis <- c(
    PFOA = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O", # Acyclic
    Benzene = "c1ccccc1",                                                # Cyclic
    Atrazine = "CCNC1=NC(=NC(=N1)Cl)NC(C)C",                             # Heterocyclic
    Empty = ""
  )

  res <- .nts_extract_murcko_scaffold_batch(smis)
  expect_s3_class(res, "data.table")
  expect_equal(nrow(res), 3L) # Empty string is filtered out

  # Check cyclic vs acyclic flags
  row_pfoa <- res[smiles == smis["PFOA"]]
  expect_false(row_pfoa$is_cyclic)

  row_atrazine <- res[smiles == smis["Atrazine"]]
  expect_true(row_atrazine$is_cyclic)
})

test_that("infer_scaffold_from_semantic_neighbors works on mock data.frame inputs", {
  # Targets (Unknowns)
  ft_df <- data.frame(
    name = c("UNK_001", "UNK_002", "UNK_003"),
    mz = c(362.9691, 428.9613, 200.0000), # UNK_001 is -CF2 from PFOA (412.966), UNK_002 is +O from PFOA
    rt = c(600, 700, 100),
    stringsAsFactors = FALSE
  )

  # Anchor (Known standard)
  anchor_df <- data.frame(
    name = "PFOA_STD",
    Compound.name = "Perfluorooctanoic acid",
    SMILES = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",
    Formula = "C8HF15O2",
    exp_mz = 412.9664,
    score_total = 0.95,
    First_Class = "Perfluoroalkyl carboxylic acids",
    stringsAsFactors = FALSE
  )

  # Synthetic 256-d embeddings with exact cosine similarities
  set.seed(42)
  base_emb <- rnorm(256)
  base_emb <- base_emb / sqrt(sum(base_emb^2))

  # UNK_001: exact cosine = 0.98 to PFOA_STD
  u1 <- rnorm(256); u1 <- u1 - sum(u1 * base_emb) * base_emb; u1 <- u1 / sqrt(sum(u1^2))
  emb_unk1 <- 0.98 * base_emb + sqrt(1 - 0.98^2) * u1

  # UNK_002: exact cosine = 0.85
  u2 <- rnorm(256); u2 <- u2 - sum(u2 * base_emb) * base_emb; u2 <- u2 / sqrt(sum(u2^2))
  emb_unk2 <- 0.85 * base_emb + sqrt(1 - 0.85^2) * u2

  # UNK_003: completely orthogonal / unrelated (cosine = 0.10)
  u3 <- rnorm(256); u3 <- u3 - sum(u3 * base_emb) * base_emb; u3 <- u3 / sqrt(sum(u3^2))
  emb_unk3 <- 0.10 * base_emb + sqrt(1 - 0.10^2) * u3

  emb_mat <- rbind(
    PFOA_STD = base_emb,
    UNK_001 = emb_unk1,
    UNK_002 = emb_unk2,
    UNK_003 = emb_unk3
  )

  # Run inference
  res <- infer_scaffold_from_semantic_neighbors(
    project = ft_df,
    anchors = anchor_df,
    embeddings = emb_mat,
    min_semantic_score = 0.70,
    top_k = 3L,
    write_to_project = FALSE,
    verbose = FALSE
  )

  expect_type(res, "list")
  expect_true(all(c("summary", "analogs") %in% names(res)))

  sum_dt <- res$summary
  expect_equal(nrow(sum_dt), 3L)

  # UNK_001 should hit PFOA_STD with -CF2
  row1 <- sum_dt[feature_id == "UNK_001"]
  expect_equal(row1$analog_count, 1L)
  expect_equal(row1$top_analog_name, "Perfluorooctanoic acid")
  expect_true(row1$top_analog_similarity >= 0.90)
  expect_true(grepl("-CF2", row1$inferred_transformation))
  expect_false(row1$is_cyclic)
  expect_true(grepl("Analog to Perfluorooctanoic acid", row1$inferred_summary))
  expect_true(grepl("-CF2", row1$inferred_summary))

  # UNK_003 has low similarity -> 0 analog hits
  row3 <- sum_dt[feature_id == "UNK_003"]
  expect_equal(row3$analog_count, 0L)
  expect_true(is.na(row3$top_analog_name))
})

test_that("infer_scaffold_from_semantic_neighbors integrates with NTSProject and accessors", {
  # Construct minimal valid NTSProject
  prj <- new("NTSProject")
  ft <- data.table::data.table(
    name = c("FT01", "FT02", "FT03"),
    mz = c(412.9664, 362.9691, 500.0000),
    rt = c(600, 520, 400)
  )
  prj@feature_table <- ft

  # Annotations: FT01 is confirmed PFOA (anchor)
  prj@db_annotations <- data.table::data.table(
    name = "FT01",
    Compound.name = "PFOA",
    SMILES = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",
    Formula = "C8HF15O2",
    score_total = 0.95,
    First_Class = "Perfluoroalkyl carboxylic acids"
  )

  # Embeddings for FT01 and FT02 (FT03 has no embedding)
  set.seed(123)
  v1 <- rnorm(256)
  v1 <- v1 / sqrt(sum(v1^2))
  uv <- rnorm(256); uv <- uv - sum(uv * v1) * v1; uv <- uv / sqrt(sum(uv^2))
  v2 <- 0.98 * v1 + sqrt(1 - 0.98^2) * uv

  emb <- rbind(FT01 = v1, FT02 = v2)
  proxy_df <- data.frame(name = c("FT01", "FT02", "FT03"), stringsAsFactors = FALSE)
  proxy_df$embedding <- rbind(FT01 = v1, FT02 = v2, FT03 = rep(0, 256))
  prj@physchem_proxy <- proxy_df

  # Run inference on NTSProject
  prj_res <- infer_scaffold_from_semantic_neighbors(
    prj,
    embeddings = emb,
    min_semantic_score = 0.70,
    top_k = 2L,
    write_to_project = TRUE,
    verbose = FALSE
  )

  expect_s4_class(prj_res, "NTSProject")
  expect_true(!is.null(attr(prj_res@physchem_proxy, "analog_scaffold_inference")))

  # Access via get_inferred_scaffolds
  sum_tbl <- get_inferred_scaffolds(prj_res, type = "summary")
  expect_s3_class(sum_tbl, "data.table")
  expect_equal(nrow(sum_tbl), 2L) # FT02 and FT03 were targets

  # FT02 should hit FT01
  row_ft02 <- sum_tbl[feature_id == "FT02"]
  expect_equal(row_ft02$top_analog_id, "FT01")
  expect_equal(row_ft02$top_analog_name, "PFOA")
  expect_true(row_ft02$top_analog_similarity > 0.95)
  expect_true(grepl("-CF2", row_ft02$inferred_transformation))

  # FT03 had no MS2 embedding
  row_ft03 <- sum_tbl[feature_id == "FT03"]
  expect_equal(row_ft03$analog_count, 0L)
  expect_true(grepl("No MS2 spectrum", row_ft03$inferred_summary))

  # Check analogs long table
  ana_tbl <- get_inferred_scaffolds(prj_res, type = "analogs")
  expect_s3_class(ana_tbl, "data.table")
  expect_equal(nrow(ana_tbl), 1L)
  expect_equal(ana_tbl$feature_id, "FT02")
  expect_equal(ana_tbl$analog_id, "FT01")

  # Log updated
  expect_true(any(grepl("infer_scaffold_from_semantic_neighbors", prj_res@log)))
})

test_that(".nts_match_transformation handles edge cases and non-scalars defensively", {
  rules <- .nts_default_environmental_transformations()
  expect_false(.nts_match_transformation(numeric(0), rules)$matched)
  expect_false(.nts_match_transformation(NA_real_, rules)$matched)
  expect_false(.nts_match_transformation(Inf, rules)$matched)
  expect_false(.nts_match_transformation(c(10, 20), rules)$matched)
  expect_false(.nts_match_transformation(NULL, rules)$matched)
})

test_that("infer_scaffold_from_semantic_neighbors handles unnamed and non-finite embeddings matrices", {
  ft_df <- data.frame(
    name = c("FT_A", "FT_B"),
    mz = c(300, 350),
    rt = c(100, 120),
    stringsAsFactors = FALSE
  )
  anchors_df <- data.frame(
    name = "FT_A",
    Compound.name = "AnchorCompound",
    SMILES = "c1ccccc1",
    exp_mz = 300,
    stringsAsFactors = FALSE
  )

  # Matrix with NULL rownames: auto-derivation gives FT_A, FT_B
  set.seed(99)
  mat <- matrix(rnorm(2 * 256), nrow = 2, ncol = 256)
  rownames(mat) <- NULL

  res <- infer_scaffold_from_semantic_neighbors(
    project = ft_df,
    anchors = anchors_df,
    embeddings = mat,
    min_semantic_score = 0.0,
    top_k = 1L,
    write_to_project = FALSE,
    verbose = FALSE
  )
  expect_equal(nrow(res$summary), 1L) # FT_B is the target
  expect_equal(res$summary$feature_id, "FT_B")
  expect_equal(res$summary$top_analog_id, "FT_A")

  # Matrix where all rows are non-finite (NA)
  mat_na <- mat
  rownames(mat_na) <- c("FT_A", "FT_B")
  mat_na[, 5] <- NA # Corrupt all rows
  expect_warning(
    infer_scaffold_from_semantic_neighbors(
      project = ft_df,
      anchors = anchors_df,
      embeddings = mat_na,
      min_semantic_score = 0.0,
      top_k = 1L,
      write_to_project = FALSE,
      verbose = FALSE
    ),
    "No valid finite MS2 embeddings rows remaining"
  )
})

test_that("infer_scaffold_from_semantic_neighbors excludes self matches and preserves anchor alignment", {
  # 3 features: FT_1, FT_2, FT_3.
  ft_df <- data.frame(
    name = c("FT_1", "FT_2", "FT_3"),
    Compound.name = c("Comp_1", "Comp_2", "Comp_3"),
    SMILES = c("c1ccccc1", "c1ncncn1", "CCCC"),
    mz = c(100, 200, 300),
    rt = c(10, 20, 30),
    stringsAsFactors = FALSE
  )

  set.seed(42)
  v1 <- rnorm(256); v1 <- v1 / sqrt(sum(v1^2))
  u <- rnorm(256); u <- u - sum(u * v1) * v1; u <- u / sqrt(sum(u^2))
  v2 <- 0.95 * v1 + sqrt(1 - 0.95^2) * u # High similarity to FT_1
  v3 <- rnorm(256); v3 <- v3 / sqrt(sum(v3^2))
  emb <- rbind(FT_1 = v1, FT_2 = v2, FT_3 = v3)

  # Explicit anchor_features in reversed order
  res <- infer_scaffold_from_semantic_neighbors(
    project = ft_df,
    target_features = "FT_1",
    anchor_features = c("FT_3", "FT_2", "FT_1"),
    embeddings = emb,
    min_semantic_score = 0.50,
    top_k = 2L,
    write_to_project = FALSE,
    verbose = FALSE
  )

  sum_dt <- res$summary
  expect_equal(nrow(sum_dt), 1L)
  # FT_1 must NOT self-match
  expect_true(sum_dt$top_analog_id != "FT_1")
  expect_equal(sum_dt$top_analog_id, "FT_2")
  expect_equal(sum_dt$top_analog_name, "Comp_2")
})
