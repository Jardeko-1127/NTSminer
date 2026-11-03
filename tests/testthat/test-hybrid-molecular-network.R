# Tests for Hybrid Molecular Network (Rule-based Edges + Deep MS2 Semantic Embeddings)
# Verifies backward compatibility, edge synthesis, filtering, collapse, and visualization.

.make_hybrid_test_proj <- function(names = c("FT01", "FT02", "FT03")) {
  d <- tempfile("nts_hybrid_test_")
  dir.create(d)
  proj <- create_nts_project(d)
  ft <- data.frame(
    name = names,
    mz = seq(300.0, by = 50.0, length.out = length(names)),
    rt = seq(100.0, by = 50.0, length.out = length(names)),
    stringsAsFactors = FALSE
  )
  proj <- set_feature_table(proj, ft)
  ms1 <- data.frame(
    name = names,
    mz = ft$mz,
    rt = ft$rt,
    annotation_level = rep("Candidate", length(names)),
    stringsAsFactors = FALSE
  )
  set_ms1_annotations(proj, ms1)
}


test_that(".nts_network_ms2_semantic_edges builds pairwise edges from embedding matrix", {
  # 3 features: F1 and F2 are highly similar (cos = 0.99), F3 is orthogonal (cos = 0)
  v1 <- rep(0, 256)
  v1[1] <- 1.0

  v2 <- rep(0, 256)
  v2[1] <- 0.99
  v2[2] <- sqrt(1 - 0.99^2)

  v3 <- rep(0, 256)
  v3[3] <- 1.0

  emb <- rbind(v1, v2, v3)
  rownames(emb) <- c("FT01", "FT02", "FT03")

  # 1. min_score = 0.70 keeps only FT01 - FT02
  e <- NTSminer:::.nts_network_ms2_semantic_edges(
    project = NULL,
    valid_ids = c("FT01", "FT02", "FT03"),
    min_score = 0.70,
    top_k = 3L,
    embeddings = emb
  )

  expect_s3_class(e, "data.table")
  expect_equal(nrow(e), 1L)
  expect_equal(e$from, "FT01")
  expect_equal(e$to, "FT02")
  expect_equal(e$match_type, "semantic")
  expect_equal(e$Interaction, "MS2_Semantic")
  expect_equal(e$edge_score, 0.99)
  expect_equal(e$semantic_similarity, 0.99)
  expect_true(grepl("0.990 semantic", e$label))

  # 2. Lower min_score to 0.50: should still connect only FT01-FT02 since FT03 has 0 similarity
  e_low <- NTSminer:::.nts_network_ms2_semantic_edges(
    project = NULL,
    valid_ids = c("FT01", "FT02", "FT03"),
    min_score = 0.50,
    top_k = 3L,
    embeddings = emb
  )
  expect_equal(nrow(e_low), 1L)

  # 3. top_k = 1 with a triad
  # Add F4 similar to F1 (cos = 0.85)
  v4 <- rep(0, 256)
  v4[1] <- 0.85
  v4[4] <- sqrt(1 - 0.85^2)
  emb4 <- rbind(emb, v4)
  rownames(emb4) <- c("FT01", "FT02", "FT03", "FT04")

  e_top1 <- NTSminer:::.nts_network_ms2_semantic_edges(
    project = NULL,
    valid_ids = c("FT01", "FT02", "FT03", "FT04"),
    min_score = 0.70,
    top_k = 1L,
    embeddings = emb4
  )
  expect_true("FT01" %in% e_top1$from || "FT01" %in% e_top1$to)

  # 4. Graceful handling of < 2 valid features
  e_empty <- NTSminer:::.nts_network_ms2_semantic_edges(
    project = NULL,
    valid_ids = c("FT01"),
    embeddings = emb
  )
  expect_equal(nrow(e_empty), 0L)
  expect_true("semantic_similarity" %in% names(e_empty))
})


test_that("get_network_data supports opt-in semantic edges and keeps backward compatibility", {
  proj <- .make_hybrid_test_proj(c("FT01", "FT02", "FT03"))

  # Create embeddings for FT01 and FT02
  v1 <- rep(0, 256)
  v1[1] <- 1.0
  v2 <- rep(0, 256)
  v2[1] <- 0.95
  v2[2] <- sqrt(1 - 0.95^2)
  v3 <- rep(0, 256)
  v3[3] <- 1.0

  emb <- rbind(v1, v2, v3)
  rownames(emb) <- c("FT01", "FT02", "FT03")

  # Case 1: Default include_semantic = FALSE (100% backward compatible)
  net_default <- get_network_data(
    proj,
    type = "all",
    ms2_edge_types = character(0),
    ms1_edge_types = character(0),
    include_semantic = FALSE
  )
  expect_equal(nrow(net_default$edges), 0L)

  # Case 2: include_semantic = TRUE with external embeddings
  net_semantic <- get_network_data(
    proj,
    type = "all",
    ms2_edge_types = character(0),
    ms1_edge_types = character(0),
    include_semantic = TRUE,
    min_semantic_score = 0.70,
    embeddings = emb
  )
  expect_equal(nrow(net_semantic$edges), 1L)
  expect_equal(net_semantic$edges$from, "FT01")
  expect_equal(net_semantic$edges$to, "FT02")
  expect_equal(net_semantic$edges$edge_category, "Semantic")
  expect_equal(net_semantic$edges$semantic_similarity, 0.95)

  # Case 3: Slice type = "semantic"
  net_slice <- get_network_data(
    proj,
    type = "semantic",
    min_semantic_score = 0.70,
    embeddings = emb,
    keep_isolated = FALSE
  )
  expect_equal(nrow(net_slice$edges), 1L)
  expect_equal(nrow(net_slice$nodes), 2L) # Only FT01 and FT02 kept when keep_isolated = FALSE
  expect_setequal(net_slice$nodes$id, c("FT01", "FT02"))
})


test_that(".nts_merged_network_edge_collapse folds semantic_similarity properly", {
  # Node pair with both Fragment match and Semantic similarity
  E <- data.table::data.table(
    from = c("FT01", "FT01"),
    to = c("FT02", "FT02"),
    edge_type = c(NA_character_, NA_character_),
    match_type = c("fragment", "semantic"),
    Interaction = c("MS2_Shared_Fragment", "MS2_Semantic"),
    label = c("2 fragments", "0.950 semantic"),
    edge_category = c("Fragment", "Semantic"),
    shared_fragments_count = c(2L, NA_integer_),
    shared_fragments_mz = c("119.0, 169.0", NA_character_),
    shared_nls_count = c(NA_integer_, NA_integer_),
    shared_nls_mz = c(NA_character_, NA_character_),
    shared_mds_count = c(NA_integer_, NA_integer_),
    shared_mds_mz = c(NA_character_, NA_character_),
    edge_score = c(0.80, 0.95),
    semantic_similarity = c(NA_real_, 0.95),
    n_steps = c(NA_integer_, NA_integer_),
    transformation = c(NA_character_, NA_character_),
    mass_difference = c(NA_real_, NA_real_),
    ms2_ladder_confirmed = c(NA, NA)
  )

  collapsed <- NTSminer:::.nts_merged_network_edge_collapse(E)
  expect_equal(nrow(collapsed), 1L)
  expect_equal(collapsed$from, "FT01")
  expect_equal(collapsed$to, "FT02")
  expect_equal(collapsed$shared_fragments_count, 2L)
  expect_equal(collapsed$shared_fragments_mz, "119.0, 169.0")
  expect_equal(collapsed$semantic_similarity, 0.95)
  expect_equal(collapsed$edge_score, 0.95) # max(0.80, 0.95)
  expect_true(grepl("MS2_Shared_Fragment", collapsed$Interaction))
  expect_true(grepl("MS2_Semantic", collapsed$Interaction))
})


test_that(".nts_edge_network_vis styles Semantic edges with purple color and dashes", {
  edges_df <- data.frame(
    from = "FT01",
    to = "FT02",
    Interaction = "MS2_Semantic",
    match_type = "semantic",
    label = "0.920 semantic",
    edge_category = "Semantic",
    semantic_similarity = 0.92,
    stringsAsFactors = FALSE
  )

  styled <- NTSminer:::.nts_edge_network_vis(edges_df, color_by_edge = TRUE)
  expect_equal(styled$color, "#9370DB") # Medium Purple
  expect_equal(styled$width, 1) # MS2 width
  expect_equal(styled$dashes[[1L]], c(5, 5))
  expect_equal(styled$arrows[[1L]], "") # Undirected look
  expect_true(grepl("Semantic similarity", styled$title))
  expect_true(grepl("0.920", styled$title))
})


test_that("generate_network_data passes semantic arguments and returns correct schema", {
  proj <- .make_hybrid_test_proj(c("FT01", "FT02"))

  v1 <- rep(0, 256); v1[1] <- 1.0
  v2 <- rep(0, 256); v2[1] <- 0.90; v2[2] <- sqrt(1 - 0.90^2)
  emb <- rbind(v1, v2)
  rownames(emb) <- c("FT01", "FT02")

  net <- generate_network_data(
    proj,
    ms1_edge_types = character(0),
    ms2_edge_types = character(0),
    include_semantic = TRUE,
    min_semantic_score = 0.80,
    embeddings = emb
  )

  expect_type(net, "list")
  expect_s3_class(net$nodes, "data.frame")
  expect_s3_class(net$edges, "data.frame")
  expect_equal(nrow(net$edges), 1L)
  expect_true("semantic_similarity" %in% names(net$edges))
  expect_equal(net$edges$semantic_similarity, 0.90)

  # Check plot_molecular_network accepts the hybrid network
  skip_if_not_installed("visNetwork")
  p <- plot_molecular_network(net, remove_isolated = FALSE)
  expect_s3_class(p, "visNetwork")
})


test_that("semantic edges survive when ms2_annotations only contains a subset of features (P1 fix)", {
  proj <- .make_hybrid_test_proj(c("FT01", "FT02", "FT03"))
  
  # Simulate ms2_annotations only recording FT01
  proj@ms2_annotations <- data.frame(
    name = "FT01",
    precursor_mz = 300.0,
    rt = 100.0,
    stringsAsFactors = FALSE
  )

  v1 <- rep(0, 256); v1[1] <- 1.0
  v2 <- rep(0, 256); v2[1] <- 0.92; v2[2] <- sqrt(1 - 0.92^2)
  v3 <- rep(0, 256); v3[1] <- 0.88; v3[3] <- sqrt(1 - 0.88^2)
  emb <- rbind(v1, v2, v3)
  rownames(emb) <- c("FT01", "FT02", "FT03")

  net <- get_network_data(
    proj,
    type = "semantic",
    min_semantic_score = 0.80,
    embeddings = emb,
    keep_isolated = FALSE
  )

  # FT02-FT03 and FT01-FT02 should both exist even though FT02/FT03 are not in ms2_annotations
  expect_true(nrow(net$edges) >= 2L)
  expect_true(any(net$edges$from == "FT02" | net$edges$to == "FT02"))
  expect_true(any(net$edges$from == "FT03" | net$edges$to == "FT03"))
})


test_that("unnamed embeddings and non-finite rows are handled defensively (P2 fixes)", {
  # 1. Unnamed matrix of matching length takes valid_ids
  v1 <- rep(0, 256); v1[1] <- 1.0
  v2 <- rep(0, 256); v2[1] <- 0.95; v2[2] <- sqrt(1 - 0.95^2)
  emb_unnamed <- rbind(v1, v2)
  rownames(emb_unnamed) <- NULL # explicitly clear rownames

  e_unnamed <- NTSminer:::.nts_network_ms2_semantic_edges(
    project = NULL,
    valid_ids = c("NodeA", "NodeB"),
    min_score = 0.80,
    embeddings = emb_unnamed
  )
  expect_equal(nrow(e_unnamed), 1L)
  expect_equal(e_unnamed$from, "NodeA")
  expect_equal(e_unnamed$to, "NodeB")

  # 2. Embeddings with Inf or NA row are safely filtered
  v_bad <- rep(NA_real_, 256)
  emb_corrupt <- rbind(v1, v2, v_bad)
  rownames(emb_corrupt) <- c("NodeA", "NodeB", "NodeC")

  e_clean <- NTSminer:::.nts_network_ms2_semantic_edges(
    project = NULL,
    valid_ids = c("NodeA", "NodeB", "NodeC"),
    min_score = 0.80,
    embeddings = emb_corrupt
  )
  expect_equal(nrow(e_clean), 1L)
  expect_equal(e_clean$from, "NodeA")
  expect_equal(e_clean$to, "NodeB")
})

