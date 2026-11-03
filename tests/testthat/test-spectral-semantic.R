# Tests for Spectral Semantic Embeddings & Matching Module
# (Spec2Vec Paradigm via Transformer Encoder)

test_that("match_spectral_semantic correctly calculates vectorized cosine similarities", {
  # Create synthetic 256-dimensional unit vectors
  # Vector 1: pure unit vector along dimension 1
  v1 <- rep(0, 256)
  v1[1] <- 1.0

  # Vector 2: pure unit vector along dimension 2 (orthogonal to v1)
  v2 <- rep(0, 256)
  v2[2] <- 1.0

  # Vector 3: 45-degree angle in plane of dims 1 and 2
  v3 <- rep(0, 256)
  v3[1] <- 1 / sqrt(2)
  v3[2] <- 1 / sqrt(2)

  # Query matrix (1 query: v1)
  query_mat <- matrix(v1, nrow = 1, ncol = 256)
  rownames(query_mat) <- "Query_01"

  # Reference matrix (3 compounds: v1, v2, v3)
  ref_mat <- rbind(v1, v2, v3)
  rownames(ref_mat) <- c("Ref_ExactMatch", "Ref_Orthogonal", "Ref_Partial")

  res <- match_spectral_semantic(query_mat, ref_mat, top_k = 3L)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3L)
  expect_equal(res$query_id, rep("Query_01", 3))
  expect_equal(res$rank, 1:3)

  # Exact match should rank 1st with cosine = 1.0
  expect_equal(res$reference_id[1], "Ref_ExactMatch")
  expect_equal(res$semantic_similarity[1], 1.0)

  # Partial match should rank 2nd with cosine = 1 / sqrt(2) = 0.7071
  expect_equal(res$reference_id[2], "Ref_Partial")
  expect_equal(res$semantic_similarity[2], round(1 / sqrt(2), 4))

  # Orthogonal match should rank 3rd with cosine = 0.0
  expect_equal(res$reference_id[3], "Ref_Orthogonal")
  expect_equal(res$semantic_similarity[3], 0.0)

  # min_similarity cutoff filtering
  res_filtered <- match_spectral_semantic(query_mat, ref_mat, top_k = 3L, min_similarity = 0.50)
  expect_equal(nrow(res_filtered), 2L)
  expect_true(all(res_filtered$semantic_similarity >= 0.50))
})


test_that("match_spectral_semantic handles data.frame input and dimension errors gracefully", {
  # 1. Invalid dimensions throw informative error
  bad_q <- matrix(runif(100), nrow = 1, ncol = 100)
  bad_r <- matrix(runif(100), nrow = 1, ncol = 100)
  expect_error(match_spectral_semantic(bad_q, bad_r), "exactly 256 dimensions")

  # 2. Data.frame with embedding matrix column
  q_mat <- matrix(0, nrow = 2, ncol = 256)
  q_mat[1, 1] <- 1.0
  q_mat[2, 2] <- 1.0
  rownames(q_mat) <- c("Q1", "Q2")

  r_mat <- matrix(0, nrow = 2, ncol = 256)
  r_mat[1, 1] <- 1.0
  r_mat[2, 3] <- 1.0
  rownames(r_mat) <- c("R1", "R2")

  df_q <- data.frame(feature_id = c("Q1", "Q2"), stringsAsFactors = FALSE)
  df_q$embedding <- q_mat

  df_r <- data.frame(feature_id = c("R1", "R2"), stringsAsFactors = FALSE)
  df_r$embedding <- r_mat

  res <- match_spectral_semantic(df_q, df_r, top_k = 1L)
  expect_equal(nrow(res), 2L)
  expect_equal(res$reference_id[res$query_id == "Q1"], "R1")
  expect_equal(res$semantic_similarity[res$query_id == "Q1"], 1.0)
})


test_that("extract_ms2_embeddings validates inputs and handles missing python gracefully", {
  fake_py <- file.path(tempdir(), "non_existent_python_binary_xyz.exe")
  pks <- matrix(c(100, 200, 50, 100), ncol = 2)

  # Non-existent Python path throws informative error
  expect_error(
    extract_ms2_embeddings(pks, precursor_mz = 412.9664, python_path = fake_py),
    "Python executable not found"
  )

  # Missing precursor_mz throws informative error
  expect_error(
    extract_ms2_embeddings(pks, precursor_mz = NULL, python_path = "some_python")
  )
})


test_that("extract_ms2_embeddings works with mock python CLI", {
  mock_py_script <- tempfile(fileext = ".bat")
  mock_r_worker <- tempfile(fileext = ".R")

  # Write mock R worker script simulating Python embed output
  cat('
    args <- commandArgs(trailingOnly = TRUE)
    in_idx <- which(args == "--input_json") + 1L
    out_idx <- which(args == "--output_json") + 1L
    in_file <- args[in_idx]
    out_file <- args[out_idx]

    inputs <- jsonlite::fromJSON(in_file, simplifyVector = FALSE)
    results <- lapply(inputs, function(it) {
      # 256-dim mock unit embedding with 1.0 at dim 1
      emb <- rep(0, 256)
      emb[1] <- 1.0
      list(
        id = it$id,
        precursor_mz = it$precursor_mz,
        n_peaks = length(it$peaks),
        embedding = as.list(emb)
      )
    })
    jsonlite::write_json(results, out_file, auto_unbox = TRUE, digits = 6)
  ', file = mock_r_worker)

  r_exe <- file.path(R.home("bin"), "Rscript.exe")
  cat(sprintf('@"%s" "%s" %%*\n', r_exe, mock_r_worker), file = mock_py_script)

  pks <- matrix(c(100.0, 200.0, 50.0, 100.0), ncol = 2)

  # Run extraction with mock script
  res_df <- extract_ms2_embeddings(
    object = pks,
    precursor_mz = 412.9664,
    python_path = mock_py_script,
    verbose = FALSE
  )

  expect_s3_class(res_df, "data.frame")
  expect_equal(nrow(res_df), 1L)
  expect_equal(res_df$feature_id, "Spectrum_0001")
  expect_equal(res_df$precursor_mz, 412.9664)
  expect_true(is.matrix(res_df$embedding))
  expect_equal(ncol(res_df$embedding), 256L)
  expect_equal(res_df$embedding[1, 1], 1.0)

  # Test as_matrix = TRUE
  res_mat <- extract_ms2_embeddings(
    object = pks,
    precursor_mz = 412.9664,
    python_path = mock_py_script,
    as_matrix = TRUE,
    verbose = FALSE
  )
  expect_true(is.matrix(res_mat))
  expect_equal(dim(res_mat), c(1L, 256L))
  expect_equal(res_mat[1, 1], 1.0)
})
