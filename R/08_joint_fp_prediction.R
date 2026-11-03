# =============================================================================
# 08_joint_fp_prediction.R
# Production R Interface for Joint 2215-bit Fingerprint & Ecotoxicity Model
# (MS2 -> Joint 2215-bit Fingerprint [ECFP2048 + MACCS167] -> Downstream Tox)
# =============================================================================

#' Load configuration and benchmark metrics of the joint fingerprint model
#'
#' Reads the bundled or user-specified JSON configuration file of the trained
#' joint 2215-bit fingerprint model (ECFP2048 + MACCS167 offset by 2048),
#' including neural architecture hyperparameters, training convergence history,
#' and independent golden benchmark metrics on 174 held-out PFAS scaffolds.
#'
#' @param config_path Optional path to `joint_fp_model_config.json`. If `NULL`,
#'   the bundled model configuration in `inst/extdata/models/joint_fp/` is resolved.
#'
#' @return An S3 object of class `nts_joint_fp_config` containing:
#'   \itemize{
#'     \item `architecture`: Network dimensions, vocabularies, and sequence lengths.
#'     \item `training_summary`: Train/validation dataset scale and loss trajectory.
#'     \item `pfas_benchmark`: Independent metrics on 174 held-out PFAS scaffolds.
#'     \item `downstream_toxicity`: Spearman rho, Pearson r, and MAE for QSAR ecotoxicity.
#'   }
#'
#' @export
load_joint_fp_model_config <- function(config_path = NULL) {
  if (is.null(config_path)) {
    paths <- .nts_resolve_joint_fp_model_paths()
    config_path <- paths$config
  }
  if (!file.exists(config_path)) {
    stop("Joint fingerprint model configuration not found: ", config_path, call. = FALSE)
  }

  raw_cfg <- jsonlite::fromJSON(config_path, simplifyVector = TRUE)
  structure(raw_cfg, class = c("nts_joint_fp_config", "list"))
}

#' @export
print.nts_joint_fp_config <- function(x, ...) {
  cat("======================================================================\n")
  cat("  NTSminer Joint 2215-bit Fingerprint & Ecotoxicity Model\n")
  cat("  (MS2 Spectrum -> ECFP2048 + MACCS167 -> Downstream QSAR Risk)\n")
  cat("======================================================================\n")
  
  arch <- x$architecture
  if (!is.null(arch)) {
    cat(sprintf("Model Type:          %s (Autoregressive Decoder)\n", arch$model_type %||% "Transformer"))
    cat(sprintf("Vocabulary Sizes:    m/z bins = %d, Fingerprint tokens = %d\n", arch$num_mz_vocabs %||% 0, arch$vocab_size %||% 0))
    cat(sprintf("Embedding / Layers:  d_model = %d, Heads = %d, Enc = %d, Dec = %d\n",
                arch$d_model %||% 0, arch$num_heads %||% 0, arch$num_encoder_layers %||% 0, arch$num_decoder_layers %||% 0))
    cat(sprintf("Sequence Constraints:Max Seq Len = %d bits, Max MS2 Peaks = %d\n", arch$max_sequence_length %||% 0, arch$max_peaks %||% 0))
  }

  tr <- x$training_summary
  if (!is.null(tr)) {
    cat("----------------------------------------------------------------------\n")
    cat(sprintf("Training Corpus:     %s spectra (%s unique scaffolds)\n",
                format(tr$training_spectra %||% 0, big.mark = ","),
                format(tr$training_scaffolds %||% 0, big.mark = ",")))
    cat(sprintf("Validation Corpus:   %s spectra (%s unique scaffolds)\n",
                format(tr$validation_spectra %||% 0, big.mark = ","),
                format(tr$validation_scaffolds %||% 0, big.mark = ",")))
    cat(sprintf("Loss / Convergence:  Train Loss = %s, Val Loss = %s\n",
                round(as.numeric(tr$train_loss %||% 0), 4),
                round(as.numeric(tr$val_loss %||% 0), 4)))
    cat(sprintf("Validation Tanimoto: %s (Held-out Evaluation)\n",
                round(as.numeric(tr$val_tanimoto %||% 0), 4)))
  }

  bench <- x$pfas_benchmark_performance
  if (!is.null(bench)) {
    cat("----------------------------------------------------------------------\n")
    cat(sprintf("PFAS Golden Benchmark (%d Scaffolds 100%% Isolated, %s Spectra):\n",
                as.integer(bench$isolated_scaffolds %||% 0L),
                format(bench$evaluated_spectra %||% 0, big.mark = ",")))
    cat(sprintf("  - Joint 2215-bit Tanimoto:  Mean = %s, Median = %s\n",
                round(as.numeric(bench$joint_2215_tanimoto_mean %||% 0), 4),
                round(as.numeric(bench$joint_2215_tanimoto_median %||% 0), 4)))
    cat(sprintf("  - MACCS167 Key Moieties:    Mean = %s\n",
                round(as.numeric(bench$maccs167_tanimoto_mean %||% 0), 4)))
    cat(sprintf("  - ECFP2048 Substructures:   Mean = %s\n",
                round(as.numeric(bench$ecfp2048_tanimoto_mean %||% 0), 4)))

    tox <- bench$downstream_toxicity_fidelity
    if (!is.null(tox)) {
      cat(sprintf("  - Downstream Ecotox Fidelity: Spearman rho = %s, Pearson r = %s, MAE = %s\n",
                  round(as.numeric(tox$spearman_rho %||% 0), 4),
                  round(as.numeric(tox$pearson_r %||% 0), 4),
                  round(as.numeric(tox$mae %||% 0), 4)))
    }
  }
  cat("======================================================================\n")
  invisible(x)
}


#' Predict joint 2215-bit fingerprints and ecotoxicity directly from MS2 spectra
#'
#' Deploys the pre-trained deep learning Transformer model to translate raw
#' MS2 fragmentation spectra and precursor m/z directly into a combined 2215-bit
#' binary molecular fingerprint (ECFP2048 and MACCS167 with bit offset 2048),
#' detects key environmental functional groups (e.g., C-F aliphatic chains,
#' carboxylic acid, sulfonates), and scores downstream ecotoxicity risk.
#'
#' This enables structural and toxicological characterization of unknown features
#' even when no reference library spectrum or confirmed SMILES structure exists.
#'
#' @param object An [NTSProject], a two-column `matrix`/`data.frame` of peaks
#'   (`mz`, `intensity`), or a list of spectra.
#' @param precursor_mz Precursor ion m/z. Required when `object` is a raw peak
#'   matrix or data frame.
#' @param peaks Two-column numeric matrix or data frame of MS2 peaks (`mz`, `intensity`).
#'   Optional if provided via `object`.
#' @param python_path Path to Python executable with PyTorch installed. If `NULL`,
#'   automatically resolved from environment variables (`NTSMINER_PYTHON`,
#'   `RETICULATE_PYTHON`) or system PATH.
#' @param model_weights Optional path to `.pt` model weights. If `NULL`, resolved
#'   from package bundled assets (`inst/extdata/models/joint_fp/`).
#' @param model_config Optional path to `.json` model configuration. If `NULL`,
#'   resolved from package bundled assets.
#' @param tox_model Optional path to downstream `tox_model.joblib`. If `NULL`,
#'   resolved from package bundled assets.
#' @param device Computing device: `"auto"` (default, uses CUDA if available),
#'   `"cuda"`, or `"cpu"`.
#' @param write_to_project For [NTSProject] input: if `TRUE` (default), merges
#'   predicted fingerprints, toxicity scores, and functional groups into
#'   `@physchem_proxy`.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param ... Additional arguments.
#'
#' @return When `object` is an [NTSProject] and `write_to_project = TRUE`, returns
#'   the updated project object. Otherwise, returns a `data.frame` containing:
#'   \itemize{
#'     \item `feature_id`: Feature identifier.
#'     \item `precursor_mz`: Precursor ion m/z.
#'     \item `n_total_bits`: Total predicted active on-bits (0..2214).
#'     \item `n_ecfp_bits`: Predicted ECFP2048 on-bits (0..2047).
#'     \item `n_maccs_bits`: Predicted MACCS167 on-bits (0..166).
#'     \item `joint_fp_bits`: Semicolon-delimited string of all active bit indices.
#'     \item `toxicity_score`: Continuous downstream ecotoxicity score (0.0 to 1.0).
#'     \item `toxicity_risk`: Categorical risk classification ("High", "Moderate", "Low", "Unknown").
#'     \item `detected_functional_groups`: Semicolon-delimited recognized chemical moieties.
#'   }
#'
#' @seealso [predict_toxicity_from_ms2], [load_joint_fp_model_config]
#'
#' @export
predict_joint_fingerprint <- function(object,
                                     precursor_mz = NULL,
                                     peaks = NULL,
                                     python_path = NULL,
                                     model_weights = NULL,
                                     model_config = NULL,
                                     tox_model = NULL,
                                     device = c("auto", "cuda", "cpu"),
                                     write_to_project = TRUE,
                                     verbose = TRUE,
                                     ...) {
  device <- match.arg(device)

  # 1. Resolve Python executable
  py_exe <- .nts_resolve_python_executable(python_path)
  if (is.null(py_exe) || !nzchar(py_exe)) {
    stop(
      "Python executable not found.\n",
      "  Please install Python (>= 3.8) with PyTorch (torch) and scikit-learn / joblib,\n",
      "  and specify `python_path` or set the environment variable `NTSMINER_PYTHON`.",
      call. = FALSE
    )
  }

  # 2. Resolve Model and Script Assets
  paths <- .nts_resolve_joint_fp_model_paths(
    weights = model_weights,
    config = model_config,
    tox = tox_model
  )
  engine_script <- .nts_resolve_joint_fp_engine_script()

  if (is.null(paths$weights) || !file.exists(paths$weights)) {
    stop("Joint fingerprint model weights not found. Please provide valid `model_weights`.", call. = FALSE)
  }
  if (is.null(paths$config) || !file.exists(paths$config)) {
    stop("Joint fingerprint model config not found. Please provide valid `model_config`.", call. = FALSE)
  }
  if (!file.exists(engine_script)) {
    stop("Python inference engine script not found at: ", engine_script, call. = FALSE)
  }

  # 3. Parse and standardize spectra input
  is_proj <- inherits(object, "NTSProject")
  spectra_list <- list()

  if (is_proj) {
    spectra_list <- .nts_extract_ms2_from_project(object)
    if (!length(spectra_list)) {
      warning("No valid MS2 spectra found in NTSProject.", call. = FALSE)
      return(object)
    }
  } else if (is.data.frame(object) || is.matrix(object) || !is.null(peaks)) {
    pk_mat <- if (!is.null(peaks)) peaks else object
    if (is.null(precursor_mz) || !is.numeric(precursor_mz) || is.na(precursor_mz)) {
      stop("`precursor_mz` must be provided for raw peak data.", call. = FALSE)
    }
    pk_df <- as.data.frame(pk_mat)
    if (ncol(pk_df) < 2L) {
      stop("Peak table must contain at least two columns: m/z and intensity.", call. = FALSE)
    }
    mz_col <- pk_df[[1L]]
    int_col <- pk_df[[2L]]
    clean_pks <- lapply(seq_along(mz_col), function(i) {
      list(as.numeric(mz_col[i]), as.numeric(int_col[i]))
    })
    spectra_list <- list(list(
      id = "Spectrum_0001",
      precursor_mz = as.numeric(precursor_mz),
      peaks = clean_pks
    ))
  } else if (is.list(object)) {
    spectra_list <- object
  } else {
    stop("`object` must be an NTSProject, a peak matrix/data.frame, or a list of spectra.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message(sprintf("Running Joint Fingerprint & Toxicity Inference on %d spectra (device: %s)...",
                    length(spectra_list), device))
  }

  # 4. Exchange data via temporary JSON files
  tmp_in <- tempfile(pattern = "joint_fp_in_", fileext = ".json")
  tmp_out <- tempfile(pattern = "joint_fp_out_", fileext = ".json")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)

  jsonlite::write_json(spectra_list, tmp_in, auto_unbox = TRUE, digits = 8)

  cmd_args <- c(
    shQuote(engine_script),
    "--input_json", shQuote(tmp_in),
    "--output_json", shQuote(tmp_out),
    "--weights_path", shQuote(paths$weights),
    "--config_path", shQuote(paths$config),
    "--device", device
  )
  if (!is.null(paths$tox) && file.exists(paths$tox)) {
    cmd_args <- c(cmd_args, "--tox_model_path", shQuote(paths$tox))
  }

  run_res <- suppressWarnings(system2(py_exe, args = cmd_args, stdout = TRUE, stderr = TRUE))
  status_code <- attr(run_res, "status") %||% 0L

  if (status_code != 0L || !file.exists(tmp_out)) {
    err_msg <- paste(run_res, collapse = "\n")
    stop(
      "Joint fingerprint inference failed (exit code ", status_code, "):\n",
      err_msg,
      call. = FALSE
    )
  }

  res_json <- jsonlite::fromJSON(tmp_out, simplifyVector = FALSE)

  # 5. Format results into clean data.frame
  rows <- lapply(res_json, function(it) {
    j_bits <- if (length(it$joint_fp_on_bits)) paste(it$joint_fp_on_bits, collapse = ";") else ""
    e_bits <- if (length(it$ecfp_on_bits)) paste(it$ecfp_on_bits, collapse = ";") else ""
    m_bits <- if (length(it$maccs_on_bits)) paste(it$maccs_on_bits, collapse = ";") else ""
    groups <- if (length(it$detected_functional_groups)) paste(it$detected_functional_groups, collapse = "; ") else "None"
    
    data.frame(
      feature_id = as.character(it$id %||% ""),
      precursor_mz = as.numeric(it$precursor_mz %||% NA_real_),
      n_total_bits = as.integer(it$n_total_bits %||% 0L),
      n_ecfp_bits = as.integer(it$n_ecfp_bits %||% 0L),
      n_maccs_bits = as.integer(it$n_maccs_bits %||% 0L),
      joint_fp_bits = j_bits,
      ecfp_bits = e_bits,
      maccs_bits = m_bits,
      toxicity_score = if (!is.null(it$toxicity_score)) as.numeric(it$toxicity_score) else NA_real_,
      toxicity_risk = as.character(it$toxicity_risk %||% "Unknown"),
      detected_functional_groups = groups,
      stringsAsFactors = FALSE
    )
  })

  df_out <- do.call(rbind, rows)
  rownames(df_out) <- NULL

  # 6. If NTSProject and write_to_project, update project slot
  if (is_proj && isTRUE(write_to_project)) {
    object <- .nts_write_joint_fp_to_project(object, df_out)
    return(object)
  }

  df_out
}


#' Quick ecotoxicity assessment directly from MS2 fragmentation
#'
#' A high-level convenience wrapper around [predict_joint_fingerprint] designed
#' specifically for ecotoxicological prioritization of non-target screening features.
#'
#' @param object An [NTSProject], peak matrix, or list of spectra.
#' @param ... Passed directly to [predict_joint_fingerprint].
#'
#' @return A `data.frame` of toxicity scores and functional group flags, or an
#'   updated [NTSProject].
#'
#' @export
predict_toxicity_from_ms2 <- function(object, ...) {
  predict_joint_fingerprint(object, ...)
}


# =============================================================================
# Spectral Semantic Embeddings (Spec2Vec Paradigm via Transformer Encoder)
# =============================================================================

#' Extract 256-dimensional semantic spectrum embeddings via Transformer Encoder
#'
#' Deploys the pre-trained Transformer Encoder to compress raw MS2 fragmentation
#' spectra and precursor m/z into 256-dimensional continuous dense embeddings.
#' These embeddings capture structural and fragmentation semantics (Spec2Vec paradigm)
#' without relying on hard m/z matching.
#'
#' @param object An [NTSProject], a two-column peak `matrix`/`data.frame` (`mz`, `intensity`),
#'   or a list of spectra.
#' @param precursor_mz Precursor ion m/z. Required when `object` is a raw peak matrix.
#' @param peaks Two-column numeric matrix or data frame of MS2 peaks (`mz`, `intensity`).
#' @param python_path Path to Python executable. If `NULL`, auto-resolved.
#' @param model_weights Optional path to `.pt` weights.
#' @param model_config Optional path to `.json` configuration.
#' @param device Computing device: `"auto"`, `"cuda"`, or `"cpu"`.
#' @param as_matrix Logical. If `TRUE`, returns an \eqn{N \times 256} numeric matrix
#'   with rownames = feature IDs. Default `FALSE`.
#' @param verbose Logical. If `TRUE`, prints progress messages.
#' @param ... Additional arguments.
#'
#' @return When `as_matrix = TRUE`, an \eqn{N \times 256} numeric matrix of L2-normalized
#'   embeddings. Otherwise, a `data.frame` containing:
#'   \itemize{
#'     \item `feature_id`: Feature identifier.
#'     \item `precursor_mz`: Precursor ion m/z.
#'     \item `embedding`: A 256-column numeric matrix of embeddings.
#'   }
#'
#' @seealso [match_spectral_semantic], [predict_joint_fingerprint]
#' @export
extract_ms2_embeddings <- function(object,
                                   precursor_mz = NULL,
                                   peaks = NULL,
                                   python_path = NULL,
                                   model_weights = NULL,
                                   model_config = NULL,
                                   device = c("auto", "cuda", "cpu"),
                                   as_matrix = FALSE,
                                   verbose = TRUE,
                                   ...) {
  device <- match.arg(device)

  # 1. Resolve Python executable
  py_exe <- .nts_resolve_python_executable(python_path)
  if (is.null(py_exe) || !nzchar(py_exe)) {
    stop(
      "Python executable not found.\n",
      "  Please install Python (>= 3.8) with PyTorch (torch),\n",
      "  and specify `python_path` or set the environment variable `NTSMINER_PYTHON`.",
      call. = FALSE
    )
  }

  # 2. Resolve Model and Script Assets
  paths <- .nts_resolve_joint_fp_model_paths(
    weights = model_weights,
    config = model_config
  )
  engine_script <- .nts_resolve_joint_fp_engine_script()

  if (is.null(paths$weights) || !file.exists(paths$weights)) {
    stop("Joint fingerprint model weights not found.", call. = FALSE)
  }
  if (is.null(paths$config) || !file.exists(paths$config)) {
    stop("Joint fingerprint model config not found.", call. = FALSE)
  }
  if (!file.exists(engine_script)) {
    stop("Python inference engine script not found at: ", engine_script, call. = FALSE)
  }

  # 3. Parse and standardize spectra input
  is_proj <- inherits(object, "NTSProject")
  spectra_list <- list()

  if (is_proj) {
    spectra_list <- .nts_extract_ms2_from_project(object)
    if (!length(spectra_list)) {
      warning("No valid MS2 spectra found in NTSProject.", call. = FALSE)
      if (isTRUE(as_matrix)) return(matrix(numeric(0), nrow = 0, ncol = 256))
      return(data.frame(feature_id = character(0), precursor_mz = numeric(0)))
    }
  } else if (is.data.frame(object) || is.matrix(object) || !is.null(peaks)) {
    pk_mat <- if (!is.null(peaks)) peaks else object
    if (is.null(precursor_mz) || !is.numeric(precursor_mz) || is.na(precursor_mz)) {
      stop("`precursor_mz` must be provided for raw peak data.", call. = FALSE)
    }
    pk_df <- as.data.frame(pk_mat)
    if (ncol(pk_df) < 2L) {
      stop("Peak table must contain at least two columns: m/z and intensity.", call. = FALSE)
    }
    mz_col <- pk_df[[1L]]
    int_col <- pk_df[[2L]]
    clean_pks <- lapply(seq_along(mz_col), function(i) {
      list(as.numeric(mz_col[i]), as.numeric(int_col[i]))
    })
    spectra_list <- list(list(
      id = "Spectrum_0001",
      precursor_mz = as.numeric(precursor_mz),
      peaks = clean_pks
    ))
  } else if (is.list(object)) {
    spectra_list <- object
  } else {
    stop("`object` must be an NTSProject, a peak matrix/data.frame, or a list of spectra.", call. = FALSE)
  }

  if (isTRUE(verbose)) {
    message(sprintf("Extracting 256-d MS2 semantic embeddings for %d spectra (device: %s)...",
                    length(spectra_list), device))
  }

  # 4. Exchange data via temporary JSON files
  tmp_in <- tempfile(pattern = "joint_embed_in_", fileext = ".json")
  tmp_out <- tempfile(pattern = "joint_embed_out_", fileext = ".json")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)

  jsonlite::write_json(spectra_list, tmp_in, auto_unbox = TRUE, digits = 8)

  cmd_args <- c(
    shQuote(engine_script),
    "--mode", "embed",
    "--input_json", shQuote(tmp_in),
    "--output_json", shQuote(tmp_out),
    "--weights_path", shQuote(paths$weights),
    "--config_path", shQuote(paths$config),
    "--device", device
  )

  run_res <- suppressWarnings(system2(py_exe, args = cmd_args, stdout = TRUE, stderr = TRUE))
  status_code <- attr(run_res, "status") %||% 0L

  if (status_code != 0L || !file.exists(tmp_out)) {
    err_msg <- paste(run_res, collapse = "\n")
    stop(
      "Embedding extraction failed (exit code ", status_code, "):\n",
      err_msg,
      call. = FALSE
    )
  }

  res_json <- jsonlite::fromJSON(tmp_out, simplifyVector = FALSE)

  # 5. Build embedding matrix
  fids <- vapply(res_json, function(it) as.character(it$id %||% ""), character(1))
  pmz <- vapply(res_json, function(it) as.numeric(it$precursor_mz %||% NA_real_), numeric(1))

  mat_rows <- lapply(res_json, function(it) as.numeric(unlist(it$embedding)))
  emb_mat <- do.call(rbind, mat_rows)
  if (is.null(emb_mat) || !is.matrix(emb_mat)) {
    emb_mat <- matrix(numeric(0), nrow = 0, ncol = 256)
  }
  rownames(emb_mat) <- fids
  colnames(emb_mat) <- paste0("dim_", seq_len(ncol(emb_mat)))

  if (isTRUE(as_matrix)) {
    return(emb_mat)
  }

  res_df <- data.frame(
    feature_id = fids,
    precursor_mz = pmz,
    stringsAsFactors = FALSE
  )
  res_df$embedding <- emb_mat
  rownames(res_df) <- NULL
  res_df
}


#' Match MS2 spectra against a reference library using semantic embedding cosine similarity
#'
#' Computes the semantic similarity between query spectrum embeddings (256-d)
#' and reference library embeddings using vectorized cosine similarity.
#' Because embeddings are L2-normalized unit vectors, cosine similarity equals
#' their inner product.
#'
#' @param query Query embeddings from [extract_ms2_embeddings], a numeric matrix
#'   (\eqn{M \times 256}), or a data.frame containing an `embedding` matrix column.
#' @param reference Reference embeddings matrix (\eqn{N \times 256}) or data.frame.
#' @param top_k Integer. Number of top semantic hits to return per query. Default `5L`.
#' @param min_similarity Numeric. Minimum cosine similarity threshold (0.0 to 1.0). Default `0.0`.
#' @param verbose Logical.
#'
#' @return A `data.frame` of top semantic matches:
#'   \itemize{
#'     \item `query_id`: Identifier of the query spectrum.
#'     \item `rank`: 1-based semantic similarity rank (1..top_k).
#'     \item `reference_id`: Identifier of the library reference compound.
#'     \item `semantic_similarity`: Cosine similarity score in range [-1.0, 1.0] (typically 0.0 to 1.0).
#'   }
#' @export
match_spectral_semantic <- function(query,
                                   reference,
                                   top_k = 5L,
                                   min_similarity = 0.0,
                                   verbose = TRUE) {
  q_mat <- .nts_as_embedding_matrix(query)
  r_mat <- .nts_as_embedding_matrix(reference)

  if (is.null(q_mat) || is.null(r_mat)) {
    stop("Both `query` and `reference` must provide 256-dimensional embedding vectors.", call. = FALSE)
  }
  if (ncol(q_mat) != 256L || ncol(r_mat) != 256L) {
    stop("Embeddings must have exactly 256 dimensions.", call. = FALSE)
  }

  q_ids <- rownames(q_mat)
  if (is.null(q_ids) || !length(q_ids)) q_ids <- paste0("Query_", seq_len(nrow(q_mat)))

  r_ids <- rownames(r_mat)
  if (is.null(r_ids) || !length(r_ids)) r_ids <- paste0("Ref_", seq_len(nrow(r_mat)))

  # Vectorized matrix product for cosine similarity
  sim_mat <- q_mat %*% t(r_mat)

  top_k <- min(max(as.integer(top_k), 1L), ncol(sim_mat))
  rows <- vector("list", nrow(sim_mat))

  for (i in seq_len(nrow(sim_mat))) {
    sims <- sim_mat[i, ]
    ord <- order(sims, decreasing = TRUE)
    sel_idx <- ord[seq_len(top_k)]
    sel_sims <- sims[sel_idx]
    keep <- sel_sims >= min_similarity

    if (any(keep)) {
      rows[[i]] <- data.frame(
        query_id = q_ids[i],
        rank = seq_len(sum(keep)),
        reference_id = r_ids[sel_idx[keep]],
        semantic_similarity = round(as.numeric(sel_sims[keep]), 4),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out) || !nrow(out)) {
    out <- data.frame(
      query_id = character(0),
      rank = integer(0),
      reference_id = character(0),
      semantic_similarity = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  rownames(out) <- NULL
  out
}


# =============================================================================
# Candidate Fingerprint Re-ranking (DeepMass-style Candidate Rescoring)
# =============================================================================

#' Tanimoto similarity between two bit index vectors
#' @keywords internal
#' @noRd
.nts_tanimoto_bits <- function(bits_pred, bits_cand) {
  if (!length(bits_pred) || !length(bits_cand)) {
    return(NA_real_)
  }
  n_inter <- length(intersect(bits_pred, bits_cand))
  n_union <- length(union(bits_pred, bits_cand))
  if (n_union == 0L) {
    return(0.0)
  }
  round(n_inter / n_union, 4)
}

#' Convert object to standard embedding matrix
#' @keywords internal
#' @noRd
.nts_as_embedding_matrix <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.matrix(x)) {
    return(x)
  }
  if (is.data.frame(x)) {
    if ("embedding" %in% names(x)) {
      emb_col <- x$embedding
      if (is.matrix(emb_col) && ncol(emb_col) == 256L) {
        if (!is.null(x$feature_id)) rownames(emb_col) <- as.character(x$feature_id)
        else if (!is.null(x$name)) rownames(emb_col) <- as.character(x$name)
        return(emb_col)
      } else if (is.list(emb_col)) {
        mat <- do.call(rbind, emb_col)
        if (is.matrix(mat) && ncol(mat) == 256L) {
          if (!is.null(x$feature_id)) rownames(mat) <- as.character(x$feature_id)
          else if (!is.null(x$name)) rownames(mat) <- as.character(x$name)
          return(mat)
        }
      }
    }
    dim_cols <- grep("^(dim_|emb_|V)[0-9]+$", names(x), value = TRUE)
    if (length(dim_cols) == 256L) {
      mat <- as.matrix(x[, dim_cols, drop = FALSE])
      if (!is.null(x$feature_id)) rownames(mat) <- as.character(x$feature_id)
      else if (!is.null(x$name)) rownames(mat) <- as.character(x$name)
      return(mat)
    }
  }
  NULL
}


#' Compute binary 2215-bit joint fingerprints from candidate SMILES
#'
#' Encodes candidate chemical structures into 2215-bit binary on-bit lists
#' (ECFP2048 + MACCS167 offset by 2048).
#'
#' @param smiles Character vector of SMILES strings.
#' @param python_path Optional Python path with RDKit installed.
#' @param verbose Logical.
#' @return A named list of integer vectors of active on-bit indices (0..2214).
#' @export
compute_candidate_joint_fp <- function(smiles,
                                      python_path = NULL,
                                      verbose = FALSE) {
  if (is.null(smiles) || !length(smiles)) return(list())
  smiles_chr <- as.character(smiles)
  valid_mask <- !is.na(smiles_chr) & nzchar(trimws(smiles_chr))
  unique_smiles <- unique(smiles_chr[valid_mask])

  if (!length(unique_smiles)) {
    res_empty <- vector("list", length(smiles_chr))
    if (!is.null(names(smiles))) {
      names(res_empty) <- names(smiles)
    } else {
      names(res_empty) <- smiles_chr
    }
    return(res_empty)
  }

  fp_dict <- list()

  # 1. Try Python RDKit via joint_fp_engine.py
  py_exe <- .nts_resolve_python_executable(python_path)
  engine_script <- .nts_resolve_joint_fp_engine_script()

  if (!is.null(py_exe) && file.exists(py_exe) && file.exists(engine_script)) {
    tmp_in <- tempfile(pattern = "joint_smiles_in_", fileext = ".json")
    tmp_out <- tempfile(pattern = "joint_smiles_out_", fileext = ".json")
    on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)

    input_list <- lapply(seq_along(unique_smiles), function(idx) {
      list(id = paste0("SMI_", idx), smiles = unique_smiles[idx])
    })
    jsonlite::write_json(input_list, tmp_in, auto_unbox = TRUE)

    cmd_args <- c(
      shQuote(engine_script),
      "--mode", "smiles_fp",
      "--input_json", shQuote(tmp_in),
      "--output_json", shQuote(tmp_out)
    )
    run_res <- suppressWarnings(system2(py_exe, args = cmd_args, stdout = TRUE, stderr = TRUE))
    status_code <- attr(run_res, "status") %||% 0L

    if (status_code == 0L && file.exists(tmp_out)) {
      res_json <- tryCatch(jsonlite::fromJSON(tmp_out, simplifyVector = FALSE), error = function(e) NULL)
      if (is.list(res_json) && length(res_json)) {
        for (item in res_json) {
          smi <- item$smiles
          bits <- as.integer(unlist(item$joint_fp_on_bits))
          if (length(bits)) {
            fp_dict[[smi]] <- bits
          }
        }
      }
    }
  }

  # 2. Try rcdk fallback for any SMILES not computed
  missing_smi <- setdiff(unique_smiles, names(fp_dict))
  if (length(missing_smi) && requireNamespace("rcdk", quietly = TRUE)) {
    for (smi in missing_smi) {
      mol <- tryCatch(rcdk::parse.smiles(smi)[[1L]], error = function(e) NULL)
      if (!is.null(mol)) {
        ecfp <- tryCatch(rcdk::get.fingerprint(mol, type = "circular"), error = function(e) NULL)
        maccs <- tryCatch(rcdk::get.fingerprint(mol, type = "maccs"), error = function(e) NULL)
        ecfp_bits <- if (!is.null(ecfp)) (as.integer(ecfp@bits) - 1L) else integer(0)
        maccs_bits <- if (!is.null(maccs)) (as.integer(maccs@bits) - 1L + 2048L) else integer(0)
        all_bits <- sort(unique(c(ecfp_bits, maccs_bits)))
        if (length(all_bits)) {
          fp_dict[[smi]] <- all_bits
        }
      }
    }
  }

  # Map back to original vector
  out <- lapply(smiles_chr, function(s) {
    if (!is.na(s) && s %in% names(fp_dict)) fp_dict[[s]] else integer(0)
  })
  if (!is.null(names(smiles))) {
    names(out) <- names(smiles)
  } else {
    names(out) <- smiles_chr
  }
  out
}


#' Rescore and re-rank annotation candidates using predicted MS2 joint fingerprints
#'
#' Evaluates the structural congruence between Transformer-predicted MS2
#' molecular fingerprints (2215-bit) and candidate chemical structures (SMILES).
#' Calculates Tanimoto similarity on the 2215-bit space and computes a blended
#' composite score to re-rank candidate hits, substantially reducing false
#' positives among isobaric and isomeric candidates.
#'
#' @param candidates A data.frame of candidates containing at least a feature
#'   column (`feature`, `name`, or `feature_id`) and a structure column (`SMILES`).
#' @param predicted_fp Predicted fingerprints from [predict_joint_fingerprint].
#'   Can be an [NTSProject] (reads from `@physchem_proxy`), a data.frame with
#'   columns `feature_id` and `joint_fp_bits` (or `MS2_Joint_FP_Bits`), or a
#'   named list of active on-bit vectors.
#' @param smiles_col Column name for SMILES in `candidates`. Default `"SMILES"`.
#' @param weight_fp Weight of fingerprint similarity in composite score (0 to 1). Default `0.30`.
#' @param score_col Existing score column to blend with. Default `"score_total"`.
#' @param feature_col Optional feature identifier column in `candidates`. If `NULL`,
#'   automatically resolved from `feature`, `name`, or `feature_id`.
#' @param python_path Optional path to Python executable for RDKit fingerprinting.
#' @param verbose Logical.
#'
#' @return A `data.frame` containing the original candidate columns plus:
#'   \itemize{
#'     \item `fp_similarity`: Tanimoto similarity between predicted and candidate fingerprints (0.0 to 1.0).
#'     \item `score_reranked`: Blended composite score.
#'     \item `rank_reranked`: 1-based rank within each feature after rescoring.
#'   }
#'   Candidates are sorted within each feature by `score_reranked` descending.
#' @export
rescore_candidates_by_joint_fp <- function(candidates,
                                          predicted_fp = NULL,
                                          smiles_col = "SMILES",
                                          weight_fp = 0.30,
                                          score_col = "score_total",
                                          feature_col = NULL,
                                          python_path = NULL,
                                          verbose = FALSE) {
  if (is.null(candidates) || !is.data.frame(candidates) || !nrow(candidates)) {
    return(candidates)
  }

  # 1. Resolve feature column
  fcol <- feature_col
  if (is.null(fcol)) {
    for (cand_name in c("feature", "name", "feature_id")) {
      if (cand_name %in% names(candidates)) {
        fcol <- cand_name
        break
      }
    }
  }
  if (is.null(fcol)) {
    stop("`candidates` must contain a feature column (`feature`, `name`, or `feature_id`).", call. = FALSE)
  }

  # 2. Resolve predicted fingerprints
  pred_bits_map <- list()
  if (!is.null(predicted_fp)) {
    if (inherits(predicted_fp, "NTSProject")) {
      proxy <- as.data.frame(get_physchem_proxy(predicted_fp), stringsAsFactors = FALSE)
      if (nrow(proxy) && "name" %in% names(proxy)) {
        fp_col_hit <- intersect(names(proxy), c("MS2_Joint_FP_Bits", "joint_fp_bits"))[1L]
        if (!is.na(fp_col_hit)) {
          for (k in seq_len(nrow(proxy))) {
            val <- proxy[[fp_col_hit]][k]
            if (!is.na(val) && nzchar(trimws(val))) {
              pred_bits_map[[as.character(proxy$name[k])]] <- as.integer(strsplit(as.character(val), ";")[[1L]])
            }
          }
        }
      }
    } else if (is.data.frame(predicted_fp)) {
      id_hit <- intersect(names(predicted_fp), c("feature_id", "name", "feature", "id"))[1L]
      fp_hit <- intersect(names(predicted_fp), c("joint_fp_bits", "MS2_Joint_FP_Bits", "joint_fp_on_bits"))[1L]
      if (!is.na(id_hit) && !is.na(fp_hit)) {
        for (k in seq_len(nrow(predicted_fp))) {
          val <- predicted_fp[[fp_hit]][[k]]
          if (is.numeric(val)) {
            pred_bits_map[[as.character(predicted_fp[[id_hit]][k])]] <- as.integer(val)
          } else if (is.character(val) && nzchar(trimws(val))) {
            pred_bits_map[[as.character(predicted_fp[[id_hit]][k])]] <- as.integer(strsplit(val, ";")[[1L]])
          }
        }
      }
    } else if (is.list(predicted_fp)) {
      for (nm in names(predicted_fp)) {
        pred_bits_map[[nm]] <- as.integer(predicted_fp[[nm]])
      }
    }
  }

  # 3. Candidate fingerprints and similarities
  n_rows <- nrow(candidates)
  fp_sims <- if ("fp_similarity" %in% names(candidates)) as.numeric(candidates$fp_similarity) else rep(NA_real_, n_rows)

  cand_bits_col <- intersect(names(candidates), c("joint_fp_bits", "candidate_fp", "fp_bits"))[1L]
  need_sim <- is.na(fp_sims)

  if (any(need_sim) && length(pred_bits_map)) {
    # If candidate bits are already stored in candidates column
    if (!is.na(cand_bits_col)) {
      for (i in which(need_sim)) {
        fid <- as.character(candidates[[fcol]][i])
        if (!is.na(fid) && fid %in% names(pred_bits_map)) {
          raw_bits <- candidates[[cand_bits_col]][[i]]
          c_bits <- if (is.numeric(raw_bits)) as.integer(raw_bits)
                    else if (is.character(raw_bits) && nzchar(trimws(raw_bits))) as.integer(strsplit(raw_bits, ";")[[1L]])
                    else integer(0)
          if (length(c_bits)) {
            p_bits <- pred_bits_map[[fid]]
            if (length(p_bits)) {
              fp_sims[i] <- .nts_tanimoto_bits(p_bits, c_bits)
            }
          }
        }
      }
      need_sim <- is.na(fp_sims)
    }

    # If still need sim and SMILES is available, compute via RDKit / rcdk
    has_smi <- smiles_col %in% names(candidates)
    if (any(need_sim) && has_smi) {
      smi_vec <- candidates[[smiles_col]][need_sim]
      cand_fp_list <- compute_candidate_joint_fp(smi_vec, python_path = python_path, verbose = verbose)
      if (length(cand_fp_list)) {
        for (i in which(need_sim)) {
          fid <- as.character(candidates[[fcol]][i])
          smi <- as.character(candidates[[smiles_col]][i])
          if (!is.na(fid) && fid %in% names(pred_bits_map) && !is.na(smi) && smi %in% names(cand_fp_list)) {
            p_bits <- pred_bits_map[[fid]]
            c_bits <- cand_fp_list[[smi]]
            if (length(p_bits) && length(c_bits)) {
              fp_sims[i] <- .nts_tanimoto_bits(p_bits, c_bits)
            }
          }
        }
      }
    }
  }

  candidates$fp_similarity <- fp_sims

  # 5. Blend score_reranked
  orig_score <- if (score_col %in% names(candidates)) as.numeric(candidates[[score_col]]) else rep(NA_real_, n_rows)
  score_reranked <- rep(NA_real_, n_rows)

  weight_fp <- min(max(as.numeric(weight_fp), 0.0), 1.0)
  for (i in seq_len(n_rows)) {
    os <- orig_score[i]
    fs <- fp_sims[i]
    if (!is.na(os) && !is.na(fs)) {
      score_reranked[i] <- (1 - weight_fp) * os + weight_fp * fs
    } else if (!is.na(fs)) {
      score_reranked[i] <- fs
    } else {
      score_reranked[i] <- os
    }
  }
  candidates$score_reranked <- round(score_reranked, 4)

  # 6. Re-order within each feature
  f_vec <- as.character(candidates[[fcol]])
  reord <- order(f_vec, -ifelse(is.na(candidates$score_reranked), -999, candidates$score_reranked),
                 -ifelse(is.na(candidates$fp_similarity), -999, candidates$fp_similarity),
                 na.last = TRUE)
  candidates <- candidates[reord, , drop = FALSE]

  # Re-compute rank_reranked within each feature
  f_sorted <- as.character(candidates[[fcol]])
  candidates$rank_reranked <- as.integer(seq_len(nrow(candidates)) - match(f_sorted, f_sorted) + 1L)

  rownames(candidates) <- NULL
  candidates
}



# =============================================================================
# Internal Helpers
# =============================================================================

#' Resolve candidate Python executable
#' @keywords internal
#' @noRd
.nts_resolve_python_executable <- function(python_path = NULL) {
  if (!is.null(python_path)) {
    if (nzchar(python_path) && file.exists(python_path)) {
      return(normalizePath(python_path, mustWork = FALSE))
    } else {
      # Explicit user path was supplied but does not exist
      return(NULL)
    }
  }

  env_py <- Sys.getenv("NTSMINER_PYTHON", "")
  if (nzchar(env_py) && file.exists(env_py)) {
    return(normalizePath(env_py, mustWork = FALSE))
  }

  ret_py <- Sys.getenv("RETICULATE_PYTHON", "")
  if (nzchar(ret_py) && file.exists(ret_py)) {
    return(normalizePath(ret_py, mustWork = FALSE))
  }

  # If reticulate is available and initialized
  if (requireNamespace("reticulate", quietly = TRUE)) {
    cfg <- tryCatch(reticulate::py_discover_config(), error = function(e) NULL)
    if (!is.null(cfg) && !is.null(cfg$python) && file.exists(cfg$python)) {
      ok <- tryCatch({
        system2(cfg$python, c("-c", "import torch"), stdout = FALSE, stderr = FALSE) == 0L
      }, error = function(e) FALSE)
      if (isTRUE(ok)) return(normalizePath(cfg$python, mustWork = FALSE))
    }
  }

  # Fallback to system PATH candidates only if torch is importable
  for (cmd in c("python", "python3")) {
    sys_p <- Sys.which(cmd)
    if (nzchar(sys_p) && file.exists(sys_p)) {
      ok <- tryCatch({
        system2(sys_p, c("-c", "import torch"), stdout = FALSE, stderr = FALSE) == 0L
      }, error = function(e) FALSE)
      if (isTRUE(ok)) return(normalizePath(sys_p, mustWork = FALSE))
    }
  }

  NULL
}

#' Resolve joint fingerprint model assets
#' @keywords internal
#' @noRd
.nts_resolve_joint_fp_model_paths <- function(weights = NULL, config = NULL, tox = NULL) {
  # 1. Weights
  w_cand <- if (!is.null(weights)) weights else c(
    .nts_extdata(c(
      file.path("models", "joint_fp", "joint_fp_model_weights_fp16.pt"),
      file.path("models", "joint_fp", "joint_fp_model_weights.pt")
    )),
    file.path("inst", "extdata", "models", "joint_fp", "joint_fp_model_weights_fp16.pt"),
    file.path("inst", "extdata", "models", "joint_fp", "joint_fp_model_weights.pt")
  )
  w_hit <- w_cand[!is.na(w_cand) & nzchar(w_cand) & file.exists(w_cand)]
  weights_path <- if (length(w_hit)) normalizePath(w_hit[1L], mustWork = FALSE) else NULL

  # 2. Config
  c_cand <- if (!is.null(config)) config else c(
    .nts_extdata(file.path("models", "joint_fp", "joint_fp_model_config.json")),
    file.path("inst", "extdata", "models", "joint_fp", "joint_fp_model_config.json")
  )
  c_hit <- c_cand[!is.na(c_cand) & nzchar(c_cand) & file.exists(c_cand)]
  config_path <- if (length(c_hit)) normalizePath(c_hit[1L], mustWork = FALSE) else NULL

  # 3. Tox model
  t_cand <- if (!is.null(tox)) tox else c(
    .nts_extdata(file.path("models", "joint_fp", "tox_model.joblib")),
    file.path("inst", "extdata", "models", "joint_fp", "tox_model.joblib")
  )
  t_hit <- t_cand[!is.na(t_cand) & nzchar(t_cand) & file.exists(t_cand)]
  tox_path <- if (length(t_hit)) normalizePath(t_hit[1L], mustWork = FALSE) else NULL

  list(
    weights = weights_path,
    config = config_path,
    tox = tox_path
  )
}

#' Resolve joint fingerprint Python engine script
#' @keywords internal
#' @noRd
.nts_resolve_joint_fp_engine_script <- function() {
  cand <- c(
    system.file("python", "joint_fp_engine.py", package = "NTSminer"),
    file.path("inst", "python", "joint_fp_engine.py"),
    file.path("python", "joint_fp_engine.py")
  )
  hit <- cand[nzchar(cand) & file.exists(cand)]
  if (length(hit)) {
    return(normalizePath(hit[1L], mustWork = FALSE))
  }
  # Fallback to the known default relative path
  normalizePath(file.path("inst", "python", "joint_fp_engine.py"), mustWork = FALSE)
}

#' Extract best MS2 spectra from an NTSProject
#' @keywords internal
#' @noRd
.nts_extract_ms2_from_project <- function(project) {
  ft <- get_feature_table(project)
  if (!nrow(ft)) return(list())

  id_col <- if ("name" %in% names(ft)) "name" else if ("feature_id" %in% names(ft)) "feature_id" else NULL
  if (is.null(id_col)) return(list())

  # Check clean_ms2 or ms2data slots
  ms2_source <- NULL
  if (.hasSlot(project, "clean_ms2") && length(project@clean_ms2)) {
    ms2_source <- project@clean_ms2
  } else if (.hasSlot(project, "ms2data") && length(project@ms2data)) {
    ms2_source <- project@ms2data
  }

  if (is.null(ms2_source) || !length(ms2_source)) {
    return(list())
  }

  out <- list()
  for (i in seq_len(nrow(ft))) {
    fid <- as.character(ft[[id_col]][i])
    pmz <- as.numeric(ft$mz[i] %||% NA_real_)
    if (is.na(pmz) || pmz <= 0) next

    spec_entry <- ms2_source[[fid]]
    if (is.null(spec_entry)) next

    # Polymorphic peak extraction
    pks <- list()
    if (is.list(spec_entry) && !is.null(spec_entry$spectrum)) {
      sp <- spec_entry$spectrum
      if ((is.data.frame(sp) || is.matrix(sp)) && nrow(sp) > 0L) {
        m_vec <- as.numeric(if ("mz" %in% names(sp)) sp$mz else sp[, 1L])
        i_vec <- as.numeric(if ("intensity" %in% names(sp)) sp$intensity else sp[, 2L])
        ok <- !is.na(m_vec) & !is.na(i_vec) & i_vec > 0
        pks <- lapply(which(ok), function(r) list(m_vec[r], i_vec[r]))
      }
    } else if (is.matrix(spec_entry) || is.data.frame(spec_entry)) {
      mat <- as.matrix(spec_entry)
      if (ncol(mat) >= 2L && nrow(mat) > 0L) {
        m_vec <- as.numeric(mat[, 1L])
        i_vec <- as.numeric(mat[, 2L])
        ok <- !is.na(m_vec) & !is.na(i_vec) & i_vec > 0
        pks <- lapply(which(ok), function(r) list(m_vec[r], i_vec[r]))
      }
    } else if (is.list(spec_entry) && !is.null(spec_entry$mz) && !is.null(spec_entry$intensity)) {
      m_vec <- as.numeric(spec_entry$mz)
      i_vec <- as.numeric(spec_entry$intensity)
      ok <- !is.na(m_vec) & !is.na(i_vec) & i_vec > 0
      pks <- lapply(which(ok), function(r) list(m_vec[r], i_vec[r]))
    } else if (is.list(spec_entry) && length(spec_entry) > 0L && is.list(spec_entry[[1L]]) && !is.null(spec_entry[[1L]]$mz)) {
      # ms2data: list of scans
      best_scan <- spec_entry[[1L]]
      m_vec <- as.numeric(best_scan$mz)
      i_vec <- as.numeric(best_scan$intensity)
      ok <- !is.na(m_vec) & !is.na(i_vec) & i_vec > 0
      pks <- lapply(which(ok), function(r) list(m_vec[r], i_vec[r]))
    }

    if (length(pks)) {
      out[[length(out) + 1L]] <- list(
        id = fid,
        precursor_mz = pmz,
        peaks = pks
      )
    }
  }

  out
}

#' Write joint fingerprint and toxicity results into NTSProject
#' @keywords internal
#' @noRd
.nts_write_joint_fp_to_project <- function(project, df_results) {
  if (!nrow(df_results)) return(project)

  proxy <- as.data.frame(get_physchem_proxy(project), stringsAsFactors = FALSE)
  if (!nrow(proxy) || !"name" %in% names(proxy)) {
    ft <- get_feature_table(project)
    proxy <- data.frame(name = as.character(ft$name), stringsAsFactors = FALSE)
  }

  # Match on name / feature_id
  m <- match(as.character(proxy$name), as.character(df_results$feature_id))

  proxy$MS2_Toxicity_Score <- df_results$toxicity_score[m]
  proxy$MS2_Toxicity_Risk <- df_results$toxicity_risk[m]
  proxy$MS2_Functional_Groups <- df_results$detected_functional_groups[m]
  proxy$MS2_Joint_FP_Bits <- df_results$joint_fp_bits[m]
  proxy$MS2_N_ECFP_Bits <- df_results$n_ecfp_bits[m]
  proxy$MS2_N_MACCS_Bits <- df_results$n_maccs_bits[m]

  project <- set_physchem_proxy(project, proxy)

  # Append log entry
  if (.hasSlot(project, "log")) {
    project@log <- .nts_log_append(
      project@log,
      "predict_joint_fingerprint",
      list(
        n_spectra = nrow(df_results),
        n_high_risk = sum(df_results$toxicity_risk == "High", na.rm = TRUE),
        n_moderate_risk = sum(df_results$toxicity_risk == "Moderate", na.rm = TRUE)
      )
    )
  }

  project
}
