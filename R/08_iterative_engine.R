## Iterative local inference helpers (decoupled radar + on-disk mini libraries)
##
## Does not modify build_structural_library or annotate_by_dual_network.
## Fingerprinting reuses `.nts_normalize_smiles_string` and `.nts_rcdk_stop`.
## Targets / queries: one SMILES at a time (`parse.smiles` + `get.fingerprint`)
## inside `tryCatch`, optional SMILES length cap, progress per row, JVM/R GC every
## 1000 fingerprint steps (`.nts_force_jvm_r_gc`). Warmup parse of "C" before Step 2.
## Parsing uses `.nts_mol_first_from_smiles_rcdk` (full SMILES then dot fragments for salts).
## Dense `target_matrix` fill (no `rbind`).
## Similarity: `calculate_similarity_batch`-style Dice via `Q %*% t(T)` + `outer()`.


#' @noRd
.nts_force_jvm_r_gc <- function() {
  if (requireNamespace("rJava", quietly = TRUE)) {
    tryCatch(
      {
        sys <- rJava::J("java/lang/System")
        rJava::.jcall(sys, "V", "gc")
      },
      error = function(e) invisible(NULL)
    )
  }
  gc(verbose = FALSE)
}


#' @noRd
.nts_java_err_is_fatal <- function(cond) {
  msg <- tryCatch(conditionMessage(cond), error = function(e) "")
  grepl(
    "OutOfMemoryError|Java heap space|InternalError|StackOverflowError|Could not create the Java Virtual Machine|java\\.lang\\.ExceptionInInitializerError",
    msg,
    ignore.case = TRUE
  )
}


#' Load precompiled target fingerprint bundle from RDS (see `search_asymmetric_smiles`).
#' @noRd
.nts_load_precompiled_target_bundle <- function(path, smiles_col, id_col, fp_matrix_nbit) {
  obj <- readRDS(path.expand(as.character(path)[1L]))
  if (!is.list(obj)) {
    stop(
      "search_asymmetric_smiles: precompiled RDS root must be a list with target_matrix and metadata.",
      call. = FALSE
    )
  }
  req <- c("target_matrix", "metadata")
  miss <- req[!req %in% names(obj)]
  if (length(miss)) {
    stop(
      sprintf("search_asymmetric_smiles: RDS must contain: %s.", paste(req, collapse = ", ")),
      call. = FALSE
    )
  }
  tm <- obj$target_matrix
  if (!is.matrix(tm)) {
    tm <- as.matrix(tm)
  }
  storage.mode(tm) <- "double"
  md <- obj$metadata
  df <- as.data.frame(md, stringsAsFactors = FALSE)
  n <- nrow(df)
  if (nrow(tm) != n) {
    stop(
      "search_asymmetric_smiles: nrow(metadata) must equal nrow(target_matrix).",
      call. = FALSE
    )
  }
  nm <- names(df)
  if (!smiles_col %in% nm) {
    stop(
      sprintf(
        "search_asymmetric_smiles: smiles_col=%s not found in RDS metadata. Columns (first 40): %s",
        smiles_col,
        paste(utils::head(nm, 40L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!id_col %in% nm) {
    stop(
      sprintf(
        "search_asymmetric_smiles: id_col=%s not found in RDS metadata. Columns (first 40): %s",
        id_col,
        paste(utils::head(nm, 40L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  L <- ncol(tm)
  if (!is.finite(L) || L < 1L) {
    stop("search_asymmetric_smiles: target_matrix must have at least one column.", call. = FALSE)
  }
  L_user <- suppressWarnings(as.integer(fp_matrix_nbit)[1L])
  if (!is.null(fp_matrix_nbit) && is.finite(L_user) && L_user > 0L && L_user != L) {
    stop(
      sprintf(
        "search_asymmetric_smiles: fp_matrix_nbit (%d) does not match precompiled matrix columns (%d).",
        L_user,
        L
      ),
      call. = FALSE
    )
  }
  ts <- obj$target_sums
  if (is.null(ts)) {
    ts <- rowSums(tm)
  } else {
    ts <- as.numeric(ts)
    if (length(ts) != n) {
      stop(
        "search_asymmetric_smiles: target_sums length must match nrow(target_matrix).",
        call. = FALSE
      )
    }
  }
  if (!is.null(obj$norm_target_smi)) {
    nts <- as.character(obj$norm_target_smi)
  } else {
    nts <- as.character(df[[smiles_col]])
  }
  if (length(nts) != n) {
    stop(
      "search_asymmetric_smiles: norm_target_smi (if present) must have length nrow(metadata).",
      call. = FALSE
    )
  }
  list(
    df = df,
    n = n,
    norm_target_smi = nts,
    target_matrix = tm,
    target_sums = ts,
    L = as.integer(L)
  )
}


#' Asymmetric structural similarity radar (M queries x N CSV targets)
#'
#' Loads a compound list with [data.table::fread] and scores **each query SMILES**
#' against **every CSV row** (M-by-N), without any N-by-N library self-comparison.
#' Target rows: normalized SMILES are fingerprinted **one at a time** with
#' `tryCatch` (skip on error; skip if `nchar(SMILES) > 500`). Molecules are resolved
#' with `.nts_mol_first_from_smiles_rcdk` (full string then dot-separated fragments),
#' matching [build_structural_library]. On-bits written into a
#' pre-allocated dense matrix (default width **307** for `substructure`, else
#' `fp_matrix_nbit` or `1024`). JVM/R GC runs every 1000 fingerprint steps. Queries use the
#' same per-string pattern. All-vs-all intersections use `Q %*% t(target_matrix)` with vectorized
#' Dice (or Tanimoto), numerically matching [fingerprint::fp.sim.matrix] on the
#' same CDK fingerprints without O(M x N) JNI pair calls. [build_structural_library]
#' still recomputes fingerprints for the mini-RDS.
#'
#' @details
#' **Precompiled RDS** must be a `list` with at least `target_matrix` (numeric
#' matrix `n` x `L`), `metadata` (`data.frame` with `n` rows, including `smiles_col`
#' and `id_col`), and optionally `target_sums` (length `n`) and `norm_target_smi`
#' (length `n`; otherwise `metadata[[smiles_col]]` is used for output SMILES).
#'
#' @param query_smiles Character vector of query SMILES (distinct queries are
#'   fingerprinted once each).
#' @param full_csv_path Path to the full structural CSV. Not read when
#'   `precompiled_rds_path` points to an existing RDS bundle.
#' @param precompiled_rds_path Optional path to an RDS file (see Details). When set
#'   and the file exists, **target** normalization and CDK fingerprinting are
#'   skipped; only query seeds are parsed in Java.
#' @param smiles_col,id_col Column names in the CSV or in RDS `metadata` (same role as
#'   [build_structural_library] `smiles_col` / `id_col` on **input** tables).
#'   Defaults follow OECD-style headers; override with e.g. `"SMILES"` and `"ID"`.
#' @param top_n Per-query hit cap before union and deduplication by row.
#' @param fingerprint_type Passed to [rcdk::get.fingerprint] `type`.
#' @param fallback_smiles_col Optional fallback SMILES column (same as [build_structural_library]).
#' @param similarity_method `"dice"` or `"tanimoto"`.
#' @param fp_matrix_nbit Fixed bit length for the dense engine (pad 0 / truncate
#'   on-bits above this length). Default `NULL`: use **307** for
#'   `fingerprint_type = "substructure"` (CDK layout), otherwise `1024`.
#' @param chunk_size Kept for API compatibility (evaluated but not used for batching);
#'   targets and queries are always fingerprinted **one SMILES at a time** with
#'   `tryCatch`. Default `1`.
#'
#' @return A `data.table` of unique hit rows: all original CSV columns for
#'   selected rows, plus canonical `SMILES` and `ID` (from `id_col`) for
#'   [build_structural_library].
#' @export
search_asymmetric_smiles <- function(query_smiles,
                                     full_csv_path = NULL,
                                     smiles_col = "QSAR_READY_SMILES",
                                     id_col = "DTXSID",
                                     top_n = 50L,
                                     fingerprint_type = "circular",
                                     fallback_smiles_col = NULL,
                                     similarity_method = "dice",
                                     fp_matrix_nbit = NULL,
                                     chunk_size = 1L,
                                     precompiled_rds_path = NULL) {
  .nts_rcdk_stop()
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }

  invisible(force(chunk_size))

  fp_type <- as.character(fingerprint_type)[1L]
  if (!nzchar(fp_type)) {
    stop("fingerprint_type must be a non-empty string", call. = FALSE)
  }

  pc_raw <- precompiled_rds_path
  pc_path <- if (!is.null(pc_raw)) {
    path.expand(as.character(pc_raw)[1L])
  } else {
    ""
  }
  use_rds <- nzchar(pc_path) && file.exists(pc_path)

  fp_type_l <- tolower(trimws(fp_type))
  fp_type_rcdk <- if (identical(fp_type_l, "substructure")) {
    "substructure"
  } else {
    fp_type
  }

  queries <- unique(trimws(as.character(query_smiles)))
  queries <- queries[!is.na(queries) & nzchar(queries)]
  if (!length(queries)) {
    return(data.table::data.table(ID = character(0), SMILES = character(0)))
  }
  n_q <- length(queries)

  if (!use_rds) {
    if (is.null(full_csv_path) || !nzchar(trimws(as.character(full_csv_path)[1L]))) {
      stop(
        "search_asymmetric_smiles: full_csv_path is required unless precompiled_rds_path is a path to an existing RDS file.",
        call. = FALSE
      )
    }
    full_csv_path <- path.expand(as.character(full_csv_path)[1L])
    .nts_check(full_csv_path, "file_exists", .arg = "full_csv_path")
    message(sprintf("  -> Loading CSV and mapping columns (%s)...", basename(full_csv_path)))
    dt <- data.table::fread(full_csv_path, na.strings = c("", "NA", "N/A"))
    if (!nrow(dt)) {
      return(data.table::data.table(ID = character(0), SMILES = character(0)))
    }
  }

  jvm_gc_every <- 1000L

  ## Warning-suppression state for the target-fingerprint loop below. It MUST be
  ## defined here (function level, before the loop): the tryCatch error handler
  ## inside that loop touches `warn_tar$n`, so a single unparsable target SMILES
  ## would otherwise raise "object 'warn_tar' not found" *inside the handler*
  ## and abort the whole search.
  warn_tar <- new.env(parent = emptyenv())
  warn_tar$n <- 0L
  warn_tar$max_detail <- 15L

  if (!use_rds) {
  nm_dt <- names(dt)
  if (!smiles_col %in% nm_dt) {
    stop(
      sprintf(
        "search_asymmetric_smiles: smiles_col=%s not found in CSV. Available columns (first 40): %s",
        smiles_col,
        paste(utils::head(nm_dt, 40L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!id_col %in% nm_dt) {
    stop(
      sprintf(
        "search_asymmetric_smiles: id_col=%s not found in CSV. Available columns (first 40): %s",
        id_col,
        paste(utils::head(nm_dt, 40L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!is.null(fallback_smiles_col)) {
    fb <- as.character(fallback_smiles_col)[1L]
    if (!nzchar(fb) || !fb %in% nm_dt) {
      stop(
        sprintf(
          "search_asymmetric_smiles: fallback_smiles_col=%s not found in CSV.",
          fallback_smiles_col
        ),
        call. = FALSE
      )
    }
    if (identical(fb, smiles_col)) {
      fallback_smiles_col <- NULL
    }
  }

  df <- as.data.frame(dt, stringsAsFactors = FALSE)
  id_vec <- as.character(df[[id_col]])
  smi_vec <- as.character(df[[smiles_col]])
  n <- length(id_vec)
  fb_vec <- NULL
  if (!is.null(fallback_smiles_col)) {
    fb_vec <- as.character(df[[fallback_smiles_col]])
    if (length(fb_vec) != n) {
      stop("Internal error: fallback_smiles_col length mismatch.", call. = FALSE)
    }
  }

  message("  -> Step 1: Normalizing targets...")
  norm_target_smi <- rep(NA_character_, n)
  use_pb_norm <- interactive() && n >= 3L
  if (use_pb_norm) {
    pb_norm <- utils::txtProgressBar(min = 0L, max = n, style = 3L)
  }
  for (j in seq_len(n)) {
    smi <- .nts_normalize_smiles_string(smi_vec[j])
    if (is.na(smi) && !is.null(fb_vec)) {
      smi <- .nts_normalize_smiles_string(fb_vec[j])
    }
    if (!is.na(smi)) {
      norm_target_smi[j] <- smi
    }
    if (use_pb_norm) {
      utils::setTxtProgressBar(pb_norm, j)
    } else if (n >= 500L && j %% 500L == 0L) {
      message(".", appendLF = FALSE)
    }
  }
  if (use_pb_norm) {
    close(pb_norm)
    cat("\n")
  } else if (n >= 500L) {
    message("")
  }

  L <- suppressWarnings(as.integer(fp_matrix_nbit)[1L])
  if (is.null(fp_matrix_nbit) || !is.finite(L) || L < 32L) {
    L <- if (identical(fp_type_l, "substructure")) {
      307L
    } else {
      1024L
    }
  }

  ix <- which(!is.na(norm_target_smi))
  if (!length(ix)) {
    warning(
      "search_asymmetric_smiles: no parsable target SMILES after normalization; returning empty table.",
      call. = FALSE,
      immediate. = TRUE
    )
    return(data.table::data.table(ID = character(0), SMILES = character(0)))
  }

  target_smiles_vector <- as.character(norm_target_smi[ix])
  n_targets <- length(target_smiles_vector)

  message("  -> Step 2: Extracting fingerprints (Individual mode with tryCatch)...")
  message(sprintf(
    "  -> (fingerprint type %s, matrix %d x %d; then %d query SMILES)",
    fp_type_rcdk,
    n,
    L,
    n_q
  ))
  target_matrix <- matrix(0L, nrow = n, ncol = L)

  n_fp_total <- n_targets + n_q
  use_pb_fp <- interactive() && n_fp_total >= 3L
  if (use_pb_fp) {
    pb_fp <- utils::txtProgressBar(min = 0L, max = n_fp_total, style = 3L)
  }

  for (k in seq_len(n_targets)) {
    smi <- target_smiles_vector[k]
    smi1 <- as.character(smi)[1L]
    nch <- suppressWarnings(nchar(smi1, allowNA = TRUE, keepNA = TRUE))
    if (length(nch) != 1L || is.na(nch)) {
      nch <- 0L
    }
    if (isTRUE(nch > 500L)) {
      warning(
        sprintf("CSV row %d: SMILES too long (>500 chars), skipping.", ix[k]),
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      tryCatch(
        {
          mol <- .nts_mol_first_from_smiles_rcdk(as.character(smi))
          if (!is.null(mol)) {
            fp <- rcdk::get.fingerprint(mol, type = fp_type_rcdk)
            bits <- methods::slot(fp, "bits")
            on_bits <- as.integer(bits)
            on_bits <- on_bits[is.finite(on_bits) & on_bits >= 1L & on_bits <= L]
            if (length(on_bits)) {
              target_matrix[ix[k], on_bits] <- 1L
            }
          }
          TRUE
        },
        error = function(e) {
          if (.nts_java_err_is_fatal(e)) {
            stop(
              "search_asymmetric_smiles: fatal JVM error while fingerprinting targets; restart R.",
              call. = FALSE
            )
          }
          warn_tar$n <- warn_tar$n + 1L
          if (warn_tar$n <= warn_tar$max_detail) {
            warning(
              sprintf("Row %d: CDK/Java error, skipping.", ix[k]),
              call. = FALSE,
              immediate. = TRUE
            )
          } else if (warn_tar$n == warn_tar$max_detail + 1L) {
            warning(
              "search_asymmetric_smiles: further per-row CDK/Java warnings suppressed (check salts / encoding).",
              call. = FALSE,
              immediate. = TRUE
            )
          }
          FALSE
        }
      )
    }
    if (use_pb_fp) {
      utils::setTxtProgressBar(pb_fp, k)
    } else if (n_fp_total >= 500L && k %% 500L == 0L) {
      message(".", appendLF = FALSE)
    }
    if (k %% jvm_gc_every == 0L) {
      .nts_force_jvm_r_gc()
    }
  }

  } else {
    message("  -> Loaded precompiled target matrix from RDS. Bypassing Target Java parsing.")
    bundle <- .nts_load_precompiled_target_bundle(pc_path, smiles_col, id_col, fp_matrix_nbit)
    df <- bundle$df
    n <- bundle$n
    norm_target_smi <- bundle$norm_target_smi
    target_matrix <- bundle$target_matrix
    target_sums <- as.numeric(bundle$target_sums)
    L <- as.integer(bundle$L)
    n_targets <- n
    n_fp_total <- n_targets + n_q
    use_pb_fp <- interactive() && n_fp_total >= 3L
    if (use_pb_fp) {
      pb_fp <- utils::txtProgressBar(min = 0L, max = n_fp_total, style = 3L)
    }
  }

  query_matrix <- matrix(0L, nrow = n_q, ncol = L)

  message("  -> Warming up CDK / JVM (parse SMILES 'C' and fingerprint)...")
  warm_ok <- tryCatch(
    {
      pm0 <- rcdk::parse.smiles("C")
      if (!is.list(pm0)) {
        pm0 <- list(pm0)
      }
      mol0 <- if (length(pm0) >= 1L) pm0[[1L]] else NULL
      if (is.null(mol0)) {
        stop("parse.smiles('C') returned no molecule", call. = FALSE)
      }
      fp0 <- rcdk::get.fingerprint(mol0, type = fp_type_rcdk)
      invisible(methods::slot(fp0, "nbit"))
      rm(pm0, mol0, fp0)
      TRUE
    },
    error = function(e) {
      if (.nts_java_err_is_fatal(e)) {
        stop(
          "search_asymmetric_smiles: fatal JVM/CDK error during warmup; restart R and try again.",
          call. = FALSE
        )
      }
      FALSE
    }
  )
  if (!isTRUE(warm_ok)) {
    stop(
      "search_asymmetric_smiles: CDK warmup failed (could not parse or fingerprint SMILES 'C'). Restart R, reload rcdk, then retry.",
      call. = FALSE
    )
  }

  warn_q <- new.env(parent = emptyenv())
  warn_q$n <- 0L
  warn_q$max_detail <- 15L

  for (i in seq_len(n_q)) {
    qsmi <- queries[i]
    q1 <- as.character(qsmi)[1L]
    nchq <- suppressWarnings(nchar(q1, allowNA = TRUE, keepNA = TRUE))
    if (length(nchq) != 1L || is.na(nchq)) {
      nchq <- 0L
    }
    if (isTRUE(nchq > 500L)) {
      warning(
        sprintf("Query %d: SMILES too long (>500 chars), skipping.", i),
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      tryCatch(
        {
          mol <- .nts_mol_first_from_smiles_rcdk(as.character(qsmi))
          if (!is.null(mol)) {
            fp <- rcdk::get.fingerprint(mol, type = fp_type_rcdk)
            bits <- methods::slot(fp, "bits")
            on_bits <- as.integer(bits)
            on_bits <- on_bits[is.finite(on_bits) & on_bits >= 1L & on_bits <= L]
            if (length(on_bits)) {
              query_matrix[i, on_bits] <- 1L
            }
          }
          TRUE
        },
        error = function(e) {
          if (.nts_java_err_is_fatal(e)) {
            stop(
              "search_asymmetric_smiles: fatal JVM error while fingerprinting queries; restart R.",
              call. = FALSE
            )
          }
          warn_q$n <- warn_q$n + 1L
          if (warn_q$n <= warn_q$max_detail) {
            warning(
              sprintf("Query row %d: CDK/Java error, skipping.", i),
              call. = FALSE,
              immediate. = TRUE
            )
          } else if (warn_q$n == warn_q$max_detail + 1L) {
            warning(
              "search_asymmetric_smiles: further query CDK/Java warnings suppressed.",
              call. = FALSE,
              immediate. = TRUE
            )
          }
          FALSE
        }
      )
    }
    if (use_pb_fp) {
      utils::setTxtProgressBar(pb_fp, n_targets + i)
    } else if (n_fp_total >= 500L && (n_targets + i) %% 500L == 0L) {
      message(".", appendLF = FALSE)
    }
    if (i %% jvm_gc_every == 0L) {
      .nts_force_jvm_r_gc()
    }
  }

  if (!use_rds && warn_tar$n > warn_tar$max_detail) {
    message(
      sprintf(
        "  -> Note: %d target row(s) hit CDK/Java errors (only first %d shown above).",
        warn_tar$n,
        warn_tar$max_detail
      )
    )
  }

  if (use_pb_fp) {
    close(pb_fp)
    cat("\n")
  } else if (n_fp_total >= 500L) {
    message("")
  }
  .nts_force_jvm_r_gc()
  if (!any(rowSums(target_matrix) > 0)) {
    warning(
      "search_asymmetric_smiles: no usable target fingerprints after extraction; returning empty table.",
      call. = FALSE,
      immediate. = TRUE
    )
    return(data.table::data.table(ID = character(0), SMILES = character(0)))
  }
  storage.mode(target_matrix) <- "double"
  if (!use_rds) {
    target_sums <- rowSums(target_matrix)
  }

  sim_meth <- tolower(as.character(similarity_method)[1L])
  tn <- suppressWarnings(as.integer(top_n)[1L])
  if (!is.finite(tn) || tn < 1L) {
    tn <- 50L
  }

  q_ok <- which(rowSums(query_matrix) > 0L)
  if (!length(q_ok)) {
    warning(
      "search_asymmetric_smiles: no valid query fingerprints; returning empty table.",
      call. = FALSE,
      immediate. = TRUE
    )
    return(data.table::data.table(ID = character(0), SMILES = character(0)))
  }

  Q <- query_matrix[q_ok, , drop = FALSE]
  storage.mode(Q) <- "double"
  M <- nrow(Q)

  n_qraw <- length(query_smiles)
  if (!is.finite(n_qraw) || n_qraw < 1L) {
    n_qraw <- 1L
  }

  message("  -> Step 3: Fast matrix calculation...")
  pick_rows <- integer(0)
  use_pb_q <- interactive() && M > 0L
  if (use_pb_q) {
    pb_q <- utils::txtProgressBar(min = 0L, max = n_qraw, style = 3L)
    on.exit(
      {
        close(pb_q)
        cat("\n")
      },
      add = TRUE
    )
  }

  sim_mat <- .nts_sim_matrix_from_bits(Q, target_matrix, method = sim_meth)

  for (i in seq_len(M)) {
    rw <- sim_mat[i, ]
    ord <- order(rw, decreasing = TRUE, na.last = TRUE)
    ord <- ord[rw[ord] > 0]
    if (length(ord)) {
      take <- ord[seq_len(min(tn, length(ord)))]
      pick_rows <- c(pick_rows, take)
    }
    if (use_pb_q) {
      utils::setTxtProgressBar(pb_q, max(1L, min(n_qraw, ceiling(i * n_qraw / M))))
    }
  }
  if (use_pb_q) {
    utils::setTxtProgressBar(pb_q, n_qraw)
  }
  pick_rows <- sort(unique(pick_rows))
  if (!length(pick_rows)) {
    return(data.table::data.table(ID = character(0), SMILES = character(0)))
  }

  out <- data.table::as.data.table(df[pick_rows, , drop = FALSE])
  out[, SMILES := norm_target_smi[pick_rows]]
  out[, ID := as.character(out[[id_col]])]
  out[]
}


#' Build a on-disk mini structural library from seed SMILES (radar + CSV + RDS)
#'
#' Runs [search_asymmetric_smiles], writes **`prefix_candidates.csv`** under
#' `output_dir`, then calls [build_structural_library] to produce
#' **`prefix_db.rds`** in the same folder. No hidden tempfiles; paths are stable
#' for auditing and for passing to [annotate_by_dual_network].
#'
#' @param initial_seeds_smiles Character vector of seed SMILES (queries for the radar).
#' @param full_csv_path Full compound CSV path. Not read when `precompiled_rds_path`
#'   points to an existing precompiled target RDS (see [search_asymmetric_smiles]).
#' @param precompiled_rds_path Optional; passed to [search_asymmetric_smiles] for
#'   instant target matrix loading (bypasses target-side Java fingerprinting).
#' @param output_dir Directory for artifacts (created if missing), default
#'   `"data/iterative_dbs"`.
#' @param prefix File name prefix (e.g. `"gen_1"`, `"gen_2"`).
#' @param top_n Per-seed radar width (see [search_asymmetric_smiles]).
#' @param fingerprint_type Fingerprinter `type` for radar and [build_structural_library].
#' @param smiles_col,id_col Passed to [search_asymmetric_smiles] (defaults match
#'   OECD lists; use `"SMILES"` / `"ID"` for generic CSVs).
#' @param fallback_smiles_col Optional fallback SMILES column on the CSV.
#' @param similarity_method `"dice"` or `"tanimoto"` for the radar matrix.
#' @param fp_matrix_nbit Passed to [search_asymmetric_smiles] (fixed bit width for
#'   the pure-R similarity engine).
#' @param chunk_size Passed to [search_asymmetric_smiles] (API compatibility).
#'
#' @details
#' Exactly one of `full_csv_path` (existing file) or `precompiled_rds_path`
#' (existing RDS bundle) must be usable; if the RDS path exists, the CSV is ignored.
#'
#' @return Character path to the written **`prefix_db.rds`** (normalized via
#'   [normalizePath] when the file exists).
#' @export
build_local_db_from_seeds <- function(initial_seeds_smiles,
                                      full_csv_path = NULL,
                                      output_dir = "data/iterative_dbs",
                                      prefix = "gen_1",
                                      top_n = 50L,
                                      fingerprint_type = "circular",
                                      smiles_col = "QSAR_READY_SMILES",
                                      id_col = "DTXSID",
                                      fallback_smiles_col = NULL,
                                      similarity_method = "dice",
                                      fp_matrix_nbit = NULL,
                                      chunk_size = 1L,
                                      precompiled_rds_path = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }

  message(sprintf("========== [build_local_db_from_seeds: %s] ==========", as.character(prefix)[1L]))

  pc_raw <- precompiled_rds_path
  pc_path <- if (!is.null(pc_raw)) {
    path.expand(as.character(pc_raw)[1L])
  } else {
    ""
  }
  use_rds <- nzchar(pc_path) && file.exists(pc_path)
  if (!use_rds) {
    if (is.null(full_csv_path) || !nzchar(trimws(as.character(full_csv_path)[1L]))) {
      stop(
        "build_local_db_from_seeds: full_csv_path is required unless precompiled_rds_path is a path to an existing RDS file.",
        call. = FALSE
      )
    }
    full_csv_path <- path.expand(as.character(full_csv_path)[1L])
    .nts_check(full_csv_path, "file_exists", .arg = "full_csv_path")
  }

  od <- path.expand(as.character(output_dir)[1L])
  if (!nzchar(od)) {
    stop("output_dir must be a non-empty path", call. = FALSE)
  }
  if (!dir.exists(od)) {
    dir.create(od, recursive = TRUE, showWarnings = FALSE)
  }

  px <- as.character(prefix)[1L]
  if (!nzchar(trimws(px))) {
    stop("prefix must be non-empty", call. = FALSE)
  }

  seeds_chr <- unique(trimws(as.character(initial_seeds_smiles)))
  seeds_chr <- seeds_chr[!is.na(seeds_chr) & nzchar(seeds_chr)]
  n_seed <- length(seeds_chr)
  message(sprintf("  [1/3] Running asymmetric radar (%d seed SMILES)...", n_seed))

  cand <- search_asymmetric_smiles(
    query_smiles = initial_seeds_smiles,
    full_csv_path = if (use_rds) NULL else full_csv_path,
    smiles_col = smiles_col,
    id_col = id_col,
    top_n = top_n,
    fingerprint_type = fingerprint_type,
    fallback_smiles_col = fallback_smiles_col,
    similarity_method = similarity_method,
    fp_matrix_nbit = fp_matrix_nbit,
    chunk_size = chunk_size,
    precompiled_rds_path = if (use_rds) pc_path else NULL
  )

  if (!nrow(cand)) {
    stop(
      "build_local_db_from_seeds: radar returned no candidates; nothing written.",
      call. = FALSE
    )
  }

  message(sprintf("  Radar done: %d candidate row(s). Writing CSV...", nrow(cand)))

  csv_path <- file.path(od, paste0(px, "_candidates.csv"))
  rds_path <- file.path(od, paste0(px, "_db.rds"))

  data.table::fwrite(cand, csv_path)

  message("  [2/3] Calling build_structural_library() to compile RDS...")

  build_structural_library(
    csv_path = csv_path,
    smiles_col = "SMILES",
    id_col = "ID",
    fallback_smiles_col = NULL,
    fingerprint_type = fingerprint_type,
    output_rds = rds_path
  )

  if (file.exists(rds_path)) {
    out_path <- normalizePath(rds_path, winslash = "/", mustWork = FALSE)
  } else {
    out_path <- rds_path
  }
  message(sprintf("  [3/3] Local database ready: %s", out_path))
  out_path
}


## Back-compat name (same as search_asymmetric_smiles); not exported.
.nts_search_asymmetric_smiles <- search_asymmetric_smiles
