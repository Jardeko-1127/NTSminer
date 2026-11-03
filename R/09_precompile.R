## Precompile large target libraries into dense fingerprint matrices (RDS).
## One-shot batch job: CSV -> rcdk fingerprints -> pure-R artifact for downstream
## asymmetric search without re-parsing targets in Java.


#' @noRd
.precompile_jvm_gc <- function() {
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


#' Precompile a target compound table into a dense fingerprint matrix (RDS)
#'
#' Reads a large structural CSV, fingerprints every surviving row with
#' [rcdk::get.fingerprint], and writes a compact RDS bundle for fast pure-R
#' loading in iterative workflows. Intended as a **one-shot** batch step to
#' avoid repeated Java target parsing.
#'
#' @param input_csv Path to the input CSV file.
#' @param output_rds Path to the output `.rds` file (parent directories are
#'   created if missing).
#' @param smi_col Column name for SMILES strings (default matches OECD-style
#'   headers).
#' @param id_col Column name for stable compound identifiers.
#' @param max_smi_len Maximum allowed SMILES length (characters); longer rows are
#'   dropped before fingerprinting.
#' @param fp_type Fingerprinter `type` passed to [rcdk::get.fingerprint].
#' @param fp_matrix_nbit Number of columns in the dense bit matrix; on-bits are
#'   kept only when `1 <= bit <= fp_matrix_nbit`.
#' @param chunk_size Number of rows parsed per `rcdk::parse.smiles` batch; each
#'   chunk ends with explicit JVM/R garbage collection.
#' @param timeout_sec Per-molecule fingerprint timeout (seconds).
#'
#' @details
#' The saved object is a `list` with elements:
#' \describe{
#'   \item{`target_matrix`}{Integer matrix `n` x `fp_matrix_nbit` (0/1).}
#'   \item{`target_sums`}{Numeric vector `rowSums(target_matrix)` (bit counts per row).}
#'   \item{`smiles`}{Character vector of SMILES used after filtering (same row
#'     order as `target_matrix`).}
#'   \item{`ids`}{Character vector of identifiers from `id_col`.}
#'   \item{`metadata`}{The filtered input table as a `data.frame` (all columns).}
#' }
#' This layout matches [search_asymmetric_smiles] `precompiled_rds_path`
#' expectations (`target_matrix`, `target_sums`, `metadata`; optional
#' `norm_target_smi` can be added separately if needed).
#'
#' @return Invisibly, the expanded `output_rds` path as a character scalar.
#' @export
precompile_target_db <- function(input_csv,
                                 output_rds,
                                 smi_col = "QSAR_READY_SMILES",
                                 id_col = "DTXSID",
                                 max_smi_len = 400L,
                                 fp_type = "circular",
                                 fp_matrix_nbit = 1024L,
                                 chunk_size = 100L,
                                 timeout_sec = 30L) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("precompile_target_db: package 'data.table' is required.", call. = FALSE)
  }
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("precompile_target_db: package 'rcdk' is required.", call. = FALSE)
  }

  input_csv <- path.expand(as.character(input_csv)[1L])
  output_rds <- path.expand(as.character(output_rds)[1L])
  .nts_check(input_csv, "file_exists", .arg = "input_csv")

  max_smi_len <- suppressWarnings(as.integer(max_smi_len)[1L])
  if (!is.finite(max_smi_len) || max_smi_len < 1L) {
    stop("precompile_target_db: max_smi_len must be a positive integer.", call. = FALSE)
  }

  fp_matrix_nbit <- suppressWarnings(as.integer(fp_matrix_nbit)[1L])
  if (!is.finite(fp_matrix_nbit) || fp_matrix_nbit < 32L) {
    stop("precompile_target_db: fp_matrix_nbit must be an integer >= 32.", call. = FALSE)
  }

  chunk_size <- max(1L, suppressWarnings(as.integer(chunk_size)[1L]))
  if (!is.finite(chunk_size)) {
    chunk_size <- 100L
  }

  fp_type_use <- as.character(fp_type)[1L]
  if (!nzchar(fp_type_use)) {
    stop("precompile_target_db: fp_type must be a non-empty string.", call. = FALSE)
  }

  message("-> [Precompile] Loading raw CSV...")
  dt <- data.table::fread(input_csv, na.strings = c("", "NA", "N/A"))
  if (!nrow(dt)) {
    stop("precompile_target_db: input CSV has zero rows.", call. = FALSE)
  }

  nm <- names(dt)
  if (!smi_col %in% nm) {
    stop(
      sprintf(
        "precompile_target_db: smiles column '%s' not found. First columns: %s",
        smi_col,
        paste(utils::head(nm, 30L), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!id_col %in% nm) {
    stop(
      sprintf(
        "precompile_target_db: id column '%s' not found. First columns: %s",
        id_col,
        paste(utils::head(nm, 30L), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  smi_raw <- trimws(as.character(dt[[smi_col]]))
  ok <- !is.na(smi_raw) & nzchar(smi_raw)
  nch <- suppressWarnings(nchar(smi_raw, type = "chars", allowNA = TRUE))
  ok <- ok & !is.na(nch) & nch <= max_smi_len
  dt <- dt[ok]
  n <- nrow(dt)
  if (!n) {
    stop(
      "precompile_target_db: no rows left after removing NA/empty/overlong SMILES.",
      call. = FALSE
    )
  }

  message(sprintf("-> [Precompile] %d row(s) after SMILES filter (max length %d).", n, max_smi_len))

  target_smiles <- trimws(as.character(dt[[smi_col]]))
  ids_chr <- as.character(dt[[id_col]])
  meta_df <- as.data.frame(dt, stringsAsFactors = FALSE)

  target_matrix <- matrix(0L, nrow = n, ncol = fp_matrix_nbit)

  message("-> [Precompile] Starting JVM (warmup parse + fingerprint)...")
  tryCatch(
    {
      pm0 <- rcdk::parse.smiles("C")
      if (!is.list(pm0)) {
        pm0 <- list(pm0)
      }
      mol0 <- if (length(pm0) >= 1L) pm0[[1L]] else NULL
      if (is.null(mol0)) {
        stop("parse.smiles('C') returned no molecule", call. = FALSE)
      }
      fp0 <- rcdk::get.fingerprint(mol0, type = fp_type_use)
      invisible(methods::slot(fp0, "nbit"))
      rm(pm0, mol0, fp0)
      TRUE
    },
    error = function(e) {
      stop(
        sprintf(
          "precompile_target_db: JVM/CDK warmup failed: %s",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
  .precompile_jvm_gc()

  message(sprintf(
    "-> [Precompile] Fingerprinting %d row(s) in chunks of %d (type '%s', %d bits)...",
    n,
    chunk_size,
    fp_type_use,
    fp_matrix_nbit
  ))

  use_pb <- interactive() && n >= 2L
  if (use_pb) {
    pb <- utils::txtProgressBar(min = 0L, max = n, style = 3L)
    on.exit(
      {
        suppressWarnings(try(close(pb), silent = TRUE))
        cat("\n")
      },
      add = TRUE
    )
  }

  chunk_starts <- seq.int(1L, n, by = chunk_size)
  for (cs in chunk_starts) {
    ce <- min(n, cs + chunk_size - 1L)
    idx <- cs:ce
    smi_chunk <- target_smiles[idx]

    mols <- tryCatch(
      rcdk::parse.smiles(smi_chunk),
      error = function(e) {
        warning(
          sprintf(
            "[Precompile] parse.smiles failed for rows [%d, %d]: %s - skipping chunk.",
            cs,
            ce,
            conditionMessage(e)
          ),
          call. = FALSE,
          immediate. = TRUE
        )
        NULL
      }
    )

    if (is.null(mols)) {
      if (use_pb) {
        utils::setTxtProgressBar(pb, ce)
      }
      .precompile_jvm_gc()
      next
    }

    if (!is.list(mols)) {
      mols <- list(mols)
    }

    n_mol <- length(mols)
    n_idx <- length(idx)
    if (n_mol != n_idx) {
      warning(
        sprintf(
          "[Precompile] parse.smiles length mismatch (mols=%d, rows=%d) for chunk [%d,%d]; aligning by minimum length.",
          n_mol,
          n_idx,
          cs,
          ce
        ),
        call. = FALSE,
        immediate. = TRUE
      )
      n_take <- min(n_mol, n_idx)
    } else {
      n_take <- n_idx
    }

    for (j in seq_len(n_take)) {
      row_i <- idx[j]
      mol <- mols[[j]]
      if (is.null(mol)) {
        if (use_pb) {
          utils::setTxtProgressBar(pb, row_i)
        }
        next
      }

      tryCatch(
        {
          .do_fp <- function() {
            fp <- rcdk::get.fingerprint(mol, type = fp_type_use)
            bits <- as.integer(methods::slot(fp, "bits"))
            bits <- bits[is.finite(bits) & bits >= 1L & bits <= fp_matrix_nbit]
            if (length(bits)) {
              target_matrix[row_i, bits] <<- 1L
            }
            TRUE
          }
          ts_val <- suppressWarnings(as.numeric(timeout_sec)[1L])
          if (is.finite(ts_val) && ts_val > 0 && requireNamespace("R.utils", quietly = TRUE)) {
            R.utils::withTimeout(.do_fp(), timeout = ts_val, onTimeout = "error")
          } else {
            .do_fp()
          }
          invisible(TRUE)
        },
        error = function(e) {
          warning(
            sprintf(
              "[Precompile] get.fingerprint failed at row %d: %s",
              row_i,
              conditionMessage(e)
            ),
            call. = FALSE,
            immediate. = TRUE
          )
          invisible(FALSE)
        }
      )

      if (use_pb) {
        utils::setTxtProgressBar(pb, row_i)
      }
    }

    rm(mols)
    .precompile_jvm_gc()
  }

  message("-> [Precompile] Computing row sums...")
  target_sums <- rowSums(target_matrix)
  storage.mode(target_matrix) <- "integer"

  out_dir <- dirname(output_rds)
  if (!nzchar(out_dir)) {
    out_dir <- "."
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  precompiled_data <- list(
    target_matrix = target_matrix,
    target_sums = target_sums,
    smiles = target_smiles,
    ids = ids_chr,
    metadata = meta_df
  )

  message(sprintf("-> [Precompile] Saving RDS to %s ...", output_rds))
  saveRDS(precompiled_data, file = output_rds, compress = "xz")
  .precompile_jvm_gc()

  message("-> [Precompile] Done.")
  invisible(output_rds)
}
