# RT prediction module ==========================================================
#
# Multi-fingerprint computation, RF/XGBoost modeling, prediction, persistence,
# and transfer learning for retention time prediction in NTS workflows.

# Fingerprint constants --------------------------------------------------------
.FP_BITS <- c(circular = 1024L, standard = 1024L, extended = 1024L,
              substructure = 307L, maccs = 166L, pubchem = 881L,
              kr = 4860L, shortestpath = 1024L)
.FP_BASE_TYPES <- c("standard", "extended", "maccs", "pubchem", "circular")
.FP_FULL_TYPES <- names(.FP_BITS)

# Format version written into newly built RT models. Bump it when the meaning of
# an existing field changes or a field becomes required. "0.1.0" is the legacy
# format still accepted by load_rt_model(): it predates `app_domain`,
# `split_strategy` / `n_groups` and `calib_shift`, all of which may be NULL there.
.NTS_RT_MODEL_VERSION <- "0.2.0"
.NTS_RT_MODEL_LEGACY_VERSIONS <- "0.1.0"


#' Resolve fingerprint types from shorthand
#' @keywords internal
#' @noRd
.nts_resolve_fp_types <- function(types) {
  types <- match.arg(types, c("base", "full", names(.FP_BITS)), several.ok = TRUE)
  out <- character(0)
  for (t in types) {
    if (t == "base") out <- union(out, .FP_BASE_TYPES)
    else if (t == "full") out <- union(out, .FP_FULL_TYPES)
    else out <- c(out, t)
  }
  unique(out)
}


# Phase 1: compute_fingerprints ------------------------------------------------

#' Compute molecular fingerprints for SMILES
#'
#' Computes multiple fingerprint types for a set of SMILES strings using rcdk.
#' Supports on-bit index output (sparse) and full 0/1 matrix output.
#' Parallel computation via BiocParallel.
#'
#' @param smiles Character vector of SMILES strings.
#' @param types Fingerprint types. Shorthand `"base"` (`standard` + `extended`
#'   + `circular` + `maccs` + `pubchem`, 4119 bits) or `"full"` (all 8 types,
#'   10310 bits), or individual: `"circular"`, `"standard"`, `"extended"`,
#'   `"substructure"`, `"maccs"`, `"pubchem"`, `"kr"`, `"shortestpath"`.
#' @param fp_mode `"onbits"` returns a data.table with list columns of on-bit
#'   indices; `"matrix"` returns a 0/1 numeric matrix.
#' @param BPPARAM A BiocParallelParam object for parallelisation.
#'
#' @return A data.table (onbits) or numeric matrix (matrix).
#' @export
compute_fingerprints <- function(smiles,
                                 types = c("base", "full", "circular",
                                           "substructure", "maccs", "pubchem",
                                           "kr", "shortestpath"),
                                 fp_mode = c("onbits", "matrix"),
                                 BPPARAM = BiocParallel::SerialParam()) {
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("rcdk is required. Install: install.packages('rcdk')", call. = FALSE)
  }
  fp_mode <- match.arg(fp_mode)
  fp_types <- .nts_resolve_fp_types(types)
  if (!length(fp_types)) {
    stop("No valid fingerprint types specified.", call. = FALSE)
  }

  smi <- as.character(smiles)
  smi <- trimws(smi)
  valid_mask <- !is.na(smi) & nzchar(smi)
  if (!any(valid_mask)) {
    stop("No valid SMILES strings provided.", call. = FALSE)
  }
  smi_clean <- smi[valid_mask]
  n <- length(smi_clean)

  nbits <- .FP_BITS[fp_types]
  total_bits <- sum(nbits)

  message(sprintf("Computing %s fingerprints for %d SMILES (%d bits) ...",
                  paste(fp_types, collapse = "+"), n, total_bits))

  fp_list <- BiocParallel::bplapply(seq_len(n), function(i) {
    s <- smi_clean[i]
    mol <- tryCatch(rcdk::parse.smiles(s)[[1]], error = function(e) NULL)
    if (is.null(mol)) {
      return(structure(list(), failed = TRUE, smiles = s))
    }
    fpl <- lapply(fp_types, function(tp) {
      tryCatch(rcdk::get.fingerprint(mol, type = tp), error = function(e) NULL)
    })
    names(fpl) <- fp_types
    fpl$failed <- FALSE
    fpl$smiles <- s
    fpl
  }, BPPARAM = BPPARAM)

  all_failed <- all(vapply(fp_list, function(x) isTRUE(x$failed), logical(1L)))
  if (all_failed) {
    stop("All SMILES failed to parse with rcdk.", call. = FALSE)
  }
  n_failed <- sum(vapply(fp_list, function(x) isTRUE(x$failed), logical(1L)))
  if (n_failed > 0) {
    warning(n_failed, " SMILES failed to parse and will have empty fingerprints.", call. = FALSE)
  }

  if (fp_mode == "onbits") {
    dt <- data.table::data.table(name = smi_clean)
    for (tp in fp_types) {
      col <- paste0("fp_", tp)
      dt[[col]] <- lapply(fp_list, function(x) {
        if (isTRUE(x$failed)) return(integer(0))
        fp <- x[[tp]]
        if (is.null(fp)) return(integer(0))
        bits <- as.integer(fp@bits)
        if (!length(bits)) return(integer(0))
        bits[bits > 0]
      })
    }
    data.table::setattr(dt, "fp_types", fp_types)
    data.table::setattr(dt, "fp_nbits", as.integer(nbits))
    return(dt)
  }

  mat <- matrix(0, nrow = n, ncol = total_bits)
  col_offset <- 0L
  fp_bounds <- list()
  for (tp in fp_types) {
    nb <- as.integer(nbits[tp])
    idx_range <- (col_offset + 1L):(col_offset + nb)
    fp_bounds[[tp]] <- c(idx_range[1L], idx_range[length(idx_range)])
    for (i in seq_len(n)) {
      if (isTRUE(fp_list[[i]]$failed)) next
      fp <- fp_list[[i]][[tp]]
      if (is.null(fp)) next
      bits <- as.integer(fp@bits)
      if (!length(bits)) next
      on <- bits[bits > 0 & bits <= nb]
      if (length(on)) mat[i, col_offset + on] <- 1
    }
    col_offset <- col_offset + nb
  }
  rownames(mat) <- smi_clean
  colnames(mat) <- unlist(lapply(fp_types, function(tp) {
    paste0(tp, "_", seq_len(nbits[tp]))
  }))
  attr(mat, "fp_types") <- fp_types
  attr(mat, "fp_bounds") <- fp_bounds
  mat
}


# Phase 2: build_rt_model ------------------------------------------------------

#' Build an RT prediction model
#'
#' Trains a random forest (ranger) or gradient boosting (xgboost) model for
#' retention time prediction from molecular fingerprints.
#'
#' @param smiles Character vector of training SMILES.
#' @param rt Numeric vector of retention times (seconds), same length as smiles.
#' @param fingerprint_types Fingerprint types (passed to [compute_fingerprints]).
#' @param method `"rf"` (ranger) or `"xgb"` (xgboost).
#' @param tune_params Named list of hyperparameter overrides.
#' @param cv_folds Number of cross-validation folds for tuning (small datasets only).
#' @param train_split Proportion of data for training (remainder = test).
#' @param group Optional grouping vector (e.g. chemical class or scaffold) the
#'   same length as the input compounds. When supplied, whole groups are held
#'   out together so structurally related analogues (homologues, isomers) do not
#'   leak between train and test. When `NULL` (default) groups are derived
#'   automatically by average-linkage clustering of the molecular fingerprints at
#'   `cluster_similarity`.
#' @param cluster_similarity Tanimoto similarity at which automatically derived
#'   groups are cut (default 0.8). Only used when `group = NULL`.
#' @param app_min_similarity Tanimoto similarity to the nearest training compound
#'   below which a prediction is flagged `"out_of_domain"` by [predict_rt()].
#' @param seed Random seed (also seeds the grouped split).
#' @param BPPARAM A BiocParallelParam for fingerprint computation.
#' @param fingerprint_matrix Optional pre-computed fingerprint matrix (rows = compounds).
#' @param descriptor_matrix Optional pre-computed descriptor matrix (rows = compounds).
#' @param verbose If TRUE, print progress messages.
#'
#' @return An `nts_rt_model` S3 object.
#' @export
build_rt_model <- function(smiles,
                           rt,
                           fingerprint_types = "base",
                           method = c("rf", "xgb"),
                           tune_params = list(),
                           cv_folds = 5,
                           train_split = 0.8,
                           group = NULL,
                           cluster_similarity = 0.8,
                           app_min_similarity = 0.6,
                           seed = 42,
                           BPPARAM = BiocParallel::SerialParam(),
                           fingerprint_matrix = NULL,
                           descriptor_matrix = NULL,
                           verbose = TRUE) {
  method <- match.arg(method)
  set.seed(seed)

  rt_vec <- suppressWarnings(as.numeric(rt))
  has_fp <- !is.null(fingerprint_matrix)
  has_desc <- !is.null(descriptor_matrix)

  if (has_fp || has_desc) {
    smi <- rep("precomputed", if (has_fp) nrow(fingerprint_matrix) else nrow(descriptor_matrix))
    n_input <- length(smi)
    .nts_rt_check_group(group, n_input, rt_vec)
    ok <- !is.na(rt_vec) & is.finite(rt_vec)
    if (!any(ok)) stop("No valid rt values.", call. = FALSE)
    rt_vec <- rt_vec[ok]
    fp_mat <- if (has_fp) fingerprint_matrix[ok, , drop = FALSE] else NULL
    desc_mat <- if (has_desc) descriptor_matrix[ok, , drop = FALSE] else NULL
    group <- if (is.null(group)) NULL else group[ok]
    n <- length(rt_vec)
    if (verbose) message("Using pre-computed matrices: ", n, " rows, ",
                         if (has_fp) ncol(fp_mat) else 0, " fp + ",
                         if (has_desc) ncol(desc_mat) else 0, " desc")
  } else {
    smi <- as.character(smiles)
    smi <- trimws(smi)
    if (length(smi) != length(rt_vec)) {
      stop("smiles and rt must have the same length.", call. = FALSE)
    }
    n_input <- length(smi)
    .nts_rt_check_group(group, n_input, rt_vec)
    ok <- !is.na(smi) & nzchar(trimws(smi)) & !is.na(rt_vec) & is.finite(rt_vec)
    if (!any(ok)) stop("No valid (smiles, rt) pairs.", call. = FALSE)
    if (verbose) message("Training with ", sum(ok), " valid compounds")
    smi <- trimws(smi[ok])
    rt_vec <- rt_vec[ok]
    group <- if (is.null(group)) NULL else group[ok]
    n <- length(smi)
    if (verbose) message("Computing fingerprints...")
    fp_mat <- compute_fingerprints(smi, types = fingerprint_types, fp_mode = "matrix", BPPARAM = BPPARAM)
    desc_mat <- NULL
  }

  # Combine fingerprints + descriptors
  if (!is.null(fp_mat) && !is.null(desc_mat)) {
    X_raw <- cbind(fp_mat, desc_mat)
  } else if (!is.null(fp_mat)) {
    X_raw <- fp_mat
  } else if (!is.null(desc_mat)) {
    X_raw <- desc_mat
  } else {
    stop("No features provided.", call. = FALSE)
  }

  fp_types <- if (!is.null(fp_mat)) attr(fp_mat, "fp_types") else character(0)
  fp_bounds <- if (!is.null(fp_mat)) attr(fp_mat, "fp_bounds") else list()
  nbits <- if (!is.null(fp_mat)) {
    vapply(fp_types, function(tp) {
      b <- fp_bounds[[tp]]; as.integer(b[2L] - b[1L] + 1L)
    }, integer(1L))
  } else integer(0)
  names(nbits) <- fp_types
  p_orig <- ncol(X_raw)
  desc_cols <- if (!is.null(desc_mat)) colnames(desc_mat) else character(0)

  # Group-aware train/test split (finding B5). A naive random split scatters
  # structural analogues across train and test, so test error is optimistic.
  # Groups are resolved (user-supplied, else fingerprint clustering) and whole
  # groups are assigned to one side. The split happens *before* the feature
  # transform is fitted, so the transform never sees the test rows (finding
  # P1-1: fitting on all rows leaked test information into the training matrix).
  grp <- .nts_rt_resolve_groups(group, n = n, fp = fp_mat,
                                cluster_similarity = cluster_similarity, seed = seed)
  split_obj <- .nts_rt_group_split(grp, n = n, train_split = train_split, seed = seed)
  train_idx <- split_obj$train_idx
  test_idx <- split_obj$test_idx
  if (verbose) message("Split strategy: ", split_obj$strategy,
                       " (", length(unique(grp)), " groups)",
                       " -> train: ", length(train_idx),
                       ", test: ", length(test_idx))

  if (verbose) message("Preprocessing features...")
  preproc <- .nts_rt_preprocess(X_raw, fit_idx = train_idx)
  fp_clean <- preproc$matrix
  p_final <- ncol(fp_clean)
  if (verbose) message("After preprocessing: ", p_final, " features (removed ", preproc$n_removed, ")")

  X_train <- fp_clean[train_idx, , drop = FALSE]
  y_train <- rt_vec[train_idx]
  X_test  <- fp_clean[test_idx, , drop = FALSE]
  y_test  <- rt_vec[test_idx]

  if (method == "rf") {
    fit <- .nts_rt_train_rf(X_train, y_train, tune_params, cv_folds, seed, verbose)
  } else {
    fit <- .nts_rt_train_xgb(X_train, y_train, tune_params, cv_folds, seed, verbose)
  }

  pred_test <- if (method == "rf") {
    as.numeric(predict(fit, X_test)$predictions)
  } else {
    as.numeric(predict(fit, xgboost::xgb.DMatrix(X_test)))
  }
  test_metrics <- .nts_rt_metrics(y_test, pred_test)
  if (verbose) message(sprintf("Test RMSE=%.1f, R2=%.3f, MAE=%.1f",
                                test_metrics$RMSE, test_metrics$R2, test_metrics$MAE))

  var_imp <- .nts_rt_var_importance(fit, method, fp_types, fp_bounds, colnames(fp_clean))

  model <- structure(list(
    model          = fit,
    method         = method,
    fp_types       = fp_types,
    fp_nbits       = nbits,
    nzv_removed    = preproc$nzv_removed,
    cor_removed    = preproc$cor_removed,
    nzv_idx        = preproc$nzv_idx,
    cor_idx        = preproc$cor_idx,
    preproc_center = preproc$center,
    preproc_scale  = preproc$scale,
    kept_colnames  = colnames(fp_clean),
    descriptor_cols = desc_cols,
    test_metrics   = test_metrics,
    var_importance = var_imp,
    training_date  = Sys.time(),
    n_training     = length(train_idx),
    p_features     = p_final,
    p_original     = p_orig,
    calib_shift    = NULL,
    calib_scale    = NULL,
    split_strategy = split_obj$strategy,
    n_groups       = length(unique(grp)),
    app_domain     = .nts_rt_build_app_domain(fp_mat, train_idx, app_min_similarity),
    model_version  = .NTS_RT_MODEL_VERSION
  ), class = "nts_rt_model")

  if (verbose) message("Model built successfully.")
  model
}


#' @export
print.nts_rt_model <- function(x, ...) {
  cat("nts_rt_model (v", x$model_version, ")\n", sep = "")
  cat("  Method:    ", x$method, "\n")
  if (length(x$fp_types)) cat("  Fingerprints:", paste(x$fp_types, collapse = ", "), "\n")
  if (length(x$descriptor_cols)) cat("  Descriptors:", length(x$descriptor_cols), "columns\n")
  cat("  Features:  ", x$p_features, " (from ", x$p_original, " raw)\n", sep = "")
  cat("  Training:  ", x$n_training, " compounds\n", sep = "")
  if (!is.null(x$split_strategy)) {
    cat("  Split:     ", x$split_strategy,
        if (!is.null(x$n_groups)) paste0(" (", x$n_groups, " groups)") else "",
        "\n", sep = "")
  }
  cat("  Test R2:   ", round(x$test_metrics$R2, 3),
      "  RMSE: ", round(x$test_metrics$RMSE, 1), " s",
      "  MAE: ", round(x$test_metrics$MAE, 1), " s\n", sep = "")
  cat("  Trained:   ", format(x$training_date), "\n", sep = "")
  if (!is.null(x$calib_shift)) {
    cat("  Calibrated: shift=", round(x$calib_shift, 1),
        ", scale=", round(x$calib_scale, 3), "\n", sep = "")
  }
  # Applicability domain: named so the user can see *how* out-of-domain is judged
  # (legacy / externally trained models carry none).
  if (!is.null(x$app_domain) && !is.null(x$app_domain$sim_threshold)) {
    cat("  App. domain: Tanimoto >= ",
        round(x$app_domain$sim_threshold, 3),
        if (!is.null(x$app_domain$fp_bin)) paste0(" (", nrow(x$app_domain$fp_bin),
                                                  " training fingerprints)") else "",
        "\n", sep = "")
  } else {
    cat("  App. domain: none (predictions are extrapolation)\n")
  }
  invisible(x)
}


# Phase 3: predict_rt + evaluate_rt_model -------------------------------------

#' Predict retention times from SMILES
#'
#' @param model An `nts_rt_model` object.
#' @param smiles Character vector of SMILES to predict.
#' @param return_interval If TRUE, return prediction intervals (lower/upper).
#' @param BPPARAM A BiocParallelParam object.
#' @param fingerprint_matrix Optional precomputed fingerprint matrix aligned to
#'   `smiles`; skips on-the-fly fingerprinting when supplied.
#' @param descriptor_matrix Optional precomputed descriptor matrix aligned to
#'   `smiles`.
#'
#' @return A data.frame with column `predicted_rt`, plus `max_train_similarity`
#'   (max Tanimoto similarity to a training compound) and `applicability`
#'   (`"in_domain"` / `"out_of_domain"`), and `lower`/`upper` if
#'   `return_interval = TRUE`. The similarity columns are `NA` when the model
#'   carries no applicability domain (e.g. models built before 0.1.20 or with no
#'   fingerprints).
#' @export
predict_rt <- function(model,
                       smiles,
                       return_interval = FALSE,
                       BPPARAM = BiocParallel::SerialParam(),
                       fingerprint_matrix = NULL,
                       descriptor_matrix = NULL) {
  if (!inherits(model, "nts_rt_model")) {
    stop("model must be an nts_rt_model object.", call. = FALSE)
  }

  smi <- trimws(as.character(smiles))
  has_desc <- !is.null(model$descriptor_cols) && length(model$descriptor_cols) > 0L
  has_fp <- !is.null(model$fp_types) && length(model$fp_types) > 0L

  if (!is.null(fingerprint_matrix)) {
    if (!is.matrix(fingerprint_matrix)) stop("fingerprint_matrix must be a matrix.", call. = FALSE)
    fp_mat <- fingerprint_matrix
  } else if (has_fp) {
    if (!length(smi) || all(is.na(smi) | !nzchar(smi))) {
      stop("No valid SMILES strings provided.", call. = FALSE)
    }
    fp_mat <- compute_fingerprints(smi, types = model$fp_types, fp_mode = "matrix", BPPARAM = BPPARAM)
  } else {
    fp_mat <- NULL
  }

  if (!is.null(descriptor_matrix)) {
    if (!is.matrix(descriptor_matrix)) stop("descriptor_matrix must be a matrix.", call. = FALSE)
    desc_mat <- descriptor_matrix
  } else if (has_desc) {
    if (!length(smi) || all(is.na(smi) | !nzchar(smi))) {
      stop("No valid SMILES strings provided.", call. = FALSE)
    }
    desc_dt <- .nts_compute_descriptors(smi, BPPARAM = BPPARAM, prune = FALSE)
    avail_desc <- intersect(model$descriptor_cols, names(desc_dt))
    if (!length(avail_desc)) stop("None of the model descriptor columns found in computed descriptors.", call. = FALSE)
    desc_mat <- as.matrix(desc_dt[, avail_desc, with = FALSE])
  } else {
    desc_mat <- NULL
  }

  # Combine fingerprints + descriptors
  if (!is.null(fp_mat) && !is.null(desc_mat)) {
    stopifnot(nrow(fp_mat) == nrow(desc_mat))
    X <- cbind(fp_mat, desc_mat)
  } else if (!is.null(fp_mat)) {
    X <- fp_mat
  } else if (!is.null(desc_mat)) {
    X <- desc_mat
  } else {
    stop("Model requires neither fingerprints nor descriptors; cannot predict.", call. = FALSE)
  }

  # Align columns to training order
  if (!is.null(model$kept_colnames)) {
    missing_cols <- setdiff(model$kept_colnames, colnames(X))
    if (length(missing_cols)) {
      # Fill with the *training centre* rather than a raw 0. The columns are
      # standardised afterwards, so the centre maps to z = 0 (the neutral value)
      # whereas a raw 0 would map to a large, arbitrary z-score and bias the
      # prediction (finding P2-2).
      ctr <- model$preproc_center
      if (!is.null(ctr) && length(ctr) == length(model$kept_colnames)) {
        names(ctr) <- model$kept_colnames
        fill_val <- suppressWarnings(as.numeric(ctr[missing_cols]))
        fill_val[!is.finite(fill_val)] <- 0
      } else {
        fill_val <- rep(0, length(missing_cols))
      }
      for (k in seq_along(missing_cols)) {
        X <- cbind(X, matrix(fill_val[k], nrow = nrow(X), ncol = 1,
                             dimnames = list(NULL, missing_cols[k])))
      }
      warning(length(missing_cols), " model feature column(s) could not be ",
              "computed (", paste(utils::head(missing_cols, 8L), collapse = ", "),
              if (length(missing_cols) > 8L) ", ..." else "",
              ") and were imputed with the training mean; predictions for these ",
              "compounds are less reliable.", call. = FALSE)
    }
    X <- X[, model$kept_colnames, drop = FALSE]
  }

  p_orig <- ncol(X)
  if (p_orig != model$p_original && is.null(model$descriptor_cols)) {
    stop("Feature dimension mismatch: expected ", model$p_original,
         " columns, got ", p_orig, ".", call. = FALSE)
  }

  # NOTE: no further pruning here. `model$kept_colnames` already holds the
  # post-preprocessing columns (see build_rt_model), so X above is already in
  # training order and width. `model$nzv_idx` / `model$cor_idx` are indices into
  # the *original* un-pruned matrix and must never be applied to X.

  if (!is.null(model$preproc_center) && length(model$preproc_center) == ncol(X)) {
    X <- scale(X, center = model$preproc_center, scale = model$preproc_scale)
  }

  pred <- if (model$method == "rf") {
    as.numeric(predict(model$model, X)$predictions)
  } else {
    as.numeric(predict(model$model, xgboost::xgb.DMatrix(X)))
  }

  if (!is.null(model$calib_shift)) {
    pred <- model$calib_shift + model$calib_scale * pred
  }

  out <- data.frame(predicted_rt = pred, stringsAsFactors = FALSE)

  # Applicability domain: max Tanimoto similarity to the training compounds.
  # A prediction far from any training chemistry is flagged out-of-domain.
  ad <- model$app_domain
  sim <- NULL
  if (!is.null(ad) && !is.null(ad$fp_bin) && !is.null(fp_mat) &&
      nrow(fp_mat) == length(pred)) {
    sim <- tryCatch(.nts_rt_tanimoto(fp_mat, ad$fp_bin), error = function(e) NULL)
  }
  if (!is.null(sim) && nrow(sim) == length(pred)) {
    max_sim <- vapply(seq_len(nrow(sim)), function(i) {
      v <- max(sim[i, ], na.rm = TRUE)
      if (is.finite(v)) v else NA_real_
    }, numeric(1L))
    thr <- suppressWarnings(as.numeric(ad$sim_threshold)[1L])
    if (!is.finite(thr)) thr <- 0.6
    out$max_train_similarity <- max_sim
    out$applicability <- ifelse(
      is.na(max_sim), NA_character_,
      ifelse(max_sim >= thr, "in_domain", "out_of_domain")
    )
  } else {
    out$max_train_similarity <- NA_real_
    out$applicability <- NA_character_
  }

  if (isTRUE(return_interval)) {
    # The interval is the model's test RMSE (±2 s), which is only calibrated for
    # compounds inside the applicability domain. Reporting it for a compound far
    # from any training chemistry would be spurious precision, so those rows get
    # NA (finding P2-5). Models without an applicability domain keep the
    # interval for every row (previous behaviour).
    se <- suppressWarnings(as.numeric(model$test_metrics$RMSE)[1L])
    if (!is.finite(se) || se < 0) se <- NA_real_
    out$lower <- if (is.na(se)) NA_real_ else pmax(0, pred - 2 * se)
    out$upper <- if (is.na(se)) NA_real_ else pred + 2 * se
    if ("applicability" %in% names(out)) {
      ood <- !is.na(out$applicability) & out$applicability == "out_of_domain"
      if (any(ood)) {
        out$lower[ood] <- NA_real_
        out$upper[ood] <- NA_real_
      }
    }
  }
  out
}


#' Evaluate an RT model on external test data
#'
#' @param model An `nts_rt_model`.
#' @param test_smiles Character vector of test SMILES.
#' @param test_rt Numeric vector of observed RTs.
#'
#' @return A list with `metrics` (named numeric), `predicted`, and `observed`.
#' @export
evaluate_rt_model <- function(model, test_smiles, test_rt) {
  if (!inherits(model, "nts_rt_model")) {
    stop("model must be an nts_rt_model object.", call. = FALSE)
  }
  obs <- suppressWarnings(as.numeric(test_rt))
  # predict_rt() silently drops rows whose SMILES is NA / empty (via
  # compute_fingerprints()'s valid_mask), so its output is shorter than the
  # input. Re-expand it against the same validity mask before comparing,
  # otherwise observed and predicted are silently misaligned.
  smi_ok <- !is.na(test_smiles) & nzchar(trimws(as.character(test_smiles)))
  pred <- predict_rt(model, test_smiles)
  pr <- pred$predicted_rt
  predicted <- rep(NA_real_, length(smi_ok))
  if (length(pr) == sum(smi_ok)) {
    predicted[smi_ok] <- pr
  } else {
    # Copying positionally here would pair each prediction with the wrong
    # compound as soon as one SMILES is invalid, silently corrupting the
    # metrics. Leave the unmatched rows NA instead and say so (finding P2-3).
    warning("evaluate_rt_model: predict_rt() returned ", length(pr),
            " predictions for ", sum(smi_ok), " valid SMILES, so they cannot be ",
            "matched to the observations; affected rows are reported as NA and ",
            "excluded from the metrics.", call. = FALSE)
  }
  ok <- !is.na(obs) & is.finite(obs) & !is.na(predicted)
  m <- .nts_rt_metrics(obs[ok], predicted[ok])
  message(sprintf("Evaluation: RMSE=%.1f, R2=%.3f, MAE=%.1f, median_err=%.1f",
                  m$RMSE, m$R2, m$MAE, m$median_error))
  list(metrics = c(RMSE = m$RMSE, R2 = m$R2, MAE = m$MAE,
                   median_error = m$median_error, mean_error = m$mean_error),
       predicted = predicted, observed = obs)
}


# Phase 4: save/load/transfer --------------------------------------------------

#' Save an RT model to disk
#' @param model An `nts_rt_model` object.
#' @param file Path to output .rds file.
#' @return Invisibly, the file path.
#' @export
save_rt_model <- function(model, file) {
  if (!inherits(model, "nts_rt_model")) stop("model must be an nts_rt_model object.", call. = FALSE)
  saveRDS(model, file)
  invisible(file)
}


#' Load an RT model from disk
#' @param file Path to .rds file.
#' @return An `nts_rt_model` object.
#' @export
load_rt_model <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  model <- readRDS(file)
  if (!inherits(model, "nts_rt_model")) stop("File does not contain an nts_rt_model object.", call. = FALSE)
  mv <- model$model_version
  if (is.null(mv) || !mv %in% c(.NTS_RT_MODEL_LEGACY_VERSIONS, .NTS_RT_MODEL_VERSION)) {
    warning("Model version ", if (is.null(mv)) "(none)" else mv,
            " may not be compatible.", call. = FALSE)
  }
  model
}


#' Transfer-learn an RT model to new data
#'
#' @param base_model An `nts_rt_model` trained on source data.
#' @param smiles Character vector of target-domain SMILES.
#' @param rt Numeric vector of target-domain RTs.
#' @param method `"shift_correct"` (default): predict on target data and fit
#'   a linear correction. `"retrain"`: train a new model from scratch.
#' @param ... Passed to `build_rt_model` when `method = "retrain"`.
#'
#' @return An updated `nts_rt_model`.
#' @export
transfer_rt_model <- function(base_model, smiles, rt,
                              method = c("shift_correct", "retrain"), ...) {
  if (!inherits(base_model, "nts_rt_model")) {
    stop("base_model must be an nts_rt_model object.", call. = FALSE)
  }
  method <- match.arg(method)
  if (method == "retrain") {
    args <- list(smiles = smiles, rt = rt,
                 fingerprint_types = base_model$fp_types,
                 method = base_model$method)
    dots <- list(...)
    for (nm in names(dots)) args[[nm]] <- dots[[nm]]
    return(do.call(build_rt_model, args))
  }
  pred <- predict_rt(base_model, smiles)
  obs <- suppressWarnings(as.numeric(rt))
  ok <- !is.na(obs) & is.finite(obs) & !is.na(pred$predicted_rt)
  if (sum(ok) < 5L) stop("Need at least 5 valid pairs for calibration; got ", sum(ok), ".", call. = FALSE)
  fit <- stats::lm(obs[ok] ~ pred$predicted_rt[ok])
  base_model$calib_shift <- as.numeric(stats::coef(fit)[1L])
  base_model$calib_scale <- as.numeric(stats::coef(fit)[2L])
  message(sprintf("Transfer calibration: obs = %.2f + %.3f * pred (R2=%.3f)",
                  base_model$calib_shift, base_model$calib_scale,
                  summary(fit)$r.squared))
  base_model
}


# ============================================================================
# Internal helpers
# ============================================================================

# --- Group-aware splitting + applicability domain (finding B5) ---------------
#
# Two leakage/misuse risks in RT modelling are addressed here:
#   * A random train/test split scatters chemical analogues (homologues,
#     isomers) across the split, so the reported test error is optimistic.
#     Groups (user-supplied or fingerprint clusters) are held out together.
#   * A model applied far outside its training chemistry gives confident but
#     meaningless RT predictions. Each prediction is scored by its Tanimoto
#     similarity to the nearest training compound and flagged in/out of domain.

# Row-aligned Tanimoto similarity between two binary fingerprint matrices
# (jaccard on the shared column set). Returns an nrow(A) x nrow(B) matrix, or
# NULL when the matrices cannot be aligned.
.nts_rt_tanimoto <- function(A, B = NULL) {
  if (is.null(B)) B <- A
  if (is.null(dim(A)) || is.null(dim(B))) return(NULL)
  if (is.null(colnames(A)) || is.null(colnames(B))) {
    if (ncol(A) != ncol(B)) return(NULL)
  } else {
    cols <- intersect(colnames(A), colnames(B))
    if (!length(cols)) return(NULL)
    A <- A[, cols, drop = FALSE]
    B <- B[, cols, drop = FALSE]
  }
  Ab <- (A > 0) * 1
  Bb <- (B > 0) * 1
  inter <- Ab %*% t(Bb)
  denom <- outer(rowSums(Ab), rowSums(Bb), "+") - inter
  sim <- inter / denom
  sim[!is.finite(sim)] <- 0
  sim
}

# Validate a user-supplied `group` vector against the input rows.
.nts_rt_check_group <- function(group, n, rt = NULL) {
  if (is.null(group)) return(invisible(NULL))
  if (length(group) != n) {
    stop("`group` must have one value per input compound (", n, " rows); got ",
         length(group), ".", call. = FALSE)
  }
  invisible(NULL)
}

# Resolve compound -> group labels: user-supplied, else average-linkage
# clustering of the fingerprints cut at `cluster_similarity`. Falls back to
# singleton groups (i.e. a plain random split) when clustering is not viable.
.nts_rt_resolve_groups <- function(group = NULL, n, fp = NULL,
                                   cluster_similarity = 0.8, seed = 42) {
  if (!is.null(group)) {
    g <- as.character(group)
    if (length(g) != n) {
      stop("`group` must have the same length as the training data (", n, ").",
           call. = FALSE)
    }
    bad <- is.na(g) | !nzchar(g)
    if (any(bad)) g[bad] <- paste0("..ungrouped..", which(bad))
    return(g)
  }
  if (is.null(fp) || nrow(fp) < 4L || ncol(fp) < 1L) {
    return(as.character(seq_len(n)))
  }
  sim_thr <- suppressWarnings(as.numeric(cluster_similarity)[1L])
  if (!is.finite(sim_thr) || sim_thr <= 0 || sim_thr >= 1) {
    return(as.character(seq_len(n)))
  }
  # Guard the O(n^2) similarity matrix for very large libraries.
  if (n > 5000L) {
    return(as.character(seq_len(n)))
  }
  sim <- .nts_rt_tanimoto(fp)
  if (is.null(sim)) return(as.character(seq_len(n)))
  sim[is.na(sim)] <- 0
  diag(sim) <- 1
  hc <- stats::hclust(stats::as.dist(1 - sim), method = "average")
  cl <- stats::cutree(hc, h = max(0, 1 - sim_thr))
  as.character(cl)
}

# Assign whole groups to train/test. Groups are added in shuffled order to the
# training side until the target size is (approximately) reached; at least one
# group is always held out for testing. When there is no real grouping this
# degrades to a seeded random split of the requested size.
.nts_rt_group_split <- function(groups, n, train_split, seed = 42) {
  p <- suppressWarnings(as.numeric(train_split)[1L])
  if (!is.finite(p) || p <= 0 || p >= 1) p <- 0.8
  target <- max(1L, min(n - 1L, round(n * p)))
  ug <- unique(groups)

  if (length(ug) <= 1L || length(ug) >= n) {
    set.seed(seed)
    train_idx <- sort(sample.int(n, target))
    return(list(
      train_idx = train_idx,
      test_idx = setdiff(seq_len(n), train_idx),
      strategy = if (length(ug) >= n) "random" else "random_single_group"
    ))
  }

  set.seed(seed)
  ord <- sample(ug)
  sizes <- vapply(ord, function(g) sum(groups == g), integer(1L))
  cum <- cumsum(sizes)
  k <- which(cum >= target)
  k <- if (length(k)) k[1L] else length(ord)
  # Never place every group on the training side.
  if (k >= length(ord)) k <- max(1L, length(ord) - 1L)
  sel <- ord[seq_len(k)]
  train_idx <- which(groups %in% sel)
  if (!length(train_idx) || length(train_idx) >= n) {
    set.seed(seed)
    train_idx <- sort(sample.int(n, target))
    return(list(train_idx = train_idx,
                test_idx = setdiff(seq_len(n), train_idx),
                strategy = "random_fallback"))
  }
  list(
    train_idx = sort(train_idx),
    test_idx = setdiff(seq_len(n), train_idx),
    strategy = "grouped"
  )
}

# Capture the training fingerprints (binary) needed to score the applicability
# of future predictions. Returns NULL when no fingerprints are available.
.nts_rt_build_app_domain <- function(fp_mat, train_idx, app_min_similarity = 0.6) {
  if (is.null(fp_mat) || !ncol(fp_mat) || !length(train_idx)) return(NULL)
  thr <- suppressWarnings(as.numeric(app_min_similarity)[1L])
  if (!is.finite(thr) || thr < 0 || thr > 1) thr <- 0.6
  list(
    fp_bin = tryCatch((fp_mat[train_idx, , drop = FALSE] > 0) * 1,
                      error = function(e) NULL),
    sim_threshold = thr
  )
}

# Screen features (near-zero variance + highly correlated) and standardise.
#
# `fit_idx` selects the rows used to *learn* the screening decisions and the
# centring/scaling constants. It must be the training rows only: fitting the
# transform on the full matrix (the historical behaviour) leaks test-set
# information into training via the column means/sds and the correlation
# structure, which makes the reported test error optimistic (finding P1-1).
# `fit_idx = NULL` keeps the old fit-on-everything behaviour.
.nts_rt_preprocess <- function(fp_mat, fit_idx = NULL) {
  p <- ncol(fp_mat)
  if (is.null(fit_idx)) fit_idx <- seq_len(nrow(fp_mat))
  fit_idx <- as.integer(fit_idx)
  fit_idx <- sort(unique(fit_idx[!is.na(fit_idx) &
                                  fit_idx >= 1L & fit_idx <= nrow(fp_mat)]))
  if (!length(fit_idx)) fit_idx <- seq_len(nrow(fp_mat))
  fit_mat <- fp_mat[fit_idx, , drop = FALSE]

  if (p < 2L) return(list(
    matrix = fp_mat, n_removed = 0L, nzv_removed = 0L, cor_removed = 0L,
    nzv_idx = integer(0), cor_idx = integer(0), center = NULL, scale = NULL,
    removed_idx = integer(0)
  ))
  col_vars <- apply(fit_mat, 2L, stats::var)
  uniq_ratio <- apply(fit_mat, 2L, function(x) max(table(x)) / length(x))
  nzv_idx <- which(col_vars < 1e-6 | uniq_ratio > 0.95)
  if (length(nzv_idx)) { mat_sub <- fit_mat[, -nzv_idx, drop = FALSE] }
  else { mat_sub <- fit_mat }

  cor_idx <- integer(0)
  if (ncol(mat_sub) > 1L) {
    cor_mat <- suppressWarnings(stats::cor(mat_sub))
    if (!all(is.na(cor_mat))) {
      cor_mat[is.na(cor_mat)] <- 0
      cor_mat[lower.tri(cor_mat, diag = TRUE)] <- 0
      hi_cor <- which(abs(cor_mat) > 0.95, arr.ind = TRUE)
      if (nrow(hi_cor)) {
        kept_cols <- setdiff(seq_len(p), nzv_idx)
        cor_idx <- kept_cols[unique(hi_cor[, 2L])]
      }
    }
  }

  all_remove <- unique(c(nzv_idx, cor_idx))
  min_keep <- 10L
  if (length(all_remove) && (p - length(all_remove)) >= min_keep) {
    keep <- setdiff(seq_len(p), all_remove)
  } else if (length(all_remove)) {
    if (length(nzv_idx) && (p - length(nzv_idx)) >= min_keep) {
      keep <- setdiff(seq_len(p), nzv_idx)
      all_remove <- nzv_idx
    } else { keep <- seq_len(p); all_remove <- integer(0) }
  } else { keep <- seq_len(p) }

  fit_clean <- fit_mat[, keep, drop = FALSE]
  cs <- scale(fit_clean, center = TRUE, scale = TRUE)
  center <- attr(cs, "scaled:center")
  scale  <- attr(cs, "scaled:scale")
  if (is.null(center)) center <- rep(0, ncol(fit_clean))
  if (is.null(scale))  scale  <- rep(1, ncol(fit_clean))
  scale[scale == 0] <- 1
  center[!is.finite(center)] <- 0
  scale[!is.finite(scale)] <- 1

  prep <- list(n_removed = length(all_remove),
               nzv_removed = length(nzv_idx), cor_removed = length(cor_idx),
               nzv_idx = as.integer(nzv_idx), cor_idx = as.integer(cor_idx),
               center = as.numeric(center), scale = as.numeric(scale),
               removed_idx = as.integer(all_remove))
  # Apply the (train-fitted) transform to every row, so callers get a matrix
  # that is already in the training column space.
  prep$matrix <- .nts_rt_apply_preproc(fp_mat, prep)
  prep
}

# Apply a transform learned by .nts_rt_preprocess() to a (new) feature matrix.
# `mat` must have the same columns as the matrix the transform was fitted on.
.nts_rt_apply_preproc <- function(mat, prep) {
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (!is.null(prep$removed_idx) && length(prep$removed_idx)) {
    rm_idx <- prep$removed_idx[prep$removed_idx >= 1L & prep$removed_idx <= ncol(mat)]
    if (length(rm_idx)) mat <- mat[, -rm_idx, drop = FALSE]
  }
  if (!is.null(prep$center) && length(prep$center) == ncol(mat)) {
    mat <- scale(mat, center = prep$center, scale = prep$scale)
  }
  mat
}


.nts_rt_train_rf <- function(X, y, tune_params, cv_folds, seed, verbose) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("ranger is required. Install: install.packages('ranger')", call. = FALSE)
  }
  p <- ncol(X)
  n <- nrow(X)
  best_mtry <- max(1L, round(sqrt(p) / 2))
  best_node <- 5L

  if (length(tune_params$mtry)) best_mtry <- as.integer(tune_params$mtry)[1L]
  if (length(tune_params$min.node.size)) best_node <- as.integer(tune_params$min.node.size)[1L]

  # For small datasets (<20K), do lightweight CV tuning
  if (n < 20000L && !length(tune_params$mtry) && !length(tune_params$min.node.size)) {
    mtry_vals <- unique(pmax(1L, round(c(sqrt(p)/2, sqrt(p)))))
    node_vals <- if (n < 100L) c(3L, 5L) else c(5L, 10L)
    fold_ids <- sample(rep(seq_len(cv_folds), length.out = n))
    best_rmse <- Inf
    for (mt in mtry_vals) {
      for (nd in node_vals) {
        cv_pred <- numeric(n)
        for (f in seq_len(cv_folds)) {
          tr <- fold_ids != f; te <- fold_ids == f
          rf <- ranger::ranger(x = X[tr, , drop = FALSE], y = y[tr],
                               mtry = mt, min.node.size = nd, num.trees = 200L,
                               seed = seed, verbose = FALSE)
          cv_pred[te] <- as.numeric(predict(rf, X[te, , drop = FALSE])$predictions)
        }
        rmse <- sqrt(mean((y - cv_pred)^2))
        if (rmse < best_rmse) { best_rmse <- rmse; best_mtry <- mt; best_node <- nd }
        if (verbose) message(sprintf("  CV mtry=%d node=%d -> RMSE=%.1f", mt, nd, rmse))
      }
    }
    if (verbose) message(sprintf("  Best: mtry=%d, min.node.size=%d, CV RMSE=%.1f",
                                  best_mtry, best_node, best_rmse))
  }

  data.frame(fold = 1L:cv_folds, RMSE = NA_real_, R2 = NA_real_, MAE = NA_real_,
             stringsAsFactors = FALSE)

  if (verbose) message(sprintf("  Training final RF: mtry=%d, node=%d, trees=500", best_mtry, best_node))
  ranger::ranger(x = X, y = y, mtry = best_mtry, min.node.size = best_node,
                 num.trees = 500L, importance = "impurity", seed = seed, verbose = verbose)
}


.nts_rt_train_xgb <- function(X, y, tune_params, cv_folds, seed, verbose) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("xgboost is required. Install: install.packages('xgboost')", call. = FALSE)
  }
  n <- nrow(X)
  dtrain <- xgboost::xgb.DMatrix(X, label = y)

  best_params <- list(nrounds = 500L, max_depth = 6L, eta = 0.05)
  if (length(tune_params$nrounds)) best_params$nrounds <- as.integer(tune_params$nrounds)
  if (length(tune_params$max_depth)) best_params$max_depth <- as.integer(tune_params$max_depth)
  if (length(tune_params$eta)) best_params$eta <- as.numeric(tune_params$eta)

  # For small datasets, do lightweight CV grid tuning
  if (n < 20000L && !length(tune_params$nrounds) && !length(tune_params$max_depth) && !length(tune_params$eta)) {
    eta_vals <- c(0.05, 0.1)
    depth_vals <- c(4L, 6L)
    fold_ids <- sample(rep(seq_len(cv_folds), length.out = n))
    best_rmse <- Inf
    for (et in eta_vals) {
      for (dp in depth_vals) {
        cv_pred <- numeric(n)
        for (f in seq_len(cv_folds)) {
          tr <- fold_ids != f; te <- fold_ids == f
          dtr <- xgboost::xgb.DMatrix(X[tr, , drop = FALSE], label = y[tr])
          dte <- xgboost::xgb.DMatrix(X[te, , drop = FALSE], label = y[te])
          bst <- xgboost::xgb.train(
            params = list(max_depth = dp, eta = et, objective = "reg:squarederror",
                          subsample = 0.8, colsample_bytree = 0.8),
            data = dtr, nrounds = 300L, verbose = 0,
            watchlist = list(test = dte), early_stopping_rounds = 20
          )
          cv_pred[te] <- as.numeric(predict(bst, dte))
        }
        rmse <- sqrt(mean((y - cv_pred)^2))
        if (rmse < best_rmse) {
          best_rmse <- rmse; best_params <- list(nrounds = 500L, max_depth = dp, eta = et)
        }
        if (verbose) message(sprintf("  CV depth=%d eta=%.2f -> RMSE=%.1f", dp, et, rmse))
      }
    }
    if (verbose) message(sprintf("  Best: depth=%d, eta=%.2f, CV RMSE=%.1f",
                                  best_params$max_depth, best_params$eta, best_rmse))
  }

  if (verbose) message(sprintf("  Training final XGB: depth=%d, eta=%.3f, nrounds=%d",
                                best_params$max_depth, best_params$eta, best_params$nrounds))
  xgboost::xgb.train(
    params = list(max_depth = best_params$max_depth, eta = best_params$eta,
                  objective = "reg:squarederror", subsample = 0.8, colsample_bytree = 0.8),
    data = dtrain, nrounds = best_params$nrounds, verbose = 0
  )
}


.nts_rt_metrics <- function(obs, pred) {
  ok <- !is.na(obs) & !is.na(pred)
  obs <- obs[ok]; pred <- pred[ok]
  if (!length(obs)) return(list(RMSE = NA_real_, R2 = NA_real_, MAE = NA_real_,
                                 median_error = NA_real_, mean_error = NA_real_))
  err <- obs - pred
  ss_res <- sum(err^2)
  ss_tot <- sum((obs - mean(obs))^2)
  list(RMSE = sqrt(mean(err^2)),
       R2 = if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_,
       MAE = mean(abs(err)),
       median_error = stats::median(abs(err)),
       mean_error = mean(err))
}


.nts_rt_var_importance <- function(fit, method, fp_types, fp_bounds, colnames_kept) {
  if (method == "rf") {
    imp <- as.numeric(ranger::importance(fit))
  } else {
    imp_mat <- xgboost::xgb.importance(model = fit)
    imp_named <- setNames(rep(0, length(colnames_kept)), colnames_kept)
    if (nrow(imp_mat)) {
      for (i in seq_len(nrow(imp_mat))) {
        feat <- imp_mat$Feature[i]
        idx <- suppressWarnings(as.integer(feat))
        if (is.na(idx)) {
          m <- match(feat, colnames_kept)
          if (!is.na(m)) idx <- m
        }
        if (!is.na(idx) && idx > 0 && idx <= length(imp_named)) {
          imp_named[idx] <- imp_mat$Gain[i]
        }
      }
    }
    imp <- as.numeric(imp_named)
  }
  n <- length(imp)
  if (!n) return(data.frame(bit_index = integer(0), type = character(0),
                             importance = numeric(0), stringsAsFactors = FALSE))
  types <- character(n)
  for (i in seq_len(n)) {
    cn <- colnames_kept[i]
    for (tp in names(fp_bounds)) {
      pat <- paste0("^", tp, "_")
      if (grepl(pat, cn)) { types[i] <- tp; break }
    }
    if (!nzchar(types[i])) types[i] <- "unknown"
  }
  data.frame(bit_index = seq_len(n), type = types, importance = imp, stringsAsFactors = FALSE)
}


# Phase 5: DB builder integration -------------------------------------------------

.nts_add_predicted_rt <- function(metadata, rt_model, smiles_col = "smiles") {
  if (!inherits(rt_model, "nts_rt_model")) {
    stop("rt_model must be an nts_rt_model object.", call. = FALSE)
  }
  if (!is.data.frame(metadata) || !nrow(metadata)) return(metadata)
  smiles_cols <- c(smiles_col, "SMILES", "smiles", "Smiles")
  found <- intersect(smiles_cols, names(metadata))
  if (!length(found)) {
    warning("No SMILES column found in metadata; cannot add predicted RT.", call. = FALSE)
    return(metadata)
  }
  sc <- found[1L]
  smi <- as.character(metadata[[sc]])
  ok <- !is.na(smi) & nzchar(trimws(smi))
  if (!sum(ok)) {
    warning("No valid SMILES in metadata; cannot add predicted RT.", call. = FALSE)
    return(metadata)
  }
  pred <- predict_rt(rt_model, smi[ok])
  metadata$predicted_rt <- NA_real_
  metadata$predicted_rt[ok] <- pred$predicted_rt
  metadata
}


# RT annotation integration for annotate_by_spectral_match ---------------------

#' Annotate RT prediction results onto db_annotations
#'
#' Uses a trained `nts_rt_model` to predict RT for each annotated compound
#' in `@db_annotations` and writes `predicted_rt` column.
#'
#' @param project An [NTSProject].
#' @param rt_model An `nts_rt_model` object.
#' @param BPPARAM A BiocParallelParam for fingerprint computation.
#'
#' @return The updated [NTSProject] invisibly.
#' @export
annotate_rt_prediction <- function(project, rt_model, BPPARAM = BiocParallel::SerialParam()) {
  if (!methods::is(project, "NTSProject")) stop("project must be an NTSProject.", call. = FALSE)
  if (!inherits(rt_model, "nts_rt_model")) stop("rt_model must be an nts_rt_model object.", call. = FALSE)

  db <- methods::slot(project, "db_annotations")
  if (!is.data.frame(db) || !nrow(db)) {
    warning("Empty db_annotations; nothing to predict.", call. = FALSE)
    return(invisible(project))
  }

  smi_col <- if ("SMILES" %in% names(db)) "SMILES"
             else if ("smiles" %in% names(db)) "smiles"
             else NULL
  if (is.null(smi_col)) {
    warning("No SMILES column in db_annotations; cannot predict RT.", call. = FALSE)
    return(invisible(project))
  }

  smi <- as.character(db[[smi_col]])
  ok <- !is.na(smi) & nzchar(trimws(smi))
  if (!any(ok)) {
    warning("No valid SMILES in db_annotations.", call. = FALSE)
    return(invisible(project))
  }

  message("Predicting RT for ", sum(ok), " / ", nrow(db), " annotations...")
  pred <- predict_rt(rt_model, smi[ok], BPPARAM = BPPARAM)
  db$predicted_rt <- NA_real_
  db$predicted_rt[ok] <- pred$predicted_rt

  # Compute rt_delta = observed_rt - predicted_rt
  feat_names <- if ("name" %in% names(db)) trimws(as.character(db$name))
                else if ("feature_id" %in% names(db)) trimws(as.character(db$feature_id))
                else rep(NA_character_, nrow(db))
  ft <- get_feature_table(project)
  if (!is.null(ft) && nrow(ft) && "rt" %in% names(ft) && "name" %in% names(ft)) {
    rt_map <- setNames(as.numeric(ft$rt), trimws(as.character(ft$name)))
    obs_rt <- rt_map[feat_names]
    db$rt_delta <- obs_rt - db$predicted_rt
  } else {
    db$rt_delta <- NA_real_
  }

  # Rank candidate annotations by |rt_delta| *within each feature*. A global
  # sort would interleave the candidates of different features, breaking the
  # feature-contiguous layout the rest of the pipeline relies on (finding P2-4).
  if (any(is.finite(db$rt_delta))) {
    db <- db[order(feat_names, abs(db$rt_delta), na.last = TRUE), , drop = FALSE]
  }

  methods::slot(project, "db_annotations") <- db
  invisible(project)
}


#' Annotate RT predictions onto @ms1_annotations
#'
#' Reads SMILES from `@ms1_annotations`, gets observed RT from `@feature_table`,
#' predicts RT via `rt_model`, computes `rt_delta = observed - predicted`,
#' and writes `predicted_rt` and `rt_delta` columns.
#'
#' Call after [merge_annotation_layers] / [apply_annotation_master_to_ms1],
#' once `@ms1_annotations` contains SMILES from spectral and network annotations.
#'
#' @param project An [NTSProject].
#' @param rt_model An `nts_rt_model` object (e.g. from [load_rt_model]).
#' @param BPPARAM A BiocParallelParam for fingerprint computation.
#'
#' @return The updated [NTSProject] invisibly.
#' @export
annotate_rt_to_ms1 <- function(project, rt_model, BPPARAM = BiocParallel::SerialParam()) {
  if (!methods::is(project, "NTSProject")) stop("project must be an NTSProject.", call. = FALSE)
  if (!inherits(rt_model, "nts_rt_model")) stop("rt_model must be an nts_rt_model object.", call. = FALSE)

  ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
  if (!nrow(ms1)) {
    warning("Empty ms1_annotations; nothing to predict.", call. = FALSE)
    return(invisible(project))
  }

  # Find SMILES column
  smi_col <- if ("master_smiles" %in% names(ms1)) "master_smiles"
             else if ("SMILES" %in% names(ms1)) "SMILES"
             else if ("smiles" %in% names(ms1)) "smiles"
             else NULL
  if (is.null(smi_col)) {
    warning("No SMILES column in ms1_annotations; run merge_annotation_layers first.", call. = FALSE)
    return(invisible(project))
  }

  smi <- as.character(ms1[[smi_col]])
  ok <- !is.na(smi) & nzchar(trimws(smi))
  if (!any(ok)) {
    warning("No valid SMILES in ms1_annotations.", call. = FALSE)
    return(invisible(project))
  }

  message("Predicting RT for ", sum(ok), " / ", nrow(ms1), " ms1 features...")
  pred <- predict_rt(rt_model, smi[ok], BPPARAM = BPPARAM)
  ms1$predicted_rt <- NA_real_
  ms1$predicted_rt[ok] <- pred$predicted_rt

  # rt_delta = observed_rt - predicted_rt
  ft <- get_feature_table(project)
  if (!is.null(ft) && nrow(ft) && "rt" %in% names(ft) && "name" %in% names(ft)) {
    rt_map <- setNames(as.numeric(ft$rt), trimws(as.character(ft$name)))
    obs_rt <- rt_map[trimws(as.character(ms1$name))]
    ms1$rt_delta <- obs_rt - ms1$predicted_rt
  } else {
    ms1$rt_delta <- NA_real_
  }

  # Sort by |rt_delta|
  if (any(is.finite(ms1$rt_delta))) {
    ms1 <- ms1[order(abs(ms1$rt_delta), na.last = TRUE), , drop = FALSE]
  }

  project <- set_ms1_annotations(project, .nts_ms1_sanitize_and_order(ms1))
  invisible(project)
}


# Descriptor computation -------------------------------------------------------

# rcdk 3.8.2 / CDK descriptors that recurse infinitely in native code and
# therefore CRASH the R process (uncatchable "segfault from C stack overflow").
# `LongestAliphaticChainDescriptor` triggers on any molecule that has an alkyl
# chain attached to an aromatic ring (e.g. toluene, benzoic acid, aspirin); a
# pure aliphatic chain or a bare aromatic ring is fine. Because the fault is a
# native segmentation fault, `tryCatch()` cannot intercept it, so the only safe
# mitigation is to drop these descriptors from the default set. Excluding one
# of ~53 descriptors has a negligible effect on the RT model.
.NTS_UNSAFE_DESCRIPTORS <- c(
  "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor"
)

#' Compute molecular descriptors using rcdk
#' @param smiles Character vector of SMILES
#' @param desc_names Character vector of descriptor names (NULL = all available)
#' @param BPPARAM BiocParallelParam for parallelisation
#' @param prune Drop all-NA / zero-variance / >50%-NA columns. Must stay `TRUE`
#'   for exploratory use, but **prediction** must pass `FALSE`: those filters are
#'   evaluated per batch, so a column that happens to be constant in one batch
#'   would be dropped and then silently imputed (with 0 when the model stores no
#'   training means) even though the model was trained on it.
#' @return data.table with SMILES + descriptor columns, NAs dropped
#' @keywords internal
#' @noRd
.nts_compute_descriptors <- function(smiles, desc_names = NULL,
                                      BPPARAM = BiocParallel::SerialParam(),
                                      prune = TRUE) {
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    stop("rcdk is required. Install: install.packages('rcdk')", call. = FALSE)
  }

  smi <- trimws(as.character(smiles))
  valid_mask <- !is.na(smi) & nzchar(smi)
  if (!any(valid_mask)) stop("No valid SMILES strings.", call. = FALSE)
  smi_clean <- smi[valid_mask]
  n <- length(smi_clean)

  if (is.null(desc_names)) desc_names <- rcdk::get.desc.names()
  # Drop descriptors known to segfault the process in native code (see above).
  unsafe <- intersect(desc_names, .NTS_UNSAFE_DESCRIPTORS)
  if (length(unsafe)) desc_names <- setdiff(desc_names, unsafe)

  message(sprintf("Computing %d descriptors for %d compounds...",
                  length(desc_names), n))

  desc_list <- BiocParallel::bplapply(seq_len(n), function(i) {
    s <- smi_clean[i]
    mol <- tryCatch(rcdk::parse.smiles(s)[[1]], error = function(e) NULL)
    if (is.null(mol)) {
      return(structure(list(smiles = s, .failed = TRUE), class = "data.frame"))
    }
    res <- tryCatch(
      suppressWarnings(rcdk::eval.desc(mol, desc_names)),
      error = function(e) NULL
    )
    if (is.null(res)) {
      return(structure(list(smiles = s, .failed = TRUE), class = "data.frame"))
    }
    res$smiles <- s
    res$.failed <- FALSE
    res
  }, BPPARAM = BPPARAM)

  all_rows <- vector("list", n)
  failed_count <- 0L
  for (i in seq_len(n)) {
    if (isTRUE(desc_list[[i]]$.failed)) {
      all_rows[[i]] <- NULL
      failed_count <- failed_count + 1L
    } else {
      all_rows[[i]] <- desc_list[[i]]
    }
  }
  if (failed_count > 0) warning(failed_count, " SMILES failed descriptor computation.", call. = FALSE)

  dt <- data.table::rbindlist(all_rows, fill = TRUE, use.names = TRUE)
  if (!nrow(dt)) stop("All SMILES failed descriptor computation.", call. = FALSE)

  dt$.failed <- NULL

  na_cols <- character(0); zv_cols <- character(0); high_na <- character(0)

  if (prune) {
    # Remove all-NA columns
    na_cols <- setdiff(names(dt), "smiles")
    na_cols <- na_cols[vapply(dt[, na_cols, with = FALSE], function(x) all(is.na(x)), logical(1L))]
    if (length(na_cols)) dt <- dt[, setdiff(names(dt), na_cols), with = FALSE]

    # Remove zero-variance columns. A column with fewer than two non-NA values has
    # no variance estimate, so keep it: a single-row batch would otherwise lose
    # every descriptor column here. Genuinely uninformative columns are caught by
    # the all-NA rule above and the >50%-NA rule below.
    desc_cols <- setdiff(names(dt), "smiles")
    zv_cols <- desc_cols[vapply(dt[, desc_cols, with = FALSE], function(x) {
      xc <- x[!is.na(x)]
      length(xc) >= 2L && stats::var(xc) < 1e-12
    }, logical(1L))]
    if (length(zv_cols)) dt <- dt[, setdiff(names(dt), zv_cols), with = FALSE]

    # Remove columns where >50% are NA
    desc_cols <- setdiff(names(dt), "smiles")
    na_thresh <- max(1L, floor(0.5 * nrow(dt)))
    high_na <- desc_cols[vapply(dt[, desc_cols, with = FALSE], function(x) {
      sum(is.na(x)) > na_thresh
    }, logical(1L))]
    if (length(high_na)) dt <- dt[, setdiff(names(dt), high_na), with = FALSE]
  }

  # Move smiles to first column
  sc <- setdiff(names(dt), "smiles")
  data.table::setcolorder(dt, c("smiles", sc))

  message(sprintf("Descriptors: %d columns for %d compounds (removed %d all-NA, %d zero-var, %d high-NA)",
                  length(sc), nrow(dt), length(na_cols), length(zv_cols), length(high_na)))
  dt[]
}
