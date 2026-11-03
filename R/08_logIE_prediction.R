# logIE prediction module =======================================================
#
# S3 class nts_logIE_model, predict_logIE, build/save/load/evaluate.
# Follows the API pattern of 08_rt_prediction.R.

# ---- S3 class -----------------------------------------------------------------

#' Create an nts_logIE_model object
#' @param model A trained xgboost or ranger model object.
#' @param method `"xgb"` or `"rf"`.
#' @param feature_set Character label for the feature combination used.
#' @param feature_cols Character vector of column names used during training.
#' @param n_features Integer number of features.
#' @param hp Character string describing hyperparameters.
#' @param train_metrics Named list with R2, RMSE, MAE.
#' @param test_metrics Named list with R2, RMSE, MAE.
#' @param cv_metrics Named list with R2, RMSE.
#' @param polarity `"positive"` or `"negative"`.
#' @param n_training,n_test Integer counts.
#' @param .best_nrounds Integer early-stopping round (XGBoost only).
#' @param .feature_centers Named numeric vector of training-set column means,
#'   used to impute required columns that a later `predict_logIE` call cannot
#'   compute. `NULL` for models built before the value was recorded.
#' @param .preproc_center,.preproc_scale Named numeric vectors of the
#'   standardisation constants learned on the *training rows* by
#'   [build_logIE_model] when `feature_screen = TRUE`. `predict_logIE` applies
#'   them so the feature space at prediction time matches training. `NULL` for
#'   externally trained models, which are used unscaled.
#' @param .split_strategy,.n_groups Character / integer describing how the
#'   train/test split was made (see [build_logIE_model]). `NULL` for externally
#'   trained models.
#' @param .training_date POSIXct training timestamp.
#'
#' @return An `nts_logIE_model` object.
#' @keywords internal
#' @noRd
new_nts_logIE_model <- function(model, method, feature_set, feature_cols,
                                 n_features, hp, train_metrics, test_metrics,
                                 cv_metrics, polarity,
                                 n_training, n_test,
                                 .best_nrounds = NULL,
                                 .feature_centers = NULL,
                                 .preproc_center = NULL,
                                 .preproc_scale = NULL,
                                 .split_strategy = NULL,
                                 .n_groups = NULL,
                                 .training_date = Sys.time()) {
  structure(list(
    model          = model,
    method         = method,
    feature_set    = feature_set,
    feature_cols   = feature_cols,
    n_features     = n_features,
    hp             = hp,
    best_nrounds   = .best_nrounds,
    train_metrics  = train_metrics,
    test_metrics   = test_metrics,
    cv_metrics     = cv_metrics,
    polarity       = match.arg(polarity, c("positive", "negative")),
    n_training     = n_training,
    n_test         = n_test,
    feature_centers = .feature_centers,
    preproc_center = .preproc_center,
    preproc_scale  = .preproc_scale,
    split_strategy = .split_strategy,
    n_groups       = .n_groups,
    training_date  = .training_date,
    version        = "0.1.0"
  ), class = "nts_logIE_model")
}


#' @export
print.nts_logIE_model <- function(x, ...) {
  cat("nts_logIE_model (v", x$version, ")\n", sep = "")
  cat("  Polarity:  ", x$polarity, "\n")
  cat("  Method:    ", x$method, "\n")
  cat("  Features:  ", x$n_features, " (", x$feature_set, ")\n", sep = "")
  cat("  Training:  ", x$n_training, " compounds\n", sep = "")
  if (!is.null(x$train_metrics)) {
    cat("  Train R2:  ", round(x$train_metrics$R2, 4),
        "  RMSE: ", round(x$train_metrics$RMSE, 4),
        "  MAE: ", round(x$train_metrics$MAE, 4), "\n", sep = "")
  }
  cat("  Test R2:   ", round(x$test_metrics$R2, 4),
      "  RMSE: ", round(x$test_metrics$RMSE, 4),
      "  MAE: ", round(x$test_metrics$MAE, 4), "\n", sep = "")
  cat("  CV R2:     ", round(x$cv_metrics$R2, 4),
      "  RMSE: ", round(x$cv_metrics$RMSE, 4), "\n", sep = "")
  if (!is.null(x$split_strategy)) {
    cat("  Split:     ", x$split_strategy,
        if (!is.null(x$n_groups)) paste0(" (", x$n_groups, " groups)") else "",
        "\n", sep = "")
  }
  if (!is.null(x$preproc_center)) {
    cat("  Scaled:     yes (", length(x$preproc_center),
        " training means applied at predict time)\n", sep = "")
  }
  cat("  Trained:   ", format(x$training_date), "\n", sep = "")
  invisible(x)
}


# ---- metadata matrix ----------------------------------------------------------

#' Assemble the metadata columns (RT / MW / Organic_Percentage) for prediction
#'
#' Builds the subset of metadata columns the model was trained on, naming each
#' column exactly as it appears in `feature_cols`. The RT column in particular
#' is not always called `"RT_avg"`: the POS in-house model was trained with a
#' column named `"RT"`. Hard-coding `"RT_avg"` there caused the model's `"RT"`
#' feature to be treated as missing and silently filled with 0.
#'
#' @param fcols Model `feature_cols`.
#' @param rt,mw,organic_pct Numeric vectors (or `NULL` to omit the column).
#' @return A matrix of the metadata columns present in `fcols`, or `NULL`.
#' @keywords internal
#' @noRd
.nts_logIE_meta_matrix <- function(fcols, rt = NULL, mw = NULL, organic_pct = NULL) {
  meta_mat <- NULL
  if (any(c("RT", "RT_avg") %in% fcols) && !is.null(rt)) {
    rt_names <- intersect(c("RT", "RT_avg"), fcols)
    if (!length(rt_names)) rt_names <- "RT_avg"
    for (rn in rt_names) {
      meta_mat <- cbind(meta_mat,
                        matrix(as.numeric(rt), ncol = 1, dimnames = list(NULL, rn)))
    }
  }
  if ("MW" %in% fcols && !is.null(mw)) {
    meta_mat <- cbind(meta_mat,
                      matrix(as.numeric(mw), ncol = 1, dimnames = list(NULL, "MW")))
  }
  if ("Organic_Percentage" %in% fcols && !is.null(organic_pct)) {
    meta_mat <- cbind(meta_mat,
                      matrix(as.numeric(organic_pct), ncol = 1,
                             dimnames = list(NULL, "Organic_Percentage")))
  }
  meta_mat
}


#' Molecular weight from SMILES for the `MW` model feature
#'
#' The bundled logIE models use `MW` (molecular weight) as a feature. MW is a
#' pure function of the structure, so it must never be imputed with 0 -- it can
#' always be computed from the SMILES. Reuses the RCDK descriptor block when the
#' model already asked for one, otherwise computes just the weight descriptor.
#'
#' @param smiles Character vector of SMILES.
#' @param desc_dt Optional data.table of already-computed RCDK descriptors
#'   (must carry `smiles` and `MW` columns).
#' @param BPPARAM A BiocParallelParam object.
#' @return Numeric vector aligned to `smiles`, `NA` where unavailable, or `NULL`
#'   when MW could not be obtained at all.
#' @keywords internal
#' @noRd
.nts_logIE_derive_mw <- function(smiles, desc_dt = NULL,
                                  BPPARAM = BiocParallel::SerialParam()) {
  smi <- trimws(as.character(smiles))
  if (!length(smi)) return(NULL)
  src <- NULL
  if (!is.null(desc_dt) && all(c("smiles", "MW") %in% names(desc_dt))) {
    src <- desc_dt
  } else if (requireNamespace("rcdk", quietly = TRUE)) {
    src <- tryCatch(
      .nts_compute_descriptors(
        smi,
        desc_names = "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor",
        BPPARAM = BPPARAM, prune = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(src) && !all(c("smiles", "MW") %in% names(src))) src <- NULL
  }
  if (is.null(src)) return(NULL)
  mw <- suppressWarnings(as.numeric(src$MW))
  idx <- match(smi, as.character(src$smiles))
  out <- rep(NA_real_, length(smi))
  hit <- !is.na(idx)
  out[hit] <- mw[idx[hit]]
  if (!any(is.finite(out))) return(NULL)
  out
}


# ---- predict_logIE ------------------------------------------------------------

#' Predict logIE from SMILES
#'
#' Computes required fingerprints and/or descriptors, then returns
#' predicted log10 ionization efficiency for each input SMILES.
#'
#' @param model An `nts_logIE_model` object.
#' @param smiles Character vector of SMILES to predict.
#' @param polarity `"auto"` (default) uses the model's stored polarity. Passing a
#'   polarity that contradicts `model$polarity` is an error: one model object
#'   serves one polarity, so load the matching model instead.
#' @param fingerprint_matrix Optional pre-computed fingerprint matrix
#'   (use when calling repeatedly on the same compounds).
#' @param descriptor_matrix Optional pre-computed descriptor matrix.
#' @param rt Optional numeric vector of retention times, used when the
#'   model includes RT as a feature. Aligned to `smiles` or to its valid subset.
#' @param mw Optional numeric vector of molecular weights. When the model uses
#'   `MW` and this is `NULL`, the value is computed from the SMILES (it is a pure
#'   function of the structure), never imputed with 0.
#' @param organic_pct Optional numeric vector of organic percentages. This is a
#'   run-level (LC gradient) value that cannot be derived from the structure, so a
#'   model that requires it must be given it explicitly.
#' @param gradient_profile Optional list with numeric vectors \code{rt} and
#'   \code{organic} (e.g. \code{list(rt = c(0, 15, 20), organic = c(5, 95, 95))})
#'   used to interpolate \code{Organic_Percentage} from \code{rt} via
#'   \code{predict_organic_percentage} when the model requires it and
#'   \code{organic_pct} is not explicitly supplied.
#' @param BPPARAM A BiocParallelParam object.
#'
#' @details When the caller cannot supply a feature column it is imputed with the
#'   training-set mean stored in `model$feature_centers`, or with 0 when the model
#'   records no means; both cases raise a warning naming the columns.
#'
#' @return A data.frame with column `predicted_logIE`.
#' @export
predict_logIE <- function(model,
                           smiles,
                           polarity = c("auto", "positive", "negative"),
                           fingerprint_matrix = NULL,
                           descriptor_matrix = NULL,
                           rt = NULL,
                           mw = NULL,
                           organic_pct = NULL,
                           gradient_profile = NULL,
                           BPPARAM = BiocParallel::SerialParam()) {
  if (!inherits(model, "nts_logIE_model")) {
    stop("model must be an nts_logIE_model object.", call. = FALSE)
  }

  polarity <- match.arg(polarity)
  if (polarity != "auto" && polarity != model$polarity) {
    # A single model object serves exactly one polarity, so a contradicting
    # `polarity` cannot be honoured. Refusing loudly beats silently returning
    # predictions from the wrong-polarity model.
    stop("Requested polarity '", polarity, "' does not match the polarity of the ",
         "supplied model ('", model$polarity, "'). Load the matching model ",
         "instead (e.g. logIE_model_inhouse_",
         if (polarity == "positive") "pos" else "neg", ".rds).", call. = FALSE)
  }

  # The model backend must be installed for prediction to work.
  backend_pkg <- if (identical(model$method, "xgb")) "xgboost" else "ranger"
  if (!requireNamespace(backend_pkg, quietly = TRUE)) {
    stop(sprintf("Model backend '%s' is required to predict with a '%s' model. Install: install.packages('%s')",
                 backend_pkg, model$method, backend_pkg), call. = FALSE)
  }

  smi <- trimws(as.character(smiles))
  fcols <- model$feature_cols

  # Determine what needs computing
  need_extended   <- any(grepl("^extended_", fcols))
  need_maccs      <- any(grepl("^maccs_", fcols))
  need_circular   <- any(grepl("^circular_", fcols))
  need_standard   <- any(grepl("^standard_", fcols))
  need_pubchem    <- any(grepl("^pubchem_", fcols))
  need_substruct  <- any(grepl("^substructure_", fcols))
  need_rcdk_desc  <- any(grepl("^rcdk_desc_", fcols))
  need_rdkit_desc <- any(grepl("^rdkit_desc_", fcols))
  need_mordred    <- any(grepl("^mordred_", fcols))
  need_rt         <- "RT_avg" %in% fcols || "RT" %in% fcols
  need_mw_col     <- "MW" %in% fcols
  need_org_pct    <- "Organic_Percentage" %in% fcols

  any_fp <- need_extended || need_maccs || need_circular || need_standard ||
            need_pubchem || need_substruct
  any_desc     <- need_rcdk_desc
  any_ext_desc <- need_rdkit_desc || need_mordred

  # Fingerprints
  if (!is.null(fingerprint_matrix)) {
    if (!is.matrix(fingerprint_matrix)) stop("fingerprint_matrix must be a matrix.", call. = FALSE)
    fp_mat <- fingerprint_matrix
  } else if (any_fp) {
    fp_types <- c(
      if (need_extended) "extended",
      if (need_maccs) "maccs",
      if (need_circular) "circular",
      if (need_standard) "standard",
      if (need_pubchem) "pubchem",
      if (need_substruct) "substructure"
    )
    if (!length(smi) || all(is.na(smi) | !nzchar(smi))) {
      stop("No valid SMILES strings provided.", call. = FALSE)
    }
    fp_mat <- compute_fingerprints(smi, types = fp_types, fp_mode = "matrix", BPPARAM = BPPARAM)
  } else {
    fp_mat <- NULL
  }

  # Descriptors (RCDK)
  desc_dt <- NULL
  if (!is.null(descriptor_matrix)) {
    if (!is.matrix(descriptor_matrix)) stop("descriptor_matrix must be a matrix.", call. = FALSE)
    desc_mat <- descriptor_matrix
  } else if (any_desc) {
    if (!length(smi) || all(is.na(smi) | !nzchar(smi))) {
      stop("No valid SMILES strings provided.", call. = FALSE)
    }
    desc_dt <- .nts_compute_descriptors(smi, BPPARAM = BPPARAM, prune = FALSE)
    desc_cols <- grep("^rcdk_desc_", fcols, value = TRUE)
    clean_cols <- sub("^rcdk_desc_", "", desc_cols)
    avail <- intersect(clean_cols, names(desc_dt))
    if (!length(avail)) stop("None of the required RCDK descriptor columns could be computed.", call. = FALSE)
    desc_mat <- as.matrix(desc_dt[, avail, with = FALSE])
    colnames(desc_mat) <- paste0("rcdk_desc_", avail)
  } else {
    desc_mat <- NULL
  }

  # External descriptors: RDKit (via reticulate)
  rdkit_mat <- NULL
  if (need_rdkit_desc && is.null(descriptor_matrix)) {
    rdkit_cols <- grep("^rdkit_desc_", fcols, value = TRUE)
    rdkit_raw  <- sub("^rdkit_desc_", "", rdkit_cols)
    rdkit_raw  <- setdiff(rdkit_raw, c("RT", "Organic_Percentage", "MW"))
    rdkit_dt   <- .nts_compute_rdkit_descriptors(smi, rdkit_raw)
    if (!is.null(rdkit_dt) && ncol(rdkit_dt) > 0L) {
      rdkit_mat <- as.matrix(rdkit_dt)
      colnames(rdkit_mat) <- paste0("rdkit_desc_", colnames(rdkit_dt))
    }
  }

  # External descriptors: Mordred (via reticulate)
  mordred_mat <- NULL
  if (need_mordred && is.null(descriptor_matrix)) {
    mordred_cols <- grep("^mordred_", fcols, value = TRUE)
    mordred_raw  <- sub("^mordred_", "", mordred_cols)
    mordred_raw  <- setdiff(mordred_raw, c("RT", "Organic_Percentage", "MW"))
    mordred_dt   <- .nts_compute_mordred_descriptors(smi, mordred_raw)
    if (!is.null(mordred_dt) && ncol(mordred_dt) > 0L) {
      mordred_mat <- as.matrix(mordred_dt)
      colnames(mordred_mat) <- paste0("mordred_", colnames(mordred_dt))
    }
  }

  # Warn if external descriptors unavailable
  ext_missing <- character(0)
  if (need_rdkit_desc && is.null(rdkit_mat) && is.null(descriptor_matrix))
    ext_missing <- c(ext_missing, "RDKit")
  if (need_mordred && is.null(mordred_mat) && is.null(descriptor_matrix))
    ext_missing <- c(ext_missing, "Mordred")
  if (length(ext_missing)) {
    warning("Model requires ", paste(ext_missing, collapse = "/"),
            " descriptors which are not available. ",
            "Install the Python packages (pip install rdkit mordred) and ",
            "reticulate R package, or supply them via descriptor_matrix. ",
            "Missing columns will be imputed from the model's training means ",
            "and predictions will be less reliable.", call. = FALSE)
  }

  # Metadata columns. Every feature block describes the same rows -- the *valid*
  # SMILES -- so a metadata vector may be aligned either to `smiles` or to that
  # valid subset.
  valid_idx <- which(!is.na(smi) & nzchar(smi))
  block_rows <- c(
    if (!is.null(fp_mat)) nrow(fp_mat),
    if (!is.null(desc_mat)) nrow(desc_mat),
    if (!is.null(rdkit_mat)) nrow(rdkit_mat),
    if (!is.null(mordred_mat)) nrow(mordred_mat)
  )
  if (length(unique(block_rows)) > 1L) {
    stop("Feature blocks have inconsistent row counts (",
         paste(unique(block_rows), collapse = " vs "),
         "); some SMILES most likely failed descriptor computation. ",
         "Remove those SMILES and retry.", call. = FALSE)
  }
  n_smi <- if (length(block_rows)) block_rows[1L] else length(smi)

  .norm_meta <- function(v, label) {
    if (is.null(v)) return(NULL)
    v <- suppressWarnings(as.numeric(v))
    if (length(v) == n_smi) return(v)
    if (length(v) == length(smi) && n_smi == length(valid_idx)) return(v[valid_idx])
    stop("`", label, "` has ", length(v), " value(s) but the feature matrix has ",
         n_smi, " row(s) (input had ", length(smi), " SMILES).", call. = FALSE)
  }
  rt          <- .norm_meta(rt, "rt")
  mw          <- .norm_meta(mw, "mw")
  # If model needs Organic_Percentage and gradient_profile is provided, interpolate from rt
  if (need_org_pct && is.null(organic_pct) && !is.null(rt) && !is.null(gradient_profile)) {
    if (is.list(gradient_profile) && all(c("rt", "organic") %in% names(gradient_profile))) {
      organic_pct <- tryCatch(
        predict_organic_percentage(rt, ref_rt = gradient_profile$rt, ref_org = gradient_profile$organic),
        error = function(e) NULL
      )
    }
  }
  organic_pct <- .norm_meta(organic_pct, "organic_pct")

  # MW is a pure function of the structure: when the model uses it and the caller
  # supplied nothing, derive it from the descriptor block (already computed for
  # both bundled models) instead of imputing 0 Da.
  if (need_mw_col && is.null(mw)) {
    mw <- .norm_meta(.nts_logIE_derive_mw(smi, desc_dt = desc_dt, BPPARAM = BPPARAM), "mw")
  }

  meta_mat <- .nts_logIE_meta_matrix(fcols, rt = rt, mw = mw, organic_pct = organic_pct)

  # Combine
  X_parts <- list()
  if (!is.null(fp_mat)) X_parts <- c(X_parts, list(fp_mat))
  if (!is.null(desc_mat)) X_parts <- c(X_parts, list(desc_mat))
  if (!is.null(rdkit_mat)) X_parts <- c(X_parts, list(rdkit_mat))
  if (!is.null(mordred_mat)) X_parts <- c(X_parts, list(mordred_mat))
  if (!is.null(meta_mat)) X_parts <- c(X_parts, list(meta_mat))
  if (!length(X_parts)) stop("No features could be assembled for prediction.", call. = FALSE)

  X <- do.call(cbind, X_parts)

  # Align to model feature order. Any column the caller could not supply is
  # imputed with the training-set mean recorded in the model -- never silently
  # with 0 (a 0 in a fingerprint bit or a descriptor is a *real* value there,
  # so 0-filling produces confident-looking but wrong predictions).
  missing_cols <- setdiff(fcols, colnames(X))
  if (length(missing_cols)) {
    centers <- model$feature_centers
    fill <- if (is.null(centers)) rep(0, length(missing_cols))
            else unname(centers[missing_cols])
    fill[is.na(fill) | !is.finite(fill)] <- 0
    for (i in seq_along(missing_cols)) {
      X <- cbind(X, matrix(fill[i], nrow = nrow(X), ncol = 1,
                           dimnames = list(NULL, missing_cols[i])))
    }
    # Metadata columns are run-level inputs the caller can always supply, unlike
    # a fingerprint bit; name the argument that would fix each one.
    meta_missing <- intersect(missing_cols, c("RT", "RT_avg", "MW", "Organic_Percentage"))
    hint <- ""
    if (length(meta_missing)) {
      arg_map <- c(RT = "rt", RT_avg = "rt", MW = "mw", Organic_Percentage = "organic_pct")
      hint <- paste0(" Supply the run-level column(s) via ",
                     paste(unique(sprintf("`%s=`", arg_map[meta_missing])), collapse = " / "),
                     " (or add them to the input table).")
    }
    warning(
      "Model requires ", length(missing_cols), " feature column(s) that could ",
      "not be computed: ",
      paste(utils::head(missing_cols, 8L), collapse = ", "),
      if (length(missing_cols) > 8L) ", ..." else "",
      ". They were imputed with ",
      if (is.null(centers)) "0 (the model stores no training means)"
      else "the training-set mean (0 where unavailable)",
      "; predictions for these compounds are less reliable.", hint,
      call. = FALSE
    )
  }
  X <- X[, fcols, drop = FALSE]

  # Models built by build_logIE_model(feature_screen = TRUE) were trained on
  # standardised features; apply the *training* constants (never refit them) so
  # the prediction space matches. Externally trained models carry no `preproc_*`
  # fields and are used as-is. Imputed columns are filled with the raw training
  # mean above, so after centring they become exactly 0 -- the neutral value.
  if (!is.null(model$preproc_center) && length(model$preproc_center) == ncol(X)) {
    X <- scale(X, center = model$preproc_center, scale = model$preproc_scale)
  }

  # Predict
  pred <- if (model$method == "rf") {
    as.numeric(predict(model$model, X)$predictions)
  } else {
    as.numeric(predict(model$model, xgboost::xgb.DMatrix(X)))
  }

  data.frame(predicted_logIE = pred, stringsAsFactors = FALSE)
}


# ---- build_logIE_model --------------------------------------------------------

#' Build a logIE prediction model
#'
#' Trains an XGBoost or Random Forest model for logIE prediction.
#' This is a standardised wrapper; the heavy HP tuning workflow is
#' performed externally (see tests/logIEpred/).
#'
#' @param smiles Character vector of SMILES.
#' @param logIE Numeric vector of observed logIE values.
#' @param fingerprint_types Fingerprint types passed to
#'   \code{compute_fingerprints}.
#' @param method `"xgb"` or `"rf"`.
#' @param tune_params Named list of hyperparameters.
#' @param train_split Proportion for training (default 0.8).
#' @param group Optional grouping vector (e.g. chemical class / homologue series)
#'   the same length as the input compounds. When supplied, whole groups are held
#'   out together. When `NULL` groups are derived automatically by clustering the
#'   fingerprints (see `cluster_similarity`), mirroring [build_rt_model].
#' @param cluster_similarity Tanimoto similarity at which automatically derived
#'   groups are cut (default 0.8). Only used when `group = NULL`.
#' @param feature_screen If `TRUE` (default), screen the features on the
#'   *training rows only* (near-zero variance + highly correlated columns) and
#'   record the standardisation constants in the model, exactly as
#'   [build_rt_model] does. This guards against the overfitting a wide,
#'   unscreened matrix invites (e.g. 4617 features on ~480 compounds) and keeps
#'   the reported test metrics honest. The screening is learned after the split,
#'   so the test rows never influence it.
#' @param seed Random seed.
#' @param BPPARAM A BiocParallelParam object.
#' @param metadata_matrix Optional pre-computed metadata matrix (e.g. carrying
#'   RT, MW, Organic_Percentage).
#' @param polarity `"positive"` (default) or `"negative"`.
#'
#' @return An `nts_logIE_model` object.
#' @export
build_logIE_model <- function(smiles,
                               logIE,
                               fingerprint_types = "base",
                               method = c("xgb", "rf"),
                               tune_params = list(),
                               train_split = 0.8,
                               group = NULL,
                               cluster_similarity = 0.8,
                               feature_screen = TRUE,
                               seed = 42,
                               BPPARAM = BiocParallel::SerialParam(),
                               fingerprint_matrix = NULL,
                               descriptor_matrix = NULL,
                               metadata_matrix = NULL,
                               polarity = c("positive", "negative")) {
  method <- match.arg(method)
  polarity <- match.arg(polarity)
  set.seed(seed)

  backend_pkg <- if (identical(method, "xgb")) "xgboost" else "ranger"
  if (!requireNamespace(backend_pkg, quietly = TRUE)) {
    stop(sprintf("Model backend '%s' is required to build a '%s' model. Install: install.packages('%s')",
                 backend_pkg, method, backend_pkg), call. = FALSE)
  }

  has_fp <- !is.null(fingerprint_matrix)
  has_desc <- !is.null(descriptor_matrix)

  if (has_fp || has_desc) {
    smi <- rep("precomputed", if (has_fp) nrow(fingerprint_matrix) else nrow(descriptor_matrix))
    y_vec <- suppressWarnings(as.numeric(logIE))
    ok <- !is.na(y_vec) & is.finite(y_vec)
    if (!any(ok)) stop("No valid logIE values.", call. = FALSE)
    y_vec <- y_vec[ok]
    fp_mat <- if (has_fp) fingerprint_matrix[ok, , drop = FALSE] else NULL
    desc_mat <- if (has_desc) descriptor_matrix[ok, , drop = FALSE] else NULL
    n <- length(y_vec)
  } else {
    smi <- trimws(as.character(smiles))
    y_vec <- suppressWarnings(as.numeric(logIE))
    if (length(smi) != length(y_vec)) {
      stop("smiles and logIE must have the same length.", call. = FALSE)
    }
    ok <- !is.na(smi) & nzchar(smi) & !is.na(y_vec) & is.finite(y_vec)
    if (!any(ok)) stop("No valid (smiles, logIE) pairs.", call. = FALSE)
    smi <- smi[ok]; y_vec <- y_vec[ok]; n <- length(smi)
    fp_mat <- compute_fingerprints(smi, types = fingerprint_types,
                                    fp_mode = "matrix", BPPARAM = BPPARAM)
    desc_mat <- NULL
  }

  # Build feature matrix
  if (!is.null(fp_mat) && !is.null(desc_mat)) {
    X <- cbind(fp_mat, desc_mat)
  } else if (!is.null(fp_mat)) {
    X <- fp_mat
  } else if (!is.null(desc_mat)) {
    X <- desc_mat
  } else {
    X <- NULL
  }
  if (!is.null(metadata_matrix)) {
    if (!is.matrix(metadata_matrix)) metadata_matrix <- as.matrix(metadata_matrix)
    X <- if (!is.null(X)) cbind(X, metadata_matrix) else metadata_matrix
  }
  if (is.null(X)) {
    stop("No features provided.", call. = FALSE)
  }

  # ---- Train/test split ------------------------------------------------------
  # Group-aware by default, mirroring build_rt_model: a random split scatters
  # structurally related analogues (homologues, isomers) across both sides, so
  # the test error is optimistic. Groups come from the caller, else from
  # average-linkage clustering of the fingerprints. The grouped splitter is only
  # used when the grouping is informative (2..n/2 distinct groups); otherwise we
  # fall back to the logIE-stratified split, which guarantees that the observed
  # logIE range is represented on both sides.
  grp <- NULL
  handled <- FALSE
  tr_idx <- integer(0); te_idx <- integer(0)
  if (!is.null(group)) {
    grp <- .nts_rt_resolve_groups(group, n = n, fp = NULL,
                                  cluster_similarity = cluster_similarity, seed = seed)
  } else {
    grp <- .nts_rt_resolve_groups(NULL, n = n, fp = fp_mat,
                                  cluster_similarity = cluster_similarity, seed = seed)
  }
  n_grp <- length(unique(grp))
  if (n_grp >= 2L && n_grp <= floor(n / 2L)) {
    sp <- .nts_rt_group_split(grp, n = n, train_split = train_split, seed = seed)
    if (length(sp$train_idx) && length(sp$test_idx)) {
      tr_idx <- sp$train_idx; te_idx <- sp$test_idx
      split_strategy <- sp$strategy
      n_groups_rep <- n_grp
      handled <- TRUE
    }
  }
  if (!handled) {
    # Stratified split on the logIE value itself.
    bins <- cut(y_vec, breaks = 10, labels = FALSE)
    bins[is.na(bins)] <- 1L
    set.seed(seed)
    for (b in seq_len(max(bins, na.rm = TRUE))) {
      bi <- which(bins == b)
      if (length(bi) > 1) {
        te_idx <- c(te_idx, sample(bi, max(1L, round(length(bi) * (1 - train_split)))))
      }
    }
    tr_idx <- setdiff(seq_len(n), te_idx)
    split_strategy <- "stratified"
    n_groups_rep <- NA_integer_
  }

  # ---- Feature screening -----------------------------------------------------
  # Fitted on the training rows only (see build_rt_model for the same rule); the
  # constants are stored so predict_logIE applies the identical transform.
  X_pre <- X
  preproc_center <- NULL
  preproc_scale <- NULL
  if (isTRUE(feature_screen)) {
    prep <- .nts_rt_preprocess(X, fit_idx = tr_idx)
    X <- prep$matrix
    preproc_center <- prep$center
    preproc_scale  <- prep$scale
    if (!is.null(preproc_center)) names(preproc_center) <- colnames(X)
    if (!is.null(preproc_scale))  names(preproc_scale)  <- colnames(X)
  }

  p_features <- ncol(X)
  feature_cols <- colnames(X)

  X_tr <- X[tr_idx, , drop = FALSE]; y_tr <- y_vec[tr_idx]
  X_te <- X[te_idx, , drop = FALSE]; y_te <- y_vec[te_idx]

  # Train
  if (method == "xgb") {
    obj <- if (!is.null(tune_params$objective)) tune_params$objective else "reg:squarederror"
    params <- tune_params
    params$objective <- obj
    if (is.null(params$nthread)) params$nthread <- 1L
    if (is.null(params$eta)) params$eta <- 0.05
    if (is.null(params$max_depth)) params$max_depth <- 6L
    dtr <- xgboost::xgb.DMatrix(X_tr, label = y_tr)
    dte <- xgboost::xgb.DMatrix(X_te, label = y_te)

    # Early stopping must not watch the held-out test set: choosing the number
    # of rounds on `dte` leaks test information into the model and makes the
    # reported test metrics optimistic. Carve an inner validation split out of
    # the *training* rows instead, and touch `dte` only once, at the end.
    n_tr <- nrow(X_tr)
    n_eval <- max(1L, floor(n_tr * 0.2))
    if (n_tr - n_eval >= 5L) {
      ev_idx <- sample.int(n_tr, n_eval)
      fit_idx <- setdiff(seq_len(n_tr), ev_idx)
      dtr_fit <- xgboost::xgb.DMatrix(X_tr[fit_idx, , drop = FALSE], label = y_tr[fit_idx])
      dva <- xgboost::xgb.DMatrix(X_tr[ev_idx, , drop = FALSE], label = y_tr[ev_idx])
      bst <- xgboost::xgb.train(
        params = params, data = dtr_fit, nrounds = 500L, verbose = 0L,
        watchlist = list(train = dtr_fit, valid = dva),
        early_stopping_rounds = 20L
      )
      best_nr <- if (is.null(bst$best_iteration) || bst$best_iteration < 1) bst$niter else bst$best_iteration
      fit <- xgboost::xgb.train(params = params, data = dtr_fit, nrounds = best_nr, verbose = 0L)
    } else {
      # Too few training rows to hold out a validation split; fit on everything
      # with the full round budget rather than watching the test set.
      best_nr <- 500L
      fit <- xgboost::xgb.train(params = params, data = dtr, nrounds = best_nr, verbose = 0L)
    }

    if (requireNamespace("xgboost", quietly = TRUE)) {
      fit <- xgboost::xgb.Booster.complete(fit)
    }
    pr_tr <- as.numeric(predict(fit, dtr))
    pr_te <- as.numeric(predict(fit, dte))
    model_obj <- fit
  } else {
    params <- c(tune_params, list(num.threads = 1L, verbose = FALSE))
    if (is.null(params$num.trees)) params$num.trees <- 500L
    if (is.null(params$mtry)) params$mtry <- max(1L, floor(p_features / 3))
    if (is.null(params$min.node.size)) params$min.node.size <- 5L

    fit <- ranger::ranger(x = X_tr, y = y_tr,
                           num.trees = params$num.trees,
                           mtry = params$mtry,
                           min.node.size = params$min.node.size,
                           num.threads = params$num.threads,
                           verbose = params$verbose)
    pr_tr <- as.numeric(predict(fit, X_tr)$predictions)
    pr_te <- as.numeric(predict(fit, X_te)$predictions)
    model_obj <- fit
    best_nr <- NULL
  }

  train_r2 <- 1 - sum((y_tr - pr_tr)^2) / sum((y_tr - mean(y_tr))^2)
  train_rmse <- sqrt(mean((y_tr - pr_tr)^2))
  train_mae <- mean(abs(y_tr - pr_tr))
  test_r2 <- 1 - sum((y_te - pr_te)^2) / sum((y_te - mean(y_te))^2)
  test_rmse <- sqrt(mean((y_te - pr_te)^2))
  test_mae <- mean(abs(y_te - pr_te))

  hp_str <- if (method == "xgb") {
    sprintf("eta=%.3f depth=%d", params$eta, params$max_depth)
  } else {
    sprintf("trees=%d mtry=%d", params$num.trees, params$mtry)
  }

  # Column means of the training rows: used by predict_logIE to impute any
  # required feature it cannot compute, instead of silently using 0.
  # NOTE: measured on the *raw* (pre-screening) scale, because predict_logIE
  # fills the missing raw values first and standardises afterwards -- a mean of
  # the already-standardised matrix would be ~0 and would not cancel out.
  feature_centers <- colMeans(X_pre[tr_idx, feature_cols, drop = FALSE], na.rm = TRUE)
  feature_centers[!is.finite(feature_centers)] <- 0

  new_nts_logIE_model(
    model         = model_obj,
    method        = method,
    feature_set   = "user",
    feature_cols  = feature_cols,
    n_features    = p_features,
    hp            = hp_str,
    train_metrics = list(R2 = train_r2, RMSE = train_rmse, MAE = train_mae),
    test_metrics  = list(R2 = test_r2, RMSE = test_rmse, MAE = test_mae),
    cv_metrics    = list(R2 = NA_real_, RMSE = NA_real_),
    polarity      = polarity,
    n_training    = length(tr_idx),
    n_test        = length(te_idx),
    .best_nrounds = best_nr,
    .feature_centers = feature_centers,
    .preproc_center  = preproc_center,
    .preproc_scale   = preproc_scale,
    .split_strategy  = split_strategy,
    .n_groups        = n_groups_rep
  )
}


# ---- save / load --------------------------------------------------------------

#' Save a logIE model to disk
#' @param model An `nts_logIE_model` object.
#' @param file Path to output .rds file.
#' @return Invisibly, the file path.
#' @export
save_logIE_model <- function(model, file) {
  if (!inherits(model, "nts_logIE_model")) {
    stop("model must be an nts_logIE_model object.", call. = FALSE)
  }
  if (identical(model$method, "xgb") && requireNamespace("xgboost", quietly = TRUE)) {
    if (inherits(model$model, "xgb.Booster")) {
      model$model <- xgboost::xgb.Booster.complete(model$model)
    }
  }
  saveRDS(model, file)
  invisible(file)
}


#' Load a logIE model from disk
#' @param file Path to .rds file.
#' @return An `nts_logIE_model` object.
#' @export
load_logIE_model <- function(file) {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  model <- readRDS(file)

  # Handle both raw model structures (from external training) and
  # nts_logIE_model S3 objects
  if (inherits(model, "nts_logIE_model")) {
    if (!model$version %in% c("0.1.0", "0.2.0")) {
      warning("Model version ", model$version, " may not be compatible.", call. = FALSE)
    }
    return(model)
  }

  # Up-convert a raw list model (from external training scripts) to S3
  if (is.list(model) && !is.null(model$model) && !is.null(model$model_type)) {
    return(.upgrade_raw_logIE_model(model))
  }

  stop("File does not contain a recognised nts_logIE_model.", call. = FALSE)
}


#' Upgrade a raw list logIE model to the S3 class
#' @param m A list with model, model_type, feature_set, etc.
#' @return An `nts_logIE_model` object.
#' @keywords internal
#' @noRd
.upgrade_raw_logIE_model <- function(m) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  new_nts_logIE_model(
    model         = m$model,
    method        = if (m$model_type == "xgboost") "xgb" else "rf",
    feature_set   = m$feature_set %||% "unknown",
    feature_cols  = m$feature_cols %||% character(0),
    n_features    = m$n_features %||% 0L,
    hp            = m$hp %||% "",
    train_metrics = m$train_metrics %||% list(R2 = NA, RMSE = NA, MAE = NA),
    test_metrics  = m$test_metrics %||% list(R2 = NA, RMSE = NA, MAE = NA),
    cv_metrics    = m$cv_metrics %||% list(R2 = NA, RMSE = NA),
    polarity      = m$polarity %||% "positive",
    n_training    = m$n_training %||% 0L,
    n_test        = m$n_test %||% 0L,
    .best_nrounds = m$best_nrounds %||% NULL,
    .feature_centers = m$feature_centers %||% NULL,
    .preproc_center = m$preproc_center %||% NULL,
    .preproc_scale  = m$preproc_scale %||% NULL,
    .split_strategy = m$split_strategy %||% NULL,
    .n_groups       = m$n_groups %||% NULL,
    .training_date = m$training_date %||% Sys.time()
  )
}


# ---- evaluate_logIE_model -----------------------------------------------------

#' Evaluate a logIE model on external test data
#'
#' @param model An `nts_logIE_model`.
#' @param test_smiles Character vector of test SMILES.
#' @param test_logIE Numeric vector of observed logIE.
#' @param ... Additional arguments passed to \code{predict_logIE}.
#'
#' @return A list with `metrics` (named numeric), `predicted`, and `observed`.
#' @export
evaluate_logIE_model <- function(model, test_smiles, test_logIE, ...) {
  if (!inherits(model, "nts_logIE_model")) {
    stop("model must be an nts_logIE_model object.", call. = FALSE)
  }
  obs <- suppressWarnings(as.numeric(test_logIE))
  # predict_logIE() drops NA / empty SMILES rows (via compute_fingerprints), so
  # its output can be shorter than the input. Re-expand it against the same
  # validity mask before comparing, otherwise predicted/observed are silently
  # misaligned.
  smi_ok <- !is.na(test_smiles) & nzchar(trimws(as.character(test_smiles)))
  pred <- predict_logIE(model, test_smiles, ...)
  pli <- pred$predicted_logIE
  predicted <- rep(NA_real_, length(smi_ok))
  if (length(pli) == sum(smi_ok)) {
    predicted[smi_ok] <- pli
  } else {
    # Copying positionally here would pair each prediction with the wrong
    # compound as soon as one SMILES is invalid, silently corrupting the
    # metrics. Leave the unmatched rows NA instead and say so (finding P2-3).
    warning("evaluate_logIE_model: predict_logIE() returned ", length(pli),
            " predictions for ", sum(smi_ok), " valid SMILES, so they cannot be ",
            "matched to the observations; affected rows are reported as NA and ",
            "excluded from the metrics.", call. = FALSE)
  }
  ok <- !is.na(obs) & is.finite(obs) & !is.na(predicted)
  yt <- obs[ok]; yp <- predicted[ok]
  ssr <- sum((yt - yp)^2); sst <- sum((yt - mean(yt))^2)
  r2 <- if (sst > 0) 1 - ssr / sst else NA
  rmse <- sqrt(mean((yt - yp)^2))
  mae <- mean(abs(yt - yp))
  list(
    metrics   = c(R2 = r2, RMSE = rmse, MAE = mae),
    predicted = predicted,
    observed  = obs
  )
}

# ---- External descriptor computation (reticulate) ------------------------------

# Helper: configure reticulate Python
.nts_configure_python <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
  if (reticulate::py_available(initialize = FALSE)) return(TRUE)
  # Resolve a Python without hard-coding a machine-specific path: an explicit
  # override wins, then `PATH`, then common conda roots derived from env vars.
  py_name <- if (.Platform$OS.type == "windows") "python.exe" else "bin/python"
  home <- Sys.getenv("USERPROFILE", unset = "")
  if (!nzchar(home)) home <- Sys.getenv("HOME", unset = "")
  conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  candidates <- c(
    Sys.getenv("NTS_PYTHON", unset = ""),
    Sys.getenv("RETICULATE_PYTHON", unset = ""),
    unname(Sys.which(c("python", "python3", "python.exe"))),
    if (nzchar(conda_prefix)) file.path(conda_prefix, py_name) else character(0),
    if (nzchar(home)) file.path(home, c("miniconda3", "anaconda3"), py_name) else character(0)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (py in candidates) {
    if (nzchar(py) && file.exists(py)) {
      tryCatch({
        reticulate::use_python(py, required = FALSE)
        if (reticulate::py_available(initialize = FALSE)) return(TRUE)
      }, error = function(e) NULL)
    }
  }
  FALSE
}

#' Compute RDKit descriptors via reticulate
#'
#' Calls RDKit Python library to compute 2D molecular descriptors.
#' Requires `rdkit` Python package installed:
#' `pip install rdkit` or `conda install -c conda-forge rdkit`.
#'
#' @param smiles Character vector of SMILES.
#' @param wanted_columns Character vector of descriptor names to return
#'   (without prefix). If NULL, returns all available.
#' @return A data.table with descriptor columns, or NULL if RDKit unavailable.
#' @keywords internal
#' @noRd
.nts_compute_rdkit_descriptors <- function(smiles, wanted_columns = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    warning("reticulate required for RDKit descriptors. Install: install.packages('reticulate')", call. = FALSE)
    return(NULL)
  }

  .nts_configure_python()

  rdkit_avail <- tryCatch({
    reticulate::py_module_available("rdkit")
  }, error = function(e) FALSE)

  if (!rdkit_avail) {
    warning("RDKit Python module not found. Install: pip install rdkit", call. = FALSE)
    return(NULL)
  }

  smi <- trimws(as.character(smiles))
  n <- length(smi)
  valid_idx <- which(!is.na(smi) & nzchar(smi))

  # Eager import
  Chem <- reticulate::import("rdkit.Chem")
  Descriptors <- reticulate::import("rdkit.Chem.Descriptors")

  # All available descriptor getter functions (name -> callable)
  getters <- list()

  # Standard descriptors from Descriptors._descList: (name, function) pairs
  if (reticulate::py_has_attr(Descriptors, "_descList")) {
    for (item in Descriptors[["_descList"]]) {
      nm <- reticulate::py_to_r(item[[1L]])
      fn <- item[[2L]]
      getters[[nm]] <- fn
    }
  } else {
    # Fallback: known descriptor list
    known <- c("MaxAbsEStateIndex", "MaxEStateIndex", "MinAbsEStateIndex",
               "MinEStateIndex", "qed", "SPS", "MolWt", "HeavyAtomMolWt",
               "ExactMolWt", "NumValenceElectrons", "MaxPartialCharge",
               "MinPartialCharge", "MaxAbsPartialCharge", "MinAbsPartialCharge",
               "FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3",
               "BCUT2D_MWHI", "BCUT2D_MWLOW", "BCUT2D_CHGHI", "BCUT2D_CHGLO",
               "BCUT2D_LOGPHI", "BCUT2D_LOGPLOW", "BCUT2D_MRHI", "BCUT2D_MRLOW",
               "AvgIpc", "BalabanJ", "BertzCT", "Chi0", "Chi0n", "Chi0v",
               "Chi1", "Chi1n", "Chi1v", "Chi2n", "Chi2v", "Chi3n", "Chi3v",
               "Chi4n", "Chi4v", "HallKierAlpha", "Ipc", "Kappa1", "Kappa2",
               "Kappa3", "LabuteASA", "PEOE_VSA1", "PEOE_VSA10",
               "PEOE_VSA11", "PEOE_VSA12", "PEOE_VSA13", "PEOE_VSA14",
               "PEOE_VSA2", "PEOE_VSA3", "PEOE_VSA4", "PEOE_VSA5",
               "PEOE_VSA6", "PEOE_VSA7", "PEOE_VSA8", "PEOE_VSA9",
               "SMR_VSA1", "SMR_VSA10", "SMR_VSA2", "SMR_VSA3", "SMR_VSA4",
               "SMR_VSA5", "SMR_VSA6", "SMR_VSA7", "SMR_VSA9",
               "SlogP_VSA1", "SlogP_VSA10", "SlogP_VSA11", "SlogP_VSA12",
               "SlogP_VSA2", "SlogP_VSA3", "SlogP_VSA4", "SlogP_VSA5",
               "SlogP_VSA6", "SlogP_VSA7", "SlogP_VSA8", "TPSA",
               "EState_VSA1", "EState_VSA10", "EState_VSA11", "EState_VSA2",
               "EState_VSA3", "EState_VSA4", "EState_VSA5", "EState_VSA6",
               "EState_VSA7", "EState_VSA8", "EState_VSA9",
               "VSA_EState1", "VSA_EState10", "VSA_EState2", "VSA_EState3",
               "VSA_EState4", "VSA_EState5", "VSA_EState6", "VSA_EState7",
               "VSA_EState8", "VSA_EState9", "FractionCSP3", "HeavyAtomCount",
               "NHOHCount", "NOCount", "NumAliphaticCarbocycles",
               "NumAliphaticHeterocycles", "NumAliphaticRings", "NumAmideBonds",
               "NumAromaticCarbocycles", "NumAromaticHeterocycles",
               "NumAromaticRings", "NumAtomStereoCenters", "NumBridgeheadAtoms",
               "NumHAcceptors", "NumHDonors", "NumHeteroatoms",
               "NumHeterocycles", "NumRotatableBonds", "NumSaturatedCarbocycles",
               "NumSaturatedHeterocycles", "NumSaturatedRings", "NumSpiroAtoms",
               "NumUnspecifiedAtomStereoCenters", "Phi", "RingCount",
               "MolLogP", "MolMR")
    for (nm in known) {
      fn <- tryCatch(Descriptors[[nm]], error = function(e) NULL)
      if (!is.null(fn)) getters[[nm]] <- fn
    }
  }

  # Fragment descriptors: fr_*
  frag_funcs <- c("fr_Al_COO", "fr_Al_OH", "fr_Al_OH_noTert", "fr_ArN",
    "fr_Ar_COO", "fr_Ar_N", "fr_Ar_NH", "fr_Ar_OH", "fr_COO", "fr_COO2",
    "fr_C_O", "fr_C_O_noCOO", "fr_C_S", "fr_Imine", "fr_NH0", "fr_NH1",
    "fr_NH2", "fr_N_O", "fr_Ndealkylation1", "fr_Ndealkylation2",
    "fr_Nhpyrrole", "fr_aldehyde", "fr_alkyl_carbamate", "fr_alkyl_halide",
    "fr_allylic_oxid", "fr_amide", "fr_amidine", "fr_aniline",
    "fr_aryl_methyl", "fr_azide", "fr_benzene", "fr_bicyclic", "fr_epoxide",
    "fr_ester", "fr_ether", "fr_furan", "fr_guanido", "fr_halogen",
    "fr_hdrzine", "fr_hdrzone", "fr_imidazole", "fr_imide", "fr_ketone",
    "fr_ketone_Topliss", "fr_lactam", "fr_lactone", "fr_methoxy",
    "fr_morpholine", "fr_nitrile", "fr_nitro", "fr_nitro_arom",
    "fr_nitro_arom_nonortho", "fr_oxazole", "fr_oxime",
    "fr_para_hydroxylation", "fr_phenol", "fr_phenol_noOrthoHbond",
    "fr_phos_acid", "fr_phos_ester", "fr_piperdine", "fr_piperzine",
    "fr_priamide", "fr_pyridine", "fr_sulfide", "fr_sulfonamd", "fr_sulfone",
    "fr_term_acetylene", "fr_tetrazole", "fr_thiazole", "fr_thiophene",
    "fr_unbrch_alkane", "fr_urea")
  for (fnm in frag_funcs) {
    fn <- tryCatch(Descriptors[[fnm]], error = function(e) NULL)
    if (!is.null(fn)) getters[[fnm]] <- fn
  }

  cat(sprintf("Computing %d RDKit descriptors for %d compounds...\n",
              length(getters), n))

  # Compute for each SMILES (only valid)
  n_getters <- length(getters)
  getter_names <- names(getters)
  result_mat <- matrix(NA_real_, nrow = n, ncol = n_getters)

  for (idx in valid_idx) {
    mol <- tryCatch(Chem$MolFromSmiles(smi[idx]), error = function(e) NULL)
    if (is.null(mol)) next
    vals <- vapply(getter_names, function(nm) {
      tryCatch(as.numeric(getters[[nm]](mol)), error = function(e) NA_real_)
    }, numeric(1L))
    result_mat[idx, ] <- vals
  }

  dt <- data.table::as.data.table(result_mat)
  data.table::setnames(dt, getter_names)

  # Filter to wanted columns if specified
  if (!is.null(wanted_columns)) {
    avail <- intersect(wanted_columns, names(dt))
    if (length(avail)) dt <- dt[, avail, with = FALSE]
  }

  dt
}


#' Compute Mordred descriptors via reticulate
#'
#' Calls Mordred Python library to compute 2D molecular descriptors.
#' Requires `mordred` Python package: `pip install mordred`.
#'
#' @param smiles Character vector of SMILES.
#' @param wanted_columns Character vector of descriptor names to return
#'   (without prefix). If NULL, returns all available.
#' @return A data.table with descriptor columns, or NULL if Mordred unavailable.
#' @keywords internal
#' @noRd
.nts_compute_mordred_descriptors <- function(smiles, wanted_columns = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    warning("reticulate required for Mordred descriptors. Install: install.packages('reticulate')", call. = FALSE)
    return(NULL)
  }

  .nts_configure_python()

  mordred_avail <- tryCatch({
    reticulate::py_module_available("mordred")
  }, error = function(e) FALSE)

  if (!mordred_avail) {
    warning("Mordred Python module not found. Install: pip install mordred", call. = FALSE)
    return(NULL)
  }

  smi <- trimws(as.character(smiles))
  n <- length(smi)

  # Import
  mordred <- reticulate::import("mordred")
  mordred_descriptors <- reticulate::import("mordred.descriptors")
  calc <- mordred$Calculator(mordred_descriptors, ignore_3D = TRUE)

  desc_names <- reticulate::py_to_r(calc$descriptors)
  cat(sprintf("Computing %d Mordred descriptors for %d compounds...\n",
              length(desc_names), n))

  # Compute
  mols <- lapply(smi, function(s) {
    tryCatch(reticulate::r_to_py(reticulate::import("rdkit.Chem")$MolFromSmiles(s)),
             error = function(e) NULL)
  })

  result_list <- lapply(mols, function(mol) {
    if (is.null(mol)) return(rep(NA_real_, length(desc_names)))
    vals <- tryCatch({
      res <- calc(mol)
      vapply(seq_along(desc_names), function(j) {
        tryCatch(as.numeric(res[[j]]), error = function(e) NA_real_)
      }, numeric(1L))
    }, error = function(e) rep(NA_real_, length(desc_names)))
    vals
  })

  dt <- data.table::as.data.table(do.call(rbind, result_list))
  data.table::setnames(dt, desc_names)

  if (!is.null(wanted_columns)) {
    avail <- intersect(wanted_columns, names(dt))
    if (length(avail)) dt <- dt[, avail, with = FALSE]
  }

  dt
}
