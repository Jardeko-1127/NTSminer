# =============================================================================
# 08_spectrum_prediction.R
# Production R Interface for FraG-FIORA Forward MS/MS Spectrum Prediction
# (Combinatorial Sub-formula Graph Cleavage + Deep Environment Propensity)
# =============================================================================

#' Predict complete high-resolution MS/MS spectra from chemical structures
#'
#' Deploys the hybrid FraG-FIORA forward MS/MS prediction engine to compute
#' exact monoisotopic fragment masses and relative intensities directly from
#' candidate SMILES strings across collision energies (NCE 10, 20, 40, or Merged).
#'
#' Combines FraGNNet's exact combinatorial bond cleavage and strict mass conservation
#' (zero unphysical hallucination) with FIORA's deep chemical environment cleavage
#' propensities, gas-phase charge localization, and collision energy modulation.
#' Automatically routes fluorinated pollutants to fine-tuned in-house PFAS model weights.
#'
#' @param smiles Character vector or single character string containing SMILES.
#' @param nce Collision energy: `"Merged"` (default), `"10"`, `"20"`, or `"40"`.
#' @param polarity Ionization polarity: `"negative"` (default, `[M-H]-`) or
#'   `"positive"` (`[M+H]+`).
#' @param model Fragmentation model parameters: `"auto"` (default), `"pfas_inhouse"`,
#'   or `"general"`. If `"auto"`, molecules containing 3 or more fluorine atoms
#'   are automatically processed with fine-tuned PFAS parameters.
#' @param min_rel_intensity Numeric. Cutoff threshold (%) to filter trace noise peaks.
#'   Default `1.0`.
#' @param top_k Optional integer. Limit output to top K most intense peaks. Default `NULL`.
#' @param python_path Optional path to Python executable with RDKit installed.
#'   If `NULL`, auto-detected from environment variables (`NTSMINER_PYTHON`, `RETICULATE_PYTHON`)
#'   or system PATH.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return When `smiles` is of length 1, returns an S3 object of class
#'   `c("nts_predicted_spectrum", "data.frame")` containing:
#'   \itemize{
#'     \item `mz`: Exact monoisotopic fragment m/z (Da).
#'     \item `intensity`: Relative intensity normalized to base peak (0.0 to 100.0).
#'     \item `formula`: Monoisotopic molecular formula of the fragment ion.
#'     \item `annotation`: In-silico cleavage mechanism or neutral loss.
#'     \item `type`: Fragment classification (`"fragment"` or `"precursor"`).
#'   }
#'   When `smiles` is a vector of length > 1, returns a named list of
#'   `nts_predicted_spectrum` objects.
#'
#' @export
predict_ms2_spectrum <- function(smiles,
                                 nce = c("Merged", "10", "20", "40"),
                                 polarity = c("negative", "positive"),
                                 model = c("auto", "pfas_inhouse", "general"),
                                 min_rel_intensity = 1.0,
                                 top_k = NULL,
                                 python_path = NULL,
                                 verbose = FALSE) {
  if (is.null(smiles) || !length(smiles)) {
    return(.nts_empty_predicted_spectrum())
  }

  nce <- if (is.numeric(nce)) as.character(nce[1L]) else match.arg(as.character(nce)[1L], c("Merged", "10", "20", "40"))
  polarity <- match.arg(polarity)
  model <- match.arg(model)

  min_rel_intensity <- suppressWarnings(as.numeric(min_rel_intensity)[1L])
  if (is.na(min_rel_intensity) || min_rel_intensity < 0) min_rel_intensity <- 1.0

  if (!is.null(top_k)) {
    top_k <- suppressWarnings(as.integer(top_k)[1L])
    if (is.na(top_k) || top_k <= 0L) top_k <- NULL
  }

  smiles_vec <- as.character(smiles)
  is_single <- (length(smiles_vec) == 1L)

  # 1. Resolve Python executable with RDKit
  py_exe <- .nts_resolve_forward_python(python_path)
  engine_script <- .nts_resolve_forward_engine_script()

  if (is.null(py_exe) || !file.exists(py_exe) || is.null(engine_script) || !file.exists(engine_script)) {
    warning("Python with RDKit is required for FraG-FIORA MS/MS forward prediction. ",
            "Please configure `python_path` or install RDKit in your Python environment. ",
            "Returning empty predicted spectrum.", call. = FALSE)
    if (is_single) {
      return(.nts_empty_predicted_spectrum(smiles_vec[1L], nce, polarity, model))
    } else {
      res_list <- lapply(smiles_vec, function(s) .nts_empty_predicted_spectrum(s, nce, polarity, model))
      names(res_list) <- names(smiles) %||% smiles_vec
      return(res_list)
    }
  }

  # 2. Build task items and resolve model weights per compound
  items <- vector("list", length(smiles_vec))
  for (i in seq_along(smiles_vec)) {
    smi <- smiles_vec[i]
    chosen_model <- model
    if (chosen_model == "auto") {
      # Count Fluorine atoms in SMILES to determine PFAS-likeness
      n_f <- nchar(gsub("[^F]", "", smi))
      chosen_model <- if (n_f >= 3L) "pfas_inhouse" else "general"
    }

    w_path <- .nts_resolve_forward_weights(chosen_model)
    cid <- if (!is.null(names(smiles))) names(smiles)[i] else paste0("SMI_", i)

    it <- list(
      id = cid,
      smiles = smi,
      nce = nce,
      polarity = polarity,
      weights = w_path %||% "",
      min_rel_intensity = min_rel_intensity,
      model = chosen_model
    )
    if (!is.null(top_k)) it$top_k <- as.integer(top_k)
    items[[i]] <- it
  }

  # 3. Execute batch prediction via temporary files
  tmp_in <- tempfile(pattern = "frag_fiora_in_", fileext = ".json")
  tmp_out <- tempfile(pattern = "frag_fiora_out_", fileext = ".json")
  on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)

  jsonlite::write_json(items, tmp_in, auto_unbox = TRUE)

  cmd_args <- c(
    shQuote(engine_script),
    "--batch_json", shQuote(tmp_in),
    "--output", shQuote(tmp_out)
  )

  run_res <- suppressWarnings(system2(py_exe, args = cmd_args, stdout = TRUE, stderr = TRUE))
  status_code <- attr(run_res, "status") %||% 0L

  if (status_code != 0L || !file.exists(tmp_out)) {
    err_msg <- paste(run_res, collapse = "\n")
    warning("FraG-FIORA MS/MS prediction failed (exit code ", status_code, "):\n", err_msg, call. = FALSE)
    if (is_single) {
      return(.nts_empty_predicted_spectrum(smiles_vec[1L], nce, polarity, model))
    } else {
      res_list <- lapply(smiles_vec, function(s) .nts_empty_predicted_spectrum(s, nce, polarity, model))
      names(res_list) <- names(smiles) %||% smiles_vec
      return(res_list)
    }
  }

  res_json <- tryCatch(jsonlite::fromJSON(tmp_out, simplifyVector = FALSE), error = function(e) NULL)
  results_list <- res_json$results %||% list()

  # 4. Parse results into S3 nts_predicted_spectrum objects
  out_spectra <- vector("list", length(smiles_vec))
  names(out_spectra) <- names(smiles) %||% smiles_vec

  for (i in seq_along(smiles_vec)) {
    smi <- smiles_vec[i]
    m_info <- items[[i]]$model

    match_entry <- NULL
    if (i <= length(results_list)) {
      match_entry <- results_list[[i]]
    }

    if (is.null(match_entry) || identical(match_entry$status, "error") || !length(match_entry$peaks)) {
      out_spectra[[i]] <- .nts_empty_predicted_spectrum(smi, nce, polarity, m_info)
    } else {
      pks <- match_entry$peaks
      mz_v <- as.numeric(vapply(pks, function(p) p$frag_mz %||% 0.0, numeric(1L)))
      int_v <- as.numeric(vapply(pks, function(p) p$intensity %||% 0.0, numeric(1L)))
      form_v <- as.character(vapply(pks, function(p) p$formula %||% "", character(1L)))
      ann_v <- as.character(vapply(pks, function(p) p$annotation %||% "", character(1L)))
      type_v <- as.character(vapply(pks, function(p) p$type %||% "fragment", character(1L)))

      df_pks <- data.frame(
        mz = mz_v,
        intensity = int_v,
        formula = form_v,
        annotation = ann_v,
        type = type_v,
        stringsAsFactors = FALSE
      )
      df_pks <- df_pks[order(df_pks$mz), , drop = FALSE]
      rownames(df_pks) <- NULL

      attr(df_pks, "smiles") <- smi
      attr(df_pks, "nce") <- nce
      attr(df_pks, "polarity") <- polarity
      attr(df_pks, "model") <- m_info
      attr(df_pks, "n_peaks") <- nrow(df_pks)
      class(df_pks) <- c("nts_predicted_spectrum", "data.frame")
      out_spectra[[i]] <- df_pks
    }
  }

  if (is_single) {
    return(out_spectra[[1L]])
  }
  out_spectra
}

#' @export
print.nts_predicted_spectrum <- function(x, ...) {
  smi <- attr(x, "smiles") %||% "<unknown>"
  nce <- attr(x, "nce") %||% "Merged"
  pol <- attr(x, "polarity") %||% "negative"
  model <- attr(x, "model") %||% "auto"
  n_pks <- nrow(x)

  cat("======================================================================\n")
  cat("  FraG-FIORA In-Silico Predicted MS/MS Spectrum (NTSminer)\n")
  cat("======================================================================\n")
  cat(sprintf("SMILES:    %s\n", smi))
  cat(sprintf("Polarity:  %s | NCE: %s | Model: %s\n", pol, nce, model))
  cat(sprintf("Peaks:     %d predicted fragment ions\n", n_pks))
  cat("----------------------------------------------------------------------\n")

  if (n_pks == 0L) {
    cat("  (No peaks predicted)\n")
  } else {
    show_n <- min(n_pks, 10L)
    sub_df <- as.data.frame(x[seq_len(show_n), c("mz", "intensity", "formula", "annotation", "type")], stringsAsFactors = FALSE)
    sub_df$mz <- sprintf("%.4f", as.numeric(sub_df$mz))
    sub_df$intensity <- sprintf("%.2f%%", as.numeric(sub_df$intensity))
    print(sub_df, row.names = FALSE)
    if (n_pks > 10L) {
      cat(sprintf("  ... and %d more peaks (omitted)\n", n_pks - 10L))
    }
  }
  cat("======================================================================\n")
  invisible(x)
}


# =============================================================================
# Candidate Rescoring by Predicted MS/MS
# =============================================================================

#' Rescore and rank candidate structures using forward MS/MS prediction
#'
#' Matches an experimental MS2 fragmentation spectrum against in-silico predicted
#' spectra for candidate chemical structures, scoring candidates by square-root cosine
#' similarity, linear cosine, and explained experimental intensity.
#'
#' @param exp_ms2 A two-column numeric matrix or data frame (`mz`, `intensity`),
#'   or a list containing `mz` and `intensity` vectors.
#' @param candidate_smiles Character vector of candidate SMILES strings.
#' @param candidate_ids Optional character vector of candidate identifiers or compound names.
#'   If `NULL`, uses names of `candidate_smiles` or sequential indices.
#' @param polarity Ionization polarity: `"negative"` (default) or `"positive"`.
#' @param nce Collision energy: `"Merged"` (default), `"10"`, `"20"`, or `"40"`.
#' @param model Model weights to use: `"auto"` (default), `"pfas_inhouse"`, or `"general"`.
#' @param mz_tol Numeric. Mass tolerance for matching experimental and predicted peaks in Da. Default `0.02`.
#' @param min_rel_intensity Minimum relative intensity cutoff (%) for predicted peaks. Default `1.0`.
#' @param weights Named numeric vector of weights for composite score calculation.
#'   Default `c(forward_cosine = 0.6, explained_intensity = 0.4)`.
#' @param python_path Optional Python path.
#' @param verbose Logical. Default `FALSE`.
#'
#' @return A `data.frame` sorted decreasingly by `score_pred`, containing:
#'   \itemize{
#'     \item `candidate_id`: Identifier of candidate.
#'     \item `smiles`: Candidate SMILES.
#'     \item `rank`: Rank among candidates (1 = best hit).
#'     \item `score_pred`: Composite prediction score (0.0 to 1.0).
#'     \item `cosine_sqrt`: Square-root weighted spectral cosine similarity.
#'     \item `cosine_linear`: Standard linear dot-product cosine similarity.
#'     \item `explained_intensity`: Proportion of total experimental intensity explained.
#'     \item `matched_peaks`: Number of matched fragment peaks.
#'     \item `total_pred_peaks`: Total number of predicted peaks.
#'   }
#' @export
rescore_candidates_by_predicted_spectrum <- function(exp_ms2,
                                                     candidate_smiles,
                                                     candidate_ids = NULL,
                                                     polarity = c("negative", "positive"),
                                                     nce = "Merged",
                                                     model = c("auto", "pfas_inhouse", "general"),
                                                     mz_tol = 0.02,
                                                     min_rel_intensity = 1.0,
                                                     weights = c(forward_cosine = 0.6, explained_intensity = 0.4),
                                                     python_path = NULL,
                                                     verbose = FALSE) {
  if (is.null(candidate_smiles) || !length(candidate_smiles)) {
    return(.nts_empty_candidate_ranking())
  }

  exp_clean <- .nts_normalize_exp_peaks(exp_ms2)
  if (!nrow(exp_clean)) {
    warning("Experimental spectrum contains no valid peaks.", call. = FALSE)
    return(.nts_empty_candidate_ranking())
  }

  smi_names <- names(candidate_smiles)
  candidate_smiles <- as.character(candidate_smiles)
  cids <- if (!is.null(candidate_ids)) {
    as.character(candidate_ids)
  } else if (!is.null(smi_names)) {
    smi_names
  } else {
    paste0("Candidate_", seq_along(candidate_smiles))
  }

  polarity <- match.arg(polarity)
  model <- match.arg(model)
  mz_tol <- suppressWarnings(as.numeric(mz_tol)[1L])
  if (is.na(mz_tol) || mz_tol <= 0) mz_tol <- 0.02

  w_cos <- as.numeric(weights["forward_cosine"] %||% 0.6)
  w_exp <- as.numeric(weights["explained_intensity"] %||% 0.4)
  w_sum <- w_cos + w_exp
  if (w_sum > 0) {
    w_cos <- w_cos / w_sum
    w_exp <- w_exp / w_sum
  } else {
    w_cos <- 0.6; w_exp <- 0.4
  }

  if (isTRUE(verbose)) {
    message(sprintf("Predicting forward MS/MS spectra for %d candidates...", length(candidate_smiles)))
  }

  # Predict batch
  pred_list <- predict_ms2_spectrum(
    smiles = candidate_smiles,
    nce = nce,
    polarity = polarity,
    model = model,
    min_rel_intensity = min_rel_intensity,
    python_path = python_path,
    verbose = verbose
  )
  if (!is.list(pred_list) || is.data.frame(pred_list)) {
    pred_list <- list(pred_list)
  }

  rows <- vector("list", length(candidate_smiles))
  for (i in seq_along(candidate_smiles)) {
    smi <- candidate_smiles[i]
    cid <- cids[i]
    pred_spec <- pred_list[[i]]

    if (is.null(pred_spec) || !nrow(pred_spec)) {
      rows[[i]] <- data.frame(
        candidate_id = cid,
        smiles = smi,
        rank = 0L,
        score_pred = 0.0,
        cosine_sqrt = 0.0,
        cosine_linear = 0.0,
        explained_intensity = 0.0,
        matched_peaks = 0L,
        total_pred_peaks = 0L,
        stringsAsFactors = FALSE
      )
      next
    }

    match_res <- .nts_match_spectral_peaks(exp_clean, pred_spec, mz_tol = mz_tol)
    comp_score <- round(w_cos * match_res$cosine_sqrt + w_exp * match_res$explained_intensity, 4)

    rows[[i]] <- data.frame(
      candidate_id = cid,
      smiles = smi,
      rank = 0L,
      score_pred = comp_score,
      cosine_sqrt = round(match_res$cosine_sqrt, 4),
      cosine_linear = round(match_res$cosine_linear, 4),
      explained_intensity = round(match_res$explained_intensity, 4),
      matched_peaks = as.integer(match_res$matched_peaks),
      total_pred_peaks = as.integer(nrow(pred_spec)),
      stringsAsFactors = FALSE
    )
  }

  df_res <- do.call(rbind, rows)
  ord <- order(df_res$score_pred, df_res$cosine_sqrt, df_res$explained_intensity, decreasing = TRUE)
  df_res <- df_res[ord, , drop = FALSE]
  df_res$rank <- seq_len(nrow(df_res))
  rownames(df_res) <- NULL
  df_res
}


# =============================================================================
# Automated In-Silico Annotation Pipeline for NTSProject
# =============================================================================

#' Annotate NTSProject features using forward MS/MS in-silico prediction
#'
#' Screens non-target features against a suspect chemical structure library,
#' performing precursor MS1 mass matching followed by automated forward MS/MS
#' spectrum prediction (FraG-FIORA) and spectral matching against experimental MS2 scans.
#'
#' Passing hits are assigned confidence level **L2b** (Probable structure - in-silico /
#' prediction library match) in accordance with Schymanski (2014) guidelines,
#' appending them to `@db_annotations` with `db_type = "prediction"`.
#'
#' @param object An [NTSProject] or a feature table `data.frame`.
#' @param suspect_table A data frame of candidate suspect compounds. If `NULL`,
#'   automatically resolves the bundled comprehensive 3-library deduplicated PFAS suspect list
#'   (`inst/extdata/tables/NTS_PFAS_suspect_3lib_dedup.csv`).
#' @param smiles_col Character. Column name in `suspect_table` containing SMILES. Default `"SMILES"`.
#' @param id_col Character. Column name for candidate ID. Default `"ID"`.
#' @param name_col Character. Column name for compound name. Default `"Name"`.
#' @param formula_col Character. Column name for molecular formula. Default `"Formula"`.
#' @param mass_col Character. Column name for monoisotopic exact mass. Default `"ExactMass"`.
#' @param target_mz_col Character. Optional column name in `suspect_table` containing verified precursor
#'   or target ion m/z (e.g. `target_mz` taking precedence over neutral exact mass, supporting
#'   in-source fragmentation and adduct overrides). Default `"target_mz"`.
#' @param rt_col Character. Optional column name in `suspect_table` containing reference retention time.
#'   Default `"RT"`.
#' @param rt_tol_s Numeric. Optional retention time tolerance in seconds. If `NULL` (default),
#'   no RT filtering is performed (relaxed mode). If specified (e.g. `150`), candidates with known
#'   RT outside this window are excluded, while candidates with unmeasured (`NA`) RT are retained.
#' @param min_ms2_peaks Integer. Minimum number of experimental MS2 fragment peaks required to
#'   evaluate a feature. Default `1L` (permits matching single prominent fragment ions such as PFBA).
#' @param mz_tol_ms1 Precursor MS1 mass tolerance in ppm. Default `10`.
#' @param mz_tol_ms2 Fragment MS2 mass tolerance in Da. Default `0.02`.
#' @param nce Collision energy for prediction: `"Merged"` (default), `"10"`, `"20"`, or `"40"`.
#' @param polarity Ionization polarity: `"negative"` (default) or `"positive"`.
#' @param model Model weights to use: `"auto"` (default), `"pfas_inhouse"`, or `"general"`.
#' @param min_cosine Minimum square-root cosine similarity threshold to retain hit. Default `0.15`.
#' @param min_explained_int Minimum explained experimental intensity threshold. Default `0.15`.
#' @param top_n_candidates Maximum number of in-silico candidates to retain per feature. Default `5L`.
#' @param skip_annotated Logical. If `TRUE` (default), skips features that already have confirmed
#'   in-house (L1) or online (L2a) annotations in `@db_annotations`.
#' @param python_path Optional Python path.
#' @param verbose Logical. Default `TRUE`.
#'
#' @return An updated [NTSProject] (or `data.frame` of annotations) with in-silico hits
#'   integrated into `@db_annotations` and run logged in `@log`.
#' @export
annotate_by_insilico_prediction <- function(object,
                                           suspect_table = NULL,
                                           smiles_col = "SMILES",
                                           id_col = "ID",
                                           name_col = "Name",
                                           formula_col = "Formula",
                                           mass_col = "ExactMass",
                                           target_mz_col = "target_mz",
                                           rt_col = "RT",
                                           rt_tol_s = NULL,
                                           min_ms2_peaks = 1L,
                                           mz_tol_ms1 = 10,
                                           mz_tol_ms2 = 0.02,
                                           nce = "Merged",
                                           polarity = c("negative", "positive"),
                                           model = c("auto", "pfas_inhouse", "general"),
                                           min_cosine = 0.15,
                                           min_explained_int = 0.15,
                                           top_n_candidates = 5L,
                                           skip_annotated = TRUE,
                                           python_path = NULL,
                                           verbose = TRUE) {
  is_proj <- inherits(object, "NTSProject")
  polarity <- match.arg(polarity)
  model <- match.arg(model)

  # 1. Resolve suspect library
  if (is.null(suspect_table)) {
    suspect_path <- .nts_extdata("tables/NTS_PFAS_suspect_3lib_dedup.csv")
    if (is.na(suspect_path) || !file.exists(suspect_path)) {
      stop("Default suspect list not found: NTS_PFAS_suspect_3lib_dedup.csv. ",
           "Please provide `suspect_table` explicitly.", call. = FALSE)
    }
    suspect_table <- utils::read.csv(suspect_path, stringsAsFactors = FALSE)
  }

  if (!is.data.frame(suspect_table) || !nrow(suspect_table)) {
    stop("`suspect_table` must be a non-empty data frame.", call. = FALSE)
  }

  req_cols <- c(smiles_col, mass_col)
  missing_cols <- setdiff(req_cols, names(suspect_table))
  if (length(missing_cols)) {
    stop("`suspect_table` is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  suspect_table <- suspect_table[!is.na(suspect_table[[mass_col]]) & suspect_table[[mass_col]] > 0, , drop = FALSE]
  suspect_table <- suspect_table[!is.na(suspect_table[[smiles_col]]) & nzchar(trimws(suspect_table[[smiles_col]])), , drop = FALSE]

  # Calculate theoretical precursor m/z for suspects
  proton_mass <- 1.007276466812
  suspect_mzs <- if (polarity == "negative") {
    suspect_table[[mass_col]] - proton_mass
  } else {
    suspect_table[[mass_col]] + proton_mass
  }

  # Canonical target_mz override when available (target_mz > ExactMass +/- proton)
  resolved_tmz_col <- if (!is.null(target_mz_col) && target_mz_col %in% names(suspect_table)) {
    target_mz_col
  } else {
    cand_tmz <- c("target_mz", "targetMz", "precursor_mz", if (polarity == "negative") "mz.neg" else "mz.pos")
    intersect(cand_tmz, names(suspect_table))[1L]
  }
  if (!is.null(resolved_tmz_col) && !is.na(resolved_tmz_col)) {
    raw_tmz <- suppressWarnings(as.numeric(suspect_table[[resolved_tmz_col]]))
    ok_tmz <- is.finite(raw_tmz) & raw_tmz > 0
    suspect_mzs[ok_tmz] <- raw_tmz[ok_tmz]
  }

  # Resolve candidate retention times when available
  suspect_rts <- rep(NA_real_, nrow(suspect_table))
  resolved_rt_col <- if (!is.null(rt_col) && rt_col %in% names(suspect_table)) {
    rt_col
  } else {
    cand_rts <- c("RT", "rt", "lib_rt", "best_rt")
    intersect(cand_rts, names(suspect_table))[1L]
  }
  if (!is.null(resolved_rt_col) && !is.na(resolved_rt_col)) {
    raw_rts <- suppressWarnings(as.numeric(suspect_table[[resolved_rt_col]]))
    suspect_rts <- raw_rts
  }

  # 2. Extract features from object
  ft <- if (is_proj) get_feature_table(object) else object
  if (!is.data.frame(ft) || !nrow(ft)) {
    if (is_proj) return(object) else return(data.frame())
  }

  feat_id_col <- if ("name" %in% names(ft)) "name" else if ("feature_id" %in% names(ft)) "feature_id" else names(ft)[1L]
  feat_mz_col <- if ("mz" %in% names(ft)) "mz" else if ("mzmed" %in% names(ft)) "mzmed" else NULL
  feat_rt_col <- if ("rt" %in% names(ft)) "rt" else if ("rtmed" %in% names(ft)) "rtmed" else NULL

  if (is.null(feat_mz_col)) {
    stop("Feature table must contain an `mz` column.", call. = FALSE)
  }

  # Exclude already confirmed features if requested
  already_annotated <- character(0)
  if (is_proj && isTRUE(skip_annotated) && .hasSlot(object, "db_annotations")) {
    db_ann <- methods::slot(object, "db_annotations")
    if (is.data.frame(db_ann) && nrow(db_ann) && "name" %in% names(db_ann)) {
      high_conf <- db_ann$db_type %in% c("inhouse", "online")
      already_annotated <- unique(as.character(db_ann$name[high_conf]))
    }
  }

  # 3. Retrieve experimental MS2 data
  ms2_source <- NULL
  if (is_proj) {
    if (.hasSlot(object, "clean_ms2") && length(object@clean_ms2)) {
      ms2_source <- object@clean_ms2
    } else if (.hasSlot(object, "ms2data") && length(object@ms2data)) {
      ms2_source <- object@ms2data
    }
  }

  n_features_checked <- 0L
  n_features_annotated <- 0L
  all_ann_rows <- list()

  for (i in seq_len(nrow(ft))) {
    fid <- as.character(ft[[feat_id_col]][i])
    if (fid %in% already_annotated) next

    fmz <- as.numeric(ft[[feat_mz_col]][i])
    frt <- if (!is.null(feat_rt_col)) as.numeric(ft[[feat_rt_col]][i]) else NA_real_

    if (is.na(fmz) || fmz <= 0) next
    n_features_checked <- n_features_checked + 1L

    # MS1 precursor match (ppm)
    ppm_diff <- abs(suspect_mzs - fmz) / fmz * 1e6
    hit_idx <- which(ppm_diff <= mz_tol_ms1)

    # Optional relaxed RT match (if rt_tol_s provided and feature has finite RT)
    if (length(hit_idx) && !is.null(rt_tol_s) && is.finite(rt_tol_s) && rt_tol_s > 0 && is.finite(frt)) {
      cand_rts <- suspect_rts[hit_idx]
      # Keep candidate if cand_rt is NA or within relaxed rt_tol_s window
      rt_pass <- is.na(cand_rts) | (abs(cand_rts - frt) <= rt_tol_s)
      hit_idx <- hit_idx[rt_pass]
    }

    if (!length(hit_idx)) next

    # Check MS2 spectrum
    exp_peaks <- NULL
    if (!is.null(ms2_source) && !is.null(ms2_source[[fid]])) {
      sp_entry <- ms2_source[[fid]]
      exp_peaks <- .nts_normalize_exp_peaks(sp_entry)
    }

    if (is.null(exp_peaks) || nrow(exp_peaks) < as.integer(min_ms2_peaks)) {
      # No MS2 spectrum available to substantiate in-silico prediction
      next
    }

    cand_sub <- suspect_table[hit_idx, , drop = FALSE]
    cand_smiles <- cand_sub[[smiles_col]]
    cand_ids <- if (id_col %in% names(cand_sub)) cand_sub[[id_col]] else paste0("SUSPECT_", hit_idx)
    cand_names <- if (name_col %in% names(cand_sub)) cand_sub[[name_col]] else cand_ids
    cand_forms <- if (formula_col %in% names(cand_sub)) cand_sub[[formula_col]] else rep("", length(hit_idx))
    cand_rts_sub <- suspect_rts[hit_idx]

    # Rescore candidates
    rescored <- rescore_candidates_by_predicted_spectrum(
      exp_ms2 = exp_peaks,
      candidate_smiles = cand_smiles,
      candidate_ids = cand_ids,
      polarity = polarity,
      nce = nce,
      model = model,
      mz_tol = mz_tol_ms2,
      python_path = python_path,
      verbose = FALSE
    )

    # Filter passing hits
    passing <- rescored[rescored$cosine_sqrt >= min_cosine | rescored$explained_intensity >= min_explained_int, , drop = FALSE]
    if (!nrow(passing)) next

    passing <- passing[seq_len(min(nrow(passing), as.integer(top_n_candidates))), , drop = FALSE]
    n_features_annotated <- n_features_annotated + 1L

    for (k in seq_len(nrow(passing))) {
      p_cid <- passing$candidate_id[k]
      p_idx <- match(p_cid, cand_ids)
      p_name <- if (!is.na(p_idx)) cand_names[p_idx] else p_cid
      p_form <- if (!is.na(p_idx)) cand_forms[p_idx] else ""
      p_smi <- passing$smiles[k]
      p_score <- passing$score_pred[k]
      p_cos <- passing$cosine_sqrt[k]
      p_exp <- passing$explained_intensity[k]
      p_mat_pks <- passing$matched_peaks[k]
      p_theo_mz <- if (!is.na(p_idx)) suspect_mzs[hit_idx[p_idx]] else fmz
      p_ppm <- abs(p_theo_mz - fmz) / fmz * 1e6
      p_lib_rt <- if (!is.na(p_idx) && is.finite(cand_rts_sub[p_idx])) cand_rts_sub[p_idx] else NA_real_
      p_rt_err <- if (is.finite(p_lib_rt) && is.finite(frt)) round(frt - p_lib_rt, 2) else NA_real_
      p_score_rt <- if (is.finite(p_rt_err)) {
        denom <- if (!is.null(rt_tol_s) && is.finite(rt_tol_s) && rt_tol_s > 0) rt_tol_s else 150.0
        round(pmax(0.0, 1.0 - (abs(p_rt_err) / denom) * 0.2), 4)
      } else {
        NA_real_
      }

      all_ann_rows[[length(all_ann_rows) + 1L]] <- data.frame(
        name = fid,
        ref_id = as.character(p_cid),
        db_type = "prediction",
        Compound.name = as.character(p_name),
        SMILES = as.character(p_smi),
        Formula = as.character(p_form),
        score_total = as.numeric(p_score),
        score_mz = round(1.0 - (p_ppm / mz_tol_ms1) * 0.2, 4),
        score_rt = p_score_rt,
        score_ms2 = as.numeric(p_cos),
        score_cosine_ms2 = as.numeric(p_cos),
        score_nl = 0.0,
        score_md = 0.0,
        score_diagnostic = 0.0,
        score_ndp = 0.0,
        ndp_raw = 0.0,
        lib_CE = as.character(nce),
        exp_CE = as.character(nce),
        best_ms2_index = 1L,
        n_ms2_subset = 1L,
        matched_nl = "",
        matched_md = "",
        matched_diagnostic = "",
        exp_mz = fmz,
        lib_mz = p_theo_mz,
        mz_error_ppm = round(p_ppm, 3),
        exp_rt = frt,
        lib_rt = p_lib_rt,
        rt_error_s = p_rt_err,
        n_matched_peaks = as.integer(p_mat_pks),
        explained_intensity = as.numeric(p_exp),
        stringsAsFactors = FALSE
      )
    }
  }

  df_new_ann <- do.call(rbind, all_ann_rows)
  if (is.null(df_new_ann)) df_new_ann <- data.frame()

  if (isTRUE(verbose)) {
    message(sprintf("In-silico forward MS/MS screening complete: %d features evaluated, %d features annotated (%d total candidate hits).",
                    n_features_checked, n_features_annotated, nrow(df_new_ann)))
  }

  if (!is_proj) {
    return(df_new_ann)
  }

  # Append to existing @db_annotations in NTSProject
  cur_ann <- if (.hasSlot(object, "db_annotations")) methods::slot(object, "db_annotations") else data.frame()
  if (nrow(df_new_ann) > 0L) {
    if (is.data.frame(cur_ann) && nrow(cur_ann) > 0L) {
      common_cols <- intersect(names(cur_ann), names(df_new_ann))
      methods::slot(object, "db_annotations") <- rbind(cur_ann[, common_cols, drop = FALSE], df_new_ann[, common_cols, drop = FALSE])
    } else {
      methods::slot(object, "db_annotations") <- df_new_ann
    }
  }

  # Append log entry
  if (.hasSlot(object, "log")) {
    object@log <- .nts_log_append(
      object@log,
      "annotate_by_insilico_prediction",
      list(
        n_features_checked = n_features_checked,
        n_features_annotated = n_features_annotated,
        n_hits_retained = nrow(df_new_ann),
        nce = nce,
        polarity = polarity
      )
    )
  }

  object
}


# =============================================================================
# In-Silico Predicted MS/MS Mirror Plotting
# =============================================================================

#' High-resolution mirror plot comparing experimental MS2 to in-silico predicted MS2
#'
#' Visualizes experimental fragmentation spectrum (top pane) vs FraG-FIORA predicted
#' spectrum (bottom pane), highlighting matched fragments and displaying similarity metrics.
#'
#' @param exp_spectrum Two-column matrix or data frame of experimental peaks (`mz`, `intensity`).
#' @param smiles SMILES string of the reference chemical structure.
#' @param nce Collision energy: `"Merged"` (default), `"10"`, `"20"`, or `"40"`.
#' @param polarity Ionization polarity: `"negative"` (default) or `"positive"`.
#' @param model Model weights: `"auto"`, `"pfas_inhouse"`, or `"general"`.
#' @param title Optional plot title.
#' @param mz_tol Peak matching tolerance in Da. Default `0.02`.
#' @param min_rel_intensity Minimum relative intensity cutoff (%) for peaks. Default `1.0`.
#' @param label_matched_mz Logical. Whether to print m/z labels on matched peaks. Default `TRUE`.
#' @param python_path Optional custom Python path.
#' @param ... Additional arguments.
#'
#' @return A `ggplot2` plot object.
#' @export
plot_predicted_mirror <- function(exp_spectrum,
                                  smiles,
                                  nce = "Merged",
                                  polarity = c("negative", "positive"),
                                  model = c("auto", "pfas_inhouse", "general"),
                                  title = NULL,
                                  mz_tol = 0.02,
                                  min_rel_intensity = 1.0,
                                  label_matched_mz = TRUE,
                                  python_path = NULL,
                                  ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_predicted_mirror()", call. = FALSE)
  }

  polarity <- match.arg(polarity)
  model <- match.arg(model)
  mz_tol <- suppressWarnings(as.numeric(mz_tol)[1L])
  if (is.na(mz_tol) || mz_tol <= 0) mz_tol <- 0.02

  # 1. Clean experimental spectrum
  exp_df <- .nts_normalize_exp_peaks(exp_spectrum, min_rel_int = min_rel_intensity)
  if (!nrow(exp_df)) {
    stop("Experimental spectrum contains no peaks above threshold.", call. = FALSE)
  }

  # 2. Predict forward spectrum
  pred_spec <- predict_ms2_spectrum(
    smiles = smiles,
    nce = nce,
    polarity = polarity,
    model = model,
    min_rel_intensity = min_rel_intensity,
    python_path = python_path
  )

  if (!nrow(pred_spec)) {
    stop("FraG-FIORA predicted no peaks for given SMILES.", call. = FALSE)
  }

  # 3. Match peaks
  match_res <- .nts_match_spectral_peaks(exp_df, pred_spec, mz_tol = mz_tol)
  cos_sqrt <- match_res$cosine_sqrt
  exp_int <- match_res$explained_intensity

  # 4. Prepare plotting data frames
  df_exp <- data.frame(
    mz = exp_df$mz,
    intensity = exp_df$intensity,
    source = "Experimental",
    matched = seq_len(nrow(exp_df)) %in% match_res$matched_exp_indices,
    stringsAsFactors = FALSE
  )

  df_pred <- data.frame(
    mz = pred_spec$mz,
    intensity = -pred_spec$intensity, # Invert y-axis for mirror display
    source = "FraG-FIORA (Predicted)",
    matched = seq_len(nrow(pred_spec)) %in% match_res$matched_pred_indices,
    stringsAsFactors = FALSE
  )

  df_all <- rbind(df_exp, df_pred)

  # Title & subtitle
  main_title <- title %||% "FraG-FIORA In-Silico MS/MS Mirror Plot"
  sub_title <- sprintf("Cosine (sqrt): %.3f | Explained Int: %.1f%% | Matched Peaks: %d / %d | NCE: %s",
                       cos_sqrt, 100 * exp_int, match_res$matched_peaks, nrow(exp_df), nce)

  p <- ggplot2::ggplot(df_all) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey40") +
    # Experimental peaks (top, red)
    ggplot2::geom_segment(
      data = df_all[df_all$source == "Experimental", , drop = FALSE],
      ggplot2::aes(x = .data[["mz"]], xend = .data[["mz"]], y = 0, yend = .data[["intensity"]]),
      colour = "#D62728",
      linewidth = 0.55
    ) +
    # Predicted peaks (bottom, blue)
    ggplot2::geom_segment(
      data = df_all[df_all$source != "Experimental", , drop = FALSE],
      ggplot2::aes(x = .data[["mz"]], xend = .data[["mz"]], y = 0, yend = .data[["intensity"]]),
      colour = "#1F77B4",
      linewidth = 0.55
    ) +
    ggplot2::scale_y_continuous(
      limits = c(-115, 115),
      breaks = c(-100, -50, 0, 50, 100),
      labels = function(x) abs(x)
    ) +
    ggplot2::labs(
      title = main_title,
      subtitle = sub_title,
      x = "m/z",
      y = "Relative Intensity (%)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9.5, colour = "grey30"),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (isTRUE(label_matched_mz)) {
    # Add text labels on matched experimental peaks
    exp_matched <- df_exp[df_exp$matched, , drop = FALSE]
    if (nrow(exp_matched) > 0L) {
      p <- p + ggplot2::geom_text(
        data = exp_matched,
        ggplot2::aes(x = .data[["mz"]], y = .data[["intensity"]] + 4, label = sprintf("%.2f", .data[["mz"]])),
        size = 2.6, colour = "#B71C1C", angle = 90, hjust = 0, vjust = 0.5
      )
    }
  }

  p
}


# =============================================================================
# Internal Helpers
# =============================================================================

#' Construct empty predicted spectrum structure
#' @keywords internal
#' @noRd
.nts_empty_predicted_spectrum <- function(smiles = "", nce = "Merged", polarity = "negative", model = "auto") {
  df <- data.frame(
    mz = numeric(0),
    intensity = numeric(0),
    formula = character(0),
    annotation = character(0),
    type = character(0),
    stringsAsFactors = FALSE
  )
  attr(df, "smiles") <- as.character(smiles %||% "")
  attr(df, "nce") <- as.character(nce %||% "Merged")
  attr(df, "polarity") <- as.character(polarity %||% "negative")
  attr(df, "model") <- as.character(model %||% "auto")
  attr(df, "n_peaks") <- 0L
  class(df) <- c("nts_predicted_spectrum", "data.frame")
  df
}

#' Construct empty candidate ranking data frame
#' @keywords internal
#' @noRd
.nts_empty_candidate_ranking <- function() {
  data.frame(
    candidate_id = character(0),
    smiles = character(0),
    rank = integer(0),
    score_pred = numeric(0),
    cosine_sqrt = numeric(0),
    cosine_linear = numeric(0),
    explained_intensity = numeric(0),
    matched_peaks = integer(0),
    total_pred_peaks = integer(0),
    stringsAsFactors = FALSE
  )
}

#' Normalize experimental peaks to base-peak 100
#' @keywords internal
#' @noRd
.nts_normalize_exp_peaks <- function(exp_ms2, min_rel_int = 0.0) {
  if (is.null(exp_ms2)) return(data.frame(mz = numeric(0), intensity = numeric(0)))

  mz_vec <- numeric(0)
  int_vec <- numeric(0)

  if (is.data.frame(exp_ms2) || is.matrix(exp_ms2)) {
    df <- as.data.frame(exp_ms2)
    if (ncol(df) >= 2L && nrow(df) > 0L) {
      mz_vec <- as.numeric(if ("mz" %in% names(df)) df$mz else df[, 1L])
      int_vec <- as.numeric(if ("intensity" %in% names(df)) df$intensity else df[, 2L])
    }
  } else if (is.list(exp_ms2)) {
    if (!is.null(exp_ms2$spectrum) && (is.data.frame(exp_ms2$spectrum) || is.matrix(exp_ms2$spectrum))) {
      return(.nts_normalize_exp_peaks(exp_ms2$spectrum, min_rel_int))
    } else if (!is.null(exp_ms2$mz) && !is.null(exp_ms2$intensity)) {
      mz_vec <- as.numeric(exp_ms2$mz)
      int_vec <- as.numeric(exp_ms2$intensity)
    } else if (length(exp_ms2) > 0L && is.list(exp_ms2[[1L]]) && !is.null(exp_ms2[[1L]]$mz)) {
      mz_vec <- as.numeric(exp_ms2[[1L]]$mz)
      int_vec <- as.numeric(exp_ms2[[1L]]$intensity)
    }
  }

  ok <- !is.na(mz_vec) & !is.na(int_vec) & is.finite(mz_vec) & is.finite(int_vec) & mz_vec > 0 & int_vec > 0
  if (!any(ok)) return(data.frame(mz = numeric(0), intensity = numeric(0)))

  mz_vec <- mz_vec[ok]
  int_vec <- int_vec[ok]

  max_int <- max(int_vec)
  rel_int <- if (max_int > 0) (int_vec / max_int) * 100.0 else rep(0.0, length(int_vec))

  keep <- rel_int >= min_rel_int
  if (!any(keep)) return(data.frame(mz = numeric(0), intensity = numeric(0)))

  res <- data.frame(
    mz = mz_vec[keep],
    intensity = rel_int[keep],
    stringsAsFactors = FALSE
  )
  res <- res[order(res$mz), , drop = FALSE]
  rownames(res) <- NULL
  res
}

#' Spectral matching between experimental and predicted peak tables
#' @keywords internal
#' @noRd
.nts_match_spectral_peaks <- function(exp_peaks, pred_peaks, mz_tol = 0.02) {
  if (!nrow(exp_peaks) || !nrow(pred_peaks)) {
    return(list(
      cosine_sqrt = 0.0,
      cosine_linear = 0.0,
      explained_intensity = 0.0,
      matched_peaks = 0L,
      matched_exp_indices = integer(0),
      matched_pred_indices = integer(0)
    ))
  }

  e_mz <- exp_peaks$mz
  e_int <- exp_peaks$intensity
  p_mz <- pred_peaks$mz
  p_int <- pred_peaks$intensity

  # Pair matching: find closest predicted peak within mz_tol for each experimental peak
  pairs <- list()
  for (i in seq_along(e_mz)) {
    diffs <- abs(p_mz - e_mz[i])
    cand_j <- which(diffs <= mz_tol)
    if (length(cand_j)) {
      best_j <- cand_j[which.min(diffs[cand_j])]
      pairs[[length(pairs) + 1L]] <- list(
        i = i,
        j = best_j,
        diff = diffs[best_j]
      )
    }
  }

  if (!length(pairs)) {
    return(list(
      cosine_sqrt = 0.0,
      cosine_linear = 0.0,
      explained_intensity = 0.0,
      matched_peaks = 0L,
      matched_exp_indices = integer(0),
      matched_pred_indices = integer(0)
    ))
  }

  # Greedy one-to-one assignment ordered by smallest mz diff
  pairs_df <- do.call(rbind, lapply(pairs, as.data.frame))
  pairs_df <- pairs_df[order(pairs_df$diff), , drop = FALSE]

  used_i <- logical(length(e_mz))
  used_j <- logical(length(p_mz))
  final_i <- integer(0)
  final_j <- integer(0)

  for (k in seq_len(nrow(pairs_df))) {
    cur_i <- pairs_df$i[k]
    cur_j <- pairs_df$j[k]
    if (!used_i[cur_i] && !used_j[cur_j]) {
      used_i[cur_i] <- TRUE
      used_j[cur_j] <- TRUE
      final_i <- c(final_i, cur_i)
      final_j <- c(final_j, cur_j)
    }
  }

  matched_e_int <- e_int[final_i]
  matched_p_int <- p_int[final_j]

  # Sqrt cosine
  num_sqrt <- sum(sqrt(matched_e_int) * sqrt(matched_p_int))
  den_sqrt <- sqrt(sum(e_int)) * sqrt(sum(p_int))
  cos_sqrt <- if (den_sqrt > 0) num_sqrt / den_sqrt else 0.0

  # Linear cosine
  num_lin <- sum(matched_e_int * matched_p_int)
  den_lin <- sqrt(sum(e_int^2)) * sqrt(sum(p_int^2))
  cos_lin <- if (den_lin > 0) num_lin / den_lin else 0.0

  # Explained intensity
  exp_int_sum <- sum(e_int)
  explained_int <- if (exp_int_sum > 0) sum(matched_e_int) / exp_int_sum else 0.0

  list(
    cosine_sqrt = min(max(as.numeric(cos_sqrt), 0.0), 1.0),
    cosine_linear = min(max(as.numeric(cos_lin), 0.0), 1.0),
    explained_intensity = min(max(as.numeric(explained_int), 0.0), 1.0),
    matched_peaks = length(final_i),
    matched_exp_indices = final_i,
    matched_pred_indices = final_j
  )
}

#' Resolve Python executable with RDKit
#' @keywords internal
#' @noRd
.nts_resolve_forward_python <- function(python_path = NULL) {
  check_cmd <- shQuote("from rdkit import Chem")
  if (!is.null(python_path)) {
    if (nzchar(python_path) && file.exists(python_path)) {
      ok <- tryCatch({
        system2(python_path, c("-c", check_cmd), stdout = FALSE, stderr = FALSE) == 0L
      }, error = function(e) FALSE)
      if (isTRUE(ok)) return(normalizePath(python_path, mustWork = FALSE))
    }
    return(NULL)
  }

  env_py <- Sys.getenv("NTSMINER_PYTHON", "")
  if (nzchar(env_py) && file.exists(env_py)) {
    return(normalizePath(env_py, mustWork = FALSE))
  }

  ret_py <- Sys.getenv("RETICULATE_PYTHON", "")
  if (nzchar(ret_py) && file.exists(ret_py)) {
    return(normalizePath(ret_py, mustWork = FALSE))
  }

  # Well-known conda environments
  common_cands <- c(
    "D:/miniconda/envs/admet_env/python.exe",
    "C:/miniconda3/envs/admet_env/python.exe",
    "C:/Anaconda3/envs/admet_env/python.exe"
  )
  check_cmd <- shQuote("from rdkit import Chem")
  for (cand in common_cands) {
    if (file.exists(cand)) {
      ok <- tryCatch({
        system2(cand, c("-c", check_cmd), stdout = FALSE, stderr = FALSE) == 0L
      }, error = function(e) FALSE)
      if (isTRUE(ok)) return(normalizePath(cand, mustWork = FALSE))
    }
  }

  # Reticulate config
  if (requireNamespace("reticulate", quietly = TRUE)) {
    cfg <- tryCatch(reticulate::py_discover_config(), error = function(e) NULL)
    if (!is.null(cfg) && !is.null(cfg$python) && file.exists(cfg$python)) {
      ok <- tryCatch({
        system2(cfg$python, c("-c", check_cmd), stdout = FALSE, stderr = FALSE) == 0L
      }, error = function(e) FALSE)
      if (isTRUE(ok)) return(normalizePath(cfg$python, mustWork = FALSE))
    }
  }

  # System PATH candidates
  for (cmd in c("python", "python3")) {
    sys_p <- Sys.which(cmd)
    if (nzchar(sys_p) && file.exists(sys_p)) {
      ok <- tryCatch({
        system2(sys_p, c("-c", check_cmd), stdout = FALSE, stderr = FALSE) == 0L
      }, error = function(e) FALSE)
      if (isTRUE(ok)) return(normalizePath(sys_p, mustWork = FALSE))
    }
  }

  NULL
}

#' Resolve FraG-FIORA engine script path
#' @keywords internal
#' @noRd
.nts_resolve_forward_engine_script <- function() {
  cand <- c(
    system.file("python", "frag_fiora_engine.py", package = "NTSminer"),
    file.path("inst", "python", "frag_fiora_engine.py"),
    file.path("python", "frag_fiora_engine.py")
  )
  hit <- cand[nzchar(cand) & file.exists(cand)]
  if (length(hit)) {
    return(normalizePath(hit[1L], mustWork = FALSE))
  }
  normalizePath(file.path("inst", "python", "frag_fiora_engine.py"), mustWork = FALSE)
}

#' Resolve model weights path
#' @keywords internal
#' @noRd
.nts_resolve_forward_weights <- function(model = c("auto", "pfas_inhouse", "general")) {
  model <- match.arg(model)
  filename <- if (model == "pfas_inhouse") "pfas_finetuned_weights.json" else "general_weights.json"

  cand <- c(
    .nts_extdata(file.path("models", "forward_ms2", filename)),
    file.path("inst", "extdata", "models", "forward_ms2", filename),
    file.path("extdata", "models", "forward_ms2", filename)
  )
  hit <- cand[!is.na(cand) & nzchar(cand) & file.exists(cand)]
  if (length(hit)) {
    return(normalizePath(hit[1L], mustWork = FALSE))
  }
  NULL
}
