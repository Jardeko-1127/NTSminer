#' Predict ADMET / toxicity endpoints (ADMETlab 3.0 and/or ADMET-AI)
#'
#' Batch prediction for environmental NTS hits (PFAS, PPCPs, pesticides, etc.)
#' with valid SMILES. Uses **ADMETlab 3.0** REST API and/or local **ADMET-AI**
#' (Chemprop-RDKit on TDC datasets). When `backend = "both"`, results are
#' merged on `SMILES` with disambiguated duplicate column names.
#'
#' @param smiles_list Character vector of SMILES strings.
#' @param backend One of `"admetlab3"`, `"admet_ai"`, or `"both"`.
#' @param timeout_sec HTTP timeout for ADMETlab requests (seconds).
#' @param admetlab_batch_size SMILES per ADMETlab request (default 5; server limit ~5 req/s).
#' @param conda_env Conda environment for ADMET-AI reticulate path (default `admet_env`).
#' @param prefer_admet_ai_cli If `TRUE`, try `admet_predict` CLI before reticulate (default `FALSE`, as in `tests/toxPi_202506065.R`).
#' @param admet_ai_python Optional explicit `python.exe` for ADMET-AI; else reads **`NTSminer_ADMET_AI_PYTHON`**.
#' @param progress Message batch progress.
#' @param ... Passed to [predict_admetlab3] or [predict_admet_ai].
#'
#' @return A `data.table` with `SMILES` and prediction columns from the selected
#'   backend(s). ADMETlab failures include `admetlab_status = "failed"`.
#'
#' @seealso [predict_episuite], [predict_properties], [predict_admetlab3],
#'   [predict_admet_ai]
#'
#' @export
predict_admet_risk <- function(smiles_list,
                              backend = c("both", "admetlab3", "admet_ai"),
                              timeout_sec = 120L,
                              admetlab_batch_size = 5L,
                              conda_env = Sys.getenv("NTSminer_ADMET_AI_ENV", "admet_env"),
                              prefer_admet_ai_cli = FALSE,
                              admet_ai_python = Sys.getenv("NTSminer_ADMET_AI_PYTHON", ""),
                              progress = TRUE,
                              ...) {
  .nts_check(NULL, "data.table_ns")
  backend <- match.arg(backend)
  sm <- unique(trimws(as.character(smiles_list)))
  sm <- sm[!is.na(sm) & nzchar(sm)]
  if (!length(sm)) {
    return(data.table::data.table(SMILES = character()))
  }

  lab <- NULL
  ai <- NULL

  if (backend %in% c("admetlab3", "both")) {
    lab <- tryCatch(
      predict_admetlab3(
        sm,
        batch_size = admetlab_batch_size,
        timeout_sec = timeout_sec,
        progress = progress,
        ...
      ),
      error = function(e) {
        warning("predict_admetlab3 failed: ", conditionMessage(e), call. = FALSE)
        data.frame(SMILES = sm, admetlab_status = "error", stringsAsFactors = FALSE)
      }
    )
  }

  if (backend %in% c("admet_ai", "both")) {
    py_path <- if (nzchar(admet_ai_python)) admet_ai_python else NULL
    ai <- tryCatch(
      predict_admet_ai(
        sm,
        conda_env = conda_env,
        prefer_cli = prefer_admet_ai_cli,
        python_path = py_path
      ),
      error = function(e) {
        warning("predict_admet_ai failed: ", conditionMessage(e), call. = FALSE)
        data.frame(
          SMILES = sm,
          admet_ai_status = "error",
          admet_ai_error = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }

  out <- if (identical(backend, "admetlab3")) {
    lab
  } else if (identical(backend, "admet_ai")) {
    ai
  } else {
    .nts_merge_admet_backends(lab, ai)
  }

  data.table::as.data.table(out)
}

# ============================================================================
# Internal-Standard Normalisation
# ============================================================================

#' Normalise sample abundances by internal standard (IS) responses
#'
#' Divides each sample column by IS response(s), then multiplies by
#' `scale_factor`. Use [match_internal_standards()] first when IS locations come
#' from a CSV (name, m/z, rt).
#'
#' @param object An [NTSProject] or a `data.frame` feature table.
#' @param is_names Character vector of IS **feature** names (`feature_table$name`).
#' @param abundance_cols Optional sample columns; auto-detected if `NULL`.
#' @param scale_factor Multiplier after division (default `1e6`).
#' @param method `"geom_mean"` (all IS combined per sample), `"nearest_rt"`
#'   (per feature: IS with closest RT; recommended for long gradients),
#'   `"loess"` (interpolate IS response along RT), or `"rt_segment"` (piecewise
#'   by RT windows; use `rt_breaks` or midpoints between IS).
#' @param is_match Optional table from [match_internal_standards()] (labels
#'   `is_anchor` with `is_name`).
#' @param rt_breaks RT boundaries in **seconds** for `method = "rt_segment"`.
#'
#' @return Updated object with `is_normalised`, `is_normalise_method`, `is_anchor`.
#' @export
#'
#' @seealso [read_internal_standards], [match_internal_standards], [run_semi_quantification]
normalize_by_IS <- function(object,
                            is_names,
                            abundance_cols = NULL,
                            scale_factor = 1e6,
                            method = c("geom_mean", "nearest_rt", "loess", "rt_segment"),
                            is_match = NULL,
                            rt_breaks = NULL) {
  method <- match.arg(method)
  if (inherits(object, "NTSProject")) {
    ft <- get_feature_table(object)
    if (!nrow(ft)) {
      stop("feature_table is empty", call. = FALSE)
    }
    out <- .normalize_by_IS_df(
      ft,
      is_names = is_names,
      abundance_cols = abundance_cols,
      scale_factor = scale_factor,
      method = method,
      is_match = is_match,
      rt_breaks = rt_breaks
    )
    object <- set_feature_table(object, out)
    methods::slot(object, "log") <- .nts_log_append(
      methods::slot(object, "log"),
      "normalize_by_IS",
      list(
        is_names = is_names,
        scale_factor = scale_factor,
        method = method
      )
    )
    return(object)
  }
  if (!is.data.frame(object) && !data.table::is.data.table(object)) {
    stop("`object` must be a feature table or NTSProject", call. = FALSE)
  }
  .normalize_by_IS_df(
    as.data.frame(object, stringsAsFactors = FALSE),
    is_names = is_names,
    abundance_cols = abundance_cols,
    scale_factor = scale_factor,
    method = method,
    is_match = is_match,
    rt_breaks = rt_breaks
  )
}


# ============================================================================
# Semi-Quantification (RIE / Concentration)
# ============================================================================

#' Sample-intensity columns for semi-quantification (exclude metadata)
#' @noRd
.smq_abundance_cols <- function(ft) {
  reserved <- c(
    .nts_feature_shape_cols,
    "name", "feature_id", "row_ID",
    "is_isf", "parent_mz", "homologous_series_id",
    "is_isotope", "parent_name",
    "isotope_shift_da", "isotope_kind", "isotope_label", "isotope_c13_n",
    "inferred_formula", "annotation_level",
    "nominal_mass", "kmd",
    "rie", "estimated_conc",
    "is_normalised", "is_normalise_method", "is_anchor",
    "SMILES", "LogP_pred", "Toxicity_Class",
    "ALogP_proxy", "TPSA", "MolecularWeight",
    "HBondDonorCount", "HBondAcceptorCount"
  )
  nm <- names(ft)
  cand <- setdiff(nm, reserved)
  cand[vapply(ft[, cand, drop = FALSE], is.numeric, logical(1L))]
}


#' Row-wise intensity for concentration (mean abundance or `intensity` column)
#' @noRd
.smq_row_intensity <- function(ft) {
  if ("intensity" %in% names(ft)) {
    return(as.numeric(ft$intensity))
  }
  ab <- .smq_abundance_cols(ft)
  if (length(ab)) {
    return(rowMeans(ft[, ab, drop = FALSE], na.rm = TRUE))
  }
  rep(NA_real_, nrow(ft))
}


#' First usable numeric column from a candidate list, restricted to a mask
#'
#' Used to hand the measured retention time (and any other run-level metadata)
#' to [predict_logIE] when the feature table carries it. Only the rows selected
#' by `mask` are returned, because `predict_logIE()` builds one matrix row per
#' *valid* SMILES.
#'
#' @param ft A feature table.
#' @param candidates Column names to try, in priority order.
#' @param mask Logical vector over `nrow(ft)`, the rows to keep.
#' @return Numeric vector (length `sum(mask)`) or `NULL` when no candidate column
#'   exists or holds no finite value.
#' @keywords internal
#' @noRd
.smq_first_numeric <- function(ft, candidates, mask) {
  hit <- candidates[candidates %in% names(ft)]
  if (!length(hit)) return(NULL)
  v <- suppressWarnings(as.numeric(ft[[hit[1L]]]))
  if (!any(is.finite(v))) return(NULL)
  v[mask]
}


#' Resolve the on-disk path of an in-house logIE model
#'
#' The shipped models live in `inst/extdata/models/logIE/`. This
#' helper searches the installed package (`system.file`) and the current working
#' directory (source-tree use: `pkgload::load_all()` / RStudio in the package
#' root) and returns the first existing model file, or `NULL` when none is found
#' (callers then fall back to the LogP mock).
#'
#' @param pol_abbr `"pos"` or `"neg"`.
#' @return Full path to the `.rds` model, or `NULL`.
#' @keywords internal
#' @noRd
.nts_resolve_logIE_model_path <- function(pol_abbr) {
  fn <- sprintf("logIE_model_inhouse_%s.rds", pol_abbr)
  cand <- c(
    .nts_extdata(c(file.path("models", "logIE", fn), file.path("logIE_models", fn), fn)),
    system.file("data", "logIE_models", fn, package = "NTSminer"),
    file.path("data", "logIE_models", fn),
    system.file("data", fn, package = "NTSminer"),
    file.path("data", fn)
  )
  cand <- cand[!is.na(cand) & nzchar(cand) & file.exists(cand)]
  if (length(cand)) cand[1L] else NULL
}


#' Resolve the on-disk path of an in-house RT model
#'
#' Mirror of `.nts_resolve_logIE_model_path()` for the retention-time models in
#' `inst/extdata/models/rt/`. Returns the first existing
#' `rt_model_inhouse_<pol_abbr>.rds`, or `NULL` when none is found.
#'
#' @param pol_abbr `"pos"` or `"neg"`.
#' @return Full path to the `.rds` model, or `NULL`.
#' @keywords internal
#' @noRd
.nts_resolve_rt_model_path <- function(pol_abbr) {
  fn <- sprintf("rt_model_inhouse_%s.rds", pol_abbr)
  cand <- c(
    .nts_extdata(c(file.path("models", "rt", fn), file.path("rt_models", fn), fn)),
    system.file("data", "rt_models", fn, package = "NTSminer"),
    file.path("data", "rt_models", fn),
    system.file("data", fn, package = "NTSminer"),
    file.path("data", fn)
  )
  cand <- cand[!is.na(cand) & nzchar(cand) & file.exists(cand)]
  if (length(cand)) cand[1L] else NULL
}


.predict_rie_df <- function(ft, default_rie = 1.0,
                              logIE_model = NULL,
                              polarity = c("auto", "positive", "negative"),
                              smiles_col = NULL) {
  polarity <- match.arg(polarity)

  # Try to resolve SMILES
  smi <- NULL
  if (!is.null(smiles_col) && smiles_col %in% names(ft)) {
    smi <- as.character(ft[[smiles_col]])
  } else {
    cand <- c("SMILES", "smiles", "Smiles", "CanonicalSMILES", "canonical_smiles")
    hit <- cand[cand %in% names(ft)]
    if (length(hit)) {
      smi <- as.character(ft[[hit[1L]]])
      smiles_col <- hit[1L]
    }
  }

  rie <- NULL
  smi_ok <- if (is.null(smi)) rep(FALSE, nrow(ft)) else
    !is.na(smi) & nzchar(trimws(as.character(smi)))

  # Try logIE model if SMILES available
  if (!is.null(smi) && any(smi_ok)) {
    if (is.null(logIE_model) && polarity != "auto") {
      pol_abbr <- if (polarity == "positive") "pos" else "neg"
      model_path <- .nts_resolve_logIE_model_path(pol_abbr)
      if (!is.null(model_path)) {
        logIE_model <- tryCatch(load_logIE_model(model_path), error = function(e) NULL)
      }
    }
    if (!is.null(logIE_model) && inherits(logIE_model, "nts_logIE_model")) {
      # The in-house logIE models use the *measured* retention time (and, for the
      # positive polarity, the organic fraction of the mobile phase) as features.
      # Pass the values from the feature table so those columns carry real data
      # instead of being imputed with 0. predict_logIE() derives MW from the
      # SMILES itself when the caller does not supply it. Both vectors must cover
      # the valid SMILES only -- predict_logIE() builds one row per valid SMILES.
      rt_vals <- .smq_first_numeric(ft, c("rt", "RT", "rt_avg", "RT_avg"), smi_ok)
      org_vals <- .smq_first_numeric(
        ft, c("Organic_Percentage", "organic_pct", "organic_percentage"), smi_ok
      )
      pred <- tryCatch(
        predict_logIE(logIE_model, smi, rt = rt_vals, organic_pct = org_vals),
        error = function(e) NULL
      )
      if (!is.null(pred)) {
        # predict_logIE() drops NA / empty SMILES rows (via compute_fingerprints),
        # so its output can be shorter than ft. Re-expand it against the same
        # validity mask instead of silently discarding every prediction and
        # falling back to the mock values.
        pli <- pred$predicted_logIE
        if (length(pli) == sum(smi_ok) && sum(smi_ok) > 0L) {
          full <- rep(NA_real_, nrow(ft))
          full[smi_ok] <- pli
          rie <- 10^full
          rie[!is.finite(rie) | rie <= 0] <- NA_real_
          rie[is.na(rie)] <- default_rie
        } else {
          warning("predict_logIE() returned ", length(pli), " predictions for ",
                  sum(smi_ok), " valid SMILES; using fallback RIE values.",
                  call. = FALSE)
        }
      }
    }
  }

  # Fallback: LogP mock
  if (is.null(rie)) {
    rie <- rep(default_rie, nrow(ft))
    logp_hits <- intersect(c("LogP", "logP", "logp", "cLogP", "log_p"), names(ft))
    if (length(logp_hits)) {
      lp <- ft[[logp_hits[1L]]]
      ok <- is.finite(lp)
      if (any(ok)) {
        lp0 <- mean(lp[ok], na.rm = TRUE)
        if (!is.finite(lp0)) lp0 <- 0
        rie[ok] <- default_rie + 0.08 * (lp[ok] - lp0)
        rie[ok] <- pmax(0.01, rie[ok])
      }
    }
  }

  ft$rie <- rie
  ft
}


#' Predict relative ionization efficiency (RIE)
#'
#' Computes RIE for each feature. When a logIE model is available and
#' SMILES are present in the feature table, RIE is derived from predicted
#' log10 ionization efficiency (`rie = 10^pred_logIE`).
#' Falls back to a LogP-based linear approximation when no model or
#' SMILES column is available.
#'
#' @param object A `data.frame`/`data.table` of features or an [NTSProject].
#' @param default_rie Positive scalar baseline RIE used for fallback or
#'   invalid predictions.
#' @param logIE_model Optional `nts_logIE_model` object. When `NULL` (default),
#'   attempts to auto-load the in-house model from
#'   `extdata/models/logIE/logIE_model_inhouse_{pos,neg}.rds` (resolved against the
#'   installed package and the source tree). Falls back to a LogP-based mock
#'   when no model can be found.
#' @param polarity `"auto"` (default) uses `"positive"`; set explicitly to
#'   `"positive"` or `"negative"` to select a specific model.
#' @param smiles_col Optional character: name of the SMILES column in the
#'   feature table. Auto-detected from common names when `NULL`.
#'
#' @return Updated table or [NTSProject] with column `rie`.
#' @export
predict_rie <- function(object,
                         default_rie = 1.0,
                         logIE_model = NULL,
                         polarity = c("auto", "positive", "negative"),
                         smiles_col = NULL) {
  polarity <- match.arg(polarity)

  # Resolve polarity from object metadata when auto
  if (polarity == "auto") {
    if (inherits(object, "NTSProject")) {
      obj_pol <- tryCatch(methods::slot(object, "polarity"), error = function(e) NULL)
      if (!is.null(obj_pol) && is.character(obj_pol) && nzchar(obj_pol)) {
        polarity <- tolower(obj_pol)
        if (!polarity %in% c("positive", "negative")) polarity <- "positive"
      }
    } else if (is.data.frame(object) || data.table::is.data.table(object)) {
      obj_pol <- attr(object, "polarity", exact = TRUE)
      if (!is.null(obj_pol) && is.character(obj_pol)) {
        polarity <- tolower(obj_pol)
        if (!polarity %in% c("positive", "negative")) polarity <- "positive"
      }
    }
  }

  if (inherits(object, "NTSProject")) {
    ft <- get_feature_table(object)
    if (!nrow(ft)) stop("feature_table is empty", call. = FALSE)
    out <- .predict_rie_df(ft, default_rie = default_rie,
                            logIE_model = logIE_model,
                            polarity = polarity,
                            smiles_col = smiles_col)
    object <- set_feature_table(object, out)
    methods::slot(object, "log") <- .nts_log_append(
      methods::slot(object, "log"), "predict_rie",
      list(default_rie = default_rie, polarity = polarity,
           logIE_model = !is.null(logIE_model))
    )
    return(object)
  }
  if (!is.data.frame(object) && !data.table::is.data.table(object)) {
    stop("`object` must be a feature table or NTSProject", call. = FALSE)
  }
  .predict_rie_df(as.data.frame(object, stringsAsFactors = FALSE),
                   default_rie = default_rie,
                   logIE_model = logIE_model,
                   polarity = polarity,
                   smiles_col = smiles_col)
}


#' @noRd
.calculate_concentration_df <- function(ft,
                                        is_amount = 100,
                                        is_intensity = 100,
                                        is_match = NULL) {
  if (!"rie" %in% names(ft)) {
    stop("Column `rie` not found; run predict_rie() first.", call. = FALSE)
  }
  int <- .smq_row_intensity(ft)
  denom <- pmax(as.numeric(ft$rie), 1e-12)
  def_amt <- as.numeric(is_amount)[1L]
  def_int <- as.numeric(is_intensity)[1L]
  if (!is.finite(def_int) || def_int <= 0) {
    stop("`is_intensity` must be positive and finite.", call. = FALSE)
  }
  if (!is.finite(def_amt) || def_amt < 0) {
    stop("`is_amount` must be finite and non-negative.", call. = FALSE)
  }
  amt_vec <- rep(def_amt, nrow(ft))
  int_vec <- rep(def_int, nrow(ft))
  if (!is.null(is_match) &&
      all(c("is_name", "mean_intensity") %in% names(is_match)) &&
      "is_anchor" %in% names(ft)) {
    amt_lk <- if ("amount" %in% names(is_match)) {
      stats::setNames(as.numeric(is_match$amount), is_match$is_name)
    } else {
      NULL
    }
    int_lk <- stats::setNames(
      as.numeric(is_match$mean_intensity),
      is_match$is_name
    )
    for (i in seq_len(nrow(ft))) {
      lab <- ft$is_anchor[i]
      if (is.na(lab) || !nzchar(lab)) {
        next
      }
      if (!is.null(amt_lk) && lab %in% names(amt_lk) && is.finite(amt_lk[[lab]])) {
        amt_vec[i] <- amt_lk[[lab]]
      }
      if (lab %in% names(int_lk) && is.finite(int_lk[[lab]]) && int_lk[[lab]] > 0) {
        int_vec[i] <- int_lk[[lab]]
      }
    }
  }
  ft$estimated_conc <- int / denom * (amt_vec / int_vec)
  ft
}


#' Estimate concentrations from abundances and RIE
#'
#' Uses `estimated_conc = intensity / rie * (is_amount / is_intensity)` where
#' `intensity` is the row mean of sample columns (or column `intensity` if
#' present). Defaults `is_amount = 100` and `is_intensity = 100` cancel so the
#' MVP reduces to `intensity / rie`; adjust both for internal-standard
#' normalization.
#'
#' @param object Feature table or [NTSProject] (must already contain `rie`).
#' @param is_amount Default known IS amount when not in `is_match$amount`.
#' @param is_intensity Default IS response when not in `is_match$mean_intensity`.
#' @param is_match Optional table from [match_internal_standards()]; when
#'   `is_anchor` is on the feature table (after [normalize_by_IS()]), per-row
#'   `amount` / `mean_intensity` are used.
#'
#' @return Updated object with column `estimated_conc`.
#' @export
#'
#' @seealso [run_semi_quantification], [normalize_by_IS]
calculate_concentration <- function(object,
                                    is_amount = 100,
                                    is_intensity = 100,
                                    is_match = NULL) {
  if (inherits(object, "NTSProject")) {
    ft <- get_feature_table(object)
    if (!nrow(ft)) {
      stop("feature_table is empty", call. = FALSE)
    }
    out <- .calculate_concentration_df(
      ft,
      is_amount = is_amount,
      is_intensity = is_intensity,
      is_match = is_match
    )
    object <- set_feature_table(object, out)
    methods::slot(object, "log") <- .nts_log_append(
      methods::slot(object, "log"),
      "calculate_concentration",
      list(is_amount = is_amount, is_intensity = is_intensity)
    )
    return(object)
  }
  if (!is.data.frame(object) && !data.table::is.data.table(object)) {
    stop("`object` must be a feature table or NTSProject", call. = FALSE)
  }
  .calculate_concentration_df(
    as.data.frame(object, stringsAsFactors = FALSE),
    is_amount = is_amount,
    is_intensity = is_intensity,
    is_match = is_match
  )
}

#' Known SMILES column names (case-insensitive first match)
#' @noRd
.tox_smiles_col <- function(names) {
  cand <- c("SMILES", "smiles", "Smiles", "InChI", "inchi")
  hit <- cand[cand %in% names]
  if (length(hit)) {
    return(hit[1L])
  }
  NULL
}


#' Pull SMILES vector aligned to feature rows
#'
#' Uses `feature_table[[smiles_col]]` or merges `db_annotations` on `name`.
#' @noRd
.tox_resolve_smiles <- function(ft, ann) {
  sc <- .tox_smiles_col(names(ft))
  smiles <- rep(NA_character_, nrow(ft))
  if (!is.null(sc)) {
    smiles <- as.character(ft[[sc]])
    return(list(smiles = smiles, smiles_col = sc))
  }
  if (nrow(ann) && "name" %in% names(ann) && "name" %in% names(ft)) {
    sca <- .tox_smiles_col(names(ann))
    if (!is.null(sca)) {
      u <- unique(ann[, c("name", sca), drop = FALSE])
      m <- match(ft$name, u$name)
      ok <- !is.na(m)
      smiles[ok] <- as.character(u[[sca]][m[ok]])
      return(list(smiles = smiles, smiles_col = sca))
    }
  }
  list(smiles = smiles, smiles_col = NULL)
}


#' Mock LogP from SMILES string length (deterministic placeholder)
#' @noRd
.tox_mock_logp <- function(smiles) {
  n <- nchar(smiles, type = "bytes", allowNA = TRUE)
  ifelse(is.na(n), NA_real_, pmin(6, pmax(-2, 0.12 * n - 1.5)))
}


#' Deterministic mock toxicity class from SMILES
#' @noRd
.tox_mock_class <- function(smiles) {
  n <- suppressWarnings(nchar(smiles, type = "bytes", allowNA = TRUE))
  ifelse(is.na(n), NA_character_, ifelse(n %% 2L == 0L, "High", "Low"))
}


#' @noRd
.nts_rcdk_physchem_descriptor_classes <- function() {
  c(
    "org.openscience.cdk.qsar.descriptors.molecular.ALogPDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.MolecularWeightDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor"
  )
}


#' @noRd
.nts_ann_resolve_smiles_physchem <- function(ann_dt) {
  n <- nrow(ann_dt)
  empty <- rep(NA_character_, n)
  col_norm <- function(nm) {
    if (!nm %in% names(ann_dt)) {
      return(empty)
    }
    vapply(as.character(ann_dt[[nm]]), .nts_normalize_smiles_string, character(1L), USE.NAMES = FALSE)
  }
  inf <- col_norm("inferred_smiles")
  s1 <- col_norm("SMILES")
  s2 <- col_norm("smiles")
  out <- inf
  miss <- is.na(out) | !nzchar(out)
  out[miss] <- s1[miss]
  miss <- is.na(out) | !nzchar(out)
  out[miss] <- s2[miss]
  out
}


#' Pure physicochemical proxy engine (rcdk descriptors)
#'
#' Receives a `data.table` with at least `name` and `smiles` (one SMILES string
#' per row). Parses each structure with **rcdk** under row-safe `tryCatch`;
#' invalid SMILES yield `NA` for all descriptor columns without aborting the
#' batch. Never mutates the input table by reference.
#'
#' @param dt `data.table` (or coercible) with columns `name`, `smiles`.
#'
#' @return `data.table` with columns `name`, `ALogP_proxy`, `TPSA`,
#'   `MolecularWeight`, `HBondDonorCount`, `HBondAcceptorCount`.
#' @noRd
.nts_engine_physchem <- function(dt) {
  .nts_check(NULL, "data.table_ns")
  x <- data.table::copy(data.table::as.data.table(dt))
  if (!all(c("name", "smiles") %in% names(x))) {
    stop(".nts_engine_physchem() requires columns `name` and `smiles`.", call. = FALSE)
  }
  n <- nrow(x)
  nm <- as.character(x$name)
  out <- data.table::data.table(
    name = nm,
    ALogP_proxy = rep(NA_real_, n),
    TPSA = rep(NA_real_, n),
    MolecularWeight = rep(NA_real_, n),
    HBondDonorCount = rep(NA_integer_, n),
    HBondAcceptorCount = rep(NA_integer_, n)
  )
  if (!n) {
    return(out[])
  }
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    warning(
      ".nts_engine_physchem: package 'rcdk' is not installed - returning NA descriptors.",
      call. = FALSE
    )
    return(out[])
  }

  sm_norm <- vapply(as.character(x$smiles), .nts_normalize_smiles_string, character(1L), USE.NAMES = FALSE)
  mols <- vector("list", n)
  parse_ok <- rep(FALSE, n)
  for (i in seq_len(n)) {
    if (is.na(sm_norm[i]) || !nzchar(sm_norm[i])) {
      next
    }
    mols[[i]] <- tryCatch(
      .nts_mol_first_from_smiles_rcdk(sm_norm[i]),
      error = function(e) NULL
    )
    parse_ok[i] <- !is.null(mols[[i]])
  }
  dclasses <- .nts_rcdk_physchem_descriptor_classes()
  idx <- which(parse_ok)
  if (!length(idx)) {
    return(out[])
  }

  al <- rep(NA_real_, n)
  tp <- rep(NA_real_, n)
  mw <- rep(NA_real_, n)
  hd <- rep(NA_integer_, n)
  ha <- rep(NA_integer_, n)

  batch <- mols[idx]
  df <- tryCatch(
    rcdk::eval.desc(batch, dclasses),
    error = function(e) NULL
  )
  if (!is.null(df) && nrow(df) == length(idx)) {
    nc <- min(5L, ncol(df))
    for (j in seq_len(nrow(df))) {
      ii <- idx[j]
      vals_num <- suppressWarnings(as.numeric(df[j, seq_len(nc), drop = TRUE]))
      if (length(vals_num) >= 5L) {
        al[ii] <- vals_num[[1L]]
        tp[ii] <- vals_num[[2L]]
        mw[ii] <- vals_num[[3L]]
        hd[ii] <- suppressWarnings(as.integer(round(vals_num[[4L]])))
        ha[ii] <- suppressWarnings(as.integer(round(vals_num[[5L]])))
      }
    }
  } else {
    for (ii in idx) {
      mol <- mols[[ii]]
      df1 <- tryCatch(
        rcdk::eval.desc(mol, dclasses),
        error = function(e) NULL
      )
      if (is.null(df1) || !nrow(df1)) {
        next
      }
      nc <- min(5L, ncol(df1))
      vals_num <- suppressWarnings(as.numeric(df1[1L, seq_len(nc), drop = TRUE]))
      if (length(vals_num) >= 5L) {
        al[ii] <- vals_num[[1L]]
        tp[ii] <- vals_num[[2L]]
        mw[ii] <- vals_num[[3L]]
        hd[ii] <- suppressWarnings(as.integer(round(vals_num[[4L]])))
        ha[ii] <- suppressWarnings(as.integer(round(vals_num[[5L]])))
      }
    }
  }

  out[, `:=`(
    ALogP_proxy = al,
    TPSA = tp,
    MolecularWeight = mw,
    HBondDonorCount = hd,
    HBondAcceptorCount = ha
  )]
  out[]
}


#' Calculate fast physicochemical / PMT proxy descriptors (rcdk)
#'
#' Builds a working table from `@ms1_annotations`, preferring
#' `inferred_smiles`, then `SMILES` / `smiles` (after `.nts_normalize_smiles_string()`).
#' Descriptor values are computed by `.nts_engine_physchem()` and merged into
#' the dedicated `@physchem_proxy` table (keyed by `name`) so `@ms1_annotations`
#' stays reserved for chromatographic / KMD / flagging semantics.
#'
#' @param project An [NTSProject] with non-empty `@ms1_annotations` and column
#'   `name` (SMILES are resolved from annotation / SMILES-like columns for CDK).
#'
#' @return Updated `project` with `ALogP_proxy`, `TPSA`, `MolecularWeight`,
#'   `HBondDonorCount`, and `HBondAcceptorCount` on `@physchem_proxy` (where
#'   rcdk/Java succeeds).
#'
#' @details
#' Requires CRAN package **rcdk** and a working Java runtime for CDK. If **rcdk**
#' is missing or `eval.desc` fails globally, the function still returns the
#' project with descriptor columns set to `NA` and emits a `warning()`.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE) && requireNamespace("rcdk", quietly = TRUE)
#' d <- tempfile("nts_pc_ex"); dir.create(d)
#' p <- NTSminer::create_nts_project(d)
#' methods::slot(p, "feature_table") <- data.frame(
#'   name = "F1", mz = 300, rt = 10, stringsAsFactors = FALSE
#' )
#' ann <- NTSminer:::.nts_init_ms1_from_ft(p@feature_table)
#' ann$SMILES <- "CCO"
#' ## `SMILES` is a physchem-tier column, so bypass the step sanitizer and store
#' ## the partitioned slot directly.
#' methods::slot(p, "ms1_annotations") <- NTSminer:::.nts_ms1_partition(ann)
#' p2 <- suppressWarnings(NTSminer::calculate_physchem_proxies(p))
#' stopifnot("ALogP_proxy" %in% names(p2@physchem_proxy))
#' @export
calculate_physchem_proxies <- function(project) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  ann <- get_ms1_annotations(project)
  .nts_check(ann, "non_empty", "name_col", .arg = "ms1_annotations")

  ann_dt <- data.table::copy(data.table::as.data.table(ann))
  ann_dt[, smiles_eff := .nts_ann_resolve_smiles_physchem(ann_dt)]
  work <- ann_dt[, .(name = as.character(name), smiles = smiles_eff)]
  phy <- .nts_engine_physchem(work)

  engine_cols <- setdiff(names(phy), "name")
  base <- unique(ann_dt[, .(name = as.character(name))], by = "name")
  merged <- merge(
    as.data.frame(base, stringsAsFactors = FALSE),
    as.data.frame(phy, stringsAsFactors = FALSE),
    by = "name",
    all.x = TRUE,
    sort = FALSE
  )
  methods::slot(project, "physchem_proxy") <- merged
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "calculate_physchem_proxies",
    list(engine_cols = engine_cols)
  )
  project
}


#' @noRd
.predict_properties_df <- function(ft,
                                   ann,
                                   episuite_wide = NULL) {
  res <- .tox_resolve_smiles(ft, ann)
  smiles <- res$smiles
  has <- !is.na(smiles) & nzchar(smiles)

  if (is.null(episuite_wide) || !nrow(episuite_wide)) {
    ft$LogP_pred <- NA_real_
    if (any(has)) {
      ft$LogP_pred[has] <- .tox_mock_logp(smiles[has])
    }
    ft$Toxicity_Class <- NA_character_
    if (any(has)) {
      ft$Toxicity_Class[has] <- .tox_mock_class(smiles[has])
    }
  } else {
    smi_col <- if ("SMILES" %in% names(episuite_wide)) "SMILES" else "smiles"
    j <- match(trimws(smiles[has]), trimws(as.character(episuite_wide[[smi_col]])))
    if ("LogKow_est" %in% names(episuite_wide)) {
      ft$LogP_pred <- NA_real_
      ft$LogP_pred[has] <- episuite_wide$LogKow_est[j]
    }
    if ("Fish_LC50" %in% names(episuite_wide)) {
      ft$Toxicity_Class <- NA_character_
      ft$Toxicity_Class[has] <- ifelse(
        is.finite(episuite_wide$Fish_LC50[j]),
        paste0("Fish_LC50=", episuite_wide$Fish_LC50[j]),
        NA_character_
      )
    }
    epi_cols <- setdiff(names(episuite_wide), smi_col)
    for (cl in epi_cols) {
      ft[[cl]] <- NA
      ft[[cl]][has] <- episuite_wide[[cl]][j]
    }
  }

  ann_out <- ann
  if (nrow(ann_out) && "name" %in% names(ann_out) && "name" %in% names(ft)) {
    addc <- unique(c(
      "LogP_pred", "Toxicity_Class",
      grep("^(LogKow|LogKoc|Biowin|LogBCF|Fish_|Daphnid_|Algae_|episuite_)", names(ft), value = TRUE)
    ))
    addc <- intersect(addc, names(ft))
    if (length(addc)) {
      sub <- unique(ft[, c("name", addc), drop = FALSE])
      ann_out <- merge(ann_out, sub, by = "name", all.x = TRUE, sort = FALSE)
    }
  }

  list(feature_table = ft, annotations = ann_out)
}


#' Predict environmental fate / ecotox (EPI Suite)
#'
#' Resolves SMILES from an [NTSProject] (`master_smiles`, `SMILES`, `library_SMILES`,
#' etc.) or accepts a character vector / `data.frame` with a `SMILES` column,
#' then calls the **EPI Suite** web API (`episuite.dev`) via [predict_episuite].
#'
#' Results are written to `@physchem_proxy` (wide EPI + ECOSAR columns) when
#' `object` is an `NTSProject`. Legacy columns `LogP_pred` (from `LogKow_est`)
#' and a short `Toxicity_Class` summary (from Fish LC50) are also filled on
#' `@feature_table` for backward compatibility.
#'
#' For CDK fast descriptors (ALogP, TPSA, ...) without network calls, use
#' [calculate_physchem_proxies].
#'
#' @param object An [NTSProject], a character vector of SMILES, or a `data.frame`
#'   with column `SMILES` (or `smiles`).
#' @param progress If `TRUE`, message EPI Suite batch progress.
#' @param write_to_physchem_proxy For `NTSProject` only: merge EPI columns into
#'   `@physchem_proxy` (default `TRUE`).
#' @param ... Passed to [predict_episuite] (`base_url`, `max_retries`, etc.).
#'
#' @return For `NTSProject`, the updated project. For SMILES input, a wide
#'   `data.frame` from [predict_episuite].
#'
#' @seealso [predict_admet_risk], [predict_episuite], [calculate_physchem_proxies]
#'
#' @export
predict_properties <- function(object,
                               progress = TRUE,
                               write_to_physchem_proxy = TRUE,
                               ...) {
  if (is.character(object)) {
    return(predict_episuite(object, progress = progress, ...))
  }
  if (is.data.frame(object) || data.table::is.data.table(object)) {
    df <- as.data.frame(object, stringsAsFactors = FALSE)
    smi <- .tox_smiles_col(names(df))
    if (is.null(smi)) {
      stop("data.frame input needs a SMILES/smiles column.", call. = FALSE)
    }
    return(predict_episuite(df[[smi]], progress = progress, ...))
  }
  if (!inherits(object, "NTSProject")) {
    stop("`object` must be NTSProject, SMILES vector, or data.frame.", call. = FALSE)
  }
  ft <- get_feature_table(object)
  if (!nrow(ft)) {
    stop("feature_table is empty", call. = FALSE)
  }
  ann <- methods::slot(object, "db_annotations")
  id_map <- .nts_toxpred_smiles_from_project(object)
  if (!nrow(id_map)) {
    warning("predict_properties: no SMILES on project; skipping EPI Suite.", call. = FALSE)
    return(object)
  }
  epi <- predict_episuite(id_map$SMILES, progress = progress, ...)
  if (isTRUE(write_to_physchem_proxy)) {
    object <- .nts_toxpred_write_physchem(object, epi, id_map)
  }
  pr <- .predict_properties_df(ft, ann, episuite_wide = epi)
  object <- set_feature_table(object, pr$feature_table)
  methods::slot(object, "db_annotations") <- pr$annotations
  methods::slot(object, "log") <- .nts_log_append(
    methods::slot(object, "log"),
    "predict_properties",
    list(n_smiles = nrow(id_map), backend = "EPI_Suite")
  )
  object
}
