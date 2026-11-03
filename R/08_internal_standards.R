# Internal-standard matching, normalisation, and semi-quantification pipeline.

#' Read internal-standard metadata from CSV
#'
#' Expected columns (flexible names): **name**, **mz** (or mass / precursor_mz),
#' **rt** (seconds, same as xcms `feature_table$rt`). Optional per-row
#' **mz_tol_ppm**, **rt_tol_sec**, **amount** (known spiked amount for
#' [calculate_concentration()]).
#'
#' @param path Path to a CSV file.
#' @param mz_tol_ppm Default m/z tolerance (ppm) when not given per row.
#' @param rt_tol_sec Default RT tolerance (seconds) when not given per row.
#'
#' @return A `data.frame` with columns `is_name`, `mz`, `rt`, and optional
#'   `mz_tol_ppm`, `rt_tol_sec`, `amount`.
#' @export
#'
#' @seealso [match_internal_standards], [normalize_by_IS], [run_semi_quantification]
read_internal_standards <- function(path,
                                    mz_tol_ppm = 10,
                                    rt_tol_sec = 45) {
  path <- path.expand(as.character(path)[1L])
  if (!file.exists(path)) {
    stop("Internal standards file not found: ", path, call. = FALSE)
  }
  raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(raw)) {
    stop("Internal standards file is empty: ", path, call. = FALSE)
  }
  nm <- names(raw)
  name_col <- .nts_resolve_csv_col(nm, c("name", "is_name", "compound", "label", "Lab_ID"))
  mz_col <- .nts_resolve_csv_col(nm, c("mz", "mass", "precursor_mz", "target_mz", "m/z"))
  rt_col <- .nts_resolve_csv_col(nm, c("rt", "retention_time", "retentiontime", "tr", "RT"))
  if (is.null(name_col) || is.null(mz_col) || is.null(rt_col)) {
    stop(
      "CSV needs name, mz (or mass), and rt columns. Found: ",
      paste(nm, collapse = ", "),
      call. = FALSE
    )
  }
  out <- data.frame(
    is_name = trimws(as.character(raw[[name_col]])),
    mz = suppressWarnings(as.numeric(raw[[mz_col]])),
    rt = suppressWarnings(as.numeric(raw[[rt_col]])),
    stringsAsFactors = FALSE
  )
  ppm_col <- .nts_resolve_csv_col(nm, c("mz_tol_ppm", "ppm_tol", "mz_tol"))
  rt_tol_col <- .nts_resolve_csv_col(nm, c("rt_tol_sec", "rt_tol", "rt_window"))
  amt_col <- .nts_resolve_csv_col(nm, c("amount", "is_amount", "spiked_amount", "conc"))
  out$mz_tol_ppm <- if (!is.null(ppm_col)) {
    suppressWarnings(as.numeric(raw[[ppm_col]]))
  } else {
    NA_real_
  }
  out$rt_tol_sec <- if (!is.null(rt_tol_col)) {
    suppressWarnings(as.numeric(raw[[rt_tol_col]]))
  } else {
    NA_real_
  }
  out$amount <- if (!is.null(amt_col)) {
    suppressWarnings(as.numeric(raw[[amt_col]]))
  } else {
    NA_real_
  }
  def_ppm <- as.numeric(mz_tol_ppm)[1L]
  def_rt <- as.numeric(rt_tol_sec)[1L]
  miss_ppm <- !is.finite(out$mz_tol_ppm)
  miss_rt <- !is.finite(out$rt_tol_sec)
  out$mz_tol_ppm[miss_ppm] <- def_ppm
  out$rt_tol_sec[miss_rt] <- def_rt
  ok <- !is.na(out$is_name) & nzchar(out$is_name) &
    is.finite(out$mz) & out$mz > 0 &
    is.finite(out$rt)
  if (!any(ok)) {
    stop("No valid internal-standard rows after parsing.", call. = FALSE)
  }
  out <- out[ok, , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' @keywords internal
#' @noRd
.nts_is_ppm_error <- function(mz_obs, mz_ref) {
  abs(mz_obs - mz_ref) / pmax(mz_ref, 1e-12) * 1e6
}


#' Match internal standards to features in a feature table
#'
#' Searches `feature_table` rows by **m/z (ppm)** and **RT (seconds)**. When
#' several features match, picks one via `pick` (`max_intensity` or `closest_rt`).
#'
#' @param object An [NTSProject] or a feature `data.frame` with `name`, `mz`, `rt`.
#' @param is_spec A `data.frame` from [read_internal_standards()] or a path to CSV.
#' @param mz_tol_ppm Default ppm tolerance (overridden by per-row values in `is_spec`).
#' @param rt_tol_sec Default RT tolerance in **seconds**.
#' @param pick `"max_intensity"` (default) or `"closest_rt"`.
#' @param abundance_cols Sample intensity columns; auto-detected if `NULL`.
#' @param on_missing `"stop"` (default) or `"warn"` when an IS has no matching feature.
#'
#' @return A `data.frame` with one row per matched IS: `is_name`, `feature_name`,
#'   `mz_ref`, `rt_ref`, `mz_obs`, `rt_obs`, `ppm_error`, `rt_error_sec`,
#'   `mean_intensity`, `n_candidates`, `match_status` (`"ok"` / `"missing"`).
#' @export
#'
#' @seealso [read_internal_standards], [normalize_by_IS], [run_semi_quantification]
match_internal_standards <- function(object,
                                       is_spec,
                                       mz_tol_ppm = 10,
                                       rt_tol_sec = 45,
                                       pick = c("max_intensity", "closest_rt"),
                                       abundance_cols = NULL,
                                       on_missing = c("stop", "warn")) {
  pick <- match.arg(pick)
  on_missing <- match.arg(on_missing)
  if (is.character(is_spec) && length(is_spec) == 1L && nzchar(is_spec)) {
    is_spec <- read_internal_standards(
      is_spec,
      mz_tol_ppm = mz_tol_ppm,
      rt_tol_sec = rt_tol_sec
    )
  }
  if (!is.data.frame(is_spec) || !nrow(is_spec)) {
    stop("`is_spec` must be a non-empty data.frame or CSV path.", call. = FALSE)
  }
  ft <- if (inherits(object, "NTSProject")) {
    get_feature_table(object)
  } else {
    as.data.frame(object, stringsAsFactors = FALSE)
  }
  if (!nrow(ft)) {
    stop("feature_table is empty.", call. = FALSE)
  }
  if (!all(c("name", "mz", "rt") %in% names(ft))) {
    stop("feature_table must contain name, mz, and rt.", call. = FALSE)
  }
  req_is <- c("is_name", "mz", "rt")
  if (!all(req_is %in% names(is_spec))) {
    stop("`is_spec` needs columns: ", paste(req_is, collapse = ", "), call. = FALSE)
  }
  if (is.null(abundance_cols)) {
    abundance_cols <- .smq_abundance_cols(ft)
  }
  if (!length(abundance_cols)) {
    stop("No numeric sample columns detected for intensity ranking.", call. = FALSE)
  }

  def_ppm <- as.numeric(mz_tol_ppm)[1L]
  def_rt <- as.numeric(rt_tol_sec)[1L]
  rows <- vector("list", nrow(is_spec))
  for (i in seq_len(nrow(is_spec))) {
    ref_name <- as.character(is_spec$is_name[i])
    ref_mz <- as.numeric(is_spec$mz[i])
    ref_rt <- as.numeric(is_spec$rt[i])
    ppm_tol <- if ("mz_tol_ppm" %in% names(is_spec) && is.finite(is_spec$mz_tol_ppm[i])) {
      is_spec$mz_tol_ppm[i]
    } else {
      def_ppm
    }
    rt_tol <- if ("rt_tol_sec" %in% names(is_spec) && is.finite(is_spec$rt_tol_sec[i])) {
      is_spec$rt_tol_sec[i]
    } else {
      def_rt
    }
    ppm_err <- .nts_is_ppm_error(as.numeric(ft$mz), ref_mz)
    rt_err <- abs(as.numeric(ft$rt) - ref_rt)
    cand <- which(is.finite(ppm_err) & ppm_err <= ppm_tol & is.finite(rt_err) & rt_err <= rt_tol)
    if (!length(cand)) {
      rows[[i]] <- data.frame(
        is_name = ref_name,
        feature_name = NA_character_,
        mz_ref = ref_mz,
        rt_ref = ref_rt,
        mz_obs = NA_real_,
        rt_obs = NA_real_,
        ppm_error = NA_real_,
        rt_error_sec = NA_real_,
        mean_intensity = NA_real_,
        amount = if ("amount" %in% names(is_spec)) is_spec$amount[i] else NA_real_,
        n_candidates = 0L,
        match_status = "missing",
        stringsAsFactors = FALSE
      )
      next
    }
    if (identical(pick, "closest_rt")) {
      j <- cand[which.min(rt_err[cand])]
    } else {
      inten <- rowMeans(ft[cand, abundance_cols, drop = FALSE], na.rm = TRUE)
      j <- cand[which.max(inten)]
    }
    rows[[i]] <- data.frame(
      is_name = ref_name,
      feature_name = as.character(ft$name[j]),
      mz_ref = ref_mz,
      rt_ref = ref_rt,
      mz_obs = as.numeric(ft$mz[j]),
      rt_obs = as.numeric(ft$rt[j]),
      ppm_error = ppm_err[j],
      rt_error_sec = rt_err[j],
      mean_intensity = rowMeans(ft[j, abundance_cols, drop = FALSE], na.rm = TRUE),
      amount = if ("amount" %in% names(is_spec)) is_spec$amount[i] else NA_real_,
      n_candidates = length(cand),
      match_status = "ok",
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  miss <- out$match_status == "missing"
  if (any(miss)) {
    msg <- paste0(
      "Internal standard(s) not matched in feature_table: ",
      paste(out$is_name[miss], collapse = ", ")
    )
    if (identical(on_missing, "stop")) {
      stop(msg, call. = FALSE)
    }
    warning(msg, call. = FALSE)
  }
  out <- out[!miss, , drop = FALSE]
  if (!nrow(out)) {
    stop("No internal standards matched.", call. = FALSE)
  }
  out
}


#' @keywords internal
#' @noRd
.nts_is_anchor_index <- function(feature_rt, is_rt) {
  if (!length(is_rt)) {
    return(rep(NA_integer_, length(feature_rt)))
  }
  vapply(
    feature_rt,
    function(rt) {
      if (!is.finite(rt)) {
        return(NA_integer_)
      }
      which.min(abs(is_rt - rt))
    },
    integer(1L)
  )
}


#' @keywords internal
#' @noRd
.nts_is_loess_predict <- function(rt_is, y_is, rt_query) {
  ok <- is.finite(rt_is) & is.finite(y_is) & y_is > 0
  if (sum(ok) < 2L) {
    return(rep(NA_real_, length(rt_query)))
  }
  rt_u <- rt_is[ok]
  y_u <- y_is[ok]
  ord <- order(rt_u)
  rt_u <- rt_u[ord]
  y_u <- y_u[ord]
  if (length(unique(rt_u)) < 2L) {
    return(rep(y_u[1L], length(rt_query)))
  }
  if (length(rt_u) >= 4L) {
    fit <- tryCatch(
      stats::loess(
        log(y_u) ~ rt_u,
        span = 0.75,
        control = stats::loess.control(surface = "direct")
      ),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      pred <- tryCatch(
        as.numeric(exp(stats::predict(fit, newdata = data.frame(rt_u = rt_query)))),
        error = function(e) NULL
      )
      if (!is.null(pred) && any(is.finite(pred))) {
        pred[!is.finite(pred) | pred <= 0] <- NA_real_
        return(pred)
      }
    }
  }
  pred <- stats::approx(
    x = rt_u,
    y = y_u,
    xout = rt_query,
    rule = 2,
    ties = base::ordered
  )$y
  pred[!is.finite(pred) | pred <= 0] <- NA_real_
  pred
}


#' @keywords internal
#' @noRd
.nts_is_rt_segment_index <- function(feature_rt, is_rt) {
  if (length(is_rt) < 2L) {
    return(rep(1L, length(feature_rt)))
  }
  ord <- order(is_rt)
  rt_s <- is_rt[ord]
  mids <- c(-Inf, (rt_s[-length(rt_s)] + rt_s[-1L]) / 2, Inf)
  findInterval(feature_rt, mids)
}


#' Extended IS normalisation (called from [normalize_by_IS])
#' @keywords internal
#' @noRd
.normalize_by_IS_df <- function(ft,
                                is_names,
                                abundance_cols = NULL,
                                scale_factor = 1e6,
                                method = c("geom_mean", "nearest_rt", "loess", "rt_segment"),
                                is_match = NULL,
                                rt_breaks = NULL) {
  method <- match.arg(method)
  stopifnot("name" %in% names(ft))
  is_names <- unique(as.character(is_names))
  is_names <- is_names[!is.na(is_names) & nzchar(is_names)]
  if (!length(is_names)) {
    stop("No internal-standard feature names provided.", call. = FALSE)
  }
  is_idx <- match(is_names, as.character(ft$name))
  if (any(is.na(is_idx))) {
    missing <- is_names[is.na(is_idx)]
    stop(
      "IS features not found in feature_table: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(abundance_cols)) {
    abundance_cols <- .smq_abundance_cols(ft)
  }
  if (!length(abundance_cols)) {
    stop("No numeric sample columns detected.", call. = FALSE)
  }

  is_mat <- as.matrix(ft[is_idx, abundance_cols, drop = FALSE])
  is_mat[is_mat <= 0 | !is.finite(is_mat)] <- NA_real_
  is_rt <- as.numeric(ft$rt[is_idx])
  names(is_rt) <- is_names

  is_label <- is_names
  if (!is.null(is_match) && all(c("feature_name", "is_name") %in% names(is_match))) {
    lab <- is_match$is_name[match(is_names, is_match$feature_name)]
    ok_lab <- !is.na(lab)
    is_label[ok_lab] <- lab[ok_lab]
  }

  sample_mat <- as.matrix(ft[, abundance_cols, drop = FALSE])
  n <- nrow(ft)
  anchor_lab <- rep(NA_character_, n)

  if (identical(method, "geom_mean")) {
    if (nrow(is_mat) == 1L) {
      divisor <- is_mat[1L, ]
    } else {
      divisor <- exp(colMeans(log(is_mat), na.rm = TRUE))
    }
    divisor[!is.finite(divisor) | divisor <= 0] <- NA_real_
    if (all(is.na(divisor))) {
      stop("All IS responses are zero or NA; cannot normalise.", call. = FALSE)
    }
    normalised <- sweep(sample_mat, 2, divisor, FUN = "/") * scale_factor
    anchor_lab[] <- paste(is_label, collapse = "+")
  } else if (identical(method, "nearest_rt")) {
    ai <- .nts_is_anchor_index(as.numeric(ft$rt), is_rt)
    if (any(is.na(ai))) {
      stop("Cannot assign IS anchor (invalid feature RT).", call. = FALSE)
    }
    normalised <- matrix(NA_real_, n, length(abundance_cols))
    colnames(normalised) <- abundance_cols
    for (j in seq_along(abundance_cols)) {
      col <- abundance_cols[j]
      div <- is_mat[ai, j]
      div[!is.finite(div) | div <= 0] <- NA_real_
      normalised[, j] <- sample_mat[, j] / div * scale_factor
    }
    anchor_lab <- is_label[ai]
  } else if (identical(method, "loess")) {
    normalised <- matrix(NA_real_, n, length(abundance_cols))
    colnames(normalised) <- abundance_cols
    ft_rt <- as.numeric(ft$rt)
    for (j in seq_along(abundance_cols)) {
      col <- abundance_cols[j]
      div <- .nts_is_loess_predict(is_rt, is_mat[, j], ft_rt)
      if (all(is.na(div))) {
        div <- if (nrow(is_mat) == 1L) {
          is_mat[1L, j]
        } else {
          exp(mean(log(is_mat[, j]), na.rm = TRUE))
        }
      }
      normalised[, j] <- sample_mat[, j] / div * scale_factor
    }
    ai <- .nts_is_anchor_index(ft_rt, is_rt)
    anchor_lab <- is_label[ai]
  } else if (identical(method, "rt_segment")) {
    if (!is.null(rt_breaks)) {
      br <- sort(as.numeric(rt_breaks))
      br <- br[is.finite(br)]
      mids <- c(-Inf, br, Inf)
      seg <- findInterval(as.numeric(ft$rt), mids)
      seg_is <- findInterval(is_rt, mids)
      # Map each feature's segment number back to the internal standard sitting
      # in that segment. `seg_is[seg]` indexed the per-IS vector by a segment
      # *number* (wrong, and out of bounds whenever #segments != #IS);
      # match(seg, seg_is) is the intended lookup. Features in a segment with no
      # IS get NA and are normalised to NA below.
      ai <- match(seg, seg_is)
    } else {
      ai <- .nts_is_rt_segment_index(as.numeric(ft$rt), is_rt)
      ord <- order(is_rt)
      ai <- ord[ai]
    }
    normalised <- matrix(NA_real_, n, length(abundance_cols))
    colnames(normalised) <- abundance_cols
    for (j in seq_along(abundance_cols)) {
      col <- abundance_cols[j]
      div <- is_mat[ai, j]
      div[!is.finite(div) | div <= 0] <- NA_real_
      normalised[, j] <- sample_mat[, j] / div * scale_factor
    }
    anchor_lab <- is_label[ai]
  }

  ft[, abundance_cols] <- as.data.frame(normalised, stringsAsFactors = FALSE)
  ft$is_normalised <- TRUE
  ft$is_normalise_method <- method
  ft$is_anchor <- anchor_lab

  message(
    "> Internal-standard normalisation (", method, "): ",
    length(is_names), " IS feature(s), ",
    length(abundance_cols), " sample column(s)."
  )
  ft
}


#' Semi-quantification pipeline: match IS -> normalise -> RIE -> concentration
#'
#' Runs, in order: [match_internal_standards()] (unless `is_match` supplied),
#' [normalize_by_IS()], [predict_rie()], and [calculate_concentration()].
#' **RIE** uses the current MVP ([predict_rie()]: LogP-based mock or
#' `default_rie`); when a logIE model is available, RIE is derived from
#' predicted ionization efficiency (`rie = 10^pred_logIE`).
#'
#' @param object An [NTSProject] or feature `data.frame`.
#' @param is_csv Path to internal-standard CSV (see [read_internal_standards()]).
#' @param is_spec Optional pre-built IS table; ignored if `is_csv` is set.
#' @param is_match Optional output from [match_internal_standards()]; skips matching.
#' @param mz_tol_ppm,rt_tol_sec,pick,on_missing Passed to matching when `is_match` is `NULL`.
#' @param normalize_method One of `"nearest_rt"` (default for PFAS gradients),
#'   `"geom_mean"`, `"loess"`, `"rt_segment"`.
#' @param rt_breaks Optional RT boundaries (seconds) for `normalize_method = "rt_segment"`.
#' @param scale_factor Passed to [normalize_by_IS()].
#' @param abundance_cols Sample columns; auto-detected if `NULL`.
#' @param predict_rie If `TRUE`, run [predict_rie()] (default `TRUE`).
#' @param default_rie Baseline RIE when no model or SMILES column exist.
#' @param logIE_model Optional `nts_logIE_model`; auto-loaded when `NULL`.
#' @param polarity `"auto"` (default) or `"positive"`/`"negative"`.
#' @param smiles_col SMILES column name; auto-detected when `NULL`.
#' @param calculate_conc If `TRUE`, run [calculate_concentration()] (default `TRUE`).
#' @param is_amount Default IS amount when not in CSV (`amount` column).
#' @param is_intensity Default IS intensity when per-IS `mean_intensity` unavailable.
#'
#' @return For [NTSProject], the updated project with a list entry
#'   `semi_quantification` in `@log` (`is_match`, method, etc.). For a plain
#'   feature table, a list with `feature_table`, `is_match`, and `log`.
#' @export
#'
#' @seealso [read_internal_standards], [match_internal_standards], [normalize_by_IS]
run_semi_quantification <- function(object,
                                    is_csv = NULL,
                                    is_spec = NULL,
                                    is_match = NULL,
                                    mz_tol_ppm = 10,
                                    rt_tol_sec = 45,
                                    pick = c("max_intensity", "closest_rt"),
                                    on_missing = c("stop", "warn"),
                                    normalize_method = c(
                                      "nearest_rt", "geom_mean", "loess", "rt_segment"
                                    ),
                                    rt_breaks = NULL,
                                    scale_factor = 1e6,
                                    abundance_cols = NULL,
                                    predict_rie = TRUE,
                                    default_rie = 1.0,
                                    logIE_model = NULL,
                                    polarity = c("auto", "positive", "negative"),
                                    smiles_col = NULL,
                                    calculate_conc = TRUE,
                                    is_amount = 100,
                                    is_intensity = 100) {
  pick <- match.arg(pick)
  on_missing <- match.arg(on_missing)
  normalize_method <- match.arg(normalize_method)
  polarity <- match.arg(polarity)

  if (is.null(is_match)) {
    spec <- if (!is.null(is_csv) && nzchar(as.character(is_csv)[1L])) {
      read_internal_standards(is_csv, mz_tol_ppm = mz_tol_ppm, rt_tol_sec = rt_tol_sec)
    } else if (!is.null(is_spec)) {
      is_spec
    } else {
      stop("Provide `is_csv`, `is_spec`, or `is_match`.", call. = FALSE)
    }
    is_match <- match_internal_standards(
      object,
      spec,
      mz_tol_ppm = mz_tol_ppm,
      rt_tol_sec = rt_tol_sec,
      pick = pick,
      abundance_cols = abundance_cols,
      on_missing = on_missing
    )
  }

  is_names <- unique(as.character(is_match$feature_name))
  is_names <- is_names[!is.na(is_names) & nzchar(is_names)]

  pipeline_log <- list(
    is_match = is_match,
    normalize_method = normalize_method,
    scale_factor = scale_factor
  )

  if (inherits(object, "NTSProject")) {
    object <- normalize_by_IS(
      object,
      is_names = is_names,
      abundance_cols = abundance_cols,
      scale_factor = scale_factor,
      method = normalize_method,
      is_match = is_match,
      rt_breaks = rt_breaks
    )
    if (isTRUE(predict_rie)) {
      object <- predict_rie(object, default_rie = default_rie,
                            logIE_model = logIE_model,
                            polarity = polarity,
                            smiles_col = smiles_col)
      pipeline_log$predict_rie <- list(
        default_rie = default_rie,
        logIE_model = !is.null(logIE_model),
        polarity = polarity
      )
    }
    if (isTRUE(calculate_conc)) {
      object <- calculate_concentration(
        object,
        is_amount = is_amount,
        is_intensity = is_intensity,
        is_match = is_match
      )
      pipeline_log$calculate_concentration <- TRUE
    }
    methods::slot(object, "log") <- .nts_log_append(
      methods::slot(object, "log"),
      "run_semi_quantification",
      pipeline_log
    )
    return(object)
  }

  if (!is.data.frame(object) && !data.table::is.data.table(object)) {
    stop("`object` must be NTSProject or a feature table.", call. = FALSE)
  }
  ft <- as.data.frame(object, stringsAsFactors = FALSE)
  ft <- .normalize_by_IS_df(
    ft,
    is_names = is_names,
    abundance_cols = abundance_cols,
    scale_factor = scale_factor,
    method = normalize_method,
    is_match = is_match,
    rt_breaks = rt_breaks
  )
  if (isTRUE(predict_rie)) {
    ft <- .predict_rie_df(ft, default_rie = default_rie,
                           logIE_model = logIE_model,
                           polarity = polarity,
                           smiles_col = smiles_col)
    pipeline_log$predict_rie <- list(
      default_rie = default_rie,
      logIE_model = !is.null(logIE_model),
      polarity = polarity
    )
  }
  if (isTRUE(calculate_conc)) {
    ft <- .calculate_concentration_df(
      ft,
      is_amount = is_amount,
      is_intensity = is_intensity,
      is_match = is_match
    )
    pipeline_log$calculate_concentration <- TRUE
  }
  list(feature_table = ft, is_match = is_match, log = pipeline_log)
}
