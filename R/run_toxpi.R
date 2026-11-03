# toxpiR integration for prepare_toxpi() output
# Ref: tests/toxPi_202506065.R

#' @keywords internal
#' @noRd
.nts_toxpi_require_toxpiR <- function(fn_label = "run_toxpi") {
  if (!requireNamespace("toxpiR", quietly = TRUE)) {
    stop(
      fn_label, ": requires package 'toxpiR'. ",
      "Install with install.packages('toxpiR').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.nts_toxpi_identity_trans_list <- function(n = 1L) {
  funcs <- stats::setNames(
    rep(list(function(x) x), n),
    paste0("f", seq_len(n))
  )
  do.call(toxpiR::TxpTransFuncList, funcs)
}

#' @keywords internal
#' @noRd
.nts_toxpi_build_model <- function(slice_names, weights = NULL) {
  .nts_toxpi_require_toxpiR()
  if (!length(slice_names)) {
    stop("No ToxPi slices to model.", call. = FALSE)
  }
  slices <- lapply(slice_names, function(nm) {
    toxpiR::TxpSlice(
      nm,
      txpTransFuncs = .nts_toxpi_identity_trans_list(1L)
    )
  })
  names(slices) <- slice_names
  txp_slices <- do.call(toxpiR::TxpSliceList, slices)
  w <- if (is.null(weights)) {
    rep(1, length(slice_names))
  } else {
    as.numeric(weights[slice_names])
  }
  if (any(!is.finite(w)) || any(w < 0)) {
    stop("Invalid slice weights for toxpiR model.", call. = FALSE)
  }
  toxpiR::TxpModel(
    txpSlices = txp_slices,
    txpWeights = w
  )
}

#' @keywords internal
#' @noRd
.nts_toxpi_validate_prep <- function(toxpi_prep) {
  if (!is.list(toxpi_prep) || !is.data.frame(toxpi_prep$data)) {
    stop(
      "`toxpi_prep` must be the list returned by prepare_toxpi() ",
      "(with $data, $slices, $weights).",
      call. = FALSE
    )
  }
  if (is.null(toxpi_prep$weights) || !length(toxpi_prep$weights)) {
    stop("`toxpi_prep$weights` is missing or empty.", call. = FALSE)
  }
  invisible(toxpi_prep)
}

#' Calculate ToxPi scores with toxpiR from prepare_toxpi() output
#'
#' Builds a [toxpiR::TxpModel] from slice-level scores in [prepare_toxpi()] output
#' (one pre-aggregated metric per slice, 0-1), then runs
#' [toxpiR::txpCalculateScores]. Diagnostic columns such as `Ready_biodeg_risk`
#' are not included in the model unless listed in `extra_slices`.
#'
#' @param toxpi_prep List from [prepare_toxpi()] (`$data`, `$slices`, `$weights`).
#' @param id_var Column in `$data` used as compound ID (default `"Name"`).
#' @param slice_names Character vector of slice columns to model; default
#'   `names(toxpi_prep$weights)`.
#' @param extra_slices Optional extra columns from `$data` to add (e.g.
#'   `"Ready_biodeg_risk"`).
#' @param weight_multiplier Scale `toxpi_prep$weights` for [toxpiR::TxpModel]
#'   (default `100` for integer-like relative weights; set `1` to use proportions).
#' @param negative_value_handling Passed to [toxpiR::txpCalculateScores]:
#'   `"keep"` or `"missing"` (GUI-like).
#' @param return_model If `TRUE`, return `list(results, model, input)`.
#'
#' @return A **toxpiR** `TxpResult` object (overall + slice scores, ranks). If
#'   `return_model = TRUE`, a list with `results`, `model`, and `input`.
#'
#' @seealso [prepare_toxpi], [plot_toxpi], [toxpiR::txpCalculateScores]
#' @export
#'
#' @examplesIf requireNamespace("toxpiR", quietly = TRUE)
#' \dontrun{
#' smiles <- c("CCO", "CC", "CCC")
#' prep <- prepare_toxpi(
#'   predict_episuite(smiles),
#'   predict_admet_risk(smiles, backend = "admet_ai"),
#'   query_smiles = smiles
#' )
#' res <- run_toxpi(prep)
#' plot_toxpi(res)
#' }
run_toxpi <- function(toxpi_prep,
                      id_var = "Name",
                      slice_names = NULL,
                      extra_slices = NULL,
                      weight_multiplier = 100,
                      negative_value_handling = c("keep", "missing"),
                      return_model = FALSE) {
  .nts_toxpi_require_toxpiR("run_toxpi")
  negative_value_handling <- match.arg(negative_value_handling)
  toxpi_prep <- .nts_toxpi_validate_prep(toxpi_prep)

  if (is.null(slice_names)) {
    slice_names <- names(toxpi_prep$weights)
  }
  if (!is.null(extra_slices) && length(extra_slices)) {
    slice_names <- unique(c(slice_names, extra_slices))
  }
  missing_cols <- setdiff(slice_names, names(toxpi_prep$data))
  if (length(missing_cols)) {
    stop(
      "Slice column(s) not found in toxpi_prep$data: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  if (!id_var %in% names(toxpi_prep$data)) {
    stop("id_var '", id_var, "' not found in toxpi_prep$data.", call. = FALSE)
  }

  w <- toxpi_prep$weights[intersect(slice_names, names(toxpi_prep$weights))]
  if (length(extra_slices)) {
    extra_w <- rep(
      mean(toxpi_prep$weights, na.rm = TRUE),
      length(extra_slices)
    )
    names(extra_w) <- extra_slices
    w <- c(w, extra_w[setdiff(extra_slices, names(w))])
  }
  w <- w[slice_names]
  if (is.finite(weight_multiplier) && weight_multiplier != 1) {
    w <- w * weight_multiplier
  }

  model <- .nts_toxpi_build_model(slice_names, weights = w)
  input <- toxpi_prep$data[, c(id_var, slice_names), drop = FALSE]

  results <- toxpiR::txpCalculateScores(
    model = model,
    input = input,
    id.var = id_var,
    negative.value.handling = negative_value_handling
  )

  if (isTRUE(return_model)) {
    return(list(results = results, model = model, input = input))
  }
  results
}

#' Plot toxpiR results (radial ToxPi charts)
#'
#' Thin wrapper around [toxpiR::plot()]. For ranked overviews, pass
#' `sort(results)` or subset with `[`.
#'
#' @param toxpi_results Object from [run_toxpi()] (toxpiR `TxpResult`).
#' @param package Graphics backend: `"gg"` (ggplot2) or `"base"`.
#' @param ... Passed to [toxpiR::plot()].
#'
#' @return Invisibly returns the plot result from toxpiR.
#' @export
#'
#' @seealso [run_toxpi], [prepare_toxpi]
plot_toxpi <- function(toxpi_results, package = c("gg", "base"), ...) {
  .nts_toxpi_require_toxpiR("plot_toxpi")
  package <- match.arg(package)
  toxpiR::plot(toxpi_results, package = package, ...)
}

#' Extract slice score matrix from a toxpiR result
#'
#' @param toxpi_results Object from [run_toxpi()].
#' @param adjusted If `TRUE` (default), weights are applied per toxpiR rules.
#'
#' @return Numeric matrix (compounds x slices).
#' @export
toxpi_slice_scores <- function(toxpi_results, adjusted = TRUE) {
  .nts_toxpi_require_toxpiR("toxpi_slice_scores")
  toxpiR::txpSliceScores(toxpi_results, adjusted = adjusted)
}

#' Extract overall ToxPi scores from a toxpiR result
#'
#' @param toxpi_results Object from [run_toxpi()].
#'
#' @return Named numeric vector of overall ToxPi scores (higher = higher priority).
#' @export
toxpi_scores <- function(toxpi_results) {
  .nts_toxpi_require_toxpiR("toxpi_scores")
  sc <- toxpiR::txpScores(toxpi_results)
  if (is.null(names(sc)) || !any(nzchar(names(sc)))) {
    ids <- methods::slot(toxpi_results, "txpIDs")
    if (length(ids) == length(sc)) {
      names(sc) <- as.character(ids)
    }
  }
  sc
}
