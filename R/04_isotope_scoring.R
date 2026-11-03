#' Score measured isotope pattern against theoretical formula
#'
#' For each parent feature in a project and each candidate formula, generates
#' the theoretical isotope envelope via Rdisop and compares it to measured
#' isotopologue intensities using cosine similarity and chi-squared goodness-of-fit.
#'
#' @param project An NTSProject with @feature_table (sample intensities) and
#'   @ms1_annotations containing `flag_isotopes` output (`parent_name`, `isotope_c13_n`,
#'   `isotope_kind`).
#' @param formulas A data.frame or character vector of candidate formulas.
#'   If data.frame, must contain `feature_name` and `formula` columns. If character,
#'   a single formula applied to all parent features.
#' @param sample_cols Character vector of sample column names for intensity
#'   extraction. If NULL, uses all numeric sample columns from @feature_table.
#' @param mz_tol_ppm m/z tolerance for matching theoretical peaks to measured
#'   isotopologues. Default 5.
#' @param min_isotopologues Minimum number of measured isotopologues required
#'   (including M+0) to attempt fitting. Default 1.
#' @param write_to_project If TRUE, writes `isotope_cosine`, `isotope_chisq`,
#'   and `isotope_formula` to @ms1_annotations. Default FALSE.
#' @param ... Reserved for future parameters.
#'
#' @return When `write_to_project = FALSE` (default), a data.table with columns:
#'   `parent_name`, `formula`, `isotope_cosine` (cosine similarity, 0-1),
#'   `isotope_chisq` (reduced chi-squared), `n_measured` (number of measured
#'   isotopologues), `n_theoretical` (number of theoretical peaks included).
#'   When `write_to_project = TRUE`, returns the updated [NTSProject] invisibly;
#'   assign back: `project <- score_isotope_fit(project, ...)`. The score table
#'   is written to `@ms1_annotations` columns `isotope_formula`, `isotope_cosine`,
#'   `isotope_chisq`.
#' @seealso [flag_isotopes], [extract_isotope_pattern_table], [calculate_exact_mass]
#' @export
score_isotope_fit <- function(project,
                               formulas,
                               sample_cols = NULL,
                               mz_tol_ppm = 5,
                               min_isotopologues = 1L,
                               write_to_project = FALSE,
                               ...) {
  .nts_check(project, "project", .arg = "project")
  if (!requireNamespace("Rdisop", quietly = TRUE)) {
    stop(
      "Rdisop is required: install.packages('BiocManager'); BiocManager::install('Rdisop')",
      call. = FALSE
    )
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }

  ft <- get_feature_table(project)
  if (!nrow(ft)) {
    stop("feature_table is empty.", call. = FALSE)
  }
  if (!"name" %in% names(ft)) {
    stop("feature_table must contain column `name`.", call. = FALSE)
  }
  if (!"mz" %in% names(ft)) {
    stop("feature_table must contain column `mz`.", call. = FALSE)
  }

  ann <- get_ms1_annotations(project)
  if (!nrow(ann)) {
    warning("ms1_annotations is empty; no isotope child data available.",
            call. = FALSE)
    return(.nts_isotope_fit_empty())
  }

  # Resolve sample columns for intensity extraction
  if (is.null(sample_cols)) {
    sample_cols <- .nts_isotope_fit_sample_cols(ft)
  } else {
    sample_cols <- intersect(as.character(sample_cols), names(ft))
  }
  if (!length(sample_cols)) {
    stop("No numeric sample columns found in feature_table.", call. = FALSE)
  }

  # Per-feature mean intensity across samples
  int_vec <- rowMeans(
    as.matrix(ft[, sample_cols, drop = FALSE]),
    na.rm = TRUE
  )
  names(int_vec) <- trimws(as.character(ft$name))

  # Parse formulas input
  formula_map <- .nts_isotope_fit_parse_formulas(formulas, ann)

  # Identify parent features (those with at least one isotope child)
  iso_meta <- .nts_isotope_fit_iso_meta(ft, ann)
  parents <- unique(iso_meta$parent_name[iso_meta$is_child])
  parents <- parents[!is.na(parents) & nzchar(parents)]

  if (!length(parents) && min_isotopologues > 1L) {
    warning("No parent features with flagged isotope children found.",
            call. = FALSE)
    return(.nts_isotope_fit_empty())
  }

  # For min_isotopologues == 1, allow parents without children (M+0 only)
  if (!length(parents)) {
    parents <- character(0)
  }
  all_parents <- unique(c(
    parents,
    if (min_isotopologues <= 1L) trimws(as.character(ann$name))
  ))
  all_parents <- all_parents[!is.na(all_parents) & nzchar(all_parents)]

  out_list <- vector("list", length(all_parents))
  out_idx <- 0L

  for (i in seq_along(all_parents)) {
    pn <- all_parents[i]
    pidx <- match(pn, trimws(as.character(ft$name)))
    if (is.na(pidx)) next

    # Build measured envelope: M+0 from parent, then child peaks sorted by m/z
    env <- .nts_isotope_fit_build_envelope(
      pn, ft, iso_meta, int_vec
    )
    if (length(env$mz) < min_isotopologues) next

    # Get candidate formulas for this feature
    fm_vec <- if (is.character(formula_map)) {
      formula_map
    } else if (is.data.frame(formula_map)) {
      fm <- formula_map$formula[formula_map$feature_name == pn]
      if (!length(fm)) next
      fm
    } else {
      character(0)
    }

    for (fm in fm_vec) {
      if (is.na(fm) || !nzchar(trimws(fm))) next
      fm <- trimws(fm)
      if (nchar(fm, type = "bytes", allowNA = TRUE) > 512L) next

      sc <- .nts_score_isotope_fit_single(
        mz_vec = env$mz,
        intensity_vec = env$intensity,
        formula = fm,
        mz_tol_ppm = mz_tol_ppm
      )

      out_idx <- out_idx + 1L
      out_list[[out_idx]] <- data.table::data.table(
        parent_name = pn,
        formula = fm,
        isotope_cosine = sc$cosine,
        isotope_chisq = sc$chisq,
        n_measured = length(env$mz),
        n_theoretical = sc$n_theoretical
      )
    }
  }

  if (out_idx == 0L) {
    return(.nts_isotope_fit_empty())
  }

  out <- data.table::rbindlist(out_list[seq_len(out_idx)])

  if (isTRUE(write_to_project) && nrow(out)) {
    project <- .nts_isotope_fit_write_to_project(project, out)
    return(invisible(project))
  }

  out[]
}


#' Score single isotope pattern against one formula
#'
#' @param mz_vec Numeric vector of measured isotopologue m/z values (length >= 1).
#' @param intensity_vec Numeric vector of measured intensities, same length as mz_vec.
#'   M+0 should be the first element. Values are normalized internally.
#' @param formula A single molecular formula string (e.g. "C8HF15O2").
#' @param mz_tol_ppm m/z tolerance for matching theoretical to measured peaks.
#' @param charge Charge state used to rescale the theoretical isotope spacing
#'   (peaks sit ~1/|charge| apart in m/z). Default 1, which assumes the measured
#'   m/z are neutral monoisotopic masses and applies no rescaling.
#'
#' @return A named list with `cosine` (0-1), `chisq` (reduced chi-squared),
#'   `n_matched` (number of theoretical peaks matched), `n_theoretical` (peaks in
#'   theoretical envelope up to the last matched measured peak).
#' @keywords internal
#' @noRd
.nts_score_isotope_fit_single <- function(mz_vec,
                                           intensity_vec,
                                           formula,
                                           mz_tol_ppm = 5,
                                           charge = 1) {
  if (missing(mz_vec) || !length(mz_vec)) {
    stop("mz_vec must be a non-empty numeric vector.", call. = FALSE)
  }
  if (missing(intensity_vec) || length(intensity_vec) != length(mz_vec)) {
    stop("intensity_vec must have the same length as mz_vec.", call. = FALSE)
  }
  if (missing(formula) || !is.character(formula) || length(formula) != 1L ||
      is.na(formula) || !nzchar(trimws(formula))) {
    stop("formula must be a single non-empty string.", call. = FALSE)
  }

  mz_vec <- suppressWarnings(as.numeric(mz_vec))
  intensity_vec <- suppressWarnings(as.numeric(intensity_vec))

  ok <- is.finite(mz_vec) & is.finite(intensity_vec) & intensity_vec > 0
  if (!any(ok)) {
    return(list(
      cosine = NA_real_, chisq = NA_real_,
      n_matched = 0L, n_theoretical = 0L
    ))
  }
  mz_vec <- mz_vec[ok]
  intensity_vec <- intensity_vec[ok]

  # Single peak (M+0 only): perfect match by definition
  if (length(mz_vec) == 1L) {
    return(list(
      cosine = 1, chisq = 0,
      n_matched = 1L, n_theoretical = 1L
    ))
  }

  if (!requireNamespace("Rdisop", quietly = TRUE)) {
    stop(
      "Rdisop is required: install.packages('BiocManager'); BiocManager::install('Rdisop')",
      call. = FALSE
    )
  }

  f <- trimws(formula)
  mol <- tryCatch(
    Rdisop::getMolecule(f),
    error = function(e) NULL
  )
  if (is.null(mol) || is.null(mol$isotopes) || !length(mol$isotopes)) {
    return(list(
      cosine = NA_real_, chisq = NA_real_,
      n_matched = 0L, n_theoretical = 0L
    ))
  }

  # Rdisop isotopes: a list where [[1]] is a 2-row matrix
  # Row 1 = masses, Row 2 = abundances (max = 1 for monoisotopic)
  iso_mat <- mol$isotopes[[1]]
  if (!is.matrix(iso_mat) || nrow(iso_mat) < 2L || !ncol(iso_mat)) {
    return(list(
      cosine = NA_real_, chisq = NA_real_,
      n_matched = 0L, n_theoretical = 0L
    ))
  }

  theo_mz <- as.numeric(iso_mat[1L, ])
  theo_ab <- as.numeric(iso_mat[2L, ])

  # Rdisop reports NEUTRAL monoisotopic masses. Isotope peaks are spaced ~1 Da
  # apart in mass but ~1/|z| apart in m/z, so rescale the spacing for multiply
  # charged ions while keeping the M+0 peak anchored. charge = 1 (the default) is
  # a no-op, so the existing neutral-mass behaviour is unchanged. The adduct
  # offset itself is NOT modelled here - a caller that passes an adduct m/z
  # (e.g. [M-H]-) must first align the M+0 peak.
  .charge <- suppressWarnings(abs(as.numeric(charge)[1L]))
  if (is.finite(.charge) && .charge > 0 && .charge != 1 && length(theo_mz) > 1L) {
    theo_mz <- theo_mz[1L] + (theo_mz - theo_mz[1L]) / .charge
  }

  # Normalize measured intensities to max = 1
  meas_norm <- intensity_vec / max(intensity_vec, na.rm = TRUE)

  # Match each theoretical peak to the closest measured peak within ppm.
  # Explicit pair vectors ensure cosine and chi-squared are computed on
  # correctly aligned (theo, meas) pairs - not on independent boolean subsets.
  n_theo <- length(theo_mz)
  pair_theo <- numeric(0)
  pair_meas <- numeric(0)
  used_meas <- rep(FALSE, length(mz_vec))

  for (ti in seq_len(n_theo)) {
    ppm_diffs <- abs(theo_mz[ti] - mz_vec) / theo_mz[ti] * 1e6
    best <- which.min(ppm_diffs)
    if (length(best) && ppm_diffs[best] <= mz_tol_ppm) {
      pair_theo <- c(pair_theo, theo_ab[ti])
      pair_meas <- c(pair_meas, meas_norm[best])
      used_meas[best] <- TRUE
    }
  }

  n_matched <- length(pair_theo)
  if (n_matched < 1L) {
    return(list(
      cosine = 0, chisq = NA_real_,
      n_matched = 0L, n_theoretical = n_theo
    ))
  }

  # Cosine similarity on paired vectors (same length, aligned)
  cosine <- sum(pair_theo * pair_meas) /
    (sqrt(sum(pair_theo^2)) * sqrt(sum(pair_meas^2)))
  cosine <- pmax(0, pmin(1, cosine))

  # Reduced chi-squared: sum((obs - exp)^2 / exp) / (n_matched - 1)
  if (n_matched > 1L) {
    chisq <- sum((pair_meas - pair_theo)^2 / pair_theo) / (n_matched - 1L)
  } else {
    chisq <- 0
  }

  list(
    cosine = cosine,
    chisq = chisq,
    n_matched = n_matched,
    n_theoretical = n_theo
  )
}


# ============================================================================
# Internal helpers for score_isotope_fit
# ============================================================================


#' @keywords internal
#' @noRd
.nts_isotope_fit_empty <- function() {
  data.table::data.table(
    parent_name = character(0),
    formula = character(0),
    isotope_cosine = numeric(0),
    isotope_chisq = numeric(0),
    n_measured = integer(0),
    n_theoretical = integer(0)
  )
}


#' Extract mean intensity per feature (for isotope envelope building)
#' @keywords internal
#' @noRd
.nts_isotope_mean_intensity <- function(project) {
  ft <- get_feature_table(project)
  sample_cols <- setdiff(names(ft), c(.nts_feature_shape_cols, "name",
    "sample", "into", "feature_id", "nominal_mass", "kmd"))
  sample_cols <- sample_cols[vapply(ft[, sample_cols, drop = FALSE], is.numeric, logical(1L))]
  int_vec <- rowMeans(as.matrix(ft[, sample_cols, drop = FALSE]), na.rm = TRUE)
  names(int_vec) <- trimws(as.character(ft$name))
  int_vec
}


#' Resolve numeric sample columns from feature_table
#' @keywords internal
#' @noRd
.nts_isotope_fit_sample_cols <- function(ft) {
  cand <- setdiff(names(ft), c(.nts_feature_shape_cols, "name",
    "sample", "into", "feature_id", "nominal_mass", "kmd"))
  cand[vapply(ft[, cand, drop = FALSE], is.numeric, logical(1L))]
}


#' Parse formulas parameter into a lookup structure
#' @keywords internal
#' @noRd
.nts_isotope_fit_parse_formulas <- function(formulas, ann) {
  if (is.character(formulas)) {
    return(trimws(formulas))
  }
  if (is.data.frame(formulas)) {
    if (!"feature_name" %in% names(formulas)) {
      stop("formulas data.frame must contain column `feature_name`.", call. = FALSE)
    }
    if (!"formula" %in% names(formulas)) {
      stop("formulas data.frame must contain column `formula`.", call. = FALSE)
    }
    return(formulas)
  }
  stop("formulas must be a character vector or data.frame.", call. = FALSE)
}


#' Build isotope metadata lookup (is_child, parent_name) from ft + ann
#' @keywords internal
#' @noRd
.nts_isotope_fit_iso_meta <- function(ft, ann) {
  n <- nrow(ft)
  out <- data.frame(
    name    = trimws(as.character(ft$name)),
    mz      = as.numeric(ft$mz),
    is_child   = rep(FALSE, n),
    parent_name = rep(NA_character_, n),
    isotope_c13_n = rep(NA_integer_, n),
    isotope_kind  = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )

  if (!nrow(ann) || !"name" %in% names(ann)) {
    return(out)
  }

  # Deduplicate annotations by name
  ann_u <- ann[!duplicated(as.character(ann$name)), , drop = FALSE]
  idx <- match(out$name, trimws(as.character(ann_u$name)))
  ok <- is.finite(idx)

  if ("is_isotope" %in% names(ann_u)) {
    out$is_child[ok] <- as.logical(ann_u$is_isotope[idx[ok]])
  }
  out$is_child[is.na(out$is_child)] <- FALSE

  if ("parent_name" %in% names(ann_u)) {
    out$parent_name[ok] <- trimws(as.character(ann_u$parent_name[idx[ok]]))
  }

  if ("isotope_c13_n" %in% names(ann_u)) {
    out$isotope_c13_n[ok] <- suppressWarnings(as.integer(ann_u$isotope_c13_n[idx[ok]]))
  }

  if ("isotope_kind" %in% names(ann_u)) {
    out$isotope_kind[ok] <- trimws(as.character(ann_u$isotope_kind[idx[ok]]))
  }

  out
}


#' Build measured isotope envelope for a parent feature
#'
#' Returns m/z and intensity vectors: M+0 first, then child peaks sorted by m/z.
#' @keywords internal
#' @noRd
.nts_isotope_fit_build_envelope <- function(parent_name, ft, iso_meta, int_vec) {
  pid <- match(parent_name, iso_meta$name)
  if (is.na(pid)) {
    return(list(mz = numeric(0), intensity = numeric(0)))
  }

  mz_vals <- iso_meta$mz[pid]
  int_vals <- int_vec[pid]

  # Find child peaks
  child_mask <- iso_meta$is_child &
    !is.na(iso_meta$parent_name) &
    iso_meta$parent_name == parent_name

  if (any(child_mask)) {
    child_mz <- iso_meta$mz[child_mask]
    child_int <- int_vec[child_mask]
    child_ord <- order(child_mz)
    mz_vals <- c(mz_vals, child_mz[child_ord])
    int_vals <- c(int_vals, child_int[child_ord])
  }

  ok <- is.finite(mz_vals) & is.finite(int_vals) & int_vals > 0
  if (!any(ok)) {
    return(list(mz = numeric(0), intensity = numeric(0)))
  }

  list(mz = mz_vals[ok], intensity = int_vals[ok])
}


#' Write isotope scoring results back to @ms1_annotations
#' @keywords internal
#' @noRd
.nts_isotope_fit_write_to_project <- function(project, out) {
  ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
  if (!nrow(ms1) || !"name" %in% names(ms1)) {
    warning("ms1_annotations empty; skipping write_to_project.", call. = FALSE)
    return(invisible(project))
  }

  # For each scored parent, keep the best (highest cosine) formula
  best <- data.table::as.data.table(out)
  if ("isotope_cosine" %in% names(best)) {
    data.table::setorder(best, -isotope_cosine, na.last = TRUE)
  }
  best <- best[!duplicated(best$parent_name)]

  ms1$isotope_formula <- NA_character_
  ms1$isotope_cosine  <- NA_real_
  ms1$isotope_chisq   <- NA_real_

  idx <- match(best$parent_name, trimws(as.character(ms1$name)))
  valid <- !is.na(idx)
  if (any(valid)) {
    ms1$isotope_formula[idx[valid]] <- as.character(best$formula[valid])
    ms1$isotope_cosine[idx[valid]]  <- as.numeric(best$isotope_cosine[valid])
    ms1$isotope_chisq[idx[valid]]   <- as.numeric(best$isotope_chisq[valid])
  }

  project <- set_ms1_annotations(project, .nts_ms1_sanitize_and_order(ms1))
  invisible(project)
}
