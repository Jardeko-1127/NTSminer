#' Validate user diagnostic ion table (`mz`, `type`, optional `name`).
#'
#' @keywords internal
#' @noRd
.nts_validate_diagnostic_ions <- function(diagnostic_ions, arg = "diagnostic_ions") {
  if (is.null(diagnostic_ions)) {
    return(NULL)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  dt <- data.table::as.data.table(diagnostic_ions)
  if (!nrow(dt)) {
    return(NULL)
  }
  req <- c("mz", "type")
  miss <- setdiff(req, names(dt))
  if (length(miss)) {
    stop(
      "`", arg, "` must have columns: ",
      paste(req, collapse = ", "),
      call. = FALSE
    )
  }
  dt[, mz := as.numeric(mz)]
  types <- tolower(trimws(as.character(dt[["type"]])))
  ok_t <- types %in% c("fragment", "nl")
  if (!all(ok_t)) {
    stop(
      "`", arg, "$type` must be \"fragment\" or \"nl\" (case-insensitive)",
      call. = FALSE
    )
  }
  dt[, type := types]
  if (!"name" %in% names(dt)) {
    dt[, name := NA_character_]
  } else {
    dt[, name := trimws(as.character(name))]
  }
  dt <- dt[is.finite(mz)]
  if (!nrow(dt)) {
    return(NULL)
  }
  dt[]
}


#' Expected fragment m/z in the MS2 spectrum from a diagnostic ion table.
#'
#' For `type == "fragment"`, uses `mz` directly. For `type == "nl"`, uses
#' `precursor_mz - mz` (neutral-loss mass stored in `mz`).
#'
#' @keywords internal
#' @noRd
.nts_diagnostic_ions_expected_mz <- function(diagnostic_ions, precursor_mz) {
  dt <- .nts_validate_diagnostic_ions(diagnostic_ions, "diagnostic_ions")
  if (is.null(dt)) {
    return(numeric(0))
  }
  pm <- as.numeric(precursor_mz)[1L]
  out <- numeric(0)
  for (i in seq_len(nrow(dt))) {
    t <- dt$type[i]
    m <- dt$mz[i]
    if (!is.finite(m)) {
      next
    }
    if (identical(t, "fragment")) {
      out <- c(out, m)
    } else if (identical(t, "nl")) {
      if (is.finite(pm)) {
        out <- c(out, pm - m)
      }
    }
  }
  unique(out[is.finite(out) & out > 0])
}


#' Score overlap of experimental m/z with expected diagnostic m/z (Da window).
#'
#' @keywords internal
#' @noRd
.nts_score_exp_vs_diagnostic_mz <- function(exp_mz, diagnostic_ions, precursor_mz, tol_da) {
  if (is.null(diagnostic_ions)) {
    return(list(score = NA_real_, matched = ""))
  }
  dmz <- .nts_diagnostic_ions_expected_mz(diagnostic_ions, precursor_mz)
  if (!length(dmz)) {
    return(list(score = NA_real_, matched = ""))
  }
  exp_mz <- as.numeric(exp_mz)
  exp_mz <- exp_mz[is.finite(exp_mz)]
  if (!length(exp_mz)) {
    return(list(score = 0, matched = ""))
  }
  tol <- as.numeric(tol_da)[1L]
  if (!is.finite(tol) || tol < 0) {
    tol <- 0.02
  }
  hit <- vapply(dmz, function(target) any(abs(exp_mz - target) <= tol), logical(1L))
  score <- sum(hit) / length(dmz)
  matched <- if (any(hit)) {
    paste(format(sort(dmz[hit]), digits = 6L, trim = TRUE), collapse = "; ")
  } else {
    ""
  }
  list(score = score, matched = matched)
}

#' MS1 candidate prefilter via interval overlap (`foverlaps`)
#'
#' Symmetric m/z and retention-time windows are built from [ppm] and `rt_tol`
#' around each feature and reference row. The m/z overlap join uses
#' [data.table::foverlaps]; RT overlap is applied in one vectorized filter (no
#' `for` / nested `lapply` over candidates).
#'
#' @param features `data.table` (or coercible) with numeric columns `mz`, `rt`.
#' @param reference `data.table` (or coercible) with numeric columns `mz`, `rt`.
#' @param mz_tol_ppm Numeric scalar, m/z tolerance in ppm (symmetric window).
#' @param rt_tol Numeric scalar, RT tolerance (same units as `rt`, symmetric).
#' @param feature_id_col Name of optional ID column in `features`; if missing,
#'   `feature_id` is added as row index.
#' @param ref_id_col Name of optional ID column in `reference`; if missing,
#'   `ref_id` is added as row index.
#'
#' @return A `data.table` of candidate pairs with m/z/RT windows and original
#'   coordinates preserved for downstream Layer-3 scoring.
#' @export
ms1_prefilter_foverlaps <- function(features,
                                    reference,
                                    mz_tol_ppm,
                                    rt_tol,
                                    feature_id_col = "feature_id",
                                    ref_id_col = "ref_id") {
  if (!data.table::is.data.table(features)) {
    features <- data.table::as.data.table(features)
  }
  if (!data.table::is.data.table(reference)) {
    reference <- data.table::as.data.table(reference)
  }
  stopifnot(
    nrow(features) > 0L,
    nrow(reference) > 0L,
    "mz" %in% names(features),
    "rt" %in% names(features),
    "mz" %in% names(reference),
    "rt" %in% names(reference),
    is.numeric(mz_tol_ppm),
    length(mz_tol_ppm) == 1L,
    mz_tol_ppm >= 0,
    is.numeric(rt_tol),
    length(rt_tol) == 1L,
    rt_tol >= 0
  )

  fdt <- data.table::copy(features)
  rdt <- data.table::copy(reference)

  if (!feature_id_col %in% names(fdt)) {
    fdt[, (feature_id_col) := .I]
  }
  if (!ref_id_col %in% names(rdt)) {
    rdt[, (ref_id_col) := .I]
  }

  ppm <- mz_tol_ppm / 1e6

  fdt[, `:=`(
    mz_lo = mz * (1 - ppm),
    mz_hi = mz * (1 + ppm),
    f_rt_lo = data.table::fifelse(is.na(rt), -Inf, rt - rt_tol),
    f_rt_hi = data.table::fifelse(is.na(rt), Inf, rt + rt_tol),
    f_mz = mz,
    f_rt = rt
  )]

  rdt[, `:=`(
    mz_lo = mz * (1 - ppm),
    mz_hi = mz * (1 + ppm),
    r_rt_lo = data.table::fifelse(is.na(rt), -Inf, rt - rt_tol),
    r_rt_hi = data.table::fifelse(is.na(rt), Inf, rt + rt_tol),
    r_mz = mz,
    r_rt = rt
  )]

  f_ov <- fdt[, c(
    "mz_lo", "mz_hi", "f_rt_lo", "f_rt_hi", "f_mz", "f_rt", feature_id_col
  ), with = FALSE]
  r_ov <- rdt[, c(
    "mz_lo", "mz_hi", "r_rt_lo", "r_rt_hi", "r_mz", "r_rt", ref_id_col
  ), with = FALSE]

  data.table::setkey(r_ov, mz_lo, mz_hi)

  hits <- data.table::foverlaps(
    f_ov,
    r_ov,
    type = "any",
    nomatch = 0L
  )

  hits <- hits[
    f_rt_hi >= r_rt_lo & f_rt_lo <= r_rt_hi
  ]

  hits[]
}


#' Match m/z features with nearest-neighbour rolling join (internal)
#'
#' Uses `data.table` rolling joins with `roll = "nearest"` on sorted m/z keys (no `%in%`, no
#' `round()`-based binning). For **difference quantities** (neutral loss, mass
#' differences, diagnostic matching), tolerance is an **absolute Da** window:
#' `abs(mz - i.mz) <= feature_tol`.
#'
#' @param exp_features,lib_features `data.table` (or coercible) with column `mz`.
#' @param feature_tol Single numeric, absolute mass tolerance in **Da**.
#'
#' @return A list with `rate_exp`, `rate_lib` in `[0, 1]` and join tables
#'   `joined_exp`, `joined_lib`.
#' @noRd
.match_features_dt <- function(exp_features, lib_features, feature_tol) {
  exp_dt <- data.table::as.data.table(exp_features)
  lib_dt <- data.table::as.data.table(lib_features)
  stopifnot(
    "mz" %in% names(exp_dt),
    "mz" %in% names(lib_dt),
    length(feature_tol) == 1L,
    is.finite(feature_tol),
    feature_tol >= 0
  )
  exp_dt <- exp_dt[is.finite(mz)]
  lib_dt <- lib_dt[is.finite(mz)]
  if (nrow(exp_dt) == 0L || nrow(lib_dt) == 0L) {
    return(list(
      rate_exp = 0,
      rate_lib = 0,
      joined_exp = data.table::data.table(),
      joined_lib = data.table::data.table()
    ))
  }

  e <- exp_dt[, list(mz)]
  l <- lib_dt[, list(mz)]
  ft <- feature_tol

  data.table::setkey(l, mz)
  j1 <- l[e, .(mz_lib = x.mz, mz_exp = i.mz), roll = "nearest",
          allow.cartesian = TRUE]
  j1[, ok_exp := abs(mz_lib - mz_exp) <= ft]

  data.table::setkey(e, mz)
  data.table::setkey(l, mz)
  j2 <- e[l, .(mz_exp = x.mz, mz_lib = i.mz), roll = "nearest",
          allow.cartesian = TRUE]
  j2[, ok_lib := abs(mz_exp - mz_lib) <= ft]

  list(
    rate_exp = mean(j1$ok_exp, na.rm = TRUE),
    rate_lib = mean(j2$ok_lib, na.rm = TRUE),
    joined_exp = j1,
    joined_lib = j2
  )
}


# ============================================================================
# Spectrum pre-processing (clean, align, NDP)
# ============================================================================

#' Normalize and denoise a peak matrix
#' @param mat Two-column numeric matrix (mz, intensity).
#' @param int_thresh Relative intensity cutoff (0-1); peaks below are removed.
#' @return Cleaned two-column matrix, or NULL if nothing survives.
#' @noRd
.clean_spectrum <- function(mat, int_thresh = 0.025) {
  if (is.null(mat) || !is.matrix(mat) || nrow(mat) == 0L) return(NULL)
  mz  <- mat[, 1L]
  int <- mat[, 2L]
  ok <- is.finite(mz) & is.finite(int) & int > 0
  if (!any(ok)) return(NULL)
  mz  <- mz[ok]
  int <- int[ok]
  max_int <- max(int)
  if (max_int <= 0) return(NULL)
  int <- int / max_int
  keep <- int >= int_thresh
  if (!any(keep)) return(NULL)
  cbind(mz = mz[keep], intensity = int[keep])
}


#' Align two peak matrices by m/z within tolerance (outer join, 0-fill)
#'
#' Produces two equal-length intensity vectors suitable for ndotproduct.
#' Greedy nearest-neighbour 1:1 matching; unmatched peaks get partner = 0.
#' @return list(exp = matrix, lib = matrix) with identical nrow.
#' @noRd
.align_spectra <- function(exp_mat, lib_mat, mz_tol = 0.01) {
  e_mz  <- exp_mat[, 1L]
  e_int <- exp_mat[, 2L]
  l_mz  <- lib_mat[, 1L]
  l_int <- lib_mat[, 2L]

  ne <- length(e_mz)
  nl <- length(l_mz)

  e_matched <- rep(FALSE, ne)
  l_matched <- rep(FALSE, nl)
  pairs_e <- integer(0)
  pairs_l <- integer(0)

  if (ne > 0L && nl > 0L) {
    diffs <- outer(e_mz, l_mz, function(a, b) abs(a - b))
    while (TRUE) {
      min_val <- min(diffs)
      if (min_val > mz_tol) break
      idx <- which(diffs == min_val, arr.ind = TRUE)[1L, ]
      ei <- idx[1L]
      li <- idx[2L]
      pairs_e <- c(pairs_e, ei)
      pairs_l <- c(pairs_l, li)
      e_matched[ei] <- TRUE
      l_matched[li] <- TRUE
      diffs[ei, ] <- Inf
      diffs[, li] <- Inf
    }
  }

  n_paired <- length(pairs_e)
  n_e_only <- sum(!e_matched)
  n_l_only <- sum(!l_matched)
  total <- n_paired + n_e_only + n_l_only

  out_mz    <- numeric(total)
  out_e_int <- numeric(total)
  out_l_int <- numeric(total)

  if (n_paired > 0L) {
    out_mz[seq_len(n_paired)]    <- (e_mz[pairs_e] + l_mz[pairs_l]) / 2
    out_e_int[seq_len(n_paired)] <- e_int[pairs_e]
    out_l_int[seq_len(n_paired)] <- l_int[pairs_l]
  }
  off <- n_paired
  if (n_e_only > 0L) {
    idx_e <- which(!e_matched)
    rng <- off + seq_along(idx_e)
    out_mz[rng]    <- e_mz[idx_e]
    out_e_int[rng] <- e_int[idx_e]
  }
  off <- off + n_e_only
  if (n_l_only > 0L) {
    idx_l <- which(!l_matched)
    rng <- off + seq_along(idx_l)
    out_mz[rng]    <- l_mz[idx_l]
    out_l_int[rng] <- l_int[idx_l]
  }

  list(
    exp = cbind(mz = out_mz, intensity = out_e_int),
    lib = cbind(mz = out_mz, intensity = out_l_int),
    n_matched = n_paired
  )
}


#' Clean, align, and compute NDP between two spectra
#' @return list(score, n_matched, n_exp, n_lib)
#' @noRd
.calculate_clean_ndp <- function(exp_mat, lib_mat,
                                 mz_tol = 0.01,
                                 int_thresh = 0.025,
                                 m = 0L, n = 0.5) {
  e <- .clean_spectrum(exp_mat, int_thresh)
  l <- .clean_spectrum(lib_mat, int_thresh)
  if (is.null(e) || is.null(l)) {
    return(list(score = 0, n_matched = 0L,
                n_exp = if (is.null(e)) 0L else nrow(e),
                n_lib = if (is.null(l)) 0L else nrow(l)))
  }
  al <- .align_spectra(e, l, mz_tol)
  ndp_val <- MsCoreUtils::ndotproduct(al$exp, al$lib, m = m, n = n, na.rm = TRUE)
  ndp_val <- max(0, min(1, ndp_val))
  list(score = ndp_val, n_matched = al$n_matched,
       n_exp = nrow(e), n_lib = nrow(l))
}


# ============================================================================
# Multidimensional scoring
# ============================================================================

#' Multidimensional MS2 score (NDP + NL + mass-diff + diagnostics)
#'
#' Pre-processes both spectra (normalize, denoise, align) before computing
#' [MsCoreUtils::ndotproduct]. Then scores neutral losses, pairwise mass
#' differences, and diagnostic ions via `.match_features_dt`.
#'
#' All sub-scores are individually in `[0, 1]`, and the final `total` is a
#' convex combination (weights sum to 1), so `total` is also `[0, 1]`.
#'
#' @param exp_peaks,ref_peaks Two-column numeric matrices (m/z, intensity).
#' @param precursor_mz,precursor_mz_ref Precursor m/z values.
#' @param feature_tol Absolute Da tolerance for NL/MD/diagnostic matching.
#' @param mz_tol_align Da tolerance for peak alignment before NDP (default
#'   `0.01`).
#' @param int_thresh Relative intensity threshold for denoising (default
#'   `0.025`, i.e. 2.5% of base peak).
#' @param weights Named numeric `c(ndp, nl, md, diag)` giving the **absolute**
#'   MS2 sub-score weights. The default `c(0.20, 0.10, 0.10, 0.10)` yields an
#'   MS2 composite with theoretical maximum `0.50`.
#' @param diagnostic_mz_lib Optional known diagnostic m/z vector (experimental
#'   spectrum). Ignored when `diagnostic_ions` is non-`NULL`.
#' @param diagnostic_ions Optional `data.frame`/`data.table` with columns
#'   `mz`, `type` (`"fragment"` or `"nl"`), and optional `name`. When `NULL`,
#'   the diagnostic sub-score uses neither this table nor library top peaks
#'   unless `use_library_top_diagnostic_peaks = TRUE`.
#' @param use_library_top_diagnostic_peaks If `TRUE` and both `diagnostic_ions`
#'   and `diagnostic_mz_lib` are unset, use the top `top_diag` library peaks as
#'   diagnostic targets (legacy behaviour). Default `FALSE`.
#' @param top_diag See `use_library_top_diagnostic_peaks`.
#' @param m,n Weighting exponents for ndotproduct (default 0, 0.5).
#'
#' @return A list containing raw component scores (`ndp`, `neutral_loss`,
#'   `mass_diff`, `diagnostic`), weighted component scores (`score_ndp`,
#'   `score_nl`, `score_md`, `score_diagnostic`), `score_ms2_composite`,
#'   `matched_nl`, `matched_md`, `matched_diagnostic`, and `weights`.
#' @export
score_spectra_multidim <- function(exp_peaks,
                                   ref_peaks,
                                   precursor_mz,
                                   precursor_mz_ref = precursor_mz,
                                   feature_tol = 0.01,
                                   mz_tol_align = 0.01,
                                   int_thresh = 0.025,
                                   weights = c(ndp = 0.20, nl = 0.10, md = 0.10, diag = 0.10),
                                   diagnostic_mz_lib = NULL,
                                   diagnostic_ions = NULL,
                                   use_library_top_diagnostic_peaks = FALSE,
                                   top_diag = 5L,
                                   m = 0L,
                                   n = 0.5) {
  .validate_peak_matrix(exp_peaks, "exp_peaks")
  .validate_peak_matrix(ref_peaks, "ref_peaks")
  stopifnot(
    length(precursor_mz) == 1L, is.finite(precursor_mz),
    length(precursor_mz_ref) == 1L, is.finite(precursor_mz_ref)
  )

  ndp_res <- .calculate_clean_ndp(
    exp_peaks, ref_peaks,
    mz_tol = mz_tol_align, int_thresh = int_thresh,
    m = m, n = n
  )
  ndp_score <- ndp_res$score

  e_clean <- .clean_spectrum(exp_peaks, int_thresh)
  l_clean <- .clean_spectrum(ref_peaks, int_thresh)
  if (is.null(e_clean)) e_clean <- exp_peaks
  if (is.null(l_clean)) l_clean <- ref_peaks

  exp_m <- e_clean[, 1L]
  ref_m <- l_clean[, 1L]

  # --- Neutral Loss ---
  exp_nl_vals <- precursor_mz - exp_m
  lib_nl_vals <- precursor_mz_ref - ref_m
  exp_nl <- data.table::data.table(mz = exp_nl_vals[exp_nl_vals > 0])
  lib_nl <- data.table::data.table(mz = lib_nl_vals[lib_nl_vals > 0])
  nl_mt <- .match_features_dt(exp_nl, lib_nl, feature_tol)
  nl_score <- (nl_mt$rate_exp + nl_mt$rate_lib) / 2

  matched_nl_str <- ""
  if (nrow(nl_mt$joined_exp) > 0L) {
    j <- nl_mt$joined_exp[ok_exp == TRUE]
    if (nrow(j)) {
      matched_nl_str <- paste(
        sprintf("%.4f", sort(unique(j$mz_exp))),
        collapse = "; "
      )
    }
  }

  # --- Mass Difference (pairwise) ---
  if (length(exp_m) >= 2L && length(ref_m) >= 2L) {
    dm_e <- .pairwise_mz_dist(exp_m)
    dm_l <- .pairwise_mz_dist(ref_m)
    md_mt <- .match_features_dt(
      data.table::data.table(mz = dm_e),
      data.table::data.table(mz = dm_l),
      feature_tol
    )
    md_score <- (md_mt$rate_exp + md_mt$rate_lib) / 2
  } else {
    md_mt <- list(rate_exp = 0, rate_lib = 0,
                  joined_exp = data.table::data.table(),
                  joined_lib = data.table::data.table())
    md_score <- 0
  }
  matched_md_str <- ""
  if (nrow(md_mt$joined_exp) > 0L) {
    j <- md_mt$joined_exp[ok_exp == TRUE]
    if (nrow(j)) {
      matched_md_str <- paste(
        sprintf("%.4f", sort(unique(j$mz_exp))),
        collapse = "; "
      )
    }
  }

  # --- Diagnostic Ions ---
  k <- max(0L, min(as.integer(top_diag)[1L], nrow(l_clean)))
  d_ions <- .nts_validate_diagnostic_ions(diagnostic_ions, "diagnostic_ions")
  if (!is.null(d_ions)) {
    diagnostic_mz_lib <- .nts_diagnostic_ions_expected_mz(d_ions, precursor_mz)
  } else if (is.null(diagnostic_mz_lib)) {
    diagnostic_mz_lib <- numeric(0)
    if (isTRUE(use_library_top_diagnostic_peaks) && k >= 1L) {
      ord <- order(l_clean[, 2L], decreasing = TRUE)[seq_len(k)]
      diagnostic_mz_lib <- l_clean[ord, 1L]
    }
  }
  if (length(diagnostic_mz_lib)) {
    d_mt <- .match_features_dt(
      data.table::data.table(mz = exp_m),
      data.table::data.table(mz = diagnostic_mz_lib),
      feature_tol
    )
    diag_score <- (d_mt$rate_exp + d_mt$rate_lib) / 2
  } else {
    d_mt <- list(rate_exp = 0, rate_lib = 0,
                 joined_exp = data.table::data.table(),
                 joined_lib = data.table::data.table())
    diag_score <- 0
  }
  matched_diag_str <- ""
  if (nrow(d_mt$joined_lib) > 0L) {
    j <- d_mt$joined_lib[ok_lib == TRUE]
    if (nrow(j)) {
      matched_diag_str <- paste(
        sprintf("%.4f", sort(unique(j$mz_lib))),
        collapse = "; "
      )
    }
  }

  # --- Weighted MS2 composite (default max = 0.50) ---
  w <- weights[c("ndp", "nl", "md", "diag")]
  w[is.na(w)] <- 0
  if (sum(w) <= 0) stop("weights must sum to a positive value", call. = FALSE)
  if (sum(w) > 1) stop("weights must sum to <= 1", call. = FALSE)

  score_ndp <- as.numeric(w["ndp"] * ndp_score)
  score_nl <- as.numeric(w["nl"] * nl_score)
  score_md <- as.numeric(w["md"] * md_score)
  score_diagnostic <- as.numeric(w["diag"] * diag_score)
  score_ms2_composite <- score_ndp + score_nl + score_md + score_diagnostic
  score_ms2_composite <- max(0, min(1, score_ms2_composite))

  list(
    total = score_ms2_composite,
    score_ms2_composite = score_ms2_composite,
    ndp = ndp_score,
    neutral_loss = nl_score,
    mass_diff = md_score,
    diagnostic = diag_score,
    score_ndp = score_ndp,
    score_nl = score_nl,
    score_md = score_md,
    score_diagnostic = score_diagnostic,
    n_matched = ndp_res$n_matched,
    matched_nl = matched_nl_str,
    matched_md = matched_md_str,
    matched_diagnostic = matched_diag_str,
    weights = w
  )
}


#' MS2 spectral similarity (normalized dot product) with alignment
#'
#' Wrapper around [MsCoreUtils::ndotproduct] that first cleans and aligns the
#' two spectra. For library-scale workflows prefer [score_spectra_multidim].
#'
#' @param exp_peaks Numeric two-column matrix (m/z, intensity).
#' @param ref_peaks Numeric two-column matrix (m/z, intensity).
#' @param mz_tol_align Absolute Da tolerance for peak alignment (default 0.01).
#' @param int_thresh Relative intensity cutoff for denoising (default 0.025).
#' @param m,n Weighting exponents for ndotproduct (default 0, 0.5).
#'
#' @return Scalar in `[0, 1]`.
#' @export
score_ms2_ndotproduct <- function(exp_peaks,
                                  ref_peaks,
                                  mz_tol_align = 0.01,
                                  int_thresh = 0.025,
                                  m = 0L,
                                  n = 0.5) {
  .validate_peak_matrix(exp_peaks, "exp_peaks")
  .validate_peak_matrix(ref_peaks, "ref_peaks")
  res <- .calculate_clean_ndp(exp_peaks, ref_peaks,
                              mz_tol = mz_tol_align,
                              int_thresh = int_thresh,
                              m = m, n = n)
  res$score
}


#' MS2 scores for many spectrum pairs (single-level `Map`)
#'
#' `MsCoreUtils::ndotproduct` is defined for one pair of matrices. Parallel
#' lists of spectra are scored with one [Map] dispatch (not nested `lapply`).
#' For multidimensional scoring use [score_spectra_multidim] row-wise via `Map`
#' at the workflow level.
#'
#' @param exp_peaks_list `list` of two-column numeric matrices.
#' @param ref_peaks_list `list` of two-column numeric matrices, same length as
#'   `exp_peaks_list`.
#' @param ... Passed to each [score_ms2_ndotproduct] call.
#'
#' @return `numeric()` vector of length `length(exp_peaks_list)`.
#' @export
score_ms2_ndotproduct_map <- function(exp_peaks_list, ref_peaks_list, ...) {
  stopifnot(
    is.list(exp_peaks_list),
    is.list(ref_peaks_list),
    length(exp_peaks_list) == length(ref_peaks_list)
  )
  if (length(exp_peaks_list) == 0L) {
    return(numeric(0L))
  }
  unlist(
    Map(
      function(e, r) score_ms2_ndotproduct(e, r, ...),
      exp_peaks_list,
      ref_peaks_list
    ),
    use.names = FALSE
  )
}


#' Batch annotation with MS1 prefilter and MS2 multidimensional scoring
#'
#' For a [NTSProject] + [databaseClass]: runs [ms1_prefilter_foverlaps] on
#' `database@spectra.info` (with `RT` normalized to `rt`), then for each hit
#' loads experimental MS2 from `@ms2data` and library MS2 from
#' `@spectra.data$Spectra.negative`,
#' calls [score_spectra_multidim], drops `score_total <= 0`, and sorts by score.
#' For `data.frame` `database`, behaviour is MS1 overlap only (unchanged).
#'
#' @param x Feature `data.frame`/`data.table` or [NTSProject].
#' @param database `data.frame` with `mz` and `rt`, [databaseClass], or
#'   [NTSDatabase] (triggers [annotate_by_spectral_match]; optional
#'   `diagnostic_ions`, `ms2_source` and `ref_spectra` via `...` are forwarded
#'   to [annotate_by_spectral_match], defaults `"ms2data"` / `"per_ce"`).
#' @param min_score Minimum cosine score for a hit to pass. Defaults to `0.5`,
#'   matching the documented examples in `README` / the workflow vignette.
#' @param mz_tol,rt_tol MS1 tolerances (`mz_tol` in **ppm**).
#' @param ion_mode Which library slot to use: `Spectra.positive` or
#'   `Spectra.negative`.
#' @param polarity Ionisation polarity used to resolve the library precursor
#'   m/z. A `target_mz` column always wins; otherwise polarity-specific
#'   `mz.neg` / `mz.pos`, then `mz`, are used, with `mass \eqn{\mp} 1.007276`
#'   only as the last resort (`"neg"` = `[M-H]-`, `"pos"` = `[M+H]+`).
#' @param feature_tol Da tolerance passed to [score_spectra_multidim].
#' @param ... Additional arguments passed to [score_spectra_multidim]
#'   (and to [annotate_by_spectral_match] on the [NTSDatabase] path, including
#'   `validate_diagnostic_fragments` / `diagnostic_fragments_csv`, which default
#'   to `TRUE` / the packaged expert table).
#'
#' @return When `x` is an [NTSProject], always the updated `x` (same S4 object
#'   class). The MS1-only and library-scoring paths return invisibly. When `x`
#'   is a plain feature table, returns the overlap table from
#'   [ms1_prefilter_foverlaps] (typically a `data.table`).
#' @export
batch_annotate <- function(x,
                           database,
                           mz_tol = 30,
                           rt_tol = 150,
                           ion_mode = c("negative", "positive"),
                           polarity = c("neg", "pos"),
                           feature_tol = 0.01,
                           min_score = 0.5,
                           ...) {
  mz_tol_ppm <- mz_tol
  ion_mode <- match.arg(ion_mode)
  polarity <- match.arg(polarity)

  if (methods::is(x, "NTSProject") && .is_reference_spectral_db(database)) {
    dots <- list(...)
    ms <- min_score
    if (is.null(ms) || !is.finite(suppressWarnings(as.numeric(ms)[1L]))) {
      ms <- 0.5
    }
    ms <- suppressWarnings(as.numeric(ms)[1L])
    di <- dots$diagnostic_ions
    ms2_src <- dots$ms2_source
    if (is.null(ms2_src)) {
      ms2_src <- "ms2data"
    }
    ref_sp <- dots$ref_spectra
    if (is.null(ref_sp)) {
      ref_sp <- "per_ce"
    }
    vdf <- dots$validate_diagnostic_fragments
    if (is.null(vdf)) {
      vdf <- TRUE
    }
    dfc <- dots$diagnostic_fragments_csv
    dps <- dots$detect_pfas_series
    if (is.null(dps)) {
      dps <- TRUE
    }
    out <- annotate_by_spectral_match(
      x,
      database,
      mz_tol_ppm = mz_tol_ppm,
      rt_tol_s = rt_tol,
      ms2_tol_da = feature_tol,
      min_score = ms,
      polarity = polarity,
      diagnostic_ions = di,
      ms2_source = ms2_src,
      ref_spectra = ref_sp,
      validate_diagnostic_fragments = vdf,
      diagnostic_fragments_csv = dfc,
      detect_pfas_series = dps,
      pfas_series_min_len = dots$pfas_series_min_len %||% 3L,
      cf2_tol_da = dots$cf2_tol_da %||% 0.01,
      max_gap_units = dots$max_gap_units %||% 3L,
      min_rel_int = dots$min_rel_int %||% 0.01,
      series_n_max = dots$series_n_max %||% 10L,
      mode = dots$pfas_series_mode %||% "loose"
    )
    if (!methods::is(out, "NTSProject")) {
      stop(
        "batch_annotate(): internal error - annotate_by_spectral_match did not return an NTSProject.",
        call. = FALSE
      )
    }
    ann <- methods::slot(out, "db_annotations")
    if (is.data.frame(ann) && nrow(ann)) {
      n0 <- nrow(ann)
      ann <- .nts_dedupe_db_annotations_batch(ann)
      if (nrow(ann) < n0) {
        message(
          "batch_annotate: deduplicated db_annotations ",
          n0, " -> ", nrow(ann), " rows (by ref_id, best score_total)"
        )
      }
      methods::slot(out, "db_annotations") <- ann
      out <- .nts_merge_db_top1_into_ms2_annotations(out)
      out <- .nts_merge_db_best_into_ms1_library_seed_cols(out)
      ## Dedup above dropped feature x compound pairs: keep @ms2_score_grid's
      ## `is_winner` flags in step with the surviving @db_annotations rows.
      out <- .nts_sync_grid_winners(out)
    }
    return(out)
  }

  if (methods::is(x, "NTSProject") && inherits(database, "databaseClass")) {
    out <- .batch_annotate_project_dbclass(
      x,
      database,
      mz_tol_ppm = mz_tol_ppm,
      rt_tol = rt_tol,
      ion_mode = ion_mode,
      polarity = polarity,
      feature_tol = feature_tol,
      min_score = min_score,
      ...
    )
    if (!methods::is(out, "NTSProject")) {
      stop(
        "batch_annotate(): internal error - expected NTSProject from databaseClass path.",
        call. = FALSE
      )
    }
    return(out)
  }

  if (methods::is(x, "NTSProject")) {
    ft <- get_feature_table(x)
    if (!nrow(ft)) {
      stop("feature_table is empty", call. = FALSE)
    }
    ref <- if (is.data.frame(database)) {
      database
    } else if (inherits(database, "databaseClass")) {
      .nts_db_ms1_table(database, polarity = polarity)
    } else {
      stop(
        "For NTSProject, `database` must be a data.frame or databaseClass",
        call. = FALSE
      )
    }
    hits <- ms1_prefilter_foverlaps(ft, ref, mz_tol_ppm = mz_tol_ppm, rt_tol = rt_tol)
    ann <- if (data.table::is.data.table(hits)) {
      as.data.frame(hits, stringsAsFactors = FALSE)
    } else {
      hits
    }
    methods::slot(x, "db_annotations") <- ann
    x <- .nts_merge_db_best_into_ms1_library_seed_cols(x)
    methods::slot(x, "log") <- .nts_log_append(
      methods::slot(x, "log"),
      "batch_annotate",
      list(mz_tol = mz_tol, rt_tol = rt_tol, mode = "ms1_only")
    )
    return(invisible(x))
  }

  ref_tbl <- if (inherits(database, "databaseClass")) {
    .nts_db_ms1_table(database, polarity = polarity)
  } else {
    database
  }
  ms1_prefilter_foverlaps(x, ref_tbl, mz_tol_ppm = mz_tol_ppm, rt_tol = rt_tol)
}


#' Deduplicate `db_annotations` after [batch_annotate]: one row per `ref_id` (best
#' `score_total`), then sort by `ref_id` ascending and `score_total` descending.
#' Rows with missing `ref_id` are kept (not collapsed with each other).
#' @keywords internal
#' @noRd
.nts_dedupe_db_annotations_batch <- function(ann) {
  if (!is.data.frame(ann) || !nrow(ann)) {
    return(ann)
  }
  id_col <- if ("ref_id" %in% names(ann)) {
    "ref_id"
  } else if ("Lab_ID" %in% names(ann)) {
    "Lab_ID"
  } else {
    return(ann)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    u <- ann
    u[[id_col]] <- trimws(as.character(u[[id_col]]))
    has_id <- !is.na(u[[id_col]]) & nzchar(u[[id_col]])
    u_ok <- u[has_id, , drop = FALSE]
    u_na <- u[!has_id, , drop = FALSE]
    if (nrow(u_ok) && "score_total" %in% names(u_ok)) {
      u_ok <- u_ok[order(u_ok$score_total, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    }
    if (nrow(u_ok)) {
      u_ok <- u_ok[!duplicated(u_ok[[id_col]]), , drop = FALSE]
    }
    u <- rbind(u_ok, u_na)
    if (nrow(u) && id_col %in% names(u)) {
      if ("score_total" %in% names(u)) {
        u <- u[order(u[[id_col]], -u$score_total, na.last = c(TRUE, TRUE)), , drop = FALSE]
      } else {
        u <- u[order(u[[id_col]], na.last = TRUE), , drop = FALSE]
      }
    }
    rownames(u) <- NULL
    return(u)
  }
  dt <- data.table::as.data.table(ann)
  dt[, (id_col) := trimws(as.character(get(id_col)))]
  ok <- dt[!is.na(get(id_col)) & nzchar(get(id_col))]
  na_rows <- dt[is.na(get(id_col)) | !nzchar(get(id_col))]
  if (nrow(ok)) {
    if ("score_total" %in% names(ok)) {
      data.table::setorder(ok, -score_total, na.last = TRUE)
    }
    ok <- unique(ok, by = id_col)
    if ("score_total" %in% names(ok)) {
      data.table::setorderv(ok, c(id_col, "score_total"), c(1L, -1L), na.last = TRUE)
    } else {
      data.table::setorderv(ok, id_col, na.last = TRUE)
    }
  }
  out <- data.table::rbindlist(list(ok, na_rows), use.names = TRUE, fill = TRUE)
  as.data.frame(out, stringsAsFactors = FALSE)
}


#' @noRd
.batch_annotate_project_dbclass <- function(project,
                                            database,
                                            mz_tol_ppm,
                                            rt_tol,
                                            ion_mode,
                                            polarity,
                                            feature_tol,
                                            ...) {
  ft <- get_feature_table(project)
  comp <- .nts_db_ms1_table(database, polarity = polarity)
  cons <- methods::slot(project, "ms2data")
  fr_raw <- methods::slot(database, "fragments.data")
  fr <- .nts_resolve_fragments(fr_raw, polarity)

  if (!nrow(ft)) {
    stop("feature_table is empty", call. = FALSE)
  }
  if (!nrow(comp)) {
    stop("Database metadata table is empty.", call. = FALSE)
  }

  fname_col <- if ("name" %in% names(ft)) {
    "name"
  } else {
    "feature_id"
  }
  hits <- ms1_prefilter_foverlaps(
    ft,
    comp,
    mz_tol_ppm = mz_tol_ppm,
    rt_tol = rt_tol,
    feature_id_col = fname_col,
    ref_id_col = "Lab.ID"
  )
  if (!nrow(hits)) {
    methods::slot(project, "db_annotations") <- data.frame()
    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "batch_annotate",
      list(mz_tol = mz_tol_ppm, rt_tol = rt_tol, n_hits = 0L)
    )
    return(invisible(project))
  }

  sp_root <- methods::slot(database, "spectra.data")[[paste0("Spectra.", ion_mode)]]
  if (is.null(sp_root)) {
    stop("spectra.data$", paste0("Spectra.", ion_mode), " is missing", call. = FALSE)
  }

  id_hit_col <- if ("Lab.ID" %in% names(hits)) {
    "Lab.ID"
  } else if ("i.Lab.ID" %in% names(hits)) {
    "i.Lab.ID"
  } else {
    "ref_id"
  }
  fname_hit_col <- if ("name" %in% names(hits)) {
    "name"
  } else {
    fname_col
  }

  rows <- vector("list", nrow(hits))
  pb <- .nts_progress_start(nrow(hits), "batch_annotate")
  on.exit(.nts_progress_done(pb), add = TRUE)
  for (k in seq_len(nrow(hits))) {
    .nts_progress_update(pb, k)
    fname <- as.character(hits[[fname_hit_col]][k])
    rid <- as.character(hits[[id_hit_col]][k])
    exp_sp <- cons[[fname]]
    if (is.null(exp_sp)) {
      next
    }
    ref_levels <- sp_root[[rid]]
    if (is.null(ref_levels)) {
      ref_levels <- sp_root[[as.character(rid)]]
    }
    if (is.null(ref_levels)) {
      next
    }
    ref_mat <- .nts_lib_consensus_matrix(ref_levels)
    if (is.null(ref_mat) || !nrow(ref_mat)) {
      next
    }

    diag_mz <- .nts_get_diag_mz(fr, rid)

    prec_e <- hits[["f_mz"]][k]
    prec_r <- hits[["r_mz"]][k]
    if (!is.finite(prec_e)) {
      prec_e <- ft$mz[match(fname, as.character(ft[[fname_col]]))]
    }
    if (!is.finite(prec_r)) {
      prec_r <- comp$mz[match(rid, as.character(comp[["Lab.ID"]]))]
    }
    if (!is.finite(prec_r) && "exact_mass" %in% names(comp)) {
      em <- comp$exact_mass[match(rid, as.character(comp[["Lab.ID"]]))]
      prec_r <- .mass_to_mz(em, polarity)
    }

    exp_subset <- .nts_exp_ms2_to_peak_matrices(exp_sp)
    if (!length(exp_subset)) {
      next
    }

    sc <- .score_ms2_subset_max(
      exp_subset = exp_subset,
      ref_mat = ref_mat,
      precursor_mz = prec_e,
      precursor_mz_ref = prec_r,
      feature_tol = feature_tol,
      diagnostic_mz_lib = diag_mz,
      ...
    )
    if (is.null(sc)) {
      next
    }

    # --- 25/25/50: MS1 m/z score (Gaussian, max 0.25) ---
    exp_mz_k <- hits[["f_mz"]][k]
    lib_mz_k <- hits[["r_mz"]][k]
    mz_err_ppm <- if (is.finite(exp_mz_k) && is.finite(lib_mz_k) && lib_mz_k > 0) {
      abs(exp_mz_k - lib_mz_k) / lib_mz_k * 1e6
    } else {
      mz_tol_ppm
    }
    score_mz <- exp(-0.5 * (mz_err_ppm / mz_tol_ppm)^2) * 0.25

    # --- 25/25/50: RT score (Gaussian, max 0.25) ---
    exp_rt_k <- hits[["f_rt"]][k]
    lib_rt_k <- hits[["r_rt"]][k]
    rt_err <- if (is.finite(exp_rt_k) && is.finite(lib_rt_k)) {
      abs(exp_rt_k - lib_rt_k)
    } else {
      rt_tol * 3
    }
    score_rt <- exp(-0.5 * (rt_err / rt_tol)^2) * 0.25

    # --- 25/25/50: MS2 composite (max 0.50) ---
    score_ms2 <- sc$score_ms2_composite

    score_total_k <- score_mz + score_rt + score_ms2

    # --- Compound name from database ---
    cname <- ""
    c_smi <- NA_character_
    rid_match <- match(rid, as.character(comp[["Lab.ID"]]))
    if (!is.na(rid_match)) {
      for (cn in c("Compound.name", "compound_name", "Name", "name", "Annotation")) {
        if (cn %in% names(comp) && !is.na(comp[[cn]][rid_match])) {
          cname <- as.character(comp[[cn]][rid_match])
          break
        }
      }
      for (cn in c("SMILES", "smiles", "Smiles")) {
        if (cn %in% names(comp) && !is.na(comp[[cn]][rid_match])) {
          c_smi <- as.character(comp[[cn]][rid_match])
          break
        }
      }
    }

    # Database provenance (inhouse / online / prediction) for this ref_id, so
    # assign_confidence_levels() can tier L1 (in-house) vs L2a (online) vs L2b
    # (prediction) directly instead of inferring it from ID patterns. Sourced
    # from the (provenance-stamped) MS1 candidate table `comp`.
    dbt_k <- if (!is.na(rid_match) && "db_type" %in% names(comp)) {
      val <- as.character(comp[["db_type"]][rid_match])
      if (length(val) && !is.na(val[1L]) && nzchar(val[1L])) val[1L] else NA_character_
    } else {
      NA_character_
    }

    rows[[k]] <- data.frame(
      name = fname,
      ref_id = rid,
      db_type = dbt_k,
      Compound.name = cname,
      SMILES = c_smi,
      score_total = score_total_k,
      score_mz = round(score_mz, 4),
      score_rt = round(score_rt, 4),
      score_ms2 = round(score_ms2, 4),
      score_ndp = round(sc$score_ndp, 4),
      score_nl = round(sc$score_nl, 4),
      score_md = round(sc$score_md, 4),
      score_diagnostic = round(sc$score_diagnostic, 4),
      ndp_raw = round(sc$ndp, 4),
      matched_nl = sc$matched_nl,
      matched_md = sc$matched_md,
      matched_diagnostic = sc$matched_diagnostic,
      best_ms2_index = sc$best_ms2_index,
      n_ms2_subset = sc$n_ms2_subset,
      mz_error_ppm = round(mz_err_ppm, 2),
      rt_error_s = round(rt_err, 1),
      exp_mz = exp_mz_k,
      lib_mz = lib_mz_k,
      exp_rt = exp_rt_k,
      lib_rt = lib_rt_k,
      n_matched_peaks = sc$n_matched,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) {
    methods::slot(project, "db_annotations") <- data.frame()
    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "batch_annotate",
      list(mz_tol = mz_tol_ppm, rt_tol = rt_tol, n_hits = 0L, note = "no_ms2_scores")
    )
    return(invisible(project))
  }

  ann <- do.call(rbind, rows)
  ann <- .nts_dedupe_db_annotations_batch(ann)
  rownames(ann) <- NULL

  methods::slot(project, "db_annotations") <- ann
  project <- .nts_merge_db_top1_into_ms2_annotations(project)
  project <- .nts_merge_db_best_into_ms1_library_seed_cols(project)
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "batch_annotate",
    list(
      mz_tol = mz_tol_ppm,
      rt_tol = rt_tol,
      ion_mode = "negative",
      n_hits_in = nrow(hits),
      n_annotated = nrow(ann)
    )
  )
  return(invisible(project))
}


#' @noRd
.nts_lib_to_peak_matrix <- function(x) {
  if (is.matrix(x) && ncol(x) >= 2L) {
    return(cbind(as.numeric(x[, 1L]), as.numeric(x[, 2L])))
  }
  if (is.data.frame(x) && ncol(x) >= 2L) {
    return(cbind(as.numeric(x[[1L]]), as.numeric(x[[2L]])))
  }
  NULL
}


#' @noRd
.nts_lib_consensus_matrix <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(.nts_lib_to_peak_matrix(x))
  }
  if (!is.list(x) || !length(x)) {
    return(NULL)
  }
  mats <- lapply(x, .nts_lib_to_peak_matrix)
  mats <- mats[!vapply(mats, is.null, logical(1L))]
  if (!length(mats)) {
    return(NULL)
  }
  if (length(mats) == 1L) {
    return(mats[[1L]])
  }
  all <- do.call(rbind, mats)
  all <- all[is.finite(all[, 1L]) & is.finite(all[, 2L]), , drop = FALSE]
  if (!nrow(all)) {
    return(NULL)
  }
  dt <- data.table::data.table(
    mz_bin = round(all[, 1L], 4L),
    intensity = all[, 2L]
  )
  agg <- dt[, .(intensity = max(intensity, na.rm = TRUE)), by = mz_bin]
  data.table::setorder(agg, mz_bin)
  as.matrix(agg[, .(mz_bin, intensity)])
}


#' Resolve fragments.data (list or data.frame) into a uniform lookup structure
#' @noRd
.nts_resolve_fragments <- function(fr_raw, polarity) {
  if (is.data.frame(fr_raw)) return(fr_raw)
  if (!is.list(fr_raw)) return(NULL)
  pol_key <- if (identical(polarity, "neg")) "NEG" else "POS"
  sub <- fr_raw[[pol_key]]
  if (is.null(sub)) sub <- fr_raw[[tolower(pol_key)]]
  if (is.null(sub)) sub <- fr_raw
  sub
}

#' Extract diagnostic m/z for a given compound ID from resolved fragments
#' @noRd
.nts_get_diag_mz <- function(fr, rid) {
  if (is.null(fr)) return(NULL)
  if (is.data.frame(fr)) {
    if (nrow(fr) && "Lab.ID" %in% names(fr) && "mz" %in% names(fr)) {
      subf <- fr[as.character(fr[["Lab.ID"]]) == rid, , drop = FALSE]
      if (nrow(subf)) return(subf$mz)
    }
    return(NULL)
  }
  if (is.list(fr)) {
    entry <- fr[[rid]]
    if (is.null(entry)) return(NULL)
    if (is.data.frame(entry) && "mz" %in% names(entry)) return(entry$mz)
    if (is.numeric(entry)) return(entry)
    return(NULL)
  }
  NULL
}

#' Exact mass -> experimental m/z with proton shift
#' @noRd
.PROTON_MASS <- 1.007276

#' @noRd
.mass_to_mz <- function(mass, polarity) {
  if (identical(polarity, "neg")) {
    mass - .PROTON_MASS
  } else {
    mass + .PROTON_MASS
  }
}


#' @noRd
.nts_db_ms1_table <- function(database, polarity = "neg") {
  if (!methods::is(database, "databaseClass")) {
    stop("`database` must be databaseClass", call. = FALSE)
  }
  sn <- methods::slotNames(database)
  info <- NULL
  if ("spectra.info" %in% sn) {
    info <- methods::slot(database, "spectra.info")
  } else if ("compounds" %in% sn) {
    info <- methods::slot(database, "compounds")
  } else {
    stop("Database object lacks the metadata table required for identification.", call. = FALSE)
  }
  if (!is.data.frame(info)) {
    info <- tryCatch(
      as.data.frame(info, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
  }
  if (!is.data.frame(info) || !nrow(info)) {
    stop("Database metadata table is missing or empty.", call. = FALSE)
  }
  dt <- data.table::as.data.table(info)
  if (!"Lab.ID" %in% names(dt)) {
    stop("Metadata table must contain Lab.ID.", call. = FALSE)
  }

  # --- Smart RT: prefer RT (uppercase) over rt when rt is all-NA ----------
  nms <- names(dt)
  if ("RT" %in% nms && "rt" %in% nms) {
    dt[, rt := ifelse(!is.na(RT), RT, rt)]
  } else if ("RT" %in% nms && !"rt" %in% nms) {
    data.table::setnames(dt, "RT", "rt")
  }
  if (!"rt" %in% names(dt)) {
    stop("Metadata table must contain an rt or RT column.", call. = FALSE)
  }

  # --- Smart m/z: target_mz > polarity mz.neg / mz.pos > mz + shift -------
  # `target_mz` (ISF / adduct / override ion set by the in-house builder)
  # takes precedence so this legacy `databaseClass` path resolves the same
  # precursor m/z as the `NTSDatabase` path (`.nts_metadata_screening_mz`)
  # instead of silently falling back to `mass -/+ 1.007276`.
  pol_col <- if (identical(polarity, "neg")) "mz.neg" else "mz.pos"
  has_target <- "target_mz" %in% nms
  has_pol <- pol_col %in% nms
  has_mz <- "mz" %in% nms
  if (!has_target && !has_pol && !has_mz) {
    stop("Metadata table must contain one of target_mz, mz, mz.neg, or mz.pos.", call. = FALSE)
  }
  dt[, calc_mz := NA_real_]
  if (has_target) {
    dt[, calc_mz := suppressWarnings(as.numeric(target_mz))]
  }
  if (has_pol) {
    dt[!is.finite(calc_mz), calc_mz := suppressWarnings(as.numeric(get(pol_col)))]
  }
  if (has_mz) {
    dt[!is.finite(calc_mz), calc_mz := .mass_to_mz(suppressWarnings(as.numeric(mz)), polarity)]
    dt[, exact_mass := mz]
  }
  dt[!is.finite(calc_mz) | calc_mz <= 0, calc_mz := NA_real_]
  dt[, mz := calc_mz]
  dt[, calc_mz := NULL]
  as.data.frame(dt)
}


#' @noRd
.nts_consensus_to_peak_matrix <- function(sp) {
  if (is.null(sp$mz) || is.null(sp$intensity)) {
    return(NULL)
  }
  cbind(as.numeric(sp$mz), as.numeric(sp$intensity))
}


#' @noRd
.nts_exp_ms2_to_peak_matrices <- function(x, labels = NULL) {
  if (is.null(x)) {
    return(list())
  }
  if (is.list(x) && length(x) && !is.null(x$mz) && !is.null(x$intensity)) {
    mat <- .nts_consensus_to_peak_matrix(x)
    if (is.null(mat) || !nrow(mat)) {
      return(list())
    }
    out <- list(mat)
    if (!is.null(labels) && length(labels)) {
      names(out) <- as.character(labels)[1L]
    }
    return(out)
  }
  if (!is.list(x) || !length(x)) {
    return(list())
  }
  mats <- lapply(x, .nts_consensus_to_peak_matrix)
  ok <- !vapply(mats, function(m) is.null(m) || !nrow(m), logical(1L))
  mats <- mats[ok]
  if (!is.null(labels) && length(labels) == length(ok)) {
    names(mats) <- as.character(labels)[ok]
  }
  mats
}


#' @noRd
.score_ms2_subset_max <- function(exp_subset,
                                  ref_mat,
                                  precursor_mz,
                                  precursor_mz_ref,
                                  feature_tol,
                                  diagnostic_mz_lib = NULL,
                                  ...) {
  if (!length(exp_subset)) {
    return(NULL)
  }
  scores <- lapply(seq_along(exp_subset), function(i) {
    sc <- score_spectra_multidim(
      exp_subset[[i]],
      ref_mat,
      precursor_mz = precursor_mz,
      precursor_mz_ref = precursor_mz_ref,
      feature_tol = feature_tol,
      diagnostic_mz_lib = diagnostic_mz_lib,
      use_library_top_diagnostic_peaks = TRUE,
      ...
    )
    sc$best_ms2_index <- i
    sc$n_ms2_subset <- length(exp_subset)
    sc
  })
  comps <- vapply(scores, function(z) z$score_ms2_composite, numeric(1L))
  if (!length(comps) || all(!is.finite(comps))) {
    return(NULL)
  }
  best_i <- which.max(comps)
  scores[[best_i]]
}


#' @noRd
.nts_progress_start <- function(total, name = "progress") {
  if (!is.numeric(total) || length(total) != 1L || !is.finite(total) || total <= 0) {
    return(list(type = "none"))
  }
  pb <- utils::txtProgressBar(min = 0, max = as.integer(total), style = 3)
  list(type = "txt", pb = pb, total = as.integer(total))
}


#' @noRd
.nts_progress_update <- function(pb, current) {
  if (is.null(pb$type) || identical(pb$type, "none")) {
    return(invisible(NULL))
  }
  current <- min(as.integer(current), pb$total)
  if (identical(pb$type, "txt")) {
    utils::setTxtProgressBar(pb$pb, value = current)
  }
  invisible(NULL)
}


#' @noRd
.nts_progress_done <- function(pb) {
  if (is.null(pb$type) || identical(pb$type, "none")) {
    return(invisible(NULL))
  }
  if (identical(pb$type, "txt")) {
    close(pb$pb)
    cat("\n")
  }
  invisible(NULL)
}


#' All unique pairwise |m/z_i - m/z_j| from fragment masses (internal)
#' @noRd
.pairwise_mz_dist <- function(mz) {
  v <- sort(unique(mz[is.finite(mz)]))
  if (length(v) < 2L) {
    return(numeric(0L))
  }
  as.numeric(stats::dist(matrix(v, ncol = 1L)))
}


.validate_peak_matrix <- function(x, arg) {
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(arg, " must be a numeric matrix", call. = FALSE)
  }
  if (ncol(x) != 2L) {
    stop(arg, " must have exactly two columns (m/z, intensity)", call. = FALSE)
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.nts_ms2_cosine <- function(mz_unk, int_unk, mz_ref, int_ref, mz_tol_da) {
  mz_unk <- as.numeric(mz_unk)
  int_unk <- as.numeric(int_unk)
  mz_ref <- as.numeric(mz_ref)
  int_ref <- as.numeric(int_ref)
  tol <- as.numeric(mz_tol_da)[1L]
  if (!is.finite(tol) || tol < 0) {
    stop("mz_tol_da must be a non-negative finite number", call. = FALSE)
  }

  ok_u <- is.finite(mz_unk) & is.finite(int_unk) & int_unk > 0
  mz_unk <- mz_unk[ok_u]
  int_unk <- int_unk[ok_u]
  ok_r <- is.finite(mz_ref) & is.finite(int_ref) & int_ref > 0
  mz_ref <- mz_ref[ok_r]
  int_ref <- int_ref[ok_r]

  if (!length(mz_unk) || !length(mz_ref)) {
    return(list(score = 0, n_matched_peaks = 0L))
  }

  nu <- sqrt(sum(int_unk^2))
  if (!is.finite(nu) || nu <= 0) {
    return(list(score = 0, n_matched_peaks = 0L))
  }
  u <- int_unk / nu

  nr <- sqrt(sum(int_ref^2))
  if (!is.finite(nr) || nr <= 0) {
    return(list(score = 0, n_matched_peaks = 0L))
  }
  r_norm <- int_ref / nr

  aligned <- vapply(
    seq_along(mz_unk),
    function(i) {
      w <- abs(mz_ref - mz_unk[i]) <= tol
      if (!any(w)) {
        return(0)
      }
      max(r_norm[w])
    },
    numeric(1L)
  )

  matched_count <- sum(aligned > 0)

  na <- sqrt(sum(aligned^2))
  if (!is.finite(na) || na <= 0) {
    return(list(score = 0, n_matched_peaks = as.integer(matched_count)))
  }
  cosine_score <- sum(u * (aligned / na))
  list(score = cosine_score, n_matched_peaks = as.integer(matched_count))
}


#' Experimental m/z values in `exp_v` that match any library value within `tol` (Da).
#' @keywords internal
#' @noRd
.nts_get_matched <- function(exp_v, lib_v, tol) {
  exp_v <- as.numeric(exp_v)
  lib_v <- as.numeric(lib_v)
  tol <- as.numeric(tol)[1L]
  exp_v <- exp_v[is.finite(exp_v)]
  lib_v <- lib_v[is.finite(lib_v)]
  if (!length(lib_v) || !length(exp_v)) {
    return(numeric(0))
  }
  exp_v[vapply(exp_v, function(x) any(abs(x - lib_v) <= tol), logical(1L))]
}


#' @keywords internal
#' @noRd
.nts_fmt_matched_vec <- function(v) {
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (!length(v)) {
    return("")
  }
  paste(format(v, digits = 8, trim = TRUE), collapse = "; ")
}


#' Experimental DF / NL / MD m/z vectors from one peak matrix (col1 = mz).
#' @keywords internal
#' @noRd
.nts_exp_ms2_derived_features <- function(u_mz, em) {
  mz <- as.numeric(em[, 1L])
  mz <- mz[is.finite(mz)]
  exp_df <- mz
  u_mz <- as.numeric(u_mz)[1L]
  exp_nl <- numeric(0)
  if (is.finite(u_mz) && length(mz)) {
    exp_nl <- u_mz - mz
    exp_nl <- exp_nl[is.finite(exp_nl) & exp_nl > 0]
  }
  exp_md <- numeric(0)
  if (length(mz) >= 2L) {
    exp_md <- as.numeric(stats::dist(mz))
  }
  list(df = exp_df, nl = exp_nl, md = exp_md)
}


#' Library feature m/z values for one compound and polarity from `features_long`.
#' @keywords internal
#' @noRd
.nts_lib_feat_mz <- function(fl, lab_id, pol_flag, feature_type) {
  ftwant <- unique(as.character(feature_type))
  v <- fl[
    Lab_ID == lab_id &
      Polarity == pol_flag &
      tolower(as.character(Feature_Type)) %in% tolower(ftwant),
    mz_value
  ]
  v <- as.numeric(v)
  unique(v[is.finite(v)])
}


#' Sub-feature scores and matched-fragment strings (DF treated as diagnostic).
#' @keywords internal
#' @noRd
.nts_subfeature_scores <- function(u_mz,
                                        em,
                                        fl,
                                        lab_id,
                                        pol_flag,
                                        ms2_tol_da) {
  der <- .nts_exp_ms2_derived_features(u_mz, em)
  lib_df <- .nts_lib_feat_mz(fl, lab_id, pol_flag, c("DF", "Diagnostic", "diagnostic"))
  lib_nl <- .nts_lib_feat_mz(fl, lab_id, pol_flag, c("NL", "Neutral_loss", "neutral_loss"))
  lib_md <- .nts_lib_feat_mz(
    fl, lab_id, pol_flag,
    c("MD", "Mass_difference", "mass_difference", "MassDiff", "mass_diff")
  )

  m_df <- .nts_get_matched(der$df, lib_df, ms2_tol_da)
  m_nl <- .nts_get_matched(der$nl, lib_nl, ms2_tol_da)
  m_md <- .nts_get_matched(der$md, lib_md, ms2_tol_da)

  score_diagnostic <- if (length(lib_df) > 0L) {
    length(m_df) / length(lib_df)
  } else {
    0
  }
  score_nl <- if (length(lib_nl) > 0L) {
    length(m_nl) / length(lib_nl)
  } else {
    0
  }
  score_md <- if (length(lib_md) > 0L) {
    length(m_md) / length(lib_md)
  } else {
    0
  }

  list(
    score_diagnostic = score_diagnostic,
    score_nl = score_nl,
    score_md = score_md,
    matched_diagnostic = .nts_fmt_matched_vec(m_df),
    matched_nl = .nts_fmt_matched_vec(m_nl),
    matched_md = .nts_fmt_matched_vec(m_md)
  )
}


#' @keywords internal
#' @noRd
.nts_meta_precursor_table <- function(metadata_tb, polarity) {
  meta <- data.table::as.data.table(metadata_tb)
  if (!nrow(meta)) {
    return(data.table::data.table(
      Lab_ID = character(0),
      prec_mz = numeric(0),
      ann_name = character(0),
      lib_rt = numeric(0),
      db_type = character(0)
    ))
  }
  if ("Lab.ID" %in% names(meta) && !"Lab_ID" %in% names(meta)) {
    data.table::setnames(meta, "Lab.ID", "Lab_ID", skip_absent = TRUE)
  }
  if (!"Lab_ID" %in% names(meta)) {
    stop("db@metadata must contain Lab_ID (or Lab.ID)", call. = FALSE)
  }
  meta[, Lab_ID := as.character(Lab_ID)]
  meta[, prec_mz := NA_real_]
  # Canonical precursor m/z resolution: target_mz > mz_neg/mz.neg >
  # mz_pos/mz_neg > Precursor_mz > mz > mass +/- 1.007276 (see
  # .nts_metadata_screening_mz in 01_db_builder.R).
  meta[, prec_mz := .nts_metadata_screening_mz(meta, polarity)]
  ann_vec <- rep(NA_character_, nrow(meta))
  for (cn in c(
    "Compound_name",
    "Compound.name",
    "compound_name",
    "Name",
    "name"
  )) {
    if (cn %in% names(meta)) {
      vv <- as.character(meta[[cn]])
      miss <- is.na(ann_vec) | !nzchar(ann_vec)
      ann_vec[miss] <- vv[miss]
    }
  }
  ann_vec[is.na(ann_vec)] <- ""
  meta[, ann_name := ann_vec]

  lib_rt <- rep(NA_real_, nrow(meta))
  for (cn in c("RT", "rt", "Rt", "retention_time", "RetentionTime")) {
    if (cn %in% names(meta)) {
      lib_rt <- as.numeric(meta[[cn]])
      break
    }
  }
  meta[, lib_rt := lib_rt]

  # Carry database provenance (inhouse / online / prediction) through to
  # db_annotations so that assign_confidence_levels can pick the right tier
  # (L1 vs L2a vs L2b) without re-guessing from the Lab_ID pattern.
  # Missing column -> NA (downstream falls back to pattern detection).
  if ("db_type" %in% names(meta)) {
    dbt <- tolower(trimws(as.character(meta$db_type)))
    dbt[!dbt %in% c("inhouse", "online", "prediction")] <- NA_character_
    meta[, db_type := dbt]
  } else {
    meta[, db_type := NA_character_]
  }

  unique(meta[, .(Lab_ID, prec_mz, ann_name, lib_rt, db_type)], by = "Lab_ID")
}


#' First usable SMILES string per `Lab_ID` from [NTSDatabase] metadata.
#' @keywords internal
#' @noRd
.nts_metadata_smiles_by_lab <- function(db) {
  meta <- data.table::as.data.table(methods::slot(db, "metadata"))
  if (!nrow(meta)) {
    return(data.table::data.table(Lab_ID = character(0), SMILES = character(0)))
  }
  if ("Lab.ID" %in% names(meta) && !"Lab_ID" %in% names(meta)) {
    data.table::setnames(meta, "Lab.ID", "Lab_ID", skip_absent = TRUE)
  }
  if (!"Lab_ID" %in% names(meta)) {
    return(data.table::data.table(Lab_ID = character(0), SMILES = character(0)))
  }
  sm_col <- NULL
  for (cn in c("SMILES", "smiles", "Smiles")) {
    if (cn %in% names(meta)) {
      sm_col <- cn
      break
    }
  }
  meta[, Lab_ID := as.character(Lab_ID)]
  if (is.null(sm_col)) {
    return(unique(meta[, .(Lab_ID, SMILES = NA_character_)], by = "Lab_ID"))
  }
  cn <- sm_col
  meta[
    ,
    .(SMILES = {
      v <- trimws(as.character(get(cn)))
      v <- v[!is.na(v) & nzchar(v)]
      if (length(v)) v[1L] else NA_character_
    }),
    by = Lab_ID
  ]
}


#' @keywords internal
#' @noRd
.nts_polarity_to_spectra_flag <- function(polarity) {
  if (identical(match.arg(polarity, c("neg", "pos")), "neg")) {
    "NEG"
  } else {
    "POS"
  }
}


#' @keywords internal
#' @noRd
.nts_nodes_df_for_annotation <- function(project) {
  ft <- get_feature_table(project)
  if (!is.data.frame(ft) || !nrow(ft)) {
    stop("feature_table is empty; cannot annotate by spectral match", call. = FALSE)
  }
  idc <- if ("name" %in% names(ft)) {
    "name"
  } else {
    "feature_id"
  }
  data.frame(
    id = as.character(ft[[idc]]),
    mz = as.numeric(ft$mz),
    rt = if ("rt" %in% names(ft)) as.numeric(ft$rt) else NA_real_,
    stringsAsFactors = FALSE
  )
}


#' @keywords internal
#' @noRd
.nts_feature_ms2_peak_list <- function(project,
                                      feature_id,
                                      ms2_source = c("clean_ms2", "ms2data"),
                                      with_labels = FALSE) {
  ms2_source <- match.arg(ms2_source)
  fid <- as.character(feature_id)[1L]
  cm <- methods::slot(project, "clean_ms2")
  if (identical(ms2_source, "ms2data")) {
    cm <- list()
  }
  if (is.list(cm) && length(cm)) {
    for (el in cm) {
      if (!is.list(el)) {
        next
      }
      eid <- el$feature_id
      if (is.null(eid)) {
        eid <- el$name
      }
      if (is.null(eid)) {
        next
      }
      if (as.character(eid) != fid) {
        next
      }
      sp <- el$spectrum
      if (is.data.frame(sp) && nrow(sp) && all(c("mz", "intensity") %in% names(sp))) {
        m <- cbind(as.numeric(sp$mz), as.numeric(sp$intensity))
        m <- m[is.finite(m[, 1L]) & is.finite(m[, 2L]) & m[, 2L] > 0, , drop = FALSE]
        if (nrow(m)) {
          lab <- if (!is.null(el$ce_levels) && length(el$ce_levels)) {
            ce_num <- suppressWarnings(as.numeric(el$ce_levels))
            ce_num <- ce_num[is.finite(ce_num)]
            if (length(ce_num)) {
              paste0("CE", paste(sort(unique(ce_num)), collapse = "/"))
            } else {
              "clean_ms2"
            }
          } else {
            "clean_ms2"
          }
          mats <- list(m)
          names(mats) <- lab
          if (isTRUE(with_labels)) {
            return(list(mats = mats, labels = lab))
          }
          return(mats)
        }
      }
    }
  }
  cons <- methods::slot(project, "ms2data")[[fid]]
  labs <- NULL
  if (is.list(cons) && length(cons) && is.null(cons$mz)) {
    labs <- vapply(cons, function(e) {
      ce <- suppressWarnings(as.numeric(e$ce)[1L])
      nm <- if (is.finite(ce)) paste0("CE", ce) else "CE"
      sid <- e$scan_idx
      if (!is.null(sid) && length(sid) && !is.na(sid[1L])) {
        nm <- paste0(nm, "#", as.character(sid[1L]))
      }
      # `scan_idx` is only unique within one source file, and a single feature
      # routinely aggregates MS2 scans from several files: qualify the label
      # with the file stem so every experimental spectrum stays unambiguous in
      # `@ms2_score_grid` and `db_annotations$exp_CE`.
      src <- e$source_file
      if (!is.null(src) && length(src) && !is.na(src[1L]) &&
          nzchar(as.character(src[1L]))) {
        stem <- sub("\\.[^.]*$", "", basename(as.character(src[1L])))
        if (nzchar(stem)) {
          nm <- paste0(nm, "@", stem)
        }
      }
      nm
    }, character(1L))
    # Last-resort uniqueness so no two experimental spectra of a feature can
    # collide in the score grid.
    labs <- make.unique(labs, sep = "~")
  }
  mats <- .nts_exp_ms2_to_peak_matrices(cons, labels = labs)
  if (isTRUE(with_labels)) {
    return(list(mats = mats, labels = names(mats)))
  }
  mats
}


#' @keywords internal
#' @noRd
.nts_merged_ref_matrix <- function(spectra_long, lab_id, pol_flag,
                                 ref_spectra = c("merged", "all")) {
  ref_spectra <- match.arg(ref_spectra)
  # Prefer the pre-computed "Merged" consensus spectrum; when a compound
  # only has per-CE rows (e.g. CE10/CE20/CE40), combine all CE levels.
  sl <- spectra_long[Lab_ID == lab_id & Polarity == pol_flag]
  if (identical(ref_spectra, "merged") && nrow(sl)) {
    sl_m <- sl[NCE == "Merged"]
    if (nrow(sl_m)) {
      sl <- sl_m
    }
  }
  if (!nrow(sl)) {
    return(NULL)
  }
  m <- cbind(as.numeric(sl$frag_mz), as.numeric(sl$intensity))
  m <- m[is.finite(m[, 1L]) & is.finite(m[, 2L]) & m[, 2L] > 0, , drop = FALSE]
  if (!nrow(m)) {
    return(NULL)
  }
  # Collapse duplicate m/z (max intensity) and sort ascending so CE-level
  # spectra combine into a single consensus reference matrix.
  if (anyDuplicated(m[, 1L])) {
    dt <- data.table::data.table(mz = m[, 1L], intensity = m[, 2L])
    dt <- dt[, .(intensity = max(intensity)), by = mz]
    m <- as.matrix(dt[order(mz)])
  } else if (is.unsorted(m[, 1L])) {
    m <- m[order(m[, 1L]), , drop = FALSE]
  }
  m
}


#' Per-collision-energy reference peak matrices for one compound
#'
#' Builds one reference peak matrix per distinct `NCE` level of a compound so
#' that MS2 matching can score the experimental spectra against each collision
#' energy separately (a `(library CE) x (experimental spectrum)` grid). Falls
#' back to the single collapsed spectrum (`.nts_merged_ref_matrix`) for the
#' legacy `"merged"` / `"all"` strategies.
#'
#' @return A named `list` of two-column peak matrices, names being the `NCE`
#'   labels (`"CE10"`, `"CE20"`, ...). `"Merged"`/`"InHouse"` rows are used only
#'   when a compound has no non-Merged CE rows.
#' @keywords internal
#' @noRd
.nts_ref_matrices_per_ce <- function(spectra_long, lab_id, pol_flag,
                                     ref_spectra = c("per_ce", "merged", "all")) {
  ref_spectra <- match.arg(ref_spectra)
  if (!identical(ref_spectra, "per_ce")) {
    m <- .nts_merged_ref_matrix(spectra_long, lab_id, pol_flag, ref_spectra)
    if (is.null(m)) {
      return(list())
    }
    out <- list(m)
    names(out) <- if (identical(ref_spectra, "merged")) "Merged" else "all"
    return(out)
  }

  sl <- data.table::as.data.table(spectra_long)
  sl <- sl[Lab_ID == lab_id & Polarity == pol_flag]
  if (!nrow(sl)) {
    return(list())
  }
  nce <- as.character(sl[["NCE"]])
  ce_vals <- unique(nce[!is.na(nce) & nzchar(nce) & nce != "Merged"])
  if (!length(ce_vals)) {
    ce_vals <- unique(nce[!is.na(nce) & nzchar(nce)])
  }
  if (!length(ce_vals)) {
    ce_vals <- "Merged"
  }

  out <- list()
  for (ce in ce_vals) {
    sub <- sl[nce == ce]
    m <- cbind(as.numeric(sub[["frag_mz"]]), as.numeric(sub[["intensity"]]))
    m <- m[is.finite(m[, 1L]) & is.finite(m[, 2L]) & m[, 2L] > 0, , drop = FALSE]
    if (!nrow(m)) {
      next
    }
    if (anyDuplicated(m[, 1L])) {
      dt <- data.table::data.table(mz = m[, 1L], intensity = m[, 2L])
      dt <- dt[, .(intensity = max(intensity)), by = mz]
      m <- as.matrix(dt[order(mz)])
    } else if (is.unsorted(m[, 1L])) {
      m <- m[order(m[, 1L]), , drop = FALSE]
    }
    out[[ce]] <- m
  }
  out
}


#' @keywords internal
#' @noRd
.nts_ms1_candidates_library <- function(nodes_dt, meta_dt, mz_tol_ppm) {
  ppm <- as.numeric(mz_tol_ppm)[1L] * 1e-6
  if (!is.finite(ppm) || ppm < 0) {
    stop("mz_tol_ppm must be non-negative finite", call. = FALSE)
  }
  nodes_dt <- data.table::as.data.table(nodes_dt)
  meta_dt <- data.table::as.data.table(meta_dt)
  nodes_dt <- nodes_dt[is.finite(mz_obs)]
  meta_dt <- meta_dt[is.finite(prec_mz)]
  if (!"rt_obs" %in% names(nodes_dt)) {
    nodes_dt[, rt_obs := NA_real_]
  } else {
    nodes_dt[, rt_obs := as.numeric(rt_obs)]
  }
  if (!nrow(nodes_dt) || !nrow(meta_dt)) {
    return(data.table::data.table(
      node_id = character(0),
      mz_obs = numeric(0),
      rt_obs = numeric(0),
      Lab_ID = character(0),
      prec_mz = numeric(0),
      ann_name = character(0),
      lib_rt = numeric(0)
    ))
  }
  nodes_dt[, k := 1L]
  meta_dt[, k := 1L]
  j <- merge(nodes_dt, meta_dt, by = "k", allow.cartesian = TRUE)
  j[, k := NULL]
  j <- j[abs(prec_mz - mz_obs) <= mz_obs * ppm]
  j[]
}


#' Assemble the `@ms2_score_grid` table from per-feature cell lists
#'
#' Binds the per-feature `(library CE) x (experimental spectrum)` cell tables
#' collected by [annotate_by_spectral_match], orders the columns and sorts rows
#' by feature / compound / CE / spectrum. Returns `data.frame()` when nothing was
#' scored (or when `return_score_grid = FALSE`).
#'
#' @param rows List of `data.table`s, one per scored `(feature, compound)` pair.
#' @return A `data.frame` (0-row when `rows` is empty).
#' @keywords internal
#' @noRd
.nts_ms2_grid_df <- function(rows) {
  if (!length(rows)) {
    return(data.frame())
  }
  dt <- data.table::rbindlist(rows, fill = TRUE)
  cols <- c(
    "name", "ref_id", "lib_CE", "exp_spectrum", "score_cosine",
    "n_matched_peaks", "passes_min_score", "is_winner", "score_ms2",
    "score_ndp", "ndp_raw", "score_nl", "score_md", "score_diagnostic"
  )
  data.table::setcolorder(dt, intersect(cols, names(dt)))
  data.table::setorderv(
    dt,
    intersect(c("name", "ref_id", "lib_CE", "exp_spectrum"), names(dt)),
    na.last = TRUE
  )
  as.data.frame(dt, stringsAsFactors = FALSE)
}


## Keep `@ms2_score_grid`'s `is_winner` flags consistent with the *current*
## `@db_annotations`. `annotate_by_spectral_match()` writes the grid and the
## full per-feature x compound annotation in one pass, but `batch_annotate()`
## afterwards *deduplicates* `@db_annotations` (one row per `ref_id`, best
## `score_total`). The dropped feature x compound pairs must lose their
## `is_winner` flag, otherwise the grid would advertise winners for rows that
## no longer exist. Only flags are cleared here -- the cells stay, so an
## offline re-tune of `min_score` still sees every scored (CE x spectrum) cell.
.nts_sync_grid_winners <- function(project) {
  g <- get_ms2_score_grid(project)
  if (!is.data.frame(g) || !nrow(g) || !"is_winner" %in% names(g)) {
    return(project)
  }
  ann <- methods::slot(project, "db_annotations")
  if (!is.data.frame(ann) || !nrow(ann) ||
      !all(c("name", "ref_id") %in% names(ann))) {
    g$is_winner <- FALSE
    return(set_ms2_score_grid(project, g))
  }
  keep <- paste(ann$name, ann$ref_id, sep = "\r")
  g$is_winner <- g$is_winner &
    (paste(g$name, g$ref_id, sep = "\r") %in% keep)
  set_ms2_score_grid(project, g)
}


#' Spectral library annotation using a flat [NTSDatabase]
#'
#' MS1 candidate search uses keyed `data.table` joins on precursor m/z
#' (`ppm`), then MS2 scoring against every library collision energy
#' (`.nts_ref_matrices_per_ce`) crossed with every experimental MS2 spectrum of
#' the feature. For each `(library CE) x (experimental spectrum)` cell both the
#' cosine (`score_cosine_ms2`) and the multidimensional composite
#' (`score_spectra_multidim`: `score_ms2` / `score_ndp` / `ndp_raw`, plus
#' `score_nl` / `score_md` / `score_diagnostic`) are computed; the cell with the
#' highest total score among those whose cosine clears `min_score` is kept and the
#' winning `lib_CE` / `exp_CE` are recorded. Each experimental spectrum is
#' labelled `CE<ce>#<scan_idx>@<file stem>` so that scans sharing a `scan_idx`
#' across different source files stay distinguishable in `exp_CE` and in
#' `@ms2_score_grid$exp_spectrum`.
#' All passing compounds are written to `db_annotations` (a feature may match
#' several library entries); the network nodes keep one row per feature. When
#' `SMILES` exists in `db@metadata`, it is copied per `Lab_ID` into
#' `db_annotations`.
#'
#' @param project An [NTSProject] with `ms2data` and/or `clean_ms2`.
#' @param db An [NTSDatabase] with `metadata` (precursor m/z per `Lab_ID`) and
#'   `spectra_long` containing per-CE (`NCE`) and/or `NCE == \"Merged\"` rows.
#' @param diagnostic_ions Optional `data.frame` with columns `mz`, `type`
#'   (`"fragment"` or `"nl"`), optional `name`. When non-`NULL`, `score_diagnostic`
#'   and `matched_diagnostic` are derived from experimental peaks vs this table
#'   (otherwise from library `features_long` diagnostic features).
#' @param mz_tol_ppm MS1 ppm window vs observed feature `mz` (also scales
#'   linear `score_mz = max(0, 1 - mz_error_ppm / mz_tol_ppm)`).
#' @param rt_tol_s RT tolerance (same units as metadata / nodes `rt`) for
#'   linear `score_rt = max(0, 1 - rt_error_s / rt_tol_s)` when both RTs exist.
#' @param ms2_tol_da Absolute m/z tolerance (Da) for fragment alignment and
#'   DF/NL/MD feature matching.
#' @param min_score Minimum cosine score for a `(library CE) x (experimental
#'   spectrum)` cell to be eligible. A candidate is kept when at least one cell
#'   clears the gate, and the reported cell is then the highest-total eligible
#'   cell (so the reported `score_cosine_ms2` is always `>= min_score`).
#' @param polarity `"neg"` or `"pos"`; selects `NEG`/`POS` rows in
#'   `spectra_long` and `mz_neg`/`mz_pos` in metadata when applicable.
#' @param ms2_source Which experimental MS2 slot to score against:
#'   `"ms2data"` (default; raw extracted spectra, maximal coverage) or
#'   `"clean_ms2"` (purified/consensus spectra; falls back to `ms2data`
#'   per feature when absent).
#' @param ref_spectra Library spectrum strategy: `"per_ce"` (default; match
#'   against each distinct `NCE` level separately, i.e. a `(library CE) x
#'   (experimental spectrum)` grid, falling back to the `"Merged"`/`"InHouse"`
#'   row only when a compound has no per-CE spectra), `"merged"` (use the
#'   pre-computed `NCE == "Merged"` consensus row, falling back to combining all
#'   per-CE rows) or `"all"` (always combine all CE rows, max intensity per m/z).
#' @param return_score_grid Logical (default `TRUE`). Record **every**
#'   `(library CE) x (experimental spectrum)` cell in `@ms2_score_grid`, not just
#'   the winner and not just cells above `min_score`. This lets a mismatched
#'   experimental MS2 subset be inspected directly and lets `min_score` be
#'   re-tuned offline. Set `FALSE` to skip the (small) extra bookkeeping and
#'   leave `@ms2_score_grid` empty.
#' @param validate_diagnostic_fragments Logical (default `TRUE`). After scoring,
#'   run the packaged diagnostic-fragment validation (see
#'   [default_diagnostic_fragments_csv]; PFAS-focused expert table): matched
#'   expert fragments are recorded in `@ms2_annotations$diagnostic_fragments`
#'   and `@ms1_annotations$annotation_level` is tagged `[MS2 Verified]` for the
#'   matching features. Uses `@ms2data`. Set `FALSE` to skip entirely.
#' @param diagnostic_fragments_csv Path to the diagnostic CSV
#'   (`fragment_formula`, `fragment_mass`) used when
#'   `validate_diagnostic_fragments = TRUE`. `NULL` (default) resolves to the
#'   packaged table; an unresolvable path disables the step silently.
#' @param detect_pfas_series Logical (default `TRUE`). After scoring, run the
#'   library-free CF2 fragment-ladder / series detection
#'   ([detect_pfas_fragment_series]) and write its columns into
#'   `@ms2_annotations`. This only **adds columns** and never changes the
#'   spectral-library hits or the confidence tiers (structural guardrail). Set
#'   `FALSE` to skip.
#' @param pfas_series_min_len,cf2_tol_da,max_gap_units,min_rel_int,series_n_max,mode
#'   Tuning for [detect_pfas_fragment_series].
#'
#' @return The updated `project` (in-place slot updates). `@db_annotations` holds
#'   the merged best-cell result per feature x compound; `@ms2_score_grid` holds
#'   the full per-cell grid (see [get_ms2_score_grid]).
#' @export
annotate_by_spectral_match <- function(project,
                                       db,
                                       mz_tol_ppm = 10,
                                       rt_tol_s = 150,
                                       ms2_tol_da = 0.02,
                                       min_score = 0.5,
                                       polarity = c("neg", "pos"),
                                       diagnostic_ions = NULL,
                                       ms2_source = c("ms2data", "clean_ms2"),
                                       ref_spectra = c("per_ce", "merged", "all"),
                                       return_score_grid = TRUE,
                                       validate_diagnostic_fragments = TRUE,
                                       diagnostic_fragments_csv = NULL,
                                       detect_pfas_series = TRUE,
                                       pfas_series_min_len = 3L,
                                       cf2_tol_da = 0.01,
                                       max_gap_units = 3L,
                                       min_rel_int = 0.01,
                                       series_n_max = 10L,
                                       mode = c("loose", "strict")) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  .nts_check(db, "reference_db", .arg = "db")

  polarity <- match.arg(polarity)
  ms2_source <- match.arg(ms2_source)
  ref_spectra <- match.arg(ref_spectra)
  return_score_grid <- isTRUE(return_score_grid)
  validate_diagnostic_fragments <- isTRUE(validate_diagnostic_fragments)
  # Default-on library-free CF2 fragment-ladder / series detection (design
  # §5.3, decision #1). It only ADDS @ms2_annotations columns and routes the
  # confidence impact through the structural guardrail, so it never changes the
  # spectral-library hits or the confidence tiers.
  detect_pfas_series <- isTRUE(detect_pfas_series)
  mode <- match.arg(mode)
  # Default-on diagnostic-fragment validation: applied on every return path so
  # the packaged expert table (PFAS-focused) always annotates `@ms2_annotations`
  # unless explicitly disabled. See .nts_apply_default_diagnostic_fragments().
  diag_validate <- function(p) {
    if (validate_diagnostic_fragments) {
      p <- .nts_apply_default_diagnostic_fragments(
        p, diagnostic_fragments_csv, mz_tol_ppm
      )
    }
    # Run after the diagnostic table so the two steps write different columns
    # and never clobber each other (design §5.3).
    if (detect_pfas_series) {
      p <- detect_pfas_fragment_series(
        p,
        cf2_tol_da = cf2_tol_da,
        min_len = pfas_series_min_len,
        max_gap_units = max_gap_units,
        min_rel_int = min_rel_int,
        series_n_max = series_n_max,
        mode = mode,
        mz_tol_ppm = mz_tol_ppm,
        polarity = polarity,
        quiet = TRUE
      )
    }
    p
  }
  pol_flag <- .nts_polarity_to_spectra_flag(polarity)
  ms2_tol_da <- as.numeric(ms2_tol_da)[1L]
  min_score <- as.numeric(min_score)[1L]
  rt_tol_s <- as.numeric(rt_tol_s)[1L]
  if (!is.finite(rt_tol_s) || rt_tol_s <= 0) {
    stop("rt_tol_s must be positive finite", call. = FALSE)
  }

  nodes_df <- .nts_nodes_df_for_annotation(project)
  nodes_dt <- data.table::as.data.table(nodes_df)
  nodes_dt[, node_id := as.character(id)]
  nodes_dt[, mz_obs := as.numeric(mz)]
  nodes_dt[, rt_obs := NA_real_]
  if ("rt" %in% names(nodes_dt)) {
    nodes_dt[, rt_obs := as.numeric(rt)]
  }

  meta_dt <- .nts_meta_precursor_table(methods::slot(db, "metadata"), polarity)
  cand <- .nts_ms1_candidates_library(nodes_dt, meta_dt, mz_tol_ppm)
  if (!nrow(cand)) {
    methods::slot(project, "db_annotations") <- data.frame()
    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "annotate_by_spectral_match",
      list(mz_tol_ppm = mz_tol_ppm, n_ms1_candidates = 0L)
    )
    return(invisible(diag_validate(project)))
  }

  sl <- data.table::as.data.table(methods::slot(db, "spectra_long"))
  if (!all(c("Lab_ID", "Polarity", "NCE", "frag_mz", "intensity") %in% names(sl))) {
    stop("db@spectra_long has unexpected schema", call. = FALSE)
  }

  fl <- data.table::as.data.table(methods::slot(db, "features_long"))
  fneed <- c("Lab_ID", "Polarity", "Feature_Type", "mz_value")
  if (nrow(fl) && !all(fneed %in% names(fl))) {
    stop(
      "db@features_long must contain columns: ",
      paste(fneed, collapse = ", "),
      call. = FALSE
    )
  }
  if (!nrow(fl)) {
    fl <- data.table::data.table(
      Lab_ID = character(0),
      Polarity = character(0),
      Feature_Type = character(0),
      mz_value = numeric(0)
    )
  }

  # Constant across the candidate loop: validated user diagnostic ions + tol.
  d_usr <- .nts_validate_diagnostic_ions(diagnostic_ions, "diagnostic_ions")
  mz_align <- suppressWarnings(max(ms2_tol_da, 0.005, na.rm = TRUE))
  if (!is.finite(mz_align)) {
    mz_align <- 0.02
  }

  score_rows <- list()
  score_grid_rows <- list()
  ri <- 0L
  pb <- .nts_progress_start(nrow(cand), "annotate_by_spectral_match")
  on.exit(.nts_progress_done(pb), add = TRUE)
  for (k in seq_len(nrow(cand))) {
    .nts_progress_update(pb, k)
    nid <- cand$node_id[k]
    lid <- as.character(cand$Lab_ID[k])
    # One reference spectrum per library collision energy (legacy strategies
    # "merged"/"all" collapse to a single reference).
    refs_ce <- .nts_ref_matrices_per_ce(sl, lid, pol_flag, ref_spectra)
    if (!length(refs_ce)) {
      next
    }
    exp_list <- .nts_feature_ms2_peak_list(project, nid, ms2_source, with_labels = TRUE)
    exp_mats <- exp_list$mats
    if (!length(exp_mats)) {
      next
    }
    exp_labs <- exp_list$labels
    # Number of experimental spectra this feature offered to the CE x spectrum
    # grid (post label de-duplication). Recorded alongside the winner so the
    # mirror plot can prove it used the same subset the scoring used.
    n_sub <- length(exp_mats)

    pm <- as.numeric(cand$prec_mz[k])
    om <- as.numeric(cand$mz_obs[k])
    if (!is.finite(om)) {
      next
    }

    # (library CE) x (experimental spectrum) grid: score every pair and keep
    # the cell with the highest multidimensional MS2 composite. score_total is
    # constant in score_mz/score_rt across cells, so its ordering is identical.
    # min_score gates each cell's cosine, so the winner is the highest-composite
    # cell *among those that clear the gate*: a compound that the legacy
    # single-reference path would have kept can never be dropped here, and every
    # reported cosine is therefore guaranteed to be >= min_score.
    best_total <- -Inf
    best <- NA_real_
    best_n <- 0L
    best_em <- NULL
    best_ce <- NA_character_
    best_explab <- NA_character_
    # Position of the winning experimental spectrum within `exp_mats`
    # (1-based). Written out as `best_ms2_index` so downstream consumers -
    # notably plot_ms2_mirror() - redraw the exact spectrum that was scored
    # instead of blindly falling back to subset #1.
    best_j <- NA_integer_
    sc_best <- NULL
    grid_rows <- list()
    gi <- 0L
    best_gi <- 0L
    for (ce_nm in names(refs_ce)) {
      ref_m <- refs_ce[[ce_nm]]
      for (j in seq_along(exp_mats)) {
        em <- exp_mats[[j]]
        res <- .nts_ms2_cosine(
          em[, 1L],
          em[, 2L],
          ref_m[, 1L],
          ref_m[, 2L],
          ms2_tol_da
        )
        sc <- res$score
        gate <- is.finite(sc) && sc >= min_score

        # The multidimensional composite is only needed for cells that clear the
        # gate (winner selection) - the grid below still records every cell's
        # cosine, so `min_score` can be re-tuned offline without re-scoring.
        sc_k <- NULL
        if (gate) {
          sc_k <- tryCatch(
            score_spectra_multidim(
              em,
              ref_m,
              precursor_mz = om,
              precursor_mz_ref = pm,
              feature_tol = ms2_tol_da,
              mz_tol_align = mz_align,
              int_thresh = 0.025,
              weights = c(ndp = 0.50, nl = 0.20, md = 0.15, diag = 0.15),
              diagnostic_ions = d_usr,
              use_library_top_diagnostic_peaks = FALSE
            ),
            error = function(e) NULL
          )
        }
        ms2_k <- if (!is.null(sc_k)) {
          suppressWarnings(as.numeric(sc_k$score_ms2_composite)[1L])
        } else {
          NA_real_
        }
        if (gate && !is.finite(ms2_k)) {
          ms2_k <- sc
        }

        exp_lab_j <- if (!is.null(exp_labs) && length(exp_labs) >= j) {
          exp_labs[j]
        } else {
          NA_character_
        }
        if (isTRUE(return_score_grid)) {
          gi <- gi + 1L
          grid_rows[[gi]] <- data.table::data.table(
            name = nid,
            ref_id = lid,
            lib_CE = as.character(ce_nm),
            exp_spectrum = as.character(exp_lab_j),
            score_cosine = sc,
            n_matched_peaks = as.integer(res$n_matched_peaks)[1L],
            passes_min_score = gate,
            score_ms2 = ms2_k,
            score_ndp = if (!is.null(sc_k)) suppressWarnings(as.numeric(sc_k$score_ndp)[1L]) else NA_real_,
            ndp_raw = if (!is.null(sc_k)) suppressWarnings(as.numeric(sc_k$ndp)[1L]) else NA_real_,
            score_nl = if (!is.null(sc_k)) suppressWarnings(as.numeric(sc_k$score_nl)[1L]) else NA_real_,
            score_md = if (!is.null(sc_k)) suppressWarnings(as.numeric(sc_k$score_md)[1L]) else NA_real_,
            score_diagnostic = if (!is.null(sc_k)) suppressWarnings(as.numeric(sc_k$score_diagnostic)[1L]) else NA_real_
          )
        }

        if (!gate) {
          next
        }
        if (is.null(best_em) || ms2_k > best_total) {
          best_total <- ms2_k
          best <- sc
          best_n <- as.integer(res$n_matched_peaks)[1L]
          best_em <- em
          best_ce <- ce_nm
          best_explab <- as.character(exp_lab_j)
          best_j <- j
          sc_best <- sc_k
          # Remember the winning cell's position in `grid_rows` (0 when the
          # grid is not being recorded): marking by index guarantees exactly one
          # `is_winner` row even when two experimental spectra share a label.
          best_gi <- gi
        }
      }
    }
    if (length(grid_rows)) {
      gdt <- data.table::rbindlist(grid_rows)
      gdt[, is_winner := FALSE]
      if (!is.null(best_em) && best_gi >= 1L && best_gi <= nrow(gdt)) {
        gdt[best_gi, is_winner := TRUE]
      }
      score_grid_rows[[length(score_grid_rows) + 1L]] <- gdt
    }
    if (is.null(best_em) || !is.finite(best)) {
      next
    }

    mz_eppm <- if (is.finite(pm) && pm > 0 && is.finite(om)) {
      abs(om - pm) / pm * 1e6
    } else {
      NA_real_
    }
    rte <- as.numeric(cand$rt_obs[k])
    rtl <- as.numeric(cand$lib_rt[k])
    rt_es <- if (is.finite(rte) && is.finite(rtl)) {
      abs(rte - rtl)
    } else {
      NA_real_
    }

    # The winning (CE x spectrum) pair drives every reported score, so cosine,
    # NL/MD/diagnostic and the multidimensional composite stay mutually
    # consistent.
    score_cosine_ms2 <- best

    sub <- .nts_subfeature_scores(
      om,
      best_em,
      fl,
      lid,
      pol_flag,
      ms2_tol_da
    )

    score_diag <- sub$score_diagnostic
    matched_diag <- sub$matched_diagnostic
    if (!is.null(d_usr)) {
      usr <- .nts_score_exp_vs_diagnostic_mz(best_em[, 1L], diagnostic_ions, om, ms2_tol_da)
      if (is.finite(usr$score)) {
        score_diag <- usr$score
        matched_diag <- usr$matched
      }
    }

    # Reuse the winning (CE x spectrum) cell's composite instead of recomputing.
    sc_multi <- sc_best

    if (!is.null(sc_multi)) {
      score_ms2 <- sc_multi$score_ms2_composite
      score_ndp <- sc_multi$score_ndp
      ndp_raw <- sc_multi$ndp
      score_nl <- sc_multi$score_nl
      score_md <- sc_multi$score_md
      sc_diag_multi <- suppressWarnings(as.numeric(sc_multi$score_diagnostic)[1L])
      if (is.finite(sc_diag_multi) && sc_diag_multi > 0) {
        score_diag <- sc_diag_multi
      }
      snl <- as.character(sc_multi$matched_nl)[1L]
      smd <- as.character(sc_multi$matched_md)[1L]
      sdi <- as.character(sc_multi$matched_diagnostic)[1L]
      matched_nl <- if (!is.na(snl) && nzchar(snl)) sc_multi$matched_nl else sub$matched_nl
      matched_md <- if (!is.na(smd) && nzchar(smd)) sc_multi$matched_md else sub$matched_md
      if (!is.na(sdi) && nzchar(sdi)) {
        matched_diag <- sc_multi$matched_diagnostic
      }
    } else {
      score_ms2 <- score_cosine_ms2
      score_ndp <- NA_real_
      ndp_raw <- NA_real_
      score_nl <- sub$score_nl
      score_md <- sub$score_md
      matched_nl <- sub$matched_nl
      matched_md <- sub$matched_md
    }

    score_mz <- if (is.finite(mz_eppm) && mz_tol_ppm > 0) {
      max(0, 1 - mz_eppm / mz_tol_ppm)
    } else {
      NA_real_
    }
    score_rt <- if (is.finite(rte) && is.finite(rtl) && is.finite(rt_es)) {
      max(0, 1 - rt_es / rt_tol_s)
    } else {
      NA_real_
    }

    if (is.finite(score_mz) && is.finite(score_rt)) {
      score_total <- 0.25 * score_mz + 0.25 * score_rt + 0.5 * score_ms2
    } else if (is.finite(score_mz)) {
      score_total <- 0.4 * score_mz + 0.6 * score_ms2
    } else {
      score_total <- score_ms2
    }

    ri <- ri + 1L
    score_rows[[ri]] <- data.table::data.table(
      node_id = nid,
      Lab_ID = lid,
      db_type = if ("db_type" %in% names(cand)) {
        as.character(cand$db_type[k])
      } else {
        NA_character_
      },
      Confidence_Score = best,
      Annotation = as.character(cand$ann_name[k]),
      exp_mz = om,
      lib_mz = pm,
      mz_error_ppm = mz_eppm,
      exp_rt = rte,
      lib_rt = rtl,
      rt_error_s = rt_es,
      n_matched_peaks = best_n,
      lib_CE = best_ce,
      exp_CE = best_explab,
      best_ms2_index = best_j,
      n_ms2_subset = n_sub,
      score_total = score_total,
      score_mz = score_mz,
      score_rt = score_rt,
      score_ms2 = score_ms2,
      score_cosine_ms2 = score_cosine_ms2,
      score_nl = score_nl,
      score_md = score_md,
      score_diagnostic = score_diag,
      matched_nl = matched_nl,
      matched_md = matched_md,
      matched_diagnostic = matched_diag,
      score_ndp = score_ndp,
      ndp_raw = ndp_raw
    )
  }

  if (!ri) {
    methods::slot(project, "db_annotations") <- data.frame()
    project <- set_ms2_score_grid(project, .nts_ms2_grid_df(score_grid_rows))
    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "annotate_by_spectral_match",
      list(
        mz_tol_ppm = mz_tol_ppm,
        n_ms1_candidates = nrow(cand),
        n_scored = 0L
      )
    )
    return(invisible(diag_validate(project)))
  }

  sc_dt <- data.table::rbindlist(score_rows)
  sc_dt[, node_id := trimws(as.character(node_id))]
  if ("score_total" %in% names(sc_dt)) {
    data.table::setorder(sc_dt, node_id, -score_total, na.last = TRUE)
  } else {
    data.table::setorder(sc_dt, node_id, -Confidence_Score, na.last = TRUE)
  }

  sm_lab <- .nts_metadata_smiles_by_lab(db)
  all_sc <- if (nrow(sm_lab)) {
    sm_lab[sc_dt, on = "Lab_ID"]
  } else {
    sc_dt[, SMILES := NA_character_]
  }

  # Keep every passing compound per feature in @db_annotations (a feature can
  # legitimately match several library entries). The network nodes below keep a
  # single row per feature (highest score_total).
  ann_df <- as.data.frame(
    all_sc[, .(
      name = node_id,
      ref_id = Lab_ID,
      # Database provenance (inhouse / online / prediction) is carried through
      # so assign_confidence_levels() can tier L1 (in-house) vs L2a (online)
      # vs L2b (prediction) instead of inferring the DB type from ID patterns.
      db_type = db_type,
      Compound.name = Annotation,
      SMILES = SMILES,
      score_total = score_total,
      score_mz = score_mz,
      score_rt = score_rt,
      score_ms2 = score_ms2,
      score_cosine_ms2 = score_cosine_ms2,
      score_nl = score_nl,
      score_md = score_md,
      score_diagnostic = score_diagnostic,
      score_ndp = score_ndp,
      ndp_raw = ndp_raw,
      lib_CE = lib_CE,
      exp_CE = exp_CE,
      best_ms2_index = best_ms2_index,
      n_ms2_subset = n_ms2_subset,
      matched_nl = matched_nl,
      matched_md = matched_md,
      matched_diagnostic = matched_diagnostic,
      exp_mz = exp_mz,
      lib_mz = lib_mz,
      mz_error_ppm = mz_error_ppm,
      exp_rt = exp_rt,
      lib_rt = lib_rt,
      rt_error_s = rt_error_s,
      n_matched_peaks = n_matched_peaks
    )],
    stringsAsFactors = FALSE
  )
  methods::slot(project, "db_annotations") <- ann_df
  project <- set_ms2_score_grid(project, .nts_ms2_grid_df(score_grid_rows))

  best <- all_sc[, .SD[1L], by = node_id]
  # (0.1.28) This step previously also stashed a legacy `@network_results` stub:
  # `nodes` = the best row per feature, `edges` = an empty from/to frame. That
  # predated the molecular-network module -- real MS2 similarity graphs come from
  # build_molecular_network() / get_network_data(), never from here, so the stub
  # only ever polluted a mined project with non-network content. It was removed
  # so `@network_results` stays empty unless a caller fills it explicitly with
  # set_network_results(). `best` is still needed for the log entry below.
  project <- .nts_merge_db_top1_into_ms2_annotations(project)
  project <- .nts_merge_db_best_into_ms1_library_seed_cols(project)
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "annotate_by_spectral_match",
    list(
      mz_tol_ppm = mz_tol_ppm,
      rt_tol_s = rt_tol_s,
      ms2_tol_da = ms2_tol_da,
      min_score = min_score,
      n_ms1_candidates = nrow(cand),
      n_nodes_annotated = nrow(best)
    )
  )
  return(invisible(diag_validate(project)))
}

# ============================================================================
# MS2 diagnostic fragment validation (data.table non-equi joins)
# ============================================================================

#' Packaged default CSV of expert fragment masses (optional legacy workflow).
#' @export
default_diagnostic_fragments_csv <- function() {
  # 走内置查找链：新布局 tables/… 优先，兼容旧布局 extdata/…；找不到返回 NA_character_
  .nts_extdata(c("tables/diagnostic_fragments.csv", "diagnostic_fragments.csv"))
}


#' Run the packaged diagnostic-fragment validation as part of annotation
#'
#' Internal default-on hook used by [annotate_by_spectral_match] (and therefore
#' by [batch_annotate] / [hierarchical_annotate]). It matches the packaged
#' expert fragment table (see [default_diagnostic_fragments_csv]; PFAS-focused,
#' Koelmel et al. 2022 + Barzen-Hanson et al. 2017) against each feature's
#' experimental MS2 (`@ms2data`) and records the matched fragment formulas in
#' `@ms2_annotations$diagnostic_fragments`, bumping
#' `@ms1_annotations$annotation_level` with `[MS2 Verified]` where applicable.
#'
#' Delegates to [validate_with_ms2] so the matching logic lives in exactly one
#' place. Silently no-ops (returning `project` unchanged) when there is no MS2
#' data or the CSV path cannot be resolved, so it is always safe to call.
#'
#' @param project An [NTSProject].
#' @param csv Diagnostic CSV path (`fragment_formula`, `fragment_mass`). `NULL`
#'   or empty resolves to [default_diagnostic_fragments_csv].
#' @param mz_tol_ppm Symmetric ppm tolerance around reference fragment masses.
#' @param int_threshold Keep fragment peaks with `intensity > int_threshold`.
#' @return The updated `project`.
#' @keywords internal
#' @noRd
.nts_apply_default_diagnostic_fragments <- function(project,
                                                    csv = NULL,
                                                    mz_tol_ppm = 10,
                                                    int_threshold = 0) {
  if (!methods::is(project, "NTSProject")) {
    return(project)
  }
  if (!length(methods::slot(project, "ms2data"))) {
    return(project)
  }
  if (is.null(csv) || !length(csv) || !nzchar(as.character(csv)[1L])) {
    csv <- default_diagnostic_fragments_csv()
  }
  csv <- as.character(csv)[1L]
  if (is.na(csv) || !nzchar(csv)) {
    return(project)
  }
  validate_with_ms2(
    project,
    diagnostic_csv_path = csv,
    mz_tol_ppm = mz_tol_ppm,
    int_threshold = int_threshold
  )
}


#' @noRd
.nts_resolve_diagnostic_csv <- function(diagnostic_csv_path) {
  if (!is.null(diagnostic_csv_path) && nzchar(as.character(diagnostic_csv_path)[1L])) {
    p <- normalizePath(as.character(diagnostic_csv_path)[1L], mustWork = TRUE)
    return(p)
  }
  stop(
    "diagnostic_csv_path is empty; pass a file path or use default_diagnostic_fragments_csv()",
    call. = FALSE
  )
}


#' @noRd
.nts_sanitize_rule_bridge_col <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^A-Za-z0-9_.]+", "_", x, perl = TRUE)
  make.names(x, unique = TRUE)
}


#' @noRd
.nts_read_diagnostic_rules_csv <- function(path) {
  path <- normalizePath(path.expand(as.character(path)[1L]), mustWork = TRUE)
  dt <- data.table::fread(path, na.strings = c("", "NA"))
  if (!nrow(dt)) {
    stop("rules_file is empty", call. = FALSE)
  }
  cn <- names(dt)
  pick <- function(w) {
    k <- which(tolower(cn) == tolower(w))
    if (!length(k)) {
      stop("rules CSV must contain column `", w, "`", call. = FALSE)
    }
    cn[k[1L]]
  }
  mc <- pick("mz")
  tc <- pick("tolerance")
  yc <- pick("type")
  nc <- pick("name")
  out <- data.table::data.table(
    rule_mz = suppressWarnings(as.numeric(dt[[mc]])),
    tol_ppm = suppressWarnings(as.numeric(dt[[tc]])),
    rule_type = tolower(trimws(as.character(dt[[yc]]))),
    rule_label = trimws(as.character(dt[[nc]]))
  )
  out <- out[is.finite(rule_mz) & is.finite(tol_ppm) & tol_ppm >= 0]
  if (!nrow(out)) {
    stop("rules CSV has no rows with finite mz and tolerance", call. = FALSE)
  }
  out[, rule_type := {
    r <- rule_type
    r[r %in% c("fragment", "frag", "f")] <- "fragment"
    r[r %in% c("neutral_loss", "nl", "neutral loss", "neutral-loss")] <- "neutral_loss"
    r[r %in% c("mass_difference", "mass_diff", "md", "mass difference")] <- "mass_difference"
    r
  }]
  bad <- !out$rule_type %in% c("fragment", "neutral_loss", "mass_difference")
  if (any(bad)) {
    stop(
      "rules CSV `type` must be fragment / neutral_loss / mass_difference (aliases allowed); bad: ",
      paste(unique(out$rule_type[bad]), collapse = ", "),
      call. = FALSE
    )
  }
  out[is.na(rule_label) | !nzchar(rule_label), rule_label := sprintf("rule_%.6f", rule_mz)]
  out[, bridge_col := vapply(rule_label, .nts_sanitize_rule_bridge_col, "")]
  dup <- duplicated(out$bridge_col) | duplicated(out$bridge_col, fromLast = TRUE)
  if (any(dup)) {
    out[dup, bridge_col := make.unique(bridge_col, sep = "_")]
  }
  out[]
}


#' @noRd
.nts_rule_mass_diff_hit_count <- function(mz_v, ref, tol_ppm) {
  u <- unique(sort(mz_v[is.finite(mz_v)]))
  if (length(u) < 2L) {
    return(0L)
  }
  ppm <- tol_ppm * 1e-6
  m <- outer(u, u, function(a, b) abs(a - b))
  md <- m[lower.tri(m, diag = FALSE)]
  md <- md[md > 1e-6 & is.finite(md)]
  if (!length(md)) {
    return(0L)
  }
  lo <- ref * (1 - ppm)
  hi <- ref * (1 + ppm)
  as.integer(sum(md >= lo & md <= hi))
}


#' @noRd
.nts_rules_match_clean_ms2_long <- function(clean_ms2, rules_dt, int_threshold) {
  thr <- suppressWarnings(as.numeric(int_threshold)[1L])
  if (!is.finite(thr)) {
    thr <- 0
  }
  rows <- list()
  ri <- 0L
  if (!is.list(clean_ms2) || !length(clean_ms2)) {
    return(data.table::data.table(
      name = character(0),
      bridge_col = character(0),
      rule_label = character(0),
      n_hits = integer(0)
    ))
  }
  fnames <- names(clean_ms2)
  if (is.null(fnames)) {
    fnames <- as.character(seq_along(clean_ms2))
  }
  for (i in seq_along(clean_ms2)) {
    fn <- fnames[i]
    el <- clean_ms2[[i]]
    mz_v <- numeric(0)
    prec <- NA_real_
    if (is.list(el)) {
      prec <- suppressWarnings(as.numeric(el$precursorMZ)[1L])
      sp <- el$spectrum
      if (is.data.frame(sp) && nrow(sp) && "mz" %in% names(sp)) {
        mz_v <- as.numeric(sp$mz)
        if ("intensity" %in% names(sp)) {
          int_v <- as.numeric(sp$intensity)
          ok <- is.finite(mz_v) & (is.na(int_v) | int_v > thr)
        } else {
          ok <- is.finite(mz_v)
        }
        mz_v <- mz_v[ok]
      }
    }
    for (r in seq_len(nrow(rules_dt))) {
      rt <- rules_dt$rule_type[r]
      ref <- rules_dt$rule_mz[r]
      tol <- rules_dt$tol_ppm[r]
      ppmf <- tol * 1e-6
      nh <- 0L
      if (rt == "fragment") {
        lo <- ref * (1 - ppmf)
        hi <- ref * (1 + ppmf)
        nh <- as.integer(sum(mz_v >= lo & mz_v <= hi, na.rm = TRUE))
      } else if (rt == "neutral_loss") {
        if (is.finite(prec) && length(mz_v)) {
          loss <- prec - mz_v
          lo <- ref * (1 - ppmf)
          hi <- ref * (1 + ppmf)
          nh <- as.integer(sum(loss >= lo & loss <= hi & is.finite(loss), na.rm = TRUE))
        }
      } else if (rt == "mass_difference") {
        nh <- .nts_rule_mass_diff_hit_count(mz_v, ref, tol)
      }
      ri <- ri + 1L
      rows[[ri]] <- data.table::data.table(
        name = as.character(fn),
        bridge_col = rules_dt$bridge_col[r],
        rule_label = rules_dt$rule_label[r],
        n_hits = nh
      )
    }
  }
  data.table::rbindlist(rows, use.names = TRUE)
}


#' @noRd
.nts_rules_bridge_wide <- function(project, long_dt) {
  ft <- get_feature_table(project)
  ft_names <- if (is.data.frame(ft) && nrow(ft) && "name" %in% names(ft)) {
    unique(as.character(ft$name))
  } else {
    character(0)
  }
  if (!nrow(long_dt)) {
    out <- data.table::data.table(name = ft_names)
    return(out[])
  }
  wide <- data.table::dcast(long_dt, name ~ bridge_col, value.var = "n_hits", fill = 0L)
  alln <- unique(c(ft_names, as.character(wide$name)))
  if (!length(alln)) {
    alln <- unique(as.character(long_dt$name))
  }
  base <- data.table::data.table(name = alln)
  merged <- wide[base, on = "name"]
  old_br <- get_ms2_annotations(project)
  if (is.data.frame(old_br) && nrow(old_br) && "name" %in% names(old_br)) {
    odt <- data.table::as.data.table(old_br)
    stat_cols <- intersect(
      c("feature_id", "mz", "rt", "mzmin", "mzmax", "n_fragments"),
      names(odt)
    )
    stat_cols <- setdiff(stat_cols, names(merged))
    if (length(stat_cols)) {
      extra <- odt[, c("name", stat_cols), with = FALSE]
      merged <- extra[merged, on = "name"]
    }
  }
  rule_cols <- setdiff(names(wide), "name")
  for (nm in rule_cols) {
    if (!nm %in% names(merged)) {
      next
    }
    x <- merged[[nm]]
    x[is.na(x)] <- 0L
    merged[[nm]] <- as.integer(round(x))
  }
  merged[]
}


#' @noRd
.nts_ms2_hits_from_diagnostic_ions <- function(ms2_long, d_usr, mz_tol_ppm, prec_by_name) {
  ppm <- as.numeric(mz_tol_ppm)[1L] * 1e-6
  hits <- list()
  hi <- 0L
  frag <- d_usr[type == "fragment"]
  if (nrow(frag)) {
    diag <- data.table::data.table(
      fragment_formula = ifelse(
        is.na(frag$name) | !nzchar(frag$name),
        sprintf("%.6f", frag$mz),
        frag$name
      ),
      fragment_mass = frag$mz
    )
    diag <- diag[is.finite(fragment_mass) & nzchar(fragment_formula)]
    if (nrow(diag)) {
      diag[, `:=`(
        mz_lo = fragment_mass * (1 - ppm),
        mz_hi = fragment_mass * (1 + ppm)
      )]
      data.table::setkey(ms2_long, mz)
      h <- ms2_long[
        diag,
        on = .(mz >= mz_lo, mz <= mz_hi),
        nomatch = 0L,
        allow.cartesian = TRUE
      ]
      if (nrow(h)) {
        hi <- hi + 1L
        hits[[hi]] <- h
      }
    }
  }
  nl <- d_usr[type == "nl"]
  if (nrow(nl) && nrow(ms2_long)) {
    nm <- as.character(ms2_long$name)
    prec <- prec_by_name[nm]
    mz_obs <- ms2_long$mz
    for (j in seq_len(nrow(nl))) {
      L <- nl$mz[j]
      lab <- if (is.na(nl$name[j]) || !nzchar(nl$name[j])) {
        sprintf("nl_%.4f", L)
      } else {
        as.character(nl$name[j])
      }
      if (!is.finite(L)) {
        next
      }
      target <- prec - L
      ok <- is.finite(target) & is.finite(mz_obs) & target > 0
      tol <- abs(target) * ppm
      tol[!is.finite(tol) | tol < 1e-6] <- 1e-6
      w <- ok & (abs(mz_obs - target) <= tol)
      if (any(w)) {
        hi <- hi + 1L
        hits[[hi]] <- data.table::copy(ms2_long)[w][, fragment_formula := lab]
      }
    }
  }
  if (!length(hits)) {
    return(ms2_long[FALSE])
  }
  data.table::rbindlist(hits, use.names = TRUE, fill = TRUE)
}


#' @noRd
.nts_ms2_list_to_long_dt <- function(ms2_list, int_threshold) {
  if (!length(ms2_list)) {
    return(data.table::data.table(name = character(0), mz = numeric(0), intensity = numeric(0)))
  }
  nms <- names(ms2_list)
  if (is.null(nms)) {
    nms <- as.character(seq_along(ms2_list))
  }
  parts <- vector("list", length(ms2_list))
  for (i in seq_along(ms2_list)) {
    nm <- nms[i]
    scans <- ms2_list[[i]]
    if (!length(scans)) {
      parts[[i]] <- data.table::data.table(name = character(0), mz = numeric(0), intensity = numeric(0))
      next
    }
    parts[[i]] <- data.table::rbindlist(
      lapply(scans, function(sp) {
        mz <- as.numeric(sp$mz)
        int <- as.numeric(sp$intensity)
        data.table::data.table(mz = mz, intensity = int)
      }),
      fill = TRUE
    )
    parts[[i]][, name := as.character(nm)]
  }
  long <- data.table::rbindlist(parts, fill = TRUE)
  thr <- as.numeric(int_threshold)[1L]
  if (!is.finite(thr)) {
    thr <- 0
  }
  long <- long[is.finite(mz) & is.finite(intensity) & intensity > thr]
  long[]
}


#' @noRd
.nts_annotation_level_after_ms2 <- function(al) {
  .al <- as.character(al)
  is_na <- is.na(.al) | !nzchar(.al)
  is_unk <- is_na |
    grepl("Isolated Unknown", .al, fixed = TRUE) |
    grepl("Unknown Series", .al, fixed = TRUE)
  has_ms2 <- grepl("[MS2 Verified]", .al, fixed = TRUE)
  out <- .al
  out[is_unk] <- "Unknown candidate [MS2 Verified]"
  rest <- !is_unk & !has_ms2
  out[rest] <- paste0(trimws(.al[rest]), " [MS2 Verified]")
  out
}


#' Validate annotations using diagnostic MS/MS fragments (expert MS2 radar)
#'
#' Validate MS2 evidence against diagnostic rules
#'
#' **Rules file mode** (`rules_file`): reads an external CSV (`mz`, `tolerance`
#' in ppm, `type`, `name`), matches each row against `@clean_ms2` spectra
#' (`fragment`, `neutral_loss`, `mass_difference`), and merges a wide table of
#' sanitized rule names (integer hit counts) into `@ms2_annotations` by
#' `name`. Requires [purify_ms2_spectra] first.
#'
#' **Legacy mode** (no `rules_file`): long-formats experimental MS2
#' (`ms2_list` / `@ms2data`), matches `diagnostic_ions` or legacy diagnostic
#' CSV (`fragment_formula`, `fragment_mass`), updates `@ms2_annotations` and
#' `@ms1_annotations$annotation_level`.
#'
#' @param project An [NTSProject] with `ms1_annotations` containing `name`.
#' @param ms2_list Named list of MS2 scan lists (each scan: `mz`, `intensity`).
#'   `NULL` uses `@ms2data`. Used only when `rules_file` is `NULL`.
#' @param diagnostic_ions `NULL` (default) skips validation. Otherwise a
#'   `data.frame` with columns `mz`, `type` (`"fragment"` or `"nl"`), optional
#'   `name`. Neutral-loss rows use `feature_table$mz` per feature as precursor.
#' @param diagnostic_csv_path CSV with at least `fragment_formula`,
#'   `fragment_mass`. Used only when `diagnostic_ions` is `NULL`. Pass
#'   [default_diagnostic_fragments_csv] to use the packaged example table.
#' @param rules_file If set, path to CSV with columns `mz`, `tolerance` (ppm),
#'   `type` (`fragment` / `neutral_loss` / `mass_difference`), `name`. Enables
#'   bridge-table mode on `@clean_ms2` (legacy arguments ignored for matching).
#' @param mz_tol_ppm Symmetric ppm tolerance around reference masses (legacy
#'   path only).
#' @param int_threshold Keep fragment peaks with `intensity > int_threshold`
#'   (filters background; default `0` drops non-positive intensities).
#'
#' @return The updated `project` (always [NTSProject]; returned invisibly).
#' @export
validate_with_ms2 <- function(project,
                              ms2_list = NULL,
                              diagnostic_ions = NULL,
                              diagnostic_csv_path = NULL,
                              rules_file = NULL,
                              mz_tol_ppm = 10,
                              int_threshold = 0) {
  if (!methods::is(project, "NTSProject")) {
    stop(
      "validate_with_ms2(): `project` must be an NTSProject object (pipeline upstream may have returned a data.table instead of the project).",
      call. = FALSE
    )
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }

  rf <- rules_file
  if (!is.null(rf) && nzchar(as.character(rf)[1L])) {
    rules_dt <- .nts_read_diagnostic_rules_csv(rf)
    clean_l <- methods::slot(project, "clean_ms2")
    if (!is.list(clean_l) || !length(clean_l)) {
      stop(
        "validate_with_ms2(rules_file=): @clean_ms2 is empty; run purify_ms2_spectra() first.",
        call. = FALSE
      )
    }
    long_h <- .nts_rules_match_clean_ms2_long(clean_l, rules_dt, int_threshold)
    project <- .nts_ensure_ms2_ann_for_ms2_list(project, clean_l)
    wide_dt <- .nts_rules_bridge_wide(project, long_h)
    w_df <- as.data.frame(wide_dt, stringsAsFactors = FALSE)
    m2_df <- as.data.frame(get_ms2_annotations(project), stringsAsFactors = FALSE)
    if (!nrow(m2_df) || !"name" %in% names(m2_df)) {
      stop("ms2_annotations is empty or lacks `name`", call. = FALSE)
    }
    only_w <- setdiff(names(w_df), names(m2_df))
    if (length(only_w)) {
      m2_df <- merge(
        m2_df,
        w_df[, c("name", only_w), drop = FALSE],
        by = "name",
        all.x = TRUE,
        sort = FALSE
      )
    }
    both <- intersect(setdiff(names(w_df), "name"), names(m2_df))
    idx <- match(m2_df$name, w_df$name)
    ok <- !is.na(idx)
    for (cl in both) {
      v <- w_df[[cl]][idx]
      if (any(ok)) {
        m2_df[[cl]][ok] <- v[ok]
      }
    }
    project <- set_ms2_annotations(project, m2_df)

    if (nrow(long_h)) {
      summ <- long_h[n_hits > 0L, .(diagnostic_fragments = paste(unique(rule_label), collapse = "; ")), by = name]
    } else {
      summ <- data.table::data.table(name = character(0), diagnostic_fragments = character(0))
    }
    m2 <- data.table::as.data.table(get_ms2_annotations(project))
    m2[, name := as.character(name)]
    if (!"diagnostic_fragments" %in% names(m2)) {
      m2[, diagnostic_fragments := NA_character_]
    }
    if (nrow(summ)) {
      m2[summ, on = "name", diagnostic_fragments := i.diagnostic_fragments]
    }
    cn <- names(m2)
    if ("diagnostic_fragments" %in% cn) {
      data.table::setcolorder(m2, c(setdiff(cn, "diagnostic_fragments"), "diagnostic_fragments"))
    }
    project <- set_ms2_annotations(project, as.data.frame(m2, stringsAsFactors = FALSE))

    ms1 <- get_ms1_annotations(project)
    if (!nrow(ms1) || !"name" %in% names(ms1)) {
      stop("ms1_annotations is empty or lacks `name`", call. = FALSE)
    }
    if (!"annotation_level" %in% names(ms1)) {
      ms1$annotation_level <- NA_character_
    }
    if (nrow(summ)) {
      hit_names <- unique(as.character(summ$name))
      w <- as.character(ms1$name) %in% hit_names
      if (any(w)) {
        ms1$annotation_level[w] <- .nts_annotation_level_after_ms2(as.character(ms1$annotation_level[w]))
      }
    }
    project <- set_ms1_annotations(project, ms1)

    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "validate_with_ms2",
      list(
        mode = "rules_file",
        rules_file = normalizePath(path.expand(as.character(rf)[1L]), mustWork = TRUE),
        n_rules = nrow(rules_dt),
        int_threshold = int_threshold
      )
    )
    message(
      "> MS2 validation (rules_file): ",
      nrow(rules_dt),
      " rules x clean_ms2 -> @ms2_annotations (",
      ncol(wide_dt) - 1L,
      " count columns)."
    )
    return(invisible(project))
  }

  if (is.null(ms2_list)) {
    ms2_list <- methods::slot(project, "ms2data")
  }
  if (!length(ms2_list)) {
    message("> MS2 validation: ms2_list / ms2data is empty - skipped.")
    return(invisible(project))
  }
  project <- .nts_ensure_ms2_ann_for_ms2_list(project, ms2_list)

  ms2_long <- .nts_ms2_list_to_long_dt(ms2_list, int_threshold)
  if (!nrow(ms2_long)) {
    message("> MS2 validation: no fragment peaks after intensity filter - skipped.")
    return(invisible(project))
  }

  ppm <- as.numeric(mz_tol_ppm)[1L] * 1e-6
  if (!is.finite(ppm) || ppm < 0) {
    stop("mz_tol_ppm must be a non-negative finite number", call. = FALSE)
  }

  d_usr <- .nts_validate_diagnostic_ions(diagnostic_ions, "diagnostic_ions")
  csv_path <- NA_character_
  if (!is.null(d_usr)) {
    ft <- get_feature_table(project)
    if (!nrow(ft) || !"name" %in% names(ft) || !"mz" %in% names(ft)) {
      stop("feature_table must contain name and mz for diagnostic_ions validation", call. = FALSE)
    }
    prec_by_name <- stats::setNames(as.numeric(ft$mz), as.character(ft$name))
    hits <- .nts_ms2_hits_from_diagnostic_ions(ms2_long, d_usr, mz_tol_ppm, prec_by_name)
  } else {
    dpath <- diagnostic_csv_path
    if (is.null(dpath) || !nzchar(as.character(dpath)[1L])) {
      message("> MS2 validation: no diagnostic_ions or diagnostic_csv_path - skipped.")
      return(invisible(project))
    }
    csv_path <- .nts_resolve_diagnostic_csv(dpath)
    diag <- data.table::fread(csv_path, na.strings = c("", "NA"))
    req <- c("fragment_formula", "fragment_mass")
    if (!all(req %in% names(diag))) {
      stop("diagnostic CSV must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
    }
    diag[, fragment_formula := trimws(as.character(fragment_formula))]
    diag[, fragment_mass := as.numeric(fragment_mass)]
    diag <- diag[is.finite(fragment_mass) & nzchar(fragment_formula)]
    if (!nrow(diag)) {
      stop("diagnostic CSV has no usable fragment rows", call. = FALSE)
    }
    diag[, `:=`(
      mz_lo = fragment_mass * (1 - ppm),
      mz_hi = fragment_mass * (1 + ppm)
    )]
    data.table::setkey(ms2_long, mz)
    hits <- ms2_long[
      diag,
      on = .(mz >= mz_lo, mz <= mz_hi),
      nomatch = 0L,
      allow.cartesian = TRUE
    ]
  }

  if (!nrow(hits)) {
    message("> MS2 validation: 0 diagnostic hits (", mz_tol_ppm, " ppm).")
    m2 <- data.table::as.data.table(get_ms2_annotations(project))
    if (!"diagnostic_fragments" %in% names(m2)) {
      m2[, diagnostic_fragments := NA_character_]
    }
    project <- set_ms2_annotations(project, as.data.frame(m2, stringsAsFactors = FALSE))
    return(invisible(project))
  }

  agg <- hits[
    ,
    .(diagnostic_fragments = paste(unique(fragment_formula), collapse = "; ")),
    by = name
  ]

  m2 <- data.table::as.data.table(get_ms2_annotations(project))
  if (!nrow(m2) || !"name" %in% names(m2)) {
    stop("ms2_annotations is empty or lacks `name`", call. = FALSE)
  }
  m2[, name := as.character(name)]

  if (!"diagnostic_fragments" %in% names(m2)) {
    m2[, diagnostic_fragments := NA_character_]
  }
  m2[agg, on = "name", diagnostic_fragments := i.diagnostic_fragments]

  cn <- names(m2)
  if ("diagnostic_fragments" %in% cn) {
    data.table::setcolorder(m2, c(setdiff(cn, "diagnostic_fragments"), "diagnostic_fragments"))
  }

  project <- set_ms2_annotations(project, as.data.frame(m2, stringsAsFactors = FALSE))

  ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
  if (!nrow(ms1) || !"name" %in% names(ms1)) {
    # No per-feature MS1 layer to update: the MS2-level `diagnostic_fragments`
    # column written above is still kept; skip the annotation_level bump so the
    # default-on hook stays safe on bare projects.
    message(
      "> MS2 validation: ms1_annotations is empty or lacks `name` - ",
      "annotation_level not updated."
    )
    return(invisible(project))
  }
  if (!"annotation_level" %in% names(ms1)) {
    ms1$annotation_level <- NA_character_
  }
  hit_names <- unique(as.character(agg$name))
  w <- as.character(ms1$name) %in% hit_names
  if (any(w)) {
    ms1$annotation_level[w] <- .nts_annotation_level_after_ms2(as.character(ms1$annotation_level[w]))
  }
  project <- set_ms1_annotations(project, ms1)

  n_feat <- nrow(agg)
  n_hit <- nrow(hits)
  message(
    "> MS2 validation: ",
    n_hit,
    " fragment-ion hits -> ",
    n_feat,
    " feature(s) with diagnostic_fragments; annotation_level updated where applicable."
  )

  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "validate_with_ms2",
    list(
      n_hits = n_hit,
      n_features_tagged = n_feat,
      mz_tol_ppm = mz_tol_ppm,
      diagnostic_csv = if (nzchar(as.character(csv_path)[1L]) && !is.na(csv_path)) csv_path else NA_character_
    )
  )
  return(invisible(project))
}

#' Hierarchical Batch Annotation across Multiple Databases
#'
#' Sequentially annotates features using a list of databases (e.g., L1, L2a, L2b).
#' Features successfully annotated at a higher confidence level are excluded from
#' searches against lower-level databases. The final `db_annotations` table 
#' includes a `Confidence_Level` column indicating the source database, which
#' is also propagated to the `annotation_level` column in `ms1_annotations`.
#'
#' @param project An \code{NTSProject}.
#' @param databases A named list of database objects (e.g., \code{list(L1 = db_inhouse, L2a = db_public)}).
#' @param mz_tol Numeric vector of m/z tolerances in ppm. If length 1, recycled for all databases.
#' @param rt_tol Numeric vector of RT tolerances in seconds. If length 1, recycled. Use high values (e.g., 9999) to ignore RT.
#' @param feature_tol Numeric vector of MS2 feature matching tolerances in Da. If length 1, recycled.
#' @param min_score Numeric vector of minimum MS2 composite scores. If length 1, recycled.
#' @param polarity Ionisation polarity (\code{"neg"} or \code{"pos"}).
#' @param ... Additional arguments passed to \code{\link{batch_annotate}}.
#'
#' @return An updated \code{NTSProject} with combined hierarchical annotations.
#' @export
hierarchical_annotate <- function(project,
                                  databases,
                                  mz_tol = 10,
                                  rt_tol = 150,
                                  feature_tol = 0.02,
                                  min_score = 0.3,
                                  polarity = c("neg", "pos"),
                                  ...) {
  polarity <- match.arg(polarity)
  
  if (!is.list(databases) || is.null(names(databases))) {
    stop("`databases` must be a named list of database objects (e.g., list(L1=db1, L2a=db2)).", call. = FALSE)
  }
  
  n_db <- length(databases)
  db_names <- names(databases)
  
  expand_param <- function(p) {
    if (length(p) == 1) rep(p, n_db) else p
  }
  mz_tols <- expand_param(mz_tol)
  rt_tols <- expand_param(rt_tol)
  feat_tols <- expand_param(feature_tol)
  min_scores <- expand_param(min_score)
  
  all_anns <- list()
  ft_full <- get_feature_table(project)
  if (!nrow(ft_full)) {
    return(project)
  }
  
  unannotated_names <- ft_full$name
  
  for (i in seq_len(n_db)) {
    lvl <- db_names[i]
    db <- databases[[i]]
    if (length(unannotated_names) == 0) break
    
    message(sprintf("\n=== Annotating Level: %s (%d features pending) ===", lvl, length(unannotated_names)))
    
    temp_proj <- copy_project(project)
    temp_proj <- set_feature_table(temp_proj, ft_full[ft_full$name %in% unannotated_names, , drop = FALSE])
    
    temp_proj <- batch_annotate(
      temp_proj,
      database = db,
      mz_tol = mz_tols[i],
      rt_tol = rt_tols[i],
      feature_tol = feat_tols[i],
      min_score = min_scores[i],
      polarity = polarity,
      ...
    )
    
    ann <- temp_proj@db_annotations
    if (is.data.frame(ann) && nrow(ann) > 0) {
      ann$Confidence_Level <- lvl
      all_anns[[length(all_anns) + 1L]] <- ann
      
      matched_names <- unique(as.character(ann$name))
      unannotated_names <- setdiff(unannotated_names, matched_names)
    }
  }
  
  message("\n=== Hierarchical Annotation Complete ===")
  if (length(all_anns) > 0) {
    master_ann <- do.call(rbind, all_anns)
    methods::slot(project, "db_annotations") <- master_ann
  } else {
    methods::slot(project, "db_annotations") <- data.frame()
  }
  
  project <- .nts_merge_db_top1_into_ms2_annotations(project)
  project <- .nts_merge_db_best_into_ms1_library_seed_cols(project)
  
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "hierarchical_annotate",
    list(levels = db_names, polarity = polarity)
  )
  
  return(invisible(project))
}
