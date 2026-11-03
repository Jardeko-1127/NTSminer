# ---------------------------------------------------------------------------
# DDA MS2 QC scoring (P4 of the DDA optimisation roadmap)
#
# `compute_ms2_qc()` decorates `@ms2_annotations` with per-feature quality
# metrics derived from `@clean_ms2` (consensus spectra) and `@ms2data` (raw
# scans + match scores written by [extract_ms2(scoring = ...)]).
# QC values are descriptive only by default -- downstream scoring is not
# auto-weighted (users can do so manually via the QC columns).
# ---------------------------------------------------------------------------

#' Default thresholds for [compute_ms2_qc()] grade assignment
#'
#' Returns a nested list defining when a feature is graded `"A"` (high) /
#' `"B"` (medium) / `"C"` (low). All A-criteria must pass; if not, any
#' B-criterion is sufficient for B; otherwise C.
#'
#' @param a_min_raw_scans Minimum `n_raw_scans` for grade A (default 3).
#' @param a_min_match_score_median Minimum median DDA `match_score` across
#'   scans for grade A (default 0.6). Requires
#'   `extract_ms2(scoring = "soft")` or `"strict"`.
#' @param a_min_spectral_entropy Minimum Shannon entropy of consensus
#'   `spectrum` (default 2.0, natural log units).
#' @param a_min_mean_purity Minimum mean isolation purity (`s_purity`) for
#'   grade A (default 0.7). `s_purity` is populated only when
#'   [extract_ms2] was run with `scoring = "soft"`/"strict" and the file exposes
#'   MS1 scans; when it is `NA` the criterion is treated as not applicable.
#' @param b_min_raw_scans Minimum `n_raw_scans` for grade B (default 1).
#' @param b_min_match_score_max Minimum max DDA `match_score` for grade B
#'   (default 0.4).
#'
#' @return Named list `list(a = list(...), b = list(...))`.
#' @export
ms2_qc_thresholds <- function(a_min_raw_scans = 3L,
                              a_min_match_score_median = 0.6,
                              a_min_spectral_entropy = 2.0,
                              a_min_mean_purity = 0.7,
                              b_min_raw_scans = 1L,
                              b_min_match_score_max = 0.4) {
  list(
    a = list(
      min_raw_scans = as.integer(a_min_raw_scans),
      min_match_score_median = as.numeric(a_min_match_score_median),
      min_spectral_entropy = as.numeric(a_min_spectral_entropy),
      min_mean_purity = as.numeric(a_min_mean_purity)
    ),
    b = list(
      min_raw_scans = as.integer(b_min_raw_scans),
      min_match_score_max = as.numeric(b_min_match_score_max)
    )
  )
}


#' Decorate `@ms2_annotations` with DDA QC metrics
#'
#' Computes per-feature quality indicators on top of the consensus spectra in
#' `@clean_ms2` and the raw scan list in `@ms2data`. Added columns:
#'
#' * `spectral_entropy` - Shannon entropy of normalised consensus
#'   `intensity` (natural log; higher = more dispersed fragment pattern).
#' * `top_peak_fraction` - `max(intensity) / sum(intensity)` of consensus.
#' * `inter_scan_similarity` - median pairwise cosine similarity between raw
#'   scans (only when `>= 2` scans exist); robust replicability proxy.
#' * `match_score_max` / `match_score_median` - propagated from per-scan
#'   `match_score` written by `extract_ms2(scoring = "soft"|"strict")`.
#'   `NA` if scoring was disabled.
#' * `mean_purity` - mean of per-scan `s_purity` (isolation purity). `NA` unless
#'   [extract_ms2] was run with scoring enabled and MS1 scans are available.
#' * `qc_grade` - `"A"` / `"B"` / `"C"` using `thresholds`.
#'
#' @param project An [NTSProject] with a non-empty `@clean_ms2`.
#' @param thresholds Output of [ms2_qc_thresholds()].
#' @param include_reproducibility Compute `inter_scan_similarity` (default
#'   `TRUE`). Set to `FALSE` for very large studies to save time.
#'
#' @return Updated `project` with the QC columns appended to
#'   `@ms2_annotations`.
#' @export
compute_ms2_qc <- function(project,
                           thresholds = ms2_qc_thresholds(),
                           include_reproducibility = TRUE) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  cm2 <- methods::slot(project, "clean_ms2")
  if (!length(cm2) || is.null(names(cm2))) {
    warning(
      "clean_ms2 is empty; run purify_ms2_spectra() first.",
      call. = FALSE
    )
    return(project)
  }
  m2d <- methods::slot(project, "ms2data")
  if (!is.list(m2d)) {
    m2d <- list()
  }

  feat_names <- names(cm2)
  n_feat <- length(feat_names)
  ent <- numeric(n_feat)
  tpf <- numeric(n_feat)
  iss <- rep(NA_real_, n_feat)
  ms_max <- rep(NA_real_, n_feat)
  ms_med <- rep(NA_real_, n_feat)
  mp <- rep(NA_real_, n_feat)
  grade <- character(n_feat)

  for (k in seq_len(n_feat)) {
    fn <- feat_names[k]
    el <- cm2[[fn]]
    sp <- el$spectrum
    ent[k] <- .nts_ms2_qc_spectral_entropy(sp)
    tpf[k] <- .nts_ms2_qc_top_peak_fraction(sp)

    scans <- m2d[[fn]]
    if (is.list(scans) && length(scans)) {
      mscores <- vapply(scans, function(s) {
        v <- suppressWarnings(as.numeric(s$match_score)[1L])
        if (length(v)) v else NA_real_
      }, numeric(1L))
      if (any(is.finite(mscores))) {
        ms_max[k] <- max(mscores, na.rm = TRUE)
        ms_med[k] <- stats::median(mscores, na.rm = TRUE)
      }
      pur <- vapply(scans, function(s) {
        v <- suppressWarnings(as.numeric(s$s_purity)[1L])
        if (length(v)) v else NA_real_
      }, numeric(1L))
      if (any(is.finite(pur))) {
        mp[k] <- mean(pur, na.rm = TRUE)
      }
      if (isTRUE(include_reproducibility) && length(scans) >= 2L) {
        iss[k] <- .nts_ms2_qc_inter_scan_similarity(scans)
      }
    }

    n_raw_val <- if (is.null(el$n_raw_scans)) {
      NA_integer_
    } else {
      suppressWarnings(as.integer(el$n_raw_scans)[1L])
    }
    grade[k] <- .nts_ms2_qc_assign_grade(
      n_raw = n_raw_val,
      ms_med = ms_med[k],
      ms_max = ms_max[k],
      entropy = ent[k],
      mean_purity = mp[k],
      th = thresholds
    )
  }

  qc_dt <- data.table::data.table(
    name = feat_names,
    spectral_entropy = ent,
    top_peak_fraction = tpf,
    inter_scan_similarity = iss,
    match_score_max = ms_max,
    match_score_median = ms_med,
    mean_purity = mp,
    qc_grade = grade
  )

  m2 <- data.table::as.data.table(get_ms2_annotations(project))
  if (!nrow(m2) || !"name" %in% names(m2)) {
    warning(
      "ms2_annotations is empty; QC metrics not attached.",
      call. = FALSE
    )
    return(project)
  }
  m2[, name := as.character(name)]
  qc_cols <- setdiff(names(qc_dt), "name")
  for (cn in qc_cols) {
    if (!cn %in% names(m2)) {
      proto <- qc_dt[[cn]]
      na_val <- if (is.character(proto)) {
        NA_character_
      } else if (is.integer(proto)) {
        NA_integer_
      } else if (is.logical(proto)) {
        NA
      } else {
        NA_real_
      }
      m2[, (cn) := na_val]
    }
  }
  idx <- match(m2$name, qc_dt$name)
  ok <- !is.na(idx)
  for (cn in qc_cols) {
    m2[ok, (cn) := qc_dt[[cn]][idx[ok]]]
  }

  project <- set_ms2_annotations(
    project,
    as.data.frame(m2, stringsAsFactors = FALSE)
  )
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "compute_ms2_qc",
    list(
      thresholds = thresholds,
      include_reproducibility = isTRUE(include_reproducibility),
      n_features = n_feat,
      grade_counts = as.list(table(factor(grade, levels = c("A", "B", "C"))))
    )
  )
  project
}


# Shannon entropy of normalised intensity vector (natural log).
# Returns 0 for empty / single-peak / non-finite inputs.
#' @noRd
.nts_ms2_qc_spectral_entropy <- function(sp) {
  if (!is.data.frame(sp) || !nrow(sp) || !"intensity" %in% names(sp)) {
    return(0)
  }
  v <- suppressWarnings(as.numeric(sp$intensity))
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) {
    return(0)
  }
  p <- v / sum(v)
  p <- p[p > 0]
  if (!length(p)) {
    return(0)
  }
  -sum(p * log(p))
}


# Top-peak intensity fraction; returns 1 for single-peak spectra, NA for empty.
#' @noRd
.nts_ms2_qc_top_peak_fraction <- function(sp) {
  if (!is.data.frame(sp) || !nrow(sp) || !"intensity" %in% names(sp)) {
    return(NA_real_)
  }
  v <- suppressWarnings(as.numeric(sp$intensity))
  v <- v[is.finite(v) & v > 0]
  if (!length(v)) {
    return(NA_real_)
  }
  max(v) / sum(v)
}


# Median pairwise cosine similarity between raw scans of one feature.
# Bins m/z at 0.01 Da to align ions before computing dot product.
#' @noRd
.nts_ms2_qc_inter_scan_similarity <- function(scans, mz_bin_da = 0.01) {
  n <- length(scans)
  if (n < 2L) {
    return(NA_real_)
  }
  vecs <- lapply(scans, function(s) {
    mz <- suppressWarnings(as.numeric(s$mz))
    int <- suppressWarnings(as.numeric(s$intensity))
    ok <- is.finite(mz) & is.finite(int) & int > 0
    if (!any(ok)) {
      return(NULL)
    }
    mz <- mz[ok]
    int <- int[ok]
    key <- round(mz / mz_bin_da) * mz_bin_da
    agg <- tapply(int, key, sum)
    list(keys = as.numeric(names(agg)), int = as.numeric(agg))
  })
  vecs <- vecs[!vapply(vecs, is.null, logical(1L))]
  if (length(vecs) < 2L) {
    return(NA_real_)
  }
  pairs <- utils::combn(seq_along(vecs), 2L, simplify = FALSE)
  sims <- vapply(pairs, function(p) {
    a <- vecs[[p[1L]]]
    b <- vecs[[p[2L]]]
    all_keys <- sort(unique(c(a$keys, b$keys)))
    va <- rep(0, length(all_keys))
    vb <- rep(0, length(all_keys))
    va[match(a$keys, all_keys)] <- a$int
    vb[match(b$keys, all_keys)] <- b$int
    na <- sqrt(sum(va * va))
    nb <- sqrt(sum(vb * vb))
    if (!is.finite(na) || !is.finite(nb) || na == 0 || nb == 0) {
      return(NA_real_)
    }
    sum(va * vb) / (na * nb)
  }, numeric(1L))
  sims <- sims[is.finite(sims)]
  if (!length(sims)) {
    return(NA_real_)
  }
  stats::median(sims)
}


# Assign A / B / C grade using thresholds; ignores NA fields gracefully.
#' @noRd
.nts_ms2_qc_assign_grade <- function(n_raw, ms_med, ms_max,
                                     entropy, mean_purity, th) {
  a <- th$a
  b <- th$b

  # Helper: returns TRUE iff the value passes (NA counts as fail).
  pass <- function(x, thr, op = ">=") {
    if (!length(thr) || !is.finite(thr)) {
      return(TRUE)
    }
    if (!length(x) || !is.finite(x)) {
      return(FALSE)
    }
    if (identical(op, ">=")) x >= thr else x > thr
  }
  # When match_score / purity weren't computed (scoring disabled), treat the
  # corresponding A-thresholds as "not applicable" rather than "fail" so
  # users running the legacy pipeline can still earn A on raw scans + entropy.
  ms_med_ok <- if (is.finite(ms_med)) pass(ms_med, a$min_match_score_median) else TRUE
  ms_max_ok_for_b <- if (is.finite(ms_max)) pass(ms_max, b$min_match_score_max) else TRUE
  purity_ok <- if (is.finite(mean_purity)) pass(mean_purity, a$min_mean_purity) else TRUE

  is_a <- pass(n_raw, a$min_raw_scans) &&
    ms_med_ok &&
    pass(entropy, a$min_spectral_entropy) &&
    purity_ok
  if (is_a) {
    return("A")
  }
  is_b <- pass(n_raw, b$min_raw_scans) && ms_max_ok_for_b
  if (is_b) {
    return("B")
  }
  "C"
}

