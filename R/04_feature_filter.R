# MS1-stage exports follow AGENTS "typical call order": optional [filter_blanks] on the
# feature-intensity matrix, then KMD (`calculate_kmd` / `filter_by_kmd`),
# heteroatom/defect helpers, isotope/ISF flags, etc.

#' Blank-aware row filter using matrix broadcasting
#'
#' Compares sample intensities to pooled blank means with a fold-change
#' threshold. **No** explicit loop over sample columns: matrices are built once
#' and compared via [base::sweep] (vectorized over rows).
#'
#' @param feature_table `data.frame` or `data.table` with numeric intensity
#'   columns.
#' @param sample_cols Character vector of column names for study samples.
#' @param blank_cols Character vector of column names for blanks / QCs used as
#'   background.
#' @param threshold Numeric multiplier: a cell must satisfy
#'   `sample > threshold * blank_mean` to count as signal (default `3`).
#'   Rows where **all** blank intensities are `NA` use `blank_mean = 0` (no
#'   detectable solvent background); any sample with finite intensity `> 0`
#'   can then satisfy the rule.
#' @param mode `"any"` keeps a row if **any** sample passes; `"all"` requires
#'   **all** samples to pass.
#'
#' @return Subset of `feature_table` (same class as input) with rows passing
#'   the blank rule.
#' @export
filter_blanks <- function(feature_table,
                          sample_cols,
                          blank_cols,
                          threshold = 3,
                          mode = c("any", "all")) {
  mode <- match.arg(mode)
  stopifnot(
    length(sample_cols) >= 1L,
    length(blank_cols) >= 1L
  )
  is_dt <- requireNamespace("data.table", quietly = TRUE) &&
    data.table::is.data.table(feature_table)
  if (is_dt) {
    sm <- as.matrix(feature_table[, sample_cols, with = FALSE])
    bm <- as.matrix(feature_table[, blank_cols, with = FALSE])
  } else {
    sm <- as.matrix(feature_table[, sample_cols, drop = FALSE])
    bm <- as.matrix(feature_table[, blank_cols, drop = FALSE])
  }
  if (!is.numeric(sm) || !is.numeric(bm)) {
    stop("sample and blank columns must be numeric", call. = FALSE)
  }
  blank_mean <- rowMeans(bm, na.rm = TRUE)
  # All-blank-NA rows: no ME signal -> treat background as 0 (do not drop feature)
  blank_mean[!is.finite(blank_mean)] <- 0
  thr_row <- threshold * blank_mean
  pass <- sweep(sm, 1L, thr_row, FUN = function(s, t) {
    ifelse(is.na(s) | is.na(t), FALSE, s > t)
  })
  if (mode == "any") {
    ok <- rowSums(pass, na.rm = TRUE) >= 1L
  } else {
    ok <- rowSums(!pass, na.rm = TRUE) == 0L & rowSums(!is.na(pass)) > 0L
  }
  feature_table[ok, , drop = FALSE]
}

#' @title KMD-based homologue filtering
#'
#' @description
#' For an [NTSProject]: calls [calculate_kmd] first to populate `KM`/`KMD`/`kmd_series_id`
#' (using `base_formula`, `max_gap`, `kmd_tol`, `mz_tol_da`, and allowed arguments in
#' `...`), then trims satellite annotation tables by Kendrick `KMD` windows (**without**
#' dropping rows from `@feature_table`). For `data.frame`/`data.table`: applies KMD
#' scaling and window filtering only (via `...$kmd_ranges` and `...$base_mz`).
#'
#' @param project An [NTSProject], or a `data.frame`/`data.table` containing `mz`.
#' @param base_formula Repeating-unit formula (default `"CF2"`).
#' @param max_gap Passed to [calculate_kmd] (DAG gap tolerance for missing homologue levels).
#' @param kmd_tol Passed to [calculate_kmd] (KMD sub-lane width).
#' @param mz_tol_da If `> 0`, DAG edges use a Da residual cutoff; `0` uses `mz_tol_ppm` only
#'   (optionally supplied in `...`).
#' @param ... `kmd_ranges` (default PFAS dual windows); other names forwarded to [calculate_kmd]
#'   (e.g. `md_range`, `mz_tol_ppm`). The table path requires `kmd_ranges`; `base_mz` is optional.
#'
#' @return Updated `project`, or a subset table with `Kendrick_Mass`/`Nominal_KM`/`KMD`.
#'
#' @export
filter_by_kmd <- function(project, base_formula = "CF2", max_gap = 2, kmd_tol = 0.005, mz_tol_da = 0.01, ...) {
  if (inherits(project, "NTSProject")) {
    return(.filter_by_kmd_ntsproject_dag(
      project,
      base_formula = base_formula,
      max_gap = max_gap,
      kmd_tol = kmd_tol,
      mz_tol_da = mz_tol_da,
      ...
    ))
  }
  dots <- list(...)
  kmd_ranges <- dots$kmd_ranges
  if (is.null(kmd_ranges)) {
    # KMD = round(kmass) - kmass is confined to [-0.5, 0.5], so a window written
    # in the [0, 1) KMD convention must be mirrored into the negative half to be
    # reachable; together the two windows cover |KMD| <= 0.15, matching the band
    # calculate_kmd() uses by default (`md_range = c(-0.15, 0.15)`).
    kmd_ranges <- list(c(0, 0.15), c(-0.15, 0))
  }
  base_mz <- dots$base_mz
  if (is.null(base_mz)) {
    base_mz <- tryCatch(
      calculate_exact_mass(as.character(base_formula)[1L]),
      error = function(e) NA_real_
    )
    if (!is.finite(base_mz)) {
      stop(sprintf("Could not parse formula or compute exact mass for: %s", base_formula), call. = FALSE)
    }
  }
  .filter_by_kmd_df_window(project, base_mz = base_mz, kmd_ranges = kmd_ranges)
}


#' @noRd
.filter_by_kmd_ntsproject_dag <- function(project, base_formula, max_gap, kmd_tol, mz_tol_da, ...) {
  dots <- list(...)
  kmd_ranges <- dots$kmd_ranges
  dots$kmd_ranges <- NULL
  base_mz_override <- if ("base_mz" %in% names(dots)) dots$base_mz else NULL
  dots$base_mz <- NULL
  calc_nm <- c(
    "md_range", "mz_tol_ppm", "rt_monotonic", "rt_max_step",
    "min_series_size", "min_length", "rt_filter", "min_spearman",
    "min_pearson_r", "min_r2", "precision", "series_digits", "series_min_n"
  )
  calc_dots <- dots[intersect(names(dots), calc_nm)]
  project <- do.call(
    calculate_kmd,
    c(
      list(
        object = project,
        base_formulas = as.character(base_formula)[1L],
        kmd_tol = kmd_tol,
        max_gap = max_gap,
        mz_tol_da = mz_tol_da
      ),
      calc_dots
    )
  )
  base_mz <- suppressWarnings(as.numeric(base_mz_override)[1L])
  if (!is.finite(base_mz) || base_mz <= 0) {
    base_mz <- tryCatch(
      calculate_exact_mass(as.character(base_formula)[1L]),
      error = function(e) NA_real_
    )
  }
  if (!is.finite(base_mz)) {
    stop(sprintf("Could not parse formula: %s", base_formula), call. = FALSE)
  }
  if (is.null(kmd_ranges)) {
    # KMD = round(kmass) - kmass is confined to [-0.5, 0.5], so a window written
    # in the [0, 1) KMD convention must be mirrored into the negative half to be
    # reachable; together the two windows cover |KMD| <= 0.15, matching the band
    # calculate_kmd() uses by default (`md_range = c(-0.15, 0.15)`).
    kmd_ranges <- list(c(0, 0.15), c(-0.15, 0))
  }
  .filter_by_kmd_ntsproject(project, base_mz = base_mz, kmd_ranges = kmd_ranges)
}


#' @noRd
.filter_by_kmd_df_window <- function(features, base_mz, kmd_ranges) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  bm <- suppressWarnings(as.numeric(base_mz)[1L])
  if (!is.finite(bm) || bm <= 0) {
    stop("`base_mz` must be a positive finite number", call. = FALSE)
  }
  dt <- data.table::as.data.table(features)
  if (!nrow(dt)) {
    dt[, Kendrick_Mass := numeric()]
    dt[, Nominal_KM := numeric()]
    dt[, KMD := numeric()]
    return(dt[])
  }

  mz_col <- NULL
  for (cand in c("mz", "MZ", "mz_value")) {
    if (cand %in% names(dt)) {
      mz_col <- cand
      break
    }
  }
  if (is.null(mz_col)) {
    stop(
      "`features` must contain one of: mz, MZ, mz_value",
      call. = FALSE
    )
  }

  nm_int <- round(bm)
  scale <- nm_int / bm

  mz_v <- suppressWarnings(as.numeric(dt[[mz_col]]))
  if (any(!is.finite(mz_v))) {
    warning("Dropping rows with non-finite m/z", call. = FALSE)
    dt <- dt[is.finite(mz_v)]
    mz_v <- mz_v[is.finite(mz_v)]
  }
  if (!nrow(dt)) {
    dt[, Kendrick_Mass := numeric()]
    dt[, Nominal_KM := numeric()]
    dt[, KMD := numeric()]
    return(dt[])
  }

  km <- mz_v * scale
  nkm <- round(km)
  kmd <- nkm - km

  dt[, Kendrick_Mass := km]
  dt[, Nominal_KM := as.numeric(nkm)]
  dt[, KMD := kmd]

  if (!is.list(kmd_ranges) || !length(kmd_ranges)) {
    stop("`kmd_ranges` must be a non-empty list of c(lo, hi) vectors", call. = FALSE)
  }

  keep <- rep(FALSE, nrow(dt))
  for (rng in kmd_ranges) {
    r <- suppressWarnings(as.numeric(rng))
    if (length(r) < 2L || !all(is.finite(r[1:2]))) {
      next
    }
    lo <- min(r[1L], r[2L])
    hi <- max(r[1L], r[2L])
    keep <- keep | (dt$KMD >= lo & dt$KMD <= hi)
  }

  dt[keep]
}


#' @noRd
.filter_by_kmd_ntsproject <- function(project, base_mz, kmd_ranges) {
  ft <- get_feature_table(project)
  if (!is.data.frame(ft) || !nrow(ft)) {
    stop("feature_table is empty", call. = FALSE)
  }
  if (!"name" %in% names(ft)) {
    stop("feature_table must contain column name", call. = FALSE)
  }
  res <- .filter_by_kmd_df_window(
    as.data.frame(ft, stringsAsFactors = FALSE),
    base_mz = base_mz,
    kmd_ranges = kmd_ranges
  )
  if (!nrow(res) || !"name" %in% names(res)) {
    kept <- character(0L)
  } else {
    kept <- unique(as.character(res$name))
    kept <- kept[!is.na(kept) & nzchar(kept)]
  }
  ms1 <- get_ms1_annotations(project)
  if (is.data.frame(ms1) && nrow(ms1) && "name" %in% names(ms1)) {
    ms1 <- ms1[as.character(ms1$name) %in% kept, , drop = FALSE]
    idx <- match(ms1$name, res$name)
    if ("Kendrick_Mass" %in% names(res)) {
      ms1$Kendrick_Mass <- ifelse(is.na(idx), NA_real_, as.numeric(res$Kendrick_Mass[idx]))
    }
    if ("Nominal_KM" %in% names(res)) {
      ms1$Nominal_KM <- ifelse(is.na(idx), NA_real_, as.numeric(res$Nominal_KM[idx]))
    }
    if ("KMD" %in% names(res)) {
      ms1$KMD <- ifelse(is.na(idx), NA_real_, as.numeric(res$KMD[idx]))
    }
    project <- set_ms1_annotations(project, .nts_ms1_sanitize_and_order(ms1))
  }
  project <- .nts_subset_project_satellites_by_name(project, kept)
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "filter_by_kmd",
    list(
      base_mz = base_mz,
      kmd_ranges = kmd_ranges,
      n_kept = length(kept),
      n_feature_table = nrow(ft)
    )
  )
  project
}


#' Mass defect class (MDC) filter (placeholder)
#'
#' Reserved for a future implementation that filters features by mass defect
#' class or related MS1 metrics. Not yet implemented.
#'
#' @param features A `data.table` or `data.frame` (passed through once implemented).
#' @param ... Reserved for future parameters (e.g. `mdc_ranges`, `polarity`).
#'
#' @return Not called; always errors.
#'
#' @export
filter_by_mdc <- function(features, ...) {
  stop(
    "filter_by_mdc() is not implemented yet; reserved for MDC-based screening.",
    call. = FALSE
  )
}

# ============================================================================
# MS1 13C isotope labelling (pre-deisotoping)
# ============================================================================

# Precise isotope mass shifts (neutral mono-isotopic substitution, Da)
.C13_MZ_DELTA <- 1.003355
.ISOTOPE_SHIFT_C13   <- 1.003355  # 13C <- 12C
.ISOTOPE_SHIFT_N15   <- 0.997035  # 15N <- 14N
.ISOTOPE_SHIFT_O18   <- 2.004246  # 18O <- 16O
.ISOTOPE_SHIFT_S34   <- 1.995796  # 34S <- 32S
.ISOTOPE_SHIFT_CL37  <- 1.997050  # 37Cl <- 35Cl
.ISOTOPE_SHIFT_BR81  <- 1.997953  # 81Br <- 79Br

# Default halogen/S isotope diffs with labels (names = isotope labels)
.ISOTOPE_HALOGEN_DEFAULTS <- c(
  "37Cl/35Cl" = 1.997050,   # 37Cl:35Cl ~ 0.32
  "81Br/79Br" = 1.997953,   # 81Br:79Br ~ 0.98
  "34S/32S"   = 1.995796    # 34S:32S ~ 0.045
)

# Minor (N/O) isotope diffs - OFF by default; opt in via
# `include_minor_isotopes = TRUE` (see flag_isotopes). 18O (2.004246) sits only
# 0.0025 Da from 13C+2 (2.006710), so it is unresolvable below ~m/z 250 at 10 ppm;
# keeping it opt-in avoids injecting that ambiguity into default runs.
.ISOTOPE_MINOR_DEFAULTS <- c(
  "15N/14N" = 0.997035,
  "18O/16O" = 2.004246
)
.ISOTOPE_MINOR_ELEMENT_BY_SHIFT <- c("0.997035" = "N", "2.004246" = "O")

# Intensity ratio validation ranges for 13C.  The upper bound is intentionally
# loose (a genuine 13C+1 peak never rivals its parent); the effective guard is
# the `cint < pint * 1.5` cap plus the 13C chain-decay rule.
.ISOTOPE_C13_RATIO_MIN <- 0.005
# Upper bound retained at 0.90: a 13C M+1/M+0 ratio above ~0.9 would require
# ~80+ carbons (natural 13C ~ 1.1% per C); typical PFAS (C4-C14) sit at
# 0.04-0.15.  An observed ratio near 0.9 therefore signals a non-monoisotopic
# "parent" (co-eluting compound / in-source fragment), so this bound is a
# useful specificity guard rather than a clip on real signal.
.ISOTOPE_C13_RATIO_MAX <- 0.90

# Exact natural abundances of the two most abundant isotopes of Cl/Br/S/N/O
# (IUPAC/CIAAW representative abundances).  Used to *compute* theoretical
# isotope envelopes instead of hard-coding intensity-ratio bands.
.ISOTOPE_ABUNDANCE <- list(
  Cl = c(major = 0.7576, minor = 0.2424),   # 35Cl / 37Cl
  Br = c(major = 0.5069, minor = 0.4931),   # 79Br / 81Br
  S  = c(major = 0.9499, minor = 0.0425),   # 32S / 34S
  N  = c(major = 0.99632, minor = 0.00368), # 14N / 15N
  O  = c(major = 0.99757, minor = 0.00205)  # 16O / 18O
)
# Highest isotope count considered when fitting an M+2 envelope.
.ISOTOPE_HALOGEN_MAX_N <- list(Cl = 4L, Br = 4L, S = 2L, N = 4L, O = 4L)
# Relative-error window for accepting an observed ratio against a theoretical
# envelope ratio (M+2 / M+0 = n * minor/major).
.ISOTOPE_ENVELOPE_REL_TOL <- 0.35
# M+2 shift -> halogen element (for envelope fitting).
.ISOTOPE_HALOGEN_ELEMENT_BY_SHIFT <- c(
  "1.997050" = "Cl",
  "1.997953" = "Br",
  "1.995796" = "S"
)

# Legacy single-halogen intensity-ratio bands (matched by shift value).
# Retained so the historical high-confidence single-halogen acceptance is
# byte-for-byte unchanged; multi-halogen ratios (2-4 Cl/Br, 2 S) are now
# *additionally* accepted through the envelope fit in
# `.nts_halogen_envelope_fit()` (previously they were silently rejected).
.ISOTOPE_HALOGEN_RATIO_TABLE <- list(
  "1.997050" = c(0.10, 0.50),   # 37Cl/35Cl
  "1.997953" = c(0.50, 1.50),   # 81Br/79Br
  "1.995796" = c(0.01, 0.10)    # 34S/32S
)

#' @noRd
.nts_build_isotope_target_diffs <- function(max_isotopes, halogen_diffs) {
  max_isotopes <- as.integer(max_isotopes)[1L]
  if (!is.finite(max_isotopes) || max_isotopes < 1L) {
    stop("max_isotopes must be a positive integer", call. = FALSE)
  }
  c13 <- seq_len(max_isotopes) * .C13_MZ_DELTA
  if (is.null(halogen_diffs)) {
    halogen_diffs <- numeric(0)
  } else {
    halogen_diffs <- as.numeric(halogen_diffs)
    halogen_diffs <- halogen_diffs[is.finite(halogen_diffs) & halogen_diffs > 0]
  }
  sort(unique(c(c13, halogen_diffs)))
}

#' @noRd
.nts_halogen_args <- function(halogen_diffs, halogen_labels) {
  if (is.null(halogen_diffs) || !length(halogen_diffs)) {
    return(list(vals = numeric(0), labels = character(0)))
  }
  vals <- suppressWarnings(as.numeric(halogen_diffs))
  if (!all(is.finite(vals)) || any(vals <= 0)) {
    stop("halogen_diffs must contain positive finite numeric values.", call. = FALSE)
  }
  n <- length(vals)
  labels <- rep(NA_character_, n)
  nm <- names(halogen_diffs)
  if (!is.null(nm) && length(nm) == n) {
    labels <- ifelse(nzchar(nm), nm, NA_character_)
  }
  if (!is.null(halogen_labels)) {
    hl <- as.character(halogen_labels)
    if (length(hl) != n) {
      stop("halogen_labels must have the same length as halogen_diffs.", call. = FALSE)
    }
    ok <- !is.na(hl) & nzchar(hl)
    labels[ok] <- hl[ok]
  }
  list(vals = vals, labels = labels)
}

#' @noRd
.nts_isotope_detail_from_shifts <- function(target_diff, max_isotopes,
                                            halogen_vals, halogen_labels,
                                            halogen_n = NULL) {
  n <- length(target_diff)
  if (!n) {
    return(data.frame(
      isotope_shift_da = numeric(0),
      isotope_kind = character(0),
      isotope_label = character(0),
      isotope_c13_n = integer(0),
      stringsAsFactors = FALSE
    ))
  }
  if (!is.null(halogen_n) && length(halogen_n) != n) {
    stop("Internal error: halogen_n length mismatch.", call. = FALSE)
  }
  c13_steps <- seq_len(max_isotopes) * .C13_MZ_DELTA
  kind <- rep(NA_character_, n)
  lab <- rep(NA_character_, n)
  ccn <- rep(NA_integer_, n)
  hv <- halogen_vals
  hn <- halogen_labels
  if (length(hv) && !length(hn)) {
    hn <- rep(NA_character_, length(hv))
  }
  if (length(hv) && length(hn) != length(hv)) {
    stop("Internal error: halogen label length mismatch.", call. = FALSE)
  }
  mv <- .ISOTOPE_MINOR_DEFAULTS
  tol_c13 <- 1e-5
  for (i in seq_len(n)) {
    td <- target_diff[i]
    if (!is.finite(td)) {
      next
    }
    matched_c13 <- FALSE
    for (j in seq_along(c13_steps)) {
      if (abs(c13_steps[j] - td) < tol_c13) {
        kind[i] <- "13C"
        ccn[i] <- as.integer(j)
        lab[i] <- paste0("13C+", j)
        matched_c13 <- TRUE
        break
      }
    }
    if (matched_c13) {
      next
    }
    if (length(hv)) {
      dh <- abs(hv - td)
      jh <- which.min(dh)
      tol_h <- max(1e-5, abs(hv[jh]) * 1e-9)
      if (dh[jh] < tol_h) {
        kind[i] <- "halogen"
        hlab <- hn[jh]
        hlab <- if (!is.na(hlab) && nzchar(hlab)) {
          hlab
        } else {
          paste0("halogen_", format(td, digits = 7, trim = TRUE))
        }
        # Multi-atom halogen envelope (n >= 2): state the count explicitly.
        if (!is.null(halogen_n) && is.finite(halogen_n[i]) && halogen_n[i] >= 2) {
          hlab <- paste0(sub("/.*$", "", hlab), "x", as.integer(halogen_n[i]))
        }
        lab[i] <- hlab
        next
      }
    }
    if (length(mv)) {
      dm <- abs(mv - td)
      jm <- which.min(dm)
      tol_m <- max(1e-5, abs(mv[jm]) * 1e-9)
      if (dm[jm] < tol_m) {
        kind[i] <- "minor"
        mlab <- names(mv)[jm]
        lab[i] <- if (!is.na(mlab) && nzchar(mlab)) {
          mlab
        } else {
          paste0("minor_", format(td, digits = 7, trim = TRUE))
        }
      }
    }
  }
  data.frame(
    isotope_shift_da = target_diff,
    isotope_kind = kind,
    isotope_label = lab,
    isotope_c13_n = ccn,
    stringsAsFactors = FALSE
  )
}

#' Classify a target mass shift as 13C, halogen, or minor (N/O) isotope
#' @noRd
.nts_classify_target_diff <- function(target_diff, max_isotopes) {
  c13_steps <- seq_len(max_isotopes) * .C13_MZ_DELTA
  for (j in seq_along(c13_steps)) {
    if (abs(c13_steps[j] - target_diff) < 1e-5) {
      return(list(kind = "13C", c13_n = j))
    }
  }
  hv <- .ISOTOPE_HALOGEN_DEFAULTS
  if (length(hv)) {
    dh <- abs(hv - target_diff)
    jh <- which.min(dh)
    tol_h <- max(1e-5, abs(hv[jh]) * 1e-9)
    if (dh[jh] < tol_h) {
      return(list(kind = "halogen", shift = unname(hv[jh]), label = names(hv)[jh]))
    }
  }
  mv <- .ISOTOPE_MINOR_DEFAULTS
  if (length(mv)) {
    dm <- abs(mv - target_diff)
    jm <- which.min(dm)
    tol_m <- max(1e-5, abs(mv[jm]) * 1e-9)
    if (dm[jm] < tol_m) {
      return(list(kind = "minor", shift = unname(mv[jm]), label = names(mv)[jm]))
    }
  }
  list(kind = NA_character_)
}


#' Theoretical isotope M+2 / M+0 envelope ratios for Cl/Br/S/N/O
#'
#' For `n` atoms of an element, the M+2 (single minor isotope) intensity
#' relative to M+0 (all major) equals `n * minor/major`, computed from exact
#' natural abundances rather than hard-coded bands.  For `n = 1` this
#' reproduces the classical \eqn{^{37}}Cl/\eqn{^{35}}Cl \eqn{\approx} 0.32,
#' \eqn{^{81}}Br/\eqn{^{79}}Br \eqn{\approx} 0.97, \eqn{^{34}}S/\eqn{^{32}}S
#' \eqn{\approx} 0.045 (and small values for N/O).
#' @noRd
.nts_halogen_envelope_table <- function(elements = c("Cl", "Br", "S")) {
  out <- list()
  for (el in elements) {
    ab <- .ISOTOPE_ABUNDANCE[[el]]
    if (is.null(ab)) next
    ratio <- unname(ab[["minor"]] / ab[["major"]])
    nmax <- .ISOTOPE_HALOGEN_MAX_N[[el]]
    if (is.null(nmax) || !is.finite(nmax) || nmax < 1L) nmax <- 1L
    ns <- seq_len(nmax)
    out[[el]] <- data.frame(
      element = el,
      n = ns,
      r_theory = ns * ratio,
      stringsAsFactors = FALSE
    )
  }
  if (!length(out)) {
    return(data.frame(element = character(0), n = integer(0),
                      r_theory = numeric(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, out)
}


#' Fit an observed child/parent intensity ratio to a theoretical isotope envelope
#'
#' Envelope-based replacement for a hard-coded ratio band: the closest
#' theoretical `M+2 / M+0` ratio (over `elements` and isotope counts 1..max_n)
#' is found and accepted when its relative error is within
#' `.ISOTOPE_ENVELOPE_REL_TOL`.  This makes the validator aware of multi-atom
#' envelopes (e.g. 2 Cl \eqn{\approx} 0.64, 4 Cl \eqn{\approx} 1.29,
#' 4 Br \eqn{\approx} 3.89) that a single-atom band silently rejects.
#'
#' @return list(accept, element, n, r_obs, r_theory, rel_err).
#' @noRd
.nts_halogen_envelope_fit <- function(cint, pint, elements = c("Cl", "Br", "S")) {
  fail <- list(
    accept = FALSE, element = NA_character_, n = NA_integer_,
    r_obs = NA_real_, r_theory = NA_real_, rel_err = NA_real_
  )
  if (length(cint) != 1L || length(pint) != 1L) return(fail)
  if (is.na(cint) || is.na(pint) || !is.finite(cint) || !is.finite(pint) || pint <= 0) {
    return(fail)
  }
  r <- cint / pint
  if (!is.finite(r) || r <= 0) return(fail)
  tab <- .nts_halogen_envelope_table(elements)
  if (!nrow(tab)) return(fail)
  rel <- abs(tab$r_theory - r) / tab$r_theory
  k <- which.min(rel)
  list(
    accept = is.finite(rel[k]) && rel[k] <= .ISOTOPE_ENVELOPE_REL_TOL,
    element = tab$element[k],
    n = as.integer(tab$n[k]),
    r_obs = r,
    r_theory = tab$r_theory[k],
    rel_err = rel[k]
  )
}


#' Validate intensity ratio (child / parent) for a given isotope shift
#'
#' Returns `TRUE` / `FALSE`, or `NA` when the ratio cannot be computed (missing
#' or non-positive parent intensity).  Callers treat `NA` as "unverified but
#' accepted" and count it (see `.flag_isotopes_df`), so an uncomputable ratio is
#' never silently indistinguishable from a validated one.
#' @noRd
.nts_validate_isotope_ratio <- function(cint, pint, target_diff, max_isotopes,
                                        include_minor_isotopes = FALSE) {
  if (is.na(cint) || is.na(pint) || !is.finite(cint) || !is.finite(pint) || pint <= 0) {
    return(NA)
  }
  ratio <- cint / pint
  cls <- .nts_classify_target_diff(target_diff, max_isotopes)
  if (isTRUE(cls$kind == "13C")) {
    n <- cls$c13_n
    # Minimum ratio decreases for higher 13C orders
    ratio_min <- if (n <= 1L) .ISOTOPE_C13_RATIO_MIN else if (n == 2L) 0.0005 else 0.00005
    return(ratio >= ratio_min && ratio <= .ISOTOPE_C13_RATIO_MAX)
  }
  if (isTRUE(cls$kind == "halogen")) {
    if (!is.null(cls$shift)) {
      key <- format(round(cls$shift, 6), nsmall = 6)
      rng <- .ISOTOPE_HALOGEN_RATIO_TABLE[[key]]
      # Legacy single-halogen band is kept verbatim (non-regressive).
      if (!is.null(rng) && ratio >= rng[1L] && ratio <= rng[2L]) {
        return(TRUE)
      }
    }
    # Multi-halogen-aware envelope fit (adds recall for n >= 2).
    return(isTRUE(.nts_halogen_envelope_fit(cint, pint, c("Cl", "Br", "S"))$accept))
  }
  if (isTRUE(cls$kind == "minor")) {
    if (!include_minor_isotopes) return(FALSE)
    return(isTRUE(.nts_halogen_envelope_fit(cint, pint, c("N", "O"))$accept))
  }
  TRUE
}


#' Enforce 13C isotope chain consistency within each parent group
#'
#' Uses target_diff to derive 13C order (no dependency on pre-computed detail columns).
#' @noRd
.nts_enforce_isotope_chain <- function(picks) {
  if (!nrow(picks)) return(picks)
  pick_dt <- data.table::as.data.table(picks)
  if (!"target_diff" %in% names(pick_dt)) return(picks)

  # Derive 13C order directly from target_diff (M+1 = 1.003355, M+2 = 2.006710, etc.)
  c13_n <- rep(NA_integer_, nrow(pick_dt))
  for (k in 1:6) {  # up to 6 13C steps
    expected <- k * .C13_MZ_DELTA
    mask <- abs(pick_dt$target_diff - expected) < 1e-5
    c13_n[mask] <- k
  }
  pick_dt$c13_n <- c13_n

  keep <- rep(TRUE, nrow(pick_dt))
  parent_names <- unique(pick_dt$parent_name[!is.na(pick_dt$parent_name)])
  for (pg in parent_names) {
    idx <- which(pick_dt$parent_name == pg)

    sub <- pick_dt[idx]
    sub$.orig_idx <- idx
    sub <- sub[!is.na(sub$c13_n) & sub$c13_n > 0L]
    if (nrow(sub) < 1L) next
    data.table::setorder(sub, c13_n)

    # All 13C children must form a complete chain starting from 1
    prev_int <- Inf
    prev_n   <- 0L
    for (r in seq_len(nrow(sub))) {
      n_val <- sub$c13_n[r]
      if (n_val != prev_n + 1L) {
        keep[sub$.orig_idx[r]] <- FALSE
        next
      }
      cint <- sub$cint[r]
      pin  <- sub$pint[r]
      if (!is.na(cint) && !is.na(pin) && cint > pin) {
        keep[sub$.orig_idx[r]] <- FALSE
        next
      }
      if (!is.na(cint) && is.finite(cint) && is.finite(prev_int) && cint > prev_int) {
        keep[sub$.orig_idx[r]] <- FALSE
        next
      }
      if (is.finite(cint)) prev_int <- cint
      prev_n <- n_val
    }
  }
  picks[keep, , drop = FALSE]
}


#' Resolve M+0 parent conflicts: strip parent role from features that are
#' themselves isotopes of a stronger peak. Does NOT reassign children (which
#' would require recomputing isotope detail columns and can produce wrong labels).
#' @noRd
.nts_resolve_parent_conflicts <- function(ft) {
  if (!nrow(ft) || !"is_isotope" %in% names(ft) || !"parent_name" %in% names(ft)) {
    return(ft)
  }
  iso_mask <- ft$is_isotope %in% TRUE & !is.na(ft$parent_name) & nzchar(ft$parent_name)
  if (!any(iso_mask)) return(ft)

  all_parents <- unique(ft$parent_name[iso_mask])
  all_parents <- all_parents[!is.na(all_parents) & nzchar(all_parents)]
  if (!length(all_parents)) return(ft)

  # Features that act as parent AND are marked as child of someone else
  dup <- all_parents[all_parents %in% ft$name[iso_mask]]
  if (!length(dup)) return(ft)

  int_all <- .nts_isotope_row_intensity(as.data.frame(ft, stringsAsFactors = FALSE))
  names(int_all) <- ft$name

  for (dp in dup) {
    claimant <- unique(ft$parent_name[ft$name == dp & iso_mask])
    claimant <- claimant[!is.na(claimant) & nzchar(claimant) & claimant != dp]
    if (!length(claimant)) next

    int_dp       <- int_all[dp]
    int_claimant <- int_all[claimant[1L]]

    if (!is.na(int_claimant) && !is.na(int_dp) && int_claimant > int_dp) {
      # Claimant is stronger: strip dp's parent role
      # Children of dp lose their parent (conservative - no reassignment)
      child_idx <- which(iso_mask & ft$parent_name == dp)
      if (length(child_idx)) {
        ft$parent_name[child_idx] <- NA_character_
        ft$is_isotope[child_idx] <- FALSE
        # Clear the whole isotope detail so no row keeps an isotope label /
        # kind / shift while `is_isotope` is FALSE (downstream code that keys
        # on `isotope_kind` would otherwise misread these as isotopes).
        if ("isotope_shift_da" %in% names(ft)) ft$isotope_shift_da[child_idx] <- NA_real_
        if ("isotope_kind" %in% names(ft)) ft$isotope_kind[child_idx] <- NA_character_
        if ("isotope_label" %in% names(ft)) ft$isotope_label[child_idx] <- NA_character_
        if ("isotope_c13_n" %in% names(ft)) ft$isotope_c13_n[child_idx] <- NA_integer_
      }
    }
  }
  ft
}


#' @noRd
.nts_isotope_row_intensity <- function(ft) {
  if ("intensity" %in% names(ft)) {
    return(as.numeric(ft[["intensity"]]))
  }
  reserved <- c(
    .nts_feature_shape_cols,
    "name", "feature_id", "row_ID",
    "is_isf", "parent_mz", "isf_annotation",
    "is_isotope", "parent_name",
    "isotope_shift_da", "isotope_kind", "isotope_label", "isotope_c13_n",
    "homologous_series_id", "kmd", "kmd_series_id",
    "standard_md", "KM", "NKM", "KMD", "mz_remainder"
  )
  ab <- setdiff(names(ft), reserved)
  ab <- ab[vapply(ft[, ab, drop = FALSE], is.numeric, logical(1L))]
  if (!length(ab)) {
    return(rep(NA_real_, nrow(ft)))
  }
  rowMeans(as.matrix(ft[, ab, drop = FALSE]), na.rm = TRUE)
}


#' Flag probable isotope peaks (M+n)
#'
#' @description
#' Marks features that are probable isotope ions relative to monoisotopic parents
#' (`parent_name`), with detail columns `isotope_shift_da`, `isotope_kind`
#' (`13C` or `halogen`), human-readable `isotope_label` (e.g. `13C+1`, `37Cl`),
#' and `isotope_c13_n` (integer order for \eqn{^{13}}C, `NA` for halogen).
#' For [NTSProject], those columns are written to \verb{@ms1_annotations}
#' (together with `is_isotope` / `parent_name`); \verb{@feature_table} stays
#' the quantitative core.
#'
#' Scans each target mass shift: `n * 1.003355` Da for
#' `n = 1, ..., max_isotopes` (\eqn{^{13}}C) plus optional `halogen_diffs` (e.g.
#' Cl/Br spacing), each within `mz_tol_ppm` and `rt_tol`. Abundance rule: child
#' mean intensity `< 1.5 *` parent when both are finite. Conflicts (one child
#' matching several parents) resolve to the **strongest parent** by intensity,
#' then smallest m/z error.
#'
#' Halogen spacing labels: use **named** `halogen_diffs` (names become defaults for
#' `isotope_label`), or pass `halogen_labels` (same length as `halogen_diffs`;
#' non-empty strings override names). Unnamed halogen hits get
#' `halogen_<shift>` style labels.
#'
#' @details
#' Uses [data.table::foverlaps] for fast interval matching on m/z.
#'
#' @param df A [NTSProject], `data.frame`, or `data.table` with `mz` and `rt`.
#' @param mz_tol_ppm m/z tolerance (ppm, from parent m/z) around each target shift.
#' @param rt_tol Retention-time tolerance (seconds).
#' @param max_isotopes Maximum 13C order to scan (default `3` => M+1..M+3 shifts).
#' @param halogen_diffs Optional numeric vector of extra mass shifts (Da), e.g.
#'   Cl/Br isotope spacing; use `numeric(0)` or `NULL` to skip.
#' @param halogen_labels Optional character vector parallel to `halogen_diffs`
#'   giving labels such as `37Cl` / `81Br`. Empty names in `halogen_diffs` are
#'   filled from here; explicit non-empty strings override `names(halogen_diffs)`.
#' @param include_minor_isotopes Logical (default `FALSE`). When `TRUE`, also
#'   scan for 15N (`0.997035` Da) and 18O (`2.004246` Da) minor-isotope shifts,
#'   which share the ~1-2 Da region with the 13C/37Cl/81Br/34S series and are
#'   therefore only searched when their exact-mass shifts are resolvable.
#'
#' @details
#' `isotope_kind` takes three values: `"13C"` for the M+1/M+2/... carbon
#' series, `"halogen"` for the 37Cl / 81Br / 34S ~2 Da shifts (validated by
#' theoretical isotope-envelope fitting, so multi-halogen patterns such as
#' 2Cl / 3Cl / 4Br are accepted), and `"minor"` for 15N / 18O (only when
#' `include_minor_isotopes = TRUE`).
#'
#' @return Table input: adds `is_isotope`, `parent_name`, `isotope_shift_da`,
#'   `isotope_kind` (`"13C"`, `"halogen"`, or `"minor"`), `isotope_label`, and
#'   `isotope_c13_n`. [NTSProject]: the isotope columns live on
#'   \verb{@ms1_annotations} only; \verb{@feature_table} width unchanged.
#' @export
flag_isotopes <- function(df,
                          mz_tol_ppm = 10,
                          rt_tol = 5,
                          max_isotopes = 3,
                          halogen_diffs = .ISOTOPE_HALOGEN_DEFAULTS,
                          halogen_labels = NULL,
                          include_minor_isotopes = FALSE) {
  if (is.null(halogen_diffs)) {
    halogen_diffs <- numeric(0)
  }
  if (!is.null(halogen_labels) && !length(halogen_diffs)) {
    stop("halogen_labels requires a non-empty halogen_diffs vector.", call. = FALSE)
  }
  include_minor_isotopes <- isTRUE(include_minor_isotopes)
  if (inherits(df, "NTSProject")) {
    proj <- df
    ft <- get_feature_table(proj)
    if (!nrow(ft)) {
      stop("feature_table is empty", call. = FALSE)
    }
    out <- .flag_isotopes_df(
      ft,
      mz_tol_ppm = mz_tol_ppm,
      rt_tol = rt_tol,
      max_isotopes = max_isotopes,
      halogen_diffs = halogen_diffs,
      halogen_labels = halogen_labels,
      include_minor_isotopes = include_minor_isotopes
    )
    iso_cols <- c(
      "is_isotope", "parent_name",
      "isotope_shift_da", "isotope_kind", "isotope_label", "isotope_c13_n"
    )
    proj <- set_feature_table(proj, out[, setdiff(names(out), iso_cols), drop = FALSE])

    ann <- get_ms1_annotations(proj)
    if (nrow(ann) && "name" %in% names(ann) && "name" %in% names(out)) {
      map <- out[, c("name", iso_cols), drop = FALSE]
      idx <- match(as.character(ann$name), as.character(map$name))
      ann$is_isotope <- ifelse(is.na(idx), FALSE, as.logical(map$is_isotope[idx]) %in% TRUE)
      ann$parent_name <- ifelse(is.na(idx), NA_character_, as.character(map$parent_name[idx]))
      ann$isotope_shift_da <- suppressWarnings(
        ifelse(is.na(idx), NA_real_, as.numeric(map$isotope_shift_da[idx]))
      )
      ann$isotope_kind <- ifelse(is.na(idx), NA_character_, as.character(map$isotope_kind[idx]))
      ann$isotope_label <- ifelse(is.na(idx), NA_character_, as.character(map$isotope_label[idx]))
      ann$isotope_c13_n <- suppressWarnings(
        ifelse(is.na(idx), NA_integer_, as.integer(map$isotope_c13_n[idx]))
      )
      proj <- set_ms1_annotations(proj, .nts_ms1_sanitize_and_order(ann))
    }

    methods::slot(proj, "log") <- .nts_log_append(
      methods::slot(proj, "log"),
      "flag_isotopes",
      list(
        mz_tol_ppm = mz_tol_ppm,
        rt_tol = rt_tol,
        max_isotopes = max_isotopes,
        halogen_diffs = as.numeric(
          if (is.null(halogen_diffs)) numeric(0) else halogen_diffs
        ),
        halogen_labels = if (is.null(halogen_labels)) {
          NULL
        } else {
          as.character(halogen_labels)
        },
        include_minor_isotopes = include_minor_isotopes
      )
    )
    return(proj)
  }

  .flag_isotopes_df(
    df,
    mz_tol_ppm = mz_tol_ppm,
    rt_tol = rt_tol,
    max_isotopes = max_isotopes,
    halogen_diffs = halogen_diffs,
    halogen_labels = halogen_labels,
    include_minor_isotopes = include_minor_isotopes
  )
}


#' @noRd
.flag_isotopes_df <- function(feature_table,
                              mz_tol_ppm = 10,
                              rt_tol = 5,
                              max_isotopes = 3,
                              halogen_diffs = .ISOTOPE_HALOGEN_DEFAULTS,
                              halogen_labels = NULL,
                              include_minor_isotopes = FALSE) {
  if (is.null(halogen_diffs)) {
    halogen_diffs <- numeric(0)
  }
  hp <- .nts_halogen_args(halogen_diffs, halogen_labels)
  if (!data.table::is.data.table(feature_table)) {
    ft <- data.table::as.data.table(feature_table)
    as_dt <- FALSE
  } else {
    ft <- data.table::copy(feature_table)
    as_dt <- TRUE
  }

  stopifnot(
    nrow(ft) >= 1L,
    "mz" %in% names(ft),
    "rt" %in% names(ft),
    is.numeric(mz_tol_ppm), length(mz_tol_ppm) == 1L, is.finite(mz_tol_ppm), mz_tol_ppm >= 0,
    is.numeric(rt_tol), length(rt_tol) == 1L, is.finite(rt_tol), rt_tol >= 0
  )

  if (!"name" %in% names(ft)) {
    ft[, name := sprintf("FT%05d", .I)]
  }

  mz_tol_ppm <- as.numeric(mz_tol_ppm)[1L]
  rt_tol <- as.numeric(rt_tol)[1L]
  include_minor_isotopes <- isTRUE(include_minor_isotopes)
  target_diffs <- .nts_build_isotope_target_diffs(max_isotopes, hp$vals)
  if (include_minor_isotopes) {
    target_diffs <- sort(unique(c(target_diffs, .ISOTOPE_MINOR_DEFAULTS)))
  }
  n_ratio_skipped <- 0L

  int_vec <- .nts_isotope_row_intensity(as.data.frame(ft, stringsAsFactors = FALSE))
  ft[, row_id := seq_len(.N)]
  ft[, pmz := as.numeric(mz)]
  ft[, prt := as.numeric(rt)]
  ft[, pint := int_vec]

  cand <- ft[, .(
    child_row_id = row_id,
    child_name = as.character(name),
    cmz = pmz,
    crt = prt,
    cint = pint,
    imz_lo = pmz,
    imz_hi = pmz
  )]
  data.table::setkey(cand, imz_lo, imz_hi)

  hit_parts <- vector("list", length(target_diffs))
  for (ki in seq_along(target_diffs)) {
    diff_val <- target_diffs[ki]
    parents <- ft[, .(
      parent_row_id = row_id,
      parent_name = as.character(name),
      pmz,
      prt,
      pint,
      mz_lo = pmz + diff_val - abs(pmz) * mz_tol_ppm / 1e6,
      mz_hi = pmz + diff_val + abs(pmz) * mz_tol_ppm / 1e6
    )]
    data.table::setkey(parents, mz_lo, mz_hi)

    hits <- data.table::foverlaps(
      cand,
      parents,
      type = "within",
      nomatch = NA,
      mult = "all"
    )

    hits <- hits[
      is.finite(parent_row_id) &
        child_row_id != parent_row_id &
        abs(crt - prt) <= rt_tol
    ]

    if (!nrow(hits)) {
      next
    }

    hits[, target_diff := diff_val]
    hits[, exp_mz := pmz + diff_val]
    hits[, mz_err := abs(cmz - exp_mz)]

    # Basic intensity cap
    hits <- hits[
      is.na(pint) | is.na(cint) | (cint < pint * 1.5)
    ]

    # Strict isotope ratio validation. An uncomputable ratio returns NA and is
    # ACCEPTED (and counted in n_ratio_skipped), exactly as the log message and
    # the `is.na(...)` branch below promise: only an explicit FALSE (a ratio
    # that was computed and failed) drops the candidate.
    if (nrow(hits)) {
      keep <- vapply(seq_len(nrow(hits)), function(i) {
        .nts_validate_isotope_ratio(
          hits$cint[i], hits$pint[i], hits$target_diff[i], max_isotopes,
          include_minor_isotopes = include_minor_isotopes
        )
      }, logical(1L))
      n_ratio_skipped <- n_ratio_skipped + sum(is.na(keep))
      hits <- hits[keep %in% TRUE | is.na(keep)]
    }

    # Envelope fit for halogen / minor-isotope hits: records the matched
    # isotope count (n) so n >= 2 can be labelled explicitly.
    if (nrow(hits)) {
      hits[, hl_n := NA_integer_]
      hits[, hl_element := NA_character_]
      kinds <- vapply(hits$target_diff, function(td) {
        k <- .nts_classify_target_diff(td, max_isotopes)$kind
        if (is.null(k) || is.na(k)) NA_character_ else as.character(k)
      }, character(1L))
      for (i in seq_along(kinds)) {
        if (is.na(kinds[i])) next
        els <- switch(kinds[i], halogen = c("Cl", "Br", "S"),
                      minor = c("N", "O"), NULL)
        if (is.null(els)) next
        f <- .nts_halogen_envelope_fit(hits$cint[i], hits$pint[i], els)
        if (is.finite(f$n)) hits$hl_n[i] <- as.integer(f$n)
        hits$hl_element[i] <- f$element
      }
    }

    if (nrow(hits)) {
      hit_parts[[ki]] <- hits
    }
  }

  hits <- data.table::rbindlist(hit_parts, use.names = TRUE, fill = TRUE)

  if (nrow(hits)) {
    if (!"hl_n" %in% names(hits)) hits[, hl_n := NA_integer_]
    data.table::setorder(hits, child_row_id, -pint, mz_err)
    picks <- hits[, .SD[1L], by = child_row_id]
  } else {
    picks <- data.table::data.table()
  }

  # Chain consistency: M+N only valid if M+(N-1) exists and intensity decays
  if (nrow(picks)) {
    picks <- .nts_enforce_isotope_chain(picks)
  }

  ft[, is_isotope := FALSE]
  ft[, parent_name := NA_character_]
  ft[, isotope_shift_da := NA_real_]
  ft[, isotope_kind := NA_character_]
  ft[, isotope_label := NA_character_]
  ft[, isotope_c13_n := NA_integer_]
  if (nrow(picks)) {
    det <- .nts_isotope_detail_from_shifts(
      picks$target_diff,
      max_isotopes,
      hp$vals,
      hp$labels,
      halogen_n = picks$hl_n
    )
    ft[picks$child_row_id, `:=`(
      is_isotope = TRUE,
      parent_name = picks$parent_name,
      isotope_shift_da = det$isotope_shift_da,
      isotope_kind = det$isotope_kind,
      isotope_label = det$isotope_label,
      isotope_c13_n = det$isotope_c13_n
    )]
  }

  # M+0 parent conflict resolution: strongest peak in RT cluster wins
  ft <- .nts_resolve_parent_conflicts(ft)

  ft[, c("row_id", "pmz", "prt", "pint") := NULL]

  if (n_ratio_skipped > 0L) {
    message(
      "- flag_isotopes: ", n_ratio_skipped,
      " candidate ratio(s) could not be computed (missing / non-positive ",
      "parent intensity) and were accepted without ratio validation."
    )
  }

  if (!as_dt) {
    return(as.data.frame(ft, stringsAsFactors = FALSE))
  }
  ft[]
}

#' PFAS / matrix-relevant neutral losses (exact mass, Da; from [calculate_exact_mass])
#'
#' Single source of truth for [recognize_isf]: returns the same loss set as
#' [flag_isf] (`.nts_default_neutral_losses_isf()` = HF, H2O, CO2, CF2, CH2,
#' SO3) so the two ISF entry points cannot drift apart.
#' @noRd
.nts_isf_nl_library_vec <- function() {
  unname(unlist(.nts_default_neutral_losses_isf(), use.names = TRUE))
}


# Shared ISF defaults - one source of truth for BOTH flag_isf and
# recognize_isf so the two entry points cannot drift apart.
.ISF_DEFAULT_RT_TOL      <- 2      # co-elution window (seconds)
.ISF_DEFAULT_CORR_THRESH <- 0.85   # minimum Pearson r across samples
.ISF_DEFAULT_LOSS_TOL    <- 0.005  # neutral-loss mass tolerance (Da)


#' Co-elution ISF recognition from correlation and neutral-loss mass rules
#'
#' Within sliding RT windows, Pearson correlations are computed **once** per
#' window via [stats::cor] on the transposed abundance matrix (no per-pair
#' correlation loops). Pairs with correlation >= `corr_thresh` and
#' |m/z_parent - m/z_fragment| matching a tabulated neutral loss within
#' `loss_tol` (default 0.005 Da) label the **lower** m/z feature as ISF and set
#' `parent_mz` to the larger m/z.
#'
#' Features already flagged as isotope peaks are excluded from candidacy (as
#' both fragment and parent) when `exclude_isotopes = TRUE` - an isotope peak is
#' the same molecule as its parent, not a distinct fragment.
#'
#' @param feature_table `data.frame`, `data.table`, or [NTSProject] (for
#'   projects, updates `@ms1_annotations` only; `@feature_table` column set is
#'   unchanged aside from any pre-existing stray ISF columns being dropped).
#' @param corr_thresh Minimum Pearson r for co-evolution (default `0.85`,
#'   shared with [flag_isf]).
#' @param rt_window Maximum RT difference for the same local window (same units
#'   as `rt`; default `2`, shared with [flag_isf]).
#' @param abundance_cols Optional character vector of sample intensity columns;
#'   if `NULL`, all numeric columns except `mz`, `rt`, `name`, `is_isf`,
#'   `parent_mz`, `homologous_series_id` are used.
#' @param loss_tol Absolute Da tolerance for neutral-loss match (default
#'   `0.005`, shared with [flag_isf]).
#' @param exclude_isotopes Logical (default `TRUE`). For a [NTSProject], features
#'   with `is_isotope` TRUE on `@ms1_annotations` are dropped from candidacy. For
#'   a plain table the same applies when the table carries an `is_isotope`
#'   column.
#'
#' @return Updated `data.table` / `data.frame` with `is_isf` and `parent_mz`, or
#'   an updated [NTSProject] with those columns on `@ms1_annotations` only.
#' @export
recognize_isf <- function(feature_table,
                          corr_thresh = .ISF_DEFAULT_CORR_THRESH,
                          rt_window = .ISF_DEFAULT_RT_TOL,
                          abundance_cols = NULL,
                          loss_tol = .ISF_DEFAULT_LOSS_TOL,
                          exclude_isotopes = TRUE) {
  exclude_isotopes <- isTRUE(exclude_isotopes)
  if (inherits(feature_table, "NTSProject")) {
    proj <- feature_table
    ft <- get_feature_table(proj)
    if (!nrow(ft)) {
      stop("feature_table is empty", call. = FALSE)
    }
    if (!"name" %in% names(ft)) {
      stop("feature_table must contain column `name` for NTSProject recognize_isf()", call. = FALSE)
    }
    iso_names <- .nts_isotope_names_from_project(proj, exclude_isotopes)
    out <- .recognize_isf_df(
      ft,
      corr_thresh = corr_thresh,
      rt_window = rt_window,
      abundance_cols = abundance_cols,
      loss_tol = loss_tol,
      isotope_names = iso_names
    )
    isf_cols <- c("is_isf", "parent_mz")
    proj <- set_feature_table(proj, out[, setdiff(names(out), isf_cols), drop = FALSE])

    ann <- get_ms1_annotations(proj)
    if (nrow(ann) && "name" %in% names(ann) && "name" %in% names(out)) {
      map <- out[, c("name", isf_cols), drop = FALSE]
      idx <- match(as.character(ann$name), as.character(map$name))
      ann$is_isf <- ifelse(is.na(idx), FALSE, as.logical(map$is_isf[idx]) %in% TRUE)
      ann$parent_mz <- ifelse(is.na(idx), NA_real_, as.numeric(map$parent_mz[idx]))
      proj <- set_ms1_annotations(proj, .nts_ms1_sanitize_and_order(ann))
    }

    methods::slot(proj, "log") <- .nts_log_append(
      methods::slot(proj, "log"),
      "recognize_isf",
      list(
        corr_thresh = corr_thresh,
        rt_window = rt_window,
        loss_tol = loss_tol,
        exclude_isotopes = exclude_isotopes
      )
    )
    return(proj)
  }
  .recognize_isf_df(
    feature_table,
    corr_thresh = corr_thresh,
    rt_window = rt_window,
    abundance_cols = abundance_cols,
    loss_tol = loss_tol
  )
}


#' Names of features flagged as isotopes (for ISF candidate exclusion)
#'
#' Returns a character vector of feature names with `is_isotope` TRUE on
#' `@ms1_annotations`, or `NULL` when `enabled` is FALSE / no such data exists.
#' @noRd
.nts_isotope_names_from_project <- function(proj, enabled = TRUE) {
  if (!isTRUE(enabled)) return(NULL)
  ann <- tryCatch(get_ms1_annotations(proj), error = function(e) NULL)
  if (is.null(ann) || !nrow(ann)) return(NULL)
  if (!all(c("name", "is_isotope") %in% names(ann))) return(NULL)
  nm <- as.character(ann$name[ann$is_isotope %in% TRUE])
  if (!length(nm)) return(NULL)
  nm
}


#' @noRd
.recognize_isf_df <- function(feature_table,
                               corr_thresh = .ISF_DEFAULT_CORR_THRESH,
                               rt_window = .ISF_DEFAULT_RT_TOL,
                               abundance_cols = NULL,
                               loss_tol = .ISF_DEFAULT_LOSS_TOL,
                               isotope_names = NULL) {
  if (!data.table::is.data.table(feature_table)) {
    ft <- data.table::as.data.table(feature_table)
    as_dt <- FALSE
  } else {
    ft <- data.table::copy(feature_table)
    as_dt <- TRUE
  }
  stopifnot(
    nrow(ft) >= 1L,
    "mz" %in% names(ft),
    "rt" %in% names(ft)
  )
  if (is.null(abundance_cols)) {
    reserved <- c(
      .nts_feature_shape_cols,
      "name", "feature_id",
      "is_isf", "parent_mz", "homologous_series_id", "is_isotope"
    )
    abundance_cols <- setdiff(names(ft), reserved)
    abundance_cols <- abundance_cols[
      vapply(ft[, abundance_cols, with = FALSE], is.numeric, logical(1L))
    ]
  }
  if (!length(abundance_cols)) {
    stop("No numeric abundance columns found; set `abundance_cols` explicitly.", call. = FALSE)
  }

  nl_lib <- .nts_isf_nl_library_vec()
  n <- nrow(ft)
  ord <- order(ft[["rt"]])
  rt_o <- ft[["rt"]][ord]
  mz_all <- ft[["mz"]]
  abund <- as.matrix(ft[, abundance_cols, with = FALSE])

  # Isotope peaks are the same molecule as their parent - exclude them from
  # candidacy in BOTH roles (fragment and parent).
  excl <- rep(FALSE, n)
  if ("is_isotope" %in% names(ft)) {
    excl <- excl | (ft[["is_isotope"]] %in% TRUE)
  }
  if (length(isotope_names) && "name" %in% names(ft)) {
    excl <- excl | (as.character(ft[["name"]]) %in% isotope_names)
  }

  pairs_acc <- list()
  pc <- 1L

  i <- 1L
  while (i <= n) {
    j <- i
    while (j < n && rt_o[j + 1L] - rt_o[i] <= rt_window) {
      j <- j + 1L
    }
    idx <- ord[seq.int(i, j)]
    idx <- idx[!excl[idx]]
    if (length(idx) >= 2L) {
      block <- abund[idx, , drop = FALSE]
      cmat <- stats::cor(t(block), use = "pairwise.complete.obs")
      cmat[is.na(cmat)] <- 0
      diag(cmat) <- 1
      nr <- length(idx)
      w <- which(upper.tri(cmat, diag = FALSE), arr.ind = TRUE)
      loc_i <- w[, 1L]
      loc_j <- w[, 2L]
      corv <- cmat[cbind(loc_i, loc_j)]
      gi <- idx[loc_i]
      gj <- idx[loc_j]
      mi <- mz_all[gi]
      mj <- mz_all[gj]
      s_row <- ifelse(mi < mj, gi, gj)
      p_row <- ifelse(mi < mj, gj, gi)
      dm <- abs(mi - mj)
      mat <- abs(outer(dm, nl_lib, "-"))
      nl_ok <- rowSums(mat <= loss_tol) >= 1L
      ok <- nl_ok & corv >= corr_thresh
      if (any(ok)) {
        pairs_acc[[pc]] <- data.table::data.table(
          s = s_row[ok],
          p = p_row[ok],
          cor = corv[ok]
        )
        pc <- pc + 1L
      }
    }
    i <- i + 1L
  }

  is_isf <- rep(FALSE, n)
  parent_mz <- rep(NA_real_, n)

  if (length(pairs_acc)) {
    pairs_all <- data.table::rbindlist(pairs_acc)
    pairs_all <- unique(pairs_all[order(-cor)], by = c("s", "p"))
    pairs_all <- unique(pairs_all[order(-cor)], by = "s")
    # Acyclic guard: a feature cannot be both a fragment and another
    # fragment's parent (in-source fragmentation is a single step). Keep the
    # parent role and drop the fragment role for any such feature.
    repeat {
      both <- intersect(pairs_all$s, pairs_all$p)
      if (!length(both)) break
      pairs_all <- pairs_all[!pairs_all$s %in% both]
    }
    is_isf[pairs_all$s] <- TRUE
    parent_mz[pairs_all$s] <- mz_all[pairs_all$p]
  }

  ft[, is_isf := is_isf]
  ft[, parent_mz := parent_mz]

  if (!as_dt) {
    return(as.data.frame(ft, stringsAsFactors = FALSE))
  }
  ft[]
}


#' Flag in-source fragments (ISF) by co-elution, neutral loss, and correlation
#'
#' Efficient ISF recognition using a three-stage pruning strategy:
#' 1. sort features by RT and only inspect local windows within `rt_tol`;
#' 2. within each window, vectorize all pairwise `delta_mz` values and match
#'    against a neutral-loss library using **absolute Da tolerance**;
#' 3. compute one Pearson correlation matrix per local window and retain pairs
#'    whose cross-sample correlation exceeds `corr_thresh`.
#'
#' If one fragment matches multiple parents, conflicts are resolved by:
#' 1. higher correlation first;
#' 2. smaller neutral-loss mass error second;
#' 3. larger parent m/z third.
#'
#' @param object A [NTSProject], `data.frame`, or `data.table`.
#' @param rt_tol Co-elution RT window in seconds (default `2`).
#' @param loss_tol_da Absolute neutral-loss tolerance in Da (default `0.005`).
#' @param corr_thresh Minimum Pearson correlation across sample abundances
#'   (default `0.85`).
#' @param neutral_losses Named list or vector of exact neutral losses. If
#'   `NULL`, uses [calculate_exact_mass] for `HF`, `H2O`, `CO2`, `CF2`, `CH2`,
#'   and `SO3`.
#' @param abundance_cols Optional abundance columns. If `NULL`, all numeric
#'   columns excluding metadata columns are used.
#' @param exclude_isotopes Logical (default `TRUE`). When `TRUE`, features
#'   already labelled as isotopes by [flag_isotopes] (`@ms1_annotations$is_isotope`)
#'   are removed from the ISF candidate set, so an isotopologue can never be
#'   reported as an in-source fragment. Parent chains are additionally forced
#'   acyclic (a fragment is never itself a parent).
#'
#' @return Updated `data.frame` / `data.table` with `is_isf`, `parent_mz`, and
#'   `isf_annotation`, or an updated [NTSProject] with those fields on
#'   `@ms1_annotations` only (`@feature_table` stays the abundance-only core).
#' @export
flag_isf <- function(object,
                     rt_tol = .ISF_DEFAULT_RT_TOL,
                     loss_tol_da = .ISF_DEFAULT_LOSS_TOL,
                     corr_thresh = .ISF_DEFAULT_CORR_THRESH,
                     neutral_losses = NULL,
                     abundance_cols = NULL,
                     exclude_isotopes = TRUE) {
  if (is.null(neutral_losses)) {
    neutral_losses <- .nts_default_neutral_losses_isf()
  }
  exclude_isotopes <- isTRUE(exclude_isotopes)
  if (inherits(object, "NTSProject")) {
    proj <- object
    ft <- get_feature_table(proj)
    if (!nrow(ft)) {
      stop("feature_table is empty", call. = FALSE)
    }
    iso_names <- .nts_isotope_names_from_project(proj, exclude_isotopes)
    out <- .flag_isf_df(
      ft,
      rt_tol = rt_tol,
      loss_tol_da = loss_tol_da,
      corr_thresh = corr_thresh,
      neutral_losses = neutral_losses,
      abundance_cols = abundance_cols,
      isotope_names = iso_names
    )
    isf_cols <- c("is_isf", "parent_mz", "isf_annotation")
    proj <- set_feature_table(proj, out[, setdiff(names(out), isf_cols), drop = FALSE])

    ann <- get_ms1_annotations(proj)
    if (nrow(ann) && "name" %in% names(ann) && "name" %in% names(out)) {
      map <- out[, c("name", isf_cols), drop = FALSE]
      idx <- match(as.character(ann$name), as.character(map$name))
      ann$is_isf <- ifelse(is.na(idx), FALSE, as.logical(map$is_isf[idx]))
      ann$parent_mz <- ifelse(is.na(idx), NA_real_, as.numeric(map$parent_mz[idx]))
      ann$isf_annotation <- ifelse(is.na(idx), NA_character_, as.character(map$isf_annotation[idx]))
      proj <- set_ms1_annotations(proj, .nts_ms1_sanitize_and_order(ann))
    }

    methods::slot(proj, "log") <- .nts_log_append(
      methods::slot(proj, "log"),
      "flag_isf",
      list(
        rt_tol = rt_tol,
        loss_tol_da = loss_tol_da,
        corr_thresh = corr_thresh,
        neutral_losses = as.list(neutral_losses),
        exclude_isotopes = exclude_isotopes
      )
    )
    return(proj)
  }

  .flag_isf_df(
    object,
    rt_tol = rt_tol,
    loss_tol_da = loss_tol_da,
    corr_thresh = corr_thresh,
    neutral_losses = neutral_losses,
    abundance_cols = abundance_cols
  )
}


#' @noRd
.nts_default_neutral_losses_isf <- function() {
  list(
    HF = calculate_exact_mass("HF"),
    H2O = calculate_exact_mass("H2O"),
    CO2 = calculate_exact_mass("CO2"),
    CF2 = calculate_exact_mass("CF2"),
    CH2 = calculate_exact_mass("CH2"),
    SO3 = calculate_exact_mass("SO3")
  )
}


#' @noRd
.flag_isf_df <- function(feature_table,
                         rt_tol = .ISF_DEFAULT_RT_TOL,
                         loss_tol_da = .ISF_DEFAULT_LOSS_TOL,
                         corr_thresh = .ISF_DEFAULT_CORR_THRESH,
                         neutral_losses,
                         abundance_cols = NULL,
                         isotope_names = NULL) {
  if (!data.table::is.data.table(feature_table)) {
    ft <- data.table::as.data.table(feature_table)
    as_dt <- FALSE
  } else {
    ft <- data.table::copy(feature_table)
    as_dt <- TRUE
  }

  stopifnot(
    nrow(ft) >= 1L,
    "mz" %in% names(ft),
    "rt" %in% names(ft),
    is.numeric(rt_tol), length(rt_tol) == 1L, is.finite(rt_tol), rt_tol >= 0,
    is.numeric(loss_tol_da), length(loss_tol_da) == 1L, is.finite(loss_tol_da), loss_tol_da >= 0,
    is.numeric(corr_thresh), length(corr_thresh) == 1L, is.finite(corr_thresh)
  )

  if (!"name" %in% names(ft)) {
    ft[, name := sprintf("FT%05d", .I)]
  }

  if (is.list(neutral_losses)) {
    nl_vals <- unlist(neutral_losses, use.names = TRUE)
  } else {
    nl_vals <- neutral_losses
  }
  nl_vals <- as.numeric(nl_vals)
  nl_nms <- names(unlist(neutral_losses, use.names = TRUE))
  if (is.null(nl_nms) || length(nl_nms) != length(nl_vals)) {
    nl_nms <- paste0("NL", seq_along(nl_vals))
  }
  keep_nl <- is.finite(nl_vals) & nzchar(nl_nms)
  nl_vals <- nl_vals[keep_nl]
  nl_nms <- nl_nms[keep_nl]
  if (!length(nl_vals)) {
    stop("neutral_losses must contain at least one finite value", call. = FALSE)
  }

  if (is.null(abundance_cols)) {
    reserved <- c(
      .nts_feature_shape_cols,
      "name", "feature_id",
      "is_isf", "parent_mz", "isf_annotation",
      "is_isotope", "parent_name",
      "isotope_shift_da", "isotope_kind", "isotope_label", "isotope_c13_n",
      "inferred_formula", "annotation_level",
      "homologous_series_id", "kmd"
    )
    abundance_cols <- setdiff(names(ft), reserved)
    abundance_cols <- abundance_cols[
      vapply(ft[, abundance_cols, with = FALSE], is.numeric, logical(1L))
    ]
  }
  if (!length(abundance_cols)) {
    stop("No numeric abundance columns found; set `abundance_cols` explicitly.", call. = FALSE)
  }

  n <- nrow(ft)
  ord <- order(ft[["rt"]], ft[["mz"]])
  rt_o <- as.numeric(ft[["rt"]][ord])
  mz_all <- as.numeric(ft[["mz"]])
  name_all <- as.character(ft[["name"]])
  abund <- as.matrix(ft[, abundance_cols, with = FALSE])

  # Isotope peaks are the same molecule as their parent - exclude them from
  # candidacy in BOTH roles (fragment and parent).
  excl <- rep(FALSE, n)
  if ("is_isotope" %in% names(ft)) {
    excl <- excl | (ft[["is_isotope"]] %in% TRUE)
  }
  if (length(isotope_names)) {
    excl <- excl | (name_all %in% isotope_names)
  }

  pairs_acc <- vector("list", n)
  pc <- 1L

  i <- 1L
  while (i <= n) {
    j <- i
    while (j < n && (rt_o[j + 1L] - rt_o[i]) <= rt_tol) {
      j <- j + 1L
    }

    idx <- ord[seq.int(i, j)]
    idx <- idx[!excl[idx]]
    if (length(idx) >= 2L) {
      nr <- length(idx)
      block_mz <- mz_all[idx]

      tri <- which(upper.tri(matrix(FALSE, nr, nr), diag = FALSE), arr.ind = TRUE)
      loc_i <- tri[, 1L]
      loc_j <- tri[, 2L]
      gi <- idx[loc_i]
      gj <- idx[loc_j]
      mi <- block_mz[loc_i]
      mj <- block_mz[loc_j]
      child <- ifelse(mi < mj, gi, gj)
      parent <- ifelse(mi < mj, gj, gi)
      child_mz <- pmin(mi, mj)
      parent_mz <- pmax(mi, mj)
      delta_mz <- parent_mz - child_mz

      err_mat <- abs(outer(delta_mz, nl_vals, "-"))
      hit_any <- rowSums(err_mat <= loss_tol_da) >= 1L

      if (any(hit_any)) {
        hit_idx <- which(hit_any)
        err_sub <- err_mat[hit_idx, , drop = FALSE]
        best_loss_col <- max.col(-err_sub, ties.method = "first")
        best_err <- err_sub[cbind(seq_along(hit_idx), best_loss_col)]
        best_loss_name <- nl_nms[best_loss_col]

        block <- abund[idx, , drop = FALSE]
        cmat <- suppressWarnings(stats::cor(t(block), use = "pairwise.complete.obs"))
        cmat[is.na(cmat)] <- 0
        diag(cmat) <- 1
        corv <- cmat[cbind(loc_i[hit_idx], loc_j[hit_idx])]

        ok <- corv >= corr_thresh
        if (any(ok)) {
          keep <- hit_idx[ok]
          pairs_acc[[pc]] <- data.table::data.table(
            child = child[keep],
            parent = parent[keep],
            cor = corv[ok],
            loss_name = best_loss_name[ok],
            loss_error = best_err[ok],
            parent_mz = parent_mz[keep]
          )
          pc <- pc + 1L
        }
      }
    }
    i <- i + 1L
  }

  is_isf <- rep(FALSE, n)
  parent_mz_out <- rep(NA_real_, n)
  isf_annotation <- rep(NA_character_, n)

  if (pc > 1L) {
    pairs_all <- data.table::rbindlist(pairs_acc[seq_len(pc - 1L)], use.names = TRUE)
    if (nrow(pairs_all)) {
      data.table::setorder(pairs_all, -cor, loss_error, -parent_mz)
      pairs_all <- unique(pairs_all, by = c("child", "parent"))
      pairs_all <- unique(pairs_all, by = "child")
      # Acyclic guard: a feature cannot be both a fragment and another
      # fragment's parent. Keep the parent role, drop the fragment role.
      repeat {
        both <- intersect(pairs_all$child, pairs_all$parent)
        if (!length(both)) break
        pairs_all <- pairs_all[!pairs_all$child %in% both]
      }

      is_isf[pairs_all$child] <- TRUE
      parent_mz_out[pairs_all$child] <- pairs_all$parent_mz
      isf_annotation[pairs_all$child] <- paste0(
        "ISF of ",
        name_all[pairs_all$parent],
        " (-",
        pairs_all$loss_name,
        ")"
      )
    }
  }

  ft[, is_isf := is_isf]
  ft[, parent_mz := parent_mz_out]
  ft[, isf_annotation := isf_annotation]

  if (!as_dt) {
    return(as.data.frame(ft, stringsAsFactors = FALSE))
  }
  ft[]
}


#' Infer parent compound candidates from in-source fragments
#'
#' When a compound undergoes in-source fragmentation (ISF), its fragment
#' *m/z* values can serve as diagnostic clues for identification.
#' `infer_isf_parents` searches a reference spectral database
#' ([NTSDatabase]) for compounds whose MS2 fragment *m/z* values match the
#' observed ISF fragment *m/z* values, and writes candidate identities to
#' both the ISF features and their parent features.
#'
#' The function groups ISF fragments by their `parent_mz`, then for each
#' fragment *m/z* performs a ppm-toleranced lookup in
#' `db@features_long`. Candidates are ranked by the number of distinct ISF
#' fragments that matched each database compound (more matching fragments
#' = stronger evidence).
#'
#' @param object An [NTSProject] with `is_isf` already populated via
#'   [flag_isf], or a `data.frame` / `data.table` containing at minimum
#'   the columns `name`, `mz`, `is_isf`, and `parent_mz`.
#' @param db An [NTSDatabase] with populated `@features_long` and
#'   `@metadata` slots.
#' @param mz_tol_ppm Numeric(1). m/z tolerance in ppm for matching ISF
#'   fragment *m/z* values to database fragment features. Default 20.
#' @param polarity `"POS"`, `"NEG"`, or `NULL`. If set, only database
#'   entries with matching polarity are considered. Default `NULL` (both
#'   polarities).
#'
#' @return An updated object of the same class as `object`, with three new
#'   columns on the annotation table:
#'   \describe{
#'     \item{`isf_candidate_ids`}{Semicolon-separated candidate compound
#'       identifiers (from `db@metadata$Lab_ID`).}
#'     \item{`isf_candidate_names`}{Semicolon-separated candidate compound
#'       names.}
#'     \item{`isf_candidate_n_fragments`}{Semicolon-separated match counts:
#'       how many distinct ISF fragments of this parent matched each
#'       candidate.}
#'   }
#'   For an [NTSProject], these are written to `@ms1_annotations`;
#'   `@feature_table` is unchanged.
#'
#' @export
infer_isf_parents <- function(object, db, mz_tol_ppm = 20, polarity = NULL) {
  if (!.is_reference_spectral_db(db)) {
    stop("`db` must be an NTSDatabase", call. = FALSE)
  }

  if (inherits(object, "NTSProject")) {
    proj <- object
    ann <- get_ms1_annotations(proj)
    if (!nrow(ann)) {
      stop("ms1_annotations is empty", call. = FALSE)
    }
    ann <- .infer_isf_parents_impl(ann, db, mz_tol_ppm, polarity)
    proj <- set_ms1_annotations(proj, .nts_ms1_sanitize_and_order(ann))
    methods::slot(proj, "log") <- .nts_log_append(
      methods::slot(proj, "log"),
      "infer_isf_parents",
      list(mz_tol_ppm = mz_tol_ppm, polarity = polarity)
    )
    return(proj)
  }

  .infer_isf_parents_impl(object, db, mz_tol_ppm, polarity)
}


#' @noRd
.infer_isf_parents_impl <- function(x, db, mz_tol_ppm, polarity) {
  if (!data.table::is.data.table(x)) {
    dt <- data.table::as.data.table(x)
    as_dt <- FALSE
  } else {
    dt <- data.table::copy(x)
    as_dt <- TRUE
  }

  dt[, isf_candidate_ids := NA_character_]
  dt[, isf_candidate_names := NA_character_]
  dt[, isf_candidate_n_fragments := NA_character_]

  need <- c("name", "mz", "is_isf", "parent_mz")
  miss <- setdiff(need, names(dt))
  if (length(miss)) {
    if (!as_dt) return(as.data.frame(dt, stringsAsFactors = FALSE))
    return(dt[])
  }

  isf <- dt[is_isf == TRUE & is.finite(parent_mz)]
  if (nrow(isf) == 0) {
    if (!as_dt) return(as.data.frame(dt, stringsAsFactors = FALSE))
    return(dt[])
  }

  feat <- data.table::copy(methods::slot(db, "features_long"))
  meta <- data.table::copy(methods::slot(db, "metadata"))

  if (nrow(feat) == 0 || nrow(meta) == 0) {
    warning("Database has empty features_long or metadata; returning without candidates.")
    if (!as_dt) return(as.data.frame(dt, stringsAsFactors = FALSE))
    return(dt[])
  }

  if (!is.null(polarity)) {
    polarity <- toupper(as.character(polarity)[1L])
    if (polarity %in% c("POS", "NEG") && "Polarity" %in% names(feat)) {
      feat <- feat[Polarity == polarity]
    }
  }

  if (nrow(feat) == 0) {
    if (!as_dt) return(as.data.frame(dt, stringsAsFactors = FALSE))
    return(dt[])
  }

  unique_parents <- unique(isf$parent_mz)
  parent_candidates <- vector("list", length(unique_parents))
  pi <- 0L

  for (pmz in unique_parents) {
    frags <- isf[parent_mz == pmz]
    frag_mz_vec <- frags$mz

    cand_list <- vector("list", length(frag_mz_vec))
    for (j in seq_along(frag_mz_vec)) {
      fmz <- frag_mz_vec[j]
      tol <- (mz_tol_ppm / 1e6) * fmz
      matched <- feat[abs(mz_value - fmz) <= tol]
      if (nrow(matched) > 0) {
        cand_list[[j]] <- unique(matched$Lab_ID)
      }
    }

    all_ids <- unlist(cand_list)
    if (length(all_ids) == 0) next

    id_counts <- sort(table(all_ids), decreasing = TRUE)
    top_ids <- names(id_counts)
    top_counts <- as.integer(id_counts)

    meta_sub <- meta[Lab_ID %in% top_ids]
    if (nrow(meta_sub) == 0) next

    name_vec <- vapply(seq_len(nrow(meta_sub)), function(k) {
      .nts_struct_lib_row_compound_name(meta_sub, k)
    }, character(1L))

    idx <- match(top_ids, meta_sub$Lab_ID)
    name_vec <- name_vec[idx]

    pi <- pi + 1L
    parent_candidates[[pi]] <- data.table::data.table(
      parent_mz = pmz,
      candidate_ids = paste(top_ids, collapse = "; "),
      candidate_names = paste(name_vec, collapse = "; "),
      candidate_n_fragments = paste(top_counts, collapse = "; ")
    )
  }

  if (pi == 0L) {
    if (!as_dt) return(as.data.frame(dt, stringsAsFactors = FALSE))
    return(dt[])
  }

  pc <- data.table::rbindlist(parent_candidates[seq_len(pi)], use.names = TRUE)

  # Write to ISF fragment rows
  isf_idx <- dt[, which(is_isf == TRUE & is.finite(parent_mz))]
  for (idx in isf_idx) {
    pmz <- dt$parent_mz[idx]
    match_row <- pc[parent_mz == pmz]
    if (nrow(match_row) > 0) {
      data.table::set(dt, idx, "isf_candidate_ids", match_row$candidate_ids[1L])
      data.table::set(dt, idx, "isf_candidate_names", match_row$candidate_names[1L])
      data.table::set(dt, idx, "isf_candidate_n_fragments", match_row$candidate_n_fragments[1L])
    }
  }

  # Write to parent feature rows (match by mz within ppm)
  for (pi2 in seq_len(nrow(pc))) {
    pmz <- pc$parent_mz[pi2]
    tol <- (mz_tol_ppm / 1e6) * pmz
    parent_idx <- dt[, which(!is_isf & abs(mz - pmz) <= tol)]
    if (length(parent_idx) > 0) {
      if (length(parent_idx) > 1L) {
        diffs <- abs(dt$mz[parent_idx] - pmz)
        parent_idx <- parent_idx[which.min(diffs)]
      }
      if (is.na(dt$isf_candidate_ids[parent_idx])) {
        data.table::set(dt, parent_idx, "isf_candidate_ids", pc$candidate_ids[pi2])
        data.table::set(dt, parent_idx, "isf_candidate_names", pc$candidate_names[pi2])
        data.table::set(dt, parent_idx, "isf_candidate_n_fragments", pc$candidate_n_fragments[pi2])
      }
    }
  }

  if (!as_dt) {
    return(as.data.frame(dt, stringsAsFactors = FALSE))
  }
  dt[]
}


#' @noRd
.nts_parse_exact_mass <- function(formula) {
  m <- calculate_exact_mass(as.character(formula)[1L])
  if (!is.finite(m) || m <= 0) {
    stop("could not parse `formula` to a positive exact mass: ", formula, call. = FALSE)
  }
  m
}


#' @noRd
.kmd_rt_audit_bad_series <- function(ft, base_nominal_span, min_spearman) {
  audit <- ft[!is.na(kmd_series_id), {
    ord <- order(mz, na.last = TRUE)
    mz_o <- as.numeric(mz[ord])
    rt_o <- as.numeric(rt[ord])
    ok <- is.finite(mz_o) & is.finite(rt_o)
    mz_o <- mz_o[ok]
    rt_o <- rt_o[ok]

    n_pt <- length(mz_o)
    mz_span <- if (n_pt) max(mz_o) - min(mz_o) else NA_real_
    span_ok <- is.finite(mz_span) && mz_span > (base_nominal_span * 0.8)
    if (n_pt < 2L || !span_ok) {
      .(
        keep = FALSE,
        spearman_r = NA_real_,
        span_ok = span_ok,
        mz_span = mz_span,
        n = n_pt
      )
    } else if (n_pt == 2L) {
      rt_diff <- rt_o[2L] - rt_o[1L]
      pair_ok <- rt_diff > -30
      .(
        keep = pair_ok,
        spearman_r = if (pair_ok) 1 else -1,
        span_ok = span_ok,
        mz_span = mz_span,
        n = n_pt
      )
    } else {
      r <- suppressWarnings(stats::cor(
        mz_o, rt_o,
        method = "spearman",
        use = "pairwise.complete.obs"
      ))
      if (!is.finite(r)) r <- 0
      .(
        keep = (r >= min_spearman),
        spearman_r = r,
        span_ok = span_ok,
        mz_span = mz_span,
        n = n_pt
      )
    }
  }, by = kmd_series_id]

  audit[keep == FALSE, as.character(kmd_series_id)]
}


#' @noRd
.kmd_rt_step_var <- function(sum_d, sum_d2, cnt) {
  if (!cnt || !is.finite(cnt) || cnt < 1L) {
    return(0)
  }
  m <- sum_d / cnt
  m2 <- sum_d2 / cnt
  v <- m2 - m * m
  if (!is.finite(v) || v < 0) 0 else v
}


#' `max_gap`: allow `steps` in `1:(max_gap+1)` repeat-unit masses between nodes
#' (missing homologues). `max_gap = 0` recovers strict one-step edges only.
#' @noRd
.kmd_edge_mass_ok <- function(mz_i, mz_j, exact_m, mz_tol_ppm, max_gap, mz_tol_da = NULL) {
  if (!is.finite(exact_m) || exact_m <= 0) {
    return(FALSE)
  }
  d <- mz_j - mz_i
  if (!is.finite(d) || d <= 0) {
    return(FALSE)
  }
  steps <- as.integer(round(d / exact_m))
  max_steps <- as.integer(max_gap) + 1L
  if (steps < 1L || steps > max_steps) {
    return(FALSE)
  }
  resid <- abs(d - as.numeric(steps) * exact_m)
  use_da <- !is.null(mz_tol_da) && length(mz_tol_da) &&
    is.finite(as.numeric(mz_tol_da)[[1L]]) && as.numeric(mz_tol_da)[[1L]] > 0
  if (use_da) {
    return(resid <= as.numeric(mz_tol_da)[[1L]])
  }
  ref <- max(abs(mz_i), abs(mz_j), 1e-12)
  resid / ref * 1e6 <= mz_tol_ppm
}


#' @noRd
.kmd_dp_better <- function(len_a, var_a, sd2_a, len_b, var_b, sd2_b) {
  if (len_a != len_b) {
    return(len_a > len_b)
  }
  if (abs(var_a - var_b) > 1e-15) {
    return(var_a < var_b)
  }
  sd2_a < sd2_b
}


#' @noRd
.kmd_extract_longest_path <- function(mz2, rt2, rid2, exact_m, mz_tol_ppm,
                                      max_gap, mz_tol_da, rt_monotonic, rt_max_step, min_series_size) {
  n <- length(mz2)
  if (n < min_series_size) {
    return(integer(0))
  }
  ord <- order(mz2, rt2, rid2, seq_len(n))
  len <- rep(1L, n)
  pred <- rep(NA_integer_, n)
  sum_d <- rep(0, n)
  sum_d2 <- rep(0, n)
  cnt_d <- rep(0L, n)

  for (pos in seq_len(n)) {
    j <- ord[pos]
    best_len <- 1L
    best_pred <- NA_integer_
    best_sd <- 0
    best_sd2 <- 0
    best_cnt <- 0L
    best_var <- 0

    if (pos > 1L) {
      for (ip in seq_len(pos - 1L)) {
        i <- ord[ip]
        if (mz2[i] >= mz2[j]) {
          next
        }
        if (!.kmd_edge_mass_ok(mz2[i], mz2[j], exact_m, mz_tol_ppm, max_gap, mz_tol_da)) {
          next
        }
        if (isTRUE(rt_monotonic) && !(rt2[j] > rt2[i])) {
          next
        }
        if (is.finite(rt_max_step) && (rt2[j] - rt2[i]) > rt_max_step) {
          next
        }

        nl <- len[i] + 1L
        nd <- rt2[j] - rt2[i]
        nsd <- sum_d[i] + nd
        nsd2 <- sum_d2[i] + nd * nd
        ncnt <- cnt_d[i] + 1L
        nv <- .kmd_rt_step_var(nsd, nsd2, ncnt)

        cur_var <- .kmd_rt_step_var(best_sd, best_sd2, best_cnt)
        if (.kmd_dp_better(nl, nv, nsd2, best_len, cur_var, best_sd2)) {
          best_len <- nl
          best_pred <- i
          best_sd <- nsd
          best_sd2 <- nsd2
          best_cnt <- ncnt
          best_var <- nv
        }
      }
    }

    len[j] <- best_len
    pred[j] <- best_pred
    sum_d[j] <- best_sd
    sum_d2[j] <- best_sd2
    cnt_d[j] <- best_cnt
  }

  elig <- which(len >= min_series_size)
  if (!length(elig)) {
    return(integer(0))
  }

  pick_var <- function(j) .kmd_rt_step_var(sum_d[j], sum_d2[j], cnt_d[j])
  jj <- elig[[1L]]
  for (k in seq_along(elig)) {
    j2 <- elig[k]
    vj <- pick_var(jj)
    v2 <- pick_var(j2)
    if (len[j2] > len[jj]) {
      jj <- j2
    } else if (len[j2] == len[jj]) {
      if (v2 < vj - 1e-15) {
        jj <- j2
      } else if (abs(v2 - vj) <= 1e-15 && rid2[j2] < rid2[jj]) {
        jj <- j2
      }
    }
  }

  path_loc <- integer(0)
  cur <- jj
  while (!is.na(cur)) {
    path_loc <- c(cur, path_loc)
    cur <- pred[cur]
  }
  path_loc
}


#' @noRd
.kmd_series_hits_for_base <- function(ft,
                                      eligible_idx,
                                      mz_num,
                                      exact_m,
                                      bf_tag,
                                      kmd_tol,
                                      mz_tol_ppm,
                                      max_gap,
                                      mz_tol_da,
                                      rt_monotonic,
                                      rt_max_step,
                                      min_series_size,
                                      rt_filter,
                                      min_spearman) {
  rt_dropped <- 0L
  nom_m <- round(exact_m)
  if (!is.finite(nom_m) || nom_m < 1) {
    stop("invalid nominal mass for base `", bf_tag, "`", call. = FALSE)
  }
  base_nominal_int <- as.integer(round(nom_m))

  mz_tol_ppm <- as.numeric(mz_tol_ppm)[1L]
  if (!is.finite(mz_tol_ppm) || mz_tol_ppm <= 0) {
    stop("mz_tol_ppm must be a positive finite number", call. = FALSE)
  }
  rt_max_step <- suppressWarnings(as.numeric(rt_max_step)[1L])
  if (!is.finite(rt_max_step) || rt_max_step <= 0) {
    rt_max_step <- Inf
  }
  max_gap <- suppressWarnings(as.integer(max_gap)[1L])
  if (!is.finite(max_gap) || max_gap < 0L) {
    stop("max_gap must be a non-negative integer", call. = FALSE)
  }

  km <- mz_num[eligible_idx] * (nom_m / exact_m)
  nkm <- round(km)
  kmd <- nkm - km
  mz_rem <- as.integer(round(mz_num[eligible_idx])) %% base_nominal_int

  wrk <- data.table::data.table(
    row_id = as.integer(eligible_idx),
    name = as.character(ft$name[eligible_idx]),
    mz_remainder = mz_rem,
    KMD = as.numeric(kmd),
    mz = as.numeric(mz_num[eligible_idx]),
    rt = as.numeric(ft$rt[eligible_idx])
  )
  data.table::setorder(wrk, mz_remainder, KMD, mz, row_id)
  wrk[, kmd_gap := c(0, diff(KMD)), by = mz_remainder]
  wrk[, sub_group := cumsum(kmd_gap > kmd_tol), by = mz_remainder]
  wrk[, n_grp := .N, by = .(mz_remainder, sub_group)]
  ok <- wrk[n_grp >= min_series_size]
  if (!nrow(ok)) {
    return(list(
      hits = data.table::data.table(name = character(0), kmd_series_id = character(0)),
      rt_dropped = 0L
    ))
  }

  series_counter <- 0L
  hits_parts <- list()

  u_grp <- unique(ok[, .(mz_remainder, sub_group)])
  data.table::setorder(u_grp, mz_remainder, sub_group)

  for (gr in seq_len(nrow(u_grp))) {
    mr <- u_grp$mz_remainder[gr]
    sg <- u_grp$sub_group[gr]
    grp <- ok[mz_remainder == mr & sub_group == sg]
    if (!nrow(grp)) {
      next
    }

    active <- rep(TRUE, nrow(grp))
    while (any(active)) {
      idx <- which(active)
      if (length(idx) < min_series_size) {
        active[] <- FALSE
        break
      }
      mz2 <- grp$mz[idx]
      rt2 <- grp$rt[idx]
      rid2 <- grp$row_id[idx]
      path_loc <- .kmd_extract_longest_path(
        mz2, rt2, rid2,
        exact_m = exact_m,
        mz_tol_ppm = mz_tol_ppm,
        max_gap = max_gap,
        mz_tol_da = mz_tol_da,
        rt_monotonic = rt_monotonic,
        rt_max_step = rt_max_step,
        min_series_size = min_series_size
      )
      if (!length(path_loc)) {
        active[] <- FALSE
        break
      }
      series_counter <- series_counter + 1L
      sid_chr <- sprintf("%s_R%02d_Series_%d", bf_tag, as.integer(mr), series_counter)
      gidx <- idx[path_loc]
      hits_parts[[length(hits_parts) + 1L]] <- data.table::data.table(
        name = grp$name[gidx],
        kmd_series_id = sid_chr
      )
      active[gidx] <- FALSE
    }
  }

  ok <- if (length(hits_parts)) {
    data.table::rbindlist(hits_parts, use.names = TRUE)
  } else {
    data.table::data.table(name = character(0), kmd_series_id = character(0))
  }

  if (isTRUE(rt_filter) && nrow(ok)) {
    tmp <- data.table::copy(ft)
    tmp[, kmd_series_id := NA_character_]
    tmp[ok, kmd_series_id := i.kmd_series_id, on = "name"]
    bad <- .kmd_rt_audit_bad_series(tmp, nom_m, min_spearman)
    rt_dropped <- length(bad)
    if (length(bad)) {
      ok <- ok[!kmd_series_id %in% bad]
    }
  }

  list(hits = ok[, .(name, kmd_series_id)], rt_dropped = rt_dropped)
}


#' @title Calculate Kendrick Mass Defect (KMD)
#'
#' @description
#' Computes Kendrick Mass (`KM`), nominal Kendrick mass (`NKM`), and Kendrick
#' Mass Defect (`KMD = NKM - KM`) from feature `mz` using the **first** entry in
#' `base_formulas`, and writes these to `@ms1_annotations` only (never to
#' `@feature_table`). Each formula in
#' `base_formulas` is parsed to an exact mass (see internal `.nts_parse_exact_mass`)
#' and integer nominal `round(exact)`; within each `mz_remainder` / `KMD` bin,
#' homologue chains are recovered as **longest paths** in a directed `m/z`-step
#' graph (each edge spans `1` to `max_gap + 1` repeat-unit mass increments within
#' `mz_tol_ppm`, optional RT monotonicity), then merged into one
#' `kmd_series_id` per feature (semicolon-separated). Rows with
#' `is_isotope == TRUE` are excluded from series assignment but still receive
#' `KM`/`NKM`/`KMD`/`standard_md`/`mz_remainder` (for the first base formula).
#'
#' @param object A [NTSProject], `data.frame`, or `data.table`.
#'   If `feature_table` contains logical `is_isotope` (from [flag_isotopes]),
#'   rows with `is_isotope == TRUE` are excluded from KMD clustering and never
#'   receive `kmd_series_id`, but still get `KM`/`NKM`/`KMD`/`mz_remainder`.
#' @param base_formulas Character vector of repeating-unit formulas, e.g.
#'   `c("CF2", "CF2O")`. Order matters: the first defines `KM`/`KMD`/`mz_remainder`.
#' @param md_range Two-element numeric vector defining the allowed standard mass
#'   defect range (`mz - round(mz)`) before KMD grouping (default
#'   `c(-0.15, 0.15)`).
#' @param kmd_tol Within each `mz_remainder` group, a new KMD sub-lane starts when
#'   `diff(KMD) > kmd_tol` after sorting by `KMD` then `m/z` (default `0.005`).
#'   Series are then found by a directed **m/z step graph** (see `mz_tol_ppm`)
#'   inside each sub-lane.
#' @param mz_tol_ppm Allowed relative error (ppm) when matching
#'   `m/z[j] - m/z[i]` to `steps *` repeat-unit exact mass (default `10`), used
#'   when `mz_tol_da` is `NULL` or non-positive.
#' @param mz_tol_da If finite and `> 0`, homologue edges use absolute Da residual
#'   `abs((mz[j]-mz[i]) - steps*exact_mass) <= mz_tol_da`; otherwise ppm rule applies.
#' @param max_gap Maximum number of missing homologue levels skipped by a single
#'   edge: allowed integer `steps` in `round((mz_j - mz_i) / exact_mass)` satisfy
#'   `1 <= steps <= max_gap + 1` (default `2`, i.e. up to three repeat units per edge).
#' @param rt_monotonic If `TRUE` (default), only allow edges `i -> j` when
#'   `rt[j] > rt[i]` (RP-LC PFAS ladders).
#' @param rt_max_step Optional maximum `rt[j] - rt[i]` per edge (seconds); `Inf`
#'   disables (default).
#' @param min_series_size Minimum chain length to assign `kmd_series_id`
#'   (default `3`).
#' @param min_length Optional alias for `min_series_size` (same meaning).
#' @param rt_filter Logical; if `TRUE`, retain only KMD series that pass a
#'   physical RT-trend audit (default `TRUE`).
#' @param min_spearman Minimum Spearman rank correlation required for `n >= 3`
#'   series in the `mz-rt` trend audit (default `0.5`).
#' @param min_pearson_r,min_r2 Backward-compatible aliases. If supplied they are
#'   converted into `min_spearman`.
#' @param precision,series_digits,series_min_n Backward-compatible aliases from
#'   the older binned KMD engine; if supplied they are translated into
#'   `kmd_tol`/`min_series_size`.
#'
#' @return Updated [NTSProject] or table with `standard_md`, `KM`, `NKM`, `KMD`,
#'   `mz_remainder`, and `kmd_series_id`.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' ft <- data.frame(
#'   name = "FT00001",
#'   mz = 300,
#'   rt = 10,
#'   stringsAsFactors = FALSE
#' )
#' out <- suppressMessages(calculate_kmd(ft, base_formulas = "CF2", min_series_size = 1L))
#' stopifnot(all(c("KM", "KMD") %in% names(out)))
#' @export
calculate_kmd <- function(object,
                          base_formulas = "CF2",
                          md_range = c(-0.15, 0.15),
                          kmd_tol = 0.005,
                          mz_tol_ppm = 10,
                          max_gap = 2,
                          mz_tol_da = NULL,
                          rt_monotonic = TRUE,
                          rt_max_step = Inf,
                          min_series_size = 3,
                          min_length = NULL,
                          rt_filter = TRUE,
                          min_spearman = 0.5,
                          min_pearson_r = NULL,
                          min_r2 = NULL,
                          precision = NULL,
                          series_digits = NULL,
                          series_min_n = NULL) {
  if (!is.null(mz_tol_da)) {
    mz_tol_da <- suppressWarnings(as.numeric(mz_tol_da)[1L])
    if (!is.finite(mz_tol_da) || mz_tol_da <= 0) {
      mz_tol_da <- NULL
    }
  }
  if (!is.null(min_length)) {
    min_series_size <- max(1L, as.integer(min_length)[1L])
  }
  if (!is.null(precision)) {
    kmd_tol <- 0.5 * 10^(-as.integer(precision)[1L])
  }
  if (!is.null(series_digits)) {
    kmd_tol <- 0.5 * 10^(-as.integer(series_digits)[1L])
  }
  if (!is.null(series_min_n)) {
    min_series_size <- series_min_n
  }
  if (!is.null(min_pearson_r)) {
    min_spearman <- as.numeric(min_pearson_r)[1L]
  }
  if (!is.null(min_r2)) {
    min_spearman <- sqrt(as.numeric(min_r2)[1L])
  }
  if (inherits(object, "NTSProject")) {
    proj <- object
    ft <- get_feature_table(proj)
    if (!nrow(ft)) {
      stop("feature_table is empty", call. = FALSE)
    }
    ms1 <- get_ms1_annotations(proj)
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("data.table is required", call. = FALSE)
    }
    work <- data.table::copy(data.table::as.data.table(ft))
    keep_nm <- intersect(names(work), c("name", "mz", "rt"))
    work <- work[, keep_nm, with = FALSE]
    if (nrow(ms1) && "name" %in% names(ms1) && "is_isotope" %in% names(ms1)) {
      iso <- ms1[, c("name", "is_isotope"), drop = FALSE]
      work <- merge(
        as.data.frame(work, stringsAsFactors = FALSE),
        iso,
        by = "name",
        all.x = TRUE,
        sort = FALSE
      )
      work <- data.table::as.data.table(work)
      work[is.na(is_isotope), is_isotope := FALSE]
    } else {
      work[, is_isotope := FALSE]
    }
    kmd_sub <- .nts_engine_kmd(
      work,
      base_formulas = base_formulas,
      md_range = md_range,
      kmd_tol = kmd_tol,
      mz_tol_ppm = mz_tol_ppm,
      max_gap = max_gap,
      mz_tol_da = mz_tol_da,
      rt_monotonic = rt_monotonic,
      rt_max_step = rt_max_step,
      min_series_size = min_series_size,
      rt_filter = rt_filter,
      min_spearman = min_spearman
    )

    ann <- get_ms1_annotations(proj)
    if (nrow(ann) && "name" %in% names(ann) && nrow(kmd_sub) && "name" %in% names(kmd_sub)) {
      ann_dt <- data.table::as.data.table(ann)
      engine_cols <- setdiff(names(kmd_sub), "name")
      drop_here <- intersect(engine_cols, names(ann_dt))
      for (cl in drop_here) {
        ann_dt[, (cl) := NULL]
      }
      merged <- kmd_sub[ann_dt, on = "name"]
      proj <- set_ms1_annotations(
        proj,
        as.data.frame(
          merged,
          stringsAsFactors = FALSE
        )
      )
    }

    proj <- set_ms1_annotations(
      proj,
      .nts_ms1_sanitize_and_order(get_ms1_annotations(proj))
    )
    methods::slot(proj, "log") <- .nts_log_append(
      methods::slot(proj, "log"),
      "calculate_kmd",
      list(
        base_formulas = as.character(base_formulas),
        md_range = md_range,
        kmd_tol = kmd_tol,
        mz_tol_ppm = mz_tol_ppm,
        max_gap = max_gap,
        mz_tol_da = mz_tol_da,
        rt_monotonic = rt_monotonic,
        rt_max_step = rt_max_step,
        min_series_size = min_series_size,
        rt_filter = rt_filter,
        min_spearman = min_spearman
      )
    )
    return(proj)
  }

  .calculate_kmd_df(
    object,
    base_formulas = base_formulas,
    md_range = md_range,
    kmd_tol = kmd_tol,
    mz_tol_ppm = mz_tol_ppm,
    max_gap = max_gap,
    mz_tol_da = mz_tol_da,
    rt_monotonic = rt_monotonic,
    rt_max_step = rt_max_step,
    min_series_size = min_series_size,
    rt_filter = rt_filter,
    min_spearman = min_spearman
  )
}


#' @noRd
.calculate_kmd_df <- function(feature_table,
                              base_formulas = "CF2",
                              md_range = c(-0.15, 0.15),
                              kmd_tol = 0.005,
                              mz_tol_ppm = 10,
                              max_gap = 2,
                              mz_tol_da = NULL,
                              rt_monotonic = TRUE,
                              rt_max_step = Inf,
                              min_series_size = 3,
                              rt_filter = TRUE,
                              min_spearman = 0.5) {
  if (!data.table::is.data.table(feature_table)) {
    ft <- data.table::copy(data.table::as.data.table(feature_table))
    as_dt <- FALSE
  } else {
    ft <- data.table::copy(feature_table)
    as_dt <- TRUE
  }

  stopifnot("mz" %in% names(ft), "rt" %in% names(ft))

  base_formulas <- as.character(base_formulas)
  if (!length(base_formulas)) {
    stop("base_formulas must be non-empty", call. = FALSE)
  }
  base_exacts <- vapply(base_formulas, .nts_parse_exact_mass, numeric(1L), USE.NAMES = TRUE)

  md_range <- sort(as.numeric(md_range))
  if (length(md_range) != 2L || any(!is.finite(md_range))) {
    stop("md_range must contain two finite numeric bounds", call. = FALSE)
  }
  kmd_tol <- as.numeric(kmd_tol)[1L]
  if (!is.finite(kmd_tol) || kmd_tol <= 0) {
    stop("kmd_tol must be a positive finite number", call. = FALSE)
  }
  mz_tol_ppm <- as.numeric(mz_tol_ppm)[1L]
  if (!is.finite(mz_tol_ppm) || mz_tol_ppm <= 0) {
    stop("mz_tol_ppm must be a positive finite number", call. = FALSE)
  }
  min_series_size <- max(1L, as.integer(min_series_size)[1L])
  max_gap <- suppressWarnings(as.integer(max_gap)[1L])
  if (!is.finite(max_gap) || max_gap < 0L) {
    stop("max_gap must be a non-negative integer", call. = FALSE)
  }
  if (!is.null(mz_tol_da)) {
    mz_tol_da <- suppressWarnings(as.numeric(mz_tol_da)[1L])
    if (!is.finite(mz_tol_da) || mz_tol_da <= 0) {
      mz_tol_da <- NULL
    }
  }
  min_spearman <- as.numeric(min_spearman)[1L]
  if (!is.finite(min_spearman) || min_spearman < -1 || min_spearman > 1) {
    stop("min_spearman must be a finite number in [-1, 1]", call. = FALSE)
  }

  if (!"name" %in% names(ft)) {
    ft[, name := sprintf("FT%05d", .I)]
  }

  mz_num <- as.numeric(ft$mz)
  standard_md <- mz_num - round(mz_num)
  bf0 <- base_formulas[[1L]]
  ex0 <- unname(base_exacts[[bf0]])
  nm0 <- round(ex0)
  if (!is.finite(nm0) || nm0 < 1) {
    stop("invalid nominal mass for first base formula `", bf0, "`", call. = FALSE)
  }
  nom0_int <- as.integer(nm0)
  km <- mz_num * (nm0 / ex0)
  nkm <- round(km)
  kmd <- nkm - km
  mz_remainder <- as.integer(round(mz_num)) %% nom0_int

  ft[, standard_md := standard_md]
  ft[, KM := km]
  ft[, NKM := nkm]
  ft[, KMD := kmd]
  ft[, mz_remainder := mz_remainder]
  ft[, kmd_series_id := NA_character_]

  iso_eligible <- rep(TRUE, nrow(ft))
  if ("is_isotope" %in% names(ft)) {
    iso <- ft$is_isotope
    iso_eligible <- is.na(iso) | !as.logical(iso)
  }

  eligible_idx <- which(
    is.finite(standard_md) &
      is.finite(mz_num) &
      standard_md >= md_range[1L] &
      standard_md <= md_range[2L] &
      iso_eligible
  )

  series_list <- vector("list", length(base_formulas))
  names(series_list) <- base_formulas
  dropped_series <- 0L

  for (bf in base_formulas) {
    pack <- .kmd_series_hits_for_base(
      ft,
      eligible_idx,
      mz_num,
      unname(base_exacts[[bf]]),
      bf,
      kmd_tol,
      mz_tol_ppm,
      max_gap,
      mz_tol_da,
      rt_monotonic,
      rt_max_step,
      min_series_size,
      rt_filter,
      min_spearman
    )
    dropped_series <- dropped_series + as.integer(pack$rt_dropped)
    series_list[[bf]] <- pack$hits
  }

  all_hits <- data.table::rbindlist(series_list, use.names = TRUE, fill = TRUE)
  if (nrow(all_hits)) {
    agg <- all_hits[, .(kmd_series_id = paste(unique(kmd_series_id), collapse = "; ")), by = name]
    ft[agg, kmd_series_id := i.kmd_series_id, on = "name"]
  }

  n_labeled <- sum(!is.na(ft$kmd_series_id))
  message(
    "> Multi-base KMD (",
    paste(base_formulas, collapse = ", "),
    "): ",
    n_labeled,
    " feature row(s) assigned kmd_series_id."
  )

  if (isTRUE(rt_filter)) {
    # Per-base series IDs are audited inside each base; here we tally dropped series (may double-count across bases).
    message(
      "[ok] RT trend filter (per base): dropped ",
      dropped_series,
      " spurious series (insufficient span or abnormal trend)."
    )
  }

  if (!as_dt) {
    return(as.data.frame(ft, stringsAsFactors = FALSE))
  }
  ft[]
}


#' Pure KMD engine facade: never mutates the caller's `dt` object. All Kendrick
#' columns are computed on an internal copy via `.calculate_kmd_df`; returns only
#' the join key `name` plus qualitative KMD columns for safe merging into
#' `@ms1_annotations`.
#'
#' @keywords internal
#' @noRd
.nts_engine_kmd <- function(dt,
                            base_formulas = "CF2",
                            md_range = c(-0.15, 0.15),
                            kmd_tol = 0.005,
                            mz_tol_ppm = 10,
                            max_gap = 2,
                            mz_tol_da = NULL,
                            rt_monotonic = TRUE,
                            rt_max_step = Inf,
                            min_series_size = 3,
                            rt_filter = TRUE,
                            min_spearman = 0.5) {
  full <- .calculate_kmd_df(
    dt,
    base_formulas = base_formulas,
    md_range = md_range,
    kmd_tol = kmd_tol,
    mz_tol_ppm = mz_tol_ppm,
    max_gap = max_gap,
    mz_tol_da = mz_tol_da,
    rt_monotonic = rt_monotonic,
    rt_max_step = rt_max_step,
    min_series_size = min_series_size,
    rt_filter = rt_filter,
    min_spearman = min_spearman
  )
  full <- data.table::as.data.table(full)
  keep <- intersect(
    c("name", "standard_md", "KM", "NKM", "KMD", "mz_remainder", "kmd_series_id"),
    names(full)
  )
  if (!length(keep) || !"name" %in% keep) {
    stop(".nts_engine_kmd(): result must contain column `name`", call. = FALSE)
  }
  full[, ..keep]
}


#' Predict organic percentage from retention time
#'
#' Uses linear interpolation (via [stats::approxfun]) to predict the organic
#' solvent percentage at a given retention time, based on a reference
#' RT-to-organic-percentage mapping from a known LC gradient.
#' The same LC gradient (e.g., RP-C18) can be reused across datasets
#' with the same chromatographic method.
#'
#' @param rt Numeric vector of retention times (seconds) to predict.
#' @param ref_rt Numeric vector of reference retention times (seconds).
#' @param ref_org Numeric vector of reference organic percentages, same
#'   length as `ref_rt`.
#' @param method Interpolation method: `"linear"` (default) or `"constant"`
#'   uses [stats::approxfun]; `"spline"` uses [stats::splinefun] with
#'   method `"monoH.FC"`.
#' @param rule Extrapolation rule for `"linear"` and `"constant"` methods:
#'   `1` returns `NA` outside the reference range; `2` (default) extends
#'   the boundary values.
#'
#' @return Numeric vector of predicted organic percentages, same length as `rt`.
#' @export
#'
#' @examples
#' ref_rt <- c(0, 150, 300, 450, 600)
#' ref_org <- c(5, 27.5, 50, 72.5, 95)
#' predict_organic_percentage(c(100, 250, 700), ref_rt, ref_org)
predict_organic_percentage <- function(rt,
                                       ref_rt,
                                       ref_org,
                                       method = c("linear", "constant", "spline"),
                                       rule = 2) {
  if (length(rt) == 0L) {
    return(numeric(0))
  }
  if (!is.numeric(rt)) {
    stop("`rt` must be numeric", call. = FALSE)
  }

  stopifnot(
    "`ref_rt` must be numeric" = is.numeric(ref_rt),
    "`ref_org` must be numeric" = is.numeric(ref_org),
    "`ref_rt` and `ref_org` must have the same length" =
      length(ref_rt) == length(ref_org),
    "`ref_rt` must not be empty" = length(ref_rt) > 0L
  )

  method <- match.arg(method)

  if (!is.numeric(rule) || length(rule) != 1L || !rule %in% c(1, 2)) {
    stop("`rule` must be 1 or 2", call. = FALSE)
  }
  rule <- as.integer(rule)

  ok <- is.finite(ref_rt) & is.finite(ref_org)
  ref_rt <- ref_rt[ok]
  ref_org <- ref_org[ok]

  if (length(ref_rt) < 1L) {
    stop("Need at least 1 finite reference point", call. = FALSE)
  }
  if (length(ref_rt) < 2L && method != "constant") {
    stop(
      "Need at least 2 finite reference points for method '", method, "'",
      call. = FALSE
    )
  }

  if (anyDuplicated(ref_rt)) {
    agg <- tapply(ref_org, ref_rt, mean, na.rm = TRUE)
    ref_rt <- as.numeric(names(agg))
    ref_org <- as.numeric(agg)
  }

  ord <- order(ref_rt)
  ref_rt <- ref_rt[ord]
  ref_org <- ref_org[ord]

  if (method == "spline") {
    fn <- tryCatch(
      stats::splinefun(ref_rt, ref_org, method = "monoH.FC"),
      error = function(e) {
        stop("spline interpolation failed: ", e$message, call. = FALSE)
      }
    )
    pred <- fn(rt)
    lo <- ref_rt[1L]
    hi <- ref_rt[length(ref_rt)]
    if (rule == 2) {
      pred[!is.na(rt) & rt < lo] <- ref_org[1L]
      pred[!is.na(rt) & rt > hi] <- ref_org[length(ref_org)]
    } else {
      pred[!is.na(rt) & (rt < lo | rt > hi)] <- NA_real_
    }
    return(as.numeric(pred))
  }

  fn <- stats::approxfun(ref_rt, ref_org, method = method, rule = rule)
  as.numeric(fn(rt))
}
