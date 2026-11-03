# ============================================================================
# MS2 CF2 fragment-ladder / homologous-series detection (library-free)
#
# Design: docs/analysis/NTSminer_MS2_CF2碎片阶梯检测设计方案_v0.1.22.md  (S1-S4)
#
# Two entry points:
#   * find_pfas_fragment_series()      -- pure, single-spectrum, side-effect free
#   * detect_pfas_fragment_series()    -- project-level, writes @ms2_annotations
#
# The neutral repeat unit CF2 has exact mass 49.9968064 Da
# (C 12.000000 + 2 x F 18.9984032).  A "ladder" is a set of fragment peaks on
# the m/z axis whose consecutive members differ by an integer multiple k x CF2.
#
# Decisions locked 2026-10-10 (design §11):
#   #1 main-pipeline switch default ON      (see 06_scoring_engine.R)
#   #2 confidence: structural guardrail only (see 07_confidence_levels.R)
#   #3 network: node attributes only         (no new edge type)
#   #4 ladder allows missing rungs           (mode = "loose" is the default)
#   #5 ship the P1-2 preset subclass rule set here
# ============================================================================

#' Neutral exact mass of the CF2 repeat unit (Da)
#'
#' `C (12.000000) + 2 x F (18.9984032) = 49.9968064`.
#' @return A numeric scalar.
#' @keywords internal
#' @noRd
.nts_pfas_cf2_da <- 49.9968064

#' Theoretical m/z of the perfluoroalkyl anion CnF(2n+1)-
#'
#' The series is strictly arithmetic with step CF2: the anchor CF3- has exact
#' m/z 68.9957596 (see the packaged `extdata/tables/diagnostic_fragments.csv`), so every higher member
#' is an exact multiple of `cf2_da` above it.  No isotope-table lookup (and
#' hence no `Rdisop` dependency) is required.
#'
#' @param n Integer vector, n >= 1 (CF3- is n = 1).
#' @param cf2_da CF2 repeat mass.
#' @return Numeric vector of monoisotopic m/z (negative mode).
#' @keywords internal
#' @noRd
.nts_pfas_cnf2n1_mz <- function(n, cf2_da = .nts_pfas_cf2_da) {
  n <- suppressWarnings(as.numeric(n))
  68.9957596 + (n - 1) * cf2_da
}

#' @keywords internal
#' @noRd
.nts_pfas_cnf2n1_label <- function(n) {
  n <- suppressWarnings(as.integer(n))
  sprintf("C%dF%d-", n, 2L * n + 1L)
}


#' Preset PFAS sub-class diagnostic rules (P1-2)
#'
#' A small, chemistry-curated rule table that maps diagnostic fragment masses
#' to PFAS sub-classes.  Schema follows the fragment/neutral-loss/type
#' convention used elsewhere in the package
#' (`R/01_subclass_fragment_catalog.R`):
#'
#' \describe{
#'   \item{`subclass`}{`"PFCA"`, `"PFSA"`, `"FTS"` (extensible).}
#'   \item{`ion_formula`}{Diagnostic ion label, e.g. `"C3F7-"`, `"SO3-"`.}
#'   \item{`ion_mz`}{Monoisotopic m/z (negative mode); `NA` for the
#'     `type = "series"` row.}
#'   \item{`polarity`}{`"neg"` / `"pos"`.}
#'   \item{`type`}{`"fragment"` or `"series"`. A `"series"` row matches when a
#'     CF2 ladder is present, letting rules require a ladder explicitly.}
#'   \item{`required`}{Logical. Every `required = TRUE` row of a sub-class must
#'     match **and** at least one row of that sub-class must match for the
#'     sub-class hint to be emitted.}
#'   \item{`label`}{Human-readable description.}
#' }
#'
#' Rules (`neg`): PFCA fires on any CnF(2n+1)- member; PFSA on any sulfonate /
#' fluorosulfonate / triflate head-group ion; FTS on a sulfonate head-group
#' **plus** a CF2 ladder (`type = "series"`, `required = TRUE`).  The
#' ether-PFAS rule set is deferred to a later stage (design §4.4).
#'
#' @param polarity Optional `"neg"` / `"pos"` filter. `NULL` (default) returns
#'   every row; `"pos"` yields a 0-row table (PFAS sub-class rules are
#'   negative-mode).
#' @param series_n_max Highest n for the CnF(2n+1)- members (default 8).
#' @return A `data.frame`.
#' @seealso [find_pfas_fragment_series]
#' @examples
#' head(pfas_diagnostic_rules())
#' @export
pfas_diagnostic_rules <- function(polarity = NULL, series_n_max = 8L) {
  n_max <- suppressWarnings(as.integer(series_n_max)[1L])
  if (!is.finite(n_max) || n_max < 1L) n_max <- 8L
  ns <- seq_len(n_max)
  cn_mz <- .nts_pfas_cnf2n1_mz(ns)
  cn_lab <- .nts_pfas_cnf2n1_label(ns)

  out <- data.frame(
    subclass = c(
      rep("PFCA", n_max),
      "PFSA", "PFSA", "PFSA",
      "FTS", "FTS", "FTS"
    ),
    ion_formula = c(
      cn_lab,
      "SO3-", "FSO3-", "CF3SO3-",
      "SO3-", "SO3H-", "CF2_ladder"
    ),
    ion_mz = c(
      cn_mz,
      79.957367, 98.9557702, 148.9525766,
      79.957367, 80.9651916, NA_real_
    ),
    polarity = "neg",
    type = c(rep("fragment", n_max + 5L), "series"),
    required = c(rep(FALSE, n_max + 5L), TRUE),
    label = c(
      cn_lab,
      "sulfonate head group", "fluorosulfonate", "triflate",
      "sulfonate head group", "sulfonic acid", "CF2 ladder"
    ),
    stringsAsFactors = FALSE
  )

  if (!is.null(polarity) && length(polarity)) {
    pol <- tolower(trimws(as.character(polarity)[1L]))
    if (pol %in% c("negative", "neg", "-")) pol <- "neg"
    if (pol %in% c("positive", "pos", "+")) pol <- "pos"
    out <- out[tolower(out$polarity) == pol, , drop = FALSE]
    rownames(out) <- NULL
  }
  out
}


#' Coerce a spectrum argument to a two-column data.frame
#' @keywords internal
#' @noRd
.nts_pfas_spectrum_df <- function(spectrum) {
  if (is.null(spectrum)) {
    return(data.frame(mz = numeric(0), intensity = numeric(0)))
  }
  if (is.data.frame(spectrum)) {
    nm <- names(spectrum)
    mz_col <- intersect(c("mz", "Mz", "MZ", "precursorMZ"), nm)
    in_col <- intersect(c("intensity", "Intensity", "int", "Int"), nm)
    if (!length(mz_col) || !length(in_col)) {
      stop(
        "find_pfas_fragment_series(): `spectrum` data.frame needs mz and intensity columns.",
        call. = FALSE
      )
    }
    return(data.frame(
      mz = as.numeric(spectrum[[mz_col[1L]]]),
      intensity = as.numeric(spectrum[[in_col[1L]]])
    ))
  }
  if (is.list(spectrum) && !is.null(spectrum$mz) && !is.null(spectrum$intensity)) {
    return(data.frame(
      mz = as.numeric(spectrum$mz),
      intensity = as.numeric(spectrum$intensity)
    ))
  }
  if (is.matrix(spectrum) && ncol(spectrum) >= 2L) {
    return(data.frame(
      mz = as.numeric(spectrum[, 1L]),
      intensity = as.numeric(spectrum[, 2L])
    ))
  }
  stop(
    "find_pfas_fragment_series(): unsupported `spectrum` (expected data.frame / ",
    "list(mz=, intensity=) / 2-column matrix).",
    call. = FALSE
  )
}


#' Longest CF2 ladder on a sorted m/z axis (DAG longest path)
#'
#' Consecutive members may differ by `k x CF2` with `k` in `[1, max_gap_units]`
#' in `"loose"` mode (missing rungs allowed, decision #4) or `k = 1` only in
#' `"strict"` mode.  Ties are broken towards the larger total span, then the
#' lower m/z start.
#'
#' @param mz Numeric vector of peak m/z (already intensity-filtered).
#' @param cf2_da,cf2_tol_da Repeat-unit mass and per-step tolerance (Da).
#' @param max_gap_units Maximum k for a single hop.
#' @param mode `"loose"` or `"strict"`.
#' @return `list(len, span, mz, mz_chr, gaps)`; `len = 0` when no chain of at
#'   least two members exists.
#' @keywords internal
#' @noRd
.nts_pfas_ladder_chain <- function(mz, cf2_da, cf2_tol_da,
                                   max_gap_units, mode = c("loose", "strict")) {
  mode <- match.arg(mode)
  empty <- list(
    len = 0L, span = 0L, mz = numeric(0), mz_chr = NA_character_,
    gaps = integer(0)
  )
  mz <- mz[is.finite(mz)]
  n <- length(mz)
  if (n < 2L) {
    return(empty)
  }
  mzs <- sort(mz)
  hi_k <- if (identical(mode, "strict")) 1L else max(1L, as.integer(max_gap_units))

  best_len  <- rep(1L, n)
  best_span <- rep(0L, n)
  prev      <- rep(NA_integer_, n)
  gap_k     <- rep(NA_integer_, n)

  for (j in seq.int(2L, n)) {
    for (i in seq.int(j - 1L, 1L)) {
      d <- mzs[j] - mzs[i]
      if (!is.finite(d) || d <= 0) next
      k <- round(d / cf2_da)
      if (k < 1L || k > hi_k) next
      if (abs(d - k * cf2_da) > cf2_tol_da) next
      cand_len  <- best_len[i] + 1L
      cand_span <- best_span[i] + k
      if (cand_len > best_len[j] ||
          (cand_len == best_len[j] && cand_span > best_span[j])) {
        best_len[j]  <- cand_len
        best_span[j] <- cand_span
        prev[j]      <- i
        gap_k[j]     <- as.integer(k)
      }
    }
  }

  if (max(best_len) < 2L) {
    return(empty)
  }
  key <- best_len * 1e6 + best_span
  j <- which.max(key)

  nodes_rev <- j
  gaps_rev  <- integer(0)
  cur <- j
  while (!is.na(prev[cur])) {
    gaps_rev  <- c(gaps_rev, gap_k[cur])
    cur       <- prev[cur]
    nodes_rev <- c(nodes_rev, cur)
  }
  nodes <- rev(nodes_rev)
  gaps  <- rev(as.integer(gaps_rev))
  members <- mzs[nodes]

  list(
    len = length(members),
    span = as.integer(sum(gaps)),
    mz = members,
    mz_chr = paste(format(members, digits = 6, trim = TRUE), collapse = "; "),
    gaps = gaps
  )
}


#' Match peaks against the packaged diagnostic-fragment expert table
#' @keywords internal
#' @noRd
.nts_pfas_diag_hits <- function(mz, polarity, diag_csv = NULL, mz_tol_ppm = 10) {
  dt <- .nts_pfas_diag_table(diag_csv)
  if (is.null(dt) || !nrow(dt)) {
    return(character(0))
  }
  pol <- .nts_pfas_norm_polarity(polarity)
  pol_col <- intersect(c("polarity", "Polarity"), names(dt))
  if (length(pol_col)) {
    keep <- tolower(trimws(as.character(dt[[pol_col[1L]]]))) == pol
    keep[is.na(keep)] <- FALSE
    dt <- dt[keep, , drop = FALSE]
  }
  if (!nrow(dt)) {
    return(character(0))
  }
  mass <- suppressWarnings(as.numeric(dt[["fragment_mass"]]))
  form <- trimws(as.character(dt[["fragment_formula"]]))
  ok <- is.finite(mass) & nzchar(form)
  dt <- dt[ok, , drop = FALSE]
  mass <- mass[ok]
  form <- form[ok]
  if (!length(mass) || !length(mz)) {
    return(character(0))
  }
  mz <- mz[is.finite(mz)]
  ppm <- mz_tol_ppm * 1e-6
  hit <- vapply(mass, function(m) {
    any(abs(mz - m) <= m * ppm)
  }, logical(1L))
  unique(form[hit])
}


#' Read / normalise a diagnostic fragment table
#' @keywords internal
#' @noRd
.nts_pfas_diag_table <- function(diag_csv = NULL) {
  if (is.data.frame(diag_csv)) {
    return(diag_csv)
  }
  path <- diag_csv
  if (is.null(path) || !length(path) || !nzchar(as.character(path)[1L])) {
    fn <- tryCatch(
      get("default_diagnostic_fragments_csv", mode = "function"),
      error = function(e) NULL
    )
    path <- if (is.null(fn)) NA_character_ else fn()
  }
  path <- as.character(path)[1L]
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  out <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  out
}


#' @keywords internal
#' @noRd
.nts_pfas_norm_polarity <- function(polarity) {
  p <- tolower(trimws(as.character(polarity)[1L]))
  if (is.na(p) || !nzchar(p)) p <- "neg"
  if (p %in% c("negative", "neg", "-", "n")) p <- "neg"
  if (p %in% c("positive", "pos", "+", "p")) p <- "pos"
  p
}


#' @keywords internal
#' @noRd
.nts_pfas_rules_table <- function(rules = NULL) {
  if (is.data.frame(rules)) {
    return(rules)
  }
  pfas_diagnostic_rules()
}


#' Apply the sub-class rule table directly to the peak list
#'
#' Rules are matched on their own `ion_mz` (ppm tolerance) rather than on the
#' diagnostic-table hits, so `pfas_diagnostic_rules()` is self-contained and
#' does not depend on the packaged CSV being resolvable.  `type = "series"`
#' rows match on the presence of a CF2 ladder instead.
#'
#' @param mz Numeric vector of (filtered) peak m/z.
#' @param ladder_present Logical; is a CF2 ladder of at least `min_len` members
#'   present?
#' @param rules Sub-class rule `data.frame`.
#' @param mz_tol_ppm Symmetric ppm tolerance.
#' @return Character vector of sub-class hints (possibly empty).
#' @keywords internal
#' @noRd
.nts_pfas_subclass_hints <- function(mz, ladder_present, rules, mz_tol_ppm = 10) {
  if (is.null(rules) || !nrow(rules) ||
      !all(c("subclass", "ion_formula") %in% names(rules))) {
    return(character(0))
  }
  mz0 <- mz[is.finite(mz)]
  ppm <- mz_tol_ppm * 1e-6
  ty <- if ("type" %in% names(rules)) {
    tolower(trimws(as.character(rules$type)))
  } else {
    rep("fragment", nrow(rules))
  }
  ty[is.na(ty) | !nzchar(ty)] <- "fragment"
  req <- if ("required" %in% names(rules)) {
    as.logical(rules$required)
  } else {
    rep(FALSE, nrow(rules))
  }
  req[is.na(req)] <- FALSE

  # A rule "matches" when its own ion_mz is present in the spectrum, or (for
  # series rows) when a ladder is present.
  rule_mz <- suppressWarnings(as.numeric(rules$ion_mz))
  label_hit <- vapply(rule_mz, function(t) {
    if (!is.finite(t) || !length(mz0)) return(FALSE)
    any(abs(mz0 - t) <= t * ppm)
  }, logical(1L))
  rule_hit <- ifelse(ty == "series", ladder_present, label_hit)

  subs <- unique(trimws(as.character(rules$subclass)))
  out <- character(0)
  for (s in subs) {
    if (is.na(s) || !nzchar(s)) next
    idx <- which(trimws(as.character(rules$subclass)) == s)
    hit_s <- rule_hit[idx]
    req_s <- req[idx]
    ty_s <- ty[idx]
    # A sub-class with fragment rows needs at least one of them to match --
    # a lone `type = "series"` row (e.g. FTS's ladder) must not fire on its own.
    frag_idx <- which(ty_s != "series")
    if (length(frag_idx) && !any(hit_s[frag_idx])) next
    if (any(req_s) && !all(hit_s[req_s])) next
    if (!any(hit_s)) next
    out <- c(out, s)
  }
  out
}


#' Detect CF2 fragment ladders and PFAS series in a single MS2 spectrum
#'
#' Library-free structural evidence for PFAS: finds the longest fragment ladder
#' whose consecutive members differ by an integer multiple of CF2
#' (49.9968064 Da), matches the perfluoroalkyl anion series
#' CnF(2n+1)-, matches the packaged diagnostic-fragment expert table, and
#' merges both into a sub-class hint and a continuous evidence score.
#'
#' Pure function: no project, no I/O side effects (the diagnostic CSV is only
#' *read*, never written).  Designed for unit testing.
#'
#' @param spectrum `data.frame(mz, intensity[, frequency])`,
#'   `list(mz =, intensity =)`, or a 2-column matrix.
#' @param precursor_mz Optional precursor m/z (recorded for context; not used
#'   to gate the ladder).
#' @param polarity `"neg"` (default) or `"pos"`.
#' @param cf2_da CF2 repeat-unit exact mass (fixed constant).
#' @param cf2_tol_da Per-step tolerance in Da.
#' @param min_len Minimum number of members for a chain to count as a ladder.
#' @param max_gap_units Maximum hop size `k` in `"loose"` mode.
#' @param min_rel_int Relative-intensity floor (vs the most intense peak).
#' @param series_n_max Highest n for the CnF(2n+1)- series.
#' @param diag_csv Diagnostic table: path, `data.frame`, or `NULL` to resolve
#'   the packaged `extdata/tables/diagnostic_fragments.csv` (skipped if unavailable).
#' @param rules Sub-class rule `data.frame`, or `NULL` for
#'   [pfas_diagnostic_rules].
#' @param mode `"loose"` (default, allows missing rungs) or `"strict"`.
#' @param mz_tol_ppm Symmetric ppm tolerance for series / diagnostic matches.
#' @param int_threshold Drop peaks with `intensity <= int_threshold`.
#' @param weights Named numeric `c(ladder, series, subclass)` (default
#'   `c(0.5, 0.3, 0.2)`; provisional, see design §7).
#' @param ... Reserved for forward compatibility (ignored).
#'
#' @return A `list` with elements `ladder` (`len`, `span`, `mz`, `mz_chr`,
#'   `gaps`), `series` (`max_n`, `hits`, `hits_chr`, `n_hits`, `gap`),
#'   `diagnostic` (`hits`, `subclass_hint`), `score`, `chain`
#'   (`from_series`, `from_ladder`) and `is_hit`.
#'
#' `chain$from_series` and `chain$from_ladder` are two *independent* chain-length
#' estimates: the perfluoroalkyl series count actually observed
#' (`max(1, series$max_n)`) and the ladder span extrapolated in CF2 units
#' (`ladder$span + 1`). They are deliberate cross-checks, **not** a
#' lower/upper interval, so either may be the larger of the two (design §4.5).
#' @seealso [detect_pfas_fragment_series], [pfas_diagnostic_rules]
#' @examples
#' sp <- data.frame(mz = c(68.9957596, 118.992566, 168.9893724),
#'                  intensity = c(100, 80, 60))
#' find_pfas_fragment_series(sp)$ladder$len   # 3
#' @export
find_pfas_fragment_series <- function(
    spectrum,
    precursor_mz = NULL,
    polarity = c("neg", "pos"),
    cf2_da = .nts_pfas_cf2_da,
    cf2_tol_da = 0.01,
    min_len = 3L,
    max_gap_units = 3L,
    min_rel_int = 0.01,
    series_n_max = 10L,
    diag_csv = NULL,
    rules = NULL,
    mode = c("loose", "strict"),
    mz_tol_ppm = 10,
    int_threshold = 0,
    weights = c(ladder = 0.5, series = 0.3, subclass = 0.2),
    ...) {
  mode <- match.arg(mode)
  pol <- .nts_pfas_norm_polarity(if (missing(polarity)) "neg" else polarity)

  cf2_da <- suppressWarnings(as.numeric(cf2_da)[1L])
  if (!is.finite(cf2_da) || cf2_da <= 0) cf2_da <- .nts_pfas_cf2_da
  cf2_tol_da <- suppressWarnings(as.numeric(cf2_tol_da)[1L])
  if (!is.finite(cf2_tol_da) || cf2_tol_da < 0) cf2_tol_da <- 0.01
  min_len <- suppressWarnings(as.integer(min_len)[1L])
  if (!is.finite(min_len) || min_len < 2L) min_len <- 3L
  max_gap_units <- suppressWarnings(as.integer(max_gap_units)[1L])
  if (!is.finite(max_gap_units) || max_gap_units < 1L) max_gap_units <- 3L
  series_n_max <- suppressWarnings(as.integer(series_n_max)[1L])
  if (!is.finite(series_n_max) || series_n_max < 1L) series_n_max <- 10L
  min_rel_int <- suppressWarnings(as.numeric(min_rel_int)[1L])
  if (!is.finite(min_rel_int) || min_rel_int < 0) min_rel_int <- 0.01
  mz_tol_ppm <- suppressWarnings(as.numeric(mz_tol_ppm)[1L])
  if (!is.finite(mz_tol_ppm) || mz_tol_ppm < 0) mz_tol_ppm <- 10
  int_threshold <- suppressWarnings(as.numeric(int_threshold)[1L])
  if (!is.finite(int_threshold)) int_threshold <- 0
  tw <- .nts_pfas_weights(weights)

  # ---- pre-processing (design §4.1) ----
  sp <- .nts_pfas_spectrum_df(spectrum)
  mz <- sp$mz
  inten <- sp$intensity
  keep <- is.finite(mz) & mz > 0 & is.finite(inten) & inten > int_threshold
  mz <- mz[keep]
  inten <- inten[keep]
  if (length(mz) && max(inten) > 0) {
    keep2 <- (inten / max(inten)) >= min_rel_int
    mz <- mz[keep2]
    inten <- inten[keep2]
  }
  if (length(mz)) {
    o <- order(mz)
    mz <- mz[o]
    inten <- inten[o]
    mz <- mz[!duplicated(mz)]
  }
  mz <- mz[is.finite(mz)]

  # ---- core ladder (design §4.2) ----
  ladder <- .nts_pfas_ladder_chain(mz, cf2_da, cf2_tol_da, max_gap_units, mode)
  ladder_present <- ladder$len >= min_len

  # ---- CnF(2n+1)- series (design §4.3) ----
  ns <- seq_len(series_n_max)
  series_mz <- .nts_pfas_cnf2n1_mz(ns, cf2_da)
  ppm <- mz_tol_ppm * 1e-6
  series_hit <- if (identical(pol, "neg") && length(mz)) {
    vapply(series_mz, function(t) any(abs(mz - t) <= t * ppm), logical(1L))
  } else {
    rep(FALSE, length(ns))
  }
  hit_n <- as.integer(ns[series_hit])
  series_max_n <- if (length(hit_n)) max(hit_n) else 0L
  series_gap <- if (length(hit_n) > 1L) {
    max(diff(sort(hit_n))) - 1L
  } else {
    0L
  }

  # ---- diagnostic table + sub-class rules (design §4.4) ----
  matched_lab <- as.character(hit_n[hit_n > 0L])
  matched_lab <- .nts_pfas_cnf2n1_label(matched_lab)
  diag_hits <- .nts_pfas_diag_hits(
    mz, pol, diag_csv = diag_csv, mz_tol_ppm = mz_tol_ppm
  )
  rule_tbl <- .nts_pfas_rules_table(rules)
  if (!is.null(rule_tbl) && length(rule_tbl) && "polarity" %in% names(rule_tbl)) {
    rp <- tolower(trimws(as.character(rule_tbl$polarity)))
    rp[is.na(rp)] <- "neg"
    rule_tbl <- rule_tbl[rp == pol, , drop = FALSE]
  }
  subclass <- .nts_pfas_subclass_hints(
    mz, ladder_present, rule_tbl, mz_tol_ppm = mz_tol_ppm
  )

  # ---- score + chain estimates (design §4.5) ----
  norm_ladder <- if (ladder_present) min(1, ladder$len / 6) else 0
  norm_series <- min(1, series_max_n / 8)
  norm_sub <- min(1, length(subclass) / 2)
  score <- tw[["ladder"]] * norm_ladder +
    tw[["series"]] * norm_series +
    tw[["subclass"]] * norm_sub
  score <- round(as.numeric(score), 4)

  # Two *independent* chain-length estimates (design §4.5), kept as cross-checks:
  # they are NOT a lower/upper interval, so either may come out larger. Hence the
  # names record where each number comes from rather than "lo"/"hi".
  chain_from_series <- max(1L, series_max_n)         # perfluoroalkyl series count
  chain_from_ladder <- as.integer(ladder$span) + 1L  # CF2 ladder span, extrapolated

  # Joint-evidence gate (design §8.7): a bare ~CF2 ladder alone is not enough to
  # raise MS2-level evidence, because 49.9968 Da is not PFAS-exclusive. Require
  # the ladder *and* a corroborating PFAS signal (the CnF(2n+1)- series or a
  # sub-class hint). This keeps `[MS2 Verified]` tagging from being flooded by
  # coincidental 50 Da spacings.
  is_hit <- isTRUE(ladder_present) &&
    (series_max_n >= 1L || length(subclass) > 0L)

  list(
    ladder = ladder,
    series = list(
      max_n = series_max_n,
      hits = hit_n,
      hits_chr = if (length(matched_lab)) paste(matched_lab, collapse = "; ") else NA_character_,
      n_hits = length(hit_n),
      gap = as.integer(series_gap)
    ),
    diagnostic = list(
      hits = diag_hits,
      subclass_hint = if (length(subclass)) paste(subclass, collapse = "; ") else NA_character_,
      subclass = subclass
    ),
    score = score,
    chain = c(from_series = as.integer(chain_from_series),
              from_ladder = as.integer(chain_from_ladder)),
    is_hit = is_hit,
    precursor_mz = if (is.null(precursor_mz)) NA_real_ else as.numeric(precursor_mz)[1L],
    polarity = pol
  )
}


#' @keywords internal
#' @noRd
.nts_pfas_weights <- function(weights = NULL) {
  w <- c(ladder = 0.5, series = 0.3, subclass = 0.2)
  if (!is.null(weights) && length(weights)) {
    nm <- names(weights)
    v <- suppressWarnings(as.numeric(weights))
    if (length(nm) == length(v) && any(nzchar(nm))) {
      for (i in seq_along(v)) {
        if (nm[i] %in% names(w) && is.finite(v[i])) w[[nm[i]]] <- v[i]
      }
    } else if (length(v) == 3L && all(is.finite(v))) {
      w <- c(ladder = v[1L], series = v[2L], subclass = v[3L])
    }
  }
  tot <- sum(w)
  if (is.finite(tot) && tot > 0) w <- w / tot
  w
}


# ---------------------------------------------------------------------------
# Project-level entry point (writes @ms2_annotations, design §5.2)
# ---------------------------------------------------------------------------

#' New `@ms2_annotations` columns produced by [detect_pfas_fragment_series]
#'
#' @return Character vector of column names.
#' @keywords internal
#' @noRd
.nts_pfas_series_columns <- function() {
  c(
    "pfas_ladder_len", "pfas_ladder_mz", "pfas_ladder_span",
    "pfas_series_max_n", "pfas_series_hits",
    "pfas_subclass_hint", "pfas_series_score",
    "pfas_chain_from_series", "pfas_chain_from_ladder"
  )
}


#' Detect PFAS fragment ladders / series across a project
#'
#' Runs [find_pfas_fragment_series] over every feature's consensus MS2 spectrum
#' (`@clean_ms2`) and writes the design §4.6 columns into `@ms2_annotations`
#' (**adds columns only; never alters existing columns**).  Features with a
#' qualifying ladder or a CnF(2n+1)- hit have `@ms1_annotations$annotation_level`
#' tagged `[MS2 Verified]` on exactly the same code path as the 0.1.21
#' diagnostic-fragment validation.
#'
#' Network integration is by **node attributes only** (decision #3): the new
#' columns flow into the network nodes through the existing merge and no edge
#' type is created.
#'
#' @param project An [NTSProject] with `@clean_ms2` (run [purify_ms2_spectra]
#'   first).
#' @param cf2_tol_da,min_len,max_gap_units,min_rel_int,series_n_max,mode,mz_tol_ppm
#'   Forwarded to [find_pfas_fragment_series].
#' @param polarity `NULL` to infer from the project (`neg`/`pos`), or an
#'   explicit `"neg"` / `"pos"`.
#' @param diag_csv,rules Forwarded to [find_pfas_fragment_series]; `NULL`
#'   resolves the packaged defaults.
#' @param tag_annotation_method Logical (default `TRUE`). When a hit feature
#'   already carries a **structural** `@db_annotations$annotation_method`
#'   (`Dual_Network` / `Propagation` / `Network_*` / `Homologue` /
#'   `Transformation`), append the `PFAS_Series` token so the value stays in the
#'   0.1.22 structural guardrail (decision #2). This is provably
#'   confidence-tier neutral (the row is structural either way); non-structural
#'   and blank rows are never touched.
#' @param quiet Logical; suppress the progress message.
#'
#' @return The updated `project` (invisibly).
#' @seealso [find_pfas_fragment_series], [validate_with_ms2]
#' @examples
#' \dontrun{
#' p <- purify_ms2_spectra(p)
#' p <- detect_pfas_fragment_series(p)
#' get_ms2_annotations(p)$pfas_ladder_len
#' }
#' @export
detect_pfas_fragment_series <- function(project,
                                        cf2_tol_da = 0.01,
                                        min_len = 3L,
                                        max_gap_units = 3L,
                                        min_rel_int = 0.01,
                                        series_n_max = 10L,
                                        mode = c("loose", "strict"),
                                        mz_tol_ppm = 10,
                                        polarity = NULL,
                                        diag_csv = NULL,
                                        rules = NULL,
                                        tag_annotation_method = TRUE,
                                        quiet = FALSE) {
  if (!methods::is(project, "NTSProject")) {
    return(project)
  }
  mode <- match.arg(mode)
  tag_annotation_method <- isTRUE(tag_annotation_method)
  quiet <- isTRUE(quiet)

  clean <- methods::slot(project, "clean_ms2")
  if (!is.list(clean) || !length(clean)) {
    if (!quiet) {
      message("> detect_pfas_fragment_series: @clean_ms2 is empty - skipped.")
    }
    return(invisible(project))
  }

  pol <- if (is.null(polarity)) .nts_project_polarity(project) else
    .nts_pfas_norm_polarity(polarity)

  # Read the diagnostic table and rules once, reuse for every spectrum.
  diag_tbl <- .nts_pfas_diag_table(diag_csv)
  rule_tbl <- .nts_pfas_rules_table(rules)

  # Make sure every clean_ms2 feature has a row in @ms2_annotations.
  project <- .nts_ensure_ms2_ann_for_ms2_list(project, clean)

  feat_names <- names(clean)
  if (is.null(feat_names)) {
    feat_names <- vapply(clean, function(z) as.character(z$feature_id)[1L], character(1L))
  }

  res <- lapply(feat_names, function(fn) {
    el <- clean[[fn]]
    spec <- if (is.list(el)) el$spectrum else el
    prec <- if (is.list(el) && !is.null(el$precursorMZ)) el$precursorMZ else NA_real_
    find_pfas_fragment_series(
      spec,
      precursor_mz = prec,
      polarity = pol,
      cf2_tol_da = cf2_tol_da,
      min_len = min_len,
      max_gap_units = max_gap_units,
      min_rel_int = min_rel_int,
      series_n_max = series_n_max,
      diag_csv = diag_tbl,
      rules = rule_tbl,
      mode = mode,
      mz_tol_ppm = mz_tol_ppm
    )
  })
  names(res) <- feat_names

  # ---- fold results into @ms2_annotations (add columns only) ----
  m2 <- as.data.frame(get_ms2_annotations(project), stringsAsFactors = FALSE)
  if (!nrow(m2) || !"name" %in% names(m2)) {
    if (!quiet) {
      message("> detect_pfas_fragment_series: @ms2_annotations has no rows - skipped.")
    }
    return(invisible(project))
  }
  m2$name <- as.character(m2$name)

  new_cols <- .nts_pfas_series_columns()
  for (cn in new_cols) {
    if (!cn %in% names(m2)) m2[[cn]] <- NA
  }

  idx <- match(m2$name, feat_names)
  pick <- function(field, mode) {
    vals <- lapply(res, field)
    out <- vector(mode, length(idx))
    na_val <- switch(mode,
      integer = NA_integer_,
      numeric = NA_real_,
      character = NA_character_,
      NA
    )
    for (k in seq_along(idx)) {
      i <- idx[k]
      out[k] <- if (is.na(i)) na_val else vals[[i]]
    }
    out
  }

  m2$pfas_ladder_len  <- pick(function(r) as.integer(r$ladder$len), "integer")
  m2$pfas_ladder_mz   <- pick(function(r) as.character(r$ladder$mz_chr), "character")
  m2$pfas_ladder_span <- pick(function(r) as.integer(r$ladder$span), "integer")
  m2$pfas_series_max_n <- pick(function(r) as.integer(r$series$max_n), "integer")
  m2$pfas_series_hits <- pick(function(r) as.character(r$series$hits_chr), "character")
  m2$pfas_subclass_hint <- pick(function(r) as.character(r$diagnostic$subclass_hint), "character")
  m2$pfas_series_score <- pick(function(r) as.numeric(r$score), "numeric")
  m2$pfas_chain_from_series <- pick(function(r) as.integer(r$chain[["from_series"]]), "integer")
  m2$pfas_chain_from_ladder <- pick(function(r) as.integer(r$chain[["from_ladder"]]), "integer")

  project <- set_ms2_annotations(project, m2)

  # ---- annotation_level tagging (same helper as 0.1.21) ----
  hit_names <- names(res)[vapply(res, function(r) isTRUE(r$is_hit), logical(1L))]
  if (length(hit_names)) {
    ms1 <- get_ms1_annotations(project)
    if (is.data.frame(ms1) && nrow(ms1) && "name" %in% names(ms1)) {
      if (!"annotation_level" %in% names(ms1)) {
        ms1$annotation_level <- NA_character_
      }
      w <- as.character(ms1$name) %in% hit_names
      if (any(w)) {
        ms1$annotation_level[w] <- .nts_annotation_level_after_ms2(
          as.character(ms1$annotation_level[w])
        )
        project <- set_ms1_annotations(project, ms1)
      }
    }
  }

  # ---- structural guardrail stamp (decision #2; tier neutral) ----
  if (tag_annotation_method && length(hit_names)) {
    project <- .nts_pfas_series_tag_annotation_method(project, hit_names)
  }

  n_lad <- sum(vapply(res, function(r) isTRUE(r$ladder$len >= as.integer(min_len)), logical(1L)))
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "detect_pfas_fragment_series",
    list(
      mode = mode,
      cf2_tol_da = cf2_tol_da,
      min_len = as.integer(min_len),
      max_gap_units = as.integer(max_gap_units),
      min_rel_int = min_rel_int,
      series_n_max = as.integer(series_n_max),
      polarity = pol,
      n_features = length(feat_names),
      n_ladder = n_lad,
      n_hit = length(hit_names)
    )
  )
  if (!quiet) {
    message(
      "> detect_pfas_fragment_series: ", length(feat_names), " spectra -> ",
      n_lad, " ladder(s), ", length(hit_names),
      " feature(s) with CF2-series evidence."
    )
  }
  invisible(project)
}


#' Append the PFAS_Series token to structural annotation_method values
#'
#' Only rows whose `annotation_method` already contains a structural token are
#' touched, so `am_structural` is TRUE before and after and the confidence tier
#' is unchanged (decision #2).  Blank / spectral-only rows are left alone.
#' @keywords internal
#' @noRd
.nts_pfas_series_tag_annotation_method <- function(project, hit_names) {
  db <- methods::slot(project, "db_annotations")
  if (!is.data.frame(db) || !nrow(db)) {
    return(project)
  }
  fcol <- if ("name" %in% names(db)) "name" else if ("feature_id" %in% names(db)) "feature_id" else NULL
  if (is.null(fcol)) {
    return(project)
  }
  if (!"annotation_method" %in% names(db)) {
    db$annotation_method <- ""
  }
  am <- trimws(as.character(db$annotation_method))
  am[is.na(am)] <- ""
  structural <- grepl("Dual_Network", am, fixed = TRUE) |
    grepl("Propagation", am, fixed = TRUE) |
    grepl("Network_", am, fixed = TRUE) |
    grepl("Homologue", am, fixed = TRUE) |
    grepl("Transformation", am, fixed = TRUE)
  already <- grepl("PFAS_Series", am, fixed = TRUE)
  feat <- trimws(as.character(db[[fcol]]))
  sel <- !is.na(feat) & feat %in% hit_names & structural & !already
  if (any(sel)) {
    db$annotation_method[sel] <- paste0(am[sel], "; PFAS_Series")
    methods::slot(project, "db_annotations") <- db
  }
  project
}


#' Infer a project's dominant polarity as "neg" / "pos"
#' @keywords internal
#' @noRd
.nts_project_polarity <- function(project) {
  ft <- tryCatch(get_feature_table(project), error = function(e) NULL)
  if (is.data.frame(ft) && nrow(ft)) {
    pcol <- intersect(c("polarity", "Polarity", "ion_mode", "mode"), names(ft))
    if (length(pcol)) {
      v <- tolower(trimws(as.character(ft[[pcol[1L]]])))
      v <- v[!is.na(v) & nzchar(v)]
      if (length(v)) {
        v[v %in% c("negative", "neg", "-")] <- "neg"
        v[v %in% c("positive", "pos", "+")] <- "pos"
        tab <- sort(table(v), decreasing = TRUE)
        if (length(tab) && names(tab)[1L] %in% c("neg", "pos")) {
          return(names(tab)[1L])
        }
      }
    }
  }
  "neg"
}
