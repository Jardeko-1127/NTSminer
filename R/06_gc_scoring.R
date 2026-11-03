# GC-EI spectral library matching - cosine similarity via OrgMassSpecR ---------

#' Score GC-EI experimental spectra against a reference library
#'
#' Compares each experimental GC-EI spectrum in an NTSproject against
#' a reference spectral library using cosine similarity via
#' [OrgMassSpecR::SpectrumSimilarity]. When \code{prefilter_threshold > 0},
#' an inverted m/z index is built for fast candidate lookup (no full-table
#' copy - only m/z + lab_id columns are extracted).
#'
#' @param project An [NTSProject] with populated \code{@gc_spectra}.
#' @param library An [NTSDatabase] containing reference spectra. The
#'   \code{spectra_long} slot must have columns \code{lab_id}, \code{mz},
#'   and \code{intensity}.
#' @param mz_tol Numeric, m/z tolerance in Da for peak matching between
#'   experimental and reference spectra. Passed as \code{t} to
#'   [OrgMassSpecR::SpectrumSimilarity]. Default 0.25.
#' @param min_intensity Numeric, minimum relative intensity (0-100 scale)
#'   for peaks to be included in matching. Passed as \code{b} to
#'   [OrgMassSpecR::SpectrumSimilarity]. Default 10.
#' @param filter_query Numeric, pre-filter query spectra before matching:
#'   peaks with relative intensity below this threshold (0-100 scale) are
#'   removed. Applied BEFORE candidate selection and cosine calculation.
#'   Default 0 (no filtering). Recommended 0.15 for nominal-mass EI data
#'   (removes sub-0.15% noise). Distinct from \code{min_intensity}, which
#'   is passed to [OrgMassSpecR::SpectrumSimilarity] for per-pair filtering.
#' @param min_score Numeric, minimum cosine similarity score (0-1) to
#'   retain a match. Default 0. Lower this to keep all results for
#'   downstream filtering.
#' @param use_deconvolved Logical, if TRUE and \code{@gc_deconvolved}
#'   is non-empty, score each deconvolved component separately instead
#'   of the raw spectra in \code{@gc_spectra}. Default FALSE.
#' @param method Character, matching method. \code{"cosine"} (default) uses
#'   [OrgMassSpecR::SpectrumSimilarity] for cosine similarity.
#'   \code{"ei_normal"} uses strict identity matching: each query peak must
#'   have a matching peak in the reference within \code{mz_tol}, and the
#'   score is the fraction of matched query intensity.
#' @param prefilter_threshold Integer, minimum number of common m/z
#'   bins required before computing cosine. Default 3 for cosine method,
#'   1 for ei_normal. Use 0 to disable pre-filtering entirely.
#'   Inverted index is built whenever this is > 0 (regardless of library size).
#' @param rt_tol Numeric, retention time tolerance in minutes. If provided
#'   together with valid RT data, only library entries whose RT is within
#'   \code{rt_tol} of the query RT are considered. Default \code{NULL}
#'   (no RT filtering).
#' @param rt_lookup Named numeric vector mapping feature names to RT
#'   values (in minutes). Required when \code{rt_tol} is set.
#' @param lib_rt_col Character, column name in \code{library@metadata}
#'   containing RT values (in minutes). Required when \code{rt_tol} is set.
#' @param cas_bridge A [data.table::data.table] with columns \code{feature}
#'   and \code{lab_id} (or \code{CAS}), specifying pre-determined
#'   query-to-library mappings that bypass the m/z index entirely.
#'   Useful for CAS-based cross-library matching. Default \code{NULL}.
#' @param max_candidates Integer, maximum number of candidate library
#'   compounds to evaluate per feature. Default \code{Inf} (no limit).
#'   Reducing this speeds up matching for large libraries but may
#'   miss the best match. Recommended: 100-200 for NIST-scale libraries.
#' @param mass_tol Numeric, molecular mass tolerance in Da for pre-filtering
#'   library candidates. When non-\code{NULL}, only library compounds whose
#'   \code{mass} (from \code{library@metadata}) is within \code{mass_tol}
#'   of the query spectrum's highest m/z are considered. This dramatically
#'   reduces candidates for large libraries (e.g. NIST23 from 346K to
#'   dozens). Set to 5 for a loose filter, 0.5 for tight. Default
#'   \code{NULL} (no mass filtering).
#' @param best_only Logical, if \code{TRUE} return only the best match
#'   per feature (by cosine). Default \code{FALSE}.
#'
#' @return A [data.table::data.table] with columns:
#'   \item{feature_id}{Feature name from \code{@gc_spectra}}
#'   \item{ref_id}{Reference compound ID from the library}
#'   \item{ref_name}{Reference compound name}
#'   \item{cosine}{Similarity score (0-1)}
#'   \item{method}{Matching method used}
#'
#' @seealso [OrgMassSpecR::SpectrumSimilarity], [annotate_gc_by_spectral_match]
#' @importFrom data.table data.table rbindlist setorder
#' @export
score_ei_spectra <- function(project,
                             library,
                             mz_tol = 0.25,
                             min_intensity = 10,
                             min_score = 0,
                             filter_query = 0,
                             use_deconvolved = FALSE,
                             method = c("cosine", "ei_normal"),
                             prefilter_threshold = NULL,
                             rt_tol = NULL,
                             rt_lookup = NULL,
                             lib_rt_col = NULL,
                             cas_bridge = NULL,
                             max_candidates = Inf,
                             mass_tol = NULL,
                             best_only = FALSE) {
  .nts_check(project, "project", .arg = "project")
  .nts_check(library, "reference_db", .arg = "library")
  if (!requireNamespace("OrgMassSpecR", quietly = TRUE)) stop("OrgMassSpecR is required", call. = FALSE)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required", call. = FALSE)

  method <- match.arg(method)

  # Default prefilter: 3 for cosine, 1 for ei_normal
  if (is.null(prefilter_threshold)) {
    prefilter_threshold <- if (method == "cosine") 3L else 1L
  }
  prefilter_threshold <- as.integer(prefilter_threshold[1L])

  gc <- project@gc_spectra
  if (!length(gc)) stop("@gc_spectra is empty; run extract_gc_spectra() first", call. = FALSE)

  if (use_deconvolved && length(project@gc_deconvolved)) {
    d <- project@gc_deconvolved
    gc <- list()
    feature_names <- character(0L)
    for (fnm in names(d)) {
      comps <- d[[fnm]]
      if (!length(comps)) {
        if (fnm %in% names(project@gc_spectra)) {
          gc[[fnm]] <- project@gc_spectra[[fnm]]
          feature_names <- c(feature_names, fnm)
        }
        next
      }
      for (ci in seq_along(comps)) {
        cid <- if (length(comps) > 1L) paste0(fnm, "_c", ci) else fnm
        gc[[cid]] <- comps[[ci]]$spectrum
        feature_names <- c(feature_names, cid)
      }
    }
    names(gc) <- feature_names
  }

  lib_spec <- library@spectra_long
  lib_meta <- library@metadata
  if (!nrow(lib_spec)) stop("library@spectra_long is empty", call. = FALSE)
  if (!all(c("lab_id", "mz", "intensity") %in% names(lib_spec))) {
    stop("library@spectra_long must have columns: lab_id, mz, intensity", call. = FALSE)
  }
  if (!"lab_id" %in% names(lib_meta)) {
    stop("library@metadata must have column: lab_id", call. = FALSE)
  }

  # Build reference spectrum index: split spectra_long by lab_id
  ref_list <- split(lib_spec, by = "lab_id", keep.by = FALSE)
  ref_name_map <- lib_meta$name
  names(ref_name_map) <- as.character(lib_meta$lab_id)

  # ---- RT lookup setup ----
  use_rt <- !is.null(rt_tol) && !is.null(rt_lookup) && !is.null(lib_rt_col)
  if (use_rt && !lib_rt_col %in% names(lib_meta)) {
    warning("lib_rt_col '", lib_rt_col, "' not found in library metadata; RT filter disabled", call. = FALSE)
    use_rt <- FALSE
  }
  lib_rt_map <- if (use_rt) setNames(lib_meta[[lib_rt_col]], lib_meta$lab_id) else NULL

  # ---- CAS-Bridge mode ----
  use_cas_bridge <- !is.null(cas_bridge) && nrow(cas_bridge) > 0L
  cas_map <- NULL
  if (use_cas_bridge) {
    if (!"feature" %in% names(cas_bridge))
      stop("cas_bridge must have a 'feature' column", call. = FALSE)
    bridge_col <- if ("lab_id" %in% names(cas_bridge)) "lab_id" else names(cas_bridge)[2L]
    cas_map <- setNames(cas_bridge[[bridge_col]], cas_bridge$feature)
  }

  # ---- Mass pre-filter setup ----
  use_mass <- !is.null(mass_tol) && "mass" %in% names(lib_meta)
  if (!is.null(mass_tol) && !use_mass) {
    warning("library@metadata has no 'mass' column; mass filter disabled", call. = FALSE)
  }
  lib_mass_map <- if (use_mass) {
    masses <- lib_meta$mass
    names(masses) <- as.character(lib_meta$lab_id)
    masses
  } else NULL

  # ---- Inverted m/z index ----
  n_lib <- length(ref_list)
  use_index <- prefilter_threshold > 0L

  if (use_index) {
    bin_factor <- max(1, floor(1 / max(mz_tol, 0.01)))
    # Lightweight: extract only the 2 columns we need (no full-table copy)
    idx_dt <- lib_spec[, .(lab_id, mz_bin = as.character(round(mz * bin_factor)))]
    mz_idx <- idx_dt[, .(lab_ids = list(unique(lab_id))), by = mz_bin]
    mz_bins_map <- setNames(mz_idx$lab_ids, mz_idx$mz_bin)
    message(sprintf("Built inverted index: %d unique m/z bins for %d library compounds",
                    nrow(mz_idx), n_lib))
  }

  feature_names <- names(gc)
  if (is.null(feature_names)) {
    feature_names <- as.character(seq_along(gc))
    names(gc) <- feature_names
  }

  results <- vector("list", length(feature_names))

  for (i in seq_along(feature_names)) {
    fnm <- feature_names[i]
    qspec <- gc[[fnm]]

    if (!is.data.frame(qspec) || !nrow(qspec) || !all(c("mz", "intensity") %in% names(qspec))) {
      next
    }

    # Pre-filter query spectrum (remove noise peaks before matching)
    if (filter_query > 0) {
      qspec <- qspec[qspec[["intensity"]] >= filter_query, , drop = FALSE]
      if (nrow(qspec) < 3) next
    }

    q_rt <- if (use_rt && fnm %in% names(rt_lookup)) rt_lookup[fnm] else NA_real_

    # ---- Candidate selection ----
    if (use_cas_bridge && fnm %in% names(cas_map)) {
      cand_ids <- as.character(cas_map[fnm])
      cand_ids <- cand_ids[!is.na(cand_ids)]
      cand_ids <- intersect(cand_ids, names(ref_list))
      if (!length(cand_ids)) next
    } else if (use_index) {
      sp_bins <- as.character(round(qspec$mz * bin_factor))
      cand_list <- unlist(lapply(sp_bins, function(b) mz_bins_map[[b]]),
                          use.names = FALSE)
      if (!length(cand_list)) next
      cand_tab <- sort(table(cand_list), decreasing = TRUE)
      cand_tab <- cand_tab[cand_tab >= prefilter_threshold]
      if (!length(cand_tab)) next
      cand_ids <- names(cand_tab)
      if (is.finite(max_candidates) && length(cand_ids) > max_candidates) {
        cand_ids <- cand_ids[seq_len(max_candidates)]
      }
    } else {
      cand_ids <- names(ref_list)
    }

    # RT pre-filter (if enabled, applied before cosine)
    if (use_rt && !is.na(q_rt) && length(cand_ids) > 1L) {
      keep <- vapply(cand_ids, function(lid) {
        l_rt <- lib_rt_map[lid]
        if (is.na(l_rt)) return(TRUE)  # no RT = keep
        abs(q_rt - l_rt) <= rt_tol
      }, logical(1L))
      cand_ids <- cand_ids[keep]
      if (!length(cand_ids)) next
    }

    # Mass pre-filter (if enabled, applied before cosine)
    if (use_mass && length(cand_ids) > 1L) {
      q_mass <- max(qspec$mz)
      if (!is.na(q_mass)) {
        keep <- vapply(cand_ids, function(lid) {
          l_mass <- lib_mass_map[lid]
          if (is.na(l_mass)) return(TRUE)
          abs(q_mass - l_mass) <= mass_tol
        }, logical(1L))
        cand_ids <- cand_ids[keep]
        if (!length(cand_ids)) next
      }
    }

    feature_results <- vector("list", length(cand_ids))
    ridx <- 0L

    for (lid in cand_ids) {
      ref_spec <- ref_list[[lid]]
      if (!is.data.frame(ref_spec) || !nrow(ref_spec)) next

      if (method == "cosine") {
        cos_val <- tryCatch(
          OrgMassSpecR::SpectrumSimilarity(
            spec.top      = ref_spec,
            spec.bottom   = qspec,
            t             = mz_tol,
            b             = min_intensity,
            print.graphic = FALSE,
            output.list   = FALSE
          ),
          error = function(e) NA_real_
        )
      } else {
        q_int_filt <- qspec$intensity[qspec$intensity >= min_intensity]
        q_mz_filt  <- qspec$mz[qspec$intensity >= min_intensity]
        if (!length(q_mz_filt)) {
          cos_val <- NA_real_
        } else {
          total_int <- sum(q_int_filt)
          matched_int <- 0
          for (pi in seq_along(q_mz_filt)) {
            if (any(abs(ref_spec$mz - q_mz_filt[pi]) <= mz_tol)) {
              matched_int <- matched_int + q_int_filt[pi]
            }
          }
          cos_val <- if (total_int > 0) matched_int / total_int else 0
        }
      }

      if (is.na(cos_val)) next

      ridx <- ridx + 1L
      feature_results[[ridx]] <- data.table::data.table(
        feature_id = fnm,
        ref_id     = lid,
        ref_name   = if (lid %in% names(ref_name_map)) ref_name_map[lid] else NA_character_,
        cosine     = cos_val,
        method     = method
      )
    }

    if (ridx > 0L) {
      results[[i]] <- data.table::rbindlist(feature_results[1L:ridx])
    }
  }

  out <- data.table::rbindlist(results)
  if (!nrow(out)) {
    return(data.table::data.table(
      feature_id = character(0L), ref_id = character(0L),
      ref_name = character(0L), cosine = numeric(0L),
      method = character(0L)
    ))
  }

  data.table::setorder(out, -cosine, na.last = TRUE)
  if (best_only) {
    out <- out[, .SD[1L], by = feature_id]
  }
  out[cosine >= min_score]
}


#' Annotate GC-EI features by spectral library matching
#'
#' Calls [score_ei_spectra] and writes the top-N matches per feature
#' to \code{project@db_annotations}.
#'
#' @param project An [NTSProject] with populated \code{@gc_spectra}.
#' @param library An [NTSDatabase] containing reference spectra.
#' @param min_score Numeric, minimum cosine similarity (0-1) for a
#'   match to be retained. Default 0.7.
#' @param top_n Integer, maximum number of top matches to retain per
#'   feature. Default 3.
#' @param ... Additional arguments passed to [score_ei_spectra].
#'
#' @return The project with \code{@db_annotations} populated.
#'
#' @seealso [score_ei_spectra], [NTSDatabase], [batch_annotate]
#' @export
annotate_gc_by_spectral_match <- function(project,
                                          library,
                                          min_score = 0.7,
                                          top_n = 3L,
                                          ...) {
  .nts_check(project, "project", .arg = "project")
  .nts_check(library, "reference_db", .arg = "library")

  scores <- score_ei_spectra(project, library, min_score = min_score, ...)

  if (!nrow(scores)) {
    warning("No GC-EI spectral matches found with min_score >= ", min_score, call. = FALSE)
    return(project)
  }

  # Take top_n per feature
  scores[, rank := seq_len(.N), by = feature_id]
  scores <- scores[rank <= top_n]
  scores[, rank := NULL]

  # Build db_annotations table
  db <- data.table::data.table(
    name            = scores$feature_id,
    ref_id          = scores$ref_id,
    Compound.name   = scores$ref_name,
    score_total     = scores$cosine,
    score_cosine    = scores$cosine,
    match_method    = "EI_cosine"
  )

  # Merge with existing db_annotations if present
  existing <- project@db_annotations
  if (is.data.frame(existing) && nrow(existing)) {
    existing <- existing[!as.character(existing$name) %in% unique(db$name), , drop = FALSE]
    for (cn in setdiff(names(db), names(existing))) {
      existing[[cn]] <- rep(NA, nrow(existing))
    }
    for (cn in setdiff(names(existing), names(db))) {
      db[[cn]] <- rep(NA, nrow(db))
    }
    db <- rbind(existing, db[, names(existing), drop = FALSE])
  }

  project@db_annotations <- as.data.frame(db, stringsAsFactors = FALSE)
  validObject(project)
  project
}
