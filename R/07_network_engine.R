#' @title Annotate transformation edges in a \code{network_data} list (read-only project)
#'
#' @description
#' Enriches \code{network_data$edges} using precursor mass differences (\code{mz1},
#' \code{mz2}) against the same rule schema as [search_transformations]
#' (\code{rule_name}, \code{mass_shift}, \code{sub_formula}, \code{add_formula}),
#' writing the ppm-matched rule name to column \code{rule_name} on edges. When
#' edges already carry \code{transformation_rule} (e.g. from
#' [search_transformations]), merges \code{sub_formula} / \code{add_formula}
#' from the matching rule rows and projects \code{transformed_formula} onto
#' target \code{nodes}. Does **not** modify \code{NTSProject} slots.
#'
#' @param project An [NTSProject] (read for optional \code{ms1_annotations} joins).
#' @param network_data Named \code{list} with \code{nodes} and \code{edges}.
#' @param rules \code{NULL} or a \code{data.table} / \code{data.frame} with columns
#'   \code{rule_name}, \code{mass_shift}, \code{sub_formula}, \code{add_formula}
#'   (same contract as [search_transformations]). \code{NULL} uses a small
#'   built-in default (oxidation / +O).
#' @param rules_csv Optional path or table forwarded to rules loading when
#'   \code{rules} is \code{NULL} (same schema as \code{rules}).
#' @param match_tol_ppm Single ppm tolerance applied to every row of \code{rules}
#'   when comparing \code{abs(mz1 - mz2)} to \code{abs(mass_shift)} (or to
#'   \code{mass_shift} when \code{abs(mass_shift)} is near zero).
#'
#' @return Enhanced \code{network_data}: \code{list(nodes =, edges =)} as
#'   \code{data.table}s with annotation columns; intermediate / join scaffold
#'   columns stripped before return.
#'
#' @export
annotate_transformations <- function(project,
                                     network_data,
                                     rules = NULL,
                                     rules_csv = NULL,
                                     match_tol_ppm = 20) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject object.", call. = FALSE)
  }
  if (!is.list(network_data) || !all(c("nodes", "edges") %in% names(network_data))) {
    stop("network_data must be a list with elements `nodes` and `edges`.", call. = FALSE)
  }
  extra <- setdiff(names(network_data), c("nodes", "edges"))
  nodes_was_df <- is.data.frame(network_data$nodes) &&
    !data.table::is.data.table(network_data$nodes)
  edges_was_df <- is.data.frame(network_data$edges) &&
    !data.table::is.data.table(network_data$edges)
  nd <- data.table::as.data.table(network_data$nodes)
  ed <- data.table::as.data.table(network_data$edges)
  if (!nrow(ed)) {
    ed[, `:=`(
      delta_mz = numeric(),
      rule_name = character()
    )]
    out <- list(nodes = nd, edges = ed)
    for (nm in extra) {
      out[[nm]] <- network_data[[nm]]
    }
    out <- .nts_network_data_strip_intermediate_cols(out, strip_edges = TRUE)
    if (nodes_was_df) {
      out$nodes <- as.data.frame(out$nodes, stringsAsFactors = FALSE)
    }
    if (edges_was_df) {
      out$edges <- as.data.frame(out$edges, stringsAsFactors = FALSE)
    }
    return(out)
  }

  if (!all(c("mz1", "mz2") %in% names(ed))) {
    if (!all(c("from", "to") %in% names(ed))) {
      stop(
        "edges must contain mz1 and mz2, or from and to with node mz for derivation.",
        call. = FALSE
      )
    }
    key <- if ("id" %in% names(nd)) "id" else if ("name" %in% names(nd)) "name" else {
      stop("nodes must contain `id` or `name` to derive mz1/mz2.", call. = FALSE)
    }
    if (!"mz" %in% names(nd)) {
      stop("nodes must contain column `mz` when edges lack mz1/mz2.", call. = FALSE)
    }
    mz_map <- stats::setNames(
      suppressWarnings(as.numeric(nd$mz)),
      as.character(nd[[key]])
    )
    ed[, mz1 := suppressWarnings(as.numeric(mz_map[as.character(from)]))]
    ed[, mz2 := suppressWarnings(as.numeric(mz_map[as.character(to)]))]
  }

  ed[, delta_mz := abs(suppressWarnings(as.numeric(mz1)) - suppressWarnings(as.numeric(mz2)))]

  if (is.null(rules)) {
    if (!is.null(rules_csv)) {
      if (is.character(rules_csv) || is.factor(rules_csv)) {
        rp <- path.expand(as.character(rules_csv)[1L])
        .nts_check(rp, "file_exists", .arg = "rules_csv")
        rules <- data.table::fread(rp, data.table = TRUE, na.strings = c("", "NA"))
      } else {
        rules <- data.table::as.data.table(rules_csv)
      }
    } else {
      rules <- .nts_default_annotate_mass_shift_rules_dt()
    }
  } else {
    rules <- data.table::as.data.table(rules)
  }

  req <- c("rule_name", "mass_shift", "sub_formula", "add_formula")
  miss <- setdiff(req, names(rules))
  if (length(miss)) {
    stop(
      "rules table missing columns: ",
      paste(miss, collapse = ", "),
      call. = FALSE
    )
  }
  rules[, rule_name := as.character(rule_name)]
  rules[, mass_shift := suppressWarnings(as.numeric(mass_shift))]
  rules[, sub_formula := as.character(sub_formula)]
  rules[, add_formula := as.character(add_formula)]
  rules <- rules[is.finite(mass_shift)]

  tol <- suppressWarnings(as.numeric(match_tol_ppm)[1L])
  if (!is.finite(tol) || tol < 0) {
    stop("match_tol_ppm must be a non-negative finite number", call. = FALSE)
  }

  ed[, rule_name := NA_character_]

  if (nrow(rules)) {
    mz1n <- suppressWarnings(as.numeric(ed$mz1))
    mz2n <- suppressWarnings(as.numeric(ed$mz2))
    ref_mz <- pmax(abs(mz1n), abs(mz2n), na.rm = TRUE)
    ref_mz[!is.finite(ref_mz) | ref_mz < 1e-12] <- 1e-6
    del <- ed$delta_mz
    for (j in seq_len(nrow(rules))) {
      sh <- rules$mass_shift[j]
      rn <- rules$rule_name[j]
      if (!is.finite(sh)) {
        next
      }
      sh_abs <- abs(sh)
      if (sh_abs < 1e-12) {
        ppm_ok <- abs(del - sh_abs) <= (tol / 1e6) * ref_mz
      } else {
        ppm_ok <- (abs(del - sh_abs) / sh_abs) * 1e6 <= tol
      }
      w <- which(
        ppm_ok &
          (is.na(ed$rule_name) | !nzchar(trimws(as.character(ed$rule_name))))
      )
      if (length(w)) {
        data.table::set(ed, i = w, j = "rule_name", value = rep(rn, length(w)))
      }
    }
  }

  if ("transformation_rule" %in% names(ed) &&
      any(!is.na(ed$transformation_rule) & nzchar(trimws(as.character(ed$transformation_rule))))) {
    keyn <- if ("id" %in% names(nd)) "id" else "name"
    nd[, feat_key := as.character(get(keyn))]

    for (cn in c("transformed_formula", "transformation_rule", "transformation_parent")) {
      if (!cn %in% names(nd)) {
        nd[, (cn) := NA_character_]
      }
    }

    trt_u <- unique(
      rules[, .(
        transformation_rule = as.character(rule_name),
        sub_formula_r = as.character(sub_formula),
        add_formula_r = as.character(add_formula)
      )],
      by = "transformation_rule"
    )
    trt_u <- trt_u[
      !is.na(transformation_rule) & nzchar(trimws(as.character(transformation_rule)))
    ]
    trt_d <- unique(
      .nts_default_transformation_rules_dt()[, .(
        transformation_rule = as.character(rule_name),
        sub_formula_r = as.character(sub_formula),
        add_formula_r = as.character(add_formula)
      )],
      by = "transformation_rule"
    )
    trt <- if (nrow(trt_u) > 0L) {
      unique(
        data.table::rbindlist(list(trt_u, trt_d), fill = TRUE),
        by = "transformation_rule"
      )
    } else {
      unique(trt_d, by = "transformation_rule")
    }
    ed <- trt[ed, on = "transformation_rule"]

    pmap <- nd[, .(
      feat_key = feat_key,
      seed_f = {
        fd <- rep(NA_character_, .N)
        if ("formula" %in% names(nd)) {
          fd <- as.character(formula)
        }
        fi <- rep(NA_character_, .N)
        if ("inferred_formula" %in% names(nd)) {
          fi <- as.character(inferred_formula)
        }
        data.table::fifelse(!is.na(fd) & nzchar(trimws(fd)), fd, fi)
      }
    )]
    pmap <- unique(pmap[!is.na(seed_f) & nzchar(trimws(seed_f))], by = "feat_key")

    ed[, parent_key := as.character(from)]
    ed <- pmap[ed, on = c(feat_key = "parent_key")]
    ed[, transformed_formula_edge := as.character(mapply(
      function(sf, su, ad) {
        .nts_apply_seed_transform(sf, su, ad)
      },
      seed_f,
      sub_formula_r,
      add_formula_r,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    ))]
    tmp_cols_ed <- c(
      "sub_formula_r", "add_formula_r", "seed_f", "parent_key", "feat_key",
      "Reaction_Name", "Trans_Type", "Type", "Mass_Shift", "Tolerance_ppm"
    )
    ex_ed <- intersect(tmp_cols_ed, names(ed))
    if (length(ex_ed) > 0L) {
      ed[, (ex_ed) := NULL]
    }

    ch <- unique(
      ed[!is.na(transformed_formula_edge) & nzchar(trimws(as.character(transformed_formula_edge))), .(
        feat_key = as.character(to),
        tf_hit = as.character(transformed_formula_edge),
        rule_hit = as.character(transformation_rule),
        par_hit = as.character(from)
      )],
      by = "feat_key"
    )
    nd <- ch[nd, on = "feat_key"]
    nd[!is.na(tf_hit) & nzchar(trimws(as.character(tf_hit))) &
      (is.na(transformed_formula) | !nzchar(trimws(as.character(transformed_formula)))),
    transformed_formula := tf_hit]
    nd[!is.na(rule_hit) & nzchar(trimws(as.character(rule_hit))) &
      (is.na(transformation_rule) | !nzchar(trimws(as.character(transformation_rule)))),
    transformation_rule := rule_hit]
    nd[!is.na(par_hit) & nzchar(trimws(as.character(par_hit))) &
      (is.na(transformation_parent) | !nzchar(trimws(as.character(transformation_parent)))),
    transformation_parent := par_hit]
    tmp_cols_nd <- c("tf_hit", "rule_hit", "par_hit")
    ex_nd <- intersect(tmp_cols_nd, names(nd))
    if (length(ex_nd) > 0L) {
      nd[, (ex_nd) := NULL]
    }
    if ("transformed_formula_edge" %in% names(ed)) {
      ed[, transformed_formula_edge := NULL]
    }
    if ("feat_key" %in% names(nd)) {
      nd[, feat_key := NULL]
    }
  }

  out <- list(nodes = nd, edges = ed)
  for (nm in extra) {
    out[[nm]] <- network_data[[nm]]
  }
  out <- .nts_network_data_strip_intermediate_cols(out, strip_edges = TRUE)
  if (nodes_was_df) {
    out$nodes <- as.data.frame(out$nodes, stringsAsFactors = FALSE)
  }
  if (edges_was_df) {
    out$edges <- as.data.frame(out$edges, stringsAsFactors = FALSE)
  }
  out
}

#' Default rules for [annotate_transformations] (same schema as [search_transformations]).
#' @noRd
.nts_default_annotate_mass_shift_rules_dt <- function() {
  data.table::data.table(
    rule_name = "Oxidation",
    mass_shift = 15.994915,
    sub_formula = NA_character_,
    add_formula = "O"
  )
}

#' Strip homologue / XCMS / join-scaffold columns from \code{network_data} tables.
#' @keywords internal
#' @noRd
.nts_network_data_strip_intermediate_cols <- function(net_data, strip_edges = FALSE) {
  if (!is.list(net_data)) {
    return(net_data)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    return(net_data)
  }
  cols_to_remove <- c(
    "inferred_name", "inference_evidence_path",
    "homologue_inference_n_paths", "homologue_inference_transformations",
    "homologue_inference_isf_pseudo", "inferred_formula", "inferred_smiles",
    "i.inferred_name",
    "exp_mz", "lib_mz", "mz_error_ppm", "exp_rt", "lib_rt", "rt_error_s",
    "n_matched_peaks", "mzmin", "mzmax", "rtmin", "rtmax",
    "library_Compound_name", "library_ref_id", "library_SMILES",
    "Reaction_Name", "Trans_Type", "Type", "Mass_Shift", "Tolerance_ppm"
  )
  # Transformation / generic network nodes: drop redundant Hill columns (edges keep
  # `transformed_formula` when present for lineage export).
  cols_to_remove_nodes <- c(cols_to_remove, "formula", "transformed_formula")
  strip_one <- function(obj, remove_names = cols_to_remove) {
    if (is.null(obj)) {
      return(obj)
    }
    was_df <- is.data.frame(obj) && !data.table::is.data.table(obj)
    dt <- data.table::as.data.table(obj)
    ex <- intersect(remove_names, colnames(dt))
    if (length(ex) > 0L) {
      dt[, (ex) := NULL]
    }
    iex <- grep("^i\\.", colnames(dt), value = TRUE)
    if (length(iex) > 0L) {
      dt[, (iex) := NULL]
    }
    if (was_df) {
      as.data.frame(dt, stringsAsFactors = FALSE)
    } else {
      dt
    }
  }
  net_data$nodes <- strip_one(net_data$nodes, cols_to_remove_nodes)
  if (isTRUE(strip_edges) && !is.null(net_data$edges)) {
    net_data$edges <- strip_one(net_data$edges, cols_to_remove)
  }
  net_data
}

# ============================================================================
# Molecular Networking (MS2 spectral similarity graph)
# ============================================================================

#' @title Pairwise MS2 scoring in a precursor m/z window
#'
#' @description
#' Avoids full O(N^2) by sorting spectra on precursor m/z and only comparing
#' pairs within `precursor_tol` Da.
#' @noRd
.mn_pairwise_scores <- function(spec_list, prec_mz, precursor_tol, min_score) {
  n <- length(spec_list)
  if (n < 2L) {
    return(data.table::data.table(i = integer(0), j = integer(0), score = numeric(0)))
  }
  ord <- order(prec_mz)
  prec_sorted <- prec_mz[ord]

  from_vec <- integer(0)
  to_vec   <- integer(0)
  sc_vec   <- numeric(0)

  for (a in seq_len(n - 1L)) {
    b <- a + 1L
    while (b <= n && (prec_sorted[b] - prec_sorted[a]) <= precursor_tol) {
      s <- tryCatch(
        .calculate_clean_ndp(spec_list[[ord[a]]], spec_list[[ord[b]]],
                             mz_tol = 0.01, int_thresh = 0.025)$score,
        error = function(e) NA_real_
      )
      if (!is.na(s) && s >= min_score) {
        from_vec <- c(from_vec, ord[a])
        to_vec   <- c(to_vec, ord[b])
        sc_vec   <- c(sc_vec, s)
      }
      b <- b + 1L
    }
  }
  data.table::data.table(i = from_vec, j = to_vec, score = sc_vec)
}

#' Build a molecular network from MS2 spectral similarity
#'
#' Computes modified cosine similarity (ndotproduct) between all consensus MS2
#' spectra whose precursor m/z differs by at most `precursor_tol` Da. Edges with
#' similarity >= `min_score` are retained. For [NTSProject] input, the graph and
#' edge table are **returned only** (nothing is written into the project).
#'
#' Top-K filtering: each node keeps at most `top_k` highest-scoring edges,
#' preventing hub-dominated topologies.
#'
#' @param project An [NTSProject] with non-empty `@clean_ms2` and
#'   `@ms2_annotations` (run [purify_ms2_spectra] first), or a named list of
#'   spectra (two-column matrices or list elements with `mz`/`intensity` and
#'   optional `precursorMz`).
#' @param precursor_tol Maximum precursor m/z difference (Da) to attempt
#'   pairwise scoring (default `1.0`).
#' @param min_score Minimum cosine similarity to create an edge (default `0.6`).
#' @param top_k Max edges per node (default `10`).
#' @param min_peaks Minimum number of fragment rows required per spectrum
#'   (after removing non-finite / non-positive intensities; default `2`).
#'
#' @return Named `list` with `igraph`, `edges` (`data.frame` of pairwise scores),
#'   and `nodes` (precursor `data.frame`: `name`, `mz`, `id`, optional `SMILES`,
#'   `Compound.name` merged from `@db_annotations` / `@ms2_annotations` /
#'   `library_SMILES` on `@ms1_annotations`) when `project` is an
#'   [NTSProject]; an [igraph::graph] only when `project` is a plain named list
#'   of spectra (legacy).
#' @export
build_molecular_network <- function(project,
                                    precursor_tol = 1.0,
                                    min_score = 0.6,
                                    top_k = 10L,
                                    min_peaks = 2L) {
  if (inherits(project, "NTSProject")) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("data.table is required", call. = FALSE)
    }
    clean <- methods::slot(project, "clean_ms2")
    m2 <- get_ms2_annotations(project)
    if (!length(clean) || is.null(names(clean))) {
      stop("clean_ms2 is empty or unnamed; run purify_ms2_spectra() first.", call. = FALSE)
    }
    if (!is.data.frame(m2) || !nrow(m2) || !"name" %in% names(m2)) {
      stop("ms2_annotations is empty or lacks `name`", call. = FALSE)
    }
    if (!"mz" %in% names(m2)) {
      stop("ms2_annotations must contain column `mz` (precursor)", call. = FALSE)
    }
    valid_ids <- intersect(names(clean), as.character(m2$name))
    valid_ids <- valid_ids[!is.na(valid_ids) & nzchar(valid_ids)]
    if (length(valid_ids) < 2L) {
      stop(
        "Fewer than 2 features after intersect(names(clean_ms2), ms2_annotations$name).",
        call. = FALSE
      )
    }
    spec_list <- .mn_spec_list_from_clean_ms2(clean, valid_ids, min_peaks = min_peaks)
    if (length(spec_list) < 2L) {
      stop(
        "Fewer than 2 spectra pass min_peaks (and non-empty consensus peaks); lower min_peaks or check clean_ms2.",
        call. = FALSE
      )
    }
    m2dt <- data.table::as.data.table(m2)
    m2dt[, name := as.character(name)]
    prec_df <- unique(
      m2dt[name %in% names(spec_list), .(name, mz = suppressWarnings(as.numeric(mz)))],
      by = "name"
    )
    if (any(!is.finite(prec_df$mz))) {
      stop("Non-finite precursor mz in ms2_annotations for some network nodes.", call. = FALSE)
    }
    res <- .build_mn_from_spectra(
      spec_list,
      precursor_tol = precursor_tol,
      min_score = min_score,
      top_k = top_k,
      min_peaks = min_peaks,
      precursor_df = prec_df,
      precursors_from_df_only = TRUE
    )
    prec_df_nodes <- .mn_enrich_nodes_with_library_seeds(project, prec_df)
    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "build_molecular_network",
      list(
        precursor_tol = precursor_tol,
        min_score = min_score,
        top_k = top_k,
        min_peaks = min_peaks,
        n_nodes = igraph::vcount(res$graph),
        n_edges = igraph::ecount(res$graph)
      )
    )
    return(invisible(list(
      igraph = res$graph,
      graph = res$graph,
      edges = res$edges,
      nodes = as.data.frame(prec_df_nodes, stringsAsFactors = FALSE)
    )))
  }
  if (is.list(project) && !inherits(project, "NTSProject")) {
    res <- .build_mn_from_spectra(
      project,
      precursor_tol = precursor_tol,
      min_score = min_score,
      top_k = top_k,
      min_peaks = min_peaks,
      precursor_df = NULL,
      precursors_from_df_only = FALSE
    )
    return(res$graph)
  }
  stop("`project` must be NTSProject or a named list of spectra", call. = FALSE)
}

# ============================================================================
# Advanced network builders (decoupled from plotting; return nodes / edges)
# ============================================================================

#' MS2 similarity network (thin wrapper)
#'
#' Calls [build_molecular_network] and returns only `list(nodes, edges)` for use with
#' [build_homologue_network] and [extract_targeted_subnetwork].
#'
#' @inheritParams build_molecular_network
#'
#' @return `list(nodes = data.frame, edges = data.frame)` with the same elements as
#'   in `build_molecular_network`.
#' @export
build_ms2_network <- function(project,
                              precursor_tol = 1.0,
                              min_score = 0.6,
                              top_k = 10L,
                              min_peaks = 2L) {
  res <- build_molecular_network(
    project,
    precursor_tol = precursor_tol,
    min_score = min_score,
    top_k = top_k,
    min_peaks = min_peaks
  )
  if (!is.list(res) || !all(c("nodes", "edges") %in% names(res))) {
    stop(
      "build_ms2_network() requires an NTSProject; for a bare spectrum list use build_molecular_network().",
      call. = FALSE
    )
  }
  list(nodes = res$nodes, edges = res$edges)
}

#' @noRd
.nts_homologue_transformation_tag <- function(mass_difference,
                                              ref_mz,
                                              target_mz_diff,
                                              tol_ppm,
                                              target_unit_name,
                                              max_gap = 3L) {
  unit <- trimws(as.character(target_unit_name)[1L])
  if (!nzchar(unit) || is.na(unit)) {
    unit <- "UNIT"
  }
  tmd0 <- suppressWarnings(as.numeric(target_mz_diff)[1L])
  .nts_check(abs(tmd0), "positive_finite", .arg = "target_mz_diff")
  tmd_mag <- abs(tmd0)
  .nts_check(tol_ppm, "positive_finite", .arg = "tol_ppm")
  tol <- suppressWarnings(as.numeric(tol_ppm)[1L])
  mg <- max(1L, as.integer(max_gap)[1L])
  if (!is.finite(mg)) {
    mg <- 1L
  }
  dm <- suppressWarnings(as.numeric(mass_difference))
  ref <- suppressWarnings(as.numeric(ref_mz))
  n <- length(dm)
  if (!length(ref)) {
    ref <- rep(NA_real_, n)
  } else if (length(ref) != n) {
    ref <- rep(ref, length.out = n)
  }
  out <- rep("Unknown", n)
  for (i in seq_len(n)) {
    d <- dm[i]
    r <- ref[i]
    if (!is.finite(d) || d == 0) {
      next
    }
    n_steps <- as.integer(round(d / tmd_mag))
    if (n_steps == 0L || abs(n_steps) > mg) {
      next
    }
    exp_dm <- as.numeric(n_steps) * tmd_mag
    tol_da <- abs(r) * (tol * 1e-6) * max(1L, abs(n_steps))
    if (!is.finite(tol_da) || tol_da < 1e-9) {
      tol_da <- 1e-6
    }
    if (abs(d - exp_dm) > tol_da) {
      next
    }
    out[i] <- sprintf("%+d%s", n_steps, unit)
  }
  out
}

#' @noRd
.nts_homologue_edges_append_reverse <- function(E,
                                                mz_map,
                                                target_mz_diff,
                                                tol_ppm,
                                                target_unit_name,
                                                max_gap = 3L) {
  if (!nrow(E) || !all(c("from", "to", "mass_difference") %in% names(E))) {
    return(E)
  }
  Er <- data.table::copy(E)
  tmp <- Er[["from"]]
  Er[, from := as.character(to)]
  Er[, to := as.character(tmp)]
  Er[, mass_difference := -suppressWarnings(as.numeric(mass_difference))]
  ref_mz <- pmax(
    abs(suppressWarnings(as.numeric(mz_map[Er[["from"]]]))),
    abs(suppressWarnings(as.numeric(mz_map[Er[["to"]]]))),
    1e-9,
    na.rm = TRUE
  )
  Er[, transformation := .nts_homologue_transformation_tag(
    mass_difference,
    ref_mz,
    target_mz_diff,
    tol_ppm,
    target_unit_name,
    max_gap = max_gap
  )]
  # Recompute the signed step count and display `label` for the reversed edge:
  # `mass_difference` was negated above, so `n_steps`/`label` must follow suit
  # (they are otherwise inherited unchanged from the forward edge, yielding
  # e.g. a reverse edge with mass_difference < 0 but n_steps > 0 -- P1-E2).
  if ("n_steps" %in% names(Er)) {
    tmd <- suppressWarnings(as.numeric(target_mz_diff)[1L])
    Er[, n_steps := {
      ns <- if (is.finite(tmd) && tmd != 0) {
        suppressWarnings(as.integer(round(mass_difference / tmd)))
      } else {
        NA_integer_
      }
      ifelse(transformation == "Unknown", NA_integer_, ns)
    }]
  }
  if ("label" %in% names(Er)) {
    Er[, label := ifelse(transformation == "Unknown", "", transformation)]
  }
  out <- data.table::rbindlist(list(E, Er), use.names = TRUE, fill = TRUE)
  unique(out, by = c("from", "to"))
}

#' Homologue / KMD-chain network (MS1 layer only)
#'
#' Builds MS1 homologue edges from [get_network_data] with `type = \"homologue\"` (unless
#' `network_data` is supplied), then optionally filters edges by target mass shift
#' `target_mz_diff` (step magnitude, Da) and `tol_ppm`. Each edge gets signed
#' `mass_difference` (`to_mz - from_mz`) and a `transformation` label (e.g. `+1CH2`,
#' `-2CF2`) derived from `target_unit_name`, `target_mz_diff`, and integer
#' `n_steps = round(mass_difference / target_mz_diff)` (gap-bridging:
#' \(|n_steps| \\le\\) `max_gap`). **No** fixed PFAS dictionary is required.
#'
#' **Pipeline:** `get_network_data(type = \"homologue\")` -> (optional) pass as
#' `network_data` here to avoid rebuilding KMD edges ->
#' [infer_homologue_structures] (requires edges with `transformation` from this
#' function, not raw KMD tables alone).
#'
#' @param project An [NTSProject].
#' @param network_data Optional `list(nodes, edges)` from a prior
#'   [get_network_data] call with `type = \"homologue\"`. When set, KMD homologue
#'   edges are **not** recomputed; only CF2/CH2 filtering and `transformation`
#'   labelling run. Arguments in `...` are then ignored.
#' @param target_mz_diff Single-step mass (Da, positive), e.g. `14.0156` (CH2), `49.9968064` (CF2).
#'   Residual between `mass_difference` and `n_steps * target_mz_diff` must fall within the ppm-derived tolerance.
#' @param target_unit_name Unit suffix written into `transformation` (e.g. `"CH2"`, `"CF2"`),
#'   combined with step count via `sprintf`, e.g. `+2CF2`, `-1CH2`.
#'   When calling [infer_homologue_structures], pass the same tokens in
#'   `transformation_units` (defaults include `CF2`/`CH2`/`O`/`H2`); unknown labels are ignored by the algebra layer.
#' @param tol_ppm Symmetric ppm tolerance: reference `max(|from_mz|,|to_mz|)`, tolerance in Da is
#'   `ref * tol_ppm * 1e-6 * max(1, |n_steps|)` to match `mass_difference` against
#'   `n_steps * target_mz_diff`.
#' @param filter_by_target_mass If `FALSE`, keep all homologue edges (still no MS2 layer),
#'   but still compute `mass_difference` and `transformation` for each edge.
#' @param bidirectional If `TRUE` (default), append reverse edges (swap
#'   `from`/`to`, negate `mass_difference`, relabel `transformation`) for
#'   undirected walks or bidirectional inference; duplicates removed with
#'   [data.table::unique]. The default is `TRUE` because downstream
#'   seed-guided expansion ([infer_homologue_structures]) otherwise only walks
#'   the stored edge direction: a seed sitting mid-chain would not reach its
#'   "upstream" ancestors, producing incomplete GEN-1/GEN-2 local libraries.
#'   Set `FALSE` only for a strictly directional store.
#' @param max_gap Maximum absolute step count allowed (integer, `>= 1`). E.g. `3L` accepts
#'   mass shifts consistent with step magnitudes 1-3 (gaps in the series may still link).
#'   Placed after `bidirectional` for backward-compatible positional argument order.
#' @param ... Forwarded to [get_network_data] (e.g. `mz_tol_ppm`, `isf_mz_tol_ppm`).
#'
#' @return `list(nodes = data.table, edges = data.table)`.
#' @export
build_homologue_network <- function(project,
                                    network_data = NULL,
                                    target_mz_diff = 49.9968064,
                                    target_unit_name = "CF2",
                                    tol_ppm = 20,
                                    filter_by_target_mass = TRUE,
                                    bidirectional = TRUE,
                                    max_gap = 3L,
                                    ...) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject.", call. = FALSE)
  }
  if (!is.null(network_data)) {
    if (!is.list(network_data) || !all(c("nodes", "edges") %in% names(network_data))) {
      stop(
        "`network_data` must be a list with elements `nodes` and `edges`.",
        call. = FALSE
      )
    }
    if (length(list(...))) {
      warning(
        "build_homologue_network: `network_data` supplied; `...` is ignored (not passed to get_network_data).",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    net <- network_data
  } else {
    net <- get_network_data(project, type = "homologue", ...)
  }
  N <- data.table::as.data.table(net$nodes)
  E <- data.table::as.data.table(net$edges)
  if (!nrow(N) || !"id" %in% names(N)) {
    return(list(nodes = N, edges = E))
  }
  if (!"mz" %in% names(N)) {
    warning(
      "build_homologue_network: nodes lack mz; skipping target_mz_diff filter and returning unfiltered network.",
      call. = FALSE,
      immediate. = TRUE
    )
    return(list(nodes = N, edges = E))
  }
  N[, id := trimws(as.character(id))]
  N[, mz := suppressWarnings(as.numeric(mz))]
  mz_map <- stats::setNames(N$mz, N$id)

  if (!nrow(E) || !isTRUE(filter_by_target_mass)) {
    if (nrow(E)) {
      if (all(c("from", "to") %in% names(E))) {
        E[, from := trimws(as.character(from))]
        E[, to := trimws(as.character(to))]
        E[, from_mz := suppressWarnings(as.numeric(mz_map[from]))]
        E[, to_mz := suppressWarnings(as.numeric(mz_map[to]))]
        E[, mass_difference := to_mz - from_mz]
        ref_mz <- pmax(abs(E$from_mz), abs(E$to_mz), 1e-9, na.rm = TRUE)
        E[, transformation := .nts_homologue_transformation_tag(
          mass_difference,
          ref_mz,
          target_mz_diff,
          tol_ppm,
          target_unit_name,
          max_gap = max_gap
        )]
        E[, c("from_mz", "to_mz") := NULL]
        if (isTRUE(bidirectional)) {
          E <- .nts_homologue_edges_append_reverse(
            E,
            mz_map,
            target_mz_diff,
            tol_ppm,
            target_unit_name,
            max_gap = max_gap
          )
        }
      }
      N <- .nts_network_prune_isolated_nodes_homologue(N, E)
    }
    return(list(nodes = N, edges = E))
  }
  if (!all(c("from", "to") %in% names(E))) {
    return(list(nodes = N, edges = E))
  }
  E[, from := trimws(as.character(from))]
  E[, to := trimws(as.character(to))]
  E[, from_mz := suppressWarnings(as.numeric(mz_map[from]))]
  E[, to_mz := suppressWarnings(as.numeric(mz_map[to]))]
  E[, mass_difference := to_mz - from_mz]
  tmd <- suppressWarnings(as.numeric(target_mz_diff)[1L])
  if (!is.finite(tmd) || tmd == 0) {
    stop("`target_mz_diff` must be finite and non-zero.", call. = FALSE)
  }
  tol <- suppressWarnings(as.numeric(tol_ppm)[1L])
  if (!is.finite(tol) || tol < 0) {
    stop("`tol_ppm` must be a non-negative finite number.", call. = FALSE)
  }
  ref_mz <- pmax(abs(E$from_mz), abs(E$to_mz), 1e-9, na.rm = TRUE)
  E[, transformation := .nts_homologue_transformation_tag(
    mass_difference,
    ref_mz,
    target_mz_diff,
    tol_ppm,
    target_unit_name,
    max_gap = max_gap
  )]
  E2 <- E[is.finite(from_mz) & is.finite(to_mz) & transformation != "Unknown"]
  E2[, c("from_mz", "to_mz") := NULL]
  if (isTRUE(bidirectional) && nrow(E2)) {
    E2 <- .nts_homologue_edges_append_reverse(
      E2,
      mz_map,
      target_mz_diff,
      tol_ppm,
      target_unit_name,
      max_gap = max_gap
    )
  }
  if (nrow(E2)) {
    N <- .nts_network_prune_isolated_nodes_homologue(N, E2)
  } else {
    N <- N[0L]
  }
  list(nodes = N, edges = E2)
}

#' @noRd
.nts_network_prune_isolated_nodes_homologue <- function(Nodes, Edges) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  if (!nrow(Edges) || !all(c("from", "to") %in% names(Edges))) {
    return(Nodes[0L])
  }
  u <- unique(c(as.character(Edges$from), as.character(Edges$to)))
  Nodes[as.character(Nodes$id) %in% u]
}

#' Extract a targeted subnetwork from network data by name
#'
#' Matches `target_name` to `nodes$id`, `nodes$name`, or `Compound.name` (if present)
#' after `trimws`, then calls [extract_subnetwork].
#'
#' @param network_data Named list with `nodes` and `edges` (e.g. from [get_network_data]).
#' @param target_name One or more node identifiers (aligned with id / name / Compound.name).
#' @param max_depth Passed to [extract_subnetwork] as `max_steps`.
#'
#' @return `list(nodes = data.table, edges = data.table)`.
#' @export
extract_targeted_subnetwork <- function(network_data,
                                        target_name,
                                        max_depth = 2L) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  if (!is.list(network_data) || !all(c("nodes", "edges") %in% names(network_data))) {
    stop("network_data must be a list with nodes and edges.", call. = FALSE)
  }
  nd <- as.data.frame(network_data$nodes, stringsAsFactors = FALSE)
  if (!nrow(nd)) {
    return(list(
      nodes = data.table::as.data.table(nd),
      edges = data.table::as.data.table(network_data$edges)[0L]
    ))
  }
  want <- unique(trimws(as.character(target_name)))
  want <- want[!is.na(want) & nzchar(want)]
  if (!length(want)) {
    stop("`target_name` resolves to an empty set.", call. = FALSE)
  }
  hit <- rep(FALSE, nrow(nd))
  if ("id" %in% names(nd)) {
    hit <- hit | (trimws(as.character(nd$id)) %in% want)
  }
  if ("name" %in% names(nd)) {
    hit <- hit | (trimws(as.character(nd$name)) %in% want)
  }
  if ("Compound.name" %in% names(nd)) {
    hit <- hit | (trimws(as.character(nd$Compound.name)) %in% want)
  }
  ids <- character(0)
  if (any(hit)) {
    if ("id" %in% names(nd)) {
      ids <- unique(as.character(nd$id[hit]))
    } else {
      ids <- unique(as.character(nd$name[hit]))
    }
  }
  if (!length(ids)) {
    warning(
      "extract_targeted_subnetwork: no nodes matched target_name among id/name/Compound.name.",
      call. = FALSE,
      immediate. = TRUE
    )
    return(list(
      nodes = data.table::as.data.table(nd)[0L],
      edges = data.table::as.data.table(network_data$edges)[0L]
    ))
  }
  extract_subnetwork(network_data, ids, max_steps = max_depth)
}

#' @noRd
.mn_spec_list_from_clean_ms2 <- function(clean, ids, min_peaks) {
  out <- list()
  min_peaks <- as.integer(min_peaks)[1L]
  if (!is.finite(min_peaks) || min_peaks < 1L) {
    min_peaks <- 1L
  }
  for (id in ids) {
    el <- clean[[id]]
    if (!is.list(el)) {
      next
    }
    spf <- el[["spectrum"]]
    if (!is.data.frame(spf) || !nrow(spf)) {
      next
    }
    mz <- suppressWarnings(as.numeric(spf$mz))
    inten <- if ("intensity" %in% names(spf)) {
      suppressWarnings(as.numeric(spf$intensity))
    } else if ("frequency" %in% names(spf)) {
      suppressWarnings(as.numeric(spf$frequency))
    } else {
      NULL
    }
    if (is.null(inten)) {
      next
    }
    ok <- is.finite(mz) & is.finite(inten) & inten > 0
    mz <- mz[ok]
    inten <- inten[ok]
    if (length(mz) < min_peaks) {
      next
    }
    out[[id]] <- cbind(mz = mz, intensity = inten)
  }
  out
}

#' Enrich molecular-network `nodes` (precursor table) with library SMILES / names
#' from `@db_annotations` (Top-1 per `name`), then `@ms2_annotations`, then
#' optional `library_SMILES` on `@ms1_annotations`. Column **`SMILES`** is
#' always uppercase for downstream seed lookup.
#' @noRd
.mn_enrich_nodes_with_library_seeds <- function(project, prec_df) {
  out <- as.data.frame(prec_df, stringsAsFactors = FALSE)
  if (!nrow(out) || !"name" %in% names(out)) {
    return(out)
  }
  out$name <- trimws(as.character(out$name))
  if (!"id" %in% names(out)) {
    out$id <- out$name
  } else {
    out$id <- trimws(as.character(out$id))
  }
  out$SMILES <- NA_character_
  out$Compound.name <- NA_character_
  if (!methods::is(project, "NTSProject")) {
    return(out)
  }
  db <- .nts_db_best_per_name(methods::slot(project, "db_annotations"))
  if (is.data.frame(db) && nrow(db) && "name" %in% names(db)) {
    db$name <- trimws(as.character(db$name))
    j <- match(out$name, db$name)
    if ("SMILES" %in% names(db)) {
      out$SMILES <- as.character(db$SMILES)[j]
    } else if ("smiles" %in% names(db)) {
      out$SMILES <- as.character(db$smiles)[j]
    }
    if ("Compound.name" %in% names(db)) {
      out$Compound.name <- as.character(db$Compound.name)[j]
    } else if ("compound_name" %in% names(db)) {
      out$Compound.name <- as.character(db$compound_name)[j]
    }
  }
  m2 <- as.data.frame(get_ms2_annotations(project), stringsAsFactors = FALSE)
  if (is.data.frame(m2) && nrow(m2) && "name" %in% names(m2)) {
    m2$name <- trimws(as.character(m2$name))
    sm2 <- rep(NA_character_, nrow(m2))
    if ("SMILES" %in% names(m2)) {
      sm2 <- trimws(as.character(m2$SMILES))
    } else if ("smiles" %in% names(m2)) {
      sm2 <- trimws(as.character(m2$smiles))
    }
    jm <- match(out$name, m2$name)
    sm2m <- sm2[jm]
    wna <- is.na(out$SMILES) | !nzchar(as.character(out$SMILES))
    okm <- wna & !is.na(sm2m) & nzchar(as.character(sm2m))
    if (any(okm, na.rm = TRUE)) {
      out$SMILES[okm] <- as.character(sm2m[okm])
    }
    cna <- is.na(out$Compound.name) | !nzchar(as.character(out$Compound.name))
    if ("Compound.name" %in% names(m2)) {
      cm2 <- trimws(as.character(m2$Compound.name)[jm])
      okc <- cna & !is.na(cm2) & nzchar(cm2)
      if (any(okc, na.rm = TRUE)) {
        out$Compound.name[okc] <- as.character(cm2[okc])
      }
    }
  }
  ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
  if (is.data.frame(ms1) && nrow(ms1) && "name" %in% names(ms1) &&
      "library_SMILES" %in% names(ms1)) {
    ms1$name <- trimws(as.character(ms1$name))
    jl <- match(out$name, ms1$name)
    lsmi <- as.character(ms1$library_SMILES)[jl]
    wls <- (is.na(out$SMILES) | !nzchar(as.character(out$SMILES))) &
      !is.na(lsmi) & nzchar(lsmi)
    if (any(wls, na.rm = TRUE)) {
      out$SMILES[wls] <- lsmi[wls]
    }
  }
  out$SMILES[is.na(out$SMILES)] <- NA_character_
  if (!"library_Compound_name" %in% names(out)) {
    out$library_Compound_name <- NA_character_
  }
  if (!"library_ref_id" %in% names(out)) {
    out$library_ref_id <- NA_character_
  }
  if (methods::is(project, "NTSProject")) {
    ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
    if (is.data.frame(ms1) && nrow(ms1) && "name" %in% names(ms1)) {
      ms1$name <- trimws(as.character(ms1$name))
      jm <- match(out$name, ms1$name)
      if ("library_Compound_name" %in% names(ms1)) {
        out$library_Compound_name <- as.character(ms1$library_Compound_name)[jm]
      }
      if ("library_ref_id" %in% names(ms1)) {
        out$library_ref_id <- as.character(ms1$library_ref_id)[jm]
      }
    }
  }
  out
}

#' Nominal exact-mass deltas (Da) for common repeat units / heteroatoms.
#' @noRd
.mn_precursor_delta_rule_table <- function() {
  data.frame(
    tag = c("CF2", "CH2", "O", "H2"),
    mass_da = c(.nts_pfas_cf2_da, 14.0156, 15.9949, 2.0156),
    stringsAsFactors = FALSE
  )
}

#' Map signed precursor mass difference to `+CF2` / `-O` / `Unknown` (ppm tol).
#' @noRd
.mn_delta_mz_to_transformation <- function(delta_mz, mz_ref, tol_ppm = 20) {
  d <- suppressWarnings(as.numeric(delta_mz)[1L])
  ref <- suppressWarnings(as.numeric(mz_ref)[1L])
  if (!is.finite(d) || d == 0) {
    return("Unknown")
  }
  if (!is.finite(ref) || ref <= 0) {
    ref <- abs(d)
  }
  tol_da <- abs(ref) * (as.numeric(tol_ppm)[1L] * 1e-6)
  if (!is.finite(tol_da) || tol_da <= 0) {
    tol_da <- 0.01
  }
  tab <- .mn_precursor_delta_rule_table()
  ad <- abs(d)
  hit <- which.min(abs(ad - tab$mass_da))
  if (!length(hit)) {
    return("Unknown")
  }
  if (abs(ad - tab$mass_da[hit]) > tol_da) {
    return("Unknown")
  }
  unit <- as.character(tab$tag[hit])
  if (!nzchar(unit)) {
    return("Unknown")
  }
  if (d > 0) {
    paste0("+", unit)
  } else {
    paste0("-", unit)
  }
}

#' @noRd
.build_mn_from_spectra <- function(spec_list,
                                   precursor_tol,
                                   min_score,
                                   top_k,
                                   min_peaks = 2L,
                                   precursor_df = NULL,
                                   precursors_from_df_only = FALSE) {
  .nts_check(NULL, "igraph_ns", .arg = "igraph")
  min_peaks <- as.integer(min_peaks)[1L]
  if (!is.finite(min_peaks) || min_peaks < 1L) {
    min_peaks <- 1L
  }
  ids <- names(spec_list)
  if (is.null(ids)) {
    ids <- as.character(seq_along(spec_list))
  }
  n <- length(spec_list)

  mat_list <- vector("list", n)
  prec_mz  <- rep(NA_real_, n)
  for (k in seq_len(n)) {
    s <- spec_list[[k]]
    if (is.matrix(s) && ncol(s) >= 2L) {
      mz <- suppressWarnings(as.numeric(s[, 1L]))
      inten <- suppressWarnings(as.numeric(s[, 2L]))
      ok <- is.finite(mz) & is.finite(inten) & inten > 0
      mz <- mz[ok]
      inten <- inten[ok]
      if (length(mz) < min_peaks) {
        mat_list[[k]] <- NULL
        next
      }
      mat_list[[k]] <- cbind(mz = mz, intensity = inten)
      prec_mz[k] <- if (!is.null(attr(s, "precursorMz"))) {
        suppressWarnings(as.numeric(attr(s, "precursorMz"))[1L])
      } else {
        NA_real_
      }
    } else if (is.list(s)) {
      m <- s$mz
      i <- s$intensity
      sp_mat <- s$spectrum
      if (!is.null(sp_mat) && is.matrix(sp_mat)) {
        mz <- suppressWarnings(as.numeric(sp_mat[, 1L]))
        inten <- suppressWarnings(as.numeric(sp_mat[, 2L]))
      } else if (!is.null(m) && !is.null(i)) {
        mz <- suppressWarnings(as.numeric(m))
        inten <- suppressWarnings(as.numeric(i))
      } else {
        mat_list[[k]] <- NULL
        next
      }
      ok <- is.finite(mz) & is.finite(inten) & inten > 0
      mz <- mz[ok]
      inten <- inten[ok]
      if (length(mz) < min_peaks) {
        mat_list[[k]] <- NULL
        next
      }
      mat_list[[k]] <- cbind(mz = mz, intensity = inten)
      prec_mz[k] <- suppressWarnings(as.numeric(
        if (!is.null(s$precursorMz)) {
          s$precursorMz[1L]
        } else if (!is.null(s$precursor)) {
          s$precursor[1L]
        } else {
          NA_real_
        }
      ))
    } else {
      mat_list[[k]] <- NULL
    }
  }

  if (isTRUE(precursors_from_df_only)) {
    if (!is.data.frame(precursor_df) || !nrow(precursor_df) ||
        !all(c("name", "mz") %in% names(precursor_df))) {
      stop(
        "precursors_from_df_only=TRUE requires precursor_df with columns name, mz",
        call. = FALSE
      )
    }
    pm <- suppressWarnings(as.numeric(precursor_df$mz)[match(ids, as.character(precursor_df$name))])
    prec_mz <- pm
  } else if (!is.null(precursor_df) && is.data.frame(precursor_df) && nrow(precursor_df) &&
             "mz" %in% names(precursor_df) && "name" %in% names(precursor_df)) {
    na_idx <- which(!is.finite(prec_mz))
    if (length(na_idx)) {
      m <- match(ids[na_idx], as.character(precursor_df$name))
      ok <- !is.na(m)
      prec_mz[na_idx[ok]] <- suppressWarnings(as.numeric(precursor_df$mz[m[ok]]))
    }
  }

  keep <- !vapply(mat_list, is.null, logical(1L)) & is.finite(prec_mz)
  edges_empty <- data.frame(
    from = character(0),
    to = character(0),
    weight = numeric(0),
    delta_mz = numeric(0),
    transformation = character(0),
    stringsAsFactors = FALSE
  )
  if (sum(keep) < 2L) {
    message("> Fewer than 2 usable spectra; returning empty network.")
    return(list(
      graph = igraph::make_empty_graph(directed = FALSE),
      edges = edges_empty
    ))
  }

  idx <- which(keep)
  sub_spec <- mat_list[idx]
  sub_prec <- prec_mz[idx]
  sub_ids  <- ids[idx]

  message("> Computing MS2 similarity network (", length(sub_spec), " nodes, precursor_tol=",
          precursor_tol, " Da, min_peaks=", min_peaks, ")...")

  edges <- .mn_pairwise_scores(sub_spec, sub_prec, precursor_tol, min_score)

  if (nrow(edges)) {
    edges[, from := sub_ids[i]]
    edges[, to   := sub_ids[j]]

    dt_long <- rbind(
      edges[, .(node = from, partner = to, score)],
      edges[, .(node = to, partner = from, score)]
    )
    dt_long <- dt_long[order(-score)]
    dt_long[, rk := seq_len(.N), by = node]
    dt_long <- dt_long[rk <= top_k]

    edge_set <- unique(
      data.table::data.table(
        a = pmin(dt_long$node, dt_long$partner),
        b = pmax(dt_long$node, dt_long$partner),
        score = dt_long$score
      )
    )

    es <- edge_set[, .(from = a, to = b, weight = score)]
    es[["from"]] <- as.character(es[["from"]])
    es[["to"]] <- as.character(es[["to"]])
    ok_e <- !is.na(es[["from"]]) & !is.na(es[["to"]]) &
      nzchar(es[["from"]]) & nzchar(es[["to"]])
    if (length(ok_e) && !all(ok_e)) {
      es <- es[ok_e]
    }
    if (!nrow(es)) {
      edges_df <- edges_empty
      g <- igraph::make_empty_graph(directed = FALSE)
      g <- igraph::add_vertices(g, length(sub_ids), name = sub_ids)
    } else {
      mf <- sub_prec[match(es[["from"]], sub_ids)]
      mt <- sub_prec[match(es[["to"]], sub_ids)]
      es[["delta_mz"]] <- mt - mf
      es[["transformation"]] <- vapply(seq_len(nrow(es)), function(k) {
        ref_mz <- mean(c(mf[k], mt[k]), na.rm = TRUE)
        .mn_delta_mz_to_transformation(es[["delta_mz"]][k], ref_mz)
      }, character(1L), USE.NAMES = FALSE)
      edges_df <- data.frame(
        from = es[["from"]],
        to = es[["to"]],
        weight = es[["weight"]],
        delta_mz = es[["delta_mz"]],
        transformation = es[["transformation"]],
        stringsAsFactors = FALSE
      )
      es_igraph <- data.frame(
        from = es[["from"]],
        to = es[["to"]],
        weight = es[["weight"]],
        delta_mz = es[["delta_mz"]],
        transformation = es[["transformation"]],
        stringsAsFactors = FALSE
      )
      g <- igraph::graph_from_data_frame(
        es_igraph,
        directed = FALSE,
        vertices = data.frame(name = sub_ids, stringsAsFactors = FALSE)
      )
    }
  } else {
    edges_df <- edges_empty
    g <- igraph::make_empty_graph(directed = FALSE)
    g <- igraph::add_vertices(g, length(sub_ids), name = sub_ids)
  }

  igraph::V(g)$precursor_mz <- sub_prec[match(igraph::V(g)$name, sub_ids)]

  message("> Molecular network built: ", igraph::vcount(g), " nodes, ",
          igraph::ecount(g), " edges")
  list(graph = g, edges = edges_df)
}

# ============================================================================
# Network-based annotation / formula propagation (internal rules)
# ============================================================================

#' Parse a signed Hill-style formula difference into a net exact mass (Da)
#'
#' Single fragment (no \code{+}/\code{-}) is passed to [calculate_exact_mass].
#' Strings like \code{"H-F+O"} are normalized to leading \code{+} then split into
#' tokens \code{[+-][A-Za-z0-9]+} (e.g. \code{+H}, \code{-F}, \code{+O}) and summed.
#'
#' @param formula_diff Character scalar.
#' @return Finite positive or negative mass, or \code{NA_real_}.
#' @noRd
.nts_formula_diff_exact_mass <- function(formula_diff) {
  fd <- trimws(as.character(formula_diff)[1L])
  if (is.na(fd) || !nzchar(fd)) {
    return(NA_real_)
  }
  if (!grepl("[+-]", fd, perl = TRUE)) {
    m <- calculate_exact_mass(fd)
    return(if (is.finite(m)) m else NA_real_)
  }
  if (!grepl("^[+-]", fd, perl = TRUE)) {
    fd <- paste0("+", fd)
  }
  toks <- regmatches(fd, gregexpr("[+-][A-Za-z0-9]+", fd, perl = TRUE))[[1L]]
  if (!length(toks)) {
    return(NA_real_)
  }
  tot <- 0
  for (t in toks) {
    sg <- substr(t, 1L, 1L)
    fm <- substring(t, 2L)
    if (!nzchar(fm)) {
      return(NA_real_)
    }
    mm <- calculate_exact_mass(fm)
    if (!is.finite(mm)) {
      return(NA_real_)
    }
    tot <- tot + if (identical(sg, "+")) mm else -mm
  }
  tot
}

#' Default transformation / homologue mass-difference rules for edge matching
#'
#' \code{formula_diff} uses optional \code{+}/\code{-} segments; \code{exact_mass} is computed
#' via [calculate_exact_mass] and `.nts_formula_diff_exact_mass`.
#'
#' @return \code{data.frame} with columns \code{rule_name}, \code{formula_diff}, \code{exact_mass}.
#' @noRd
.nts_default_transformations <- function() {
  rules <- data.frame(
    rule_name = c(
      "CF2_addition",
      "CH2_addition",
      "Oxygenation",
      "Defluorination_Hydroxylation"
    ),
    formula_diff = c("CF2", "CH2", "O", "H-F+O"),
    stringsAsFactors = FALSE
  )
  rules$exact_mass <- vapply(
    rules$formula_diff,
    function(x) {
      m <- .nts_formula_diff_exact_mass(x)
      if (!is.finite(m)) {
        NA_real_
      } else {
        m
      }
    },
    numeric(1L),
    USE.NAMES = FALSE
  )
  rules
}

#' @title Network-based formula propagation on masses
#'
#' @description
#' Breadth-first search from seed nodes with non-empty `formula` on `network_data$nodes`,
#' along `edges`, matching candidate mass differences (precursor \code{m/z} deltas and
#' optional \code{shared_mds_mz} values) to `rules$exact_mass` within ppm tolerance, and
#' applying `rules$formula_diff` to infer `propagated_formula` on unknown nodes.
#'
#' @param network_data Named list containing `nodes` and `edges` (e.g. outputs of
#'   [get_network_data] or [generate_network_data]).
#' @param rules Rule table with columns `rule_name`, `formula_diff`, `exact_mass`; if `NULL`,
#'   `.nts_default_transformations` is used.
#' @param mz_tol_ppm ppm tolerance for mass-difference matching.
#' @param max_steps Maximum BFS depth (edge count from seeds).
#'
#' @return Named list: `nodes` (`data.table` with added `propagated_formula`,
#'   `propagation_rule`, `propagation_depth`), `edges` (copy of input edges), and
#'   `stats` (e.g. `n_seeds`, `n_propagated`).
#' @export
propagate_network_annotations <- function(network_data,
                                          rules = NULL,
                                          mz_tol_ppm = 10,
                                          max_steps = 3L) {
  if (!is.list(network_data) || !all(c("nodes", "edges") %in% names(network_data))) {
    stop("network_data must be a list with elements `nodes` and `edges`.", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  N <- data.table::as.data.table(network_data$nodes)
  E <- data.table::as.data.table(network_data$edges)
  if (!nrow(N)) {
    return(list(
      nodes = N,
      edges = if (nrow(E)) data.table::copy(E) else E,
      stats = list(n_seeds = 0L, n_propagated = 0L, n_queue_peak = 0L)
    ))
  }
  if (!"id" %in% names(N) && "name" %in% names(N)) {
    N[, id := as.character(name)]
  }
  if (!"id" %in% names(N)) {
    stop("nodes must contain column `id` (or `name` to alias id).", call. = FALSE)
  }
  if (!"mz" %in% names(N)) {
    stop("nodes must contain numeric column `mz`.", call. = FALSE)
  }
  if (!"formula" %in% names(N)) {
    N[, formula := NA_character_]
  }
  N[, id := as.character(id)]
  N[, mz := suppressWarnings(as.numeric(mz))]
  N[, formula_chr := ifelse(
    is.na(formula) | !nzchar(trimws(as.character(formula))),
    NA_character_,
    trimws(as.character(formula))
  )]

  if (is.null(rules)) {
    rules <- .nts_default_transformations()
  }
  if (!is.data.frame(rules)) {
    stop("rules must be a data.frame.", call. = FALSE)
  }
  rules <- as.data.frame(rules, stringsAsFactors = FALSE)
  required_cols <- c("rule_name", "formula_diff", "exact_mass")
  if (!all(required_cols %in% colnames(rules))) {
    stop("rules must have columns: rule_name, formula_diff, exact_mass.", call. = FALSE)
  }

  max_steps <- max(0L, as.integer(max_steps)[1L])
  ppm <- suppressWarnings(as.numeric(mz_tol_ppm)[1L])
  if (!is.finite(ppm) || ppm < 0) {
    stop("mz_tol_ppm must be non-negative finite.", call. = FALSE)
  }

  if (!nrow(E) || !all(c("from", "to") %in% names(E))) {
    N[, propagated_formula := NA_character_]
    N[, propagation_rule := NA_character_]
    N[, propagation_depth := NA_integer_]
    N[, formula_chr := NULL]
    n_seed <- sum(!is.na(N$formula_chr) & nzchar(N$formula_chr), na.rm = TRUE)
    return(list(
      nodes = N,
      edges = if (nrow(E)) data.table::copy(E) else E,
      stats = list(
        n_seeds = as.integer(n_seed),
        n_propagated = 0L,
        n_queue_peak = 0L
      )
    ))
  }
  E[, from := as.character(from)]
  E[, to := as.character(to)]

  adj <- data.table::rbindlist(
    list(
      E[, .(from, to)],
      E[, .(from = to, to = from)]
    ),
    use.names = TRUE
  )
  adj <- unique(adj, by = c("from", "to"))
  data.table::setkey(adj, from)

  E_lookup <- data.table::copy(E)
  E_lookup[, ek := paste(pmin(from, to), pmax(from, to), sep = "\t")]
  data.table::setkey(E_lookup, ek)

  N[, propagated_formula := NA_character_]
  N[, propagation_rule := NA_character_]
  N[, propagation_depth := NA_integer_]

  id_row <- stats::setNames(seq_len(nrow(N)), as.character(N$id))
  mz_id <- N$mz
  names(mz_id) <- as.character(N$id)

  seed_ids <- N$id[!is.na(N$formula_chr) & nzchar(N$formula_chr)]
  n_seed <- length(seed_ids)
  for (sid in seed_ids) {
    ir <- id_row[[as.character(sid)]]
    N[ir, `:=`(
      propagated_formula = formula_chr[ir],
      propagation_depth = 0L
    )]
  }

  Q <- vector("list", max(64L, nrow(N) * 4L))
  q_hi <- 0L
  q_push <- function(id, dep, fr) {
    q_hi <<- q_hi + 1L
    if (q_hi > length(Q)) {
      length(Q) <<- max(length(Q) * 2L, q_hi)
    }
    Q[[q_hi]] <<- list(id = id, dep = dep, fr = fr)
  }

  for (sid in seed_ids) {
    ir <- id_row[[as.character(sid)]]
    q_push(sid, 0L, as.character(N$formula_chr[ir]))
  }

  q_lo <- 1L
  q_peak <- q_hi
  n_prop <- 0L

  while (q_lo <= q_hi) {
    cur <- Q[[q_lo]]
    q_lo <- q_lo + 1L
    u <- cur$id
    du <- cur$dep
    fu <- cur$fr
    if (is.na(fu) || !nzchar(fu)) {
      next
    }
    if (du >= max_steps) {
      next
    }
    mu <- mz_id[[as.character(u)]]
    if (!is.finite(mu)) {
      next
    }
    nb <- unique(as.character(adj[from == u, to]))
    if (!length(nb)) {
      next
    }
    for (v in nb) {
      vc <- as.character(v)
      vr <- id_row[[vc]]
      if (is.na(vr)) {
        next
      }
      pv <- N$propagated_formula[vr]
      if (!is.na(pv) && nzchar(pv)) {
        next
      }
      mv <- mz_id[[vc]]
      if (!is.finite(mv)) {
        next
      }
      ek <- paste(sort(c(u, vc)), collapse = "\t")
      er <- E_lookup[list(ek), nomatch = NULL]
      er1 <- if (nrow(er) > 0L) er[1L] else NULL
      cands <- .nts_collect_signed_md_candidates(mu, mv, er1)
      if (!length(cands)) {
        next
      }
      ref <- max(abs(mu), abs(mv), 1e-6, na.rm = TRUE)
      best <- .nts_best_mass_rule_match(cands, rules, ref, ppm)
      if (is.null(best)) {
        next
      }
      fd <- as.character(rules$formula_diff[best$i])
      fv <- .nts_apply_formula_diff_direction(fu, fd, best$dir)
      if (is.na(fv) || !nzchar(fv)) {
        next
      }
      N[vr, `:=`(
        propagated_formula = fv,
        propagation_rule = as.character(rules$rule_name[best$i]),
        propagation_depth = as.integer(du + 1L)
      )]
      n_prop <- n_prop + 1L
      q_push(vc, du + 1L, fv)
    }
    q_peak <- max(q_peak, q_hi - q_lo + 1L, na.rm = TRUE)
  }

  N[, formula_chr := NULL]
  list(
    nodes = N,
    edges = data.table::copy(E),
    stats = list(
      n_seeds = as.integer(n_seed),
      n_propagated = as.integer(n_prop),
      n_queue_peak = as.integer(q_peak)
    )
  )
}

#' @noRd
.nts_collect_signed_md_candidates <- function(mz_u, mz_v, edge_row) {
  su <- suppressWarnings(as.numeric(mz_u)[1L])
  sv <- suppressWarnings(as.numeric(mz_v)[1L])
  if (!is.finite(su) || !is.finite(sv)) {
    return(numeric(0))
  }
  d0 <- sv - su
  sg <- if (!is.finite(d0) || d0 == 0) 1 else sign(d0)
  out <- d0
  if (!is.null(edge_row) && nrow(edge_row) && "shared_mds_mz" %in% names(edge_row)) {
    s <- as.character(edge_row[["shared_mds_mz"]])[1L]
    if (!is.na(s) && nzchar(s)) {
      s_norm <- gsub("\\s*\\+\\s*", ";", s)
      parts <- trimws(strsplit(s_norm, "[;,\n\t|]+", perl = TRUE)[[1L]])
      parts <- parts[vapply(parts, nzchar, logical(1L))]
      for (p in parts) {
        z <- suppressWarnings(as.numeric(p))
        if (is.finite(z)) {
          out <- c(out, sg * abs(z))
        }
      }
    }
  }
  unique(out[is.finite(out)])
}

#' @noRd
.nts_best_mass_rule_match <- function(candidates, rules, ref_mz, ppm_tol) {
  best <- NULL
  for (d in candidates) {
    if (!is.finite(d)) {
      next
    }
    for (i in seq_len(nrow(rules))) {
      m <- suppressWarnings(as.numeric(rules$exact_mass[i]))
      if (!is.finite(m) || m == 0) {
        next
      }
      e1 <- abs(d - m) / ref_mz * 1e6
      e2 <- abs(d + m) / ref_mz * 1e6
      e <- min(e1, e2)
      dir <- if (e1 <= e2) 1L else -1L
      if (e <= ppm_tol && (is.null(best) || e < best$err)) {
        best <- list(i = i, err = e, dir = dir)
      }
    }
  }
  best
}

#' @noRd
.nts_apply_formula_diff_direction <- function(formula, formula_diff, direction) {
  formula <- as.character(formula)[1L]
  direction <- as.integer(sign(direction[1L]))
  if (!direction %in% c(-1L, 1L)) {
    return(NA_character_)
  }
  if (is.na(formula) || !nzchar(trimws(formula))) {
    return(NA_character_)
  }
  fd <- trimws(as.character(formula_diff)[1L])
  if (is.na(fd) || !nzchar(fd)) {
    return(NA_character_)
  }
  if (!grepl("[+-]", fd, perl = TRUE)) {
    return(adjust_homologue_formula(formula, direction, fd))
  }
  fd2 <- if (!grepl("^[+-]", fd, perl = TRUE)) {
    paste0("+", fd)
  } else {
    fd
  }
  toks <- regmatches(fd2, gregexpr("[+-][A-Za-z0-9]+", fd2, perl = TRUE))[[1L]]
  if (!length(toks)) {
    return(NA_character_)
  }
  fcnt <- .parse_formula_element_counts(formula)
  if (!length(names(fcnt))) {
    return(NA_character_)
  }
  net <- structure(as.integer(fcnt), names = names(fcnt))
  for (t in toks) {
    s0 <- if (substr(t, 1L, 1L) == "+") {
      1L
    } else {
      -1L
    }
    s <- as.integer(s0 * direction)
    fm <- substring(t, 2L)
    if (!nzchar(fm)) {
      return(NA_character_)
    }
    bc <- .parse_formula_element_counts(fm)
    if (!length(names(bc))) {
      return(NA_character_)
    }
    for (e in names(bc)) {
      ad <- as.integer(bc[[e]]) * s
      if (e %in% names(net)) {
        net[[e]] <- net[[e]] + ad
      } else {
        net[e] <- ad
      }
    }
  }
  if (any(net < 0L, na.rm = TRUE)) {
    return(NA_character_)
  }
  net <- net[!is.na(net) & net > 0L]
  if (!length(net)) {
    return(NA_character_)
  }
  .format_sorted_formula(net)
}

#' Propagate annotations through molecular network (neighbour voting)
#'
#' For each unannotated node that has at least one annotated neighbour, assigns
#' the most frequent annotation label among its direct neighbours (1-hop). Ties
#' are broken by highest average edge weight. Optionally propagates up to
#' `max_hops` layers.
#'
#' The labels written are **inferred by neighbour voting**, not spectral library
#' matches. They are therefore persisted on `@ms1_annotations` (as a dedicated
#' `propagated_<label_col>` column) and **never** injected into
#' `@db_annotations`, which is the spectral-library layer that drives confidence
#' grading and the physchem proxy. Mixing inferred votes there would create
#' phantom "Database Confirmed" rows.
#'
#' @param project An [NTSProject] with `@db_annotations` (read for labels).
#' @param graph An [igraph::graph] from [build_molecular_network] (return value
#'   `$igraph` or `$graph`).
#' @param label_col Column in `@db_annotations` containing the label to propagate
#'   (default `"ref_id"`).
#' @param max_hops Maximum propagation depth (default `2`).
#' @param min_neighbours Minimum annotated neighbours required to propagate
#'   (default `2`).
#'
#' @return Updated [NTSProject] with propagated labels in `@ms1_annotations`
#'   (column `propagated_<label_col>`).
#' @export
propagate_annotations <- function(project,
                                  graph,
                                  label_col = "ref_id",
                                  max_hops = 2L,
                                  min_neighbours = 2L) {
  stopifnot(inherits(project, "NTSProject"))
  .nts_check(NULL, "igraph_ns", .arg = "igraph")

  g <- graph
  if (is.null(g) || !inherits(g, "igraph") || igraph::vcount(g) == 0L) {
    stop("`graph` must be a non-empty igraph (see build_molecular_network()$igraph).", call. = FALSE)
  }

  ann <- methods::slot(project, "db_annotations")
  ft  <- get_feature_table(project)
  if (!nrow(ann) || !label_col %in% names(ann)) {
    warning("No annotations with column '", label_col, "'; nothing to propagate.",
            call. = FALSE)
    return(project)
  }

  join_key <- intersect(c("name", "feature_id"), intersect(names(ann), names(ft)))
  if (!length(join_key)) {
    warning("Cannot link annotations to features (no shared key column).", call. = FALSE)
    return(project)
  }
  key <- join_key[1L]

  node_names <- igraph::V(g)$name
  ann_map <- stats::setNames(as.character(ann[[label_col]]), as.character(ann[[key]]))
  ann_map <- ann_map[!is.na(ann_map) & nzchar(ann_map)]

  labels <- rep(NA_character_, length(node_names))
  names(labels) <- node_names
  hit <- intersect(node_names, names(ann_map))
  labels[hit] <- ann_map[hit]

  prop_col <- paste0("propagated_", label_col)
  propagated <- rep(NA_character_, length(node_names))
  names(propagated) <- node_names

  weights <- if ("weight" %in% igraph::edge_attr_names(g)) {
    igraph::E(g)$weight
  } else if ("edge_score" %in% igraph::edge_attr_names(g)) {
    igraph::E(g)$edge_score
  } else {
    NULL
  }

  for (hop in seq_len(max_hops)) {
    unknown <- node_names[is.na(labels) & is.na(propagated)]
    if (!length(unknown)) break

    labeled_nodes <- node_names[!is.na(labels) & nzchar(labels)]
    if (!length(labeled_nodes)) break

    adj_list <- igraph::adjacent_vertices(g, labeled_nodes, mode = "all")
    adj_names <- unique(unlist(lapply(adj_list, function(v) v$name), use.names = FALSE))
    cand_unknown <- intersect(unknown, adj_names)
    if (!length(cand_unknown)) break

    changed <- 0L
    for (u in cand_unknown) {
      nb <- igraph::neighbors(g, u)$name
      nb_labels <- labels[nb]
      nb_labels <- nb_labels[!is.na(nb_labels) & nzchar(nb_labels)]
      if (length(nb_labels) < min_neighbours) next

      tbl <- sort(table(nb_labels), decreasing = TRUE)
      winner <- names(tbl)[1L]

      if (length(tbl) > 1L && tbl[1L] == tbl[2L]) {
        tied <- names(tbl)[tbl == tbl[1L]]
        if (!is.null(weights) && is.numeric(weights)) {
          avg_w <- vapply(tied, function(lab) {
            nb_with_lab <- nb[nb_labels == lab]
            eids <- igraph::get.edge.ids(g, c(rbind(u, nb_with_lab)))
            val <- mean(weights[eids[eids > 0]], na.rm = TRUE)
            if (is.finite(val)) val else 0
          }, numeric(1L))
          idx <- which.max(avg_w)
          winner <- if (length(idx) && is.finite(idx) && idx > 0) tied[idx] else tied[1L]
        } else {
          winner <- tied[1L]
        }
      }
      if (length(winner) == 1L && !is.na(winner) && nzchar(winner)) {
        propagated[u] <- winner
        changed <- changed + 1L
      }
    }
    labels[!is.na(propagated) & is.na(labels)] <- propagated[!is.na(propagated) & is.na(labels)]
    if (changed == 0L) break
    message("> Label propagation hop ", hop, ": +", changed, " node(s)")
  }

  # Inferred labels are persisted on @ms1_annotations (joined by `name`), NOT on
  # @db_annotations (see the roxygen note above): db_annotations is the
  # spectral-library layer and must not receive voted/inferred rows.
  n_prop <- sum(!is.na(propagated))
  n_written <- 0L
  if (n_prop > 0L) {
    ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
    if (nrow(ms1) && "name" %in% names(ms1)) {
      if (!prop_col %in% names(ms1)) {
        ms1[[prop_col]] <- NA_character_
      }
      idx <- match(
        trimws(as.character(ms1$name)),
        trimws(as.character(node_names))
      )
      hit <- !is.na(idx)
      if (any(hit)) {
        ms1[[prop_col]][hit] <- as.character(propagated[idx[hit]])
        n_written <- sum(hit)
        project <- set_ms1_annotations(
          project,
          .nts_ms1_sanitize_and_order(ms1)
        )
      }
    }
  }

  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "propagate_annotations",
    list(label_col = label_col, max_hops = max_hops,
         min_neighbours = min_neighbours, n_propagated = n_prop,
         n_written_ms1 = n_written)
  )
  message("> Label propagation finished: ", n_prop,
          " node(s) carry propagated labels (", n_written, " written to ms1)")
  project
}

# ============================================================================
# KMD Homologous Series
# ============================================================================

#' Resolve Kendrick repeating unit (exact / integer nominal)
#' @noRd
.kmd_base <- function(base) {
  b <- trimws(as.character(base)[1L])
  if (is.na(b) || !nzchar(b)) {
    stop("`base` must be a non-empty Hill-style formula string.", call. = FALSE)
  }
  ex <- suppressWarnings(tryCatch(
    calculate_exact_mass(b),
    error = function(e) NA_real_
  ))
  if (!is.finite(ex) || ex <= 0) {
    stop("Could not compute exact mass for base: ", b, call. = FALSE)
  }
  list(exact = ex, nominal = as.integer(round(ex)))
}

#' KMD-based homologous series tagging
#'
#' Computes Kendrick mass (`KM = mz * nominal/exact`) and Kendrick mass defect
#' `KMD = KM - floor(KM)`, rounds features into `kmd_bin = round(KMD / hs_tol)`
#' strips, and assigns a common `homologous_series_id` to bins with at least
#' `n_min` members. Also adds `nominal_mass` (integer, `round(mz)`) and `kmd`
#' for interpretation.
#'
#' @param feature_table `data.frame`, `data.table`, or [NTSProject].
#' @param base Repeating unit: `"CF2"` (default) or `"CH2"` (exact masses per
#'   PRD).
#' @param hs_tol Width of KMD discretisation; bins have width ~`hs_tol`.
#' @param n_min Minimum features in a bin to call it a homologous series.
#'
#' @return Table with `nominal_mass`, `kmd`, `homologous_series_id`, or updated
#'   [NTSProject].
#' @export
filter_homologues <- function(feature_table,
                              base = "CF2",
                              hs_tol = 0.005,
                              n_min = 3) {
  if (inherits(feature_table, "NTSProject")) {
    proj <- feature_table
    ft <- get_feature_table(proj)
    if (!nrow(ft)) {
      stop("feature_table is empty", call. = FALSE)
    }
    out <- .filter_homologues_df(ft, base = base, hs_tol = hs_tol, n_min = n_min)
    proj <- set_feature_table(proj, out)
    methods::slot(proj, "log") <- .nts_log_append(
      methods::slot(proj, "log"),
      "filter_homologues",
      list(base = base, hs_tol = hs_tol, n_min = n_min)
    )
    return(proj)
  }
  .filter_homologues_df(feature_table, base = base, hs_tol = hs_tol, n_min = n_min)
}

#' @noRd
.filter_homologues_df <- function(feature_table,
                                   base = "CF2",
                                   hs_tol = 0.005,
                                   n_min = 3) {
  bu <- .kmd_base(base)
  exact_b <- bu$exact
  nominal_b <- bu$nominal

  if (!data.table::is.data.table(feature_table)) {
    ft <- data.table::as.data.table(feature_table)
    as_dt <- FALSE
  } else {
    ft <- data.table::copy(feature_table)
    as_dt <- TRUE
  }
  stopifnot("mz" %in% names(ft))
  n_min <- max(1L, as.integer(n_min)[1L])
  hs_tol <- as.numeric(hs_tol)[1L]
  if (!is.finite(hs_tol) || hs_tol <= 0) {
    stop("`hs_tol` must be positive and finite.", call. = FALSE)
  }

  ft[, `:=`(
    nominal_mass = as.integer(round(mz)),
    KM = mz * (nominal_b / exact_b)
  )]
  ft[, kmd := KM - floor(KM)]
  ft[, KM := NULL]
  ft[, kmd_bin := round(kmd / hs_tol)]

  ok_bins <- ft[, .N, by = kmd_bin][N >= n_min, kmd_bin]
  ft[, homologous_series_id := NA_integer_]
  if (length(ok_bins)) {
    id_map <- data.table::data.table(
      kmd_bin = ok_bins,
      homologous_series_id = seq_along(ok_bins)
    )
    ft[id_map, homologous_series_id := i.homologous_series_id, on = "kmd_bin"]
  }
  ft[, kmd_bin := NULL]

  if (!as_dt) {
    return(as.data.frame(ft, stringsAsFactors = FALSE))
  }
  ft[]
}

# ============================================================================
# KMD-anchor homologue targeting and formula inference (data.table + base R)
# ============================================================================

#' @noRd
.parse_formula_element_counts <- function(formula) {
  f <- gsub("\\s+", "", as.character(formula)[1L])
  if (is.na(f) || !nzchar(f)) {
    return(structure(integer(0), names = character(0)))
  }
  m <- gregexpr("([A-Z][a-z]*)(\\d*)", f, perl = TRUE)
  mt <- m[[1L]]
  if (any(mt < 0L)) {
    return(structure(integer(0), names = character(0)))
  }
  tok <- regmatches(f, m)[[1L]]
  el <- character(0)
  ct <- integer(0)
  for (t in tok) {
    ml <- regexec("^([A-Z][a-z]*)(\\d*)$", t, perl = TRUE)
    parts <- regmatches(t, ml)[[1L]]
    if (length(parts) < 2L) {
      return(structure(integer(0), names = character(0)))
    }
    elem <- parts[[2L]]
    ns <- parts[[3L]]
    cnt <- if (nzchar(ns)) suppressWarnings(as.integer(ns)) else 1L
    if (is.na(cnt) || cnt < 1L) {
      return(structure(integer(0), names = character(0)))
    }
    idx <- match(elem, el)
    if (is.na(idx)) {
      el <- c(el, elem)
      ct <- c(ct, cnt)
    } else {
      ct[idx] <- ct[idx] + cnt
    }
  }
  structure(ct, names = el)
}

#' @noRd
.hill_element_order <- function(elems) {
  u <- sort(unique(elems))
  has_c <- "C" %in% u
  has_h <- "H" %in% u
  if (has_c) {
    rest <- sort(setdiff(u, c("C", "H")))
    if (has_h) {
      c("C", "H", rest)
    } else {
      c("C", rest)
    }
  } else {
    u
  }
}

#' @noRd
.format_sorted_formula <- function(named_counts) {
  if (!length(named_counts)) {
    return(NA_character_)
  }
  nn <- as.integer(named_counts)
  names(nn) <- names(named_counts)
  nn <- nn[!is.na(nn) & nn > 0L]
  if (!length(nn)) {
    return(NA_character_)
  }
  o <- .hill_element_order(names(nn))
  parts <- vapply(o, function(e) {
    n <- as.integer(nn[e])
    if (is.na(n) || n <= 0L) {
      return(NA_character_)
    }
    if (n == 1L) e else paste0(e, n)
  }, character(1L))
  if (any(is.na(parts))) {
    return(NA_character_)
  }
  paste(parts, collapse = "")
}

#' @title Adjust a molecular formula by adding repeat units
#'
#' @description
#' Parses `formula` and `base_formula` with Base R (`gregexpr` + `regmatches`),
#' adds `(count in base) * steps` to each element, and returns a Hill-system
#' string (C first, then H if present, then other elements alphabetically;
#' stoichiometry 1 is implicit). Missing counts
#' after an element symbol are treated as `1`. Returns `NA_character_` if
#' inputs are invalid or any resulting stoichiometry is non-positive.
#'
#' @param formula Character molecular formula.
#' @param steps Integer repeat-unit steps (may be negative).
#' @param base_formula Character formula of one repeat unit (e.g. `"CF2"`,
#'   `"CF2O"`).
#'
#' @return Adjusted formula string, or `NA_character_`.
#' @aliases adjust_homologue_formula_v
#' @export
adjust_homologue_formula <- function(formula, steps, base_formula) {
  formula <- as.character(formula)[1L]
  steps <- suppressWarnings(as.integer(round(as.numeric(steps)[1L])))
  base_formula <- as.character(base_formula)[1L]
  if (is.na(formula) || !nzchar(formula) || is.na(steps) || is.na(base_formula) || !nzchar(base_formula)) {
    return(NA_character_)
  }
  fcnt <- .parse_formula_element_counts(formula)
  bcnt <- .parse_formula_element_counts(base_formula)
  if (!length(names(fcnt)) || !length(names(bcnt))) {
    return(NA_character_)
  }
  elems <- unique(c(names(fcnt), names(bcnt)))
  outv <- integer(length(elems))
  names(outv) <- elems
  for (e in elems) {
    fe <- if (e %in% names(fcnt)) as.integer(fcnt[e]) else 0L
    be <- if (e %in% names(bcnt)) as.integer(bcnt[e]) else 0L
    outv[e] <- fe + be * steps
  }
  if (any(outv < 0L | is.na(outv))) {
    return(NA_character_)
  }
  outv <- outv[outv > 0L]
  if (!length(outv)) {
    return(NA_character_)
  }
  .format_sorted_formula(outv)
}

#' @export
#' @rdname adjust_homologue_formula
adjust_homologue_formula_v <- Vectorize(
  adjust_homologue_formula,
  vectorize.args = c("formula", "steps", "base_formula")
)

#' @title Adjust formula by CF2 repeat units
#'
#' @description
#' Thin wrapper around [adjust_homologue_formula] with `base_formula = "CF2"`.
#' Output element order is alphabetical (e.g. `C9F19`), not necessarily
#' the input order.
#'
#' @param formula Character molecular formula (single string per call).
#' @param n_cf2 Integer number of `CF2` units to add (may be negative).
#'
#' @return Adjusted formula string, or `NA_character_`.
#' @export
adjust_pfas_formula <- function(formula, n_cf2) {
  adjust_homologue_formula(formula, n_cf2, "CF2")
}

#' @export
#' @rdname adjust_pfas_formula
adjust_pfas_formula_v <- Vectorize(adjust_pfas_formula, vectorize.args = c("formula", "n_cf2"))

#' @noRd
.homologue_series_base_formula <- function(kid) {
  kid <- as.character(kid)[1L]
  if (is.na(kid) || !nzchar(kid)) {
    return(NA_character_)
  }
  m <- regexec("^(.+)_R[0-9]{2}_Series_[0-9]+$", kid, perl = TRUE)
  p <- regmatches(kid, m)[[1L]]
  if (length(p) >= 2L && nzchar(p[[2L]])) {
    p[[2L]]
  } else {
    NA_character_
  }
}

#' @noRd
.homologue_resolve_base_and_denom <- function(kid_part) {
  b <- .homologue_series_base_formula(kid_part)
  if (is.na(b) || !nzchar(b)) {
    b <- "CF2"
  }
  w <- suppressWarnings(tryCatch(
    as.integer(round(.nts_parse_exact_mass(b))),
    error = function(e) NA_integer_
  ))
  if (is.na(w) || w < 1L) {
    w <- as.integer(round(calculate_exact_mass("CF2")))
  }
  list(base = b, denom = w)
}

#' Default repeat-unit tokens for [infer_homologue_structures] algebraic names.
#' @noRd
.nts_default_homologue_transformation_units <- function() {
  c("CF2", "CH2", "O", "H2")
}

#' Signed integer CF2 step count on one edge label (e.g. `+2CF2` -> 2, `-CF2` -> -1).
#' Non-CF2 / unknown labels return 0.
#' @noRd
.nts_homologue_edge_cf2_step_count <- function(transformation) {
  tr <- trimws(as.character(transformation)[1L])
  if (is.na(tr) || !nzchar(tr)) {
    return(0L)
  }
  tru <- toupper(tr)
  if (tru == "UNKNOWN") {
    return(0L)
  }
  if (!grepl("CF2", tru, fixed = TRUE)) {
    return(0L)
  }
  op_edge <- NA_character_
  n_mult <- NA_integer_
  if (grepl("^[+-]\\s*\\d", tru, perl = TRUE)) {
    m1 <- regexec("^([+-])\\s*(\\d+)\\s*CF2\\s*$", tru, perl = TRUE)
    r1 <- regmatches(tru, m1)[[1L]]
    if (length(r1) >= 3L) {
      op_edge <- r1[[2L]]
      n_mult <- suppressWarnings(as.integer(r1[[3L]]))
    }
  }
  if (is.na(n_mult) || n_mult < 1L || is.na(op_edge) || !nzchar(op_edge)) {
    m2 <- regexec("^([+-])\\s*CF2\\s*$", tru, perl = TRUE)
    r2 <- regmatches(tru, m2)[[1L]]
    if (length(r2) >= 2L) {
      op_edge <- r2[[2L]]
      n_mult <- 1L
    }
  }
  if (is.na(n_mult) || n_mult < 1L || is.na(op_edge) || !nzchar(op_edge)) {
    return(0L)
  }
  if (op_edge == "+") {
    as.integer(n_mult)
  } else if (op_edge == "-") {
    -as.integer(n_mult)
  } else {
    0L
  }
}

#' Algebraic homologue name: stacks `PFOA + 1 CF2` + `+1CF2` -> `PFOA + 2 CF2`
#' when the parent matches the counted pattern; otherwise appends
#' `parent " +/- " n " " unit` or `paste0(parent, " ", transformation)` so the
#' function **never** returns `NA` for a non-empty parent (except both inputs empty).
#' @noRd
.nts_algebraic_transform <- function(parent_name,
                                       transformation,
                                       transformation_units = NULL) {
  parent_name <- trimws(as.character(parent_name)[1L])
  transformation <- trimws(as.character(transformation)[1L])

  .bad_lbl <- function(x) {
    x <- trimws(as.character(x)[1L])
    is.na(x) | !nzchar(x) | toupper(x) %in% c("NA", "NAN", "NULL")
  }

  .fb <- function() {
    pn <- parent_name
    tr <- transformation
    if (.bad_lbl(pn) && .bad_lbl(tr)) {
      return(NA_character_)
    }
    if (.bad_lbl(pn)) {
      return(tr)
    }
    if (.bad_lbl(tr)) {
      return(pn)
    }
    paste0(pn, " ", tr)
  }

  if (.bad_lbl(transformation)) {
    if (!.bad_lbl(parent_name)) {
      return(parent_name)
    }
    return(NA_character_)
  }
  if (toupper(transformation) == "UNKNOWN") {
    if (!.bad_lbl(parent_name)) {
      return(parent_name)
    }
    return(NA_character_)
  }
  if (.bad_lbl(parent_name)) {
    return(transformation)
  }

  units_allowed <- transformation_units
  if (is.null(units_allowed)) {
    units_allowed <- .nts_default_homologue_transformation_units()
  } else {
    units_allowed <- unique(as.character(units_allowed))
    units_allowed <- units_allowed[!is.na(units_allowed) & nzchar(units_allowed)]
    if (!length(units_allowed)) {
      units_allowed <- .nts_default_homologue_transformation_units()
    }
  }

  op_edge <- NA_character_
  unit_raw <- NA_character_
  n_mult <- NA_integer_
  if (grepl("^[+-]\\s*\\d", transformation, perl = TRUE)) {
    m1 <- regexec("^([+-])\\s*(\\d+)\\s*([A-Za-z0-9]+)\\s*$", transformation, perl = TRUE)
    r1 <- regmatches(transformation, m1)[[1L]]
    if (length(r1) >= 4L) {
      op_edge <- r1[[2L]]
      n_mult <- suppressWarnings(as.integer(r1[[3L]]))
      unit_raw <- r1[[4L]]
    }
  }
  if (is.na(n_mult) || n_mult < 1L || is.na(op_edge) || !nzchar(op_edge)) {
    m2 <- regexec("^([+-])\\s*([A-Za-z0-9]+)\\s*$", transformation, perl = TRUE)
    r2 <- regmatches(transformation, m2)[[1L]]
    if (length(r2) >= 3L) {
      op_edge <- r2[[2L]]
      unit_raw <- r2[[3L]]
      n_mult <- 1L
    }
  }
  if (is.na(n_mult) || n_mult < 1L || is.na(op_edge) || !nzchar(op_edge) ||
    is.na(unit_raw) || !nzchar(unit_raw)) {
    return(.fb())
  }

  unit_canon <- unit_raw
  if (!unit_raw %in% units_allowed) {
    u2 <- toupper(unit_raw)
    hit <- units_allowed[toupper(units_allowed) == u2]
    if (length(hit)) {
      unit_canon <- hit[[1L]]
    }
  }

  step_edge <- if (op_edge == "+") {
    as.integer(n_mult)
  } else {
    -as.integer(n_mult)
  }
  sig_display <- if (op_edge == "+") "+" else "-"

  cnt <- 0L
  base <- parent_name
  stacked_ok <- FALSE
  pat <- paste0(
    "^(.+?)\\s+([+-])\\s*(\\d+)\\s+(",
    paste(units_allowed, collapse = "|"),
    ")\\s*$"
  )
  mm <- regexec(pat, base, perl = TRUE)
  parts <- regmatches(base, mm)[[1L]]
  if (length(parts) >= 5L) {
    u0 <- parts[[5L]]
    u0c <- u0
    if (!u0 %in% units_allowed) {
      u2 <- toupper(u0)
      hit0 <- units_allowed[toupper(units_allowed) == u2]
      if (length(hit0)) {
        u0c <- hit0[[1L]]
      }
    }
    if (u0c == unit_canon) {
      op0 <- parts[[3L]]
      n0 <- suppressWarnings(as.integer(parts[[4L]]))
      if (!is.na(n0)) {
        cnt <- if (op0 == "+") n0 else -as.integer(n0)
        base <- trimws(parts[[2L]])
        if (nzchar(base)) {
          stacked_ok <- TRUE
        }
      }
    }
  }

  if (stacked_ok) {
    cnt2 <- as.integer(cnt + step_edge)
    if (cnt2 == 0L) {
      return(trimws(base))
    }
    if (cnt2 > 0L) {
      return(sprintf("%s + %d %s", trimws(base), cnt2, unit_canon))
    }
    return(sprintf("%s - %d %s", trimws(base), as.integer(abs(cnt2)), unit_canon))
  }

  paste0(parent_name, " ", sig_display, " ", as.integer(n_mult), " ", unit_canon)
}

#' @noRd
.nts_homologue_edge_is_isf <- function(transformation) {
  tr <- toupper(trimws(as.character(transformation)))
  is.na(tr) | !nzchar(tr) | tr == "ISF_LINK" | grepl("^ISF_", tr, perl = TRUE) |
    grepl("_ISF$", tr, perl = TRUE)
}

#' @noRd
.infer_homologue_structures_network <- function(network_data,
                                                 seeds,
                                                 max_steps = 3L,
                                                 transformation_units = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  if (!is.list(network_data) || !all(c("nodes", "edges") %in% names(network_data))) {
    stop("network_data must be a list with elements `nodes` and `edges`.", call. = FALSE)
  }
  max_steps <- max(0L, as.integer(max_steps)[1L])
  seeds <- unique(as.character(seeds))
  seeds <- seeds[!is.na(seeds) & nzchar(seeds)]
  if (!length(seeds)) {
    stop("`seeds` must contain at least one non-empty feature id.", call. = FALSE)
  }

  N <- data.table::as.data.table(network_data$nodes)
  E <- data.table::as.data.table(network_data$edges)
  Eout <- if (nrow(E)) data.table::copy(E) else E

  if (!nrow(N)) {
    return(list(
      nodes = N,
      edges = Eout,
      inference = data.table::data.table(),
      stats = list(n_seeds = 0L, n_steps = 0L, n_inferred = 0L)
    ))
  }
  if (!"id" %in% names(N) && "name" %in% names(N)) {
    N[, id := as.character(name)]
  }
  if (!"id" %in% names(N)) {
    stop("nodes must contain column `id` (or `name` to alias `id`).", call. = FALSE)
  }
  N[, id := as.character(id)]
  if (!"formula" %in% names(N)) {
    N[, formula := NA_character_]
  } else {
    N[, formula := as.character(formula)]
  }
  if (!".nts_seed_label" %in% names(N)) {
    N[, .nts_seed_label := NA_character_]
  } else {
    N[, .nts_seed_label := trimws(as.character(.nts_seed_label))]
  }

  if (!nrow(E) || !all(c("from", "to") %in% names(E))) {
    N2 <- data.table::copy(N)
    if (".nts_seed_label" %in% names(N2)) {
      N2[, .nts_seed_label := NULL]
    }
    return(list(
      nodes = N2,
      edges = Eout,
      inference = data.table::data.table(),
      stats = list(n_seeds = 0L, n_steps = 0L, n_inferred = 0L, n_waves = 0L, note = "no edges")
    ))
  }
  E[, from := as.character(from)]
  E[, to := as.character(to)]
  if (!"transformation" %in% names(E)) {
    E[, transformation := NA_character_]
  }
  if (!"smiles_op_payload" %in% names(E)) {
    E[, smiles_op_payload := NA_character_]
  }

  sh <- N[id %chin% seeds]
  if (nrow(sh)) {
    sh[, `:=`(
      seed_root = id,
      evidence_path = paste0("SEED:", id),
      depth = 0L
    )]
    sh[, sx := trimws(as.character(.nts_seed_label))]
    sh[, ok_sx := !is.na(sx) & nzchar(sx) & !toupper(sx) %in% c("NA", "NAN", "NULL")]
    if ("library_Compound_name" %in% names(sh)) {
      sh[, alt_lib := trimws(as.character(library_Compound_name))]
    } else {
      sh[, alt_lib := NA_character_]
    }
    sh[, ok_alt_lib := !is.na(alt_lib) & nzchar(alt_lib) & !toupper(alt_lib) %in% c("NA", "NAN", "NULL")]
    if ("Compound.name" %in% names(sh)) {
      sh[, alt_cmp := trimws(as.character(Compound.name))]
    } else {
      sh[, alt_cmp := NA_character_]
    }
    sh[, ok_alt_cmp := !is.na(alt_cmp) & nzchar(alt_cmp) & !toupper(alt_cmp) %in% c("NA", "NAN", "NULL")]
    sh[, label := data.table::fifelse(
      ok_sx,
      sx,
      data.table::fifelse(
        ok_alt_lib,
        alt_lib,
        data.table::fifelse(ok_alt_cmp, alt_cmp, as.character(id))
      )
    )]
    sh[, c("sx", "ok_sx", "alt_lib", "alt_cmp", "ok_alt_lib", "ok_alt_cmp") := NULL]
    sh[, cum_cf2 := 0L]
  }
  seed_hit <- sh
  if (!nrow(seed_hit)) {
    warning(
      "No `seeds` matched any `nodes$id` (or `name` when aliased to `id`); returning empty inference.",
      call. = FALSE
    )
    N2 <- data.table::copy(N)
    if (".nts_seed_label" %in% names(N2)) {
      N2[, .nts_seed_label := NULL]
    }
    return(list(
      nodes = N2,
      edges = Eout,
      inference = data.table::data.table(),
      stats = list(
        n_seeds = 0L,
        n_steps = max_steps,
        n_inferred = 0L,
        n_waves = 0L,
        note = "no_seed_match"
      )
    ))
  }

  # BFS: `state` always includes seed rows so long bridges from seeds stay
  # eligible every wave; `annotated_to` blocks re-targeting nodes already
  # resolved as `to` (seeds + newly inferred).
  state <- seed_hit[!is.na(label) & nzchar(label), .(id, label, seed_root, evidence_path, depth, cum_cf2)]
  annotated_to <- unique(state$id)
  history_nodes <- list()
  step <- 0L

  while (TRUE) {
    if (max_steps > 0L && step >= max_steps) {
      break
    }
    step <- step + 1L

    state_r <- data.table::copy(state)
    data.table::setnames(state_r, "id", "sid", skip_absent = TRUE)
    # Rename the parent's `label` *before* the join: homologue edges from
    # get_network_data() / build_homologue_network() always carry a `label`
    # column (the transformation tag). After `E[state_r, ...]` the incoming
    # (i) `label` would be shadowed to `i.label`, so `parent_label := label`
    # would read the edge tag instead of the parent compound name (P0-E1).
    data.table::setnames(state_r, "label", "parent_label", skip_absent = TRUE)
    M <- E[state_r, on = c(from = "sid"), nomatch = NA, allow.cartesian = TRUE]
    M <- M[!(to %chin% annotated_to)]
    if (!nrow(M)) {
      break
    }

    tr <- trimws(as.character(M$transformation))
    tru <- toupper(tr)
    unk <- is.na(tr) | !nzchar(tr) | tru == "UNKNOWN"
    isf <- .nts_homologue_edge_is_isf(M$transformation)
    M[, inferred_name := NA_character_]
    ok <- which(!unk & !isf)
    if (length(ok)) {
      inn <- rep(NA_character_, nrow(M))
      inn[ok] <- as.character(mapply(
        function(p, t) {
          .nts_algebraic_transform(
            p,
            t,
            transformation_units = transformation_units
          )
        },
        as.character(M[["parent_label"]][ok]),
        as.character(M[["transformation"]][ok]),
        USE.NAMES = FALSE
      ))
      M[, inferred_name := inn]
    }
    dcf2 <- vapply(seq_len(nrow(M)), function(i) {
      .nts_homologue_edge_cf2_step_count(M[["transformation"]][i])
    }, integer(1L), USE.NAMES = FALSE)
    M[, dcf2 := as.integer(dcf2)]
    if (!"cum_cf2" %in% names(M)) {
      M[, cum_cf2 := 0L]
    }
    M[, cum_cf2 := as.integer(cum_cf2) + as.integer(dcf2)]
    M[, dcf2 := NULL]
    M[, isf_pseudo := isf | unk]
    M[isf_pseudo == TRUE, inferred_name := NA_character_]
    M[, evidence_path := paste0(
      evidence_path,
      " -> [",
      ifelse(is.na(transformation), "?", as.character(transformation)),
      "] -> ",
      as.character(to)
    )]

    rec <- M[, .(
      name = as.character(to),
      inferred_name = as.character(inferred_name),
      evidence_path = as.character(evidence_path),
      seed_root = as.character(seed_root),
      step = as.integer(step),
      transformation = as.character(transformation),
      isf_pseudo = as.logical(isf_pseudo),
      cum_cf2 = as.integer(cum_cf2)
    )]
    history_nodes[[length(history_nodes) + 1L]] <- rec

    okm <- M[
      isf_pseudo == FALSE &
        !is.na(inferred_name) &
        nzchar(trimws(as.character(inferred_name))),
    ]
    if (!nrow(okm)) {
      break
    }

    new_state <- okm[, .(
      id = as.character(to),
      label = as.character(inferred_name),
      seed_root = as.character(seed_root),
      evidence_path = as.character(evidence_path),
      depth = as.integer(step),
      cum_cf2 = as.integer(cum_cf2)
    )]
    state <- unique(
      data.table::rbindlist(list(state, new_state), use.names = TRUE, fill = TRUE),
      by = c("id", "label")
    )
    new_tos <- unique(as.character(okm$to))
    annotated_to <- unique(c(annotated_to, new_tos))
  }

  parts <- history_nodes[vapply(history_nodes, function(z) {
    data.table::is.data.table(z) && nrow(z) > 0L
  }, logical(1L))]
  if (!length(parts)) {
    H <- data.table::data.table()
  } else {
    H <- data.table::rbindlist(parts, fill = TRUE)
  }
  if (!nrow(H)) {
    N2 <- data.table::copy(N)
    if (".nts_seed_label" %in% names(N2)) {
      N2[, .nts_seed_label := NULL]
    }
    return(list(
      nodes = N2,
      edges = Eout,
      inference = data.table::data.table(),
      stats = list(
        n_seeds = nrow(seed_hit),
        n_steps = max_steps,
        n_inferred = 0L,
        n_waves = 0L
      )
    ))
  }

  agg <- H[, .(
    evidence_path = paste(unique(evidence_path), collapse = " | "),
    n_paths = .N,
    transformation = paste(unique(transformation), collapse = " | "),
    isf_pseudo = any(isf_pseudo, na.rm = TRUE),
    seed_root = data.table::first(as.character(seed_root)),
    cum_cf2 = as.integer(cum_cf2[which.max(step)])
  ), by = .(name, inferred_name)]
  agg[, has_nm := !is.na(inferred_name) & nzchar(trimws(as.character(inferred_name)))]
  data.table::setorder(agg, name, -n_paths, -has_nm)
  agg[, has_nm := NULL]
  win <- agg[, .SD[1L], by = name]

  N2 <- data.table::copy(N)
  win_small <- win[, .(
    name,
    inferred_name,
    seed_root,
    cum_cf2,
    inference_evidence_path = evidence_path,
    homologue_inference_n_paths = n_paths,
    homologue_inference_transformations = transformation,
    homologue_inference_isf_pseudo = isf_pseudo
  )]
  fseed <- unique(
    N[, .(seed_root = as.character(id), anchor_seed_formula = as.character(formula))],
    by = "seed_root"
  )
  win_small <- fseed[win_small, on = "seed_root"]
  win_small[, inferred_formula := {
    af <- anchor_seed_formula
    k <- cum_cf2
    out <- rep(NA_character_, length(af))
    for (ii in seq_along(out)) {
      a <- af[ii]
      kk <- k[ii]
      if (is.na(kk) || !is.finite(kk)) {
        next
      }
      if (is.na(a) || !nzchar(trimws(as.character(a)))) {
        next
      }
      out[ii] <- tryCatch(
        as.character(adjust_homologue_formula(as.character(a)[1L], as.integer(kk)[1L], "CF2"))[1L],
        error = function(e) NA_character_
      )
    }
    out
  }]
  win_small[, inferred_smiles := NA_character_]
  win_small[, c("anchor_seed_formula", "seed_root", "cum_cf2") := NULL]
  if (!"name" %in% names(N2)) {
    N2[, name := id]
  }
  N2 <- win_small[N2, on = "name"]
  if (".nts_seed_label" %in% names(N2)) {
    N2[, .nts_seed_label := NULL]
  }

  list(
    nodes = N2,
    edges = Eout,
    inference = win,
    raw_inference = H,
    stats = list(
      n_seeds = nrow(seed_hit),
      n_steps = max_steps,
      n_inferred = nrow(win),
      n_raw_rows = nrow(H),
      n_waves = as.integer(length(parts))
    )
  )
}

#' Assemble enhanced `network_data` (nodes + edges + passthrough slots) for plotting.
#'
#' @keywords internal
#' @noRd
.nts_homologue_enrich_network_data <- function(project, network_data, seeds_ch, inf_out) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  extra <- setdiff(names(network_data), c("nodes", "edges"))
  nds <- data.table::as.data.table(inf_out$nodes)
  eds <- data.table::as.data.table(inf_out$edges)

  infer_cols <- c(
    "inferred_name", "inferred_formula", "inferred_smiles",
    "inference_evidence_path", "homologue_inference_n_paths",
    "homologue_inference_transformations", "homologue_inference_isf_pseudo"
  )
  for (cn in infer_cols) {
    if (!cn %in% names(nds)) {
      if (cn == "homologue_inference_n_paths") {
        nds[, (cn) := NA_integer_]
      } else if (cn == "homologue_inference_isf_pseudo") {
        nds[, (cn) := NA]
      } else {
        nds[, (cn) := NA_character_]
      }
    }
  }

  nds[, homologue_steps := as.character(homologue_inference_transformations)]
  nds[, evidence_path := as.character(inference_evidence_path)]

  id_chr <- as.character(nds$id)
  nds[, evidence_source := NA_character_]
  seeds_ch <- unique(as.character(seeds_ch))
  seeds_ch <- seeds_ch[!is.na(seeds_ch) & nzchar(seeds_ch)]
  if (length(seeds_ch)) {
    nds[id_chr %chin% seeds_ch, evidence_source := "Seed"]
  }
  inf_mask <- !is.na(nds[["inferred_name"]]) &
    nzchar(trimws(as.character(nds[["inferred_name"]])))
  nds[inf_mask & (is.na(evidence_source) | !nzchar(as.character(evidence_source))), evidence_source := "Inferred"]

  nds[, homologue_label := NA_character_]
  nds[evidence_source == "Inferred" & inf_mask, homologue_label := trimws(as.character(inferred_name))]

  libv <- if ("library_Compound_name" %in% names(nds)) {
    trimws(as.character(nds[["library_Compound_name"]]))
  } else {
    rep(NA_character_, nrow(nds))
  }
  cmpv <- if ("Compound.name" %in% names(nds)) {
    trimws(as.character(nds[["Compound.name"]]))
  } else {
    rep(NA_character_, nrow(nds))
  }
  seed_idx <- !is.na(nds[["evidence_source"]]) & nds[["evidence_source"]] == "Seed"
  if (any(seed_idx, na.rm = TRUE)) {
    hl_seed <- ifelse(
      !is.na(libv) & nzchar(libv) & !toupper(libv) %in% c("NA", "NAN", "NULL"),
      libv,
      ifelse(
        !is.na(cmpv) & nzchar(cmpv) & !toupper(cmpv) %in% c("NA", "NAN", "NULL"),
        cmpv,
        id_chr
      )
    )
    ii <- which(seed_idx)
    nds[ii, homologue_label := hl_seed[ii]]
  }

  annf <- tryCatch(get_ms1_annotations(project), error = function(e) NULL)
  if (is.data.frame(annf) && nrow(annf) && "name" %in% names(annf)) {
    nm <- as.character(annf$name)
    ix <- match(id_chr, nm)
    sm_cols <- intersect(c("library_SMILES", "SMILES", "smiles"), names(annf))
    for (sc in sm_cols) {
      src <- as.character(annf[[sc]])[ix]
      miss <- is.na(nds[["inferred_smiles"]]) | !nzchar(trimws(as.character(nds[["inferred_smiles"]])))
      w <- which(miss & !is.na(ix) & !is.na(src) & nzchar(trimws(as.character(src))))
      if (length(w)) {
        data.table::set(nds, i = w, j = "inferred_smiles", value = src[w])
      }
    }
  }

  if (".nts_seed_label" %in% names(nds)) {
    nds[, .nts_seed_label := NULL]
  }

  out <- list(
    nodes = as.data.frame(nds, stringsAsFactors = FALSE),
    edges = as.data.frame(eds, stringsAsFactors = FALSE)
  )
  for (nm in extra) {
    out[[nm]] <- network_data[[nm]]
  }
  out
}

#' @keywords internal
#' @noRd
.nts_homologue_strip_scaffold_node_cols <- function(net_data) {
  .nts_network_data_strip_intermediate_cols(net_data, strip_edges = FALSE)
}

#' @keywords internal
#' @noRd
.infer_homologue_structures_project <- function(project) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject object.", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  ann_raw <- get_ms1_annotations(project)
  if (!nrow(ann_raw) || !"name" %in% names(ann_raw)) {
    stop("ms1_annotations is empty or lacks `name`", call. = FALSE)
  }

  ft <- get_feature_table(project)
  if (!nrow(ft) || !"name" %in% names(ft) || !"mz" %in% names(ft)) {
    stop("feature_table must be non-empty with `name` and `mz`", call. = FALSE)
  }

  ann <- data.table::as.data.table(ann_raw)
  db0 <- .nts_db_best_per_name(methods::slot(project, "db_annotations"))
  if (is.data.frame(db0) && nrow(db0)) {
    ddt <- data.table::as.data.table(db0)
    ann <- ddt[ann, on = "name"]
  }
  had_mz <- "mz" %in% names(ann_raw)
  ftz <- unique(data.table::as.data.table(ft[, c("name", "mz"), drop = FALSE]), by = "name")
  data.table::setnames(ftz, c("name", "ft_mz"))
  ann <- ftz[ann, on = "name"]
  if (had_mz) {
    ann[!is.finite(mz), mz := ft_mz]
    ann[, ft_mz := NULL]
  } else {
    data.table::setnames(ann, "ft_mz", "mz")
  }

  if (!"formula" %in% names(ann)) {
    ann[, formula := NA_character_]
  }
  if (!"score_total" %in% names(ann)) {
    ann[, score_total := NA_real_]
  }
  if (!"Compound.name" %in% names(ann)) {
    ann[, Compound.name := NA_character_]
  }
  if (!"kmd_series_id" %in% names(ann)) {
    ann[, kmd_series_id := NA_character_]
  }

  ann[, row_id := .I]
  ann[, kid_raw := as.character(kmd_series_id)]

  bycols <- setdiff(names(ann), "kmd_series_id")
  long <- ann[, {
    kid0 <- as.character(kmd_series_id)
    if (is.na(kid0) || !nzchar(trimws(kid0))) {
      .(kid_part = NA_character_, series_idx = NA_integer_)
    } else {
      parts <- trimws(strsplit(kid0, ";\\s*", perl = TRUE)[[1L]])
      parts <- parts[nzchar(parts) & !(toupper(parts) %in% c("NA", "NAN"))]
      if (!length(parts)) {
        .(kid_part = NA_character_, series_idx = NA_integer_)
      } else {
        .(kid_part = parts, series_idx = seq_along(parts))
      }
    }
  }, by = bycols]

  long[, isolated := is.na(kid_part) | !nzchar(as.character(kid_part))]
  long[, has_form := !is.na(formula) & nzchar(as.character(formula))]
  long[, st := suppressWarnings(as.numeric(score_total))]
  long[is.na(st), st := -Inf]

  bd <- lapply(long$kid_part, function(k) {
    if (is.na(k) || !nzchar(k)) {
      list(base = "CF2", denom = as.integer(round(calculate_exact_mass("CF2"))))
    } else {
      .homologue_resolve_base_and_denom(k)
    }
  })
  long[, base_for_step := vapply(bd, function(z) z$base, character(1L))]
  long[, step_denom := vapply(bd, function(z) z$denom, integer(1L))]

  cand <- long[!isolated & (has_form == TRUE)]
  anchors <- if (nrow(cand)) {
    pick <- cand[, .SD[which.max(st)], by = kid_part]
    pick[, .(
      kid_part = kid_part,
      anchor_feature = as.character(name),
      anchor_mz = as.numeric(mz),
      anchor_formula = as.character(formula),
      anchor_name = {
        cn <- as.character(Compound.name)
        ifelse(is.na(cn) | !nzchar(cn), as.character(name), cn)
      }
    )]
  } else {
    data.table::data.table(
      kid_part = character(0),
      anchor_feature = character(0),
      anchor_mz = numeric(0),
      anchor_formula = character(0),
      anchor_name = character(0)
    )
  }

  if (nrow(anchors)) {
    long <- merge(long, anchors, by = "kid_part", all.x = TRUE, sort = FALSE)
  } else {
    long[, anchor_feature := NA_character_]
    long[, anchor_mz := NA_real_]
    long[, anchor_formula := NA_character_]
    long[, anchor_name := NA_character_]
  }

  long[, series_steps := suppressWarnings(round((as.numeric(mz) - as.numeric(anchor_mz)) / as.numeric(step_denom)))]
  long[is.na(series_steps), series_steps := NA_integer_]

  long[, is_anchor := !is.na(anchor_feature) & nzchar(as.character(anchor_feature)) &
    as.character(name) == as.character(anchor_feature)]

  long[, inferred_formula := NA_character_]
  long[is_anchor == TRUE, inferred_formula := as.character(formula)]
  infer_mask <- !long[["isolated"]] &
    !is.na(long[["anchor_feature"]]) &
    nzchar(as.character(long[["anchor_feature"]])) &
    long[["is_anchor"]] %in% FALSE
  infer_idx <- which(infer_mask)
  if (length(infer_idx)) {
    vinf <- rep(NA_character_, length(infer_idx))
    for (jj in seq_along(infer_idx)) {
      i <- infer_idx[jj]
      af <- long[["anchor_formula"]][i]
      st <- long[["series_steps"]][i]
      bf <- long[["base_for_step"]][i]
      if (is.na(af) || !nzchar(trimws(as.character(af)))) {
        next
      }
      if (is.na(bf) || !nzchar(trimws(as.character(bf)))) {
        next
      }
      st1 <- suppressWarnings(as.integer(round(as.numeric(st)[1L])))
      if (!is.finite(st1)) {
        next
      }
      vinf[jj] <- tryCatch(
        {
          out1 <- adjust_homologue_formula(
            as.character(af)[1L],
            st1,
            as.character(bf)[1L]
          )
          as.character(out1)[1L]
        },
        error = function(e) NA_character_
      )
    }
    long[infer_idx, inferred_formula := vinf]
  }

  long[, annotation_level := data.table::fcase(
    isolated, "Isolated Unknown",
    is.na(anchor_feature) | !nzchar(as.character(anchor_feature)), "Unknown Series (No Anchor)",
    is_anchor, "Database Confirmed (Anchor)",
    default = sprintf("Inferred from %s (%+d %s)", anchor_name, series_steps, base_for_step)
  )]

  long[, has_inf := !is.na(inferred_formula) & nzchar(inferred_formula)]
  data.table::setorder(long, row_id, -has_inf, series_idx)
  collapsed <- unique(long, by = "row_id")
  collapsed[, kmd_series_id := kid_raw]
  collapsed[, c(
    "kid_raw", "kid_part", "series_idx", "base_for_step", "step_denom",
    "has_inf", "isolated", "has_form", "st",
    "anchor_feature", "anchor_mz", "anchor_formula", "anchor_name",
    "is_anchor", "series_steps", "row_id"
  ) := NULL]

  out_df <- as.data.frame(collapsed, stringsAsFactors = FALSE)
  db_only <- c(
    "ref_id", "score_total", "score_mz", "score_rt", "score_ms2", "score_ndp",
    "score_nl", "score_md", "score_diagnostic", "ndp_raw", "matched_nl", "matched_md",
    "matched_diagnostic", "best_ms2_index", "n_ms2_subset", "mz_error_ppm", "rt_error_s",
    "exp_mz", "lib_mz", "exp_rt", "lib_rt", "n_matched_peaks", "Compound.name",
    "formula", "SMILES", "smiles"
  )
  ms1 <- get_ms1_annotations(project)
  carry <- setdiff(names(out_df), db_only)
  carry <- carry[carry != "name"]
  idx <- match(as.character(ms1$name), as.character(out_df$name))
  for (cl in carry) {
    if (!cl %in% names(out_df)) next
    ms1[[cl]] <- out_df[[cl]][idx]
  }
  # Confidence wiring (G): a KMD-anchor *derived* row (annotation_level
  # "Inferred from <anchor>") carries no `annotation_method`, so
  # assign_confidence_levels() would never place it in a tier. Tag those rows
  # "Homologue" (recognised as L3c by the ms1 `annotation_method` branch) while
  # leaving spectral/DB rows and any pre-existing method untouched.
  lvl_ms1 <- trimws(as.character(out_df$annotation_level[idx]))
  is_derived <- !is.na(lvl_ms1) & startsWith(lvl_ms1, "Inferred from")
  if (any(is_derived)) {
    if (!"annotation_method" %in% names(ms1)) {
      ms1$annotation_method <- NA_character_
    }
    cur_am <- trimws(as.character(ms1$annotation_method))
    fill_am <- is_derived & (is.na(cur_am) | !nzchar(cur_am))
    ms1$annotation_method[fill_am] <- "Homologue"
  }
  project <- set_ms1_annotations(project, .nts_ms1_sanitize_and_order(ms1))

  tab <- table(out_df$annotation_level)
  msg <- paste0(
    "> Homologue inference finished - ",
    paste(names(tab), tab, sep = ": ", collapse = " | ")
  )
  message(msg)

  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "infer_homologue_from_kmd_anchors",
    list(counts = as.list(tab))
  )
  project
}

#' KMD-anchor homologue formula inference (no network required)
#'
#' Unpacks semicolon-separated `kmd_series_id` values (multi-base KMD from
#' [calculate_kmd]), selects database anchors, and infers formulas with
#' [adjust_homologue_formula]. For **SMILES propagation on a molecular
#' network graph**, use [infer_homologue_structures] with explicit
#' `network_data`.
#'
#' @param project An [NTSProject].
#'
#' @return Updated `project`.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' d <- tempfile("nts_inf"); dir.create(d)
#' p <- NTSminer::create_nts_project(d)
#' methods::slot(p, "feature_table") <- data.frame(
#'   name = "F1", mz = 400, rt = 1, stringsAsFactors = FALSE
#' )
#' p <- NTSminer::set_ms1_annotations(p, NTSminer:::.nts_init_ms1_from_ft(p@feature_table))
#' ann <- NTSminer::get_ms1_annotations(p)
#' ann$kmd_series_id <- "S1"
#' p <- NTSminer::set_ms1_annotations(p, ann)
#' methods::slot(p, "db_annotations") <- data.frame(
#'   name = "F1", score_total = 0.9, formula = "C8F17",
#'   Compound.name = "A", stringsAsFactors = FALSE
#' )
#' p2 <- suppressMessages(infer_homologue_from_kmd_anchors(p))
#' stopifnot("inferred_formula" %in% names(NTSminer::get_ms1_annotations(p2)))
#' @export
infer_homologue_from_kmd_anchors <- function(project) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject object.", call. = FALSE)
  }
  .infer_homologue_structures_project(project)
}

#' @noRd
.nts_infer_homologue_build_seed_labels <- function(ndf, seed_col) {
  n <- nrow(ndf)
  lab <- rep(NA_character_, n)
  if (!is.null(seed_col) && length(seed_col) && nzchar(trimws(as.character(seed_col[1L])))) {
    cn <- trimws(as.character(seed_col[1L]))
    if (cn %in% names(ndf)) {
      lab <- trimws(as.character(ndf[[cn]]))
    }
  }
  miss <- is.na(lab) | !nzchar(lab)
  if (any(miss)) {
    for (cn in c(
      "library_Compound_name", "Compound.name",
      "library_ref_id", "ref_id", "Library_ID"
    )) {
      if (!cn %in% names(ndf)) {
        next
      }
      fill <- trimws(as.character(ndf[[cn]]))
      m2 <- is.na(lab) | !nzchar(lab)
      ok <- m2 & !is.na(fill) & nzchar(fill)
      if (any(ok)) {
        lab[ok] <- fill[ok]
      }
    }
  }
  lab
}

#' Label network edges with signed target-unit transformations from node `mz`
#'
#' Fills an edge's `transformation` label from the **node `mz` difference** when
#' the edge carries no usable label (absent / `NA` / `"Unknown"`). This lets the
#' algebraic derivation ([infer_homologue_structures]) run on **any** network
#' whose nodes expose `mz` - e.g. the MS2 fragment / neutral-loss /
#' mass-difference *feature networks* - not only the KMD homologue network whose
#' labels come from [build_homologue_network]. Existing non-empty labels are
#' preserved, so this is a pure fallback.
#'
#' @param edges `data.frame` with `from`/`to` node ids.
#' @param nodes `data.frame` with `id` (or `name`) and numeric `mz`.
#' @param target_mz_diff Positive single-step mass (Da), e.g. `49.9968064` (CF2).
#' @param tol_ppm ppm tolerance against `n_steps * target_mz_diff`.
#' @param max_gap Maximum absolute step count (integer `>= 1`).
#' @param target_unit_name Unit token written into `transformation` (e.g. `CF2`).
#' @return `edges` with `transformation` / `n_steps` / `mass_difference` filled.
#' @keywords internal
#' @noRd
.nts_homologue_tag_edges_from_dmz <- function(edges, nodes,
                                              target_mz_diff = 49.9968064,
                                              tol_ppm = 20,
                                              max_gap = 3L,
                                              target_unit_name = "CF2") {
  if (!is.data.frame(edges) || !nrow(edges) ||
      !all(c("from", "to") %in% names(edges))) {
    return(edges)
  }
  ndf <- as.data.frame(nodes, stringsAsFactors = FALSE)
  idc <- if ("id" %in% names(ndf)) {
    "id"
  } else if ("name" %in% names(ndf)) {
    "name"
  } else {
    NULL
  }
  if (is.null(idc) || !"mz" %in% names(ndf)) {
    return(edges)
  }
  mz_map <- suppressWarnings(as.numeric(ndf$mz))
  names(mz_map) <- trimws(as.character(ndf[[idc]]))
  edf <- as.data.frame(edges, stringsAsFactors = FALSE)
  mz_f <- unname(mz_map[trimws(as.character(edf$from))])
  mz_t <- unname(mz_map[trimws(as.character(edf$to))])
  dm <- mz_t - mz_f
  tmd <- suppressWarnings(as.numeric(target_mz_diff)[1L])
  if (!is.finite(tmd) || tmd <= 0) {
    return(edges)
  }
  mg <- max(1L, as.integer(max_gap)[1L])
  if (!is.finite(mg)) mg <- 1L
  ns <- suppressWarnings(as.integer(round(dm / tmd)))
  tol <- suppressWarnings(as.numeric(tol_ppm)[1L])
  if (!is.finite(tol) || tol < 0) tol <- 20
  tol_da <- abs(pmax(abs(mz_f), abs(mz_t))) * (tol * 1e-6) * pmax(1L, abs(ns))
  ok <- !is.na(ns) & ns != 0L & abs(ns) <= mg & abs(dm - ns * tmd) <= tol_da
  unit <- trimws(as.character(target_unit_name)[1L])
  if (is.na(unit) || !nzchar(unit)) unit <- "UNIT"
  lab <- ifelse(ok, sprintf("%+d%s", ns, unit), "Unknown")
  cur <- if ("transformation" %in% names(edf)) {
    trimws(as.character(edf$transformation))
  } else {
    rep(NA_character_, nrow(edf))
  }
  keep <- !is.na(cur) & nzchar(cur) & cur != "Unknown"
  lab[keep] <- cur[keep]
  edf$transformation <- lab
  if (!"mass_difference" %in% names(edf)) {
    edf$mass_difference <- rep(NA_real_, nrow(edf))
  } else {
    edf$mass_difference <- suppressWarnings(as.numeric(edf$mass_difference))
  }
  fill_md <- is.na(edf$mass_difference) & !is.na(dm)
  edf$mass_difference[fill_md] <- dm[fill_md]
  if (!"n_steps" %in% names(edf)) {
    edf$n_steps <- rep(NA_integer_, nrow(edf))
  }
  edf$n_steps[ok] <- ns[ok]
  edf
}

#' Network homologue **algebraic** inference on `network_data` (enhanced list)
#'
#' @description
#' Breadth-first propagation uses **seed compound names / ref ids** on nodes
#' and **signed** edge labels such as `+CF2` / `+1CF2` / `+2CF2` /
#' `-3CH2` (from [build_molecular_network] or [build_homologue_network]
#' `transformation`). The function **does not** modify the [NTSProject]; it
#' returns a copy of `network_data` with homologue columns merged into
#' `nodes` for use with [plot_molecular_network].
#'
#' When the edges carry no `transformation` label (e.g. a raw MS2 feature
#' network from [generate_network_data]), set `derive_from_dmz = TRUE` to
#' derive signed labels from the **node `mz` difference** against
#' `target_mz_diff` (see `.nts_homologue_tag_edges_from_dmz`); the graph is then
#' walked exactly as for the KMD homologue network. Default `FALSE` keeps the
#' previous behaviour untouched.
#'
#' `network_data` must be produced upstream - this function **never**
#' recomputes MS2 similarity internally. For homologue **name** propagation, pass
#' the output of [build_homologue_network] (edges need `transformation`), not only
#' raw [get_network_data] with `type = \"homologue\"`.
#'
#' When `seeds` is `NULL` or empty, seeds default to every `nodes$id` with a
#' non-empty seed label (see `seed_col` and fallback columns).
#'
#' @param project An [NTSProject] with `@ms1_annotations$name` aligned to
#'   `network_data$nodes$id` (or `nodes$name` when aliased to `id`). Read only.
#' @param network_data **Required** named `list` with `nodes` and `edges`
#'   (`data.frame` or `data.table`). Typically from [get_network_data],
#'   [generate_network_data], or [build_homologue_network]; must not be `NULL`.
#' @param seeds Optional character vector of seed feature ids.
#' @param max_steps Maximum number of BFS **waves** (outer iterations): each
#'   wave expands **all** edges whose `from` is any known `(id, label)` state
#'   (seeds plus previously inferred nodes) into `to` nodes not yet annotated.
#'   Set `0L` for no cap (run until a wave adds no new target). Default `3L`.
#' @param seed_col Optional single column name on `nodes` to use as the primary
#'   seed label; if `NULL`, falls back to `library_Compound_name`, then
#'   `Compound.name`, then `library_ref_id` / `ref_id`.
#' @param transformation_units Character vector of allowed edge repeat tokens
#'   (e.g. `CF2`, `CH2`; numeric multipliers in edges such as `+2CF2` use the
#'   same token) for algebraic name propagation; `NULL` uses the built-in
#'   default `c("CF2", "CH2", "O", "H2")`.
#' @param derive_from_dmz If `TRUE`, fill missing edge `transformation` labels
#'   from the node `mz` difference against `target_mz_diff` (via
#'   `.nts_homologue_tag_edges_from_dmz`) so the derivation runs on any network
#'   whose nodes expose `mz` (e.g. MS2 feature networks). Default `FALSE`.
#' @param target_mz_diff Positive single-step mass (Da) used when
#'   `derive_from_dmz = TRUE` (default `49.9968064`, CF2).
#' @param dmz_tol_ppm ppm tolerance for mass-difference matching when
#'   `derive_from_dmz = TRUE` (default `20`).
#' @param max_gap Maximum absolute step count when `derive_from_dmz = TRUE`
#'   (integer `>= 1`, default `3L`).
#' @param target_unit_name Unit token written into derived `transformation`
#'   labels when `derive_from_dmz = TRUE` (default `"CF2"`).
#'
#' @return A **named list** (enhanced `network_data`): same structure as the
#'   input (including any extra list elements), with `edges` unchanged and
#'   `nodes` carrying display-oriented columns `homologue_label` (for
#'   [plot_molecular_network]), `homologue_steps`, `evidence_path`, and
#'   `evidence_source` (`"Seed"` / `"Inferred"`). Intermediate propagation
#'   columns (`inferred_name`, `inferred_formula`, `inferred_smiles`,
#'   `inference_evidence_path`, `homologue_inference_*`, join artifacts such as
#'   `i.inferred_name`), plus redundant XCMS / spectral-library match columns
#'   (`exp_mz`, `lib_mz`, `mz_error_ppm`, `exp_rt`, `lib_rt`, `rt_error_s`,
#'   `n_matched_peaks`, `mzmin`, `mzmax`, `rtmin`, `rtmax`, `library_Compound_name`,
#'   `library_ref_id`, `library_SMILES`) are removed before return when present.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' d <- tempfile("nts_net"); dir.create(d)
#' p <- NTSminer::create_nts_project(d)
#' methods::slot(p, "feature_table") <- data.frame(
#'   name = c("A", "B"), mz = c(300, 350.002), rt = c(1, 2), stringsAsFactors = FALSE
#' )
#' p <- NTSminer::set_ms1_annotations(p, NTSminer:::.nts_init_ms1_from_ft(p@feature_table))
#' net <- list(
#'   nodes = data.frame(
#'     id = c("A", "B"), name = c("A", "B"),
#'     library_Compound_name = c("PFOA", NA_character_),
#'     formula = c("C8F17", NA_character_),
#'     stringsAsFactors = FALSE
#'   ),
#'   edges = data.frame(
#'     from = "A", to = "B", transformation = "+CF2", weight = 0.9,
#'     stringsAsFactors = FALSE
#'   )
#' )
#' net2 <- infer_homologue_structures(p, net, seeds = "A", max_steps = 1L)
#' nd <- net2$nodes
#' stopifnot(nd$homologue_label[nd$id == "B"] == "PFOA + 1 CF2")
#' stopifnot(!"inferred_formula" %in% names(nd))
#' @export
infer_homologue_structures <- function(project,
                                       network_data,
                                       seeds = NULL,
                                       max_steps = 3L,
                                       seed_col = NULL,
                                       transformation_units = NULL,
                                       derive_from_dmz = FALSE,
                                       target_mz_diff = 49.9968064,
                                       dmz_tol_ppm = 20,
                                       max_gap = 3L,
                                       target_unit_name = "CF2") {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject object.", call. = FALSE)
  }
  if (missing(network_data) || is.null(network_data)) {
    stop(
      paste0(
        "infer_homologue_structures requires a precomputed `network_data` list. ",
        "For homologue name propagation use build_homologue_network(proj, ...) first; ",
        "for MS2 graphs use generate_network_data(proj). ",
        "Then: infer_homologue_structures(project, network_data, ...)."
      ),
      call. = FALSE
    )
  }
  if (!is.list(network_data) || !all(c("nodes", "edges") %in% names(network_data))) {
    stop(
      "network_data must be a named list containing `nodes` and `edges` (data.frames or data.tables).",
      call. = FALSE
    )
  }
  nd <- network_data$nodes
  ed <- network_data$edges
  if (is.null(nd) || is.null(ed)) {
    stop("network_data$nodes and network_data$edges must be non-NULL.", call. = FALSE)
  }
  if ((!is.data.frame(nd) && !data.table::is.data.table(nd)) ||
      (!is.data.frame(ed) && !data.table::is.data.table(ed))) {
    stop(
      "network_data$nodes and network_data$edges must be data.frames or data.tables.",
      call. = FALSE
    )
  }
  ndf <- as.data.frame(nd, stringsAsFactors = FALSE)
  edf <- as.data.frame(ed, stringsAsFactors = FALSE)
  if (!nrow(ndf)) {
    stop("network_data$nodes is empty.", call. = FALSE)
  }
  if (!"id" %in% names(ndf) && "name" %in% names(ndf)) {
    ndf$id <- trimws(as.character(ndf$name))
  } else if ("id" %in% names(ndf)) {
    ndf$id <- trimws(as.character(ndf$id))
  } else {
    stop("network_data$nodes must contain column `id` or `name`.", call. = FALSE)
  }

  ndf$.nts_seed_label <- .nts_infer_homologue_build_seed_labels(ndf, seed_col)

  seeds_ch <- unique(as.character(seeds))
  seeds_ch <- seeds_ch[!is.na(seeds_ch) & nzchar(seeds_ch)]
  if (length(seeds_ch)) {
    id0 <- trimws(as.character(ndf$id))
    seeds_ch <- seeds_ch[seeds_ch %in% id0]
    if (!length(seeds_ch)) {
      warning(
        "infer_homologue_structures: explicit `seeds` do not intersect nodes$id; skipped algebraic propagation.",
        call. = FALSE
      )
      return(.nts_homologue_strip_scaffold_node_cols(
        .nts_homologue_enrich_network_data(
          project,
          network_data,
          character(0L),
          list(nodes = ndf, edges = edf)
        )
      ))
    }
  } else {
    sv <- trimws(as.character(ndf$.nts_seed_label))
    if (!any(!is.na(sv) & nzchar(sv))) {
      warning(
        "infer_homologue_structures: no usable seed labels on nodes (library_Compound_name / ref_id, etc.) and no explicit seeds; skipped algebraic propagation.",
        call. = FALSE
      )
      return(.nts_homologue_strip_scaffold_node_cols(
        .nts_homologue_enrich_network_data(
          project,
          network_data,
          character(0L),
          list(nodes = ndf, edges = edf)
        )
      ))
    }
    seeds_ch <- unique(as.character(ndf$id)[!is.na(sv) & nzchar(sv)])
  }

  net_list <- list(nodes = ndf, edges = edf)
  if (isTRUE(derive_from_dmz)) {
    net_list$edges <- .nts_homologue_tag_edges_from_dmz(
      net_list$edges,
      ndf,
      target_mz_diff = target_mz_diff,
      tol_ppm = dmz_tol_ppm,
      max_gap = max_gap,
      target_unit_name = target_unit_name
    )
  }
  out <- .infer_homologue_structures_network(
    net_list,
    seeds = seeds_ch,
    max_steps = max_steps,
    transformation_units = transformation_units
  )
  .nts_homologue_strip_scaffold_node_cols(
    .nts_homologue_enrich_network_data(project, network_data, seeds_ch, out)
  )
}

# ============================================================================
# NAP-style structural inference: mass-window filtering + fingerprint Tanimoto (BFS)
# ============================================================================

#' @noRd
.NTS_PROTON_MASS_NAP <- 1.007276

#' @noRd
.nts_theory_mz_from_exact <- function(exact_mass, polarity, adduct_mass_offset = NULL) {
  em <- suppressWarnings(as.numeric(exact_mass))
  if (!is.null(adduct_mass_offset)) {
    off <- suppressWarnings(as.numeric(adduct_mass_offset)[1L])
    if (!is.finite(off)) {
      return(rep(NA_real_, length(em)))
    }
    return(em + off)
  }
  pol <- tolower(trimws(as.character(polarity))[1L])
  if (pol %in% c("pos", "positive", "+")) {
    return(em + .NTS_PROTON_MASS_NAP)
  }
  em - .NTS_PROTON_MASS_NAP
}

#' @noRd
.nts_prep_structural_db_nap <- function(structural_db, polarity, adduct_mass_offset) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  db <- if (data.table::is.data.table(structural_db)) {
    data.table::copy(structural_db)
  } else {
    data.table::as.data.table(structural_db)
  }
  nm <- names(db)
  fcol <- nm[match("formula", tolower(nm), nomatch = 0L)][1L]
  mcol <- nm[match("exact_mass", tolower(nm), nomatch = 0L)][1L]
  scol <- nm[match("smiles", tolower(nm), nomatch = 0L)][1L]
  if (is.na(fcol) || is.na(mcol) || is.na(scol)) {
    stop(
      "structural_db must contain columns formula, exact_mass, SMILES (case-insensitive).",
      call. = FALSE
    )
  }
  data.table::setnames(db, old = c(fcol, mcol, scol), new = c("formula", "exact_mass", "SMILES"))
  db[, exact_mass := suppressWarnings(as.numeric(exact_mass))]
  db[, formula := as.character(formula)]
  db[, SMILES := trimws(as.character(SMILES))]
  db <- db[is.finite(exact_mass) & nzchar(SMILES)]
  db[, ion_mz := .nts_theory_mz_from_exact(exact_mass, polarity, adduct_mass_offset)]
  db <- db[is.finite(ion_mz)]
  data.table::setorder(db, ion_mz)
  db[]
}

#' @noRd
.nts_ms1_merge_nap_structure_inference <- function(project, inference_dt) {
  if (!methods::is(project, "NTSProject")) {
    stop("internal: NTSProject required", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  if (!data.table::is.data.table(inference_dt) || !nrow(inference_dt) || !"name" %in% names(inference_dt)) {
    return(project)
  }
  ann_df <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
  if (!nrow(ann_df) || !"name" %in% names(ann_df)) {
    stop("ms1_annotations is empty or lacks `name`.", call. = FALSE)
  }
  inf <- data.table::as.data.table(inference_dt)
  inf_map <- inf[, .(
    name = as.character(name),
    nap_structure_smiles = as.character(nap_structure_smiles),
    nap_structure_formula = as.character(nap_structure_formula),
    nap_structure_tanimoto = suppressWarnings(as.numeric(nap_structure_tanimoto)),
    nap_structure_seed = as.character(nap_structure_seed),
    nap_structure_depth = suppressWarnings(as.integer(nap_structure_depth))
  )]
  drop <- intersect(
    c(
      "nap_structure_smiles",
      "nap_structure_formula",
      "nap_structure_tanimoto",
      "nap_structure_seed",
      "nap_structure_depth"
    ),
    names(ann_df)
  )
  if (length(drop)) {
    ann_df <- ann_df[, setdiff(names(ann_df), drop), drop = FALSE]
  }
  merged <- merge(
    ann_df,
    as.data.frame(inf_map, stringsAsFactors = FALSE),
    by = "name",
    all.x = TRUE,
    sort = FALSE
  )
  project <- set_ms1_annotations(project, .nts_ms1_sanitize_and_order(merged))
  project
}

#' @noRd
.nts_infer_structure_node_smiles <- function(ndf, seed_smiles_col) {
  smi <- rep(NA_character_, nrow(ndf))
  if (!is.null(seed_smiles_col) && nzchar(seed_smiles_col) && seed_smiles_col %in% names(ndf)) {
    v <- trimws(as.character(ndf[[seed_smiles_col]]))
    ok <- !is.na(v) & nzchar(v)
    smi[ok] <- v[ok]
    return(smi)
  }
  for (cn in c("library_SMILES", "library_smiles", "SMILES", "smiles", "structure_SMILES")) {
    if (!cn %in% names(ndf)) {
      next
    }
    v <- trimws(as.character(ndf[[cn]]))
    fill <- !is.na(v) & nzchar(v) & (is.na(smi) | !nzchar(smi))
    smi[fill] <- v[fill]
  }
  smi
}

#' Network inference by mass window and structural fingerprint similarity (NAP-style)
#'
#' During BFS propagation, **fingerprints are not computed for the full library**: for each
#' neighbour, observed m/z is matched against `structural_db` within a ppm mass window, then
#' Tanimoto similarity is computed between candidates and the parent SMILES using
#' `rcdk` + `fingerprint::fp.sim.matrix`; above-threshold hits inherit structure and seed the
#' next wave. Default negative-ion mode uses `[M-H]-`: `ion_mz = exact_mass - 1.007276`.
#'
#' @param project [NTSProject] (writes `@ms1_annotations`, may back-fill `mz`).
#' @param network_data Named list with `nodes` and `edges`.
#' @param structural_db `data.table` with `formula`, `exact_mass`, and `SMILES` (column names case-insensitive).
#' @param polarity `"negative"` (default) or `"positive"`; or override with `adduct_mass_offset`.
#' @param adduct_mass_offset If set, theoretical ion m/z = `exact_mass + offset` (custom adducts).
#' @param mass_ppm ppm window matching neighbour m/z to library theoretical ion m/z (default `5`).
#' @param min_tanimoto Minimum Tanimoto to accept an inference (default `0.6`).
#' @param max_steps Maximum BFS depth (edge count).
#' @param seeds Explicit seed `id` values; `NULL` uses all nodes with non-empty SMILES.
#' @param seed_smiles_col Preferred SMILES column name.
#' @param fp_type `type` argument passed to [rcdk::get.fingerprint].
#' @param similarity_method Similarity metric for candidate ranking (default
#'   `"tanimoto"`).
#' @param timeout_sec Per-molecule fingerprint timeout (seconds); slow
#'   molecules are skipped.
#' @param write_to_ms1 If `TRUE`, merge `nap_structure_*` columns into `@ms1_annotations` by `name`.
#' @param return_inference_only If `TRUE`, do not write back `project`; return inference table only.
#'
#' @return **invisibly** `list(inference = data.table, stats = list)`. When
#'   `return_inference_only = FALSE` and `write_to_ms1 = TRUE`, `project` is updated in place.
#'
#' @details
#' Invalid SMILES are skipped via `tryCatch` in parse and fingerprint steps. Large candidate
#' sets can be slow on the JVM fingerprint path; narrow `structural_db` or reduce graph degree.
#'
#' @export
infer_by_structure <- function(project,
                               network_data,
                               structural_db,
                               polarity = c("negative", "positive"),
                               adduct_mass_offset = NULL,
                               mass_ppm = 5,
                               min_tanimoto = 0.6,
                               max_steps = 20L,
                               seeds = NULL,
                               seed_smiles_col = NULL,
                               fp_type = "circular",
                               similarity_method = "tanimoto",
                               timeout_sec = 30L,
                               write_to_ms1 = TRUE,
                               return_inference_only = FALSE) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  polarity <- match.arg(polarity)
  if (missing(network_data) || is.null(network_data)) {
    stop("`network_data` is required (e.g. from build_ms2_network() or get_network_data()).", call. = FALSE)
  }
  if (!is.list(network_data) || !all(c("nodes", "edges") %in% names(network_data))) {
    stop("network_data must contain `nodes` and `edges`.", call. = FALSE)
  }

  nd <- as.data.frame(network_data$nodes, stringsAsFactors = FALSE)
  ed <- as.data.frame(network_data$edges, stringsAsFactors = FALSE)
  if (!nrow(nd)) {
    stop("network_data$nodes is empty.", call. = FALSE)
  }
  if (!"id" %in% names(nd) && "name" %in% names(nd)) {
    nd$id <- trimws(as.character(nd$name))
  } else if ("id" %in% names(nd)) {
    nd$id <- trimws(as.character(nd$id))
  } else {
    stop("network_data$nodes must have `id` or `name`.", call. = FALSE)
  }
  if (!"name" %in% names(nd)) {
    nd$name <- nd$id
  } else {
    nd$name <- trimws(as.character(nd$name))
  }

  nd <- as.data.frame(nd, stringsAsFactors = FALSE)
  if (!"mz" %in% names(nd)) {
    nd$mz <- NA_real_
  }
  nd$mz <- suppressWarnings(as.numeric(nd$mz))
  if (any(!is.finite(nd$mz))) {
    ft <- get_feature_table(project)
    if (is.data.frame(ft) && nrow(ft) && all(c("name", "mz") %in% names(ft))) {
      ftx <- unique(data.table::as.data.table(ft[, c("name", "mz")]), by = "name")
      ndt <- data.table::as.data.table(nd)
      ndt <- merge(ndt, ftx, by = "name", all.x = TRUE, suffixes = c("", ".ft"))
      ndt[!is.finite(mz), mz := suppressWarnings(as.numeric(mz.ft))]
      ndt[, mz.ft := NULL]
      nd <- as.data.frame(ndt, stringsAsFactors = FALSE)
      nd$mz <- suppressWarnings(as.numeric(nd$mz))
    }
  }

  db <- .nts_prep_structural_db_nap(structural_db, polarity = polarity, adduct_mass_offset = adduct_mass_offset)
  if (!nrow(db)) {
    warning("infer_by_structure: structural_db empty after validation; nothing to do.", call. = FALSE, immediate. = TRUE)
    out <- list(
      inference = data.table::data.table(),
      stats = list(n_seeds = 0L, n_assigned = 0L, n_db = 0L)
    )
    return(invisible(out))
  }
  if (nrow(db) > 50000L) {
    warning(
      "infer_by_structure: structural_db has ",
      nrow(db),
      " rows; fingerprint comparisons may be slow.",
      call. = FALSE,
      immediate. = TRUE
    )
  }

  nd$`.nap_seed_smiles` <- .nts_infer_structure_node_smiles(nd, seed_smiles_col)
  id_all <- as.character(nd$id)
  names_m <- as.character(nd$name)

  seeds_ch <- unique(as.character(seeds))
  seeds_ch <- seeds_ch[!is.na(seeds_ch) & nzchar(seeds_ch)]
  if (length(seeds_ch)) {
    seeds_ch <- seeds_ch[seeds_ch %in% id_all]
  } else {
    seeds_ch <- id_all[!is.na(nd$`.nap_seed_smiles`) & nzchar(nd$`.nap_seed_smiles`)]
    seeds_ch <- unique(seeds_ch)
  }
  if (!length(seeds_ch)) {
    warning(
      "infer_by_structure: no seeds with SMILES (set seeds / seed_smiles_col or annotate nodes).",
      call. = FALSE,
      immediate. = TRUE
    )
    out <- list(
      inference = data.table::data.table(),
      stats = list(n_seeds = 0L, n_assigned = 0L, n_db = nrow(db))
    )
    return(invisible(out))
  }

  max_steps <- suppressWarnings(as.integer(max_steps)[1L])
  if (!is.finite(max_steps) || max_steps < 1L) {
    stop("`max_steps` must be a positive integer.", call. = FALSE)
  }
  mass_ppm <- suppressWarnings(as.numeric(mass_ppm)[1L])
  if (!is.finite(mass_ppm) || mass_ppm <= 0) {
    stop("`mass_ppm` must be positive.", call. = FALSE)
  }
  min_tanimoto <- suppressWarnings(as.numeric(min_tanimoto)[1L])
  if (!is.finite(min_tanimoto)) {
    stop("`min_tanimoto` must be finite.", call. = FALSE)
  }

  ed <- data.table::as.data.table(ed)
  if (nrow(ed) && all(c("from", "to") %in% names(ed))) {
    neigh <- data.table::rbindlist(
      list(
        ed[, .(id = trimws(as.character(from)), nb = trimws(as.character(to)))],
        ed[, .(id = trimws(as.character(to)), nb = trimws(as.character(from)))]
      ),
      use.names = TRUE
    )
    neigh <- unique(neigh, by = c("id", "nb"))
    data.table::setkey(neigh, id)
  } else {
    neigh <- data.table::data.table(id = character(0), nb = character(0))
    data.table::setkey(neigh, id)
  }

  assigned_smi <- new.env(parent = emptyenv(), hash = TRUE)
  assigned_score <- new.env(parent = emptyenv(), hash = TRUE)
  assigned_seed <- new.env(parent = emptyenv(), hash = TRUE)
  assigned_depth <- new.env(parent = emptyenv(), hash = TRUE)
  assigned_formula <- new.env(parent = emptyenv(), hash = TRUE)
  fp_cache <- new.env(parent = emptyenv(), hash = TRUE)

  queue_id <- character(0)
  queue_hop <- integer(0)

  push_queue <- function(idv, hopv) {
    idv <- as.character(idv)
    queue_id <<- c(queue_id, idv)
    queue_hop <<- c(queue_hop, hopv)
  }

  for (sid in seeds_ch) {
    sid <- as.character(sid)[1L]
    i_sid <- match(sid, id_all, nomatch = NA_integer_)
    if (is.na(i_sid)) {
      next
    }
    ssmi <- nd$`.nap_seed_smiles`[i_sid]
    if (is.null(ssmi) || length(ssmi) < 1L || is.na(ssmi[1L]) || !nzchar(as.character(ssmi[1L]))) {
      next
    }
    ssmi <- as.character(ssmi[1L])
    assigned_smi[[sid]] <- ssmi
    assigned_score[[sid]] <- 1.0
    assigned_seed[[sid]] <- sid
    assigned_depth[[sid]] <- 0L
    assigned_formula[[sid]] <- NA_character_
    fp <- .nts_fp_object_one(ssmi, fp_type, timeout_sec = timeout_sec)
    if (!is.null(fp)) {
      fp_cache[[sid]] <- fp
    }
    push_queue(sid, 0L)
  }

  if (!length(queue_id)) {
    warning("infer_by_structure: no valid seed SMILES after filtering.", call. = FALSE, immediate. = TRUE)
    out <- list(
      inference = data.table::data.table(),
      stats = list(n_seeds = 0L, n_assigned = 0L, n_db = nrow(db))
    )
    return(invisible(out))
  }

  qh <- 1L
  while (qh <= length(queue_id)) {
    parent <- as.character(queue_id[qh])[1L]
    hop_parent <- queue_hop[qh]
    qh <- qh + 1L
    if (hop_parent >= max_steps) {
      next
    }
    ps <- assigned_smi[[parent]]
    if (is.null(ps) || !nzchar(ps)) {
      next
    }
    fp_parent <- fp_cache[[parent]]
    if (is.null(fp_parent)) {
      fp_parent <- .nts_fp_object_one(ps, fp_type, timeout_sec = timeout_sec)
      if (!is.null(fp_parent)) {
        fp_cache[[parent]] <- fp_parent
      }
    }

    nbs <- if (nrow(neigh)) {
      unique(as.character(neigh[list(parent), nomatch = NULL]$nb))
    } else {
      character(0)
    }
    if (!length(nbs)) {
      next
    }

    for (nb in nbs) {
      nb <- as.character(nb)[1L]
      i_nb <- match(nb, id_all, nomatch = NA_integer_)
      if (is.na(i_nb)) {
        next
      }
      qmz <- suppressWarnings(as.numeric(nd$mz[i_nb]))
      if (!is.finite(qmz)) {
        next
      }
      ref_mz <- max(abs(qmz), 1e-9)
      win <- (mass_ppm / 1e6) * ref_mz
      cand <- db[data.table::between(ion_mz, qmz - win, qmz + win)]
      if (!nrow(cand)) {
        next
      }

      best_sc <- -Inf
      best_row <- NULL
      for (k in seq_len(nrow(cand))) {
        sm_c <- cand$SMILES[k]
        ps_ch <- trimws(as.character(ps))
        sm_ch <- trimws(as.character(sm_c))
        if (nzchar(ps_ch) && nzchar(sm_ch) && identical(ps_ch, sm_ch)) {
          sc <- 1.0
        } else if (is.null(fp_parent)) {
          sc <- NA_real_
        } else {
          fp_c <- .nts_fp_object_one(sm_c, fp_type, timeout_sec = timeout_sec)
          sc <- .nts_fp_similarity_pair(fp_parent, fp_c, method = similarity_method)
        }
        if (is.finite(sc) && sc > best_sc) {
          best_sc <- sc
          best_row <- k
        }
      }
      if (!is.finite(best_sc) || best_sc < min_tanimoto) {
        next
      }

      prev <- assigned_score[[nb]]
      if (!is.null(prev) && is.finite(prev) && prev >= best_sc) {
        next
      }

      assigned_smi[[nb]] <- cand$SMILES[best_row]
      assigned_score[[nb]] <- best_sc
      assigned_seed[[nb]] <- assigned_seed[[parent]]
      assigned_formula[[nb]] <- as.character(cand$formula[best_row])
      new_depth <- hop_parent + 1L
      assigned_depth[[nb]] <- new_depth
      fp_nb <- .nts_fp_object_one(cand$SMILES[best_row], fp_type, timeout_sec = timeout_sec)
      if (!is.null(fp_nb)) {
        fp_cache[[nb]] <- fp_nb
      }
      push_queue(nb, new_depth)
    }
  }

  out_ids <- id_all[vapply(id_all, function(z) !is.null(assigned_smi[[z]]), logical(1L))]
  out_ids <- setdiff(out_ids, seeds_ch)
  inf_rows <- vector("list", length(seeds_ch) + length(out_ids))
  ii <- 0L
  add_row <- function(fid) {
    ii <<- ii + 1L
    idx <- match(fid, id_all)
    nm <- if (!is.na(idx)) names_m[idx] else fid
    fm <- assigned_formula[[fid]]
    if (is.null(fm)) {
      fm <- NA_character_
    }
    inf_rows[[ii]] <<- data.table::data.table(
      name = nm,
      nap_structure_smiles = assigned_smi[[fid]],
      nap_structure_formula = as.character(fm),
      nap_structure_tanimoto = as.numeric(assigned_score[[fid]]),
      nap_structure_seed = as.character(assigned_seed[[fid]]),
      nap_structure_depth = as.integer(assigned_depth[[fid]])
    )
  }
  for (sid in seeds_ch) {
    if (!is.null(assigned_smi[[sid]])) {
      add_row(sid)
    }
  }
  for (fid in out_ids) {
    add_row(fid)
  }
  inference <- if (ii > 0L) {
    data.table::rbindlist(inf_rows[seq_len(ii)], fill = TRUE)
  } else {
    data.table::data.table()
  }
  if (nrow(inference)) {
    data.table::setorder(inference, -nap_structure_tanimoto)
    inference <- inference[, .SD[1L], by = name]
  }

  st <- list(
    n_seeds = length(seeds_ch),
    n_assigned = nrow(inference),
    n_db = nrow(db)
  )

  if (isTRUE(return_inference_only)) {
    return(invisible(list(inference = inference, stats = st, project = project)))
  }
  if (isTRUE(write_to_ms1) && nrow(inference)) {
    project <- .nts_ms1_merge_nap_structure_inference(project, inference)
  }
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "infer_by_structure",
    st
  )
  invisible(list(inference = inference, stats = st, project = project))
}

# ============================================================================
# Environmental transformation / degradate m/z radar (data.table non-equi joins)
# ============================================================================

#' @noRd
.nts_apply_seed_transform <- function(seed_f, sub_f, add_f) {
  seed_f <- as.character(seed_f)[1L]
  if (is.na(seed_f) || !nzchar(seed_f)) {
    return(NA_character_)
  }
  f <- seed_f
  sub_f <- sub_f[1L]
  add_f <- add_f[1L]
  if (!is.na(sub_f) && nzchar(as.character(sub_f))) {
    f <- adjust_homologue_formula(f, -1L, as.character(sub_f))
  }
  if (is.na(f) || !nzchar(f)) {
    return(NA_character_)
  }
  if (!is.na(add_f) && nzchar(as.character(add_f))) {
    f <- adjust_homologue_formula(f, 1L, as.character(add_f))
  }
  f
}

#' Default rule table for [search_transformations] (mass_shift / sub_formula /
#' add_formula schema - distinct from `.nts_default_transformations` used by
#' [propagate_network_annotations]).
#' @noRd
.nts_default_transformation_rules_dt <- function() {
  data.table::data.table(
    rule_name = c(
      "Decarboxylation",
      "Reductive_Defluorination",
      "Hydroxylation",
      "Ether_Cleavage"
    ),
    mass_shift = c(-43.989829, -17.990578, 15.994915, -115.989829),
    sub_formula = c("CO2", "F", NA_character_, "C2F4O"),
    add_formula = c(NA_character_, "H", "O", NA_character_)
  )
}

#' @title Search for environmental transformation products (network_data pipeline)
#'
#' @description
#' Builds Cartesian combinations of seed features (rows with database or
#' inferred Hill formula, and/or non-empty `SMILES`, `library_SMILES` from the
#' merged feature table, or `Annotation`) and transformation rules, predicts
#' product m/z, and matches unknown targets (no such anchors, `is_isotope`
#' false) via `data.table` non-equi join on m/z. With `all_vs_all = TRUE`, or
#' when no anchors exist, uses every MS1 row as a seed. For each target with
#' multiple hits, the match with smallest absolute m/z error is kept.
#' **Does not** modify the [NTSProject]; returns `list(nodes, edges)` suitable
#' for [annotate_transformations] / visualization.
#'
#' @param project An [NTSProject] with `feature_table` and `ms1_annotations`
#'   linked by `name` (optional `db_annotations` merged for `formula` / scores).
#' @param rules Rule `data.table` / `data.frame` with columns `rule_name`,
#'   `mass_shift`, `sub_formula`, `add_formula`. `NULL` uses
#'   `.nts_default_transformation_rules_dt()`.
#' @param rules_dt Deprecated alias for `rules` (positional second argument).
#' @param mz_tol_ppm Symmetric ppm tolerance around predicted product m/z.
#' @param all_vs_all If `TRUE`, every MS1 row is used as a transformation seed
#'   (Cartesian with rules). If `FALSE`, seeds are rows with database/inferred
#'   formula and/or non-empty `SMILES`, `library_SMILES`, or `Annotation`. When
#'   no such anchors exist, automatically falls back to the same all-row seeding
#'   with an informative message.
#'
#' @return `list(nodes = data.table, edges = data.table)` - `nodes` mirrors
#'   MS1 feature rows with optional `transformed_formula`,
#'   `transformation_rule`, `transformation_parent`, `annotation_level` on hits;
#'   `edges` holds one row per accepted hit with `from` = seed id, `to` = target
#'   id, and `transformation_rule` (rule name).
#' @export
search_transformations <- function(project, rules = NULL, rules_dt = NULL, mz_tol_ppm = 5,
                                   all_vs_all = FALSE) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  rules_in <- if (!is.null(rules)) rules else rules_dt
  if (is.null(rules_in)) {
    rules_in <- .nts_default_transformation_rules_dt()
  }
  if (!data.table::is.data.table(rules_in)) {
    rules_in <- data.table::as.data.table(rules_in)
  }
  req <- c("rule_name", "mass_shift", "sub_formula", "add_formula")
  if (!all(req %in% names(rules_in))) {
    stop("`rules` must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
  }

  edges_empty <- data.table::data.table(
    from = character(0L),
    to = character(0L),
    transformation_rule = character(0L),
    transformation_parent = character(0L),
    seed_mz = numeric(0L),
    pred_mz = numeric(0L),
    target_mz = numeric(0L),
    abs_mz_error = numeric(0L),
    mass_shift = numeric(0L),
    sub_formula = character(0L),
    add_formula = character(0L),
    transformed_formula = character(0L),
    annotation_level = character(0L)
  )

  ann_raw <- get_ms1_annotations(project)
  ft <- get_feature_table(project)
  if (!nrow(ann_raw) || !"name" %in% names(ann_raw)) {
    message("> transformation search: ms1_annotations empty or lacks name - skipped.")
    return(list(nodes = data.table::as.data.table(ann_raw), edges = edges_empty))
  }
  if (!nrow(ft) || !"name" %in% names(ft) || !"mz" %in% names(ft)) {
    stop("feature_table must be non-empty with `name` and `mz`", call. = FALSE)
  }

  ann <- data.table::as.data.table(ann_raw)
  db0 <- .nts_db_best_per_name(methods::slot(project, "db_annotations"))
  if (is.data.frame(db0) && nrow(db0)) {
    ddt <- data.table::as.data.table(db0)
    ann <- ddt[ann, on = "name"]
  }
  had_mz <- "mz" %in% names(ann_raw)
  ftz <- unique(data.table::as.data.table(ft[, c("name", "mz"), drop = FALSE]), by = "name")
  data.table::setnames(ftz, c("name", "ft_mz"))
  ann <- ftz[ann, on = "name"]
  if (had_mz) {
    ann[!is.finite(mz), mz := ft_mz]
    ann[, ft_mz := NULL]
  } else {
    data.table::setnames(ann, "ft_mz", "mz")
  }

  ft_struct_cols <- intersect(c("library_SMILES", "SMILES", "Annotation"), names(ft))
  if (length(ft_struct_cols) > 0L) {
    ftx <- unique(
      data.table::as.data.table(ft[, c("name", ft_struct_cols), drop = FALSE]),
      by = "name"
    )
    ann <- merge(ann, ftx, by = "name", all.x = TRUE, suffixes = c("", ".ft"))
    for (cn in ft_struct_cols) {
      sc <- paste0(cn, ".ft")
      if (sc %in% names(ann)) {
        if (cn %in% names(ann) && cn != sc) {
          ann[, (cn) := data.table::fcoalesce(as.character(get(cn)), as.character(get(sc)))]
          ann[, (sc) := NULL]
        } else {
          data.table::setnames(ann, old = sc, new = cn)
        }
      }
    }
  }

  if (!"formula" %in% names(ann)) {
    ann[, formula := NA_character_]
  }
  if (!"inferred_formula" %in% names(ann)) {
    ann[, inferred_formula := NA_character_]
  }
  if (!"is_isotope" %in% names(ann)) {
    ann[, is_isotope := FALSE]
  }
  if (!"annotation_level" %in% names(ann)) {
    ann[, annotation_level := NA_character_]
  }

  ann[, has_db := !is.na(formula) & nzchar(as.character(formula))]
  ann[, has_inf := !is.na(inferred_formula) & nzchar(as.character(inferred_formula))]
  ann[, has_struct := FALSE]
  for (col in c("SMILES", "library_SMILES", "Annotation")) {
    if (col %in% names(ann)) {
      ann[, has_struct := has_struct | (!is.na(get(col)) & nzchar(trimws(as.character(get(col)))))]
    }
  }
  ann[, seed_anchor := has_db | has_inf | has_struct]
  iso <- as.logical(ann$is_isotope)
  iso[is.na(iso)] <- FALSE

  seeds <- ann[seed_anchor == TRUE][
    ,
    .(
      seed_parent = as.character(name),
      seed_mz = as.numeric(mz),
      seed_formula = {
        fd <- as.character(formula)
        fi <- as.character(inferred_formula)
        data.table::fifelse(
          has_db & nzchar(fd),
          fd,
          data.table::fifelse(has_inf & nzchar(fi), fi, NA_character_)
        )
      },
      only_struct = has_struct & !has_db & !has_inf
    )
  ]
  seeds <- seeds[is.finite(seed_mz) & (nzchar(as.character(seed_formula)) | only_struct == TRUE)]
  if ("only_struct" %in% names(seeds)) {
    seeds[, only_struct := NULL]
  }

  if (isTRUE(all_vs_all) || nrow(seeds) == 0L) {
    message(
      "\u25b6 \u5f00\u542f\u5168\u5c40\u65e0\u5dee\u522b\u5339\u914d\u6a21\u5f0f (all_vs_all = TRUE \u6216\u672a\u53d1\u73b0\u951a\u70b9)\u3002"
    )
    seeds <- ann[
      ,
      .(
        seed_parent = as.character(name),
        seed_mz = as.numeric(mz),
        seed_formula = {
          fd <- as.character(formula)
          fi <- as.character(inferred_formula)
          data.table::fifelse(
            has_db & nzchar(fd),
            fd,
            data.table::fifelse(has_inf & nzchar(fi), fi, NA_character_)
          )
        }
      )
    ]
    seeds <- seeds[is.finite(seed_mz)]
  } else {
    message(sprintf(
      "\u25b6 \u627e\u5230 %d \u4e2a\u5b9e\u8bc1\u951a\u70b9 (Seeds) \u4f5c\u4e3a\u8f6c\u5316\u7f51\u7edc\u8d77\u70b9\u3002",
      nrow(seeds)
    ))
  }

  targets <- ann[
    !seed_anchor & !iso,
    .(name = as.character(name), mz = as.numeric(mz))
  ]
  targets <- targets[is.finite(mz)]

  finalize_nodes <- function(ann_dt, edges_dt) {
    if ("has_db" %in% names(ann_dt)) {
      ann_dt[, has_db := NULL]
    }
    if ("has_inf" %in% names(ann_dt)) {
      ann_dt[, has_inf := NULL]
    }
    if ("has_struct" %in% names(ann_dt)) {
      ann_dt[, has_struct := NULL]
    }
    if ("seed_anchor" %in% names(ann_dt)) {
      ann_dt[, seed_anchor := NULL]
    }
    if (!"transformed_formula" %in% names(ann_dt)) {
      ann_dt[, transformed_formula := NA_character_]
    }
    if (!"transformation_rule" %in% names(ann_dt)) {
      ann_dt[, transformation_rule := NA_character_]
    }
    if (!"transformation_parent" %in% names(ann_dt)) {
      ann_dt[, transformation_parent := NA_character_]
    }
    .nts_network_data_strip_intermediate_cols(
      list(nodes = ann_dt, edges = edges_dt),
      strip_edges = TRUE
    )
  }

  if (!nrow(seeds) || !nrow(targets)) {
    message(
      "> transformation search: seeds = ",
      nrow(seeds),
      ", targets = ",
      nrow(targets),
      " - no Cartesian match pipeline run."
    )
    return(finalize_nodes(data.table::copy(ann), edges_empty))
  }

  seeds[, join_k := 1L]
  rules_cj <- data.table::copy(rules_in)
  rules_cj[, join_k := 1L]
  expanded <- merge(seeds, rules_cj, by = "join_k", allow.cartesian = TRUE)
  expanded[, join_k := NULL]
  expanded[, pred_mz := seed_mz + as.numeric(mass_shift)]
  tol <- as.numeric(mz_tol_ppm)[1L] * 1e-6
  if (!is.finite(tol) || tol < 0) {
    stop("mz_tol_ppm must be a non-negative finite number", call. = FALSE)
  }
  expanded[, `:=`(
    mz_lower = pred_mz * (1 - tol),
    mz_upper = pred_mz * (1 + tol)
  )]

  data.table::setkey(targets, mz)
  hits <- targets[
    expanded,
    on = .(mz >= mz_lower, mz <= mz_upper),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]

  if (!nrow(hits)) {
    message("> transformation search: 0 hits (m/z within ", mz_tol_ppm, " ppm).")
    return(finalize_nodes(data.table::copy(ann), edges_empty))
  }

  hits[, err := abs(as.numeric(mz) - as.numeric(pred_mz))]
  data.table::setorder(hits, name, err)
  best <- hits[, .SD[1L], by = name]

  best[, transformed_formula := as.character(mapply(
    .nts_apply_seed_transform,
    seed_formula,
    sub_formula,
    add_formula,
    SIMPLIFY = TRUE,
    USE.NAMES = FALSE
  ))]
  best[, transformation_rule := as.character(rule_name)]
  best[, transformation_parent := as.character(seed_parent)]
  best[, annotation_level := paste0(
    "Transformation: ",
    transformation_rule,
    " from ",
    transformation_parent
  )]

  ann[, c("has_db", "has_inf") := NULL]

  if (!"transformed_formula" %in% names(ann)) {
    ann[, transformed_formula := NA_character_]
  }
  if (!"transformation_rule" %in% names(ann)) {
    ann[, transformation_rule := NA_character_]
  }
  if (!"transformation_parent" %in% names(ann)) {
    ann[, transformation_parent := NA_character_]
  }

  ann[
    best,
    on = "name",
    `:=`(
      transformed_formula = i.transformed_formula,
      transformation_rule = i.transformation_rule,
      transformation_parent = i.transformation_parent,
      annotation_level = i.annotation_level
    )
  ]

  n_hit <- nrow(best)
  n_ok <- sum(!is.na(best$transformed_formula) & nzchar(as.character(best$transformed_formula)))
  message(
    "> transformation search: ",
    n_hit,
    " target(s) matched; ",
    n_ok,
    " with valid Hill formula after sub/add."
  )

  edges_dt <- best[, .(
    from = as.character(seed_parent),
    to = as.character(name),
    transformation_rule = as.character(transformation_rule),
    transformation_parent = as.character(transformation_parent),
    seed_mz = as.numeric(seed_mz),
    pred_mz = as.numeric(pred_mz),
    target_mz = as.numeric(mz),
    abs_mz_error = as.numeric(err),
    mass_shift = as.numeric(mass_shift),
    sub_formula = as.character(sub_formula),
    add_formula = as.character(add_formula),
    transformed_formula = as.character(transformed_formula),
    annotation_level = as.character(annotation_level)
  )]

  finalize_nodes(ann, edges_dt)
}
