#' Collapse character vectors for merged annotation tables
#' @keywords internal
#' @noRd
.nts_collapse_semi <- function(x, sep = "; ") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) {
    return(NA_character_)
  }
  paste(x, collapse = sep)
}


#' Long annotation rows from [NTSProject] spectral `db_annotations`
#' @keywords internal
#' @noRd
.nts_annotation_layer_from_db <- function(project) {
  empty <- data.frame(
    feature = character(0),
    ref_id = character(0),
    compound_name = character(0),
    smiles = character(0),
    score = numeric(0),
    propagation_level = integer(0),
    annotation_method = character(0),
    stringsAsFactors = FALSE
  )
  if (!methods::is(project, "NTSProject")) {
    return(empty)
  }
  db <- .nts_db_best_per_name(.nts_get_project_db_annotations(project))
  if (!is.data.frame(db) || !nrow(db)) {
    return(empty)
  }
  feat_col <- if ("name" %in% names(db)) {
    "name"
  } else if ("feature_id" %in% names(db)) {
    "feature_id"
  } else {
    return(empty)
  }
  feat <- trimws(as.character(db[[feat_col]]))
  rid <- if ("ref_id" %in% names(db)) {
    trimws(as.character(db$ref_id))
  } else if ("Lab.ID" %in% names(db)) {
    trimws(as.character(db[["Lab.ID"]]))
  } else {
    NA_character_
  }
  cname <- rep(NA_character_, nrow(db))
  for (col in c("Compound.name", "compound_name")) {
    if (col %in% names(db)) {
      v <- trimws(as.character(db[[col]]))
      miss <- is.na(cname) | !nzchar(cname)
      cname[miss] <- v[miss]
    }
  }
  smi <- rep(NA_character_, nrow(db))
  for (col in c("SMILES", "smiles")) {
    if (col %in% names(db)) {
      v <- trimws(as.character(db[[col]]))
      miss <- is.na(smi) | !nzchar(smi)
      smi[miss] <- v[miss]
    }
  }
  sc <- if ("score_total" %in% names(db)) {
    suppressWarnings(as.numeric(db$score_total))
  } else {
    rep(NA_real_, nrow(db))
  }
  ok <- !is.na(feat) & nzchar(feat) & !is.na(rid) & nzchar(rid)
  if (!any(ok)) {
    return(empty)
  }
  data.frame(
    feature = feat[ok],
    ref_id = rid[ok],
    compound_name = cname[ok],
    smiles = smi[ok],
    score = sc[ok],
    propagation_level = NA_integer_,
    annotation_method = "Spectral_DB",
    stringsAsFactors = FALSE
  )
}


#' Long annotation rows from annotated `network_data` (e.g. [annotate_by_dual_network])
#' @keywords internal
#' @noRd
.nts_annotation_layer_from_net <- function(net, method_label) {
  empty <- data.frame(
    feature = character(0),
    ref_id = character(0),
    compound_name = character(0),
    smiles = character(0),
    score = numeric(0),
    propagation_level = integer(0),
    annotation_method = character(0),
    stringsAsFactors = FALSE
  )
  if (is.null(net) || !is.list(net) || is.null(net$nodes) || !is.data.frame(net$nodes)) {
    return(empty)
  }
  nd <- as.data.frame(net$nodes, stringsAsFactors = FALSE)
  if (!nrow(nd)) {
    return(empty)
  }
  feat <- if ("id" %in% names(nd)) {
    trimws(as.character(nd$id))
  } else if ("name" %in% names(nd)) {
    trimws(as.character(nd$name))
  } else {
    return(empty)
  }
  rid <- rep(NA_character_, nrow(nd))
  if ("propagated_ID" %in% names(nd)) {
    rid <- trimws(as.character(nd$propagated_ID))
  }
  for (col in c("ref_id", "Library_ID")) {
    if (!col %in% names(nd)) {
      next
    }
    v <- trimws(as.character(nd[[col]]))
    miss <- is.na(rid) | !nzchar(rid)
    rid[miss] <- v[miss]
  }
  cname <- rep(NA_character_, nrow(nd))
  for (col in c("Compound.name", "Annotation")) {
    if (!col %in% names(nd)) {
      next
    }
    v <- trimws(as.character(nd[[col]]))
    miss <- is.na(cname) | !nzchar(cname)
    cname[miss] <- v[miss]
  }
  smi <- rep(NA_character_, nrow(nd))
  for (col in c("propagated_SMILES", "SMILES", "smiles")) {
    if (!col %in% names(nd)) {
      next
    }
    v <- trimws(as.character(nd[[col]]))
    miss <- is.na(smi) | !nzchar(smi)
    smi[miss] <- v[miss]
  }
  sc <- if ("propagated_similarity" %in% names(nd)) {
    suppressWarnings(as.numeric(nd$propagated_similarity))
  } else if ("score_total" %in% names(nd)) {
    suppressWarnings(as.numeric(nd$score_total))
  } else {
    rep(NA_real_, nrow(nd))
  }
  pl <- if ("Propagation_Level" %in% names(nd)) {
    suppressWarnings(as.integer(nd$Propagation_Level))
  } else {
    rep(NA_integer_, nrow(nd))
  }
  # Homologue network nodes carry `evidence_source` (Seed / Inferred) and their
  # `Propagation_Level` column, when present, is entirely NA. Fall back to
  # evidence_source so count_annotation_layers() can distinguish level-0 (seed)
  # from propagated rows (P2-E5).
  if (identical(method_label, "Homologue") && all(is.na(pl)) &&
      "evidence_source" %in% names(nd)) {
    es <- trimws(as.character(nd$evidence_source))
    pl <- ifelse(!is.na(es) & es == "Seed", 0L,
                 ifelse(!is.na(es) & es == "Inferred", 1L, NA_integer_))
  }
  if (identical(method_label, "Homologue") && "homologue_label" %in% names(nd)) {
    hl <- trimws(as.character(nd$homologue_label))
    miss <- is.na(cname) | !nzchar(cname)
    cname[miss] <- hl[miss]
  }
  ok <- !is.na(feat) & nzchar(feat) & !is.na(rid) & nzchar(rid)
  if (identical(method_label, "Homologue")) {
    ok <- ok | (
      !is.na(feat) & nzchar(feat) &
        !is.na(cname) & nzchar(cname)
    )
  }
  if (!any(ok)) {
    return(empty)
  }
  data.frame(
    feature = feat[ok],
    ref_id = rid[ok],
    compound_name = cname[ok],
    smiles = smi[ok],
    score = sc[ok],
    propagation_level = pl[ok],
    annotation_method = method_label,
    stringsAsFactors = FALSE
  )
}


#' Merge spectral and network annotation layers into one feature table
#'
#' Combines Top-1 spectral hits from `@db_annotations` ([batch_annotate] /
#' [annotate_by_spectral_match]) with optional annotated molecular networks
#' (outputs of [generate_network_data] + [annotate_by_dual_network], or
#' [infer_homologue_structures] when passed as `net_homologue_ann`).
#'
#' For each feature, values from multiple sources are de-duplicated and joined
#' with `"; "` in wide columns. `annotation_method` lists contributing sources
#' (e.g. `Spectral_DB; Network_DF; Network_MD`). Per-method `ref_id_*` columns
#' are added for auditing disagreements between layers.
#'
#' @param project An [NTSProject] (used for `db_annotations` and, when
#'   `include_all_features = TRUE`, the feature table).
#' @param net_frag_ann Optional `list(nodes, edges)` after fragment-network
#'   propagation. Labelled `Network_DF` in `annotation_method`.
#' @param net_nl_ann Optional neutral-loss network; labelled `Network_NL`.
#' @param net_md_ann Optional mass-difference network; labelled `Network_MD`.
#' @param net_homologue_ann Optional homologue / KMD network (e.g. from
#'   [infer_homologue_structures]); uses `ref_id`, `Library_ID`, `homologue_label`,
#'   and `Compound.name`. Labelled `Homologue`.
#' @param include_all_features If `TRUE` (default), all feature names from the
#'   project feature table are retained; features without any annotation have
#'   `NA` in merged columns.
#' @param method_sep Separator between methods or values when collapsing (default
#'   `"; "`).
#' @param write_to_ms1 If `TRUE`, after building the wide table, merge columns
#'   onto `@ms1_annotations` and attach the updated [NTSProject] as attribute
#'   `"project"` on the returned table. You **must** reassign:
#'   `project <- attr(annotation_master, "project")` (R copies S4 objects on
#'   slot replacement inside functions). Prefer
#'   `project <- apply_annotation_master_to_ms1(project, annotation_master)` when
#'   you already hold the wide table.
#'
#' @return A `data.frame` with columns `feature`, `annotation_method`, `ref_id`,
#'   `compound_name`, `smiles`, `score`, `propagation_level`, plus one
#'   `ref_id_<Method>` column per contributing method.
#'
#' @seealso [count_annotation_layers] for per-layer feature counts,
#'   [apply_annotation_master_to_ms1] to persist the wide table on the project.
#'
#' @export
merge_annotation_layers <- function(project,
                                    net_frag_ann = NULL,
                                    net_nl_ann = NULL,
                                    net_md_ann = NULL,
                                    net_homologue_ann = NULL,
                                    include_all_features = TRUE,
                                    method_sep = "; ",
                                    write_to_ms1 = FALSE) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject.", call. = FALSE)
  }
  method_sep <- as.character(method_sep)[1L]
  if (!nzchar(method_sep)) {
    stop("`method_sep` must be non-empty.", call. = FALSE)
  }

  long <- rbind(
    .nts_annotation_layer_from_db(project),
    .nts_annotation_layer_from_net(net_frag_ann, "Network_DF"),
    .nts_annotation_layer_from_net(net_nl_ann, "Network_NL"),
    .nts_annotation_layer_from_net(net_md_ann, "Network_MD"),
    .nts_annotation_layer_from_net(net_homologue_ann, "Homologue")
  )
  if (!nrow(long)) {
    warning("merge_annotation_layers: no annotation rows found.", call. = FALSE)
    out <- data.frame(
      feature = character(0),
      annotation_method = character(0),
      ref_id = character(0),
      compound_name = character(0),
      smiles = character(0),
      score = numeric(0),
      propagation_level = character(0),
      stringsAsFactors = FALSE
    )
    if (isTRUE(include_all_features)) {
      ft <- get_feature_table(project)
      if (nrow(ft) && "name" %in% names(ft)) {
        out <- data.frame(
          feature = trimws(as.character(ft$name)),
          stringsAsFactors = FALSE
        )
      }
    }
    if (isTRUE(write_to_ms1) && nrow(out)) {
      project <- apply_annotation_master_to_ms1(project, annotation_master = out)
      attr(out, "project") <- project
    }
    return(out)
  }

  split_cols <- c(
    "ref_id", "compound_name", "smiles", "score", "propagation_level",
    "annotation_method"
  )
  agg <- lapply(split_cols, function(col) {
    v <- tapply(long[[col]], long$feature, .nts_collapse_semi, sep = method_sep)
    data.frame(feature = names(v), val = unname(v), stringsAsFactors = FALSE)
  })
  names(agg) <- split_cols
  out <- agg[[1L]][, "feature", drop = FALSE]
  for (col in split_cols) {
    out[[col]] <- agg[[col]]$val
  }

  methods_used <- unique(as.character(long$annotation_method))
  for (m in methods_used) {
    sub <- long[long$annotation_method == m, , drop = FALSE]
    sub_agg <- tapply(sub$ref_id, sub$feature, .nts_collapse_semi, sep = method_sep)
    col_nm <- paste0("ref_id_", gsub("[^A-Za-z0-9]+", "_", m))
    ref_df <- data.frame(
      feature = names(sub_agg),
      val = unname(sub_agg),
      stringsAsFactors = FALSE
    )
    names(ref_df)[2L] <- col_nm
    out <- merge(out, ref_df, by = "feature", all.x = TRUE)
  }

  if (isTRUE(include_all_features)) {
    ft <- get_feature_table(project)
    if (nrow(ft) && "name" %in% names(ft)) {
      all_f <- data.frame(
        feature = unique(trimws(as.character(ft$name))),
        stringsAsFactors = FALSE
      )
      out <- merge(all_f, out, by = "feature", all.x = TRUE)
    }
  }

  out <- out[order(out$feature), , drop = FALSE]
  rownames(out) <- NULL
  if (isTRUE(write_to_ms1)) {
    project <- apply_annotation_master_to_ms1(project, annotation_master = out)
    attr(out, "project") <- project
  }
  out
}


#' Merge `annotation_master` columns into `@ms1_annotations`
#'
#' Left-joins the wide table from [merge_annotation_layers] (or an equivalent
#' `data.frame`) onto `@ms1_annotations` by `feature` / `name`. Column `smiles`
#' is stored as `master_smiles` because bare `SMILES` / `smiles` are stripped
#' from the MS1 tier. Per-layer `ref_id_*` audit columns are copied as-is.
#'
#' @param project An [NTSProject] with non-empty `@ms1_annotations$name`.
#' @param annotation_master Optional wide table. If `NULL`, built by calling
#'   [merge_annotation_layers] with `...`.
#' @param overwrite If `TRUE` (default), replace existing merged columns; if
#'   `FALSE`, only fill `NA` or empty character values.
#' @param ... Passed to [merge_annotation_layers] when `annotation_master` is
#'   `NULL`.
#'
#' @return The updated [NTSProject] (always assign back:
#'   `project <- apply_annotation_master_to_ms1(project, ...)`).
#'
#' @seealso [merge_annotation_layers]
#'
#' @export
apply_annotation_master_to_ms1 <- function(project,
                                           annotation_master = NULL,
                                           overwrite = TRUE,
                                           ...) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject.", call. = FALSE)
  }
  if (is.null(annotation_master)) {
    annotation_master <- merge_annotation_layers(project, ...)
  }
  if (!is.data.frame(annotation_master) || !nrow(annotation_master)) {
    warning(
      "apply_annotation_master_to_ms1: empty annotation_master - ms1_annotations unchanged.",
      call. = FALSE
    )
    return(invisible(project))
  }
  feat_col <- if ("feature" %in% names(annotation_master)) {
    "feature"
  } else if ("name" %in% names(annotation_master)) {
    "name"
  } else {
    stop(
      "annotation_master must contain column `feature` or `name`.",
      call. = FALSE
    )
  }

  ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
  if (!nrow(ms1) || !"name" %in% names(ms1)) {
    stop("ms1_annotations is empty or lacks `name`.", call. = FALSE)
  }

  master <- as.data.frame(annotation_master, stringsAsFactors = FALSE)
  master$.nts_join_name <- trimws(as.character(master[[feat_col]]))
  master <- master[!is.na(master$.nts_join_name) & nzchar(master$.nts_join_name), , drop = FALSE]
  if (!nrow(master)) {
    warning(
      "apply_annotation_master_to_ms1: no valid feature keys - ms1_annotations unchanged.",
      call. = FALSE
    )
    return(invisible(project))
  }

  rename_map <- c(
    smiles = "master_smiles",
    score = "master_score",
    propagation_level = "master_propagation_level"
  )
  skip_cols <- c(feat_col, ".nts_join_name")
  bring <- setdiff(names(master), skip_cols)
  for (src in names(rename_map)) {
    if (src %in% names(master)) {
      names(master)[names(master) == src] <- rename_map[[src]]
    }
  }
  bring <- setdiff(names(master), c(".nts_join_name", skip_cols))
  if (!length(bring)) {
    warning(
      "apply_annotation_master_to_ms1: no columns to merge - ms1_annotations unchanged.",
      call. = FALSE
    )
    return(invisible(project))
  }

  ms1_names <- trimws(as.character(ms1$name))
  m_names <- master$.nts_join_name
  idx <- match(ms1_names, m_names)
  for (cl in bring) {
    newv <- master[[cl]][idx]
    if (!cl %in% names(ms1) || isTRUE(overwrite)) {
      ms1[[cl]] <- newv
      next
    }
    cur <- ms1[[cl]]
    if (is.character(cur)) {
      empty <- is.na(cur) | !nzchar(trimws(as.character(cur)))
      cur[empty] <- newv[empty]
      ms1[[cl]] <- cur
    } else {
      empty <- is.na(cur)
      cur[empty] <- newv[empty]
      ms1[[cl]] <- cur
    }
  }

  project <- set_ms1_annotations(project, .nts_ms1_sanitize_and_order(ms1))
  invisible(project)
}


#' First non-empty value across candidate columns
#'
#' Returns, per row of `nd`, the first value found among `candidates` (in order)
#' that is non-missing / non-empty. Used to read a node column that different
#' network producers name differently (e.g. `propagated_ID` vs `ref_id`).
#'
#' @param nd A `data.frame` of network nodes.
#' @param candidates Character vector of candidate column names, most preferred
#'   first.
#' @param mode `"chr"` (default) returns character; `"num"` coerces with
#'   [as.numeric].
#' @return A vector of length `nrow(nd)`.
#' @keywords internal
#' @noRd
.nts_net_col_first <- function(nd, candidates, mode = c("chr", "num")) {
  mode <- match.arg(mode)
  n <- nrow(nd)
  if (mode == "num") {
    out <- rep(NA_real_, n)
    for (col in candidates) {
      if (!col %in% names(nd)) next
      v <- suppressWarnings(as.numeric(nd[[col]]))
      miss <- is.na(out) & !is.na(v)
      out[miss] <- v[miss]
    }
  } else {
    out <- rep(NA_character_, n)
    for (col in candidates) {
      if (!col %in% names(nd)) next
      v <- trimws(as.character(nd[[col]]))
      miss <- (is.na(out) | !nzchar(out)) & !is.na(v) & nzchar(v)
      out[miss] <- v[miss]
    }
  }
  out
}


#' Node-level results of an annotated network, prefixed for one MS2 edge layer
#'
#' Projects the per-feature **propagated** identity of an annotated feature
#' network (output of [annotate_by_dual_network] / [propagate_network_annotations],
#' both of which write `propagated_ID` / `propagated_SMILES` only for nodes whose
#' identity was *inferred* by network propagation, and leave seeds at
#' `Propagation_Level == 0` with `propagated_ID == NA`) into five
#' `net_<prefix>_*` columns keyed by feature `name`: `_ref_id` (`propagated_ID`),
#' `_smiles` (`propagated_SMILES`), `_score` (`propagated_similarity`), `_hops`
#' (`Propagation_Level`) and `_seed` (`propagated_from_node`, the donor node).
#' The seed's own library identity is deliberately NOT captured here: a direct
#' library hit is not network-propagation evidence (it already lives on
#' `@db_annotations` as L1/L2), so only genuinely propagated nodes are returned.
#'
#' @param net Annotated `list(nodes, edges)` or `NULL`.
#' @param prefix Layer prefix (`"df"` / `"nl"` / `"md"`).
#' @return A 5+1-column `data.frame` (`name` + the five `net_<prefix>_*`), 0-row
#'   when `net` is `NULL` or carries no propagated (`propagated_ID`) identity.
#' @keywords internal
#' @noRd
.nts_net_node_layer_prefixed <- function(net, prefix) {
  nms <- paste0("net_", prefix, c("_ref_id", "_smiles", "_score", "_hops", "_seed"))
  empty <- data.frame(name = character(0), stringsAsFactors = FALSE)
  empty[[nms[1L]]] <- character(0)
  empty[[nms[2L]]] <- character(0)
  empty[[nms[3L]]] <- numeric(0)
  empty[[nms[4L]]] <- integer(0)
  empty[[nms[5L]]] <- character(0)
  if (is.null(net) || !is.list(net) || is.null(net$nodes) || !is.data.frame(net$nodes)) {
    return(empty)
  }
  nd <- as.data.frame(net$nodes, stringsAsFactors = FALSE)
  if (!nrow(nd)) {
    return(empty)
  }
  nm <- if ("id" %in% names(nd)) {
    trimws(as.character(nd$id))
  } else if ("name" %in% names(nd)) {
    trimws(as.character(nd$name))
  } else {
    return(empty)
  }
  ## Only the *propagated* identity (not the seed's own Library_ID/ref_id).
  ref <- .nts_net_col_first(nd, c("propagated_ID"))
  smi <- .nts_net_col_first(nd, c("propagated_SMILES", "propagated_smiles"))
  scr <- .nts_net_col_first(nd, c("propagated_similarity", "score_total"), mode = "num")
  hop <- .nts_net_col_first(nd, c("Propagation_Level", "propagation_level"), mode = "num")
  sed <- .nts_net_col_first(nd, c("propagated_from_node", "seed"))
  ok <- !is.na(nm) & nzchar(nm) &
    ((!is.na(ref) & nzchar(ref)) | (!is.na(smi) & nzchar(smi)))
  if (!any(ok)) {
    return(empty)
  }
  out <- data.frame(name = nm[ok], stringsAsFactors = FALSE)
  out[[nms[1L]]] <- ref[ok]
  out[[nms[2L]]] <- smi[ok]
  out[[nms[3L]]] <- scr[ok]
  out[[nms[4L]]] <- as.integer(hop[ok])
  out[[nms[5L]]] <- sed[ok]
  out
}


#' Left-join MS2 network columns onto `@ms2_annotations` by `name`
#' @keywords internal
#' @noRd
.nts_ms2_write_net_columns <- function(project, out) {
  ms2 <- as.data.frame(get_ms2_annotations(project), stringsAsFactors = FALSE)
  if (!nrow(ms2) || !"name" %in% names(ms2)) {
    return(project)
  }
  idx <- match(trimws(as.character(ms2$name)), out$name)
  for (cl in setdiff(names(out), "name")) {
    ms2[[cl]] <- out[[cl]][idx]
  }
  set_ms2_annotations(project, ms2)
}


#' Merge MS2 feature-network layers into `@ms2_annotations`
#'
#' Split-landing counterpart of [merge_annotation_layers] for the **MS2 tier**:
#' projects the per-feature propagated identity of the DFI / neutral-loss /
#' mass-difference feature networks (output of [annotate_by_dual_network],
#' [propagate_network_annotations]) into `net_df_*` / `net_nl_*` / `net_md_*`
#' columns and, when `write_back = TRUE`, left-joins them onto
#' `@ms2_annotations` by `name`. These are MS2-derived structural evidence and
#' belong on the MS2 tier (the homologue/MS1 evidence goes to
#' `@ms1_annotations$network` via [merge_ms1_homologue_layers]).
#'
#' Two summary columns are added: `net_ms2_methods` lists which edge types
#' supplied a `ref_id` (`"DFI"` / `"NL"` / `"MD"`, `"; "`-joined), and
#' `net_ms2_agreement` is `1` when all present edge types agree on the same
#' `ref_id`, `0` when they disagree, `NA` when none is present.
#'
#' @param project An [NTSProject].
#' @param net_frag_ann,net_nl_ann,net_md_ann Annotated feature-network
#'   `list(nodes, edges)` for the DFI / NL / MD layer, or `NULL`.
#' @param write_back If `TRUE`, merge the columns onto `@ms2_annotations` and
#'   attach the updated project as attribute `"project"` on the returned table
#'   (reassign: `project <- attr(out, "project")`). Default `FALSE` only returns
#'   the table (no slot write).
#'
#' @return A `data.frame` with `name` plus the `net_df_*` / `net_nl_*` /
#'   `net_md_*` / `net_ms2_methods` / `net_ms2_agreement` columns.
#'
#' @seealso [merge_ms1_homologue_layers], [merge_annotation_layers],
#'   [count_annotation_layers]
#' @export
merge_ms2_network_layers <- function(project,
                                     net_frag_ann = NULL,
                                     net_nl_ann = NULL,
                                     net_md_ann = NULL,
                                     write_back = FALSE) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject.", call. = FALSE)
  }
  layers <- list(
    df = .nts_net_node_layer_prefixed(net_frag_ann, "df"),
    nl = .nts_net_node_layer_prefixed(net_nl_ann, "nl"),
    md = .nts_net_node_layer_prefixed(net_md_ann, "md")
  )
  all_names <- unique(unlist(lapply(layers, function(z) z$name), use.names = FALSE))
  out <- data.frame(name = as.character(all_names), stringsAsFactors = FALSE)
  for (pref in c("df", "nl", "md")) {
    z <- layers[[pref]]
    idx <- match(out$name, z$name)
    for (cl in setdiff(names(z), "name")) {
      out[[cl]] <- z[[cl]][idx]
    }
  }
  lab <- c("DFI", "NL", "MD")
  ref_cols <- paste0("net_", c("df", "nl", "md"), "_ref_id")
  present <- vapply(ref_cols, function(cc) {
    v <- out[[cc]]
    !is.na(v) & nzchar(v)
  }, logical(nrow(out)))
  if (nrow(out) && !is.matrix(present)) {
    present <- matrix(present, nrow = nrow(out))
  }
  out$net_ms2_methods <- if (!nrow(out)) {
    character(0)
  } else {
    vapply(seq_len(nrow(out)), function(i) {
      p <- lab[present[i, ]]
      if (!length(p)) NA_character_ else paste(p, collapse = ";")
    }, character(1L))
  }
  out$net_ms2_agreement <- if (!nrow(out)) {
    integer(0)
  } else {
    vapply(seq_len(nrow(out)), function(i) {
      ids <- vapply(ref_cols, function(cc) out[[cc]][i], character(1L))
      ids <- unique(ids[!is.na(ids) & nzchar(ids)])
      if (!length(ids)) NA_integer_ else as.integer(length(ids) == 1L)
    }, integer(1L))
  }
  out <- out[order(out$name), , drop = FALSE]
  rownames(out) <- NULL
  if (isTRUE(write_back)) {
    project <- .nts_ms2_write_net_columns(project, out)
    attr(out, "project") <- project
  }
  out
}


#' Merge homologue-network annotations into the MS1 `network` step
#'
#' Split-landing counterpart of [merge_annotation_layers] for the **MS1 tier**:
#' projects the per-feature homologue result of an enriched homologue network
#' (output of [infer_homologue_structures] / [build_homologue_network]) into
#' `homologue_*` columns and, when `write_back = TRUE`, writes them into the
#' `network` step of `@ms1_annotations` (registered in `.nts_ms1_groups()`).
#'
#' Columns: `homologue_label` (`"PFOA + 1 CF2"`), `homologue_steps` (signed CF2
#' step distance to the seed), `homologue_series_id` (`kmd_series_id` when
#' present), `homologue_seed` (seed id, when carried), `homologue_ref_id`
#' (inferred compound `ref_id` / `Library_ID`), `homologue_smiles` (inferred
#' SMILES; stored prefixed because bare `smiles` is stripped from the MS1 tier)
#' and `homologue_level` (`0` seed / `1` inferred, from `evidence_source`).
#'
#' @param project An [NTSProject].
#' @param net_homologue_ann Enriched homologue network `list(nodes, edges)`
#'   (e.g. from [infer_homologue_structures]) or `NULL`.
#' @param write_back If `TRUE`, write the columns into `@ms1_annotations$network`
#'   and attach the updated project as attribute `"project"` on the returned
#'   table (reassign: `project <- attr(out, "project")`). Default `FALSE` only
#'   returns the table (no slot write).
#'
#' @return A `data.frame` with `name` plus the `homologue_*` columns.
#'
#' @seealso [merge_ms2_network_layers], [merge_annotation_layers],
#'   [infer_homologue_structures]
#' @export
merge_ms1_homologue_layers <- function(project,
                                       net_homologue_ann = NULL,
                                       write_back = FALSE) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject.", call. = FALSE)
  }
  cols <- c(
    "homologue_label", "homologue_steps", "homologue_series_id",
    "homologue_seed", "homologue_ref_id", "homologue_smiles", "homologue_level"
  )
  out <- data.frame(name = character(0), stringsAsFactors = FALSE)
  for (cl in cols) {
    out[[cl]] <- if (cl %in% c("homologue_steps", "homologue_level")) integer(0) else character(0)
  }
  if (!is.null(net_homologue_ann) && is.list(net_homologue_ann) &&
      !is.null(net_homologue_ann$nodes) && is.data.frame(net_homologue_ann$nodes) &&
      nrow(net_homologue_ann$nodes)) {
    nd <- as.data.frame(net_homologue_ann$nodes, stringsAsFactors = FALSE)
    nm <- if ("id" %in% names(nd)) {
      trimws(as.character(nd$id))
    } else if ("name" %in% names(nd)) {
      trimws(as.character(nd$name))
    } else {
      rep(NA_character_, nrow(nd))
    }
    lab <- .nts_net_col_first(nd, c("homologue_label"))
    stp <- .nts_net_col_first(nd, c("homologue_steps", "n_steps"), mode = "num")
    sid <- .nts_net_col_first(nd, c("kmd_series_id", "homologue_series_id"))
    sed <- .nts_net_col_first(nd, c("homologue_seed", "seed"))
    ref <- .nts_net_col_first(nd, c("ref_id", "Library_ID", "propagated_ID"))
    smi <- .nts_net_col_first(nd, c("SMILES", "smiles", "propagated_SMILES"))
    if ("evidence_source" %in% names(nd)) {
      es <- trimws(as.character(nd$evidence_source))
      lvl <- ifelse(!is.na(es) & es == "Seed", 0L,
             ifelse(!is.na(es) & es == "Inferred", 1L, NA_integer_))
    } else {
      lvl <- as.integer(.nts_net_col_first(
        nd, c("Propagation_Level", "propagation_level"), mode = "num"
      ))
    }
    ok <- !is.na(nm) & nzchar(nm) &
      ((!is.na(lab) & nzchar(lab)) | (!is.na(ref) & nzchar(ref)))
    if (any(ok)) {
      out <- data.frame(name = nm[ok], stringsAsFactors = FALSE)
      out$homologue_label <- lab[ok]
      out$homologue_steps <- as.integer(stp[ok])
      out$homologue_series_id <- sid[ok]
      out$homologue_seed <- sed[ok]
      out$homologue_ref_id <- ref[ok]
      out$homologue_smiles <- smi[ok]
      out$homologue_level <- as.integer(lvl[ok])
    }
  }
  out <- out[order(out$name), , drop = FALSE]
  rownames(out) <- NULL
  if (isTRUE(write_back) && nrow(out) && any(nzchar(out$name))) {
    project <- set_ms1_annotations(project, out, step = "network")
    attr(out, "project") <- project
  }
  out
}


#' Count annotated features per annotation layer
#'
#' Summarises how many unique features have a non-empty compound id in each
#' spectral / network layer passed to [merge_annotation_layers].
#'
#' @param project An [NTSProject].
#' @param net_frag_ann,net_nl_ann,net_md_ann,net_homologue_ann Optional annotated
#'   network lists (same as [merge_annotation_layers]).
#'
#' @return A `data.frame` with columns `method`, `n_features` (unique features
#'   with `ref_id`), `n_level0` (network seeds only), `n_propagated` (network
#'   `Propagation_Level > 0` only). Spectral rows have `NA` for level columns.
#'
#' @export
count_annotation_layers <- function(project,
                                    net_frag_ann = NULL,
                                    net_nl_ann = NULL,
                                    net_md_ann = NULL,
                                    net_homologue_ann = NULL) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject.", call. = FALSE)
  }

  .count_one <- function(long, label) {
    if (!nrow(long)) {
      return(data.frame(
        method = label,
        n_features = 0L,
        n_level0 = 0L,
        n_propagated = 0L,
        stringsAsFactors = FALSE
      ))
    }
    feat <- unique(long$feature)
    pl <- long$propagation_level
    data.frame(
      method = label,
      n_features = length(feat),
      n_level0 = length(unique(long$feature[!is.na(pl) & pl == 0L])),
      n_propagated = length(unique(long$feature[!is.na(pl) & pl > 0L])),
      stringsAsFactors = FALSE
    )
  }

  layers <- list(
    Spectral_DB = .nts_annotation_layer_from_db(project),
    Network_DF = .nts_annotation_layer_from_net(net_frag_ann, "Network_DF"),
    Network_NL = .nts_annotation_layer_from_net(net_nl_ann, "Network_NL"),
    Network_MD = .nts_annotation_layer_from_net(net_md_ann, "Network_MD"),
    Homologue = .nts_annotation_layer_from_net(net_homologue_ann, "Homologue")
  )
  out <- do.call(rbind, lapply(names(layers), function(nm) {
    row <- .count_one(layers[[nm]], nm)
    if (nm == "Spectral_DB") {
      row$n_level0 <- NA_integer_
      row$n_propagated <- NA_integer_
    }
    row
  }))
  rownames(out) <- NULL
  out
}
