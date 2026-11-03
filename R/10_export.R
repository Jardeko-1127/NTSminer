#' Null-coalescing operator (internal)
#' @noRd
`%||%` <- function(a, b) {
  if (missing(a) || is.null(a)) {
    return(b)
  }
  if (length(a) == 1L && is.na(a)) {
    return(b)
  }
  a
}


#' Numeric abundance columns for export (exclude metadata)
#' @noRd
.nts_abundance_colnames <- function(ft) {
  reserved <- c(
    .nts_feature_shape_cols,
    "name", "feature_id", "row_ID",
    "is_isf", "parent_mz", "homologous_series_id",
    "nominal_mass", "kmd"
  )
  nm <- names(ft)
  cand <- setdiff(nm, reserved)
  cand[vapply(ft[, cand, drop = FALSE], is.numeric, logical(1L))]
}


#' Export GNPS FBMN inputs (quant table + MGF)
#'
#' Writes `quantification_table.csv` and `spectra.mgf` in the layout the GNPS
#' Feature-Based Molecular Networking (FBMN) workflow expects:
#' the quant table starts with the literal columns `row ID`, `row m/z` and
#' `row retention time` (**minutes**, converted from the seconds held in
#' `feature_table$rt`) followed by one column per sample; each MGF entry carries
#' `FEATURE_ID=` (matching the quant table `row ID`) rather than `TITLE=`.
#' `spectra.mgf` is built from `@ms2data`. Each MS2 entry should be a `list`
#' with optional keys `feature_id` / `id` (matches `feature_table$name` if
#' present), `precursorMz` / `precursor`, optional `rt`, and `mz` + `intensity`
#' or two-column `spectrum` matrix.
#'
#' @param project An [NTSProject] with non-empty `@feature_table`.
#' @param out_dir Output directory (created if needed).
#'
#' @return Invisibly, a character vector of written file paths.
#' @export
export_to_gnps <- function(project, out_dir) {
  .nts_check(project, "project", .arg = "project")
  ft <- get_feature_table(project)
  .nts_check(ft, "non_empty", .arg = "feature_table")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  abund_cols <- .nts_abundance_colnames(ft)
  ft_export <- ft
  if (!"name" %in% names(ft_export)) {
    ft_export$row_ID <- sprintf("F%06d", seq_len(nrow(ft_export)))
  }
  ids <- as.character(if ("name" %in% names(ft)) ft_export[["name"]] else ft_export[["row_ID"]])
  if (!"mz" %in% names(ft_export)) {
    stop("feature_table must contain an `mz` column for GNPS export.", call. = FALSE)
  }
  # GNPS FBMN requires the literal column names `row ID`, `row m/z` and
  # `row retention time`, with retention time in **minutes**
  # (`feature_table$rt` is in seconds -- see 03_feature_extract.R).
  rt_min <- if ("rt" %in% names(ft_export)) {
    as.numeric(ft_export[["rt"]]) / 60
  } else {
    rep(NA_real_, nrow(ft_export))
  }
  qtab <- data.frame(
    ids, as.numeric(ft_export[["mz"]]), rt_min,
    stringsAsFactors = FALSE, check.names = FALSE
  )
  names(qtab) <- c("row ID", "row m/z", "row retention time")
  if (length(abund_cols)) {
    ab <- ft_export[, abund_cols, drop = FALSE]
    for (cn in abund_cols) {
      qtab[[cn]] <- ab[[cn]]
    }
  }
  qpath <- file.path(out_dir, "quantification_table.csv")
  write.csv(qtab, qpath, row.names = FALSE)

  sp <- methods::slot(project, "ms2data")
  mgf_path <- file.path(out_dir, "spectra.mgf")
  if (!length(sp)) {
    warning("ms2data is empty; wrote quantification table only.", call. = FALSE)
    return(invisible(qpath))
  }
  lines_out <- character(0)
  for (k in seq_along(sp)) {
    s <- sp[[k]]
    if (is.null(s$mz) && is.matrix(s$spectrum)) {
      s$mz <- s$spectrum[, 1L]
      s$intensity <- s$spectrum[, 2L]
    }
    pm <- s$precursorMz %||% s$precursor %||% NA_real_
    tit <- s$feature_id %||% s$id %||% paste0("sp_", k)
    rt <- s$rt %||% NA_real_
    mz <- s$mz
    int <- s$intensity
    if (length(mz) != length(int) || !length(mz)) {
      next
    }
    ord <- order(mz)
    mz <- mz[ord]
    int <- int[ord]
    blk <- c(
      "BEGIN IONS",
      # FBMN matches spectra to the quant table through `FEATURE_ID=` (the
      # value must equal the `row ID` column above), not `TITLE=`.
      paste0("FEATURE_ID=", tit),
      paste0("PEPMASS=", pm),
      if (is.finite(rt)) paste0("RTINSECONDS=", rt) else character(0),
      "MSLEVEL=2",
      "CHARGE=1+",
      sprintf("%.6f %.6f", mz, int),
      "END IONS"
    )
    lines_out <- c(lines_out, blk)
  }
  writeLines(lines_out, mgf_path)
  invisible(c(qpath, mgf_path))
}


#' Merge annotation scores onto features (first match column)
#' @noRd
.nts_annotation_scores <- function(ft, ann, score_column = NULL) {
  if (!nrow(ann)) {
    return(rep(NA_real_, nrow(ft)))
  }
  jn <- intersect(names(ft), names(ann))
  jn <- jn[jn %in% c("name", "feature_id", "ref_id", "mz")]
  if (!length(jn)) {
    return(rep(NA_real_, nrow(ft)))
  }
  key <- jn[1L]
  if (is.null(score_column)) {
    cn <- names(ann)
    hit <- cn[grepl("(score|total|ndp|similarity)", cn, ignore.case = TRUE)]
    score_column <- if (length(hit)) hit[1L] else NULL
  }
  if (is.null(score_column) || !score_column %in% names(ann)) {
    return(rep(NA_real_, nrow(ft)))
  }
  u <- unique(ann[, c(key, score_column), drop = FALSE])
  m <- merge(
    ft[, key, drop = FALSE],
    u,
    by = key,
    all.x = TRUE,
    sort = FALSE
  )
  vals <- m[[score_column]]
  as.numeric(vals)
}


#' Export SIRIUS `.ms` files for low-confidence, non-ISF features
#'
#' One file per feature: `compound_<id>.ms` with `>compound`, `>parentmass`,
#' `>ms2` blocks. Candidates: `is_isf` is `FALSE` or missing, and annotation
#' score is `NA` or **<=** `score_threshold`.
#'
#' @param project An [NTSProject] with `@ms2data` MS2 lists (named by feature
#'   id or index string).
#' @param out_dir Output directory.
#' @param score_threshold Maximum score still exported as low-confidence.
#' @param score_column Column in `@db_annotations`; auto-detected if `NULL`.
#'
#' @return Invisibly, paths to written `.ms` files.
#' @export
export_to_sirius <- function(project,
                             out_dir,
                             score_threshold = 0.3,
                             score_column = NULL) {
  .nts_check(project, "project", .arg = "project")
  ft <- get_feature_table(project)
  sp <- methods::slot(project, "ms2data")
  ann <- methods::slot(project, "db_annotations")
  .nts_check(ft, "non_empty", .arg = "feature_table")
  out_dir <- normalizePath(out_dir, mustWork = FALSE)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  scores <- .nts_annotation_scores(ft, ann, score_column)
  ms1 <- get_ms1_annotations(project)
  isf <- if (nrow(ms1) && "is_isf" %in% names(ms1)) {
    ii <- match(as.character(ft$name), as.character(ms1$name))
    v <- ms1$is_isf[ii]
    v[is.na(v)] <- FALSE
    as.logical(v)
  } else if ("is_isf" %in% names(ft)) {
    ft$is_isf
  } else {
    rep(FALSE, nrow(ft))
  }
  isf[is.na(isf)] <- FALSE
  low <- is.na(scores) | scores <= score_threshold
  take <- !isf & low
  if (!any(take)) {
    warning("No features matched export filters.", call. = FALSE)
    return(invisible(character(0)))
  }

  idvec <- if ("name" %in% names(ft)) {
    ft$name
  } else {
    sprintf("F%06d", seq_len(nrow(ft)))
  }
  paths <- character(0)
  for (i in which(take)) {
    sid <- as.character(idvec[i])
    sp_i <- sp[[sid]]
    if (is.null(sp_i)) {
      sp_i <- sp[[i]]
    }
    if (is.null(sp_i)) {
      next
    }
    if (is.null(sp_i$mz) && is.matrix(sp_i$spectrum)) {
      sp_i$mz <- sp_i$spectrum[, 1L]
      sp_i$intensity <- sp_i$spectrum[, 2L]
    }
    pm <- sp_i$precursorMz %||% sp_i$precursor %||% ft$mz[i]
    mz <- sp_i$mz
    int <- sp_i$intensity
    if (length(mz) != length(int) || !length(mz)) {
      next
    }
    fn <- file.path(out_dir, paste0("compound_", sid, ".ms"))
    ord <- order(mz)
    mz <- mz[ord]
    int <- int[ord]
    lines <- c(
      paste0(">compound ", sid),
      paste0(">parentmass ", pm),
      ">ms2",
      sprintf("%.6f\t%.6f", mz, int)
    )
    writeLines(lines, fn)
    paths <- c(paths, fn)
  }
  invisible(paths)
}

# ============================================================================
# Global results export (master wide table + network + audit log)
# ============================================================================


#' Export master wide table (quantitative + MS1 qualitative)
#'
#' Read-only: copies `@feature_table` and `@ms1_annotations` into local
#' `data.table` objects, **never** assigns into the [NTSProject] slots.
#' Left-joins qualitative columns onto every feature row on `name` (all
#' quantitative rows retained). Writes `master_results.csv` via
#' [data.table::fwrite].
#'
#' @param project An [NTSProject] with non-empty `@feature_table` and column
#'   `name`.
#' @param out_dir Output directory (created with `recursive = TRUE` if missing).
#'
#' @return Invisibly, the absolute path to `master_results.csv`.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' d <- tempfile("nts_em_ex"); dir.create(d)
#' p <- NTSminer::create_nts_project(d)
#' methods::slot(p, "feature_table") <- data.frame(
#'   name = c("F1", "F2"), mz = 1:2, rt = 1:2, s1 = 10:11, stringsAsFactors = FALSE
#' )
#' methods::slot(p, "ms1_annotations") <- data.frame(
#'   name = "F1", KMD = 0.1, stringsAsFactors = FALSE
#' )
#' od <- tempfile("nts_em_out"); dir.create(od)
#' mp <- NTSminer::export_master_results(p, out_dir = od)
#' stopifnot(file.exists(mp))
#' @export
export_master_results <- function(project, out_dir = "results") {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  ft <- get_feature_table(project)
  .nts_check(ft, "non_empty", "name_col", .arg = "feature_table")

  out_dir0 <- path.expand(out_dir)
  if (!dir.exists(out_dir0)) {
    dir.create(out_dir0, recursive = TRUE, showWarnings = FALSE)
  }
  out_dir_abs <- normalizePath(out_dir0, winslash = "/", mustWork = TRUE)

  ft_dt <- data.table::copy(data.table::as.data.table(ft))
  ms1_raw <- get_ms1_annotations(project)
  ms1 <- data.table::as.data.table(ms1_raw)

  if (!nrow(ms1)) {
    merged <- ft_dt
  } else if (!"name" %in% names(ms1)) {
    warning(
      "ms1_annotations lacks column `name`; exporting quantitative feature_table only.",
      call. = FALSE
    )
    merged <- ft_dt
  } else {
    ms1 <- data.table::copy(ms1)
    ## Drop EVERY column shared with ft_dt (not just mz/rt) before the join.
    ## Otherwise data.table right-joins the remaining overlaps as `i.`-prefixed
    ## duplicate columns (i.mzmin, i.rtmax, ...). ft_dt wins on overlap, as mz/rt
    ## already did.
    dupe_meta <- setdiff(intersect(names(ms1), names(ft_dt)), "name")
    if (length(dupe_meta)) {
      ms1[, (dupe_meta) := NULL]
    }
    merged <- ms1[ft_dt, on = "name"]
  }

  pc <- methods::slot(project, "physchem_proxy")
  if (is.data.frame(pc) && nrow(pc) && "name" %in% names(pc)) {
    merged <- .nts_left_join_add_cols(merged, data.table::as.data.table(pc))
  }

  master_path <- file.path(out_dir_abs, "master_results.csv")
  data.table::fwrite(merged, master_path, bom = TRUE, quote = "auto")
  invisible(master_path)
}


#' Export Cytoscape / Gephi-style network tables
#'
#' Writes `cytoscape_nodes.csv` and `cytoscape_edges.csv` from a pre-built
#' `network_data` object (e.g. [get_network_data] or [generate_network_data]).
#' **Does not** accept an [NTSProject] or recompute MS2 similarity.
#'
#' @param network_data Named `list` with `nodes` and `edges` (`data.frame` or
#'   `data.table`).
#' @param out_dir Output directory (created with `recursive = TRUE` if missing).
#'
#' @return Invisibly, a character vector: paths to nodes and edges CSV files.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' net <- list(
#'   nodes = data.frame(id = "a", mz = 100, stringsAsFactors = FALSE),
#'   edges = data.frame(from = "a", to = "a", stringsAsFactors = FALSE)
#' )
#' od <- tempfile("nts_cy_out"); dir.create(od)
#' paths <- NTSminer::export_cytoscape_network(net, out_dir = od)
#' stopifnot(length(paths) == 2L, all(file.exists(paths)))
#' @export
export_cytoscape_network <- function(network_data, out_dir = "results") {
  .nts_check(NULL, "data.table_ns")
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

  out_dir0 <- path.expand(out_dir)
  if (!dir.exists(out_dir0)) {
    dir.create(out_dir0, recursive = TRUE, showWarnings = FALSE)
  }
  out_dir_abs <- normalizePath(out_dir0, winslash = "/", mustWork = TRUE)

  ndt <- data.table::as.data.table(nd)
  edt <- data.table::as.data.table(ed)
  np <- file.path(out_dir_abs, "cytoscape_nodes.csv")
  ep <- file.path(out_dir_abs, "cytoscape_edges.csv")
  data.table::fwrite(ndt, np, bom = TRUE, quote = "auto")
  data.table::fwrite(edt, ep, bom = TRUE, quote = "auto")
  invisible(c(np, ep))
}


#' @noRd
.nts_format_log_entry <- function(entry) {
  if (!is.list(entry)) {
    return(as.character(entry)[1L])
  }
  ts <- entry[["timestamp"]]
  st <- entry[["step"]]
  pr <- entry[["params"]]
  ts_chr <- if (!inherits(ts, c("POSIXct", "POSIXt"))) {
    as.character(ts)[1L]
  } else {
    format(ts, usetz = TRUE)
  }
  st_chr <- as.character(st)[1L]
  if (is.null(st_chr) || !nzchar(st_chr)) {
    st_chr <- "<unknown_step>"
  }
  pr_str <- ""
  if (is.list(pr) && length(pr)) {
    nm <- names(pr)
    if (is.null(nm)) {
      nm <- rep("", length(pr))
    }
    parts <- vapply(seq_along(pr), function(i) {
      v <- pr[[i]]
      lab <- if (nzchar(nm[i])) paste0(nm[i], "=") else ""
      paste0(lab, paste(deparse(v, width.cutoff = 500L), collapse = " "))
    }, character(1L))
    pr_str <- paste(parts, collapse = "; ")
  }
  paste0("[", ts_chr, "] ", st_chr, if (nzchar(pr_str)) paste0(" | ", pr_str) else "")
}


#' @noRd
.nts_left_join_add_cols <- function(base_dt, right_dt, drop_meta = c("mz", "rt")) {
  if (is.null(right_dt) || !nrow(right_dt) || !"name" %in% names(right_dt)) {
    return(base_dt)
  }
  add <- setdiff(names(right_dt), c("name", drop_meta, names(base_dt)))
  if (!length(add)) {
    return(base_dt)
  }
  sub <- right_dt[, c("name", add), with = FALSE]
  sub[base_dt, on = "name"]
}


#' Export assembled project results for Excel / downstream stats
#'
#' Builds a **master feature table** by left-joining (on `name`)
#' `ms1_annotations` (base), `feature_table` (sample intensities),
#' best-ranked `db_annotations` row per feature, and `ms2_annotations`.
#' Duplicate metadata columns (`mz`, `rt`) from satellite tables are dropped
#' so `mz` / `rt` stay aligned with the MS1 layer. When [compute_ms2_qc()] has
#' run, its per-feature DDA quality columns (`spectral_entropy`,
#' `top_peak_fraction`, `inter_scan_similarity`, `match_score_max`,
#' `match_score_median`, `mean_purity`, `qc_grade`) are carried into the master
#' table as descriptive columns. Optional network tables and a text processing
#' log are also written.
#'
#' **Network tables caliber.** When the merged network is non-empty, the
#' function writes `Network_Nodes.csv` / `Network_Edges.csv` built with
#' `get_network_data(project, type = "all")` -- i.e. the **path-A** merged graph
#' (KMD homologue + transformation + ISF + MS2 fragment / neutral-loss /
#' mass-difference layers), with **default scale guards** (`top_k = Inf`,
#' `max_edges = Inf`, `min_degree = 0L`). The homologue edges here are the
#' path-A KMD edges: **directionally stored** (consecutive members within a KMD
#' series, `from` -> `to`) and **without** the `target_mz_diff` filter /
#' bidirectional augmentation applied by [build_homologue_network]. To obtain the
#' filtered, bidirectional homologue network (or the algebraic
#' [infer_homologue_structures] view), call `build_homologue_network()` /
#' `infer_homologue_structures()` explicitly. A project with no purified MS2
#' degrades to an MS1-only network and emits a warning; it does not fail.
#'
#' @param project An [NTSProject].
#' @param output_dir Directory to create (if needed) and write CSV/TXT into
#'   (default `"NTS_Results"`).
#'
#' @return Invisibly returns `project`.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' d <- tempfile("nts_epr"); dir.create(d)
#' p <- NTSminer::create_nts_project(d)
#' methods::slot(p, "feature_table") <- data.frame(
#'   name = "F1", mz = 100, rt = 1, s1 = 5, stringsAsFactors = FALSE
#' )
#' methods::slot(p, "ms1_annotations") <- data.frame(
#'   name = "F1", mz = 100, rt = 1, KMD = 0.1, stringsAsFactors = FALSE
#' )
#' od <- tempfile("nts_epr_out"); dir.create(od)
#' NTSminer::export_project_results(p, output_dir = od)
#' stopifnot(file.exists(file.path(od, "Master_Feature_Table.csv")))
#' @export
export_project_results <- function(project, output_dir = "NTS_Results") {
  .nts_check(project, "data.table_ns", "project", .arg = "project")

  out_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)

  ms1 <- data.table::as.data.table(get_ms1_annotations(project))
  if (!nrow(ms1) || !"name" %in% names(ms1)) {
    warning(
      "ms1_annotations is empty or lacks `name`; Master_Feature_Table.csv may be empty.",
      call. = FALSE
    )
  }

  ft <- data.table::as.data.table(get_feature_table(project))
  ft_j <- NULL
  if (nrow(ft) && "name" %in% names(ft)) {
    ft_nc <- setdiff(names(ft), c("name", intersect(c("mz", "rt"), names(ft))))
    ## Also drop the ft columns that ms1 already carries, so the join below does
    ## not emit `i.`-prefixed duplicates (i.mzmin, i.rtmax, ...). ms1 wins on
    ## overlap, matching the mz/rt handling above.
    ft_nc <- setdiff(ft_nc, names(ms1))
    if (length(ft_nc)) {
      ft_j <- ft[, c("name", ft_nc), with = FALSE]
    }
  }

  out <- if (!is.null(ft_j)) {
    ft_j[ms1, on = "name"]
  } else {
    data.table::copy(ms1)
  }

  db <- .nts_db_best_per_name(methods::slot(project, "db_annotations"))
  if (is.data.frame(db) && nrow(db)) {
    out <- .nts_left_join_add_cols(out, data.table::as.data.table(db))
  }

  m2 <- data.table::as.data.table(get_ms2_annotations(project))
  if (nrow(m2) && "name" %in% names(m2)) {
    out <- .nts_left_join_add_cols(out, m2)
  }

  master_path <- file.path(out_dir, "Master_Feature_Table.csv")
  data.table::fwrite(out, master_path, bom = TRUE)

  net_files <- character(0L)
  nd_full <- tryCatch(
    get_network_data(project, type = "all"),
    error = function(e) NULL
  )
  nd <- NULL
  ed <- NULL
  if (is.list(nd_full)) {
    nd <- nd_full$nodes
    ed <- nd_full$edges
  }
  nr <- get_network_results(project)
  if (is.list(nr)) {
    if (!is.data.frame(nd) || !nrow(nd)) {
      if (is.data.frame(nr$nodes) && nrow(nr$nodes)) {
        nd <- nr$nodes
      }
    }
    if (!is.data.frame(ed) || !nrow(ed)) {
      if (is.data.frame(nr$edges) && nrow(nr$edges)) {
        ed <- nr$edges
      }
    }
  }
  if (is.data.frame(nd) && nrow(nd)) {
    np <- file.path(out_dir, "Network_Nodes.csv")
    data.table::fwrite(as.data.frame(nd, stringsAsFactors = FALSE), np, bom = TRUE)
    net_files <- c(net_files, basename(np))
  }
  if (is.data.frame(ed) && nrow(ed)) {
    ep <- file.path(out_dir, "Network_Edges.csv")
    data.table::fwrite(as.data.frame(ed, stringsAsFactors = FALSE), ep, bom = TRUE)
    net_files <- c(net_files, basename(ep))
  }

  log <- methods::slot(project, "log")
  log_path <- file.path(out_dir, "Processing_Log.txt")
  if (length(log)) {
    lines <- vapply(log, .nts_format_log_entry, character(1L), USE.NAMES = FALSE)
  } else {
    lines <- "(no log entries)"
  }
  writeLines(lines, log_path, useBytes = FALSE)

  msg <- paste0(
    "\n",
    strrep("=", 60L),
    "\n  NTSminer - project results exported\n",
    strrep("=", 60L),
    "\n  Output directory: ", out_dir, "\n",
    "  - Master_Feature_Table.csv (", nrow(out), " rows x ", ncol(out), " columns)\n",
    "  - Processing_Log.txt\n"
  )
  if (length(net_files)) {
    msg <- paste0(msg, "  - ", paste(net_files, collapse = "\n  - "), "\n")
  }
  msg <- paste0(msg, strrep("=", 60L), "\n")
  message(msg)

  invisible(project)
}
