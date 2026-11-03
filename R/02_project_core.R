#' Create an empty [NTSProject]
#'
#' @param dir Path to a directory of raw MS files (must exist).
#'
#' @return A valid [NTSProject] instance.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' d <- tempfile("nts_ex"); dir.create(d)
#' p <- create_nts_project(d)
#' stopifnot(methods::is(p, "NTSProject"))
#' @export
create_nts_project <- function(dir) {
  dir <- normalizePath(dir, mustWork = TRUE)
  methods::new(
    "NTSProject",
    raw_data_path = dir,
    feature_table = data.frame(),
    ms1_annotations = list(),
    xdata = NULL,
    ms2data = list(),
    ms2_annotations = list(),
    clean_ms2 = list(),
    db_annotations = data.frame(),
    physchem_proxy = data.frame(),
    sirius_fingerprints = data.table::data.table(),
    network_results = list(),
    gc_spectra = list(),
    gc_deconvolved = list(),
    log = list()
  )
}

#' Deep copy an [NTSProject] (tabular isolation)
#'
#' Returns a new S4 instance whose `data.table`-backed or `data.frame` tabular
#' slots are duplicated with [data.table::copy] (or an equivalent deep copy for
#' plain `data.frame`s) so downstream `:=` pipelines cannot mutate the original
#' object by reference. Nested lists under `@ms2data` / `@clean_ms2` are walked
#' recursively to duplicate embedded tabular objects; `@xdata` is copied by
#' shallow reference when a lossless structural duplicate is not attempted here.
#'
#' @param project An [NTSProject].
#'
#' @return A new [NTSProject] with the same slot names as `project`.
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' d <- tempfile("nts_cp"); dir.create(d)
#' p <- create_nts_project(d)
#' p2 <- copy_project(p)
#' stopifnot(methods::is(p2, "NTSProject"))
#' @export
copy_project <- function(project) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  sn <- methods::slotNames(project)
  args <- vector("list", length(sn))
  names(args) <- sn
  for (nm in sn) {
    val <- if (identical(nm, "feature_table")) {
      get_feature_table(project)
    } else if (identical(nm, "network_results")) {
      get_network_results(project)
    } else {
      ## NOTE: for `ms1_annotations` do NOT use get_ms1_annotations(): with
      ## step = NULL it returns a *flattened data.frame*, which would replace the
      ## slot's canonical named list of sub-tables. Fall through to the raw slot
      ## so .nts_copy_project_slot() deep-copies the list (same as ms2_annotations).
      methods::slot(project, nm)
    }
    args[[nm]] <- .nts_copy_project_slot(val, nm)
  }
  do.call(methods::new, c(list(Class = "NTSProject"), args))
}

#' @keywords internal
#' @noRd
.nts_copy_tabular_df_like <- function(x) {
  if (!is.data.frame(x) && !data.table::is.data.table(x)) {
    stop("expected a data.frame or data.table", call. = FALSE)
  }
  dt <- data.table::copy(data.table::as.data.table(x))
  if (data.table::is.data.table(x)) {
    return(dt)
  }
  as.data.frame(dt, stringsAsFactors = FALSE)
}

#' @keywords internal
#' @noRd
.nts_copy_nested_r <- function(x, depth = 0L, max_depth = 32L) {
  if (!length(x)) {
    return(x)
  }
  if (depth > max_depth) {
    return(x)
  }
  if (inherits(x, "igraph")) {
    return(x)
  }
  if (isS4(x) && !data.table::is.data.table(x)) {
    return(x)
  }
  if (data.table::is.data.table(x)) {
    return(data.table::copy(x))
  }
  if (is.data.frame(x)) {
    return(as.data.frame(
      data.table::copy(data.table::as.data.table(x)),
      stringsAsFactors = FALSE
    ))
  }
  if (is.list(x)) {
    return(lapply(x, .nts_copy_nested_r, depth = depth + 1L, max_depth = max_depth))
  }
  x
}

#' @keywords internal
#' @noRd
.nts_copy_project_slot <- function(val, slot_name) {
  if (slot_name %in% c(
    "feature_table", "db_annotations", "physchem_proxy"
  )) {
    return(.nts_copy_tabular_df_like(val))
  }
  ## `ms1_annotations` / `ms2_annotations` are named lists of sub-tables;
  ## deep-copy every piece so a `data.table`-backed one cannot mutate by ref.
  if (slot_name %in% c(
    "ms1_annotations", "ms2_annotations", "ms2data", "clean_ms2", "gc_spectra",
    "gc_deconvolved", "network_results", "log"
  )) {
    return(.nts_copy_nested_r(val))
  }
  val
}

#' @keywords internal
#' @noRd
.nts_log_append <- function(log, step, params = list()) {
  c(
    log,
    list(list(
      timestamp = Sys.time(),
      step = step,
      params = params
    ))
  )
}

#' @keywords internal
#' @noRd
.nts_xcms_feature_meta_cols <- function() {
  c(
    "name", "mz", "rt",
    "mzmed", "mzmin", "mzmax",
    "rtmed", "rtmin", "rtmax",
    "npeaks", "peakidx", "row_ID", "row_id", "feature_id",
    "is_filled", "chrom_peak_id"
  )
}

#' @keywords internal
#' @noRd
.nts_feature_intensity_colnames <- function(ft) {
  if (!nrow(ft) || !length(names(ft))) {
    return(character(0L))
  }
  meta <- .nts_xcms_feature_meta_cols()
  setdiff(names(ft), meta)
}

#' Round per-sample intensity columns of a feature table to integers
#'
#' The stored quantitative core (`@feature_table`) keeps peak areas / response
#' intensities as whole numbers. Only the sample (abundance) columns are rounded;
#' `mz`, `rt` and any other meta column stay numeric. `NA` is preserved, and
#' values are kept as `double` when they would overflow the integer range.
#'
#' @param ft A `data.frame` / `data.table` shaped like `@feature_table`.
#' @return `ft` with integer abundance columns.
#' @keywords internal
#' @noRd
.nts_round_feature_intensities <- function(ft) {
  if (!is.data.frame(ft) || !nrow(ft) || !length(names(ft))) {
    return(ft)
  }
  int_cols <- .nts_feature_intensity_colnames(ft)
  int_cols <- int_cols[
    vapply(int_cols, function(nm) is.numeric(ft[[nm]]), logical(1L))
  ]
  if (!length(int_cols)) {
    return(ft)
  }
  for (nm in int_cols) {
    v <- ft[[nm]]
    if (is.integer(v)) {
      next
    }
    r <- suppressWarnings(round(as.numeric(v)))
    ft[[nm]] <- if (any(abs(r) > .Machine$integer.max, na.rm = TRUE)) {
      r
    } else {
      suppressWarnings(as.integer(r))
    }
  }
  ft
}

#' @keywords internal
#' @noRd
.nts_feature_table_core <- function(ft) {
  if (!nrow(ft)) {
    return(ft)
  }
  if (!"name" %in% names(ft)) {
    ft$name <- sprintf("FT%05d", seq_len(nrow(ft)))
  }
  if (!all(c("mz", "rt") %in% names(ft))) {
    stop("feature_table must contain columns mz and rt", call. = FALSE)
  }
  ann_strip <- c(
    "annotation_level",
    "kmd_series_id",
    "KM",
    "NKM",
    "KMD",
    "standard_md",
    "mz_remainder",
    "Kendrick_Mass",
    "Nominal_KM"
  )
  bad <- names(ft)[
    names(ft) %in% ann_strip | grepl("^kmd_", names(ft), ignore.case = TRUE)
  ]
  if (length(bad)) {
    ft <- ft[, setdiff(names(ft), bad), drop = FALSE]
  }
  ## drop xcms / data.table artefacts and non-quantitative channel tags
  noise <- intersect(
    c("X1", "ms_level"),
    names(ft)
  )
  noise <- c(noise, grep("^X[0-9]+$", names(ft), value = TRUE))
  noise <- unique(noise)
  if (length(noise)) {
    ft <- ft[, setdiff(names(ft), noise), drop = FALSE]
  }
  ab <- .nts_feature_intensity_colnames(ft)
  ab <- ab[vapply(ab, function(nm) is.numeric(ft[[nm]]), logical(1L))]
  ft[, unique(c("name", "mz", "rt", ab)), drop = FALSE]
}

#' Column-to-step map for the `@ms1_annotations` named list
#'
#' `@ms1_annotations` is stored as a **named list** of sub-tables, one per
#' processing step, rather than a single wide `data.frame`. This function is
#' the single source of truth for which column belongs to which step. Every
#' sub-table keeps its own `name` key, so the pieces re-join losslessly by
#' feature and stay aligned with `@feature_table$name`.
#'
#' @return An ordered named `list` of character vectors:
#'   \describe{
#'     \item{base}{Feature identity and bounds (`name`, `mzmin`, `mzmax`,
#'       `rtmin`, `rtmax`, `mz`, `rt`).}
#'     \item{eic}{EIC shape metrics from [extract_features] (`.nts_ms1_eic_shape_colnames()`).}
#'     \item{library}{Reference-library identity columns (`library_*`).}
#'     \item{isotopes}{[flag_isotopes] output.}
#'     \item{isf}{[flag_isf] / [recognize_isf] output.}
#'     \item{kmd}{[calculate_kmd] output (KM / NKM / KMD / homologue series).}
#'     \item{adduct}{Adduct assignment (`adduct_type`, `neutral_mass`).}
#'     \item{master}{`annotation_master` merge output (`master_*`).}
#'     \item{network}{Structural / homologue propagation (`nap_structure_*`,
#'       `inferred_formula`, `inferred_name`).}
#'     \item{annotation}{Final confidence and RT-prediction columns.}
#'     \item{extra}{Any column that matches no step above. Never dropped.}
#'   }
#'
#' This function is the **full** column-to-step schema map. The stored slot,
#' however, only keeps the steps that actually carry data (see
#' `.nts_ms1_partition`), so `names(@ms1_annotations)` lists just the
#' populated steps rather than all eleven.
#' @keywords internal
#' @noRd
.nts_ms1_groups <- function() {
  list(
    base = c("name", "mzmin", "mzmax", "rtmin", "rtmax", "mz", "rt"),
    eic = .nts_ms1_eic_shape_colnames(),
    library = c("library_Compound_name", "library_ref_id", "library_SMILES"),
    isotopes = c(
      "is_isotope", "parent_name", "isotope_shift_da", "isotope_kind",
      "isotope_label", "isotope_c13_n", "isotope_formula", "isotope_cosine"
    ),
    isf = c(
      "is_isf", "parent_mz", "isf_annotation", "isf_candidate_ids",
      "isf_candidate_names", "isf_candidate_n_fragments"
    ),
    kmd = c(
      "standard_md", "KM", "NKM", "KMD", "mz_remainder", "kmd_series_id",
      "Kendrick_Mass", "Nominal_KM"
    ),
    adduct = c("adduct_type", "neutral_mass"),
    master = c("master_smiles", "master_score", "master_propagation_level"),
    network = c(
      "inferred_formula", "inferred_name", "nap_structure_smiles",
      "nap_structure_formula", "nap_structure_seed", "nap_structure_tanimoto",
      "nap_structure_depth",
      "homologue_label", "homologue_steps", "homologue_series_id",
      "homologue_seed", "homologue_ref_id", "homologue_smiles",
      "homologue_level"
    ),
    annotation = c(
      "annotation_level", "confidence_level", "rt_confidence", "predicted_rt"
    ),
    extra = character(0)
  )
}

#' Sub-table names carried by `@ms2_annotations`
#'
#' `@ms2_annotations` is a named `list` mirroring `@ms1_annotations`: at most two
#' sub-tables, `annotation` (the per-feature MS2 evidence table) and
#' `score_grid` (the full per-cell MS2 score grid from
#' [annotate_by_spectral_match]). Only populated sub-tables are stored.
#' @keywords internal
#' @noRd
.nts_ms2_steps <- function() c("annotation", "score_grid")

#' Canonicalise an `@ms2_annotations` value into its stored named-list form
#'
#' Accepts a flat `data.frame` (legacy slot shape, treated as `annotation`), a
#' canonical list, or any list carrying `annotation` / `score_grid` elements.
#' Sub-tables with no rows are dropped so `names()` lists only the populated
#' pieces, mirroring `.nts_ms1_partition()`.
#' @keywords internal
#' @noRd
.nts_ms2_sanitize <- function(value) {
  out <- list()
  if (is.data.frame(value)) {
    if (nrow(value)) {
      out$annotation <- as.data.frame(value, stringsAsFactors = FALSE)
    }
    return(out)
  }
  if (is.list(value)) {
    a <- value[["annotation"]]
    g <- value[["score_grid"]]
    if (is.data.frame(a) && nrow(a)) {
      out$annotation <- as.data.frame(a, stringsAsFactors = FALSE)
    }
    if (is.data.frame(g) && nrow(g)) {
      out$score_grid <- as.data.frame(g, stringsAsFactors = FALSE)
    }
  }
  out
}

#' Read `@ms2_annotations` in its canonical named-list form
#'
#' Handles both the current named-list slot and the pre-refactor layout (a flat
#' `data.frame` in `@ms2_annotations` plus a separate `@ms2_score_grid` slot), so
#' projects saved before the merge keep reading correctly. The stray
#' `ms2_score_grid` attribute (if any) is folded into the list as `score_grid`.
#' @param x An [NTSProject].
#' @return A named `list` with any of `annotation` / `score_grid`.
#' @keywords internal
#' @noRd
.nts_ms2_read_list <- function(x) {
  raw <- tryCatch(methods::slot(x, "ms2_annotations"), error = function(e) NULL)
  out <- .nts_ms2_sanitize(raw)
  if (!length(out) && is.list(raw) && length(raw) &&
      !any(vapply(raw, is.data.frame, logical(1L))) &&
      all(vapply(raw, function(z) is.atomic(z) && !is.list(z), logical(1L))) &&
      length(unique(vapply(raw, length, integer(1L)))) == 1L &&
      "name" %in% names(raw)) {
    ## Legacy flat `data.frame` coerced to a bare list by S4 slot assignment:
    ## every element is a same-length atomic vector and `name` is a key.
    out <- .nts_ms2_sanitize(as.data.frame(raw, stringsAsFactors = FALSE))
  }
  if (is.null(out$score_grid)) {
    legacy <- attr(x, "ms2_score_grid")
    if (is.data.frame(legacy) && nrow(legacy)) {
      out$score_grid <- as.data.frame(legacy, stringsAsFactors = FALSE)
    }
  }
  out
}

#' Does one `@ms1_annotations` sub-table actually carry data?
#'
#' A step is worth storing only when it holds at least one column beyond the
#' shared `name` key **and** at least one of those columns has a non-`NA` value.
#' This is what lets `.nts_ms1_partition` drop steps that exist only as empty
#' scaffolding - e.g. `adduct_type` / `neutral_mass` are created as all-`NA` by
#' `.nts_ms1_sanitize_and_order` but stay `NA` until adduct assignment actually
#' runs, so the `adduct` step should not appear in the slot until then.
#'
#' @param z A candidate sub-table (`data.frame`).
#' @return `TRUE` when `z` carries data, else `FALSE`.
#' @keywords internal
#' @noRd
.nts_ms1_step_has_data <- function(z) {
  if (!is.data.frame(z) || !nrow(z) || !ncol(z)) {
    return(FALSE)
  }
  cols <- setdiff(names(z), "name")
  if (!length(cols)) {
    return(FALSE)
  }
  any(vapply(
    cols,
    function(cl) {
      v <- z[[cl]]
      is.atomic(v) && any(!is.na(v))
    },
    logical(1L)
  ))
}

#' Split a flat `@ms1_annotations` view into the named-list slot form
#'
#' Every data-carrying sub-table keeps a `name` key (added when absent) so that
#' `.nts_ms1_flatten` can re-join the pieces without loss. A step is stored
#' **only when it carries data** (see `.nts_ms1_step_has_data`): steps whose
#' mapped columns are absent or entirely `NA` are dropped, so
#' `names(@ms1_annotations)` mirrors the steps that were actually populated
#' instead of the full schema. Columns that match no step are parked in `extra`
#' (kept whenever it holds data). Dropping a step loses nothing recoverable -
#' the canonical flat columns are re-created by `.nts_ms1_sanitize_and_order` on
#' the next write, and [get_ms1_annotations] returns the same canonical flat
#' view.
#'
#' @param flat A `data.frame` shaped like the flat `@ms1_annotations` view.
#' @return A named `list` of sub-tables, one element per **data-carrying**
#'   `.nts_ms1_groups()` entry (empty `list` when `flat` has no rows).
#' @keywords internal
#' @noRd
.nts_ms1_partition <- function(flat) {
  gm <- .nts_ms1_groups()
  if (!is.data.frame(flat) || !nrow(flat) || !length(names(flat))) {
    return(list())
  }
  nm <- names(flat)
  assigned <- unique(unlist(gm, use.names = FALSE))
  extra_cols <- setdiff(nm, assigned)
  out <- vector("list", length(gm))
  names(out) <- names(gm)
  for (g in names(gm)) {
    cols <- if (identical(g, "extra")) extra_cols else gm[[g]]
    cols <- intersect(cols, nm)
    if (!length(cols)) {
      out[[g]] <- data.frame()
      next
    }
    keep <- if ("name" %in% cols) cols else unique(c("name", cols))
    out[[g]] <- flat[, keep, drop = FALSE]
  }
  ## Keep only the steps that actually carry data, so the slot's names mirror
  ## the populated processing steps rather than the full schema.
  keep <- vapply(out, .nts_ms1_step_has_data, logical(1L))
  out[keep]
}

#' Re-join the `@ms1_annotations` named list into a flat `data.frame`
#'
#' Left-joins every sub-table onto `base` (or the first sub-table carrying
#' `name`), preserving base row order and column order:
#' `base`, `eic`, `library`, `isotopes`, `isf`, `kmd`, `adduct`, `master`,
#' `network`, `annotation`, `extra`. Only the steps actually present in the slot
#' contribute columns, so a slot that stores just the populated steps (see
#' `.nts_ms1_partition`) flattens cleanly; missing steps simply add nothing.
#'
#' Accepts three shapes: the canonical named `list` of sub-tables, a plain
#' `data.frame` (legacy flat slot), and the *mangled* bare list produced when a
#' `data.frame` is pushed straight into the `"list"` slot with `@<-` or
#' [methods::slot]`<-` (S4 coerces it to a list of equal-length columns). The
#' last case is recovered so pre-refactor scripts keep reading correctly.
#'
#' @param x The `@ms1_annotations` slot.
#' @return A `data.frame`; 0-row if `x` holds no data.
#' @keywords internal
#' @noRd
.nts_ms1_flatten <- function(x) {
  if (is.data.frame(x)) {
    return(x)
  }
  if (!is.list(x) || !length(x)) {
    return(data.frame())
  }
  ok <- vapply(
    x,
    function(z) is.data.frame(z) && nrow(z) && "name" %in% names(z),
    logical(1L)
  )
  tabs <- x[ok]
  if (!length(tabs)) {
    ## Flat `data.frame` coerced to a bare list by S4 slot assignment:
    ## every element is a same-length atomic vector and `name` is a key.
    bare <- !any(vapply(x, is.data.frame, logical(1L))) &&
      all(vapply(x, function(z) is.atomic(z) && !is.list(z), logical(1L))) &&
      length(unique(vapply(x, length, integer(1L)))) == 1L
    if (bare && "name" %in% names(x)) {
      return(as.data.frame(x, stringsAsFactors = FALSE))
    }
    return(data.frame())
  }
  gm <- names(.nts_ms1_groups())
  nm <- names(tabs)
  nm <- c(intersect(gm, nm), setdiff(nm, gm))
  tabs <- tabs[nm]
  out <- as.data.frame(tabs[[1L]], stringsAsFactors = FALSE)
  for (k in nm[-1L]) {
    t <- as.data.frame(tabs[[k]], stringsAsFactors = FALSE)
    add <- setdiff(names(t), names(out))
    if (!length(add)) {
      next
    }
    pos <- match(as.character(out$name), as.character(t$name))
    for (cl in add) {
      out[[cl]] <- t[[cl]][pos]
    }
  }
  out
}

#' Column names that must never live on `@ms1_annotations` (legacy / wrong tier).
#' @keywords internal
#' @noRd
.nts_ms1_forbidden_names <- function() {
  c(
    "ALogP_proxy", "TPSA", "MolecularWeight", "HBondDonorCount", "HBondAcceptorCount",
    "formula", "SMILES", "smiles", "inferred_smiles",
    "score_cosine_ms2", "f_mz", "r_mz"
  )
}

#' Canonical (flat) `@ms1_annotations` column order: bounds before mz/rt,
#' EIC shape metrics, other columns alphabetically, flag / KMD tail.
#' @keywords internal
#' @noRd
.nts_ms1_sanitize_and_order <- function(ms1) {
  if (!is.data.frame(ms1)) {
    ## Accept the slot's named-list form (or any list of sub-tables) and fall
    ## back to the flat view; flat output is what every legacy caller expects.
    ms1 <- .nts_ms1_flatten(ms1)
  }
  if (!is.data.frame(ms1) || !nrow(ms1)) {
    return(ms1)
  }
  forb <- intersect(.nts_ms1_forbidden_names(), names(ms1))
  if (length(forb)) {
    ms1[forb] <- NULL
  }
  head_c <- c("name", "mzmin", "mzmax", "rtmin", "rtmax", "mz", "rt")
  tail_c <- c(
    "is_isotope", "parent_name",
    "isotope_shift_da", "isotope_kind", "isotope_label", "isotope_c13_n",
    "is_isf", "parent_mz", "isf_annotation",
    "isf_candidate_ids", "isf_candidate_names",
    "isf_candidate_n_fragments",
    "adduct_type", "neutral_mass",
    "standard_md", "KM", "NKM", "KMD", "mz_remainder", "kmd_series_id",
    "inferred_formula", "inferred_name", "annotation_level"
  )
  for (cl in head_c) {
    if (!cl %in% names(ms1)) {
      if (cl == "name") {
        ms1[[cl]] <- character(nrow(ms1))
      } else {
        ms1[[cl]] <- rep(NA_real_, nrow(ms1))
      }
    }
  }
  for (cl in tail_c) {
    if (!cl %in% names(ms1)) {
      if (cl %in% c("is_isotope", "is_isf")) {
        ms1[[cl]] <- rep(FALSE, nrow(ms1))
      } else if (cl %in% c(
        "parent_mz", "standard_md", "KM", "NKM", "KMD", "mz_remainder",
        "isotope_shift_da", "neutral_mass"
      )) {
        ms1[[cl]] <- rep(NA_real_, nrow(ms1))
      } else if (cl == "isotope_c13_n") {
        ms1[[cl]] <- rep(NA_integer_, nrow(ms1))
      } else {
        ms1[[cl]] <- rep(NA_character_, nrow(ms1))
      }
    }
  }
  eic_c <- .nts_ms1_eic_shape_colnames()
  rest <- setdiff(names(ms1), c(head_c, tail_c, eic_c))
  rest <- rest[order(rest)]
  ord <- c(
    intersect(head_c, names(ms1)),
    intersect(eic_c, names(ms1)),
    rest,
    intersect(tail_c, names(ms1))
  )
  ord <- ord[!duplicated(ord)]
  ms1[, ord, drop = FALSE]
}

#' @keywords internal
#' @noRd
.nts_init_ms1_from_ft <- function(ft, bounds_df = NULL) {
  if (!nrow(ft)) {
    return(data.frame())
  }
  ft <- .nts_feature_table_core(ft)
  n <- nrow(ft)
  mzmin <- rep(NA_real_, n)
  mzmax <- rep(NA_real_, n)
  rtmin <- rep(NA_real_, n)
  rtmax <- rep(NA_real_, n)
  if (is.data.frame(bounds_df) && nrow(bounds_df) == n) {
    if ("mzmin" %in% names(bounds_df)) {
      mzmin <- suppressWarnings(as.numeric(bounds_df$mzmin))
    }
    if ("mzmax" %in% names(bounds_df)) {
      mzmax <- suppressWarnings(as.numeric(bounds_df$mzmax))
    }
    if ("rtmin" %in% names(bounds_df)) {
      rtmin <- suppressWarnings(as.numeric(bounds_df$rtmin))
    }
    if ("rtmax" %in% names(bounds_df)) {
      rtmax <- suppressWarnings(as.numeric(bounds_df$rtmax))
    }
  }
  base <- data.frame(
    name = ft$name,
    mzmin = mzmin,
    mzmax = mzmax,
    rtmin = rtmin,
    rtmax = rtmax,
    mz = ft$mz,
    rt = ft$rt,
    is_isotope = rep(FALSE, n),
    parent_name = rep(NA_character_, n),
    isotope_shift_da = rep(NA_real_, n),
    isotope_kind = rep(NA_character_, n),
    isotope_label = rep(NA_character_, n),
    isotope_c13_n = rep(NA_integer_, n),
    is_isf = rep(FALSE, n),
    parent_mz = rep(NA_real_, n),
    isf_annotation = rep(NA_character_, n),
    isf_candidate_ids = rep(NA_character_, n),
    isf_candidate_names = rep(NA_character_, n),
    isf_candidate_n_fragments = rep(NA_character_, n),
    adduct_type = rep(NA_character_, n),
    neutral_mass = rep(NA_real_, n),
    standard_md = rep(NA_real_, n),
    KM = rep(NA_real_, n),
    NKM = rep(NA_real_, n),
    KMD = rep(NA_real_, n),
    mz_remainder = rep(NA_real_, n),
    kmd_series_id = rep(NA_character_, n),
    inferred_formula = rep(NA_character_, n),
    annotation_level = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
  base <- base[!duplicated(as.character(base$name)), , drop = FALSE]
  .nts_ms1_sanitize_and_order(base)
}

#' @keywords internal
#' @noRd
.nts_db_best_per_name <- function(db) {
  if (!is.data.frame(db) || !nrow(db) || !"name" %in% names(db)) {
    return(data.frame())
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    u <- db
    if ("score_total" %in% names(u)) {
      u <- u[order(u$score_total, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    }
    u <- u[!duplicated(as.character(u$name)), , drop = FALSE]
    rownames(u) <- NULL
    return(u)
  }
  dt <- data.table::as.data.table(db)
  if ("score_total" %in% names(dt)) {
    data.table::setorder(dt, -score_total, na.last = TRUE)
  }
  out <- unique(dt, by = "name")
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Safe `db_annotations` data.frame from an [NTSProject] or a list-like object.
#' @keywords internal
#' @noRd
.nts_get_project_db_annotations <- function(project) {
  db_raw <- NULL
  if (methods::is(project, "NTSProject")) {
    db_raw <- methods::slot(project, "db_annotations")
  } else if (is.list(project) && !is.null(project$db_annotations)) {
    db_raw <- project$db_annotations
  }
  if (is.null(db_raw)) {
    return(data.frame())
  }
  if (!is.data.frame(db_raw)) {
    db_raw <- tryCatch(
      as.data.frame(db_raw, stringsAsFactors = FALSE),
      error = function(e) data.frame()
    )
  }
  db_raw
}

#' @keywords internal
#' @noRd
.nts_ensure_ms2_ann_for_ms2_list <- function(project, ms2_list) {
  if (!length(ms2_list)) {
    return(project)
  }
  nms <- names(ms2_list)
  if (is.null(nms)) {
    return(project)
  }
  ft <- get_feature_table(project)
  cur <- get_ms2_annotations(project)
  existing <- if (nrow(cur) && "name" %in% names(cur)) {
    as.character(cur$name)
  } else {
    character(0L)
  }
  need <- nms[!nms %in% existing]
  if (!length(need)) {
    return(project)
  }
  if (!nrow(ft) || !all(c("name", "mz", "rt") %in% names(ft))) {
    stop("feature_table must contain name, mz, rt to expand ms2_annotations", call. = FALSE)
  }
  ft_names <- as.character(ft$name)
  matched <- need[need %in% ft_names]
  unmapped <- setdiff(need, ft_names)
  if (length(unmapped)) {
    warning(
      "ms2_annotations: ", length(unmapped),
      " feature(s) not found in feature_table and will be skipped: ",
      paste(utils::head(unmapped, 10), collapse = ", "),
      if (length(unmapped) > 10) " ..." else "",
      call. = FALSE
    )
  }
  if (!length(matched)) {
    return(project)
  }
  rows <- lapply(matched, function(nm) {
    m <- match(nm, ft_names)
    data.frame(
      name = nm,
      mz = ft$mz[m],
      rt = ft$rt[m],
      stringsAsFactors = FALSE
    )
  })
  add <- do.call(rbind, rows)
  if (nrow(cur)) {
    for (cn in setdiff(names(cur), names(add))) {
      add[[cn]] <- NA
    }
    for (cn in setdiff(names(add), names(cur))) {
      cur[[cn]] <- rep(NA, nrow(cur))
    }
    cur <- rbind(cur, add[, names(cur), drop = FALSE])
  } else {
    cur <- add
  }
  project <- set_ms2_annotations(project, cur)
  project
}

#' Build `@ms2_annotations` from `@ms1_annotations` + `feature_table` (mz/rt
#' fallback) + raw MS2 peak counts.  Features whose names cannot be
#' resolved to finite mz/rt in either source are **excluded with a warning**
#' rather than stopping the pipeline.
#' Build `@ms2_annotations` from `@ms1_annotations` + `feature_table` (mz/rt
#' fallback) + raw MS2 peak counts.  Features whose names cannot be
#' resolved to finite mz/rt in either source are **excluded with a warning**
#' rather than stopping the pipeline.
#' Build `@ms2_annotations` from `@ms1_annotations` + `feature_table` (mz/rt
#' fallback) + raw MS2 peak counts.  Features whose names cannot be
#' resolved to finite mz/rt in either source are **excluded with a warning**
#' rather than stopping the pipeline.
#' Build `@ms2_annotations` from `@ms1_annotations` + `feature_table` (mz/rt
#' fallback) + raw MS2 peak counts.  Features whose names cannot be
#' resolved to finite mz/rt in either source are **excluded with a warning**
#' rather than stopping the pipeline.
#' @keywords internal
#' @noRd
.nts_build_ms2_annotations_from_ms2data <- function(project, ms2_named_list) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  if (!length(ms2_named_list)) {
    project <- set_ms2_annotations(project, data.frame())
    return(project)
  }
  nms <- names(ms2_named_list)
  if (is.null(nms)) {
    stop("ms2data list must be named by feature id", call. = FALSE)
  }
  n_frag <- vapply(ms2_named_list, function(lst) {
    if (!length(lst)) {
      return(0L)
    }
    as.integer(sum(vapply(lst, function(s) {
      if (!is.list(s) || is.null(s$mz)) {
        return(0L)
      }
      length(s$mz)
    }, integer(1L))))
  }, integer(1L))
  seed <- data.table::data.table(name = nms, n_fragments = n_frag)

  ms1 <- data.table::as.data.table(get_ms1_annotations(project))
  if (nrow(ms1) && "name" %in% names(ms1)) {
    ms1[, name := as.character(name)]
    cols_keep <- intersect(c("name", "mz", "rt"), names(ms1))
    oms1 <- unique(ms1[, cols_keep, with = FALSE], by = "name")
    for (cn in c("mz", "rt")) {
      if (!cn %in% names(oms1)) {
        oms1[, (cn) := NA_real_]
      }
    }
    ann <- oms1[seed, on = "name"]
  } else {
    ann <- seed
    ann[, `:=`(mz = NA_real_, rt = NA_real_)]
  }

  ft <- data.table::as.data.table(get_feature_table(project))
  if (nrow(ft) && all(c("name", "mz", "rt") %in% names(ft))) {
    ftu <- unique(
      ft[, .(name = as.character(name), mz_ft = as.numeric(mz), rt_ft = as.numeric(rt))],
      by = "name"
    )
    ann <- ftu[, .(name, mz_ft, rt_ft)][ann, on = "name"]
    ann[is.na(mz) | !is.finite(mz), mz := mz_ft]
    ann[is.na(rt) | !is.finite(rt), rt := rt_ft]
    ann[, `:=`(mz_ft = NULL, rt_ft = NULL)]
  }

  missing_idx <- which(!is.finite(ann$mz) | !is.finite(ann$rt))
  if (length(missing_idx)) {
    missing_names <- ann$name[missing_idx]
    warning(
      "extract_ms2: ", length(missing_names),
      " MS2 feature(s) have no mz/rt in ms1_annotations or feature_table ",
      "and will be excluded from @ms2_annotations: ",
      paste(utils::head(unique(missing_names), 10), collapse = ", "),
      if (length(missing_names) > 10) " ..." else "",
      call. = FALSE
    )
    ann <- ann[-missing_idx]
  }

  project <- set_ms2_annotations(project, as.data.frame(ann, stringsAsFactors = FALSE))
  project
}

#' Update `n_fragments` (and CE summary) in `@ms2_annotations` from
#' consensus `@clean_ms2` spectra.
#' @keywords internal
#' @noRd
.nts_update_ms2_annotations_n_fragments <- function(project) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  cl <- methods::slot(project, "clean_ms2")
  if (!length(cl) || is.null(names(cl))) {
    return(project)
  }
  m2 <- data.table::as.data.table(get_ms2_annotations(project))
  if (!nrow(m2) || !"name" %in% names(m2)) {
    return(project)
  }
  m2[, name := as.character(name)]
  nmap <- vapply(names(cl), function(nm) {
    sp <- cl[[nm]]$spectrum
    if (is.data.frame(sp)) {
      as.integer(nrow(sp))
    } else {
      0L
    }
  }, integer(1L))
  ce_levels_map <- lapply(cl, function(el) {
    suppressWarnings(as.numeric(el$ce_levels))
  })
  n_ce_vec <- vapply(ce_levels_map, function(v) {
    sum(is.finite(v))
  }, integer(1L))
  ce_min_vec <- vapply(ce_levels_map, function(v) {
    v2 <- v[is.finite(v)]
    if (length(v2)) min(v2) else NA_real_
  }, numeric(1L))
  ce_max_vec <- vapply(ce_levels_map, function(v) {
    v2 <- v[is.finite(v)]
    if (length(v2)) max(v2) else NA_real_
  }, numeric(1L))
  idx <- match(m2$name, names(nmap))
  ok <- !is.na(idx)
  if (!"n_fragments" %in% names(m2)) {
    m2[, n_fragments := NA_integer_]
  }
  if (!"n_ce_levels" %in% names(m2)) {
    m2[, n_ce_levels := NA_integer_]
  }
  if (!"ce_min" %in% names(m2)) {
    m2[, ce_min := NA_real_]
  }
  if (!"ce_max" %in% names(m2)) {
    m2[, ce_max := NA_real_]
  }
  m2[ok, n_fragments := nmap[idx[ok]]]
  m2[ok, n_ce_levels := n_ce_vec[idx[ok]]]
  m2[ok, ce_min := ce_min_vec[idx[ok]]]
  m2[ok, ce_max := ce_max_vec[idx[ok]]]
  project <- set_ms2_annotations(project, as.data.frame(m2, stringsAsFactors = FALSE))
  .nts_ms2_strip_mz_bounds(project)
}

#' Remove legacy m/z interval columns from `@ms2_annotations` (bounds live on MS1).
#' @keywords internal
#' @noRd
.nts_ms2_strip_mz_bounds <- function(project) {
  m2 <- get_ms2_annotations(project)
  if (!is.data.frame(m2) || !nrow(m2)) {
    return(project)
  }
  drop <- intersect(c("mzmin", "mzmax"), names(m2))
  if (length(drop)) {
    m2 <- m2[, setdiff(names(m2), drop), drop = FALSE]
    project <- set_ms2_annotations(project, m2)
  }
  project
}

#' Merge best `db_annotations` library fields into `@ms2_annotations` by `name`.
#' @keywords internal
#' @noRd
.nts_merge_db_top1_into_ms2_annotations <- function(project) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  db <- methods::slot(project, "db_annotations")
  m2 <- get_ms2_annotations(project)
  if (!is.data.frame(m2) || !nrow(m2) || !"name" %in% names(m2)) {
    return(project)
  }
  if (!is.data.frame(db) || !nrow(db) || !"name" %in% names(db)) {
    return(project)
  }
  best <- .nts_db_best_per_name(db)
  want <- c(
    "name", "ref_id", "Compound.name", "SMILES",
    "matched_nl", "matched_md", "matched_diagnostic"
  )
  cols <- intersect(want, names(best))
  if (length(cols) < 2L) {
    return(project)
  }
  bdt <- data.table::as.data.table(best[, cols, drop = FALSE])
  m2dt <- data.table::as.data.table(m2)
  for (cn in setdiff(names(bdt), "name")) {
    if (cn %in% names(m2dt)) {
      m2dt[, (cn) := NULL]
    }
  }
  out <- bdt[m2dt, on = "name"]
  project <- set_ms2_annotations(project, as.data.frame(out, stringsAsFactors = FALSE))
  project
}

#' Drop features with no usable MS2 scans (empty list or only zero-peak scans)
#' @keywords internal
#' @noRd
.nts_ms2data_drop_empty <- function(L) {
  if (!length(L)) {
    return(L)
  }
  nms <- names(L)
  if (is.null(nms)) {
    nms <- as.character(seq_along(L))
    names(L) <- nms
  }
  keep <- vapply(L, function(el) {
    if (!length(el)) {
      return(FALSE)
    }
    any(vapply(el, function(s) {
      is.list(s) && length(s$mz) > 0L &&
        any(is.finite(suppressWarnings(as.numeric(s$mz))))
    }, logical(1L)))
  }, logical(1L))
  L[keep]
}

#' Subset satellite slots by feature name; never modifies `feature_table`
#' @keywords internal
#' @noRd
.nts_subset_project_satellites_by_name <- function(project, keep_names) {
  keep_names <- unique(as.character(keep_names))
  keep_names <- keep_names[!is.na(keep_names) & nzchar(keep_names)]
  ms1 <- get_ms1_annotations(project)
  if (is.data.frame(ms1) && nrow(ms1) && "name" %in% names(ms1)) {
    project <- set_ms1_annotations(
      project,
      ms1[
        as.character(ms1$name) %in% keep_names,
        ,
        drop = FALSE
      ]
    )
  }
  sp <- methods::slot(project, "ms2data")
  if (length(sp) && !is.null(names(sp))) {
    methods::slot(project, "ms2data") <- sp[names(sp) %in% keep_names]
  }
  cl <- methods::slot(project, "clean_ms2")
  if (length(cl) && !is.null(names(cl))) {
    methods::slot(project, "clean_ms2") <- cl[names(cl) %in% keep_names]
  }
  m2 <- get_ms2_annotations(project)
  if (is.data.frame(m2) && nrow(m2) && "name" %in% names(m2)) {
    project <- set_ms2_annotations(
      project,
      m2[
        as.character(m2$name) %in% keep_names,
        ,
        drop = FALSE
      ]
    )
  }
  db <- methods::slot(project, "db_annotations")
  if (is.data.frame(db) && nrow(db) && "name" %in% names(db)) {
    methods::slot(project, "db_annotations") <- db[
      as.character(db$name) %in% keep_names,
      ,
      drop = FALSE
    ]
  }
  pc <- methods::slot(project, "physchem_proxy")
  if (is.data.frame(pc) && nrow(pc) && "name" %in% names(pc)) {
    methods::slot(project, "physchem_proxy") <- pc[
      as.character(pc$name) %in% keep_names,
      ,
      drop = FALSE
    ]
  }
  project
}

#' @keywords internal
#' @noRd
.nts_merge_layers_for_network <- function(project) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  ms1 <- data.table::as.data.table(get_ms1_annotations(project))
  if (!nrow(ms1) || !"name" %in% names(ms1)) {
    return(ms1[0L])
  }
  db <- .nts_db_best_per_name(methods::slot(project, "db_annotations"))
  if (is.data.frame(db) && nrow(db)) {
    ddt <- data.table::as.data.table(db)
    ms1 <- ddt[ms1, on = "name"]
  }
  m2 <- data.table::as.data.table(get_ms2_annotations(project))
  if (nrow(m2) && "name" %in% names(m2)) {
    m2_only <- setdiff(names(m2), c("name", "mz", "rt"))
    if (length(m2_only)) {
      m2_sub <- m2[, c("name", m2_only), with = FALSE]
      ms1 <- m2_sub[ms1, on = "name"]
    }
  }
  ms1[]
}

#' Merge Top-1 `db_annotations` (and optional `ms2_annotations`) seed fields onto
#' `@ms1_annotations` by `name` using columns **not** stripped by
#' `.nts_ms1_sanitize_and_order` (`library_*`), so downstream network / SMILES
#' tooling can read stable per-feature library hits without polluting forbidden
#' core names.
#' @keywords internal
#' @noRd
.nts_merge_db_best_into_ms1_library_seed_cols <- function(project) {
  if (!methods::is(project, "NTSProject")) {
    return(project)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    return(project)
  }
  ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
  if (!is.data.frame(ms1) || !nrow(ms1) || !"name" %in% names(ms1)) {
    return(project)
  }
  db <- .nts_db_best_per_name(methods::slot(project, "db_annotations"))
  if (!is.data.frame(db) || !nrow(db) || !"name" %in% names(db)) {
    return(project)
  }
  db$name <- trimws(as.character(db$name))
  smi_db <- rep(NA_character_, nrow(db))
  if ("SMILES" %in% names(db)) {
    smi_db <- trimws(as.character(db$SMILES))
  } else if ("smiles" %in% names(db)) {
    smi_db <- trimws(as.character(db$smiles))
  }
  cname_db <- rep(NA_character_, nrow(db))
  if ("Compound.name" %in% names(db)) {
    cname_db <- trimws(as.character(db$Compound.name))
  } else if ("compound_name" %in% names(db)) {
    cname_db <- trimws(as.character(db$compound_name))
  }
  ref_db <- rep(NA_character_, nrow(db))
  if ("ref_id" %in% names(db)) {
    ref_db <- trimws(as.character(db$ref_id))
  } else if ("Lab.ID" %in% names(db)) {
    ref_db <- trimws(as.character(db[["Lab.ID"]]))
  }
  lvl_db <- rep(NA_character_, nrow(db))
  if ("Confidence_Level" %in% names(db)) {
    lvl_db <- trimws(as.character(db$Confidence_Level))
  }
  lib <- data.frame(
    name = db$name,
    library_SMILES = smi_db,
    library_Compound_name = cname_db,
    library_ref_id = ref_db,
    annotation_level = lvl_db,
    stringsAsFactors = FALSE
  )
  lib <- lib[!duplicated(lib$name), , drop = FALSE]
  idx <- match(trimws(as.character(ms1$name)), lib$name)
  ms1$library_SMILES <- lib$library_SMILES[idx]
  ms1$library_Compound_name <- lib$library_Compound_name[idx]
  ms1$library_ref_id <- lib$library_ref_id[idx]

  alvl_fill <- lib$annotation_level[idx]
  w_lvl <- which(!is.na(alvl_fill) & nzchar(alvl_fill))
  if (length(w_lvl)) {
    ms1$annotation_level[w_lvl] <- alvl_fill[w_lvl]
  }

  m2 <- as.data.frame(get_ms2_annotations(project), stringsAsFactors = FALSE)
  if (is.data.frame(m2) && nrow(m2) && "name" %in% names(m2)) {
    m2$name <- trimws(as.character(m2$name))
    smi_m2 <- rep(NA_character_, nrow(m2))
    if ("SMILES" %in% names(m2)) {
      smi_m2 <- trimws(as.character(m2$SMILES))
    } else if ("smiles" %in% names(m2)) {
      smi_m2 <- trimws(as.character(m2$smiles))
    }
    m2s <- data.frame(name = m2$name, m2_lib_smiles = smi_m2, stringsAsFactors = FALSE)
    m2s <- m2s[!duplicated(m2s$name), , drop = FALSE]
    j2 <- match(trimws(as.character(ms1$name)), m2s$name)
    fill <- m2s$m2_lib_smiles[j2]
    w <- which(is.na(ms1$library_SMILES) | !nzchar(as.character(ms1$library_SMILES)))
    if (length(w)) {
      fw <- fill[w]
      ok <- !is.na(fw) & nzchar(as.character(fw))
      if (any(ok)) {
        ms1$library_SMILES[w[ok]] <- as.character(fw[ok])
      }
    }
  }

  project <- set_ms1_annotations(project, .nts_ms1_sanitize_and_order(ms1))
  project
}

#' @keywords internal
#' @noRd
.nts_ft_ms1_view <- function(project) {
  ft <- get_feature_table(project)
  ms1 <- get_ms1_annotations(project)
  if (!nrow(ft)) {
    return(ft)
  }
  df <- ft
  if (nrow(ms1) && "name" %in% names(ms1)) {
    idx <- match(as.character(df$name), as.character(ms1$name))
    for (cl in setdiff(names(ms1), "name")) {
      df[[cl]] <- ms1[[cl]][idx]
    }
  }
  df
}

#' @keywords internal
#' @noRd
.nts_list_ms_files <- function(path,
                               pattern = "\\.[Mm][Zz][Mm][Ll]$|\\.[Mm][Zz][Xx][Mm][Ll]$") {
  files <- list.files(path, pattern = pattern, full.names = TRUE, recursive = TRUE)
  sort(files)
}
