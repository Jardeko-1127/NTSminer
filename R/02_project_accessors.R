#' GC-EI full-scan spectra per feature
#'
#' @param x An [NTSProject].
#' @return The `gc_spectra` slot (named `list`).
#' @export
get_gc_spectra <- function(x) {
  .nts_check(x, "project", .arg = "x")
  x@gc_spectra
}

#' @param value Named list of per-feature `data.frame(mz, intensity)`.
#' @rdname get_gc_spectra
#' @export
set_gc_spectra <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  x@gc_spectra <- value
  validObject(x)
  x
}

#' GC-EI deconvolution results
#'
#' @param x An [NTSProject].
#' @return The `gc_deconvolved` slot (named `list` of component lists).
#' @export
get_gc_deconvolved <- function(x) {
  .nts_check(x, "project", .arg = "x")
  x@gc_deconvolved
}

#' @param value Named list of deconvolution components.
#' @rdname get_gc_deconvolved
#' @export
set_gc_deconvolved <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  x@gc_deconvolved <- value
  validObject(x)
  x
}

# --- Getters (NTSProject facade) -------------------------------------------

#' Feature table (immutable abundance core)
#'
#' @param x An [NTSProject].
#' @return The `feature_table` slot (`data.frame` or `data.table`).
#' @export
get_feature_table <- function(x) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "feature_table")
}

#' MS1 qualitative annotations
#'
#' The `ms1_annotations` slot is a **named list of sub-tables**, one per
#' processing step that actually carries data (see [NTSProject] and
#' `.nts_ms1_groups()`). Steps whose columns are absent or entirely `NA` are
#' not stored, so `names()` lists only the populated steps. This getter
#' returns, by default, the merged **flat** `data.frame` view (one row per
#' feature, columns from every step) so existing code keeps working. Pass
#' `step` to read one step's own sub-table instead, which is never polluted by
#' columns from the other steps.
#'
#' @param x An [NTSProject].
#' @param step Optional step name. One of `"base"`, `"eic"`, `"library"`,
#'   `"isotopes"`, `"isf"`, `"kmd"`, `"adduct"`, `"master"`, `"network"`,
#'   `"annotation"`, `"extra"`. `NULL` (default) returns the merged flat view.
#' @return With `step = NULL`, a `data.frame` (0-row when the slot is empty).
#'   With `step` set, that step's sub-table (`data.frame`, possibly 0-row when
#'   the step is valid but carries no data), or `NULL` when the slot itself is
#'   empty.
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' d <- tempfile("nts_ms1"); dir.create(d)
#' p <- create_nts_project(d)
#' p <- set_ms1_annotations(
#'   p,
#'   data.frame(name = "F1", mz = 400, rt = 60, is_isf = TRUE)
#' )
#' names(get_ms1_annotations(p))           # flat merged view
#' get_ms1_annotations(p, "isf")           # just the flag_isf sub-table
#' @export
get_ms1_annotations <- function(x, step = NULL) {
  .nts_check(x, "project", .arg = "x")
  slot <- methods::slot(x, "ms1_annotations")
  if (is.null(step)) {
    ## The slot stores only the steps that carry data (see
    ## [.nts_ms1_partition]); re-apply the canonical schema so the merged flat
    ## view keeps its full legacy column set regardless of which steps ran.
    return(.nts_ms1_sanitize_and_order(.nts_ms1_flatten(slot)))
  }
  step <- as.character(step)[[1L]]
  gm <- .nts_ms1_groups()
  if (!step %in% names(gm)) {
    stop(
      "Unknown `ms1_annotations` step: '", step, "'. Valid steps: ",
      paste(names(gm), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.list(slot) || !length(slot)) {
    return(NULL)
  }
  ## Canonical named-list slot: read the step straight out when present. A
  ## valid step that carries no data was pruned from the slot, so hand back an
  ## empty sub-table rather than `NULL`.
  if (any(vapply(slot, is.data.frame, logical(1L)))) {
    if (step %in% names(slot)) {
      return(slot[[step]])
    }
    return(data.frame())
  }
  ## Legacy slot shapes (flat data.frame, or one mangled into a bare list).
  sub <- .nts_ms1_partition(.nts_ms1_flatten(slot))
  if (is.null(sub[[step]])) data.frame() else sub[[step]]
}

#' MS2 annotations (per-feature view / sub-tables)
#'
#' `@ms2_annotations` is a named `list` mirroring `@ms1_annotations`, holding at
#' most two populated sub-tables: `annotation` (the per-feature MS2 evidence
#' table) and `score_grid` (the full per-cell score grid, also read with
#' [get_ms2_score_grid]). `names()` lists just the populated pieces. This getter
#' returns, by default, the per-feature `annotation` table so existing callers
#' keep working; pass `step` to read one sub-table on its own.
#'
#' @param x An [NTSProject].
#' @param step Optional sub-table name, `"annotation"` or `"score_grid"`. `NULL`
#'   (default) returns the per-feature `annotation` view.
#' @return With `step = NULL` (or `"annotation"`), the per-feature MS2
#'   `data.frame` (0-row when empty). With `step = "score_grid"`, the score-grid
#'   `data.frame` (0-row when empty).
#' @export
get_ms2_annotations <- function(x, step = NULL) {
  .nts_check(x, "project", .arg = "x")
  lst <- .nts_ms2_read_list(x)
  if (is.null(step)) {
    return(if (is.null(lst$annotation)) data.frame() else lst$annotation)
  }
  step <- as.character(step)[[1L]]
  if (!step %in% .nts_ms2_steps()) {
    stop(
      "Unknown `ms2_annotations` step: '", step, "'. Valid steps: ",
      paste(.nts_ms2_steps(), collapse = ", "), ".",
      call. = FALSE
    )
  }
  sub <- lst[[step]]
  if (is.null(sub)) data.frame() else sub
}

#' Cached network results (optional list)
#'
#' @param x An [NTSProject].
#' @return The `network_results` slot.
#' @export
get_network_results <- function(x) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "network_results")
}

#' Library / score hits on an [NTSProject] (`db_annotations` slot)
#'
#' @param x An [NTSProject].
#' @return The `db_annotations` slot (`data.frame`).
#' @export
get_db_annotations <- function(x) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "db_annotations")
}

#' Full MS2 score grid (per feature x compound x CE x spectrum)
#'
#' The `score_grid` element of the `@ms2_annotations` list holds every scored MS2
#' cell produced by [annotate_by_spectral_match], **including cells below
#' `min_score`**. Where `@db_annotations` keeps only the merged result (the best
#' passing cell per feature x compound), this keeps each `(library CE) x
#' (experimental spectrum)` subset separately, so you can see how every
#' experimental MS2 scan matched and re-tune `min_score` offline. Filter with
#' `passes_min_score`, or read `is_winner` for the cell behind each
#' `db_annotations` row.
#'
#' @param x An [NTSProject].
#' @return The `score_grid` sub-table (`data.frame`; 0-row before annotation).
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' d <- tempfile("nts_grid"); dir.create(d)
#' p <- create_nts_project(d)
#' nrow(get_ms2_score_grid(p))   # 0 before annotate_by_spectral_match()
#' @export
get_ms2_score_grid <- function(x) {
  .nts_check(x, "project", .arg = "x")
  # Projects saved before this refactor kept the grid in a standalone
  # `@ms2_score_grid` slot (missing on yet older `.rds`). `.nts_ms2_read_list()`
  # folds either shape into the `score_grid` element and returns an empty grid
  # when nothing has been scored, so legacy objects keep reading.
  sub <- .nts_ms2_read_list(x)$score_grid
  if (is.null(sub)) data.frame() else sub
}

#' @rdname get_ms2_score_grid
#' @param value New `score_grid` `data.frame`.
#' @return `x` with the grid stored under `@ms2_annotations$score_grid`.
#' @export
set_ms2_score_grid <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  set_ms2_annotations(x, value, step = "score_grid")
}

#' Cached CDK-style physchem proxies
#'
#' @param x An [NTSProject].
#' @return The `physchem_proxy` slot (`data.frame`).
#' @export
get_physchem_proxy <- function(x) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "physchem_proxy")
}

#' Experimental MS2 scan lists per feature
#'
#' @param x An [NTSProject].
#' @return The `ms2data` slot (named `list`).
#' @export
get_ms2data <- function(x) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "ms2data")
}

#' Purified MS2 spectra per feature
#'
#' @param x An [NTSProject].
#' @return The `clean_ms2` slot (named `list`).
#' @export
get_clean_ms2 <- function(x) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "clean_ms2")
}

#' Pipeline step log
#'
#' @param x An [NTSProject].
#' @return The `log` slot (`list`).
#' @export
get_project_log <- function(x) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "log")
}

# --- Setters (full slot replacement) --------------------------------------

#' @rdname get_feature_table
#' @param value New `feature_table` object. Per-sample intensity columns are
#'   rounded to whole numbers on assignment (see
#'   `.nts_round_feature_intensities`); `mz`/`rt` are left numeric.
#' @return `x` with updated slot.
#' @export
set_feature_table <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "feature_table") <- .nts_round_feature_intensities(value)
  validObject(x)
  x
}

#' @rdname get_ms1_annotations
#' @param value New `ms1_annotations` content. With `step = NULL`, a flat
#'   `data.frame` (one row per feature) or a full named list of sub-tables;
#'   either is canonicalised, split per step and stored. With `step` set, a
#'   `data.frame` for that one step, keyed by `name`, whose columns replace the
#'   current values of those columns (rows are matched by `name`). Only steps
#'   that carry data are kept in the stored list (see `.nts_ms1_partition`).
#' @return `x` with an updated `ms1_annotations` **list** slot.
#' @export
set_ms1_annotations <- function(x, value, step = NULL) {
  .nts_check(x, "project", .arg = "x")

  if (is.null(step)) {
    methods::slot(x, "ms1_annotations") <-
      .nts_ms1_partition(.nts_ms1_sanitize_and_order(value))
    validObject(x)
    return(x)
  }

  step <- as.character(step)[[1L]]
  gm <- .nts_ms1_groups()
  if (!step %in% names(gm)) {
    stop(
      "Unknown `ms1_annotations` step: '", step, "'. Valid steps: ",
      paste(names(gm), collapse = ", "), ".",
      call. = FALSE
    )
  }
  value <- .nts_ms1_flatten(value)
  if (!is.data.frame(value) || !"name" %in% names(value)) {
    stop(
      "`value` for `step = \"", step,
      "\" must be a `data.frame` with a `name` column.",
      call. = FALSE
    )
  }
  cols <- if (identical(step, "extra")) {
    setdiff(names(value), c("name", unique(unlist(gm, use.names = FALSE))))
  } else {
    intersect(gm[[step]], names(value))
  }
  sub <- value[, unique(c("name", cols)), drop = FALSE]

  cur <- .nts_ms1_flatten(methods::slot(x, "ms1_annotations"))
  if (!nrow(cur) || !"name" %in% names(cur)) {
    flat <- sub
  } else {
    flat <- cur
    idx <- match(as.character(flat$name), as.character(sub$name))
    for (cl in cols) {
      flat[[cl]] <- sub[[cl]][idx]
    }
  }
  methods::slot(x, "ms1_annotations") <-
    .nts_ms1_partition(.nts_ms1_sanitize_and_order(flat))
  validObject(x)
  x
}

#' @rdname get_ms2_annotations
#' @param value New `ms2_annotations` content. With `step = NULL`, either a
#'   flat per-feature `data.frame` (stored as the `annotation` sub-table, the
#'   existing `score_grid` kept) or a full named `list` of sub-tables. With
#'   `step = "annotation"` / `"score_grid"`, a `data.frame` for that sub-table.
#' @return `x` with an updated `ms2_annotations` **list** slot.
#' @export
set_ms2_annotations <- function(x, value, step = NULL) {
  .nts_check(x, "project", .arg = "x")
  lst <- .nts_ms2_read_list(x)
  if (is.null(step)) {
    if (is.data.frame(value)) {
      lst <- .nts_ms2_sanitize(
        list(annotation = value, score_grid = lst$score_grid)
      )
    } else {
      lst <- .nts_ms2_sanitize(value)
    }
  } else {
    step <- as.character(step)[[1L]]
    if (!step %in% .nts_ms2_steps()) {
      stop(
        "Unknown `ms2_annotations` step: '", step, "'. Valid steps: ",
        paste(.nts_ms2_steps(), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (!is.data.frame(value)) {
      stop(
        "`value` for `step = \"", step, "\" must be a `data.frame`.",
        call. = FALSE
      )
    }
    lst[[step]] <- as.data.frame(value, stringsAsFactors = FALSE)
    lst <- .nts_ms2_sanitize(lst)
  }
  methods::slot(x, "ms2_annotations") <- lst
  ## Migrate a pre-refactor object: drop the stray legacy slot attribute now
  ## that its value lives in the list.
  if (!is.null(attr(x, "ms2_score_grid"))) {
    attr(x, "ms2_score_grid") <- NULL
  }
  validObject(x)
  x
}

#' @rdname get_network_results
#' @param value New `network_results` object.
#' @return `x` with updated slot.
#' @export
set_network_results <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "network_results") <- value
  validObject(x)
  x
}

#' @rdname get_db_annotations
#' @param value New `db_annotations` `data.frame`.
#' @export
set_db_annotations <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "db_annotations") <- value
  validObject(x)
  x
}

#' @rdname get_physchem_proxy
#' @param value New `physchem_proxy` `data.frame`.
#' @export
set_physchem_proxy <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "physchem_proxy") <- value
  validObject(x)
  x
}

#' @rdname get_ms2data
#' @param value New `ms2data` named `list`.
#' @export
set_ms2data <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "ms2data") <- value
  validObject(x)
  x
}

#' @rdname get_clean_ms2
#' @param value New `clean_ms2` named `list`.
#' @export
set_clean_ms2 <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "clean_ms2") <- value
  validObject(x)
  x
}

#' @rdname get_project_log
#' @param value New `log` `list`.
#' @export
set_project_log <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  methods::slot(x, "log") <- value
  validObject(x)
  x
}

#' CSI:FingerID fingerprint table
#'
#' @param x An [NTSProject].
#' @param as_matrix Logical. If `TRUE`, returns a numeric matrix with
#'   rows named by `feature | mz=M/Z | RT=seconds` and columns named
#'   `FP_1` ... `FP_N`. Only features with non-empty fingerprint vectors
#'   are included. Default `FALSE` returns the raw [data.table::data.table].
#' @return The `sirius_fingerprints` slot: a [data.table::data.table]
#'   by default, or a numeric matrix when `as_matrix = TRUE`.
#' @export
get_sirius_fingerprints <- function(x, as_matrix = FALSE) {
  .nts_check(x, "project", .arg = "x")
  fp <- methods::slot(x, "sirius_fingerprints")

  if (!as_matrix) return(fp)

  # Build fingerprint matrix with informative row/column names
  ok <- lengths(fp$fingerprint) > 0L
  if (!any(ok)) {
    warning("No non-empty fingerprint vectors in @sirius_fingerprints.",
            call. = FALSE)
    return(matrix(numeric(0), nrow = 0L, ncol = 0L))
  }

  ms1 <- get_ms1_annotations(x)
  if (is.data.frame(ms1) && nrow(ms1) && all(c("mz", "rt") %in% names(ms1))) {
    mz_val <- ms1$mz[match(fp$name[ok], ms1$name)]
    rt_val <- ms1$rt[match(fp$name[ok], ms1$name)]
    row_labels <- paste0(
      fp$name[ok], " | mz=", round(mz_val, 4L), " | RT=", round(rt_val, 1L)
    )
  } else {
    row_labels <- fp$name[ok]
  }

  fmat <- do.call(rbind, fp$fingerprint[ok])
  rownames(fmat) <- row_labels
  colnames(fmat) <- paste0("FP_", seq_len(ncol(fmat)))

  fmat
}

#' @rdname get_sirius_fingerprints
#' @param value New `sirius_fingerprints` `data.table` with required columns:
#'   `name`, `adduct`, `molecular_formula`, `rank`, `sirius_score`,
#'   `zodiac_score`, `confidence`, `explained_peaks`, `fingerprint` (list),
#'   `fp_version`.
#' @export
set_sirius_fingerprints <- function(x, value) {
  .nts_check(x, "project", .arg = "x")
  .validate_sirius_fingerprints(value)
  # The slot is declared `data.table` (see 00_AllClasses.R) while the validator
  # above deliberately accepts a plain `data.frame`. Without this coercion a
  # data.frame assignment dies inside `checkSlotAssignment()` with a confusing
  # "is(value, \"data.table\") is not TRUE" error -- which is exactly what the
  # 0-row path of compute_sirius_fingerprints() used to hit whenever SIRIUS
  # produced no results. Coerce here so the accessor honours the validator.
  if (!is.null(value) && is.data.frame(value) &&
      !data.table::is.data.table(value)) {
    value <- data.table::as.data.table(value)
  }
  methods::slot(x, "sirius_fingerprints") <- value
  validObject(x)
  x
}

# --- Sirius fingerprints validator (internal) ---

.REQUIRED_SIRIUS_COLS <- c(
  "name", "adduct", "molecular_formula", "rank",
  "sirius_score", "zodiac_score", "confidence",
  "explained_peaks", "fingerprint", "fp_version",
  "canopus_npc", "canopus_superclass", "canopus_class", "canopus_subclass"
)

.validate_sirius_fingerprints <- function(value) {
  # Allow empty (0-row) table
  if (is.null(value)) return(invisible(TRUE))
  if (!data.table::is.data.table(value) && !is.data.frame(value)) {
    stop("`value` must be a data.table or data.frame.", call. = FALSE)
  }
  if (nrow(value) == 0L) return(invisible(TRUE))

  missing <- setdiff(.REQUIRED_SIRIUS_COLS, names(value))
  if (length(missing)) {
    stop(
      "Missing required columns in sirius_fingerprints: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Validate fingerprint is a list-column
  if (!is.list(value$fingerprint)) {
    stop(
      "`fingerprint` column must be a list-column (each element a numeric vector).",
      call. = FALSE
    )
  }

  invisible(TRUE)
}
