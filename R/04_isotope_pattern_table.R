# Isotope envelope table (M, M+1, ...) export ----------------------------------------------


#' @keywords internal
#' @noRd
.nts_pattern_intensity_vec <- function(ft, sample) {
  ab <- .nts_feature_intensity_colnames(ft)
  ab <- ab[vapply(ab, function(nm) is.numeric(ft[[nm]]), logical(1L))]
  s0 <- sample
  if (is.null(s0) || (is.character(s0) && !nzchar(paste0(s0, collapse = "")))) {
    if (!length(ab)) {
      return(rep(NA_real_, nrow(ft)))
    }
    return(rowMeans(as.matrix(ft[, ab, drop = FALSE]), na.rm = TRUE))
  }
  sig <- as.character(s0)[1L]
  if (!sig %in% names(ft)) {
    stop("Unknown `sample` column in feature_table: ", sig, call. = FALSE)
  }
  suppressWarnings(as.numeric(ft[[sig]]))
}


#' @keywords internal
#' @noRd
.nts_join_isotope_meta_to_ft <- function(ft, ann) {
  meta_cols <- c(
    "is_isotope", "parent_name", "isotope_shift_da", "isotope_kind",
    "isotope_label", "isotope_c13_n"
  )
  n <- nrow(ft)
  side <- data.frame(
    is_isotope = rep(FALSE, n),
    parent_name = rep(NA_character_, n),
    isotope_shift_da = rep(NA_real_, n),
    isotope_kind = rep(NA_character_, n),
    isotope_label = rep(NA_character_, n),
    isotope_c13_n = rep(NA_integer_, n),
    stringsAsFactors = FALSE
  )
  if (!nrow(ann) || !"name" %in% names(ann)) {
    return(side)
  }
  ann2 <- ann[!duplicated(as.character(ann$name)), , drop = FALSE]
  idx <- match(as.character(ft$name), as.character(ann2$name))
  ok <- is.finite(idx)
  for (cl in meta_cols) {
    if (!cl %in% names(ann2)) {
      next
    }
    side[[cl]][ok] <- ann2[[cl]][idx[ok]]
  }
  iso <- side$is_isotope
  iso[is.na(iso)] <- FALSE
  side$is_isotope <- as.logical(iso)
  side
}


#' @keywords internal
#' @noRd
.nts_isotope_peak_sort_key <- function(peak) {
  p <- as.character(peak)
  n <- length(p)
  v <- rep(1e6, n)
  if (!n) {
    return(numeric(0))
  }
  m <- p == "M"
  v[m] <- 0
  g <- grepl("^M\\+[0-9]+$", p) & !m
  if (any(g)) {
    v[g] <- as.numeric(sub("^M\\+", "", p[g]))
  }
  v
}


#' Extract isotope pattern table from a project
#'
#' Returns a **long** table: one row per peak (monoisotopic \code{M}, each
#' \eqn{^{13}}C satellite \code{M+1}, \code{M+2}, ... up to \code{max_c13_order},
#' and optionally halogen satellites from [flag_isotopes]). Columns are
#' suitable for checking m/z, intensity, and annotation consistency in a
#' spreadsheet or \code{View()}.
#'
#' @param project An [NTSProject] after [flag_isotopes] and [extract_features].
#' @param sample `NULL` (default): row-wise mean across numeric sample columns in
#'   `feature_table`; or one column name.
#' @param max_c13_order Maximum \eqn{^{13}}C order (default `3` -> through `M+3`).
#' @param include_halogen If `TRUE`, include rows with `isotope_kind == "halogen"`
#'   (`peak` = label, e.g. `37Cl`).
#' @param require_isotope_child If `TRUE` (default), omit parents with no flagged
#'   isotope child.
#'
#' @return A `data.frame` with columns `parent_name`, `peak`, `feature_name`,
#'   `mz`, `rt`, `intensity`, `isotope_kind`, `isotope_label`,
#'   `isotope_shift_da`, `isotope_c13_n`, `intensity_sample`.
#'
#' @seealso [flag_isotopes]
#' @export
extract_isotope_pattern_table <- function(project,
                                          sample = NULL,
                                          max_c13_order = 3L,
                                          include_halogen = TRUE,
                                          require_isotope_child = TRUE) {
  .nts_check(project, "project", .arg = "project")
  max_c13_order <- as.integer(max_c13_order)[1L]
  if (!is.finite(max_c13_order) || max_c13_order < 1L) {
    stop("max_c13_order must be a positive integer.", call. = FALSE)
  }

  ft <- get_feature_table(project)
  ann <- get_ms1_annotations(project)
  if (!nrow(ft)) {
    stop("feature_table is empty.", call. = FALSE)
  }
  if (!all(c("name", "mz", "rt") %in% names(ft))) {
    stop("feature_table must contain columns name, mz, rt.", call. = FALSE)
  }

  int_vec <- .nts_pattern_intensity_vec(ft, sample)
  sample_lab <- if (is.null(sample) || (is.character(sample) &&
    !nzchar(paste0(sample, collapse = "")))) {
    NA_character_
  } else {
    as.character(sample)[1L]
  }

  side <- .nts_join_isotope_meta_to_ft(ft, ann)
  iso_mask <- side$is_isotope %in% TRUE
  parent_ids <- unique(as.character(side$parent_name[iso_mask]))
  parent_ids <- parent_ids[!is.na(parent_ids) & nzchar(parent_ids)]

  empty_iso <- function() {
    data.frame(
      parent_name = character(0),
      peak = character(0),
      feature_name = character(0),
      mz = numeric(0),
      rt = numeric(0),
      intensity = numeric(0),
      isotope_kind = character(0),
      isotope_label = character(0),
      isotope_shift_da = numeric(0),
      isotope_c13_n = integer(0),
      intensity_sample = character(0),
      stringsAsFactors = FALSE
    )
  }

  if (!length(parent_ids)) {
    warning(
      "No isotope-flagged features (is_isotope=TRUE); run flag_isotopes() first.",
      call. = FALSE
    )
    return(empty_iso())
  }

  long_parts <- vector("list", length(parent_ids))
  for (pi in seq_along(parent_ids)) {
    pid <- parent_ids[pi]
    iM <- match(pid, as.character(ft$name))[1L]
    if (is.na(iM)) {
      next
    }
    ch_idx <- which(iso_mask & as.character(side$parent_name) == pid)
    if (!length(ch_idx) && isTRUE(require_isotope_child)) {
      next
    }

    rows <- list()
    rows[[1L]] <- data.frame(
      parent_name = pid,
      peak = "M",
      feature_name = pid,
      mz = as.numeric(ft$mz[iM]),
      rt = as.numeric(ft$rt[iM]),
      intensity = int_vec[iM],
      isotope_kind = NA_character_,
      isotope_label = NA_character_,
      isotope_shift_da = NA_real_,
      isotope_c13_n = NA_integer_,
      intensity_sample = sample_lab,
      stringsAsFactors = FALSE
    )
    next_row <- 2L
    done <- "M"

    if (length(ch_idx)) {
      for (j in ch_idx) {
        fn <- as.character(ft$name[j])
        kind <- side$isotope_kind[j]
        lbl <- side$isotope_label[j]
        c13n <- side$isotope_c13_n[j]
        sh <- side$isotope_shift_da[j]
        pk <- NA_character_
        if (!is.na(kind) && kind == "13C" &&
            is.finite(c13n) && c13n >= 1L && c13n <= max_c13_order) {
          pk <- paste0("M+", c13n)
        } else if (isTRUE(include_halogen) && !is.na(kind) && kind == "halogen") {
          pk <- if (!is.na(lbl) && nzchar(lbl)) as.character(lbl) else "halogen"
        }
        if (is.na(pk) || pk %in% done) {
          next
        }
        done <- c(done, pk)
        rows[[next_row]] <- data.frame(
          parent_name = pid,
          peak = pk,
          feature_name = fn,
          mz = as.numeric(ft$mz[j]),
          rt = as.numeric(ft$rt[j]),
          intensity = int_vec[j],
          isotope_kind = as.character(kind),
          isotope_label = as.character(lbl),
          isotope_shift_da = suppressWarnings(as.numeric(sh)),
          isotope_c13_n = if (is.finite(c13n)) as.integer(c13n) else NA_integer_,
          intensity_sample = sample_lab,
          stringsAsFactors = FALSE
        )
        next_row <- next_row + 1L
      }
    }
    long_parts[[pi]] <- if (next_row > 2L || !isTRUE(require_isotope_child)) {
      do.call(rbind, rows[seq_len(next_row - 1L)])
    } else {
      NULL
    }
  }

  long_df <- long_parts[!vapply(long_parts, is.null, logical(1L))]
  if (!length(long_df)) {
    warning("No isotope patterns could be assembled.", call. = FALSE)
    return(empty_iso())
  }
  long_df <- do.call(rbind, long_df)
  rownames(long_df) <- NULL
  ord <- order(
    long_df$parent_name,
    .nts_isotope_peak_sort_key(long_df$peak),
    long_df$peak
  )
  long_df <- long_df[ord, , drop = FALSE]
  rownames(long_df) <- NULL
  long_df
}
