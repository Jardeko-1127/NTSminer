#' Summarize shared fragment / neutral-loss features by metadata subclass
#'
#' Pipeline: (1) read `Lab_ID` and subclass from `@metadata`; (2) attach subclass
#' to `@features_long` by `Lab_ID` (after `trimws`); (3) cluster m/z per
#' `(subclass, Feature_Type)` using **`mz_tol_ppm`** or **`mz_digits`** (or
#' **`mz_resolutions`** for multi-resolution sweep); (4) compute **`prevalence`**
#' = `n_present / n_lab_class` where **`n_lab_class`** follows
#' **`prevalence_denom`** (**`metadata`** vs compounds with features).
#'
#' @param db An [NTSDatabase].
#' @param class_col Metadata column for subclass labels (default `"Second_Class"`).
#' @param subclass_filter `NULL` or character vector of subclass **values** to keep.
#' @param lab_id_col Metadata column matching `features_long$Lab_ID` (default `"Lab_ID"`).
#' @param polarity `NULL` or `"POS"` / `"NEG"` to filter `features_long`.
#' @param feature_types Character vector of `Feature_Type` (default NL, DF, MD).
#' @param min_prevalence Keep rows with `prevalence >=` this (default `0.5`).
#' @param min_n_lab_class Minimum distinct `Lab_ID` per subclass in **metadata** (default `2`).
#' @param min_frequency Optional minimum `Frequency` in `features_long`.
#' @param prevalence_denom `metadata` (default) or `with_features` for prevalence denominator.
#' @param mz_tol_ppm `NULL` uses `mz_digits` only; default `10` for consecutive ppm clustering.
#'   Ignored when `mz_resolutions` is provided.
#' @param mz_tol_da Optional absolute gap (Da) in ppm clustering. Ignored when
#'   `mz_resolutions` is provided.
#' @param mz_digits Used when `mz_tol_ppm` is `NULL`; partial named vector merges
#'   with defaults. Ignored when `mz_resolutions` is provided.
#' @param mz_resolutions `NULL` or a **named list** of resolution specifications
#'   for multi-resolution sweep. Each element must be a list with **one** of:
#'   - `list(ppm = <numeric>)` - ppm consecutive clustering (optional `da = ...`)
#'   - `list(digits = <named/int vector>)` - round-based clustering
#'   Example: `list("digits3" = list(digits = c(DF=3, NL=3, MD=3)),
#'                   "digits2" = list(digits = c(DF=2, NL=2, MD=2)),
#'                   "ppm50"   = list(ppm = 50))`.
#'   When provided, `mz_tol_ppm`, `mz_tol_da`, `mz_digits` are **ignored** and
#'   output tables gain a **`resolution_label`** column.
#' @param quiet If `TRUE`, suppress messages.
#'
#' @return `invisible(list(by_feature, by_feature_candidates, class_counts,
#'   subclass_coverage, params))`. **`subclass_coverage`** has **`n_lab_metadata`**,
#'   **`n_lab_with_features`**, **`max_prevalence`** (best pre-cutoff prevalence per subclass).
#'   When `mz_resolutions` is used, the first two data frames include a
#'   **`resolution_label`** column.
#'
#' @seealso [plot_subclass_fragment_features]
#' @export
summarize_subclass_fragment_features <- function(
    db,
    class_col = "Second_Class",
    subclass_filter = NULL,
    lab_id_col = "Lab_ID",
    polarity = NULL,
    feature_types = c("NL", "DF", "MD"),
    min_prevalence = 0.5,
    min_n_lab_class = 2L,
    min_frequency = NULL,
    prevalence_denom = c("metadata", "with_features"),
    mz_tol_ppm = 10,
    mz_tol_da = NULL,
    mz_digits = NULL,
    mz_resolutions = NULL,
    quiet = FALSE) {
  if (!methods::is(db, "NTSDatabase")) {
    stop("`db` must be an NTSDatabase", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  class_col <- as.character(class_col)[1L]
  lab_id_col <- as.character(lab_id_col)[1L]
  mp <- suppressWarnings(as.numeric(min_prevalence)[1L])
  if (!is.finite(mp) || mp <= 0 || mp > 1) {
    stop("min_prevalence must be in (0, 1]", call. = FALSE)
  }
  min_n_lab_class <- suppressWarnings(as.integer(min_n_lab_class)[1L])
  if (!is.finite(min_n_lab_class) || min_n_lab_class < 1L) {
    stop("min_n_lab_class must be a positive integer", call. = FALSE)
  }
  quiet <- isTRUE(quiet)
  prevalence_denom <- match.arg(prevalence_denom, c("metadata", "with_features"))
  use_multi_res <- !is.null(mz_resolutions)

  md <- data.table::as.data.table(methods::slot(db, "metadata"))
  fl <- data.table::as.data.table(methods::slot(db, "features_long"))

  if (!lab_id_col %in% names(md)) {
    stop("metadata has no column `", lab_id_col, "` (see lab_id_col).", call. = FALSE)
  }
  if (!class_col %in% names(md)) {
    stop(
      "metadata has no column named `", class_col, "`.\n",
      "  `class_col` is the column **name** holding subclasses per row (e.g. Second_Class).\n",
      "  Use `subclass_filter = \"PFCAs\"` to restrict **values**.\n",
      "  Available: ", paste(head(names(md), 30L), collapse = ", "),
      if (ncol(md) > 30L) ", ..." else "",
      call. = FALSE
    )
  }

  empty_tbl <- function() {
    data.table::data.table(
      subclass = character(),
      Feature_Type = character(),
      mz_key = numeric(),
      mz_value = numeric(),
      n_lab_class = integer(),
      n_present = integer(),
      n_fragment_rows = integer(),
      prevalence = numeric(),
      median_frequency = numeric(),
      mean_frequency = numeric()
    )
  }

  if (!nrow(fl)) {
    if (!quiet) message("summarize_subclass_fragment_features: @features_long is empty.")
    cc <- unique(md[
      ,
      .(
        subclass = trimws(as.character(.SD[[1L]])),
        Lab_ID = trimws(as.character(.SD[[2L]]))
      ),
      .SDcols = c(class_col, lab_id_col)
    ])
    cc <- cc[!is.na(subclass) & nzchar(subclass)]
    nlab <- if (nrow(cc)) {
      cc[, .(n_lab_class = data.table::uniqueN(Lab_ID)), by = subclass]
    } else {
      data.table::data.table(subclass = character(), n_lab_class = integer())
    }
    ns0 <- nlab[n_lab_class >= min_n_lab_class]
    scov <- if (nrow(ns0)) {
      o <- data.table::copy(ns0)
      data.table::setnames(o, "n_lab_class", "n_lab_metadata")
      o[, `:=`(n_lab_with_features = 0L, max_prevalence = NA_real_)]
      o
    } else {
      data.table::data.table(
        subclass = character(),
        n_lab_metadata = integer(),
        n_lab_with_features = integer(),
        max_prevalence = numeric()
      )
    }
    return(invisible(list(
      by_feature = empty_tbl(),
      by_feature_candidates = NULL,
      class_counts = nlab,
      subclass_coverage = scov,
      params = list(
        class_col = class_col,
        subclass_filter = subclass_filter,
        lab_id_col = lab_id_col,
        polarity = polarity,
        feature_types = feature_types,
        min_prevalence = mp,
        min_n_lab_class = min_n_lab_class,
        min_frequency = min_frequency,
        prevalence_denom = prevalence_denom,
        mz_tol_ppm = mz_tol_ppm,
        mz_tol_da = mz_tol_da,
        use_ppm_clustering = !is.null(mz_tol_ppm),
        mz_digits = if (is.null(mz_tol_ppm)) {
          if (is.null(mz_digits)) .nts_default_features_long_mz_digits() else mz_digits
        } else {
          NULL
        }
      )
    )))
  }

  nmap <- c(
    lab_id       = "Lab_ID",
    Lab_ID       = "Lab_ID",
    polarity     = "Polarity",
    Polarity     = "Polarity",
    feature_type = "Feature_Type",
    Feature_Type = "Feature_Type",
    mz           = "mz_value",
    mz_value     = "mz_value",
    frequency    = "Frequency",
    Frequency    = "Frequency"
  )
  oldn <- names(fl)
  rename_mask <- tolower(oldn) %in% tolower(names(nmap)) & oldn != nmap[oldn]
  if (any(rename_mask)) {
    data.table::setnames(fl, oldn[rename_mask], nmap[oldn[rename_mask]])
  }
  fneed <- c("Lab_ID", "Polarity", "Feature_Type", "mz_value", "Frequency")
  miss <- setdiff(fneed, names(fl))
  if (length(miss)) {
    stop(
      "@features_long missing columns: ", paste(miss, collapse = ", "),
      ". Available: ", paste(utils::head(names(fl), 20L), collapse = ", "),
      call. = FALSE
    )
  }

  md2 <- md[
    ,
    .(
      Lab_ID = trimws(as.character(.SD[[1L]])),
      subclass = trimws(as.character(.SD[[2L]]))
    ),
    .SDcols = c(lab_id_col, class_col)
  ]
  md2 <- md2[!is.na(subclass) & nzchar(subclass)]
  if (!is.null(subclass_filter) && length(subclass_filter)) {
    sf <- unique(trimws(as.character(subclass_filter)))
    sf <- sf[!is.na(sf) & nzchar(sf)]
    if (length(sf)) md2 <- md2[subclass %in% sf]
  }
  if (!nrow(md2)) {
    stop(
      "No rows left after subclass filters in metadata$", class_col, ".",
      call. = FALSE
    )
  }

  nlab_sub <- md2[, .(n_lab_class = data.table::uniqueN(Lab_ID)), by = subclass]
  nlab_sub <- nlab_sub[n_lab_class >= min_n_lab_class]
  if (!nrow(nlab_sub)) {
    stop(
      "Every subclass has fewer than min_n_lab_class (", min_n_lab_class,
      ") lab entries in metadata after filtering.",
      call. = FALSE
    )
  }

  fl2 <- data.table::copy(fl)
  fl2[, Lab_ID := trimws(as.character(Lab_ID))]
  fl2[, Polarity := toupper(as.character(Polarity))]
  fl2[, Feature_Type := toupper(as.character(Feature_Type))]
  ft_ok <- toupper(as.character(feature_types))
  fl2 <- fl2[Feature_Type %in% ft_ok]
  if (!nrow(fl2) && !quiet) {
    message("summarize_subclass_fragment_features: no rows after feature_types filter.")
  }

  if (!is.null(polarity)) {
    pol <- toupper(as.character(polarity)[1L])
    fl2 <- fl2[Polarity == pol]
  }

  fl2[, mz_value := as.numeric(as.character(mz_value))]
  fl2[, Frequency := as.numeric(as.character(Frequency))]
  fl2 <- fl2[is.finite(mz_value)]
  if (!is.null(min_frequency)) {
    th <- suppressWarnings(as.numeric(min_frequency)[1L])
    if (is.finite(th)) {
      fl2 <- fl2[is.finite(Frequency) & Frequency >= th]
    }
  }

  fl2 <- md2[fl2, on = "Lab_ID", nomatch = 0L]
  if (!nrow(fl2)) {
    if (!quiet) {
      message("summarize_subclass_fragment_features: no features_long rows join metadata.")
    }
    scov <- data.table::copy(nlab_sub)
    data.table::setnames(scov, "n_lab_class", "n_lab_metadata")
    scov[, `:=`(n_lab_with_features = 0L, max_prevalence = NA_real_)]
    return(invisible(list(
      by_feature = empty_tbl(),
      by_feature_candidates = NULL,
      class_counts = nlab_sub[],
      subclass_coverage = scov,
      params = list(
        class_col = class_col,
        subclass_filter = subclass_filter,
        lab_id_col = lab_id_col,
        polarity = polarity,
        feature_types = ft_ok,
        min_prevalence = mp,
        min_n_lab_class = min_n_lab_class,
        min_frequency = min_frequency,
        prevalence_denom = prevalence_denom,
        mz_resolutions = mz_resolutions,
        use_ppm_clustering = if (use_multi_res) NA else !is.null(mz_tol_ppm),
        mz_tol_ppm = if (use_multi_res) NULL else mz_tol_ppm,
        mz_tol_da = if (use_multi_res) NULL else mz_tol_da,
        mz_digits = if (use_multi_res) NULL else .nts_resolve_features_long_mz_digits(mz_digits)
      )
    )))
  }

  # ---- resolution loop helper: cluster mz_key, aggregate, filter ----
  .run_one_resolution <- function(fl_res, res_label,
                                   ppm = NULL, da = NULL,
                                   digits_spec = NULL,
                                   nlab_sub, mp, quiet) {
    if (!is.null(ppm)) {
      fl_res[
        ,
        cluster_id := {
          if (.N <= 1L) 1L else {
            o <- order(mz_value)
            cid <- .nts_cluster_ids_consecutive_ppm(mz_value[o], ppm, da)
            cid[order(o)]
          }
        },
        by = .(subclass, Feature_Type)
      ]
      fl_res[, mz_key := stats::median(mz_value), by = .(subclass, Feature_Type, cluster_id)]
      fl_res[, cluster_id := NULL]
    } else {
      stopifnot(!is.null(digits_spec))
      ndig <- names(digits_spec)
      dig_vec <- setNames(digits_spec, toupper(names(digits_spec)))
      fl_res[, mz_key := round(mz_value, {
        ft <- toupper(Feature_Type[1L])
        d <- dig_vec[[ft]]
        if (is.null(d) || !is.finite(d)) 5L else as.integer(d)
      }), by = Feature_Type]
    }

    frag_n <- fl_res[, .(n_fragment_rows = .N), by = .(subclass, Feature_Type, mz_key)]
    feat_lab <- fl_res[, .(n_lab_with_features = data.table::uniqueN(Lab_ID)), by = subclass]

    scov <- data.table::copy(nlab_sub)
    data.table::setnames(scov, "n_lab_class", "n_lab_metadata")
    scov[feat_lab, n_lab_with_features := i.n_lab_with_features, on = "subclass"]
    scov[is.na(n_lab_with_features), n_lab_with_features := 0L]
    scov[, max_prevalence := NA_real_]

    lab_feat <- fl_res[
      ,
      .(
        freq_lab = {
          m <- suppressWarnings(max(Frequency, na.rm = TRUE))
          if (!is.finite(m)) NA_real_ else m
        },
        mz_rep = mean(mz_value, na.rm = TRUE)
      ),
      by = .(Lab_ID, subclass, Feature_Type, mz_key)
    ]
    lab_feat <- lab_feat[is.finite(mz_rep)]
    lab_feat <- lab_feat[subclass %in% nlab_sub$subclass]

    mz_mean <- lab_feat[
      ,
      .(mz_value = stats::median(mz_rep, na.rm = TRUE)),
      by = .(subclass, Feature_Type, mz_key)
    ]

    agg <- lab_feat[
      ,
      .(
        n_present = data.table::uniqueN(Lab_ID),
        median_frequency = suppressWarnings(
          stats::median(freq_lab[is.finite(freq_lab)], na.rm = TRUE)
        ),
        mean_frequency = mean(freq_lab[is.finite(freq_lab)], na.rm = TRUE)
      ),
      by = .(subclass, Feature_Type, mz_key)
    ]
    agg <- nlab_sub[agg, on = "subclass"]
    if (identical(prevalence_denom, "with_features")) {
      n_feat_labs <- fl_res[, .(n_lab_class = data.table::uniqueN(Lab_ID)), by = subclass]
      agg[n_feat_labs, n_lab_class := i.n_lab_class, on = "subclass"]
    }
    agg[is.na(n_lab_class) | n_lab_class <= 0L, n_lab_class := 1L]
    agg[, prevalence := n_present / n_lab_class]
    agg <- mz_mean[agg, on = c("subclass", "Feature_Type", "mz_key")]
    agg[frag_n, n_fragment_rows := i.n_fragment_rows,
      on = c("subclass", "Feature_Type", "mz_key")]
    data.table::setorder(agg, subclass, Feature_Type, -prevalence, mz_key)

    agg_all <- agg[]
    if (nrow(agg_all)) {
      smx <- agg_all[, .(max_prevalence = max(prevalence, na.rm = TRUE)), by = subclass]
      scov[smx, max_prevalence := i.max_prevalence, on = "subclass"]
    }

    agg_out <- agg[prevalence >= mp]
    cand <- NULL
    if (!nrow(agg_out)) {
      cand <- head(agg_all, 50L)
    }

    list(by_feature = agg_out[], candidates = cand, subclass_coverage = scov[])
  }

  # ---- multi-resolution sweep ----
  if (use_multi_res) {
    res_specs <- .nts_parse_mz_resolutions(mz_resolutions)
    all_bf <- list()
    all_cand <- list()
    all_scov <- list()
    res_labels <- character()

    for (ri in seq_along(res_specs)) {
      rlabel <- names(res_specs)[ri]
      spec <- res_specs[[ri]]
      fl_res <- data.table::copy(fl2)

      if (!is.null(spec$ppm)) {
        r <- .run_one_resolution(fl_res, rlabel,
          ppm = spec$ppm, da = if (!is.null(spec$da)) spec$da else NULL,
          digits_spec = NULL, nlab_sub = nlab_sub,
          mp = mp, quiet = quiet)
      } else if (!is.null(spec$digits)) {
        r <- .run_one_resolution(fl_res, rlabel,
          ppm = NULL, da = NULL,
          digits_spec = spec$digits,
          nlab_sub = nlab_sub,
          mp = mp, quiet = quiet)
      } else {
        next
      }

      r$by_feature[, resolution_label := rlabel]
      r$subclass_coverage[, resolution_label := rlabel]
      all_bf[[ri]] <- r$by_feature
      all_cand[[ri]] <- r$candidates
      all_scov[[ri]] <- r$subclass_coverage
      res_labels <- c(res_labels, rlabel)
    }

    bf <- data.table::rbindlist(all_bf, use.names = TRUE, fill = TRUE)
    scov <- data.table::rbindlist(all_scov, use.names = TRUE, fill = TRUE)
    cand <- NULL
    cand_list <- all_cand[!vapply(all_cand, is.null, logical(1L))]
    if (length(cand_list)) {
      cand <- data.table::rbindlist(lapply(cand_list, function(cd) {
        cd[, resolution_label := "?"]
        cd
      }), use.names = TRUE, fill = TRUE)
    }

    if (!nrow(bf) && !quiet) {
      message(
        "summarize_subclass_fragment_features: no feature reaches min_prevalence = ",
        mp, " at any resolution. See $by_feature_candidates."
      )
    }

    return(invisible(list(
      by_feature = bf,
      by_feature_candidates = cand,
      class_counts = nlab_sub[],
      subclass_coverage = scov,
      params = list(
        class_col = class_col,
        subclass_filter = subclass_filter,
        lab_id_col = lab_id_col,
        polarity = polarity,
        feature_types = ft_ok,
        min_prevalence = mp,
        min_n_lab_class = min_n_lab_class,
        min_frequency = min_frequency,
        prevalence_denom = prevalence_denom,
        mz_resolutions = res_labels
      )
    )))
  }

  # ---- single-resolution (backward compatible) ----
  use_ppm <- !is.null(mz_tol_ppm)
  ppm_use <- if (use_ppm) {
    p <- suppressWarnings(as.numeric(mz_tol_ppm)[1L])
    if (!is.finite(p) || p < 0) 10 else p
  } else {
    NULL
  }
  da_use <- if (!is.null(mz_tol_da)) {
    d <- suppressWarnings(as.numeric(mz_tol_da)[1L])
    if (!is.finite(d) || d < 0) NULL else d
  } else {
    NULL
  }
  dig_map <- if (!use_ppm) {
    .nts_resolve_features_long_mz_digits(mz_digits)
  } else {
    NULL
  }

  if (use_ppm) {
    fl2[
      ,
      cluster_id := {
        if (.N <= 1L) return(1L)
        o <- order(mz_value)
        cid <- .nts_cluster_ids_consecutive_ppm(mz_value[o], ppm_use, da_use)
        cid[order(o)]
      },
      by = .(subclass, Feature_Type)
    ]
    fl2[, mz_key := stats::median(mz_value), by = .(subclass, Feature_Type, cluster_id)]
    fl2[, cluster_id := NULL]
  } else {
    ndig <- names(dig_map)
    dig_vec <- setNames(dig_map, toupper(names(dig_map)))
    fl2[, mz_key := round(mz_value, {
      ft <- toupper(Feature_Type[1L])
      d <- dig_vec[[ft]]
      if (is.null(d) || !is.finite(d)) 5L else as.integer(d)
    }), by = Feature_Type]
  }

  frag_n <- fl2[, .(n_fragment_rows = .N), by = .(subclass, Feature_Type, mz_key)]

  feat_lab <- fl2[, .(n_lab_with_features = data.table::uniqueN(Lab_ID)), by = subclass]
  subclass_coverage <- data.table::copy(nlab_sub)
  data.table::setnames(subclass_coverage, "n_lab_class", "n_lab_metadata")
  subclass_coverage[feat_lab, n_lab_with_features := i.n_lab_with_features,
    on = "subclass"]
  subclass_coverage[is.na(n_lab_with_features), n_lab_with_features := 0L]
  subclass_coverage[, max_prevalence := NA_real_]

  if (!quiet) {
    no_feat <- subclass_coverage[n_lab_metadata > 0L & n_lab_with_features == 0L]
    if (nrow(no_feat)) {
      message(
        "summarize_subclass_fragment_features: ",
        nrow(no_feat), " subclass(es) have metadata but no features_long after join/filters. ",
        "Example: ", paste(utils::head(no_feat$subclass, 3L), collapse = ", "), "."
      )
    }
  }

  lab_feat <- fl2[
    ,
    .(
      freq_lab = {
        m <- suppressWarnings(max(Frequency, na.rm = TRUE))
        if (!is.finite(m)) NA_real_ else m
      },
      mz_rep = mean(mz_value, na.rm = TRUE)
    ),
    by = .(Lab_ID, subclass, Feature_Type, mz_key)
  ]
  lab_feat <- lab_feat[is.finite(mz_rep)]
  lab_feat <- lab_feat[subclass %in% nlab_sub$subclass]

  mz_mean <- lab_feat[
    ,
    .(mz_value = stats::median(mz_rep, na.rm = TRUE)),
    by = .(subclass, Feature_Type, mz_key)
  ]

  agg <- lab_feat[
    ,
    .(
      n_present = data.table::uniqueN(Lab_ID),
      median_frequency = suppressWarnings(
        stats::median(freq_lab[is.finite(freq_lab)], na.rm = TRUE)
      ),
      mean_frequency = mean(freq_lab[is.finite(freq_lab)], na.rm = TRUE)
    ),
    by = .(subclass, Feature_Type, mz_key)
  ]
  agg <- nlab_sub[agg, on = "subclass"]
  if (identical(prevalence_denom, "with_features")) {
    n_feat_labs <- fl2[, .(n_lab_class = data.table::uniqueN(Lab_ID)), by = subclass]
    agg[n_feat_labs, n_lab_class := i.n_lab_class, on = "subclass"]
  }
  agg[is.na(n_lab_class) | n_lab_class <= 0L, n_lab_class := 1L]
  agg[, prevalence := n_present / n_lab_class]
  agg <- mz_mean[agg, on = c("subclass", "Feature_Type", "mz_key")]
  agg[frag_n, n_fragment_rows := i.n_fragment_rows,
    on = c("subclass", "Feature_Type", "mz_key")]
  data.table::setorder(agg, subclass, Feature_Type, -prevalence, mz_key)

  agg_all <- agg[]
  if (nrow(agg_all)) {
    smx <- agg_all[, .(max_prevalence = max(prevalence, na.rm = TRUE)), by = subclass]
    subclass_coverage[smx, max_prevalence := i.max_prevalence, on = "subclass"]
  }

  agg_out <- agg[prevalence >= mp]

  if (!quiet && nrow(agg_all)) {
    prev_drop <- setdiff(unique(agg_all$subclass), unique(agg_out$subclass))
    prev_drop <- intersect(
      prev_drop,
      subclass_coverage[n_lab_with_features > 0L]$subclass
    )
    if (length(prev_drop)) {
      mx_sub <- agg_all[
        subclass %in% prev_drop,
        .(max_prev = max(prevalence, na.rm = TRUE)),
        by = subclass
      ]
      message(
        "summarize_subclass_fragment_features: min_prevalence removed all rows for: ",
        paste(utils::head(mx_sub$subclass, 6L), collapse = ", "),
        ". Try prevalence_denom = \"with_features\" or lower min_prevalence."
      )
    }
  }

  cand <- NULL
  if (!nrow(agg_out)) {
    mx <- if (nrow(agg_all)) max(agg_all$prevalence, na.rm = TRUE) else NA_real_
    if (!quiet) {
      message(
        "summarize_subclass_fragment_features: no feature reaches min_prevalence = ",
        mp, ". Max prevalence = ",
        if (is.finite(mx)) format(mx, digits = 3) else "NA",
        ". See $by_feature_candidates."
      )
    }
    cand <- head(agg_all, 50L)
  }

  invisible(list(
    by_feature = agg_out[],
    by_feature_candidates = cand,
    class_counts = nlab_sub[],
    subclass_coverage = subclass_coverage[],
    params = list(
      class_col = class_col,
      subclass_filter = subclass_filter,
      lab_id_col = lab_id_col,
      polarity = polarity,
      feature_types = ft_ok,
      min_prevalence = mp,
      min_n_lab_class = min_n_lab_class,
      min_frequency = min_frequency,
      prevalence_denom = prevalence_denom,
      use_ppm_clustering = use_ppm,
      mz_tol_ppm = if (use_ppm) ppm_use else NULL,
      mz_tol_da = da_use,
      mz_digits = dig_map
    )
  ))
}


#' @keywords internal
#' @noRd
.nts_parse_mz_resolutions <- function(x) {
  stopifnot(is.list(x), !is.null(names(x)), all(nzchar(names(x))))
  out <- list()
  for (i in seq_along(x)) {
    label <- names(x)[i]
    spec <- x[[i]]
    if (!is.list(spec)) spec <- as.list(spec)
    nms <- tolower(names(spec))

    if ("ppm" %in% nms) {
      p <- suppressWarnings(as.numeric(spec[[match("ppm", nms)]][1L]))
      if (!is.finite(p) || p < 0) stop("ppm must be >= 0 in mz_resolutions[[\"", label, "\"]]", call. = FALSE)
      da_val <- NULL
      if ("da" %in% nms) {
        d <- suppressWarnings(as.numeric(spec[[match("da", nms)]][1L]))
        if (is.finite(d) && d >= 0) da_val <- d
      }
      out[[label]] <- list(ppm = p, da = da_val)
    } else if ("digits" %in% nms) {
      dig <- spec[[match("digits", nms)]]
      if (is.numeric(dig) && length(dig) >= 1L) {
        if (is.null(names(dig))) {
          d <- as.integer(dig[1L])
          out[[label]] <- list(digits = c(DF = d, NL = d, MD = d))
        } else {
          resolved <- .nts_resolve_features_long_mz_digits(dig)
          out[[label]] <- list(digits = resolved)
        }
      } else {
        stop("digits in mz_resolutions[[\"", label, "\"]] must be numeric", call. = FALSE)
      }
    } else {
      stop("Each element of mz_resolutions must contain `ppm` or `digits`. Problem: \"", label, "\"", call. = FALSE)
    }
  }
  out
}


#' @keywords internal
#' @noRd
.nts_cluster_ids_consecutive_ppm <- function(mz_sorted, ppm, mz_tol_da = NULL) {
  n <- length(mz_sorted)
  if (n == 0L) return(integer())
  if (n == 1L) return(1L)
  ppm <- suppressWarnings(as.numeric(ppm)[1L])
  if (!is.finite(ppm) || ppm < 0) ppm <- 10
  da <- NULL
  if (!is.null(mz_tol_da)) {
    da <- suppressWarnings(as.numeric(mz_tol_da)[1L])
    if (!is.finite(da) || da < 0) da <- NULL
  }
  cid <- integer(n)
  cid[1L] <- 1L
  kid <- 1L
  for (i in 2L:n) {
    dm <- mz_sorted[i] - mz_sorted[i - 1L]
    ref <- mz_sorted[i - 1L]
    within <- FALSE
    if (is.finite(ref) && ref > 0) {
      within <- (dm / ref) * 1e6 <= ppm
    }
    if (!within && !is.null(da) && is.finite(dm)) {
      within <- dm <= da
    }
    if (within) {
      cid[i] <- cid[i - 1L]
    } else {
      kid <- kid + 1L
      cid[i] <- kid
    }
  }
  cid
}


#' Heatmap of subclass vs fragment prevalence (or median_frequency)
#'
#' @param summary Return of [summarize_subclass_fragment_features], or a `data.table`
#'   with columns like `$by_feature`.
#' @param value `prevalence` or `median_frequency`.
#' @param max_features Max y-axis categories.
#' @param feature_select `union_per_subclass` or `global_max`.
#' @param subclasses Optional character vector; ranking uses these subclasses only.
#' @param subclasses_ignore_case Logical.
#' @param label_digits Rounding for m/z labels / collapse (default 3).
#' @param collapse_mz_label Merge same rounded m/z per subclass for display.
#' @param sort_y `mz_asc`, `mz_desc`, or `selection`.
#' @param title Optional plot title.
#' @param resolution_label Optional character - select one resolution from a
#'   multi-resolution summary (specific value of the `resolution_label` column).
#'   When `NULL` (default) and the table has a `resolution_label` column, the
#'   first resolution is used.
#' @param print_plot Print ggplot in interactive use.
#'
#' @return ggplot2 object (invisibly).
#' @importFrom ggplot2 aes ggplot geom_tile scale_fill_gradient labs theme theme_bw
#'   element_text coord_cartesian
#' @export
plot_subclass_fragment_features <- function(
    summary,
    value = c("prevalence", "median_frequency"),
    max_features = 50L,
    feature_select = c("union_per_subclass", "global_max"),
    subclasses = NULL,
    subclasses_ignore_case = FALSE,
    label_digits = 3L,
    collapse_mz_label = TRUE,
    sort_y = c("mz_asc", "mz_desc", "selection"),
    title = NULL,
    resolution_label = NULL,
    print_plot = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for plot_subclass_fragment_features", call. = FALSE)
  }
  value <- match.arg(value)
  feature_select <- match.arg(feature_select)
  sort_y <- match.arg(sort_y)
  subclasses_ignore_case <- isTRUE(subclasses_ignore_case)
  max_features <- suppressWarnings(as.integer(max_features)[1L])
  if (!is.finite(max_features) || max_features < 1L) max_features <- 50L

  if (is.list(summary) && !inherits(summary, "data.frame") &&
      "by_feature" %in% names(summary)) {
    dt <- data.table::as.data.table(summary$by_feature)
    if (!nrow(dt)) {
      cand <- summary$by_feature_candidates
      if (!is.null(cand) && nrow(data.table::as.data.table(cand)) > 0L) {
        message(
          "plot_subclass_fragment_features: by_feature empty; plotting by_feature_candidates."
        )
        dt <- data.table::as.data.table(cand)
      }
    }
    # If multi-resolution, filter to one resolution
    if ("resolution_label" %in% names(dt)) {
      avail <- unique(dt$resolution_label)
      pick <- if (is.null(resolution_label) || !nzchar(as.character(resolution_label)[1L])) {
        avail[1L]
      } else {
        as.character(resolution_label)[1L]
      }
      if (!pick %in% avail) {
        stop("`resolution_label` = \"", pick, "\" not found. Available: ",
             paste(avail, collapse = ", "), ".", call. = FALSE)
      }
      dt <- dt[resolution_label == pick]
      if (!nrow(dt)) {
        stop("No data for `resolution_label` = \"", pick, "\".", call. = FALSE)
      }
      dt[, resolution_label := NULL]
    }
  } else {
    dt <- data.table::as.data.table(summary)
  }
  if (!nrow(dt)) {
    stop(
      "No rows to plot. Lower min_prevalence in summarize_subclass_fragment_features.",
      call. = FALSE
    )
  }
  req <- c("subclass", "Feature_Type", "mz_key", value)
  miss <- setdiff(req, names(dt))
  if (length(miss)) {
    stop("Table missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  dt <- dt[is.finite(dt[[value]])]
  if (!nrow(dt)) stop("No finite values for `", value, "`.", call. = FALSE)

  dt[, subclass := trimws(as.character(subclass))]

  digits <- suppressWarnings(as.integer(label_digits)[1L])
  if (!is.finite(digits) || digits < 0L) digits <- 3L

  vn <- value
  dt[, fill_plot := .SD[[1L]], .SDcols = vn]
  sc_rank <- unique(trimws(as.character(subclasses)))
  sc_rank <- sc_rank[!is.na(sc_rank) & nzchar(sc_rank)]

  if (length(sc_rank)) {
    if (subclasses_ignore_case) {
      sc_key <- toupper(sc_rank)
      dt_rank <- dt[toupper(subclass) %chin% sc_key]
    } else {
      dt_rank <- dt[subclass %chin% sc_rank]
    }
    if (!nrow(dt_rank)) {
      lev <- sort(unique(dt$subclass))
      nlev <- length(lev)
      lev_q <- if (nlev <= 25L) lev else c(head(lev, 22L), paste0("(+", nlev - 22L, " more)"))

      cov_hints <- character()
      summ_l <- is.list(summary) && !inherits(summary, "data.frame") &&
        ("by_feature" %in% names(summary)) && ("subclass_coverage" %in% names(summary))
      if (summ_l) {
        scov <- data.table::as.data.table(summary[["subclass_coverage"]])
        mp_thresh <- NA_real_
        if (is.list(summary$params) && !is.null(summary$params$min_prevalence)) {
          mp_thresh <- suppressWarnings(as.numeric(summary$params$min_prevalence)[1L])
        }
        if (nrow(scov) && "subclass" %in% names(scov)) {
          for (s in sc_rank) {
            rw <- scov[trimws(as.character(subclass)) == s]
            if (!nrow(rw) && subclasses_ignore_case) {
              rw <- scov[toupper(trimws(as.character(subclass))) == toupper(trimws(s))]
            }
            if (nrow(rw) && rw$n_lab_with_features[1L] > 0L) {
              mxp <- if ("max_prevalence" %in% names(rw)) {
                suppressWarnings(as.numeric(rw$max_prevalence[1L]))
              } else {
                NA_real_
              }
              if (!is.finite(mxp)) mxp <- NA_real_
              line <- paste0(
                "`", s, "`: n_lab_with_features=", rw$n_lab_with_features[1L],
                " but absent from $by_feature (all fragments below min_prevalence)."
              )
              if (is.finite(mxp)) {
                line <- paste0(line, " Best prevalence~", format(mxp, digits = 3), ".")
              }
              if (is.finite(mxp) && is.finite(mp_thresh)) {
                sug <- max(0.05, min(mxp * 0.995, mp_thresh - 1e-6))
                line <- paste0(
                  line, " Try min_prevalence <= ", format(sug, digits = 3), "."
                )
              }
              cov_hints <- c(cov_hints, line)
            }
          }
        }
      }

      if (length(cov_hints)) {
        ic_part <- paste0(
          "\n",
          paste0("  - ", cov_hints, collapse = "\n"),
          "\n  (PFCA vs PFCAs are different labels if both exist in metadata.)"
        )
      } else {
        sug_lines <- character(length(sc_rank))
        for (i in seq_along(sc_rank)) {
          ag <- agrep(sc_rank[i], lev, ignore.case = TRUE, max.distance = 0.25, value = TRUE)
          sug_lines[i] <- if (length(ag)) ag[[1L]] else NA_character_
        }
        sug_ok <- sug_lines[!is.na(sug_lines) & nzchar(sug_lines)]
        sug_part <- if (length(sug_ok)) {
          paste0(" Did you mean: ", paste(unique(sug_ok), collapse = ", "), "?")
        } else {
          ""
        }
        ic_part <- if (!subclasses_ignore_case) {
          paste0(" Use subclasses_ignore_case=TRUE for case-insensitive match.", sug_part)
        } else {
          sug_part
        }
      }
      stop(
        "No rows match subclasses (", paste(sc_rank, collapse = ", "), ").\n",
        "  Present: ", paste(lev_q, collapse = ", "), ".",
        ic_part,
        call. = FALSE
      )
    }
  } else {
    dt_rank <- dt
  }

  if (feature_select == "global_max") {
    feat_score <- dt_rank[, .(sc = max(fill_plot, na.rm = TRUE)), by = .(Feature_Type, mz_key)]
    data.table::setorder(feat_score, -sc, Feature_Type, mz_key)
    topk <- head(feat_score, max_features)
    topk[, feat_ord := seq_len(.N)]
    topk <- topk[, .(Feature_Type, mz_key, feat_ord)]
  } else {
    n_sc <- data.table::uniqueN(dt_rank$subclass)
    k_each <- max(1L, as.integer(ceiling(max_features / max(1L, n_sc))))
    feat_union <- data.table::data.table(Feature_Type = character(), mz_key = numeric())
    repeat {
      keys_per <- dt_rank[
        ,
        head(.SD[order(-fill_plot, Feature_Type, mz_key)], k_each),
        by = subclass
      ]
      feat_union <- unique(keys_per[, .(Feature_Type, mz_key)])
      if (nrow(feat_union) <= max_features || k_each <= 1L) break
      k_each <- k_each - 1L
    }
    glob_sc_full <- dt[, .(sc = max(fill_plot, na.rm = TRUE)), by = .(Feature_Type, mz_key)]
    topk <- glob_sc_full[feat_union, on = c("Feature_Type", "mz_key"), nomatch = NULL]
    if (!nrow(topk)) {
      stop("Internal error: empty union feature table.", call. = FALSE)
    }
    data.table::setorder(topk, -sc, Feature_Type, mz_key)
    if (nrow(topk) > max_features) topk <- head(topk, max_features)
    topk[, feat_ord := seq_len(.N)]
    topk <- topk[, .(Feature_Type, mz_key, feat_ord)]
  }

  dt <- dt[topk, on = c("Feature_Type", "mz_key"), nomatch = NULL]

  if (isTRUE(collapse_mz_label)) {
    dt[, mz_bin := round(mz_key, digits)]
    dt <- dt[
      ,
      .(
        fill_plot = max(fill_plot, na.rm = TRUE),
        mz_key = stats::median(mz_key),
        feat_ord = min(feat_ord, na.rm = TRUE)
      ),
      by = .(subclass, Feature_Type, mz_bin)
    ]
    dt[, feat_label := paste0(
      Feature_Type, " ",
      formatC(mz_bin, format = "f", digits = digits)
    )]
    if (sort_y == "selection") {
      ord_tab <- dt[, .(sc = max(fill_plot, na.rm = TRUE), fo = min(feat_ord)), by = feat_label]
      data.table::setorder(ord_tab, fo, -sc, feat_label)
      lvl <- unique(ord_tab$feat_label)
      dt[, feat_label := factor(feat_label, levels = rev(lvl))]
    } else {
      uniq <- unique(dt[, .(feat_label, mz_plot = mz_bin, ft = Feature_Type)])
      if (sort_y == "mz_asc") {
        data.table::setorder(uniq, mz_plot, ft, feat_label)
      } else {
        data.table::setorder(uniq, -mz_plot, ft, feat_label)
      }
      lvl <- uniq$feat_label[!duplicated(uniq$feat_label)]
      dt[, feat_label := factor(feat_label, levels = lvl)]
    }
    dt[, mz_bin := NULL]
  } else {
    dt[, base_label := paste0(
      Feature_Type, " ",
      formatC(mz_key, format = "f", digits = digits)
    )]
    dt[, n_mz_in_base := uniqueN(mz_key), by = base_label]
    dt[, feat_label := base_label]
    dt[n_mz_in_base > 1L, feat_label := paste0(
      base_label, " (", formatC(mz_key, format = "fg", digits = max(8L, digits + 2L)), ")"
    )]
    dt[, c("base_label", "n_mz_in_base") := NULL]

    ul <- unique(dt[
      ,
      .(
        feat_label = data.table::first(feat_label),
        feat_ord = data.table::first(feat_ord),
        mz_key = data.table::first(mz_key),
        Feature_Type = data.table::first(Feature_Type)
      ),
      by = .(Feature_Type, mz_key)
    ])
    if (sort_y == "selection") {
      data.table::setorder(ul, feat_ord, Feature_Type, mz_key)
      lvl <- unique(ul$feat_label)
      dt[, feat_label := factor(feat_label, levels = rev(lvl))]
    } else if (sort_y == "mz_asc") {
      data.table::setorder(ul, mz_key, Feature_Type)
      lvl <- unique(ul$feat_label)
      dt[, feat_label := factor(feat_label, levels = lvl)]
    } else {
      data.table::setorder(ul, -mz_key, Feature_Type)
      lvl <- unique(ul$feat_label)
      dt[, feat_label := factor(feat_label, levels = lvl)]
    }
  }

  ttl <- if (is.null(title) || !nzchar(as.character(title)[1L])) {
    paste0("Subclass fragment features (", vn, ")")
  } else {
    as.character(title)[1L]
  }

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = subclass, y = feat_label, fill = fill_plot)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey90") +
    ggplot2::labs(x = "Subclass", y = "Feature (type m/z)", fill = vn, title = ttl) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8)
    ) +
    ggplot2::coord_cartesian(expand = FALSE)

  if (isTRUE(print_plot)) print(p)
  invisible(p)
}
