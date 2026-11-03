#' Assign Schymanski confidence levels to annotated features
#'
#' @description
#' Implements the Schymanski (2014) confidence level scheme for non-target
#' screening (NTS) annotations. Each feature in `@ms1_annotations` receives a
#' `confidence_level` label (L1-L5) driven by the **provenance of the reference
#' library** behind its best hit:
#'
#' - **L1 (Confirmed Structure):** In-house reference-standard library match.
#' - **L2a (Probable - Online library):** Online spectral-library match.
#' - **L2b (Probable - Predicted library):** In-silico / prediction library match.
#' - **L3 (Tentative Candidates):** Sub-divided into 3a (Network propagation + RT),
#'   3b (Network propagation, no RT), 3c (Homologue), 3d (Transformation),
#'   3e (Multi-method consensus).
#' - **L4 (Unique Formula):** Isotope pattern / formula evidence.
#' - **L5 (Exact Mass):** All remaining features.
#'
#' The library tiers are assigned **purely by `db_type`** — in-house -> L1,
#' online -> L2a, prediction -> L2b — and are **not** gated on retention-time
#' agreement or MS2 cosine similarity. A feature whose best match came from an
#' in-house library is L1 even without RT confirmation; an online match is L2a
#' even when no measured RT is available or the MS2 cosine is low. Only genuine
#' library matches are tiered this way: `@db_annotations` rows whose
#' `annotation_method` marks a *structural* result (network propagation,
#' homologue series, transformation) keep their L3/L4 tier even if they happen
#' to carry a `db_type` value. RT agreement still splits the
#' network-propagation tier (L3a vs L3b). Each feature receives the **highest**
#' applicable level when multiple evidence sources are available. DB type is
#' read from `@db_annotations$db_type` (or auto-detected from the `Library_ID` /
#' `ref_id` pattern).
#'
#' @param project An [NTSProject] with populated `@db_annotations` and
#'   `@ms1_annotations`.
#' @param ms2_cosine_threshold **Deprecated / unused.** Kept so existing calls
#'   keep working. The online tier (L2a) is no longer gated on MS2 cosine
#'   similarity; the library tiers are assigned purely by `db_type`.
#' @param propagation_max_level Maximum BFS propagation level for Level 3a/3b.
#'   Default `3`.
#' @param rt_tol_sec RT tolerance in seconds used to split the network
#'   propagation tier: a propagated feature earns L3a ("with RT") only when the
#'   observed-vs-library RT error is within this tolerance, otherwise L3b
#'   ("no RT"). RT is never used as a hard candidate filter, and no longer
#'   affects the library tiers (L1 / L2a / L2b). Default `30`.
#' @param ... Reserved for future parameters.
#'
#' @return Updated [NTSProject] with `confidence_level` and `rt_confidence`
#'   columns appended to `@ms1_annotations`.
#'
#' @seealso [merge_annotation_layers], [count_annotation_layers],
#'   [annotate_by_dual_network], [build_homologue_network],
#'   [annotate_transformations]
#'
#' @export
assign_confidence_levels <- function(project,
                                     ms2_cosine_threshold = 0.7,
                                     propagation_max_level = 3L,
                                     rt_tol_sec = 30,
                                     ...) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject.", call. = FALSE)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }

  ms1 <- as.data.frame(get_ms1_annotations(project), stringsAsFactors = FALSE)
  if (!nrow(ms1)) {
    ft <- get_feature_table(project)
    if (!nrow(ft) || !"name" %in% names(ft)) {
      warning("Empty ms1_annotations and feature_table; returning project unchanged.",
              call. = FALSE)
      return(invisible(project))
    }
    ms1 <- as.data.frame(ft, stringsAsFactors = FALSE)
  }
  if (!"name" %in% names(ms1)) {
    stop("ms1_annotations must contain column `name`.", call. = FALSE)
  }

  db <- .nts_get_project_db_annotations(project)

  # Per-feature evidence summary with the 10-level priority scheme
  dt <- .nts_confidence_build_feature_table(
    ms1 = ms1,
    db = db,
    ms2_cosine_threshold = ms2_cosine_threshold,
    propagation_max_level = propagation_max_level,
    rt_tol_sec = rt_tol_sec,
    project = project
  )

  dt <- .nts_confidence_assign_levels(dt)

  ms1 <- .nts_confidence_write_to_ms1(ms1, dt)

  project <- set_ms1_annotations(project, .nts_ms1_sanitize_and_order(ms1))
  invisible(project)
}
# ============================================================================
# Internal helpers for assign_confidence_levels
# ============================================================================

# 10-level priority map: lower = higher confidence
# The three library tiers are keyed on db_type (provenance) only:
#  1 = L1  (in-house library)
#  2 = L2a (online library)
#  3 = L2b (prediction library)
#  4 = L3a (propagation + RT)
#  5 = L3b (propagation - RT)
#  6 = L3c (homologue)
#  7 = L3d (transformation)
#  8 = L3e (multi-method)
#  9 = L4  (isotope)
# 10 = L5  (rest)
.nts_confidence_priority <- c(
  L1  = 1L,
  L2a = 2L,
  L2b = 3L,
  L3a = 4L,
  L3b = 5L,
  L3c = 6L,
  L3d = 7L,
  L3e = 8L,
  L4  = 9L,
  L5  = 10L
)

.nts_confidence_label <- function(priority) {
  factor(priority,
    levels = 1L:10L,
    labels = c("L1", "L2a", "L2b", "L3a", "L3b", "L3c", "L3d", "L3e",
               "L4", "L5")
  )
}


#' Detect database type from db_type column or Library_ID pattern
#'
#' Resolution order per row:
#' 1. an explicit `db_type` column (values other than inhouse/online/prediction
#'    and `NA` are treated as unresolved);
#' 2. `Lab_ID` pattern matching for any still-unresolved row;
#' 3. default to `"online"`.
#'
#' When the `db_type` column is entirely absent from `db`, a one-shot warning is
#' emitted (once per R session): a library that is in fact in-house but carries
#' no `db_type` and has non-matching IDs would otherwise be silently demoted to
#' the online tiers, making L1/L2 unreachable. Set the provenance with
#' [set_db_type] (or rebuild the library) to silence it.
#' @keywords internal
#' @noRd
.nts_confidence_db_type <- function(db) {
  n <- nrow(db)
  if (!n) return(rep("online", 0L))

  out <- rep(NA_character_, n)

  # 1) explicit db_type column (non-empty, recognised values only)
  if ("db_type" %in% names(db)) {
    dt <- tolower(trimws(as.character(db$db_type)))
    dt[!dt %in% c("inhouse", "online", "prediction")] <- NA_character_
    out <- dt
  } else {
    # Column absent entirely -> cannot know provenance. Warn once per session:
    # an in-house library without db_type silently becomes "online".
    if (!isTRUE(getOption("NTSminer.db_type_warned", FALSE))) {
      warning(
        "db_annotations has no `db_type` column; database provenance is being ",
        "inferred from `Lab_ID` patterns only and unmatched rows default to ",
        "'online'. If this is an in-house standard library, set `db_type` ",
        "once with set_db_type() (or rebuild it) so L1/L2 stay reachable.",
        call. = FALSE
      )
      options(NTSminer.db_type_warned = TRUE)
    }
  }

  # 2) pattern match for rows still unresolved
  needs <- is.na(out)
  if (any(needs)) {
    lib_col <- if ("Library_ID" %in% names(db)) "Library_ID"
               else if ("ref_id" %in% names(db)) "ref_id"
               else if ("Lab.ID" %in% names(db)) "Lab.ID"
               else if ("Lab_ID" %in% names(db)) "Lab_ID"
               else NULL
    if (!is.null(lib_col)) {
      ids <- tolower(as.character(db[[lib_col]]))
      ids[is.na(ids)] <- ""
      out[needs & grepl("inhouse|in_house|nts_inhouse", ids)] <- "inhouse"
      out[needs & grepl("predicted|prediction|pred_", ids)] <- "prediction"
    }
  }

  # 3) default remainder to online
  out[is.na(out)] <- "online"
  out
}


#' Detect whether a db_annotation has RT *agreement* evidence
#'
#' "Has RT" means the library RT is present AND the observed-vs-library RT
#' error is within `rt_tol_sec`. Merely having a library RT value is NOT
#' enough: the `L3a` ("propagation + RT") tier relies on RT *corroborating* the
#' identification, so a match whose RT is far off tolerance must fall to `L3b`.
#' (The library tiers L1/L2a/L2b are provenance-only and no longer consult RT.
#' Candidate generation itself is m/z-only; RT is used as a scoring/confirmation
#' axis, never a hard candidate filter.)
#' @keywords internal
#' @noRd
.nts_confidence_has_rt_vec <- function(db, rt_tol_sec = 30) {
  n <- nrow(db)
  if (!n) return(rep(FALSE, 0L))
  if (!"lib_rt" %in% names(db)) return(rep(FALSE, n))
  # Coerce element-wise: `is.numeric()` on the whole column returns a single
  # logical, so a character `lib_rt` (common after MSP/CSV import) would make
  # every row "no RT". Parse per element instead so numeric-looking strings
  # still count as RT evidence.
  lib_rt_num <- suppressWarnings(as.numeric(db$lib_rt))
  has_lib_rt <- !is.na(lib_rt_num) & is.finite(lib_rt_num)
  # Without a measured RT error we cannot confirm agreement -> not "has RT".
  if (!"rt_error_s" %in% names(db)) {
    return(rep(FALSE, n))
  }
  rt_err <- suppressWarnings(as.numeric(db$rt_error_s))
  tol <- suppressWarnings(as.numeric(rt_tol_sec)[1L])
  if (!is.finite(tol) || tol <= 0) tol <- 30
  has_lib_rt & !is.na(rt_err) & is.finite(rt_err) & rt_err <= tol
}


#' Build per-feature evidence summary with the 10-level priority scheme
#' @keywords internal
#' @noRd
.nts_confidence_build_feature_table <- function(ms1,
                                                db,
                                                ms2_cosine_threshold,
                                                propagation_max_level,
                                                rt_tol_sec = 30,
                                                project = NULL) {
  n <- nrow(ms1)
  dt <- data.table::data.table(
    name      = trimws(as.character(ms1$name)),
    priority  = rep(10L, n)   # default L5
  )

  # ---- MS1-derived evidence ----

  # L3c: Homologue series (priority 6)
  if ("homologous_series_id" %in% names(ms1)) {
    hs <- trimws(as.character(ms1$homologous_series_id))
    hs_mask <- !is.na(hs) & nzchar(hs) & hs != "NA"
    dt[hs_mask, priority := pmin(priority, 6L)]
  }
  if ("kmd_series_id" %in% names(ms1)) {
    ks <- trimws(as.character(ms1$kmd_series_id))
    ks_mask <- !is.na(ks) & nzchar(ks) & ks != "NA"
    dt[ks_mask, priority := pmin(priority, 6L)]
  }
  # L3c: Homologue provenance written to the merged ms1 layer. The network path
  # reaches confidence via merge_annotation_layers(write_to_ms1 = TRUE) /
  # apply_annotation_master_to_ms1(), which store the contributing sources in
  # `ms1$annotation_method` (e.g. "Spectral_DB; Homologue"). The db-side L3c
  # branch further below reads `@db_annotations$annotation_method`, which the
  # homologue network layer never writes -- so without this ms1-side read the
  # network propagation could not reach L3c at all (P1-E4).
  if ("annotation_method" %in% names(ms1)) {
    am_ms1 <- trimws(as.character(ms1$annotation_method))
    am_ms1[is.na(am_ms1)] <- ""
    hm_ms1 <- nzchar(am_ms1) & grepl("Homologue", am_ms1, fixed = TRUE)
    if (any(hm_ms1)) {
      dt[hm_ms1, priority := pmin(priority, 6L)]
    }
  }
  # L3c via the split homologue layer (Plan B). merge_ms1_homologue_layers()
  # writes `homologue_ref_id` / `homologue_label` into the `network` step of
  # @ms1_annotations (registered in .nts_ms1_groups()); they surface in the flat
  # ms1 view here. Either column being non-empty means the feature carries a
  # homologue-network identity -> L3c (priority 6).
  for (hcol in c("homologue_ref_id", "homologue_label")) {
    if (hcol %in% names(ms1)) {
      hv <- trimws(as.character(ms1[[hcol]]))
      hv_mask <- !is.na(hv) & nzchar(hv) & hv != "NA"
      if (any(hv_mask)) {
        dt[hv_mask, priority := pmin(priority, 6L)]
      }
    }
  }

  # L3d: Transformation evidence (priority 7)
  if ("transformed_formula" %in% names(ms1)) {
    tf <- trimws(as.character(ms1$transformed_formula))
    tf_mask <- !is.na(tf) & nzchar(tf) & tf != "NA"
    dt[tf_mask, priority := pmin(priority, 7L)]
  }

  # L4: Isotope pattern (priority 9)
  has_iso_score <- FALSE
  if (is.data.frame(db) && nrow(db) && "isotope_cosine" %in% names(db)) {
    feat_col_db <- if ("name" %in% names(db)) "name"
                   else if ("feature_id" %in% names(db)) "feature_id"
    if (!is.null(feat_col_db)) {
      db_feat_iso <- trimws(as.character(db[[feat_col_db]]))
      iso_cos <- suppressWarnings(as.numeric(db$isotope_cosine))
      iso_valid <- !is.na(db_feat_iso) & nzchar(db_feat_iso) &
        !is.na(iso_cos) & is.finite(iso_cos) & iso_cos >= 0.9
      if (any(iso_valid)) {
        l4_feats <- unique(db_feat_iso[iso_valid])
        dt[name %in% l4_feats, priority := pmin(priority, 9L)]
        has_iso_score <- TRUE
      }
    }
  }
  if (!has_iso_score &&
      "isotope_kind" %in% names(ms1) && "isotope_c13_n" %in% names(ms1)) {
    ik <- ms1$isotope_kind
    ic <- ms1$isotope_c13_n
    l4_mask <- !is.na(ik) & nzchar(as.character(ik)) &
      !is.na(ic) & is.finite(as.numeric(ic))
    dt[l4_mask, priority := pmin(priority, 9L)]
  }

  # ---- DB-derived evidence ----
  # Per-feature RT-agreement map (populated from the db rows below). The MS2-tier
  # network branch at the end reuses it to split L3a (RT agreement) from L3b.
  has_rt_by_name <- NULL
  if (is.data.frame(db) && nrow(db)) {
    feat_col <- if ("name" %in% names(db)) "name"
                else if ("feature_id" %in% names(db)) "feature_id"
                else NULL

    if (!is.null(feat_col)) {
      db_feat <- trimws(as.character(db[[feat_col]]))
      db_valid <- !is.na(db_feat) & nzchar(db_feat)

      db_type  <- .nts_confidence_db_type(db)
      has_rt   <- .nts_confidence_has_rt_vec(db, rt_tol_sec = rt_tol_sec)
      has_rt_by_name <- tapply(has_rt, db_feat, function(v) any(v, na.rm = TRUE))

      # `annotation_method` separates a genuine spectral-library match (blank, or
      # "Spectral_DB") from a *structural* annotation produced by the network /
      # series engines (network propagation, homologue series, transformation).
      # Structural rows can carry a `db_type` value only incidentally, so they
      # must NOT be re-labelled as library hits -- otherwise e.g. a propagated
      # feature would surface as L2a. A row keeps a library tier only when it is
      # not structural-only.
      am <- if ("annotation_method" %in% names(db)) {
        trimws(as.character(db$annotation_method))
      } else {
        rep(NA_character_, nrow(db))
      }
      am[is.na(am)] <- ""
      am_structural <- grepl("Dual_Network", am, fixed = TRUE) |
        grepl("Propagation",    am, fixed = TRUE) |
        grepl("Network_",       am, fixed = TRUE) |
        grepl("Homologue",      am, fixed = TRUE) |
        grepl("Transformation", am, fixed = TRUE) |
        grepl("PFAS_Series",    am, fixed = TRUE)
      am_spectral <- grepl("Spectral_DB", am, fixed = TRUE)
      lib_row     <- db_valid & !(am_structural & !am_spectral)

      # ---- Library tiers: provenance ONLY ----------------------------------
      # in-house -> L1, online -> L2a, prediction -> L2b. Assigned purely by
      # `db_type`; RT agreement and MS2 cosine are deliberately NOT gates here
      # (RT only splits the propagation tier L3a/L3b below).
      for (fn in unique(db_feat[lib_row])) {
        types <- unique(db_type[lib_row & db_feat == fn])
        for (tp in types) {
          if (tp == "inhouse") {
            dt[name == fn, priority := pmin(priority, 1L)]   # L1
          } else if (tp == "online") {
            dt[name == fn, priority := pmin(priority, 2L)]   # L2a
          } else if (tp == "prediction") {
            dt[name == fn, priority := pmin(priority, 3L)]   # L2b
          }
        }
      }

      # L3a / L3b: structural network propagation (split by RT). The token set
      # here MUST mirror `am_structural` above: the merge layer emits
      # `annotation_method = "Network_DF"/"Network_NL"/"Network_MD"`, while the
      # network engines emit `Dual_Network`/`Propagation`. Omitting "Network_"
      # let a `Network_DF` structural row slip past every tier and fall to L5
      # (P1-D2).
      dn_mask <- db_valid & !is.na(am) & nzchar(am) &
        (grepl("Dual_Network", am, fixed = TRUE) |
          grepl("Propagation", am, fixed = TRUE) |
          grepl("Network_",    am, fixed = TRUE))
      if (any(dn_mask)) {
        pl <- rep(NA_integer_, nrow(db))
        if ("propagation_level" %in% names(db)) {
          pl <- suppressWarnings(as.integer(db$propagation_level))
        } else if ("Propagation_Level" %in% names(db)) {
          pl <- suppressWarnings(as.integer(db[["Propagation_Level"]]))
        }
        dn_valid <- dn_mask & (!is.na(pl) & pl >= 1L & pl <= propagation_max_level)
        if (any(dn_valid)) {
          dn_feats <- unique(db_feat[dn_valid])
          for (fn in dn_feats) {
            fn_row <- dn_valid & db_feat == fn
            fn_rt <- any(has_rt[fn_row], na.rm = TRUE)
            if (fn_rt) dt[name == fn, priority := pmin(priority, 4L)]  # L3a
            else        dt[name == fn, priority := pmin(priority, 5L)]  # L3b
          }
        }
      }

      # L3c: Homologue from db (priority 6)
      hm_mask <- db_valid & !is.na(am) & nzchar(am) &
        grepl("Homologue", am, fixed = TRUE)
      if (any(hm_mask)) {
        l3c_feats <- unique(db_feat[hm_mask])
        dt[name %in% l3c_feats, priority := pmin(priority, 6L)]
      }

      # L3d: Transformation from db (priority 7)
      tr_mask <- db_valid & !is.na(am) & nzchar(am) &
        grepl("Transformation", am, fixed = TRUE)
      if (any(tr_mask)) {
        l3d_feats <- unique(db_feat[tr_mask])
        dt[name %in% l3d_feats, priority := pmin(priority, 7L)]
      }
      if ("rule_name" %in% names(db)) {
        rn <- trimws(as.character(db$rule_name))
        rn_mask <- db_valid & !is.na(rn) & nzchar(rn) & rn != "NA"
        if (any(rn_mask)) {
          l3d_feats2 <- unique(db_feat[rn_mask])
          dt[name %in% l3d_feats2, priority := pmin(priority, 7L)]
        }
      }

      # L3e: Multi-method consensus (priority 8)
      mm_mask <- db_valid & !is.na(am) & nzchar(am) &
        grepl(";", am, fixed = TRUE)
      if (any(mm_mask)) {
        l3e_feats <- unique(db_feat[mm_mask])
        dt[name %in% l3e_feats, priority := pmin(priority, 8L)]
      }
    }
  }

  # Also read ms1 annotation_level for propagation L3 evidence. The dual-network
  # engine writes the spelled-out form "Network Propagated: <id> (Dice x.xxx)"
  # (see annotate_by_dual_network), whereas the merged ms1 layer uses
  # "Network_DF/NL/MD"/"Dual_Network". Matching on the "Network" stem covers all
  # of them; the old "Network_" (underscore) token missed "Network Propagated"
  # entirely, so that string never reached L3b (P1-D1).
  if ("annotation_level" %in% names(ms1)) {
    al <- trimws(as.character(ms1$annotation_level))
    al[is.na(al)] <- ""
    dn_al_mask <- nzchar(al) & (grepl("Network", al, fixed = TRUE) |
      grepl("Propagation", al, fixed = TRUE))
    if (any(dn_al_mask)) {
      # Conservative: no RT info from ms1 alone -> L3b
      dt[dn_al_mask, priority := pmin(priority, 5L)]
    }
  }

  # Structural provenance written by merge_annotation_layers(write_to_ms1 = TRUE)
  # lands in `ms1$annotation_method` ("Network_DF"/"Network_NL"/"Network_MD") +
  # `ms1$master_propagation_level` -- NOT in `ms1$annotation_level`, nor in
  # `@db_annotations$annotation_method` that the db-side branch above reads. Read
  # the ms1 method column here so the merge -> write-back -> confidence main path
  # can actually reach L3 (P1-D2). ms1 carries no RT-agreement axis, so split
  # conservatively to L3b (priority 5); a level present but out of the accepted
  # domain is rejected, mirroring the db-side `propagation_max_level` guard.
  if ("annotation_method" %in% names(ms1)) {
    am_ms1n <- trimws(as.character(ms1$annotation_method))
    am_ms1n[is.na(am_ms1n)] <- ""
    net_ms1 <- nzchar(am_ms1n) & (grepl("Network_", am_ms1n, fixed = TRUE) |
      grepl("Dual_Network", am_ms1n, fixed = TRUE) |
      grepl("Propagation",    am_ms1n, fixed = TRUE))
    if (any(net_ms1)) {
      pl_ms1 <- if ("master_propagation_level" %in% names(ms1)) {
        suppressWarnings(as.integer(ms1$master_propagation_level))
      } else {
        rep(NA_integer_, nrow(ms1))
      }
      pl_ok <- is.na(pl_ms1) |
        (pl_ms1 >= 1L & pl_ms1 <= propagation_max_level)
      net_valid <- net_ms1 & pl_ok
      if (any(net_valid)) {
        dt[net_valid, priority := pmin(priority, 5L)]
      }
    }
  }

  # Structural network evidence split onto the MS2 tier (Plan B). The split write
  # path merge_ms2_network_layers() lands `net_df_ref_id` / `net_nl_ref_id` /
  # `net_md_ref_id` on @ms2_annotations. Any of them non-empty means the feature
  # carries an MS2-derived structural identity -> L3b (priority 5); when the
  # feature ALSO has RT agreement (reused from the db side, `has_rt_by_name`),
  # promote to L3a (priority 4). Read via the existing get_ms2_annotations().
  if (!is.null(project) && methods::is(project, "NTSProject")) {
    ms2 <- tryCatch(
      as.data.frame(get_ms2_annotations(project), stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.data.frame(ms2) && nrow(ms2) && "name" %in% names(ms2)) {
      ref_cols <- intersect(
        c("net_df_ref_id", "net_nl_ref_id", "net_md_ref_id"),
        names(ms2)
      )
      if (length(ref_cols)) {
        nm2 <- trimws(as.character(ms2$name))
        present <- rep(FALSE, nrow(ms2))
        for (cc in ref_cols) {
          v <- trimws(as.character(ms2[[cc]]))
          present <- present | (!is.na(v) & nzchar(v))
        }
        net_feats <- unique(nm2[present])
        net_feats <- net_feats[!is.na(net_feats) & nzchar(net_feats)]
        if (length(net_feats)) {
          dt[name %in% net_feats, priority := pmin(priority, 5L)]  # L3b
          if (!is.null(has_rt_by_name) && length(has_rt_by_name)) {
            rt_feats <- names(has_rt_by_name)[has_rt_by_name]
            rt_hit <- intersect(net_feats, rt_feats)
            if (length(rt_hit)) {
              dt[name %in% rt_hit, priority := pmin(priority, 4L)]  # L3a
            }
          }
        }
      }
    }
  }

  dt[]
}


#' Assign Schymanski confidence level from priority
#' @keywords internal
#' @noRd
.nts_confidence_assign_levels <- function(dt) {
  dt[, confidence_level := as.character(.nts_confidence_label(priority))]
  dt[]
}


#' Write confidence_level and rt_confidence back to ms1_annotations
#' @keywords internal
#' @noRd
.nts_confidence_write_to_ms1 <- function(ms1, dt) {
  ms1$confidence_level <- NA_character_
  ms1$rt_confidence   <- "unavailable"

  idx <- match(as.character(ms1$name), as.character(dt$name))
  valid <- !is.na(idx)
  if (any(valid)) {
    ms1$confidence_level[valid] <- dt$confidence_level[idx[valid]]
  }

  ms1
}
