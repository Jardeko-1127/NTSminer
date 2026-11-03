# =============================================================================
# 07_annotation_rank.R —— 注释合并 / 优先级排序（只读派生）
# -----------------------------------------------------------------------------
# 目的
#   把一条龙跑出来的**四路注释证据**合并成「一张候选排序表 + 每特征一个最优注释」：
#     库命中(L1/L2a/L2b) + MetFrag 预测(L2c) + 网络/同系物传播(L3) + 自证结构(阶梯/诊断离子)
#
# 设计（与置信度 10 档体系保持一致）
#   主键 = **识别来源**的置信度档；★ 看的是「怎么识别出来的」，不是「化合物属于哪个库」：
#   in-house 化合物若只被 MetFrag 预测命中（库并未命中该特征），档位是 L2c 而不是 L1。
#     L1 > L2a > L2b > L2c > L3a > L3b > L3c > L3d > L3e > L4 > L5
#   （与 .nts_confidence_priority 的相对顺序一致，唯一新增的是在 L2b 与 L3a 之间插 L2c。）
#
#   叠加 = **独立证据族**有限晋升（半档 "*"），**不做无界投票** —— 来源之间高度相关：
#   MetFrag 的候选库往往就是库并集、网络传播身份又全部源自库命中种子。
#     F1 实测参考谱 : in-house(L1) / online(L2a)
#     F2 预测       : insilico 预测谱(L2b) / MetFrag(L2c)
#     F3 自证结构   : CF2 阶梯 / 诊断离子（特征级，只做平票与子类对照，不用于识别）
#     F4 关系推断   : 同系物 ±kCF2 锚定 / 特征网络传播身份（ev_tier 0 不计，见下）
#
#   复核清单（merge_conflicts）四条 crisp 判据：
#     mf_vs_lib_differs      实测谱库(inhouse/online)首选 与 MetFrag Top-1 不同化合物
#     mf_vs_insilico_differs 仅预测库时，insilico 谱库首选 与 MetFrag Top-1 不同化合物
#     net_ident_differs      网络传播身份 != 最优注释，且最优档不是实测谱档
#     tied_alternatives      同一最优档并列多个化合物（证据分不开）
#   conflict_status: none / resolved_by_lib（最终档 L1/L2a，实测谱已定案）/ open（需人工看）
#
# 只读：不写工程、不写 ckpt、不写 @network_results。
# =============================================================================


# ---------------------------------------------------------------------------
# 内部小工具
# ---------------------------------------------------------------------------

#' Missing-or-empty test for rank helpers
#' @param x Any vector.
#' @return Logical vector.
#' @keywords internal
#' @noRd
.nts_rank_blank <- function(x) {
  is.na(x) | !nzchar(trimws(as.character(x)))
}

#' @keywords internal
#' @noRd
.nts_rank_num <- function(x) suppressWarnings(as.numeric(x))

#' @keywords internal
#' @noRd
.nts_rank_chr <- function(x) as.character(x)

#' Fetch a column, or a length-`nrow` default when absent
#' @keywords internal
#' @noRd
.nts_rank_col <- function(df, nm, default = NA) {
  if (is.null(df) || !nm %in% names(df)) {
    return(rep(default, nrow(df)))
  }
  df[[nm]]
}

#' Resolve the feature-key column of a table
#'
#' Tables in the pipeline key on `feature`, `name` or `feature_id` depending on
#' the producer; `NTSminer`'s own exports use `name`.
#' @return The column name, or `NULL` when none is present.
#' @keywords internal
#' @noRd
.nts_rank_feature_col <- function(df) {
  for (nm in c("feature", "name", "feature_id")) {
    if (nm %in% names(df)) {
      return(nm)
    }
  }
  NULL
}

#' Confidence ladder used for ranking (ascending = better)
#' @keywords internal
#' @noRd
.nts_rank_ladder <- function() {
  c(
    L1 = 1, L2a = 2, L2b = 3, L2c = 4, L3a = 5, L3b = 6,
    L3c = 7, L3d = 8, L3e = 9, L4 = 10, L5 = 11
  )
}

#' Library set of a `ref_id` (from its `<set>_<lab_id>` prefix)
#' @return `"inhouse"` / `"online"` / `"insilico"`, else `NA`.
#' @keywords internal
#' @noRd
.nts_rank_lib_of <- function(ref_id) {
  s <- sub("_.*$", "", .nts_rank_chr(ref_id))
  s[!grepl("^(inhouse|online|insilico)$", s)] <- NA_character_
  s
}

#' Library set of a candidate, falling back to the `db_type` column
#' @keywords internal
#' @noRd
.nts_rank_lib_of2 <- function(ref_id, db_type) {
  s <- .nts_rank_lib_of(ref_id)
  d <- tolower(.nts_rank_chr(db_type))
  d[!grepl("^(inhouse|online|insilico)$", d)] <- NA_character_
  ifelse(is.na(s), d, s)
}

#' Base confidence level of a library-identification set
#'
#' `inhouse` -> `L1`, `online` -> `L2a`, `insilico` -> `L2b` (matches
#' `.nts_confidence_priority`); anything else -> `NA`.
#' @keywords internal
#' @noRd
.nts_rank_base_level <- function(set) {
  out <- rep(NA_character_, length(set))
  out[set == "inhouse"] <- "L1"
  out[set == "online"] <- "L2a"
  out[set == "insilico"] <- "L2b"
  out
}

#' Read a CSV written by the pipeline (tolerates a UTF-8 BOM)
#' @return A `data.frame`, or `NULL` when the file is absent / unreadable.
#' @keywords internal
#' @noRd
.nts_rank_read_csv <- function(path) {
  if (is.null(path) || length(path) != 1L || is.na(path) || !file.exists(path)) {
    return(NULL)
  }
  tryCatch(
    utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE,
                    fileEncoding = "UTF-8-BOM"),
    error = function(e) {
      utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    }
  )
}


# ---------------------------------------------------------------------------
# 核心：纯数据框
# ---------------------------------------------------------------------------

#' Merge and rank annotation evidence into one best call per feature
#'
#' Pure, side-effect-free core of the annotation merge: takes the candidate
#' long tables produced by a pipeline run in memory and returns four
#' `data.frame`s. Use [merge_annotation_sources] to run it straight off a
#' pipeline product directory.
#'
#' Ranking key is the **confidence level of the identification source**, not the
#' library the compound happens to live in:
#' `L1` (in-house hit) > `L2a` (online hit) > `L2b` (predicted-spectrum hit) >
#' `L2c` (MetFrag only, no library hit) > `L3b`/`L3c` (network propagation /
#' homologue). `L2c` is introduced by this function; the packaged confidence
#' ladder jumps from `L2b` to `L3a`.
#'
#' Corroboration is a **bounded half-rung** marker, never additive voting: a
#' candidate supported by two or more *independent evidence families* gets a
#' trailing `"*"` on `decision_label` (same level, never promoted past it).
#' The families are F1 measured reference spectra (`inhouse` / `online`),
#' F2 prediction (`insilico` hits, MetFrag), and F4 relationship inference
#' (homologue / feature-network propagated identity, or `ev_tier` 1-2 in the
#' candidate table). `ev_tier == 0` is deliberately excluded: it only means
#' "this candidate is itself an L1 hit", which is not independent support.
#'
#' @param candidates Library-hit candidate long table (one row per feature x
#'   candidate). Needs a feature column (`feature` / `name` / `feature_id`) and
#'   `ref_id`; recognises `candidate` (or `Compound.name`), `db_type`,
#'   `score_total`, `score_ms2`, `ev_tier`, `network_support`, `SMILES`.
#'   Column names come from the `ref_id` prefixes, so `ref_id` should look like
#'   `<set>_<lab_id>` with `<set>` in `inhouse` / `online` / `insilico`.
#' @param metfrag_hits Optional MetFrag candidate long table (same feature
#'   column plus `identifier` — the library `ref_id` — and `rank`,
#'   `score_total`, `is_top1`, `Compound.name`).
#' @param network_layers Optional long table of network-propagated identities:
#'   feature column + `ref_id` + `net_kind` (`"DFI"` / `"NL"` / `"MD"` /
#'   `"HOMO"`; default `"NET"`). Homologue identities must be tagged `"HOMO"`
#'   so they can land on `L3c`.
#' @param feature_meta Optional feature-level table (feature column plus
#'   `mz`, `rt`, `confidence_level`, `pfas_ladder_len`, `pfas_subclass_hint`),
#'   used for reported coordinates and the self-evidence structure flags.
#' @param suspect Optional `data.frame` mapping compound ids to names /
#'   InChIKeys (columns `ID` / `Identifier` / `ref_id` / `Lab_ID` plus `Name`
#'   and optionally `InChIKey`). Used to fill blank candidate names and to add
#'   the `inchikey` column.
#' @param verbose If `TRUE`, report progress and the level tally via
#'   [message].
#'
#' @return A named `list` of four `data.frame`s:
#'   \describe{
#'     \item{`candidates`}{long table, one row per feature x candidate, ranked
#'       (`rank_in_feature`, `is_top`, `base_level`, `decision_tier`,
#'       `decision_label`, `n_families`, `family_F1/F2/F4`, `flags`, scores).}
#'     \item{`final`}{one row per feature: the winning annotation plus
#'       `tied_alternatives`, `review_reason`, `conflict_status` and the
#'       per-family best (`fam1_best`, `fam2_best`, `fam4_best`,
#'       `lib_best_ref`, `mf_top1_ref`, `net_ident_ref`, `mf_vs_lib`).}
#'     \item{`conflicts`}{the subset of `final` with a non-empty
#'       `review_reason` (the human review list).}
#'     \item{`summary`}{counts by `decision_label`
#'       (`final_label`, `n_features`).}
#'   }
#'
#' @seealso [merge_annotation_sources] to run this off a pipeline run
#'   directory; [merge_annotation_layers] / [count_annotation_layers] for the
#'   project-slot annotation layers.
#'
#' @examples
#' cand <- data.frame(
#'   feature = c("F1", "F1", "F2"),
#'   ref_id = c("inhouse_L01", "insilico_L07", "online_L03"),
#'   candidate = c("PFOA", "PFHxA isomer", "PFOS"),
#'   db_type = c("inhouse", "prediction", "online"),
#'   score_total = c(0.88, 0.55, 0.77),
#'   stringsAsFactors = FALSE
#' )
#' mf <- data.frame(
#'   name = c("F1", "F2"),
#'   identifier = c("inhouse_L01", "insilico_L09"),
#'   Compound.name = c("PFOA", "PFOS (predicted)"),
#'   score_total = c(120, 90),
#'   rank = c(1L, 1L),
#'   is_top1 = c(TRUE, TRUE),
#'   stringsAsFactors = FALSE
#' )
#' net <- data.frame(
#'   feature = c("F1", "F2"),
#'   ref_id = c("insilico_L07", "insilico_L09"),
#'   net_kind = c("DFI", "HOMO"),
#'   stringsAsFactors = FALSE
#' )
#' out <- rank_annotation_candidates(cand, metfrag_hits = mf, network_layers = net)
#' out$final[, c("name", "candidate", "decision_label", "n_families",
#'               "conflict_status", "review_reason")]
#'
#' @export
rank_annotation_candidates <- function(candidates,
                                       metfrag_hits = NULL,
                                       network_layers = NULL,
                                       feature_meta = NULL,
                                       suspect = NULL,
                                       predicted_fp = NULL,
                                       weight_fp = 0.25,
                                       verbose = TRUE) {
  if (is.null(candidates) || !is.data.frame(candidates)) {
    stop("`candidates` must be a data.frame (one row per feature x candidate).",
         call. = FALSE)
  }
  if (!nrow(candidates)) {
    stop("`candidates` has no rows.", call. = FALSE)
  }
  cand_fcol <- .nts_rank_feature_col(candidates)
  if (is.null(cand_fcol)) {
    stop("`candidates` needs a feature column: `feature`, `name` or `feature_id`.",
         call. = FALSE)
  }
  if (!"ref_id" %in% names(candidates)) {
    stop("`candidates` needs a `ref_id` column.", call. = FALSE)
  }
  say <- function(...) if (isTRUE(verbose)) message(...)

  # Optional Candidate Re-ranking using predicted joint fingerprints
  if (!is.null(predicted_fp) && !"fp_similarity" %in% names(candidates)) {
    candidates <- rescore_candidates_by_joint_fp(
      candidates,
      predicted_fp = predicted_fp,
      weight_fp = weight_fp,
      verbose = verbose
    )
  }

  ## ------------------------------------------------------- 1. 特征级信息表
  if (!is.null(feature_meta) && is.data.frame(feature_meta)) {
    fm_col <- .nts_rank_feature_col(feature_meta)
    if (is.null(fm_col)) {
      stop("`feature_meta` needs a feature column: `feature`, `name` or `feature_id`.",
           call. = FALSE)
    }
    feat <- data.frame(
      name = .nts_rank_chr(feature_meta[[fm_col]]),
      mz = .nts_rank_num(.nts_rank_col(feature_meta, "mz")),
      rt = .nts_rank_num(.nts_rank_col(feature_meta, "rt")),
      confidence_level = .nts_rank_chr(.nts_rank_col(feature_meta, "confidence_level")),
      pfas_ladder_len = .nts_rank_num(.nts_rank_col(feature_meta, "pfas_ladder_len")),
      pfas_subclass_hint = .nts_rank_chr(.nts_rank_col(feature_meta, "pfas_subclass_hint")),
      stringsAsFactors = FALSE
    )
  } else {
    feat <- data.frame(
      name = character(), mz = numeric(), rt = numeric(),
      confidence_level = character(), pfas_ladder_len = numeric(),
      pfas_subclass_hint = character(), stringsAsFactors = FALSE
    )
  }

  ## ------------------------------------------------------- 2. 候选长表并集
  # ident_set = 识别来源库（定档用）；comp_db = 化合物自身所属库（仅 provenance 展示）
  CAND_COLS <- c(
    "feature", "ref_id", "candidate", "src_lib", "src_mf", "src_net",
    "ident_set", "comp_db", "score_lib", "score_lib_ms2",
    "mf_rank", "mf_score", "mf_top1", "ev_tier", "network_support",
    "net_kind", "SMILES", "db_type_raw", "fp_similarity", "score_reranked"
  )
  parts <- list()

  lib_name <- .nts_rank_col(candidates, "candidate")
  if (all(.nts_rank_blank(lib_name)) && "Compound.name" %in% names(candidates)) {
    lib_name <- candidates[["Compound.name"]]
  }
  parts[[length(parts) + 1L]] <- data.frame(
    feature = .nts_rank_chr(candidates[[cand_fcol]]),
    ref_id = .nts_rank_chr(candidates$ref_id),
    candidate = .nts_rank_chr(lib_name),
    src_lib = TRUE, src_mf = FALSE, src_net = FALSE,
    ident_set = .nts_rank_lib_of2(candidates$ref_id,
                                  .nts_rank_col(candidates, "db_type")),
    comp_db = .nts_rank_lib_of(candidates$ref_id),
    score_lib = .nts_rank_num(.nts_rank_col(candidates, "score_total")),
    score_lib_ms2 = .nts_rank_num(.nts_rank_col(candidates, "score_ms2")),
    mf_rank = NA_real_, mf_score = NA_real_, mf_top1 = FALSE,
    ev_tier = .nts_rank_num(.nts_rank_col(candidates, "ev_tier")),
    network_support = .nts_rank_chr(.nts_rank_col(candidates, "network_support")),
    net_kind = NA_character_,
    SMILES = .nts_rank_chr(.nts_rank_col(candidates, "SMILES")),
    db_type_raw = .nts_rank_chr(.nts_rank_col(candidates, "db_type")),
    fp_similarity = .nts_rank_num(.nts_rank_col(candidates, "fp_similarity")),
    score_reranked = .nts_rank_num(.nts_rank_col(candidates, "score_reranked")),
    stringsAsFactors = FALSE
  )

  if (!is.null(metfrag_hits)) {
    if (!is.data.frame(metfrag_hits)) {
      stop("`metfrag_hits` must be a data.frame or NULL.", call. = FALSE)
    }
    if (nrow(metfrag_hits)) {
      mf_col <- .nts_rank_feature_col(metfrag_hits)
      if (is.null(mf_col)) {
        stop("`metfrag_hits` needs a feature column: `feature`, `name` or `feature_id`.",
             call. = FALSE)
      }
      rid <- .nts_rank_chr(.nts_rank_col(metfrag_hits, "identifier"))
      if (all(.nts_rank_blank(rid))) {
        # `identifier` 缺列 / 全空时退回 `ref_id`（有些导出只留一列）
        rid <- .nts_rank_chr(.nts_rank_col(metfrag_hits, "ref_id"))
      }
      mf_feat <- .nts_rank_chr(metfrag_hits[[mf_col]])
      ok <- !.nts_rank_blank(rid) & !.nts_rank_blank(mf_feat) &
        !.nts_rank_blank(.nts_rank_chr(.nts_rank_col(metfrag_hits, "Compound.name")))
      if (any(ok)) {
        mf2 <- metfrag_hits[ok, , drop = FALSE]
        rid <- rid[ok]
        parts[[length(parts) + 1L]] <- data.frame(
          feature = .nts_rank_chr(mf2[[mf_col]]),
          ref_id = rid,
          candidate = .nts_rank_chr(.nts_rank_col(mf2, "Compound.name")),
          src_lib = FALSE, src_mf = TRUE, src_net = FALSE,
          ident_set = NA_character_, comp_db = .nts_rank_lib_of(rid),
          score_lib = NA_real_, score_lib_ms2 = NA_real_,
          mf_rank = .nts_rank_num(.nts_rank_col(mf2, "rank")),
          mf_score = .nts_rank_num(.nts_rank_col(mf2, "score_total")),
          mf_top1 = toupper(.nts_rank_chr(.nts_rank_col(mf2, "is_top1"))) %in%
            c("TRUE", "T", "1"),
          ev_tier = NA_real_, network_support = NA_character_,
          net_kind = NA_character_,
          SMILES = .nts_rank_chr(.nts_rank_col(mf2, "SMILES")),
          db_type_raw = "prediction",
          fp_similarity = .nts_rank_num(.nts_rank_col(mf2, "fp_similarity")),
          score_reranked = .nts_rank_num(.nts_rank_col(mf2, "score_reranked")),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  net_only <- NULL
  if (!is.null(network_layers)) {
    if (!is.data.frame(network_layers)) {
      stop("`network_layers` must be a data.frame or NULL.", call. = FALSE)
    }
    if (nrow(network_layers)) {
      nl_col <- .nts_rank_feature_col(network_layers)
      if (is.null(nl_col)) {
        stop("`network_layers` needs a feature column: `feature`, `name` or `feature_id`.",
             call. = FALSE)
      }
      if (!"ref_id" %in% names(network_layers)) {
        stop("`network_layers` needs a `ref_id` column.", call. = FALSE)
      }
      kind <- .nts_rank_chr(.nts_rank_col(network_layers, "net_kind"))
      kind[.nts_rank_blank(kind)] <- "NET"
      nr <- data.frame(
        feature = .nts_rank_chr(network_layers[[nl_col]]),
        ref_id = .nts_rank_chr(network_layers$ref_id),
        net_kind = kind,
        stringsAsFactors = FALSE
      )
      nr <- nr[!.nts_rank_blank(nr$feature) & !.nts_rank_blank(nr$ref_id), ,
               drop = FALSE]
      if (nrow(nr)) {
        net_only <- nr
      }
    }
  }

  cand <- do.call(rbind, parts)
  cand <- cand[!.nts_rank_blank(cand$feature), , drop = FALSE]
  if (!nrow(cand)) {
    stop("No candidate rows with a non-empty feature id.", call. = FALSE)
  }

  # 同一 (feature, ref_id) 合并
  fnb <- function(x) {
    x <- x[!.nts_rank_blank(x)]
    if (length(x)) x[1] else NA_character_
  }
  fnn <- function(x) {
    x <- x[!is.na(x)]
    if (length(x)) x[1] else NA_real_
  }
  agg <- split(seq_len(nrow(cand)), paste(cand$feature, cand$ref_id, sep = "\r"))
  rows <- lapply(agg, function(ix) {
    d <- cand[ix, , drop = FALSE]
    data.frame(
      feature = d$feature[1], ref_id = d$ref_id[1],
      candidate = fnb(d$candidate),
      src_lib = any(d$src_lib), src_mf = any(d$src_mf), src_net = FALSE,
      ident_set = fnb(d$ident_set), comp_db = fnb(d$comp_db),
      score_lib = fnn(d$score_lib), score_lib_ms2 = fnn(d$score_lib_ms2),
      mf_rank = fnn(d$mf_rank), mf_score = fnn(d$mf_score),
      mf_top1 = any(d$mf_top1), ev_tier = fnn(d$ev_tier),
      network_support = fnb(d$network_support), net_kind = fnb(d$net_kind),
      SMILES = fnb(d$SMILES), db_type_raw = fnb(d$db_type_raw),
      fp_similarity = fnn(d$fp_similarity), score_reranked = fnn(d$score_reranked),
      stringsAsFactors = FALSE
    )
  })
  cand <- do.call(rbind, rows)
  rownames(cand) <- NULL

  # 仅网络传播身份：已存在则补 src_net / net_kind，否则作为新候选补进
  if (!is.null(net_only)) {
    key <- function(f, r) paste(f, r, sep = "\r")
    add <- list()
    for (i in seq_len(nrow(net_only))) {
      j <- match(key(net_only$feature[i], net_only$ref_id[i]),
                 key(cand$feature, cand$ref_id))
      if (!is.na(j)) {
        cand$src_net[j] <- TRUE
        cand$net_kind[j] <- if (is.na(cand$net_kind[j]) || !nzchar(cand$net_kind[j])) {
          net_only$net_kind[i]
        } else {
          paste(unique(c(strsplit(cand$net_kind[j], ";")[[1]],
                         net_only$net_kind[i])), collapse = ";")
        }
      } else {
        add[[length(add) + 1L]] <- data.frame(
          feature = net_only$feature[i], ref_id = net_only$ref_id[i],
          candidate = NA_character_,
          src_lib = FALSE, src_mf = FALSE, src_net = TRUE,
          ident_set = NA_character_, comp_db = .nts_rank_lib_of(net_only$ref_id[i]),
          score_lib = NA_real_, score_lib_ms2 = NA_real_, mf_rank = NA_real_,
          mf_score = NA_real_, mf_top1 = FALSE, ev_tier = NA_real_,
          network_support = paste0("propagated(", net_only$net_kind[i], ")"),
          net_kind = net_only$net_kind[i], SMILES = NA_character_,
          db_type_raw = NA_character_,
          fp_similarity = NA_real_, score_reranked = NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(add)) {
      cand <- rbind(cand, do.call(rbind, add))
    }
  }
  if (!nrow(cand)) {
    stop("No candidates left after merging the supplied sources.", call. = FALSE)
  }

  ## ---------------------------------------------------- 3. 名称 / InChIKey 补全
  namemap <- NULL
  ikmap <- NULL
  if (!is.null(suspect) && is.data.frame(suspect) && nrow(suspect)) {
    id_col <- intersect(c("ID", "Identifier", "ref_id", "Lab_ID"), names(suspect))
    nm_col <- intersect(c("Name", "Compound.name", "compound_name"), names(suspect))
    if (length(id_col) && length(nm_col)) {
      keys <- .nts_rank_chr(suspect[[id_col[1]]])
      namemap <- stats::setNames(.nts_rank_chr(suspect[[nm_col[1]]]), keys)
      ik_col <- intersect(c("InChIKey", "INCHIKEY", "inchikey"), names(suspect))
      if (length(ik_col)) {
        ikmap <- stats::setNames(.nts_rank_chr(suspect[[ik_col[1]]]), keys)
      }
    }
  }
  fill_name <- function(rid, cur) {
    out <- .nts_rank_chr(cur)
    blank <- .nts_rank_blank(out)
    if (!is.null(namemap)) {
      idx <- match(rid, names(namemap))
      hit <- blank & !is.na(idx)
      out[hit] <- namemap[idx[hit]]
    }
    blank <- .nts_rank_blank(out)
    out[blank] <- rid[blank]
    out
  }
  cand$candidate <- fill_name(cand$ref_id, cand$candidate)
  cand$inchikey <- if (!is.null(ikmap)) unname(ikmap[cand$ref_id]) else NA_character_

  ## ----------------------------------------- 4. 基线档 + 证据族 + 半档晋升
  lv <- rep(NA_character_, nrow(cand))
  lv[cand$src_lib] <- .nts_rank_base_level(cand$ident_set[cand$src_lib])
  lv[cand$src_mf & is.na(lv)] <- "L2c"
  only_net <- is.na(lv) & !cand$src_lib & !cand$src_mf
  lv[only_net] <- ifelse(cand$net_kind[only_net] == "HOMO", "L3c", "L3b")
  cand$base_level <- lv

  cand$family_F1 <- cand$src_lib & cand$ident_set %in% c("inhouse", "online")
  cand$family_F2 <- cand$src_mf | (cand$src_lib & cand$ident_set == "insilico")
  # F4 = 来自「别的特征 / 别的方法」的关系证据（ev_tier 1 = 锚定 L1 种子，2 = 传播身份 == 本候选）。
  # ev_tier == 0 只表示「它本身就是 L1 命中」，不是独立佐证，不能计入 F4。
  cand$family_F4 <- (!is.na(cand$ev_tier) & cand$ev_tier %in% c(1, 2)) | cand$src_net
  #「身份级」网络证据只有 ev_tier == 2 或仅网络传播身份；ev_tier == 1 是「与 L1 种子有 ±kCF2 关系」，
  # 是关系而非竞争身份，不参与冲突判定。
  cand$family_F4_ident <- (!is.na(cand$ev_tier) & cand$ev_tier == 2) | cand$src_net
  cand$n_families <- as.integer(cand$family_F1) + as.integer(cand$family_F2) +
    as.integer(cand$family_F4)

  # MetFrag 得分 == 0 表示库内无匹配碎片，区分度为 0，不参与晋升（只当线索）
  mf_valid <- cand$src_mf & cand$mf_top1 & !is.na(cand$mf_score) & cand$mf_score > 0
  tier <- rep(NA_integer_, nrow(cand))
  tier[lv == "L1"] <- 1L
  tier[lv == "L2a"] <- 2L
  tier[lv == "L2b"] <- ifelse(mf_valid[lv == "L2b"], 3L, 4L)
  tier[lv == "L2c"] <- ifelse(cand$family_F4[lv == "L2c"], 5L, 6L)
  tier[lv == "L3a"] <- 7L
  tier[lv == "L3b"] <- ifelse(mf_valid[lv == "L3b"], 6L, 8L)
  tier[lv == "L3c"] <- 9L
  tier[lv %in% c("L3d", "L3e")] <- 10L
  tier[lv == "L4"] <- 11L
  tier[lv == "L5"] <- 12L
  tier[is.na(tier)] <- 12L
  cand$decision_tier <- as.integer(tier)

  cand$decision_label <- paste0(cand$base_level, ifelse(cand$n_families >= 2L, "*", ""))

  ident_cross <- !is.na(cand$ident_set) &
    cand$ident_set %in% c("inhouse", "online") &
    !is.na(cand$comp_db) & cand$comp_db == "insilico"
  cand$flags <- vapply(seq_len(nrow(cand)), function(i) {
    f <- character(0)
    if (cand$src_lib[i] && !is.na(cand$ident_set[i])) {
      f <- c(f, paste0("LIB:", cand$ident_set[i]))
    }
    if (cand$src_mf[i]) {
      f <- c(f, if (isTRUE(cand$mf_top1[i])) "MF_TOP1" else "MF_CAND")
    }
    if (cand$src_net[i]) {
      f <- c(f, paste0("NET:", cand$net_kind[i]))
    }
    if (!is.na(cand$ev_tier[i]) && cand$ev_tier[i] >= 1) {
      f <- c(f, paste0("NET_TIER", as.integer(cand$ev_tier[i])))
    }
    if (ident_cross[i]) {
      f <- c(f, "ID_ACROSS_LIB")
    }
    if (!is.na(cand$fp_similarity[i]) && cand$fp_similarity[i] >= 0.70) {
      f <- c(f, "FP_MATCH")
    }
    paste(unique(f), collapse = ";")
  }, character(1))

  ## --------------------------------------------------- 5. 特征级自证结构信号
  mi <- match(cand$feature, feat$name)
  cand$flag_ladder <- ifelse(
    is.na(mi), FALSE,
    !is.na(feat$pfas_ladder_len[mi]) & feat$pfas_ladder_len[mi] >= 3
  )
  cand$flag_subclass <- ifelse(is.na(mi), NA_character_, feat$pfas_subclass_hint[mi])

  ## -------------------------------------- 6. 每特征排序（tier -> 族数 -> 分数 -> 指纹相似度）
  cand$score_for_rank <- ifelse(
    !is.na(cand$score_reranked), cand$score_reranked,
    ifelse(!is.na(cand$score_lib), cand$score_lib,
           ifelse(!is.na(cand$mf_score), cand$mf_score, NA_real_))
  )
  fp_sim_order <- if ("fp_similarity" %in% names(cand)) {
    v <- cand$fp_similarity
    v[is.na(v)] <- -1.0
    -v
  } else 0
  cand <- cand[order(cand$feature, cand$decision_tier, -cand$n_families,
                     -cand$score_for_rank, fp_sim_order, cand$ref_id, na.last = TRUE), ,
               drop = FALSE]
  fkey <- as.character(cand$feature)
  cand$rank_in_feature <- as.integer(seq_len(nrow(cand)) - match(fkey, fkey) + 1L)
  cand$is_top <- cand$rank_in_feature == 1L
  cand$n_candidates <- as.integer(table(fkey)[fkey])

  ## -------------------------------------- 7. 各证据族首选 + 冲突标记
  feature_list <- unique(cand$feature)
  fam_pick <- function(col) {
    vapply(feature_list, function(f) {
      d <- cand[cand$feature == f & cand[[col]], , drop = FALSE]
      if (!nrow(d)) {
        return(NA_character_)
      }
      d <- d[order(d$decision_tier, -d$n_families, -d$score_for_rank,
                   na.last = TRUE), , drop = FALSE]
      d$ref_id[1]
    }, character(1))
  }
  f1p <- fam_pick("family_F1")
  f2p <- fam_pick("family_F2")
  f4p <- fam_pick("family_F4")

  # 库首选分两个 scope：measured（只用 inhouse/online 实测参考谱）用于与 MetFrag 对照，
  # any（含 insilico）仅展示 —— 拿预测库首选去比 MetFrag 属同族内部噪声。
  best_by_score <- function(f, mask) {
    d <- cand[cand$feature == f & mask, , drop = FALSE]
    if (!nrow(d)) {
      return(NA_character_)
    }
    d <- d[order(-d$score_lib, d$decision_tier, na.last = TRUE), , drop = FALSE]
    d$ref_id[1]
  }
  lib_best_meas <- vapply(
    feature_list, best_by_score, character(1),
    mask = cand$src_lib & cand$ident_set %in% c("inhouse", "online")
  )
  lib_best_any <- vapply(
    feature_list, best_by_score, character(1), mask = cand$src_lib
  )
  mf_top <- vapply(feature_list, function(f) {
    d <- cand[cand$feature == f & cand$src_mf & cand$mf_top1, , drop = FALSE]
    if (!nrow(d)) {
      return(NA_character_)
    }
    d <- d[order(d$mf_rank, na.last = TRUE), , drop = FALSE]
    d$ref_id[1]
  }, character(1))
  lib_best_ref <- ifelse(is.na(lib_best_meas), lib_best_any, lib_best_meas)
  lib_best_scope <- ifelse(
    !is.na(lib_best_meas), "measured",
    ifelse(!is.na(lib_best_any), "insilico", "none")
  )
  mf_vs_lib <- ifelse(
    is.na(lib_best_ref) & is.na(mf_top), "na",
    ifelse(is.na(lib_best_ref), "mf_only",
           ifelse(is.na(mf_top), "lib_only",
                  ifelse(lib_best_ref == mf_top, "agree", "differs")))
  )
  net_ident <- vapply(feature_list, function(f) {
    d <- cand[cand$feature == f & cand$family_F4_ident, , drop = FALSE]
    if (!nrow(d)) {
      return(NA_character_)
    }
    d <- d[order(d$decision_tier, na.last = TRUE), , drop = FALSE]
    d$ref_id[1]
  }, character(1))
  famtbl <- data.frame(
    name = feature_list, fam1_best = f1p, fam2_best = f2p, fam4_best = f4p,
    lib_best_ref = lib_best_ref, lib_best_any = lib_best_any,
    lib_best_scope = lib_best_scope, mf_top1_ref = mf_top, mf_vs_lib = mf_vs_lib,
    net_ident_ref = net_ident, stringsAsFactors = FALSE
  )

  ## ------------------------------------------------------------ 8. 候选排序表
  out_cols <- c(
    "feature", "rank_in_feature", "is_top", "ref_id", "candidate", "inchikey",
    "ident_set", "comp_db", "base_level", "decision_tier", "decision_label",
    "n_families", "family_F1", "family_F2", "family_F4",
    "score_lib", "score_lib_ms2", "mf_rank", "mf_score", "mf_top1",
    "fp_similarity", "score_reranked",
    "ev_tier", "net_kind", "network_support", "flags",
    "flag_ladder", "flag_subclass", "n_candidates", "SMILES"
  )
  cr <- cand[, out_cols, drop = FALSE]

  ## ------------------------------------------------------------ 9. 每特征最优
  top <- cand[cand$is_top, , drop = FALSE]
  ti <- match(top$feature, feature_list)
  top$fam1_best <- famtbl$fam1_best[ti]
  top$fam2_best <- famtbl$fam2_best[ti]
  top$fam4_best <- famtbl$fam4_best[ti]
  top$lib_best_ref <- famtbl$lib_best_ref[ti]
  top$lib_best_any <- famtbl$lib_best_any[ti]
  top$lib_best_scope <- famtbl$lib_best_scope[ti]
  top$mf_top1_ref <- famtbl$mf_top1_ref[ti]
  top$mf_vs_lib <- famtbl$mf_vs_lib[ti]
  top$net_ident_ref <- famtbl$net_ident_ref[ti]
  top$tied_alternatives <- vapply(seq_len(nrow(top)), function(i) {
    d <- cand[
      cand$feature == top$feature[i] &
        cand$decision_tier == top$decision_tier[i] &
        cand$ref_id != top$ref_id[i], , drop = FALSE
    ]
    if (!nrow(d)) NA_character_ else paste(unique(d$candidate), collapse = " | ")
  }, character(1))
  top$review_reason <- vapply(seq_len(nrow(top)), function(i) {
    r <- character(0)
    if (identical(top$mf_vs_lib[i], "differs")) {
      if (identical(top$lib_best_scope[i], "measured")) {
        r <- c(r, "mf_vs_lib_differs")
      } else if (top$decision_tier[i] >= 3L) {
        r <- c(r, "mf_vs_insilico_differs")
      }
    }
    if (!is.na(top$tied_alternatives[i])) {
      r <- c(r, "tied_alternatives")
    }
    if (top$decision_tier[i] >= 3L && !is.na(top$net_ident_ref[i]) &&
        top$net_ident_ref[i] != top$ref_id[i]) {
      r <- c(r, "net_ident_differs")
    }
    paste(r, collapse = ";")
  }, character(1))
  top$conflict_status <- ifelse(
    !nzchar(top$review_reason), "none",
    ifelse(top$decision_tier <= 2L, "resolved_by_lib", "open")
  )
  fi <- match(top$feature, feat$name)
  top$mz <- feat$mz[fi]
  top$rt <- feat$rt[fi]
  top$report_confidence <- feat$confidence_level[fi]

  fa_cols <- c(
    "feature", "mz", "rt", "candidate", "ref_id", "inchikey", "base_level",
    "decision_tier", "decision_label", "n_families", "flags",
    "fp_similarity", "score_reranked",
    "tied_alternatives", "conflict_status", "review_reason",
    "fam1_best", "fam2_best", "fam4_best",
    "lib_best_ref", "lib_best_any", "lib_best_scope",
    "mf_top1_ref", "mf_vs_lib", "net_ident_ref",
    "score_lib", "mf_rank", "mf_score", "mf_top1",
    "ev_tier", "network_support", "flag_ladder", "flag_subclass",
    "n_candidates", "report_confidence"
  )
  fa <- top[, fa_cols, drop = FALSE]
  fa <- fa[order(fa$decision_tier, -fa$n_families, fa$feature), , drop = FALSE]
  names(fa)[names(fa) == "feature"] <- "name"

  ## -------------------------------------------------------------- 10. 冲突表
  cf <- fa[nzchar(fa$review_reason), , drop = FALSE]
  cf <- cf[order(cf$decision_tier, cf$name), , drop = FALSE]

  ## ---------------------------------------------------------------- 11. 汇总
  sm <- as.data.frame(table(fa$decision_label), stringsAsFactors = FALSE)
  names(sm) <- c("final_label", "n_features")
  sm$ord <- .nts_rank_ladder()[sub("\\*$", "", sm$final_label)]
  sm <- sm[order(sm$ord), c("final_label", "n_features")]
  rownames(sm) <- NULL

  say(sprintf("[rank] candidates %d / features %d", nrow(cr), length(feature_list)))
  say(sprintf(
    "[rank] best per feature %d | tied %d | conflict none/resolved/open = %d/%d/%d | review list %d",
    nrow(fa), sum(!is.na(fa$tied_alternatives)),
    sum(fa$conflict_status == "none"), sum(fa$conflict_status == "resolved_by_lib"),
    sum(fa$conflict_status == "open"), nrow(cf)
  ))
  if (verbose) {
    message("[rank] by final label:\n",
            paste(utils::capture.output(print(sm, row.names = FALSE)), collapse = "\n"))
    message(sprintf("[rank] MetFrag-only (no library hit) features: %d",
                    sum(fa$base_level == "L2c")))
  }

  list(candidates = cr, final = fa, conflicts = cf, summary = sm)
}


# ---------------------------------------------------------------------------
# 读目录包装
# ---------------------------------------------------------------------------

#' Merge and rank annotations straight off a pipeline run directory
#'
#' Convenience wrapper around [rank_annotation_candidates]: reads the standard
#' products of a pipeline run and (by default) writes the four merge tables into
#' `out_dir`. Nothing upstream is modified — this is a read-only derived step
#' (it never touches the project, checkpoints, or `@network_results`).
#'
#' Inputs are looked up as:
#' \describe{
#'   \item{candidates}{`network/candidate_ranking.csv`; falls back to
#'     `candidates_all_sources.csv` filtered to `source == "library"` when the
#'     network step did not run.}
#'   \item{metfrag_hits}{`metfrag/metfrag_candidates.csv` (optional).}
#'   \item{network_layers}{`multilib/<set>/ms2_annotations_with_qc.csv`
#'     (`net_df_ref_id` / `net_nl_ref_id` / `net_md_ref_id`) plus
#'     `ms1_annotations_with_confidence.csv` (`homologue_ref_id`, tagged
#'     `"HOMO"`). The `<set>` sub-directory defaults to the one containing
#'     `combined`, else the first that has MS2 annotations.}
#'   \item{feature_meta}{`annotation_master.csv` (optional).}
#'   \item{suspect}{`suspect_csv` when given; otherwise a
#'     `NTS_PFAS_suspect_3lib_dedup.csv` next to `run_dir` or the packaged copy
#'     (`bundled_nts_file()`), plus `multilib/per_library_hits_long.csv` as a
#'     name fallback.}
#' }
#'
#' @param run_dir Directory holding the pipeline products (read-only).
#' @param out_dir Where to write the four CSVs. Default
#'   `<run_dir>/merged`. Created when missing.
#' @param suspect_csv Optional id -> name / InChIKey table (CSV) used to fill
#'   blank candidate names and the `inchikey` column.
#' @param multilib_set Optional name of the `multilib/<set>` sub-directory to
#'   read network identities from (default: the `combined` set, else the first
#'   one that has MS2 annotations).
#' @param write_files If `TRUE` (default), write
#'   `candidates_ranked.csv` / `final_annotation.csv` / `merge_conflicts.csv` /
#'   `merge_summary.csv` into `out_dir`.
#' @param verbose If `TRUE`, report progress.
#'
#' @return Invisibly, the same named list of four `data.frame`s as
#'   [rank_annotation_candidates].
#'
#' @seealso [rank_annotation_candidates], [merge_annotation_layers]
#'
#' @examples
#' \dontrun{
#' merge_annotation_sources("path/to/run_dir")          # -> run_dir/merged/*.csv
#' out <- merge_annotation_sources("path/to/run_dir", write_files = FALSE)
#' head(out$final)
#' }
#'
#' @export
merge_annotation_sources <- function(run_dir,
                                     out_dir = file.path(run_dir, "merged"),
                                     suspect_csv = NULL,
                                     multilib_set = NULL,
                                     predicted_fp = NULL,
                                     weight_fp = 0.25,
                                     write_files = TRUE,
                                     verbose = TRUE) {
  if (!is.character(run_dir) || length(run_dir) != 1L || is.na(run_dir) ||
      !dir.exists(run_dir)) {
    stop("`run_dir` must be the path of an existing directory.", call. = FALSE)
  }
  say <- function(...) if (isTRUE(verbose)) message(...)

  ## --------------------------------------------------------- 1. 库命中候选
  cand_path <- file.path(run_dir, "network", "candidate_ranking.csv")
  candidates <- .nts_rank_read_csv(cand_path)
  if (is.null(candidates)) {
    ca <- .nts_rank_read_csv(file.path(run_dir, "candidates_all_sources.csv"))
    if (!is.null(ca) && nrow(ca)) {
      ca <- ca[.nts_rank_chr(.nts_rank_col(ca, "source")) == "library", ,
               drop = FALSE]
      ca_col <- .nts_rank_feature_col(ca)
      if (nrow(ca) && !is.null(ca_col) && "ref_id" %in% names(ca)) {
        say("[merge] network/candidate_ranking.csv not found; ",
            "falling back to candidates_all_sources.csv (library rows only)")
        candidates <- data.frame(
          feature = .nts_rank_chr(ca[[ca_col]]),
          ref_id = .nts_rank_chr(ca$ref_id),
          candidate = .nts_rank_chr(.nts_rank_col(ca, "candidate")),
          db_type = .nts_rank_chr(.nts_rank_col(ca, "db_type")),
          score_total = .nts_rank_num(.nts_rank_col(ca, "score_total")),
          score_ms2 = .nts_rank_num(.nts_rank_col(ca, "score_ms2")),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (is.null(candidates) || !nrow(candidates)) {
    stop("No usable candidate table in `run_dir` (looked for ",
         "network/candidate_ranking.csv and candidates_all_sources.csv).",
         call. = FALSE)
  }

  ## ------------------------------------------------------------- 2. MetFrag
  metfrag_hits <- .nts_rank_read_csv(
    file.path(run_dir, "metfrag", "metfrag_candidates.csv")
  )

  ## --------------------------------------------------------- 3. 特征级信息
  feature_meta <- .nts_rank_read_csv(file.path(run_dir, "annotation_master.csv"))

  ## ------------------------------------------------- 4. 网络传播身份（长表）
  mdir <- file.path(run_dir, "multilib")
  sets <- if (dir.exists(mdir)) {
    list.dirs(mdir, recursive = FALSE)
  } else {
    character(0)
  }
  set_dir <- NULL
  if (!is.null(multilib_set)) {
    cand_dir <- file.path(mdir, multilib_set)
    if (dir.exists(cand_dir)) {
      set_dir <- cand_dir
    }
  }
  if (is.null(set_dir) && length(sets)) {
    has_m2 <- sets[file.exists(file.path(sets, "ms2_annotations_with_qc.csv"))]
    pick <- has_m2[grepl("combined", basename(has_m2))]
    set_dir <- if (length(pick)) pick[1] else if (length(has_m2)) has_m2[1] else NULL
  }
  layers <- list()
  if (!is.null(set_dir)) {
    ms2ann <- .nts_rank_read_csv(file.path(set_dir, "ms2_annotations_with_qc.csv"))
    if (!is.null(ms2ann) && nrow(ms2ann)) {
      nm2 <- .nts_rank_chr(ms2ann[[.nts_rank_feature_col(ms2ann)]])
      for (k in list(c("DFI", "net_df_ref_id"), c("NL", "net_nl_ref_id"),
                     c("MD", "net_md_ref_id"))) {
        if (!k[2] %in% names(ms2ann)) {
          next
        }
        v <- .nts_rank_chr(ms2ann[[k[2]]])
        keep <- !.nts_rank_blank(v)
        if (any(keep)) {
          layers[[length(layers) + 1L]] <- data.frame(
            feature = nm2[keep], ref_id = v[keep], net_kind = k[1],
            stringsAsFactors = FALSE
          )
        }
      }
    }
    ms1ann <- .nts_rank_read_csv(
      file.path(set_dir, "ms1_annotations_with_confidence.csv")
    )
    if (!is.null(ms1ann) && nrow(ms1ann) && "homologue_ref_id" %in% names(ms1ann)) {
      v <- .nts_rank_chr(ms1ann$homologue_ref_id)
      keep <- !.nts_rank_blank(v)
      if (any(keep)) {
        layers[[length(layers) + 1L]] <- data.frame(
          feature = .nts_rank_chr(ms1ann[[.nts_rank_feature_col(ms1ann)]])[keep],
          ref_id = v[keep], net_kind = "HOMO", stringsAsFactors = FALSE
        )
      }
    }
  }
  network_layers <- if (length(layers)) do.call(rbind, layers) else NULL

  ## ------------------------------------------------------- 5. 名称 / InChIKey
  suspect <- NULL
  if (is.null(suspect_csv)) {
    cand_paths <- c(
      file.path(run_dir, "NTS_PFAS_suspect_3lib_dedup.csv"),
      tryCatch(bundled_nts_file("extdata/tables/NTS_PFAS_suspect_3lib_dedup.csv",
                                mustWork = FALSE),
               error = function(e) NA_character_)
    )
    hit <- cand_paths[!is.na(cand_paths) & file.exists(cand_paths)]
    if (length(hit)) {
      suspect_csv <- hit[1]
    }
  }
  sus_parts <- list()
  sus_raw <- .nts_rank_read_csv(suspect_csv)
  if (!is.null(sus_raw) && nrow(sus_raw)) {
    idc <- intersect(c("ID", "Identifier", "ref_id", "Lab_ID"), names(sus_raw))
    nmc <- intersect(c("Name", "Compound.name", "compound_name"), names(sus_raw))
    if (length(idc) && length(nmc)) {
      d <- data.frame(
        ID = .nts_rank_chr(sus_raw[[idc[1]]]),
        Name = .nts_rank_chr(sus_raw[[nmc[1]]]),
        stringsAsFactors = FALSE
      )
      ikc <- intersect(c("InChIKey", "INCHIKEY", "inchikey"), names(sus_raw))
      if (length(ikc)) {
        d$InChIKey <- .nts_rank_chr(sus_raw[[ikc[1]]])
      }
      sus_parts[[length(sus_parts) + 1L]] <- d
    }
  }
  # 库命中长表作为名称兜底（放在后面：同键时上面的表优先）
  pl <- .nts_rank_read_csv(file.path(run_dir, "multilib",
                                     "per_library_hits_long.csv"))
  if (!is.null(pl) && all(c("ref_id", "compound") %in% names(pl))) {
    sus_parts[[length(sus_parts) + 1L]] <- data.frame(
      ID = .nts_rank_chr(pl$ref_id),
      Name = .nts_rank_chr(pl$compound),
      stringsAsFactors = FALSE
    )
  }
  suspect <- if (length(sus_parts)) {
    do.call(rbind, lapply(sus_parts, function(d) {
      if (!"InChIKey" %in% names(d)) {
        d$InChIKey <- NA_character_
      }
      d[, c("ID", "Name", "InChIKey"), drop = FALSE]
    }))
  } else NULL

  ## --------------------------------------------------------------- 6. 计算
  res <- rank_annotation_candidates(
    candidates = candidates,
    metfrag_hits = metfrag_hits,
    network_layers = network_layers,
    feature_meta = feature_meta,
    suspect = suspect,
    predicted_fp = predicted_fp,
    weight_fp = weight_fp,
    verbose = verbose
  )

  ## --------------------------------------------------------------- 7. 落盘
  if (isTRUE(write_files)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    files <- c(
      candidates_ranked.csv = "candidates",
      final_annotation.csv = "final",
      merge_conflicts.csv = "conflicts",
      merge_summary.csv = "summary"
    )
    for (fn in names(files)) {
      utils::write.csv(res[[files[[fn]]]], file.path(out_dir, fn),
                       row.names = FALSE, na = "")
    }
    say("[merge] wrote ", length(files), " tables to ", out_dir)
  }
  invisible(res)
}
