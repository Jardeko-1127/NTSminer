## Molecular network data construction (nodes / edges)
##
## Exported: [get_network_data], [generate_network_data], [extract_subnetwork],
##   [filter_network].
## MS2 list flattening / scan consensus: `05_ms2_extract.R`.

# ============================================================================
# Molecular network data assembly (nodes / edges for visNetwork; data.table)
# ============================================================================

#' @noRd
.nts_expand_kmd_series_long <- function(ann_dt) {
  if (!nrow(ann_dt) || !"kmd_series_id" %in% names(ann_dt)) {
    return(data.table::data.table(
      name = character(0),
      mz = numeric(0),
      kid = character(0)
    ))
  }
  x <- ann_dt[!is.na(kmd_series_id) & nzchar(as.character(kmd_series_id))]
  if (!nrow(x)) {
    return(data.table::data.table(
      name = character(0),
      mz = numeric(0),
      kid = character(0)
    ))
  }
  bycols <- setdiff(names(x), "kmd_series_id")
  x[, {
    parts <- trimws(strsplit(as.character(kmd_series_id), ";\\s*", perl = TRUE)[[1L]])
    parts <- parts[nzchar(parts) & !(toupper(parts) %in% c("NA", "NAN"))]
    if (!length(parts)) {
      .(kid = character(0))
    } else {
      .(kid = parts)
    }
  }, by = bycols]
}

#' @noRd
.nts_ms2_edge_empty_dt <- function() {
  data.table::data.table(
    from = character(0),
    to = character(0),
    match_type = character(0),
    Interaction = character(0),
    label = character(0),
    shared_fragments_count = integer(0),
    shared_fragments_mz = character(0),
    shared_nls_count = integer(0),
    shared_nls_mz = character(0),
    shared_mds_count = integer(0),
    shared_mds_mz = character(0),
    edge_score = numeric(0),
    semantic_similarity = numeric(0)
  )
}


#' Normalised MS2 edge weight: shared items over the smaller feature's peak set.
#'
#' `n_hit / min(n_a, n_b)` in `(0, 1]`; larger = stronger shared evidence. Counts
#' come from the post-filter consensus peaks used to build the edge, so a
#' 2-shared-of-2 spectrum scores higher than 2-shared-of-40.
#' @noRd
.nts_ms2_edge_score <- function(n_hit, n_a, n_b) {
  denom <- pmax(1, pmin(as.numeric(n_a), as.numeric(n_b)))
  s <- as.numeric(n_hit) / denom
  s[!is.finite(s)] <- NA_real_
  pmin(1, pmax(0, s))
}


#' @noRd
.nts_network_ms2_fragment_edges <- function(clean_dt, mz_tol_ppm, min_k) {
  if (!nrow(clean_dt)) {
    return(.nts_ms2_edge_empty_dt())
  }
  ppm <- as.numeric(mz_tol_ppm)[1L] * 1e-6
  u <- unique(clean_dt[is.finite(mz), .(name, mz)], by = c("name", "mz"))
  if (!nrow(u)) {
    return(.nts_ms2_edge_empty_dt())
  }
  fr_cnt <- u[, .(n_peaks = .N), by = .(name)]
  u[, `:=`(mz_lo = mz * (1 - ppm), mz_hi = mz * (1 + ppm))]
  L <- u[, .(name_a = name, mz_a = mz, mz_lo, mz_hi)]
  R <- u[, .(name_b = name, mz_b = mz, mz_b_val = mz)]
  data.table::setkey(R, mz_b)
  hits_fr <- R[
    L,
    on = .(mz_b >= mz_lo, mz_b <= mz_hi),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]
  hits_fr <- hits_fr[name_a < name_b]
  if (!nrow(hits_fr)) {
    return(.nts_ms2_edge_empty_dt())
  }
  if (min_k > 1L) {
    cand_pairs <- hits_fr[, .N, by = .(name_a, name_b)][N >= min_k, .(name_a, name_b)]
    if (!nrow(cand_pairs)) {
      return(.nts_ms2_edge_empty_dt())
    }
    hits_fr <- hits_fr[cand_pairs, on = .(name_a, name_b), nomatch = 0L]
  }
  # 对称共享计数：只数「左向」匹配键会因多个近端峰落进同一对端峰而重复计数
  # 进而越 min(n_a, n_b)（并放宽 min_k 闸门、将 edge_score 虚高钳到 1）。取两侧
  # 去重匹配键数的较小者；键列表同取较小侧，使长度与 n_hit 一致。
  frag_n <- hits_fr[
    ,
    {
      ua <- sort(unique(mz_a[is.finite(mz_a)]))
      ub <- sort(unique(mz_b_val[is.finite(mz_b_val)]))
      keys <- if (length(ua) <= length(ub)) ua else ub
      .(
        n_hit = min(length(ua), length(ub)),
        shared_fragments_mz = if (length(keys)) {
          paste(format(keys, digits = 8L, trim = TRUE), collapse = "; ")
        } else {
          ""
        }
      )
    },
    by = .(name_a, name_b)
  ]
  frag_n <- frag_n[n_hit >= min_k]
  if (!nrow(frag_n)) {
    return(.nts_ms2_edge_empty_dt())
  }
  frag_n[, n_a := fr_cnt$n_peaks[match(name_a, fr_cnt$name)]]
  frag_n[, n_b := fr_cnt$n_peaks[match(name_b, fr_cnt$name)]]
  frag_n[, edge_score := .nts_ms2_edge_score(n_hit, n_a, n_b)]
  frag_n[
    ,
    .(
      from = name_a,
      to = name_b,
      Interaction = "MS2_Shared_Fragment",
      label = paste0(n_hit, " shared m/z"),
      shared_fragments_count = as.integer(n_hit),
      shared_fragments_mz = shared_fragments_mz,
      shared_nls_count = NA_integer_,
      shared_nls_mz = NA_character_,
      shared_mds_count = NA_integer_,
      shared_mds_mz = NA_character_,
      edge_score = edge_score
    )
  ]
}


#' @noRd
.nts_network_ms2_nl_edges <- function(clean_dt, mz_tol_ppm, min_k) {
  if (!nrow(clean_dt) || !"neutral_loss" %in% names(clean_dt)) {
    return(.nts_ms2_edge_empty_dt())
  }
  ppm <- as.numeric(mz_tol_ppm)[1L] * 1e-6
  loss_u <- unique(
    clean_dt[is.finite(neutral_loss), .(name, loss = neutral_loss)],
    by = c("name", "loss")
  )
  if (!nrow(loss_u)) {
    return(.nts_ms2_edge_empty_dt())
  }
  nl_cnt <- loss_u[, .(n_peaks = .N), by = .(name)]
  loss_u[, `:=`(mz_lo = loss * (1 - ppm), mz_hi = loss * (1 + ppm))]
  Ll <- loss_u[, .(name_a = name, loss_a = loss, mz_lo, mz_hi)]
  Rl <- loss_u[, .(name_b = name, loss_b = loss, loss_b_val = loss)]
  data.table::setkey(Rl, loss_b)
  hits_ls <- Rl[
    Ll,
    on = .(loss_b >= mz_lo, loss_b <= mz_hi),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]
  hits_ls <- hits_ls[name_a < name_b]
  if (!nrow(hits_ls)) {
    return(.nts_ms2_edge_empty_dt())
  }
  if (min_k > 1L) {
    cand_pairs <- hits_ls[, .N, by = .(name_a, name_b)][N >= min_k, .(name_a, name_b)]
    if (!nrow(cand_pairs)) {
      return(.nts_ms2_edge_empty_dt())
    }
    hits_ls <- hits_ls[cand_pairs, on = .(name_a, name_b), nomatch = 0L]
  }
  # 对称共享计数（同 DFI）：min(两侧去重匹配键数)，键列表取较小侧。
  loss_n <- hits_ls[
    ,
    {
      ua <- sort(unique(loss_a[is.finite(loss_a)]))
      ub <- sort(unique(loss_b_val[is.finite(loss_b_val)]))
      keys <- if (length(ua) <= length(ub)) ua else ub
      .(
        n_hit = min(length(ua), length(ub)),
        shared_nls_mz = if (length(keys)) {
          paste(format(keys, digits = 8L, trim = TRUE), collapse = "; ")
        } else {
          ""
        }
      )
    },
    by = .(name_a, name_b)
  ]
  loss_n <- loss_n[n_hit >= min_k]
  if (!nrow(loss_n)) {
    return(.nts_ms2_edge_empty_dt())
  }
  loss_n[, n_a := nl_cnt$n_peaks[match(name_a, nl_cnt$name)]]
  loss_n[, n_b := nl_cnt$n_peaks[match(name_b, nl_cnt$name)]]
  loss_n[, edge_score := .nts_ms2_edge_score(n_hit, n_a, n_b)]
  loss_n[
    ,
    .(
      from = name_a,
      to = name_b,
      Interaction = "MS2_Shared_Neutral_Loss",
      label = paste0(n_hit, " shared losses"),
      shared_fragments_count = NA_integer_,
      shared_fragments_mz = NA_character_,
      shared_nls_count = as.integer(n_hit),
      shared_nls_mz = shared_nls_mz,
      shared_mds_count = NA_integer_,
      shared_mds_mz = NA_character_,
      edge_score = edge_score
    )
  ]
}


#' Mass-difference (MD) MS2 edges: two features share intra-spectrum pairwise
#' |Δmz| values (>= `md_min_da`). NOTE on `edge_score` denominator: unlike DFI /
#' NL, `n_a` / `n_b` here count *intra-spectrum pairwise differences*
#' (k(k-1)/2-scale for k peaks), not peaks, so the same `edge_score` formula
#' yields systematically smaller values for MD than for DFI / NL. Compare
#' `edge_score` only within an edge layer, not across layers.
#' @noRd
.nts_network_ms2_massdiff_edges <- function(clean_dt, mz_tol_ppm, min_k, md_min_da = 1.0) {
  if (!nrow(clean_dt)) {
    return(.nts_ms2_edge_empty_dt())
  }
  md_floor <- suppressWarnings(as.numeric(md_min_da)[1L])
  if (!is.finite(md_floor)) {
    md_floor <- 1.0
  }
  ppm <- as.numeric(mz_tol_ppm)[1L] * 1e-6
  diff_long <- clean_dt[
    ,
    {
      u <- unique(mz[is.finite(mz)])
      u <- sort(u)
      if (length(u) < 2L) {
        .(mass_diff = NA_real_)
      } else {
        m <- outer(u, u, function(a, b) abs(a - b))
        md <- m[lower.tri(m, diag = FALSE)]
        md <- md[md >= md_floor & is.finite(md)]
        if (!length(md)) {
          .(mass_diff = NA_real_)
        } else {
          .(mass_diff = md)
        }
      }
    },
    by = name
  ]
  diff_long <- diff_long[is.finite(mass_diff)]
  if (!nrow(diff_long)) {
    return(.nts_ms2_edge_empty_dt())
  }
  du <- unique(diff_long[, .(name, mass_diff)], by = c("name", "mass_diff"))
  md_cnt <- du[, .(n_peaks = .N), by = .(name)]
  du[, `:=`(md_lo = mass_diff * (1 - ppm), md_hi = mass_diff * (1 + ppm))]
  Ld <- du[, .(name_a = name, md_a = mass_diff, md_lo, md_hi)]
  Rd <- du[, .(name_b = name, md_b = mass_diff, md_b_val = mass_diff)]
  data.table::setkey(Rd, md_b)
  hits_md <- Rd[
    Ld,
    on = .(md_b >= md_lo, md_b <= md_hi),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]
  hits_md <- hits_md[name_a < name_b]
  if (!nrow(hits_md)) {
    return(.nts_ms2_edge_empty_dt())
  }
  if (min_k > 1L) {
    cand_pairs <- hits_md[, .N, by = .(name_a, name_b)][N >= min_k, .(name_a, name_b)]
    if (!nrow(cand_pairs)) {
      return(.nts_ms2_edge_empty_dt())
    }
    hits_md <- hits_md[cand_pairs, on = .(name_a, name_b), nomatch = 0L]
  }
  # 对称共享计数（同 DFI/NL）：min(两侧去重匹配键数)，键列表取较小侧。
  md_n <- hits_md[
    ,
    {
      ua <- sort(unique(md_a[is.finite(md_a)]))
      ub <- sort(unique(md_b_val[is.finite(md_b_val)]))
      keys <- if (length(ua) <= length(ub)) ua else ub
      .(
        n_hit = min(length(ua), length(ub)),
        shared_mds_mz = if (length(keys)) {
          paste(format(keys, digits = 8L, trim = TRUE), collapse = "; ")
        } else {
          ""
        }
      )
    },
    by = .(name_a, name_b)
  ]
  md_n <- md_n[n_hit >= min_k]
  if (!nrow(md_n)) {
    return(.nts_ms2_edge_empty_dt())
  }
  md_n[, n_a := md_cnt$n_peaks[match(name_a, md_cnt$name)]]
  md_n[, n_b := md_cnt$n_peaks[match(name_b, md_cnt$name)]]
  md_n[, edge_score := .nts_ms2_edge_score(n_hit, n_a, n_b)]
  md_n[
    ,
    .(
      from = name_a,
      to = name_b,
      Interaction = "MS2_Shared_Mass_Diff",
      label = paste0(n_hit, " shared mass diffs"),
      shared_fragments_count = NA_integer_,
      shared_fragments_mz = NA_character_,
      shared_nls_count = NA_integer_,
      shared_nls_mz = NA_character_,
      shared_mds_count = as.integer(n_hit),
      shared_mds_mz = shared_mds_mz,
      edge_score = edge_score
    )
  ]
}


#' Aligned ndotproduct after shifting the partner spectrum by `shift` Da.
#' @noRd
.nts_shift_cosine <- function(A, B, shift, mz_tol = 0.01, int_thresh = 0.025) {
  if (is.finite(shift) && shift != 0) {
    B <- cbind(mz = B[, 1L] + shift, intensity = B[, 2L])
  }
  .calculate_clean_ndp(A, B, mz_tol = mz_tol, int_thresh = int_thresh)$score
}


#' Cosine / modified-cosine MS2 similarity edges (path B folded into path A).
#'
#' Plain cosine (`modified = FALSE`) is ndotproduct over the aligned consensus
#' peak sets. Modified cosine (`modified = TRUE`) additionally allows the whole
#' partner spectrum to be shifted by the precursor m/z difference (bounded by
#' `shift_max_da`), capturing "A's fragment = B's fragment + delta" analog
#' relations DFI / NL / MD cannot express. Pairs are scored only when their
#' precursor m/z differ by at most `precursor_tol` (plain) or `shift_max_da`
#' (modified). `edge_score` is the (max) similarity; `shared_*` columns are `NA`.
#' @noRd
.nts_network_ms2_cosine_edges <- function(clean_dt,
                                          prec_dt,
                                          min_score = 0.6,
                                          precursor_tol = 1.0,
                                          modified = FALSE,
                                          shift_max_da = 500,
                                          mz_tol = 0.01,
                                          int_thresh = 0.025,
                                          min_peaks = 2L) {
  if (!nrow(clean_dt) || !nrow(prec_dt)) {
    return(.nts_ms2_edge_empty_dt())
  }
  sp <- clean_dt[is.finite(mz) & is.finite(intensity) & intensity > 0]
  if (!nrow(sp)) {
    return(.nts_ms2_edge_empty_dt())
  }
  ids <- unique(as.character(sp$name))
  specs <- lapply(ids, function(nm) {
    d <- sp[name == nm]
    m <- cbind(mz = as.numeric(d$mz), intensity = as.numeric(d$intensity))
    m <- m[order(m[, 1L]), , drop = FALSE]
    rownames(m) <- NULL
    m
  })
  names(specs) <- ids
  npeak <- vapply(specs, nrow, integer(1L))
  mp <- suppressWarnings(as.integer(min_peaks)[1L])
  if (!is.finite(mp) || mp < 1L) mp <- 1L
  keep <- npeak >= mp
  ids <- ids[keep]
  specs <- specs[keep]
  if (length(ids) < 2L) {
    return(.nts_ms2_edge_empty_dt())
  }
  pv <- suppressWarnings(as.numeric(prec_dt$prec_mz[match(ids, as.character(prec_dt$name))]))
  ok <- is.finite(pv)
  ids <- ids[ok]
  specs <- specs[ok]
  pv <- pv[ok]
  if (length(ids) < 2L) {
    return(.nts_ms2_edge_empty_dt())
  }
  prec_tol <- suppressWarnings(as.numeric(precursor_tol)[1L])
  if (!is.finite(prec_tol) || prec_tol <= 0) prec_tol <- 1.0
  win <- prec_tol
  sh_max <- suppressWarnings(as.numeric(shift_max_da)[1L])
  if (!is.finite(sh_max) || sh_max < 0) sh_max <- 500
  if (isTRUE(modified)) win <- sh_max
  minsc <- suppressWarnings(as.numeric(min_score)[1L])
  if (!is.finite(minsc)) minsc <- 0.6
  o <- order(pv)
  ids <- ids[o]
  specs <- specs[o]
  pv <- pv[o]
  N <- length(ids)
  fr <- character(0)
  to <- character(0)
  sc <- numeric(0)
  for (a in seq_len(N - 1L)) {
    b <- a + 1L
    while (b <= N && (pv[b] - pv[a]) <= win) {
      d <- pv[a] - pv[b]
      shifts <- 0
      if (isTRUE(modified) && is.finite(d) && d != 0 && abs(d) <= sh_max) {
        shifts <- c(0, d, -d)
      }
      best <- -Inf
      for (s in shifts) {
        v <- tryCatch(
          .nts_shift_cosine(specs[[a]], specs[[b]], s, mz_tol = mz_tol, int_thresh = int_thresh),
          error = function(e) NA_real_
        )
        if (is.finite(v) && v > best) best <- v
      }
      if (is.finite(best) && best >= minsc) {
        fr <- c(fr, ids[a])
        to <- c(to, ids[b])
        sc <- c(sc, best)
      }
      b <- b + 1L
    }
  }
  if (!length(fr)) {
    return(.nts_ms2_edge_empty_dt())
  }
  lo <- pmin(fr, to)
  hi <- pmax(fr, to)
  data.table::data.table(
    from = lo,
    to = hi,
    Interaction = if (isTRUE(modified)) "MS2_Modified_Cosine" else "MS2_Cosine",
    label = paste0(
      sprintf("%.3f", sc),
      if (isTRUE(modified)) " modified cosine" else " cosine"
    ),
    shared_fragments_count = NA_integer_,
    shared_fragments_mz = NA_character_,
    shared_nls_count = NA_integer_,
    shared_nls_mz = NA_character_,
    shared_mds_count = NA_integer_,
    shared_mds_mz = NA_character_,
    edge_score = sc
  )
}


#' Compute MS2 spectral semantic similarity edges using embeddings (256-d or DreaMS/Spec2Vec)
#' @noRd
.nts_network_ms2_semantic_edges <- function(project,
                                            valid_ids,
                                            min_score = 0.70,
                                            top_k = 3L,
                                            embeddings = NULL,
                                            python_path = NULL) {
  emb_mat <- NULL
  if (!is.null(embeddings)) {
    emb_mat <- .nts_as_embedding_matrix(embeddings)
    if (is.null(emb_mat) && is.matrix(embeddings)) {
      emb_mat <- embeddings
    }
  }

  if (is.null(emb_mat) && inherits(project, "NTSProject")) {
    proxy <- get_physchem_proxy(project)
    if (is.data.frame(proxy) && nrow(proxy)) {
      emb_mat <- .nts_as_embedding_matrix(proxy)
    }
    if (is.null(emb_mat)) {
      clean_l <- methods::slot(project, "clean_ms2")
      if (is.list(clean_l) && length(clean_l)) {
        emb_res <- tryCatch(
          extract_ms2_embeddings(project, python_path = python_path, write_proxy = FALSE, verbose = FALSE),
          error = function(e) NULL
        )
        if (!is.null(emb_res)) {
          emb_mat <- .nts_as_embedding_matrix(emb_res)
          if (is.null(emb_mat) && is.matrix(emb_res)) {
            emb_mat <- emb_res
          }
        }
      }
    }
  }

  if (is.null(emb_mat) || !is.matrix(emb_mat) || nrow(emb_mat) < 2L) {
    return(.nts_ms2_edge_empty_dt())
  }

  fids <- rownames(emb_mat)
  if (is.null(fids)) {
    if (nrow(emb_mat) == length(valid_ids)) {
      fids <- as.character(valid_ids)
    } else {
      fids <- paste0("FT_", seq_len(nrow(emb_mat)))
    }
    rownames(emb_mat) <- fids
  }
  keep_mask <- fids %in% as.character(valid_ids)
  emb_mat <- emb_mat[keep_mask, , drop = FALSE]
  if (nrow(emb_mat) < 2L) {
    return(.nts_ms2_edge_empty_dt())
  }

  # Drop non-finite rows (NA/NaN/Inf) defensively
  finite_rows <- rowSums(is.finite(emb_mat)) == ncol(emb_mat)
  emb_mat <- emb_mat[finite_rows, , drop = FALSE]
  if (nrow(emb_mat) < 2L) {
    return(.nts_ms2_edge_empty_dt())
  }

  # L2-normalize rows
  row_norms <- sqrt(rowSums(emb_mat^2))
  row_norms[row_norms < 1e-12] <- 1.0
  emb_mat <- emb_mat / row_norms

  # Dot product equals cosine similarity for L2-normalized vectors
  sim_mat <- emb_mat %*% t(emb_mat)
  diag(sim_mat) <- 0

  min_s <- suppressWarnings(as.numeric(min_score)[1L])
  if (!is.finite(min_s) || min_s < 0) min_s <- 0.70
  tk <- if (is.finite(top_k)) max(as.integer(top_k), 1L) else ncol(sim_mat)

  n_nodes <- nrow(sim_mat)
  node_names <- rownames(sim_mat)
  edges_list <- vector("list", n_nodes)

  for (i in seq_len(n_nodes)) {
    sims <- sim_mat[i, ]
    v_idx <- which(sims >= min_s)
    if (length(v_idx)) {
      ord <- order(sims[v_idx], decreasing = TRUE)
      sel_j <- v_idx[ord[seq_len(min(length(ord), tk))]]
      edges_list[[i]] <- data.table::data.table(
        from = node_names[i],
        to = node_names[sel_j],
        sim = as.numeric(sims[sel_j])
      )
    }
  }

  all_e <- data.table::rbindlist(edges_list[!vapply(edges_list, is.null, logical(1L))], use.names = TRUE)
  if (!nrow(all_e)) {
    return(.nts_ms2_edge_empty_dt())
  }

  # Canonicalize undirected pairs: from < to
  all_e[, c("u", "v") := .(pmin(from, to), pmax(from, to))]
  all_e <- all_e[u != v]
  if (!nrow(all_e)) {
    return(.nts_ms2_edge_empty_dt())
  }
  # Keep maximum similarity between u and v
  all_e <- all_e[, .(sim = max(sim)), by = .(u, v)]

  data.table::data.table(
    from = all_e$u,
    to = all_e$v,
    match_type = "semantic",
    Interaction = "MS2_Semantic",
    label = sprintf("%.3f semantic", all_e$sim),
    shared_fragments_count = NA_integer_,
    shared_fragments_mz = NA_character_,
    shared_nls_count = NA_integer_,
    shared_nls_mz = NA_character_,
    shared_mds_count = NA_integer_,
    shared_mds_mz = NA_character_,
    edge_score = round(all_e$sim, 4L),
    semantic_similarity = round(all_e$sim, 4L)
  )
}



#' ISF edges from `is_isf` / `parent_mz` on the node table (sourced from
#' `@ms1_annotations` via `.nts_network_nodes_dt`, not `@feature_table`).
#' @noRd
.nts_build_isf_edges_dt <- function(nodes_dt, mz_tol_ppm = 20) {
  empty <- data.table::data.table(
    from = character(0),
    to = character(0),
    edge_type = character(0),
    label = character(0)
  )
  if (!nrow(nodes_dt)) {
    return(empty)
  }
  req <- c("name", "mz", "is_isf", "parent_mz")
  if (!all(req %in% names(nodes_dt))) {
    return(empty)
  }
  ppm <- suppressWarnings(as.numeric(mz_tol_ppm)[1L]) * 1e-6
  if (!is.finite(ppm) || ppm <= 0) {
    ppm <- 20 * 1e-6
  }
  dt <- data.table::copy(nodes_dt)
  dt[, `:=`(name = as.character(name), mz = as.numeric(mz))]
  ch <- dt[is_isf %in% TRUE & is.finite(parent_mz)]
  if (!nrow(ch)) {
    return(empty)
  }
  ch <- ch[, .(to = name, pc = as.numeric(parent_mz))]
  ch <- ch[is.finite(pc)]
  if (!nrow(ch)) {
    return(empty)
  }
  ch[, tol := pmax(abs(pc) * ppm, 1e-6)]
  ch[, `:=`(lo = pc - tol, hi = pc + tol)]
  par <- dt[, .(from = name, mz = mz)]
  hits <- par[
    ch,
    on = .(mz >= lo, mz <= hi),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]
  hits <- hits[from != to]
  if (!nrow(hits)) {
    return(empty)
  }
  hits[, d := abs(mz - pc)]
  hits <- hits[, .SD[which.min(d)], by = .(to)]
  out <- hits[
    ,
    .(
      from = as.character(from),
      to = as.character(to),
      edge_type = "ISF",
      label = "ISF"
    )
  ]
  unique(out, by = c("from", "to"))
}


#' @param ms1_layers Subset of `c("homologue", "transformation", "isf")` to
#'   compute. Empty vector returns an empty edge table without touching KMD / ISF
#'   machinery.
#' @param cf2_ladder_boost When `TRUE` (default), a Homologue edge whose two
#'   endpoints each carry an MS2 CF2 fragment ladder (`pfas_ladder_len >= 3` in
#'   the node table, from [detect_pfas_fragment_series]) is flagged
#'   `ms2_ladder_confirmed = TRUE` and given `edge_score = 1`. This **reinforces**
#'   the MS1 homologue edge with MS2 evidence without inventing new feature-level
#'   edges. No-op when the node table lacks `pfas_ladder_len`.
#' @noRd
.nts_ms1_layer_edges_dt <- function(nodes_dt,
                                    valid_ids,
                                    isf_mz_tol_ppm = 20,
                                    ms1_layers = c("homologue", "transformation", "isf"),
                                    cf2_ladder_boost = TRUE) {
  ms1_empty <- data.table::data.table(
    from = character(0),
    to = character(0),
    edge_type = character(0),
    label = character(0),
    mass_difference = numeric(0),
    n_steps = integer(0),
    transformation = character(0),
    edge_score = numeric(0),
    ms2_ladder_confirmed = logical(0)
  )
  layers <- .nts_normalize_ms1_edge_types(ms1_layers)
  if (!length(layers)) {
    return(ms1_empty)
  }
  ladder_len <- if ("pfas_ladder_len" %in% names(nodes_dt)) {
    v <- suppressWarnings(as.numeric(nodes_dt$pfas_ladder_len))
    names(v) <- as.character(nodes_dt$name)
    v
  } else {
    NULL
  }
  edge_list <- list()
  if ("homologue" %in% layers) {
    kmd_long <- .nts_expand_kmd_series_long(nodes_dt)
    if (nrow(kmd_long) && "mz" %in% names(kmd_long)) {
      data.table::setorder(kmd_long, kid, mz, name)
      kmd_long[, from := data.table::shift(name, n = 1L, type = "lag"), by = kid]
      kmd_long[, mz_from := data.table::shift(mz, n = 1L, type = "lag"), by = kid]
      kmd_long[, mass_difference := as.numeric(mz) - as.numeric(mz_from)]
      ek <- kmd_long[
        !is.na(from) & from != name,
        .(
          from = as.character(from),
          to = as.character(name),
          edge_type = "Homologue",
          mass_difference = as.numeric(mass_difference),
          ref_mz = pmax(abs(as.numeric(mz_from)), abs(as.numeric(mz)), na.rm = TRUE),
          label = ""
        )
      ]
      if (nrow(ek)) {
        # CF2 step count + signed transformation label (`+2CF2` etc.), using the
        # exact CF2 constant and gap tolerance (P1-1).
        ek[, transformation := .nts_homologue_transformation_tag(
          mass_difference,
          ref_mz,
          .nts_pfas_cf2_da,
          20,
          "CF2",
          max_gap = 3L
        )]
        ek[, n_steps := {
          ns <- suppressWarnings(as.integer(round(mass_difference / .nts_pfas_cf2_da)))
          ifelse(transformation == "Unknown", NA_integer_, ns)
        }]
        ek[, label := ifelse(transformation == "Unknown", "", transformation)]
        ek[, ref_mz := NULL]
        # (S3-a) Reinforce these MS1 edges with per-feature MS2 CF2-ladder evidence.
        ek[, ms2_ladder_confirmed := FALSE]
        ek[, edge_score := NA_real_]
        if (isTRUE(cf2_ladder_boost) && !is.null(ladder_len)) {
          la <- ladder_len[as.character(ek$from)]
          lb <- ladder_len[as.character(ek$to)]
          ok <- is.finite(la) & is.finite(lb) & la >= 3 & lb >= 3
          if (any(ok)) {
            ek[ok, ms2_ladder_confirmed := TRUE]
            ek[ok, edge_score := 1.0]
          }
        }
        edge_list[[length(edge_list) + 1L]] <- ek
      }
    }
  }
  if ("transformation" %in% layers) {
    et <- nodes_dt[
      !is.na(transformation_parent) & nzchar(as.character(transformation_parent)),
      .(
        from = as.character(transformation_parent),
        to = as.character(name),
        edge_type = "Transformation",
        label = {
          r <- as.character(transformation_rule)
          ifelse(is.na(r) | !nzchar(r), "", r)
        }
      )
    ]
    if (nrow(et)) {
      edge_list[[length(edge_list) + 1L]] <- et
    }
  }
  if ("isf" %in% layers) {
    ei <- .nts_build_isf_edges_dt(nodes_dt, isf_mz_tol_ppm)
    if (nrow(ei)) {
      edge_list[[length(edge_list) + 1L]] <- ei
    }
  }
  if (!length(edge_list)) {
    return(ms1_empty)
  }
  E <- data.table::rbindlist(edge_list, use.names = TRUE, fill = TRUE)
  E <- E[, .SD[1L], by = .(from, to, edge_type)]
  keep <- E$from %in% valid_ids & E$to %in% valid_ids
  E[keep]
}


#' Drop join debris and redundant scoring columns from network nodes.
#'
#' Removes \code{i.*} join suffix columns, \code{f_mz}, **all** \code{score_*}
#' columns (including \code{score_cosine_ms2}), and \code{inferred_formula} after
#' formula has been merged. The pruned node table therefore carries **no**
#' scoring columns at all; read scores from \code{@db_annotations} /
#' \code{@ms2_score_grid} instead.
#' @noRd
.nts_network_nodes_prune_columns <- function(nodes_dt) {
  if (!nrow(nodes_dt) || !requireNamespace("data.table", quietly = TRUE)) {
    return(nodes_dt)
  }
  DT <- data.table::as.data.table(nodes_dt)
  nms <- names(DT)
  rmv <- character(0)
  rmv <- c(rmv, nms[grepl("^i\\.", nms)])
  rmv <- c(rmv, nms[nms == "f_mz"])
  
  # Remove all score_* columns including score_cosine_ms2
  rmv <- c(rmv, nms[grepl("^score_", nms)])
  
  # Remove specific columns requested by the user
  requested_removals <- c(
    "ndp_raw", "exp_mz", "lib_mz", "mz_error_ppm",
    "exp_rt", "lib_rt", "rt_error_s", "n_matched_peaks",
    "mzmin", "mzmax", "rtmin", "rtmax",
    "library_Compound_name", "library_ref_id", "library_SMILES"
  )
  rmv <- c(rmv, intersect(nms, requested_removals))
  
  rmv <- c(rmv, intersect(nms, "inferred_formula"))
  rmv <- unique(rmv[nzchar(rmv)])
  if (length(rmv)) {
    hit <- intersect(rmv, names(DT))
    if (length(hit)) {
      DT[, (hit) := NULL]
    }
  }
  DT[]
}


#' @noRd
.nts_network_nodes_dt <- function(project) {
  ann <- .nts_merge_layers_for_network(project)
  ft <- data.table::as.data.table(get_feature_table(project))
  if (!nrow(ann) || !"name" %in% names(ann)) {
    return(list(nodes_dt = ann[0L], ok = FALSE, err = NULL))
  }
  # Distinguish an empty feature_table from a table that is present but missing
  # columns: both used to raise the same "must contain name, mz, rt" message,
  # which misdirected users of an empty project (P3-1).
  if (!nrow(ft)) {
    return(list(
      nodes_dt = ann[0L], ok = FALSE,
      err = "feature_table is empty (no rows); network construction needs at least one feature."
    ))
  }
  if (!all(c("name", "mz", "rt") %in% names(ft))) {
    return(list(nodes_dt = ann[0L], ok = FALSE, err = "feature_table must contain name, mz, rt"))
  }
  ftz <- unique(ft[, .(name, mz = as.numeric(mz), rt = as.numeric(rt))], by = "name")
  
  # Ensure all features in clean_ms2 are represented in ann
  ms2_names <- names(methods::slot(project, "clean_ms2"))
  all_names <- unique(c(ann$name, ms2_names))
  
  if (length(all_names) > nrow(ann)) {
    # Add missing MS2 names to ann
    missing_names <- setdiff(ms2_names, ann$name)
    if (length(missing_names)) {
      missing_dt <- data.table::data.table(name = missing_names)
      ann <- data.table::rbindlist(list(ann, missing_dt), fill = TRUE)
    }
  }

  ann <- ftz[ann, on = "name"]
  if (!"annotation_level" %in% names(ann)) {
    ann[, annotation_level := NA_character_]
  }
  if (!"diagnostic_fragments" %in% names(ann)) {
    ann[, diagnostic_fragments := NA_character_]
  }
  if (!"formula" %in% names(ann)) {
    ann[, formula := NA_character_]
  }
  if (!"inferred_formula" %in% names(ann)) {
    ann[, inferred_formula := NA_character_]
  }
  if (!"kmd_series_id" %in% names(ann)) {
    ann[, kmd_series_id := NA_character_]
  }
  if (!"transformation_parent" %in% names(ann)) {
    ann[, transformation_parent := NA_character_]
  }
  if (!"transformation_rule" %in% names(ann)) {
    ann[, transformation_rule := NA_character_]
  }
  ms2_names <- names(methods::slot(project, "clean_ms2"))
  noise <- grepl("Isolated Unknown", as.character(ann$annotation_level), fixed = TRUE) &
    (is.na(ann$diagnostic_fragments) | !nzchar(as.character(ann$diagnostic_fragments))) &
    !(ann$name %in% ms2_names) # Do NOT filter out nodes that have purified MS2
  nodes_dt <- ann[!noise]
  if (!nrow(nodes_dt)) {
    return(list(nodes_dt = nodes_dt, ok = FALSE, err = NULL))
  }
  nodes_dt[
    ,
    formula_out := {
      fd <- as.character(formula)
      fi <- as.character(inferred_formula)
      out <- ifelse(!is.na(fd) & nzchar(fd), fd, fi)
      ifelse(is.na(out), NA_character_, out)
    }
  ]
  nodes_dt[, id := as.character(name)]
  if (!"is_isf" %in% names(nodes_dt)) {
    nodes_dt[, is_isf := FALSE]
  }
  if (!"mz_remainder" %in% names(nodes_dt)) {
    nodes_dt[, mz_remainder := NA_character_]
  }
  if (!"Propagation_Level" %in% names(nodes_dt)) {
    nodes_dt[, Propagation_Level := NA_integer_]
  }
  if (!"n_fragments" %in% names(nodes_dt)) {
    nodes_dt[, n_fragments := NA_real_]
  }
  if (!"parent_mz" %in% names(nodes_dt)) {
    nodes_dt[, parent_mz := NA_real_]
  }
  if ("formula" %in% names(nodes_dt)) {
    nodes_dt[
      is.na(formula) | !nzchar(trimws(as.character(formula))),
      formula := as.character(formula_out)
    ]
  } else {
    nodes_dt[, formula := as.character(formula_out)]
  }
  nodes_dt[, formula_out := NULL]
  db_raw <- .nts_get_project_db_annotations(project)
  db_df <- .nts_db_best_per_name(db_raw)
  if (is.data.frame(db_df) && nrow(db_df)) {
    db_dt <- data.table::as.data.table(db_df)
    if ("name" %in% names(db_dt)) {
      db_dt[, name := as.character(name)]
    } else if ("feature_id" %in% names(db_dt)) {
      db_dt[, name := as.character(feature_id)]
    }
    if ("name" %in% names(db_dt) && nrow(db_dt)) {
      db_dt <- unique(db_dt, by = "name")
      db_cols <- setdiff(names(db_dt), "name")
      db_cols <- setdiff(db_cols, intersect(c("mz", "rt"), names(db_dt)))
      overlap <- intersect(db_cols, names(nodes_dt))
      overlap <- setdiff(overlap, c("mz", "rt"))
      if (length(overlap)) {
        nodes_dt[, (overlap) := NULL]
      }
      db_join <- db_dt[, c("name", db_cols), with = FALSE]
      nodes_dt <- db_join[nodes_dt, on = "name"]
    }
  }
  if ("ref_id" %in% names(nodes_dt) && !"Library_ID" %in% names(nodes_dt)) {
    nodes_dt[, Library_ID := as.character(ref_id)]
  }
  if ("Compound.name" %in% names(nodes_dt) && !"Annotation" %in% names(nodes_dt)) {
    nodes_dt[, Annotation := as.character(Compound.name)]
  }
  smi_candidates <- intersect(
    c("SMILES", "smiles", "QSAR_READY_SMILES", "Smiles", "canonical_smiles"),
    names(nodes_dt)
  )
  if (length(smi_candidates)) {
    cn <- smi_candidates[1L]
    smi_vec <- as.character(nodes_dt[[cn]])
    nodes_dt[, SMILES := smi_vec]
  } else {
    nodes_dt[, SMILES := NA_character_]
  }
  nodes_dt <- .nts_network_nodes_prune_columns(nodes_dt)
  list(nodes_dt = nodes_dt, ok = TRUE, err = NULL)
}


#' @noRd
.nts_ms2_layer_edges_dt <- function(project,
                                    valid_ids,
                                    ms2_edge_types,
                                    mz_tol_ppm,
                                    min_shared_ms2,
                                    max_bg_freq,
                                    min_frag_freq = 0.5,
                                    md_min_da = 1.0,
                                    precursor_tol = 1.0,
                                    cosine_min_score = 0.6,
                                    shift_max_da = 500,
                                    cosine_mz_tol = 0.01,
                                    cosine_min_peaks = 2L,
                                    min_semantic_score = 0.70,
                                    semantic_top_k = 3L,
                                    embeddings = NULL,
                                    python_path = NULL) {
  types <- unique(as.character(ms2_edge_types))
  need_ms2 <- any(types %in% c(
    "fragment", "neutral_loss", "mass_difference", "cosine", "modified_cosine", "semantic"
  ))
  empty <- .nts_ms2_edge_empty_dt()
  if (!need_ms2) {
    return(empty)
  }

  need_peaks <- any(types %in% c(
    "fragment", "neutral_loss", "mass_difference", "cosine", "modified_cosine"
  ))

  ms2_anchor_names <- character(0)
  if (inherits(project, "NTSProject")) {
    m2df <- tryCatch(get_ms2_annotations(project), error = function(e) NULL)
    if (is.data.frame(m2df) && nrow(m2df) && "name" %in% names(m2df)) {
      ms2_anchor_names <- unique(as.character(m2df$name))
      ms2_anchor_names <- ms2_anchor_names[!is.na(ms2_anchor_names) & nzchar(ms2_anchor_names)]
    }
  }

  parts_ms2 <- list()

  if (need_peaks) {
    .nts_check(mz_tol_ppm, "positive_finite", .arg = "mz_tol_ppm")
    if (!length(ms2_anchor_names)) {
      stop(
        "MS2 network edges require non-empty @ms2_annotations with column `name` ",
        "(seeded by purify_ms2_spectra() or extract_ms2()).",
        call. = FALSE
      )
    }
    clean_l <- methods::slot(project, "clean_ms2")
    if (!is.list(clean_l) || !length(clean_l)) {
      stop(
        "No purified MS2 found; run purify_ms2_spectra() to build clean_ms2 first.",
        call. = FALSE
      )
    }
    nt <- .nts_clean_ms2_network_tables(clean_l, ms2_anchor_names = ms2_anchor_names)
    clean_dt <- nt$clean_dt
    prec_dt <- nt$prec_dt
    if (!nrow(clean_dt)) {
      stop(
        "No purified MS2 peaks for features listed in @ms2_annotations$name ",
        "(check clean_ms2 names / feature_id alignment).",
        call. = FALSE
      )
    }
    clean_dt <- clean_dt[!is.na(name) & nzchar(as.character(name))]
    clean_dt <- clean_dt[is.finite(mz) & is.finite(intensity)]
    mf <- min_frag_freq[1L]
    fr <- if ("frequency" %in% names(clean_dt)) {
      suppressWarnings(as.numeric(clean_dt$frequency))
    } else {
      NULL
    }
    skip_frag_freq <- FALSE
    if (length(fr) && any(is.finite(fr))) {
      rfin <- fr[is.finite(fr)]
      skip_frag_freq <- min(rfin, na.rm = TRUE) >= mf - 1e-12
    }
    if (
      !is.null(fr) &&
        !skip_frag_freq &&
        is.numeric(mf) &&
        is.finite(mf)
    ) {
      clean_dt <- clean_dt[is.na(fr) | fr >= mf]
    }
    if (!nrow(clean_dt)) {
      stop(
        "All consensus peaks were removed by min_frag_freq / NA handling; ",
        "lower min_frag_freq or inspect spectrum `frequency` in @clean_ms2.",
        call. = FALSE
      )
    }
    clean_dt <- .nts_clean_ms2_filter_bg_mz(clean_dt, max_bg_freq)
    if (!nrow(clean_dt)) {
      stop(
        "All peaks removed by ubiquitous-ion filter (max_bg_freq); ",
        "raise max_bg_freq or disable with NULL.",
        call. = FALSE
      )
    }
    clean_nl <- .nts_ms2_clean_with_precursor(clean_dt, prec_dt)
    if (!is.null(max_bg_freq) && is.numeric(max_bg_freq) && is.finite(max_bg_freq[1L]) &&
        max_bg_freq[1L] > 0 && max_bg_freq[1L] < 1 && nrow(clean_nl) && "neutral_loss" %in% names(clean_nl)) {
      mbf_nl <- max_bg_freq[1L]
      clean_nl[, nl_k := round(as.numeric(neutral_loss), 2L)]
      n_feat_nl <- clean_nl[, data.table::uniqueN(name)]
      if (n_feat_nl >= 3L) {
        ub_nl <- clean_nl[, .(nf = data.table::uniqueN(name)), by = nl_k]
        bad_nl <- ub_nl[nf / n_feat_nl > mbf_nl, nl_k]
        if (length(bad_nl)) {
          clean_nl <- clean_nl[!nl_k %in% bad_nl]
        }
      }
      clean_nl[, nl_k := NULL]
    }
    min_k <- suppressWarnings(as.integer(min_shared_ms2)[1L])
    if (length(min_k) != 1L || is.na(min_k) || min_k < 1L) {
      stop(
        "`min_shared_ms2` must be a single positive integer (got ",
        paste(deparse(min_shared_ms2), collapse = " "), ").",
        call. = FALSE
      )
    }

    if ("fragment" %in% types) {
      parts_ms2[[length(parts_ms2) + 1L]] <- .nts_network_ms2_fragment_edges(
        clean_dt,
        mz_tol_ppm,
        min_k
      )
    }
    if ("neutral_loss" %in% types) {
      parts_ms2[[length(parts_ms2) + 1L]] <- .nts_network_ms2_nl_edges(
        clean_nl,
        mz_tol_ppm,
        min_k
      )
    }
    if ("mass_difference" %in% types) {
      parts_ms2[[length(parts_ms2) + 1L]] <- .nts_network_ms2_massdiff_edges(
        clean_dt,
        mz_tol_ppm,
        min_k,
        md_min_da = md_min_da
      )
    }
    if ("cosine" %in% types) {
      parts_ms2[[length(parts_ms2) + 1L]] <- .nts_network_ms2_cosine_edges(
        clean_dt,
        prec_dt,
        min_score = cosine_min_score,
        precursor_tol = precursor_tol,
        modified = FALSE,
        shift_max_da = shift_max_da,
        mz_tol = cosine_mz_tol,
        min_peaks = cosine_min_peaks
      )
    }
    if ("modified_cosine" %in% types) {
      parts_ms2[[length(parts_ms2) + 1L]] <- .nts_network_ms2_cosine_edges(
        clean_dt,
        prec_dt,
        min_score = cosine_min_score,
        precursor_tol = precursor_tol,
        modified = TRUE,
        shift_max_da = shift_max_da,
        mz_tol = cosine_mz_tol,
        min_peaks = cosine_min_peaks
      )
    }
  }

  if ("semantic" %in% types) {
    parts_ms2[[length(parts_ms2) + 1L]] <- .nts_network_ms2_semantic_edges(
      project,
      valid_ids,
      min_score = min_semantic_score,
      top_k = semantic_top_k,
      embeddings = embeddings,
      python_path = python_path
    )
  }

  em <- if (!length(parts_ms2)) {
    empty
  } else {
    data.table::rbindlist(parts_ms2, use.names = TRUE, fill = TRUE)
  }
  if (!nrow(em)) {
    return(empty)
  }
  em[, match_type := NA_character_]
  em[Interaction == "MS2_Shared_Fragment", match_type := "fragment"]
  em[Interaction == "MS2_Shared_Neutral_Loss", match_type := "neutral_loss"]
  em[Interaction == "MS2_Shared_Mass_Diff", match_type := "mass_difference"]
  em[Interaction == "MS2_Cosine", match_type := "cosine"]
  em[Interaction == "MS2_Modified_Cosine", match_type := "modified_cosine"]
  em[Interaction == "MS2_Semantic", match_type := "semantic"]

  valid_ms2 <- if (length(ms2_anchor_names)) intersect(as.character(valid_ids), ms2_anchor_names) else as.character(valid_ids)
  is_sem <- em$Interaction == "MS2_Semantic"
  keep <- ifelse(
    is_sem,
    em$from %in% as.character(valid_ids) & em$to %in% as.character(valid_ids),
    em$from %in% valid_ms2 & em$to %in% valid_ms2
  )
  em[keep]
}


#' Normalise edge-layer tokens for `.nts_build_network_edges`.
#'
#' Canonical MS2 layer ids are `fragment` (DFI, alias `dfi`), `neutral_loss`
#' (aliases `nl`), `mass_difference` (aliases `md`, `massdiff`, `mass_diff`);
#' MS1 ids (`homologue`, `transformation`, `isf`) pass through unchanged.
#' Matching is case-insensitive on the alias table.
#' @noRd
.nts_normalize_network_edge_layers <- function(edge_layers) {
  x <- unique(as.character(edge_layers))
  x <- x[!is.na(x) & nzchar(trimws(x))]
  x <- trimws(x)
  lx <- tolower(x)
  x[lx %in% c("dfi", "fragment")] <- "fragment"
  x[lx %in% c("nl", "neutral_loss", "neutralloss", "neutral loss")] <- "neutral_loss"
  x[lx %in% c(
    "md", "massdiff", "mass_diff", "mass difference", "mass_difference",
    "massdifference"
  )] <- "mass_difference"
  x[lx %in% c("cos", "cosine", "ndp", "ndotproduct")] <- "cosine"
  x[lx %in% c(
    "modcos", "modified_cosine", "modified cosine", "modifiedcosine",
    "modified-cosine"
  )] <- "modified_cosine"
  x[lx %in% c("semantic", "spec2vec", "embedding", "embeddings", "spectral_semantic")] <- "semantic"
  unique(x)
}


#' Normalise MS1 network layer ids for `.nts_ms1_layer_edges_dt`.
#'
#' Canonical ids: `homologue` (aliases `kmd`, `hom`, `homolog`, `series`),
#' `transformation` (aliases `trans`, `transform`), `isf`. Case-insensitive.
#' @noRd
.nts_normalize_ms1_edge_types <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }
  v <- unique(as.character(x))
  v <- v[!is.na(v) & nzchar(trimws(v))]
  v <- trimws(v)
  lv <- tolower(v)
  v[lv %in% c("kmd", "hom", "homologue", "homolog", "series")] <- "homologue"
  v[lv %in% c("trans", "transform", "transformation")] <- "transformation"
  v[lv == "isf"] <- "isf"
  allowed <- c("homologue", "transformation", "isf")
  unique(v[v %in% allowed])
}


#' `NULL` or length-0 input yields no MS1 layers for `type = "all"`.
#' @noRd
.nts_resolve_ms1_edge_types_arg <- function(ms1_edge_types) {
  if (is.null(ms1_edge_types) || !length(ms1_edge_types)) {
    return(character(0))
  }
  .nts_normalize_ms1_edge_types(ms1_edge_types)
}


#' Internal factory: MS1 + MS2 network edge tables from canonical layer ids.
#'
#' @param edge_layers Character vector including any of `"homologue"`,
#'   `"transformation"`, `"isf"`, `"fragment"`, `"neutral_loss"`,
#'   `"mass_difference"` (aliases `"nl"`, `"massdiff"`, `"mass_diff"`),
#'   or `"semantic"`.
#'
#' @return `list(ms1, ms2)` of [data.table::data.table].
#' @noRd
.nts_build_network_edges <- function(project,
                                      nodes_dt,
                                      valid_ids,
                                      edge_layers,
                                      mz_tol_ppm = 10,
                                      min_shared_ms2 = 2,
                                      max_bg_freq = 0.20,
                                      min_frag_freq = 0.5,
                                      isf_mz_tol_ppm = 20,
                                      md_min_da = 1.0,
                                      cf2_ladder_boost = TRUE,
                                      precursor_tol = 1.0,
                                      cosine_min_score = 0.6,
                                      shift_max_da = 500,
                                      cosine_mz_tol = 0.01,
                                      cosine_min_peaks = 2L,
                                      min_semantic_score = 0.70,
                                      semantic_top_k = 3L,
                                      embeddings = NULL,
                                      python_path = NULL) {
  el <- .nts_normalize_network_edge_layers(edge_layers)
  ms1_codes <- c("homologue", "transformation", "isf")
  ms2_codes <- c(
    "fragment", "neutral_loss", "mass_difference", "cosine", "modified_cosine", "semantic"
  )
  ms1_sub <- intersect(el, ms1_codes)
  ms2_sub <- intersect(el, ms2_codes)

  ms1_empty <- data.table::data.table(
    from = character(0),
    to = character(0),
    edge_type = character(0),
    label = character(0)
  )
  out_ms1 <- if (length(ms1_sub)) {
    .nts_ms1_layer_edges_dt(
      nodes_dt,
      valid_ids,
      isf_mz_tol_ppm = isf_mz_tol_ppm,
      ms1_layers = ms1_sub,
      cf2_ladder_boost = cf2_ladder_boost
    )
  } else {
    ms1_empty
  }

  out_ms2 <- if (length(ms2_sub)) {
    .nts_ms2_layer_edges_dt(
      project,
      valid_ids,
      ms2_sub,
      mz_tol_ppm,
      min_shared_ms2,
      max_bg_freq,
      min_frag_freq = min_frag_freq,
      md_min_da = md_min_da,
      precursor_tol = precursor_tol,
      cosine_min_score = cosine_min_score,
      shift_max_da = shift_max_da,
      cosine_mz_tol = cosine_mz_tol,
      cosine_min_peaks = cosine_min_peaks,
      min_semantic_score = min_semantic_score,
      semantic_top_k = semantic_top_k,
      embeddings = embeddings,
      python_path = python_path
    )
  } else {
    .nts_ms2_edge_empty_dt()
  }

  list(ms1 = out_ms1, ms2 = out_ms2)
}


#' @noRd
.nts_merged_network_edge_collapse <- function(edges_dt) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    return(edges_dt)
  }
  E <- data.table::as.data.table(edges_dt)
  if (!nrow(E) || !all(c("from", "to") %in% names(E))) {
    return(edges_dt)
  }
  .pu <- function(x) {
    u <- unique(as.character(x))
    u <- u[!is.na(u) & nzchar(u)]
    if (!length(u)) {
      ""
    } else {
      paste(u, collapse = " + ")
    }
  }
  .fi <- function(x) {
    x <- suppressWarnings(as.integer(x))
    x <- stats::na.omit(x)
    if (!length(x)) {
      NA_integer_
    } else {
      as.integer(x[1L])
    }
  }
  .fc <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(x)]
    if (!length(x)) {
      NA_character_
    } else {
      x[1L]
    }
  }
  .fn <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else x[1L]
  }
  .fx <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    if (!length(x)) NA_real_ else max(x)
  }
  .fl <- function(x) {
    x <- x[!is.na(x)]
    if (!length(x)) NA else any(as.logical(x))
  }
  E[, `:=`(from = as.character(from), to = as.character(to))]
  E[, `:=`(from = pmin(from, to), to = pmax(from, to))]
  need <- c(
    "edge_type",
    "match_type",
    "Interaction",
    "label",
    "shared_fragments_count",
    "shared_fragments_mz",
    "shared_nls_count",
    "shared_nls_mz",
    "shared_mds_count",
    "shared_mds_mz",
    "edge_score",
    "semantic_similarity",
    "n_steps",
    "transformation",
    "mass_difference",
    "ms2_ladder_confirmed"
  )
  miss_nm <- setdiff(need, names(E))
  if (length(miss_nm)) {
    mz_cols <- miss_nm[grepl("_mz$", miss_nm)]
    int_cols <- miss_nm[miss_nm %in% c(
      "shared_fragments_count", "shared_nls_count", "shared_mds_count", "n_steps"
    )]
    num_cols <- miss_nm[miss_nm %in% c("edge_score", "mass_difference", "semantic_similarity")]
    lgl_cols <- miss_nm[miss_nm %in% c("ms2_ladder_confirmed")]
    oth <- setdiff(miss_nm, c(mz_cols, int_cols, num_cols, lgl_cols))
    if (length(mz_cols)) {
      E[, (mz_cols) := NA_character_]
    }
    if (length(int_cols)) {
      E[, (int_cols) := NA_integer_]
    }
    if (length(num_cols)) {
      E[, (num_cols) := NA_real_]
    }
    if (length(lgl_cols)) {
      E[, (lgl_cols) := NA]
    }
    if (length(oth)) {
      E[, (oth) := NA_character_]
    }
  }
  if (!"edge_category" %in% names(E)) {
    E[, edge_category := NA_character_]
  }
  dups <- anyDuplicated(E, by = c("from", "to"))
  if (dups == 0L) {
    return(E)
  }
  is_dup <- duplicated(E, by = c("from", "to")) | duplicated(E, by = c("from", "to"), fromLast = TRUE)
  E_single <- E[!is_dup]
  E_coll <- E[
    is_dup,
    .(
      edge_type = .pu(edge_type),
      match_type = .pu(match_type),
      Interaction = .pu(Interaction),
      label = .pu(label),
      edge_category = .pu(edge_category),
      shared_fragments_count = .fi(shared_fragments_count),
      shared_fragments_mz = .fc(shared_fragments_mz),
      shared_nls_count = .fi(shared_nls_count),
      shared_nls_mz = .fc(shared_nls_mz),
      shared_mds_count = .fi(shared_mds_count),
      shared_mds_mz = .fc(shared_mds_mz),
      edge_score = .fx(edge_score),
      semantic_similarity = .fx(semantic_similarity),
      n_steps = .fi(n_steps),
      transformation = .fc(transformation),
      mass_difference = .fn(mass_difference),
      ms2_ladder_confirmed = .fl(ms2_ladder_confirmed)
    ),
    by = .(from, to)
  ]
  data.table::rbindlist(list(E_single, E_coll), use.names = TRUE, fill = TRUE)
}


#' @noRd
.nts_network_assign_edge_category <- function(E) {
  if (!nrow(E)) {
    return(E)
  }
  E <- data.table::as.data.table(E)
  if (!"edge_category" %in% names(E)) {
    E[, edge_category := NA_character_]
  }
  if ("edge_type" %in% names(E)) {
    et <- as.character(E$edge_type)
    E[et == "Homologue", edge_category := "Homologue"]
    E[et == "Transformation", edge_category := "Transformation"]
    E[et == "ISF", edge_category := "ISF"]
  }
  if ("match_type" %in% names(E)) {
    mt <- as.character(E$match_type)
    E[!is.na(mt) & mt == "fragment", edge_category := "Fragment"]
    E[!is.na(mt) & mt == "neutral_loss", edge_category := "NeutralLoss"]
    E[!is.na(mt) & mt == "mass_difference", edge_category := "MassDifference"]
    E[!is.na(mt) & mt == "cosine", edge_category := "Cosine"]
    E[!is.na(mt) & mt == "modified_cosine", edge_category := "ModifiedCosine"]
    E[!is.na(mt) & mt == "semantic", edge_category := "Semantic"]
  }
  if (!"edge_score" %in% names(E)) {
    E[, edge_score := NA_real_]
  }
  if ("match_type" %in% names(E)) {
    okm <- !is.na(E$match_type) & nzchar(as.character(E$match_type))
    if (any(okm)) {
      # MS2 count edges already carry a normalised edge_score from the builders;
      # only fill rows that have no finite score (fallback = 1).
      cur <- suppressWarnings(as.numeric(E$edge_score))
      fill <- okm & !is.finite(cur)
      if (any(fill)) {
        E[fill, edge_score := 1.0]
      }
    }
  }
  E[]
}


#' Filter edges by physical slice: matches \code{Interaction} tokens (and
#' \code{edge_category} after collapse) so merged Homologue+Fragment rows stay
#' visible when the slice matches either token.
#' @noRd
.nts_network_edges_filter_by_type <- function(E, type) {
  type <- match.arg(
    type,
    c(
      "all", "homologue", "transformation", "fragment", "neutral_loss",
      "mass_diff", "isf", "cosine", "modified_cosine", "semantic"
    )
  )
  if (type == "all" || !nrow(E)) {
    return(E)
  }
  E <- data.table::as.data.table(E)
  intr <- if ("Interaction" %in% names(E)) as.character(E$Interaction) else rep("", nrow(E))
  intr[is.na(intr)] <- ""
  ec <- if ("edge_category" %in% names(E)) as.character(E$edge_category) else rep("", nrow(E))
  ec[is.na(ec)] <- ""
  hay <- paste(intr, ec)
  modcos <- grepl("Modified.?Cosine", hay, ignore.case = TRUE)
  keep <- switch(type,
    homologue = grepl("Homologue", hay, ignore.case = TRUE),
    transformation = grepl("Transformation", hay, ignore.case = TRUE),
    fragment = grepl("Fragment", hay, ignore.case = TRUE),
    neutral_loss = grepl("Neutral_Loss|NeutralLoss", hay, ignore.case = TRUE),
    mass_diff = grepl("Mass_Diff|MassDifference", hay, ignore.case = TRUE),
    isf = grepl("ISF", hay, ignore.case = TRUE),
    cosine = grepl("Cosine", hay, ignore.case = TRUE) & !modcos,
    modified_cosine = modcos,
    semantic = grepl("Semantic", hay, ignore.case = TRUE)
  )
  E[keep]
}


#' Fill empty `edge_type` / `label` on merged edges (MS2-only graphs, collapse).
#' @noRd
.nts_network_edges_finalize_display <- function(E) {
  if (!nrow(E) || !requireNamespace("data.table", quietly = TRUE)) {
    return(E)
  }
  E <- data.table::as.data.table(E)
  if ("match_type" %in% names(E)) {
    w <- (is.na(E$edge_type) | !nzchar(as.character(E$edge_type))) &
      !is.na(E$match_type) & nzchar(as.character(E$match_type))
    if (any(w)) {
      if (!"edge_type" %in% names(E)) {
        E[, edge_type := NA_character_]
      }
      E[w, edge_type := as.character(match_type)]
    }
  }
  if (!"edge_type" %in% names(E)) {
    E[, edge_type := ""]
  }
  E[is.na(edge_type) | !nzchar(as.character(edge_type)), edge_type := ""]
  if ("label" %in% names(E) && "Interaction" %in% names(E)) {
    E[is.na(label) | !nzchar(as.character(label)), label := as.character(Interaction)]
  }
  if ("label" %in% names(E)) {
    E[is.na(label) | !nzchar(as.character(label)), label := ""]
  } else {
    E[, label := ""]
  }
  E[]
}


#' @noRd
.nts_network_prune_isolated_nodes <- function(Nodes, Edges) {
  if (!nrow(Nodes)) {
    return(Nodes)
  }
  if (!nrow(Edges) || !all(c("from", "to") %in% names(Edges))) {
    return(Nodes[0L])
  }
  u <- unique(c(as.character(Edges$from), as.character(Edges$to)))
  Nodes[as.character(Nodes$id) %in% u]
}


#' Scale guards for a collapsed edge table (defaults are no-ops).
#'
#' `top_k` keeps an edge when it ranks among the strongest `top_k` edges of *at
#' least one* of its two endpoints (a union rule, not a strict per-node cap) --
#' see the note below. `max_edges` caps the total (strongest first). Row order is
#' preserved except when `max_edges` truncates. Missing / non-finite `edge_score`
#' sorts last.
#' @noRd
.nts_network_apply_edge_limits <- function(E, top_k = Inf, max_edges = Inf) {
  if (!nrow(E) || !requireNamespace("data.table", quietly = TRUE)) {
    return(E)
  }
  E <- data.table::as.data.table(E)
  sc <- if ("edge_score" %in% names(E)) suppressWarnings(as.numeric(E$edge_score)) else rep(NA_real_, nrow(E))
  sc[!is.finite(sc)] <- -Inf
  tk <- suppressWarnings(as.numeric(top_k)[1L])
  if (is.finite(tk) && tk >= 1) {
    # Long form: one (node, edge_index, score) row per endpoint, then rank each
    # edge within every node it touches (strongest score first) and keep an edge
    # if it survives on at least one endpoint. This is a *union* rule: a single
    # node may still end up with more than `top_k` edges when several of its
    # edges rank top-`top_k` at their other endpoint. It is retained (rather than
    # a strict per-node cap) because a hard cap fragments the graph; the roxygen
    # for `top_k` documents the union semantics accordingly (P2-1).
    # Plain vectors avoid data.table NSE.
    idx <- seq_len(nrow(E))
    nd <- c(as.character(E$from), as.character(E$to))
    ei <- c(idx, idx)
    ss <- c(sc, sc)
    o <- order(nd, -ss, ei)
    nd <- nd[o]
    ei <- ei[o]
    rk <- as.integer(stats::ave(seq_along(nd), nd, FUN = seq_along))
    keep_i <- sort(unique(ei[rk <= tk]))
    E <- E[keep_i]
    sc <- sc[keep_i]
  }
  me <- suppressWarnings(as.numeric(max_edges)[1L])
  if (is.finite(me) && me >= 0 && nrow(E) > me) {
    ord <- order(-sc, seq_len(nrow(E)))
    keep_i2 <- sort(ord[seq_len(as.integer(me))])
    E <- E[keep_i2]
  }
  E[]
}



# Exported network builders: [generate_network_data] wraps [get_network_data]
# (AGENTS step 7). `get_network_data` appears above so its implementation stays
# in scope before the thin wrapper.


#' @title Extract network layer tables (decoupled from plotting)
#'
#' @description
#' Returns `nodes` and `edges` as **data.table** without mutating the project.
#' Nodes carry evidence-based `group`, optional `SMILES`,
#' and a slimmed annotation view (join debris and most `score_*` columns
#' removed; `score_cosine_ms2` kept when present). Edges carry `edge_category` (`Homologue`,
#' `ISF`, `Fragment`, `NeutralLoss`, `MassDifference`) plus shared MS2 evidence
#' columns. MS1 and MS2 edge sets are merged for `type = "all"`, then parallel
#' edges collapsed; subset `type` filters by `Interaction` / `edge_category`
#' tokens and drops isolated nodes.
#'
#' @param project An [NTSProject].
#' @param type Network slice: `"all"` (default), `"homologue"`,
#'   `"transformation"`, `"isf"`, or one MS2 layer (`fragment`, `neutral_loss`,
#'   `mass_diff`, `cosine`, `modified_cosine`, `semantic`). Slices filter on `Interaction` /
#'   `edge_category` tokens; `"all"` merges MS1 + requested MS2 layers then
#'   collapses parallel edges (preserving shared MS2 evidence columns).
#' @param ms2_edge_types Passed through `.nts_build_network_edges` to the MS2
#'   branch when MS2 edges are built (ignored when `type` requests a single MS2
#'   layer 'then only that layer is computed). May include `"cosine"` (plain
#'   ndotproduct) and `"modified_cosine"` (shift-aware) similarity layers; both
#'   are opt-in and off by default so `type = "all"` behaviour is unchanged.
#' @param ms1_edge_types MS1 edge layers for `type = "all"`: subset of
#'   `c("homologue", "transformation", "isf")`. `NULL` or `character(0)` builds
#'   **no** KMD / transformation / ISF edges (pure MS2 graphs when combined with
#'   `ms2_edge_types`). Slices `type = "homologue"` / `"transformation"` / `"isf"`
#'   always build only that MS1 layer regardless of this argument.
#' @inheritParams generate_network_data
#' @param isf_mz_tol_ppm m/z tolerance (ppm) when linking [flag_isf] rows to
#'   parents for `ISF` edges.
#' @param md_min_da Minimum intra-spectrum fragment spacing (Da) kept when
#'   building `mass_difference` (MD) edges (default `1.0`, was hard-coded).
#' @param min_ms2_score Optional numeric floor on `edge_score` for MS2 edges.
#'   `edge_score` is now a **normalised** weight in `(0, 1]`
#'   (`n_shared / min(peaks_a, peaks_b)` for DFI/NL/MD, `1` for reinforced
#'   homologue edges). Default `NULL` keeps every edge (no filtering); set a
#'   number (e.g. `0.5`) to drop weak MS2 edges.
#' @param top_k,max_edges,min_degree Optional scale guards (defaults are no-ops,
#'   preserving the pre-0.1.25 behaviour). `top_k` keeps an edge when it ranks
#'   among the strongest `top_k` edges of **at least one** of its endpoints (a
#'   *union* rule -- an individual node may therefore keep more than `top_k`
#'   edges); `max_edges` caps the total edge count (strongest first);
#'   `min_degree` peels nodes whose degree is below this threshold and
#'   re-trims the edges to the survivors (an iterative k-core, so the returned
#'   `nodes` / `edges` stay self-consistent).
#' @param cf2_ladder_boost When `TRUE` (default), a Homologue edge whose two
#'   endpoints each carry an MS2 CF2 fragment ladder (`pfas_ladder_len >= 3`) is
#'   flagged `ms2_ladder_confirmed = TRUE` and given `edge_score = 1`. Reinforces
#'   MS1 homologue edges with MS2 evidence; no new feature-level edges.
#' @param precursor_tol Precursor m/z window (Da) within which **plain cosine**
#'   pairs are scored (default `1.0`). Ignored for MS2 layers other than
#'   `"cosine"` / `"modified_cosine"`.
#' @param cosine_min_score Minimum cosine similarity for a `"cosine"` /
#'   `"modified_cosine"` edge (default `0.6`). Unlike DFI / NL / MD these layers
#'   do not use `min_shared_ms2`; `edge_score` is the similarity itself.
#' @param shift_max_da Maximum precursor m/z difference (Da) at which a
#'   **modified-cosine** shift is attempted (default `500`, GNPS-aligned). A
#'   partner pair further apart than this in precursor m/z is not scored.
#' @param cosine_mz_tol Fragment m/z tolerance (Da) for aligning peaks inside the
#'   cosine computation (default `0.01`).
#' @param cosine_min_peaks Minimum consensus peaks a spectrum must have to
#'   participate in cosine / modified-cosine edges (default `2`).
#' @param include_semantic Logical (default `FALSE`). When `TRUE`, builds MS2
#'   spectral semantic similarity edges from 256-d embeddings (Transformer or
#'   external model such as DreaMS/Spec2Vec).
#' @param min_semantic_score Minimum cosine similarity for a semantic edge (default `0.70`).
#' @param semantic_top_k Maximum number of strongest semantic neighbors per node (default `3L`).
#' @param embeddings Optional precomputed embedding matrix or data.frame. When `NULL`
#'   and `include_semantic = TRUE`, retrieved from `@physchem_proxy$ms2_embeddings`
#'   or extracted via [extract_ms2_embeddings].
#' @param python_path Optional path to Python executable for embedding extraction.
#'
#' @return `list(nodes = data.table, edges = data.table)`.
#'
#' @importFrom data.table := as.data.table copy data.table is.data.table rbindlist
#' @export
get_network_data <- function(project,
                             type = c(
                               "all",
                               "homologue",
                               "transformation",
                               "fragment",
                               "neutral_loss",
                               "mass_diff",
                               "isf",
                               "cosine",
                               "modified_cosine",
                               "semantic"
                             ),
                             ms2_edge_types = c("fragment", "neutral_loss", "mass_difference"),
                             ms1_edge_types = c("homologue", "transformation", "isf"),
                             mz_tol_ppm = 10,
                             min_shared_ms2 = 2,
                             max_bg_freq = NULL,
                             min_frag_freq = 0.5,
                             max_frag_freq = NULL,
                             isf_mz_tol_ppm = 20,
                             md_min_da = 1.0,
                             min_ms2_score = NULL,
                             keep_isolated = TRUE,
                             top_k = Inf,
                             max_edges = Inf,
                             min_degree = 0L,
                             cf2_ladder_boost = TRUE,
                             precursor_tol = 1.0,
                             cosine_min_score = 0.6,
                             shift_max_da = 500,
                             cosine_mz_tol = 0.01,
                             cosine_min_peaks = 2L,
                             include_semantic = FALSE,
                             min_semantic_score = 0.70,
                             semantic_top_k = 3L,
                             embeddings = NULL,
                             python_path = NULL) {
  type <- match.arg(type)
  if (!is.null(max_frag_freq)) {
    warning(
      "`max_frag_freq` has been renamed to `min_frag_freq` ",
      "(it keeps consensus peaks with frequency >= the value).",
      call. = FALSE
    )
    min_frag_freq <- max_frag_freq
  }
  .nts_check(project, "data.table_ns", "project", .arg = "project")

  ms1_for_all <- .nts_resolve_ms1_edge_types_arg(ms1_edge_types)
  ms1_layers <- if (identical(type, "all")) {
    ms1_for_all
  } else if (identical(type, "homologue")) {
    c("homologue")
  } else if (identical(type, "transformation")) {
    c("transformation")
  } else if (identical(type, "isf")) {
    c("isf")
  } else {
    character(0)
  }
  build_ms1 <- length(ms1_layers) > 0L
  build_ms2 <- if (identical(type, "all")) {
    length(ms2_edge_types) > 0L || isTRUE(include_semantic)
  } else {
    type %in% c("fragment", "neutral_loss", "mass_diff", "cosine", "modified_cosine", "semantic")
  }

  m2_types <- if (type == "fragment") {
    "fragment"
  } else if (type == "neutral_loss") {
    "neutral_loss"
  } else if (type == "mass_diff") {
    "mass_difference"
  } else if (type == "cosine") {
    "cosine"
  } else if (type == "modified_cosine") {
    "modified_cosine"
  } else if (type == "semantic") {
    "semantic"
  } else {
    cur_m2 <- ms2_edge_types
    if (isTRUE(include_semantic) && !"semantic" %in% cur_m2) {
      cur_m2 <- c(cur_m2, "semantic")
    }
    cur_m2
  }

  # MS2 layers requested: if peak-based layers are requested but the project has no
  # purified MS2, degrade gracefully or stop.
  if (build_ms2) {
    clean_l0 <- methods::slot(project, "clean_ms2")
    peak_m2_types <- c("fragment", "neutral_loss", "mass_difference", "cosine", "modified_cosine")
    need_peak_types <- any(m2_types %in% peak_m2_types)
    has_clean <- is.list(clean_l0) && length(clean_l0) > 0L
    if (need_peak_types && !has_clean) {
      if (build_ms1 || ("semantic" %in% m2_types)) {
        warning(
          "No purified MS2 (@clean_ms2 is empty): skipping peak-based MS2 edges. ",
          "Run purify_ms2_spectra() for MS2 peak edges, or pass ",
          "ms2_edge_types = character(0) to silence this warning.",
          call. = FALSE
        )
        m2_types <- setdiff(m2_types, peak_m2_types)
        if (!length(m2_types)) {
          build_ms2 <- FALSE
        }
      } else {
        stop(
          "No purified MS2 found and no MS1 edge layers requested; ",
          "run purify_ms2_spectra() or request MS1 layers (e.g. type = \"homologue\").",
          call. = FALSE
        )
      }
    }
  }

  nb <- .nts_network_nodes_dt(project)
  if (!isTRUE(nb$ok)) {
    if (!is.null(nb$err)) {
      stop(nb$err, call. = FALSE)
    }
    return(list(
      nodes = data.table::data.table(),
      edges = data.table::data.table()
    ))
  }
  Nodes <- data.table::copy(nb$nodes_dt)
  
  # If we are only building MS2 edges (no MS1 layers), restrict nodes 
  # to the universe of features that actually have purified MS2 spectra or embeddings.
  if (!build_ms1 && build_ms2) {
    ms2_names <- names(methods::slot(project, "clean_ms2"))
    emb_names <- if (!is.null(embeddings)) {
      if (is.matrix(embeddings) && !is.null(rownames(embeddings))) rownames(embeddings)
      else if (is.data.frame(embeddings)) {
        if (!is.null(embeddings$feature_id)) as.character(embeddings$feature_id)
        else if (!is.null(embeddings$name)) as.character(embeddings$name)
        else character(0)
      } else character(0)
    } else character(0)
    allowed_ids <- unique(c(ms2_names, emb_names))
    if (length(allowed_ids)) {
      Nodes <- Nodes[id %in% allowed_ids]
    }
  }

  valid_ids <- unique(Nodes$id)
  grp <- .nts_network_nodes_plot_group(project, Nodes)
  if (length(grp) == nrow(Nodes)) {
    Nodes[, group := grp]
  }

  edge_layers <- character(0)
  if (build_ms1) {
    edge_layers <- c(edge_layers, ms1_layers)
  }
  if (build_ms2) {
    edge_layers <- c(edge_layers, m2_types)
  }
  eb <- .nts_build_network_edges(
    project,
    Nodes,
    valid_ids,
    edge_layers,
    mz_tol_ppm = mz_tol_ppm,
    min_shared_ms2 = min_shared_ms2,
    max_bg_freq = max_bg_freq,
    min_frag_freq = min_frag_freq,
    isf_mz_tol_ppm = isf_mz_tol_ppm,
    md_min_da = md_min_da,
    cf2_ladder_boost = cf2_ladder_boost,
    precursor_tol = precursor_tol,
    cosine_min_score = cosine_min_score,
    shift_max_da = shift_max_da,
    cosine_mz_tol = cosine_mz_tol,
    cosine_min_peaks = cosine_min_peaks,
    min_semantic_score = min_semantic_score,
    semantic_top_k = semantic_top_k,
    embeddings = embeddings,
    python_path = python_path
  )
  e_ms1 <- eb$ms1
  e_ms2 <- eb$ms2
  if (nrow(e_ms1)) {
    e_ms1[, match_type := NA_character_]
    e_ms1[, Interaction := NA_character_]
    e_ms1[edge_type == "Homologue", Interaction := "KMD_Homologue"]
    e_ms1[edge_type == "Transformation", Interaction := "Transformation"]
    e_ms1[edge_type == "ISF", Interaction := "ISF_Link"]
    e_ms1 <- .nts_network_assign_edge_category(e_ms1)
  }
  if (nrow(e_ms2)) {
    if (!"edge_type" %in% names(e_ms2)) {
      e_ms2[, edge_type := NA_character_]
    }
    e_ms2 <- .nts_network_assign_edge_category(e_ms2)
    if (
      !is.null(min_ms2_score) &&
        "edge_score" %in% names(e_ms2) &&
        any(is.finite(e_ms2$edge_score))
    ) {
      msc <- suppressWarnings(as.numeric(min_ms2_score)[1L])
      if (is.finite(msc)) {
        e_ms2 <- e_ms2[is.na(edge_score) | edge_score >= msc]
      }
    }
  }

  e_merged <- data.table::rbindlist(list(e_ms1, e_ms2), use.names = TRUE, fill = TRUE)
  if (nrow(e_merged)) {
    e_merged <- .nts_merged_network_edge_collapse(e_merged)
    e_merged <- .nts_network_edges_filter_by_type(e_merged, type)
    if (nrow(e_merged)) {
      e_merged <- .nts_network_edges_finalize_display(e_merged)
      e_merged <- .nts_network_apply_edge_limits(
        e_merged,
        top_k = top_k,
        max_edges = max_edges
      )
      mdg <- suppressWarnings(as.integer(min_degree)[1L])
      if (is.finite(mdg) && mdg > 0 && nrow(e_merged)) {
        # Iterative k-core peel: repeatedly drop nodes whose degree falls below
        # `min_degree` and re-trim edges to the surviving endpoints, until the
        # node set is stable. This keeps the returned (nodes, edges) pair a
        # self-consistent subgraph. Previously only `Nodes` was filtered and
        # `e_merged` was returned whole, so edges still referenced dropped ids
        # (dangling edges) and the node/edge invariants broke (P1-1).
        fr <- as.character(e_merged$from)
        to <- as.character(e_merged$to)
        keep_ids <- as.character(Nodes$id)
        repeat {
          active <- fr %in% keep_ids & to %in% keep_ids
          deg <- table(c(fr[active], to[active]))
          new_keep <- names(deg)[deg >= mdg]
          if (identical(sort(new_keep), sort(keep_ids))) break
          keep_ids <- new_keep
        }
        e_merged <- e_merged[fr %in% keep_ids & to %in% keep_ids]
        Nodes <- Nodes[as.character(id) %in% keep_ids]
      } else if (!keep_isolated) {
        Nodes <- .nts_network_prune_isolated_nodes(Nodes, e_merged)
      }
    } else if (!identical(type, "all") && !keep_isolated) {
      Nodes <- Nodes[0L]
    }
  }

  list(nodes = Nodes, edges = e_merged)
}


#' @title Build merged MS1+MS2 network tables for dual-network workflows
#'
#' @description
#' Generates visNetwork-style `nodes` and `edges` tables from an [NTSProject]
#' via [get_network_data]: KMD homologue chains, MS1 transformations/ISF, and
#' optional MS2 fragment / neutral-loss / mass-difference layers using
#' `@clean_ms2` from [purify_ms2_spectra]. Does **not** mutate the project
#' object; returns plain `data.frame`s suitable for [annotate_by_dual_network]
#' and [plot_molecular_network].
#'
#' @param project An [NTSProject].
#' @param ms2_edge_types Which MS2-derived edge layers to build: `"fragment"`
#'   (shared fragment m/z), `"neutral_loss"` (shared neutral losses vs
#'   precursor), `"mass_difference"` (shared intra-spectrum mass differences),
#'   plus the opt-in similarity layers `"cosine"` (plain ndotproduct) and
#'   `"modified_cosine"` (shift-aware). If any of these are selected,
#'   `@clean_ms2` must be non-empty.
#' @param ms1_edge_types MS1-derived edge layers for the merged graph: subset of
#'   `c("homologue", "transformation", "isf")`. `NULL` or `character(0)` builds
#'   **no** MS1 edges (pure MS2 graph when combined with non-empty
#'   `ms2_edge_types`). Default keeps all three layers.
#' @param mz_tol_ppm Ppm tolerance for fragment, loss, and mass-difference
#'   matching windows.
#' @param min_shared_ms2 Minimum number of shared items (per edge type) required
#'   to draw an MS2 edge between two features.
#' @param max_bg_freq Optional numeric in `(0, 1)`. When the number of features
#'   with purified MS2 is at least 3, drop m/z bins present in more than this
#'   fraction of features (ubiquitous / background ions) before MS2 edge
#'   building. Default `NULL` disables filtering.
#' @param min_frag_freq Numeric (default `0.5`). Per-consensus-peak **floor** on
#'   `spectrum$frequency` in `@clean_ms2` before MS2 fragment / NL / mass-diff
#'   mining: keep peaks with `frequency >= min_frag_freq`, or with missing /
#'   non-finite `frequency` (retained). Low-recurrence sporadic ions are dropped
#'   so networking favours reproducible (e.g. diagnostic) fragments. If every
#'   finite `frequency` already meets the floor, the filter is skipped. Use
#'   `NULL` or a non-finite value to disable the filter. (Formerly named
#'   `max_frag_freq`; the old name is still accepted with a deprecation warning.)
#' @param max_frag_freq Deprecated alias for `min_frag_freq`. Accepted for
#'   back-compat; when supplied it overrides `min_frag_freq` and emits a
#'   deprecation warning.
#' @param isf_mz_tol_ppm,min_ms2_score Passed to [get_network_data] when building
#'   the merged MS1+MS2 graph (ISF linking tolerance and optional MS2 edge score
#'   floor). `min_ms2_score` defaults to `NULL` (no filtering).
#' @param md_min_da,top_k,max_edges,min_degree,cf2_ladder_boost Passed to
#'   [get_network_data] (MD minimum spacing, per-node / total edge caps, minimum
#'   node degree, and CF2-ladder reinforcement of homologue edges). All defaults
#'   are no-ops except `cf2_ladder_boost`.
#' @param precursor_tol,cosine_min_score,shift_max_da,cosine_mz_tol,cosine_min_peaks
#'   Passed to [get_network_data] for the opt-in `"cosine"` / `"modified_cosine"`
#'   MS2 similarity layers (precursor window, similarity floor, max shift, peak
#'   m/z tolerance, minimum peaks). Inert unless those layer ids are requested in
#'   `ms2_edge_types`.
#' @param include_semantic,min_semantic_score,semantic_top_k,embeddings,python_path
#'   Passed to [get_network_data] for MS2 spectral semantic similarity edges
#'   (opt-in flag, similarity floor, top-k neighbors cap, precomputed embeddings,
#'   and Python interpreter path).
#' @param keep_isolated Logical (default `TRUE`). If `TRUE`, keeps nodes that
#'   have no edges. For molecular networks, it is recommended to keep these
#'   to see all features with purified MS2.
#'
#' @importFrom data.table := as.data.table copy data.table is.data.table rbindlist
#'
#' @return Named list: `nodes` and `edges` as plain `data.frame`s. `nodes` keeps
#'   merged MS1/MS2 annotation columns, best-per-feature rows from
#'   `@db_annotations` (e.g. `ref_id`, `Compound.name`; `score_cosine_ms2` kept
#'   when present, other `score_*` dropped), optional aliases `Library_ID` and
#'   `Annotation` for [annotate_by_dual_network] seeds, plus `id`, `mz`, `rt`,
#'   `formula`. MS2-derived `edges` may include
#'   `shared_fragments_*`, `shared_nls_*`, and `shared_mds_*` (mass-difference
#'   MS2 edges) for tooltips and interpretation, plus the normalised `edge_score`
#'   and (MS1) `mass_difference` / `n_steps` / `transformation` /
#'   `ms2_ladder_confirmed` columns.
#' @export
generate_network_data <- function(project,
                                  ms2_edge_types = c("fragment", "neutral_loss", "mass_difference"),
                                  ms1_edge_types = c("homologue", "transformation", "isf"),
                                  mz_tol_ppm = 10,
                                  min_shared_ms2 = 2,
                                  max_bg_freq = NULL,
                                  min_frag_freq = 0.5,
                                  max_frag_freq = NULL,
                                  isf_mz_tol_ppm = 20,
                                  md_min_da = 1.0,
                                  min_ms2_score = NULL,
                                  keep_isolated = TRUE,
                                  top_k = Inf,
                                  max_edges = Inf,
                                  min_degree = 0L,
                                  cf2_ladder_boost = TRUE,
                                  precursor_tol = 1.0,
                                  cosine_min_score = 0.6,
                                  shift_max_da = 500,
                                  cosine_mz_tol = 0.01,
                                  cosine_min_peaks = 2L,
                                  include_semantic = FALSE,
                                  min_semantic_score = 0.70,
                                  semantic_top_k = 3L,
                                  embeddings = NULL,
                                  python_path = NULL) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")

  nd <- tryCatch(
    get_network_data(
      project,
      type = "all",
      ms2_edge_types = ms2_edge_types,
      ms1_edge_types = ms1_edge_types,
      mz_tol_ppm = mz_tol_ppm,
      min_shared_ms2 = min_shared_ms2,
      max_bg_freq = max_bg_freq,
      min_frag_freq = min_frag_freq,
      max_frag_freq = max_frag_freq,
      isf_mz_tol_ppm = isf_mz_tol_ppm,
      md_min_da = md_min_da,
      min_ms2_score = min_ms2_score,
      keep_isolated = keep_isolated,
      top_k = top_k,
      max_edges = max_edges,
      min_degree = min_degree,
      cf2_ladder_boost = cf2_ladder_boost,
      precursor_tol = precursor_tol,
      cosine_min_score = cosine_min_score,
      shift_max_da = shift_max_da,
      cosine_mz_tol = cosine_mz_tol,
      cosine_min_peaks = cosine_min_peaks,
      include_semantic = include_semantic,
      min_semantic_score = min_semantic_score,
      semantic_top_k = semantic_top_k,
      embeddings = embeddings,
      python_path = python_path
    ),
    error = function(e) e
  )
  if (inherits(nd, "error")) {
    stop(nd$message, call. = FALSE)
  }
  Nodes <- nd$nodes
  Edges <- nd$edges
  if (!nrow(Nodes)) {
    return(list(nodes = data.frame(), edges = data.frame()))
  }

  if (!nrow(Edges)) {
    Edges <- data.table::data.table(
      from = character(0),
      to = character(0),
      edge_type = character(0),
      match_type = character(0),
      Interaction = character(0),
      label = character(0),
      edge_category = character(0),
      shared_fragments_count = integer(0),
      shared_fragments_mz = character(0),
      shared_nls_count = integer(0),
      shared_nls_mz = character(0),
      shared_mds_count = integer(0),
      shared_mds_mz = character(0),
      edge_score = numeric(0),
      semantic_similarity = numeric(0),
      n_steps = integer(0),
      transformation = character(0),
      mass_difference = numeric(0),
      ms2_ladder_confirmed = logical(0)
    )
  }

  list(
    nodes = as.data.frame(Nodes, stringsAsFactors = FALSE),
    edges = as.data.frame(Edges, stringsAsFactors = FALSE)
  )
}


#' Tag ego-centre rows on a subnet node table (`group = "Subnet_Root"`, colors).
#' @noRd
.nts_subnetwork_tag_root_nodes <- function(kn, target_ok) {
  if (!nrow(kn) || !length(target_ok)) {
    return(kn)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    return(kn)
  }
  kn <- data.table::as.data.table(kn)
  roots <- unique(as.character(target_ok[!is.na(target_ok) & nzchar(target_ok)]))
  if (!length(roots)) {
    return(kn)
  }
  kn[, id := as.character(id)]
  hit <- kn$id %in% roots
  if ("name" %in% names(kn)) {
    hit <- hit | (as.character(kn$name) %in% roots)
  }
  if (!any(hit)) {
    return(kn)
  }
  if (!"group" %in% names(kn)) {
    kn[, group := NA_character_]
  }
  kn[hit, group := "Subnet_Root"]
  kn[hit, `:=`(color.background = "#FF8C00", color.border = "#CC7000")]
  if (!"color" %in% names(kn)) {
    kn[, color := NA_character_]
  }
  kn[hit, color := "#FF8C00"]
  kn[]
}


#' @noRd
.nts_subnetwork_format_shared_mz_line_html <- function(s) {
  s <- trimws(as.character(s)[1L])
  if (is.na(s) || !nzchar(s)) {
    return("")
  }
  bits <- trimws(strsplit(s, ";\\s*", perl = TRUE)[[1L]])
  bits <- bits[!is.na(bits) & nzchar(bits)]
  if (!length(bits)) {
    return("")
  }
  paste0("<b>Shared m/z:</b> ", paste(bits, collapse = ", "), "<br>")
}


#' @noRd
.nts_subnetwork_edge_is_fragment <- function(ec, mt, intr) {
  (!is.na(ec) & ec == "Fragment") |
    (!is.na(mt) & grepl("^fragment$", mt, ignore.case = TRUE)) |
    grepl("MS2_Shared_Fragment", intr, fixed = TRUE)
}


#' Enrich `edges$title` (fragment shared m/z summary) and clear `edges$label` for visNetwork.
#' @noRd
.nts_subnetwork_enrich_edges_visual <- function(E) {
  if (!nrow(E)) {
    return(E)
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    return(E)
  }
  E <- data.table::as.data.table(E)
  if (!"title" %in% names(E)) {
    E[, title := NA_character_]
  }
  E[, title := ifelse(is.na(title), "", as.character(title))]
  if (!"label" %in% names(E)) {
    E[, label := NA_character_]
  }
  E[, lab0 := ifelse(is.na(label), "", as.character(label))]

  intr <- if ("Interaction" %in% names(E)) {
    as.character(E$Interaction)
  } else {
    rep("", nrow(E))
  }
  ec <- if ("edge_category" %in% names(E)) {
    as.character(E$edge_category)
  } else {
    rep(NA_character_, nrow(E))
  }
  mt <- if ("match_type" %in% names(E)) {
    as.character(E$match_type)
  } else {
    rep(NA_character_, nrow(E))
  }

  is_frag <- .nts_subnetwork_edge_is_fragment(ec, mt, intr)
  is_frag[is.na(is_frag)] <- FALSE

  sfc <- rep(NA_integer_, nrow(E))
  if ("shared_fragments_count" %in% names(E)) {
    sfc <- suppressWarnings(as.integer(E$shared_fragments_count))
  }

  frag_src <- rep(NA_character_, nrow(E))
  if ("matching_mz_list" %in% names(E)) {
    frag_src <- as.character(E$matching_mz_list)
  }
  if ("shared_fragments_mz" %in% names(E)) {
    sf <- as.character(E$shared_fragments_mz)
    frag_src <- data.table::fifelse(!is.na(sf) & nzchar(sf), sf, frag_src)
  }

  ttl <- E$title
  lab <- E$lab0
  for (i in seq_len(nrow(E))) {
    t_i <- ttl[i]
    if (identical(is_frag[i], TRUE)) {
      has_desc <- grepl("shared\\s+m/z", t_i, ignore.case = TRUE, perl = TRUE)
      if (!has_desc) {
        ni <- sfc[i]
        fsi <- frag_src[i]
        if (!is.finite(ni) || is.na(ni) || ni < 1L) {
          bits <- if (!is.na(fsi) && nzchar(fsi)) {
            trimws(strsplit(fsi, ";\\s*", perl = TRUE)[[1L]])
          } else {
            character(0)
          }
          bits <- bits[!is.na(bits) & nzchar(bits)]
          ni <- length(bits)
        }
        line_n <- if (is.finite(ni) && !is.na(ni) && ni >= 1L) {
          paste0("<b>", ni, " shared m/z</b><br>")
        } else {
          ""
        }
        fz <- .nts_subnetwork_format_shared_mz_line_html(fsi)
        if (!nzchar(fz) && nzchar(line_n)) {
          fz <- "<b>Shared m/z:</b> <i>(not recorded)</i><br>"
        }
        block <- paste0(line_n, fz)
        t_i <- if (nzchar(trimws(t_i))) paste0(block, t_i) else block
      }
    }
    li <- lab[i]
    if (nzchar(li)) {
      if (!grepl(li, t_i, fixed = TRUE)) {
        sep <- if (nzchar(trimws(t_i))) "<br>" else ""
        t_i <- paste0(t_i, sep, "<b>Label</b><br>", li)
      }
    }
    ttl[i] <- t_i
  }
  E[, title := ttl]
  E[, label := ""]
  E[, lab0 := NULL]
  E[]
}


#' Extract ego-network (igraph `make_ego_graph`)
#'
#' Uses [igraph::make_ego_graph] on an undirected graph from `edges`, then keeps
#' the union of induced subgraphs centred on each valid `target_id`.
#' Resolved ego centres are tagged on `nodes` with `group = "Subnet_Root"` and
#' orange `color` / `color.background` so [plot_molecular_network] can render
#' them without extra metadata.
#'
#' @param network_list List with `nodes` and `edges` (`data.frame` or
#'   `data.table`; `nodes` must have column `id`, `edges` columns `from`, `to`).
#' @param target_id Character vector of node `id` values (ego centres).
#' @param max_steps Radius `order` passed to [igraph::make_ego_graph] (`0` =
#'   centres only, no neighbour edges).
#'
#' @return `list(nodes, edges, subnet_meta)` where `subnet_meta$subnet_roots`
#'   repeats the resolved ego-centre node ids (`target_id` entries that lie in
#'   the graph) for downstream plotting (e.g. [plot_molecular_network]).
#'   Subnet `edges` have `label` cleared to `""` and `title` enriched (fragment
#'   edges: shared m/z count and list when missing) so [plot_molecular_network]
#'   can render without edge captions while keeping hover text.
#' @export
extract_subnetwork <- function(network_list, target_id, max_steps = 3L) {
  .nts_check(NULL, "data.table_ns", "igraph_ns", .arg = "igraph")
  if (!is.list(network_list) || !all(c("nodes", "edges") %in% names(network_list))) {
    stop("network_list must be a list with components nodes and edges", call. = FALSE)
  }
  nodes_df <- network_list$nodes
  edges_df <- network_list$edges
  if (!is.data.frame(nodes_df) || !"id" %in% names(nodes_df)) {
    stop("network_list$nodes must be a data.frame-like object with column id", call. = FALSE)
  }
  if (!is.data.frame(edges_df)) {
    edges_df <- data.frame(
      from = character(0),
      to = character(0),
      stringsAsFactors = FALSE
    )
  }
  if (!all(c("from", "to") %in% names(edges_df))) {
    stop("network_list$edges must contain columns from and to", call. = FALSE)
  }

  target_id <- unique(as.character(target_id))
  target_id <- target_id[!is.na(target_id) & nzchar(target_id)]
  valid_ids <- unique(as.character(nodes_df$id))
  ok <- target_id %in% valid_ids
  if (length(target_id) && !all(ok)) {
    bad <- target_id[!ok]
    warning(
      sprintf(
        "extract_subnetwork: ignoring %d target_id(s) not in nodes$id: %s",
        length(bad),
        paste(utils::head(bad, 10L), collapse = ", ")
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  target_ok <- target_id[ok]
  k <- suppressWarnings(as.integer(max_steps)[1L])
  if (!is.finite(k) || k < 0L) {
    stop("max_steps must be a non-negative integer", call. = FALSE)
  }

  empty_nodes <- data.table::as.data.table(nodes_df)[0L]
  empty_edges <- data.table::as.data.table(edges_df)[0L]
  meta_roots <- function(roots) {
    list(subnet_roots = unique(as.character(roots[!is.na(roots) & nzchar(roots)])))
  }
  if (!length(target_ok)) {
    return(list(
      nodes = empty_nodes,
      edges = empty_edges,
      subnet_meta = meta_roots(character(0))
    ))
  }

  v_all <- unique(as.character(nodes_df$id))
  if (!nrow(edges_df)) {
    N0 <- data.table::as.data.table(nodes_df)
    N0[, id := as.character(id)]
    kn <- .nts_subnetwork_tag_root_nodes(N0[id %in% target_ok], target_ok)
    return(list(
      nodes = kn,
      edges = empty_edges,
      subnet_meta = meta_roots(target_ok)
    ))
  }

  el <- data.table::as.data.table(edges_df)[, .(from = as.character(from), to = as.character(to))]
  g <- igraph::graph_from_data_frame(
    d = as.data.frame(el, stringsAsFactors = FALSE),
    directed = FALSE,
    vertices = v_all
  )
  vnames <- igraph::V(g)$name
  vhit <- match(target_ok, vnames)
  vhit <- vhit[is.finite(vhit)]
  if (!length(vhit)) {
    return(list(
      nodes = empty_nodes,
      edges = empty_edges,
      subnet_meta = meta_roots(target_ok)
    ))
  }

  ego_l <- igraph::make_ego_graph(g, order = k, nodes = vhit, mode = "all", mindist = 0L)
  vu <- unique(unlist(lapply(ego_l, function(h) igraph::V(h)$name)))
  if (!length(vu)) {
    return(list(
      nodes = empty_nodes,
      edges = empty_edges,
      subnet_meta = meta_roots(target_ok)
    ))
  }

  N <- data.table::as.data.table(nodes_df)
  N[, id := as.character(id)]
  kn <- .nts_subnetwork_tag_root_nodes(N[id %in% vu], target_ok)
  E <- data.table::as.data.table(edges_df)
  E[, `:=`(from = as.character(from), to = as.character(to))]
  ke <- E[from %in% vu & to %in% vu]
  ke <- .nts_subnetwork_enrich_edges_visual(ke)
  list(
    nodes = kn,
    edges = ke,
    subnet_meta = meta_roots(target_ok)
  )
}


#' Filter network by seed nodes and hop radius (BFS)
#'
#' Extracts an induced subgraph around one or more seeds by expanding along
#' undirected edges for `max_step` rounds (breadth-first), then keeps all edges
#' whose endpoints lie in the accumulated node set. Seeds are resolved against
#' `nodes$id` and, when present, `nodes$name`.
#'
#' @param network_data List with `nodes` and `edges` (e.g. from
#'   [generate_network_data] or [get_network_data]).
#' @param nodes Character vector of seed identifiers (feature `id` or `name`).
#' @param max_step Non-negative integer number of BFS expansion rounds from the
#'   current node set. Use `0` to keep only seeds (and edges with both ends in
#'   the seed set).
#'
#' @return `list(nodes = ..., edges = ...)` using the same table classes as the
#'   input (`data.frame` vs [data.table::data.table]) for compatibility with
#'   [plot_molecular_network]. Subnet `edges` use the same visual sanitization as
#'   [extract_subnetwork] (`label` cleared, `title` enriched for hover).
#'
#' @export
filter_network <- function(network_data, nodes, max_step = 1L) {
  if (!is.list(network_data) || !all(c("nodes", "edges") %in% names(network_data))) {
    stop("network_data must be a list with components nodes and edges", call. = FALSE)
  }
  nodes_df <- network_data$nodes
  edges_df <- network_data$edges
  if (!is.data.frame(nodes_df) || !"id" %in% names(nodes_df)) {
    stop("network_data$nodes must be a data.frame-like object with column id", call. = FALSE)
  }
  if (!is.data.frame(edges_df)) {
    edges_df <- data.frame(
      from = character(0),
      to = character(0),
      stringsAsFactors = FALSE
    )
  }
  if (!all(c("from", "to") %in% names(edges_df))) {
    stop("network_data$edges must contain columns from and to", call. = FALSE)
  }

  use_df <- !data.table::is.data.table(nodes_df)
  Ndt <- data.table::as.data.table(nodes_df)
  Edt <- data.table::as.data.table(edges_df)
  Ndt[, id := as.character(id)]
  Edt[, `:=`(from = as.character(from), to = as.character(to))]

  empty_out <- function(roots = character(0)) {
    n0 <- Ndt[0L]
    e0 <- Edt[0L]
    if (use_df) {
      n0 <- as.data.frame(n0, stringsAsFactors = FALSE)
      e0 <- as.data.frame(e0, stringsAsFactors = FALSE)
    }
    list(
      nodes = n0,
      edges = e0,
      subnet_meta = list(subnet_roots = unique(as.character(roots[!is.na(roots) & nzchar(roots)])))
    )
  }

  req <- unique(as.character(nodes))
  req <- req[!is.na(req) & nzchar(req)]
  if (!length(req)) {
    warning(
      "filter_network: `nodes` is empty; returning empty nodes and edges.",
      call. = FALSE,
      immediate. = TRUE
    )
    return(empty_out(character(0)))
  }

  idc <- Ndt$id
  has_name <- "name" %in% names(Ndt)
  nm <- if (has_name) as.character(Ndt$name) else NULL
  hits <- character(0)
  for (r in req) {
    if (r %in% idc) {
      hits <- c(hits, r)
    } else if (has_name) {
      w <- which(nm == r)
      if (length(w)) {
        hits <- c(hits, idc[w])
      }
    }
  }
  hits <- unique(hits)

  miss <- req[!vapply(req, function(r) {
    r %in% idc || (has_name && r %in% nm)
  }, logical(1L))]
  if (length(miss)) {
    warning(
      sprintf(
        "filter_network: requested node(s) not found in network_data$nodes (id or name): %s",
        paste(utils::head(miss, 15L), collapse = ", ")
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  if (!length(hits)) {
    warning(
      "filter_network: no valid seed nodes; returning empty nodes and edges.",
      call. = FALSE,
      immediate. = TRUE
    )
    return(empty_out(character(0)))
  }

  k <- suppressWarnings(as.integer(max_step)[1L])
  if (!is.finite(k) || k < 0L) {
    stop("max_step must be a non-negative integer", call. = FALSE)
  }

  V <- unique(as.character(hits))
  if (nrow(Edt) && k > 0L) {
    for (.step in seq_len(k)) {
      rel <- Edt[from %in% V | to %in% V]
      if (!nrow(rel)) {
        break
      }
      V <- unique(c(V, as.character(rel$from), as.character(rel$to)))
    }
  }

  kn <- .nts_subnetwork_tag_root_nodes(Ndt[id %in% V], unique(as.character(hits)))
  ke <- Edt[from %in% V & to %in% V]
  ke <- .nts_subnetwork_enrich_edges_visual(ke)
  if (use_df) {
    kn <- as.data.frame(kn, stringsAsFactors = FALSE)
    ke <- as.data.frame(ke, stringsAsFactors = FALSE)
  }
  list(
    nodes = kn,
    edges = ke,
    subnet_meta = list(subnet_roots = unique(as.character(hits)))
  )
}
