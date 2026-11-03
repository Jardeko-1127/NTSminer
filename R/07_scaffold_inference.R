# =============================================================================
# 07_scaffold_inference.R —— 无库未知物类似物与化学骨架自动推导 (DreaMS 范式)
# -----------------------------------------------------------------------------
# 目的
#   针对非靶向质谱中海量查无此人、无直接谱库命中的未知污染物（Level 4 / Level 5），
#   借助 Transformer / DreaMS / Spec2Vec 256 维质谱语义空间的高维连续相似度：
#     1. 极速检索其 Top-K 语义最近邻已知类似物（高置信度同批特征或参考库化合物）；
#     2. 计算精确质量差（Δm/z），比对通用环境化学转化字典（+O, -CO2, ±CF2, ±CH2, ±SO3 等）；
#     3. 提取类似物的 Bemis-Murcko Scaffold（环系母核）与 ClassyFire / 化学分类信息；
#     4. 综合推导出无库未知物的母核类别、类似物关系与结构修饰标签。
# =============================================================================

#' Default environmental transformation and mass shift dictionary
#'
#' Provides an extensive, domain-general dictionary of common chemical
#' modifications, environmental transformations, homologue series steps,
#' and phase I/II metabolic conjugations.
#'
#' @return A [data.table::data.table] with columns:
#'   \itemize{
#'     \item \code{rule_name}: Descriptive name of the transformation.
#'     \item \code{mass_shift}: Exact monoisotopic mass difference in Daltons.
#'     \item \code{delta_formula}: Chemical formula difference.
#'     \item \code{category}: Transformation category.
#'   }
#' @export
default_transformation_rules <- function() {
  .nts_default_environmental_transformations()
}

#' @keywords internal
#' @noRd
.nts_default_environmental_transformations <- function() {
  data.table::data.table(
    rule_name = c(
      "Isomer / Identical mass",
      "+O (Hydroxylation)",
      "+2O (Dihydroxylation)",
      "-O (Deoxygenation)",
      "-CO2 (Decarboxylation)",
      "+CH2 (Alkyl chain homologue)",
      "+2CH2 (Ethyl chain homologue)",
      "-CH2 (Demethylation)",
      "+CF2 (Perfluoroalkyl chain homologue)",
      "+2CF2 (Perfluoroalkyl 2-step)",
      "-CF2 (Defluorination step)",
      "-2CF2 (Defluorination 2-step)",
      "+SO3 (Sulfonation)",
      "-SO3 (Desulfonation)",
      "+H2O (Hydration)",
      "-H2O (Dehydration)",
      "-H2 (Desaturation / Oxidation)",
      "+H2 (Hydrogenation / Reduction)",
      "+Cl - H (Chlorination)",
      "-Cl + H (Dechlorination)",
      "+Br - H (Bromination)",
      "-Br + H (Debromination)",
      "-F + OH (F Hydrolysis)",
      "+C2F4O (Perfluoroether addition)",
      "-C2F4O (Perfluoroether cleavage)",
      "+C3F6O (HFPO addition)",
      "+Glucuronide (+C6H8O6)",
      "+Acetyl (+C2H2O)"
    ),
    mass_shift = c(
      0.000000,
      15.994915,
      31.989829,
      -15.994915,
      -43.989829,
      14.015650,
      28.031300,
      -14.015650,
      49.996806,
      99.993613,
      -49.996806,
      -99.993613,
      79.956815,
      -79.956815,
      18.010565,
      -18.010565,
      -2.015650,
      2.015650,
      33.961028,
      -33.961028,
      77.910512,
      -77.910512,
      -1.995663,
      115.988527,
      -115.988527,
      165.985333,
      176.032089,
      42.010565
    ),
    delta_formula = c(
      "0",
      "+O",
      "+O2",
      "-O",
      "-CO2",
      "+CH2",
      "+C2H4",
      "-CH2",
      "+CF2",
      "+C2F4",
      "-CF2",
      "-C2F4",
      "+SO3",
      "-SO3",
      "+H2O",
      "-H2O",
      "-H2",
      "+H2",
      "+Cl-H",
      "-Cl+H",
      "+Br-H",
      "-Br+H",
      "-F+OH",
      "+C2F4O",
      "-C2F4O",
      "+C3F6O",
      "+C6H8O6",
      "+C2H2O"
    ),
    category = c(
      "Isomer",
      "Phase_I_Metabolism",
      "Phase_I_Metabolism",
      "Phase_I_Metabolism",
      "Environmental_Transformation",
      "Homologue_Chain",
      "Homologue_Chain",
      "Homologue_Chain",
      "Homologue_Chain",
      "Homologue_Chain",
      "Homologue_Chain",
      "Homologue_Chain",
      "Conjugation",
      "Environmental_Transformation",
      "Hydration_Dehydration",
      "Hydration_Dehydration",
      "Redox",
      "Redox",
      "Halogenation",
      "Halogenation",
      "Halogenation",
      "Halogenation",
      "Hydrolysis",
      "Ether_Modification",
      "Ether_Modification",
      "Polyether_Chain",
      "Phase_II_Conjugation",
      "Phase_II_Conjugation"
    )
  )
}

#' @keywords internal
#' @noRd
.nts_match_transformation <- function(delta_mz, rules_dt = NULL, mz_tol_ppm = 15.0) {
  if (is.null(rules_dt) || !nrow(rules_dt)) {
    rules_dt <- .nts_default_environmental_transformations()
  }
  if (is.null(delta_mz) || length(delta_mz) != 1L || !is.numeric(delta_mz) || !is.finite(delta_mz)) {
    return(list(
      matched = FALSE,
      rule_name = NA_character_,
      delta_formula = NA_character_,
      delta_ppm = NA_real_,
      category = NA_character_
    ))
  }

  # 1. Check for isomer / identical mass
  if (abs(delta_mz) < 0.005) {
    return(list(
      matched = TRUE,
      rule_name = "Isomer / Identical mass",
      delta_formula = "0",
      delta_ppm = round(abs(delta_mz) / 100 * 1e6, 2L),
      category = "Isomer"
    ))
  }

  diffs <- abs(delta_mz - rules_dt$mass_shift)
  # Calculate ppm error relative to mass shift (or 100 Da floor for small shifts)
  ppm_errors <- (diffs / pmax(100.0, abs(rules_dt$mass_shift))) * 1e6
  
  # Allow match if ppm <= mz_tol_ppm OR absolute difference <= 0.005 Da
  match_mask <- (ppm_errors <= mz_tol_ppm) | (diffs <= 0.005)
  
  if (any(match_mask)) {
    best_idx <- which(match_mask)[which.min(diffs[match_mask])]
    return(list(
      matched = TRUE,
      rule_name = as.character(rules_dt$rule_name[best_idx]),
      delta_formula = as.character(rules_dt$delta_formula[best_idx]),
      delta_ppm = round(ppm_errors[best_idx], 2L),
      category = as.character(rules_dt$category[best_idx])
    ))
  }

  list(
    matched = FALSE,
    rule_name = sprintf("%+.4f Da", delta_mz),
    delta_formula = "",
    delta_ppm = NA_real_,
    category = "Unclassified"
  )
}

#' @keywords internal
#' @noRd
.nts_extract_murcko_scaffold_batch <- function(smiles_vec, python_path = NULL) {
  smiles_vec <- as.character(smiles_vec)
  u_smiles <- unique(smiles_vec[!is.na(smiles_vec) & nzchar(trimws(smiles_vec))])
  
  empty_res <- data.table::data.table(
    smiles = character(0L),
    scaffold_smiles = character(0L),
    generic_scaffold = character(0L),
    is_cyclic = logical(0L),
    num_rings = integer(0L)
  )
  if (!length(u_smiles)) {
    return(empty_res)
  }

  # 1. Try Python RDKit batch execution via joint_fp_engine.py --mode scaffold
  py_exe <- .nts_resolve_python_executable(python_path)
  engine_script <- .nts_resolve_joint_fp_engine_script()

  has_py <- !is.null(py_exe) && is.character(py_exe) && nzchar(py_exe) && file.exists(py_exe)
  has_eng <- !is.null(engine_script) && is.character(engine_script) && nzchar(engine_script) && file.exists(engine_script)

  if (isTRUE(has_py) && isTRUE(has_eng)) {
    tmp_in <- tempfile(pattern = "scaffold_in_", fileext = ".json")
    tmp_out <- tempfile(pattern = "scaffold_out_", fileext = ".json")
    on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)

    input_list <- lapply(seq_along(u_smiles), function(i) {
      list(id = sprintf("SMI_%d", i), smiles = u_smiles[i])
    })
    jsonlite::write_json(input_list, tmp_in, auto_unbox = TRUE)

    cmd_args <- c(
      shQuote(engine_script),
      "--mode", "scaffold",
      "--input_json", shQuote(tmp_in),
      "--output_json", shQuote(tmp_out)
    )
    run_res <- suppressWarnings(system2(py_exe, args = cmd_args, stdout = TRUE, stderr = TRUE))
    if (file.exists(tmp_out) && file.info(tmp_out)$size > 5L) {
      py_data <- tryCatch(jsonlite::fromJSON(tmp_out), error = function(e) NULL)
      if (is.data.frame(py_data) && "scaffold_smiles" %in% names(py_data)) {
        return(data.table::data.table(
          smiles = as.character(py_data$smiles),
          scaffold_smiles = as.character(py_data$scaffold_smiles),
          generic_scaffold = as.character(py_data$generic_scaffold),
          is_cyclic = as.logical(py_data$is_cyclic),
          num_rings = as.integer(py_data$num_rings)
        ))
      }
    }
  }

  # 2. Fallback: Pure R regex ring detection (graceful degradation)
  has_digits <- grepl("[0-9]", u_smiles)
  has_aromatic <- grepl("[cbnops]", u_smiles)
  is_cyc <- has_digits | has_aromatic

  data.table::data.table(
    smiles = u_smiles,
    scaffold_smiles = rep(NA_character_, length(u_smiles)),
    generic_scaffold = rep(NA_character_, length(u_smiles)),
    is_cyclic = is_cyc,
    num_rings = ifelse(is_cyc, 1L, 0L)
  )
}

#' Infer Chemical Scaffold and Analogs for Unannotated Features from Semantic Neighbors
#'
#' @description
#' Implements the DreaMS / Spec2Vec analog and scaffold retrieval paradigm for non-targeted
#' screening. For unannotated features (or features with low confidence levels, such as Level 4/5),
#' this function retrieves their nearest semantic neighbors (cosine similarity of 256-d MS2
#' embeddings \eqn{\ge} \code{min_semantic_score}) that possess verified or putative chemical annotations.
#'
#' From these reference analogs, it extracts:
#' \enumerate{
#'   \item The top reference analogs (name, SMILES, formula, similarity score).
#'   \item Putative mass shifts and transformations (\eqn{\Delta m/z}, matched against a dictionary of common environmental chemical transformations).
#'   \item Bemis-Murcko scaffolds (core cyclic/aromatic frameworks extracted via RDKit).
#'   \item ClassyFire / chemical taxonomy classes and subclasses (if present).
#'   \item A human-readable synthesis summary and inferred analog classification.
#' }
#'
#' @param project An \code{NTSProject} object, or a list/data.frame of query features.
#' @param network_data Optional molecular network list with \code{nodes} and \code{edges}
#'   (e.g., from [get_network_data()]). If \code{NULL}, semantic similarities are computed directly
#'   from 256-d MS2 embeddings via BLAS matrix multiplication.
#' @param reference_db Optional reference database for searching analogs beyond the current
#'   project (e.g., \code{"inhouse"}, \code{"online"}, an \code{NTSDatabase} object, or path to an RDS file).
#'   Defaults to \code{NULL} (within-project search).
#' @param target_features Optional character vector of feature names/IDs to analyze. If \code{NULL} (default),
#'   automatically selects unannotated features (features without verified L1/L2a/L2b annotations).
#' @param anchor_features Optional character vector of feature names/IDs to use as reference seeds.
#'   Defaults to features with confirmed SMILES and names in \code{@db_annotations}.
#' @param anchors Optional \code{data.frame} of custom anchor compounds with columns
#'   \code{name} / \code{anchor_id}, \code{SMILES}, and optional \code{mz}, \code{Formula}, \code{Class}.
#' @param min_semantic_score Minimum cosine similarity threshold for semantic analogs (default 0.70).
#' @param top_k Maximum number of top semantic analogs to retain per feature (default 3L).
#' @param transformation_rules Optional table or list of mass difference transformation rules.
#'   Defaults to [default_transformation_rules()].
#' @param mz_tol_ppm Mass tolerance in ppm for matching transformation rules (default 15.0 ppm).
#' @param embeddings Optional precomputed 256-dimensional MS2 embeddings matrix. If \code{NULL},
#'   retrieved from \code{project@physchem_proxy$ms2_embeddings} or computed on-the-fly.
#' @param python_path Optional custom path to Python executable for RDKit-based Murcko scaffold extraction.
#' @param write_to_project Logical; if \code{TRUE} (default for \code{NTSProject}), writes the
#'   resulting inference tables into \code{project@physchem_proxy} and logs the operation.
#' @param verbose Logical; print progress messages (default \code{TRUE}).
#'
#' @return If \code{project} is an \code{NTSProject} and \code{write_to_project = TRUE}, returns the
#'   modified \code{NTSProject}. Otherwise, returns a list containing:
#'   \itemize{
#'     \item \code{summary}: A [data.table::data.table] with one row per target feature containing
#'           top analog match, \eqn{\Delta m/z}, transformation rule, Murcko scaffold, and concise summary text.
#'     \item \code{analogs}: A [data.table::data.table] containing the detailed long list of all Top-K analog hits.
#'   }
#' @export
infer_scaffold_from_semantic_neighbors <- function(project,
                                                   network_data = NULL,
                                                   reference_db = NULL,
                                                   target_features = NULL,
                                                   anchor_features = NULL,
                                                   anchors = NULL,
                                                   min_semantic_score = 0.70,
                                                   top_k = 3L,
                                                   transformation_rules = NULL,
                                                   mz_tol_ppm = 15.0,
                                                   embeddings = NULL,
                                                   python_path = NULL,
                                                   write_to_project = TRUE,
                                                   verbose = TRUE) {
  is_proj <- methods::is(project, "NTSProject")

  # 1. Resolve feature table and metadata
  ft_df <- if (is_proj) {
    ft_raw <- get_feature_table(project)
    if (is.null(ft_raw) || !nrow(ft_raw)) {
      if (isTRUE(verbose)) message("Feature table is empty; cannot infer scaffolds.")
      return(if (is_proj && write_to_project) project else list(summary = data.table::data.table(), analogs = data.table::data.table()))
    }
    data.frame(
      name = as.character(ft_raw$name),
      mz = as.numeric(ft_raw$mz),
      rt = as.numeric(ft_raw$rt),
      stringsAsFactors = FALSE
    )
  } else if (is.data.frame(project)) {
    as.data.frame(project)
  } else if (is.list(project) && "name" %in% names(project)) {
    as.data.frame(project)
  } else {
    stop("`project` must be an NTSProject, data.frame, or list containing feature metadata.", call. = FALSE)
  }

  name_col <- if ("name" %in% names(ft_df)) "name" else if ("feature_id" %in% names(ft_df)) "feature_id" else names(ft_df)[1L]
  ft_df$feature_id <- as.character(ft_df[[name_col]])
  ft_df$mz <- as.numeric(ft_df$mz %||% NA_real_)
  ft_df$rt <- as.numeric(ft_df$rt %||% NA_real_)

  # 2. Resolve transformation rules
  rules_dt <- if (!is.null(transformation_rules)) {
    data.table::as.data.table(transformation_rules)
  } else {
    .nts_default_environmental_transformations()
  }

  # 3. Resolve Embeddings matrix
  emb_mat <- if (!is.null(embeddings)) {
    .nts_as_embedding_matrix(embeddings)
  } else if (is_proj) {
    proxy <- get_physchem_proxy(project)
    m <- if (is.data.frame(proxy) && nrow(proxy)) .nts_as_embedding_matrix(proxy) else NULL
    if (is.null(m)) {
      m <- tryCatch(
        extract_ms2_embeddings(project, python_path = python_path, verbose = verbose),
        error = function(e) NULL
      )
    }
    m
  } else {
    NULL
  }

  if (is.null(emb_mat) && !is.null(embeddings) && is.matrix(embeddings)) {
    emb_mat <- embeddings
  }

  if (is.null(emb_mat) || !is.matrix(emb_mat) || nrow(emb_mat) < 1L) {
    warning("No valid MS2 embeddings available; cannot perform semantic analog retrieval.", call. = FALSE)
    empty_out <- list(summary = data.table::data.table(), analogs = data.table::data.table())
    return(if (is_proj && write_to_project) project else empty_out)
  }

  # Auto-derive row names if omitted (e.g. anonymous embeddings matrix)
  if (is.null(rownames(emb_mat))) {
    if (nrow(emb_mat) == nrow(ft_df)) {
      rownames(emb_mat) <- as.character(ft_df$feature_id)
    } else {
      rownames(emb_mat) <- paste0("FT_", seq_len(nrow(emb_mat)))
    }
  }

  # Defensively prune non-finite rows (NA/NaN/Inf)
  finite_rows <- rowSums(is.finite(emb_mat)) == ncol(emb_mat)
  if (!all(finite_rows)) {
    emb_mat <- emb_mat[finite_rows, , drop = FALSE]
  }
  if (nrow(emb_mat) < 1L) {
    warning("No valid finite MS2 embeddings rows remaining after filtering non-finite values.", call. = FALSE)
    empty_out <- list(summary = data.table::data.table(), analogs = data.table::data.table())
    return(if (is_proj && write_to_project) project else empty_out)
  }

  # Ensure embeddings are L2 normalized
  r_norms <- sqrt(rowSums(emb_mat^2))
  r_norms[r_norms < 1e-12] <- 1.0
  emb_mat <- emb_mat / r_norms

  # 4. Resolve Anchors (Confirmed reference seeds with SMILES)
  anchors_list <- list()

  # 4a. From NTSProject db_annotations
  if (is_proj) {
    db_raw <- as.data.frame(project@db_annotations)
    if (nrow(db_raw) > 0L && "SMILES" %in% names(db_raw)) {
      valid_smi_mask <- !is.na(db_raw$SMILES) & nzchar(trimws(as.character(db_raw$SMILES)))
      if (any(valid_smi_mask)) {
        db_valid <- db_raw[valid_smi_mask, , drop = FALSE]
        if ("score_total" %in% names(db_valid)) {
          db_valid <- db_valid[order(db_valid$score_total, decreasing = TRUE), ]
        }
        db_dedup <- db_valid[!duplicated(db_valid$name), ]
        
        mz_map <- stats::setNames(ft_df$mz, ft_df$feature_id)
        
        anchors_list[[length(anchors_list) + 1L]] <- data.table::data.table(
          anchor_id = as.character(db_dedup$name),
          anchor_name = as.character(db_dedup$Compound.name %||% db_dedup$name),
          anchor_smiles = as.character(db_dedup$SMILES),
          anchor_formula = as.character(db_dedup$Formula %||% NA_character_),
          anchor_mz = as.numeric(mz_map[as.character(db_dedup$name)] %||% db_dedup$exp_mz %||% NA_real_),
          chemical_class = as.character(db_dedup$First_Class %||% db_dedup$Class %||% NA_character_),
          source = "project_anchor"
        )
      }
    }
  }

  # 4b. From custom anchors data.frame
  if (!is.null(anchors) && is.data.frame(anchors) && nrow(anchors)) {
    a_id_col <- if ("anchor_id" %in% names(anchors)) "anchor_id" else if ("name" %in% names(anchors)) "name" else names(anchors)[1L]
    a_name_col <- if ("Compound.name" %in% names(anchors)) "Compound.name" else if ("anchor_name" %in% names(anchors)) "anchor_name" else a_id_col
    a_smi_col <- if ("SMILES" %in% names(anchors)) "SMILES" else if ("anchor_smiles" %in% names(anchors)) "anchor_smiles" else NA_character_
    a_mz_col <- if ("exp_mz" %in% names(anchors)) "exp_mz" else if ("mz" %in% names(anchors)) "mz" else NA_character_
    
    anchors_list[[length(anchors_list) + 1L]] <- data.table::data.table(
      anchor_id = as.character(anchors[[a_id_col]]),
      anchor_name = as.character(anchors[[a_name_col]]),
      anchor_smiles = if (!is.na(a_smi_col)) as.character(anchors[[a_smi_col]]) else rep(NA_character_, nrow(anchors)),
      anchor_formula = as.character(anchors$Formula %||% NA_character_),
      anchor_mz = if (!is.na(a_mz_col)) as.numeric(anchors[[a_mz_col]]) else rep(NA_real_, nrow(anchors)),
      chemical_class = as.character(anchors$First_Class %||% anchors$Class %||% NA_character_),
      source = "custom_anchor"
    )
  }

  # 4c. From external reference_db
  if (!is.null(reference_db)) {
    ref_obj <- if (is.character(reference_db) && length(reference_db) == 1L && reference_db %in% c("inhouse", "online", "prediction")) {
      load_builtin_pfas_library(reference_db)
    } else if (is.character(reference_db) && file.exists(reference_db)) {
      readRDS(reference_db)
    } else {
      reference_db
    }

    if (methods::is(ref_obj, "NTSDatabase")) {
      meta <- as.data.frame(ref_obj@metadata)
      if (nrow(meta) > 0L && "SMILES" %in% names(meta)) {
        smi_ok <- !is.na(meta$SMILES) & nzchar(trimws(as.character(meta$SMILES)))
        if (any(smi_ok)) {
          meta_val <- meta[smi_ok, ]
          anchors_list[[length(anchors_list) + 1L]] <- data.table::data.table(
            anchor_id = as.character(meta_val$Lab_ID %||% meta_val$ID %||% meta_val$name),
            anchor_name = as.character(meta_val$name %||% meta_val$Compound.name %||% meta_val$Lab_ID),
            anchor_smiles = as.character(meta_val$SMILES),
            anchor_formula = as.character(meta_val$Formula %||% NA_character_),
            anchor_mz = as.numeric(meta_val$target_mz %||% meta_val$mass %||% NA_real_),
            chemical_class = as.character(meta_val$First_Class %||% meta_val$Class %||% NA_character_),
            source = "reference_library"
          )
        }
      }
    }
  }

  anchors_dt <- if (length(anchors_list)) {
    data.table::rbindlist(anchors_list, fill = TRUE, use.names = TRUE)
  } else {
    data.table::data.table()
  }

  # 4d. Fallback if project is data.frame and anchor_features is passed
  if (!nrow(anchors_dt) && !is.null(anchor_features)) {
    af_ids <- intersect(as.character(anchor_features), rownames(emb_mat))
    if (length(af_ids)) {
      m_idx <- match(af_ids, ft_df$feature_id)
      af_rows <- ft_df[m_idx, , drop = FALSE]
      anchors_dt <- data.table::data.table(
        anchor_id = af_ids,
        anchor_name = if (nrow(af_rows) && "Compound.name" %in% names(af_rows)) af_rows$Compound.name else af_ids,
        anchor_smiles = if (nrow(af_rows) && "SMILES" %in% names(af_rows)) af_rows$SMILES else rep(NA_character_, length(af_ids)),
        anchor_formula = if (nrow(af_rows) && "Formula" %in% names(af_rows)) af_rows$Formula else rep(NA_character_, length(af_ids)),
        anchor_mz = if (nrow(af_rows) && "mz" %in% names(af_rows)) af_rows$mz else rep(NA_real_, length(af_ids)),
        chemical_class = if (nrow(af_rows) && "First_Class" %in% names(af_rows)) af_rows$First_Class else if (nrow(af_rows) && "Class" %in% names(af_rows)) af_rows$Class else rep(NA_character_, length(af_ids)),
        source = "custom_anchor"
      )
    }
  }

  if (nrow(anchors_dt)) {
    anchors_dt <- unique(anchors_dt, by = "anchor_id")
  }

  # Filter anchor_features if explicitly specified
  if (!is.null(anchor_features) && length(anchor_features) && nrow(anchors_dt)) {
    anchors_dt <- anchors_dt[anchor_id %in% as.character(anchor_features)]
  }

  # Only retain anchors that have corresponding embedding vectors
  common_anchor_ids <- intersect(anchors_dt$anchor_id, rownames(emb_mat))
  anchors_dt <- anchors_dt[anchor_id %in% common_anchor_ids]

  if (!nrow(anchors_dt)) {
    if (isTRUE(verbose)) message("No anchor features with valid MS2 embeddings found.")
    empty_out <- list(summary = data.table::data.table(), analogs = data.table::data.table())
    return(if (is_proj && write_to_project) project else empty_out)
  }

  # 5. Determine Target Features (Unannotated features)
  if (is.null(target_features)) {
    target_ids <- setdiff(ft_df$feature_id, anchors_dt$anchor_id)
  } else {
    target_ids <- intersect(as.character(target_features), ft_df$feature_id)
  }

  # Only keep targets with MS2 embeddings
  target_ids_with_ms2 <- intersect(target_ids, rownames(emb_mat))
  target_ids_no_ms2 <- setdiff(target_ids, rownames(emb_mat))

  if (isTRUE(verbose)) {
    message(sprintf("Running Scaffold & Analog Inference on %d target features (%d with MS2, %d MS1-only) against %d anchors (min_score=%.2f, top_k=%d)...",
                    length(target_ids), length(target_ids_with_ms2), length(target_ids_no_ms2), nrow(anchors_dt), min_semantic_score, top_k))
  }

  # 6. Extract Murcko Scaffolds for all unique anchor SMILES
  valid_anchor_smiles <- anchors_dt$anchor_smiles[!is.na(anchors_dt$anchor_smiles) & nzchar(trimws(anchors_dt$anchor_smiles))]
  scaffolds_map_dt <- if (length(valid_anchor_smiles)) {
    .nts_extract_murcko_scaffold_batch(valid_anchor_smiles, python_path = python_path)
  } else {
    data.table::data.table(smiles = character(0L), scaffold_smiles = character(0L), generic_scaffold = character(0L), is_cyclic = logical(0L), num_rings = integer(0L))
  }

  if (nrow(scaffolds_map_dt)) {
    anchors_dt <- scaffolds_map_dt[anchors_dt, on = "smiles == anchor_smiles"]
    data.table::setnames(anchors_dt, "smiles", "anchor_smiles")
  } else {
    anchors_dt[, `:=`(scaffold_smiles = NA_character_, generic_scaffold = NA_character_, is_cyclic = NA, num_rings = 0L)]
  }

  # 7. Compute Semantic Similarities via Matrix Dot Product
  min_s <- suppressWarnings(as.numeric(min_semantic_score)[1L])
  if (!is.finite(min_s) || min_s < 0) min_s <- 0.70
  tk <- if (is.finite(top_k)) max(as.integer(top_k), 1L) else 3L

  mz_lookup <- stats::setNames(ft_df$mz, ft_df$feature_id)
  rt_lookup <- stats::setNames(ft_df$rt, ft_df$feature_id)

  analogs_hit_list <- vector("list", length(target_ids_with_ms2))
  summary_row_list <- vector("list", length(target_ids))

  if (length(target_ids_with_ms2)) {
    target_emb <- emb_mat[target_ids_with_ms2, , drop = FALSE]
    anchor_emb <- emb_mat[anchors_dt$anchor_id, , drop = FALSE]

    # BLAS matrix multiplication: target x anchor
    sim_matrix <- target_emb %*% t(anchor_emb)

    for (i in seq_along(target_ids_with_ms2)) {
      tid <- target_ids_with_ms2[i]
      tmz <- mz_lookup[[tid]] %||% NA_real_
      trt <- rt_lookup[[tid]] %||% NA_real_

      sims <- sim_matrix[i, ]
      sims <- pmin(pmax(sims, -1.0), 1.0)
      hit_idx <- which(sims >= min_s & anchors_dt$anchor_id != tid)

      if (length(hit_idx)) {
        ord <- order(sims[hit_idx], decreasing = TRUE)
        sel_hits <- hit_idx[ord[seq_len(min(length(ord), tk))]]

        cur_analogs <- vector("list", length(sel_hits))
        for (h in seq_along(sel_hits)) {
          aid_idx <- sel_hits[h]
          anc_row <- anchors_dt[aid_idx]
          sim_val <- round(as.numeric(sims[aid_idx]), 4L)
          dmz <- if (is.finite(tmz) && is.finite(anc_row$anchor_mz)) tmz - anc_row$anchor_mz else NA_real_

          trans_info <- .nts_match_transformation(dmz, rules_dt = rules_dt, mz_tol_ppm = mz_tol_ppm)

          cur_analogs[[h]] <- data.table::data.table(
            feature_id = tid,
            rank = as.integer(h),
            analog_id = anc_row$anchor_id,
            analog_name = anc_row$anchor_name,
            analog_smiles = anc_row$anchor_smiles,
            analog_formula = anc_row$anchor_formula,
            analog_mz = anc_row$anchor_mz,
            similarity = sim_val,
            delta_mz = if (is.finite(dmz)) round(dmz, 4L) else NA_real_,
            delta_ppm = trans_info$delta_ppm,
            transformation = trans_info$rule_name,
            transformation_formula = trans_info$delta_formula,
            transformation_category = trans_info$category,
            murcko_scaffold = anc_row$scaffold_smiles,
            generic_scaffold = anc_row$generic_scaffold,
            is_cyclic = anc_row$is_cyclic,
            chemical_class = anc_row$chemical_class,
            source = anc_row$source
          )
        }
        analogs_hit_list[[i]] <- data.table::rbindlist(cur_analogs, fill = TRUE)

        # Summary for feature with hits
        top_h <- cur_analogs[[1L]]
        
        # Format human-readable summary text
        scaff_desc <- if (isTRUE(top_h$is_cyclic) && !is.na(top_h$murcko_scaffold) && nzchar(top_h$murcko_scaffold)) {
          sprintf("Scaffold: %s; ", top_h$murcko_scaffold)
        } else if (identical(top_h$is_cyclic, FALSE)) {
          "Acyclic; "
        } else {
          ""
        }
        class_desc <- if (!is.na(top_h$chemical_class) && nzchar(top_h$chemical_class)) {
          sprintf("Class: %s", top_h$chemical_class)
        } else {
          ""
        }
        shift_desc <- if (!is.na(top_h$transformation) && nzchar(top_h$transformation)) {
          sprintf("shift: %s (Δmz=%+.4f%s)",
                  top_h$transformation,
                  top_h$delta_mz,
                  if (!is.na(top_h$delta_ppm)) sprintf(", %.1f ppm", top_h$delta_ppm) else "")
        } else {
          "shift: Unknown"
        }

        summary_text <- sprintf("Analog to %s (sim=%.3f); %s; %s%s",
                                top_h$analog_name,
                                top_h$similarity,
                                shift_desc,
                                scaff_desc,
                                class_desc)

        summary_row_list[[i]] <- data.table::data.table(
          feature_id = tid,
          precursor_mz = tmz,
          rt = trt,
          analog_count = as.integer(length(sel_hits)),
          top_analog_id = top_h$analog_id,
          top_analog_name = top_h$analog_name,
          top_analog_smiles = top_h$analog_smiles,
          top_analog_similarity = top_h$similarity,
          delta_mz = top_h$delta_mz,
          inferred_transformation = top_h$transformation,
          murcko_scaffold = top_h$murcko_scaffold,
          is_cyclic = top_h$is_cyclic,
          chemical_class = top_h$chemical_class,
          inferred_summary = trimws(summary_text)
        )
      } else {
        # Feature with MS2 but no hits >= min_s
        summary_row_list[[i]] <- data.table::data.table(
          feature_id = tid,
          precursor_mz = tmz,
          rt = trt,
          analog_count = 0L,
          top_analog_id = NA_character_,
          top_analog_name = NA_character_,
          top_analog_smiles = NA_character_,
          top_analog_similarity = NA_real_,
          delta_mz = NA_real_,
          inferred_transformation = NA_character_,
          murcko_scaffold = NA_character_,
          is_cyclic = NA,
          chemical_class = NA_character_,
          inferred_summary = sprintf("No semantic analogs found above threshold (min_score=%.2f)", min_s)
        )
      }
    }
  }

  # Append MS1-only features if any
  if (length(target_ids_no_ms2)) {
    offset <- length(target_ids_with_ms2)
    for (j in seq_along(target_ids_no_ms2)) {
      tid <- target_ids_no_ms2[j]
      tmz <- mz_lookup[[tid]] %||% NA_real_
      trt <- rt_lookup[[tid]] %||% NA_real_
      summary_row_list[[offset + j]] <- data.table::data.table(
        feature_id = tid,
        precursor_mz = tmz,
        rt = trt,
        analog_count = 0L,
        top_analog_id = NA_character_,
        top_analog_name = NA_character_,
        top_analog_smiles = NA_character_,
        top_analog_similarity = NA_real_,
        delta_mz = NA_real_,
        inferred_transformation = NA_character_,
        murcko_scaffold = NA_character_,
        is_cyclic = NA,
        chemical_class = NA_character_,
        inferred_summary = "No MS2 spectrum available for semantic embedding"
      )
    }
  }

  summary_dt <- data.table::rbindlist(summary_row_list[!vapply(summary_row_list, is.null, logical(1L))], fill = TRUE)
  analogs_dt <- data.table::rbindlist(analogs_hit_list[!vapply(analogs_hit_list, is.null, logical(1L))], fill = TRUE)

  if (isTRUE(verbose)) {
    n_mapped <- sum(summary_dt$analog_count > 0L)
    message(sprintf("Inference complete: %d / %d target features successfully mapped to semantic analogs (%.1f%%).",
                    n_mapped, nrow(summary_dt), (n_mapped / max(1L, nrow(summary_dt))) * 100))
  }

  # 8. Project Write-back
  if (is_proj && isTRUE(write_to_project)) {
    proxy <- as.data.frame(get_physchem_proxy(project), stringsAsFactors = FALSE)
    if (!nrow(proxy) || !"name" %in% names(proxy)) {
      ft <- get_feature_table(project)
      proxy <- data.frame(name = as.character(ft$name), stringsAsFactors = FALSE)
    }

    # Match summary_dt into proxy columns
    m <- match(as.character(proxy$name), as.character(summary_dt$feature_id))
    proxy$Inferred_Analog_Name <- summary_dt$top_analog_name[m]
    proxy$Inferred_Analog_SMILES <- summary_dt$top_analog_smiles[m]
    proxy$Inferred_Analog_Sim <- summary_dt$top_analog_similarity[m]
    proxy$Inferred_Transformation <- summary_dt$inferred_transformation[m]
    proxy$Inferred_Murcko_Scaffold <- summary_dt$murcko_scaffold[m]
    proxy$Inferred_Scaffold_Summary <- summary_dt$inferred_summary[m]

    # Store full tables as attribute on physchem_proxy
    attr(proxy, "analog_scaffold_inference") <- list(
      summary = summary_dt,
      analogs = analogs_dt
    )

    project <- set_physchem_proxy(project, proxy)
    
    # Update log
    log_msg <- sprintf("infer_scaffold_from_semantic_neighbors: evaluated %d features, %d mapped to semantic analogs (min_score=%.2f, top_k=%d)",
                       nrow(summary_dt), sum(summary_dt$analog_count > 0L), min_s, tk)
    if (.hasSlot(project, "log")) {
      project@log <- .nts_log_append(
        project@log,
        "infer_scaffold_from_semantic_neighbors",
        list(
          message = log_msg,
          n_evaluated = nrow(summary_dt),
          n_mapped = sum(summary_dt$analog_count > 0L),
          min_score = min_s,
          top_k = tk
        )
      )
    }

    return(project)
  }

  list(summary = summary_dt, analogs = analogs_dt)
}

#' Retrieve Inferred Scaffolds and Analogs from NTSProject
#'
#' Extracts the scaffold and analog inference tables stored in an [NTSProject].
#'
#' @param project An [NTSProject] object.
#' @param type One of \code{"summary"} (one row per target feature) or \code{"analogs"}
#'   (long table of all top-K hits).
#'
#' @return A [data.table::data.table] of inferred scaffolds or analog matches.
#' @export
get_inferred_scaffolds <- function(project, type = c("summary", "analogs")) {
  if (!methods::is(project, "NTSProject")) {
    stop("`project` must be an NTSProject object.", call. = FALSE)
  }
  type <- match.arg(type)
  proxy <- get_physchem_proxy(project)
  inf <- attr(proxy, "analog_scaffold_inference")
  if (is.null(inf) || !is.list(inf)) {
    if (type == "summary" && is.data.frame(proxy) && "Inferred_Analog_Name" %in% names(proxy)) {
      return(data.table::as.data.table(proxy))
    }
    return(data.table::data.table())
  }
  if (type == "summary") {
    inf$summary %||% data.table::data.table()
  } else {
    inf$analogs %||% data.table::data.table()
  }
}
