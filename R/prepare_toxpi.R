# ToxPi input preparation from EPI Suite + ADMET-AI predictions
# Ref: tests/toxPi_202506065.R

# ---- constants ----------------------------------------------------------------

.nts_toxpi_scenario_weights <- function(scenario, slice_names) {
  scenario <- match.arg(scenario, c("environmental", "drug_safety", "pesticide", "equal"))
  base <- list(
    environmental = c(
      Persistence = 0.25, BioDegradation = 0.15, AquaticTox = 0.20,
      SystemicTox = 0.20, Genotox = 0.10, Tox21 = 0.10
    ),
    drug_safety = c(
      Persistence = 0.05, BioDegradation = 0.05, AquaticTox = 0.05,
      SystemicTox = 0.35, Genotox = 0.20, Tox21 = 0.15
    ),
    pesticide = c(
      Persistence = 0.20, BioDegradation = 0.20, AquaticTox = 0.25,
      SystemicTox = 0.15, Genotox = 0.10, Tox21 = 0.10
    )
  )
  if (identical(scenario, "equal")) {
    w <- rep(1 / length(slice_names), length(slice_names))
    names(w) <- slice_names
    return(w)
  }
  w <- base[[scenario]]
  missing <- setdiff(slice_names, names(w))
  if (length(missing)) {
    extra <- rep(0, length(missing))
    names(extra) <- missing
    w <- c(w, extra)
  }
  w <- w[slice_names]
  if ("Tox21_NR" %in% slice_names && "Tox21_SR" %in% slice_names && "Tox21" %in% names(base[[scenario]])) {
    tox21_w <- base[[scenario]]["Tox21"]
    w["Tox21_NR"] <- tox21_w / 2
    w["Tox21_SR"] <- tox21_w / 2
    w <- w[names(w) != "Tox21"]
  }
  w / sum(w)
}

.nts_toxpi_tox21_nr_cols <- function() {
  c(
    "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase",
    "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma"
  )
}

.nts_toxpi_tox21_sr_cols <- function() {
  c("SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53")
}

.nts_toxpi_required_epi <- function() {
  c(
    "LogKow_est", "LogKoc_est", "LogKoa_est", "LogBCF", "LogBAF",
    paste0("Biowin", seq_len(7L)),
    "Fish_LC50", "Daphnid_LC50", "Algae_LC50",
    "Fish_chV", "Daphnid_chV", "Algae_chV"
  )
}

.nts_toxpi_required_admet <- function() {
  c(
    "LD50_Zhu", "DILI", "ClinTox", "hERG",
    "AMES", "Carcinogens_Lagunin"
  )
}

# ---- internal helpers ---------------------------------------------------------

#' @keywords internal
#' @noRd
.nts_toxpi_msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    message(...)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.nts_toxpi_warn <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    warning(..., call. = FALSE)
  }
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.nts_toxpi_resolve_col <- function(df, candidates) {
  nms <- names(df)
  for (cand in candidates) {
    if (cand %in% nms) {
      return(cand)
    }
    alt <- make.names(cand)
    if (alt %in% nms) {
      return(alt)
    }
  }
  NA_character_
}

#' @keywords internal
#' @noRd
.nts_toxpi_inchikey_from_mol <- function(mol) {
  if (is.null(mol)) {
    return(NA_character_)
  }
  if (!requireNamespace("rJava", quietly = TRUE)) {
    return(NA_character_)
  }
  ac <- if (inherits(mol, "jobjRef")) {
    mol
  } else if (is.object(mol) && !is.null(mol@jobj)) {
    mol@jobj
  } else {
    return(NA_character_)
  }
  tryCatch({
    factory <- rJava::.jcall(
      "org/openscience/cdk/inchi/InChIGeneratorFactory",
      "Lorg/openscience/cdk/inchi/InChIGeneratorFactory;",
      "getInstance"
    )
    gen <- rJava::.jcall(
      factory,
      "Lorg/openscience/cdk/inchi/InChIGenerator;",
      "getInChIGenerator",
      ac
    )
    key <- rJava::.jcall(gen, "Ljava/lang/String;", "getInchiKey")
    as.character(key)[1L]
  }, error = function(e) {
    # Legacy CDK (< 2.8): getInChIKey(); newer CDK uses getInchiKey() above.
    tryCatch({
      gen <- rJava::J("org.openscience.cdk.inchi.InChIGeneratorFactory")$getInstance()$getInChIGenerator(ac)
      as.character(gen$getInChIKey())[1L]
    }, error = function(e2) NA_character_)
  })
}

#' @keywords internal
#' @noRd
.nts_toxpi_canonical_smiles_one <- function(smi) {
  smi <- trimws(as.character(smi)[1L])
  if (is.na(smi) || !nzchar(smi)) {
    return(NA_character_)
  }
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    return(smi)
  }
  mol <- .nts_mol_first_from_smiles_rcdk(smi)
  if (is.null(mol)) {
    return(smi)
  }
  out <- tryCatch({
    if (requireNamespace("rJava", quietly = TRUE)) {
      sg <- rJava::J("org/openscience/cdk/smiles/SmilesGenerator")$getInstance()
      cs <- sg$create(mol)
      as.character(cs)[1L]
    } else {
      NA_character_
    }
  }, error = function(e) NA_character_)
  if (!is.na(out) && nzchar(out)) {
    return(out)
  }
  if ("get.smiles" %in% getNamespaceExports("rcdk")) {
    out2 <- tryCatch(as.character(rcdk::get.smiles(mol)[1L]), error = function(e) NA_character_)
    if (!is.na(out2) && nzchar(out2)) {
      return(out2)
    }
  }
  smi
}

#' Canonicalize SMILES for cross-engine alignment
#'
#' Normalizes SMILES strings before [predict_episuite()] and
#' [predict_admet_risk()] so EPI / ADMET tables can be merged reliably in
#' [prepare_toxpi()]. Uses **rcdk** (CDK `SmilesGenerator`) by default; optional
#' RDKit Python via **reticulate** when `method = "rdkit_python"`.
#'
#' @param smiles_list Character vector of SMILES.
#' @param method `"auto"` (rcdk then trimmed input), `"rcdk"`, or `"rdkit_python"`.
#'
#' @return Character vector of canonical SMILES (same length as input; failed parse -> `NA`).
#' @export
#'
#' @seealso [prepare_toxpi]
canonicalize_smiles <- function(smiles_list,
                                method = c("auto", "rcdk", "rdkit_python")) {
  method <- match.arg(method)
  smi <- trimws(as.character(smiles_list))
  if (!length(smi)) {
    return(character())
  }

  canon <- vapply(smi, function(s) {
    if (is.na(s) || !nzchar(s)) {
      return(NA_character_)
    }
    if (method %in% c("auto", "rcdk")) {
      out <- .nts_toxpi_canonical_smiles_one(s)
      if (!is.na(out) && nzchar(out)) {
        return(out)
      }
    }
    if (method %in% c("auto", "rdkit_python") && requireNamespace("reticulate", quietly = TRUE)) {
      out <- tryCatch({
        rdkit <- reticulate::import("rdkit.Chem", delay_load = TRUE)
        mol <- rdkit$MolFromSmiles(s)
        if (is.null(mol)) {
          return(NA_character_)
        }
        as.character(rdkit$MolToSmiles(mol))
      }, error = function(e) NA_character_)
      if (!is.na(out) && nzchar(out)) {
        return(out)
      }
    }
    if (method == "rdkit_python") {
      return(NA_character_)
    }
    s
  }, character(1L), USE.NAMES = FALSE)

  stats::setNames(canon, names(smiles_list))
}

#' @keywords internal
#' @noRd
.nts_toxpi_inchikeys <- function(smiles) {
  smiles <- trimws(as.character(smiles))
  vapply(
    smiles,
    function(s) {
      if (is.na(s) || !nzchar(s)) {
        return(NA_character_)
      }
      mol <- if (requireNamespace("rcdk", quietly = TRUE)) {
        .nts_mol_first_from_smiles_rcdk(s)
      } else {
        NULL
      }
      .nts_toxpi_inchikey_from_mol(mol)
    },
    character(1),
    USE.NAMES = FALSE
  )
}

#' @keywords internal
#' @noRd
.nts_toxpi_structural_keys <- function(smiles) {
  smiles <- trimws(as.character(smiles))
  inchikeys <- .nts_toxpi_inchikeys(smiles)
  canon <- vapply(smiles, .nts_toxpi_canonical_smiles_one, character(1L), USE.NAMES = FALSE)
  keys <- inchikeys
  miss <- is.na(keys) | !nzchar(keys)
  keys[miss] <- canon[miss]
  miss2 <- is.na(keys) | !nzchar(keys)
  keys[miss2] <- smiles[miss2]
  keys
}

#' @keywords internal
#' @noRd
.nts_toxpi_pick_query_smiles <- function(epi, admet, query_smiles = NULL) {
  if (!is.null(query_smiles)) {
    return(as.character(query_smiles))
  }
  if ("query_SMILES" %in% names(epi)) {
    return(as.character(epi$query_SMILES))
  }
  if ("query_SMILES" %in% names(admet)) {
    return(as.character(admet$query_SMILES))
  }
  as.character(admet$SMILES)
}

#' @keywords internal
#' @noRd
.nts_toxpi_merge_by_order <- function(epi, admet, query_smiles = NULL, verbose = TRUE) {
  ne <- nrow(epi)
  na <- nrow(admet)
  if (ne != na) {
    stop(
      "Row-order alignment requires equal nrow (epi=", ne, ", admet=", na, "). ",
      "Use align_by = 'structure' or pass `query_smiles` of length ",
      max(ne, na), ".",
      call. = FALSE
    )
  }
  query_smiles <- .nts_toxpi_pick_query_smiles(epi, admet, query_smiles)
  if (length(query_smiles) != ne) {
    stop("`query_smiles` length (", length(query_smiles), ") must match nrow (", ne, ").", call. = FALSE)
  }

  sm_epi <- trimws(as.character(epi$SMILES))
  sm_admet <- trimws(as.character(admet$SMILES))
  mismatch <- !is.na(sm_epi) & !is.na(sm_admet) & sm_epi != sm_admet
  if (any(mismatch)) {
    .nts_toxpi_warn(
      sum(mismatch), " row(s): EPI returned different SMILES than ADMET; merged by row order. ",
      "Run canonicalize_smiles() before predict_*(), or pass `query_smiles`.",
      verbose = verbose
    )
  }

  admet_drop <- setdiff(names(admet), c("SMILES", "query_SMILES"))
  merged <- cbind(
    data.frame(
      query_SMILES = query_smiles,
      SMILES = query_smiles,
      SMILES_epi = sm_epi,
      SMILES_admet = sm_admet,
      stringsAsFactors = FALSE
    ),
    epi[, setdiff(names(epi), c("SMILES", "query_SMILES")), drop = FALSE],
    admet[, admet_drop, drop = FALSE]
  )
  merged
}

#' @keywords internal
#' @noRd
.nts_toxpi_merge_by_query_smiles <- function(epi, admet, verbose = TRUE) {
  if (!"query_SMILES" %in% names(epi) && !"query_SMILES" %in% names(admet)) {
    stop("align_by = 'query_smiles' requires query_SMILES column or `query_smiles` argument.", call. = FALSE)
  }
  if (!"query_SMILES" %in% names(epi)) {
    epi$query_SMILES <- admet$query_SMILES
  }
  if (!"query_SMILES" %in% names(admet)) {
    admet$query_SMILES <- epi$query_SMILES
  }
  epi$._toxpi_key <- trimws(as.character(epi$query_SMILES))
  admet$._toxpi_key <- trimws(as.character(admet$query_SMILES))
  merged <- merge(epi, admet, by = "._toxpi_key", all = FALSE, suffixes = c("_epi", "_admet"))
  if (nrow(merged) < max(nrow(epi), nrow(admet))) {
    .nts_toxpi_warn(
      "query_SMILES merge: ", max(nrow(epi), nrow(admet)) - nrow(merged),
      " compound(s) dropped (no EPI+ADMET pair).",
      verbose = verbose
    )
  }
  merged$query_SMILES <- merged$._toxpi_key
  merged$SMILES <- merged$query_SMILES
  merged
}

#' @keywords internal
#' @noRd
.nts_toxpi_merge_by_structure <- function(epi, admet, verbose = TRUE) {
  epi$._toxpi_key <- .nts_toxpi_structural_keys(epi$SMILES)
  admet$._toxpi_key <- .nts_toxpi_structural_keys(admet$SMILES)

  dup <- c(
    epi$._toxpi_key[duplicated(epi$._toxpi_key) & !is.na(epi$._toxpi_key)],
    admet$._toxpi_key[duplicated(admet$._toxpi_key) & !is.na(admet$._toxpi_key)]
  )
  dup <- unique(dup[nzchar(dup)])
  if (length(dup)) {
    stop(
      "Structural merge: duplicate join keys (", length(dup), "). ",
      "Use align_by = 'order' with equal-length tables or pass `query_smiles`.",
      call. = FALSE
    )
  }

  merged <- merge(epi, admet, by = "._toxpi_key", all = FALSE, suffixes = c("_epi", "_admet"))
  n_epi <- nrow(epi)
  n_admet <- nrow(admet)
  if (nrow(merged) < max(n_epi, n_admet)) {
    .nts_toxpi_warn(
      "Structural merge: ", max(n_epi, n_admet) - nrow(merged),
      " compound(s) lacked an EPI+ADMET pair.",
      verbose = verbose
    )
  }

  sm_epi <- if ("SMILES_epi" %in% names(merged)) merged$SMILES_epi else merged$SMILES
  sm_admet <- if ("SMILES_admet" %in% names(merged)) merged$SMILES_admet else merged$SMILES
  merged$query_SMILES <- if ("query_SMILES_epi" %in% names(merged)) {
    merged$query_SMILES_epi
  } else if ("query_SMILES_admet" %in% names(merged)) {
    merged$query_SMILES_admet
  } else {
    sm_admet
  }
  merged$SMILES <- ifelse(
    !is.na(sm_admet) & nzchar(as.character(sm_admet)),
    sm_admet,
    sm_epi
  )
  merged
}

#' @keywords internal
#' @noRd
.nts_toxpi_merge_inputs <- function(epi_data,
                                   admet_data,
                                   align_by = c("auto", "order", "structure", "query_smiles"),
                                   query_smiles = NULL,
                                   verbose = TRUE) {
  align_by <- match.arg(align_by)
  epi <- as.data.frame(epi_data, stringsAsFactors = FALSE)
  admet <- as.data.frame(admet_data, stringsAsFactors = FALSE)
  if (!"SMILES" %in% names(epi)) {
    stop("epi_data must contain a SMILES column.", call. = FALSE)
  }
  if (!"SMILES" %in% names(admet)) {
    stop("admet_data must contain a SMILES column.", call. = FALSE)
  }

  if (!is.null(query_smiles)) {
    qs <- as.character(query_smiles)
    if (length(qs) != nrow(epi) && length(qs) != nrow(admet) && length(qs) != max(nrow(epi), nrow(admet))) {
      stop(
        "`query_smiles` length (", length(qs), ") must match epi/admet nrow.",
        call. = FALSE
      )
    }
    if (length(qs) == nrow(epi)) {
      epi$query_SMILES <- qs
    }
    if (length(qs) == nrow(admet)) {
      admet$query_SMILES <- qs
    }
  }

  if (identical(align_by, "auto")) {
    align_by <- if ("query_SMILES" %in% names(epi) && "query_SMILES" %in% names(admet)) {
      "query_smiles"
    } else if (nrow(epi) == nrow(admet)) {
      "order"
    } else {
      "structure"
    }
    .nts_toxpi_msg("Merge strategy: align_by = '", align_by, "'.", verbose = verbose)
  }

  merged <- switch(
    align_by,
    order = .nts_toxpi_merge_by_order(epi, admet, query_smiles = query_smiles, verbose = verbose),
    query_smiles = .nts_toxpi_merge_by_query_smiles(epi, admet, verbose = verbose),
    structure = .nts_toxpi_merge_by_structure(epi, admet, verbose = verbose)
  )

  n_merged <- nrow(merged)
  if (n_merged > max(nrow(epi), nrow(admet))) {
    stop(
      "Merge produced ", n_merged, " rows (> input). ",
      "Use canonicalize_smiles() before predict or align_by = 'order'.",
      call. = FALSE
    )
  }

  admet_cols <- c(
    "LD50_Zhu", "DILI", "ClinTox", "hERG", "AMES", "Carcinogens_Lagunin",
    .nts_toxpi_tox21_nr_cols(), .nts_toxpi_tox21_sr_cols()
  )
  ld50_col <- .nts_toxpi_resolve_col(merged, "LD50_Zhu")
  fish_col <- .nts_toxpi_resolve_col(merged, "Fish_LC50")
  if (!is.na(ld50_col) && !is.na(fish_col)) {
    admet_missing <- is.na(merged[[ld50_col]])
    epi_missing <- is.na(merged[[fish_col]])
    if (any(admet_missing & !epi_missing) || any(epi_missing & !admet_missing)) {
      stop(
        "Misaligned merge detected: ",
        sum(admet_missing & !epi_missing), " row(s) EPI-only, ",
        sum(epi_missing & !admet_missing), " ADMET-only. ",
        "Use canonicalize_smiles() before predict_*() or align_by = 'order'.",
        call. = FALSE
      )
    }
  }

  merged
}

#' @keywords internal
#' @noRd
.nts_toxpi_short_smiles_label <- function(x, max_len = 40L) {
  x <- trimws(as.character(x))
  long <- !is.na(x) & nzchar(x) & nchar(x) > max_len
  if (any(long)) {
    x[long] <- paste0(substr(x[long], 1L, max_len), "...")
  }
  x
}

#' Per-row compound labels for ToxPi / toxpiR (coalesce id columns, then SMILES)
#' @keywords internal
#' @noRd
.nts_toxpi_coalesce_names <- function(merged, compound_ids = NULL) {
  n <- nrow(merged)
  if (!n) {
    return(character())
  }
  if (!is.null(compound_ids)) {
    lab <- trimws(as.character(compound_ids))
    if (length(lab) != n) {
      stop(
        "`compound_ids` length (", length(lab), ") must match merged rows (", n, ").",
        call. = FALSE
      )
    }
    return(lab)
  }
  candidates <- c(
    "compound_id", "compound_id_epi", "compound_name", "compound_name_epi",
    "Lab_ID", "Lab_ID_epi", "lab_id", "lab_id_epi",
    "episuite_Name", "episuite_Name_epi", "Name", "Name_epi",
    "name", "name_epi"
  )
  out <- rep(NA_character_, n)
  for (col in candidates) {
    if (!col %in% names(merged)) {
      next
    }
    v <- trimws(as.character(merged[[col]]))
    fill <- (is.na(out) | !nzchar(out)) & !is.na(v) & nzchar(v)
    if (any(fill)) {
      out[fill] <- v[fill]
    }
  }
  fill <- is.na(out) | !nzchar(out)
  if (any(fill) && "query_SMILES" %in% names(merged)) {
    v <- .nts_toxpi_short_smiles_label(merged$query_SMILES)
    ok <- fill & !is.na(v) & nzchar(v)
    out[ok] <- v[ok]
  }
  fill <- is.na(out) | !nzchar(out)
  if (any(fill) && "SMILES" %in% names(merged)) {
    v <- .nts_toxpi_short_smiles_label(merged$SMILES)
    ok <- fill & !is.na(v) & nzchar(v)
    out[ok] <- v[ok]
  }
  fill <- is.na(out) | !nzchar(out)
  if (any(fill)) {
    out[fill] <- paste0("CMP_", sprintf("%03d", seq_len(sum(fill))))
  }
  out
}

#' Attach user compound ids to EPI/ADMET tables before merge (by row order)
#' @keywords internal
#' @noRd
.nts_toxpi_attach_compound_ids <- function(epi, admet, compound_ids, query_smiles = NULL) {
  cid <- trimws(as.character(compound_ids))
  ne <- nrow(epi)
  na <- nrow(admet)
  nq <- if (is.null(query_smiles)) 0L else length(query_smiles)
  if (nq > 0L && length(cid) == nq) {
    if (ne == nq) {
      epi$compound_id <- cid
    }
    if (na == nq) {
      admet$compound_id <- cid
    }
    if (!"compound_id" %in% names(epi) && !"compound_id" %in% names(admet)) {
      stop(
        "`compound_ids` length matches `query_smiles` but not nrow(epi) or nrow(admet).",
        call. = FALSE
      )
    }
    return(list(epi = epi, admet = admet))
  }
  if (length(cid) == ne) {
    epi$compound_id <- cid
    return(list(epi = epi, admet = admet))
  }
  if (length(cid) == na) {
    admet$compound_id <- cid
    return(list(epi = epi, admet = admet))
  }
  stop(
    "`compound_ids` length (", length(cid), ") must match nrow(epi) (", ne,
    "), nrow(admet) (", na, "), or length(query_smiles) (", nq, ").",
    call. = FALSE
  )
}

#' @keywords internal
#' @noRd
.nts_toxpi_norm_biowin <- function(vals, col_name, verbose = TRUE) {
  vals <- as.numeric(vals)
  col_name <- as.character(col_name)[1L]
  if (identical(col_name, "Biowin3")) {
    bad <- !is.na(vals) & (vals < 1 | vals > 5)
    if (any(bad)) {
      .nts_toxpi_warn(
        col_name, ": timeframe Biowin outside [1, 5]; clipped before scoring.",
        verbose = verbose
      )
      vals[bad] <- pmin(pmax(vals[bad], 1), 5)
    }
    return(.nts_biowin3_to_risk(vals))
  }
  bad <- !is.na(vals) & (vals < 0 | vals > 1)
  if (any(bad)) {
    .nts_toxpi_warn(
      col_name, ": probability Biowin outside [0, 1]; clipped before scoring.",
      verbose = verbose
    )
  }
  .nts_biowin_persistence_risk(vals, col_name)
}

#' @keywords internal
#' @noRd
.nts_toxpi_normalize_minmax <- function(x, direction = c("positive", "negative")) {
  direction <- match.arg(direction)
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1L]) || !is.finite(rng[2L]) || diff(rng) == 0) {
    return(rep(0.5, length(x)))
  }
  score <- (x - rng[1L]) / (rng[2L] - rng[1L])
  if (identical(direction, "negative")) {
    score <- 1 - score
  }
  score
}

#' @keywords internal
#' @noRd
.nts_toxpi_aquatic_log <- function(x) {
  x <- as.numeric(x)
  bad <- !is.na(x) & x <= 0
  if (any(bad)) {
    x[bad] <- NA_real_
  }
  -log10(pmax(x, 1e-6, na.rm = FALSE))
}

#' @keywords internal
#' @noRd
.nts_toxpi_clip01 <- function(x) {
  x <- as.numeric(x)
  x[!is.na(x) & x < 0] <- 0
  x[!is.na(x) & x > 1] <- 1
  x
}

#' @keywords internal
#' @noRd
.nts_toxpi_prep_lc50 <- function(x, col_name, verbose = TRUE) {
  x <- as.numeric(x)
  bad <- !is.na(x) & x <= 0
  if (any(bad)) {
    .nts_toxpi_warn(col_name, ": replaced ", sum(bad), " non-positive value(s) with NA.", verbose = verbose)
    x[bad] <- NA_real_
  }
  high <- !is.na(x) & x > 1000
  if (any(high)) {
    .nts_toxpi_msg(col_name, ": ", sum(high), " value(s) > 1000 mg/L (low acute concern).", verbose = verbose)
  }
  .nts_toxpi_aquatic_log(x)
}

#' @keywords internal
#' @noRd
.nts_toxpi_slice_rowmean <- function(norm_mat) {
  if (!is.matrix(norm_mat)) {
    norm_mat <- as.matrix(norm_mat)
  }
  if (!ncol(norm_mat)) {
    return(rep(NA_real_, nrow(norm_mat)))
  }
  apply(norm_mat, 1L, function(row) {
    if (all(is.na(row))) {
      NA_real_
    } else {
      mean(row, na.rm = TRUE)
    }
  })
}

#' @keywords internal
#' @noRd
.nts_toxpi_bioaccum_col <- function(df, bcf_baf = c("baf", "bcf", "max")) {
  bcf_baf <- match.arg(bcf_baf)
  bcf_col <- .nts_toxpi_resolve_col(df, "LogBCF")
  baf_col <- .nts_toxpi_resolve_col(df, "LogBAF")
  if (identical(bcf_baf, "bcf") && !is.na(bcf_col)) {
    return(df[[bcf_col]])
  }
  if (identical(bcf_baf, "baf") && !is.na(baf_col)) {
    return(df[[baf_col]])
  }
  bcf_vals <- if (!is.na(bcf_col)) df[[bcf_col]] else rep(NA_real_, nrow(df))
  baf_vals <- if (!is.na(baf_col)) df[[baf_col]] else rep(NA_real_, nrow(df))
  pmax(bcf_vals, baf_vals, na.rm = TRUE)
}

#' @keywords internal
#' @noRd
.nts_toxpi_write_csv <- function(data, weights, outfile) {
  slice_cols <- setdiff(names(data), c("Name", "SMILES", "episuite_Name"))
  header <- c("Name", slice_cols)
  weight_row <- c("Weights", as.character(unname(weights[slice_cols])))
  mat <- data[, c("Name", slice_cols), drop = FALSE]
  lines <- c(
    paste(header, collapse = ","),
    paste(weight_row, collapse = ","),
    apply(mat, 1L, function(row) paste(row, collapse = ","))
  )
  writeLines(lines, outfile, useBytes = TRUE)
  invisible(outfile)
}

# ---- exported -----------------------------------------------------------------

#' Validate EPI Suite and ADMET-AI tables before ToxPi preparation
#'
#' Checks that required endpoint columns are present (with optional `make.names`
#' aliases for ADMET Tox21 hyphenated names). Does not run normalization.
#'
#' @param epi_data Output of [predict_episuite()] or compatible `data.frame`.
#' @param admet_data Output of [predict_admet_risk()] with `backend = "admet_ai"`.
#' @param tox21_agg Tox21 aggregation mode passed to [prepare_toxpi()].
#'
#' @return Invisibly `TRUE` if all required columns resolve; otherwise `stop()`.
#' @export
validate_toxpi_input <- function(epi_data,
                                 admet_data,
                                 tox21_agg = c("max", "mean", "split")) {
  tox21_agg <- match.arg(tox21_agg)
  epi <- as.data.frame(epi_data, stringsAsFactors = FALSE)
  admet <- as.data.frame(admet_data, stringsAsFactors = FALSE)

  if (!"SMILES" %in% names(epi)) {
    stop("epi_data missing column: SMILES", call. = FALSE)
  }
  if (!"SMILES" %in% names(admet)) {
    stop("admet_data missing column: SMILES", call. = FALSE)
  }

  epi_missing <- character()
  for (col in .nts_toxpi_required_epi()) {
    if (is.na(.nts_toxpi_resolve_col(epi, col))) {
      if (identical(col, "LogKow_est") &&
          !is.na(.nts_toxpi_resolve_col(epi, "LogP_pred"))) {
        next
      }
      epi_missing <- c(epi_missing, col)
    }
  }

  admet_missing <- character()
  for (col in .nts_toxpi_required_admet()) {
    if (is.na(.nts_toxpi_resolve_col(admet, col))) {
      admet_missing <- c(admet_missing, col)
    }
  }

  tox21_cols <- c(.nts_toxpi_tox21_nr_cols(), .nts_toxpi_tox21_sr_cols())
  tox21_missing <- vapply(
    tox21_cols,
    function(col) is.na(.nts_toxpi_resolve_col(admet, col)),
    logical(1L)
  )
  if (any(tox21_missing)) {
    admet_missing <- c(admet_missing, tox21_cols[tox21_missing])
  }

  if (length(epi_missing) || length(admet_missing)) {
    stop(
      "Missing required columns - EPI: ",
      if (length(epi_missing)) paste(epi_missing, collapse = ", ") else "(none)",
      "; ADMET: ",
      if (length(admet_missing)) paste(admet_missing, collapse = ", ") else "(none)",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' Prepare ToxPi slice scores from EPI Suite and ADMET-AI predictions
#'
#' Merges [predict_episuite()] and [predict_admet_risk()] (`backend = "admet_ai"`)
#' outputs, normalizes endpoints within the batch, aggregates six ToxPi slices
#' (Persistence, BioDegradation, AquaticTox, SystemicTox, Genotox, Tox21), and
#' BioDegradation uses EPI BIOWIN 1-7 (Biowin3 timeframe bins; others `1-p`).
#' `Ready_biodeg_risk` is a diagnostic column on `$data` (raw Biowin3 >= 3 & Biowin5 >= 0.5),
#' not averaged into the BioDegradation slice. See [analyze_episuite_biodegradation()].
#'
#' Join uses row order when `nrow(epi) == nrow(admet)` (default **`align_by = "auto"`**),
#' or structural keys (InChIKey / canonical SMILES) otherwise. **Recommend**
#' [canonicalize_smiles()] before [predict_episuite()] / [predict_admet_risk()].
#'
#' **SystemicTox** scores from ADMET-AI may sit in a narrow range for PFAS
#' (training-domain limitation); interpret comparatively, not as absolute risk.
#'
#' @param epi_data EPI Suite wide table (`predict_episuite()`).
#' @param admet_data ADMET-AI wide table (`predict_admet_risk(..., backend = "admet_ai")`).
#' @param query_smiles Optional original query SMILES (same order as prediction input).
#' @param compound_ids Optional compound labels (`Lab_ID`, `name`, etc.) in the same order
#'   as `query_smiles` or `nrow(epi_data)`; used when `episuite_Name` is missing or all `NA`.
#' @param align_by Merge strategy: `"auto"`, `"order"`, `"structure"`, or `"query_smiles"`.
#' @param scenario Weight preset: `"environmental"`, `"drug_safety"`, `"pesticide"`, or `"equal"`.
#' @param tox21_agg Tox21 aggregation: `"max"`, `"mean"`, or `"split"` (NR / SR sub-slices).
#' @param bcf_baf Bioaccumulation column for Persistence: `"baf"`, `"bcf"`, or `"max"`.
#' @param output_format `"toxpiR"` (list) or `"csv"` (also writes when `outfile` set).
#' @param outfile Optional path for desktop ToxPi CSV (Name + slice scores + weight row).
#' @param verbose If `FALSE`, suppress `message()` / `warning()` from this function.
#'
#' @return A list with:
#'   \describe{
#'     \item{`data`}{`data.frame` of slice scores (0-1) plus `Name` / `SMILES`.}
#'     \item{`slices`}{Named list mapping slice names to contributing normalized column names.}
#'     \item{`weights`}{Named numeric vector of slice weights (sum 1).}
#'   }
#'
#' @seealso [validate_toxpi_input], [canonicalize_smiles], [analyze_episuite_biodegradation],
#'   [predict_episuite], [predict_admet_risk]
#' @export
prepare_toxpi <- function(epi_data,
                          admet_data,
                          query_smiles = NULL,
                          compound_ids = NULL,
                          align_by = c("auto", "order", "structure", "query_smiles"),
                          scenario = c("environmental", "drug_safety", "pesticide", "equal"),
                          tox21_agg = c("max", "mean", "split"),
                          bcf_baf = c("baf", "bcf", "max"),
                          output_format = c("toxpiR", "csv"),
                          outfile = NULL,
                          verbose = TRUE) {
  scenario <- match.arg(scenario)
  tox21_agg <- match.arg(tox21_agg)
  bcf_baf <- match.arg(bcf_baf)
  output_format <- match.arg(output_format)
  align_by <- match.arg(align_by)

  validate_toxpi_input(epi_data, admet_data, tox21_agg = tox21_agg)

  if (!is.null(compound_ids)) {
    attached <- .nts_toxpi_attach_compound_ids(
      as.data.frame(epi_data, stringsAsFactors = FALSE),
      as.data.frame(admet_data, stringsAsFactors = FALSE),
      compound_ids,
      query_smiles = query_smiles
    )
    epi_data <- attached$epi
    admet_data <- attached$admet
  }

  merged <- .nts_toxpi_merge_inputs(
    epi_data,
    admet_data,
    align_by = align_by,
    query_smiles = query_smiles,
    verbose = verbose
  )
  n <- nrow(merged)
  if (n < 3L) {
    .nts_toxpi_warn(
      "Fewer than 3 compounds (n=", n, "); batch min-max normalization may be unstable.",
      verbose = verbose
    )
  }

  name_col <- .nts_toxpi_coalesce_names(merged, compound_ids = NULL)

  slice_defs <- list()
  slice_scores <- data.frame(
    Name = name_col,
    SMILES = as.character(merged$SMILES),
    stringsAsFactors = FALSE
  )

  # --- Persistence ---
  logp_col <- .nts_toxpi_resolve_col(merged, c("LogKow_est", "LogP_pred"))
  koc_col <- .nts_toxpi_resolve_col(merged, "LogKoc_est")
  koa_col <- .nts_toxpi_resolve_col(merged, "LogKoa_est")
  bioaccum <- .nts_toxpi_bioaccum_col(merged, bcf_baf = bcf_baf)
  pers_mat <- cbind(
    LogP = .nts_toxpi_normalize_minmax(merged[[logp_col]], "positive"),
    LogKoc = .nts_toxpi_normalize_minmax(merged[[koc_col]], "positive"),
    Bioaccum = .nts_toxpi_normalize_minmax(bioaccum, "positive"),
    LogKoa = .nts_toxpi_normalize_minmax(merged[[koa_col]], "positive")
  )
  colnames(pers_mat) <- c("LogP_norm", "LogKoc_norm", "Bioaccum_norm", "LogKoa_norm")
  slice_scores$Persistence <- .nts_toxpi_slice_rowmean(pers_mat)
  slice_defs$Persistence <- colnames(pers_mat)
  .nts_toxpi_msg(
    "Persistence score range: ",
    paste(round(range(slice_scores$Persistence, na.rm = TRUE), 3), collapse = " - "),
    verbose = verbose
  )

  # --- BioDegradation (Biowin1-7 only; Ready_biodeg_risk is diagnostic, not in slice mean) ---
  biowin_cols <- paste0("Biowin", seq_len(7L))
  biowin_cols <- biowin_cols[!vapply(biowin_cols, function(c) is.na(.nts_toxpi_resolve_col(merged, c)), logical(1L))]
  bio_mat <- do.call(cbind, lapply(biowin_cols, function(col) {
    resolved <- .nts_toxpi_resolve_col(merged, col)
    .nts_toxpi_norm_biowin(merged[[resolved]], col, verbose = verbose)
  }))
  if (!is.null(bio_mat) && ncol(bio_mat)) {
    colnames(bio_mat) <- paste0(biowin_cols, "_norm")
    slice_scores$BioDegradation <- .nts_toxpi_slice_rowmean(bio_mat)
    slice_defs$BioDegradation <- colnames(bio_mat)
    .nts_toxpi_msg(
      "BioDegradation score range: ",
      paste(round(range(slice_scores$BioDegradation, na.rm = TRUE), 3), collapse = " - "),
      verbose = verbose
    )
  }

  b3_col <- .nts_toxpi_resolve_col(merged, "Biowin3")
  b5_col <- .nts_toxpi_resolve_col(merged, "Biowin5")
  if (!is.na(b3_col) && !is.na(b5_col)) {
    slice_scores$Ready_biodeg_risk <- mapply(
      .nts_biowin_ready_biodeg_risk,
      merged[[b3_col]],
      merged[[b5_col]],
      USE.NAMES = FALSE
    )
  }

  # --- AquaticTox ---
  lc50_cols <- c(
    "Fish_LC50", "Daphnid_LC50", "Algae_LC50",
    "Fish_chV", "Daphnid_chV", "Algae_chV"
  )
  aq_mat <- do.call(cbind, lapply(lc50_cols, function(col) {
    resolved <- .nts_toxpi_resolve_col(merged, col)
    transformed <- .nts_toxpi_prep_lc50(merged[[resolved]], resolved, verbose = verbose)
    .nts_toxpi_normalize_minmax(transformed, "positive")
  }))
  colnames(aq_mat) <- paste0(lc50_cols, "_norm")
  slice_scores$AquaticTox <- .nts_toxpi_slice_rowmean(aq_mat)
  slice_defs$AquaticTox <- colnames(aq_mat)
  .nts_toxpi_msg(
    "AquaticTox score range: ",
    paste(round(range(slice_scores$AquaticTox, na.rm = TRUE), 3), collapse = " - "),
    verbose = verbose
  )

  # --- SystemicTox ---
  ld50_col <- .nts_toxpi_resolve_col(merged, "LD50_Zhu")
  sys_cols <- c("DILI", "ClinTox", "hERG")
  sys_mat <- cbind(
    LD50_Zhu_norm = .nts_toxpi_normalize_minmax(merged[[ld50_col]], "negative")
  )
  for (col in sys_cols) {
    resolved <- .nts_toxpi_resolve_col(merged, col)
    clipped <- .nts_toxpi_clip01(merged[[resolved]])
    out_rng <- clipped[!is.na(clipped)] < 0 | clipped[!is.na(clipped)] > 1
    if (length(out_rng) && any(out_rng)) {
      .nts_toxpi_warn(resolved, ": ADMET probability outside [0, 1]; clipped.", verbose = verbose)
    }
    sys_mat <- cbind(
      sys_mat,
      .nts_toxpi_normalize_minmax(clipped, "positive")
    )
    colnames(sys_mat)[ncol(sys_mat)] <- paste0(col, "_norm")
  }
  slice_scores$SystemicTox <- .nts_toxpi_slice_rowmean(sys_mat)
  slice_defs$SystemicTox <- colnames(sys_mat)
  .nts_toxpi_msg(
    "SystemicTox score range: ",
    paste(round(range(slice_scores$SystemicTox, na.rm = TRUE), 3), collapse = " - "),
    verbose = verbose
  )

  # --- Genotox ---
  ames_col <- .nts_toxpi_resolve_col(merged, "AMES")
  carc_col <- .nts_toxpi_resolve_col(merged, "Carcinogens_Lagunin")
  geno_mat <- cbind(
    AMES_norm = .nts_toxpi_clip01(merged[[ames_col]]),
    Carcinogens_Lagunin_norm = .nts_toxpi_clip01(merged[[carc_col]])
  )
  slice_scores$Genotox <- .nts_toxpi_slice_rowmean(geno_mat)
  slice_defs$Genotox <- colnames(geno_mat)
  .nts_toxpi_msg(
    "Genotox score range: ",
    paste(round(range(slice_scores$Genotox, na.rm = TRUE), 3), collapse = " - "),
    verbose = verbose
  )

  # --- Tox21 ---
  nr_cols <- vapply(
    .nts_toxpi_tox21_nr_cols(),
    function(col) .nts_toxpi_resolve_col(merged, col),
    character(1L)
  )
  sr_cols <- vapply(
    .nts_toxpi_tox21_sr_cols(),
    function(col) .nts_toxpi_resolve_col(merged, col),
    character(1L)
  )
  .nts_toxpi_rowmax <- function(mat) {
    apply(mat, 1L, function(row) {
      if (all(is.na(row))) {
        NA_real_
      } else {
        max(row, na.rm = TRUE)
      }
    })
  }
  nr_mat <- do.call(cbind, lapply(nr_cols, function(col) {
    .nts_toxpi_normalize_minmax(.nts_toxpi_clip01(merged[[col]]), "positive")
  }))
  sr_mat <- do.call(cbind, lapply(sr_cols, function(col) {
    .nts_toxpi_normalize_minmax(.nts_toxpi_clip01(merged[[col]]), "positive")
  }))
  if (!is.null(nr_mat)) {
    colnames(nr_mat) <- paste0(nr_cols, "_norm")
  }
  if (!is.null(sr_mat)) {
    colnames(sr_mat) <- paste0(sr_cols, "_norm")
  }

  if (identical(tox21_agg, "split")) {
    slice_scores$Tox21_NR <- .nts_toxpi_rowmax(nr_mat)
    slice_scores$Tox21_SR <- .nts_toxpi_rowmax(sr_mat)
    slice_defs$Tox21_NR <- colnames(nr_mat)
    slice_defs$Tox21_SR <- colnames(sr_mat)
    .nts_toxpi_msg(
      "Tox21_NR / Tox21_SR score ranges: ",
      paste(round(range(slice_scores$Tox21_NR, na.rm = TRUE), 3), collapse = " - "),
      " / ",
      paste(round(range(slice_scores$Tox21_SR, na.rm = TRUE), 3), collapse = " - "),
      verbose = verbose
    )
  } else {
    all_mat <- cbind(nr_mat, sr_mat)
    if (identical(tox21_agg, "max")) {
      slice_scores$Tox21 <- .nts_toxpi_rowmax(all_mat)
    } else {
      slice_scores$Tox21 <- .nts_toxpi_slice_rowmean(all_mat)
    }
    slice_defs$Tox21 <- colnames(all_mat)
    .nts_toxpi_msg(
      "Tox21 score range: ",
      paste(round(range(slice_scores$Tox21, na.rm = TRUE), 3), collapse = " - "),
      verbose = verbose
    )
  }

  slice_names <- setdiff(names(slice_scores), c("Name", "SMILES"))
  weights <- .nts_toxpi_scenario_weights(scenario, slice_names)
  rownames(slice_scores) <- NULL

  out <- list(
    data = slice_scores,
    slices = slice_defs,
    weights = weights
  )

  if (!is.null(outfile) || identical(output_format, "csv")) {
    path <- outfile
    if (is.null(path)) {
      stop("output_format = 'csv' requires a non-NULL `outfile`.", call. = FALSE)
    }
    .nts_toxpi_write_csv(slice_scores, weights, path)
    out$outfile <- path
    .nts_toxpi_msg("Wrote ToxPi CSV: ", path, verbose = verbose)
  }

  out
}
