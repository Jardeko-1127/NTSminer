# External toxicity / fate prediction engines (EPI Suite, ADMETlab 3.0, ADMET-AI)
# Ref: tests/toxPi_202506065.R, tests/toxpred/

#' Default EPI Suite submit API base URL
#' @keywords internal
#' @noRd
.nts_episuite_default_base_url <- function() {
  "https://episuite.dev/EpiWebSuite/api/submit"
}

#' Default ADMETlab 3.0 API base URL
#' @keywords internal
#' @noRd
.nts_admetlab_default_base <- function() {
  "https://admetlab3.scbdd.com"
}

#' Random delay between API calls (rate limiting)
#' @keywords internal
#' @noRd
.nts_toxpred_sleep_jitter <- function(range_sec = c(1.5, 3.0)) {
  if (length(range_sec) < 2L) {
    return(invisible(NULL))
  }
  Sys.sleep(stats::runif(1L, min = range_sec[1L], max = range_sec[2L]))
  invisible(NULL)
}

#' @keywords internal
#' @noRd
.nts_require_httr_json <- function(fn_label = "toxpred") {
  if (!requireNamespace("httr", quietly = TRUE) ||
      !requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      fn_label, ": requires packages `httr` and `jsonlite`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# ---- EPI Suite ----------------------------------------------------------------

#' Call EPI Suite submit API for one SMILES
#' @keywords internal
#' @noRd
.nts_episuite_call <- function(smiles,
                                base_url = .nts_episuite_default_base_url(),
                                max_retries = 3L,
                                delay_sec = c(1.5, 3.0),
                                timeout_sec = 90) {
  .nts_require_httr_json("EPI Suite")
  smiles <- trimws(as.character(smiles)[1L])
  if (is.na(smiles) || !nzchar(smiles)) {
    return(NULL)
  }
  url <- paste0(base_url, "?smiles=", utils::URLencode(smiles, reserved = TRUE))
  for (retry in seq_len(max_retries)) {
    if (retry > 1L) {
      .nts_toxpred_sleep_jitter(delay_sec * 1.5)
    } else {
      .nts_toxpred_sleep_jitter(delay_sec)
    }
    response <- tryCatch(
      httr::GET(
        url,
        httr::add_headers(Accept = "application/json"),
        httr::timeout(timeout_sec)
      ),
      error = function(e) {
        message("EPI Suite request error: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(response)) {
      next
    }
    if (httr::status_code(response) == 200L) {
      return(httr::content(response, as = "parsed", type = "application/json", encoding = "UTF-8"))
    }
    message(
      "EPI Suite HTTP ", httr::status_code(response),
      " for SMILES: ", substr(smiles, 1L, 50L)
    )
  }
  NULL
}

#' Flatten one EPI Suite JSON response to a wide data.frame row
#' @keywords internal
#' @noRd
.nts_episuite_row_from_response <- function(api_response, smiles_query = NA_character_) {
  if (is.null(api_response)) {
    return(data.frame(SMILES = smiles_query, stringsAsFactors = FALSE))
  }
  chem <- api_response$chemicalProperties
  basic <- data.frame(
    episuite_Name = if (is.null(chem$name) || !nzchar(chem$name)) "Unknown" else chem$name,
    episuite_CAS = if (is.null(chem$cas)) NA_character_ else as.character(chem$cas),
    SMILES = if (is.null(chem$smiles)) smiles_query else chem$smiles,
    episuite_Weight = chem$molecularWeight,
    episuite_Formula = chem$molecularFormula,
    stringsAsFactors = FALSE
  )
  phys <- data.frame(
    LogKow_est = api_response$logKow$estimatedValue$value,
    LogKoc_est = api_response$logKoc$estimatedValue$value,
    MP_est_C = api_response$meltingPoint$estimatedValue$value,
    BP_est_C = api_response$boilingPoint$estimatedValue$value,
    VP_est_mmHg = api_response$vaporPressure$estimatedValue$value,
    WS_est_mgL = api_response$waterSolubilityFromLogKow$estimatedValue$value,
    HLC_est = api_response$henrysLawConstant$estimatedValue$value,
    LogKoa_est = api_response$logKoa$estimatedValue$value,
    stringsAsFactors = FALSE
  )
  bio_data <- api_response$biodegradationRate$models
  biowin <- rep(NA_real_, 7L)
  if (!is.null(bio_data)) {
    vals <- vapply(bio_data, function(x) x$value, numeric(1L))
    biowin[seq_len(min(7L, length(vals)))] <- vals[seq_len(min(7L, length(vals)))]
  }
  names(biowin) <- paste0("Biowin", seq_len(7L))
  biocon <- api_response$bioconcentration
  bcf <- c(LogBCF = NA_real_, LogBAF = NA_real_)
  if (!is.null(biocon)) {
    bcf["LogBCF"] <- biocon$logBioconcentrationFactor
    bcf["LogBAF"] <- biocon$logBioaccumulationFactor
  }
  ecosar <- c(
    Fish_LC50 = NA_real_, Daphnid_LC50 = NA_real_, Algae_LC50 = NA_real_,
    Fish_chV = NA_real_, Daphnid_chV = NA_real_, Algae_chV = NA_real_
  )
  if (!is.null(api_response$ecosar) && !is.null(api_response$ecosar$modelResults)) {
    er <- api_response$ecosar$modelResults
    if (length(er) >= 6L) {
      conc <- vapply(er, function(x) x$concentration, numeric(1L))
      ecosar <- conc[seq_len(6L)]
      names(ecosar) <- c(
        "Fish_LC50", "Daphnid_LC50", "Algae_LC50",
        "Fish_chV", "Daphnid_chV", "Algae_chV"
      )
    }
  }
  out <- cbind(
    basic,
    phys,
    as.data.frame(t(biowin), stringsAsFactors = FALSE),
    as.data.frame(t(bcf), stringsAsFactors = FALSE),
    as.data.frame(t(ecosar), stringsAsFactors = FALSE)
  )
  out$LogP_pred <- out$LogKow_est
  rownames(out) <- NULL
  out
}

#' Classify EPI BIOWIN timeframe value (Biowin3 only)
#' @keywords internal
#' @noRd
.nts_biowin_timeframe_label <- function(value) {
  v <- as.numeric(value)[1L]
  if (is.na(v)) {
    return(NA_character_)
  }
  if (v >= 5) {
    "Hours"
  } else if (v >= 4) {
    "Days"
  } else if (v >= 3) {
    "Weeks"
  } else if (v >= 2) {
    "Months"
  } else {
    "Longer"
  }
}

#' Classify a single EPI BIOWIN model result
#' @keywords internal
#' @noRd
.nts_biowin_classify_model <- function(model_name, value) {
  model_name <- as.character(model_name)[1L]
  value <- as.numeric(value)[1L]
  if (is.na(model_name) || !nzchar(model_name) || is.na(value)) {
    return(list(
      Classification = "Unknown",
      Criterion = "No specific criterion available"
    ))
  }
  if (grepl("Linear|Non-Linear|Anaerobic", model_name, ignore.case = TRUE) &&
      !grepl("MITI", model_name, ignore.case = TRUE)) {
    return(list(
      Classification = if (value >= 0.5) {
        "Biodegrades fast"
      } else {
        "Does not biodegrade fast"
      },
      Criterion = paste(
        "Greater than or equal to 0.5 indicates 'Biodegrades fast';",
        "Less than 0.5 indicates 'Does not biodegrade fast'"
      )
    ))
  }
  if (grepl("MITI", model_name, ignore.case = TRUE)) {
    return(list(
      Classification = if (value >= 0.5) {
        "Readily degradable"
      } else {
        "Not readily degradable"
      },
      Criterion = paste(
        "Greater than or equal to 0.5 indicates 'Readily degradable';",
        "Less than 0.5 indicates 'Not readily degradable'"
      )
    ))
  }
  if (grepl("Timeframe", model_name, ignore.case = TRUE)) {
    tf <- .nts_biowin_timeframe_label(value)
    return(list(
      Classification = paste0("About ", tf),
      Criterion = "5.00 -> Hours; 4.00 -> Days; 3.00 -> Weeks; 2.00 -> Months; 1.00 -> Longer"
    ))
  }
  list(
    Classification = "Unknown",
    Criterion = "No specific criterion available"
  )
}

#' Ready biodegradability from Biowin3 timeframe + Biowin5 probability
#' @keywords internal
#' @noRd
.nts_biowin_ready_biodegradability <- function(biowin3, biowin5) {
  b3 <- as.numeric(biowin3)[1L]
  b5 <- as.numeric(biowin5)[1L]
  if (is.na(b3) || is.na(b5)) {
    return(list(
      prediction = NA_character_,
      biowin3_timeframe = if (is.na(b3)) NA_character_ else .nts_biowin_timeframe_label(b3),
      biowin3_meets = NA,
      biowin5_meets = NA,
      ready_biodegradable = NA
    ))
  }
  meets3 <- !is.na(b3) && b3 >= 3
  meets5 <- !is.na(b5) && b5 >= 0.5
  ready <- isTRUE(meets3) && isTRUE(meets5)
  list(
    prediction = if (ready) "YES" else "NO",
    biowin3_timeframe = if (is.na(b3)) NA_character_ else .nts_biowin_timeframe_label(b3),
    biowin3_meets = meets3,
    biowin5_meets = meets5,
    ready_biodegradable = ready
  )
}

#' Biowin3 ultimate timeframe (1-5) to persistence risk score
#' @keywords internal
#' @noRd
.nts_biowin3_to_risk <- function(value) {
  v <- as.numeric(value)
  out <- rep(NA_real_, length(v))
  ok <- !is.na(v)
  out[ok & v >= 5] <- 0
  out[ok & v >= 4 & v < 5] <- 0.25
  out[ok & v >= 3 & v < 4] <- 0.5
  out[ok & v >= 2 & v < 3] <- 0.75
  out[ok & v < 2] <- 1
  out
}

#' BIOWIN persistence-risk score from wide-table column (EPI standard order)
#' @keywords internal
#' @noRd
.nts_biowin_persistence_risk <- function(value, biowin_col) {
  col <- as.character(biowin_col)[1L]
  if (identical(col, "Biowin3")) {
    bad <- !is.na(value) & (value < 1 | value > 5)
    if (any(bad)) {
      value <- as.numeric(value)
      value[bad] <- pmin(pmax(value[bad], 1), 5)
    }
    return(.nts_biowin3_to_risk(value))
  }
  p <- as.numeric(value)
  p[!is.na(p) & p < 0] <- 0
  p[!is.na(p) & p > 1] <- 1
  1 - p
}

#' Persistence-risk score for Ready biodegradability (0 = YES, 1 = NO)
#' @keywords internal
#' @noRd
.nts_biowin_ready_biodeg_risk <- function(biowin3, biowin5) {
  rd <- .nts_biowin_ready_biodegradability(biowin3, biowin5)
  if (is.na(rd$prediction)) {
    return(NA_real_)
  }
  if (isTRUE(rd$ready_biodegradable)) {
    0
  } else {
    1
  }
}

#' Analyze EPI Suite BIOWIN / ready biodegradability results
#'
#' Parses `biodegradationRate$models` from an EPI Suite API response, classifies
#' each BIOWIN model (probability vs timeframe), and applies the EPI **Ready
#' Biodegradability** rule (Biowin3 weeks-or-faster AND Biowin5 >= 0.5).
#'
#' @param api_response Parsed JSON from [predict_episuite()] / EPI submit API.
#'
#' @return A list with:
#'   \describe{
#'     \item{`bio_results`}{Per-model `data.frame` (`Name`, `Model`, `Value`, `Classification`, `Criterion`).}
#'     \item{`prediction_result`}{Ready biodegradability summary `data.frame` or `"Insufficient data"`.}
#'   }
#'
#' @seealso [predict_episuite], [prepare_toxpi]
#' @export
analyze_episuite_biodegradation <- function(api_response) {
  if (is.null(api_response) || is.null(api_response$biodegradationRate$models)) {
    return(list(
      bio_results = data.frame(
        Name = paste0("Biowin", seq_len(7L)),
        Model = NA_character_,
        Value = NA_real_,
        Classification = NA_character_,
        Criterion = NA_character_,
        stringsAsFactors = FALSE
      ),
      prediction_result = "Insufficient data to determine"
    ))
  }

  bio_data <- api_response$biodegradationRate$models
  bio_results <- data.frame(
    Name = paste0("Biowin", seq_len(length(bio_data))),
    Model = vapply(bio_data, function(x) as.character(x$name), character(1L)),
    Value = vapply(bio_data, function(x) as.numeric(x$value), numeric(1L)),
    stringsAsFactors = FALSE
  )

  cls <- lapply(seq_len(nrow(bio_results)), function(i) {
    .nts_biowin_classify_model(bio_results$Model[i], bio_results$Value[i])
  })
  bio_results$Classification <- vapply(cls, function(x) x$Classification, character(1L))
  bio_results$Criterion <- vapply(cls, function(x) x$Criterion, character(1L))

  b3_idx <- grepl("Ultimate Biodegradation Timeframe", bio_results$Model, ignore.case = TRUE)
  b5_idx <- grepl("MITI Linear Model Prediction", bio_results$Model, ignore.case = TRUE)
  b3_val <- bio_results$Value[b3_idx][1L]
  b5_val <- bio_results$Value[b5_idx][1L]

  prediction_result <- if (length(b3_val) && !is.na(b3_val) && length(b5_val) && !is.na(b5_val)) {
    rd <- .nts_biowin_ready_biodegradability(b3_val, b5_val)
    data.frame(
      Results = paste0("Ready Biodegradability Prediction: ", rd$prediction),
      Criterion = paste(
        "Criteria for the YES or NO prediction: If the Biowin3 (ultimate survey model)",
        "result is 'weeks' or faster (i.e. 'days', 'days to weeks', or 'weeks') AND",
        "the Biowin5 (MITI linear model) probability is >= 0.5, then the prediction is",
        "YES (readily biodegradable). Otherwise NO."
      ),
      stringsAsFactors = FALSE
    )
  } else {
    "Insufficient data to determine"
  }

  list(
    bio_results = bio_results,
    prediction_result = prediction_result
  )
}

#' Predict environmental fate / ecotox via EPI Suite (batch)
#'
#' Calls the public EPI Suite web API (`episuite.dev`) for each unique SMILES.
#' Suitable for PFAS, PPCPs, pesticides, and other environmental contaminants
#' when a valid SMILES is available.
#'
#' @param smiles_list Character vector of SMILES strings.
#' @param base_url EPI Suite submit API base (default episuite.dev).
#' @param max_retries Retries per SMILES on transient HTTP errors.
#' @param delay_sec Random delay range (seconds) between requests.
#' @param timeout_sec HTTP timeout per request.
#' @param progress If `TRUE`, message progress.
#'
#' @return A `data.frame` with one row per SMILES (wide EPI + ECOSAR columns).
#'   Failed requests return a row with `SMILES` only and `NA` elsewhere.
#'
#' @seealso [predict_properties]
#' @export
predict_episuite <- function(smiles_list,
                             base_url = .nts_episuite_default_base_url(),
                             max_retries = 3L,
                             delay_sec = c(1.5, 3.0),
                             timeout_sec = 90,
                             progress = TRUE) {
  smiles <- unique(trimws(as.character(smiles_list)))
  smiles <- smiles[!is.na(smiles) & nzchar(smiles)]
  if (!length(smiles)) {
    return(data.frame())
  }
  rows <- vector("list", length(smiles))
  for (i in seq_along(smiles)) {
    if (isTRUE(progress)) {
      message("[EPI Suite ", i, "/", length(smiles), "] ", substr(smiles[i], 1L, 50L))
    }
    api <- .nts_episuite_call(
      smiles[i],
      base_url = base_url,
      max_retries = max_retries,
      delay_sec = delay_sec,
      timeout_sec = timeout_sec
    )
    rows[[i]] <- .nts_episuite_row_from_response(api, smiles_query = smiles[i])
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# ---- ADMETlab 3.0 -------------------------------------------------------------

#' @keywords internal
#' @noRd
.nts_admetlab_post <- function(url, body_list, timeout_sec = 120L) {
  .nts_require_httr_json("ADMETlab")
  body_raw <- jsonlite::toJSON(body_list, auto_unbox = TRUE, null = "null")
  tryCatch(
    httr::POST(
      url,
      body = body_raw,
      httr::content_type_json(),
      httr::timeout(timeout_sec)
    ),
    error = function(e) {
      message("ADMETlab POST error (", url, "): ", conditionMessage(e))
      NULL
    }
  )
}

#' @keywords internal
#' @noRd
.nts_admetlab_wash <- function(smiles, base = .nts_admetlab_default_base()) {
  url <- paste0(sub("/$", "", base), "/api/washmol")
  resp <- .nts_admetlab_post(url, list(smiles = as.character(smiles)))
  if (is.null(resp) || httr::status_code(resp) != 200L) {
    return(as.character(smiles))
  }
  parsed <- httr::content(resp, as = "parsed", type = "application/json", encoding = "UTF-8")
  washed <- parsed$smiles
  if (is.null(washed)) {
    washed <- parsed$data
  }
  if (is.character(washed) && length(washed) == length(smiles)) {
    return(washed)
  }
  as.character(smiles)
}

#' @keywords internal
#' @noRd
.nts_admetlab_parse <- function(parsed, smiles_query) {
  if (is.data.frame(parsed)) {
    if (!"SMILES" %in% names(parsed) && !is.null(rownames(parsed))) {
      parsed$SMILES <- rownames(parsed)
    }
    return(as.data.frame(parsed, stringsAsFactors = FALSE))
  }
  if (is.list(parsed) && !is.null(parsed$data) && is.data.frame(parsed$data)) {
    return(as.data.frame(parsed$data, stringsAsFactors = FALSE))
  }
  if (is.list(parsed) && length(parsed) == length(smiles_query)) {
    rows <- lapply(seq_along(smiles_query), function(i) {
      item <- parsed[[i]]
      if (is.null(item)) {
        return(data.frame(SMILES = smiles_query[i], stringsAsFactors = FALSE))
      }
      if (is.list(item) && !is.null(item$smiles)) {
        item$SMILES <- item$smiles
      }
      as.data.frame(item, stringsAsFactors = FALSE)
    })
    return(do.call(rbind, rows))
  }
  data.frame(SMILES = smiles_query, stringsAsFactors = FALSE)
}

#' Predict ADMET endpoints via ADMETlab 3.0 REST API
#'
#' @param smiles_list Character vector of SMILES.
#' @param base ADMETlab 3.0 server base URL.
#' @param endpoints API paths tried in order (default `/api/admet`, `/api/single/admet`, `/api/predict`).
#' @param batch_size SMILES per HTTP request (respect server rate limits; default 5).
#' @param delay_sec Pause between batches (seconds).
#' @param wash_first If `TRUE`, standardize SMILES via `/api/washmol` first.
#' @param timeout_sec HTTP timeout.
#' @param progress Message batch progress.
#'
#' @return Wide `data.frame`; failed batches have `admetlab_status = "failed"`.
#'
#' @seealso [predict_admet_risk]
#' @export
predict_admetlab3 <- function(smiles_list,
                              base = .nts_admetlab_default_base(),
                              endpoints = c("/api/admet", "/api/single/admet", "/api/predict"),
                              batch_size = 5L,
                              delay_sec = 0.25,
                              wash_first = TRUE,
                              timeout_sec = 120L,
                              progress = TRUE) {
  smiles <- unique(trimws(as.character(smiles_list)))
  smiles <- smiles[!is.na(smiles) & nzchar(smiles)]
  if (!length(smiles)) {
    return(data.frame())
  }
  batch_size <- max(1L, as.integer(batch_size)[1L])
  chunks <- split(smiles, ceiling(seq_along(smiles) / batch_size))

  .batch_one <- function(chunk) {
    sm <- as.character(chunk)
    if (isTRUE(wash_first)) {
      sm <- .nts_admetlab_wash(sm, base = base)
    }
    payloads <- list(
      list(smiles = sm),
      list(smi = sm),
      list(smiles = sm, model = "DMPNN")
    )
    for (ep in endpoints) {
      url <- paste0(sub("/$", "", base), ep)
      for (pl in payloads) {
        if (length(sm) == 1L && identical(ep, "/api/single/admet")) {
          pl <- list(smi = sm[[1L]])
        }
        resp <- .nts_admetlab_post(url, pl, timeout_sec = timeout_sec)
        if (is.null(resp)) {
          next
        }
        if (httr::status_code(resp) == 200L) {
          parsed <- httr::content(resp, as = "parsed", type = "application/json", encoding = "UTF-8")
          out <- .nts_admetlab_parse(parsed, sm)
          out$admetlab_backend <- "ADMETlab3"
          attr(out, "endpoint") <- ep
          return(out)
        }
      }
    }
    data.frame(SMILES = sm, admetlab_status = "failed", admetlab_backend = "ADMETlab3", stringsAsFactors = FALSE)
  }

  parts <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    if (isTRUE(progress)) {
      message("[ADMETlab3 ", i, "/", length(chunks), "] n=", length(chunks[[i]]))
    }
    parts[[i]] <- .batch_one(chunks[[i]])
    if (i < length(chunks)) {
      Sys.sleep(delay_sec)
    }
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# ---- ADMET-AI -----------------------------------------------------------------

#' @keywords internal
#' @noRd
.nts_admet_ai_cli_available <- function() {
  nzchar(Sys.which("admet_predict"))
}

#' Candidate python.exe paths for ADMET-AI conda env
#' @keywords internal
#' @noRd
.nts_admet_ai_python_candidates <- function(conda_env = "admet_env") {
  envs <- unique(c(
    conda_env,
    Sys.getenv("NTSminer_ADMET_AI_ENV", unset = "admet_env"),
    Sys.getenv("TOXPRED_ADMET_AI_ENV", unset = "")
  ))
  envs <- envs[nzchar(envs)]
  py_name <- if (.Platform$OS.type == "windows") "python.exe" else "bin/python"
  roots <- unique(c(
    Sys.getenv("CONDA_PREFIX", unset = NA_character_),
    dirname(Sys.getenv("CONDA_EXE", unset = NA_character_)),
    file.path(Sys.getenv("USERPROFILE", unset = "~"), "anaconda3"),
    file.path(Sys.getenv("USERPROFILE", unset = "~"), "miniconda3"),
    file.path(Sys.getenv("USERPROFILE", unset = "~"), ".conda"),
    file.path(Sys.getenv("LOCALAPPDATA", unset = ""), "miniconda3"),
    # System-wide install only; user-specific drive letters belong in the
    # CONDA_PREFIX / NTSminer_ADMET_AI_* env vars handled above.
    "C:/ProgramData/anaconda3"
  ))
  roots <- roots[!is.na(roots) & nzchar(roots)]
  c(
    Sys.getenv("NTSminer_ADMET_AI_PYTHON", unset = NA_character_),
    Sys.getenv("TOXPRED_ADMET_AI_PYTHON", unset = NA_character_),
    Sys.getenv("RETICULATE_PYTHON", unset = NA_character_),
    unlist(lapply(envs, function(env) {
      unlist(lapply(roots, function(root) {
        c(
          file.path(root, "envs", env, py_name),
          if (identical(basename(root), env)) file.path(root, py_name) else character()
        )
      }), use.names = FALSE)
    }), use.names = FALSE)
  )
}

#' Activate Python for ADMET-AI (explicit path, env var, or conda env)
#' @keywords internal
#' @noRd
.nts_admet_ai_configure_python <- function(conda_env = "admet_env", python_path = NULL) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("ADMET-AI requires the `reticulate` package.", call. = FALSE)
  }
  pick <- function(py) {
    if (!nzchar(py) || !file.exists(py)) {
      return(FALSE)
    }
    tryCatch({
      reticulate::use_python(py, required = TRUE)
      TRUE
    }, error = function(e) FALSE)
  }
  if (!is.null(python_path) && nzchar(python_path)) {
    if (!pick(python_path)) {
      stop("python_path not usable: ", python_path, call. = FALSE)
    }
    return(invisible(python_path))
  }
  for (py in unique(.nts_admet_ai_python_candidates(conda_env))) {
    if (pick(py)) {
      return(invisible(py))
    }
  }
  tryCatch({
    reticulate::use_condaenv(conda_env, required = TRUE)
    invisible(Sys.getenv("RETICULATE_PYTHON", unset = NA_character_))
  }, error = function(e) {
    stop(
      "ADMET-AI Python not found. As in tests/toxPi_202506065.R, use conda env `",
      conda_env, "` or set Sys.setenv(NTSminer_ADMET_AI_PYTHON = ",
      "'<path/to/admet_env/python.exe>'). Original error: ",
      conditionMessage(e),
      call. = FALSE
    )
  })
}

#' @keywords internal
#' @noRd
.nts_admet_ai_is_py <- function(x) {
  cls <- class(x)
  any(grepl("^python\\.", cls)) || inherits(x, "python.builtin.object")
}

#' Coerce ADMET-AI predict output to R (skip if already native)
#' @keywords internal
#' @noRd
.nts_admet_ai_to_r <- function(x) {
  if (.nts_admet_ai_is_py(x)) {
    return(reticulate::py_to_r(x))
  }
  x
}

#' Parse one ADMETModel$predict() result (ref: tests/toxPi_202506065.R)
#' @keywords internal
#' @noRd
.nts_admet_ai_row_from_predict <- function(result, smiles) {
  if (.nts_admet_ai_is_py(result) && reticulate::py_has_attr(result, "reset_index")) {
    out <- .nts_admet_ai_to_r(result$reset_index())
    if ("index" %in% names(out) && !"SMILES" %in% names(out)) {
      names(out)[names(out) == "index"] <- "SMILES"
    }
    if (!"SMILES" %in% names(out)) {
      out$SMILES <- smiles
    }
    return(out)
  }

  r <- .nts_admet_ai_to_r(result)
  if (is.data.frame(r)) {
    r$SMILES <- smiles
    return(r)
  }
  if (is.list(r)) {
    mat <- tryCatch(
      do.call(cbind, r),
      error = function(e) NULL
    )
    if (!is.null(mat)) {
      out <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
      out$SMILES <- smiles
      return(out)
    }
  }
  out <- as.data.frame(t(unlist(r)), stringsAsFactors = FALSE, check.names = FALSE)
  out$SMILES <- smiles
  out
}

#' Predict ADMET via local ADMET-AI (reticulate, optional CLI)
#'
#' Default path matches **`tests/toxPi_202506065.R`**: `reticulate::use_condaenv("admet_env")`
#' + `admet_ai$ADMETModel()$predict()`. Set **`NTSminer_ADMET_AI_PYTHON`** to the env's
#' `python.exe` if `use_condaenv()` fails (common on older conda builds).
#'
#' @param smiles_list Character vector of SMILES.
#' @param conda_env Conda environment name (default `admet_env`).
#' @param prefer_cli If `TRUE`, try `admet_predict` CLI first; on failure fall back to reticulate.
#' @param python_path Optional explicit path to `python.exe` (overrides auto-discovery).
#'
#' @return Wide `data.frame` with `SMILES` and TDC ADMET columns.
#'
#' @seealso [predict_admet_risk]
#' @export
predict_admet_ai <- function(smiles_list,
                             conda_env = Sys.getenv("NTSminer_ADMET_AI_ENV", "admet_env"),
                             prefer_cli = FALSE,
                             python_path = NULL) {
  smiles <- unique(trimws(as.character(smiles_list)))
  smiles <- smiles[!is.na(smiles) & nzchar(smiles)]
  if (!length(smiles)) {
    return(data.frame())
  }
  tbl <- data.frame(SMILES = smiles, stringsAsFactors = FALSE)

  if (isTRUE(prefer_cli) && .nts_admet_ai_cli_available()) {
    cli_out <- tryCatch({
      tmp_in <- tempfile(fileext = ".csv")
      tmp_out <- tempfile(fileext = ".csv")
      on.exit(unlink(c(tmp_in, tmp_out)), add = TRUE)
      utils::write.csv(tbl, tmp_in, row.names = FALSE)
      cmd <- sprintf(
        "admet_predict --data_path %s --save_path %s --smiles_column SMILES",
        shQuote(tmp_in, type = "cmd"),
        shQuote(tmp_out, type = "cmd")
      )
      status <- system(cmd)
      if (status != 0L || !file.exists(tmp_out)) {
        stop("admet_predict failed (exit ", status, ").", call. = FALSE)
      }
      out <- utils::read.csv(tmp_out, stringsAsFactors = FALSE, check.names = FALSE)
      if (!"SMILES" %in% names(out) && "smiles" %in% names(out)) {
        names(out)[names(out) == "smiles"] <- "SMILES"
      }
      out$admet_ai_backend <- "ADMET-AI"
      out
    }, error = function(e) {
      message("ADMET-AI CLI failed: ", conditionMessage(e), "; using reticulate.")
      NULL
    })
    if (!is.null(cli_out)) {
      return(cli_out)
    }
  }

  .nts_admet_ai_configure_python(conda_env = conda_env, python_path = python_path)
  admet_mod <- reticulate::import("admet_ai")
  model <- admet_mod$ADMETModel()
  rows <- lapply(smiles, function(s) {
    pred <- model$predict(smiles = s)
    .nts_admet_ai_row_from_predict(pred, s)
  })
  out <- do.call(rbind, rows)
  out$admet_ai_backend <- "ADMET-AI"
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.nts_merge_admet_backends <- function(lab, ai) {
  if (is.null(lab) || !nrow(lab)) {
    return(ai)
  }
  if (is.null(ai) || !nrow(ai)) {
    return(lab)
  }
  smi_lab <- if ("SMILES" %in% names(lab)) "SMILES" else if ("smiles" %in% names(lab)) "smiles" else NULL
  smi_ai <- if ("SMILES" %in% names(ai)) "SMILES" else if ("smiles" %in% names(ai)) "smiles" else NULL
  if (is.null(smi_lab) || is.null(smi_ai)) {
    return(cbind(lab, ai))
  }
  if (smi_lab != "SMILES") {
    names(lab)[names(lab) == smi_lab] <- "SMILES"
  }
  if (smi_ai != "SMILES") {
    names(ai)[names(ai) == smi_ai] <- "SMILES"
  }
  dup <- intersect(setdiff(names(lab), "SMILES"), setdiff(names(ai), "SMILES"))
  if (length(dup)) {
    names(lab)[match(dup, names(lab))] <- paste0(dup, "_admetlab3")
    names(ai)[match(dup, names(ai))] <- paste0(dup, "_admetai")
  }
  merge(lab, ai, by = "SMILES", all = TRUE, sort = FALSE)
}

#' Resolve SMILES from NTSProject for tox prediction
#' @keywords internal
#' @noRd
.nts_toxpred_smiles_from_project <- function(project) {
  ann <- get_ms1_annotations(project)
  ft <- get_feature_table(project)
  if (!nrow(ann) || !"name" %in% names(ann)) {
    stop("ms1_annotations is empty or lacks `name`.", call. = FALSE)
  }
  smi_cols <- c(
    "master_smiles", "SMILES", "smiles", "library_SMILES",
    "inferred_smiles"
  )
  smi <- rep(NA_character_, nrow(ann))
  for (cl in smi_cols) {
    if (!cl %in% names(ann)) {
      next
    }
    v <- trimws(as.character(ann[[cl]]))
    miss <- is.na(smi) | !nzchar(smi)
    smi[miss] <- v[miss]
  }
  out <- data.frame(
    name = as.character(ann$name),
    SMILES = smi,
    stringsAsFactors = FALSE
  )
  if (nrow(ft) && "name" %in% names(ft)) {
    out <- merge(
      unique(ft[, intersect(names(ft), c("name", "mz", "rt")), drop = FALSE]),
      out,
      by = "name",
      all.x = TRUE,
      sort = FALSE
    )
  }
  out <- out[!is.na(out$SMILES) & nzchar(out$SMILES), , drop = FALSE]
  out[!duplicated(out$SMILES), , drop = FALSE]
}

#' Merge wide tox table onto project `@physchem_proxy` by SMILES / name
#' @keywords internal
#' @noRd
.nts_toxpred_write_physchem <- function(project, tox_wide, id_map) {
  if (!nrow(tox_wide)) {
    return(project)
  }
  smi_col <- if ("SMILES" %in% names(tox_wide)) "SMILES" else "smiles"
  proxy <- as.data.frame(get_physchem_proxy(project), stringsAsFactors = FALSE)
  if (!nrow(proxy) || !"name" %in% names(proxy)) {
    proxy <- unique(id_map[, c("name", "SMILES"), drop = FALSE])
  }
  idx <- match(trimws(as.character(proxy$name)), trimws(as.character(id_map$name)))
  proxy_smiles <- id_map$SMILES[idx]
  j <- match(trimws(as.character(tox_wide[[smi_col]])), trimws(as.character(proxy_smiles)))
  bring <- setdiff(names(tox_wide), smi_col)
  for (cl in bring) {
    proxy[[cl]] <- tox_wide[[cl]][j]
  }
  if ("LogKow_est" %in% names(proxy)) {
    proxy$LogP_pred <- proxy$LogKow_est
  }
  project <- set_physchem_proxy(project, proxy)
  project
}
