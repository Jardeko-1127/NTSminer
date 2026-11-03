# ClassyFire (ChemOnt) structure classification via SMILES
# API: http://classyfire.wishartlab.com/

#' Default ClassyFire API base URL
#' @keywords internal
#' @noRd
.nts_classyfire_default_base <- function() {
  url <- Sys.getenv("NTSminer_CLASSYFIRE_BASE_URL", unset = "http://classyfire.wishartlab.com")
  if (!nzchar(url)) {
    url <- "http://classyfire.wishartlab.com"
  }
  sub("/$", "", url)
}

#' Use system HTTP(S)_PROXY for ClassyFire (default: direct, no proxy)
#' @keywords internal
#' @noRd
.nts_classyfire_use_system_proxy <- function() {
  identical(toupper(Sys.getenv("NTSminer_CLASSYFIRE_USE_PROXY", unset = "FALSE")), "TRUE")
}

#' httr config: bypass broken proxies unless NTSminer_CLASSYFIRE_USE_PROXY=TRUE
#' @keywords internal
#' @noRd
.nts_classyfire_httr_config <- function() {
  if (isTRUE(.nts_classyfire_use_system_proxy())) {
    return(httr::config())
  }
  httr::config(proxy = "")
}

#' @keywords internal
#' @noRd
.nts_classyfire_do_get <- function(url, timeout_sec = 60) {
  httr::GET(
    url = url,
    httr::accept_json(),
    httr::timeout(timeout_sec),
    .nts_classyfire_httr_config()
  )
}

#' @keywords internal
#' @noRd
.nts_classyfire_do_post <- function(url, body_list, timeout_sec = 90) {
  body_raw <- jsonlite::toJSON(body_list, auto_unbox = TRUE, null = "null")
  httr::POST(
    url = url,
    body = body_raw,
    httr::content_type_json(),
    httr::accept_json(),
    httr::timeout(timeout_sec),
    .nts_classyfire_httr_config()
  )
}

#' @keywords internal
#' @noRd
.nts_classyfire_submit_failed <- function(last_error = NULL,
                                        http_status = NA_integer_,
                                        response_body = NULL) {
  structure(
    list(
      ok = FALSE,
      last_error = last_error,
      http_status = http_status,
      response_body = response_body
    ),
    class = "nts_classyfire_submit_failed"
  )
}

#' @keywords internal
#' @noRd
.nts_classyfire_stop_submit_failed <- function(fail = NULL,
                                              last_error = NULL,
                                              http_status = NA_integer_,
                                              response_body = NULL) {
  if (inherits(fail, "nts_classyfire_submit_failed")) {
    last_error <- fail$last_error
    http_status <- fail$http_status
    response_body <- fail$response_body
  }
  err <- if (is.null(last_error)) "" else conditionMessage(last_error)
  proxy_env <- paste(
    c(
      Sys.getenv("http_proxy", unset = ""),
      Sys.getenv("https_proxy", unset = ""),
      Sys.getenv("all_proxy", unset = ""),
      Sys.getenv("HTTP_PROXY", unset = ""),
      Sys.getenv("HTTPS_PROXY", unset = "")
    ),
    collapse = " "
  )
  msg <- "ClassyFire submit failed (no query id)."
  if (is.finite(http_status)) {
    msg <- paste0(msg, "\n  HTTP status: ", http_status)
  }
  if (!is.null(response_body) && nzchar(response_body)) {
    rb <- substr(as.character(response_body), 1L, 400L)
    msg <- paste0(msg, "\n  Response: ", rb)
  }
  if (is.finite(http_status) && http_status >= 500L) {
    msg <- paste0(
      msg,
      "\n  Server error (5xx): retry later, or use prefer_inchikey_lookup = TRUE ",
      "to fetch existing classifications by InChIKey without POST."
    )
  }
  if (nzchar(err)) {
    msg <- paste0(msg, "\n  Last error: ", err)
  }
  if (nzchar(proxy_env) || grepl("proxy", err, ignore.case = TRUE)) {
    msg <- paste0(
      msg,
      "\n  A system HTTP proxy is set but ClassyFire requests use a direct connection by default.",
      "\n  In R, try:\n",
      "    Sys.unsetenv(c('http_proxy','https_proxy','all_proxy','HTTP_PROXY','HTTPS_PROXY'))\n",
      "    cf <- classify_classyfire(...)\n",
      "  Or set NTSminer_CLASSYFIRE_BASE_URL if you need a mirror (must speak ClassyFire API).",
      "\n  To force the system proxy: Sys.setenv(NTSminer_CLASSYFIRE_USE_PROXY = 'TRUE')"
    )
  }
  stop(msg, call. = FALSE)
}

#' @keywords internal
#' @noRd
.nts_classyfire_level_name <- function(node) {
  if (is.null(node) || !length(node)) {
    return(NA_character_)
  }
  if (is.list(node) && !is.null(node$name) && nzchar(as.character(node$name)[1L])) {
    return(as.character(node$name)[1L])
  }
  NA_character_
}

#' Flatten one ClassyFire entity JSON object to a data.frame row
#' @keywords internal
#' @noRd
.nts_classyfire_entity_to_row <- function(ent) {
  if (!is.list(ent) || !length(ent)) {
    return(NULL)
  }
  data.frame(
    identifier = if (!is.null(ent$identifier)) as.character(ent$identifier)[1L] else NA_character_,
    SMILES = if (!is.null(ent$smiles)) as.character(ent$smiles)[1L] else NA_character_,
    InChIKey = if (!is.null(ent$inchikey)) as.character(ent$inchikey)[1L] else NA_character_,
    Kingdom = .nts_classyfire_level_name(ent$kingdom),
    Superclass = .nts_classyfire_level_name(ent$superclass),
    Class = .nts_classyfire_level_name(ent$class),
    Subclass = .nts_classyfire_level_name(ent$subclass),
    Direct_parent = .nts_classyfire_level_name(ent$direct_parent),
    classification_version = if (!is.null(ent$classification_version)) {
      as.character(ent$classification_version)[1L]
    } else {
      NA_character_
    },
    classyfire_status = "classified",
    cf_status = "found",
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.nts_classyfire_lookup_row_not_found <- function(identifier,
                                                 smiles = NA_character_,
                                                 inchikey = NA_character_) {
  data.frame(
    identifier = identifier,
    SMILES = smiles,
    InChIKey = inchikey,
    Kingdom = NA_character_,
    Superclass = NA_character_,
    Class = NA_character_,
    Subclass = NA_character_,
    Direct_parent = NA_character_,
    classification_version = NA_character_,
    classyfire_status = "not_found",
    cf_status = "not_found",
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.nts_classyfire_parse_query_json <- function(parsed, input_ids = NULL) {
  .nts_require_httr_json("classyfire")
  rows <- list()
  if (!is.null(parsed$entities) && length(parsed$entities)) {
    for (ent in parsed$entities) {
      row <- .nts_classyfire_entity_to_row(ent)
      if (!is.null(row)) {
        rows[[length(rows) + 1L]] <- row
      }
    }
  }
  if (length(rows)) {
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
  } else {
    out <- data.frame(
      identifier = character(),
      SMILES = character(),
      InChIKey = character(),
      Kingdom = character(),
      Superclass = character(),
      Class = character(),
      Subclass = character(),
      Direct_parent = character(),
      classification_version = character(),
      classyfire_status = character(),
      cf_status = character(),
      stringsAsFactors = FALSE
    )
  }
  if (!is.null(input_ids) && length(input_ids)) {
    miss <- setdiff(input_ids, out$identifier)
    if (length(miss)) {
      extra <- lapply(miss, function(id) {
        data.frame(
          identifier = id,
          SMILES = NA_character_,
          InChIKey = NA_character_,
          Kingdom = NA_character_,
          Superclass = NA_character_,
          Class = NA_character_,
          Subclass = NA_character_,
          Direct_parent = NA_character_,
          classification_version = NA_character_,
          classyfire_status = "unclassified",
          cf_status = "not_found",
          stringsAsFactors = FALSE
        )
      })
      out <- rbind(out, do.call(rbind, extra))
      rownames(out) <- NULL
    }
  }
  out
}

#' Build tab-separated query_input for ClassyFire (id \\t SMILES per line)
#' @keywords internal
#' @noRd
.nts_classyfire_query_input <- function(smiles, ids = NULL) {
  smiles <- trimws(as.character(smiles))
  ok <- !is.na(smiles) & nzchar(smiles)
  smiles <- smiles[ok]
  if (!length(smiles)) {
    stop("No valid SMILES to submit to ClassyFire.", call. = FALSE)
  }
  if (is.null(ids)) {
    ids <- paste0("CMP_", seq_along(smiles))
  } else {
    ids <- as.character(ids)[ok]
    if (length(ids) != length(smiles)) {
      stop("`compound_ids` length must match non-empty SMILES.", call. = FALSE)
    }
  }
  paste(paste(ids, smiles, sep = "\t"), collapse = "\n")
}

#' POST a structure query to ClassyFire
#' @keywords internal
#' @noRd
.nts_classyfire_submit <- function(smiles,
                                  ids = NULL,
                                  label = "ntsminer",
                                  base_url = .nts_classyfire_default_base(),
                                  timeout_sec = 90) {
  .nts_require_httr_json("classyfire")
  query_input <- .nts_classyfire_query_input(smiles, ids = ids)
  body <- list(
    label = as.character(label)[1L],
    query_input = query_input,
    query_type = "STRUCTURE"
  )
  base <- sub("/$", "", base_url)
  url <- paste0(base, "/queries")
  last_err <- NULL
  resp <- tryCatch(
    .nts_classyfire_do_post(url, body_list = body, timeout_sec = timeout_sec),
    error = function(e) {
      last_err <<- conditionMessage(e)
      message("ClassyFire POST error: ", last_err)
      NULL
    }
  )
  if (is.null(resp) || httr::status_code(resp) != 200L) {
    code <- if (!is.null(resp)) httr::status_code(resp) else NA_integer_
    resp_body <- if (!is.null(resp)) {
      tryCatch(httr::content(resp, as = "text", encoding = "UTF-8"), error = function(e) NULL)
    } else {
      NULL
    }
    message("ClassyFire submit HTTP ", code)
    return(.nts_classyfire_submit_failed(
      last_error = last_err,
      http_status = code,
      response_body = resp_body
    ))
  }
  httr::content(resp, as = "parsed", type = "application/json", encoding = "UTF-8")
}

#' Submit one SMILES at a time (fallback when batch POST returns 5xx)
#' @keywords internal
#' @noRd
.nts_classyfire_submit_one_by_one <- function(smiles,
                                              ids,
                                              label = "ntsminer",
                                              base_url = .nts_classyfire_default_base(),
                                              timeout_sec = 90,
                                              progress = TRUE) {
  rows <- vector("list", length(smiles))
  for (i in seq_along(smiles)) {
    if (isTRUE(progress)) {
      message("[ClassyFire] single submit ", i, "/", length(smiles))
    }
    post <- .nts_classyfire_submit(
      smiles[i],
      ids = ids[i],
      label = paste0(label, "_", i),
      base_url = base_url,
      timeout_sec = timeout_sec
    )
    if (inherits(post, "nts_classyfire_submit_failed") || is.null(post$id)) {
      next
    }
    parsed <- .nts_classyfire_wait_query(
      post$id,
      base_url = base_url,
      timeout_sec = timeout_sec,
      progress = FALSE
    )
    part <- .nts_classyfire_parse_query_json(parsed, input_ids = ids[i])
    if (nrow(part)) {
      rows[[i]] <- part
    }
    Sys.sleep(6)
  }
  if (!length(rows)) {
    return(NULL)
  }
  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (!length(rows)) {
    return(NULL)
  }
  do.call(rbind, rows)
}

#' GET /entities for one structure; retry with canonical SMILES if first key misses
#' @keywords internal
#' @noRd
.nts_classyfire_entity_lookup_one <- function(smiles,
                                             inchi_key,
                                             base_url,
                                             progress = TRUE) {
  key <- trimws(as.character(inchi_key)[1L])
  ent <- if (nzchar(key)) {
    .nts_classyfire_entity_by_inchikey(key, base_url = base_url)
  } else {
    NULL
  }
  if (!is.null(ent)) {
    return(list(entity = ent, inchikey = key, smiles_used = smiles))
  }
  smi <- trimws(as.character(smiles)[1L])
  if (is.na(smi) || !nzchar(smi)) {
    return(NULL)
  }
  canon <- tryCatch(
    canonicalize_smiles(smi),
    error = function(e) NA_character_
  )
  if (is.na(canon) || !nzchar(canon) || identical(canon, smi)) {
    return(NULL)
  }
  key2 <- tryCatch(
    .nts_toxpi_inchikeys(canon)[1L],
    error = function(e) NA_character_
  )
  if (is.na(key2) || !nzchar(key2) || identical(key2, key)) {
    return(NULL)
  }
  if (isTRUE(progress)) {
    message(
      "[ClassyFire] retry canonical SMILES -> InChIKey ", key2
    )
  }
  ent2 <- .nts_classyfire_entity_by_inchikey(key2, base_url = base_url)
  if (is.null(ent2)) {
    return(NULL)
  }
  list(entity = ent2, inchikey = key2, smiles_used = canon)
}

#' Try ClassyFire /entities by InChIKey (rcdk) before POST submit
#' @keywords internal
#' @noRd
.nts_classyfire_lookup_by_inchikey <- function(smiles, ids, base_url, progress = TRUE) {
  n <- length(smiles)
  if (isTRUE(progress)) {
    message("[ClassyFire] InChIKey lookup called, n=", n)
  }
  if (!requireNamespace("rcdk", quietly = TRUE)) {
    if (isTRUE(progress)) {
      message("[ClassyFire] InChIKey lookup skipped: package 'rcdk' not available")
    }
    return(data.frame())
  }
  keys <- tryCatch(
    .nts_toxpi_inchikeys(smiles),
    error = function(e) {
      if (isTRUE(progress)) {
        message("[ClassyFire] InChIKey generation failed: ", conditionMessage(e))
      }
      rep(NA_character_, n)
    }
  )
  if (length(keys) != n) {
    keys <- rep(NA_character_, n)
  }
  if (isTRUE(progress)) {
    message("[ClassyFire] InChIKeys: ", paste(keys, collapse = ", "))
  }
  rows <- vector("list", n)
  for (i in seq_len(n)) {
    if (is.na(keys[i]) || !nzchar(keys[i])) {
      rows[[i]] <- .nts_classyfire_lookup_row_not_found(ids[i], smiles[i], NA_character_)
      next
    }
    hit1 <- .nts_classyfire_entity_lookup_one(
      smiles[i],
      keys[i],
      base_url = base_url,
      progress = FALSE
    )
    if (isTRUE(progress)) {
      message(
        "[ClassyFire] GET /entities/", keys[i],
        ": ", if (is.null(hit1)) "not found" else "ok"
      )
    }
    if (is.null(hit1)) {
      rows[[i]] <- .nts_classyfire_lookup_row_not_found(ids[i], smiles[i], keys[i])
      next
    }
    row <- .nts_classyfire_entity_to_row(hit1$entity)
    if (is.null(row)) {
      rows[[i]] <- .nts_classyfire_lookup_row_not_found(ids[i], smiles[i], hit1$inchikey)
      next
    }
    row$identifier <- ids[i]
    row$SMILES <- smiles[i]
    row$InChIKey <- hit1$inchikey
    row$cf_status <- "found"
    rows[[i]] <- row
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Poll ClassyFire until query completes
#' @keywords internal
#' @noRd
.nts_classyfire_wait_query <- function(query_id,
                                      base_url = .nts_classyfire_default_base(),
                                      timeout_sec = 300,
                                      poll_sec = 2,
                                      progress = TRUE) {
  .nts_require_httr_json("classyfire")
  if (!is.finite(query_id)) {
    stop("Invalid ClassyFire query id.", call. = FALSE)
  }
  url <- paste0(sub("/$", "", base_url), "/queries/", query_id, ".json")
  deadline <- Sys.time() + as.numeric(timeout_sec)
  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    resp <- tryCatch(
      .nts_classyfire_do_get(url, timeout_sec = 60),
      error = function(e) NULL
    )
    if (is.null(resp) || httr::status_code(resp) != 200L) {
      if (Sys.time() > deadline) {
        stop("ClassyFire poll timeout (HTTP error).", call. = FALSE)
      }
      Sys.sleep(poll_sec)
      next
    }
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    parsed <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
    st_raw <- parsed$status
    st <- tolower(as.character(if (is.null(st_raw)) "" else st_raw))
    ent <- parsed$entities
    n_ent <- if (is.null(ent)) 0L else length(ent)
    if (isTRUE(progress) && attempt %% 5L == 1L) {
      message("[ClassyFire] query ", query_id, " status=", parsed$status,
              " entities=", n_ent)
    }
    if (n_ent > 0L || st %in% c("done", "complete", "finished")) {
      return(parsed)
    }
    if (st %in% c("failed", "error")) {
      stop("ClassyFire query failed (status: ", parsed$status, ").", call. = FALSE)
    }
    if (Sys.time() > deadline) {
      stop(
        "ClassyFire poll timeout after ", timeout_sec, " s (last status: ",
        parsed$status, ").",
        call. = FALSE
      )
    }
    Sys.sleep(poll_sec)
  }
}

#' GET existing ClassyFire entity by InChIKey
#' @keywords internal
#' @noRd
.nts_classyfire_entity_by_inchikey <- function(inchi_key,
                                                base_url = .nts_classyfire_default_base(),
                                                timeout_sec = 60) {
  .nts_require_httr_json("classyfire")
  key <- gsub("^InChIKey=", "", trimws(as.character(inchi_key)[1L]))
  if (!nzchar(key)) {
    return(NULL)
  }
  url <- paste0(sub("/$", "", base_url), "/entities/", key, ".json")
  resp <- tryCatch(
    .nts_classyfire_do_get(url, timeout_sec = timeout_sec),
    error = function(e) NULL
  )
  if (is.null(resp) || httr::status_code(resp) != 200L) {
    return(NULL)
  }
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  if (!nzchar(trimws(txt)) || identical(trimws(txt), "{}")) {
    return(NULL)
  }
  parsed <- jsonlite::fromJSON(txt, simplifyVector = FALSE)
  if (!is.list(parsed) || !length(parsed)) {
    return(NULL)
  }
  parsed
}

#' Classify compounds with ClassyFire (ChemOnt) from SMILES
#'
#' Submits structures to the public [ClassyFire](http://classyfire.wishartlab.com/)
#' API and returns ChemOnt **Kingdom / Superclass / Class / Subclass** and
#' **direct parent** in a wide `data.frame`. Use [canonicalize_smiles()] first
#' when SMILES may differ in aromaticity or salt form.
#'
#' @param smiles_list Character vector of SMILES, or a `data.frame` with a
#'   `SMILES` / `smiles` column.
#' @param compound_ids Optional identifiers (same order as SMILES); default
#'   `CMP_001`, ... Tab-separated with SMILES for the API.
#' @param label Query label sent to ClassyFire (default `"ntsminer"`).
#' @param canonicalize If `TRUE`, run [canonicalize_smiles()] before submit.
#' @param base_url API root (default Wishart lab server).
#' @param timeout_sec Max seconds to wait for query completion.
#' @param poll_sec Seconds between status polls.
#' @param progress If `TRUE`, show poll messages.
#' @param prefer_inchikey_lookup If `TRUE` (default), try **GET** `/entities/{InChIKey}`
#'   via **rcdk** before POST (no rate limit on cache hits).
#' @param fallback_single If batch POST fails with HTTP 5xx, submit one SMILES per query
#'   (slower; respects ~12 POST/min API limit).
#'
#' @section Network / proxy:
#' By default, ClassyFire HTTP calls **do not** use `http_proxy` / `https_proxy`
#' (avoids broken local proxy DNS). Set `Sys.setenv(NTSminer_CLASSYFIRE_USE_PROXY = "TRUE")`
#' to use the system proxy. Optional `NTSminer_CLASSYFIRE_BASE_URL` overrides the API host.
#'
#' @return `data.frame` with columns `identifier`, `SMILES`, `InChIKey`,
#'   `Kingdom`, `Superclass`, `Class`, `Subclass`, `Direct_parent`,
#'   `classification_version`, `classyfire_status` (`classified` / `unclassified` / `not_found`),
#'   `cf_status` (`found` / `not_found` for InChIKey cache hits).
#'
#' @seealso [get_classyfire_classification], [canonicalize_smiles]
#' @export
#'
#' @examples
#' \dontrun{
#' classify_classyfire(c("CCO", "CC"))
#' }
classify_classyfire <- function(smiles_list,
                               compound_ids = NULL,
                               label = "ntsminer",
                               canonicalize = TRUE,
                               base_url = .nts_classyfire_default_base(),
                               timeout_sec = 300,
                               poll_sec = 2,
                               progress = TRUE,
                               prefer_inchikey_lookup = TRUE,
                               fallback_single = TRUE) {
  if (is.data.frame(smiles_list)) {
    df <- as.data.frame(smiles_list, stringsAsFactors = FALSE)
    smi_col <- intersect(c("SMILES", "smiles", "master_smiles"), names(df))[1L]
    if (is.na(smi_col)) {
      stop("smiles_list data.frame needs a SMILES column.", call. = FALSE)
    }
    smiles <- df[[smi_col]]
    if (is.null(compound_ids)) {
      for (id_col in c("name", "Name", "Lab_ID", "identifier", "compound_id")) {
        if (id_col %in% names(df)) {
          compound_ids <- df[[id_col]]
          break
        }
      }
    }
  } else {
    nm <- names(smiles_list)
    smiles <- trimws(as.character(smiles_list))
    if (is.null(compound_ids) && !is.null(nm) && any(nzchar(nm))) {
      compound_ids <- nm
    }
  }
  ok <- !is.na(smiles) & nzchar(smiles)
  if (!any(ok)) {
    return(data.frame())
  }
  smiles <- smiles[ok]
  if (!is.null(compound_ids)) {
    compound_ids <- as.character(compound_ids)[ok]
  }

  if (isTRUE(canonicalize)) {
    smiles <- canonicalize_smiles(smiles)
  }

  ids <- compound_ids
  if (is.null(ids)) {
    ids <- paste0("CMP_", seq_along(smiles))
  }

  out <- NULL
  pending <- seq_along(smiles)

  if (isTRUE(prefer_inchikey_lookup)) {
    hit <- tryCatch(
      .nts_classyfire_lookup_by_inchikey(
        smiles,
        ids,
        base_url = base_url,
        progress = progress
      ),
      error = function(e) {
        warning(
          "[ClassyFire] InChIKey lookup error: ",
          conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
    if (!is.null(hit) && nrow(hit) > 0L) {
      out <- hit
      found_ids <- hit$identifier[!is.na(hit$cf_status) & hit$cf_status == "found"]
      pending <- which(!(ids %in% found_ids))
      if (isTRUE(progress)) {
        message(
          "[ClassyFire] InChIKey lookup: ",
          length(found_ids), " / ", length(smiles), " found"
        )
      }
    } else if (isTRUE(progress)) {
      message("[ClassyFire] InChIKey lookup: 0 / ", length(smiles), " found")
    }
  }

  if (length(pending)) {
    sm_p <- smiles[pending]
    id_p <- ids[pending]
    if (isTRUE(progress)) {
      message("[ClassyFire] POST submit ", length(sm_p), " structure(s)...")
    }

    post <- .nts_classyfire_submit(
      sm_p,
      ids = id_p,
      label = label,
      base_url = base_url
    )

    post_out <- NULL
    if (inherits(post, "nts_classyfire_submit_failed") || is.null(post$id)) {
      if (isTRUE(fallback_single)) {
        if (isTRUE(progress)) {
          message("[ClassyFire] batch POST failed; trying one compound per query...")
        }
        post_out <- .nts_classyfire_submit_one_by_one(
          sm_p,
          ids = id_p,
          label = label,
          base_url = base_url,
          timeout_sec = timeout_sec,
          progress = progress
        )
        if (is.null(post_out) || !nrow(post_out)) {
          .nts_classyfire_stop_submit_failed(
            fail = if (inherits(post, "nts_classyfire_submit_failed")) post else NULL
          )
        }
      } else {
        warning(
          "ClassyFire POST failed and fallback_single = FALSE. ",
          "Returning partial results (InChIKey lookup hits only).",
          call. = FALSE
        )
      }
    } else {
      parsed <- .nts_classyfire_wait_query(
        query_id = post$id,
        base_url = base_url,
        timeout_sec = timeout_sec,
        poll_sec = poll_sec,
        progress = progress
      )
      post_out <- .nts_classyfire_parse_query_json(parsed, input_ids = id_p)
    }

    if (!is.null(post_out) && nrow(post_out)) {
      out <- if (is.null(out) || !nrow(out)) post_out else rbind(out, post_out)
      rownames(out) <- NULL
    }
  }

  base <- .nts_classyfire_parse_query_json(list(entities = list()), input_ids = ids)
  base$SMILES <- smiles
  if (!is.null(out) && nrow(out)) {
    for (i in seq_len(nrow(out))) {
      j <- match(out$identifier[i], base$identifier)
      if (!is.na(j)) {
        base[j, names(out)] <- out[i, ]
        base$SMILES[j] <- ifelse(
          is.na(out$SMILES[i]) || !nzchar(out$SMILES[i]),
          smiles[j],
          out$SMILES[i]
        )
      }
    }
  }
  base[match(ids, base$identifier), , drop = FALSE]
}

#' Fetch existing ClassyFire classification by InChIKey
#'
#' Uses the `/entities/{InChIKey}.json` endpoint (no new structure submission).
#' Useful when structures were classified previously and indexed by InChIKey.
#'
#' @param inchi_keys Character vector of InChIKeys (with or without `InChIKey=` prefix).
#' @param base_url API root.
#' @param progress If `TRUE`, message per key.
#'
#' @return Same column layout as [classify_classyfire()]; failed keys get
#'   `classyfire_status = "not_found"`.
#' @export
#'
#' @seealso [classify_classyfire]
get_classyfire_classification <- function(inchi_keys,
                                            base_url = .nts_classyfire_default_base(),
                                            progress = TRUE) {
  keys <- unique(trimws(as.character(inchi_keys)))
  keys <- keys[!is.na(keys) & nzchar(keys)]
  if (!length(keys)) {
    return(data.frame())
  }
  rows <- vector("list", length(keys))
  for (i in seq_along(keys)) {
    if (isTRUE(progress)) {
      message("[ClassyFire] entity ", i, "/", length(keys))
    }
    ent <- .nts_classyfire_entity_by_inchikey(keys[i], base_url = base_url)
    if (!is.null(ent)) {
      row <- .nts_classyfire_entity_to_row(ent)
      if (!is.null(row)) {
        row$identifier <- keys[i]
        rows[[i]] <- row
        next
      }
    }
    rows[[i]] <- .nts_classyfire_lookup_row_not_found(
      identifier = keys[i],
      smiles = NA_character_,
      inchikey = keys[i]
    )
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
