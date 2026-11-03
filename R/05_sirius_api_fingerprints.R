#' Extract CSI:FingerID fingerprints via SIRIUS REST API (fallback)
#'
#' Connects to a running SIRIUS 6.x REST API service, locates previously
#' computed project files (created by \code{compute_sirius_fingerprints()}),
#' and retrieves the numeric fingerprint vectors for each feature's top-ranked
#' formula candidate.  The resulting vectors are written into the
#' \code{fingerprint} list-column of \verb{@sirius_fingerprints}.
#'
#' @details
#' Since v0.1.12, \code{compute_sirius_fingerprints()} attempts to export
#' fingerprint vectors directly via CLI (\code{export --type=fingerprints}).
#' If that succeeds, this function is \emph{not} needed.  Use this function
#' as a fallback when CLI export is unsupported by the SIRIUS version and
#' fingerprint columns are empty after \code{compute_sirius_fingerprints()}.
#'
#' This function automates the REST API process:
#'
#' \enumerate{
#'   \item Locates \code{.sirius} project directories (via
#'     \code{sirius_projects_dir}, \code{output_dir}, or by probing known
#'     temporary locations).
#'   \item Stops any existing SIRIUS background service, then starts a fresh
#'     headless service with a workspace that contains all discovered
#'     \code{.sirius} projects (or re-uses one already running at
#'     \code{api_base}).
#'   \item For each visible project, queries aligned features and their
#'     top-ranked formula candidates (with \code{optFields = "predictedFingerprint"}).
#'   \item Maps each feature name to its fingerprint vector and populates
#'     the \code{fingerprint} list-column.
#' }
#'
#' @param project An \code{\link{NTSProject}} whose \verb{@sirius_fingerprints}
#'   slot contains the formula-metadata rows produced by
#'   \code{compute_sirius_fingerprints()}.
#' @param sirius_projects_dir Directory containing one or more \code{.sirius}
#'   project subdirectories.  If \code{NULL} (default) and \code{output_dir}
#'   is given, the function searches there; otherwise it probes common
#'   temporary locations.
#' @param output_dir The same \code{output_dir} that was passed to
#'   \code{compute_sirius_fingerprints()}.  If non-\code{NULL}, the function
#'   will search recursively for \code{.sirius} directories inside it.
#' @param sirius_exe Path to the SIRIUS executable. Required if
#'   \code{start_service = TRUE}.  Default \code{NULL} reads environment
#'   variable \code{NTSminer_SIRIUS_EXE} / \code{SIRIUS_EXE}.
#' @param port TCP port for the headless service.  If \code{NULL} (default)
#'   the function picks a random free port automatically, avoiding conflicts
#'   from \code{TIME_WAIT} on Windows.
#' @param api_base URL of an already-running SIRIUS REST API service
#'   (e.g. \code{"http://localhost:52423"}).  If \code{NULL}, a new service
#'   is started automatically.
#' @param start_service If \code{TRUE} (default when \code{api_base = NULL}),
#'   the function starts its own temporary service. Set to \code{FALSE} to
#'   use an existing service via \code{api_base}.
#' @param keep_service If \code{TRUE}, the temporary service is \emph{not}
#'   stopped after extraction. Only relevant when \code{start_service = TRUE}.
#' @param top_n_rank Only extract fingerprints for formula candidates whose
#'   \code{rank <= top_n_rank} (default \code{1}, i.e. the top hit only).
#' @param timeout_sec Max seconds to wait for the service to become ready
#'   (default \code{180}).  On Windows the first headless start may be slow
#'   because SIRIUS loads all \code{.sirius} projects and initializes its
#'   embedded Web server.
#'
#' @return The updated \code{project} with fingerprint vectors populated in
#'   \verb{@sirius_fingerprints}.
#'
#' @section Prerequisites:
#' \itemize{
#'   \item \code{compute_sirius_fingerprints()} must have been called with
#'     \code{keep_intermediates = TRUE} so that the per-batch \code{.sirius}
#'     project directories still exist on disk.
#'   \item The SIRIUS subscription must have CSI:FingerID / Fingerprint
#'     prediction enabled (a valid login is required).
#'   \item Packages \code{httr} and \code{jsonlite} must be installed.
#' }
#'
#' @examplesIf requireNamespace("httr", quietly = TRUE) && requireNamespace("jsonlite", quietly = TRUE)
#' \dontrun{
#' p <- compute_sirius_fingerprints(p, keep_intermediates = TRUE)
#' p <- extract_sirius_fingerprints_api(p,
#'   output_dir = "./tests/sirius_test",
#'   sirius_exe = "/opt/sirius/bin/sirius")
#' fp <- get_sirius_fingerprints(p)
#' str(fp$fingerprint[1:3])   # now contains numeric vectors
#' }
#' @export
extract_sirius_fingerprints_api <- function(project,
                                            sirius_projects_dir = NULL,
                                            output_dir = NULL,
                                            sirius_exe = NULL,
                                            port = NULL,
                                            api_base = NULL,
                                            start_service = is.null(api_base),
                                            keep_service = FALSE,
                                            top_n_rank = 1L,
                                            timeout_sec = 180) {

  # --- Dependency check ---
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required. Please install it.", call. = FALSE)
  }

  # local NULL-coalesce
  `%or%` <- function(x, y) if (is.null(x)) y else x

  # --- Project check ---
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  fp_table <- methods::slot(project, "sirius_fingerprints")
  if (is.null(fp_table) || nrow(fp_table) == 0L) {
    stop("No sirius_fingerprints found. Run compute_sirius_fingerprints() first.", call. = FALSE)
  }

  has_fp <- vapply(fp_table$fingerprint, function(x) length(x) > 0L, logical(1L))
  n_missing <- sum(!has_fp)
  message(">> ", n_missing, " / ", nrow(fp_table), " rows have missing fingerprint vectors.")
  if (!n_missing) {
    message(">> All fingerprint vectors already present. Nothing to do.")
    return(invisible(project))
  }

  # --- Resolve SIRIUS executable ---
  if (start_service) {
    if (is.null(sirius_exe) || !nzchar(sirius_exe)) {
      sirius_exe <- Sys.getenv("NTSminer_SIRIUS_EXE", Sys.getenv("SIRIUS_EXE", ""))
    }
    if (!nzchar(sirius_exe)) {
      stop(
        "SIRIUS path not set.\n",
        "  - If SIRIUS is already running (e.g. GUI is open), use:\n",
        "      extract_sirius_fingerprints_api(..., api_base = \"http://localhost:52423\", start_service = FALSE)\n",
        "  - Otherwise, pass `sirius_exe =` or set env NTSminer_SIRIUS_EXE (or SIRIUS_EXE).",
        call. = FALSE)
    }
    if (!file.exists(sirius_exe)) {
      stop("SIRIUS executable not found at: ", sirius_exe, call. = FALSE)
    }
  }

  # --- Set up API connection ---
  if (start_service) {
    # --- Pre-flight: verify SIRIUS login ---
    login_check <- tryCatch(
      system2(sirius_exe, c("login", "--show"),
              stdout = TRUE, stderr = FALSE, timeout = 30),
      error = function(e) NULL
    )
    if (is.null(login_check) || !any(grepl("Logged in as", login_check))) {
      warning("SIRIUS does not appear to be logged in. ",
              "Run '", sirius_exe, " login' first, or the service may fail.")
    } else {
      message(">> SIRIUS login confirmed.")
    }

    # Find all .sirius project files
    output_dir_abs <- if (!is.null(output_dir) && dir.exists(output_dir)) {
      normalizePath(output_dir, mustWork = TRUE)
    } else output_dir
    sirius_dirs <- .find_sirius_project_dirs(sirius_projects_dir, output_dir_abs)
    if (!length(sirius_dirs)) {
      stop(
        "No .sirius project directories found.\n",
        "  - Pass sirius_projects_dir = <dir> containing batch_*.sirius folders\n",
        "  - Or pass output_dir = <same path as compute_sirius_fingerprints output_dir>\n",
        "  - Make sure keep_intermediates = TRUE was used in compute_sirius_fingerprints.",
        call. = FALSE
      )
    }
    message(">> Found ", length(sirius_dirs), " .sirius project files.")

    # Resolve port
    if (is.null(port) || is.na(port) || port <= 0L) {
      port <- .find_free_port()
    }
    message(">> Using port ", port, " for headless SIRIUS service.")

    # Register cleanup
    if (!keep_service) {
      on.exit({
        .kill_sirius_service(port)
      }, add = TRUE)
    }

    # Build feature map by iterating over each project file
    # SIRIUS 6.x only supports one --project argument per service instance
    feature_map <- list()
    proj_count <- 0L

    for (j in seq_along(sirius_dirs)) {
      proj_file <- sirius_dirs[j]
      proj_name <- gsub("\\.sirius$", "", basename(proj_file))
      message(sprintf(">> [%d/%d] Opening project: %s", j, length(sirius_dirs), proj_name))

      # Kill any previous service
      .kill_sirius_service(port)

      # Start service with this single project
      api_url <- .start_sirius_service(sirius_exe, dirname(proj_file),
                                        port, timeout_sec,
                                        project_files = proj_file)
      if (is.null(api_url)) {
        warning("Failed to start service for project ", proj_name, call. = FALSE)
        next
      }

      # Get aligned features for this project
      af_path <- paste0("/api/projects/", .url_enc(proj_name), "/aligned-features")
      af_list <- .api_get_json(api_url, af_path)
      if (is.null(af_list) || !length(af_list)) {
        warning("No aligned features in project ", proj_name, call. = FALSE)
        next
      }

      for (af in af_list) {
        fname <- trimws(af$name)
        if (is.null(fname) || !nzchar(fname)) next
        if (is.null(feature_map[[fname]])) {
          feature_map[[fname]] <- list(
            projectId = proj_name,
            alignedFeatureId = af$alignedFeatureId
          )
        }
      }
      proj_count <- proj_count + 1L
      message(sprintf(">>   %d features mapped from this project", length(af_list)))
    }

    if (!length(feature_map)) {
      stop("Could not map any feature names to alignedFeatureIds. ",
           "The .sirius projects may not contain matching data.", call. = FALSE)
    }
    message(">> Mapped ", length(feature_map), " feature(s) across ", proj_count, " project(s).")

    # Start final service with first project for fingerprint extraction
    .kill_sirius_service(port)
    api_url <- .start_sirius_service(sirius_exe, dirname(sirius_dirs[1]),
                                      port, timeout_sec,
                                      project_files = sirius_dirs[1])
    proj_list <- list(list(projectId = gsub("\\.sirius$", "", basename(sirius_dirs[1]))))
  } else {
    api_url <- sub("/+$", "", api_base)
    message(">> Using existing SIRIUS API at ", api_url,
            " - projects must already be visible to that service.")

    # --- Bypass proxy for localhost only ---
    old_no_proxy <- Sys.getenv("no_proxy", unset = "")
    new_val <- if (nzchar(old_no_proxy)) {
      paste0(old_no_proxy, ",localhost,127.0.0.1")
    } else {
      "localhost,127.0.0.1"
    }
    Sys.setenv(no_proxy = new_val)
    on.exit({
      if (nzchar(old_no_proxy)) Sys.setenv(no_proxy = old_no_proxy)
      else Sys.unsetenv("no_proxy")
    }, add = TRUE)

    # --- List projects ---
    proj_list <- .api_get_json(api_url, "/api/projects")
    if (!length(proj_list)) {
      stop("No projects visible via the API. Ensure .sirius files are accessible.", call. = FALSE)
    }
    message(">> API reports ", length(proj_list), " project(s).")

    # --- Build feature name mapping ---
    feature_map <- .build_feature_map(proj_list, api_url)
    if (!length(feature_map)) {
      stop("Could not map any feature names to alignedFeatureIds. ",
           "The .sirius projects may not contain matching data.", call. = FALSE)
    }
    message(">> Mapped ", length(feature_map), " feature(s) to alignedFeatureId.")
  }

  message(">> SIRIUS REST API at ", api_url)

  # --- Bypass proxy for localhost only ---
  old_no_proxy <- Sys.getenv("no_proxy", unset = "")
  new_val <- if (nzchar(old_no_proxy)) {
    paste0(old_no_proxy, ",localhost,127.0.0.1")
  } else {
    "localhost,127.0.0.1"
  }
  Sys.setenv(no_proxy = new_val)
  on.exit({
    if (nzchar(old_no_proxy)) Sys.setenv(no_proxy = old_no_proxy)
    else Sys.unsetenv("no_proxy")
  }, add = TRUE)

  # --- Extract fingerprints (per-project loop) ---
  n_total <- nrow(fp_table)
  n_ok <- sum(has_fp)
  failed_names <- character()

  pb <- NULL
  pb <- utils::txtProgressBar(min = 0, max = n_total, style = 3)

  # Group features by project for efficient extraction
  if (start_service) {
    # Iterate over each project file, restarting service with --project <file>
    for (j in seq_along(sirius_dirs)) {
      proj_file <- sirius_dirs[j]
      proj_name <- gsub("\\.sirius$", "", basename(proj_file))

      # Find rows in fp_table whose feature names belong to this project
      proj_features <- names(feature_map)[
        vapply(feature_map, function(v) identical(v$projectId, proj_name), logical(1L))
      ]
      if (!length(proj_features)) next

      # Restart service with this project
      .kill_sirius_service(port)
      api_url <- .start_sirius_service(sirius_exe, dirname(proj_file),
                                        port, timeout_sec,
                                        project_files = proj_file)
      if (is.null(api_url)) next

      message(sprintf(">> [%d/%d] Extracting from project: %s (%d features)",
                      j, length(sirius_dirs), proj_name, length(proj_features)))

      for (fn in proj_features) {
        utils::setTxtProgressBar(pb, n_ok)
        # Skip only if every candidate row for this feature already has a vector
        # (the target row is not known until the API answer arrives).
        rows_fn <- which(fp_table$name == fn)
        if (!length(rows_fn)) next
        if (!any(lengths(fp_table$fingerprint[rows_fn]) == 0L)) next

        entry <- feature_map[[fn]]
        if (is.null(entry)) {
          failed_names <- c(failed_names, fn)
          next
        }

        pid  <- entry$projectId
        afid <- entry$alignedFeatureId

        # Fetch formulas with predictedFingerprint optField
        form_path <- paste0(
          "/api/projects/", .url_enc(pid),
          "/aligned-features/", .url_enc(afid),
          "/formulas?optFields=predictedFingerprint"
        )
        formulas <- .api_get_json(api_url, form_path)
        if (is.null(formulas) || !length(formulas)) {
          failed_names <- c(failed_names, fn)
          next
        }

        # Pick the first candidate within top_n_rank
        top <- NULL
        for (fc in formulas) {
          rk <- fc$rank %or% NA_integer_
          if (is.na(rk)) next
          if (rk <= top_n_rank) { top <- fc; break }
        }
        if (is.null(top)) {
          failed_names <- c(failed_names, fn)
          next
        }

        fp_vec <- top$predictedFingerprint
        if (is.null(fp_vec) || !length(fp_vec)) {
          fid <- top$formulaId %or% NULL
          if (!is.null(fid)) {
            fp_path <- paste0(
              "/api/projects/", .url_enc(pid),
              "/aligned-features/", .url_enc(afid),
              "/formulas/", .url_enc(fid),
              "/fingerprint"
            )
            fp_vec <- .api_get_json(api_url, fp_path)
          }
        }

        # Map the returned candidate to its exact row (rank/formula/adduct)
        cand <- .nts_api_pick_fp_row(fp_table, fn, top)
        if (!length(cand)) {
          failed_names <- c(failed_names, fn)
          next
        }
        i <- cand[1L]

        if (length(fp_vec) > 0L) {
          fp_vals <- as.numeric(fp_vec)
          if (anyNA(fp_vals)) next
          data.table::set(fp_table, i, "fingerprint", list(list(fp_vals)))
          n_ok <- n_ok + 1L
        } else {
          failed_names <- c(failed_names, fn)
        }
      }
    }
  } else {
    # Using existing service: single pass over features. Kept as its own loop
    # (NOT merged with the per-project loop above); only the row-fill logic is
    # shared via .nts_api_pick_fp_row(). One API call per distinct feature: the
    # returned candidate is written to the row it actually belongs to, so a
    # feature's other candidate rows are left empty rather than duplicated.
    first_of_feature <- !duplicated(fp_table$name)
    for (i in seq_len(n_total)) {
      utils::setTxtProgressBar(pb, i)
      if (!first_of_feature[i]) next
      fname <- fp_table$name[i]

      rows_fn <- which(fp_table$name == fname)
      if (!any(lengths(fp_table$fingerprint[rows_fn]) == 0L)) next

      entry <- feature_map[[fname]]
      if (is.null(entry)) {
        failed_names <- c(failed_names, fname)
        next
      }

      pid  <- entry$projectId
      afid <- entry$alignedFeatureId

      form_path <- paste0(
        "/api/projects/", .url_enc(pid),
        "/aligned-features/", .url_enc(afid),
        "/formulas?optFields=predictedFingerprint"
      )
      formulas <- .api_get_json(api_url, form_path)
      if (is.null(formulas) || !length(formulas)) {
        failed_names <- c(failed_names, fname)
        next
      }

      top <- NULL
      for (fc in formulas) {
        rk <- fc$rank %or% NA_integer_
        if (is.na(rk)) next
        if (rk <= top_n_rank) { top <- fc; break }
      }
      if (is.null(top)) {
        failed_names <- c(failed_names, fname)
        next
      }

      fp_vec <- top$predictedFingerprint
      if (is.null(fp_vec) || !length(fp_vec)) {
        fid <- top$formulaId %or% NULL
        if (!is.null(fid)) {
          fp_path <- paste0(
            "/api/projects/", .url_enc(pid),
            "/aligned-features/", .url_enc(afid),
            "/formulas/", .url_enc(fid),
            "/fingerprint"
          )
          fp_vec <- .api_get_json(api_url, fp_path)
        }
      }

      # Map the returned candidate to its exact row (rank/formula/adduct)
      cand <- .nts_api_pick_fp_row(fp_table, fname, top)
      if (!length(cand)) {
        failed_names <- c(failed_names, fname)
        next
      }

      if (length(fp_vec) > 0L) {
        fp_vals <- as.numeric(fp_vec)
        if (anyNA(fp_vals)) {
          warning("Fingerprint vector for feature '", fname,
                  "' contains NA after conversion; skipping.", call. = FALSE)
          failed_names <- c(failed_names, fname)
          next
        }
        data.table::set(fp_table, cand[1L], "fingerprint", list(list(fp_vals)))
        n_ok <- n_ok + 1L
      } else {
        failed_names <- c(failed_names, fname)
      }
    }
  }

  if (!is.null(pb)) close(pb)

  # --- Update project ---
  # Ensure fp_table stays as data.table (list-column mods may coerce)
  fp_table <- data.table::as.data.table(fp_table)
  methods::slot(project, "sirius_fingerprints") <- fp_table
  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "extract_sirius_fingerprints_api",
    list(api_url = api_url,
         n_projects = length(proj_list),
         n_features_mapped = length(feature_map),
         n_fingerprints_extracted = n_ok)
  )
  if (length(failed_names)) {
    message(">> ", length(failed_names), " feature(s) could not be resolved: ",
            paste(utils::head(failed_names, 5), collapse = ", "),
            if (length(failed_names) > 5) " ..." else "")
  }
  message("\n>> Done. ", n_ok, " / ", n_total, " rows have fingerprint vectors.")
  invisible(project)
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' @keywords internal
#' @noRd
.nts_api_pick_fp_row <- function(fp_table, fn, top) {
  # Map the API-returned formula candidate back to the exact row of `fp_table`
  # it belongs to. Matching on `name` alone would write a rank-2 fingerprint
  # into the rank-1 row whenever the caller's `top_n_rank` selection lands on a
  # rank > 1 (and would leave the other candidate rows permanently empty).
  #
  # Candidate identity fields from the SIRIUS 6 REST /formulas endpoint:
  #   top$rank, top$molecularFormula, top$adduct, top$formulaId
  # (`%or%` is local to extract_sirius_fingerprints_api(); use a self-contained
  #  NULL/length guard here so this helper is usable on its own.)
  pick1 <- function(x, default) if (is.null(x) || !length(x)) default else x[[1L]]
  target_rank  <- as.integer(pick1(top$rank, NA_integer_))
  target_frmla <- as.character(pick1(top$molecularFormula, NA_character_))
  target_add   <- as.character(pick1(top$adduct, NA_character_))

  cand <- which(fp_table$name == fn)
  if (!length(cand)) return(integer())

  if (!is.na(target_rank) && target_rank %in% fp_table$rank[cand]) {
    cand <- cand[fp_table$rank[cand] == target_rank]
  }
  if (!is.na(target_frmla)) {
    hit <- cand[!is.na(fp_table$molecular_formula[cand]) &
                fp_table$molecular_formula[cand] == target_frmla]
    if (length(hit)) cand <- hit
  }
  if (!is.na(target_add)) {
    hit <- cand[!is.na(fp_table$adduct[cand]) & fp_table$adduct[cand] == target_add]
    if (length(hit)) cand <- hit
  }

  cand <- cand[lengths(fp_table$fingerprint[cand]) == 0L]   # only fill empty rows
  if (length(cand) > 1L) {
    warning("Fingerprint for feature '", fn, "' (rank ", target_rank,
            ") matches ", length(cand), " rows; filling the first only.", call. = FALSE)
  }
  cand
}

#' @keywords internal
.find_sirius_project_dirs <- function(sirius_projects_dir = NULL, output_dir = NULL) {
  candidates <- character()

  if (!is.null(sirius_projects_dir) && dir.exists(sirius_projects_dir)) {
    subdirs <- list.dirs(sirius_projects_dir, recursive = TRUE, full.names = TRUE)
    candidates <- grep("\\.sirius$", subdirs, value = TRUE)
    # also check files (SIRIUS may create file-placeholders instead of dirs)
    f <- list.files(sirius_projects_dir, pattern = "\\.sirius$",
                    recursive = TRUE, full.names = TRUE)
    candidates <- unique(c(candidates, f))
  }

  if (!length(candidates) && !is.null(output_dir) && dir.exists(output_dir)) {
    search_dirs <- unique(c(output_dir, dirname(output_dir),
                            list.dirs(output_dir, recursive = TRUE)))
    for (d in search_dirs) {
      sub <- list.dirs(d, recursive = FALSE, full.names = TRUE)
      candidates <- c(candidates, grep("\\.sirius$", sub, value = TRUE))
      # also check files
      f <- list.files(d, pattern = "\\.sirius$", full.names = TRUE)
      candidates <- c(candidates, f)
    }
    candidates <- unique(candidates)
  }

  if (!length(candidates)) {
    tmp_root <- tempdir()
    candidates <- grep("sirius.*\\.sirius$",
                       list.dirs(tmp_root, recursive = FALSE), value = TRUE)
  }

  # Filter: keep only items that are either directories or non-empty files
  candidates <- candidates[
    vapply(candidates, function(x) dir.exists(x) || (file.exists(x) && file.size(x) > 0), logical(1L))
  ]

  unique(normalizePath(candidates, mustWork = FALSE))
}


#' @keywords internal
.start_sirius_service <- function(sirius_exe, workspace, port, timeout_sec,
                                   project_files = NULL) {
  # Resolve port: use a random free port if not specified (avoids TIME_WAIT)
  if (is.null(port) || is.na(port) || port <= 0L) {
    port <- .find_free_port()
  }
  .kill_all_sirius_services(port)

  base_url <- paste0("http://localhost:", port)

  # Build argument list with --project flags for each .sirius file
  args <- c("--workspace", workspace, "service", "-p", as.character(port), "--headless")

  # Prepend --project flags so SIRIUS opens the projects at startup
  # (paths must NOT be shQuote'd - the sprintf('\"%s\"') below handles quoting)
  if (length(project_files)) {
    proj_args <- c()
    for (pf in project_files) {
      proj_args <- c(proj_args, "--project", pf)
    }
    args <- c(proj_args, args)
  }

  # --- Launch via PowerShell Start-Process (CMD START is unreliable) ---
  if (.Platform$OS.type == "windows") {
    ps_arg_list <- paste(
      sprintf('"%s"', args),
      collapse = ","
    )
    ps_cmd <- sprintf(
      'Start-Process -FilePath "%s" -ArgumentList %s -NoNewWindow -PassThru',
      sirius_exe, ps_arg_list
    )
    system2("powershell", c("-NoProfile", "-NonInteractive", "-Command", shQuote(ps_cmd)),
            wait = FALSE, stdout = FALSE, stderr = FALSE)
  } else {
    system2(sirius_exe, args,
            wait = FALSE, stdout = FALSE, stderr = FALSE)
  }

  # --- Wait loop with periodic status output ---
  deadline <- Sys.time() + timeout_sec
  last_alive_report <- 0L

  while (Sys.time() < deadline) {
    pid_now <- .find_pid_on_port(port)

    elapsed <- as.integer(difftime(Sys.time(), deadline - timeout_sec, units = "secs"))
    if (elapsed - last_alive_report >= 30 && elapsed >= 30) {
      last_alive_report <- elapsed
      if (!is.null(pid_now) && !is.na(pid_now) && pid_now > 0) {
        message(sprintf("  [%ds] waiting for SIRIUS (PID %d on port %d)...",
                        elapsed, pid_now, port))
      } else {
        message(sprintf("  [%ds] no process on port %d yet...", elapsed, port))
      }
    }

    # Try API endpoint
    resp <- tryCatch(
      httr::GET(paste0(base_url, "/api/projects"),
                httr::accept("application/json"),
                httr::timeout(5)),
      error = function(e) NULL
    )
    if (!is.null(resp) && resp$status_code == 200L) {
      message(">> SIRIUS service ready (port ", port, ")")
      return(base_url)
    }

    Sys.sleep(3)
  }

  # --- Headless failed - clean up orphan ---
  stuck_pid <- .find_pid_on_port(port)
  if (!is.null(stuck_pid) && !is.na(stuck_pid) && stuck_pid > 0) {
    system2("taskkill", c("/f", "/pid", as.character(stuck_pid)),
            stdout = FALSE, stderr = FALSE)
    Sys.sleep(1)
  }

  # --- Fallback 1: auto-probe GUI API on port 52423 ---
  gui_url <- "http://localhost:52423"
  gui_resp <- tryCatch(
    httr::GET(paste0(gui_url, "/api/projects"),
              httr::accept("application/json"), httr::timeout(5)),
    error = function(e) NULL
  )
  if (!is.null(gui_resp) && gui_resp$status_code == 200L) {
    # Auto-copy projects to SIRIUS GUI workspace
    gui_ws <- file.path(Sys.getenv("USERPROFILE", "~"), ".sirius-6.3")
    n_copied <- 0L
    if (dir.exists(gui_ws)) {
      for (sd in list.dirs(workspace, recursive = FALSE, full.names = TRUE)) {
        if (grepl("\\.sirius$", sd)) {
          target <- file.path(gui_ws, basename(sd))
          if (!dir.exists(target)) {
            file.copy(sd, gui_ws, recursive = TRUE)
            n_copied <- n_copied + 1L
          }
        }
      }
    }
    message(sprintf(
      "\nAuto-detected SIRIUS GUI API at %s (%d project(s) copied to GUI ws).\nRe-trying with start_service = FALSE...\n",
      gui_url, n_copied))
    return(gui_url)
  }

  # --- Fallback 2: copy to GUI workspace and tell user what to do ---
  gui_ws <- file.path(Sys.getenv("USERPROFILE", "~"), ".sirius-6.3")
  n_copied <- 0L
  if (dir.exists(gui_ws)) {
    for (sd in list.dirs(workspace, recursive = FALSE, full.names = TRUE)) {
      if (grepl("\\.sirius$", sd)) {
        target <- file.path(gui_ws, basename(sd))
        if (!dir.exists(target)) {
          file.copy(sd, gui_ws, recursive = TRUE)
          n_copied <- n_copied + 1L
        }
      }
    }
  }

  # List workspace (diagnostic)
  stop(.headless_instructions_msg(timeout_sec, workspace, gui_ws, n_copied), call. = FALSE)
}


#' @keywords internal
.headless_instructions_msg <- function(timeout_sec, workspace, gui_ws, n_copied) {
  base_url <- "http://localhost:52423"
  msg <- sprintf(
    "SIRIUS headless service could not become ready within %d sec.\n\n",
    timeout_sec)
  msg <- paste0(msg,
    "This is a known issue on some Windows + SIRIUS 6.x combinations -\n",
    "the JVM binds the port but the embedded Web server never completes\n",
    "initialization. This is not an R code problem, but we offer two paths:\n\n")

  if (n_copied > 0L) {
    msg <- paste0(msg,
      sprintf(".sirius projects (%d) have been auto-copied to:\n  %s\n",
              n_copied, gui_ws),
      "The GUI will pick them up automatically.\n\n")
  }

  msg <- paste0(msg,
    "=== Option A (recommended, no headless) ===\n",
    "  1. Close any existing SIRIUS GUI.\n",
    "  2. Launch SIRIUS GUI (via Start Menu or sirius.exe).\n",
    "  3. Wait for it to fully load.\n",
    "  4. In R:\n",
    "       extract_sirius_fingerprints_api(project,\n",
    sprintf("         api_base = \"%s\", start_service = FALSE)\n\n", base_url),
    "=== Option B (retry with longer wait) ===\n",
    "  Specify a longer timeout_sec:\n",
    "       extract_sirius_fingerprints_api(project, ..., timeout_sec = 600)\n")

  msg
}


#' @keywords internal
.find_free_port <- function(lower = 10000L, upper = 60000L) {
  for (attempt in 1:10) {
    candidate <- sample(lower:upper, 1L)
    raw <- system2("netstat", "-ano", stdout = TRUE, stderr = FALSE)
    pattern <- sprintf(":%d\\b", candidate)
    if (!any(grepl(pattern, raw))) {
      return(candidate)
    }
  }
  # Fallback
  8899L
}


#' @keywords internal
.kill_all_sirius_services <- function(port = NULL) {
  if (.Platform$OS.type == "windows") {
    # Kill specific java process listening on the given port
    if (!is.null(port)) {
      pid <- .find_pid_on_port(port)
      if (length(pid) && !is.na(pid) && pid > 0) {
        system2("taskkill", c("/f", "/pid", as.character(pid)),
                stdout = FALSE, stderr = FALSE)
        Sys.sleep(1)
        # Verify port is released
        deadline <- Sys.time() + 5
        while (Sys.time() < deadline) {
          pid2 <- .find_pid_on_port(port)
          if (is.null(pid2) || is.na(pid2) || pid2 <= 0) break
          Sys.sleep(0.3)
        }
        return(invisible())
      }
    }
    # If no specific port given, find and kill any headless SIRIUS
    raw <- system2("netstat", "-ano", stdout = TRUE, stderr = FALSE)
    lines <- grep("LISTENING", raw, value = TRUE)
    if (!length(lines)) return(invisible())
    parts <- strsplit(trimws(lines), "\\s+")
    for (p in parts) {
      pid <- suppressWarnings(as.integer(p[length(p)]))
      if (!is.na(pid) && pid > 0) {
        # Only kill if the process image is java (targeted)
        task_info <- tryCatch(
          system2("tasklist", c("/fi", sprintf("PID eq %d", pid), "/fo", "csv", "/nh"),
                  stdout = TRUE, stderr = FALSE),
          error = function(e) NULL
        )
        if (length(task_info) && any(grepl("java", task_info, ignore.case = TRUE))) {
          system2("taskkill", c("/f", "/pid", as.character(pid)),
                  stdout = FALSE, stderr = FALSE)
        }
      }
    }
    Sys.sleep(1)
  } else {
    system2("pkill", "-f \"sirius.*service.*headless\"",
            stdout = FALSE, stderr = FALSE)
  }
  invisible()
}


#' @keywords internal
.find_pid_on_port <- function(port) {
  tryCatch({
    raw <- system2("netstat", "-ano", stdout = TRUE, stderr = FALSE)
    lines <- grep(sprintf(":%d\\b", port), raw, value = TRUE)
    if (length(lines)) {
      parts <- strsplit(trimws(lines[1]), "\\s+")[[1]]
      suppressWarnings(as.integer(parts[length(parts)]))
    } else NULL
  }, error = function(e) NULL)
}


#' @keywords internal
.kill_sirius_service <- function(port) {
  if (.Platform$OS.type == "windows") {
    pid <- .find_pid_on_port(port)
    if (length(pid) && !is.na(pid) && pid > 0) {
      system2("taskkill", c("/f", "/pid", as.character(pid)),
              stdout = FALSE, stderr = FALSE)
      Sys.sleep(1)
    } else {
      .kill_all_sirius_services()
    }
  } else {
    system2("pkill", "-f \"sirius.*service.*headless\"",
            stdout = FALSE, stderr = FALSE)
  }
  invisible()
}


#' @keywords internal
.api_get_json <- function(base_url, path) {
  url <- paste0(base_url, path)
  resp <- tryCatch(
    httr::GET(url, httr::accept("application/json"), httr::timeout(30)),
    error = function(e) {
      warning("HTTP GET ", url, " failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(resp)) return(NULL)
  if (resp$status_code != 200L) {
    warning("GET ", url, " returned status ", resp$status_code, call. = FALSE)
    return(NULL)
  }
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  if (is.null(txt) || !nzchar(txt)) return(NULL)
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}


#' @keywords internal
.build_feature_map <- function(proj_list, api_url) {
  result <- list()
  for (pj in proj_list) {
    pid <- pj$projectId
    af_path <- paste0("/api/projects/", .url_enc(pid), "/aligned-features")
    af_list <- .api_get_json(api_url, af_path)
    if (is.null(af_list)) next
    for (af in af_list) {
      fname <- trimws(af$name)
      if (is.null(fname) || !nzchar(fname)) next
      if (is.null(result[[fname]])) {
        result[[fname]] <- list(projectId = pid,
                                alignedFeatureId = af$alignedFeatureId)
      }
    }
  }
  result
}


#' @keywords internal
.url_enc <- function(x) {
  gsub(" ", "%20", x, fixed = TRUE)
}
