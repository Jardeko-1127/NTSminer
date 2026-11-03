#' Export NTSProject MS2 to SIRIUS .ms format
#'
#' Converts purified MS2 spectra from [purify_ms2_spectra] into a single .ms file
#' compatible with SIRIUS 5.x.
#'
#' @param project An [NTSProject] with populated `@clean_ms2`.
#' @param file Path to the output `.ms` file.
#' @param ionization Ionization mode (e.g., `[M-H]-`, `[M+H]+`, `[M+Cl]-`).
#' @param ms1_mz_tol m/z tolerance for MS1 precursor matching (default 10 ppm).
#'
#' @return Invisibly, the path to the exported file.
#' @export
export_to_sirius_ms <- function(project, file, ionization = "[M-H]-", ms1_mz_tol = 10) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  ms2_list <- methods::slot(project, "clean_ms2")
  if (!length(ms2_list)) {
    stop("No purified MS2 data found in project. Run purify_ms2_spectra() first.", call. = FALSE)
  }

  con <- file(file, open = "wt")
  on.exit(close(con))

  for (fid in names(ms2_list)) {
    spec <- ms2_list[[fid]]
    cat(">compound", fid, "\n", file = con)
    cat(">feature_id", fid, "\n", file = con)   # SIRIUS 6 uses >feature_id for mapping
    cat(">parentmass", spec$precursorMZ, "\n", file = con) # SIRIUS 6 requires parentmass
    cat(">ionization", ionization, "\n", file = con) 
    cat(">ms1\n", file = con)
    cat(spec$precursorMZ, "100\n", file = con)
    cat(">ms2\n", file = con)
    peaks <- spec$spectrum
    if (nrow(peaks) > 0) {
      for (i in seq_len(nrow(peaks))) {
        cat(peaks$mz[i], peaks$intensity[i], "\n", file = con)
      }
    }
    cat("\n", file = con)
  }

  message(">> Exported ", length(ms2_list), " features to ", file)
  return(invisible(file))
}

# =========================================================================
# Internal: write raw ms2data to .ms format (multiple >ms2 blocks)
# =========================================================================

#' @keywords internal
.write_ms2data_ms <- function(ms2_raw, file, ionization = "[M-H]-") {
  # ms2_raw is a named list of features, each being a list of scan records.
  # Each scan record: list(feature_id, id, scan_idx, precursorMz, rt,
  #                        mz = numeric, intensity = numeric, ...)
  con <- file(file, open = "wt")
  on.exit(close(con))

  for (fid in names(ms2_raw)) {
    scans <- ms2_raw[[fid]]
    if (!length(scans)) next

    # Use the precursor mz from the first scan
    pm <- scans[[1L]]$precursorMz
    if (!is.finite(pm)) next

    cat(">compound", fid, "\n", file = con)
    cat(">feature_id", fid, "\n", file = con)
    cat(">parentmass", pm, "\n", file = con)
    cat(">ionization", ionization, "\n", file = con)
    cat(">ms1\n", file = con)
    cat(pm, "100\n", file = con)

    for (sc in scans) {
      mz <- sc$mz
      int <- sc$intensity
      if (!length(mz) || !length(int)) next
      # keep only finite peaks
      ok <- is.finite(mz) & is.finite(int) & int > 0
      if (!any(ok)) next
      mz <- mz[ok]
      int <- int[ok]
      cat(">ms2\n", file = con)
      for (j in seq_along(mz)) {
        cat(mz[j], int[j], "\n", file = con)
      }
    }
    cat("\n", file = con)
  }

  n_feat <- length(ms2_raw)
  n_scan <- sum(vapply(ms2_raw, length, integer(1L)))
  message(">> Exported ", n_feat, " features (", n_scan, " scans) to ", file)
  invisible(file)
}

# =========================================================================
# Internal: SIRIUS output-dir lifecycle + batch-size / column helpers
# =========================================================================

#' @keywords internal
#' @noRd
.nts_sirius_output_marker <- function(output_dir) {
  # Sentinel written when NTSminer creates an output directory. Its presence is
  # the only thing that authorises recursive deletion of that directory, so a
  # user directory that merely happens to share the name is never wiped.
  file.path(output_dir, ".nts_sirius_output")
}

#' @keywords internal
#' @noRd
.nts_sirius_dir_is_ours <- function(output_dir) {
  file.exists(.nts_sirius_output_marker(output_dir))
}

#' @keywords internal
#' @noRd
.nts_sirius_mark_output_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  writeLines("NTSminer SIRIUS output directory",
             file.path(output_dir, ".nts_sirius_output"))
  invisible(output_dir)
}

#' @keywords internal
#' @noRd
.nts_sirius_prepare_output_dir <- function(output_dir, overwrite = FALSE) {
  # Guarded lifecycle for a SIRIUS output directory. The directory is deleted
  # ONLY when it carries the marker written by NTSminer, or when the caller
  # explicitly opted in with `overwrite = TRUE`; otherwise a pre-existing
  # directory is left untouched and an error is raised. (An earlier version
  # unconditionally `unlink(recursive = TRUE)`-ed whatever path it was given.)
  if (dir.exists(output_dir) || file.exists(output_dir)) {
    if (!.nts_sirius_dir_is_ours(output_dir) && !isTRUE(overwrite)) {
      stop("`output_dir` already exists and was not created by NTSminer:\n  ",
           output_dir,
           "\n  Refusing to delete it. Pass `overwrite = TRUE` if that is really intended.",
           call. = FALSE)
    }
    if (!.nts_sirius_dir_is_ours(output_dir)) {
      warning("Overwriting a directory that was not created by NTSminer: ",
              output_dir, call. = FALSE)
    }
    message(">> Cleaning existing NTSminer SIRIUS project space: ", output_dir)
    unlink(output_dir, recursive = TRUE)
  }
  .nts_sirius_mark_output_dir(output_dir)
  invisible(output_dir)
}

#' @keywords internal
#' @noRd
.nts_sirius_normalize_elements <- function(x) {
  # Normalise the `elements` argument to the string passed to SIRIUS
  # `--elements-considered`. The value is an **element-symbol string** and is
  # case-sensitive ("Cl", "Br", "Se" are two-letter symbols), so the case is
  # preserved verbatim - uppercasing it would turn "SClBrSeI" into "SCLBRSEI"
  # and break the two-letter symbols. C/H/N/O/P are always considered by SIRIUS
  # and are stripped ONLY when the caller supplies them as separate tokens
  # (`"C,H,N,O,P,S"`); stripping them out of a concatenated string such as
  # `"SClBrSeI"` cannot be done reliably. Returns NULL when nothing is left
  # (the caller then omits the flag and SIRIUS uses its own default).
  #
  # NOTE: the exact accepted format of `--elements-considered` (concatenated
  # symbols vs. comma-separated list) is still to be confirmed against a real
  # CLI run; the default `"SClBrSeI"` reproduces the known-good, correctly cased
  # format used by the previous default `"CHNOPFClBr"`.
  if (is.null(x) || !length(x)) return(NULL)
  x <- trimws(as.character(x))
  x <- x[nzchar(x)]
  if (!length(x)) return(NULL)
  if (length(x) == 1L && grepl("[,;[:space:]]", x)) {
    x <- unlist(strsplit(x, "[,;[:space:]]+"), use.names = FALSE)
    x <- x[nzchar(x)]
  }
  if (length(x) > 1L) x <- setdiff(x, c("C", "H", "N", "O", "P"))
  if (!length(x)) return(NULL)
  paste(x, collapse = "")
}

#' @keywords internal
#' @noRd
.nts_sirius_elements_flag <- function(elements) {
  # Build the `--elements-considered` fragment (possibly empty) for a SIRIUS
  # command line. Warns once when the set omits S, because sulfur-containing
  # formulas (PFAS sulfonates / sulfonamides) will then not be found.
  els <- .nts_sirius_normalize_elements(elements)
  if (is.null(els)) return(character(0))
  if (!grepl("S", els, fixed = TRUE)) {
    warning("`elements` does not include S; sulfur-containing formulas ",
            "(PFAS sulfonates/sulfonamides) may not be detected.", call. = FALSE)
  }
  c("--elements-considered", els)
}

#' @keywords internal
#' @noRd
.nts_check_batch_size <- function(batch_size, n_total) {
  # `ceiling(n_total / batch_size)` used to propagate whatever the caller passed:
  # batch_size = 0 -> Inf -> seq_len(Inf) error; batch_size = NA -> NA batches;
  # a vector -> recycled silently. Validate and coerce to a single integer.
  if (is.null(batch_size)) return(as.integer(n_total))
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      !is.finite(batch_size) || batch_size < 1) {
    stop("`batch_size` must be NULL or a single finite number >= 1; got: ",
         paste(utils::capture.output(str(batch_size)), collapse = " "),
         call. = FALSE)
  }
  as.integer(floor(batch_size))
}

#' @keywords internal
#' @noRd
.nts_canopus_npc_columns <- function(cn) {
  # Resolve CANOPUS NPC classification columns across SIRIUS versions.
  # ALWAYS returns the same 6 fields (unused ones are NA_character_), so callers
  # never hit `if (is.na(NULL))`.
  #   mode       = "multi" (SIRIUS 6: NPC#pathway/superclass/class)
  #              | "single" (older SIRIUS: one bare `npc` column)
  #              | "none"
  #   across     = column name for the whole string (`canopus_npc`) in "single" mode
  #   pathway / superclass / class / subclass = column names, or NA_character_
  #
  # Matching is EXACT (case-insensitive), never `grep("npc", ...)`: a fuzzy match
  # would hit `npc#pathway probability` first and silently mis-assign levels.
  lc <- tolower(cn)
  pick <- function(p) { i <- which(lc == p); if (length(i)) cn[i[1L]] else NA_character_ }
  none <- list(mode = "none", across = NA_character_, pathway = NA_character_,
               superclass = NA_character_, class = NA_character_, subclass = NA_character_)
  if (any(lc == "npc")) {
    none$mode   <- "single"
    none$across <- cn[which(lc == "npc")[1L]]
    return(none)
  }
  pw <- pick("npc#pathway"); sc <- pick("npc#superclass")
  cl <- pick("npc#class");   sb <- pick("npc#subclass")
  if (!is.na(pw) || !is.na(sc) || !is.na(cl) || !is.na(sb)) {
    none$mode <- "multi"; none$pathway <- pw
    none$superclass <- sc; none$class <- cl; none$subclass <- sb
    return(none)
  }
  none
}

#' @keywords internal
#' @noRd
.nts_sirius_repair_id <- function(x, orig_feat_names) {
  # SIRIUS 6 prefixes mappingFeatureId with batch/tmp noise. Recover the true
  # feature name: exact match, else longest suffix at an underscore boundary
  # (the boundary prevents "FT01" from matching "batch_FT011").
  if (x %in% orig_feat_names) return(x)
  pattern <- paste0("(^|_)(", paste(orig_feat_names, collapse = "|"), ")$")
  m <- regexec(pattern, x, perl = TRUE)
  caps <- regmatches(x, m)[[1L]]
  if (length(caps) >= 3L && nzchar(caps[3L])) {
    return(caps[3L])
  }
  x   # last resort: keep as-is
}

#' @keywords internal
#' @noRd
.nts_sirius_apply_canopus_npc <- function(canopus_dt, npc_cols) {
  # Shared by compute_sirius_fingerprints() (inline) and parse_canopus_results().
  # Fills canopus_npc / canopus_superclass / canopus_class / canopus_subclass.
  # Returns the (mutated) data.table.
  nr <- nrow(canopus_dt)
  none_warning <- function() {
    warning(
      "CANOPUS table found but no NPC classification columns (looked for ",
      "`npc`, `npc#pathway`, `npc#superclass`, `npc#class`, `npc#subclass`). ",
      "First columns: ", paste(utils::head(names(canopus_dt), 20L), collapse = ", "),
      call. = FALSE
    )
    canopus_dt[, canopus_npc        := rep(NA_character_, nr)]
    canopus_dt[, canopus_superclass := rep(NA_character_, nr)]
    canopus_dt[, canopus_class      := rep(NA_character_, nr)]
    canopus_dt[, canopus_subclass   := rep(NA_character_, nr)]
    canopus_dt
  }

  if (identical(npc_cols$mode, "multi")) {
    gv <- function(col) if (!is.na(col)) as.character(canopus_dt[[col]]) else rep(NA_character_, nr)
    pw <- gv(npc_cols$pathway); sc <- gv(npc_cols$superclass)
    cl <- gv(npc_cols$class);   sb <- gv(npc_cols$subclass)
    lv <- vapply(seq_len(nr), function(k) {
      v <- c(pw[k], sc[k], cl[k])
      v <- v[!is.na(v) & nzchar(trimws(v))]
      if (length(v)) paste(trimws(v), collapse = "; ") else NA_character_
    }, character(1L))
    canopus_dt[, canopus_npc        := lv]
    canopus_dt[, canopus_superclass := sc]
    canopus_dt[, canopus_class      := cl]
    canopus_dt[, canopus_subclass   := sb]   # SIRIUS 6 has no NPC#subclass -> NA
    return(canopus_dt)
  }

  if (identical(npc_cols$mode, "single")) {
    raw_val   <- as.character(canopus_dt[[npc_cols$across]])
    npc_split <- strsplit(raw_val, "[\\:/;]\\s*")
    canopus_dt[, canopus_npc        := raw_val]
    canopus_dt[, canopus_superclass := vapply(npc_split, function(x) if (length(x) >= 1L) trimws(x[1L]) else NA_character_, character(1L))]
    canopus_dt[, canopus_class      := vapply(npc_split, function(x) if (length(x) >= 2L) trimws(x[2L]) else NA_character_, character(1L))]
    canopus_dt[, canopus_subclass   := vapply(npc_split, function(x) if (length(x) >= 3L) trimws(x[3L]) else NA_character_, character(1L))]
    return(canopus_dt)
  }

  none_warning()
}

#' @keywords internal
#' @noRd
.nts_sirius_merge_canopus <- function(out, canopus_dt, orig_feat_names) {
  # Merge CANOPUS classifications into `out` keyed by feature identity, not just
  # by (formula, adduct, rank): two different features that happen to share a
  # formula must not borrow each other's classification.
  npc_cols <- .nts_canopus_npc_columns(names(canopus_dt))
  canopus_dt <- .nts_sirius_apply_canopus_npc(canopus_dt, npc_cols)

  id_candidates <- c("mappingfeatureid", "mappingFeatureId")
  id_col <- id_candidates[id_candidates %in% names(canopus_dt)][1L]
  if (!is.na(id_col)) {
    canopus_dt[, name := vapply(
      trimws(as.character(canopus_dt[[id_col]])),
      .nts_sirius_repair_id, character(1L), USE.NAMES = FALSE,
      orig_feat_names = orig_feat_names)]
    join_by <- c("name", "molecular_formula", "adduct", "rank")
  } else {
    warning(
      "CANOPUS summary has no `mappingFeatureId`; falling back to a ",
      "(molecular_formula, adduct, rank) join. Classification may cross-talk ",
      "between different features that share a formula.", call. = FALSE
    )
    join_by <- c("molecular_formula", "adduct", "rank")
  }

  # merge() cannot match NA keys, so an all-NA rank column would match nothing.
  # Sentinel -1L keeps those rows joinable (A-P2-2).
  if ("rank" %in% join_by) {
    if (all(is.na(canopus_dt$rank))) canopus_dt[is.na(rank), rank := -1L]
    if ("rank" %in% names(out) && all(is.na(out$rank))) out[is.na(rank), rank := -1L]
  }

  keep <- intersect(names(canopus_dt),
    unique(c(join_by, "canopus_npc", "canopus_superclass",
             "canopus_class", "canopus_subclass")))
  canopus_dt <- canopus_dt[, ..keep]
  canopus_dt <- unique(canopus_dt, by = join_by)

  n_before <- nrow(out)
  out <- merge(out, canopus_dt, by = join_by, all.x = TRUE, sort = FALSE)
  if (nrow(out) != n_before) {
    warning("CANOPUS merge changed row count (", n_before, " -> ", nrow(out),
            "); this indicates a non-unique CANOPUS key.", call. = FALSE)
  }
  out
}

#' @keywords internal
#' @noRd
.nts_system2_capture <- function(exe, args) {
  # Thin, mockable seam around base::system2: returns the exit status together
  # with the captured stdout/stderr. Kept separate so the failure-handling path
  # can be exercised by tests without spawning a real process.
  res <- suppressWarnings(system2(exe, args = args, stdout = TRUE, stderr = TRUE))
  status <- attr(res, "status")
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = as.character(res))
}

#' @keywords internal
#' @noRd
.nts_sirius_run_cli <- function(exe, args, log_file = NULL, cli_log = NULL,
                                quiet = FALSE, echo = TRUE,
                                on_fail = c("stop", "warn", "none"),
                                label = NULL) {
  # Run one SIRIUS command line. The exit status is captured (never discarded:
  # an earlier version passed `stdout = FALSE` and threw the status away, so a
  # failed invocation was indistinguishable from a successful one). `log_file`
  # receives the full stdout/stderr; `cli_log` receives an appended command line.
  # On failure the log path is included in the stop()/warning() message.
  on_fail <- match.arg(on_fail)
  if (!is.null(cli_log)) {
    cat(sprintf("[%s] %s %s\n", format(Sys.time()), exe, paste(args, collapse = " ")),
        file = cli_log, append = TRUE)
  }
  r <- .nts_system2_capture(exe, args)
  if (!is.null(log_file)) writeLines(r$output, log_file)
  if (isTRUE(echo) && !isTRUE(quiet) && length(r$output)) cat(r$output, sep = "\n")
  if (r$status != 0L) {
    msg <- paste0(
      "SIRIUS failed with exit status ", r$status,
      if (!is.null(label)) paste0(" ", label) else "", ".",
      "\n  Command: ", exe, " ", paste(args, collapse = " "),
      if (!is.null(log_file)) paste0("\n  Output : ", log_file) else ""
    )
    if (on_fail == "stop") stop(msg, call. = FALSE)
    if (on_fail == "warn") warning(msg, call. = FALSE)
  }
  invisible(r)
}

#' @keywords internal
#' @noRd
.nts_sirius_read_canopus <- function(files) {
  # Read a set of per-batch `canopus_formula_summary.tsv` files. Returns a list
  # of data.frames (lower-cased column names, as the downstream code expects).
  # A file that exists but has ZERO data rows is reported explicitly: an empty
  # CANOPUS table used to disappear silently, which is indistinguishable from
  # "CANOPUS not run".
  out <- list()
  empty_files <- character(0)
  for (f in files) {
    if (!file.exists(f)) next
    tbl <- tryCatch(
      utils::read.table(f, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE, quote = "\"",
                        comment.char = "", check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(tbl) || !nrow(tbl)) {
      empty_files <- c(empty_files, f)
      next
    }
    names(tbl) <- tolower(names(tbl))
    out[[length(out) + 1L]] <- tbl
  }
  if (length(empty_files)) {
    warning("CANOPUS summary file(s) contained no data rows (header only): ",
            paste(basename(empty_files), collapse = ", "), call. = FALSE)
  }
  out
}

#' @keywords internal
#' @noRd
.nts_sirius_col_map <- function(cn) {
  # Resolve the score columns of a SIRIUS `formula_identifications` table by
  # exact name first, then a conservative fuzzy fallback. `explained_peaks`
  # deliberately maps to `numExplainedPeaks` (an integer COUNT); a fuzzy
  # `grep("explained", ...)` would wrongly pick `explainedIntensity` (a 0-1
  # fraction). When no confidence column exists (always the case for SIRIUS 6
  # formula-level output) a single warning is emitted here.
  col_map <- list(
    sirius_score    = grep("^SiriusScore$", cn, value = TRUE),
    zodiac_score    = grep("^ZodiacScore$", cn, value = TRUE),
    confidence      = grep("^Confidence$|^confidence$", cn, value = TRUE),
    explained_peaks = grep("^numExplainedPeaks$", cn, value = TRUE),
    rank            = grep("^rank$|^Rank$", cn, value = TRUE)
  )
  if (!length(col_map$sirius_score)) {
    col_map$sirius_score <- grep("sirius.*score", cn, ignore.case = TRUE, value = TRUE)
    if (length(col_map$sirius_score)) {
      preferred <- grep("normaliz", col_map$sirius_score, ignore.case = TRUE,
                        invert = TRUE, value = TRUE)
      col_map$sirius_score <- if (length(preferred)) preferred[1L] else col_map$sirius_score[1L]
    }
  }
  if (!length(col_map$zodiac_score)) {
    col_map$zodiac_score <- grep("zodiac", cn, ignore.case = TRUE, value = TRUE)
  }
  if (!length(col_map$confidence)) {
    # Fuzzy fallback, but prefer ConfidenceScoreApproximate over any other
    # "confiden..." column that may appear in future SIRIUS formula tables.
    cand <- grep("confiden", cn, ignore.case = TRUE, value = TRUE)
    if (length(cand)) {
      approx <- grep("^ConfidenceScoreApproximate$", cand, ignore.case = TRUE, value = TRUE)
      col_map$confidence <- if (length(approx)) approx[1L] else cand[1L]
    }
  }
  if (!length(col_map$explained_peaks)) {
    # NOTE: deliberately NOT `grep("explained", ...)` - that would match
    # `explainedIntensity` (a 0-1 fraction), not the peak count. Prefer
    # `numExplainedPeaks`; if it is absent, leave explained_peaks as NA.
    col_map$explained_peaks <- grep("numExplainedPeaks", cn, ignore.case = TRUE, value = TRUE)
  }
  if (!length(col_map$rank)) {
    col_map$rank <- grep("rank", cn, ignore.case = TRUE, value = TRUE)
  }
  if (!length(col_map$confidence)) {
    warning("No `confidence` column in formula_identifications.tsv; `confidence` ",
            "will be NA. SIRIUS 6 exposes ConfidenceScoreExact/Approximate only ",
            "in structure_identifications.tsv.", call. = FALSE)
  }
  col_map
}

#' Run SIRIUS with Progress Feedback
#'
#' Executes SIRIUS for fingerprinting. If `batch_size` is specified, the input 
#' features are processed in chunks to provide a progress bar.
#'
#' @param project An [NTSProject].
#' @param output_dir Directory to store SIRIUS results.
#' @param ionization Ionization mode for the batch (default `[M-H]-`).
#' @param elements Set of elements SIRIUS is allowed to *consider* when
#'   searching formulas, passed through as `--elements-considered`. C, H, N, O
#'   and P are always considered by SIRIUS and are stripped from this value.
#'   Default `"SClBrSeI"` - sulfur must stay in the set for PFAS sulfonates /
#'   sulfonamides to be found. Pass `NULL` (or a set that reduces to nothing) to
#'   omit the flag and use the SIRIUS default. Values may be given as a single
#'   string or a character vector; separators (`,` `;` whitespace) are ignored.
#' @param batch_size Maximum features per SIRIUS invocation (default `20`); a
#'   progress bar is shown across batches.
#' @param subset Vector of feature names to process, or an integer count (e.g., `10` for first 10).
#'   Default `NULL` (all features).
#' @param quiet If `TRUE`, suppress SIRIUS console output. Default `FALSE`.
#' @param sirius_exe Path to the SIRIUS executable. Default `NULL` reads
#'   environment variable `NTSminer_SIRIUS_EXE` (or `SIRIUS_EXE`); you can
#'   also pass an explicit path.
#' @param overwrite If `FALSE` (default), refuse to delete an existing
#'   `output_dir` that NTSminer did not create (tracked by the marker file
#'   `.nts_sirius_output`). Set `TRUE` to force replacement - the existing
#'   directory is then deleted with a warning.
#' @param ... Additional command-line flags for SIRIUS (e.g., `--ppm-max 20`).
#'
#' @export
run_sirius <- function(project,
                       output_dir,
                       ionization = "[M-H]-",
                       # Rare elements SIRIUS may auto-detect. C/H/N/O/P are
                       # always considered; S must stay in this set for PFAS
                       # sulfonates/sulfonamides to be detected.
                       elements = "SClBrSeI",
                       batch_size = 20,
                       subset = NULL,
                       quiet = FALSE, # set FALSE to see SIRIUS live errors
                       sirius_exe = NULL,
                       overwrite = FALSE,
                       ...) {

  .nts_check(project, "project", .arg = "project")
  if (is.null(sirius_exe) || !nzchar(sirius_exe)) {
    sirius_exe <- Sys.getenv("NTSminer_SIRIUS_EXE", Sys.getenv("SIRIUS_EXE", ""))
  }
  if (!nzchar(sirius_exe)) {
    stop(
      "SIRIUS path not set: pass `sirius_exe = ` or set env `NTSminer_SIRIUS_EXE` ",
      "(or `SIRIUS_EXE`).",
      call. = FALSE
    )
  }
  ms2_list <- methods::slot(project, "clean_ms2")
  feat_names <- names(ms2_list)
  
  # Handle subsetting
  if (!is.null(subset)) {
    if (is.numeric(subset)) {
      feat_names <- head(feat_names, subset)
    } else {
      feat_names <- intersect(feat_names, subset)
    }
  }
  
  n_total <- length(feat_names)
  
  if (!file.exists(sirius_exe)) stop("SIRIUS not found at: ", sirius_exe)
  
  # normalise the suffix and convert to an absolute path
  if (!grepl("\\.sirius$", output_dir)) {
    output_dir <- paste0(output_dir, ".sirius")
  }
  
  # resolve the absolute path up front to avoid repeated normalize() inside the loop
  # mustWork = FALSE allows a path that does not exist yet
  abs_output_dir <- normalizePath(output_dir, mustWork = FALSE)

  # Clear the old project space for a clean run - but only if NTSminer created
  # it (marker file present) or the caller explicitly opted in with `overwrite`.
  .nts_sirius_prepare_output_dir(abs_output_dir, overwrite = overwrite)

  # Handle batching
  batch_size <- .nts_check_batch_size(batch_size, n_total)
  n_batches <- ceiling(n_total / batch_size)

  # Resolve `elements` once (empty when the flag should be omitted)
  elements_flag <- .nts_sirius_elements_flag(elements)
  
  # Progress Bar Setup
  pb <- NULL
  is_shiny <- !is.null(shiny::getDefaultReactiveDomain())
  if (!is_shiny) {
    pb <- utils::txtProgressBar(min = 0, max = n_batches, style = 3)
  }

  for (i in seq_len(n_batches)) {
    idx_start <- (i - 1) * batch_size + 1
    idx_end <- min(i * batch_size, n_total)
    batch_feats <- feat_names[idx_start:idx_end]
    
    msg <- sprintf("> [%d/%d] Batch %d: %s", idx_end, n_total, i, 
                   paste(head(batch_feats, 3), collapse = ", "))
    if (length(batch_feats) > 3) msg <- paste0(msg, "...")
    
    if (is_shiny) {
      shiny::incProgress(1/n_batches, detail = msg)
    } else {
      cat("\n", msg, "\n")
      utils::setTxtProgressBar(pb, i)
    }

    # 1. Create temporary .ms file
    tmp_ms <- tempfile(pattern = paste0("sirius_batch_", i, "_"), fileext = ".ms")
    sub_project <- project
    methods::slot(sub_project, "clean_ms2") <- ms2_list[batch_feats]
    export_to_sirius_ms(sub_project, tmp_ms, ionization = ionization)

    # 2. Run SIRIUS (Compute)
    extra_args <- list(...)
    args <- c(
      "-i", shQuote(normalizePath(tmp_ms)),
      "-o", shQuote(abs_output_dir),
      as.character(unlist(extra_args)),
      "formula",
      elements_flag,
      "fingerprint"
    )
    
    # Capture BOTH the console output and the exit status. Previously the status
    # was discarded and the output thrown away when `quiet = TRUE`, so a failed
    # SIRIUS invocation was indistinguishable from a successful one.
    cli_log   <- file.path(abs_output_dir, "sirius_cli.log")
    batch_log <- file.path(abs_output_dir, sprintf("sirius_batch_%d.log", i))
    .nts_sirius_run_cli(
      sirius_exe, args,
      log_file = batch_log, cli_log = cli_log,
      quiet = quiet, echo = !quiet, on_fail = "stop",
      label = paste0("on batch ", i, " (", length(batch_feats), " feature(s))")
    )

    # 3. Export fingerprints (only supported on SIRIUS 5.x; 6.x no-op)
    export_dir <- gsub("\\.sirius$", "_export", abs_output_dir)
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)

    if (file.exists(abs_output_dir)) {
        message(">> Exporting fingerprints to: ", export_dir)
        fp_ok <- .nts_sirius_export_fingerprints(sirius_exe, abs_output_dir, export_dir)
        if (is.null(fp_ok)) {
          message(">> Fingerprint CLI export not available; ",
                  "use parse_sirius_fingerprints() on the _export directory ",
                  "or extract_sirius_fingerprints_api().")
        }
    }
    
    unlink(tmp_ms)
  }

  if (!is_shiny) close(pb)
  message("\n>> SIRIUS processing complete. Results exported to: ", gsub("\\.sirius$", "_export", abs_output_dir))
}

# =========================================================================
# Internal: CLI-based fingerprint export (SIRIUS 6.x / 5.x)
# =========================================================================

#' @keywords internal
.nts_sirius_export_fingerprints <- function(sirius_exe, project_dir, output_dir) {
  if (!dir.exists(project_dir)) return(NULL)

  # Method A: SIRIUS 6.x syntax
  args6 <- c(
    "--project", shQuote(normalizePath(project_dir)),
    "export",
    "--type=fingerprints",
    "--output", shQuote(normalizePath(output_dir))
  )
  status6 <- suppressWarnings(system2(sirius_exe, args = args6,
    stdout = FALSE, stderr = FALSE))
  if (identical(status6, 0L)) {
    fpt_files <- list.files(output_dir, pattern = "\\.fpt$",
      recursive = TRUE, full.names = TRUE)
    if (length(fpt_files)) return(output_dir)
  }

  # Method B: SIRIUS 5.x legacy syntax
  args5 <- c(
    "-p", shQuote(normalizePath(project_dir)),
    "export",
    "--type=fingerprints",
    "--output", shQuote(normalizePath(output_dir))
  )
  status5 <- suppressWarnings(system2(sirius_exe, args = args5,
    stdout = FALSE, stderr = FALSE))
  if (identical(status5, 0L)) {
    fpt_files <- list.files(output_dir, pattern = "\\.fpt$",
      recursive = TRUE, full.names = TRUE)
    if (length(fpt_files)) return(output_dir)
  }

  # Method C: Look for pre-existing .fpt files inside project dir
  fpt_in_proj <- list.files(project_dir, pattern = "\\.fpt$",
    recursive = TRUE, full.names = TRUE)
  if (length(fpt_in_proj)) {
    for (f in fpt_in_proj) {
      file.copy(f, file.path(output_dir, basename(f)), overwrite = TRUE)
    }
    return(output_dir)
  }

  NULL
}

#' @keywords internal
.nts_parse_fpt_files <- function(export_dir, feat_names) {
  fpt_files <- list.files(export_dir, pattern = "\\.fpt$",
    recursive = TRUE, full.names = TRUE)
  if (!length(fpt_files)) return(list())

  fp_map <- list()
  for (f in fpt_files) {
    parts <- strsplit(f, "[/\\\\]")[[1L]]
    comp_id <- NA_character_
    for (p in parts) {
      if (p %in% feat_names) {
        comp_id <- p
        break
      }
    }
    if (is.na(comp_id)) next
    vals <- tryCatch(scan(f, quiet = TRUE), error = function(e) NULL)
    if (length(vals) && is.null(fp_map[[comp_id]])) {
      fp_map[[comp_id]] <- vals
    }
  }
  fp_map
}

#' Parse SIRIUS .fpt fingerprints into @sirius_fingerprints (legacy)
#'
#' Reads CSI:FingerID \code{.fpt} files produced by [run_sirius] and stores
#' the fingerprint vectors into the standard \code{@sirius_fingerprints}
#' slot.
#'
#' @details
#' This function works with the legacy \code{run_sirius()} workflow.
#' For the modern SIRIUS 6.x toolchain, prefer [compute_sirius_fingerprints]
#' which now auto-exports fingerprints via CLI. If CLI export is unavailable,
#' use [extract_sirius_fingerprints_api] as fallback.
#'
#' \code{.fpt} files only contain raw fingerprint vectors \emph{without}
#' formula metadata (adduct, molecular_formula, scores). All metadata
#' columns will be \code{NA}.
#'
#' @param project An [NTSProject].
#' @param sirius_out_dir The output directory from [run_sirius] (where
#'   \code{.fpt} files were exported).
#' @return The updated \code{project} with \code{@sirius_fingerprints}
#'   populated (fingerprint vectors only; metadata columns are \code{NA}).
#' @export
parse_sirius_fingerprints <- function(project, sirius_out_dir) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")

  feat_names <- names(methods::slot(project, "clean_ms2"))
  n_total <- length(feat_names)

  # Reuse the shared .fpt parser
  fp_list <- .nts_parse_fpt_files(sirius_out_dir, feat_names)

  if (!length(fp_list)) {
    warning("No .fpt files found in ", sirius_out_dir,
      ". Possible reasons: SIRIUS failed to find formula candidates for any feature.")
    return(project)
  }

  message(">> Found ", length(fp_list), " fingerprint vectors. Mapping to features...")

  n_success <- length(fp_list)
  n_fail <- n_total - n_success

  # Build data.table matching @sirius_fingerprints schema
  # .fpt files only contain raw fingerprints -- formula metadata is NA
  if (n_success > 0L) {
    out_dt <- data.table::data.table(
      name               = names(fp_list),
      adduct             = NA_character_,
      molecular_formula  = NA_character_,
      rank               = 1L,
      sirius_score       = NA_real_,
      zodiac_score       = NA_real_,
      confidence         = NA_real_,
      explained_peaks    = NA_integer_,
      fingerprint        = unname(fp_list),
      fp_version         = "",
      canopus_npc        = NA_character_,
      canopus_superclass = NA_character_,
      canopus_class      = NA_character_,
      canopus_subclass   = NA_character_
    )
  } else {
    out_dt <- data.table::data.table(
      name = character(0),
      adduct = character(0),
      molecular_formula = character(0),
      rank = integer(0),
      sirius_score = numeric(0),
      zodiac_score = numeric(0),
      confidence = numeric(0),
      explained_peaks = integer(0),
      fingerprint = list(),
      fp_version = character(0),
      canopus_npc = character(0),
      canopus_superclass = character(0),
      canopus_class = character(0),
      canopus_subclass = character(0)
    )
  }

  project <- set_sirius_fingerprints(project, out_dt)

  message("-------------------------------------------")
  message(">> SIRIUS Fingerprint Summary (.fpt -> @sirius_fingerprints):")
  message(">> Total Features:      ", n_total)
  message(">> Successfully Parsed: ", n_success)
  message(">> Failed / No Formula: ", n_fail)
  if (n_fail > 0) message(">> (Note: Some compounds were skipped because SIRIUS couldn't predict a formula)")
  message(">> NOTE: .fpt files only contain fingerprint vectors. Formula metadata (adduct, scores, etc.) is NA.")
  message(">>       For full metadata, use compute_sirius_fingerprints() + extract_sirius_fingerprints_api().")
  message("-------------------------------------------")

  return(project)
}


# =========================================================================
# High-level: compute_sirius_fingerprints
# =========================================================================

#' Compute CSI:FingerID fingerprints via SIRIUS
#'
#' End-to-end wrapper that exports MS2 spectra from the project to SIRIUS
#' \code{.ms} format, runs the SIRIUS toolchain (\code{formulas +
#' fingerprints + canopus (optional) + write-summaries}), parses the resulting formula table,
#' and stores the structured result into \code{@sirius_fingerprints}.
#'
#' Fingerprint vectors are not exported by the CLI toolchain and must be
#' obtained separately via \code{extract_sirius_fingerprints_api()}, which
#' starts a SIRIUS REST API service (requires \code{sirius login} first).
#' See \code{?extract_sirius_fingerprints_api} for details.
#'
#' By default (\code{ms2_source = "clean_ms2"}), the purified/merged spectra
#' from \code{@clean_ms2} are used (one \code{>ms2} block per feature).
#' When \code{ms2_source = "ms2data"}, \emph{all} raw matched scans from
#' \code{@ms2data} are written (multiple \code{>ms2} blocks per feature),
#' letting SIRIUS perform its own scan-level selection and merging.
#'
#' @param project An [NTSProject].
#' @param ms2_source Which MS2 slot to use: \code{"clean_ms2"} (default,
#'   purified merged spectra) or \code{"ms2data"} (all raw matched scans).
#' @param output_dir Directory for SIRIUS project workspace and any
#'   exported files. If \code{NULL} (default), a temporary path is used
#'   and cleaned up unless \code{keep_intermediates = TRUE}.
#' @param ionization Ionization string passed to SIRIUS
#'   (e.g. \code{"[M+H]+"}, \code{"[M-H]-"}). Default \code{"[M-H]-"}.
#' @param elements Set of elements SIRIUS is allowed to *consider* when
#'   searching formulas, passed through as \code{--elements-considered}. C, H,
#'   N, O and P are always considered by SIRIUS and are stripped from this
#'   value. Default \code{"SClBrSeI"} - sulfur must stay in the set for PFAS
#'   sulfonates / sulfonamides to be found. Pass \code{NULL} (or a set that
#'   reduces to nothing) to omit the flag and use the SIRIUS default.
#' @param sirius_exe Path to the SIRIUS executable. If \code{NULL}
#'   (default), reads environment variables \code{NTSminer_SIRIUS_EXE}
#'   then \code{SIRIUS_EXE}.
#' @param batch_size Maximum features per SIRIUS invocation
#'   (default 20).  Used to keep individual runs fast and provide
#'   incremental progress feedback.
#' @param keep_intermediates If \code{TRUE}, the temporary \code{.ms}
#'   files and per-batch SIRIUS project directories are not deleted;
#'   the final combined output directory is never deleted regardless.
#' @param run_canopus If \code{TRUE}, adds CANOPUS NPC classification to the
#'   SIRIUS toolchain (\code{canopus} step before \code{write-summaries}).
#'   NPC pathway classifications are parsed from batch output and
#'   merged into \code{@sirius_fingerprints}. Default \code{FALSE}.
#' @param overwrite If \code{FALSE} (default), refuse to write into an existing
#'   \code{output_dir} that NTSminer did not create (tracked by the marker file
#'   \code{.nts_sirius_output}). Set \code{TRUE} to reuse/replace such a
#'   directory - a warning is emitted.
#' @param ... Additional command-line flags for the SIRIUS \code{formulas}
#'   step (e.g. \code{--ppm-max 10}, \code{--candidates 10}).
#'
#' @return The updated \code{project} with \code{@sirius_fingerprints}
#'   populated.
#'
#' @section Slot schema (\code{@sirius_fingerprints}):
#' \code{data.frame} with one row per feature per formula candidate:
#' \describe{
#'   \item{name}{Feature name - primary key into \code{@clean_ms2}.}
#'   \item{adduct}{Adduct used (e.g. \code{[M-H]-}).}
#'   \item{molecular_formula}{Predicted neutral formula.}
#'   \item{rank}{Rank of this formula for this feature (1 = top).}
#'   \item{sirius_score}{SIRIUS fragmentation-tree score.}
#'   \item{zodiac_score}{ZODIAC score if computed, else \code{NA}.}
#'   \item{confidence}{CSI:FingerID confidence in \verb{[0,1]}. SIRIUS 6 writes
#'     \code{ConfidenceScoreExact}/\code{ConfidenceScoreApproximate} only in
#'     \code{structure_identifications.tsv}, not in
#'     \code{formula_identifications.tsv}, so this is \code{NA} for SIRIUS 6
#'     formula-level output.}
#'   \item{explained_peaks}{Integer count of MS2 peaks explained by this formula
#'     (from the \code{numExplainedPeaks} column - \emph{not}
#'     \code{explainedIntensity}, which is a 0-1 fraction).}
#'   \item{fingerprint}{List-column of numeric CSI:FingerID fingerprint
#'     vectors (length varies by SIRIUS version; typical length 7921 or
#'     similar).}
#'   \item{fp_version}{Fingerprint version string (e.g. \code{"v3"}).}
#'   \item{canopus_npc}{CANOPUS NPC pathway classification string (e.g.
#'     ``Amino acids and peptides''; in SIRIUS 6 the pathway, superclass and
#'     class columns are joined with ``; ''). \code{NA} if CANOPUS was not run
#'     or produced no result.}
#'   \item{canopus_superclass}{Top-level CANOPUS superclass (parsed from NPC).}
#'   \item{canopus_class}{Mid-level CANOPUS class (parsed from NPC).}
#'   \item{canopus_subclass}{Finest CANOPUS subclass. SIRIUS 6 does not emit an
#'     \code{NPC#subclass} column, so this is \code{NA} for SIRIUS 6 output
#'     - this is expected, not a parsing failure.}
#' }
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' \dontrun{
#' p <- create_nts_project("raw_data")
#' p <- extract_features(p, param = XcmsParam())
#' p <- extract_ms2(p, "ms2_data")
#' p <- compute_sirius_fingerprints(p, ms2_source = "ms2data",
#'   sirius_exe = "/opt/sirius/bin/sirius")
#' fp <- get_sirius_fingerprints(p)
#' head(fp[, c("name", "molecular_formula", "rank", "sirius_score")])
#' }
#' @export
compute_sirius_fingerprints <- function(project,
                                        ms2_source = c("clean_ms2", "ms2data"),
                                        output_dir = NULL,
                                        ionization = "[M-H]-",
                                        # Rare elements SIRIUS may auto-detect.
                                        # C/H/N/O/P are always considered; S must
                                        # stay in this set for PFAS sulfonates /
                                        # sulfonamides to be detected.
                                        elements = "SClBrSeI",
                                        sirius_exe = NULL,
                                        batch_size = 20,
                                        keep_intermediates = FALSE,
                                        run_canopus = FALSE,
                                        overwrite = FALSE,
                                        ...) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")
  ms2_source <- match.arg(ms2_source)

  # --- resolve SIRIUS path ---
  if (is.null(sirius_exe) || !nzchar(sirius_exe)) {
    sirius_exe <- Sys.getenv("NTSminer_SIRIUS_EXE", Sys.getenv("SIRIUS_EXE", ""))
  }
  if (!nzchar(sirius_exe)) {
    stop(
      "SIRIUS path not set: pass `sirius_exe = ` or set env NTSminer_SIRIUS_EXE ",
      "(or SIRIUS_EXE).",
      call. = FALSE
    )
  }
  if (!file.exists(sirius_exe)) {
    stop("SIRIUS executable not found at: ", sirius_exe, call. = FALSE)
  }

  # --- resolve MS2 source ---
  if (identical(ms2_source, "clean_ms2")) {
    ms2_list <- methods::slot(project, "clean_ms2")
    if (!length(ms2_list)) {
      stop("No purified MS2 data found (@clean_ms2 empty). Run purify_ms2_spectra() first, ",
           "or set ms2_source = \"ms2data\".", call. = FALSE)
    }
    feat_names <- names(ms2_list)
  } else {
    ms2_raw <- methods::slot(project, "ms2data")
    if (!length(ms2_raw)) {
      stop("No raw MS2 data found (@ms2data empty). Run extract_ms2() first.", call. = FALSE)
    }
    feat_names <- names(ms2_raw)
  }
  n_total <- length(feat_names)

  # --- resolve output directory ---
  if (is.null(output_dir)) {
    output_dir <- tempfile("nts_sirius_")
    on.exit(
      if (!keep_intermediates && dir.exists(output_dir)) {
        unlink(output_dir, recursive = TRUE)
      },
      add = TRUE
    )
  }
  output_dir_abs <- normalizePath(output_dir, mustWork = FALSE)
  .nts_sirius_prepare_output_dir(output_dir_abs, overwrite = overwrite)
  output_dir <- normalizePath(output_dir_abs, mustWork = TRUE)

  batch_size <- .nts_check_batch_size(batch_size, n_total)
  n_batches <- ceiling(n_total / batch_size)

  # Resolve `elements` once (empty when the flag should be omitted)
  elements_flag <- .nts_sirius_elements_flag(elements)

  # Per-batch output dirs to avoid overwriting
  batch_out_dirs <- character(n_batches)
  batch_status   <- integer(n_batches)   # exit status per batch (0 = ok)

  pb <- NULL
  is_shiny <- !is.null(shiny::getDefaultReactiveDomain())
  if (!is_shiny) {
    pb <- utils::txtProgressBar(min = 0, max = n_batches, style = 3)
  }

  for (i in seq_len(n_batches)) {
    idx_start <- (i - 1L) * batch_size + 1L
    idx_end   <- min(i * batch_size, n_total)
    batch_feats <- feat_names[idx_start:idx_end]

    msg <- sprintf(" [%d/%d] Batch %d: %s", idx_end, n_total, i,
                   paste(head(batch_feats, 3), collapse = ", "))
    if (length(batch_feats) > 3) msg <- paste0(msg, " ...")

    if (is_shiny) {
      shiny::incProgress(1 / n_batches, detail = msg)
    } else {
      cat("\n", msg, "\n")
      utils::setTxtProgressBar(pb, i)
    }

    # --- 1. write batch .ms file ---
    tmp_ms <- tempfile(pattern = paste0("sirius_batch_", i, "_"), fileext = ".ms")
    if (identical(ms2_source, "clean_ms2")) {
      # Single >ms2 per feature from clean_ms2
      sub_project <- project
      methods::slot(sub_project, "clean_ms2") <- ms2_list[batch_feats]
      export_to_sirius_ms(sub_project, tmp_ms, ionization = ionization)
    } else {
      # Multiple >ms2 blocks per feature from raw ms2data
      ms2_raw <- methods::slot(project, "ms2data")
      .write_ms2data_ms(ms2_raw[batch_feats], tmp_ms, ionization = ionization)
    }

    # --- 2. SIRIUS project directory (per batch, alongside batch output) ---
    batch_proj <- file.path(output_dir, paste0("sirius_proj_", i, "_batch.sirius"))

    # Per-batch output subdirectory - each batch writes its own summary file
    batch_out <- file.path(output_dir, paste0("batch_", i))
    dir.create(batch_out, recursive = TRUE, showWarnings = FALSE)
    batch_out_dirs[i] <- batch_out

    extra_args <- list(...)
    toolchain_args <- c(
      "--input", shQuote(normalizePath(tmp_ms)),
      "--project", shQuote(batch_proj),
      as.character(unlist(extra_args)),
      "formulas",
      elements_flag,
      "fingerprints",
      if (run_canopus) "canopus",
      "write-summaries",
      "--output", shQuote(batch_out)
    )

    # run SIRIUS for this batch -- capture output + exit status so failures are
    # never silently swallowed (previously `stdout = FALSE` discarded both).
    cli_log   <- file.path(output_dir, "sirius_cli.log")
    batch_log <- file.path(output_dir, sprintf("sirius_batch_%d.log", i))
    r_cli <- .nts_sirius_run_cli(
      sirius_exe, toolchain_args,
      log_file = batch_log, cli_log = cli_log,
      quiet = !is_shiny, echo = !is_shiny, on_fail = "warn",
      label = paste0("for batch ", i)
    )
    batch_status[i] <- r_cli$status  # kept for the post-loop failure summary

    # Optionally try CLI fingerprint export (only supported in SIRIUS 5.x)
    if (keep_intermediates && dir.exists(batch_proj)) {
      fp_export_dir <- file.path(output_dir, paste0("fp_export_", i))
      dir.create(fp_export_dir, recursive = TRUE, showWarnings = FALSE)
      fp_tried <- .nts_sirius_export_fingerprints(sirius_exe, batch_proj, fp_export_dir)
      if (!is.null(fp_tried)) {
        message(">> Batch ", i, ": CLI fingerprint export succeeded")
      }
    }

    # Cleanup
    unlink(tmp_ms)
    if (!keep_intermediates) {
      unlink(batch_proj, recursive = TRUE)
    }
  }

  if (!is_shiny) close(pb)
  cat("\n")

  # ======================================================================
  # 3. Merge per-batch CSV files + parse
  # ======================================================================

  all_meta <- list()
  for (i in seq_len(n_batches)) {
    bd <- batch_out_dirs[i]
    f <- file.path(bd, "formula_identifications.csv")
    if (!file.exists(f)) {
      f <- file.path(bd, "formula_identifications.tsv")
      sep_val2 <- "\t"
    } else {
      sep_val2 <- ","
    }
    if (!file.exists(f)) next
    tbl <- tryCatch(
      utils::read.table(f, header = TRUE, sep = sep_val2,
                        stringsAsFactors = FALSE, quote = "\"",
                        comment.char = "", check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(tbl) || !nrow(tbl)) next
    all_meta[[length(all_meta) + 1L]] <- tbl
  }

  n_failed  <- sum(batch_status != 0L)
  n_with_tsv <- length(all_meta)

  # A batch that exits 0 but writes no summary table produced nothing (SIRIUS
  # found no formula candidate). Counting only the exit code would report such a
  # run as fully successful.
  has_tsv <- vapply(batch_out_dirs, function(bd) {
    isTRUE(file.exists(file.path(bd, "formula_identifications.tsv")) ||
             file.exists(file.path(bd, "formula_identifications.csv")))
  }, logical(1L))
  empty_batches <- which(batch_status == 0L & !has_tsv)
  n_empty <- length(empty_batches)

  # Every batch failed -> abort loudly. (Guarded on n_batches > 0: an empty
  # project produces zero batches and must not look like a total failure.)
  if (n_batches > 0L && n_failed == n_batches) {
    stop("SIRIUS failed for all ", n_batches, " batch(es) (exit status: ",
         paste(unique(batch_status), collapse = ", "), "). Aborting; see logs in:\n  ",
         output_dir,
         call. = FALSE)
  }

  if (!length(all_meta)) {
    if (n_failed > 0L) {
      # Partial failure + nothing parsed: this used to return a silent empty
      # table, indistinguishable from "SIRIUS ran fine but found nothing".
      warning(n_failed, " of ", n_batches, " SIRIUS batch(es) failed (non-zero exit ",
              "status) and NO formula_identifications table was produced. ",
              "Returning an empty fingerprint table. See ",
              file.path(output_dir, "sirius_batch_*.log"), call. = FALSE)
    } else {
      warning("No formula_identifications summary found in any batch. ",
              "SIRIUS may have produced no results.", call. = FALSE)
    }
    # Must be a `data.table`: @sirius_fingerprints is declared `data.table`, and
    # a plain 0-row data.frame used to blow up in set_sirius_fingerprints() with
    # "is(value, \"data.table\") is not TRUE" -- masking the real warning above.
    empty_out <- data.table::data.table(
      name = character(0),
      adduct = character(0),
      molecular_formula = character(0),
      rank = integer(0),
      sirius_score = numeric(0),
      zodiac_score = numeric(0),
      confidence = numeric(0),
      explained_peaks = integer(0),
      fingerprint = list(),
      fp_version = character(0),
      canopus_npc = character(0),
      canopus_superclass = character(0),
      canopus_class = character(0),
      canopus_subclass = character(0)
    )
    project <- set_sirius_fingerprints(project, empty_out)
    methods::slot(project, "log") <- .nts_log_append(
      methods::slot(project, "log"),
      "compute_sirius_fingerprints",
      list(output_dir = output_dir,
           n_features = n_total,
           n_batches = n_batches,
           n_failed_batches = n_failed,
           n_empty_batches = n_empty,
           n_matched = 0L)
    )
    return(invisible(project))
  }

  if (n_failed > 0L) {
    warning(n_failed, " of ", n_batches, " SIRIUS batch(es) failed (non-zero exit ",
            "status); results are PARTIAL - the features in those batches have no ",
            "formula candidates. See ",
            file.path(output_dir, "sirius_batch_*.log"), call. = FALSE)
  }
  if (n_empty > 0L) {
    warning(n_empty, " of ", n_batches, " SIRIUS batch(es) exited 0 but produced no ",
            "formula_identifications table (no formula candidate found for those ",
            "features). batch(es): ", paste(utils::head(empty_batches, 10L), collapse = ", "),
            if (n_empty > 10L) " ..." else "", call. = FALSE)
  }

  # `rbindlist(fill = TRUE)` tolerates batches whose column sets differ (SIRIUS
  # version / flag differences). `do.call(rbind, ...)` would error out with
  # "numbers of columns of arguments do not match" and, worse, silently shift
  # values if a later batch gained a column.
  meta_table <- data.table::rbindlist(all_meta, fill = TRUE)

  # ======================================================================
  # 3b. Parse CANOPUS classifications (per-batch canopus_formula_summary.tsv)
  # ======================================================================
  all_canopus <- .nts_sirius_read_canopus(
    file.path(batch_out_dirs, "canopus_formula_summary.tsv")
  )

  # Rename the mapping column to "name" if needed
  # SIRIUS 6 may have an empty "name" column - prefer mappingFeatureId
  if ("mappingFeatureId" %in% names(meta_table)) {
    meta_table$name <- trimws(as.character(meta_table$mappingFeatureId))
  } else if ("name" %in% names(meta_table) && any(nzchar(trimws(meta_table$name)))) {
    meta_table$name <- trimws(as.character(meta_table$name))
  } else {
    id_col <- grep("^alignedFeatureId$|^featureId$|^compoundId$|^id$",
                   names(meta_table), value = TRUE)[1L]
    if (is.na(id_col)) id_col <- names(meta_table)[1L]
    names(meta_table)[names(meta_table) == id_col] <- "name"
    meta_table$name <- as.character(meta_table$name)
  }

  # Repair mappingFeatureId - SIRIUS 6 prefixes it with batch/tmp noise.
  # The true feature name is the suffix after the last underscore before "FT" or the
  # original name that appears in names(clean_ms2).
  # Strategy: if mappingFeatureId ends with an exact feat_names entry, use that;
  # otherwise try the longest suffix match.
  orig_feat_names <- names(methods::slot(project, "clean_ms2"))
  meta_table$name <- vapply(meta_table$name, .nts_sirius_repair_id,
                            character(1L), USE.NAMES = FALSE,
                            orig_feat_names = orig_feat_names)

  # `rbindlist(fill = TRUE)` leaves a column absent entirely when no batch had
  # it, and `meta_table[[col]]` then returns NULL instead of NA. Materialise the
  # columns this function relies on so downstream indexing stays length-safe.
  # NOTE: this must run *after* the name-resolution block above - inventing an
  # all-NA `mappingFeatureId` earlier would make that block always take the
  # `mappingFeatureId` branch and lose the fallbacks.
  for (cc in c("name", "adduct", "molecularFormula", "molecular_formula",
               "precursorFormula", "numExplainedPeaks", "explainedIntensity",
               "SiriusScore", "ZodiacScore", "rank", "Rank",
               "mappingFeatureId", "alignedFeatureId", "compoundId")) {
    if (!cc %in% names(meta_table)) meta_table[[cc]] <- NA
  }

  # --- build the output data.frame ---
  # We iterate over meta_table rows; for each row we try to attach its
  # fingerprint vector. The meta table may have one row per formula candidate
  # per feature (rank column distinguishes them).

  # Fingerprint vectors are populated from CLI export (--type=fingerprints)
  # collected in batch_fp above. For SIRIUS versions that don't support CLI
  # export, the fingerprint column will be empty and
  # extract_sirius_fingerprints_api() can be used as fallback.

  # Column detection - use exact SIRIUS 6 canonical names with fallback to fuzzy
  message(">> CSV column names: ", paste(names(meta_table), collapse = " | "))
  cn <- names(meta_table)

  # Resolve score columns (exact names first, then conservative fuzzy fallback).
  # `.nts_sirius_col_map()` also emits the one-off "no confidence column" warning.
  col_map <- .nts_sirius_col_map(cn)

  has_sirius   <- if (length(col_map$sirius_score)) col_map$sirius_score[1L] else NA_character_
  has_zodiac   <- if (length(col_map$zodiac_score)) col_map$zodiac_score[1L] else NA_character_
  has_conf     <- if (length(col_map$confidence)) col_map$confidence[1L] else NA_character_
  has_expl     <- if (length(col_map$explained_peaks)) col_map$explained_peaks[1L] else NA_character_
  rank_col     <- if (length(col_map$rank)) col_map$rank[1L] else NA_character_
  has_rank     <- !is.na(rank_col)

  # grab adduct and formula columns
  adduct_col <- grep("^adduct$", names(meta_table), value = TRUE, ignore.case = TRUE)[1L]
  formula_col <- grep("^molecularFormula$|^molecular_formula$|^formula$",
                      names(meta_table), value = TRUE, ignore.case = TRUE)[1L]

  # Validate that critical columns actually exist
  if (is.na(adduct_col) || is.null(meta_table[[adduct_col]])) {
    stop("Cannot find adduct column in SIRIUS formula_identifications output. ",
         "Available columns: ", paste(cn, collapse = ", "), call. = FALSE)
  }
  if (is.na(formula_col) || is.null(meta_table[[formula_col]])) {
    stop("Cannot find molecularFormula column in SIRIUS formula_identifications output. ",
         "Available columns: ", paste(cn, collapse = ", "), call. = FALSE)
  }

  out <- data.table::data.table(
    name = meta_table$name,
    adduct = as.character(meta_table[[adduct_col]]),
    molecular_formula = as.character(meta_table[[formula_col]]),
    rank = if (has_rank) as.integer(meta_table[[rank_col]]) else NA_integer_,
    sirius_score = if (!is.na(has_sirius)) as.numeric(meta_table[[has_sirius]]) else NA_real_,
    zodiac_score = if (!is.na(has_zodiac)) as.numeric(meta_table[[has_zodiac]]) else NA_real_,
    confidence = if (!is.na(has_conf)) as.numeric(meta_table[[has_conf]]) else NA_real_,
    explained_peaks = if (!is.na(has_expl)) as.integer(meta_table[[has_expl]]) else NA_integer_
  )

  # --- Merge CANOPUS NPC classifications ---
  if (length(all_canopus)) {
    canopus_dt <- data.table::rbindlist(all_canopus, fill = TRUE)
    # Standardize column names expected from SIRIUS CANOPUS output
    if ("formularank" %in% names(canopus_dt)) {
      data.table::setnames(canopus_dt, "formularank", "rank")
    }
    if ("molecularformula" %in% names(canopus_dt)) {
      data.table::setnames(canopus_dt, "molecularformula", "molecular_formula")
    }
    if (!all(c("molecular_formula", "adduct", "rank") %in% names(canopus_dt))) {
      warning("CANOPUS summary is missing one of molecular_formula/adduct/rank; ",
              "cannot merge. Columns: ",
              paste(utils::head(names(canopus_dt), 20L), collapse = ", "), call. = FALSE)
      out[, canopus_npc := NA_character_]
      out[, canopus_superclass := NA_character_]
      out[, canopus_class := NA_character_]
      out[, canopus_subclass := NA_character_]
    } else {
      out <- .nts_sirius_merge_canopus(out, canopus_dt, orig_feat_names)
    }
  } else {
    if (run_canopus) {
      warning("CANOPUS was enabled but no canopus_formula_summary.tsv found in any batch output.",
              call. = FALSE)
    }
    out[, canopus_npc := NA_character_]
    out[, canopus_superclass := NA_character_]
    out[, canopus_class := NA_character_]
    out[, canopus_subclass := NA_character_]
  }

  # One-off visibility line: make it obvious which score columns are actually
  # populated (a silent all-NA column is the failure mode this guards against).
  message(sprintf(
    ">> sirius_score/zodiac_score/confidence/explained_peaks: NA counts = %d/%d/%d/%d (of %d rows)",
    sum(is.na(out$sirius_score)), sum(is.na(out$zodiac_score)),
    sum(is.na(out$confidence)), sum(is.na(out$explained_peaks)), nrow(out)))

  # Attach fingerprint list-column
  # SIRIUS 6.x does NOT export fingerprint vectors via write-summaries.
  # Fingerprints are stored in the .sirius project file and accessible
  # through the REST API only. Use extract_sirius_fingerprints_api()
  # after ensuring sirius login has been completed.
  fp_col <- vector("list", nrow(out))
  fpv_col <- character(nrow(out))
  for (k in seq_len(nrow(out))) {
    fp_col[[k]] <- numeric(0)
  }
  fpv_col[] <- ""
  out[, fingerprint := fp_col]
  out[, fp_version := fpv_col]

  # keep as data.table for the slot
  project <- set_sirius_fingerprints(project, out)

  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "compute_sirius_fingerprints",
    list(
      output_dir = output_dir,
      n_features = n_total,
      n_batches = n_batches,
      n_failed_batches = n_failed,
      n_empty_batches = n_empty,
      n_candidates = nrow(out)
    )
  )

  message(">> compute_sirius_fingerprints: ",
          nrow(out), " candidate rows (", n_total, " features)" )
  message(">> Fingerprint vectors not populated. To fill them:")
  message(">>   1. Run: sirius login  (in terminal)")
  message(">>   2. Run: extract_sirius_fingerprints_api(project, ...)")
  message(">> Fingerprints stored in @sirius_fingerprints. Use get_sirius_fingerprints() to access.")

  invisible(project)
}


# =========================================================================
# CANOPUS standalone parser
# =========================================================================

#' Parse CANOPUS classification results
#'
#' Reads CANOPUS NPC pathway classification output from SIRIUS batch
#' directories and merges the results into \code{@sirius_fingerprints}.
#'
#' Use this function when CANOPUS was already run (via
#' \code{compute_sirius_fingerprints(run_canopus = TRUE)}) but you want
#' to re-parse the classification results without re-running SIRIUS,
#' or when CANOPUS was run externally against the same batch output
#' directories.
#'
#' Expected file per batch: \code{canopus_formula_summary.tsv}. The NPC
#' classification columns are resolved by exact name and may take two shapes:
#' SIRIUS 6 writes separate \code{NPC#pathway} / \code{NPC#superclass} /
#' \code{NPC#class} columns (concatenated into \code{canopus_npc} with
#' \code{"; "}); older SIRIUS writes a single \code{NPC} column holding a
#' \code{"a: b: c"} string, which is split into the three levels using
#' \verb{":"} / \verb{";"} / \verb{"/"} as delimiters. SIRIUS 6 has no
#' \code{NPC#subclass}, so \code{canopus_subclass} is \code{NA} there by
#' design. If no NPC column can be found, the four columns are filled with
#' \code{NA} and a warning lists the columns actually present.
#'
#' Merging is keyed on the feature identity (\code{mappingFeatureId}, repaired
#' to the \code{@clean_ms2} name) together with
#' \code{(molecular_formula, adduct, rank)}, so two features that share a
#' formula do not inherit each other's classification. If the CANOPUS table
#' carries no \code{mappingFeatureId}, a three-column key is used and a warning
#' is emitted. Rows without a match receive \code{NA} in all four columns.
#'
#' @param project An [NTSProject] with existing \code{@sirius_fingerprints}
#'   (typically from \code{compute_sirius_fingerprints()}).
#' @param output_dir The \code{output_dir} used in
#'   \code{compute_sirius_fingerprints()}, containing per-batch
#'   subdirectories (\code{batch_1}, \code{batch_2}, ...) with
#'   \code{canopus_formula_summary.tsv} files.
#'
#' @return The updated \code{project} with CANOPUS columns populated
#'   in \code{@sirius_fingerprints}.
#'
#' @section Slot columns added/updated:
#' \describe{
#'   \item{canopus_npc}{Raw CANOPUS NPC pathway string (character).}
#'   \item{canopus_superclass}{First level of NPC pathway (character).}
#'   \item{canopus_class}{Second level of NPC pathway (character).}
#'   \item{canopus_subclass}{Third level of NPC pathway (character).}
#' }
#'
#' @examplesIf requireNamespace("data.table", quietly = TRUE)
#' \dontrun{
#' p <- create_nts_project("raw_data")
#' p <- extract_features(p, param = XcmsParam())
#' p <- extract_ms2(p, "ms2_data")
#' p <- compute_sirius_fingerprints(p, output_dir = "sirius_out")
#' p <- parse_canopus_results(p, output_dir = "sirius_out")
#' head(get_sirius_fingerprints(p)[, c("name", "canopus_npc", "canopus_class")])
#' }
#' @export
parse_canopus_results <- function(project, output_dir) {
  .nts_check(project, "data.table_ns", "project", .arg = "project")

  fp <- methods::slot(project, "sirius_fingerprints")
  if (!data.table::is.data.table(fp) || !nrow(fp)) {
    stop(
      "@sirius_fingerprints is empty. Run compute_sirius_fingerprints() first.",
      call. = FALSE
    )
  }

  if (!dir.exists(output_dir)) {
    stop("output_dir does not exist: ", output_dir, call. = FALSE)
  }

  # Discover batch subdirectories
  batch_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = TRUE)
  batch_dirs <- batch_dirs[grepl("batch_\\d+$", basename(batch_dirs))]
  if (!length(batch_dirs)) {
    stop(
      "No batch_* subdirectories found in output_dir. ",
      "Ensure output_dir matches the one used in compute_sirius_fingerprints().",
      call. = FALSE
    )
  }

  # Parse CANOPUS TSV from each batch
  all_canopus <- .nts_sirius_read_canopus(
    file.path(batch_dirs, "canopus_formula_summary.tsv")
  )

  if (!length(all_canopus)) {
    warning(
      "No canopus_formula_summary.tsv found in any batch subdirectory. ",
      "Ensure CANOPUS was run (compute_sirius_fingerprints(run_canopus = TRUE)).",
      call. = FALSE
    )
    fp[, canopus_npc := NA_character_]
    fp[, canopus_superclass := NA_character_]
    fp[, canopus_class := NA_character_]
    fp[, canopus_subclass := NA_character_]
    project <- set_sirius_fingerprints(project, fp)
    return(invisible(project))
  }

  canopus_dt <- data.table::rbindlist(all_canopus, fill = TRUE)
  message(">> Parsed CANOPUS results: ", nrow(canopus_dt), " rows from ",
          length(all_canopus), " batch(es)")

  # Standardize key column names
  if ("formularank" %in% names(canopus_dt)) {
    data.table::setnames(canopus_dt, "formularank", "rank")
  }
  if ("molecularformula" %in% names(canopus_dt)) {
    data.table::setnames(canopus_dt, "molecularformula", "molecular_formula")
  }

  if (!all(c("molecular_formula", "adduct", "rank") %in% names(canopus_dt))) {
    warning("CANOPUS TSV missing key columns (molecular_formula/adduct/rank); cannot merge. Columns: ",
            paste(utils::head(names(canopus_dt), 20L), collapse = ", "),
            call. = FALSE)
    fp[, canopus_npc := NA_character_]
    fp[, canopus_superclass := NA_character_]
    fp[, canopus_class := NA_character_]
    fp[, canopus_subclass := NA_character_]
    project <- set_sirius_fingerprints(project, fp)
    return(invisible(project))
  }

  # Same NPC resolution / merge as compute_sirius_fingerprints() (shared helper)
  orig_feat_names <- names(methods::slot(project, "clean_ms2"))

  # Remove existing CANOPUS columns before merge (to avoid duplicates / stale
  # classifications surviving a re-parse).
  existing_canopus_cols <- intersect(
    c("canopus_npc", "canopus_superclass", "canopus_class", "canopus_subclass"),
    names(fp)
  )
  if (length(existing_canopus_cols)) {
    fp[, (existing_canopus_cols) := NULL]
  }

  fp <- .nts_sirius_merge_canopus(fp, canopus_dt, orig_feat_names)
  message(">> Merged CANOPUS NPC classifications: ",
          sum(!is.na(fp$canopus_npc)), " of ", nrow(fp), " rows classified.")

  project <- set_sirius_fingerprints(project, fp)

  methods::slot(project, "log") <- .nts_log_append(
    methods::slot(project, "log"),
    "parse_canopus_results",
    list(
      output_dir = output_dir,
      n_canopus_rows = nrow(canopus_dt),
      n_classified = sum(!is.na(fp$canopus_npc))
    )
  )

  invisible(project)
}
