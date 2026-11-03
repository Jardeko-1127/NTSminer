#' Robust BiocParallel backend for xcms / heavy jobs
#'
#' On **RStudio** (`RSTUDIO=1`) or **Windows**, uses [BiocParallel::SnowParam]
#' with explicit `tasks` chunking to reduce memory spikes; on Unix non-RStudio
#' interactive shells, uses [BiocParallel::MulticoreParam] with matching
#' `tasks`.
#'
#' @param threads Integer, number of parallel workers (coerced to `>= 1`).
#'
#' @return A [BiocParallel::BiocParallelParam] object.
#' @export
get_robust_bpparam <- function(threads = 1L) {
  threads <- max(1L, as.integer(threads)[1L])
  if (threads == 1L) {
    return(BiocParallel::SerialParam())
  }
  tasks <- max(threads * 4L, 8L)
  use_snow <- Sys.getenv("RSTUDIO") == "1" || .Platform$OS.type == "windows"
  if (use_snow) {
    BiocParallel::SnowParam(
      workers = threads,
      type = "SOCK",
      tasks = tasks,
      exportglobals = FALSE,
      exportvariables = FALSE,
      log = FALSE
    )
  } else {
    BiocParallel::MulticoreParam(
      workers = threads,
      tasks = tasks,
      exportglobals = FALSE,
      exportvariables = FALSE,
      log = FALSE
    )
  }
}


#' Assert NTSminer is library-installed for SnowParam / MulticoreParam workers.
#' @param internal_names Character vector of required unexported symbols.
#' @keywords internal
#' @noRd
.nts_assert_parallel_ntsminer <- function(internal_names = character(0)) {
  installed <- tryCatch(
    utils::packageDescription("NTSminer", lib.loc = .libPaths()[1L]),
    error = function(e) NULL
  )
  if (is.null(installed)) {
    stop(
      "Parallel workers require NTSminer installed in the R library ",
      "(remotes::install_local()). ",
      "devtools::load_all() / pkgload does not propagate to SnowParam workers on Windows.",
      call. = FALSE
    )
  }
  if (!requireNamespace("NTSminer", quietly = TRUE)) {
    stop("Failed to load installed NTSminer for parallel workers.", call. = FALSE)
  }
  ns <- asNamespace("NTSminer")
  if (length(internal_names)) {
    missing <- internal_names[!vapply(
      internal_names,
      exists,
      logical(1L),
      envir = ns,
      inherits = FALSE
    )]
    if (length(missing)) {
      stop(
        "Installed NTSminer is missing ",
        paste(missing, collapse = ", "),
        ". Reinstall with remotes::install_local().",
        call. = FALSE
      )
    }
  }
  loaded_ver <- tryCatch(
    as.character(utils::packageVersion("NTSminer")),
    error = function(e) NA_character_
  )
  installed_ver <- as.character(installed$Version)
  if (!is.na(loaded_ver) && loaded_ver != installed_ver) {
    warning(
      "Loaded NTSminer (", loaded_ver, ") differs from installed (",
      installed_ver, "). Parallel workers use the installed copy; ",
      "run remotes::install_local() and restart R.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


#' Call an unexported NTSminer helper from a BiocParallel worker.
#' @keywords internal
#' @noRd
.nts_parallel_call_internal <- function(name, ...) {
  .nts_assert_parallel_ntsminer(name)
  fn <- getFromNamespace(name, "NTSminer")
  fn(...)
}


#' Monoisotopic masses of the most abundant isotope (Da)
#'
#' Built-in fallback used by [calculate_exact_mass] when the optional `Rdisop`
#' package is not installed, so exact-mass (KMD / neutral-loss / homologue)
#' logic never depends on an optional package being present.
#'
#' @keywords internal
#' @noRd
.NTS_MONOISOTOPIC_MASS <- c(
  H = 1.00782503207, D = 2.0141017778, He = 4.00260325415,
  Li = 7.01600455, B = 11.0093054, C = 12.00000000000,
  N = 14.0030740048, O = 15.9949146196, F = 18.99840322,
  Na = 22.9897692809, Mg = 23.98504170, Al = 26.98153863,
  Si = 27.9769265325, P = 30.97376163, S = 31.97207100,
  Cl = 34.96885268, Ar = 39.9623831225, K = 38.96370668,
  Ca = 39.96259098, Cr = 51.9405075, Mn = 54.9380451,
  Fe = 55.9349375, Co = 58.9331950, Ni = 57.9353429,
  Cu = 62.9295975, Zn = 63.9291422, As = 74.9215965,
  Se = 79.9165213, Br = 78.9183371, Mo = 97.9054048,
  Ag = 106.9050970, Cd = 113.9033585, Sn = 119.9021947,
  Sb = 120.9038157, I = 126.9044730, Ba = 137.9052472,
  Hg = 201.9706430, Pb = 207.9766521
)

#' Parse a Hill-style formula into element counts
#'
#' Handles element symbols, optional counts, and one or more levels of
#' parenthesised groups with an optional multiplier, e.g.
#' `"C3HF5O2"`, `"C24H19F34N2O8PS2"`, `"(CF2)3"`.
#'
#' @param formula Formula string.
#' @return A named integer vector of element counts, or `NULL` when the string
#'   cannot be parsed.
#' @keywords internal
#' @noRd
.nts_formula_counts <- function(formula) {
  f <- gsub("[[:space:]]+", "", as.character(formula)[1L])
  if (is.na(f) || !nzchar(f)) return(NULL)

  rec <- function(s) {
    out <- integer(0)
    i <- 1L
    n <- nchar(s)
    while (i <= n) {
      ch <- substr(s, i, i)
      if (identical(ch, "(")) {
        depth <- 0L
        j <- i
        repeat {
          cj <- substr(s, j, j)
          if (identical(cj, "(")) depth <- depth + 1L
          if (identical(cj, ")")) {
            depth <- depth - 1L
            if (depth == 0L) break
          }
          j <- j + 1L
          if (j > n) return(NULL)
        }
        sub_cnt <- rec(substr(s, i + 1L, j - 1L))
        if (is.null(sub_cnt)) return(NULL)
        i <- j + 1L
        num <- ""
        while (i <= n && grepl("[0-9]", substr(s, i, i))) {
          num <- paste0(num, substr(s, i, i)); i <- i + 1L
        }
        mult <- if (nzchar(num)) as.integer(num) else 1L
        for (nm in names(sub_cnt)) {
          add <- sub_cnt[[nm]] * mult
          out[nm] <- (if (nm %in% names(out)) out[[nm]] else 0L) + add
        }
      } else if (grepl("[A-Z]", ch)) {
        sym <- ch
        i <- i + 1L
        while (i <= n && grepl("[a-z]", substr(s, i, i))) {
          sym <- paste0(sym, substr(s, i, i)); i <- i + 1L
        }
        num <- ""
        while (i <= n && grepl("[0-9]", substr(s, i, i))) {
          num <- paste0(num, substr(s, i, i)); i <- i + 1L
        }
        cnt <- if (nzchar(num)) as.integer(num) else 1L
        out[sym] <- (if (sym %in% names(out)) out[[sym]] else 0L) + cnt
      } else {
        return(NULL)
      }
    }
    out
  }

  res <- rec(f)
  if (!length(res)) NULL else res
}

# one-time notice when the built-in mass table substitutes for Rdisop
.nts_extmass_env <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.nts_note_builtin_mass <- function() {
  if (isTRUE(.nts_extmass_env$noted)) return(invisible(NULL))
  .nts_extmass_env$noted <- TRUE
  message(
    "NTSminer: package 'Rdisop' is not installed -- exact masses are computed ",
    "from NTSminer's built-in monoisotopic mass table. Results agree with ",
    "Rdisop to <1e-6 Da; install Rdisop (BiocManager::install(\"Rdisop\")) for ",
    "bit-identical values to earlier releases."
  )
  invisible(NULL)
}

#' Exact monoisotopic mass from a Hill-style formula
#'
#' Returns the monoisotopic exact mass of a Hill-style formula, used by the KMD
#' and neutral-loss (ISF) logic instead of hard-coded mass constants.
#'
#' When the optional [Rdisop] package is installed it is consulted first, so
#' values stay bit-identical to earlier NTSminer releases. When it is absent the
#' mass is computed from NTSminer's built-in monoisotopic mass table
#' (`.NTS_MONOISOTOPIC_MASS`), so KMD / ISF / homologue stepping keep working on
#' a minimal install instead of returning `NA_real_`.
#'
#' @param formula Target formula string, e.g. `"CF2"`, `"CH2"`, `"H2O"`.
#'
#' @return A non-negative finite scalar; `NA_real_` with `warning` if the
#'   formula cannot be parsed by either engine.
#' @export
calculate_exact_mass <- function(formula) {
  if (missing(formula)) {
    return(NA_real_)
  }
  f <- tryCatch(as.character(formula)[1L], error = function(e) NA_character_)
  if (is.na(f) || !nzchar(trimws(f))) {
    return(NA_real_)
  }
  f <- trimws(f)
  if (nchar(f, type = "bytes", allowNA = TRUE) > 512L) {
    warning("Formula string too long; skipping parse.", call. = FALSE)
    return(NA_real_)
  }

  # 1) Rdisop first (keeps historical values bit-identical)
  if (requireNamespace("Rdisop", quietly = TRUE)) {
    out <- tryCatch(
      {
        mol <- Rdisop::getMolecule(f)
        em <- mol$exactmass
        if (is.null(em) || !length(em)) NA_real_ else suppressWarnings(as.numeric(em[1L]))
      },
      error = function(e) NA_real_
    )
    if (is.finite(out) && out > 0) {
      return(out)
    }
  } else {
    .nts_note_builtin_mass()
  }

  # 2) built-in monoisotopic mass table
  cnt <- .nts_formula_counts(f)
  if (!is.null(cnt)) {
    unknown <- setdiff(names(cnt), names(.NTS_MONOISOTOPIC_MASS))
    if (!length(unknown)) {
      return(as.numeric(sum(cnt * .NTS_MONOISOTOPIC_MASS[names(cnt)])))
    }
    warning(sprintf(
      "No monoisotopic mass for element(s) %s in formula '%s'.",
      paste(unknown, collapse = ", "), f
    ), call. = FALSE)
    return(NA_real_)
  }

  warning(sprintf("Could not parse formula: %s", f), call. = FALSE)
  NA_real_
}


#' @keywords internal
#' @noRd
.nts_check_flag <- function(x, flag, arg) {
  flag <- trimws(as.character(flag)[1L])
  switch(
    flag,
    "data.table_ns" = {
      if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("data.table is required", call. = FALSE)
      }
    },
    "project" = ,
    "nts_project" = {
      if (!methods::is(x, "NTSProject")) {
        stop(sprintf("`%s` must be an NTSProject", arg), call. = FALSE)
      }
    },
    "ntsdatabase" = ,
    "reference_db" = {
      if (!methods::is(x, "NTSDatabase")) {
        stop(sprintf("`%s` must be an NTSDatabase object", arg), call. = FALSE)
      }
    },
    "non_empty" = {
      if (is.null(x)) {
        stop(sprintf("`%s` must be non-NULL", arg), call. = FALSE)
      }
      if (is.data.frame(x) || data.table::is.data.table(x)) {
        if (!nrow(x)) {
          stop(sprintf("`%s` must be non-empty", arg), call. = FALSE)
        }
      } else if (is.list(x)) {
        if (!length(x)) {
          stop(sprintf("`%s` must be a non-empty list", arg), call. = FALSE)
        }
      } else {
        stop(
          sprintf(
            "`.nts_check(..., \"non_empty\")` for `%s`: expected data.frame, data.table, or list",
            arg
          ),
          call. = FALSE
        )
      }
    },
    "name_col" = {
      if (!is.data.frame(x) && !data.table::is.data.table(x)) {
        stop(sprintf("`%s` must be a data.frame or data.table", arg), call. = FALSE)
      }
      if (!"name" %in% names(x)) {
        stop(sprintf("`%s` must contain column `name`", arg), call. = FALSE)
      }
    },
    "ms2" = {
      if (!methods::is(x, "NTSProject")) {
        stop(sprintf("`%s` must be an NTSProject", arg), call. = FALSE)
      }
      m2 <- get_ms2_annotations(x)
      if (!is.data.frame(m2) || !nrow(m2) || !"name" %in% names(m2)) {
        stop(
          sprintf("`%s@ms2_annotations` must be a non-empty data.frame with column `name`", arg),
          call. = FALSE
        )
      }
    },
    "file_exists" = {
      p <- as.character(x)
      if (!length(p)) {
        stop(sprintf("`%s` must be a character file path", arg), call. = FALSE)
      }
      p <- unique(trimws(p))
      p <- p[nzchar(p) & !is.na(p)]
      if (!length(p)) {
        stop(sprintf("`%s` must be a non-blank file path", arg), call. = FALSE)
      }
      miss <- unique(p[!file.exists(p)])
      if (length(miss)) {
        stop(
          sprintf(
            "`%s` must exist on disk; missing: %s",
            arg,
            paste(utils::head(miss, 5L), collapse = "; ")
          ),
          call. = FALSE
        )
      }
    },
    "igraph_ns" = {
      if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("Package 'igraph' is required.", call. = FALSE)
      }
    },
    "positive_finite" = {
      if (is.null(x)) {
        stop(sprintf("`%s` must be non-NULL", arg), call. = FALSE)
      }
      v <- suppressWarnings(as.numeric(x)[1L])
      if (!length(x)) {
        stop(sprintf("`%s` must have length >= 1", arg), call. = FALSE)
      }
      if (!is.finite(v) || v <= 0) {
        stop(sprintf("`%s` must be a finite number > 0", arg), call. = FALSE)
      }
    },
    stop(sprintf("Unknown `.nts_check()` flag: %s", flag), call. = FALSE)
  )
}


#' Standard parameter checks (internal)
#'
#' Runs a sequence of guards identified by `flag` strings. Use at the top of
#' exported helpers to keep `if` / `stop` noise out of the main body.
#'
#' @param x Object to validate (ignored for namespace-only flags such as
#'   `"data.table_ns"` or `"igraph_ns"`).
#' @param ... Character strings: `"project"`, `"nts_project"`, `"ntsdatabase"`,
#'   `"reference_db"`, `"data.table_ns"`, `"igraph_ns"`, `"file_exists"`,
#'   `"positive_finite"`, `"non_empty"`, `"name_col"`, `"ms2"`.
#'   Multiple flags may be passed; `"data.table_ns"` is always evaluated first
#'   when present.
#' @param .arg Label used in error messages (without backticks in the template).
#'
#' @return `x` invisibly.
#' @keywords internal
#' @noRd
.nts_check <- function(x, ..., .arg = "argument") {
  fl <- unlist(list(...), recursive = TRUE, use.names = FALSE)
  fl <- as.character(fl)
  fl <- fl[!is.na(fl) & nzchar(fl)]
  if (!length(fl)) {
    return(invisible(x))
  }
  arg <- as.character(.arg)[1L]
  prio <- c(
    "data.table_ns",
    "igraph_ns",
    "project",
    "nts_project",
    "ntsdatabase",
    "reference_db",
    "file_exists",
    "positive_finite",
    "non_empty",
    "name_col",
    "ms2"
  )
  ord <- match(fl, prio, nomatch = 999L)
  fl <- fl[order(ord, seq_along(fl))]
  for (f in fl) {
    .nts_check_flag(x, f, arg)
  }
  invisible(x)
}

#' Preset adduct ion lists for database matching
#'
#' Convenience constants for the `adducts` parameter in
#' [batch_annotate] and [annotate_by_spectral_match].
#' Use the `*_ALL` variants to include all common adducts,
#' or pick individual tokens for a custom combination.
#'
#' @format A character vector of adduct tokens.
#' @name adduct_presets
#' @aliases ADDUCTS_POS ADDUCTS_POS_ALL ADDUCTS_NEG ADDUCTS_NEG_ALL
#'
#' @examples
#' # Match only with [M+H]+ (default for positive mode)
#' ADDUCTS_POS
#'
#' # Match with all common positive-mode adducts
#' ADDUCTS_POS_ALL
#'
#' # Match with [M-H]- and [M+Cl]- in negative mode
#' ADDUCTS_NEG_ALL
#'
#' # Custom combination: only M+H and M+Na
#' c("M+H", "M+Na")
NULL

#' @rdname adduct_presets
#' @export
ADDUCTS_POS <- c("M+H")

#' @rdname adduct_presets
#' @export
ADDUCTS_POS_ALL <- c("M+H", "M+Na", "M+K", "M+NH4")

#' @rdname adduct_presets
#' @export
ADDUCTS_NEG <- c("M-H")

#' @rdname adduct_presets
#' @export
ADDUCTS_NEG_ALL <- c("M-H", "M+Cl", "M-H-H2O")


#' Resolve an external executable across platforms
#'
#' Uniform resolution so callers stop hard-coding machine-specific Windows
#' paths. Lookup order: `envvar` (when set and non-empty), then `Sys.which(name)`
#' (PATH), then the first existing entry of `candidates`.
#'
#' @param name Bare executable name looked up on `PATH` (e.g. `"msconvert"`).
#' @param candidates Character vector of full paths to probe, in order.
#' @param envvar Optional environment-variable name whose value, when set,
#'   is used verbatim (lets users point at a custom install without editing code).
#' @return A length-1 character path/command, or `NULL` when nothing resolves.
#' @keywords internal
#' @noRd
.nts_resolve_executable <- function(name, candidates = character(0),
                                    envvar = NULL) {
  if (!is.null(envvar)) {
    ev <- trimws(Sys.getenv(envvar, unset = ""))
    if (nzchar(ev)) return(ev)
  }
  hit <- Sys.which(name)
  if (nzchar(hit)) return(unname(hit))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  for (p in candidates) {
    if (file.exists(p)) return(p)
  }
  NULL
}
