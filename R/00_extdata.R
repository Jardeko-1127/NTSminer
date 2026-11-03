# =============================================================================
# 内置数据（inst/extdata）统一查找链 + 常用 PFAS 谱库 + 合并助手
# -----------------------------------------------------------------------------
# 背景：包内自带的表 / 模型 / 谱库散布在 inst/extdata 与（被 .Rbuildignore 排除的）
#   顶层 data/ 里，历史上代码用 system.file("data", ...) 去取，装上包后恒取不到。
#   这里提供**单一查找链**，让所有「默认指向包内文件」的解析点都走同一入口，
#   并在源码树（pkgload::load_all / setwd(包根)）下也能正确回退。
#
# 约定（相对 <pkg>/extdata 的路径）：
#   libraries/   常用 PFAS 谱库（内置）
#   tables/      专家诊断碎片表 / 三库并集去重 suspect 清单
#                （NTS_PFAS_suspect_3lib_dedup.csv：inhouse+online+insilico 按
#                 InChIKey 去重的化合物表，供两个网络传播与 MetFrag LocalCSV 复用）
#   models/      logIE / RT 模型
# =============================================================================


#' Resolve a bundled `inst/extdata` file across install / source-tree layouts
#'
#' NTSminer ships a handful of tables, models and reference libraries under
#' `inst/extdata`. When the package is *installed*, these live at
#' `<pkg>/extdata/<rel>`; when running from the **source tree** (RStudio,
#' `pkgload::load_all()`), they are under `<pkg-root>/inst/extdata/<rel>`, and a
#' bare `<cwd>/extdata/<rel>` is also accepted. `system.file("data", ...)` (the
#' historical call site) never resolves because the top-level `data/` directory
#' is excluded from the build by `.Rbuildignore`, so everything routes through
#' this helper instead.
#'
#' @param rel Character vector of candidate paths relative to `extdata/`. The
#'   **first** one that exists is returned, allowing new ("tables/diagnostic_fragments.csv")
#'   and legacy ("diagnostic_fragments.csv") layouts to coexist.
#' @param package Package name (default `"NTSminer"`).
#' @param mustWork If `TRUE`, raise an informative error when nothing is found;
#'   otherwise return `NA_character_` so callers can degrade gracefully.
#' @return Normalised absolute path (length 1), or `NA_character_`.
#' @keywords internal
#' @noRd
.nts_extdata <- function(rel, package = "NTSminer", mustWork = FALSE) {
  rel <- as.character(rel)
  if (!length(rel)) {
    if (mustWork) stop("`.nts_extdata()` needs at least one candidate path.", call. = FALSE)
    return(NA_character_)
  }
  roots <- c(
    system.file("extdata", package = package),
    file.path("inst", "extdata"),
    file.path("extdata")
  )
  roots <- roots[nzchar(roots)]
  cand <- as.vector(outer(roots, rel, function(root, r) file.path(root, r)))
  cand <- cand[nzchar(cand) & file.exists(cand)]
  if (!length(cand)) {
    if (mustWork) {
      stop(sprintf(
        paste0("Bundled data not found: extdata/%s\n",
               "  expected (installed): <library>/%s/extdata/%s\n",
               "  expected (source tree): <pkg-root>/inst/extdata/%s\n",
               "  If running from source, run from the package root or use pkgload::load_all()."),
        rel[1L], package, rel[1L], rel[1L]), call. = FALSE)
    }
    return(NA_character_)
  }
  normalizePath(cand[1L], winslash = "/", mustWork = TRUE)
}


# -----------------------------------------------------------------------------
# 常用 PFAS 谱库（内置）
# -----------------------------------------------------------------------------

# 内置库的固定文件名（相对 extdata）。db_type 由库内逐行 db_type 决定置信度档：
#   inhouse -> L1 ; online -> L2a ; prediction -> L2b
.NTS_BUILTIN_LIBRARIES <- c(
  inhouse    = "libraries/NTS_InHouse_PFAS.rds",
  online     = "libraries/NTS_Online_PFAS_OECD.rds",
  prediction = "libraries/NTS_Pred_PFAS_OECD.rds"
)


#' Paths to the PFAS libraries bundled with NTSminer
#'
#' NTSminer ships three commonly used PFAS reference libraries so a fresh
#' install can annotate out of the box: an in-house standard library
#' (`inhouse`, confidence L1), an OECD online library (`online`, L2a) and an
#' OECD in-silico/prediction library (`prediction`, L2b). The very large
#' community libraries (MoNA-unified ~15 MB, spectraverse ~600 MB) are **not**
#' bundled.
#'
#' @param which Character subset of `c("inhouse", "online", "prediction")`;
#'   `NULL` (default) returns all three.
#' @return A named character vector of absolute `.rds` paths.
#' @seealso [load_builtin_pfas_library()], [merge_nts_libraries()]
#' @examples
#' \donttest{
#' nts_builtin_library_paths()
#' }
#' @export
nts_builtin_library_paths <- function(which = NULL) {
  all <- .NTS_BUILTIN_LIBRARIES
  if (is.null(which)) which <- names(all)
  which <- match.arg(which, names(all), several.ok = TRUE)
  paths <- vapply(all[which], function(r) .nts_extdata(r), character(1L))
  if (anyNA(paths)) {
    stop(sprintf(
      paste0("Bundled reference librar(y/ies) not found: %s\n",
             "  These libraries ship with the package under <library>/NTSminer/extdata/libraries/.\n",
             "  If running from source, check inst/extdata/libraries/ for the matching .rds files."),
      paste(which[is.na(paths)], collapse = ", ")), call. = FALSE)
  }
  paths
}


#' Load the bundled PFAS library (optionally combining the three)
#'
#' @param which Character subset of `c("inhouse", "online", "prediction")`;
#'   `NULL` (default) loads all three.
#' @param combine Logical. When `TRUE` (default) and more than one library is
#'   requested, the selected libraries are merged into a **single**
#'   [NTSDatabase] via [merge_nts_libraries()], so one call to
#'   `batch_annotate()` can reach L1/L2a/L2b at once (each row keeps its own
#'   `db_type`). When exactly one library is requested, that library object is
#'   returned unchanged (no `Lab_ID` prefixing). When `FALSE`, a named `list`
#'   of [NTSDatabase]s is returned.
#' @return An [NTSDatabase]; a named `list` of [NTSDatabase]s when
#'   `combine = FALSE` and `which` has length > 1.
#' @seealso [nts_builtin_library_paths()], [merge_nts_libraries()]
#' @examples
#' \donttest{
#' db <- load_builtin_pfas_library("inhouse")
#' db_all <- load_builtin_pfas_library()          # merged (L1 + L2a + L2b)
#' lst <- load_builtin_pfas_library(combine = FALSE)
#' }
#' @export
load_builtin_pfas_library <- function(which = NULL, combine = TRUE) {
  paths <- nts_builtin_library_paths(which)
  dbs <- lapply(paths, readRDS)
  names(dbs) <- names(paths)
  for (nm in names(dbs)) {
    if (!methods::is(dbs[[nm]], "NTSDatabase")) {
      stop("bundled library is not an NTSDatabase object: ", nm, " (", paths[[nm]], ")", call. = FALSE)
    }
  }
  if (length(dbs) == 1L) return(dbs[[1L]])          # 单库：原样返回，不加前缀
  if (combine) return(merge_nts_libraries(dbs))     # 多库：默认并成一个
  dbs                                               # 多库且 combine=FALSE：返回 named list
}


#' Merge several NTSDatabase reference libraries into one
#'
#' Combines a **named** list of [NTSDatabase] objects into a single database
#' whose `metadata`/`spectra_long` carry each library's rows verbatim, so the
#' per-row `db_type` survives and `batch_annotate()` can distinguish L1 / L2a /
#' L2b hits in one pass. `Lab_ID` values are prefixed with the list-name
#' (`"<name>_"`) to avoid cross-library collisions; `spectra_long$frag_mz` is
#' renamed to `mz` (the input expected by [build_nts_database()]).
#'
#' @param dbs A **named** list of [NTSDatabase] objects (e.g.
#'   `list(inhouse = db1, online = db2)`), typically from
#'   [load_builtin_pfas_library()] with `combine = FALSE`.
#' @return A single merged [NTSDatabase].
#' @seealso [load_builtin_pfas_library()]
#' @examples
#' \donttest{
#' dbs <- load_builtin_pfas_library(combine = FALSE)
#' db <- merge_nts_libraries(dbs)
#' }
#' @export
merge_nts_libraries <- function(dbs) {
  if (!is.list(dbs) || !length(dbs)) stop("`dbs` must be a non-empty list.", call. = FALSE)
  if (is.null(names(dbs)) || any(!nzchar(names(dbs)))) {
    stop("`dbs` must be a *named* list (names are used as Lab_ID prefixes).", call. = FALSE)
  }
  for (nm in names(dbs)) {
    if (!methods::is(dbs[[nm]], "NTSDatabase")) {
      stop("Element '", nm, "' is not an NTSDatabase.", call. = FALSE)
    }
  }
  md_l <- list(); sp_l <- list()
  for (nm in names(dbs)) {
    md <- data.table::as.data.table(methods::slot(dbs[[nm]], "metadata"))
    sp <- data.table::as.data.table(methods::slot(dbs[[nm]], "spectra_long"))
    if (!"mz" %in% names(sp) && "frag_mz" %in% names(sp)) {
      data.table::setnames(sp, "frag_mz", "mz")
    }
    if ("INCHIKEY" %in% names(md) && !"InChIKey" %in% names(md)) {
      data.table::setnames(md, "INCHIKEY", "InChIKey")
    }
    pre <- paste0(nm, "_")
    if ("Lab_ID" %in% names(md)) md[, Lab_ID := paste0(pre, Lab_ID)]
    if ("Lab_ID" %in% names(sp)) sp[, Lab_ID := paste0(pre, Lab_ID)]
    md_l[[nm]] <- md; sp_l[[nm]] <- sp
  }
  md_all <- data.table::rbindlist(md_l, use.names = TRUE, fill = TRUE)
  sp_all <- data.table::rbindlist(sp_l, use.names = TRUE, fill = TRUE)
  # db_type = NA：metadata 已有逐行 db_type 列时会被保留（不覆盖）
  build_nts_database(
    metadata = md_all, spectra_long = sp_all,
    noise_threshold = 0.01, top_k = 30L,
    add_merged_layer = TRUE, db_type = NA_character_,
    version_str = "1.0.0"
  )
}


# -----------------------------------------------------------------------------
# 内置包资源（谱库、表和模型）
# -----------------------------------------------------------------------------

#' Resolve a bundled file across install / source-tree layouts
#'
#' Thin (exported) wrapper around the internal lookup chain used by every
#' "defaults to a file shipped inside the package" code path. Useful when a
#' user needs the on-disk path of a bundled table, model or library — for
#' example `bundled_nts_file("extdata/tables/diagnostic_fragments.csv")`.
#'
#' @param rel Path relative to the package, e.g.
#'   `"extdata/tables/diagnostic_fragments.csv"`. A
#'   leading `"extdata/"` (or `"inst/extdata/"`) is optional.
#' @param mustWork If `TRUE` (default), error when nothing is found; otherwise
#'   return `NA_character_`.
#' @return Normalised absolute path (length 1), or `NA_character_`.
#' @seealso [nts_builtin_library_paths()], [load_builtin_pfas_library()]
#' @examples
#' bundled_nts_file("extdata/tables/diagnostic_fragments.csv")
#' @export
bundled_nts_file <- function(rel, mustWork = TRUE) {
  rel <- as.character(rel)
  if (length(rel) == 1L) {
    rel <- sub("^(inst/)?extdata/", "", rel)
  }
  .nts_extdata(rel, mustWork = mustWork)
}


# Package data paths are resolved through this shared lookup chain; raw sample
# data are supplied by the user and are not bundled with the package.
