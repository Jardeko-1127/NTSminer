#' NTSminer project container (S4)
#'
#' Normalized layout: lean `feature_table`, MS1 mining in `ms1_annotations`,
#' MS2 evidence in `ms2data` and `ms2_annotations`, library hits in
#' `db_annotations`, optional `physchem_proxy`, optional `network_results`.
#' Primary MS2 similarity graphs still come from [build_molecular_network] return
#' values; `network_results` may cache tables for export only.
#'
#' @slot raw_data_path Directory containing mzML/mzXML (recursive scan supported).
#' @slot feature_table `name`, `mz`, `rt`, and per-sample intensity columns only
#'   (XCMS row count preserved; no KMD or `annotation_level` here).
#' @slot ms1_annotations Named `list` of MS1 sub-tables, one per **populated**
#'   processing step, instead of a single wide `data.frame`. Each element is a
#'   `data.frame` carrying its own `name` key (aligned with
#'   `feature_table$name`), so the pieces re-join losslessly. A step is stored
#'   only when it actually carries data (at least one non-`name` column with a
#'   non-`NA` value), so `names()` lists just the steps that ran; the possible
#'   steps are `base` (identity, bounds, `mz`, `rt`), `eic` (EIC shape metrics),
#'   `library` (`library_*` columns), `isotopes` ([flag_isotopes]),
#'   `isf` ([flag_isf] / [recognize_isf]), `kmd` ([calculate_kmd]),
#'   `adduct`, `master` (annotation-master merge), `network`
#'   (homologue / structural propagation), `annotation` (confidence and
#'   RT-prediction columns) and `extra` (any column that matches no step, so
#'   nothing is silently dropped). Use
#'   [get_ms1_annotations]`(x)` for the merged flat view,
#'   [get_ms1_annotations]`(x, "kmd")` for one step, and
#'   [set_ms1_annotations] for assignment. See `.nts_ms1_groups()` for the
#'   authoritative column-to-step map.
#' @slot xdata `XCMSnExp` result from [extract_features] (carries aligned RT).
#' @slot ms2data Named list of experimental MS2 scan lists per feature (see
#'   [export_to_gnps]).
#' @slot ms2_annotations Named `list` of MS2 sub-tables, mirroring
#'   `ms1_annotations`. Two possible elements, stored only when populated: the
#'   `annotation` table - one row per feature with MS2 evidence (`mz`, `rt`,
#'   `n_fragments` from [extract_ms2] / [purify_ms2_spectra], optional
#'   `diagnostic_fragments`, library Top-1 merge columns from
#'   [annotate_by_spectral_match], and rule-bridge columns from
#'   [validate_with_ms2]) - and the `score_grid` table (see below). Use
#'   [get_ms2_annotations]`(x)` for the merged per-feature view,
#'   [get_ms2_annotations]`(x, "score_grid")` for one sub-table, and
#'   [get_ms2_score_grid] for the grid specifically. The `score_grid` element
#'   is the full MS2 score grid from [annotate_by_spectral_match], one row per
#'   `(feature x library compound x library CE x experimental spectrum)` cell -
#'   including cells that did **not** clear `min_score`, so the gate can be
#'   re-tuned offline without re-scoring. Columns: `name`, `ref_id`, `lib_CE`,
#'   `exp_spectrum`, `score_cosine` (raw dot-product cosine, the value
#'   `min_score` gates), `n_matched_peaks`, `passes_min_score`, `score_ms2`
#'   (the multidimensional composite used to rank passing cells; `NA` when the
#'   cell failed the gate), `score_ndp`, `score_nl`, `score_md`,
#'   `score_diagnostic`, and `is_winner` (the cell behind the reported
#'   `db_annotations` row).
#' @slot clean_ms2 Named list of per-feature purified MS2 (see [purify_ms2_spectra]):
#'   each element has `feature_id`, `precursorMZ`, `rt`, `n_raw_scans`, and
#'   `spectrum` (`data.frame` with `mz`, `intensity`, `frequency`).
#' @slot db_annotations Library matching / multidimensional scoring table
#'   (possibly multiple rows per feature).
#' @slot physchem_proxy Optional `data.frame` keyed by `name` with CDK-style
#'   descriptors from [calculate_physchem_proxies] (kept **out** of
#'   `ms1_annotations` to avoid mixing qualitative tiers).
#' @slot sirius_fingerprints Structured [data.table::data.table] from
#'   [compute_sirius_fingerprints]: one row per feature per formula candidate,
#'   with `name`, `adduct`, `molecular_formula`, `rank`, `sirius_score`,
#'   `zodiac_score`, `confidence`, `explained_peaks`, a nested `fingerprint`
#'   list-column (numeric vector from CSI:FingerID), and `fp_version`.
#'   The setter validates required columns before assignment.
#' @slot network_results Optional cached `list` (e.g. `nodes` / `edges`) for
#'   export or dual-network workflows; not required for core mining steps.
#' @slot gc_spectra Named list of per-feature GC-EI full-scan spectra (see
#'   [extract_gc_spectra]). Each element is a `data.frame` with columns
#'   `mz` and `intensity`.
#' @slot gc_deconvolved Named list of deconvolution results (see
#'   [deconvolve_gc_spectra]). Names correspond to feature_table$name.
#'   Each element is a list of components, each with \code{$spectrum}
#'   (data.frame: mz, intensity) and \code{$apex_rt} (numeric).
#' @slot log List of `list(timestamp, step, params)`.
#'
#' @name NTSProject-class
#' @aliases NTSProject
#' @exportClass NTSProject
methods::setClass(
  "NTSProject",
  slots = c(
    raw_data_path = "character",
    feature_table = "data.frame",
    ms1_annotations = "list",
    xdata = "ANY",
    ms2data = "list",
    ms2_annotations = "list",
    clean_ms2 = "list",
    db_annotations = "data.frame",
    physchem_proxy = "data.frame",
    sirius_fingerprints = "data.table",
    network_results = "list",
    gc_spectra = "list",
    gc_deconvolved = "list",
    log = "list"
  ),
  prototype = list(
    # Explicit empty defaults so a bare `new("NTSProject")` is already valid:
    # without this, the `data.table`-typed `sirius_fingerprints` slot defaults
    # to a plain `list`, which R's own slot-class check rejects. Mirrors
    # [create_nts_project].
    raw_data_path = character(0),
    feature_table = data.frame(),
    ms1_annotations = list(),
    xdata = NULL,
    ms2data = list(),
    ms2_annotations = list(),
    clean_ms2 = list(),
    db_annotations = data.frame(),
    physchem_proxy = data.frame(),
    sirius_fingerprints = data.table::data.table(),
    network_results = list(),
    gc_spectra = list(),
    gc_deconvolved = list(),
    log = list()
  ),
  validity = function(object) {
    # Cheap, permissive structural checks: the slot classes already constrain
    # assignment, but an explicit validity method makes `validObject()` a
    # meaningful, uniform post-set check (mirrors NTSDatabase) and catches
    # corrupt/legacy objects loaded from disk.
    if (length(object@raw_data_path) > 1L) {
      return("slot 'raw_data_path' must be a character vector of length <= 1")
    }
    df_slots <- c(
      "feature_table", "db_annotations", "physchem_proxy"
    )
    for (nm in df_slots) {
      v <- methods::slot(object, nm)
      if (!is.null(v) && !is.data.frame(v)) {
        return(sprintf("slot '%s' must be a data.frame", nm))
      }
    }
    list_slots <- c(
      "ms1_annotations", "ms2_annotations", "ms2data", "clean_ms2",
      "network_results", "gc_spectra", "gc_deconvolved", "log"
    )
    for (nm in list_slots) {
      v <- methods::slot(object, nm)
      if (!is.null(v) && !is.list(v)) {
        return(sprintf("slot '%s' must be a list", nm))
      }
    }
    sf <- methods::slot(object, "sirius_fingerprints")
    if (!is.null(sf) && !is.data.frame(sf) && !is.list(sf)) {
      return("slot 'sirius_fingerprints' must be a data.table / data.frame")
    }
    TRUE
  }
)

#' Reference library container (S4)
#'
#' Compatible with the original `databaseClass` from `05_classdatbase02.R`.
#'
#' @slot database.info `list` with metadata (e.g. `Version`, `RT`).
#' @slot spectra.info `data.frame` with compound metadata (`Lab.ID`, `mz`,
#'   `RT`, `Compound.name`, `Formula`, `mz.pos`, `mz.neg`, ...).
#' @slot spectra.data `list` with `Spectra.positive` and/or `Spectra.negative`.
#'   Each is a named list keyed by `Lab.ID`, containing sub-lists of NCE-level
#'   `data.frame(mz, intensity)`.
#' @slot fragments.data `list` with `NEG` and/or `POS`, each a named list of
#'   diagnostic fragment tables keyed by `Lab.ID`.
#'
#' @name databaseClass-class
#' @aliases databaseClass
#' @exportClass databaseClass
methods::setClass(
  "databaseClass",
  representation(
    database.info = "list",
    spectra.info = "data.frame",
    spectra.data = "list",
    fragments.data = "list"
  ),
  prototype = list(
    database.info = list(),
    spectra.info = data.frame(matrix(nrow = 0, ncol = 0), stringsAsFactors = FALSE),
    spectra.data = list(),
    fragments.data = list()
  )
)

#' Flat reference spectral library (S4)
#'
#' Row-oriented replacement for deeply nested lists in [databaseClass]. All
#' spectral and fragment content lives in long `data.table` tables for fast
#' keyed joins and lower memory overhead than `list` trees.
#'
#' @slot metadata One row per compound / reference entry (e.g. `Lab.ID`, `mz`,
#'   `RT`, formula and names). Typically derived from legacy `spectra.info`.
#' @slot spectra_long One row per MS2 peak: keys such as `polarity`, `lab_id`,
#'   `nce`, plus `mz`, `intensity` (exact column set defined by the ETL).
#' @slot features_long One row per aggregated fragment feature: keys such as
#'   `polarity`, `lab_id`, `feature_type`, plus `mz`, `frequency` (exact
#'   columns defined by the ETL).
#' @slot version Character scalar identifying the in-memory schema (e.g.
#'   `"1.0.0"`).
#'
#' @name NTSDatabase-class
#' @aliases NTSDatabase-class NTSDatabase
#' @exportClass NTSDatabase
methods::setClass(
  "NTSDatabase",
  slots = c(
    metadata = "data.table",
    spectra_long = "data.table",
    features_long = "data.table",
    version = "character"
  ),
  prototype = list(
    metadata = data.table::data.table(),
    spectra_long = data.table::data.table(),
    features_long = data.table::data.table(),
    version = "0.1.0"
  ),
  validity = function(object) {
    if (length(object@version) != 1L) {
      return("slot 'version' must be a character vector of length 1")
    }
    if (!data.table::is.data.table(object@metadata)) {
      return("slot 'metadata' must be a data.table")
    }
    if (!data.table::is.data.table(object@spectra_long)) {
      return("slot 'spectra_long' must be a data.table")
    }
    if (!data.table::is.data.table(object@features_long)) {
      return("slot 'features_long' must be a data.table")
    }
    TRUE
  }
)


#' @keywords internal
#' @noRd
.is_reference_spectral_db <- function(x) {
  methods::is(x, "NTSDatabase")
}
