#' MS1 annotation columns for EIC-shape proxies from [xcms::chromPeaks].
#' @keywords internal
#' @noRd
.nts_ms1_eic_shape_colnames <- function() {
  c(
    "apex_pos_ratio_median",
    "into_maxo_ratio_median",
    "peakwidth_median_s",
    "sn_median"
  )
}
