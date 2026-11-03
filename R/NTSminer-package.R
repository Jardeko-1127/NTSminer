#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @import methods
#' @import ggplot2
#' @import igraph
#' @importFrom BiocParallel bpparam register SnowParam MulticoreParam bplapply bpnworkers
#' @importFrom MsCoreUtils ndotproduct
#' @importFrom MSnbase readMSData precursorMz rtime mz intensity spectra fileNames fData
#' @importFrom xcms CentWaveParam findChromPeaks adjustRtime ObiwarpParam groupChromPeaks PeakDensityParam fillChromPeaks featureDefinitions featureValues chromatogram
#' @importFrom stats cor dist setNames predict median quantile sd
#' @importFrom utils write.csv txtProgressBar setTxtProgressBar head tail write.table getFromNamespace
#' @importFrom grDevices pdf dev.off
#' @importFrom visNetwork visEdges visEvents visLegend visNetwork visNodes visOptions visPhysics
#' @importFrom rcdk parse.smiles get.exact.mass get.fingerprint eval.desc
#' @importFrom fingerprint distance
NULL
