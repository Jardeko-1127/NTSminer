#' Create an in-house database with interactive review
#'
#' @param cmpd_list Compound metadata list
#' @param mzml_file Path to mzML files
#' @param polarity Polarity (e.g. "NEG" or "POS")
#' @param mz_tol_ppm m/z tolerance in ppm
#' @param mz_tol_da m/z tolerance in Da
#' @param rt_gap_tol RT gap tolerance
#' @param interactive_review Whether to pause and ask for user confirmation
#' @param noise_threshold Noise threshold for database building
#' @param top_k Maximum number of fragments to retain
#' 
#' @return An NTSDatabase object
#' @export
create_inhouse_database <- function(cmpd_list, mzml_file, polarity = "NEG", 
                                    mz_tol_ppm = 30, mz_tol_da = 0.02, 
                                    rt_gap_tol = 9999, interactive_review = TRUE,
                                    noise_threshold = 0.01, top_k = 30L) {
  
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table is required", call. = FALSE)
  }
  
  # Step 1: Extract multi-dimensional data
  message(">> Step 1: Extracting MS2 spectra...")
  hits <- extract_mix_standards(
    cmpd_list = cmpd_list,
    mzml_file = mzml_file,
    polarity = polarity,
    mz_tol_ppm = mz_tol_ppm,
    mz_tol_da = mz_tol_da, 
    rt_gap_tol = rt_gap_tol
  )
  
  # Step 2: Interactive review checkpoint
  valid_hits <- hits[!is.na(matched_rt)]
  if (nrow(valid_hits) == 0) {
    stop("No standard spectra extracted.")
  }
  
  message(sprintf("\nTarget: Successfully extracted %d compounds.", nrow(valid_hits)))
  message("Measured retention time shifts (for reference):")
  
  name_col <- if ("name" %in% names(valid_hits)) {
    "name"
  } else if ("Name" %in% names(valid_hits)) {
    "Name"
  } else {
    "Lab_ID"
  }
  
  cols_to_show <- intersect(c(name_col, "target_rt", "matched_rt", "rt_drift"), names(valid_hits))
  print(valid_hits[, ..cols_to_show])
  
  # Step 3: Flatten ammo and force normalization
  message(">> Step 3: Flattening and normalizing MS2 spectra...")
  max_chunks <- 0L
  for (i in seq_len(nrow(valid_hits))) {
    sp <- valid_hits$ms2_spectra[[i]]
    if (length(sp)) {
      max_chunks <- max_chunks + as.integer(length(sp))
    }
  }
  if (max_chunks < 1L) {
    max_chunks <- 1L
  }
  sp_list <- vector("list", max_chunks)
  k <- 0L
  for (i in seq_len(nrow(valid_hits))) {
    lab_id <- valid_hits$Lab_ID[i]
    specs <- valid_hits$ms2_spectra[[i]]
    if (length(specs) > 0) {
      for (nce in names(specs)) {
        dt <- data.table::as.data.table(specs[[nce]])
        if (nrow(dt) > 0 && "intensity" %in% names(dt)) {
          max_int <- max(dt$intensity, na.rm = TRUE)
          if (is.finite(max_int) && max_int > 0) {
            dt[, intensity := (intensity / max_int) * 1000]
            dt <- dt[intensity >= 5]
            if (nrow(dt) > 0) {
              dt[, Lab_ID := lab_id]
              dt[, NCE := nce]
              k <- k + 1L
              sp_list[[k]] <- dt
            }
          }
        }
      }
    }
  }
  sp_list <- sp_list[seq_len(k)]
  
  if (length(sp_list) == 0) {
    stop("No valid MS2 fragments after normalization and filtering!")
  }
  
  sp_long <- data.table::rbindlist(sp_list, use.names = TRUE, fill = TRUE)
  
  if ("NCE" %in% names(sp_long)) {
    nce_vals <- unique(sp_long$NCE)
    message(sprintf("[DIAG] Detected NCE values: %s", paste(nce_vals, collapse = ", ")))
  }
  
  # Step 4: Compile S4 database
  message(">> Step 4: Compiling NTSDatabase...")
  db_inhouse <- build_nts_database(
    metadata = valid_hits,
    spectra_long = sp_long,
    noise_threshold = noise_threshold,
    top_k = top_k,
    db_type = "inhouse"
  )
  
  # Step 5: High-level metadata cleaning (Safe mode)
  message(">> Step 5: Cleaning metadata...")
  md <- methods::slot(db_inhouse, "metadata")
  
  if ("matched_rt" %in% names(md)) {
    md[, RT := matched_rt]
  }
  
  cols_to_remove <- c("target_rt", "ms2_spectra", "diagnostic_log", "matched_rt", "matched_intensity", "rt_drift")
  keep_cols <- setdiff(names(md), cols_to_remove)
  md <- md[, .SD, .SDcols = keep_cols]
  
  methods::slot(db_inhouse, "metadata") <- md
  
  message(">> Done! Database created successfully.")
  return(db_inhouse)
}
