# Global variable declarations for R CMD check -------------------------------
#
# data.table and other NSE (non-standard evaluation) code refer to columns by
# bare name (e.g. `dt[, .(mz = mz_value)]`, `dt[order(rt)]`, `dt[, ..keep]`).
# R CMD check cannot see these as defined symbols and reports them as
# "no visible binding for global variable". Declaring them here suppresses
# those false positives -- it creates nothing and changes no behaviour.
#
# Generated from the R CMD check "R code for possible problems" NOTE.
utils::globalVariables(c(
  "..cols_to_show", "..comp_cols", "..keep", "..keep_cols", "..merge_cols",
  ".nts_seed_label", "ADDUCT", "Annotation", "CAS", "CAS_from_cm",
  "CAS_standalone", "CHARGE", "Charge", "Collision_energy",
  "Confidence_Score", "ExactMass", "FORMULA", "Feature_Type", "Formula",
  "Frequency", "ID", "INCHIKEY", "INSTRUMENT_TYPE", "InChIKey",
  "Instrument_type", "Kendrick_Mass", "Lab_ID", "Library_ID", "N", "NCE",
  "NIST_number", "Nominal_KM", "Polarity", "Precursor_type",
  "Propagation_Level", "PubChem_CID", "PubChem_line", "RT", "SMILES", "a",
  "add_formula", "add_formula_r", "alt_cmp", "alt_lib", "anchor_feature",
  "anchor_formula", "anchor_mz", "anchor_name", "anchor_seed_formula",
  "ann_name", "apex_rt", "b", "base_for_step", "base_label", "basepeak_int",
  "basepeak_mz", "block_id", "bp_bin", "bridge_col", "calc_mz",
  "canopus_class", "canopus_npc", "canopus_subclass", "canopus_superclass",
  "cas", "child_row_id", "cint", "cluster", "cluster_id", "cmz", "color",
  "compound", "compound_id", "compound_key", "confidence_level",
  "correct_Polarity", "cosine", "crt", "cum_cf2", "d", "delta", "depth",
  "diagnostic_log", "edge_category", "edge_score", "err", "evidence_source",
  "exact_mass", "exp_CE", "exp_mz", "exp_rt", "f_rt_hi", "f_rt_lo",
  "feat_key", "feat_label", "feat_ord", "feature_id", "feature_name",
  "filename", "fill_plot", "final_Pol", "fingerprint", "fo", "fp_version",
  "frag_mz", "frag_mz_1", "frag_mz_2", "fragment_formula", "fragment_mass",
  "freq_lab", "frequency", "from_mz", "ft", "ft_mz", "group", "has_db",
  "has_form", "has_inf", "has_struct", "hi", "hl_element", "hl_n",
  "homologue_label", "homologue_steps", "i.annotation_level",
  "i.diagnostic_fragments", "i.homologous_series_id", "i.max_prevalence",
  "i.mz", "i.n_fragment_rows", "i.n_lab_class", "i.n_lab_with_features",
  "i.transformation_parent", "i.transformation_rule",
  "i.transformed_formula", "id_chr", "imax", "imz_hi", "imz_lo", "int",
  "int_max", "int_rel", "ion_mode", "ion_mz", "ion_rank", "is_anchor",
  "is_winner", "isf_candidate_ids", "isf_candidate_n_fragments",
  "isf_candidate_names", "isolated", "isotope_cosine", "j", "join_k", "k",
  "keep", "kid", "kid_raw", "kmd_gap", "lab0", "lab_id", "lib_CE", "lib_mz",
  "lib_rt", "library_Compound_name", "line", "lo", "loss", "loss_a",
  "loss_b", "loss_error", "mass_diff", "mass_difference", "mass_shift",
  "matched_diagnostic", "matched_intensity", "matched_rt", "max_int",
  "max_prevalence", "max_rk", "md", "md_a", "md_b", "md_hi", "md_lo", "mx",
  "mz.ft", "mz_a", "mz_b", "mz_err", "mz_error_ppm", "mz_exp", "mz_ft",
  "mz_group", "mz_hi", "mz_k", "mz_key", "mz_label", "mz_lib", "mz_lo",
  "mz_lower", "mz_max", "mz_min", "mz_neg", "mz_obs", "mz_plot", "mz_pos",
  "mz_rep", "mz_tol_abs", "mz_upper", "mz_value", "n_fragment_rows",
  "n_fragments", "n_grp", "n_hit", "n_hits", "n_ions", "n_lab_class",
  "n_lab_metadata", "n_lab_with_features", "n_matched_peaks",
  "n_mz_in_base", "n_present", "n_raw_scans", "n_scans", "name_a", "name_b",
  "nap_structure_depth", "nap_structure_formula", "nap_structure_seed",
  "nap_structure_smiles", "nap_structure_tanimoto", "nce_level",
  "network_propagation", "neutral_loss", "neutral_mass_raw", "new_Lab_ID",
  "nf", "node_id", "npc", "occurrence", "ok_alt_cmp", "ok_alt_lib",
  "ok_exp", "ok_lib", "ok_sx", "only_struct", "par_hit", "parent_key",
  "parent_row_id", "pc", "pcid", "peak_data", "peak_ratio", "pint", "pl_i",
  "pmz", "ppm", "prec_mass", "prec_mz", "prec_mz_work", "precursorMZ",
  "pred_mz", "prevalence", "priority", "propagated_ID", "propagated_SMILES",
  "propagated_from_node", "propagated_mcs_score", "propagated_similarity",
  "prt", "pw", "r_rt_hi", "r_rt_lo", "raw_freq", "ref_id", "rel_int",
  "relative", "resolution_label", "retentionTime", "rt_bin", "rt_drift",
  "rt_error_s", "rt_ft", "rt_hi", "rt_lo", "rt_max", "rt_min", "rt_obs",
  "rule_hit", "rule_label", "rule_mz", "rule_type", "sc", "scan_id",
  "scan_order", "score", "score_diagnostic", "seed_anchor", "seed_f",
  "seed_formula", "seed_mz", "seed_parent", "seqNum", "series_steps",
  "shared_fragments_count", "shared_fragments_mz", "shared_mds_count",
  "shared_mds_mz", "shared_nls_count", "shared_nls_mz", "smi", "step_denom",
  "sub_formula", "sub_formula_r", "sub_group", "subclass", "sx",
  "target_rt", "tf_hit", "tic", "tic_apex", "title", "to_mz", "tol",
  "tol_ppm", "transformation_parent", "transformation_rule",
  "transformed_formula", "transformed_formula_edge", "type", "x.mz",
  # --- network layer (0.1.25): derived edge columns --------------------------
  "n_steps", "n_a", "n_b", "mz_from", "transformation", "ms2_ladder_confirmed",
  # --- 0.1.30 R CMD check NOTE: remaining data.table NSE columns -------------
  "NCE_Count", "best_ms2_index", "db_type", "gid", "loss_b_val", "md_b_val",
  "mz_b_val", "n_ms2_subset", "purity", "ref_mz",
  # --- 0.1.37 network layer: hybrid semantic network NSE columns -------------
  "semantic_similarity", "sim", "u", "v",
  # --- 0.1.38 scaffold inference layer: semantic analog & scaffold retrieval --
  "anchor_id", "anchor_smiles", "similarity", "analog_count", "generic_scaffold",
  "scaffold_smiles", "is_cyclic", "num_rings",
  # --- 0.1.39 forward spectrum prediction layer: FraG-FIORA in-silico MS/MS ---
  "candidate_id", "score_pred", "cosine_sqrt", "cosine_linear", "explained_intensity",
  "matched_peaks", "total_pred_peaks"
))


# Shared feature_table metadata columns ----------------------------------------
#
# Hand-written (not generated by R CMD check). These are the position/shape
# columns of an `@feature_table` (see the column set built by
# 03_feature_extract.R): xcms medians/min/max in Da for m/z and *seconds* for
# RT, plus the peak count. They are numeric, so any "take the numeric columns"
# heuristic mistakes them for sample abundances unless excluded explicitly.
# Use this constant wherever a metadata blocklist is needed so the five
# pickers in the package cannot drift apart again.
.nts_feature_shape_cols <- c(
  "mz", "mzmed", "mzmin", "mzmax",
  "rt", "rtmed", "rtmin", "rtmax",
  "npeaks"
)
