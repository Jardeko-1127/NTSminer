# ---------------------------------------------------------------------------
# Annotation merge / priority ranking (0.1.32)
#   rank_annotation_candidates() : pure, in-memory core
#   merge_annotation_sources()   : read-only wrapper over a pipeline run dir
# ---------------------------------------------------------------------------

# Small helper: a candidate long table with only the columns we care about.
cand_df <- function(feature, ref_id, candidate = ref_id, db_type = NA_character_,
                    score_total = NA_real_, score_ms2 = NA_real_,
                    ev_tier = NA_real_, network_support = NA_character_,
                    SMILES = NA_character_) {
  data.frame(
    feature = feature, ref_id = ref_id, candidate = candidate,
    db_type = db_type, score_total = score_total, score_ms2 = score_ms2,
    ev_tier = ev_tier, network_support = network_support, SMILES = SMILES,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# 1. base level = confidence of the identification source
# ---------------------------------------------------------------------------

test_that("base level follows the identification-source library", {
  d <- rbind(
    cand_df("F1", "inhouse_L01", "PFOA", score_total = 0.90),
    cand_df("F2", "online_L02", "PFOS", score_total = 0.80),
    cand_df("F3", "insilico_L03", "PFHxA", score_total = 0.70)
  )
  out <- rank_annotation_candidates(d, verbose = FALSE)
  f <- out$final
  expect_equal(f$base_level[match(c("F1", "F2", "F3"), f$name)],
               c("L1", "L2a", "L2b"))
  expect_equal(f$decision_label, c("L1", "L2a", "L2b"))
  expect_equal(f$n_families, c(1L, 1L, 1L))
  expect_true(all(f$conflict_status == "none"))
})

test_that("db_type fills in for ref_ids without a known library prefix", {
  d <- cand_df("F1", "L001", "X", db_type = "online", score_total = 0.5)
  out <- rank_annotation_candidates(d, verbose = FALSE)
  expect_equal(out$final$base_level, "L2a")
  expect_true(is.na(out$candidates$comp_db))   # own library unknown, only display
  expect_equal(out$candidates$ident_set, "online")
})

# ---------------------------------------------------------------------------
# 2. corroboration is a bounded half-rung, never additive voting
# ---------------------------------------------------------------------------

test_that("two independent families add a bounded '*' rung", {
  d <- rbind(
    cand_df("F1", "inhouse_L01", "PFOA", score_total = 0.90),
    cand_df("F1", "insilico_L07", "PFOA-pred", score_total = 0.50)
  )
  mf <- data.frame(name = "F1", identifier = "inhouse_L01",
                   Compound.name = "PFOA", score_total = 120, rank = 1L,
                   is_top1 = TRUE, stringsAsFactors = FALSE)
  net <- data.frame(feature = "F1", ref_id = "inhouse_L01", net_kind = "DFI",
                    stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, metfrag_hits = mf,
                                    network_layers = net, verbose = FALSE)
  f <- out$final
  expect_equal(f$ref_id, "inhouse_L01")
  expect_equal(f$decision_label, "L1*")
  expect_equal(f$n_families, 3L)       # F1 (library) + F2 (MetFrag) + F4 (network)
  expect_equal(f$decision_tier, 1L)    # the '*' never promotes past the level
})

test_that("ev_tier == 0 is not independent support", {
  d <- rbind(
    cand_df("F1", "inhouse_L01", "PFOA", score_total = 0.90, ev_tier = 0),
    cand_df("F1", "online_L02", "PFOS", score_total = 0.70, ev_tier = 2)
  )
  out <- rank_annotation_candidates(d, verbose = FALSE)
  cr <- out$candidates
  r0 <- cr[cr$ref_id == "inhouse_L01", ]
  r2 <- cr[cr$ref_id == "online_L02", ]
  expect_false(r0$family_F4)
  expect_equal(r0$n_families, 1L)
  expect_equal(r0$decision_label, "L1")   # no '*': "it is itself an L1 hit" is not support
  expect_true(r2$family_F4)
  expect_equal(r2$n_families, 2L)
  expect_equal(r2$decision_label, "L2a*")
})

# ---------------------------------------------------------------------------
# 3. MetFrag: its own L2c rung, and promotion only when it really scores
# ---------------------------------------------------------------------------

test_that("MetFrag top-1 promotes a predicted-spectrum hit, and only when it scores", {
  d <- cand_df("F1", "insilico_L01", "pred", score_total = 0.5)
  mf_ok <- data.frame(name = "F1", identifier = "insilico_L01",
                      Compound.name = "pred", score_total = 100, rank = 1L,
                      is_top1 = TRUE, stringsAsFactors = FALSE)
  a <- rank_annotation_candidates(d, mf_ok, verbose = FALSE)
  b <- rank_annotation_candidates(d, transform(mf_ok, score_total = 0),
                                  verbose = FALSE)
  cc <- rank_annotation_candidates(d, transform(mf_ok, is_top1 = FALSE),
                                   verbose = FALSE)
  expect_equal(a$final$decision_tier, 3L)   # valid top-1 lifts L2b
  expect_equal(b$final$decision_tier, 4L)   # score 0 -> no discrimination, no lift
  expect_equal(cc$final$decision_tier, 4L)  # not top-1 -> no lift
  expect_equal(b$final$base_level, "L2b")
})

test_that("network-only candidates land on L3b / L3c by relation kind", {
  d <- cand_df("F1", "insilico_L01", "seed", score_total = 0.9)
  net <- data.frame(feature = c("F1", "F1"),
                    ref_id = c("insilico_L99", "insilico_L98"),
                    net_kind = c("DFI", "HOMO"), stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, network_layers = net, verbose = FALSE)
  cr <- out$candidates
  expect_equal(cr$base_level[cr$ref_id == "insilico_L99"], "L3b")
  expect_equal(cr$base_level[cr$ref_id == "insilico_L98"], "L3c")
  expect_equal(cr$decision_tier[cr$ref_id == "insilico_L99"], 8L)
  expect_equal(cr$decision_tier[cr$ref_id == "insilico_L98"], 9L)
  expect_true(all(cr$src_net[cr$ref_id %in% c("insilico_L99", "insilico_L98")]))
  expect_true(all(cr$n_families[cr$ref_id %in% c("insilico_L99", "insilico_L98")] == 1L))
})

test_that("a MetFrag-only candidate is L2c and steps up with relationship support", {
  d <- cand_df("F1", "insilico_L01", "seed", score_total = 0.9)
  mf <- data.frame(name = c("F1", "F1"),
                   identifier = c("insilico_L97", "insilico_L96"),
                   Compound.name = c("mf_rel", "mf_plain"),
                   score_total = c(44, 30), rank = c(1L, 1L),
                   is_top1 = c(TRUE, FALSE), stringsAsFactors = FALSE)
  net <- data.frame(feature = "F1", ref_id = "insilico_L97", net_kind = "MD",
                    stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, metfrag_hits = mf,
                                    network_layers = net, verbose = FALSE)
  cr <- out$candidates
  expect_equal(cr$base_level[cr$ref_id == "insilico_L97"], "L2c")
  expect_equal(cr$decision_tier[cr$ref_id == "insilico_L97"], 5L)
  expect_equal(cr$base_level[cr$ref_id == "insilico_L96"], "L2c")
  expect_equal(cr$decision_tier[cr$ref_id == "insilico_L96"], 6L)
})

# ---------------------------------------------------------------------------
# 4. ordering inside a feature
# ---------------------------------------------------------------------------

test_that("candidates are ordered by tier, not by library score", {
  d <- rbind(
    cand_df("F1", "insilico_L01", "pred-b", score_total = 0.95),
    cand_df("F1", "online_L02", "online-a", score_total = 0.40),
    cand_df("F1", "inhouse_L03", "inh", score_total = 0.10)
  )
  out <- rank_annotation_candidates(d, verbose = FALSE)
  expect_equal(out$final$ref_id, "inhouse_L03")   # L1 wins despite the lowest score
  cr <- out$candidates
  expect_equal(cr$ref_id[order(cr$rank_in_feature)],
               c("inhouse_L03", "online_L02", "insilico_L01"))
  expect_equal(cr$rank_in_feature, 1:3)
  expect_equal(cr$n_candidates, c(3L, 3L, 3L))
  expect_true(cr$is_top[1])
})

# ---------------------------------------------------------------------------
# 5. the four review criteria + conflict_status
# ---------------------------------------------------------------------------

test_that("tied candidates on the same tier are flagged for review", {
  d <- rbind(
    cand_df("F1", "online_L01", "PFOS-A", score_total = 0.70),
    cand_df("F1", "online_L02", "PFOS-B", score_total = 0.69)
  )
  out <- rank_annotation_candidates(d, verbose = FALSE)
  f <- out$final
  expect_equal(f$decision_label, "L2a")
  expect_match(f$tied_alternatives, "PFOS-")
  expect_match(f$review_reason, "tied_alternatives")
  expect_equal(f$conflict_status, "resolved_by_lib")  # a measured library decides
  expect_equal(nrow(out$conflicts), 1L)
  expect_equal(out$conflicts$name, "F1")
})

test_that("MetFrag disagreeing with a measured library hit is a resolvable conflict", {
  d <- cand_df("F1", "inhouse_L01", "PFOA", score_total = 0.9)
  mf <- data.frame(name = "F1", identifier = "insilico_L02",
                   Compound.name = "PFNA", score_total = 99, rank = 1L,
                   is_top1 = TRUE, stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, metfrag_hits = mf, verbose = FALSE)
  f <- out$final
  expect_equal(f$ref_id, "inhouse_L01")
  expect_equal(f$mf_vs_lib, "differs")
  expect_equal(f$lib_best_scope, "measured")
  expect_match(f$review_reason, "mf_vs_lib_differs")
  expect_equal(f$conflict_status, "resolved_by_lib")
})

test_that("a propagated identity that loses to a non-measured best is an open conflict", {
  d <- cand_df("F1", "insilico_L01", "pred", score_total = 0.9)
  net <- data.frame(feature = "F1", ref_id = "insilico_L77", net_kind = "DFI",
                    stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, network_layers = net, verbose = FALSE)
  f <- out$final
  expect_equal(f$ref_id, "insilico_L01")
  expect_equal(f$net_ident_ref, "insilico_L77")
  expect_match(f$review_reason, "net_ident_differs")
  expect_equal(f$conflict_status, "open")
})

test_that("MetFrag vs an insilico-only library is a prediction-level conflict", {
  d <- cand_df("F1", "insilico_L01", "pred-a", score_total = 0.9)
  mf <- data.frame(name = "F1", identifier = "insilico_L02",
                   Compound.name = "pred-b", score_total = 99, rank = 1L,
                   is_top1 = TRUE, stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, metfrag_hits = mf, verbose = FALSE)
  f <- out$final
  expect_equal(f$lib_best_scope, "insilico")
  expect_equal(f$mf_vs_lib, "differs")
  # predicted-vs-predicted disagreement is the softer criterion: "insilico"
  # (not "measured") sends it to mf_vs_insilico_differs, which stays open
  expect_match(f$review_reason, "mf_vs_insilico_differs")
  expect_false(grepl("mf_vs_lib_differs", f$review_reason))
  expect_equal(f$conflict_status, "open")
})

# ---------------------------------------------------------------------------
# 6. naming / feature-key handling
# ---------------------------------------------------------------------------

test_that("blank candidate names are filled from the suspect table, else the ref_id", {
  d <- cand_df("F1", "insilico_L01", candidate = "", score_total = 0.5)
  sus <- data.frame(ID = "insilico_L01", Name = "PFHxA (pred)",
                    InChIKey = "ABC", stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, suspect = sus, verbose = FALSE)
  expect_equal(out$final$candidate, "PFHxA (pred)")
  expect_equal(out$final$inchikey, "ABC")

  out2 <- rank_annotation_candidates(d, verbose = FALSE)
  expect_equal(out2$final$candidate, "insilico_L01")
  expect_true(is.na(out2$final$inchikey))
})

test_that("Compound.name is used when the candidate column is blank", {
  d <- data.frame(feature = "F1", ref_id = "inhouse_L01", candidate = "",
                  Compound.name = "PFOA", score_total = 0.9,
                  stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, verbose = FALSE)
  expect_equal(out$final$candidate, "PFOA")
})

test_that("`name` and `feature_id` are accepted as the feature column", {
  d1 <- data.frame(name = "F1", ref_id = "inhouse_L01", candidate = "PFOA",
                   score_total = 0.9, stringsAsFactors = FALSE)
  out1 <- rank_annotation_candidates(d1, verbose = FALSE)
  expect_equal(out1$final$name, "F1")
  expect_equal(out1$final$base_level, "L1")

  d2 <- data.frame(feature_id = "F2", ref_id = "online_L01", candidate = "PFOS",
                   score_total = 0.8, stringsAsFactors = FALSE)
  out2 <- rank_annotation_candidates(d2, verbose = FALSE)
  expect_equal(out2$final$name, "F2")
  expect_equal(out2$final$base_level, "L2a")
})

test_that("feature_meta supplies coordinates, report confidence and structure flags", {
  d <- cand_df("F1", "inhouse_L01", "PFOA", score_total = 0.9)
  fm <- data.frame(name = "F1", mz = 413.0, rt = 300.5,
                   confidence_level = "L1", pfas_ladder_len = 4,
                   pfas_subclass_hint = "PFCAs", stringsAsFactors = FALSE)
  out <- rank_annotation_candidates(d, feature_meta = fm, verbose = FALSE)
  f <- out$final
  expect_equal(f$mz, 413.0)
  expect_equal(f$rt, 300.5)
  expect_equal(f$report_confidence, "L1")
  expect_true(f$flag_ladder)
  expect_equal(f$flag_subclass, "PFCAs")
})

# ---------------------------------------------------------------------------
# 7. output shape, summary, idempotence, input validation
# ---------------------------------------------------------------------------

test_that("summary counts features per final label and the run is idempotent", {
  d <- rbind(
    cand_df("F1", "inhouse_L01", "A", score_total = 0.9),
    cand_df("F2", "online_L01", "B", score_total = 0.8)
  )
  a <- rank_annotation_candidates(d, verbose = FALSE)
  b <- rank_annotation_candidates(d, verbose = FALSE)
  expect_identical(a, b)
  expect_named(a, c("candidates", "final", "conflicts", "summary"),
               ignore.order = TRUE)
  expect_equal(a$summary$final_label, c("L1", "L2a"))
  expect_equal(sum(a$summary$n_features), 2L)
  expect_equal(nrow(a$final), 2L)
  expect_true(all(c("name", "decision_label", "review_reason") %in% names(a$final)))
})

test_that("bad input is rejected with a clear error", {
  expect_error(rank_annotation_candidates(NULL, verbose = FALSE), "data\\.frame")
  expect_error(rank_annotation_candidates(data.frame(), verbose = FALSE), "no rows")
  expect_error(
    rank_annotation_candidates(
      data.frame(feature = "F1", score_total = 1), verbose = FALSE
    ),
    "ref_id"
  )
  expect_error(
    rank_annotation_candidates(data.frame(x = 1), verbose = FALSE),
    "feature column"
  )
  expect_error(
    rank_annotation_candidates(cand_df("", "inhouse_L01", "A"), verbose = FALSE),
    "non-empty feature id"
  )
})

# ---------------------------------------------------------------------------
# 8. the run_dir wrapper
# ---------------------------------------------------------------------------

write_run_dir <- function(root) {
  dir.create(file.path(root, "network"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "metfrag"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "multilib", "combined"), recursive = TRUE,
             showWarnings = FALSE)
  utils::write.csv(
    data.frame(feature = c("F1", "F2"), ref_id = c("inhouse_L01", "online_L02"),
               candidate = c("PFOA", "PFOS"), db_type = c("inhouse", "online"),
               score_total = c(0.9, 0.8), score_ms2 = c(0.95, 0.88),
               stringsAsFactors = FALSE),
    file.path(root, "network", "candidate_ranking.csv"), row.names = FALSE
  )
  utils::write.csv(
    data.frame(name = "F1", identifier = "insilico_L09", Compound.name = "PFHxA",
               score_total = 55, rank = 1L, is_top1 = TRUE,
               stringsAsFactors = FALSE),
    file.path(root, "metfrag", "metfrag_candidates.csv"), row.names = FALSE
  )
  utils::write.csv(
    data.frame(name = c("F1", "F2"), mz = c(413.0, 499.9), rt = c(300.5, 420.0),
               confidence_level = c("L1", "L2a"), stringsAsFactors = FALSE),
    file.path(root, "annotation_master.csv"), row.names = FALSE
  )
  utils::write.csv(
    data.frame(name = c("F1", "F2"), net_df_ref_id = c("insilico_L88", NA),
               stringsAsFactors = FALSE),
    file.path(root, "multilib", "combined", "ms2_annotations_with_qc.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(name = "F2", homologue_ref_id = "insilico_L87",
               stringsAsFactors = FALSE),
    file.path(root, "multilib", "combined", "ms1_annotations_with_confidence.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(ID = c("inhouse_L01", "online_L02", "insilico_L88", "insilico_L87"),
               Name = c("PFOA", "PFOS", "NFDHA", "GenX"),
               InChIKey = c("K1", "K2", "K3", "K4"), stringsAsFactors = FALSE),
    file.path(root, "NTS_PFAS_suspect_3lib_dedup.csv"), row.names = FALSE
  )
  invisible(root)
}

test_that("merge_annotation_sources reads a run dir and writes the four tables", {
  root <- tempfile("nts_rank_run")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  write_run_dir(root)

  res <- merge_annotation_sources(root, verbose = FALSE)
  merged <- file.path(root, "merged")
  expect_true(all(file.exists(file.path(merged, c(
    "candidates_ranked.csv", "final_annotation.csv",
    "merge_conflicts.csv", "merge_summary.csv"
  )))))

  f <- res$final
  expect_equal(sort(f$name), c("F1", "F2"))
  expect_equal(f$candidate[f$name == "F1"], "PFOA")
  expect_equal(f$base_level[f$name == "F1"], "L1")
  expect_equal(f$base_level[f$name == "F2"], "L2a")
  expect_equal(f$mz[f$name == "F1"], 413.0)

  # network identities come from the `combined` multilib set: MS2 edge + homologue
  expect_equal(f$net_ident_ref[f$name == "F1"], "insilico_L88")
  expect_equal(f$net_ident_ref[f$name == "F2"], "insilico_L87")
  # ...and the suspect table names them (network identities have no library hit)
  expect_equal(f$inchikey[f$name == "F1"], "K1")
  cr <- res$candidates
  expect_equal(cr$candidate[cr$ref_id == "insilico_L88"], "NFDHA")
  expect_equal(cr$base_level[cr$ref_id == "insilico_L87"], "L3c")

  back <- utils::read.csv(file.path(merged, "final_annotation.csv"),
                          stringsAsFactors = FALSE)
  expect_equal(nrow(back), nrow(f))
  expect_equal(back$name, f$name)
  expect_equal(back$decision_label, f$decision_label)
})

test_that("write_files = FALSE leaves the run directory byte-identical", {
  root <- tempfile("nts_rank_ro")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  write_run_dir(root)
  before <- tools::md5sum(list.files(root, recursive = TRUE, full.names = TRUE))

  res <- merge_annotation_sources(root, out_dir = file.path(root, "out2"),
                                  write_files = FALSE, verbose = FALSE)
  expect_false(dir.exists(file.path(root, "out2")))
  after <- tools::md5sum(list.files(root, recursive = TRUE, full.names = TRUE))
  expect_identical(unname(before), unname(after))
  expect_s3_class(res$final, "data.frame")
})

test_that("the wrapper falls back to candidates_all_sources.csv (library rows only)", {
  root <- tempfile("nts_rank_fb")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(root, recursive = TRUE)
  utils::write.csv(
    data.frame(name = c("F1", "F1"), source = c("library", "metfrag"),
               ref_id = c("inhouse_L01", NA), candidate = c("PFOA", "PFHxA"),
               db_type = c("inhouse", NA), score_total = c(0.9, 50),
               stringsAsFactors = FALSE),
    file.path(root, "candidates_all_sources.csv"), row.names = FALSE
  )
  res <- merge_annotation_sources(root, write_files = FALSE, verbose = FALSE)
  expect_equal(res$final$name, "F1")
  expect_equal(res$final$candidate, "PFOA")
})

test_that("the wrapper rejects a missing run dir or an empty one", {
  expect_error(merge_annotation_sources(file.path(tempfile("nope"), "x"),
                                        verbose = FALSE),
               "existing directory")
  root <- tempfile("nts_rank_empty")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(root, recursive = TRUE)
  expect_error(merge_annotation_sources(root, verbose = FALSE),
               "No usable candidate")
})
