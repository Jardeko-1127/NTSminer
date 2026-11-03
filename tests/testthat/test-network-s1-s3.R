# Tests for the molecular-network S1-S3 hardening (feature network, path A).
# Design: docs/analysis/NTSminer_分子网络设计方案_v0.1.24.md
# Scope: S1 (vocabulary unify + edge_score normalisation), S2 (scale guards +
# homologue step/label + md_min_da), S3 (cf2_ladder_boost). All additive: no
# project (@prj) slot is written by any of these helpers.

CF2 <- 49.9968064

# ---------------------------------------------------------------------------
# S1-b: normalized edge_score
# ---------------------------------------------------------------------------

test_that(".nts_ms2_edge_score is n_hit / min(n_a, n_b), clamped to (0, 1]", {
  # 2 shared of a 2-peak spectrum vs a 40-peak spectrum -> saturated by the
  # smaller spectrum, so the reproducible pair scores 1.
  expect_equal(.nts_ms2_edge_score(2, 2, 40), 1)
  expect_equal(.nts_ms2_edge_score(2, 40, 2), 1)
  # 2 shared of two 40-peak spectra -> weak
  expect_equal(.nts_ms2_edge_score(2, 40, 40), 0.05)
  # never exceeds 1 even if n_hit > min(n_a, n_b) (guards against count drift)
  expect_equal(.nts_ms2_edge_score(5, 3, 3), 1)
  # vectorised
  expect_equal(.nts_ms2_edge_score(c(1, 2), c(4, 4), c(4, 4)), c(0.25, 0.5))
  # non-finite numerators -> NA; degenerate denominators clamp to 1 (never NaN/Inf)
  expect_true(is.na(.nts_ms2_edge_score(NA, 2, 3)))
  expect_equal(.nts_ms2_edge_score(1, 0, 0), 1)
  expect_false(any(is.nan(.nts_ms2_edge_score(c(1, 2, 3), c(0, 2, NA), c(2, 2, 2)))))
})

# ---------------------------------------------------------------------------
# S1-a / S1-b: vocabulary aliases + category split
# ---------------------------------------------------------------------------

test_that(".nts_normalize_network_edge_layers maps aliases case-insensitively", {
  expect_equal(
    .nts_normalize_network_edge_layers(c("DFI", "NL", "MD")),
    c("fragment", "neutral_loss", "mass_difference")
  )
  expect_equal(
    .nts_normalize_network_edge_layers(c("fragment", "Neutral Loss", "mass_diff")),
    c("fragment", "neutral_loss", "mass_difference")
  )
  # dedup after normalisation
  expect_equal(
    .nts_normalize_network_edge_layers(c("nl", "neutral_loss")),
    "neutral_loss"
  )
})

test_that(".nts_normalize_ms1_edge_types maps aliases and drops unknowns", {
  expect_equal(
    .nts_normalize_ms1_edge_types(c("KMD", "hom", "series", "trans", "transform", "isf")),
    c("homologue", "transformation", "isf")
  )
  expect_equal(.nts_normalize_ms1_edge_types("bogus"), character(0))
  expect_equal(.nts_normalize_ms1_edge_types(NULL), character(0))
})

test_that(".nts_network_assign_edge_category splits Transformation from Homologue", {
  E <- data.table::data.table(
    from = c("A", "A", "A"),
    to = c("B", "B", "B"),
    edge_type = c("Homologue", "Transformation", "ISF")
  )
  out <- .nts_network_assign_edge_category(E)
  expect_equal(out$edge_category, c("Homologue", "Transformation", "ISF"))
})

test_that(".nts_network_edges_filter_by_type keeps slices disjoint", {
  E <- data.table::data.table(
    from = c("A", "A", "A", "A"),
    to = c("B", "B", "B", "B"),
    edge_category = c("Homologue", "Transformation", "ISF", "Fragment"),
    Interaction = c(NA, NA, NA, "MS2_Shared_Fragment")
  )
  expect_equal(nrow(.nts_network_edges_filter_by_type(E, "homologue")), 1L)
  expect_equal(nrow(.nts_network_edges_filter_by_type(E, "transformation")), 1L)
  expect_equal(nrow(.nts_network_edges_filter_by_type(E, "isf")), 1L)
  expect_equal(nrow(.nts_network_edges_filter_by_type(E, "fragment")), 1L)
  expect_equal(nrow(.nts_network_edges_filter_by_type(E, "all")), 4L)
})

# ---------------------------------------------------------------------------
# S2-a: scale guards (top_k / max_edges)
# ---------------------------------------------------------------------------

test_that(".nts_network_apply_edge_limits is a no-op by default", {
  E <- data.table::data.table(
    from = c("A", "A", "B", "C", "D"),
    to = c("B", "C", "C", "D", "E"),
    edge_score = c(0.9, 0.5, 0.7, 0.2, 0.1)
  )
  out <- .nts_network_apply_edge_limits(E)
  expect_equal(nrow(out), 5L)
  expect_equal(out$from, E$from)
  expect_equal(out$edge_score, E$edge_score)
})

test_that(".nts_network_apply_edge_limits keeps top_k strongest edges per node", {
  E <- data.table::data.table(
    from = c("A", "A", "B", "C", "D"),
    to = c("B", "C", "C", "D", "E"),
    edge_score = c(0.9, 0.5, 0.7, 0.2, 0.1)
  )
  out <- .nts_network_apply_edge_limits(E, top_k = 1)
  # A keeps A-B (0.9); C keeps B-C (0.7); D keeps C-D (0.2); D-E survives because
  # E's only edge is D-E. A-C (0.5) is bottom-ranked for both its endpoints.
  expect_setequal(paste(out$from, out$to), c("A B", "B C", "C D", "D E"))
})

test_that(".nts_network_apply_edge_limits caps total at max_edges by score", {
  E <- data.table::data.table(
    from = c("A", "A", "B", "C", "D"),
    to = c("B", "C", "C", "D", "E"),
    edge_score = c(0.9, 0.5, 0.7, 0.2, 0.1)
  )
  out <- .nts_network_apply_edge_limits(E, max_edges = 2)
  expect_equal(nrow(out), 2L)
  expect_setequal(paste(out$from, out$to), c("A B", "B C"))
  # finite edge_score sorts ahead of NA / non-finite
  E2 <- data.table::data.table(
    from = c("A", "B"), to = c("B", "C"),
    edge_score = c(NA_real_, 0.3)
  )
  out2 <- .nts_network_apply_edge_limits(E2, max_edges = 1)
  expect_equal(paste(out2$from, out2$to), "B C")
})

# ---------------------------------------------------------------------------
# S2-b: homologue step count + signed label
# ---------------------------------------------------------------------------

test_that(".nts_ms1_layer_edges_dt adds CF2 step count + signed label", {
  nodes <- data.table::data.table(
    name = c("A", "B", "C"),
    mz = c(100, 100 + CF2, 100 + 2 * CF2),
    kmd_series_id = c("S1", "S1", "S1"),
    pfas_ladder_len = c(4, 4, 4)
  )
  e <- .nts_ms1_layer_edges_dt(
    nodes, valid_ids = nodes$name,
    ms1_layers = "homologue", cf2_ladder_boost = FALSE
  )
  expect_true(all(c("mass_difference", "n_steps", "transformation") %in% names(e)))
  expect_equal(e$n_steps, c(1L, 1L))
  expect_equal(e$transformation, c("+1CF2", "+1CF2"))
  expect_equal(e$label, c("+1CF2", "+1CF2"))
  expect_equal(e$mass_difference, c(CF2, CF2), tolerance = 1e-6)
})

test_that(".nts_ms1_layer_edges_dt handles a 2-step gap within max_gap", {
  nodes <- data.table::data.table(
    name = c("A", "C"),
    mz = c(100, 100 + 2 * CF2),
    kmd_series_id = c("S1", "S1")
  )
  e <- .nts_ms1_layer_edges_dt(
    nodes, valid_ids = nodes$name,
    ms1_layers = "homologue", cf2_ladder_boost = FALSE
  )
  expect_equal(e$n_steps, 2L)
  expect_equal(e$transformation, "+2CF2")
})

# ---------------------------------------------------------------------------
# S3-a: cf2_ladder_boost reinforces (never invents) MS1 homologue edges
# ---------------------------------------------------------------------------

test_that("cf2_ladder_boost marks ladder-confirmed edges and scores them 1", {
  nodes <- data.table::data.table(
    name = c("A", "B", "C"),
    mz = c(100, 100 + CF2, 100 + 2 * CF2),
    kmd_series_id = c("S1", "S1", "S1"),
    pfas_ladder_len = c(4, 4, 2) # A-B both ladders strong; C below threshold
  )
  # edges oriented along the series: A->B and B->C
  e <- .nts_ms1_layer_edges_dt(
    nodes, valid_ids = nodes$name,
    ms1_layers = "homologue", cf2_ladder_boost = TRUE
  )
  ab <- e[from == "A" & to == "B"]
  bc <- e[from == "B" & to == "C"]
  expect_true(isTRUE(ab$ms2_ladder_confirmed))
  expect_equal(ab$edge_score, 1)
  expect_false(isTRUE(bc$ms2_ladder_confirmed))
  expect_true(is.na(bc$edge_score))
  # no new feature-level edges are invented: still exactly the series-adjacent pairs
  expect_equal(nrow(e), 2L)
})

test_that("cf2_ladder_boost = FALSE leaves MS1 edges unscored", {
  nodes <- data.table::data.table(
    name = c("A", "B"),
    mz = c(100, 100 + CF2),
    kmd_series_id = c("S1", "S1"),
    pfas_ladder_len = c(5, 5)
  )
  e <- .nts_ms1_layer_edges_dt(
    nodes, valid_ids = nodes$name,
    ms1_layers = "homologue", cf2_ladder_boost = FALSE
  )
  expect_false(isTRUE(e$ms2_ladder_confirmed))
  expect_true(is.na(e$edge_score))
})

test_that("cf2_ladder_boost is a no-op when the node table has no pfas_ladder_len", {
  nodes <- data.table::data.table(
    name = c("A", "B"),
    mz = c(100, 100 + CF2),
    kmd_series_id = c("S1", "S1")
  )
  e <- .nts_ms1_layer_edges_dt(
    nodes, valid_ids = nodes$name,
    ms1_layers = "homologue", cf2_ladder_boost = TRUE
  )
  expect_false(isTRUE(e$ms2_ladder_confirmed))
  expect_true(is.na(e$edge_score))
})

# ---------------------------------------------------------------------------
# S2-c: md_min_da parameterises the intra-spectrum spacing floor
# ---------------------------------------------------------------------------

test_that(".nts_network_ms2_massdiff_edges honours md_min_da", {
  clean <- data.table::data.table(
    name = c("A", "A", "B", "B"),
    mz = c(50, 50.4, 80, 80.4) # only a 0.4 Da spacing in each spectrum
  )
  # default floor 1.0 Da -> no spacing survives -> no edge
  e0 <- .nts_network_ms2_massdiff_edges(clean, mz_tol_ppm = 10, min_k = 1)
  expect_equal(nrow(e0), 0L)
  # floor lowered -> the shared 0.4 Da spacing connects A and B
  e1 <- .nts_network_ms2_massdiff_edges(clean, mz_tol_ppm = 10, min_k = 1, md_min_da = 0.3)
  expect_equal(nrow(e1), 1L)
  expect_equal(e1$Interaction, "MS2_Shared_Mass_Diff")
  expect_equal(e1$edge_score, 1) # 1 shared of 1 spacing in each spectrum
})

# ---------------------------------------------------------------------------
# S1-b: collapse aggregation of the new derived columns
# ---------------------------------------------------------------------------

test_that(".nts_merged_network_edge_collapse folds derived edge columns", {
  E <- data.table::data.table(
    from = c("A", "A"),
    to = c("B", "B"),
    edge_type = c("Homologue", "Homologue"),
    label = c("+2CF2", "+1CF2"),
    edge_score = c(NA_real_, 1),
    n_steps = c(2L, NA_integer_),
    transformation = c("+2CF2", "+1CF2"),
    mass_difference = c(99.9936128, 49.9968064),
    ms2_ladder_confirmed = c(FALSE, TRUE)
  )
  out <- .nts_merged_network_edge_collapse(E)
  expect_equal(nrow(out), 1L)
  expect_equal(out$edge_score, 1) # max finite
  expect_equal(out$n_steps, 2L) # first finite
  expect_equal(out$transformation, "+2CF2") # first non-empty
  expect_equal(out$mass_difference, 99.9936128, tolerance = 1e-6) # first finite
  expect_true(isTRUE(out$ms2_ladder_confirmed)) # any TRUE
})
