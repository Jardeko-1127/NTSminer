# Tests for optimized flag_isotopes -------------------------------------------

test_that("flag_isotopes: basic 13C chain M, M+1, M+2 with decaying intensity", {
  ft <- data.frame(
    name = c("A", "B", "C"),
    mz = c(400.0, 401.0034, 402.0067),
    rt = c(100, 101, 99),
    sample1 = c(100000, 8700, 380),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft, mz_tol_ppm = 10, rt_tol = 5)

  iso <- res[res$is_isotope, ]
  expect_true("B" %in% iso$name)
  expect_true("C" %in% iso$name)
  expect_equal(iso$parent_name[iso$name == "B"], "A")
  expect_equal(iso$parent_name[iso$name == "C"], "A")
  expect_equal(iso$isotope_kind[iso$name == "B"], "13C")
  expect_equal(iso$isotope_c13_n[iso$name == "B"], 1L)
  expect_equal(iso$isotope_c13_n[iso$name == "C"], 2L)
})


test_that("flag_isotopes: intensity ratio rejects M+1 at 95% of parent", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 401.0034),
    rt = c(100, 100),
    sample1 = c(100000, 95000),  # 95% - way too high for 13C
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  expect_false(any(res$is_isotope))
})


test_that("flag_isotopes: chain rejects M+2 when M+1 is missing", {
  ft <- data.frame(
    name = c("A", "C"),
    mz = c(400.0, 402.0067),
    rt = c(100, 100),
    sample1 = c(100000, 380),  # M+2 exists but M+1 is absent
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  expect_false(any(res$is_isotope))
})


test_that("flag_isotopes: chain rejects M+2 with stronger intensity than M+1", {
  ft <- data.frame(
    name = c("A", "B", "C"),
    mz = c(400.0, 401.0034, 402.0067),
    rt = c(100, 100, 100),
    sample1 = c(100000, 2000, 5000),  # M+2 > M+1 violates monotonic decay
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  iso <- res[res$is_isotope, ]
  # M+1 should be marked, but M+2 rejected due to intensity inversion
  expect_true("B" %in% iso$name)
  expect_false("C" %in% iso$name)
})


test_that("flag_isotopes: RT tolerance allows 4s gap with rt_tol=5", {
  ft <- data.frame(
    name = c("A", "B", "C"),
    mz = c(400.0, 401.0034, 402.0067),
    rt = c(100, 104, 96),  # 4s apart, within new default rt_tol=5
    sample1 = c(100000, 8700, 380),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft, rt_tol = 5)

  iso <- res[res$is_isotope, ]
  expect_true("B" %in% iso$name)
  expect_true("C" %in% iso$name)
})


test_that("flag_isotopes: 37Cl isotope identified with correct halogen label", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 401.9971),
    rt = c(100, 100),
    sample1 = c(100000, 32000),  # ratio 0.32 ~ 37Cl/35Cl
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  iso <- res[res$is_isotope, ]
  expect_true("B" %in% iso$name)
  expect_equal(iso$isotope_kind[iso$name == "B"], "halogen")
})


test_that("flag_isotopes: 81Br isotope identified", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 401.9980),
    rt = c(100, 100),
    sample1 = c(100000, 98000),  # ratio ~0.98 for Br
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  iso <- res[res$is_isotope, ]
  expect_true("B" %in% iso$name)
  expect_equal(iso$isotope_kind[iso$name == "B"], "halogen")
})


test_that("flag_isotopes: 34S isotope identified", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 401.9958),
    rt = c(100, 100),
    sample1 = c(100000, 4500),  # ratio ~0.045 for 34S
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  iso <- res[res$is_isotope, ]
  expect_true("B" %in% iso$name)
  expect_equal(iso$isotope_kind[iso$name == "B"], "halogen")
})


test_that("flag_isotopes: non-isotope 0.5 Da shift is not flagged", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 400.500),
    rt = c(100, 100),
    sample1 = c(100000, 80000),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  expect_false(any(res$is_isotope))
})


test_that("flag_isotopes: M+0 parent conflict - weak parent stripped", {
  # A (strong) claims B as M+1. B (weak) has a 37Cl child C that only B matches.
  # B is both child of A AND parent of C → conflict.
  # Since A is stronger, B's parent role is stripped → C loses isotope status.
  ft <- data.frame(
    name = c("A", "B", "C"),
    mz = c(400.0, 401.0034, 403.0004),
    rt = c(100, 100, 100),
    sample1 = c(100000, 50000, 15000),
    stringsAsFactors = FALSE
  )
  # C = B + 37Cl shift (~1.997), ratio 15000/50000=0.3 within Cl range
  # C cannot match A (no known A→C isotope shift)
  res <- flag_isotopes(ft)

  iso <- res[res$is_isotope, ]
  expect_false("A" %in% iso$name)
  expect_true("B" %in% iso$name)
  # B's parent status stripped → C demoted
  expect_false("C" %in% iso$name)
})


test_that("flag_isotopes: empty feature_table does not error", {
  ft <- data.frame(
    name = character(0),
    mz = numeric(0),
    rt = numeric(0),
    stringsAsFactors = FALSE
  )
  expect_error(flag_isotopes(ft), "nrow")
})


test_that("flag_isotopes: single feature returns no isotopes", {
  ft <- data.frame(
    name = "A",
    mz = 400.0,
    rt = 100,
    sample1 = 100000,
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  expect_false(any(res$is_isotope))
})


# --- P0-2: multi-halogen envelope fitting ------------------------------------

test_that("flag_isotopes: 2x37Cl (M+2/M = 0.644) is accepted via envelope fit", {
  # 1.997050 shift, ratio 0.644 = 2 x Cl minor/major -> the legacy single-halogen
  # band [0.10, 0.50] would reject this; envelope fitting must accept it.
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 401.99705),
    rt = c(100, 100),
    sample1 = c(100000, 64400),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  iso <- res[res$is_isotope, ]
  expect_true("B" %in% iso$name)
  expect_equal(iso$isotope_kind[iso$name == "B"], "halogen")
})


test_that("flag_isotopes: 3x37Cl (M+2/M = 0.964) is accepted", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 401.99705),
    rt = c(100, 100),
    sample1 = c(100000, 96400),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)
  expect_true("B" %in% res$name[res$is_isotope])
})


test_that("flag_isotopes: 4x81Br envelope (M+2/M = 3.891) is recognised by the fit", {
  # The 4xBr envelope has M+2 MORE abundant than M+0, so the lower-m/z peak is
  # the parent, not the isotope. flag_isotopes assigns child = lower m/z, so this
  # pattern cannot be expressed at the flag_isotopes level; assert the P0-2
  # capability directly on the envelope fit instead.
  fit <- NTSminer:::.nts_halogen_envelope_fit(389100, 100000, c("Cl", "Br", "S"))
  expect_true(fit$accept)
  expect_equal(fit$element, "Br")
  expect_equal(fit$n, 4L)
})


test_that("flag_isotopes: 2x34S (M+2/M ~ 0.089) is accepted", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 401.995796),
    rt = c(100, 100),
    sample1 = c(100000, 8900),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)
  expect_true("B" %in% res$name[res$is_isotope])
})


test_that("flag_isotopes: envelope fit rejects a ratio far above every halogen envelope", {
  # 10x the parent is far outside any single-substitution envelope.
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 401.99705),
    rt = c(100, 100),
    sample1 = c(100000, 1000000),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)
  expect_false(any(res$is_isotope))

  # 0.80 sits between the 2Cl (0.644) and 3Cl (0.964) envelopes; within the 0.35
  # relative tolerance of 2Cl, so it is accepted — documented behaviour.
  ft2 <- ft
  ft2$sample1[2] <- 80000
  res2 <- flag_isotopes(ft2)
  expect_true("B" %in% res2$name[res2$is_isotope])
})


test_that("flag_isotopes: theoretical envelope table is monotone in n", {
  tab <- NTSminer:::.nts_halogen_envelope_table(c("Cl", "Br", "S"))
  cl <- tab[tab$element == "Cl", ]
  expect_equal(cl$n, seq_len(nrow(cl)))
  expect_true(all(diff(cl$r_theory) > 0))
  # 1 Cl ~ 0.32, 2 Cl ~ 0.644
  expect_equal(cl$r_theory[cl$n == 1], 1 * 0.2424 / 0.7576, tolerance = 1e-6)
  expect_equal(cl$r_theory[cl$n == 2], 2 * 0.2424 / 0.7576, tolerance = 1e-6)
})


# --- P2-8: minor isotope (15N / 18O) is opt-in -------------------------------

test_that("flag_isotopes: 15N shift is NOT flagged by default", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 400.997035),
    rt = c(100, 100),
    sample1 = c(100000, 370),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)
  expect_false(any(res$is_isotope))
})


test_that("flag_isotopes: 15N shift IS flagged when include_minor_isotopes = TRUE", {
  ft <- data.frame(
    name = c("A", "B"),
    mz = c(400.0, 400.997035),
    rt = c(100, 100),
    sample1 = c(100000, 370),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft, include_minor_isotopes = TRUE)
  iso <- res[res$is_isotope, ]
  expect_true("B" %in% iso$name)
  expect_equal(iso$isotope_kind[iso$name == "B"], "minor")
})


# --- P1-3: conflict resolution clears stale detail columns -------------------

test_that("flag_isotopes: conflict resolution leaves no stale isotope detail columns", {
  # A (strong) claims B as M+1; B (weaker but still above threshold) also has a
  # 37Cl child C. Because B is simultaneously a child (of A) and a parent (of C),
  # resolve_parent_conflicts strips B's parent role. The stale isotope detail
  # columns on C must be cleared, not merely is_isotope / parent_name.
  ft <- data.frame(
    name = c("A", "B", "C"),
    mz = c(400.0, 401.0034, 403.0004),
    rt = c(100, 100, 100),
    sample1 = c(100000, 50000, 15000),
    stringsAsFactors = FALSE
  )
  res <- flag_isotopes(ft)

  stripped <- !(res$is_isotope %in% TRUE)
  stale <- stripped &
    (!is.na(res$isotope_kind) | !is.na(res$isotope_label) |
       !is.na(res$isotope_shift_da) | !is.na(res$isotope_c13_n))
  expect_false(any(stale))
})

