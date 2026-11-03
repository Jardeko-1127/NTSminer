# Tests for in-source fragment (ISF) recognition ------------------------------
#
# Covers the P0-1 (isotope exclusion), P1-5 (acyclic parent chains) and P1-6
# (unified defaults) fixes to flag_isf() / recognize_isf().
#
# .flag_isf_df() returns `is_isf`, `parent_mz` and `isf_annotation`; the parent
# is identified by its m/z (parent_mz), not by name.

.CF2 <- 49.9968064


test_that("ISF defaults are unified between flag_isf and recognize_isf", {
  expect_equal(NTSminer:::.ISF_DEFAULT_RT_TOL, 2)
  expect_equal(NTSminer:::.ISF_DEFAULT_CORR_THRESH, 0.85)
  expect_equal(NTSminer:::.ISF_DEFAULT_LOSS_TOL, 0.005)

  # Both exported entry points must default to the shared constants (the
  # formals are symbols referring to the constants, so evaluate them).
  expect_equal(eval(formals(flag_isf)$rt_tol), NTSminer:::.ISF_DEFAULT_RT_TOL)
  expect_equal(eval(formals(flag_isf)$corr_thresh), NTSminer:::.ISF_DEFAULT_CORR_THRESH)
  expect_equal(eval(formals(flag_isf)$loss_tol_da), NTSminer:::.ISF_DEFAULT_LOSS_TOL)
  expect_equal(eval(formals(recognize_isf)$rt_window), NTSminer:::.ISF_DEFAULT_RT_TOL)
  expect_equal(eval(formals(recognize_isf)$corr_thresh), NTSminer:::.ISF_DEFAULT_CORR_THRESH)
  expect_equal(eval(formals(recognize_isf)$loss_tol), NTSminer:::.ISF_DEFAULT_LOSS_TOL)
})


test_that("flag_isf detects a CF2 neutral-loss pair with correlated abundance", {
  ft <- data.frame(
    name = c("P", "F"),
    mz = c(400.0, 400.0 - .CF2),
    rt = c(100, 100),
    s1 = c(1000, 300),
    s2 = c(2000, 600),
    s3 = c(4000, 1200),
    stringsAsFactors = FALSE
  )
  res <- NTSminer:::.flag_isf_df(ft, neutral_losses = c(CF2 = .CF2))

  expect_true(isTRUE(res$is_isf[res$name == "F"]))
  expect_equal(res$parent_mz[res$name == "F"], 400.0)
  expect_false(isTRUE(res$is_isf[res$name == "P"]))
})


test_that("flag_isf excludes isotope peaks from the fragment role (P0-1)", {
  ft <- data.frame(
    name = c("P", "F"),
    mz = c(400.0, 400.0 - .CF2),
    rt = c(100, 100),
    s1 = c(1000, 300),
    s2 = c(2000, 600),
    s3 = c(4000, 1200),
    is_isotope = c(FALSE, TRUE),  # F is a genuine fragment BUT labelled isotope
    stringsAsFactors = FALSE
  )
  res_excl <- NTSminer:::.flag_isf_df(ft, neutral_losses = c(CF2 = .CF2))
  expect_false(isTRUE(res_excl$is_isf[res_excl$name == "F"]))

  # Without the isotope mask the same peak is flagged (demonstrates the fix).
  ft_no_iso <- ft
  ft_no_iso$is_isotope <- FALSE
  res_all <- NTSminer:::.flag_isf_df(ft_no_iso, neutral_losses = c(CF2 = .CF2))
  expect_true(isTRUE(res_all$is_isf[res_all$name == "F"]))
})


test_that("flag_isf excludes isotope peaks from the parent role (P0-1)", {
  ft <- data.frame(
    name = c("P", "F"),
    mz = c(400.0, 400.0 - .CF2),
    rt = c(100, 100),
    s1 = c(1000, 300),
    s2 = c(2000, 600),
    s3 = c(4000, 1200),
    is_isotope = c(TRUE, FALSE),  # P (the parent) is labelled an isotope
    stringsAsFactors = FALSE
  )
  res <- NTSminer:::.flag_isf_df(ft, neutral_losses = c(CF2 = .CF2))
  expect_false(any(res$is_isf %in% TRUE))
})


test_that("flag_isf enforces acyclic parent chains (P1-5)", {
  # P -> F1 -> F2, each a CF2 loss. F1 is both a child (of P) and a parent (of
  # F2); the acyclic guard must demote F1's parent role so no feature is ever a
  # fragment and a parent at the same time.
  ft <- data.frame(
    name = c("P", "F1", "F2"),
    mz = c(400.0, 400.0 - .CF2, 400.0 - 2 * .CF2),
    rt = c(100, 100, 100),
    s1 = c(1000, 800, 600),
    s2 = c(2000, 1600, 1200),
    s3 = c(4000, 3200, 2400),
    stringsAsFactors = FALSE
  )
  res <- NTSminer:::.flag_isf_df(ft, neutral_losses = c(CF2 = .CF2))

  kids_mz <- res$mz[res$is_isf %in% TRUE]
  pars_mz <- res$parent_mz[res$is_isf %in% TRUE]
  expect_length(intersect(kids_mz, pars_mz), 0L)
})


test_that("recognize_isf and flag_isf agree on a simple fragment pair", {
  ft <- data.frame(
    name = c("P", "F"),
    mz = c(400.0, 400.0 - .CF2),
    rt = c(100, 100),
    s1 = c(1000, 300),
    s2 = c(2000, 600),
    s3 = c(4000, 1200),
    stringsAsFactors = FALSE
  )
  a <- NTSminer:::.flag_isf_df(ft, neutral_losses = c(CF2 = .CF2))
  b <- NTSminer:::.recognize_isf_df(ft, loss_tol = 0.005, rt_window = 2, corr_thresh = 0.85)

  expect_true(isTRUE(a$is_isf[a$name == "F"]))
  expect_true(isTRUE(b$is_isf[b$name == "F"]))
  expect_equal(
    a$parent_mz[a$name == "F"],
    b$parent_mz[b$name == "F"]
  )
})


test_that("flag_isf: empty feature table errors cleanly", {
  ft <- data.frame(name = character(0), mz = numeric(0), rt = numeric(0))
  expect_error(NTSminer:::.flag_isf_df(ft, neutral_losses = c(CF2 = .CF2)), "nrow")
})
