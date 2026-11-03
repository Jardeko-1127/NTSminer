# Regression tests for the rcdk descriptor computation path.
#
# rcdk 3.8.2 / CDK exposes `LongestAliphaticChainDescriptor`, which recurses
# infinitely in native code on any molecule that has an alkyl chain attached to
# an aromatic ring (toluene, benzoic acid, aspirin, ...). The result is a hard
# "segfault from C stack overflow" that kills the R process and cannot be caught
# by `tryCatch()`. `.nts_compute_descriptors()` must therefore drop it from the
# default descriptor set.

test_that(".NTS_UNSAFE_DESCRIPTORS contains the crash-prone descriptor", {
  expect_true(
    "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor" %in%
      NTSminer:::.NTS_UNSAFE_DESCRIPTORS
  )
})

test_that(".nts_compute_descriptors does not crash on aromatic-chain molecules", {
  skip_if_not_installed("rcdk")
  # Previously: "segfault from C stack overflow" (process death).
  dt <- NTSminer:::.nts_compute_descriptors(
    c("c1ccccc1", "CC(=O)OC1=CC=CC=C1C(=O)O", "CCCC")
  )
  expect_true(nrow(dt) >= 1L)
  expect_true("smiles" %in% names(dt))
  # The offending descriptor's column must not be present.
  expect_false(any(grepl("nAtomLAC", names(dt))))
})

test_that(".nts_compute_descriptors honours an explicit desc_names that excludes the unsafe ones", {
  skip_if_not_installed("rcdk")
  safe <- setdiff(rcdk::get.desc.names(), NTSminer:::.NTS_UNSAFE_DESCRIPTORS)
  dt <- NTSminer:::.nts_compute_descriptors("CC(=O)OC1=CC=CC=C1C(=O)O", desc_names = safe)
  expect_true(nrow(dt) == 1L)
})
