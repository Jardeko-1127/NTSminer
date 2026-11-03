test_that("annotate_by_dual_network syncs ref_id, Compound.name, SMILES on propagated nodes", {
  skip_if_not_installed("rcdk")
  skip_if_not_installed("fingerprint")
  skip_if_not_installed("data.table")
  has_rcdk <- suppressWarnings(tryCatch(
    !is.null(rcdk::parse.smiles("C")[[1L]]),
    error = function(e) FALSE
  ))
  if (!isTRUE(has_rcdk)) {
    skip("rcdk / CDK not available (JVM)")
  }

  tmpd <- tempfile("nts_id_sync")
  dir.create(tmpd)
  csvf <- file.path(tmpd, "lib.csv")
  writeLines(
    c("ID,Compound.name,SMILES", "P1,Perfluoro A,CC", "ancX,Other,CCC"),
    csvf
  )
  rdsf <- file.path(tmpd, "lib.rds")
  build_structural_library(
    csvf,
    smiles_col = "SMILES",
    id_col = "ID",
    output_rds = rdsf
  )
  lib <- readRDS(rdsf)
  lib_df <- as.data.frame(lib, stringsAsFactors = FALSE)
  em_p <- as.numeric(lib_df[lib_df[["ID"]] == "P1", "ExactMass", drop = TRUE][1L])
  proton <- 1.007276
  mz_u <- as.numeric(em_p - proton)

  net <- list(
    nodes = data.frame(
      id = c("anc1", "u1"),
      mz = c(30.0, mz_u),
      rt = c(1, 2),
      ref_id = c("SEED_REF", NA_character_),
      `Compound.name` = c("Seed compound", NA_character_),
      SMILES = c(NA_character_, NA_character_),
      annotation_level = c("Database Confirmed (Anchor)", "Unknown candidate"),
      Library_ID = c("P1", NA_character_),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(
      from = "anc1",
      to = "u1",
      Interaction = "MS2_Shared_Fragment",
      stringsAsFactors = FALSE
    ),
    polarity = "negative"
  )

  out <- annotate_by_dual_network(
    net,
    rdsf,
    mz_tol_ppm = 50,
    min_similarity = 0.2,
    similarity_method = "dice",
    polarity = "auto"
  )
  u1 <- out$nodes[out$nodes$id == "u1", , drop = FALSE]
  anc <- out$nodes[out$nodes$id == "anc1", , drop = FALSE]

  expect_equal(as.character(u1$ref_id[[1L]]), "P1")
  expect_equal(as.character(u1$Compound.name[[1L]]), "Perfluoro A")
  expect_equal(as.character(u1$SMILES[[1L]]), "CC")
  expect_equal(as.character(u1$propagated_ID[[1L]]), "P1")
  expect_equal(as.character(anc$ref_id[[1L]]), "SEED_REF")
  expect_equal(as.character(anc$Compound.name[[1L]]), "Seed compound")
})
