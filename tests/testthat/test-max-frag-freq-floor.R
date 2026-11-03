test_that("min_frag_freq keeps peaks with frequency >= threshold (floor)", {
  skip_if_not_installed("data.table")
  d <- tempfile("nts_mff")
  dir.create(d)
  p <- NTSminer::create_nts_project(d)
  p <- NTSminer::set_feature_table(p, data.frame(
    name = c("P", "C"),
    mz = c(400, 356.01),
    rt = c(1, 2),
    stringsAsFactors = FALSE
  ))
  p <- NTSminer::set_ms1_annotations(p, data.frame(
    name = c("P", "C"),
    mz = c(400, 356.01),
    rt = c(1, 2),
    annotation_level = c("Known", "Known"),
    stringsAsFactors = FALSE
  ))
  p <- NTSminer::set_ms2_annotations(p, data.frame(
    name = c("P", "C"),
    mz = c(400, 356.01),
    rt = c(1, 2),
    stringsAsFactors = FALSE
  ))
  cl <- list(
    P = list(
      feature_id = "P",
      precursorMZ = 400,
      rt = 1,
      spectrum = data.frame(
        mz = c(50, 60),
        intensity = c(1, 1),
        frequency = c(0.9, 0.9),
        stringsAsFactors = FALSE
      )
    ),
    C = list(
      feature_id = "C",
      precursorMZ = 356.01,
      rt = 2,
      spectrum = data.frame(
        mz = c(50, 60),
        intensity = c(1, 1),
        frequency = c(0.9, 0.9),
        stringsAsFactors = FALSE
      )
    )
  )
  methods::slot(p, "clean_ms2") <- cl
  g0 <- NTSminer::generate_network_data(
    p,
    ms2_edge_types = "fragment",
    ms1_edge_types = NULL,
    min_shared_ms2 = 2L,
    min_frag_freq = 0.8
  )
  expect_true(any(g0$edges$Interaction == "MS2_Shared_Fragment"))

  # Deprecated alias `max_frag_freq` still works but warns (renamed to
  # `min_frag_freq`; it has always been a *floor*, not a cap).
  expect_warning(
    g_alias <- NTSminer::generate_network_data(
      p,
      ms2_edge_types = "fragment",
      ms1_edge_types = NULL,
      min_shared_ms2 = 2L,
      max_frag_freq = 0.8
    ),
    "max_frag_freq.*renamed.*min_frag_freq"
  )
  expect_setequal(g_alias$edges$from, g0$edges$from)
  expect_setequal(g_alias$edges$to, g0$edges$to)

  cl_low <- list(
    P = list(
      feature_id = "P",
      precursorMZ = 400,
      rt = 1,
      spectrum = data.frame(
        mz = c(50, 60),
        intensity = c(1, 1),
        frequency = c(0.2, 0.2),
        stringsAsFactors = FALSE
      )
    ),
    C = list(
      feature_id = "C",
      precursorMZ = 356.01,
      rt = 2,
      spectrum = data.frame(
        mz = c(50, 60),
        intensity = c(1, 1),
        frequency = c(0.2, 0.2),
        stringsAsFactors = FALSE
      )
    )
  )
  methods::slot(p, "clean_ms2") <- cl_low
  expect_error(
    NTSminer::generate_network_data(
      p,
      ms2_edge_types = "fragment",
      ms1_edge_types = NULL,
      min_shared_ms2 = 2L,
      min_frag_freq = 0.8
    ),
    "min_frag_freq"
  )
})
