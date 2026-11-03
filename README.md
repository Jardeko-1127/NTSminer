# NTSminer

NTSminer is an R package for non-target screening (NTS) of complex environmental samples using LC–MS(/MS). It supports feature extraction and curation, spectral annotation, and **multi-layer molecular network construction and annotation propagation**. Network layers can represent CF2 homologue relationships, diagnostic-fragment similarity, neutral losses, exact-mass differences, and MS2 semantic similarity.

- **Version:** 0.1.41 (see [`DESCRIPTION`](DESCRIPTION)); release history is recorded in [`NEWS.md`](NEWS.md).
- **License:** MIT (see [`LICENSE`](LICENSE)).

## Installation

Install the package from GitHub. NTSminer depends on R and several Bioconductor packages; install the required dependencies for your platform before or during package installation.

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("Jardeko-1127/NTSminer", dependencies = TRUE)
```

To load the package source during development:

```r
pkgload::load_all()
```

## Working with your data

The repository does not include raw demonstration mzML files. Supply paths to your own MS1 and MS2 data. The package includes commonly used PFAS reference libraries, a diagnostic-fragment table, and model assets.

```r
ms1_dir <- "path/to/MS1"
ms2_dir <- "path/to/MS2"

# Load the bundled in-house, OECD online, and OECD in-silico PFAS libraries.
# By default, the three libraries are combined into one NTSDatabase object.
db <- load_builtin_pfas_library()

diagnostic_table <- bundled_nts_file(
  "extdata/tables/diagnostic_fragments.csv"
)
```

A typical LC–MS/MS workflow includes feature extraction, isotope and in-source-fragmentation checks, CF2-based KMD analysis, MS2 extraction and cleanup, diagnostic-fragment validation, library matching, and network construction. See the [workflow vignette](vignettes/ntsminer-workflow.Rmd) for function examples. Raw MS1/MS2 data and project-specific results are supplied by the user.

## Multi-layer molecular networks

NTSminer can construct and inspect multiple complementary network layers, including:

- **CF2 homologue layer:** connects features consistent with a CF2 mass increment and homologue-series constraints.
- **Diagnostic-fragment layer:** links features that share characteristic MS2 fragments.
- **Neutral-loss layer:** links features with shared neutral losses.
- **Mass-difference layer:** links features by exact precursor-mass differences.
- **MS2 semantic-similarity layer:** connects spectra with similar learned embeddings.

Layers can be examined separately or combined for evidence integration and annotation propagation. Structural evidence and spectral similarity provide complementary paths for organizing and prioritizing candidates; propagated annotations should be interpreted with their supporting evidence and confidence levels.

## Candidate ranking and spectral models

The package includes tools for library matching, candidate merging and ranking, forward MS2 spectrum prediction, joint-fingerprint prediction and similarity scoring, and semantic spectrum search. Prediction features use the model assets shipped under `inst/extdata/models/`; some workflows require an available Python environment and optional dependencies.

## Semi-quantification and risk assessment

NTSminer includes internal-standard matching and response normalization, relative ionization efficiency (RIE) prediction, concentration estimation, EPI Suite and ADMET risk workflows, and ToxPi scoring and visualization.

### logIE and RIE example

```r
model <- load_logIE_model(
  bundled_nts_file("extdata/models/logIE/logIE_model_inhouse_pos.rds")
)
pred <- predict_logIE(smiles_vec, model = model, polarity = "positive")
rie <- predict_rie(feature_table, logIE_model = model, polarity = "positive")
```

### ADMET and ToxPi example

These functions require the relevant backend services or optional packages. For example:

```r
smiles_tbl <- utils::read.csv(
  "path/to/compound_smiles.csv",
  stringsAsFactors = FALSE
)
smiles_canon <- canonicalize_smiles(smiles_tbl$SMILES)
epi <- predict_episuite(smiles_canon)
admet <- predict_admet_risk(smiles_canon, backend = "admet_ai")
tox_input <- prepare_toxpi(
  epi,
  admet,
  query_smiles = smiles_canon,
  compound_ids = smiles_tbl$name
)
tox_scores <- run_toxpi(tox_input, weight_multiplier = 1)
plot_toxpi(tox_scores)
```

## Shiny interface and package checks

```r
NTSminer::launch_ntsminer()

devtools::check()
devtools::test()
```

`devtools::check()` and `devtools::test()` require the package development dependencies. Optional services and tools such as Python-based models, EPI Suite, ADMET backends, and ToxPi have their own installation requirements.
