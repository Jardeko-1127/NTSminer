# NTSminer 0.1.41 submission release

This repository is a curated source snapshot prepared from the NTSminer development tree for manuscript-related public release.

## Included

- Core R package source, generated help pages, namespace, and workflow vignette.
- Runtime PFAS reference libraries and diagnostic/suspect tables.
- Inference assets needed by packaged prediction features (forward-MS2, joint fingerprint, RT, and logIE).
- Package-level `testthat` tests and the Sirius parsing fixture.

## Excluded

- Raw MS1/MS2 demonstration mzML files and demo pipeline scripts.
- Development-only training datasets, model training/evaluation outputs, and one-off analysis scripts.
- Historical/duplicate PFAS database exports and example internal-standard CSVs.
- Raw-data integration runs and research project outputs.

## Release notes

- Package version is `0.1.41`, read from `DESCRIPTION`.
- The source archive was built with R 4.4.1 using `R CMD build --no-build-vignettes --no-manual`.
- The package build completed. The test suite was not run as part of this packaging task.
- The largest included model weight file is below GitHub's 100 MB per-file limit; the generated source tarball is not tracked in Git.
- No GitHub remote was configured or authenticated when this snapshot was prepared. The repository is ready for connection to the user's chosen GitHub repository after the target URL and GitHub authentication are available.
