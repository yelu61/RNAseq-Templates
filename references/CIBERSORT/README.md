# CIBERSORT resources

This directory contains the standalone CIBERSORT implementation and signature matrices used by `RNAseq_lib/tme_utils.R`.

## Files

- `CIBERSORT.R` — Standalone source for `cibersort()`, `CoreAlg()`, and `doPerm()`. Sourced by `run_native_cibersort()`.
- `LM22.txt` — Human peripheral blood leukocyte signature matrix (22 cell types).
- `cibersort_mouse_22.csv` — Mouse immune-cell signature matrix (25 cell types; MGI gene symbols).

## Usage

See `run_native_cibersort()` in `RNAseq_lib/tme_utils.R`. The wrapper loads the signatures, prepares the expression matrix, and calls the CIBERSORT core algorithm.

## Dependencies

`CIBERSORT.R` requires the R packages `e1071`, `preprocessCore`, `future`, `furrr`, and `purrr`. These should be installed via `install_dependencies.R`.

## License

The CIBERSORT package is licensed under AGPL-3 (see the original package `DESCRIPTION`). The signature matrices are from the original CIBERSORT publication and mouse signature resources; use them in accordance with their respective licenses.
