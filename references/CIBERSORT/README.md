# CIBERSORT resources

This directory contains the standalone CIBERSORT implementation and signature matrices used by `RNAseq_lib/tme_utils.R`.

The native implementation runs locally and **does not depend on IOBR**. IOBR's
CIBERSORT integration is a separate execution path; enabling or installing it
is not required to call this bundled implementation.

## Files

- `CIBERSORT.R` — Standalone source for `cibersort()`, `CoreAlg()`, and `doPerm()`. Sourced by `run_native_cibersort()`.
- `LM22.txt` — Human peripheral blood leukocyte signature matrix (22 cell types).
- `cibersort_mouse_22.csv` — Mouse immune-cell signature matrix (25 cell types; MGI gene symbols).

With `perm = 0`, permutation P values are unavailable and exported as `NA`;
they must not be interpreted as a significance test. Single-sample matrices
retain their sample and cell-type dimensions.

## Usage

See `run_native_cibersort()` in `RNAseq_lib/tme_utils.R`. The wrapper loads the signatures, prepares the expression matrix, and calls the CIBERSORT core algorithm.

Match the expression identifiers to the selected reference:

| Native input | Signature | Ortholog conversion before this call |
|---|---|---|
| Mouse expression with MGI symbols | `cibersort_mouse_22.csv` | Do not convert to human symbols |
| Human expression with HGNC symbols | `LM22.txt` | Not needed |

Pass the mouse signature explicitly: the wrapper's default reference is human
LM22. The native engine matches symbols case-sensitively. Expression scale,
reference coverage and suitability for the tissue still require validation.

The TME runner now retains a native-species matrix for this engine and builds
human orthologs separately only for human-reference methods. The incorrect
v0.11.0 mouse-runner outputs must be recomputed; see the
[freeze audit](../FREEZE_AUDIT.md).

`run_native_cibersort(is_log = FALSE)` declares linear abundance. With
`is_log = TRUE`, the wrapper first reverses log2(TPM+1). It then calls the bundled
engine with `mixture_scale = "linear"`, so low-valued linear data are never
exponentiated a second time. Direct `cibersort()` callers may explicitly set
`mixture_scale = "linear"` or `"log2"` (plain log2 abundance); the default
`"auto"` preserves the historical max<50 heuristic for compatibility. New code
should always declare its scale. Custom scripts lacking this argument are
rejected by the wrapper rather than silently reinterpreting their input.

The wrapper records reference overlap and stops when matched genes do not
outnumber reference cell types. This is a minimal fitting guard, not a biological
validation threshold. Human/mouse native regression tests use their actual
bundled references without IOBR. Automatic native-vs-IOBR comparison is limited
to human inputs with identical verified reference matrices and QN settings;
mouse 25-class and human 22-class fractions are not equivalent categories.

## Dependencies

`CIBERSORT.R` requires the R packages `e1071`, `preprocessCore`, `future`, `furrr`, and `purrr`. These should be installed via `install_dependencies.R`.

## License

The CIBERSORT package is licensed under AGPL-3 (see the original package `DESCRIPTION`). The signature matrices are from the original CIBERSORT publication and mouse signature resources; use them in accordance with their respective licenses.
