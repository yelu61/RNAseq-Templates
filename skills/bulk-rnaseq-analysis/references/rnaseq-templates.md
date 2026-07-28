# RNAseq-Templates Backend

## Repository contract

Validate the root by checking for:

- `RNAseq_lib/`
- `notebooks/`
- `references/TEMPLATE_SELECTION.md`
- `examples/run_demo_smoke_test.R`

Read the repository's `README.md` and
`references/TEMPLATE_SELECTION.md` before selecting a template. The repository
currently targets downstream analysis from expression matrices, not FASTQ.

## Template selection

| Question | Entry point |
| --- | --- |
| Standard group comparison | `notebooks/RNAseq_General.ipynb` |
| limma-voom or batch-heavy contrast | `notebooks/RNAseq_limma_voom_Template.ipynb` |
| Longitudinal/time-course | `notebooks/RNAseq_TimeCourse_Template.ipynb` |
| Immune/stromal deconvolution | `notebooks/RNAseq_TME_Deconvolution_Template.ipynb` |
| Co-expression modules | `notebooks/RNAseq_WGCNA_Template.ipynb` |
| Lightweight TCGA/GEO teaching workflow | `notebooks/RNAseq_TCGA_GEO_Template.ipynb` |
| Headless standard pipeline | `templates/General/run_analysis.R` |

For advanced TCGA/TARGET/GTEx analysis, route to the TCGA toolkit instead of
expanding the lightweight teaching notebook.

## Standard workflow

1. Identify raw counts versus normalized expression and inspect metadata.
2. Confirm sample identifiers match exactly and groups have adequate size.
3. Copy the selected notebook or General CLI template into the analysis
   project; keep `RNAseq_lib/` resolvable.
4. Edit the parameter block/configuration before changing analysis code.
5. Run input validation before fitting a model.
6. Execute and retain the native output structure.
7. Record `sessionInfo.txt`, parameters, comparisons, and warnings.

For the General CLI runner:

```bash
Rscript templates/General/run_analysis.R
```

Use a project-local copy of `config.R`, `run_analysis.R`, and
`visualize_results.R` when the analysis must not write into the template
repository.

## Input-scale rules

- DESeq2: integer-like raw counts.
- limma-voom: raw counts with library-size modelling.
- WGCNA/time-course clustering: variance-stabilized or otherwise appropriate
  normalized expression.
- TME: TPM, `log2(TPM + 1)`, or raw counts plus gene lengths as documented.
- GEO: inspect platform annotation and preprocessing; do not assume downloaded
  values are raw counts.

## Validation

Run focused unit tests after helper-library changes:

```bash
Rscript tests/testthat.R
```

Validate notebook structure:

```bash
Rscript tests/validate_notebooks.R
```

Run the bundled smoke workflow when dependencies and network/cache conditions
allow:

```bash
Rscript examples/run_demo_smoke_test.R
```
