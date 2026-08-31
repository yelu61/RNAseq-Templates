# RNAseq_TME — command-line template

> **Maintenance correction:** this checkout separates native-species and
> human-reference matrices, fixes explicit linear input to native CIBERSORT,
> and uses the actual ssGSEA algorithm. Mouse native-CIBERSORT and files labelled
> `ssGSEA_*` from v0.11.0 must be regenerated in a new run. Local regression tests
> cover human/mouse native execution; they do not validate every IOBR method or
> biological interpretation. See the [freeze audit](../../references/FREEZE_AUDIT.md).

The bundled native CIBERSORT engine itself does **not** require IOBR. Its
mouse reference is `cibersort_mouse_22.csv`, and its human reference is `LM22.txt`.
See [standalone CIBERSORT usage](../../references/CIBERSORT/README.md).

Non-interactive tumor-microenvironment (TME) deconvolution pipeline. It shares
the analysis intent and helper library of `notebooks/RNAseq_TME_Deconvolution_Template.ipynb`;
the production runner is authoritative: expression →
ESTIMATE / IOBR / native CIBERSORT → ssGSEA immune-signature scoring.

## Files

| File | Purpose |
|------|---------|
| `config.R` | **The only file you edit.** Paths, input mode, samples/groups, method switches. |
| `run_analysis.R` | Full pipeline: expression → TPM → ESTIMATE / IOBR / native CIBERSORT → ssGSEA → plots. |
| `visualize_results.R` | Cheap re-plotting from saved results (no deconvolution/ssGSEA recompute). |

## Usage

```bash
# 1. Copy this folder into your project, then edit config.R.
# 2. Run the standard analysis:
Rscript run_analysis.R                    # uses ./config.R
Rscript run_analysis.R path/to/config.R   # uses a different config

# 3. (Optional) iterate on figures without re-running the pipeline:
Rscript visualize_results.R
```

The runner preserves the invocation directory. Relative inputs and `OUTDIR`
resolve from that run root, while a relative config argument is resolved before
loading. For production, invoke it from a fresh `analysis/runs/<run_id>/`.

`RNAseq_lib` is located automatically (env `RNASEQ_LIB_DIR`, a `./RNAseq_lib`
next to the script, the repo root via `rprojroot`, or `../RNAseq_lib`). The
bundled CIBERSORT references are located via `rprojroot`
(`references/CIBERSORT/`) unless `CIBERSORT_SCRIPT` / `CIBERSORT_SIGNATURE` are
set explicitly in `config.R`.

## Input modes

- **`raw_counts` (default, recommended).** A featureCounts-like table with gene
  lengths. Prefer the counting annotation's featureCounts `Length` or an
  equivalent exonic length via `GENE_LENGTH_COLUMN`; a genomic gene start/end
  span is not generally the same length. The coordinate fallback is appropriate
  only when that span matches the counted feature. Counts are converted to TPM
  with `build_tme_tpm()` over the full counted feature universe before ID mapping
  or filtering; stable feature IDs must be unique.
- **`expression`.** A pre-computed TPM or log2(TPM+1) matrix (`EXPR_UNIT`).
  VST/rlog is rejected — it cannot be converted back to a linear scale and is
  invalid for ESTIMATE/CIBERSORT/EPIC/xCell.

For human-signature methods, this runner converts mouse symbols to human
orthologs using the installed `babelgene` data; mouse Ensembl-to-symbol
conversion uses the input `gene_name` annotation when present, otherwise the
local annotation package. General's optional TME block uses the same preparation
functions. An explicit `ORTHOLOG_CACHE` (`TME_ORTHOLOG_CACHE` in General) retains
the legacy biomaRt/cache route. Check mapping losses and signature overlap
before interpreting scores. Human-symbol input can run the native
ESTIMATE/ssGSEA path offline with its packages installed. IOBR additionally
needs method reference bundles, and methods such as xCell may perform online
lookups even after reference caching.

## Species and symbol routing: required contract

The production runner and TME notebook implement these separate branches.

Distinguish within-species ID annotation (Ensembl → MGI/HGNC symbol) from
cross-species ortholog mapping (mouse MGI → human HGNC). Changing symbol case
does not perform ortholog mapping.

1. Prepare the annotation table and confirm species, feature IDs and length
   definition. Keep original IDs and input files unchanged.
2. Compute TPM from original counts and matching **same-species** feature
   lengths before renaming/collapsing expression rows. Existing TPM skips this
   step; explicitly declared log2(TPM+1) is converted back to linear TPM.
   VST/rlog is not an alternative source of TPM. Do not restrict input to DEGs.
3. If needed, convert Ensembl IDs to symbols **within the original species**.
   Already-valid MGI/HGNC symbols need no such conversion. The mapping table
   records original IDs, chosen native/human symbols, unmapped IDs and source.
   The local annotation helper selects the first mapping; duplicate symbols
   retain the row with highest mean linear expression. Review these choices.
4. Branch by the reference used by each method. Keep an original-species
   expression matrix; create a separate human-ortholog matrix only when a
   selected method requires a human reference. Do not overwrite the original.
5. Check actual reference-gene overlap and mapping losses before fitting or
   scoring, then record the matrix/reference species, versions and checksums
   with each result.

| Analysis branch | Expression / reference | Cross-species mapping |
|---|---|---|
| Human native CIBERSORT | HGNC / `LM22.txt` | None |
| Mouse native CIBERSORT | MGI / `cibersort_mouse_22.csv` | None; preserve mouse symbols |
| Human ESTIMATE / configured human-reference IOBR methods | HGNC / method-specific human reference | None |
| Mouse data using those human-reference methods | Separate HGNC ortholog matrix / human reference | Only after the native-species branch is preserved, before reference matching |
| Current bundled 28 immune gene sets | HGNC / human-symbol gene sets | Mouse data require ortholog mapping for this gene-set collection |
| A separately validated mouse gene-set collection | MGI / mouse gene sets | None; this is not the current bundled collection |

Human-reference analysis of mouse data is an ortholog-based extrapolation.
Mapping success alone does not validate the human reference in mouse tissue.
Record the limitation and coverage; do not treat scores as measured cell counts.

Both TME entry points now compute TPM before DEG/biotype filtering and share
native-symbol/ortholog preparation. General retains its final selected sample
set and reads the original input afresh. Match the feature universe, length
column, annotation, sample set and ortholog-cache choice to obtain the same
inputs; General still has no native-CIBERSORT branch. See
[babelgene's documented offline mapping](https://cran.r-project.org/web/packages/babelgene/vignettes/babelgene-intro.html).

Both runners and the TME notebook call the shared `run_tme_ssgsea()`, which
constructs `ssgseaParam(normalize = TRUE)` and records gene-set coverage. Its
scores are enrichment scores, not cell proportions. Previous v0.11.0 files
used `gsvaParam()` and require recomputation; renaming them is not sufficient.
The bundled collection still uses human symbols. See the
[GSVA method documentation](https://bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html).

Automatic native-vs-IOBR comparison is skipped for mouse projects: the mouse reference has 25
cell categories, whereas human LM22 has 22 with different definitions. A
predefined broad-cell mapping may support exploratory trend comparison, not
quantitative equivalence. Human native/IOBR results with the same LM22 and
preprocessing can be checked for implementation agreement; this is not
independent biological validation.

## Output layout

Outputs use the numbered run layout (set `OUTDIR <- "."` in config.R to write
these at the project root; set a path to nest the run under a directory):

```
0-Config/analysis_config_used.R     frozen config snapshot (reproducibility)
4-TME/                              TPM_matrix.csv, ESTIMATE_scores.csv,
                                    TPM_native_symbols.csv, TPM_human_symbols.csv,
                                    TME_gene_mapping.csv, TME_input_coverage.csv,
                                    IOBR_<method>.csv + IOBR_TME_combined.csv,
                                    CIBERSORT_native_*, ssGSEA_*, all PDFs,
                                    tme_results.Rdata (for visualize_results.R)
Analysis_summary.txt                text summary of the run
sessionInfo.txt                     package versions (reproducibility)
run_manifest.csv                    inputs/config/code signatures + lifecycle
```

## Optional methods and graceful degradation

ESTIMATE, IOBR, and native CIBERSORT are **optional**. Each is gated by a
`RUN_*` switch in `config.R`. A missing package or bundled native-CIBERSORT
resource disables that module with a clear message:

- `RUN_ESTIMATE` needs the `estimate` package.
- `RUN_IOBR` needs the `IOBR` package. On first use, recent IOBR releases may
  also fetch method reference bundles (`common_genes`, `TRef`, `xCell.data`,
  `lm22`). Populate that cache while online before an offline production run.
  If `RUN_IOBR=TRUE` and every requested method fails, the runner stops rather
  than silently reporting a complete deconvolution.
- `RUN_CIBERSORT` needs the bundled `references/CIBERSORT/CIBERSORT.R` script,
  a signature matrix (auto LM22 / `cibersort_mouse_22.csv`), and the
  `e1071`/`preprocessCore`/`future`/`furrr`/`purrr` packages.

ssGSEA (sections 8–9) needs `GSVA`; it too degrades gracefully if absent.

The bundled CI smoke test is deterministic and offline: it exercises native
ESTIMATE plus ssGSEA with `RUN_IOBR=FALSE`. IOBR should be validated separately
on a machine where its method data are already cached.

For production, follow the [fresh-run example](../../BEST_PRACTICES.md#1--reproducible-run).
Point expression/metadata paths to the previous run or original inputs using
absolute paths; do not copy old `1-DEG/` outputs into the new run root.

The runner rejects existing native output artifacts before writing a new
configuration snapshot. Completed runs also contain `run_inputs.csv` and
`run_runtime.csv`; unknown references make a run ineligible for duplicate
grouping. These files record provenance, not a restorable environment.

## Notes

- Both the runner and notebook prepare plotting metadata before native
  CIBERSORT, including runs that disable all human-reference methods.
- `SPECIES` selects the native reference and preserves matching native symbols;
  human-reference methods receive a separate matrix only when requested.
- `RUN_CIBERSORT_COMPARISON` (default TRUE) emits concordance tables and PDFs
  only for human data after checking the actual native and IOBR reference
  matrices and QN settings match, and the input maximum is at least 50 to avoid
  IOBR's implicit low-value anti-log transform. Unavailable reference evidence disables the
  comparison; matching method names alone is insufficient.
- The shared HTML report (`reports/analysis_report.qmd`) is oriented to the
  General (DESeq2) layout, so `GENERATE_HTML_REPORT` defaults to `FALSE`; the
  wiring is present if you add a TME-specific report template.
- `visualize_results.R` reloads `4-TME/tme_results.Rdata` and re-plots TME
  barplots/boxplots/heatmaps, the ESTIMATE boxplot, and the ssGSEA heatmap
  without re-running any deconvolution.
