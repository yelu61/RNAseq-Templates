# Troubleshooting Guide

This guide lists common issues and how to fix them.

## Installation and Environment

### Q: `Rscript install_dependencies.R` reports missing packages or install failures

1. Make sure you are running the command from the repository root (where `install_dependencies.R` is located).
2. If a Bioconductor package fails, open RStudio Console and run:
   ```r
   if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install("包名")
   ```
3. If a CRAN package fails, run:
   ```r
   install.packages("包名")
   ```
4. For `IOBR` specifically, it installs from GitHub; ensure `remotes` or `devtools` is installed first.

### Q: Notebook reports `could not find function "..."`

- Make sure the notebook's `LIB_DIR` points to the correct `RNAseq_lib/` directory.
- Re-run the library-loading cell (cell-4 in most notebooks).
- Ensure all helpers are sourced, including new files such as `batch_utils.R` or `design_utils.R`.

## Input Files

### Q: Error "Input file not found"

- Check that `INPUT_FILE` path is correct. If the file is in `./0-Data/`, use `./0-Data/your_file.tsv`.
- In RStudio, use the Files pane to verify the file exists and the name is spelled correctly.

### Q: "Count columns contain non-numeric or missing values"

- Your count matrix must contain only integer raw counts in sample columns.
- Do not use TPM, FPKM, normalized values, or empty cells.
- Excel files sometimes convert gene IDs to dates; save as TSV/CSV if possible.

### Q: "SAMPLE_NAMES length must equal count columns"

- The number of sample names in `SAMPLE_NAMES` must exactly match the number of count columns.
- Check that `COUNT_COLS` or auto-detected columns are not including annotation columns such as `gene_biotype`.

## Experimental Design

### Q: "COMPARISONS contains group names not present in GROUP_LEVELS"

- The second and third values in each `COMPARISONS` entry must match names in `GROUP_LEVELS` exactly (case-sensitive).
- Example: if `GROUP_LEVELS = c("Control", "Treatment")`, a comparison should be `c("Treatment_vs_Control", "Treatment", "Control")`.

### Q: Paired design error "Some pair_id groups do not contain all conditions"

- Each `PAIR_ID` must appear once for every condition involved in the comparison.
- Check that `PAIR_ID` is the same length as `SAMPLE_NAMES` and correctly paired.

## DEG and Enrichment

### Q: 0 DEGs at all thresholds

- Verify group assignments are correct.
- Check that `GROUP_LEVELS` order matches your intended control/treatment direction.
- Try a looser threshold (e.g., `p_cutoff = 0.10`, `log2fc = 0.5`).
- For very small sample sizes, DESeq2 may report all NA p-values; check the warning messages.

### Q: "GO ORA skipped: only X mapped Entrez IDs"

- The input gene list is too small for ORA. This is informational; the pipeline continues with other thresholds/comparisons.
- If your species is not human/mouse, the Entrez mapping database (`org.Hs.eg.db` / `org.Mm.eg.db`) will not work; you need a custom mapping.

### Q: GSEA returns very few or no terms

- Make sure the rank column (`GSEA_RANK_COLUMN`) exists in the result table. The default `"stat"` is the DESeq2 Wald statistic.
- Check that enough genes map to Entrez IDs.

## Visualization

### Q: Some PDF plots are missing

- Missing plots usually mean there were no significant terms/genes to visualize. Check the console messages.
- Verify that `3-Visualization/` directory was created and you have write permission.

### Q: Batch PCA not generated

- `BATCH_VECTOR` must be non-NULL and have the same length as `SAMPLE_NAMES` / `ncol(vsd)`.
- Check that the batch-effect diagnostics cell ran after the PCA cell.

### Q: PDF plots fail with "failed to load cairo DLL" / X11 errors

- This is an environment issue (missing X11/cairo), not an analysis error. The analysis tables (CSV) are still produced.
- On macOS install XQuartz; on Linux install `libx11`/`cairo` dev packages. The demo smoke tests only require PDF figures when `capabilities("cairo")` is available.

### Q: HTML report is skipped

- Report generation needs the **quarto** CLI (https://quarto.org) plus the R `quarto` package, or the **rmarkdown** package. If neither is present the notebook prints a message and continues.
- The report only assembles already-saved CSV/PDF outputs, so re-running just the report cell after installing quarto/rmarkdown is enough.

## TCGA / GEO

### Q: Clinical KM / survival plots fail with "object of type 'symbol' is not subsettable"

- This was a bug in `plot_km_by_group_pdf()` under R ≥ 4.x where the survfit formula was captured as a symbol. Update `RNAseq_lib/survival_utils.R` to the latest version — the formula is now inlined into the survfit call.

## TME deconvolution

### Q: ESTIMATE fails with "找不到对象 'common_genes'" (object 'common_genes' not found)

- `requireNamespace("estimate")` does not load the package's lazy-data objects, but `filterCommonGenes()` / `estimateScore()` use them directly. The template now calls `utils::data("common_genes", package = "estimate")` and `utils::data("SI_geneset", package = "estimate")` first. Update the notebook / template if you have an older copy.

### Q: ssGSEA fails with "object 'immune_gene_sets' not found"

- Older TME templates referenced an undefined `immune_gene_sets` variable. It is now a built-in constant in `RNAseq_lib/tme_utils.R` (28 immune signatures, Charoentong et al. 2017). Source `tme_utils.R` before the ssGSEA cell.



### Q: `GENE_ID_MAP_FILE` not found

- In the updated notebook, `GENE_ID_MAP_FILE` defaults to `NULL` and the ID map is derived automatically from `SummarizedExperiment::rowData(se)`.
- Only set it to a file path if you want to use a custom mapping.

### Q: GEO download fails

- Make sure `GEOquery` is installed via `install_dependencies.R`.
- Some GEO series are very large and require a stable network connection; consider downloading manually and using local-file mode.

## Still stuck?

Prepare the following information when asking for help:

1. Which step failed (parameter cell, DESeq2, enrichment, visualization?).
2. The full error message (red text in RStudio).
3. Your parameter configuration (screenshot or copied text).
4. The `sessionInfo.txt` file generated by the notebook.
5. The `Analysis_summary.txt` if it was produced.
