# RNAseq-Templates

Notebook-first bulk RNA-seq analysis templates with a lightweight shared R helper library.

This repository is designed for practical project work: copy a notebook into a project, edit the parameter cell, and keep reusable validation, DEG thresholding, enrichment, and publication-style plotting code in `RNAseq_lib/`.

## Contents

```text
notebooks/
  RNAseq_General.ipynb
  RNAseq_TME_Deconvolution_Template.ipynb
  RNAseq_WGCNA_Template.ipynb
RNAseq_lib/
  io_utils.R
  deg_utils.R
  enrichment_utils.R
  plot_utils.R
examples/
  sample_config_template.R
```

`RNAseq_TCGA_GEO_Template.ipynb` is intentionally not included in this first repository version. Public data mining can later become an optional module or a separate `PublicData-Templates` repository.

## Recommended Usage

1. Copy a notebook from `notebooks/` into a concrete project folder.
2. Keep access to `RNAseq_lib/` by either:
   - copying the `RNAseq_lib/` directory next to the project notebook, or
   - adjusting `LIB_DIR` in the notebook to point back to this repository.
3. Edit only the parameter cell first: species, input paths, sample names, groups, comparisons, thresholds, and optional gene sets.
4. Run all cells in order.

The general notebook writes analysis outputs to:

```text
1-DEG/
2-GSEA/
3-Visualization/
Analysis_summary.txt
sessionInfo.txt
```

All major visualizations are saved as editable PDF files under `3-Visualization/`. Threshold-specific ORA figures are written to `3-Visualization/<threshold>/`.

## General Bulk RNA-seq Workflow

`RNAseq_General.ipynb` includes:

- raw count loading and validation
- protein-coding filtering when annotation is available
- duplicate gene symbol handling
- low-count filtering based on the smallest group size
- DESeq2 differential expression with `lfcShrink(type = "ashr")`
- full all-gene DESeq2 result tables preserving raw and shrunken log2FC
- multi-threshold DEG outputs with `THRESHOLD_GRID`
- PCA and sample distance QC
- volcano plots and annotated DEG heatmaps
- GO/KEGG ORA with expressed-gene background, including dotplot, barplot, UP/DOWN bidirectional barplot, and optional cnetplot/emapplot
- GSEA GO/KEGG from full ranked gene lists, including dotplot, NES barplot, ridgeplot, and running enrichment curves
- GSVA for custom gene sets with heatmap and pairwise boxplot statistics
- single-gene expression plots with pairwise adjusted P values and mean-difference effect labels
- optional DoRothEA/VIPER TF activity analysis

## Multi-threshold DEG + ORA

The general notebook supports threshold grids such as:

```r
DEG_PVALUE_COLUMN <- "padj" # "padj" recommended; use "pvalue" only for exploratory screening
DEG_LFC_COLUMN <- "log2FoldChange_shrunken" # DEG thresholding and volcano plots
GSEA_RANK_COLUMN <- "stat" # recommended for preranked GSEA
THRESHOLD_GRID <- data.frame(
  name     = c("strict", "standard", "loose"),
  p_cutoff = c(0.01, 0.05, 0.10),
  log2fc   = c(1.5, 1.0, 0.5),
  stringsAsFactors = FALSE
)
DEFAULT_THRESHOLD <- "standard"
```

Multi-threshold output applies to:

- full all-gene result tables under `1-DEG/all_genes/` are saved once per comparison and are not threshold-filtered
- DEG result tables
- DEG count summaries
- GO ORA
- KEGG ORA
- ORA summary plots
- threshold-specific ORA PDFs under `3-Visualization/<threshold>/`

GSEA is intentionally not repeated by DEG threshold because it uses the full ranked gene list. The default rank is the DESeq2 Wald statistic (`stat`), while DEG calls and volcano plots default to shrunken log2FC.

## Topic Templates

- `RNAseq_TME_Deconvolution_Template.ipynb`: ssGSEA immune/stromal signatures, optional ESTIMATE, optional CIBERSORT skeleton, group comparisons, heatmaps.
- `RNAseq_WGCNA_Template.ipynb`: expression filtering, sample QC, soft-threshold selection, module detection, module-trait correlation, hub gene export.

## GitHub Sync

After creating an empty GitHub repository, connect and push:

```bash
git remote add origin <GITHUB_REPO_URL>
git branch -M main
git push -u origin main
```

## Notes

- This is not a formal R package. The helper library is deliberately lightweight so notebooks remain readable.
- For publication figures, edit the saved PDFs from `3-Visualization/`.
- For project-specific models with batch, paired design, or covariates, edit `DESIGN_FORMULA`, for example `~ batch + condition`.
