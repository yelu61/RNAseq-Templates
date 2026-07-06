# RNAseq-Templates

Notebook-first bulk RNA-seq analysis templates with a lightweight shared R helper library.

This repository is designed for practical project work: copy a notebook into a project, edit the parameter cell, and keep reusable validation, DEG thresholding, enrichment, and publication-style plotting code in `RNAseq_lib/`.

## Contents

```text
notebooks/
  RNAseq_General.ipynb
  RNAseq_TME_Deconvolution_Template.ipynb
  RNAseq_WGCNA_Template.ipynb
  RNAseq_TCGA_GEO_Template.ipynb
  RNAseq_limma_voom_Template.ipynb
  RNAseq_TimeCourse_Template.ipynb
RNAseq_lib/
  io_utils.R
  deg_utils.R
  enrichment_utils.R
  plot_utils.R
  tme_utils.R
  tcga_utils.R
  limma_voom_utils.R
  timecourse_utils.R
  survival_utils.R
examples/
  sample_config_template.R
  demo_data/                    # demo count table and metadata
  demo_RNAseq_General/          # pre-configured notebook ready to run
GETTING_STARTED.md              # zero-programming-experience guide
install_dependencies.R
```

`RNAseq_TCGA_GEO_Template.ipynb` is now included. It supports both TCGAbiolinks-based TCGA download and local GEO/TCGA expression matrices, with Tumor vs Normal DEG, single-gene expression, survival KM curves, and ORA/GSEA. Public data mining can later be split into a separate `PublicData-Templates` repository if it grows further.

## Quick Start (No Programming Experience Needed)

1. Install **R** and **RStudio** ([instructions in GETTING_STARTED.md](GETTING_STARTED.md)).
2. Install all analysis dependencies:
   ```r
   Rscript install_dependencies.R
   ```
3. Open the pre-configured demo notebook:
   ```
   examples/demo_RNAseq_General/RNAseq_General.ipynb
   ```
4. Click **Cell → Run All** and wait 2–5 minutes.
5. Check the generated `1-DEG/`, `2-GSEA/`, and `3-Visualization/` folders.

For a step-by-step guide, see [GETTING_STARTED.md](GETTING_STARTED.md).

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

`RNAseq_General.ipynb` also exports reusable intermediate files for topic notebooks:

```text
1-DEG/vsd_matrix.csv
1-DEG/colData.csv
```

The WGCNA and TME templates default to these exports, so the usual workflow is: run the general notebook first, then copy/run a topic notebook in the same project folder.

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
- key-gene and custom gene-set heatmaps
- GO/KEGG ORA with expressed-gene background, including dotplot, barplot, UP/DOWN bidirectional barplot, and optional cnetplot/emapplot
- GSEA GO/KEGG from full ranked gene lists, including dotplot, journal-style NES barplot, ridgeplot, and top activated/suppressed running enrichment curves
- GSVA for custom gene sets with heatmap, combined boxplot, and per-signature violin/box/jitter PDFs
- single-gene expression plots with mean bar, SEM, sample points, and layered pairwise P values
- optional DoRothEA/VIPER TF activity analysis
- DEG set overlap visualization (UpSet plot + Jaccard heatmap) across comparisons and thresholds

## Multi-threshold DEG + ORA

The general notebook supports threshold grids such as:

```r
DEG_PVALUE_COLUMN <- "padj" # "padj" recommended; use "pvalue" only for exploratory screening
DEG_LFC_COLUMN <- "log2FoldChange_shrunken" # DEG statistics/volcano/heatmap/ORA
GSEA_RANK_COLUMN <- "stat" # recommended for preranked GSEA; independent of DEG_LFC_COLUMN
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
- DEG result tables using the selected `DEG_LFC_COLUMN`
- DEG count summaries using the selected `DEG_LFC_COLUMN`
- raw-vs-shrunken DEG count comparison in `1-DEG/DEG_lfc_strategy_summary.csv`
- GO ORA
- KEGG ORA
- ORA summary plots
- threshold-specific ORA PDFs under `3-Visualization/<threshold>/`

GSEA is intentionally not repeated by DEG threshold because it uses the full ranked gene list. The default rank is the DESeq2 Wald statistic (`stat`), while DEG calls, volcano plots, DEG heatmaps, ORA, and compareCluster follow the selected `DEG_LFC_COLUMN`.

## Topic Templates

- `RNAseq_TME_Deconvolution_Template.ipynb`: ssGSEA immune/stromal signatures, native ESTIMATE, **IOBR multi-algorithm deconvolution (CIBERSORT/EPIC/xCell/ESTIMATE)**, optional native CIBERSORT, group comparisons, and heatmaps.
- `RNAseq_WGCNA_Template.ipynb`: expression filtering, sample QC, soft-threshold selection, module detection, module-trait correlation, hub gene export.
- `RNAseq_TCGA_GEO_Template.ipynb`: TCGA data download via TCGAbiolinks or local GEO matrices, counts/TPM preparation, Tumor vs Normal DESeq2 DEG, single-gene expression, **KM survival (median/quartile), univariate and multivariate Cox regression**, ORA/GSEA.
- `RNAseq_limma_voom_Template.ipynb`: `edgeR` + `limma-voom` alternative DEG workflow with optional batch correction and multi-contrast analysis.
- `RNAseq_TimeCourse_Template.ipynb`: Mfuzz soft clustering of time-series expression, trend plots, cluster heatmaps, and per-cluster ORA.

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
- Run `Rscript install_dependencies.R` once to install all required CRAN, Bioconductor, and GitHub packages.
