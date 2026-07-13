# RNAseq_lib

Lightweight R helper scripts shared by the RNA-seq notebook templates.

Source order used by the notebooks:

```r
source("RNAseq_lib/plot_utils.R")
source("RNAseq_lib/io_utils.R")
source("RNAseq_lib/deg_utils.R")
source("RNAseq_lib/enrichment_utils.R")
source("RNAseq_lib/batch_utils.R")
source("RNAseq_lib/design_utils.R")
source("RNAseq_lib/data_utils.R")
source("RNAseq_lib/tme_utils.R")
source("RNAseq_lib/tcga_utils.R")
source("RNAseq_lib/limma_voom_utils.R")
source("RNAseq_lib/timecourse_utils.R")
source("RNAseq_lib/survival_utils.R")
```

The scripts intentionally stay small and dependency-light instead of becoming a formal R package. This keeps each notebook readable while avoiding copy-pasted validation, DEG thresholding, enrichment, and PDF plotting code. Run `Rscript install_dependencies.R` from the repository root to install all dependencies.

## New helpers (publication-grade enrichment visualization)

- `plot_utils.R`
  - `default_enrichment_themes()` – cancer/immunology GO-BP theme dictionary.
  - `match_enrichment_themes()` – assign terms to themes by regex.
  - `prepare_theme_dotplot_df()` – select top terms per theme.
- `plot_gsea_term_figure_pdf()` – publication-quality single-term GSEA running-enrichment figure (NES-aware colors, clean title/subtitle with NES/p/adj p, panel styling).
  - `plot_gsea_term_figures_from_df()` – batch wrapper to generate single-term figures for a data frame of selected terms.
- `enrichment_utils.R`
  - `enrich_result_to_df()` – safely convert clusterProfiler results to data frames.
  - `build_multi_comparison_enrich_df()` – combine enrichment results across comparisons.
- `design_utils.R`
  - `validate_batch_design()` – warn when `BATCH_VECTOR` is provided but the DESeq2 design formula does not adjust for batch.
