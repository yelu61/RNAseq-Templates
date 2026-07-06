# RNAseq_lib

Lightweight R helper scripts shared by the RNA-seq notebook templates.

Source order used by the notebooks:

```r
source("RNAseq_lib/plot_utils.R")
source("RNAseq_lib/io_utils.R")
source("RNAseq_lib/deg_utils.R")
source("RNAseq_lib/enrichment_utils.R")
source("RNAseq_lib/tme_utils.R")
source("RNAseq_lib/tcga_utils.R")
source("RNAseq_lib/limma_voom_utils.R")
source("RNAseq_lib/timecourse_utils.R")
source("RNAseq_lib/survival_utils.R")
```

The scripts intentionally stay small and dependency-light instead of becoming a formal R package. This keeps each notebook readable while avoiding copy-pasted validation, DEG thresholding, enrichment, and PDF plotting code. Run `Rscript install_dependencies.R` from the repository root to install all dependencies.
