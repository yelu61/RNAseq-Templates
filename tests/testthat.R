library(testthat)

# Locate the helper library relative to the repository root.
# When running tests interactively, fall back to the current working directory.
repo_root <- tryCatch(
  rprojroot::find_root(rprojroot::is_git_root),
  error = function(e) getwd()
)
setwd(repo_root)
lib_dir <- file.path(repo_root, "RNAseq_lib")

source(file.path(lib_dir, "plot_utils.R"))
source(file.path(lib_dir, "io_utils.R"))
source(file.path(lib_dir, "deg_utils.R"))
source(file.path(lib_dir, "enrichment_utils.R"))
source(file.path(lib_dir, "design_utils.R"))
source(file.path(lib_dir, "tme_utils.R"))
source(file.path(lib_dir, "data_utils.R"))
source(file.path(lib_dir, "limma_voom_utils.R"))
source(file.path(lib_dir, "batch_utils.R"))
source(file.path(lib_dir, "survival_utils.R"))
source(file.path(lib_dir, "timecourse_utils.R"))
source(file.path(lib_dir, "report_utils.R"))
source(file.path(lib_dir, "geo_utils.R"))
source(file.path(lib_dir, "tcga_utils.R"))

# Run all tests under tests/testthat/
test_dir(file.path(repo_root, "tests", "testthat"))
