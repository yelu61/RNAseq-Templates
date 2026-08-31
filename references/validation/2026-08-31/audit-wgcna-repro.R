# Minimal reproduction of the v0.11.0 WGCNA runner's hub-export defect.
# No repository files or full analysis outputs are written.
suppressPackageStartupMessages(library(WGCNA))
cat(R.version.string, "\nWGCNA", as.character(packageVersion("WGCNA")), "\n")
set.seed(7)
base <- rnorm(30)
datExpr <- sapply(1:60, function(i) base + rnorm(30, sd = 0.15))
colnames(datExpr) <- paste0("g", seq_len(ncol(datExpr)))
rownames(datExpr) <- paste0("s", seq_len(nrow(datExpr)))

# Same network options and subsequent operations as the production runner;
# power is fixed here to avoid unrelated soft-threshold selection.
net <- blockwiseModules(
  datExpr, power = 6, networkType = "signed", TOMType = "signed",
  minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE,
  verbose = 0
)
moduleColors <- labels2colors(net$colors)
MEs <- orderMEs(net$MEs)
mod <- unique(moduleColors[moduleColors != "grey"])[1]
mod_genes <- colnames(datExpr)[moduleColors == mod]
ME <- MEs[[paste0("ME", mod)]]
kME <- cor(datExpr[, mod_genes, drop = FALSE], ME, use = "p")
hub <- data.frame(gene = mod_genes, module = mod, kME = as.numeric(kME))

cat("Eigengene columns:", paste(names(MEs), collapse = ","), "\n")
cat("Module color:", mod, "\n")
cat("Requested eigengene column:", paste0("ME", mod), "\n")
cat("Eigengene lookup returned NULL:", is.null(ME), "\n")
cat("Expected hub rows:", length(mod_genes), "\n")
cat("Actual hub rows:", nrow(hub), "\n")
cat("Duplicated gene rows:", sum(duplicated(hub$gene)), "\n")
cat("Computed correlation dimensions:", paste(dim(kME), collapse = " x "), "\n")
stopifnot(is.null(ME), nrow(hub) == length(mod_genes)^2,
          anyDuplicated(hub$gene) > 0)
cat("REPRODUCED: hub export contains gene-gene correlations, not kME.\n")
