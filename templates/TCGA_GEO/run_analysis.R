#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_TCGA_GEO — non-interactive analysis pipeline
# =============================================================================
# Runs the SAME workflow as notebooks/RNAseq_TCGA_GEO_Template.ipynb, driven by
# a config file instead of a parameter cell. Intended use:
#
#   1. Copy templates/TCGA_GEO/ (config.R + run_analysis.R) into your project.
#   2. Edit config.R (data source, paths, thresholds, survival genes).
#   3. Run:  Rscript run_analysis.R
#      or:   Rscript run_analysis.R path/to/other_config.R
#
# Data acquisition is 3-way (config.R picks exactly one):
#   DOWNLOAD_FROM_GDC = TRUE   -> TCGAbiolinks GDC download (needs network)
#   else DOWNLOAD_FROM_GEO = TRUE -> GEOquery SeriesMatrix (needs network)
#   else                          -> LOCAL_* files (fully offline; used by the demo)
#
# LIB_DIR resolution: if ./RNAseq_lib exists next to this script it is used;
# otherwise the repository root is located via rprojroot and its RNAseq_lib is
# used. Set the environment variable RNASEQ_LIB_DIR to override explicitly.
# =============================================================================

options(stringsAsFactors = FALSE)

# ---- Locate this script and the config file ----------------------------------
invocation_dir <- normalizePath(getwd(), mustWork = TRUE)
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", cmd_args, value = TRUE))
file_arg <- gsub("~\\+~", " ", file_arg)
script_dir <- if (length(file_arg) > 0 && nzchar(file_arg)) dirname(normalizePath(file_arg[1])) else invocation_dir

user_args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(user_args) >= 1) normalizePath(path.expand(user_args[1]), mustWork = FALSE) else file.path(script_dir, "config.R")
if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path,
       "\nProvide one as: Rscript run_analysis.R path/to/config.R")
}
config_path <- normalizePath(config_path, mustWork = TRUE)

cat("========================================\n")
cat("RNAseq_TCGA_GEO — run_analysis.R\n")
cat("Working dir :", getwd(), "\n")
cat("Config file :", config_path, "\n")
cat("========================================\n\n")

# Source the config into the global environment so all parameters are available.
source(config_path, local = globalenv())

# ---- Resolve RNAseq_lib -------------------------------------------------------
lib_dir <- Sys.getenv("RNASEQ_LIB_DIR", unset = NA_character_)
if (is.na(lib_dir) || !dir.exists(lib_dir)) {
  if (dir.exists(file.path(invocation_dir, "RNAseq_lib"))) {
    lib_dir <- file.path(invocation_dir, "RNAseq_lib")
  } else if (dir.exists(file.path(script_dir, "RNAseq_lib"))) {
    lib_dir <- file.path(script_dir, "RNAseq_lib")
  } else {
    repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root, path = script_dir), error = function(e) NA_character_)
    if (!is.na(repo_root) && dir.exists(file.path(repo_root, "RNAseq_lib"))) {
      lib_dir <- file.path(repo_root, "RNAseq_lib")
    } else if (dir.exists(file.path(dirname(script_dir), "RNAseq_lib"))) {
      lib_dir <- file.path(dirname(script_dir), "RNAseq_lib")
    } else {
      stop("Could not locate RNAseq_lib. Set RNASEQ_LIB_DIR or run from a project ",
           "that has RNAseq_lib/ alongside run_analysis.R.")
    }
  }
}
cat("RNAseq_lib  :", normalizePath(lib_dir), "\n\n")

# ---- Load libraries -----------------------------------------------------------
species_org_db <- if (SPECIES == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
if (!requireNamespace(species_org_db, quietly = TRUE)) {
  stop("Annotation package '", species_org_db, "' is required. ",
       "Install with BiocManager::install('", species_org_db, "')")
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(SummarizedExperiment)
  library(DESeq2)
  library(clusterProfiler)
  library(species_org_db, character.only = TRUE)
  library(survival)
  library(survminer)
})
species_org   <- get(species_org_db)
organism_code <- if (SPECIES == "human") "hsa" else "mmu"

for (f in c("plot_utils.R", "deg_utils.R", "enrichment_utils.R", "tcga_utils.R",
            "survival_utils.R", "geo_utils.R", "data_utils.R", "report_utils.R", "run_utils.R")) {
  source(file.path(lib_dir, f))
}
theme_set(theme_publication())
initialize_run_lifecycle(globalenv())

# ---- Output directories + config snapshot ------------------------------------
for (d in c("0-Config", "1-DEG", "2-GSEA", "3-Visualization", "6-Survival")) {
  dir.create(file.path(OUTDIR, d), showWarnings = FALSE, recursive = TRUE)
}
deg_dir  <- file.path(OUTDIR, "1-DEG")
gsea_dir <- file.path(OUTDIR, "2-GSEA")
viz_dir  <- file.path(OUTDIR, "3-Visualization")
surv_dir <- file.path(OUTDIR, "6-Survival")

config_objects <- c(
  "SPECIES", "DOWNLOAD_FROM_GDC", "DOWNLOAD_FROM_GEO", "GEO_ACCESSION",
  "TCGA_PROJECT", "TCGA_DATA_CATEGORY", "TCGA_DATA_TYPE", "TCGA_WORKFLOW",
  "GDC_COUNTS_ASSAY", "GDC_TPM_ASSAY", "GENE_ID_MAP_FILE",
  "LOCAL_COUNTS_FILE", "LOCAL_TPM_FILE", "LOCAL_CLINICAL_FILE", "LOCAL_GENE_COLUMN",
  "MIN_COUNT_PER_SAMPLE_FRAC", "MIN_COUNT", "TUMOR_NORMAL_DESIGN",
  "DEG_LFC_CUTOFF", "DEG_PADJ_CUTOFF",
  "GENES_FOR_SURVIVAL", "CLINICAL_VARS_FOR_KM", "TIME_UNIT",
  "OUTDIR", "GENERATE_HTML_REPORT", "REPORT_TITLE", "RUN_ROLE",
  "PARENT_RUN_ID", "RUN_CHANGE_NOTE", "RUN_RETENTION"
)
config_lines <- c(
  "# RNAseq_TCGA_GEO analysis configuration snapshot",
  paste0("# Saved: ", Sys.time()),
  unlist(lapply(config_objects, function(obj) {
    if (exists(obj, inherits = TRUE)) c(paste0("\n", obj, " <- "), capture.output(dput(get(obj, inherits = TRUE)))) else NULL
  }))
)
writeLines(config_lines, file.path(OUTDIR, "0-Config", "analysis_config_used.R"))

# =============================================================================
# 3. Data acquisition (GDC download / GEO download / local files)
# =============================================================================
if (isTRUE(DOWNLOAD_FROM_GDC)) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Package 'TCGAbiolinks' is required for DOWNLOAD_FROM_GDC = TRUE. ",
         "Install with BiocManager::install('TCGAbiolinks'). This mode also needs network access.")
  }
  message("Downloading from GDC via TCGAbiolinks (network required): ", TCGA_PROJECT)
  dir.create(file.path(OUTDIR, "GDCdata"), showWarnings = FALSE, recursive = TRUE)
  query <- build_tcga_query(
    project = TCGA_PROJECT,
    data.category = TCGA_DATA_CATEGORY,
    data.type = TCGA_DATA_TYPE,
    workflow.type = TCGA_WORKFLOW
  )
  TCGAbiolinks::GDCdownload(query, directory = file.path(OUTDIR, "GDCdata"))
  se <- TCGAbiolinks::GDCprepare(query, save = TRUE,
                                  save.filename = file.path(OUTDIR, paste0(TCGA_PROJECT, "_mRNA.Rdata")),
                                  directory = file.path(OUTDIR, "GDCdata"))
  assays_list <- extract_tcga_assays(se)
  cat("Available assays:", paste(names(assays_list), collapse = ", "), "\n")

  # Discover count and TPM assays by name instead of hard-coded index.
  assay_names <- names(assays_list)
  counts_assay <- if (!is.null(GDC_COUNTS_ASSAY)) GDC_COUNTS_ASSAY else {
    candidates <- assay_names[grepl("count|unstrand", assay_names, ignore.case = TRUE)]
    if (length(candidates) == 0) stop("No count assay found in SummarizedExperiment. Available: ", paste(assay_names, collapse = ", "))
    candidates[1]
  }
  tpm_assay <- if (!is.null(GDC_TPM_ASSAY)) GDC_TPM_ASSAY else {
    candidates <- assay_names[grepl("tpm", assay_names, ignore.case = TRUE)]
    if (length(candidates) == 0) stop("No TPM assay found in SummarizedExperiment. Available: ", paste(assay_names, collapse = ", "))
    candidates[1]
  }
  if (!counts_assay %in% assay_names) stop("Requested count assay not found: ", counts_assay)
  if (!tpm_assay %in% assay_names) stop("Requested TPM assay not found: ", tpm_assay)
  cat("Using count assay:", counts_assay, "; TPM assay:", tpm_assay, "\n")

  # Map ENSEMBL IDs to gene symbols and de-duplicate.
  if (is.null(GENE_ID_MAP_FILE)) {
    id_map <- build_id_map_from_se(se, id_col = "gene_id", symbol_col = "gene_name", type_col = "gene_type")
    cat("Derived gene ID map from SummarizedExperiment rowData:", nrow(id_map), "rows\n")
  } else {
    id_map <- read.delim(GENE_ID_MAP_FILE, stringsAsFactors = FALSE)
  }
  counts_raw   <- symbolize_and_dedup(assays_list[[counts_assay]], id_map, id_col = "gene_id", symbol_col = "gene_name")
  tpm_raw      <- symbolize_and_dedup(assays_list[[tpm_assay]], id_map, id_col = "gene_id", symbol_col = "gene_name")
  clinical_raw <- extract_tcga_clinical(se)
} else if (isTRUE(DOWNLOAD_FROM_GEO)) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package 'GEOquery' is required for DOWNLOAD_FROM_GEO = TRUE. ",
         "Install with BiocManager::install('GEOquery'). This mode also needs network access.")
  }
  message("Downloading GEO SeriesMatrix (network required): ", GEO_ACCESSION)
  gse <- download_geo_series_matrix(GEO_ACCESSION, destdir = file.path(OUTDIR, "0-Data"))
  geo <- parse_geo_series_matrix(gse)
  counts_raw   <- prepare_geo_counts(geo$expr, geo$feature)
  clinical_raw <- geo$pdata
  tpm_raw      <- NULL
  cat("GEO assay:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
} else {
  message("Local-file mode (offline).")
  counts_raw <- read_expression_matrix(LOCAL_COUNTS_FILE, gene_column = LOCAL_GENE_COLUMN)
  tpm_raw <- if (!is.null(LOCAL_TPM_FILE) && file.exists(LOCAL_TPM_FILE)) {
    read_expression_matrix(LOCAL_TPM_FILE, gene_column = LOCAL_GENE_COLUMN)
  } else NULL
  if (!is.null(tpm_raw)) {
    validate_expression_contract(tpm_raw, expected = "tpm")
    validate_samples_match(colnames(tpm_raw), colnames(counts_raw), context = "TPM vs counts")
  }
  clinical_raw <- read_metadata(LOCAL_CLINICAL_FILE, sample_column = "barcode", required_columns = "barcode")
}

# DESeq2 is valid only for raw integer counts. GEO SeriesMatrix values are often
# normalized microarray/log expression and will fail here; use limma for those data.
tryCatch(validate_count_matrix(counts_raw), error = function(e) {
  stop("Differential-expression input is not raw integer counts: ", conditionMessage(e),
       " For normalized GEO SeriesMatrix data, use a limma workflow or obtain raw RNA-seq counts.")
})
if (!is.null(tpm_raw)) validate_expression_contract(tpm_raw, expected = "tpm")
cat("Counts:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
if (!is.null(tpm_raw)) cat("TPM:", nrow(tpm_raw), "genes x", ncol(tpm_raw), "samples\n")

# =============================================================================
# 4. Clinical cleanup and sample filtering
# =============================================================================
colnames(clinical_raw) <- make.names(colnames(clinical_raw), unique = TRUE)

# Infer tumor/normal grouping from the TCGA barcode when one is available.
if ("barcode" %in% colnames(clinical_raw)) {
  clinical_raw$tissue_type <- infer_tcga_tumor_normal(clinical_raw$barcode)
}

validate_samples_match(colnames(counts_raw), clinical_raw$barcode, context = "counts vs clinical")

# Keep only samples present in counts.
common_samples <- intersect(colnames(counts_raw), clinical_raw$barcode)
counts_raw <- counts_raw[, common_samples, drop = FALSE]
if (!is.null(tpm_raw)) tpm_raw <- tpm_raw[, common_samples, drop = FALSE]
clinical <- clinical_raw[match(common_samples, clinical_raw$barcode), ]

has_condition <- "tissue_type" %in% colnames(clinical)
if (has_condition) {
  clinical$condition <- factor(clinical$tissue_type, levels = c("Normal", "Tumor"))
}
if (isTRUE(TUMOR_NORMAL_DESIGN)) {
  if (!has_condition) {
    stop("TUMOR_NORMAL_DESIGN = TRUE requires a tissue_type derivable from TCGA barcodes.")
  }
  if (anyNA(clinical$condition) || !all(c("Normal", "Tumor") %in% clinical$condition)) {
    stop("Tumor/Normal design requires every retained sample to have a valid tissue type and both groups to be present.")
  }
}

write.csv(clinical, file.path(surv_dir, "clinical_clean.csv"), row.names = FALSE)
if (has_condition) print(table(clinical$condition, useNA = "ifany"))

# =============================================================================
# 5. DESeq2 Tumor vs Normal DEG
# =============================================================================
res_df <- NULL
sig_genes <- character(0)
dds <- NULL
if (isTRUE(TUMOR_NORMAL_DESIGN)) {
  keep_genes <- rowSums(counts_raw > MIN_COUNT) >= ceiling(MIN_COUNT_PER_SAMPLE_FRAC * ncol(counts_raw))
  counts_filt <- as.matrix(counts_raw[keep_genes, ])
  mode(counts_filt) <- "numeric"

  colData <- data.frame(row.names = colnames(counts_filt), condition = clinical$condition)
  dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = colData, design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("condition", "Tumor", "Normal"))
  res_df <- as.data.frame(res) |>
    rownames_to_column("gene_name") |>
    arrange(padj)
  write.csv(res_df, file.path(deg_dir, "DESeq2_Tumor_vs_Normal.csv"), row.names = FALSE)

  sig_genes <- res_df |>
    filter(!is.na(padj), padj < DEG_PADJ_CUTOFF, abs(log2FoldChange) > DEG_LFC_CUTOFF) |>
    pull(gene_name)
  cat("Significant DEGs:", length(sig_genes), "\n")
}

# =============================================================================
# 6. DEG Visualization
# =============================================================================
if (isTRUE(TUMOR_NORMAL_DESIGN) && !is.null(res_df)) {
  plot_volcano_pdf(
    res_df, comp_name = "Tumor_vs_Normal",
    pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
    filename = file.path(viz_dir, "Volcano_Tumor_vs_Normal.pdf"),
    pvalue_column = "padj", lfc_column = "log2FoldChange"
  )
}

# =============================================================================
# 7. Single-gene expression + Survival (KM by median/quantile, clinical KM, Cox)
# =============================================================================
# Log-scale expression for sample comparison: log2(TPM + 1) when TPM is present,
# otherwise VST derived from the validated raw counts (never log2(raw counts + 1),
# which leaves library size as a confounder).
if (!is.null(tpm_raw)) {
  expr_log <- log2(as.matrix(tpm_raw) + 1)
} else {
  dds_for_expression <- if (!is.null(dds)) dds else {
    DESeqDataSetFromMatrix(as.matrix(counts_raw),
      colData = data.frame(row.names = colnames(counts_raw)), design = ~ 1)
  }
  expr_log <- assay(vst(dds_for_expression, blind = !isTRUE(TUMOR_NORMAL_DESIGN)))
}

# Survival requires the TCGA clinical columns; skip gracefully when they are
# absent (e.g. a GEO SeriesMatrix phenotype table without follow-up fields).
surv_required <- c("barcode", "vital_status", "days_to_death", "days_to_last_follow_up")
surv_df <- NULL
if (all(surv_required %in% colnames(clinical))) {
  surv_df <- prepare_tcga_survival(clinical)
  cat("Survival samples:", nrow(surv_df), " | Events:", sum(surv_df$status), "\n")
} else {
  message("Survival skipped: clinical table lacks columns ",
          paste(setdiff(surv_required, colnames(clinical)), collapse = ", "))
}

# Per-gene expression boxplot + a SINGLE median-split KM + a quartile-split KM.
# (The notebook computed the median-split KM twice, in sections 7 and 7.5; this
# runner keeps one median-split path and adds the quartile split alongside it.)
for (gene in GENES_FOR_SURVIVAL) {
  if (!gene %in% rownames(expr_log)) {
    warning("Gene not found in expression matrix: ", gene)
    next
  }
  if (has_condition) {
    plot_tcga_gene_boxplot_pdf(
      expr_log, gene, clinical$condition,
      filename = file.path(viz_dir, paste0("Expression_boxplot_", gene, ".pdf")),
      title = paste(gene, "expression")
    )
  }
  if (!is.null(surv_df)) {
    surv_df[[gene]] <- as.numeric(expr_log[gene, surv_df$barcode])
    plot_km_by_median_pdf(
      surv_df, value_col = gene,
      filename = file.path(surv_dir, paste0("KM_", gene, ".pdf")),
      title = paste(gene, "survival"),
      time_unit = TIME_UNIT
    )
    surv_df[[paste0(gene, "_quartile")]] <- stratify_by_quantile(surv_df[[gene]], n_groups = 4)
    plot_km_by_group_pdf(
      surv_df, group_col = paste0(gene, "_quartile"),
      filename = file.path(surv_dir, paste0("KM_quartile_", gene, ".pdf")),
      title = paste(gene, "quartile survival"),
      time_unit = TIME_UNIT
    )
  }
}

# KM by clinical variables (categorical or median-split numeric).
if (!is.null(surv_df) && length(CLINICAL_VARS_FOR_KM) > 0) {
  run_clinical_km(
    surv_df,
    clinical_df = clinical,
    var_cols = CLINICAL_VARS_FOR_KM,
    outdir = surv_dir,
    time_unit = TIME_UNIT
  )
}

# Univariate + multivariate Cox for the genes of interest.
uni_cox <- NULL
multi_cox <- NULL
if (!is.null(surv_df)) {
  surv_vars <- intersect(GENES_FOR_SURVIVAL, colnames(surv_df))
  if (length(surv_vars) > 0) {
    uni_cox <- run_univariate_cox(surv_df, vars = surv_vars)
    if (!is.null(uni_cox)) {
      write.csv(uni_cox, file.path(surv_dir, "univariate_Cox.csv"), row.names = FALSE)
      plot_cox_forest_pdf(
        uni_cox,
        filename = file.path(surv_dir, "univariate_Cox_forest.pdf"),
        title = "Univariate Cox Regression"
      )
    }
  }
  if (length(surv_vars) >= 2) {
    multi_cox <- run_multivariate_cox(surv_df, vars = surv_vars)
    if (!is.null(multi_cox)) {
      write.csv(multi_cox, file.path(surv_dir, "multivariate_Cox.csv"), row.names = FALSE)
      plot_cox_forest_pdf(
        multi_cox,
        filename = file.path(surv_dir, "multivariate_Cox_forest.pdf"),
        title = "Multivariate Cox Regression"
      )
    }
  }
}

# =============================================================================
# 8. ORA and GSEA (+ theme maps + single-term GSEA figures)
# =============================================================================
ego <- NULL
ekegg <- NULL
gsea_go <- NULL
gsea_kegg <- NULL
if (isTRUE(TUMOR_NORMAL_DESIGN) && !is.null(res_df)) {
  # ORA
  deg_list <- genes_for_enrichment(
    res_df, pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
    pvalue_column = "padj", lfc_column = "log2FoldChange"
  )
  universe <- map_symbols_to_entrez(rownames(counts_raw), species_org)
  ego   <- run_go_ora(deg_list$sig, org_db = species_org, universe = universe$ENTREZID)
  ekegg <- run_kegg_ora(deg_list$sig, org_db = species_org, universe = universe$ENTREZID, organism = organism_code)
  if (!is.null(ego)) {
    write.csv(as.data.frame(ego), file.path(gsea_dir, "GO_ORA_Tumor_vs_Normal.csv"), row.names = FALSE)
    plot_enrich_suite_pdf(ego, file.path(viz_dir, "GO_ORA_Tumor_vs_Normal"), "GO ORA")
  }
  if (!is.null(ekegg)) {
    write.csv(as.data.frame(ekegg), file.path(gsea_dir, "KEGG_ORA_Tumor_vs_Normal.csv"), row.names = FALSE)
    plot_enrich_suite_pdf(ekegg, file.path(viz_dir, "KEGG_ORA_Tumor_vs_Normal"), "KEGG ORA")
  }

  # GSEA
  ranked        <- ranked_gene_list(res_df, rank_column = "stat")
  entrez_ranked <- make_entrez_ranked_list(ranked, species_org)
  gsea_go   <- run_go_gsea(entrez_ranked, org_db = species_org)
  gsea_kegg <- run_kegg_gsea(entrez_ranked, organism = organism_code)
  if (!is.null(gsea_go)) {
    write.csv(as.data.frame(gsea_go), file.path(gsea_dir, "GO_GSEA_Tumor_vs_Normal.csv"), row.names = FALSE)
    plot_gsea_suite_pdf(gsea_go, file.path(viz_dir, "GO_GSEA_Tumor_vs_Normal"), "GO GSEA")
  }
  if (!is.null(gsea_kegg)) {
    write.csv(as.data.frame(gsea_kegg), file.path(gsea_dir, "KEGG_GSEA_Tumor_vs_Normal.csv"), row.names = FALSE)
    plot_gsea_suite_pdf(gsea_kegg, file.path(viz_dir, "KEGG_GSEA_Tumor_vs_Normal"), "KEGG GSEA")
  }

  # ---- 8.1 Theme dot-heatmaps -----------------------------------------------
  theme_outdir <- file.path(viz_dir, "ThemeEnrichment")
  dir.create(theme_outdir, showWarnings = FALSE, recursive = TRUE)
  go_ora_map    <- list()
  kegg_ora_map  <- list()
  go_gsea_map   <- list()
  kegg_gsea_map <- list()
  if (!is.null(ego))       go_ora_map[["Tumor_vs_Normal"]]    <- ego
  if (!is.null(ekegg))     kegg_ora_map[["Tumor_vs_Normal"]]  <- ekegg
  if (!is.null(gsea_go))   go_gsea_map[["Tumor_vs_Normal"]]   <- gsea_go
  if (!is.null(gsea_kegg)) kegg_gsea_map[["Tumor_vs_Normal"]] <- gsea_kegg

  theme_defs <- default_enrichment_themes()
  if (length(go_ora_map) > 0) {
    plot_theme_dotheatmap_from_results(go_ora_map, file.path(theme_outdir, "Theme_dotheatmap_GO_ORA.pdf"),
      title = "GO ORA Biological Themes", subtitle = "GO-BP ORA | Tumor vs Normal", theme_defs = theme_defs, ontology_filter = "BP")
  }
  if (length(kegg_ora_map) > 0) {
    plot_theme_dotheatmap_from_results(kegg_ora_map, file.path(theme_outdir, "Theme_dotheatmap_KEGG_ORA.pdf"),
      title = "KEGG ORA Pathway Themes", subtitle = "KEGG ORA | Tumor vs Normal", theme_defs = theme_defs, ontology_filter = NULL)
  }
  if (length(go_gsea_map) > 0) {
    plot_theme_dotheatmap_from_results(go_gsea_map, file.path(theme_outdir, "Theme_dotheatmap_GO_GSEA.pdf"),
      title = "GO GSEA Biological Themes", subtitle = "GO-BP GSEA | Tumor vs Normal", theme_defs = theme_defs, ontology_filter = "BP")
  }
  if (length(kegg_gsea_map) > 0) {
    plot_theme_dotheatmap_from_results(kegg_gsea_map, file.path(theme_outdir, "Theme_dotheatmap_KEGG_GSEA.pdf"),
      title = "KEGG GSEA Pathway Themes", subtitle = "KEGG GSEA | Tumor vs Normal", theme_defs = theme_defs, ontology_filter = NULL)
  }

  # ---- 8.2 Single-term GSEA figures ------------------------------------------
  single_term_outdir <- file.path(theme_outdir, "single_term_gsea")
  dir.create(single_term_outdir, showWarnings = FALSE, recursive = TRUE)
  if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
    ggo_df <- as.data.frame(gsea_go)
    ggo_df <- ggo_df[order(ggo_df$p.adjust, -abs(ggo_df$NES)), ]
    top_terms <- rbind(utils::head(ggo_df[ggo_df$NES > 0, ], 3),
                       utils::head(ggo_df[ggo_df$NES < 0, ], 3))
    plot_gsea_term_figures_from_df(gsea_go, top_terms,
      outdir = file.path(single_term_outdir, "GO_Tumor_vs_Normal"),
      contrast_label = "Tumor vs Normal", prefix = "gseaplot2_GO")
  }
  if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
    gkegg_df <- as.data.frame(gsea_kegg)
    gkegg_df <- gkegg_df[order(gkegg_df$p.adjust, -abs(gkegg_df$NES)), ]
    top_terms <- rbind(utils::head(gkegg_df[gkegg_df$NES > 0, ], 3),
                       utils::head(gkegg_df[gkegg_df$NES < 0, ], 3))
    plot_gsea_term_figures_from_df(gsea_kegg, top_terms,
      outdir = file.path(single_term_outdir, "KEGG_Tumor_vs_Normal"),
      contrast_label = "Tumor vs Normal", prefix = "gseaplot2_KEGG")
  }
}

# =============================================================================
# 9. Save results, summary, session info
# =============================================================================
# Targeted save so visualize_results.R can re-plot without re-running
# DESeq2 / GSEA / survival.
save(res_df, sig_genes, counts_raw, expr_log, clinical, surv_df,
     ego, ekegg, gsea_go, gsea_kegg, uni_cox, multi_cox,
     GENES_FOR_SURVIVAL, CLINICAL_VARS_FOR_KM, TIME_UNIT,
     SPECIES, TUMOR_NORMAL_DESIGN, DEG_PADJ_CUTOFF, DEG_LFC_CUTOFF,
     file = file.path(deg_dir, "TCGA_GEO_results.Rdata"))

summary_report <- paste0(
  "========================================\n",
  "TCGA/GEO Tumor-vs-Normal Analysis Summary\n",
  "========================================\n\n",
  "1. Data Overview\n",
  "   - Species: ", SPECIES, "\n",
  "   - Samples retained: ", ncol(counts_raw), "\n",
  "   - Genes (rows): ", nrow(counts_raw), "\n",
  "   - TPM available: ", ifelse(is.null(tpm_raw), "no", "yes"), "\n",
  "   - Tumor/Normal design: ", ifelse(isTRUE(TUMOR_NORMAL_DESIGN), "yes", "no"), "\n\n"
)
if (isTRUE(TUMOR_NORMAL_DESIGN) && !is.null(res_df)) {
  n_up   <- sum(res_df$padj < DEG_PADJ_CUTOFF & res_df$log2FoldChange >  DEG_LFC_CUTOFF, na.rm = TRUE)
  n_down <- sum(res_df$padj < DEG_PADJ_CUTOFF & res_df$log2FoldChange < -DEG_LFC_CUTOFF, na.rm = TRUE)
  summary_report <- paste0(summary_report,
    "2. Differential Expression (padj < ", DEG_PADJ_CUTOFF,
    " & |log2FC| > ", DEG_LFC_CUTOFF, ")\n",
    "   - Tumor_vs_Normal: ", n_up + n_down, " DEGs (Up: ", n_up, ", Down: ", n_down, ")\n\n")
}
if (!is.null(surv_df)) {
  summary_report <- paste0(summary_report,
    "3. Survival\n",
    "   - Survival samples: ", nrow(surv_df), " (events: ", sum(surv_df$status), ")\n",
    "   - Genes tested: ", paste(intersect(GENES_FOR_SURVIVAL, colnames(surv_df)), collapse = ", "), "\n\n")
}
summary_report <- paste0(summary_report,
  "4. Output Files\n",
  "   - ", file.path(OUTDIR, "0-Config", "analysis_config_used.R"), "\n",
  "   - ", deg_dir, " (DESeq2 table, results Rdata)\n",
  "   - ", gsea_dir, " (GO/KEGG ORA + GSEA tables)\n",
  "   - ", viz_dir, " (volcano, expression boxplots, ORA/GSEA suites, ThemeEnrichment/)\n",
  "   - ", surv_dir, " (clinical_clean.csv, KM + Cox forest figures, Cox/clinical-KM tables)\n",
  "\n========================================\n")
cat(summary_report)
writeLines(summary_report, file.path(OUTDIR, "Analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

# =============================================================================
# 10. HTML report (optional)
# =============================================================================
if (isTRUE(GENERATE_HTML_REPORT)) {
  has_quarto    <- requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))
  has_rmarkdown <- requireNamespace("rmarkdown", quietly = TRUE)
  if (!has_quarto && !has_rmarkdown) {
    message("Skipping HTML report: install the quarto CLI (+ R 'quarto' package) or the 'rmarkdown' package.")
  } else {
    report_path <- tryCatch(
      render_analysis_report(outdir = OUTDIR, report_file = file.path(OUTDIR, "RNAseq_report.html"),
                             params = list(title = REPORT_TITLE, author = Sys.info()[["user"]])),
      error = function(e) { message("HTML report failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(report_path)) cat("HTML report written to:", report_path, "\n")
  }
}

manifest_input <- if (isTRUE(DOWNLOAD_FROM_GDC) || isTRUE(DOWNLOAD_FROM_GEO)) {
  NA_character_
} else {
  LOCAL_COUNTS_FILE
}
write_template_run_manifest(
  run_dir = OUTDIR, config_path = config_path, input_file = manifest_input,
  config_objects = config_objects, lib_dir = lib_dir, runner_file = file_arg,
  envir = globalenv()
)

cat("\n========================================\n")
cat("RNAseq_TCGA_GEO analysis COMPLETE.\n")
cat("Outputs under:", normalizePath(OUTDIR), "\n")
cat("========================================\n")
