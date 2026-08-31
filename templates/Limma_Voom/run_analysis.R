#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_Limma_Voom — non-interactive analysis pipeline
# =============================================================================
# Runs the SAME workflow as notebooks/RNAseq_limma_voom_Template.ipynb, driven by
# a config file instead of a parameter cell. Intended use:
#
#   1. Copy templates/Limma_Voom/ (config.R + run_analysis.R) into your project.
#   2. Edit config.R (paths, samples, groups, comparisons, thresholds).
#   3. Run:  Rscript run_analysis.R
#      or:   Rscript run_analysis.R path/to/other_config.R
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
cat("RNAseq_Limma_Voom — run_analysis.R\n")
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

source(file.path(lib_dir, "run_utils.R"))
assert_fresh_run_dir(OUTDIR)

# ---- Load libraries -----------------------------------------------------------
species_org_db <- if (SPECIES == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
if (!requireNamespace(species_org_db, quietly = TRUE)) {
  stop("Annotation package '", species_org_db, "' is required. ",
       "Install with BiocManager::install('", species_org_db, "')")
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(sva)
  library(clusterProfiler)
  library(species_org_db, character.only = TRUE)
})
species_org   <- get(species_org_db)
organism_code <- if (SPECIES == "human") "hsa" else "mmu"

for (f in c("plot_utils.R", "io_utils.R", "data_utils.R", "deg_utils.R",
            "enrichment_utils.R", "limma_voom_utils.R", "report_utils.R")) {
  source(file.path(lib_dir, f))
}
theme_set(theme_publication())
initialize_run_lifecycle(globalenv())

# ---- Output directories + config snapshot ------------------------------------
for (d in c("0-Config", "1-DEG", "2-GSEA", "3-Visualization")) {
  dir.create(file.path(OUTDIR, d), showWarnings = FALSE, recursive = TRUE)
}

config_objects <- c(
  "SPECIES", "INPUT_FILE", "INPUT_FORMAT", "GENE_NAME_COL", "BIOTYPE_COL",
  "BIOTYPE_FILTER", "COUNT_COLS", "SAMPLE_NAMES", "GROUPS", "GROUP_LEVELS",
  "BATCH_VECTOR", "COMPARISONS", "DEG_PADJ_CUTOFF", "DEG_LFC_CUTOFF", "MIN_COUNT",
  "MIN_SAMPLE_FRAC", "OUTDIR", "GENERATE_HTML_REPORT", "REPORT_TITLE",
  "RUN_ROLE", "PARENT_RUN_ID", "RUN_CHANGE_NOTE", "RUN_RETENTION"
)
config_lines <- c(
  "# RNAseq_Limma_Voom analysis configuration snapshot",
  paste0("# Saved: ", Sys.time()),
  unlist(lapply(config_objects, function(obj) {
    if (exists(obj, inherits = TRUE)) c(paste0("\n", obj, " <- "), capture.output(dput(get(obj, inherits = TRUE)))) else NULL
  }))
)
writeLines(config_lines, file.path(OUTDIR, "0-Config", "analysis_config_used.R"))

# =============================================================================
# 3. Load counts and build design matrix
# =============================================================================
rawcount <- read_count_table(INPUT_FILE, INPUT_FORMAT)
if (!is.null(BIOTYPE_COL) && BIOTYPE_COL %in% colnames(rawcount)) {
  rawcount <- rawcount[rawcount[[BIOTYPE_COL]] == BIOTYPE_FILTER, ]
}
count_col_names <- detect_count_columns(rawcount, GENE_NAME_COL, COUNT_COLS)
validate_sample_design(SAMPLE_NAMES, GROUPS, GROUP_LEVELS, COMPARISONS, count_col_names)

countData <- build_count_matrix(rawcount, GENE_NAME_COL, count_col_names, SAMPLE_NAMES)
validate_count_matrix(countData)

group <- factor(GROUPS, levels = GROUP_LEVELS)
design <- make_group_design(group, batch = BATCH_VECTOR)
if (!is.null(BATCH_VECTOR)) cat("Batch included as a covariate in the limma model.\n")

# =============================================================================
# 4. limma-voom pipeline
# =============================================================================
dge <- prepare_dge_for_voom(countData, group = group,
                            min_counts_per_sample = MIN_COUNT,
                            min_sample_frac = MIN_SAMPLE_FRAC)
v <- run_voom(dge, design = design,
              plot_file = file.path(OUTDIR, "3-Visualization", "voom_mean_variance_trend.pdf"))

res_list    <- run_limma_contrasts(v, design, COMPARISONS)
deg_summary <- write_limma_results(res_list, outdir = file.path(OUTDIR, "1-DEG"))
print(deg_summary)

# =============================================================================
# 5. Visualization
# =============================================================================
group_colors <- make_group_colors(GROUP_LEVELS)

for (comp_name in names(res_list)) {
  plot_volcano_pdf(
    res_list[[comp_name]], comp_name = comp_name,
    pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
    filename = file.path(OUTDIR, "3-Visualization", paste0("Volcano_", comp_name, ".pdf")),
    pvalue_column = "padj", lfc_column = "log2FoldChange"
  )
}
plot_deg_summary_pdf(deg_summary,
                     filename = file.path(OUTDIR, "3-Visualization", "DEG_summary.pdf"))

# Top DEG heatmap
all_sig_genes <- unique(unlist(lapply(res_list, function(res) {
  res |>
    dplyr::filter(!is.na(padj), padj < DEG_PADJ_CUTOFF, abs(log2FoldChange) > DEG_LFC_CUTOFF) |>
    dplyr::pull(gene_name)
})))
if (length(all_sig_genes) >= 2) {
  mat <- v$E[intersect(all_sig_genes, rownames(v$E)), ]
  plot_expression_heatmap_pdf(
    mat,
    filename = file.path(OUTDIR, "3-Visualization", "DEG_heatmap.pdf"),
    title = "limma-voom DEGs",
    group = GROUPS, group_levels = GROUP_LEVELS, group_colors = group_colors,
    width = mm_to_in(183), height = mm_to_in(247)
  )
}

# =============================================================================
# 6. ORA and GSEA
# =============================================================================
universe <- map_symbols_to_entrez(rownames(countData), species_org)$ENTREZID

go_ora_map   <- list()
kegg_ora_map <- list()
go_gsea_map  <- list()
kegg_gsea_map <- list()

for (comp_name in names(res_list)) {
  genes <- genes_for_enrichment(
    res_list[[comp_name]],
    pvalue_thresh = DEG_PADJ_CUTOFF, log2fc_thresh = DEG_LFC_CUTOFF,
    pvalue_column = "padj", lfc_column = "log2FoldChange"
  )
  ego   <- run_go_ora(genes$sig, org_db = species_org, universe = universe)
  ekegg <- run_kegg_ora(genes$sig, org_db = species_org, universe = universe, organism = organism_code)
  if (!is.null(ego)) {
    write.csv(as.data.frame(ego), file.path(OUTDIR, "2-GSEA", paste0("GO_ORA_", comp_name, ".csv")), row.names = FALSE)
    plot_enrich_suite_pdf(ego, file.path(OUTDIR, "3-Visualization", paste0("GO_ORA_", comp_name)), paste("GO ORA", comp_name))
    go_ora_map[[comp_name]] <- ego
  }
  if (!is.null(ekegg)) {
    write.csv(as.data.frame(ekegg), file.path(OUTDIR, "2-GSEA", paste0("KEGG_ORA_", comp_name, ".csv")), row.names = FALSE)
    plot_enrich_suite_pdf(ekegg, file.path(OUTDIR, "3-Visualization", paste0("KEGG_ORA_", comp_name)), paste("KEGG ORA", comp_name))
    kegg_ora_map[[comp_name]] <- ekegg
  }

  ranked        <- ranked_gene_list_limma(res_list[[comp_name]], rank_column = "t")
  entrez_ranked <- make_entrez_ranked_list(ranked, species_org)
  gsea_go   <- run_go_gsea(entrez_ranked, org_db = species_org)
  gsea_kegg <- run_kegg_gsea(entrez_ranked, organism = organism_code)
  if (!is.null(gsea_go)) {
    write_gsea_tables(gsea_go, file.path(OUTDIR, "2-GSEA", paste0("GO_GSEA_", comp_name, ".csv")))
    plot_gsea_suite_pdf(gsea_go, file.path(OUTDIR, "3-Visualization", paste0("GO_GSEA_", comp_name)), paste("GO GSEA", comp_name))
    go_gsea_map[[comp_name]] <- gsea_go
  }
  if (!is.null(gsea_kegg)) {
    write_gsea_tables(gsea_kegg, file.path(OUTDIR, "2-GSEA", paste0("KEGG_GSEA_", comp_name, ".csv")))
    plot_gsea_suite_pdf(gsea_kegg, file.path(OUTDIR, "3-Visualization", paste0("KEGG_GSEA_", comp_name)), paste("KEGG GSEA", comp_name))
    kegg_gsea_map[[comp_name]] <- gsea_kegg
  }
}

# ---- 6.1 Theme dot-heatmaps ---------------------------------------------------
theme_outdir <- file.path(OUTDIR, "3-Visualization", "ThemeEnrichment")
dir.create(theme_outdir, showWarnings = FALSE, recursive = TRUE)
theme_defs <- default_enrichment_themes()

if (length(go_ora_map) > 0) {
  plot_theme_dotheatmap_from_results(go_ora_map, file.path(theme_outdir, "Theme_dotheatmap_GO_ORA.pdf"),
    title = "GO ORA Biological Themes", subtitle = paste("GO-BP ORA | padj <", DEG_PADJ_CUTOFF, "& |log2FC| >", DEG_LFC_CUTOFF),
    theme_defs = theme_defs, ontology_filter = "BP")
}
if (length(kegg_ora_map) > 0) {
  plot_theme_dotheatmap_from_results(kegg_ora_map, file.path(theme_outdir, "Theme_dotheatmap_KEGG_ORA.pdf"),
    title = "KEGG ORA Pathway Themes", subtitle = paste("KEGG ORA | padj <", DEG_PADJ_CUTOFF, "& |log2FC| >", DEG_LFC_CUTOFF),
    theme_defs = theme_defs, ontology_filter = NULL)
}
if (length(go_gsea_map) > 0) {
  plot_theme_dotheatmap_from_results(go_gsea_map, file.path(theme_outdir, "Theme_dotheatmap_GO_GSEA.pdf"),
    title = "GO GSEA Biological Themes", subtitle = "GO-BP GSEA", theme_defs = theme_defs, ontology_filter = "BP")
}
if (length(kegg_gsea_map) > 0) {
  plot_theme_dotheatmap_from_results(kegg_gsea_map, file.path(theme_outdir, "Theme_dotheatmap_KEGG_GSEA.pdf"),
    title = "KEGG GSEA Pathway Themes", subtitle = "KEGG GSEA", theme_defs = theme_defs, ontology_filter = NULL)
}

# ---- 6.2 Single-term GSEA figures ---------------------------------------------
single_term_outdir <- file.path(OUTDIR, "3-Visualization", "ThemeEnrichment", "single_term_gsea")
dir.create(single_term_outdir, showWarnings = FALSE, recursive = TRUE)

for (comp_name in names(res_list)) {
  ggo   <- go_gsea_map[[comp_name]]
  gkegg <- kegg_gsea_map[[comp_name]]

  if (!is.null(ggo)) {
    ggo_df <- significant_gsea_terms(ggo)
    ggo_df <- ggo_df[order(ggo_df$p.adjust, -abs(ggo_df$NES)), ]
    top_terms <- rbind(utils::head(ggo_df[ggo_df$NES > 0, ], 3),
                       utils::head(ggo_df[ggo_df$NES < 0, ], 3))
    plot_gsea_term_figures_from_df(ggo, top_terms,
      outdir = file.path(single_term_outdir, paste0("GO_", comp_name)),
      contrast_label = comp_name, prefix = "gseaplot2_GO")
  }
  if (!is.null(gkegg)) {
    gkegg_df <- significant_gsea_terms(gkegg)
    gkegg_df <- gkegg_df[order(gkegg_df$p.adjust, -abs(gkegg_df$NES)), ]
    top_terms <- rbind(utils::head(gkegg_df[gkegg_df$NES > 0, ], 3),
                       utils::head(gkegg_df[gkegg_df$NES < 0, ], 3))
    plot_gsea_term_figures_from_df(gkegg, top_terms,
      outdir = file.path(single_term_outdir, paste0("KEGG_", comp_name)),
      contrast_label = comp_name, prefix = "gseaplot2_KEGG")
  }
}

# =============================================================================
# 7. Save results, summary, session info
# =============================================================================
# Targeted save so visualize_results.R can re-plot without re-running voom/GSEA.
save(res_list, deg_summary, v, dge, countData, group, GROUPS, GROUP_LEVELS,
     COMPARISONS, SPECIES, DEG_PADJ_CUTOFF, DEG_LFC_CUTOFF,
     go_ora_map, kegg_ora_map, go_gsea_map, kegg_gsea_map,
     file = file.path(OUTDIR, "1-DEG", "limma_voom_results.Rdata"))

summary_report <- paste0(
  "========================================\n",
  "limma-voom Analysis Summary\n",
  "========================================\n\n",
  "1. Data Overview\n",
  "   - Species: ", SPECIES, "\n",
  "   - Samples: ", ncol(countData), "\n",
  "   - Groups: ", paste(GROUP_LEVELS, collapse = ", "), "\n",
  "   - Batch covariate: ", ifelse(is.null(BATCH_VECTOR), "none", "yes"), "\n",
  "   - Genes after filtering: ", nrow(dge), "\n\n",
  "2. Differential Expression (padj < ", DEG_PADJ_CUTOFF,
  " & |log2FC| > ", DEG_LFC_CUTOFF, ")\n"
)
for (comp_name in names(res_list)) {
  res <- res_list[[comp_name]]
  n_up   <- sum(res$padj < DEG_PADJ_CUTOFF & res$log2FoldChange >  DEG_LFC_CUTOFF, na.rm = TRUE)
  n_down <- sum(res$padj < DEG_PADJ_CUTOFF & res$log2FoldChange < -DEG_LFC_CUTOFF, na.rm = TRUE)
  summary_report <- paste0(summary_report, "   - ", comp_name, ": ",
                           n_up + n_down, " DEGs (Up: ", n_up, ", Down: ", n_down, ")\n")
}
summary_report <- paste0(summary_report, "\n3. Output Files\n",
                         "   - ", file.path(OUTDIR, "0-Config", "analysis_config_used.R"), "\n",
                         "   - ", file.path(OUTDIR, "1-DEG"), " (per-comparison DEG tables, results Rdata)\n",
                         "   - ", file.path(OUTDIR, "2-GSEA"), " (GO/KEGG ORA + GSEA tables)\n",
                         "   - ", file.path(OUTDIR, "3-Visualization"), " (volcano, heatmap, ORA/GSEA, theme maps)\n",
                         "\n========================================\n")
cat(summary_report)
writeLines(summary_report, file.path(OUTDIR, "Analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

# =============================================================================
# 8. HTML report (optional)
# =============================================================================
if (isTRUE(GENERATE_HTML_REPORT)) {
  has_quarto    <- requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))
  has_rmarkdown <- requireNamespace("rmarkdown", quietly = TRUE)
  if (!has_quarto && !has_rmarkdown) {
    message("Skipping HTML report: install the quarto CLI (+ R 'quarto' package) or the 'rmarkdown' package.")
  } else {
    report_path <- tryCatch(
      render_analysis_report(outdir = OUTDIR, report_file = file.path(OUTDIR, "RNAseq_report.html"),
                             template = file.path(dirname(normalizePath(lib_dir)), "reports", "analysis_report.qmd"),
                             params = list(title = REPORT_TITLE, author = Sys.info()[["user"]])),
      error = function(e) { message("HTML report failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(report_path)) cat("HTML report written to:", report_path, "\n")
  }
}

write_template_run_manifest(
  run_dir = OUTDIR, config_path = config_path, input_file = INPUT_FILE,
  input_files = c(annotation_db = AnnotationDbi::dbfile(species_org),
                  kegg_reference = NA_character_),
  config_objects = config_objects, lib_dir = lib_dir, runner_file = file_arg,
  envir = globalenv()
)

cat("\n========================================\n")
cat("RNAseq_Limma_Voom analysis COMPLETE.\n")
cat("Outputs under:", normalizePath(OUTDIR), "\n")
cat("========================================\n")
