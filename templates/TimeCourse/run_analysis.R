#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_TimeCourse — non-interactive analysis pipeline
# =============================================================================
# Runs the SAME workflow as notebooks/RNAseq_TimeCourse_Template.ipynb, driven by
# a config file instead of a parameter cell. Intended use:
#
#   1. Copy templates/TimeCourse/ (config.R + run_analysis.R) into your project.
#   2. Edit config.R (paths, time column, Mfuzz + DEG settings).
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
cat("RNAseq_TimeCourse — run_analysis.R\n")
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
# BUG FIX (vs the notebook): select the OrgDb from SPECIES instead of the
# hard-coded org.Hs.eg.db the notebook used for the Mfuzz cluster ORA.
species_org_db <- if (SPECIES == "human") "org.Hs.eg.db" else "org.Mm.eg.db"
if (!requireNamespace(species_org_db, quietly = TRUE)) {
  stop("Annotation package '", species_org_db, "' is required. ",
       "Install with BiocManager::install('", species_org_db, "')")
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(DESeq2)
  library(species_org_db, character.only = TRUE)
})
org_db <- get(species_org_db)

for (f in c("plot_utils.R", "io_utils.R", "data_utils.R", "deg_utils.R",
            "enrichment_utils.R", "timecourse_utils.R", "report_utils.R", "run_utils.R")) {
  source(file.path(lib_dir, f))
}
theme_set(theme_publication())
initialize_run_lifecycle(globalenv())

# Mfuzz clustering is optional and has a heavy dependency. Degrade gracefully when
# the package is missing so aggregation + time-point DEG still run (mirrors the
# RUN_TME lazy pattern in templates/General/run_analysis.R).
if (isTRUE(RUN_MFUZZ) && !requireNamespace("Mfuzz", quietly = TRUE)) {
  message("RUN_MFUZZ is TRUE but the 'Mfuzz' package is not installed; skipping Mfuzz clustering.")
  RUN_MFUZZ <- FALSE
}
if (isTRUE(RUN_MFUZZ)) {
  suppressPackageStartupMessages({
    library(Mfuzz)
    library(ComplexHeatmap)
    library(circlize)
  })
}

# ---- Output directories + config snapshot ------------------------------------
for (d in c("0-Config", "1-DEG_Timepoint", "3-Visualization", "5-TimeCourse")) {
  dir.create(file.path(OUTDIR, d), showWarnings = FALSE, recursive = TRUE)
}

config_objects <- c(
  "SPECIES", "EXPR_FILE", "META_FILE", "GENE_COLUMN", "SAMPLE_COLUMN",
  "TIME_COLUMN", "GROUP_COLUMN", "TIME_LEVELS", "RUN_MFUZZ", "MFUZZ_N_CLUSTERS",
  "MFUZZ_MIN_ACORE", "MFUZZ_SEED", "RAW_COUNTS_FILE", "COUNT_META_FILE",
  "COUNT_GENE_COL", "COUNT_SAMPLE_COL", "COUNT_BIOTYPE_COL", "COUNT_BIOTYPE_FILTER",
  "SUBJECT_COL", "RUN_TIMEPOINT_DEG", "BASELINE_TIME", "DEG_PADJ_CUTOFF",
  "DEG_LFC_CUTOFF", "MIN_COUNT", "OUTDIR", "GENERATE_HTML_REPORT", "REPORT_TITLE",
  "RUN_ROLE", "PARENT_RUN_ID", "RUN_CHANGE_NOTE", "RUN_RETENTION"
)
config_lines <- c(
  "# RNAseq_TimeCourse analysis configuration snapshot",
  paste0("# Saved: ", Sys.time()),
  unlist(lapply(config_objects, function(obj) {
    if (exists(obj, inherits = TRUE)) c(paste0("\n", obj, " <- "), capture.output(dput(get(obj, inherits = TRUE)))) else NULL
  }))
)
writeLines(config_lines, file.path(OUTDIR, "0-Config", "analysis_config_used.R"))

# Mfuzz-specific outputs live under 5-TimeCourse/ to keep the numbered areas clean.
tc_dir <- file.path(OUTDIR, "5-TimeCourse")

# =============================================================================
# 3. Load expression and metadata
# =============================================================================
expr <- read_expression_matrix(EXPR_FILE, gene_column = GENE_COLUMN)

required_meta_cols <- c(SAMPLE_COLUMN, TIME_COLUMN)
if (!is.null(GROUP_COLUMN)) required_meta_cols <- c(required_meta_cols, GROUP_COLUMN)
meta <- read_metadata(
  META_FILE,
  sample_column = SAMPLE_COLUMN,
  required_columns = required_meta_cols,
  time_column = TIME_COLUMN,
  time_levels = TIME_LEVELS,
  group_column = GROUP_COLUMN
)
if (is.null(TIME_LEVELS)) TIME_LEVELS <- levels(meta[[TIME_COLUMN]])

validate_samples_match(colnames(expr), meta[[SAMPLE_COLUMN]], strict_order = TRUE)
expr <- expr[, meta[[SAMPLE_COLUMN]], drop = FALSE]

# Mfuzz/clustering requires normalized expression derived from raw counts.
scale_info <- validate_expression_contract(expr, expected = "vst")

cat("Expression:", nrow(expr), "genes x", ncol(expr), "samples\n")
print(table(meta[[TIME_COLUMN]], useNA = "ifany"))

# =============================================================================
# 4. Aggregate expression by time point
# =============================================================================
expr_mean <- aggregate_expr_by_group(expr, meta[[TIME_COLUMN]])
expr_mean <- expr_mean[, TIME_LEVELS, drop = FALSE]  # ensure order
write.csv(expr_mean, file.path(tc_dir, "mean_expression_by_time.csv"))
cat("Aggregated expression:", nrow(expr_mean), "genes x", ncol(expr_mean), "time points\n")

# =============================================================================
# 5. Mfuzz time-course soft clustering
# =============================================================================
cluster_df  <- NULL
mfuzz_result <- NULL
eset        <- NULL
ora_results <- list()
if (isTRUE(RUN_MFUZZ)) {
  eset <- prepare_mfuzz_eset(expr_mean)
  set.seed(MFUZZ_SEED)  # honor MFUZZ_SEED for reproducible clustering
  mfuzz_result <- run_mfuzz(eset, n_clusters = MFUZZ_N_CLUSTERS, seed = MFUZZ_SEED)
  cluster_df <- extract_mfuzz_clusters(mfuzz_result, eset = eset, min_acore = MFUZZ_MIN_ACORE)
  write_mfuzz_cluster_table(cluster_df, file.path(tc_dir, "mfuzz_clusters.csv"))

  cat("Mfuzz cluster sizes:\n")
  print(summarize_mfuzz_clusters(cluster_df))

  # Trend plots
  plot_mfuzz_trends_pdf(
    eset, mfuzz_result,
    filename = file.path(tc_dir, "mfuzz_trends.pdf"),
    time_labels = TIME_LEVELS,
    width = mm_to_in(183), height = 6.4
  )

  # Heatmap of core genes by cluster
  core_df <- cluster_df[cluster_df$core_gene, ]
  if (nrow(core_df) >= 2) {
    group_colors <- make_group_colors(TIME_LEVELS)
    plot_timecourse_heatmap_pdf(
      expr, core_df,
      group_vec = meta[[TIME_COLUMN]],
      group_levels = TIME_LEVELS,
      group_colors = group_colors,
      filename = file.path(tc_dir, "mfuzz_core_heatmap.pdf"),
      width = mm_to_in(183), height = mm_to_in(247)
    )
  }

  # ORA per cluster (species-aware OrgDb; see BUG FIX above)
  universe <- map_symbols_to_entrez(rownames(expr), org_db)$ENTREZID
  ora_results <- run_mfuzz_cluster_ora(cluster_df, org_db = org_db, universe = universe)
  for (cl_name in names(ora_results)) {
    prefix <- file.path(tc_dir, paste0("GO_ORA_", cl_name))
    write.csv(as.data.frame(ora_results[[cl_name]]), paste0(prefix, ".csv"), row.names = FALSE)
    plot_enrich_suite_pdf(ora_results[[cl_name]], prefix, cl_name)
  }

  # ---- 5.1 Publication-grade theme dot-heatmap ------------------------------
  theme_outdir <- file.path(OUTDIR, "3-Visualization", "ThemeEnrichment")
  dir.create(theme_outdir, showWarnings = FALSE, recursive = TRUE)
  theme_defs <- default_enrichment_themes()

  # WORKAROUND: the notebook's drop_empty uses sapply(), which returns list() on an
  # empty input and crashes `lst[list()]`. vapply(logical(1)) is robust when no
  # cluster produced a significant ORA result (common for sparse/simulated data).
  drop_empty <- function(lst) lst[vapply(lst, function(x) !is.null(x) && nrow(as.data.frame(x)) > 0, logical(1))]
  ora_map <- drop_empty(ora_results)

  if (length(ora_map) > 0) {
    plot_theme_dotheatmap_from_results(
      ora_map,
      filename = file.path(theme_outdir, "Theme_dotheatmap_GO_ORA_mfuzz_clusters.pdf"),
      title = "GO ORA Biological Themes (Mfuzz Clusters)",
      subtitle = "GO-BP ORA | top terms per theme per cluster",
      theme_defs = theme_defs,
      ontology_filter = "BP",
      top_n = 6
    )
  }
} else {
  cat("Mfuzz clustering skipped (RUN_MFUZZ = FALSE).\n")
}

# =============================================================================
# 6. Time-point vs baseline DEG (optional)
# =============================================================================
tp_res_list <- list()
tp_summary  <- NULL
baseline_time <- NA_character_
if (isTRUE(RUN_TIMEPOINT_DEG)) {
  if (!file.exists(RAW_COUNTS_FILE)) {
    stop("RAW_COUNTS_FILE not found: ", RAW_COUNTS_FILE,
         "\nProvide raw integer counts with columns for gene_name and samples, ",
         "or set RUN_TIMEPOINT_DEG <- FALSE.")
  }

  rawcount <- read_count_table(RAW_COUNTS_FILE, "tsv")
  if (!is.null(COUNT_BIOTYPE_COL) && COUNT_BIOTYPE_COL %in% colnames(rawcount)) {
    rawcount <- rawcount[rawcount[[COUNT_BIOTYPE_COL]] == COUNT_BIOTYPE_FILTER, ]
  }
  count_col_names <- detect_count_columns(rawcount, COUNT_GENE_COL, NULL)
  count_meta <- read.csv(COUNT_META_FILE, check.names = FALSE)
  count_samples <- intersect(count_col_names, count_meta[[COUNT_SAMPLE_COL]])
  if (length(count_samples) == 0) {
    stop("No matching samples between raw counts and count metadata.")
  }
  rawcount <- rawcount[, c(COUNT_GENE_COL, count_samples), drop = FALSE]
  count_meta <- count_meta[match(count_samples, count_meta[[COUNT_SAMPLE_COL]]), ]

  countData <- build_count_matrix(
    rawcount, COUNT_GENE_COL, count_samples, count_meta[[COUNT_SAMPLE_COL]],
    duplicate_report_file = file.path(OUTDIR, "1-DEG_Timepoint", "Duplicated_gene_symbols.csv")
  )
  countData <- filter_low_count_genes(countData, count_meta[[TIME_COLUMN]], MIN_COUNT)$count_data

  col_data <- data.frame(
    sample = count_meta[[COUNT_SAMPLE_COL]],
    stringsAsFactors = FALSE
  )
  col_data[[TIME_COLUMN]] <- count_meta[[TIME_COLUMN]]
  # Only include the group/condition column in the design when it actually varies
  # across samples; a single-level condition breaks the DESeq2 design.
  condition_col <- NULL
  if (!is.null(GROUP_COLUMN) && GROUP_COLUMN %in% colnames(count_meta) &&
      length(unique(count_meta[[GROUP_COLUMN]])) > 1) {
    col_data[[GROUP_COLUMN]] <- count_meta[[GROUP_COLUMN]]
    condition_col <- GROUP_COLUMN
  }
  if (!is.null(SUBJECT_COL) && SUBJECT_COL %in% colnames(count_meta)) {
    col_data[[SUBJECT_COL]] <- count_meta[[SUBJECT_COL]]
  }
  rownames(col_data) <- col_data$sample

  if (is.null(BASELINE_TIME)) {
    baseline_time <- TIME_LEVELS[1]
  } else {
    baseline_time <- BASELINE_TIME
  }

  cat("Running time-point vs baseline DEG with baseline:", baseline_time, "\n")
  tp_res_list <- run_timepoint_vs_baseline_deseq2(
    count_data = countData,
    col_data = col_data,
    time_col = TIME_COLUMN,
    baseline_time = baseline_time,
    condition_col = condition_col,
    subject_col = SUBJECT_COL,
    alpha = max(DEG_PADJ_CUTOFF, 0.05)
  )

  tp_deg_dir <- file.path(OUTDIR, "1-DEG_Timepoint")
  tp_summary <- write_timepoint_deg_results(
    tp_res_list,
    outdir = tp_deg_dir,
    pvalue_column = "padj",
    lfc_column = "log2FoldChange_shrunken"
  )

  plot_timepoint_deg_summary_pdf(
    tp_summary,
    filename = file.path(OUTDIR, "3-Visualization", "Timepoint_DEG_summary.pdf")
  )

  for (comp_name in names(tp_res_list)) {
    plot_volcano_pdf(
      tp_res_list[[comp_name]],
      comp_name = comp_name,
      pvalue_thresh = DEG_PADJ_CUTOFF,
      log2fc_thresh = DEG_LFC_CUTOFF,
      pvalue_column = "padj",
      lfc_column = "log2FoldChange_shrunken",
      filename = file.path(OUTDIR, "3-Visualization", paste0("Volcano_", comp_name, ".pdf"))
    )
  }

  cat("Time-point vs baseline DEG complete. Results saved to", tp_deg_dir, "\n")
} else {
  cat("Time-point DEG skipped (RUN_TIMEPOINT_DEG = FALSE).\n")
}

# =============================================================================
# 7. Save results, summary, session info
# =============================================================================
# Targeted save so visualize_results.R can re-plot trends/heatmaps/theme maps
# without re-running Mfuzz / DESeq2.
save(expr_mean, cluster_df, mfuzz_result, eset, ora_results,
     tp_res_list, tp_summary,
     TIME_LEVELS, SPECIES, RUN_MFUZZ, RUN_TIMEPOINT_DEG,
     MFUZZ_N_CLUSTERS, MFUZZ_MIN_ACORE, MFUZZ_SEED,
     baseline_time, DEG_PADJ_CUTOFF, DEG_LFC_CUTOFF,
     file = file.path(tc_dir, "timecourse_results.Rdata"))

summary_report <- paste0(
  "========================================\n",
  "Time-Course (Mfuzz) Analysis Summary\n",
  "========================================\n\n",
  "1. Data Overview\n",
  "   - Species: ", SPECIES, "\n",
  "   - Samples: ", ncol(expr), "\n",
  "   - Genes: ", nrow(expr), "\n",
  "   - Time points: ", paste(TIME_LEVELS, collapse = ", "), "\n\n"
)
if (isTRUE(RUN_MFUZZ) && !is.null(cluster_df)) {
  summary_report <- paste0(summary_report,
    "2. Mfuzz Soft Clustering\n",
    "   - Clusters: ", MFUZZ_N_CLUSTERS, " (seed ", MFUZZ_SEED, ")\n",
    "   - Clustered genes: ", nrow(cluster_df), "\n",
    "   - Core genes (acore >= ", MFUZZ_MIN_ACORE, "): ", sum(cluster_df$core_gene), "\n",
    "   - Clusters with ORA: ", length(ora_results), "\n\n")
} else {
  summary_report <- paste0(summary_report, "2. Mfuzz Soft Clustering\n   - skipped\n\n")
}
if (isTRUE(RUN_TIMEPOINT_DEG) && !is.null(tp_summary)) {
  summary_report <- paste0(summary_report,
    "3. Time-point vs Baseline DEG (baseline = ", baseline_time,
    "; padj < ", DEG_PADJ_CUTOFF, " & |log2FC| > ", DEG_LFC_CUTOFF, ")\n")
  for (i in seq_len(nrow(tp_summary))) {
    summary_report <- paste0(summary_report, "   - ", tp_summary$Comparison[i], ": ",
                             tp_summary$Total[i], " DEGs (Up: ", tp_summary$UP[i],
                             ", Down: ", tp_summary$DOWN[i], ")\n")
  }
  summary_report <- paste0(summary_report, "\n")
} else {
  summary_report <- paste0(summary_report, "3. Time-point vs Baseline DEG\n   - skipped\n\n")
}
summary_report <- paste0(summary_report,
  "4. Output Files\n",
  "   - ", file.path(OUTDIR, "0-Config", "analysis_config_used.R"), "\n",
  "   - ", file.path(OUTDIR, "1-DEG_Timepoint"), " (time-point DEG tables)\n",
  "   - ", file.path(OUTDIR, "3-Visualization"), " (volcano, DEG summary, theme maps)\n",
  "   - ", file.path(OUTDIR, "5-TimeCourse"), " (mean expression, Mfuzz clusters/trends, ORA, results Rdata)\n",
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
                             params = list(title = REPORT_TITLE, author = Sys.info()[["user"]])),
      error = function(e) { message("HTML report failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(report_path)) cat("HTML report written to:", report_path, "\n")
  }
}

write_template_run_manifest(
  run_dir = OUTDIR, config_path = config_path, input_file = EXPR_FILE,
  config_objects = config_objects, lib_dir = lib_dir, runner_file = file_arg,
  envir = globalenv()
)

cat("\n========================================\n")
cat("RNAseq_TimeCourse analysis COMPLETE.\n")
cat("Outputs under:", normalizePath(OUTDIR), "\n")
cat("========================================\n")
