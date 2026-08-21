#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_TME — non-interactive TME deconvolution pipeline
# =============================================================================
# Runs the SAME workflow as notebooks/RNAseq_TME_Deconvolution_Template.ipynb,
# driven by a config file instead of a parameter cell. Intended use:
#
#   1. Copy templates/TME/ (config.R + run_analysis.R) into your project.
#   2. Edit config.R (paths, input mode, samples/groups, method switches).
#   3. Run:  Rscript run_analysis.R
#      or:   Rscript run_analysis.R path/to/other_config.R
#
# LIB_DIR resolution: if ./RNAseq_lib exists next to this script it is used;
# otherwise the repository root is located via rprojroot and its RNAseq_lib is
# used. Set the environment variable RNASEQ_LIB_DIR to override explicitly.
#
# The bundled CIBERSORT references are located via rprojroot
# (references/CIBERSORT/) unless CIBERSORT_SCRIPT / CIBERSORT_SIGNATURE are set.
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
cat("RNAseq_TME — run_analysis.R\n")
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

# Project root used to locate bundled resources (references/CIBERSORT/).
repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root, path = script_dir),
                      error = function(e) invocation_dir)

# ---- Optional-module graceful degradation --------------------------------------
# ESTIMATE / IOBR / native CIBERSORT are optional: when the required package or
# bundled script is absent the sub-switch is downgraded to FALSE with a message
# rather than aborting the run (parent+sub-switch pattern from templates/General).
if (isTRUE(RUN_ESTIMATE) && !requireNamespace("estimate", quietly = TRUE)) {
  message("RUN_ESTIMATE is TRUE but the 'estimate' package is not installed; skipping ESTIMATE.")
  RUN_ESTIMATE <- FALSE
}
if (isTRUE(RUN_IOBR) && !requireNamespace("IOBR", quietly = TRUE)) {
  message("RUN_IOBR is TRUE but the 'IOBR' package is not installed; skipping IOBR deconvolution.")
  RUN_IOBR <- FALSE
}
if (isTRUE(RUN_CIBERSORT)) {
  # Resolve auto paths if not provided by the user.
  if (is.null(CIBERSORT_SCRIPT)) {
    CIBERSORT_SCRIPT <- file.path(repo_root, "references", "CIBERSORT", "CIBERSORT.R")
  }
  if (is.null(CIBERSORT_SIGNATURE)) {
    sig_name <- if (SPECIES == "mouse") "cibersort_mouse_22.csv" else "LM22.txt"
    CIBERSORT_SIGNATURE <- file.path(repo_root, "references", "CIBERSORT", sig_name)
  }
  cibersort_pkgs <- c("e1071", "preprocessCore", "future", "furrr", "purrr")
  cibersort_missing <- cibersort_pkgs[!vapply(cibersort_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (!file.exists(CIBERSORT_SCRIPT)) {
    message("RUN_CIBERSORT is TRUE but the CIBERSORT script was not found at ", CIBERSORT_SCRIPT,
            "; skipping native CIBERSORT.")
    RUN_CIBERSORT <- FALSE
  } else if (!file.exists(CIBERSORT_SIGNATURE)) {
    message("RUN_CIBERSORT is TRUE but the signature file was not found at ", CIBERSORT_SIGNATURE,
            "; skipping native CIBERSORT.")
    RUN_CIBERSORT <- FALSE
  } else if (length(cibersort_missing) > 0) {
    message("RUN_CIBERSORT is TRUE but required packages are missing (",
            paste(cibersort_missing, collapse = ", "), "); skipping native CIBERSORT.")
    RUN_CIBERSORT <- FALSE
  }
}
# ssGSEA (sections 8-9) needs GSVA; degrade gracefully if it is unavailable.
RUN_SSGSEA <- TRUE
if (!requireNamespace("GSVA", quietly = TRUE)) {
  message("The 'GSVA' package is not installed; skipping ssGSEA immune-signature scoring.")
  RUN_SSGSEA <- FALSE
}

# ---- Load libraries -----------------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggpubr)
  library(corrplot)
  library(limma)
})
if (isTRUE(RUN_SSGSEA)) suppressPackageStartupMessages(library(GSVA))

for (f in c("plot_utils.R", "io_utils.R", "tme_utils.R", "data_utils.R", "report_utils.R")) {
  source(file.path(lib_dir, f))
}
theme_set(theme_publication())

# ---- Output directories + config snapshot ------------------------------------
CONFIG_DIR <- file.path(OUTDIR, "0-Config")
TME_DIR    <- file.path(OUTDIR, "4-TME")
for (d in c(CONFIG_DIR, TME_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

config_objects <- c(
  "INPUT_MODE", "RAW_COUNTS_FILE", "RAW_COUNTS_FORMAT", "EXPR_FILE", "EXPR_UNIT",
  "META_FILE", "GENE_COLUMN", "GENE_LENGTH_COLUMN", "GENE_LENGTH_UNIT",
  "GENE_START_COL", "GENE_END_COL", "SAMPLE_COLUMN", "GROUP_COLUMN", "GROUP_LEVELS",
  "SPECIES", "GROUP_COLORS", "RUN_ESTIMATE", "RUN_IOBR", "IOBR_METHODS", "IOBR_PERM",
  "IOBR_ARRAYS", "RUN_CIBERSORT", "CIBERSORT_SCRIPT", "CIBERSORT_SIGNATURE",
  "CIBERSORT_PERM", "CIBERSORT_QN", "RUN_CIBERSORT_COMPARISON", "OUTDIR",
  "GENERATE_HTML_REPORT", "REPORT_TITLE"
)
config_lines <- c(
  "# RNAseq_TME analysis configuration snapshot",
  paste0("# Saved: ", Sys.time()),
  unlist(lapply(config_objects, function(obj) {
    if (exists(obj, inherits = TRUE)) c(paste0("\n", obj, " <- "), capture.output(dput(get(obj, inherits = TRUE)))) else NULL
  }))
)
writeLines(config_lines, file.path(CONFIG_DIR, "analysis_config_used.R"))

# =============================================================================
# 3. Load expression and metadata
# =============================================================================
meta <- read_metadata(
  META_FILE,
  sample_column = SAMPLE_COLUMN,
  required_columns = GROUP_COLUMN,
  group_column = GROUP_COLUMN,
  group_levels = GROUP_LEVELS
)
if (is.null(GROUP_LEVELS)) GROUP_LEVELS <- levels(meta[[GROUP_COLUMN]])

# Build a method-valid expression matrix. Raw counts are converted to TPM;
# expression mode only accepts TPM or log2(TPM+1), never VST/rlog.
if (INPUT_MODE == "raw_counts") {
  if (!file.exists(RAW_COUNTS_FILE)) stop("RAW_COUNTS_FILE not found: ", RAW_COUNTS_FILE)
  message("Computing TPM from raw integer counts for TME deconvolution.")

  raw_annot <- read_count_table(RAW_COUNTS_FILE, input_format = RAW_COUNTS_FORMAT)
  if (is.null(GENE_COLUMN)) GENE_COLUMN <- colnames(raw_annot)[1]
  if (!GENE_COLUMN %in% colnames(raw_annot)) stop("GENE_COLUMN not found: ", GENE_COLUMN)
  count_col_names <- detect_count_columns(raw_annot, GENE_COLUMN, annotation_cols = NULL)

  # Build raw count matrix with Ensembl/symbol row names
  counts_mat <- as.matrix(raw_annot[, count_col_names, drop = FALSE])
  mode(counts_mat) <- "numeric"
  rownames(counts_mat) <- as.character(raw_annot[[GENE_COLUMN]])
  validate_count_matrix(counts_mat)

  # Extract gene lengths and compute TPM
  gene_lengths_kb <- extract_gene_lengths(
    raw_annot,
    id_col = GENE_COLUMN,
    length_col = GENE_LENGTH_COLUMN,
    start_col = GENE_START_COL,
    end_col = GENE_END_COL,
    length_unit = GENE_LENGTH_UNIT
  )
  expr_tpm <- counts_to_tpm(counts_mat, gene_lengths_kb)

  # Subset to common samples
  common_samples <- intersect(colnames(expr_tpm), meta[[SAMPLE_COLUMN]])
  if (length(common_samples) == 0) {
    stop("No common samples between RAW_COUNTS_FILE and META_FILE.")
  }
  expr_tpm <- expr_tpm[, common_samples, drop = FALSE]
  meta <- meta[match(common_samples, meta[[SAMPLE_COLUMN]]), ]

  validate_expression_contract(expr_tpm, expected = "tpm")
  expr_tme_input <- as.data.frame(expr_tpm, check.names = FALSE)
  expr_is_log2_tpm <- FALSE

  # Save TPM matrix for reference
  write.csv(expr_tpm, file.path(TME_DIR, "TPM_matrix.csv"))
  cat("TPM matrix saved to", file.path(TME_DIR, "TPM_matrix.csv"), "\n")
} else {
  if (!file.exists(EXPR_FILE)) stop("EXPR_FILE not found: ", EXPR_FILE)
  expr_mat <- read_expression_matrix(EXPR_FILE, gene_column = GENE_COLUMN)

  common_samples <- intersect(colnames(expr_mat), meta[[SAMPLE_COLUMN]])
  if (length(common_samples) == 0) {
    stop("No common samples between EXPR_FILE and META_FILE.")
  }
  expr_mat <- expr_mat[, common_samples, drop = FALSE]
  meta <- meta[match(common_samples, meta[[SAMPLE_COLUMN]]), ]

  validate_expression_contract(expr_mat, expected = EXPR_UNIT)
  expr_is_log2_tpm <- identical(EXPR_UNIT, "log2_tpm")
  expr_tme_input <- as.data.frame(expr_mat, check.names = FALSE)
}

validate_samples_match(colnames(expr_tme_input), meta[[SAMPLE_COLUMN]], strict_order = TRUE)
cat("Expression:", nrow(expr_tme_input), "genes x", ncol(expr_tme_input), "samples\n")
print(table(meta[[GROUP_COLUMN]], useNA = "ifany"))

# Convert mouse gene symbols to HGNC for TME methods (IOBR/ESTIMATE/CIBERSORT/
# EPIC/xCell use human signatures). For mouse Ensembl IDs, convert to MGI first
# then to HGNC. babelgene/biomaRt conversion needs network; the human-symbol path
# does not.
input_id_type <- detect_gene_id_type(rownames(expr_tme_input))
if (SPECIES == "mouse") {
  expr_symbol <- if (input_id_type == "ensembl") {
    convert_expression_rownames(expr_tme_input, species = "mouse", target = "symbol")
  } else expr_tme_input
  expr_tme <- convert_expression_rownames(expr_symbol, species = "mouse", target = "human_symbol")
} else if (input_id_type == "ensembl") {
  expr_tme <- convert_expression_rownames(expr_tme_input, species = "human", target = "symbol")
} else {
  expr_tme <- deduplicate_expression_by_symbol(expr_tme_input)
}

# Only log2(TPM+1) has a valid inverse. VST/rlog is rejected above.
expr_tme <- prepare_tme_expression(expr_tme, is_log = expr_is_log2_tpm, species = "human", verbose = TRUE)
validate_tme_input(expr_tme)
cat("TME expression matrix:", nrow(expr_tme), "genes x", ncol(expr_tme), "samples\n")

# Resolve group colors for consistent plotting
if (is.null(GROUP_COLORS) || !all(GROUP_LEVELS %in% names(GROUP_COLORS))) {
  group_colors <- make_group_colors(GROUP_LEVELS)
} else {
  group_colors <- GROUP_COLORS[GROUP_LEVELS]
}

# NOTE (bug fix): the notebook defines group_df_for_plot only in a LATER
# visualization section, but the native-CIBERSORT section references it earlier,
# so a top-to-bottom run with RUN_CIBERSORT = TRUE fails. Define it here, BEFORE
# any section that uses it.
group_df_for_plot <- meta[, c(SAMPLE_COLUMN, GROUP_COLUMN), drop = FALSE]

# =============================================================================
# 4. ESTIMATE score (native implementation)
# =============================================================================
estimate_scores <- NULL
if (isTRUE(RUN_ESTIMATE)) {
  # ESTIMATE relies on the package lazy-data `common_genes`; requireNamespace()
  # alone does not resolve it, so load it explicitly before calling the API.
  utils::data("common_genes", package = "estimate", envir = environment())
  utils::data("SI_geneset", package = "estimate", envir = environment())
  estimate_df <- data.frame(NAME = rownames(expr_tme), Description = NA, expr_tme, check.names = FALSE)
  write.table(estimate_df, file.path(TME_DIR, "estimate_input.gct"), sep = "\t", quote = FALSE, row.names = FALSE)
  estimate::filterCommonGenes(input.f = file.path(TME_DIR, "estimate_input.gct"),
                              output.f = file.path(TME_DIR, "estimate_common_genes.gct"), id = "GeneSymbol")
  estimate::estimateScore(file.path(TME_DIR, "estimate_common_genes.gct"),
                          file.path(TME_DIR, "estimate_scores.gct"), platform = "illumina")
  estimate_scores <- read.table(file.path(TME_DIR, "estimate_scores.gct"), skip = 2, header = TRUE, sep = "\t", check.names = FALSE)
  rownames(estimate_scores) <- estimate_scores$NAME
  estimate_scores <- as.data.frame(t(estimate_scores[, -(1:2)]))
  estimate_scores[[SAMPLE_COLUMN]] <- rownames(estimate_scores)
  write.csv(estimate_scores, file.path(TME_DIR, "ESTIMATE_scores.csv"), row.names = FALSE)
  cat("ESTIMATE scores saved to", file.path(TME_DIR, "ESTIMATE_scores.csv"), "\n")
}

# =============================================================================
# 5. IOBR multi-algorithm TME deconvolution
# =============================================================================
iobr_results <- list()
tme_combined <- NULL
if (isTRUE(RUN_IOBR)) {
  iobr_results <- run_iobr_deconvolution(
    expr_tme,
    methods = IOBR_METHODS,
    perm = IOBR_PERM,
    arrays = IOBR_ARRAYS,
    id_column = SAMPLE_COLUMN
  )
  # Save individual results
  for (method in names(iobr_results)) {
    write.csv(iobr_results[[method]], file.path(TME_DIR, paste0("IOBR_", method, ".csv")), row.names = FALSE)
  }
  # Combine all results
  if (length(iobr_results) >= 2) {
    tme_combined <- combine_tme_results(iobr_results, id_column = SAMPLE_COLUMN)
    write.csv(tme_combined, file.path(TME_DIR, "IOBR_TME_combined.csv"), row.names = FALSE)
    cat("Combined TME table:", nrow(tme_combined), "samples x", ncol(tme_combined), "features\n")
  }
}

# =============================================================================
# 6. Native CIBERSORT (optional)
# =============================================================================
native_cibersort <- NULL
if (isTRUE(RUN_CIBERSORT)) {
  native_cibersort <- run_native_cibersort(
    expr_tme,
    signature_file = CIBERSORT_SIGNATURE,
    cibersort_script = CIBERSORT_SCRIPT,
    is_log = FALSE,                        # expr_tme is already de-logged by prepare_tme_expression()
    perm = CIBERSORT_PERM,
    QN = CIBERSORT_QN,
    id_column = SAMPLE_COLUMN,
    verbose = TRUE
  )
  write.csv(native_cibersort, file.path(TME_DIR, "CIBERSORT_native_results.csv"), row.names = FALSE)
  cat("Native CIBERSORT results saved to", file.path(TME_DIR, "CIBERSORT_native_results.csv"), "\n")

  # Native CIBERSORT stacked barplot + boxplot + per-cell-type plots
  native_cib_long <- melt_tme_results(native_cibersort, id_column = SAMPLE_COLUMN,
                                       group_df = group_df_for_plot,
                                       sample_col = SAMPLE_COLUMN, group_col = GROUP_COLUMN)
  native_cib_long <- native_cib_long |> dplyr::filter(!grepl("P-value|Correlation|RMSE", .data$cell_type))

  cib_bar_size <- calc_tme_barplot_size(n_samples = length(unique(native_cib_long[[SAMPLE_COLUMN]])),
                                        n_celltypes = length(unique(native_cib_long$cell_type)))
  plot_tme_barplot_pdf(native_cib_long, group_col = GROUP_COLUMN, sample_col = SAMPLE_COLUMN,
                       filename = file.path(TME_DIR, "CIBERSORT_native_barplot.pdf"),
                       title = "Native CIBERSORT Cell Fractions",
                       width = cib_bar_size["width"], height = cib_bar_size["height"])

  cib_box_size <- calc_tme_boxplot_size(n_celltypes = length(unique(native_cib_long$cell_type)))
  plot_tme_boxplot_pdf(native_cib_long, group_col = GROUP_COLUMN, value_col = "fraction",
                       filename = file.path(TME_DIR, "CIBERSORT_native_boxplot.pdf"),
                       title = "Native CIBERSORT Cell Fraction by Group",
                       width = cib_box_size["width"], height = cib_box_size["height"],
                       group_colors = group_colors)

  plot_tme_per_celltype_pdf(
    native_cib_long,
    group_col = GROUP_COLUMN, value_col = "fraction",
    filename_prefix = file.path(TME_DIR, "CIBERSORT_native"),
    title_prefix = "Native CIBERSORT",
    group_colors = group_colors
  )

  # Broad-category aggregation
  native_cib_cat <- aggregate_tme_by_category(
    native_cibersort,
    id_column = SAMPLE_COLUMN,
    category_map = get_cibersort_category_map(SPECIES),
    method = "sum"
  )
  write.csv(native_cib_cat$wide, file.path(TME_DIR, "CIBERSORT_native_category_results.csv"), row.names = FALSE)

  native_cib_cat_long <- native_cib_cat$long |>
    dplyr::rename(cell_type = category, fraction = value)
  if (!is.null(group_df_for_plot)) {
    native_cib_cat_long <- native_cib_cat_long |>
      dplyr::left_join(group_df_for_plot, by = SAMPLE_COLUMN)
  }

  cib_cat_bar_size <- calc_tme_barplot_size(
    n_samples = length(unique(native_cib_cat_long[[SAMPLE_COLUMN]])),
    n_celltypes = length(unique(native_cib_cat_long$cell_type))
  )
  plot_tme_barplot_pdf(native_cib_cat_long, group_col = GROUP_COLUMN, sample_col = SAMPLE_COLUMN,
                       filename = file.path(TME_DIR, "CIBERSORT_native_category_barplot.pdf"),
                       title = "Native CIBERSORT Broad Categories",
                       width = cib_cat_bar_size["width"], height = cib_cat_bar_size["height"])

  cib_cat_box_size <- calc_tme_boxplot_size(n_celltypes = length(unique(native_cib_cat_long$cell_type)))
  plot_tme_boxplot_pdf(native_cib_cat_long, group_col = GROUP_COLUMN, value_col = "fraction",
                       filename = file.path(TME_DIR, "CIBERSORT_native_category_boxplot.pdf"),
                       title = "Native CIBERSORT Broad Categories by Group",
                       width = cib_cat_box_size["width"], height = cib_cat_box_size["height"],
                       group_colors = group_colors)

  plot_tme_per_celltype_pdf(
    native_cib_cat_long,
    group_col = GROUP_COLUMN, value_col = "fraction",
    filename_prefix = file.path(TME_DIR, "CIBERSORT_native_category"),
    title_prefix = "Native CIBERSORT Category",
    group_colors = group_colors
  )
}

# Native vs IOBR CIBERSORT comparison
if (isTRUE(RUN_CIBERSORT) && isTRUE(RUN_IOBR) && "cibersort" %in% names(iobr_results) &&
    isTRUE(RUN_CIBERSORT_COMPARISON) && !is.null(native_cibersort)) {
  # compare_native_iobr_cibersort() strips IOBR's "_CIBERSORT" column suffix and
  # otherwise normalizes cell-type names, so the raw IOBR table can be passed in.
  cmp <- compare_native_iobr_cibersort(
    native_cibersort,
    iobr_results[["cibersort"]],
    id_column = SAMPLE_COLUMN,
    method = "pearson"
  )
  write.csv(cmp$summary, file.path(TME_DIR, "CIBERSORT_native_vs_IOBR_summary.csv"), row.names = FALSE)
  write.csv(cmp$long, file.path(TME_DIR, "CIBERSORT_native_vs_IOBR_long.csv"), row.names = FALSE)
  cat("Native vs IOBR CIBERSORT summary saved to", file.path(TME_DIR, "CIBERSORT_native_vs_IOBR_summary.csv"), "\n")

  n_common_cells <- length(unique(cmp$long$cell_type))
  cmp_width <- min(20, max(10, 2.8 * ceiling(sqrt(n_common_cells))))
  cmp_height <- min(18, max(8, 2.8 * ceiling(n_common_cells / ceiling(sqrt(n_common_cells)))))

  plot_cibersort_correlation_pdf(
    cmp$long,
    filename = file.path(TME_DIR, "CIBERSORT_native_vs_IOBR_correlation.pdf"),
    title = "Native vs IOBR CIBERSORT Fractions",
    width = cmp_width, height = cmp_height
  )
  plot_cibersort_difference_pdf(
    cmp$long,
    filename = file.path(TME_DIR, "CIBERSORT_native_vs_IOBR_difference.pdf"),
    title = "Native - IOBR CIBERSORT Difference",
    width = cmp_width, height = cmp_height
  )
}

# =============================================================================
# 7. TME visualization (IOBR)
# =============================================================================
if (isTRUE(RUN_IOBR) && "estimate" %in% names(iobr_results)) {
  est_long <- melt_estimate_scores(iobr_results[["estimate"]], id_column = SAMPLE_COLUMN,
                                   group_df = group_df_for_plot,
                                   sample_col = SAMPLE_COLUMN, group_col = GROUP_COLUMN)
  plot_estimate_boxplot_pdf(est_long, group_col = GROUP_COLUMN,
    filename = file.path(TME_DIR, "IOBR_ESTIMATE_scores_boxplot.pdf"),
    title = "ESTIMATE Scores by Group", group_colors = group_colors,
    save_individual = TRUE, individual_prefix = file.path(TME_DIR, "IOBR_ESTIMATE"))
  plot_tme_heatmap_pdf(iobr_results[["estimate"]], meta, group_col = GROUP_COLUMN,
    sample_col = SAMPLE_COLUMN, group_colors = group_colors,
    filename = file.path(TME_DIR, "IOBR_ESTIMATE_heatmap.pdf"), title = "IOBR ESTIMATE Scores")
}

# Fraction-based methods: complete composition and group-comparison outputs.
for (method_name in intersect(c("cibersort", "epic"), names(iobr_results))) {
  method_label <- toupper(method_name)
  method_long <- melt_tme_results(iobr_results[[method_name]], id_column = SAMPLE_COLUMN,
    group_df = group_df_for_plot, sample_col = SAMPLE_COLUMN, group_col = GROUP_COLUMN) |>
    dplyr::filter(!grepl("P-value|Correlation|RMSE", .data$cell_type))
  bar_size <- calc_tme_barplot_size(length(unique(method_long[[SAMPLE_COLUMN]])), length(unique(method_long$cell_type)))
  box_size <- calc_tme_boxplot_size(length(unique(method_long$cell_type)))
  prefix <- file.path(TME_DIR, paste0("IOBR_", method_label))
  plot_tme_barplot_pdf(method_long, group_col = GROUP_COLUMN, sample_col = SAMPLE_COLUMN,
    filename = paste0(prefix, "_barplot.pdf"), title = paste(method_label, "Cell Fractions"),
    width = bar_size["width"], height = bar_size["height"])
  plot_tme_boxplot_pdf(method_long, group_col = GROUP_COLUMN, value_col = "fraction",
    filename = paste0(prefix, "_boxplot.pdf"), title = paste(method_label, "Cell Fractions by Group"),
    width = box_size["width"], height = box_size["height"], group_colors = group_colors)
  plot_tme_per_celltype_pdf(method_long, group_col = GROUP_COLUMN, value_col = "fraction",
    filename_prefix = prefix, title_prefix = method_label, group_colors = group_colors)
}

# Broad CIBERSORT categories use the project-maintained manual mapping.
# Workaround (library limitation): aggregate_tme_by_category() defaults to
# suffix_strip = "_xCell$", but IOBR cibersort columns carry a "_CIBERSORT"
# suffix, so the suffix must be passed explicitly or no cell types match.
if (isTRUE(RUN_IOBR) && "cibersort" %in% names(iobr_results)) {
  cib_cat <- aggregate_tme_by_category(iobr_results[["cibersort"]], id_column = SAMPLE_COLUMN,
    category_map = get_cibersort_category_map(SPECIES), method = "sum",
    suffix_strip = "_CIBERSORT$")
  write.csv(cib_cat$wide, file.path(TME_DIR, "IOBR_CIBERSORT_category_results.csv"), row.names = FALSE)
}

if (isTRUE(RUN_IOBR) && "xcell" %in% names(iobr_results)) {
  plot_tme_heatmap_pdf(iobr_results[["xcell"]], meta, group_col = GROUP_COLUMN,
    sample_col = SAMPLE_COLUMN, group_colors = group_colors,
    filename = file.path(TME_DIR, "IOBR_xCell_heatmap.pdf"), title = "xCell Scores", width = 10, height = 12)
  xcell_cat <- aggregate_tme_by_category(iobr_results[["xcell"]], id_column = SAMPLE_COLUMN,
    category_map = get_xcell_category_map(), method = "mean")
  write.csv(xcell_cat$wide, file.path(TME_DIR, "IOBR_xCell_category_results.csv"), row.names = FALSE)
  plot_tme_heatmap_pdf(xcell_cat$wide, meta, group_col = GROUP_COLUMN, sample_col = SAMPLE_COLUMN,
    group_colors = group_colors, filename = file.path(TME_DIR, "IOBR_xCell_category_heatmap.pdf"),
    title = "xCell Broad Categories", width = 8, height = 7)
}

# =============================================================================
# 8. ssGSEA immune signature scoring
# =============================================================================
# The bundled immune signatures are human gene symbols. Use the same prepared
# expression matrix that goes into TME methods (already non-log TPM/HGNC for
# mouse, or the loaded normalized matrix for human).
ssgsea_scores <- NULL
gs <- list()
if (isTRUE(RUN_SSGSEA)) {
  expr_for_ssgsea <- as.data.frame(expr_tme, check.names = FALSE)

  row_upper <- toupper(rownames(expr_for_ssgsea))
  upper_to_real <- setNames(rownames(expr_for_ssgsea), row_upper)
  gs <- lapply(immune_gene_sets, function(x) unique(upper_to_real[intersect(toupper(x), names(upper_to_real))]))
  gs <- gs[lengths(gs) >= 2]

  if (length(gs) == 0) {
    message("No immune signature genes matched the expression matrix; skipping ssGSEA. ",
            "Check that row names are human gene symbols (or that mouse-to-human conversion succeeded).")
    RUN_SSGSEA <- FALSE
  } else {
    params <- gsvaParam(as.matrix(expr_for_ssgsea), gs, kcdf = "Gaussian", minSize = 2, maxSize = Inf)
    ssgsea_scores <- gsva(params, verbose = FALSE)
    write.csv(ssgsea_scores, file.path(TME_DIR, "ssGSEA_immune_scores.csv"))
    cat("ssGSEA immune scores saved to", file.path(TME_DIR, "ssGSEA_immune_scores.csv"), "\n")
  }
}

# =============================================================================
# 9. ssGSEA group comparison and heatmap
# =============================================================================
if (isTRUE(RUN_SSGSEA) && !is.null(ssgsea_scores)) {
  score_df <- as.data.frame(t(ssgsea_scores)) %>% rownames_to_column(SAMPLE_COLUMN) %>% left_join(meta, by = SAMPLE_COLUMN)
  score_long <- score_df %>% pivot_longer(cols = names(gs), names_to = "signature", values_to = "score")

  score_long[[GROUP_COLUMN]] <- factor(score_long[[GROUP_COLUMN]], levels = GROUP_LEVELS)
  score_comparisons <- if (length(GROUP_LEVELS) >= 2) combn(GROUP_LEVELS, 2, simplify = FALSE) else list()
  plot_group_boxplot_pdf(
    score_long,
    value_col = "score",
    group_col = GROUP_COLUMN,
    facet_col = "signature",
    comparisons = score_comparisons,
    method = "t.test",
    title = "ssGSEA Immune Signature Scores",
    ylab = "ssGSEA score",
    group_colors = group_colors,
    filename = file.path(TME_DIR, "ssGSEA_group_boxplot.pdf"),
    width = 12,
    height = 8
  )

  ann <- data.frame(Group = as.character(meta[[GROUP_COLUMN]]))
  rownames(ann) <- meta[[SAMPLE_COLUMN]]
  ann$Group <- factor(ann$Group, levels = GROUP_LEVELS)
  annotation_colors <- list(Group = group_colors)

  # Order samples by group, then cluster within each group, so replicates of the
  # same condition appear together while preserving within-group structure.
  order_grouped_columns_ssgsea <- function(m, group_vec, group_levels = NULL) {
    if (is.null(group_levels)) group_levels <- unique(group_vec)
    ordered_cols <- character(0)
    for (g in group_levels) {
      idx <- which(group_vec == g)
      if (length(idx) == 0) next
      if (length(idx) == 1) {
        ordered_cols <- c(ordered_cols, colnames(m)[idx])
      } else {
        sub_m <- m[, idx, drop = FALSE]
        sd_rows <- apply(sub_m, 1, stats::sd, na.rm = TRUE)
        sub_m_var <- sub_m[sd_rows > 0 | is.na(sd_rows), , drop = FALSE]
        if (ncol(sub_m_var) >= 2 && nrow(sub_m_var) >= 2) {
          d <- stats::dist(t(sub_m_var))
          hc <- stats::hclust(d)
          ordered_cols <- c(ordered_cols, colnames(sub_m)[hc$order])
        } else {
          ordered_cols <- c(ordered_cols, colnames(sub_m))
        }
      }
    }
    ordered_cols
  }

  ssgsea_col_order <- order_grouped_columns_ssgsea(ssgsea_scores, ann$Group, GROUP_LEVELS)
  ssgsea_scores_plot <- ssgsea_scores[, ssgsea_col_order, drop = FALSE]
  ann_plot <- ann[ssgsea_col_order, , drop = FALSE]

  pheatmap(ssgsea_scores_plot, annotation_col = ann_plot, annotation_colors = annotation_colors,
           scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
           filename = file.path(TME_DIR, "ssGSEA_heatmap.pdf"), width = 8, height = 7)
}

# =============================================================================
# 10. Save results, summary, session info
# =============================================================================
# Targeted save so visualize_results.R can re-plot without re-running deconvolution.
save(expr_tme, expr_tme_input, iobr_results, tme_combined, estimate_scores,
     native_cibersort, ssgsea_scores, gs, meta, group_df_for_plot, group_colors,
     SAMPLE_COLUMN, GROUP_COLUMN, GROUP_LEVELS, SPECIES, GROUP_COLORS,
     RUN_ESTIMATE, RUN_IOBR, RUN_CIBERSORT, RUN_SSGSEA,
     file = file.path(TME_DIR, "tme_results.Rdata"))

summary_report <- paste0(
  "========================================\n",
  "TME Deconvolution Analysis Summary\n",
  "========================================\n\n",
  "1. Data Overview\n",
  "   - Input mode: ", INPUT_MODE, "\n",
  "   - Species: ", SPECIES, "\n",
  "   - Samples: ", ncol(expr_tme), "\n",
  "   - Groups: ", paste(GROUP_LEVELS, collapse = ", "), "\n",
  "   - TME genes after prep: ", nrow(expr_tme), "\n\n",
  "2. Methods Run\n",
  "   - ESTIMATE (native): ", ifelse(isTRUE(RUN_ESTIMATE) && !is.null(estimate_scores), "yes", "no"), "\n",
  "   - IOBR (", ifelse(length(iobr_results) > 0, paste(names(iobr_results), collapse = ", "), "none"), ")\n",
  "   - Native CIBERSORT: ", ifelse(!is.null(native_cibersort), "yes", "no"), "\n",
  "   - ssGSEA signatures scored: ", ifelse(!is.null(ssgsea_scores), nrow(ssgsea_scores), 0), "\n\n",
  "3. Output Files\n",
  "   - ", file.path(CONFIG_DIR, "analysis_config_used.R"), "\n",
  "   - ", TME_DIR, " (TPM matrix, ESTIMATE/IOBR/CIBERSORT tables, ssGSEA, all PDFs, tme_results.Rdata)\n",
  "\n========================================\n"
)
cat(summary_report)
writeLines(summary_report, file.path(OUTDIR, "Analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))

# =============================================================================
# 11. HTML report (optional)
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

cat("\n========================================\n")
cat("RNAseq_TME analysis COMPLETE.\n")
cat("Outputs under:", normalizePath(OUTDIR), "\n")
cat("========================================\n")
