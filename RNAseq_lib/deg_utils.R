# Differential expression helpers for bulk RNA-seq templates.

validate_threshold_grid <- function(threshold_grid, default_threshold = NULL) {
  if ("padj" %in% colnames(threshold_grid) && !"p_cutoff" %in% colnames(threshold_grid)) {
    threshold_grid$p_cutoff <- threshold_grid$padj
  }
  required_cols <- c("name", "p_cutoff", "log2fc")
  missing_cols <- setdiff(required_cols, colnames(threshold_grid))
  if (length(missing_cols) > 0) {
    stop("THRESHOLD_GRID is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (any(duplicated(threshold_grid$name))) {
    stop("THRESHOLD_GRID$name must be unique.")
  }
  if (!is.null(default_threshold) && !default_threshold %in% threshold_grid$name) {
    stop("DEFAULT_THRESHOLD must be one of THRESHOLD_GRID$name.")
  }
  threshold_grid
}

run_deseq2_model <- function(count_data, col_data, design_formula) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = design_formula
  )
  dds <- DESeq2::DESeq(dds)
  dds
}

extract_deseq2_results <- function(dds, comparisons, condition_col = "condition", alpha = 0.05, shrink_type = "ashr") {
  res_list <- list()
  for (comp in comparisons) {
    comp_name <- comp[1]
    treat <- comp[2]
    ctrl <- comp[3]
    message("Processing DESeq2 contrast: ", comp_name, " (", treat, " vs ", ctrl, ")")
    res <- DESeq2::results(dds, contrast = c(condition_col, treat, ctrl), alpha = alpha)
    if (!is.null(shrink_type)) {
      res <- DESeq2::lfcShrink(dds, contrast = c(condition_col, treat, ctrl), res = res, type = shrink_type)
    }
    res_df <- as.data.frame(res)
    res_df$gene_name <- rownames(res_df)
    res_list[[comp_name]] <- res_df
  }
  res_list
}

validate_pvalue_column <- function(res_df, pvalue_column = "padj") {
  if (!pvalue_column %in% colnames(res_df)) {
    stop("P-value column not found in DESeq2 result: ", pvalue_column)
  }
  if (!pvalue_column %in% c("padj", "pvalue")) {
    warning("Using non-standard p-value column: ", pvalue_column)
  }
  invisible(TRUE)
}

mark_deg_by_threshold <- function(res_df, pvalue_thresh, log2fc_thresh, pvalue_column = "padj") {
  validate_pvalue_column(res_df, pvalue_column)
  out <- res_df
  pvals <- out[[pvalue_column]]
  out$significance <- ifelse(
    !is.na(pvals) & pvals < pvalue_thresh & abs(out$log2FoldChange) > log2fc_thresh,
    ifelse(out$log2FoldChange > 0, "Up", "Down"),
    "Not_Sig"
  )
  out
}

build_deg_threshold_sets <- function(res_list, threshold_grid, pvalue_column = "padj") {
  deg_by_threshold <- list()
  for (i in seq_len(nrow(threshold_grid))) {
    th <- threshold_grid[i, ]
    deg_by_threshold[[th$name]] <- lapply(
      res_list,
      mark_deg_by_threshold,
      pvalue_thresh = th$p_cutoff,
      log2fc_thresh = th$log2fc,
      pvalue_column = pvalue_column
    )
  }
  deg_by_threshold
}

summarize_deg_thresholds <- function(deg_by_threshold) {
  do.call(rbind, lapply(names(deg_by_threshold), function(th_name) {
    do.call(rbind, lapply(names(deg_by_threshold[[th_name]]), function(comp_name) {
      res <- deg_by_threshold[[th_name]][[comp_name]]
      data.frame(
        Threshold = th_name,
        Comparison = comp_name,
        UP = sum(res$significance == "Up", na.rm = TRUE),
        DOWN = sum(res$significance == "Down", na.rm = TRUE),
        Total = sum(res$significance != "Not_Sig", na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

get_threshold_result <- function(deg_by_threshold, threshold_name) {
  if (!threshold_name %in% names(deg_by_threshold)) {
    stop("Threshold not found: ", threshold_name)
  }
  deg_by_threshold[[threshold_name]]
}

write_deg_threshold_outputs <- function(deg_by_threshold, outdir = "1-DEG") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (th_name in names(deg_by_threshold)) {
    th_dir <- file.path(outdir, th_name)
    dir.create(th_dir, showWarnings = FALSE, recursive = TRUE)
    for (comp_name in names(deg_by_threshold[[th_name]])) {
      utils::write.csv(
        deg_by_threshold[[th_name]][[comp_name]],
        file.path(th_dir, paste0("DEG_results_", comp_name, ".csv")),
        row.names = FALSE
      )
    }
  }
  summary_df <- summarize_deg_thresholds(deg_by_threshold)
  utils::write.csv(summary_df, file.path(outdir, "DEG_threshold_summary.csv"), row.names = FALSE)
  summary_df
}

genes_for_enrichment <- function(res_df, pvalue_thresh, log2fc_thresh, pvalue_column = "padj") {
  validate_pvalue_column(res_df, pvalue_column)
  sig <- res_df |>
    dplyr::filter(!is.na(.data[[pvalue_column]]), .data[[pvalue_column]] < pvalue_thresh, abs(log2FoldChange) > log2fc_thresh) |>
    dplyr::pull(gene_name)
  up <- res_df |>
    dplyr::filter(!is.na(.data[[pvalue_column]]), .data[[pvalue_column]] < pvalue_thresh, log2FoldChange > log2fc_thresh) |>
    dplyr::pull(gene_name)
  down <- res_df |>
    dplyr::filter(!is.na(.data[[pvalue_column]]), .data[[pvalue_column]] < pvalue_thresh, log2FoldChange < -log2fc_thresh) |>
    dplyr::pull(gene_name)
  list(sig = sig, up = up, down = down)
}

ranked_gene_list <- function(res_df) {
  gene_list <- res_df$log2FoldChange
  names(gene_list) <- res_df$gene_name
  gene_list <- gene_list[!is.na(gene_list)]
  sort(gene_list, decreasing = TRUE)
}
