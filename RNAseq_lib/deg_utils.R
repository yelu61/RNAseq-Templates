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

    raw_res <- DESeq2::results(dds, contrast = c(condition_col, treat, ctrl), alpha = alpha)
    raw_df <- as.data.frame(raw_res)

    out <- data.frame(
      baseMean = raw_df$baseMean,
      log2FoldChange_raw = raw_df$log2FoldChange,
      lfcSE_raw = raw_df$lfcSE,
      stat = raw_df$stat,
      pvalue = raw_df$pvalue,
      padj = raw_df$padj,
      stringsAsFactors = FALSE
    )

    if (!is.null(shrink_type)) {
      shrink_res <- DESeq2::lfcShrink(dds, contrast = c(condition_col, treat, ctrl), res = raw_res, type = shrink_type)
      shrink_df <- as.data.frame(shrink_res)
      out$log2FoldChange_shrunken <- shrink_df$log2FoldChange
      out$lfcSE_shrunken <- shrink_df$lfcSE
    } else {
      out$log2FoldChange_shrunken <- out$log2FoldChange_raw
      out$lfcSE_shrunken <- out$lfcSE_raw
    }

    out$log2FoldChange <- out$log2FoldChange_shrunken
    out$lfcSE <- out$lfcSE_shrunken
    out$gene_name <- rownames(raw_df)
    res_list[[comp_name]] <- out
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

validate_lfc_column <- function(res_df, lfc_column = "log2FoldChange") {
  if (!lfc_column %in% colnames(res_df)) {
    stop("Log2 fold-change column not found in DESeq2 result: ", lfc_column)
  }
  invisible(TRUE)
}

validate_numeric_column <- function(res_df, column, label = "numeric column") {
  if (!column %in% colnames(res_df)) {
    stop(label, " not found in DESeq2 result: ", column)
  }
  if (!is.numeric(res_df[[column]])) {
    stop(label, " must be numeric: ", column)
  }
  invisible(TRUE)
}

mark_deg_by_threshold <- function(res_df, pvalue_thresh, log2fc_thresh, pvalue_column = "padj", lfc_column = "log2FoldChange") {
  validate_pvalue_column(res_df, pvalue_column)
  validate_lfc_column(res_df, lfc_column)
  out <- res_df
  pvals <- out[[pvalue_column]]
  lfcs <- out[[lfc_column]]
  out$significance <- ifelse(
    !is.na(pvals) & !is.na(lfcs) & pvals < pvalue_thresh & abs(lfcs) > log2fc_thresh,
    ifelse(lfcs > 0, "Up", "Down"),
    "Not_Sig"
  )
  out
}

build_deg_threshold_sets <- function(res_list, threshold_grid, pvalue_column = "padj", lfc_column = "log2FoldChange") {
  deg_by_threshold <- list()
  for (i in seq_len(nrow(threshold_grid))) {
    th <- threshold_grid[i, ]
    deg_by_threshold[[th$name]] <- lapply(
      res_list,
      mark_deg_by_threshold,
      pvalue_thresh = th$p_cutoff,
      log2fc_thresh = th$log2fc,
      pvalue_column = pvalue_column,
      lfc_column = lfc_column
    )
  }
  deg_by_threshold
}

summarize_deg_thresholds <- function(deg_by_threshold, threshold_grid = NULL, pvalue_column = NULL, lfc_column = NULL) {
  do.call(rbind, lapply(names(deg_by_threshold), function(th_name) {
    th <- NULL
    if (!is.null(threshold_grid)) {
      th <- threshold_grid[threshold_grid$name == th_name, , drop = FALSE]
      if (nrow(th) == 0) th <- NULL
    }
    do.call(rbind, lapply(names(deg_by_threshold[[th_name]]), function(comp_name) {
      res <- deg_by_threshold[[th_name]][[comp_name]]
      row <- data.frame(
        Threshold = th_name,
        Comparison = comp_name,
        UP = sum(res$significance == "Up", na.rm = TRUE),
        DOWN = sum(res$significance == "Down", na.rm = TRUE),
        Total = sum(res$significance != "Not_Sig", na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      if (!is.null(th)) {
        row$pvalue_column <- pvalue_column
        row$p_cutoff <- th$p_cutoff
        row$lfc_column <- lfc_column
        row$log2fc_cutoff <- th$log2fc
        row$definition <- paste0(pvalue_column, " < ", th$p_cutoff, " & abs(", lfc_column, ") > ", th$log2fc)
      }
      row
    }))
  }))
}

get_threshold_result <- function(deg_by_threshold, threshold_name) {
  if (!threshold_name %in% names(deg_by_threshold)) {
    stop("Threshold not found: ", threshold_name)
  }
  deg_by_threshold[[threshold_name]]
}

write_all_gene_deg_results <- function(res_list, outdir = "1-DEG") {
  all_gene_dir <- file.path(outdir, "all_genes")
  dir.create(all_gene_dir, showWarnings = FALSE, recursive = TRUE)
  for (comp_name in names(res_list)) {
    res <- res_list[[comp_name]]
    preferred_cols <- c(
      "gene_name", "baseMean",
      "log2FoldChange_raw", "lfcSE_raw",
      "log2FoldChange_shrunken", "lfcSE_shrunken",
      "log2FoldChange", "lfcSE",
      "stat", "pvalue", "padj"
    )
    ordered_cols <- c(intersect(preferred_cols, colnames(res)), setdiff(colnames(res), preferred_cols))
    utils::write.csv(
      res[, ordered_cols, drop = FALSE],
      file.path(all_gene_dir, paste0("DESeq2_all_genes_", comp_name, ".csv")),
      row.names = FALSE
    )
  }
  invisible(all_gene_dir)
}

write_deg_threshold_outputs <- function(deg_by_threshold, outdir = "1-DEG", threshold_grid = NULL, pvalue_column = NULL, lfc_column = NULL) {
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
  summary_df <- summarize_deg_thresholds(deg_by_threshold, threshold_grid, pvalue_column, lfc_column)
  utils::write.csv(summary_df, file.path(outdir, "DEG_threshold_summary.csv"), row.names = FALSE)
  summary_df
}

summarize_lfc_strategies <- function(res_list, threshold_grid, pvalue_column = "padj",
                                     lfc_columns = c("log2FoldChange_raw", "log2FoldChange_shrunken")) {
  rows <- list()
  for (comp_name in names(res_list)) {
    res <- res_list[[comp_name]]
    validate_pvalue_column(res, pvalue_column)
    for (lfc_column in intersect(lfc_columns, colnames(res))) {
      validate_lfc_column(res, lfc_column)
      for (i in seq_len(nrow(threshold_grid))) {
        th <- threshold_grid[i, ]
        marked <- mark_deg_by_threshold(
          res,
          pvalue_thresh = th$p_cutoff,
          log2fc_thresh = th$log2fc,
          pvalue_column = pvalue_column,
          lfc_column = lfc_column
        )
        rows[[length(rows) + 1]] <- data.frame(
          Threshold = th$name,
          Comparison = comp_name,
          pvalue_column = pvalue_column,
          p_cutoff = th$p_cutoff,
          lfc_column = lfc_column,
          lfc_strategy = sub("^log2FoldChange_", "", lfc_column),
          log2fc_cutoff = th$log2fc,
          UP = sum(marked$significance == "Up", na.rm = TRUE),
          DOWN = sum(marked$significance == "Down", na.rm = TRUE),
          Total = sum(marked$significance != "Not_Sig", na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  dplyr::bind_rows(rows)
}

write_lfc_strategy_summary <- function(res_list, threshold_grid, outdir = "1-DEG",
                                       pvalue_column = "padj",
                                       lfc_columns = c("log2FoldChange_raw", "log2FoldChange_shrunken")) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  summary_df <- summarize_lfc_strategies(res_list, threshold_grid, pvalue_column, lfc_columns)
  utils::write.csv(summary_df, file.path(outdir, "DEG_lfc_strategy_summary.csv"), row.names = FALSE)
  summary_df
}

summarize_deg_diagnostics <- function(res_list, threshold_grid, pvalue_columns = c("pvalue", "padj"), lfc_columns = c("log2FoldChange_raw", "log2FoldChange_shrunken")) {
  rows <- list()
  for (comp_name in names(res_list)) {
    res <- res_list[[comp_name]]
    row <- data.frame(
      Comparison = comp_name,
      tested_genes = nrow(res),
      stringsAsFactors = FALSE
    )

    for (pcol in intersect(pvalue_columns, colnames(res))) {
      pvals <- res[[pcol]]
      row[[paste0(pcol, "_non_na")]] <- sum(!is.na(pvals))
      row[[paste0(pcol, "_lt_0.05")]] <- sum(!is.na(pvals) & pvals < 0.05)
    }

    for (lfc_col in intersect(lfc_columns, colnames(res))) {
      lfcs <- res[[lfc_col]]
      lfc_label <- sub("^log2FoldChange_", "", lfc_col)
      for (cutoff in c(0.5, 1, 1.5, 2)) {
        cutoff_label <- gsub("\\.", "_", as.character(cutoff))
        row[[paste0("abs_", lfc_label, "_lfc_gt_", cutoff_label)]] <- sum(!is.na(lfcs) & abs(lfcs) > cutoff)
      }
    }

    for (i in seq_len(nrow(threshold_grid))) {
      th <- threshold_grid[i, ]
      for (pcol in intersect(pvalue_columns, colnames(res))) {
        pvals <- res[[pcol]]
        for (lfc_col in intersect(lfc_columns, colnames(res))) {
          lfcs <- res[[lfc_col]]
          lfc_label <- sub("^log2FoldChange_", "", lfc_col)
          prefix <- paste(th$name, pcol, lfc_label, sep = "_")
          sig <- !is.na(pvals) & !is.na(lfcs) & pvals < th$p_cutoff & abs(lfcs) > th$log2fc
          row[[paste0(prefix, "_total")]] <- sum(sig)
          row[[paste0(prefix, "_up")]] <- sum(sig & lfcs > 0)
          row[[paste0(prefix, "_down")]] <- sum(sig & lfcs < 0)
        }
      }
    }

    rows[[length(rows) + 1]] <- row
  }
  dplyr::bind_rows(rows)
}

write_deg_diagnostic_summary <- function(res_list, threshold_grid, outdir = "1-DEG", pvalue_columns = c("pvalue", "padj"), lfc_columns = c("log2FoldChange_raw", "log2FoldChange_shrunken")) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  summary_df <- summarize_deg_diagnostics(res_list, threshold_grid, pvalue_columns, lfc_columns)
  utils::write.csv(summary_df, file.path(outdir, "DEG_diagnostic_summary.csv"), row.names = FALSE)
  summary_df
}

genes_for_enrichment <- function(res_df, pvalue_thresh, log2fc_thresh, pvalue_column = "padj", lfc_column = "log2FoldChange") {
  validate_pvalue_column(res_df, pvalue_column)
  validate_lfc_column(res_df, lfc_column)
  sig <- res_df |>
    dplyr::filter(!is.na(.data[[pvalue_column]]), !is.na(.data[[lfc_column]]), .data[[pvalue_column]] < pvalue_thresh, abs(.data[[lfc_column]]) > log2fc_thresh) |>
    dplyr::pull(gene_name)
  up <- res_df |>
    dplyr::filter(!is.na(.data[[pvalue_column]]), !is.na(.data[[lfc_column]]), .data[[pvalue_column]] < pvalue_thresh, .data[[lfc_column]] > log2fc_thresh) |>
    dplyr::pull(gene_name)
  down <- res_df |>
    dplyr::filter(!is.na(.data[[pvalue_column]]), !is.na(.data[[lfc_column]]), .data[[pvalue_column]] < pvalue_thresh, .data[[lfc_column]] < -log2fc_thresh) |>
    dplyr::pull(gene_name)
  list(sig = sig, up = up, down = down)
}

ranked_gene_list <- function(res_df, rank_column = "stat") {
  validate_numeric_column(res_df, rank_column, label = "GSEA rank column")
  gene_list <- res_df[[rank_column]]
  names(gene_list) <- res_df$gene_name
  gene_list <- gene_list[!is.na(gene_list)]
  sort(gene_list, decreasing = TRUE)
}
