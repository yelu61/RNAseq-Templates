# Batch effect diagnostics for bulk RNA-seq templates.

# Plot PCA colored by a batch vector and save as PDF.
plot_pca_by_batch_pdf <- function(vsd, batch_vec, filename,
                                    group_colors = NULL,
                                    width = 7, height = 6) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for PCA.")
  }
  batch_vec <- align_to_sample_order(batch_vec, colnames(vsd), "BATCH_VECTOR")
  if (length(batch_vec) != ncol(vsd)) {
    stop("batch_vec length (", length(batch_vec), ") must match vsd samples (", ncol(vsd), ").")
  }

  col_data <- SummarizedExperiment::colData(vsd)
  col_data$batch <- factor(batch_vec)
  SummarizedExperiment::colData(vsd) <- col_data

  pca_data <- DESeq2::plotPCA(vsd, intgroup = "batch", returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  if (is.null(group_colors)) {
    batch_levels <- levels(pca_data$batch)
    group_colors <- make_group_colors(batch_levels)
  }

  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = batch, fill = batch)) +
    ggplot2::geom_point(size = 4, alpha = 0.85) +
    ggplot2::xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggplot2::ggtitle("PCA of Gene Expression Profiles (colored by batch)") +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_fill_manual(values = group_colors) +
    theme_publication(base_size = 8) +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1)

  save_pdf_plot(p, filename, width = width, height = height)
  invisible(p)
}

# Summarize percentage of variance explained by batch for top PCs.
# Returns a data frame with columns PC, R2_batch, R2_condition (if condition supplied).
summarize_pve_by_batch <- function(vsd, batch_vec, condition_vec = NULL, ntop = 500, npcs = 10) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required.")
  }
  batch_vec <- align_to_sample_order(batch_vec, colnames(vsd), "BATCH_VECTOR")
  if (!is.null(condition_vec)) {
    condition_vec <- align_to_sample_order(condition_vec, colnames(vsd), "condition_vec")
  }
  if (length(batch_vec) != ncol(vsd)) {
    stop("batch_vec length must match vsd samples.")
  }

  # Compute PCA on top variable genes (matching DESeq2::plotPCA logic)
  mat <- SummarizedExperiment::assay(vsd)
  rv <- matrixStats::rowVars(mat)
  top_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, nrow(mat)))]
  pca <- stats::prcomp(t(mat[top_genes, , drop = FALSE]))
  pcs <- pca$x[, seq_len(min(npcs, ncol(pca$x))), drop = FALSE]

  get_r2 <- function(pc, factor_vec) {
    if (length(unique(factor_vec)) < 2) return(NA_real_)
    fit <- summary(stats::aov(pc ~ factor_vec))
    fit[[1]][["Sum Sq"]][1] / sum(fit[[1]][["Sum Sq"]])
  }

  rows <- lapply(seq_len(ncol(pcs)), function(i) {
    row <- data.frame(
      PC = paste0("PC", i),
      R2_batch = get_r2(pcs[, i], factor(batch_vec)),
      stringsAsFactors = FALSE
    )
    if (!is.null(condition_vec) && length(condition_vec) == ncol(vsd)) {
      row$R2_condition <- get_r2(pcs[, i], factor(condition_vec))
    }
    row
  })
  do.call(rbind, rows)
}

# Plot batch vs condition PVE as a grouped barplot.
plot_batch_pve_pdf <- function(pve_df, filename,
                                width = 8, height = 5) {
  if (is.null(pve_df) || nrow(pve_df) == 0) return(invisible(NULL))
  value_cols <- grep("^R2_", colnames(pve_df), value = TRUE)
  if (length(value_cols) == 0) return(invisible(NULL))

  df <- tidyr::pivot_longer(
    pve_df,
    cols = dplyr::all_of(value_cols),
    names_to = "source",
    values_to = "R2"
  )
  df$source <- gsub("^R2_", "", df$source)
  df$PC <- factor(df$PC, levels = unique(df$PC))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC, y = R2, fill = source)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = c("batch" = "#E07B54", "condition" = "#6F6F6F")) +
    ggplot2::labs(x = NULL, y = "Variance explained", title = "Variance explained by batch and condition") +
    theme_publication(base_size = 8) +
    ggplot2::theme(legend.position = "top", legend.title = ggplot2::element_blank())
  save_pdf_plot(p, filename, width = width, height = height)
  invisible(p)
}
