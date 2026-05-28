# Plotting helpers for bulk RNA-seq templates. All helpers save editable PDFs.

theme_publication <- function(base_size = 12, base_family = "Arial") {
  ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.text = ggplot2::element_text(color = "black", size = base_size),
      axis.title = ggplot2::element_text(color = "black", size = base_size + 2, face = "bold"),
      plot.title = ggplot2::element_text(color = "black", size = base_size + 4, face = "bold", hjust = 0.5),
      legend.background = ggplot2::element_rect(fill = NA, color = NA),
      legend.key = ggplot2::element_rect(fill = NA, color = NA),
      legend.text = ggplot2::element_text(size = base_size),
      legend.title = ggplot2::element_text(size = base_size, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "#E8E8E8", color = "black"),
      strip.text = ggplot2::element_text(face = "bold", size = base_size)
    )
}

make_group_colors <- function(group_levels) {
  n <- length(group_levels)
  if (n == 2) {
    colors <- c("#6F6F6F", "#E07B54")
  } else if (n == 3) {
    colors <- c("#6F6F6F", "#E07B54", "#6D65A3")
  } else if (n <= 9) {
    colors <- RColorBrewer::brewer.pal(n, "Set1")
  } else {
    colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n)
  }
  stats::setNames(colors, group_levels)
}

save_pdf_plot <- function(plot, filename, width = 7, height = 6) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, device = grDevices::cairo_pdf)
  invisible(filename)
}

plot_pca_pdf <- function(vsd, group_levels, group_colors, filename, intgroup = "condition") {
  pca_data <- DESeq2::plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  pca_data[[intgroup]] <- factor(pca_data[[intgroup]], levels = group_levels)
  percent_var <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = .data[[intgroup]], fill = .data[[intgroup]])) +
    ggplot2::geom_point(size = 4, alpha = 0.85) +
    ggplot2::xlab(paste0("PC1: ", percent_var[1], "% variance")) +
    ggplot2::ylab(paste0("PC2: ", percent_var[2], "% variance")) +
    ggplot2::ggtitle("PCA of Gene Expression Profiles") +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::scale_fill_manual(values = group_colors) +
    theme_publication(base_size = 14) +
    ggplot2::theme(legend.position = "right", aspect.ratio = 1)

  if (all(table(pca_data[[intgroup]]) >= 3)) {
    p <- p + ggplot2::stat_ellipse(type = "norm", level = 0.95, geom = "polygon", alpha = 0.12, linewidth = 0.5)
  }
  save_pdf_plot(p, filename, width = 7, height = 6)
  p
}

plot_sample_distance_pdf <- function(vsd, filename) {
  sample_dists <- stats::dist(t(SummarizedExperiment::assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)
  rownames(sample_dist_matrix) <- colnames(vsd)
  colnames(sample_dist_matrix) <- colnames(vsd)
  colors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  p <- pheatmap::pheatmap(
    sample_dist_matrix,
    clustering_distance_rows = sample_dists,
    clustering_distance_cols = sample_dists,
    col = colors,
    main = "Sample-to-Sample Distance Heatmap",
    display_numbers = TRUE,
    number_format = "%.1f",
    fontsize_number = 8
  )
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(filename, width = 7, height = 6)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grDevices::dev.off()
  invisible(p)
}

plot_deg_summary_pdf <- function(deg_summary, filename) {
  deg_long <- tidyr::pivot_longer(deg_summary, cols = c("UP", "DOWN"), names_to = "Change", values_to = "Number")
  p <- ggplot2::ggplot(deg_long, ggplot2::aes(x = Comparison, y = Number, fill = Change)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = Number), position = ggplot2::position_dodge(width = 0.7), vjust = -0.5, size = 3.5, fontface = "bold") +
    ggplot2::facet_wrap(~ Threshold, scales = "free_x") +
    ggplot2::scale_fill_manual(values = c("UP" = "#d6604d", "DOWN" = "#4393c3")) +
    ggplot2::labs(x = NULL, y = "Number of DEGs", title = "DEG Statistics Across Thresholds") +
    theme_publication(base_size = 10) +
    ggplot2::theme(legend.position = "top", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_pdf_plot(p, filename, width = max(8, 3 * length(unique(deg_summary$Threshold))), height = 5)
  p
}

plot_volcano_pdf <- function(res_df, comp_name, padj_thresh, log2fc_thresh, filename) {
  res <- res_df[!is.na(res_df$padj), ]
  top_genes <- res |>
    dplyr::filter(padj < padj_thresh, abs(log2FoldChange) > log2fc_thresh) |>
    dplyr::arrange(padj) |>
    utils::head(15) |>
    dplyr::pull(gene_name)

  keyvals <- rep("darkgrey", nrow(res))
  names(keyvals) <- rep("Not-significant", nrow(res))
  keyvals[res$log2FoldChange > log2fc_thresh & res$padj < padj_thresh] <- "#d6604d"
  names(keyvals)[res$log2FoldChange > log2fc_thresh & res$padj < padj_thresh] <- "Up-regulated"
  keyvals[res$log2FoldChange < -log2fc_thresh & res$padj < padj_thresh] <- "#4393c3"
  names(keyvals)[res$log2FoldChange < -log2fc_thresh & res$padj < padj_thresh] <- "Down-regulated"
  max_fc <- max(abs(res$log2FoldChange), na.rm = TRUE)

  p <- EnhancedVolcano::EnhancedVolcano(
    res,
    lab = res$gene_name,
    x = "log2FoldChange",
    y = "padj",
    selectLab = top_genes,
    max.overlaps = 30,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    boxedLabels = TRUE,
    xlim = c(-max_fc * 1.1, max_fc * 1.1),
    title = comp_name,
    pCutoff = padj_thresh,
    FCcutoff = log2fc_thresh,
    xlab = bquote(~Log[2]~"fold change"),
    ylab = bquote(~-Log[10]~italic(p-adjusted)),
    pointSize = 1.5,
    labSize = 3.5,
    colCustom = keyvals,
    colAlpha = 4/5,
    legendPosition = "bottom"
  )
  save_pdf_plot(p, filename, width = 8, height = 10)
  p
}

plot_enrich_dotplot <- function(enrich_result, filename, title, show_category = 15, width = 8, height = 10) {
  p <- enrichplot::dotplot(enrich_result, showCategory = show_category, title = title) +
    theme_publication(base_size = 10)
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_gsea_suite_pdf <- function(gsea_result, prefix, title_prefix, show_category = 15) {
  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) return(invisible(NULL))
  plot_enrich_dotplot(gsea_result, paste0(prefix, "_dotplot.pdf"), paste(title_prefix, "Dotplot"), show_category = show_category)
  try({
    p_ridge <- enrichplot::ridgeplot(gsea_result, showCategory = show_category) +
      ggplot2::labs(title = paste(title_prefix, "Ridgeplot")) +
      theme_publication(base_size = 10)
    save_pdf_plot(p_ridge, paste0(prefix, "_ridgeplot.pdf"), width = 9, height = 10)
  }, silent = TRUE)
  top_ids <- head(as.data.frame(gsea_result)$ID, 3)
  if (length(top_ids) > 0) {
    p_running <- enrichplot::gseaplot2(gsea_result, geneSetID = top_ids, title = paste(title_prefix, "Top terms"))
    save_pdf_plot(p_running, paste0(prefix, "_running.pdf"), width = 9, height = 7)
  }
  invisible(TRUE)
}
