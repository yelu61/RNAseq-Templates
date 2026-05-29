# Plotting helpers for bulk RNA-seq templates. All helpers save editable PDFs.

theme_publication <- function(base_size = 12, base_family = "sans") {
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

wrap_term_labels <- function(x, width = 45) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

parse_ratio_numeric <- function(x) {
  vapply(strsplit(as.character(x), "/"), function(z) {
    if (length(z) != 2) return(NA_real_)
    numerator <- suppressWarnings(as.numeric(z[1]))
    denominator <- suppressWarnings(as.numeric(z[2]))
    if (!is.finite(numerator) || !is.finite(denominator) || denominator == 0) return(NA_real_)
    numerator / denominator
  }, numeric(1))
}

zscore_rows <- function(mat, cap = 2) {
  mat <- as.matrix(mat)
  z <- t(scale(t(mat)))
  z[is.na(z)] <- 0
  if (!is.null(cap)) {
    z[z > cap] <- cap
    z[z < -cap] <- -cap
  }
  z
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

plot_sample_qc_pdf <- function(sample_qc, filename, group_colors = NULL) {
  qc_long <- tidyr::pivot_longer(
    sample_qc,
    cols = c("library_size", "detected_genes", "zero_fraction", "median_sample_correlation"),
    names_to = "metric",
    values_to = "value"
  )
  qc_long$sample <- factor(qc_long$sample, levels = sample_qc$sample)

  p <- ggplot2::ggplot(qc_long, ggplot2::aes(x = sample, y = value, fill = group)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = NULL, y = NULL, title = "Sample-level QC Metrics") +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  if (!is.null(group_colors) && "group" %in% colnames(sample_qc)) {
    p <- p + ggplot2::scale_fill_manual(values = group_colors, na.value = "#999999")
  }
  save_pdf_plot(p, filename, width = max(8, 0.55 * nrow(sample_qc)), height = 10)
  p
}

plot_sample_correlation_pdf <- function(sample_qc, filename) {
  sample_cor <- attr(sample_qc, "sample_cor")
  if (is.null(sample_cor)) return(invisible(NULL))
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(filename, width = 7, height = 6)
  pheatmap::pheatmap(
    sample_cor,
    color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(255),
    main = "Sample Spearman Correlation"
  )
  grDevices::dev.off()
  invisible(sample_cor)
}

plot_count_distribution_pdf <- function(count_data, filename, max_points_per_sample = 5000) {
  count_df <- as.data.frame(count_data, check.names = FALSE)
  count_df$gene <- rownames(count_df)
  count_long <- tidyr::pivot_longer(count_df, cols = -gene, names_to = "sample", values_to = "count")
  count_long$log10_count <- log10(count_long$count + 1)
  if (!is.null(max_points_per_sample)) {
    count_long <- count_long |>
      dplyr::group_by(sample) |>
      dplyr::slice_sample(n = max_points_per_sample) |>
      dplyr::ungroup()
  }
  p <- ggplot2::ggplot(count_long, ggplot2::aes(x = sample, y = log10_count)) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.65, fill = "#B7D7E8", color = "black") +
    ggplot2::labs(x = NULL, y = "log10(count + 1)", title = "Raw Count Distribution by Sample") +
    theme_publication(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_pdf_plot(p, filename, width = max(8, 0.55 * ncol(count_data)), height = 5)
  p
}

plot_filter_retention_pdf <- function(retention_df, filename) {
  p <- ggplot2::ggplot(retention_df, ggplot2::aes(x = sample, y = retained_count_fraction)) +
    ggplot2::geom_col(fill = "#7AA6C2", width = 0.75) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "Retained count fraction", title = "Library Size Retained After Gene Filtering") +
    theme_publication(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_pdf_plot(p, filename, width = max(8, 0.55 * nrow(retention_df)), height = 5)
  p
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

plot_lfc_strategy_summary_pdf <- function(lfc_strategy_summary, filename) {
  if (is.null(lfc_strategy_summary) || nrow(lfc_strategy_summary) == 0) return(invisible(NULL))
  strategy_levels <- unique(lfc_strategy_summary$lfc_strategy)
  lfc_strategy_summary$lfc_strategy <- factor(lfc_strategy_summary$lfc_strategy, levels = strategy_levels)
  p <- ggplot2::ggplot(
    lfc_strategy_summary,
    ggplot2::aes(x = Comparison, y = Total, fill = lfc_strategy)
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75), width = 0.68) +
    ggplot2::geom_text(
      ggplot2::aes(label = Total),
      position = ggplot2::position_dodge(width = 0.75),
      vjust = -0.45,
      size = 3.1,
      fontface = "bold"
    ) +
    ggplot2::facet_wrap(~ Threshold, scales = "free_x") +
    ggplot2::scale_fill_manual(values = c("raw" = "#7AA6C2", "shrunken" = "#D9896A"), drop = FALSE) +
    ggplot2::labs(x = NULL, y = "Number of DEGs", fill = "LFC strategy", title = "Raw vs Shrunken LFC DEG Counts") +
    theme_publication(base_size = 10) +
    ggplot2::theme(legend.position = "top", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_pdf_plot(p, filename, width = max(8, 3 * length(unique(lfc_strategy_summary$Threshold))), height = 5)
  p
}

plot_volcano_pdf <- function(res_df, comp_name, pvalue_thresh, log2fc_thresh, filename, pvalue_column = "padj", lfc_column = "log2FoldChange") {
  validate_pvalue_column(res_df, pvalue_column)
  validate_lfc_column(res_df, lfc_column)
  res <- res_df[!is.na(res_df[[pvalue_column]]) & !is.na(res_df[[lfc_column]]), ]
  top_genes <- res |>
    dplyr::filter(.data[[pvalue_column]] < pvalue_thresh, abs(.data[[lfc_column]]) > log2fc_thresh) |>
    dplyr::arrange(.data[[pvalue_column]]) |>
    utils::head(15) |>
    dplyr::pull(gene_name)

  keyvals <- rep("darkgrey", nrow(res))
  names(keyvals) <- rep("Not-significant", nrow(res))
  keyvals[res[[lfc_column]] > log2fc_thresh & res[[pvalue_column]] < pvalue_thresh] <- "#d6604d"
  names(keyvals)[res[[lfc_column]] > log2fc_thresh & res[[pvalue_column]] < pvalue_thresh] <- "Up-regulated"
  keyvals[res[[lfc_column]] < -log2fc_thresh & res[[pvalue_column]] < pvalue_thresh] <- "#4393c3"
  names(keyvals)[res[[lfc_column]] < -log2fc_thresh & res[[pvalue_column]] < pvalue_thresh] <- "Down-regulated"
  max_fc <- max(abs(res[[lfc_column]]), na.rm = TRUE)

  p <- EnhancedVolcano::EnhancedVolcano(
    res,
    lab = res$gene_name,
    x = lfc_column,
    y = pvalue_column,
    selectLab = top_genes,
    max.overlaps = 30,
    drawConnectors = TRUE,
    widthConnectors = 0.75,
    boxedLabels = TRUE,
    xlim = c(-max_fc * 1.1, max_fc * 1.1),
    title = comp_name,
    pCutoff = pvalue_thresh,
    FCcutoff = log2fc_thresh,
    xlab = bquote(~Log[2]~"fold change"),
    ylab = bquote(~-Log[10]~.(pvalue_column)),
    pointSize = 1.5,
    labSize = 3.5,
    colCustom = keyvals,
    colAlpha = 4/5,
    legendPosition = "bottom"
  )
  save_pdf_plot(p, filename, width = 8, height = 10)
  p
}

plot_expression_heatmap_pdf <- function(mat, filename, title = NULL, group = NULL, group_levels = NULL,
                                        group_colors = NULL, scale_rows = TRUE, z_cap = 2,
                                        show_row_names = TRUE, show_column_names = FALSE,
                                        row_font_size = 7, column_font_size = 8,
                                        cluster_rows = TRUE, cluster_columns = TRUE,
                                        column_split = TRUE, width = 8, height = 10,
                                        heatmap_name = "Z-score") {
  mat <- as.matrix(mat)
  if (nrow(mat) == 0 || ncol(mat) == 0) return(invisible(NULL))
  plot_mat <- if (scale_rows) zscore_rows(mat, cap = z_cap) else mat
  col_breaks <- if (scale_rows) c(-z_cap, 0, z_cap) else stats::quantile(plot_mat, c(0.02, 0.5, 0.98), na.rm = TRUE)
  if (length(unique(col_breaks)) < 3) {
    center <- stats::median(plot_mat, na.rm = TRUE)
    span <- max(abs(plot_mat - center), na.rm = TRUE)
    if (!is.finite(span) || span == 0) span <- 1
    col_breaks <- c(center - span, center, center + span)
  }
  col_fun <- circlize::colorRamp2(col_breaks, c("#2166ac", "white", "#b2182b"))

  top_annotation <- NULL
  split_vec <- NULL
  if (!is.null(group)) {
    group <- factor(group, levels = group_levels %||% unique(group))
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Group = group,
      col = if (!is.null(group_colors)) list(Group = group_colors) else NULL,
      annotation_name_side = "left"
    )
    if (isTRUE(column_split)) split_vec <- group
  }

  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(filename, width = width, height = height)
  ht <- ComplexHeatmap::Heatmap(
    plot_mat,
    name = heatmap_name,
    col = col_fun,
    top_annotation = top_annotation,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    column_split = split_vec,
    cluster_column_slices = FALSE,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_gp = grid::gpar(fontsize = row_font_size),
    column_names_gp = grid::gpar(fontsize = column_font_size),
    column_title = title,
    use_raster = nrow(plot_mat) > 500,
    show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
          title = expression(~'Z-score of VST'),
          title_gp = gpar(col = "black", cex = 0.75),
          title_position = "leftcenter-rot"
        )
  )
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grDevices::dev.off()
  invisible(ht)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

format_p_for_label <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p)))
}

pairwise_effect_table <- function(data, value_col, group_col, comparisons = NULL,
                                  facet_col = NULL, method = "t.test", p_adjust_method = "BH") {
  df <- data[!is.na(data[[value_col]]) & !is.na(data[[group_col]]), , drop = FALSE]
  if (is.null(comparisons)) {
    comparisons <- utils::combn(levels(factor(df[[group_col]])), 2, simplify = FALSE)
  }
  facets <- if (is.null(facet_col)) NA_character_ else unique(as.character(df[[facet_col]]))
  rows <- list()
  for (facet in facets) {
    sub_df <- if (is.null(facet_col)) df else df[as.character(df[[facet_col]]) == facet, , drop = FALSE]
    value_range <- range(sub_df[[value_col]], na.rm = TRUE)
    value_span <- diff(value_range)
    if (!is.finite(value_span) || value_span == 0) value_span <- max(abs(value_range), na.rm = TRUE) * 0.1 + 1
    for (i in seq_along(comparisons)) {
      pair <- comparisons[[i]]
      x1 <- sub_df[sub_df[[group_col]] == pair[1], value_col, drop = TRUE]
      x2 <- sub_df[sub_df[[group_col]] == pair[2], value_col, drop = TRUE]
      if (length(stats::na.omit(x1)) < 2 || length(stats::na.omit(x2)) < 2) next
      p_val <- tryCatch({
        if (method == "wilcox.test") stats::wilcox.test(x1, x2)$p.value else stats::t.test(x1, x2)$p.value
      }, error = function(e) NA_real_)
      rows[[length(rows) + 1]] <- data.frame(
        group1 = pair[1],
        group2 = pair[2],
        p = p_val,
        mean1 = mean(x1, na.rm = TRUE),
        mean2 = mean(x2, na.rm = TRUE),
        delta = mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE),
        y.position = max(value_range, na.rm = TRUE) + value_span * (0.12 + 0.12 * (i - 1)),
        stringsAsFactors = FALSE
      )
      if (!is.null(facet_col)) rows[[length(rows)]][[facet_col]] <- facet
    }
  }
  if (length(rows) == 0) return(data.frame())
  out <- dplyr::bind_rows(rows)
  if (!is.null(facet_col) && facet_col %in% colnames(out)) {
    out <- out |>
      dplyr::group_by(.data[[facet_col]]) |>
      dplyr::mutate(p.adj = stats::p.adjust(p, method = p_adjust_method)) |>
      dplyr::ungroup()
  } else {
    out$p.adj <- stats::p.adjust(out$p, method = p_adjust_method)
  }
  out$label <- paste0("p.adj=", format_p_for_label(out$p.adj), "\nDelta=", sprintf("%.2f", out$delta))
  out
}

plot_group_boxplot_pdf <- function(data, value_col, group_col, filename, facet_col = NULL,
                                   comparisons = NULL, method = "t.test", p_adjust_method = "BH",
                                   title = NULL, ylab = NULL, group_colors = NULL,
                                   show_ns = TRUE, width = 7, height = 6) {
  plot_data <- data
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = levels(factor(plot_data[[group_col]])))
  stat_tbl <- pairwise_effect_table(
    plot_data, value_col = value_col, group_col = group_col,
    comparisons = comparisons, facet_col = facet_col,
    method = method, p_adjust_method = p_adjust_method
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.88, linewidth = 0.45) +
    ggplot2::geom_jitter(width = 0.12, size = 2, color = "black", shape = 21, fill = "white", alpha = 0.9) +
    ggplot2::labs(title = title, x = NULL, y = ylab %||% value_col) +
    theme_publication(base_size = 11) +
    ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(face = "bold"))
  if (!is.null(group_colors)) p <- p + ggplot2::scale_fill_manual(values = group_colors)
  if (!is.null(facet_col)) p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)), scales = "free_y")
  if (nrow(stat_tbl) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      stat_tbl,
      label = "label",
      hide.ns = !show_ns,
      tip.length = 0.01,
      bracket.size = 0.35,
      size = 2.8
    ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.28)))
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_group_violin_boxplot_pdf <- function(data, value_col, group_col, filename,
                                          comparisons = NULL, method = "t.test",
                                          p_adjust_method = "BH", title = NULL,
                                          ylab = NULL, group_colors = NULL,
                                          show_ns = FALSE, width = 5.8, height = 6) {
  plot_data <- data[!is.na(data[[value_col]]) & !is.na(data[[group_col]]), , drop = FALSE]
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = levels(factor(data[[group_col]])))
  stat_tbl <- pairwise_effect_table(
    plot_data, value_col = value_col, group_col = group_col,
    comparisons = comparisons, method = method, p_adjust_method = p_adjust_method
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    ggplot2::geom_violin(trim = FALSE, width = 0.88, alpha = 0.34, color = NA) +
    ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white", alpha = 0.82, linewidth = 0.42, color = "#2F2F2F") +
    ggplot2::stat_summary(fun = stats::median, geom = "point", shape = 95, size = 7, color = "#2F2F2F") +
    ggplot2::geom_jitter(width = 0.08, size = 2.2, shape = 21, color = "#222222", stroke = 0.25, alpha = 0.9) +
    ggplot2::labs(title = title, x = NULL, y = ylab %||% value_col) +
    theme_publication(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      panel.grid.major.y = ggplot2::element_line(color = "#E7E7E7", linewidth = 0.25),
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, vjust = 1)
    )
  if (!is.null(group_colors)) p <- p + ggplot2::scale_fill_manual(values = group_colors)
  if (nrow(stat_tbl) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      stat_tbl,
      label = "label",
      hide.ns = !show_ns,
      tip.length = 0.01,
      bracket.size = 0.35,
      size = 3
    ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.07, 0.30)))
  } else {
    p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.07, 0.12)))
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

sem_upper <- function(x) mean(x, na.rm = TRUE) + stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
sem_lower <- function(x) mean(x, na.rm = TRUE) - stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

plot_group_bar_sem_pdf <- function(data, value_col, group_col, filename,
                                   comparisons = NULL, method = "t.test",
                                   title = NULL, ylab = NULL, group_colors = NULL,
                                   show_ns = FALSE, width = 5.5, height = 6,
                                   y_from_zero = TRUE) {
  plot_data <- data[!is.na(data[[value_col]]) & !is.na(data[[group_col]]), , drop = FALSE]
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = levels(factor(data[[group_col]])))
  if (is.null(comparisons)) {
    comparisons <- utils::combn(levels(plot_data[[group_col]]), 2, simplify = FALSE)
  }

  stat_tbl <- tryCatch(
    ggpubr::compare_means(
      stats::as.formula(paste(value_col, "~", group_col)),
      data = plot_data,
      method = method,
      comparisons = comparisons
    ),
    error = function(e) data.frame()
  )
  if (nrow(stat_tbl) > 0) {
    ymax <- max(plot_data[[value_col]], na.rm = TRUE)
    ymin <- min(plot_data[[value_col]], na.rm = TRUE)
    span <- ymax - ymin
    if (!is.finite(span) || span == 0) span <- max(abs(plot_data[[value_col]]), na.rm = TRUE) * 0.1 + 1
    stat_tbl$y.position <- ymax + span * (0.10 + (seq_len(nrow(stat_tbl)) - 1) * 0.12)
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    ggplot2::stat_summary(geom = "bar", fun = mean, width = 0.6, color = "white", alpha = 0.95) +
    ggplot2::stat_summary(geom = "errorbar", fun.min = sem_lower, fun.max = sem_upper, width = 0.2, linewidth = 0.45, color = "black") +
    ggplot2::geom_jitter(width = 0.12, size = 2.1, color = "black", alpha = 0.85, shape = 21, fill = "white") +
    ggplot2::labs(title = title, x = NULL, y = ylab %||% value_col) +
    ggplot2::theme_test(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      axis.text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 15, face = "bold")
    )
  if (!is.null(group_colors)) p <- p + ggplot2::scale_fill_manual(values = group_colors)
  if (nrow(stat_tbl) > 0) {
    p <- p + ggpubr::stat_pvalue_manual(
      stat_tbl,
      label = "p.format",
      y.position = "y.position",
      hide.ns = !show_ns,
      tip.length = 0.01,
      bracket.size = 0.38
    )
  }
  if (isTRUE(y_from_zero)) {
    p <- p + ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.15)))
  } else {
    p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18)))
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

prepare_enrich_df <- function(enrich_result, show_category = 15) {
  if (is.null(enrich_result)) return(data.frame())
  df <- as.data.frame(enrich_result)
  if (nrow(df) == 0) return(df)
  df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
  df <- utils::head(df, show_category)
  df$Description_wrapped <- factor(wrap_term_labels(df$Description), levels = rev(wrap_term_labels(df$Description)))
  if ("GeneRatio" %in% colnames(df)) {
    df$GeneRatioNumeric <- parse_ratio_numeric(df$GeneRatio)
  } else {
    df$GeneRatioNumeric <- df$Count / max(df$Count, na.rm = TRUE)
  }
  if (!"FoldEnrichment" %in% colnames(df)) {
    if (all(c("GeneRatio", "BgRatio") %in% colnames(df))) {
      bg_ratio <- parse_ratio_numeric(df$BgRatio)
      df$FoldEnrichment <- df$GeneRatioNumeric / bg_ratio
    } else {
      df$FoldEnrichment <- df$GeneRatioNumeric
    }
  }
  df$log10_padj <- -log10(pmax(df$p.adjust, .Machine$double.xmin))
  df$log10_pvalue <- -log10(pmax(df$pvalue, .Machine$double.xmin))
  df
}

plot_enrich_dotplot <- function(enrich_result, filename, title, show_category = 15, width = 8, height = 10) {
  p <- enrichplot::dotplot(enrich_result, showCategory = show_category, title = title) +
    ggplot2::scale_y_discrete(labels = function(x) wrap_term_labels(x, width = 48)) +
    theme_publication(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_enrich_barplot_pdf <- function(enrich_result, filename, title, show_category = 20, width = 8.5, height = 6,
                                    ontology = NULL, fill_by = c("pvalue", "padj")) {
  fill_by <- match.arg(fill_by)
  df <- prepare_enrich_df(enrich_result, show_category)
  if (!is.null(ontology) && "ONTOLOGY" %in% colnames(df)) df <- df[df$ONTOLOGY == ontology, , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))
  df <- df[order(df$FoldEnrichment, decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, show_category)
  df$Term <- factor(wrap_term_labels(df$Description, width = 55), levels = rev(wrap_term_labels(df$Description, width = 55)))
  max_x <- max(df$FoldEnrichment, na.rm = TRUE)
  label_x <- max_x * 0.015
  df$TextX <- label_x
  fill_col <- if (fill_by == "pvalue") "log10_pvalue" else "log10_padj"
  fill_name <- if (fill_by == "pvalue") "-log10(P)" else "-log10(adj. P)"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = FoldEnrichment, y = Term, fill = .data[[fill_col]])) +
    ggplot2::geom_col(width = 0.78, alpha = 0.78, color = "white", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(x = TextX, label = wrap_term_labels(Description, width = 58)),
      hjust = 0, size = 4.1, color = "#222222", lineheight = 0.92
    ) +
    ggplot2::scale_fill_distiller(palette = "RdPu", direction = 1, name = fill_name) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = "Fold Enrichment", y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 12, base_family = "Times") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 11, color = "black"),
      plot.title = ggplot2::element_text(size = 15, hjust = 0.5, face = "bold"),
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_enrich_bidirectional_barplot_pdf <- function(up_result, down_result, filename, title,
                                                  show_category = 15, width = 16, height = 10,
                                                  ontology = NULL, fill_by = c("pvalue", "padj")) {
  fill_by <- match.arg(fill_by)
  up_df <- prepare_enrich_df(up_result, show_category)
  down_df <- prepare_enrich_df(down_result, show_category)
  if (!is.null(ontology) && "ONTOLOGY" %in% colnames(up_df)) up_df <- up_df[up_df$ONTOLOGY == ontology, , drop = FALSE]
  if (!is.null(ontology) && "ONTOLOGY" %in% colnames(down_df)) down_df <- down_df[down_df$ONTOLOGY == ontology, , drop = FALSE]
  if (nrow(up_df) > 0 && "FoldEnrichment" %in% colnames(up_df)) {
    up_df <- utils::head(up_df[order(up_df$FoldEnrichment, decreasing = TRUE), , drop = FALSE], show_category)
  }
  if (nrow(down_df) > 0 && "FoldEnrichment" %in% colnames(down_df)) {
    down_df <- utils::head(down_df[order(down_df$FoldEnrichment, decreasing = TRUE), , drop = FALSE], show_category)
  }
  if (nrow(up_df) > 0) up_df$Direction <- "UP"
  if (nrow(down_df) > 0) down_df$Direction <- "DOWN"
  df <- dplyr::bind_rows(up_df, down_df)
  if (nrow(df) == 0) return(invisible(NULL))
  df$SignedFoldEnrichment <- ifelse(df$Direction == "UP", df$FoldEnrichment, -df$FoldEnrichment)
  df <- df[order(df$Direction, abs(df$FoldEnrichment), decreasing = TRUE), , drop = FALSE]
  df$Term <- factor(wrap_term_labels(df$Description, width = 52), levels = unique(wrap_term_labels(df$Description[order(df$SignedFoldEnrichment)], width = 52)))
  max_x <- max(abs(df$SignedFoldEnrichment), na.rm = TRUE)
  df$TextX <- ifelse(df$Direction == "UP", -max_x * 0.025, max_x * 0.025)
  df$TextHjust <- ifelse(df$Direction == "UP", 1, 0)
  fill_col <- if (fill_by == "pvalue") "log10_pvalue" else "log10_padj"
  fill_name <- if (fill_by == "pvalue") "-log10(P)" else "-log10(adj. P)"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = SignedFoldEnrichment, y = Term, fill = .data[[fill_col]])) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "#333333") +
    ggplot2::geom_col(width = 0.78, alpha = 0.78, color = "white", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(x = TextX, label = wrap_term_labels(Description, width = 52), hjust = TextHjust),
      size = 4.1, color = "#222222", lineheight = 0.92
    ) +
    ggplot2::scale_fill_distiller(palette = "PuOr", direction = -1, name = fill_name) +
    ggplot2::scale_x_continuous(
      labels = function(x) abs(x),
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::labs(x = "Fold Enrichment", y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 12, base_family = "Times") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 11, color = "black"),
      plot.title = ggplot2::element_text(size = 15, hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_enrich_network_pdf <- function(enrich_result, prefix, show_category = 8) {
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) return(invisible(NULL))
  try({
    p_cnet <- enrichplot::cnetplot(enrich_result, showCategory = show_category, circular = FALSE, colorEdge = TRUE) +
      ggplot2::ggtitle("Gene-term network") +
      theme_publication(base_size = 9)
    save_pdf_plot(p_cnet, paste0(prefix, "_cnetplot.pdf"), width = 10, height = 8)
  }, silent = TRUE)
  try({
    sim_obj <- enrichplot::pairwise_termsim(enrich_result)
    p_emap <- enrichplot::emapplot(sim_obj, showCategory = show_category) +
      ggplot2::ggtitle("Term similarity network") +
      theme_publication(base_size = 9)
    save_pdf_plot(p_emap, paste0(prefix, "_emapplot.pdf"), width = 9, height = 8)
  }, silent = TRUE)
  invisible(TRUE)
}

plot_enrich_suite_pdf <- function(enrich_result, prefix, title_prefix, show_category = 15) {
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) return(invisible(NULL))
  plot_enrich_dotplot(enrich_result, paste0(prefix, "_dotplot.pdf"), paste(title_prefix, "Dotplot"), show_category = show_category)
  plot_enrich_barplot_pdf(enrich_result, paste0(prefix, "_barplot.pdf"), paste(title_prefix, "Top terms"), show_category = show_category)
  plot_enrich_network_pdf(enrich_result, prefix, show_category = min(8, show_category))
  invisible(TRUE)
}

plot_comparecluster_dotplot_pdf <- function(compare_result, filename, title, show_category = 10,
                                            include_all = TRUE, width = 10, height = 12) {
  if (is.null(compare_result) || nrow(as.data.frame(compare_result)) == 0) return(invisible(NULL))
  p <- enrichplot::dotplot(compare_result, showCategory = show_category, includeAll = include_all, title = title) +
    ggplot2::scale_y_discrete(labels = function(x) wrap_term_labels(x, width = 48)) +
    theme_publication(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8.5))
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_gsea_nes_barplot_pdf <- function(gsea_result, filename, title, show_category = 20, width = 9.5, height = 8) {
  df <- as.data.frame(gsea_result)
  if (nrow(df) == 0 || !"NES" %in% colnames(df)) return(invisible(NULL))
  if (!"p.adjust" %in% colnames(df)) df$p.adjust <- df$pvalue
  df <- df[order(df$p.adjust, -abs(df$NES)), , drop = FALSE]
  df <- utils::head(df, show_category)
  df$Direction <- ifelse(df$NES >= 0, "Activated", "Suppressed")
  df$log10_padj <- -log10(pmax(df$p.adjust, .Machine$double.xmin))
  wrapped_terms <- wrap_term_labels(df$Description, width = 54)
  df$Description_wrapped <- factor(wrapped_terms, levels = wrapped_terms[order(df$NES)])
  max_nes <- max(abs(df$NES), na.rm = TRUE)
  df$TextX <- ifelse(df$NES >= 0, max_nes * 0.035, -max_nes * 0.035)
  df$TextHjust <- ifelse(df$NES >= 0, 0, 1)
  df$PointX <- df$NES + ifelse(df$NES >= 0, max_nes * 0.035, -max_nes * 0.035)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = Description_wrapped, fill = Direction)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "#2F2F2F") +
    ggplot2::geom_col(width = 0.76, alpha = 0.88, color = "white", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(x = TextX, label = wrap_term_labels(Description, width = 54), hjust = TextHjust),
      size = 3.8, color = "#1F1F1F", lineheight = 0.9
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = PointX, size = log10_padj),
      shape = 21, fill = "white", color = "#222222", stroke = 0.25,
      inherit.aes = TRUE
    ) +
    ggplot2::scale_fill_manual(values = c("Activated" = "#C8473E", "Suppressed" = "#3778A8")) +
    ggplot2::scale_size_continuous(name = "-log10(FDR)", range = c(2.2, 5.2)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.14, 0.14))) +
    ggplot2::labs(x = "Normalized enrichment score (NES)", y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 12, base_family = "Times") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 11, color = "black"),
      plot.title = ggplot2::element_text(size = 15, hjust = 0.5, face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 11, face = "bold"),
      legend.text = ggplot2::element_text(size = 10)
    )
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_gsea_suite_pdf <- function(gsea_result, prefix, title_prefix, show_category = 15) {
  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) return(invisible(NULL))
  plot_enrich_dotplot(gsea_result, paste0(prefix, "_dotplot.pdf"), paste(title_prefix, "Dotplot"), show_category = show_category)
  plot_gsea_nes_barplot_pdf(gsea_result, paste0(prefix, "_NES_barplot.pdf"), paste(title_prefix, "NES"), show_category = show_category)
  try({
    p_ridge <- enrichplot::ridgeplot(gsea_result, showCategory = show_category) +
      ggplot2::labs(title = paste(title_prefix, "Ridgeplot")) +
      ggplot2::scale_y_discrete(labels = function(x) wrap_term_labels(x, width = 46)) +
      ggplot2::theme_classic(base_size = 11, base_family = "Times") +
      ggplot2::theme(
        axis.title = ggplot2::element_text(size = 12, face = "bold"),
        axis.text = ggplot2::element_text(color = "black"),
        plot.title = ggplot2::element_text(size = 15, hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
    save_pdf_plot(p_ridge, paste0(prefix, "_ridgeplot.pdf"), width = 9.5, height = 10)
  }, silent = TRUE)
  gsea_df <- as.data.frame(gsea_result)
  if ("p.adjust" %in% colnames(gsea_df)) {
    gsea_df <- gsea_df[order(gsea_df$p.adjust, -abs(gsea_df$NES)), , drop = FALSE]
  }
  top_ids <- head(gsea_df$ID, 3)
  if (length(top_ids) > 0) {
    p_running <- enrichplot::gseaplot2(gsea_result, geneSetID = top_ids, title = paste(title_prefix, "Top terms")) +
      ggplot2::theme_classic(base_size = 11, base_family = "Times") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"))
    save_pdf_plot(p_running, paste0(prefix, "_running.pdf"), width = 9, height = 7)
  }
  top_up <- head(gsea_df$ID[gsea_df$NES > 0], 3)
  top_down <- head(gsea_df$ID[gsea_df$NES < 0], 3)
  if (length(top_up) > 0) {
    p_up <- enrichplot::gseaplot2(gsea_result, geneSetID = top_up, title = paste(title_prefix, "Top activated terms")) +
      ggplot2::theme_classic(base_size = 11, base_family = "Times") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"))
    save_pdf_plot(p_up, paste0(prefix, "_running_top_activated.pdf"), width = 9, height = 7)
  }
  if (length(top_down) > 0) {
    p_down <- enrichplot::gseaplot2(gsea_result, geneSetID = top_down, title = paste(title_prefix, "Top suppressed terms")) +
      ggplot2::theme_classic(base_size = 11, base_family = "Times") +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 14, hjust = 0.5, face = "bold"))
    save_pdf_plot(p_down, paste0(prefix, "_running_top_suppressed.pdf"), width = 9, height = 7)
  }
  invisible(TRUE)
}
