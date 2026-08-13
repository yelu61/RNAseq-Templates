# Plotting helpers for bulk RNA-seq templates. All helpers save editable PDFs.
# For color conventions and figure sizing see references/VISUALIZATION_STYLE_GUIDE.md.

# Return TRUE only when a PNG device can actually be opened. capabilities("cairo")
# can report TRUE even when the cairo/X11 shared libraries are broken (common on
# headless or partially-configured systems), so probe by opening a real device.
.has_working_cairo <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    cache <<- isTRUE(tryCatch({
      if (!capabilities("cairo")) stop("no cairo")
      tmp <- tempfile(fileext = ".png")
      grDevices::png(tmp, type = "cairo")
      grDevices::dev.off()
      unlink(tmp)
      TRUE
    }, error = function(e) FALSE, warning = function(w) {
      # A warning while opening the device (e.g. "failed to load cairo DLL")
      # means cairo is not actually usable.
      FALSE
    }))
    cache
  }
})

# Close any graphics devices left open by an earlier plotting call. Some plotting
# paths (notably EnhancedVolcano/ggrepel via ggsave) can leak an open device; a
# later function that opens its own device and draws with ComplexHeatmap then
# fails at dev.off() with "写入失败"/"write failed". Closing strays first keeps
# each helper self-contained and order-independent.
.close_leaked_devices <- function() {
  while (grDevices::dev.cur() > 1) {
    try(grDevices::dev.off(), silent = TRUE)
  }
  invisible(NULL)
}

# Standard palette constants used across templates.
RNAseq_PALETTE <- list(
  up_regulated    = "#d6604d",
  down_regulated  = "#4393c3",
  not_significant = "#999999",
  group_two       = c("#6F6F6F", "#E07B54"),
  group_three     = c("#6F6F6F", "#E07B54", "#6D65A3"),
  heatmap_up      = "#b2182b",
  heatmap_down    = "#2166ac",
  heatmap_mid     = "#f7f7f7",
  gsea_activated  = "#C8473E",
  gsea_suppressed = "#3778A8"
)

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
  if (is.null(plot)) {
    message("Skipping PDF: plot object is NULL for ", filename)
    return(invisible(NULL))
  }
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  # Write to a sibling temporary file first. If a plotting method errors while
  # the device is open (a common failure mode for enrichplot/patchwork objects),
  # no zero-content PDF is promoted to the canonical output path.
  tmp <- tempfile(
    pattern = paste0(".", basename(filename), "."),
    tmpdir = dirname(filename), fileext = ".pdf"
  )
  on.exit(unlink(tmp), add = TRUE)
  ggplot2::ggsave(tmp, plot = plot, width = width, height = height, device = device)
  size <- file.info(tmp)$size
  if (!file.exists(tmp) || is.na(size) || size < 1500) {
    stop("Refusing to save an empty or incomplete PDF: ", filename, call. = FALSE)
  }
  if (!file.rename(tmp, filename)) {
    if (!file.copy(tmp, filename, overwrite = TRUE)) {
      stop("Could not promote validated PDF to: ", filename, call. = FALSE)
    }
  }
  invisible(filename)
}

wrap_term_labels <- function(x, width = 45, max_lines = 2) {
  vapply(x, function(s) {
    s <- trimws(gsub("[[:space:]]+", " ", as.character(s)))
    lines <- strwrap(s, width = width)
    if (length(lines) > max_lines) {
      lines <- lines[seq_len(max_lines)]
      lines[max_lines] <- paste0(sub("[[:punct:][:space:]]+$", "", lines[max_lines]), "…")
    }
    paste(lines, collapse = "\n")
  }, character(1))
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
  diag(sample_cor) <- 1
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

plot_deg_summary_pdf <- function(deg_summary, filename, threshold_col = "Threshold") {
  has_threshold <- threshold_col %in% colnames(deg_summary)
  deg_long <- tidyr::pivot_longer(deg_summary, cols = c("UP", "DOWN"), names_to = "Change", values_to = "Number")
  p <- ggplot2::ggplot(deg_long, ggplot2::aes(x = Comparison, y = Number, fill = Change)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = Number), position = ggplot2::position_dodge(width = 0.7), vjust = -0.5, size = 3.5, fontface = "bold") +
    ggplot2::scale_fill_manual(values = c("UP" = "#d6604d", "DOWN" = "#4393c3")) +
    ggplot2::labs(x = NULL, y = "Number of DEGs", title = "DEG Statistics") +
    theme_publication(base_size = 10) +
    ggplot2::theme(legend.position = "top", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  if (has_threshold) {
    p <- p + ggplot2::facet_wrap(~ .data[[threshold_col]], scales = "free_x")
  }
  save_pdf_plot(p, filename, width = max(8, ifelse(has_threshold, 3 * length(unique(deg_summary[[threshold_col]])), 6)), height = 5)
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
  .close_leaked_devices()
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
    # Rasterization goes through ImageMagick, which needs a working graphics
    # device. Disable it when no usable cairo device exists (e.g. minimal or
    # headless systems) so the heatmap still renders instead of erroring.
    use_raster = nrow(plot_mat) > 500 && .has_working_cairo(),
    show_heatmap_legend = TRUE,
        heatmap_legend_param = list(
          title = heatmap_name,
          title_gp = grid::gpar(col = "black", cex = 0.75),
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

# Compact significance stars for an adjusted p-value vector (ns/*/**/***). Used
# by the group plot helpers when `label_style = "stars"` is requested, which is
# easier to read than the full "p.adj=...\nDelta=..." label on small or
# many-comparison facets.
format_p_stars <- function(p) {
  ifelse(is.na(p), "ns",
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
      ifelse(p < 0.05, "*", "ns"))))
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

# Validate a caller-supplied pairwise table and add plot positions when absent.
# This lets figures reuse the exact multiplicity correction reported in a result
# table (for example, one global BH pass across a pathway panel) instead of
# silently recomputing a different per-facet correction.
.prepare_plot_stat_table <- function(stat_table, data, value_col, facet_col = NULL,
                                     label_style = c("full", "stars")) {
  label_style <- match.arg(label_style)
  stat_tbl <- as.data.frame(stat_table)
  if (nrow(stat_tbl) == 0) return(stat_tbl)
  required <- c("group1", "group2")
  if (!is.null(facet_col)) required <- c(required, facet_col)
  missing_cols <- setdiff(required, colnames(stat_tbl))
  if (length(missing_cols) > 0) {
    stop("`stat_table` is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  p_col <- if ("p.adj" %in% colnames(stat_tbl)) "p.adj" else if ("padj" %in% colnames(stat_tbl)) "padj" else NULL
  if (is.null(p_col) && !"label" %in% colnames(stat_tbl)) {
    stop("`stat_table` must contain `p.adj`, `padj`, or a preformatted `label` column.")
  }
  if (!is.null(p_col)) {
    if (label_style == "stars") {
      stat_tbl$label <- format_p_stars(stat_tbl[[p_col]])
    } else if (!"label" %in% colnames(stat_tbl)) {
      delta_text <- if ("delta" %in% colnames(stat_tbl)) {
        paste0("\nDelta=", sprintf("%.2f", stat_tbl$delta))
      } else ""
      stat_tbl$label <- paste0("p.adj=", format_p_for_label(stat_tbl[[p_col]]), delta_text)
    }
  }
  if (!"y.position" %in% colnames(stat_tbl)) {
    stat_tbl$y.position <- NA_real_
    facets <- if (is.null(facet_col)) NA_character_ else unique(as.character(stat_tbl[[facet_col]]))
    for (facet in facets) {
      stat_idx <- if (is.null(facet_col)) seq_len(nrow(stat_tbl)) else which(as.character(stat_tbl[[facet_col]]) == facet)
      plot_values <- if (is.null(facet_col)) data[[value_col]] else data[as.character(data[[facet_col]]) == facet, value_col]
      plot_values <- plot_values[is.finite(plot_values)]
      if (length(plot_values) == 0) next
      value_range <- range(plot_values)
      value_span <- diff(value_range)
      if (!is.finite(value_span) || value_span == 0) value_span <- max(abs(value_range)) * 0.1 + 1
      stat_tbl$y.position[stat_idx] <- max(value_range) + value_span * (0.12 + 0.12 * (seq_along(stat_idx) - 1))
    }
  }
  stat_tbl
}

plot_group_boxplot_pdf <- function(data, value_col, group_col, filename, facet_col = NULL,
                                   comparisons = NULL, method = "t.test", p_adjust_method = "BH",
                                   title = NULL, ylab = NULL, group_colors = NULL,
                                   show_ns = TRUE, width = 7, height = 6,
                                   label_style = c("full", "stars"),
                                   stat_table = NULL) {
  label_style <- match.arg(label_style)
  plot_data <- data
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = levels(factor(plot_data[[group_col]])))
  stat_tbl <- if (is.null(stat_table)) {
    pairwise_effect_table(
      plot_data, value_col = value_col, group_col = group_col,
      comparisons = comparisons, facet_col = facet_col,
      method = method, p_adjust_method = p_adjust_method
    )
  } else {
    .prepare_plot_stat_table(stat_table, plot_data, value_col, facet_col, label_style)
  }
  if (is.null(stat_table) && label_style == "stars" && nrow(stat_tbl) > 0) {
    stat_tbl$label <- format_p_stars(stat_tbl$p.adj)
  }

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
                                          show_ns = FALSE, width = 5.8, height = 6,
                                          label_style = c("full", "stars"),
                                          stat_table = NULL) {
  label_style <- match.arg(label_style)
  plot_data <- data[!is.na(data[[value_col]]) & !is.na(data[[group_col]]), , drop = FALSE]
  plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = levels(factor(data[[group_col]])))
  stat_tbl <- if (is.null(stat_table)) {
    pairwise_effect_table(
      plot_data, value_col = value_col, group_col = group_col,
      comparisons = comparisons, method = method, p_adjust_method = p_adjust_method
    )
  } else {
    .prepare_plot_stat_table(stat_table, plot_data, value_col, NULL, label_style)
  }
  if (is.null(stat_table) && label_style == "stars" && nrow(stat_tbl) > 0) {
    stat_tbl$label <- format_p_stars(stat_tbl$p.adj)
  }

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
                                   p_adjust_method = "BH",
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
      comparisons = comparisons,
      p.adjust.method = p_adjust_method
    ),
    error = function(e) data.frame()
  )
  if (nrow(stat_tbl) > 0) {
    if (!"p.adj" %in% colnames(stat_tbl)) {
      stat_tbl$p.adj <- stats::p.adjust(stat_tbl$p, method = p_adjust_method)
    }
    stat_tbl$p.adj.format <- format_p_for_label(stat_tbl$p.adj)
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
      label = "p.adj.format",
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
  df <- prepare_enrich_df(enrich_result, show_category)
  if (nrow(df) == 0) {
    message("Skipping enrichment dotplot: no significant terms for ", title)
    return(invisible(NULL))
  }
  labels <- wrap_term_labels(df$Description, width = 42, max_lines = 2)
  df$Term <- factor(labels, levels = rev(labels))
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = GeneRatioNumeric, y = Term, size = Count, color = log10_padj)
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_viridis_c(option = "magma", direction = -1, name = "-log10(FDR)") +
    ggplot2::scale_size_continuous(name = "Gene count", range = c(2.5, 7)) +
    ggplot2::labs(x = "Gene ratio", y = NULL, title = title) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9, lineheight = 0.9),
      panel.grid.major.y = ggplot2::element_line(color = "#ECECEC", linewidth = 0.3),
      legend.position = "right"
    )
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_enrich_barplot_pdf <- function(enrich_result, filename, title, show_category = 20, width = 8.5, height = 6,
                                    ontology = NULL, fill_by = c("pvalue", "padj")) {
  fill_by <- match.arg(fill_by)
  df <- prepare_enrich_df(enrich_result, show_category)
  if (!is.null(ontology) && "ONTOLOGY" %in% colnames(df)) df <- df[df$ONTOLOGY == ontology, , drop = FALSE]
  if (nrow(df) == 0) {
    message("Skipping barplot: no significant enrichment terms for ", title)
    return(invisible(NULL))
  }
  df <- df[order(df$FoldEnrichment, decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, show_category)
  labels <- wrap_term_labels(df$Description, width = 42, max_lines = 2)
  df$Term <- factor(labels, levels = rev(labels))
  fill_col <- if (fill_by == "pvalue") "log10_pvalue" else "log10_padj"
  fill_name <- if (fill_by == "pvalue") "-log10(P)" else "-log10(adj. P)"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = FoldEnrichment, y = Term, fill = .data[[fill_col]])) +
    ggplot2::geom_col(width = 0.78, alpha = 0.78, color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_distiller(palette = "RdPu", direction = 1, name = fill_name) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = "Fold Enrichment", y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 12, base_family = "Times") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9, color = "black", lineheight = 0.9),
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
  # Saved CSVs can infer unused columns such as geneID with different types
  # when one direction is empty/sparse. Bind only the fields used by this plot
  # so a presentation-only refresh does not depend on irrelevant CSV typing.
  plot_cols <- c("Description", "FoldEnrichment", "log10_pvalue", "log10_padj", "Direction")
  if (nrow(up_df) > 0) up_df <- up_df[, plot_cols, drop = FALSE]
  if (nrow(down_df) > 0) down_df <- down_df[, plot_cols, drop = FALSE]
  df <- dplyr::bind_rows(up_df, down_df)
  if (nrow(df) == 0) {
    message("Skipping bidirectional barplot: no significant enrichment terms for ", title)
    return(invisible(NULL))
  }
  df$SignedFoldEnrichment <- ifelse(df$Direction == "UP", df$FoldEnrichment, -df$FoldEnrichment)
  df <- df[order(df$Direction, abs(df$FoldEnrichment), decreasing = TRUE), , drop = FALSE]
  labels <- wrap_term_labels(df$Description, width = 42, max_lines = 2)
  df$Term <- factor(labels, levels = unique(labels[order(df$SignedFoldEnrichment)]))
  fill_col <- if (fill_by == "pvalue") "log10_pvalue" else "log10_padj"
  fill_name <- if (fill_by == "pvalue") "-log10(P)" else "-log10(adj. P)"

  p <- ggplot2::ggplot(df, ggplot2::aes(x = SignedFoldEnrichment, y = Term, fill = .data[[fill_col]])) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "#333333") +
    ggplot2::geom_col(width = 0.78, alpha = 0.78, color = "white", linewidth = 0.25) +
    ggplot2::scale_fill_distiller(palette = "PuOr", direction = -1, name = fill_name) +
    ggplot2::scale_x_continuous(
      labels = function(x) abs(x),
      expand = ggplot2::expansion(mult = c(0.08, 0.08))
    ) +
    ggplot2::labs(x = "Fold Enrichment", y = NULL, title = title) +
    ggplot2::theme_classic(base_size = 12, base_family = "Times") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9, color = "black", lineheight = 0.9),
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
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
    message("Skipping network plot: no significant enrichment terms for ", prefix)
    return(invisible(NULL))
  }
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
  if (is.null(enrich_result) || nrow(as.data.frame(enrich_result)) == 0) {
    message("Skipping enrichment plot suite: no significant terms for ", title_prefix)
    return(invisible(NULL))
  }
  plot_enrich_dotplot(enrich_result, paste0(prefix, "_dotplot.pdf"), paste(title_prefix, "Dotplot"), show_category = show_category)
  plot_enrich_barplot_pdf(enrich_result, paste0(prefix, "_barplot.pdf"), paste(title_prefix, "Top terms"), show_category = show_category)
  plot_enrich_network_pdf(enrich_result, prefix, show_category = min(8, show_category))
  invisible(TRUE)
}

plot_comparecluster_dotplot_pdf <- function(compare_result, filename, title, show_category = 10,
                                            include_all = TRUE, width = 10, height = 12) {
  if (is.null(compare_result) || nrow(as.data.frame(compare_result)) == 0) {
    message("Skipping compareCluster dotplot: no significant terms for ", title)
    return(invisible(NULL))
  }
  p <- enrichplot::dotplot(compare_result, showCategory = show_category, includeAll = include_all, title = title) +
    ggplot2::scale_y_discrete(labels = function(x) wrap_term_labels(x, width = 48)) +
    theme_publication(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8.5))
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

plot_gsea_nes_barplot_pdf <- function(gsea_result, filename, title, show_category = 20, width = 9.5, height = 8) {
  df <- as.data.frame(gsea_result)
  if (nrow(df) == 0 || !"NES" %in% colnames(df)) {
    message("Skipping GSEA NES barplot: no significant terms for ", title)
    return(invisible(NULL))
  }
  if (!"p.adjust" %in% colnames(df)) df$p.adjust <- df$pvalue
  df <- df[order(df$p.adjust, -abs(df$NES)), , drop = FALSE]
  df <- utils::head(df, show_category)
  # KEGG results can carry NA or duplicated Descriptions (the online KEGG map
  # supplies IDs without names). Fill missing labels from ID and deduplicate so
  # the factor levels below are unique; otherwise factor() errors with
  # "factor level is duplicated".
  desc <- as.character(df$Description)
  desc[is.na(desc) | desc == ""] <- as.character(df$ID)[is.na(desc) | desc == ""]
  df$Description <- make.unique(desc, sep = " ")
  df$Direction <- ifelse(df$NES >= 0, "Activated", "Suppressed")
  df$log10_padj <- -log10(pmax(df$p.adjust, .Machine$double.xmin))
  wrapped_terms <- wrap_term_labels(df$Description, width = 44, max_lines = 2)
  df$Description_wrapped <- factor(wrapped_terms, levels = wrapped_terms[order(df$NES)])
  max_nes <- max(abs(df$NES), na.rm = TRUE)
  df$PointX <- df$NES + ifelse(df$NES >= 0, max_nes * 0.035, -max_nes * 0.035)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = Description_wrapped, fill = Direction)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, color = "#2F2F2F") +
    ggplot2::geom_col(width = 0.76, alpha = 0.88, color = "white", linewidth = 0.25) +
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
      axis.text.y = ggplot2::element_text(size = 9, color = "black", lineheight = 0.9),
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
  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) {
    message("Skipping GSEA plot suite: no significant terms for ", title_prefix)
    return(invisible(NULL))
  }
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
  message(
    "GSEA overview saved for ", title_prefix,
    "; create one running curve per explicitly selected term with ",
    "plot_gsea_term_figures_from_df()."
  )
  invisible(TRUE)
}

# Publication-quality single-term GSEA running-enrichment figure.
# Refines enrichplot::gseaplot2 output with NES-aware colors, clean title/subtitle,
# and panel-level styling suitable for manuscripts.
plot_gsea_term_figure_pdf <- function(gsea_result, gene_set_id,
                                       filename = NULL,
                                       title = NULL,
                                       subtitle = NULL,
                                       contrast_label = NULL,
                                       width = 10, height = 7.5,
                                       panel_heights = c(1.5, 0.25, 0.4),
                                       pvalue_table = FALSE,
                                       base_family = "Times") {
  if (is.null(gsea_result) || nrow(as.data.frame(gsea_result)) == 0) {
    message("Skipping single-term GSEA figure: no significant terms")
    return(invisible(NULL))
  }
  gsea_df <- as.data.frame(gsea_result)
  if (is.character(gene_set_id)) {
    term_stats <- gsea_df[gsea_df$ID == gene_set_id, , drop = FALSE]
  } else if (is.numeric(gene_set_id)) {
    term_stats <- gsea_df[gene_set_id, , drop = FALSE]
  } else {
    stop("gene_set_id must be a term ID string or row index.")
  }
  if (nrow(term_stats) == 0) {
    warning("gene_set_id not found in gsea_result: ", gene_set_id)
    return(invisible(NULL))
  }

  nes <- term_stats$NES[[1]]
  pval <- term_stats$pvalue[[1]]
  padj <- term_stats$p.adjust[[1]]
  term_desc <- term_stats$Description[[1]]

  # Synthetic or sparse data can yield NA NES/p-values for a term (e.g. unbalanced
  # gene-level statistics). Such terms can't be plotted meaningfully; skip instead
  # of erroring on `if (nes >= 0)` or downstream gseaplot2 calls.
  if (is.na(nes) || length(nes) == 0) {
    message("Skipping single-term GSEA figure (NA NES) for: ", term_desc)
    return(invisible(NULL))
  }

  fmt_p <- function(x) {
    x <- as.numeric(x)
    ifelse(is.na(x), "NA", ifelse(x < 0.001, sprintf("%.2e", x), sprintf("%.4f", x)))
  }
  stats_text <- sprintf("NES = %.3f | p = %s | adj p = %s", nes, fmt_p(pval), fmt_p(padj))

  colors_pos <- list(es = "#b2182b", heatmap_low = "#f4a582", heatmap_high = "#ca0020")
  colors_neg <- list(es = "#2166ac", heatmap_low = "#92c5de", heatmap_high = "#0571b0")
  heatmap_mid <- "#f7f7f7"
  cols <- if (nes >= 0) colors_pos else colors_neg

  title_text <- if (!is.null(title)) title else stringr::str_wrap(term_desc, width = 80)
  if (is.null(subtitle)) {
    subtitle_parts <- stats_text
    if (!is.null(contrast_label)) subtitle_parts <- paste0(contrast_label, "\n", subtitle_parts)
    subtitle_text <- subtitle_parts
  } else {
    subtitle_text <- subtitle
  }

  p <- enrichplot::gseaplot2(gsea_result, geneSetID = gene_set_id, pvalue_table = pvalue_table, color = cols$es, title = "")

  p <- lapply(p, function(x) {
    x + ggplot2::theme(
      text = ggplot2::element_text(family = base_family),
      axis.text = ggplot2::element_text(size = 10, color = "grey20", family = base_family),
      axis.title = ggplot2::element_text(size = 11, color = "grey20", family = base_family),
      axis.ticks = ggplot2::element_line(color = "grey40"),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
  })

  panel1_yrange <- range(ggplot2::ggplot_build(p[[1]])$data[[1]]$y, na.rm = TRUE)
  panel1_pad <- diff(panel1_yrange) * 0.05
  p[[1]] <- p[[1]] +
    ggplot2::scale_y_continuous(
      limits = c(panel1_yrange[1] - panel1_pad, panel1_yrange[2] + panel1_pad),
      breaks = pretty(c(panel1_yrange[1] - panel1_pad, panel1_yrange[2] + panel1_pad), n = 5),
      expand = c(0, 0)
    ) +
    ggplot2::labs(y = "Running enrichment score") +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 6, 12))

  p[[2]] <- p[[2]] +
    ggplot2::scale_fill_gradient2(low = cols$heatmap_low, mid = heatmap_mid, high = cols$heatmap_high, midpoint = 0) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  p[[3]] <- p[[3]] +
    ggplot2::scale_y_continuous(breaks = c(-10, -5, 0, 5, 10)) +
    ggplot2::coord_cartesian(ylim = c(-10, 10)) +
    ggplot2::labs(y = "Ranking metric", x = "Rank in ordered gene list") +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      axis.title.y = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(6, 20, 12, 12)
    )

  p[[1]] <- p[[1]] +
    ggplot2::labs(title = title_text, subtitle = subtitle_text) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 17, face = "bold", hjust = 0.5, color = "black", margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = 12, hjust = 0.5, color = "grey45", margin = ggplot2::margin(b = 8), lineheight = 1.25)
    )

  attr(p, "params") <- list(ncol = 1, heights = panel_heights)
  class(p) <- c("gglist", "list")

  if (!is.null(filename)) {
    dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
    device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggplot2::ggsave(filename, plot = p, width = width, height = height, device = device, limitsize = FALSE)
  }
  p
}

# Plot publication-grade single-term GSEA figures for selected terms.
# term_df must contain columns: ID (or Description), and optionally a grouping column.
plot_gsea_term_figures_from_df <- function(gsea_result, term_df,
                                           outdir,
                                           contrast_label = NULL,
                                           id_col = "ID",
                                           width = 10, height = 7.5,
                                           prefix = "gseaplot2") {
  if (nrow(term_df) == 0) return(invisible(NULL))
  if (!id_col %in% colnames(term_df) && !"Description" %in% colnames(term_df)) {
    stop("term_df must contain 'ID' or 'Description' column")
  }
  if (!id_col %in% colnames(term_df)) id_col <- "Description"
  gsea_df <- as.data.frame(gsea_result)

  for (i in seq_len(nrow(term_df))) {
    term_id <- term_df[[id_col]][[i]]
    if (id_col == "Description") {
      term_id <- gsea_df$ID[match(term_id, gsea_df$Description)]
    }
    if (is.na(term_id)) next
    term_desc <- gsea_df$Description[match(term_id, gsea_df$ID)]
    safe_name <- gsub("[^a-zA-Z0-9_-]", "_", tolower(term_desc))
    filename <- file.path(outdir, paste0(prefix, "_", safe_name, ".pdf"))
    plot_gsea_term_figure_pdf(
      gsea_result, term_id,
      filename = filename,
      contrast_label = contrast_label,
      width = width, height = height
    )
  }
  invisible(TRUE)
}

# Publication-quality enrichment theme dot-heatmap.
# Inspired by theme-based GO ORA dot-heatmaps used in manuscript figures.

# Default theme dictionary for common cancer / immunology GO-BP themes.
default_enrichment_themes <- function() {
  data.frame(
    figure_group = c(
      rep("Immune & Inflammation", 8),
      rep("Cellular Processes", 6),
      rep("Metabolism & Bioenergetics", 3),
      rep("Development & Tissue", 3),
      rep("Genome Regulation", 3)
    ),
    theme = c(
      "Interferon Response", "Antigen Presentation", "Cytotoxicity / Killing",
      "T/NK Activation", "Cytokine / Chemokine", "Inflammation",
      "Leukocyte Migration / Adhesion", "TNF-NFkappaB Signaling",
      "Apoptosis", "Cell Cycle / Proliferation", "DNA Damage Response",
      "Barrier / Junction", "Focal Adhesion / ECM", "Synapse / Neural",
      "Glycolysis", "OXPHOS", "Lipid Metabolism",
      "Angiogenesis", "Wnt Signaling", "EMT / Migration",
      "Chromatin / Epigenetic", "Transcription / RNA processing", "Ribosome / Translation"
    ),
    keywords = c(
      "interferon|response to interferon-alpha|response to interferon-gamma|interferon-mediated signaling",
      "antigen processing and presentation|mhc class|mhc protein complex|peptide antigen assembly",
      "natural killer cell mediated|cell killing|leukocyte mediated cytotoxicity|cytotoxicity",
      "t cell activation|leukocyte activation|lymphocyte activation|adaptive immune response|natural killer cell activation",
      "cytokine|cytokine-mediated signaling|chemokine|interleukin production",
      "inflammatory response|regulation of inflammatory response|acute inflammatory response",
      "leukocyte cell-cell adhesion|leukocyte migration|regulation of leukocyte migration|monocyte migration|granulocyte migration",
      "tumor necrosis factor|nf-kappab|canonical nf-kappab|i-kappab kinase",
      "apoptotic|intrinsic apoptotic|execution phase of apoptosis|regulation of apoptotic",
      "cell cycle|mitotic|proliferation|dna replication|chromosome segregation",
      "dna damage|dna repair|double-strand break|cellular response to dna damage",
      "tight junction|adherens junction|cell-cell junction|cell junction organization",
      "focal adhesion|extracellular matrix|ecm|cell adhesion",
      "synapse|synaptic|neurotransmitter|axon|dendrite",
      "glycolytic|glycolysis|glucose metabolic|hexose metabolic",
      "oxidative phosphorylation|atp synthesis|aerobic respiration|electron transport chain|respiratory electron transport",
      "lipid metabolic|fatty acid|cholesterol|lipid biosynthetic",
      "angiogenesis|vasculature development|blood vessel development",
      "wnt signaling|canonical wnt|beta-catenin",
      "epithelial to mesenchymal|mesenchymal|cell migration|motility",
      "chromatin organization|heterochromatin|histone|epigenetic|nucleosome",
      "transcription|rna processing|rna splicing|gene expression",
      "cytoplasmic translation|mitochondrial translation|ribosome|translation"
    ),
    stringsAsFactors = FALSE
  )
}

# Shorten GO/KEGG term labels for publication figures.
clean_term_label <- function(x) {
  x <-
    stringr::str_replace_all(x, stringr::regex("positive regulation of", ignore_case = TRUE), "+reg") |>
    stringr::str_replace_all(stringr::regex("negative regulation of", ignore_case = TRUE), "-reg") |>
    stringr::str_replace_all(stringr::regex("regulation of", ignore_case = TRUE), "reg") |>
    stringr::str_replace_all(stringr::regex("antigen processing and presentation", ignore_case = TRUE), "Ag processing & presentation") |>
    stringr::str_replace_all(stringr::regex("natural killer cell mediated", ignore_case = TRUE), "NK-mediated") |>
    stringr::str_replace_all(stringr::regex("leukocyte mediated cytotoxicity", ignore_case = TRUE), "leukocyte cytotoxicity") |>
    stringr::str_replace_all(stringr::regex("cytokine-mediated signaling pathway", ignore_case = TRUE), "cytokine-mediated signaling") |>
    stringr::str_replace_all(stringr::regex("tumor necrosis factor", ignore_case = TRUE), "TNF") |>
    stringr::str_replace_all(stringr::regex("canonical nf-kappab signal transduction", ignore_case = TRUE), "canonical NF-kappaB signaling") |>
    stringr::str_replace_all(stringr::regex("oxidative phosphorylation", ignore_case = TRUE), "OXPHOS") |>
    stringr::str_replace_all(stringr::regex("cytoplasmic translation", ignore_case = TRUE), "cyto translation") |>
    stringr::str_replace_all(stringr::regex("mitochondrial translation", ignore_case = TRUE), "mito translation") |>
    stringr::str_replace_all(stringr::regex("extracellular matrix", ignore_case = TRUE), "ECM") |>
    stringr::str_replace_all(stringr::regex("response to", ignore_case = TRUE), "resp. to") |>
    stringr::str_squish()
  x
}

# Match enrichment terms to a theme dictionary using case-insensitive regex.
match_enrichment_themes <- function(df, theme_defs, description_col = "Description") {
  if (nrow(df) == 0) return(df)
  if (!description_col %in% colnames(df)) {
    stop("Description column not found: ", description_col)
  }
  desc_lower <- tolower(as.character(df[[description_col]]))
  hits <- list()
  for (i in seq_len(nrow(theme_defs))) {
    pattern <- theme_defs$keywords[[i]]
    idx <- which(stringr::str_detect(desc_lower, stringr::regex(pattern, ignore_case = TRUE)))
    if (length(idx) == 0) next
    hit_df <- df[idx, , drop = FALSE]
    hit_df$figure_group <- theme_defs$figure_group[[i]]
    hit_df$theme <- theme_defs$theme[[i]]
    hits[[length(hits) + 1]] <- hit_df
  }
  if (length(hits) == 0) return(df[0, , drop = FALSE])
  dplyr::bind_rows(hits)
}

# Prepare a data frame for theme dot-heatmap: pick top terms per theme by p.adjust and Count.
prepare_theme_dotplot_df <- function(theme_hits,
                                      comparison_col = "comparison",
                                      top_n = 6,
                                      p_value_col = "p.adjust",
                                      count_col = "Count") {
  if (nrow(theme_hits) == 0) return(theme_hits)
  required <- c("figure_group", "theme", "Description", comparison_col, p_value_col, count_col)
  missing <- setdiff(required, colnames(theme_hits))
  if (length(missing) > 0) {
    stop("theme_hits is missing columns: ", paste(missing, collapse = ", "))
  }
  theme_hits$neg_log10_padj <- -log10(pmax(as.numeric(theme_hits[[p_value_col]]), .Machine$double.xmin))

  top_terms <- theme_hits |>
    dplyr::group_by(.data$figure_group, .data$theme, .data$Description) |>
    dplyr::summarise(
      best_padj = min(as.numeric(.data[[p_value_col]]), na.rm = TRUE),
      max_count = max(as.numeric(.data[[count_col]]), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(.data$figure_group, .data$theme, .data$best_padj, dplyr::desc(.data$max_count)) |>
    dplyr::group_by(.data$figure_group, .data$theme) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()

  plot_df <- theme_hits |>
    dplyr::semi_join(top_terms, by = c("figure_group", "theme", "Description")) |>
    dplyr::left_join(
      top_terms |> dplyr::select("figure_group", "theme", "Description", "best_padj"),
      by = c("figure_group", "theme", "Description")
    )

  plot_df$figure_group <- factor(plot_df$figure_group, levels = unique(plot_df$figure_group))
  plot_df$theme <- factor(plot_df$theme, levels = unique(plot_df$theme))
  if (is.character(plot_df[[comparison_col]])) {
    plot_df[[comparison_col]] <- factor(plot_df[[comparison_col]], levels = unique(plot_df[[comparison_col]]))
  }
  plot_df$label <- clean_term_label(plot_df$Description)
  plot_df$label <- forcats::fct_reorder(plot_df$label, plot_df$best_padj, .desc = TRUE)
  plot_df
}

# Draw a publication-quality theme dot-heatmap and save as PDF.
plot_theme_dotheatmap_pdf <- function(plot_df,
                                       filename,
                                       title = "Enrichment Theme Dot-heatmap",
                                       subtitle = NULL,
                                       comparison_col = "comparison",
                                       count_col = "Count",
                                       width = 12,
                                       height = NULL,
                                       max_width = 30,
                                       max_height = 30,
                                       fill_name = expression(-log[10]("adj. P")),
                                       size_name = "Gene count") {
  if (nrow(plot_df) == 0) {
    message("Skipping theme dot-heatmap: no plot data for ", title)
    return(invisible(NULL))
  }
  use_group_in_strip <- dplyr::n_distinct(plot_df$figure_group) > 1
  plot_df$facet_label <- if (use_group_in_strip) {
    paste(plot_df$figure_group, plot_df$theme, sep = ": ")
  } else {
    as.character(plot_df$theme)
  }
  plot_df$facet_label <- factor(plot_df$facet_label, levels = unique(plot_df$facet_label))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[comparison_col]], y = .data$label)) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data[[count_col]], fill = .data$neg_log10_padj),
      shape = 21, color = "grey20", stroke = 0.25, alpha = 0.92
    ) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::scale_fill_gradientn(
      colors = c("#E8EDF1", "#94C5B2", "#3E8E7E", "#B74949"),
      name = fill_name
    ) +
    ggplot2::scale_size_continuous(range = c(1.8, 8.5), name = size_name) +
    ggplot2::facet_grid(
      facet_label ~ ., scales = "free_y", space = "free_y", switch = "y"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.background.y = ggplot2::element_rect(fill = "#464C55", color = NA),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold", color = "white", size = 9.5),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 9, hjust = 0),
      axis.ticks.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 10, color = "grey35"),
      plot.margin = ggplot2::margin(8, 12, 8, 12)
    )

  if (is.null(height)) {
    height <- max(8, 0.32 * dplyr::n_distinct(plot_df$label) + 0.6 * dplyr::n_distinct(plot_df$facet_label))
  }
  if (!is.numeric(width) || length(width) != 1 || !is.finite(width) || width <= 0 ||
      !is.numeric(height) || length(height) != 1 || !is.finite(height) || height <= 0) {
    stop("Theme dot-heatmap dimensions must be finite positive numbers.")
  }
  if (width > max_width || height > max_height) {
    warning(sprintf(
      "Theme dot-heatmap dimensions capped from %.1f x %.1f to %.1f x %.1f inches; reduce `top_n` for a less dense figure.",
      width, height, min(width, max_width), min(height, max_height)
    ))
    width <- min(width, max_width)
    height <- min(height, max_height)
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

# Wrapper: build a theme dot-heatmap directly from enrichment result objects or data frames.
plot_theme_dotheatmap_from_results <- function(result_map,
                                              filename,
                                              title = "Enrichment Theme Dot-heatmap",
                                              subtitle = NULL,
                                              theme_defs = default_enrichment_themes(),
                                              ontology_filter = "BP",
                                              top_n = 6,
                                              comparison_col = "comparison",
                                              width = 12,
                                              height = NULL) {
  df <- build_multi_comparison_enrich_df(
    result_map,
    comparison_col = comparison_col,
    ontology_filter = ontology_filter
  )
  if (nrow(df) == 0) {
    warning("No enrichment terms to plot.")
    return(invisible(NULL))
  }
  themed <- match_enrichment_themes(df, theme_defs)
  if (nrow(themed) == 0) {
    warning("No terms matched the theme dictionary.")
    return(invisible(NULL))
  }
  plot_df <- prepare_theme_dotplot_df(themed, comparison_col = comparison_col, top_n = top_n)
  plot_theme_dotheatmap_pdf(
    plot_df, filename, title = title, subtitle = subtitle,
    comparison_col = comparison_col, width = width, height = height
  )
}
