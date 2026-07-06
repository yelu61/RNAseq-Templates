# Tumor microenvironment (TME) deconvolution helpers for bulk RNA-seq templates.

# Reverse log-transformed expression if the input is log-scale.
# IOBR/ESTIMATE/CIBERSORT usually expect non-log TPM-like values.
undo_log_expr <- function(expr, is_log = TRUE, log_base = 2) {
  if (!isTRUE(is_log)) return(expr)
  if (log_base == 2) {
    2^expr - 1
  } else {
    log_base^expr - 1
  }
}

# Run IOBR deconvolution for one or more methods.
# Methods supported by IOBR (typical): "cibersort", "epic", "xcell", "estimate".
run_iobr_deconvolution <- function(expr, methods = c("cibersort", "epic", "xcell", "estimate"),
                                   perm = 1000, arrays = FALSE, id_column = "ID") {
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    stop("Package 'IOBR' is required for IOBR-based deconvolution.")
  }

  results <- list()
  for (method in methods) {
    message("Running IOBR deconvolution: ", method)
    res <- tryCatch({
      IOBR::deconvo_tme(eset = expr, method = method, arrays = arrays, perm = perm)
    }, error = function(e) {
      warning("IOBR method ", method, " failed: ", conditionMessage(e))
      NULL
    })
    if (!is.null(res) && nrow(res) > 0) {
      # IOBR returns the sample ID column with varying names (ID, Samples, etc.)
      id_col <- intersect(c("ID", "Sample", "Samples", "sample"), colnames(res))[1]
      if (is.na(id_col)) id_col <- colnames(res)[1]
      colnames(res)[colnames(res) == id_col] <- id_column
      results[[method]] <- res
    }
  }
  if (length(results) == 0) {
    stop("All IOBR deconvolution methods failed.")
  }
  results
}

# Combine IOBR result tables by the common ID column.
combine_tme_results <- function(result_list, id_column = "ID") {
  if (length(result_list) == 0) return(NULL)
  combined <- result_list[[1]]
  if (length(result_list) > 1) {
    for (i in seq(2, length(result_list))) {
      combined <- dplyr::inner_join(combined, result_list[[i]], by = id_column)
    }
  }
  combined
}

# Convert a wide IOBR result table to long format for ggplot.
melt_tme_results <- function(tme_df, id_column = "ID", group_df = NULL,
                             sample_col = "sample", group_col = "condition") {
  value_cols <- setdiff(colnames(tme_df), id_column)
  long <- tme_df |>
    tidyr::pivot_longer(cols = dplyr::all_of(value_cols), names_to = "cell_type", values_to = "fraction") |>
    dplyr::rename(!!sample_col := .data[[id_column]]) |>
    dplyr::mutate(fraction = as.numeric(.data$fraction))
  if (!is.null(group_df) && all(c(sample_col, group_col) %in% colnames(group_df))) {
    long <- long |>
      dplyr::left_join(group_df[, c(sample_col, group_col), drop = FALSE], by = sample_col)
  }
  long
}

# Plot stacked barplot of cell fraction across samples.
plot_tme_barplot_pdf <- function(long_df, group_col = "condition", sample_col = "sample",
                                  filename, title = "TME Cell Fractions", width = 12, height = 7) {
  if (nrow(long_df) == 0 || !all(c(sample_col, "cell_type", "fraction") %in% colnames(long_df))) {
    return(invisible(NULL))
  }
  n_celltypes <- length(unique(long_df$cell_type))
  fill_colors <- if (n_celltypes <= 12) {
    RColorBrewer::brewer.pal(n_celltypes, "Paired")
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_celltypes)
  }

  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = .data[[sample_col]], y = .data$fraction, fill = .data$cell_type)) +
    ggplot2::geom_col(position = "stack", width = 0.85, color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = setNames(fill_colors, unique(long_df$cell_type))) +
    ggplot2::labs(x = NULL, y = "Fraction / Score", title = title, fill = "Cell type") +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "right"
    )
  if (!is.null(group_col) && group_col %in% colnames(long_df)) {
    p <- p + ggplot2::facet_grid(~ .data[[group_col]], scales = "free_x", space = "free_x")
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

# Plot box/violin plots for each cell type with pairwise statistics.
plot_tme_boxplot_pdf <- function(long_df, group_col = "condition", value_col = "fraction",
                                  filename, title = "TME Cell Fraction by Group",
                                  comparisons = NULL, method = "t.test", show_ns = FALSE,
                                  width = 14, height = 10) {
  if (nrow(long_df) == 0 || !all(c(group_col, value_col, "cell_type") %in% colnames(long_df))) {
    return(invisible(NULL))
  }
  plot_data <- long_df
  plot_data[[group_col]] <- factor(plot_data[[group_col]])

  if (is.null(comparisons)) {
    group_levels <- levels(plot_data[[group_col]])
    if (length(group_levels) >= 2) {
      comparisons <- utils::combn(group_levels, 2, simplify = FALSE)
    }
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
    ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", alpha = 0.85, linewidth = 0.4, color = "#2F2F2F") +
    ggplot2::geom_jitter(width = 0.1, size = 1.8, shape = 21, color = "#222222", stroke = 0.25, alpha = 0.9) +
    ggplot2::facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
    ggplot2::labs(x = NULL, y = value_col, title = title) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, size = 8),
      legend.position = "none",
      strip.text = ggplot2::element_text(size = 9, face = "bold")
    )

  if (length(comparisons) > 0) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = method,
      hide.ns = !show_ns,
      tip.length = 0.015,
      bracket.size = 0.35,
      size = 2.8
    )
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

# Extract ESTIMATE-like scores from an IOBR estimate result and make them long-format.
melt_estimate_scores <- function(estimate_df, id_column = "ID", group_df = NULL,
                                 sample_col = "sample", group_col = "condition") {
  score_cols <- intersect(c("StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity"), colnames(estimate_df))
  if (length(score_cols) == 0) return(NULL)
  long <- estimate_df |>
    dplyr::rename(!!sample_col := .data[[id_column]]) |>
    tidyr::pivot_longer(cols = dplyr::all_of(score_cols), names_to = "score_type", values_to = "score") |>
    dplyr::mutate(score = as.numeric(.data$score))
  if (!is.null(group_df) && all(c(sample_col, group_col) %in% colnames(group_df))) {
    long <- long |>
      dplyr::left_join(group_df[, c(sample_col, group_col), drop = FALSE], by = sample_col)
  }
  long
}

# Plot ESTIMATE-style scores across groups.
plot_estimate_boxplot_pdf <- function(long_df, group_col = "condition", filename,
                                      title = "ESTIMATE Scores by Group", method = "t.test",
                                      show_ns = FALSE, width = 10, height = 6) {
  if (is.null(long_df) || nrow(long_df) == 0) return(invisible(NULL))
  plot_data <- long_df
  plot_data[[group_col]] <- factor(plot_data[[group_col]])
  group_levels <- levels(plot_data[[group_col]])
  comparisons <- if (length(group_levels) >= 2) utils::combn(group_levels, 2, simplify = FALSE) else list()

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[group_col]], y = .data$score, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.12, size = 2, shape = 21, color = "#222222", stroke = 0.25, fill = "white", alpha = 0.9) +
    ggplot2::facet_wrap(~ score_type, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Score", title = title) +
    theme_publication(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      legend.position = "none"
    )
  if (length(comparisons) > 0) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = method,
      hide.ns = !show_ns,
      tip.length = 0.015,
      bracket.size = 0.4,
      size = 3
    )
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}
