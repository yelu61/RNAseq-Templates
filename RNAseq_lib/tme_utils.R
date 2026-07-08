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

# Validate that an expression matrix is suitable for TME deconvolution.
# Checks: numeric, non-negative, rownames present.
# Warns if rownames look like Ensembl IDs rather than gene symbols.
validate_tme_input <- function(expr) {
  if (is.null(expr)) stop("Expression matrix is NULL.")
  if (!is.matrix(expr) && !is.data.frame(expr)) {
    stop("Expression matrix must be a matrix or data.frame.")
  }

  rnames <- rownames(expr)
  if (is.null(rnames) || any(rnames == "" | is.na(rnames))) {
    stop("Expression matrix must have non-empty rownames (gene symbols).")
  }

  expr <- as.data.frame(expr, check.names = FALSE)
  if (!all(vapply(expr, is.numeric, logical(1)))) {
    stop("Expression matrix must contain numeric values only.")
  }
  if (anyNA(expr)) stop("Expression matrix contains missing values.")
  if (any(expr < 0, na.rm = TRUE)) stop("Expression matrix contains negative values.")

  looks_ensembl <- any(grepl("^(ENS|ENSMUS|ENSMUST|ENSG)", rnames, ignore.case = TRUE))
  if (looks_ensembl) {
    warning("Rownames look like Ensembl IDs. TME deconvolution methods expect gene symbols (e.g. HGNC or MGI symbols).")
  }

  invisible(TRUE)
}

# Deduplicate an expression matrix by gene symbol, keeping the row with the
# highest average expression. Returns a matrix with unique rownames.
deduplicate_expression_by_symbol <- function(expr, symbol_col = NULL) {
  expr <- as.data.frame(expr, check.names = FALSE)
  if (!is.null(symbol_col)) {
    if (!symbol_col %in% colnames(expr)) {
      stop("symbol_col not found in expression data: ", symbol_col)
    }
    symbols <- as.character(expr[[symbol_col]])
    expr <- expr[, setdiff(colnames(expr), symbol_col), drop = FALSE]
  } else {
    symbols <- rownames(expr)
  }
  if (is.null(symbols) || any(symbols == "" | is.na(symbols))) {
    stop("Gene symbols are missing or empty; cannot deduplicate.")
  }

  mean_expr <- rowMeans(expr, na.rm = TRUE)
  ord <- order(mean_expr, decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]
  symbols <- symbols[ord]
  keep <- !duplicated(symbols)
  expr <- expr[keep, , drop = FALSE]
  rownames(expr) <- symbols[keep]
  expr
}

# Convert a mouse gene-symbol expression matrix to human orthologs using biomaRt.
# The input should have MGI symbols as rownames. Output has HGNC symbols as rownames.
convert_mouse_symbols_to_human <- function(expr,
                                           mouse_attr = "mgi_symbol",
                                           human_attr = "hgnc_symbol",
                                           host = "https://dec2021.archive.ensembl.org/",
                                           verbose = TRUE) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required for mouse-to-human conversion.\n",
         "Run: BiocManager::install('biomaRt')")
  }

  expr <- as.data.frame(expr, check.names = FALSE)
  validate_tme_input(expr)

  mouse_symbols <- rownames(expr)
  mouse_symbols <- mouse_symbols[mouse_symbols != "" & !is.na(mouse_symbols)]

  # Connect to Ensembl marts. Archive host avoids dataset version drift.
  human <- tryCatch(
    biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host),
    error = function(e) stop("Failed to connect to human Ensembl mart: ", conditionMessage(e))
  )
  mouse <- tryCatch(
    biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host),
    error = function(e) stop("Failed to connect to mouse Ensembl mart: ", conditionMessage(e))
  )

  lds <- tryCatch(
    biomaRt::getLDS(
      attributes  = mouse_attr,
      filters     = mouse_attr,
      mart        = mouse,
      values      = mouse_symbols,
      attributesL = human_attr,
      martL       = human,
      uniqueRows  = TRUE
    ),
    error = function(e) stop("biomaRt::getLDS failed: ", conditionMessage(e))
  )

  if (nrow(lds) == 0) {
    stop("No mouse-to-human ortholog mappings were found. Check that rownames are MGI gene symbols.")
  }

  mgi_col <- grep("^MGI\\.symbol$|^mgi_symbol$", colnames(lds), value = TRUE)[1]
  hgnc_col <- grep("^HGNC\\.symbol$|^hgnc_symbol$", colnames(lds), value = TRUE)[1]
  if (is.na(mgi_col) || is.na(hgnc_col)) {
    stop("Unexpected columns returned by biomaRt::getLDS: ", paste(colnames(lds), collapse = ", "))
  }

  # Build a clean mapping table
  mapping <- data.frame(
    mgi_symbol = toupper(as.character(lds[[mgi_col]])),
    hgnc_symbol = as.character(lds[[hgnc_col]]),
    stringsAsFactors = FALSE
  )
  mapping <- mapping[mapping$hgnc_symbol != "" & !is.na(mapping$hgnc_symbol), ]

  if (nrow(mapping) == 0) {
    stop("All mouse-to-human mappings returned empty HGNC symbols.")
  }

  # Work with upper-case MGI symbols for matching
  expr_upper <- expr
  rownames(expr_upper) <- toupper(rownames(expr_upper))

  expr_mapped <- merge(
    mapping,
    data.frame(mgi_symbol = rownames(expr_upper), stringsAsFactors = FALSE),
    by = "mgi_symbol",
    all.y = FALSE
  )

  # Join expression values
  expr_with_human <- cbind(
    expr_upper[expr_mapped$mgi_symbol, , drop = FALSE],
    hgnc_symbol = expr_mapped$hgnc_symbol
  )

  # Deduplicate on human symbol, keeping the mouse row with highest average expression
  expr_human <- deduplicate_expression_by_symbol(expr_with_human, symbol_col = "hgnc_symbol")

  if (verbose) {
    message(
      "Mouse-to-human conversion:",
      " input genes = ", length(mouse_symbols),
      "; mapped to human = ", nrow(mapping),
      "; unique HGNC genes after dedup = ", nrow(expr_human),
      "; unmapped = ", length(setdiff(mouse_symbols, mapping$mgi_symbol))
    )
  }

  expr_human
}

# One-stop preparation for TME deconvolution.
# - Reverses log transformation if requested.
# - If species == "mouse", converts MGI symbols to HGNC symbols via biomaRt.
prepare_tme_expression <- function(expr,
                                   is_log = TRUE,
                                   species = c("human", "mouse"),
                                   log_base = 2,
                                   verbose = TRUE) {
  species <- match.arg(species)
  expr <- as.data.frame(expr, check.names = FALSE)
  validate_tme_input(expr)

  expr <- undo_log_expr(expr, is_log = is_log, log_base = log_base)
  expr <- as.data.frame(expr, check.names = FALSE)
  # undo_log_expr may return a matrix; ensure rownames survive the conversion
  if (is.null(rownames(expr))) rownames(expr) <- rownames(as.data.frame(expr, check.names = FALSE))

  if (species == "mouse") {
    if (verbose) message("Converting mouse expression matrix to human orthologs for TME deconvolution...")
    expr <- convert_mouse_symbols_to_human(expr, verbose = verbose)
  } else {
    if (verbose) message("Using human expression matrix for TME deconvolution.")
    expr <- deduplicate_expression_by_symbol(expr)
  }

  validate_tme_input(expr)
  expr
}

# Calculate a reasonable size for a TME stacked barplot.
# Width scales with the number of samples; height is modest because cell types
# are stacked.
calc_tme_barplot_size <- function(n_samples,
                                   n_celltypes = NULL,
                                   base_width = 0.55,
                                   min_width = 6,
                                   max_width = 24,
                                   height = 7) {
  width <- max(min_width, min(max_width, base_width * n_samples + 2))
  c(width = width, height = height)
}

# Calculate a reasonable size for a faceted TME boxplot.
# Width/height scale with the number of facets (cell types) arranged in `ncol`
# columns.
calc_tme_boxplot_size <- function(n_celltypes,
                                   ncol = 4,
                                   base_width = 2.4,
                                   base_height = 2.2,
                                   min_width = 8,
                                   max_width = 22,
                                   min_height = 6,
                                   max_height = 20) {
  nrows <- ceiling(n_celltypes / ncol)
  width <- max(min_width, min(max_width, base_width * min(ncol, n_celltypes) + 1.5))
  height <- max(min_height, min(max_height, base_height * nrows + 1.5))
  c(width = width, height = height)
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
                                  filename, title = "TME Cell Fractions",
                                  width = 12, height = 7,
                                  group_colors = NULL) {
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
    if (!is.null(group_colors)) {
      p <- p + ggplot2::theme(
        strip.background = ggplot2::element_rect(fill = "white", color = "black"),
        strip.text = ggplot2::element_text(color = "black", face = "bold", size = 10)
      )
    }
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}

# Plot box/violin plots for each cell type with pairwise statistics.
plot_tme_boxplot_pdf <- function(long_df, group_col = "condition", value_col = "fraction",
                                  filename, title = "TME Cell Fraction by Group",
                                  comparisons = NULL, method = "t.test", show_ns = FALSE,
                                  width = 14, height = 10,
                                  group_colors = NULL) {
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
  if (!is.null(group_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = group_colors)
  }

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

# Plot one cell type per PDF using the publication violin/box style.
plot_tme_per_celltype_pdf <- function(long_df, group_col = "condition", value_col = "fraction",
                                       celltype_col = "cell_type", sample_col = "sample",
                                       filename_prefix, title_prefix = "TME",
                                       group_colors = NULL, method = "t.test",
                                       width = 5.5, height = 6,
                                       show_ns = FALSE) {
  if (nrow(long_df) == 0 || !all(c(group_col, value_col, celltype_col) %in% colnames(long_df))) {
    return(invisible(NULL))
  }
  plot_data <- long_df[!is.na(long_df[[value_col]]) & !is.na(long_df[[group_col]]), , drop = FALSE]
  plot_data[[group_col]] <- factor(plot_data[[group_col]])
  celltypes <- unique(as.character(plot_data[[celltype_col]]))

  if (is.null(group_colors)) {
    group_colors <- make_group_colors(levels(plot_data[[group_col]]))
  }
  comparisons <- if (length(levels(plot_data[[group_col]])) >= 2) {
    utils::combn(levels(plot_data[[group_col]]), 2, simplify = FALSE)
  } else {
    list()
  }

  for (ct in celltypes) {
    sub <- plot_data[plot_data[[celltype_col]] == ct, , drop = FALSE]
    if (nrow(sub) == 0) next

    stat_tbl <- tryCatch(
      pairwise_effect_table(
        sub, value_col = value_col, group_col = group_col,
        comparisons = comparisons, method = method, p_adjust_method = "BH"
      ),
      error = function(e) data.frame()
    )

    safe_title <- gsub("[[:space:]]+", "_", ct)
    safe_title <- gsub("[^[:alnum:]_-]", "", safe_title)
    out_file <- paste0(filename_prefix, "_", safe_title, ".pdf")

    p <- ggplot2::ggplot(sub, ggplot2::aes(x = .data[[group_col]], y = .data[[value_col]], fill = .data[[group_col]])) +
      ggplot2::geom_violin(trim = FALSE, width = 0.88, alpha = 0.34, color = NA) +
      ggplot2::geom_boxplot(width = 0.22, outlier.shape = NA, fill = "white", alpha = 0.82, linewidth = 0.42, color = "#2F2F2F") +
      ggplot2::stat_summary(fun = stats::median, geom = "point", shape = 95, size = 7, color = "#2F2F2F") +
      ggplot2::geom_jitter(width = 0.08, size = 2.2, shape = 21, color = "#222222", stroke = 0.25, alpha = 0.9) +
      ggplot2::labs(title = paste(title_prefix, ct), x = NULL, y = value_col) +
      ggplot2::scale_fill_manual(values = group_colors) +
      theme_publication(base_size = 12) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
        panel.grid.major.y = ggplot2::element_line(color = "#E7E7E7", linewidth = 0.25),
        axis.text.x = ggplot2::element_text(angle = 25, hjust = 1, vjust = 1)
      )

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
    save_pdf_plot(p, out_file, width = width, height = height)
  }
  invisible(celltypes)
}

# Plot a TME score heatmap (xCell or ESTIMATE) with consistent group annotation.
plot_tme_heatmap_pdf <- function(tme_df, meta,
                                  group_col = "condition", sample_col = "sample",
                                  group_colors = NULL,
                                  filename, title = "TME Scores",
                                  width = 10, height = 12,
                                  scale = "row") {
  if (nrow(tme_df) == 0 || !sample_col %in% colnames(tme_df)) return(invisible(NULL))
  mat <- as.data.frame(tme_df)
  rownames(mat) <- as.character(mat[[sample_col]])
  mat <- mat[, setdiff(colnames(mat), sample_col), drop = FALSE]
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"

  samples <- intersect(rownames(mat), meta[[sample_col]])
  if (length(samples) == 0) return(invisible(NULL))
  mat <- mat[samples, , drop = FALSE]

  ann <- data.frame(
    Group = as.character(meta[match(samples, meta[[sample_col]]), group_col]),
    stringsAsFactors = FALSE
  )
  rownames(ann) <- samples
  ann$Group <- factor(ann$Group, levels = unique(as.character(ann$Group)))

  annotation_colors <- NULL
  if (!is.null(group_colors)) {
    annotation_colors <- list(Group = group_colors[levels(ann$Group)])
  }

  pheatmap::pheatmap(t(mat),
                     annotation_col = ann,
                     annotation_colors = annotation_colors,
                     scale = scale,
                     cluster_cols = TRUE,
                     cluster_rows = TRUE,
                     main = title,
                     filename = filename,
                     width = width,
                     height = height)
  invisible(filename)
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
                                      show_ns = FALSE, width = 10, height = 6,
                                      group_colors = NULL) {
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
  if (!is.null(group_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = group_colors)
  }
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
