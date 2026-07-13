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
                                   is_log = FALSE,
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
                                  group_colors = NULL, p_adjust_method = "BH") {
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
      p.adjust.method = p_adjust_method,
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
# Samples are ordered by group, then clustered within each group, so replicates
# of the same condition appear together while still preserving within-group
# hierarchical structure.
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

  # Build a group-ordered, within-group clustered column order.
  # Do this by hierarchical clustering within each group slice and concatenating.
  order_grouped_columns <- function(m, group_vec, group_levels = NULL) {
    if (is.null(group_levels)) group_levels <- unique(group_vec)
    ordered_cols <- character(0)
    for (g in group_levels) {
      idx <- which(group_vec == g)
      if (length(idx) == 0) next
      if (length(idx) == 1) {
        ordered_cols <- c(ordered_cols, colnames(m)[idx])
      } else {
        sub_m <- m[, idx, drop = FALSE]
        # Remove constant rows before clustering to avoid dist failures
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

  group_vec <- ann$Group
  names(group_vec) <- rownames(ann)
  col_order <- order_grouped_columns(t(mat), group_vec, levels(ann$Group))
  mat <- mat[col_order, , drop = FALSE]
  ann <- ann[col_order, , drop = FALSE]

  pheatmap::pheatmap(t(mat),
                     annotation_col = ann,
                     annotation_colors = annotation_colors,
                     scale = scale,
                     cluster_cols = FALSE,
                     cluster_rows = TRUE,
                     main = title,
                     filename = filename,
                     width = width,
                     height = height)
  invisible(filename)
}

# Extract ESTIMATE-like scores from an IOBR estimate result and make them long-format.
# IOBR appends the method name as a suffix (e.g. "StromalScore_estimate"); native
# ESTIMATE does not. Accept both naming styles and strip any method suffix.
melt_estimate_scores <- function(estimate_df, id_column = "ID", group_df = NULL,
                                 sample_col = "sample", group_col = "condition") {
  score_cols <- intersect(c(
    "StromalScore", "ImmuneScore", "ESTIMATEScore", "TumorPurity",
    "StromalScore_estimate", "ImmuneScore_estimate",
    "ESTIMATEScore_estimate", "TumorPurity_estimate"
  ), colnames(estimate_df))
  if (length(score_cols) == 0) return(NULL)
  long <- estimate_df |>
    dplyr::rename(!!sample_col := dplyr::all_of(id_column)) |>
    tidyr::pivot_longer(cols = dplyr::all_of(score_cols), names_to = "score_type", values_to = "score") |>
    dplyr::mutate(
      score = as.numeric(.data$score),
      score_type = gsub("_estimate$", "", .data$score_type)
    )
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
                                      group_colors = NULL, ncol = 4,
                                      save_individual = FALSE, individual_prefix = NULL,
                                      p_adjust_method = "BH") {
  if (is.null(long_df) || nrow(long_df) == 0) return(invisible(NULL))
  plot_data <- long_df
  plot_data[[group_col]] <- factor(plot_data[[group_col]])
  group_levels <- levels(plot_data[[group_col]])
  comparisons <- if (length(group_levels) >= 2) utils::combn(group_levels, 2, simplify = FALSE) else list()

  if (is.null(ncol)) ncol <- min(4, length(unique(plot_data$score_type)))

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[group_col]], y = .data$score, fill = .data[[group_col]])) +
    ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, linewidth = 0.4) +
    ggplot2::geom_jitter(width = 0.12, size = 2, shape = 21, color = "#222222", stroke = 0.25, fill = "white", alpha = 0.9) +
    ggplot2::facet_wrap(~ score_type, scales = "free_y", ncol = ncol) +
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
      p.adjust.method = p_adjust_method,
      hide.ns = !show_ns,
      tip.length = 0.015,
      bracket.size = 0.4,
      size = 3
    )
  }
  save_pdf_plot(p, filename, width = width, height = height)

  if (isTRUE(save_individual) && !is.null(individual_prefix)) {
    score_types <- unique(as.character(plot_data$score_type))
    for (st in score_types) {
      sub <- plot_data[plot_data$score_type == st, , drop = FALSE]
      sub_comparisons <- if (length(group_levels) >= 2) utils::combn(group_levels, 2, simplify = FALSE) else list()
      p_st <- ggplot2::ggplot(sub, ggplot2::aes(x = .data[[group_col]], y = .data$score, fill = .data[[group_col]])) +
        ggplot2::geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85, linewidth = 0.4) +
        ggplot2::geom_jitter(width = 0.12, size = 2.5, shape = 21, color = "#222222", stroke = 0.25, fill = "white", alpha = 0.9) +
        ggplot2::labs(x = NULL, y = "Score", title = paste(title, "-", st)) +
        theme_publication(base_size = 12) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
          legend.position = "none"
        )
      if (!is.null(group_colors)) {
        p_st <- p_st + ggplot2::scale_fill_manual(values = group_colors)
      }
      if (length(sub_comparisons) > 0) {
        p_st <- p_st + ggpubr::stat_compare_means(
          comparisons = sub_comparisons,
          method = method,
          p.adjust.method = p_adjust_method,
          hide.ns = !show_ns,
          tip.length = 0.015,
          bracket.size = 0.4,
          size = 3
        )
      }
      out_file <- paste0(individual_prefix, "_", gsub("[^[:alnum:]_-]", "_", st), ".pdf")
      save_pdf_plot(p_st, out_file, width = 5.5, height = 6)
    }
  }

  p
}

# -----------------------------------------------------------------------------
# Native CIBERSORT wrapper + native-vs-IOBR comparison
# -----------------------------------------------------------------------------

# Default paths for the bundled CIBERSORT resources.
# These locate files from the project git root when possible, falling back to
# a path relative to the current working directory.
.default_cibersort_script <- function() {
  root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root), error = function(e) NULL)
  candidates <- c(
    if (!is.null(root)) file.path(root, "references", "CIBERSORT", "CIBERSORT.R") else NULL,
    file.path("references", "CIBERSORT", "CIBERSORT.R")
  )
  for (cand in candidates) {
    if (file.exists(cand)) return(cand)
  }
  candidates[length(candidates)]
}
.default_cibersort_signature <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  sig_name <- if (species == "human") "LM22.txt" else "cibersort_mouse_22.csv"
  root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root), error = function(e) NULL)
  candidates <- c(
    if (!is.null(root)) file.path(root, "references", "CIBERSORT", sig_name) else NULL,
    file.path("references", "CIBERSORT", sig_name)
  )
  for (cand in candidates) {
    if (file.exists(cand)) return(cand)
  }
  candidates[length(candidates)]
}

# Read a CIBERSORT signature file into a matrix.
# Supports the bundled LM22.txt (tab-delimited) and cibersort_mouse_22.csv.
read_cibersort_signature <- function(signature_file) {
  if (is.matrix(signature_file) || is.data.frame(signature_file)) {
    return(as.matrix(signature_file))
  }
  if (!is.character(signature_file) || length(signature_file) != 1) {
    stop("signature_file must be a file path or a matrix/data.frame.")
  }
  if (!file.exists(signature_file)) {
    stop("CIBERSORT signature file not found: ", signature_file)
  }

  ext <- tolower(tools::file_ext(signature_file))
  if (ext == "csv") {
    sig <- utils::read.csv(signature_file, row.names = 1, check.names = FALSE)
  } else {
    sig <- utils::read.delim(signature_file, row.names = 1, check.names = FALSE)
  }
  sig <- as.matrix(sig)
  if (!is.numeric(sig)) {
    stop("Signature matrix must be numeric after reading.")
  }
  sig
}

# Run the bundled native CIBERSORT implementation on a prepared expression
# matrix. The expression matrix should have gene symbols as rownames and
# non-log normalized values (TPM/FPKM-like). Use prepare_tme_expression() first
# if the input is log-scaled and/or mouse data that needs ortholog conversion.
run_native_cibersort <- function(expr,
                                  signature_file = .default_cibersort_signature("human"),
                                  cibersort_script = .default_cibersort_script(),
                                  is_log = FALSE,
                                  perm = 1000,
                                  QN = FALSE,
                                  id_column = "ID",
                                  verbose = TRUE) {
  if (!file.exists(cibersort_script)) {
    stop("CIBERSORT script not found: ", cibersort_script,
         "\nMake sure references/CIBERSORT/CIBERSORT.R exists.")
  }

  expr <- as.data.frame(expr, check.names = FALSE)
  validate_tme_input(expr)

  # Native CIBERSORT expects non-log data. It performs its own anti-log when
  # max(expr) < 50, but being explicit here avoids surprises.
  expr <- undo_log_expr(expr, is_log = is_log, log_base = 2)
  expr <- as.data.frame(expr, check.names = FALSE)
  validate_tme_input(expr)
  expr <- deduplicate_expression_by_symbol(expr)

  if (verbose) {
    message(
      "Running native CIBERSORT: ", ncol(expr), " samples; ",
      nrow(expr), " genes; perm = ", perm, "; QN = ", QN
    )
  }

  # The bundled CIBERSORT.R needs e1071/preprocessCore/future/furrr/purrr.
  required_pkgs <- c("e1071", "preprocessCore", "future", "furrr", "purrr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Native CIBERSORT requires packages: ", paste(missing_pkgs, collapse = ", "),
         "\nInstall via BiocManager/install.packages as needed.")
  }

  sig <- read_cibersort_signature(signature_file)
  if (verbose) {
    message("Signature matrix: ", nrow(sig), " genes x ", ncol(sig), " cell types")
  }

  # Source in a private environment to avoid polluting the global namespace,
  # then call cibersort().
  cibersort_env <- new.env(parent = .BaseNamespaceEnv)
  source(cibersort_script, local = cibersort_env)
  if (!exists("cibersort", envir = cibersort_env)) {
    stop("cibersort() function was not found in ", cibersort_script)
  }
  cibersort_fn <- get("cibersort", envir = cibersort_env)

  res_mat <- cibersort_fn(sig_matrix = sig, mixture_file = as.matrix(expr),
                          perm = perm, QN = QN)

  # Convert to tidy data.frame
  res_df <- as.data.frame(res_mat, stringsAsFactors = FALSE)
  res_df[[id_column]] <- rownames(res_df)
  res_df <- res_df[, c(id_column, setdiff(colnames(res_df), id_column)), drop = FALSE]

  # Ensure numeric columns are numeric
  numeric_cols <- setdiff(colnames(res_df), id_column)
  for (col in numeric_cols) {
    res_df[[col]] <- as.numeric(res_df[[col]])
  }

  if (verbose) {
    message("Native CIBERSORT completed for ", nrow(res_df), " samples.")
  }
  res_df
}

# Align native and IOBR CIBERSORT result tables and compute per-cell-type
# concordance metrics (Pearson/Spearman correlation, RMSE, mean absolute
# difference, paired t-test p-value).
compare_native_iobr_cibersort <- function(native_df,
                                          iobr_df,
                                          id_column = "ID",
                                          method = c("pearson", "spearman")) {
  method <- match.arg(method)

  if (!id_column %in% colnames(native_df) || !id_column %in% colnames(iobr_df)) {
    stop("id_column '", id_column, "' must be present in both native_df and iobr_df.")
  }

  native_df[[id_column]] <- as.character(native_df[[id_column]])
  iobr_df[[id_column]] <- as.character(iobr_df[[id_column]])

  # IOBR may rename the sample ID column; we keep the caller-supplied id_column.
  common_samples <- intersect(native_df[[id_column]], iobr_df[[id_column]])
  if (length(common_samples) == 0) {
    stop("No common sample IDs between native and IOBR CIBERSORT results.")
  }

  native_sub <- native_df[native_df[[id_column]] %in% common_samples, , drop = FALSE]
  iobr_sub <- iobr_df[iobr_df[[id_column]] %in% common_samples, , drop = FALSE]

  rownames(native_sub) <- native_sub[[id_column]]
  rownames(iobr_sub) <- iobr_sub[[id_column]]
  native_sub <- native_sub[common_samples, , drop = FALSE]
  iobr_sub <- iobr_sub[common_samples, , drop = FALSE]

  # Cell-type columns: everything except the ID and the diagnostic columns
  native_cell_cols <- setdiff(colnames(native_sub), c(id_column, "P-value", "Correlation", "RMSE"))
  iobr_cell_cols <- setdiff(colnames(iobr_sub), id_column)

  # Try to match cell-type names. IOBR names may contain dots, spaces, or
  # slightly different capitalization. Use a relaxed matching that tolerates
  # differences in punctuation/whitespace.
  .normalize_celltype_name <- function(x) {
    x <- tolower(as.character(x))
    x <- gsub("[._]+", " ", x)
    x <- gsub("\\s+", " ", x)
    x <- trimws(x)
    x
  }
  native_norm <- .normalize_celltype_name(native_cell_cols)
  iobr_norm <- .normalize_celltype_name(iobr_cell_cols)

  common_cells <- intersect(native_norm, iobr_norm)
  if (length(common_cells) == 0) {
    stop("No common cell types between native and IOBR CIBERSORT results.")
  }

  summary_list <- list()
  for (cell in common_cells) {
    native_col <- native_cell_cols[which(native_norm == cell)[1]]
    iobr_col <- iobr_cell_cols[which(iobr_norm == cell)[1]]
    x <- as.numeric(native_sub[[native_col]])
    y <- as.numeric(iobr_sub[[iobr_col]])

    valid <- stats::complete.cases(x, y)
    x_valid <- x[valid]
    y_valid <- y[valid]

    n <- length(x_valid)
    corr <- if (n >= 2) stats::cor(x_valid, y_valid, method = method) else NA_real_
    rmse <- if (n >= 1) sqrt(mean((x_valid - y_valid)^2, na.rm = TRUE)) else NA_real_
    mae <- if (n >= 1) mean(abs(x_valid - y_valid), na.rm = TRUE) else NA_real_
    ttest <- tryCatch(
      stats::t.test(x_valid, y_valid, paired = TRUE),
      error = function(e) NULL
    )

    summary_list[[cell]] <- data.frame(
      cell_type = cell,
      n = n,
      native_mean = mean(x_valid, na.rm = TRUE),
      iobr_mean = mean(y_valid, na.rm = TRUE),
      correlation = corr,
      rmse = rmse,
      mae = mae,
      paired_t_pvalue = if (!is.null(ttest)) ttest$p.value else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL

  # Also return the aligned long-format data for plotting
  long_list <- list()
  for (cell in common_cells) {
    native_col <- native_cell_cols[which(native_norm == cell)[1]]
    iobr_col <- iobr_cell_cols[which(iobr_norm == cell)[1]]
    long_list[[cell]] <- data.frame(
      sample = common_samples,
      cell_type = cell,
      native = as.numeric(native_sub[[native_col]]),
      iobr = as.numeric(iobr_sub[[iobr_col]]),
      difference = as.numeric(native_sub[[native_col]]) - as.numeric(iobr_sub[[iobr_col]]),
      stringsAsFactors = FALSE
    )
  }
  long_df <- do.call(rbind, long_list)

  list(summary = summary_df, long = long_df)
}

# Faceted scatter plot of native vs IOBR CIBERSORT fractions, one panel per
# matched cell type. Optionally overlays the overall Pearson correlation.
plot_cibersort_correlation_pdf <- function(long_df,
                                           filename,
                                           title = "Native vs IOBR CIBERSORT",
                                           width = 14,
                                           height = 12) {
  if (nrow(long_df) == 0 || !all(c("native", "iobr", "cell_type") %in% colnames(long_df))) {
    return(invisible(NULL))
  }

  # per-cell-type correlation for facet labels
  corr_df <- long_df |>
    dplyr::group_by(.data$cell_type) |>
    dplyr::summarise(
      r = stats::cor(.data$native, .data$iobr, use = "pairwise.complete.obs"),
      .groups = "drop"
    ) |>
    dplyr::mutate(label = sprintf("r = %.2f", .data$r))

  plot_data <- long_df |>
    dplyr::left_join(corr_df, by = "cell_type")

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = .data$native, y = .data$iobr)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#888888") +
    ggplot2::geom_point(shape = 21, size = 2, alpha = 0.8, fill = "#4A90A4", color = "#222222", stroke = 0.25) +
    ggplot2::geom_text(
      data = corr_df,
      ggplot2::aes(x = Inf, y = -Inf, label = .data$label),
      hjust = 1.1, vjust = -0.4, size = 3, color = "#333333", inherit.aes = FALSE
    ) +
    ggplot2::facet_wrap(~ cell_type, scales = "free", ncol = 5) +
    ggplot2::labs(
      x = "Native CIBERSORT fraction",
      y = "IOBR CIBERSORT fraction",
      title = title
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1,
      strip.text = ggplot2::element_text(size = 8, face = "bold")
    )

  save_pdf_plot(p, filename, width = width, height = height)
  p
}

# Bland-Altman-style difference plot: (native - IOBR) vs mean, faceted by cell
# type. Helps identify cell types where the two methods disagree systematically.
plot_cibersort_difference_pdf <- function(long_df,
                                          filename,
                                          title = "Native - IOBR CIBERSORT Difference",
                                          width = 14,
                                          height = 12) {
  if (nrow(long_df) == 0 || !all(c("native", "iobr", "cell_type") %in% colnames(long_df))) {
    return(invisible(NULL))
  }

  plot_data <- long_df |>
    dplyr::mutate(mean = (.data$native + .data$iobr) / 2,
                  diff = .data$native - .data$iobr)

  # limits of agreement per cell type
  loa_df <- plot_data |>
    dplyr::group_by(.data$cell_type) |>
    dplyr::summarise(
      mean_diff = mean(.data$diff, na.rm = TRUE),
      sd_diff = stats::sd(.data$diff, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      upper_loa = .data$mean_diff + 1.96 * .data$sd_diff,
      lower_loa = .data$mean_diff - 1.96 * .data$sd_diff
    )

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = .data$mean, y = .data$diff)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "#444444") +
    ggplot2::geom_point(shape = 21, size = 2, alpha = 0.8, fill = "#B85C38", color = "#222222", stroke = 0.25) +
    ggplot2::geom_hline(
      data = loa_df,
      ggplot2::aes(yintercept = .data$mean_diff),
      linetype = "solid", color = "#4A90A4", linewidth = 0.6
    ) +
    ggplot2::geom_hline(
      data = loa_df,
      ggplot2::aes(yintercept = .data$upper_loa),
      linetype = "dotted", color = "#4A90A4", linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      data = loa_df,
      ggplot2::aes(yintercept = .data$lower_loa),
      linetype = "dotted", color = "#4A90A4", linewidth = 0.4
    ) +
    ggplot2::facet_wrap(~ cell_type, scales = "free", ncol = 5) +
    ggplot2::labs(
      x = expression((Native + IOBR) / 2),
      y = expression(Native - IOBR),
      title = title
    ) +
    theme_publication(base_size = 10) +
    ggplot2::theme(
      aspect.ratio = 1,
      strip.text = ggplot2::element_text(size = 8, face = "bold")
    )

  save_pdf_plot(p, filename, width = width, height = height)
  p
}

# -----------------------------------------------------------------------------
# Broad-category aggregation of TME deconvolution results
# -----------------------------------------------------------------------------

# Return a mapping from CIBERSORT fine-grained cell types to broad categories.
# Supports the bundled human LM22 signature and the mouse 25-cell signature.
get_cibersort_category_map <- function(species = c("human", "mouse")) {
  species <- match.arg(species)

  if (species == "human") {
    mapping <- c(
      "B cells naive" = "B cells",
      "B cells memory" = "B cells",
      "Plasma cells" = "B cells",
      "T cells CD8" = "CD8 T cells",
      "T cells CD4 naive" = "CD4 T cells",
      "T cells CD4 memory resting" = "CD4 T cells",
      "T cells CD4 memory activated" = "CD4 T cells",
      "T cells CD4 follicular helper" = "CD4 T cells",
      "T cells follicular helper" = "CD4 T cells",
      "T cells regulatory (Tregs)" = "CD4 T cells",
      "T cells gamma delta" = "Other T cells",
      "NK cells resting" = "NK cells",
      "NK cells activated" = "NK cells",
      "Monocytes" = "Monocytes",
      "Macrophages M0" = "Macrophages",
      "Macrophages M1" = "Macrophages",
      "Macrophages M2" = "Macrophages",
      "Dendritic cells resting" = "Dendritic cells",
      "Dendritic cells activated" = "Dendritic cells",
      "Mast cells resting" = "Mast cells",
      "Mast cells activated" = "Mast cells",
      "Eosinophils" = "Eosinophils",
      "Neutrophils" = "Neutrophils"
    )
    # IOBR returns LM22 columns with underscores instead of spaces; keep both
    # the canonical spaced names and the IOBR-style underscored names.
    mapping_underscore <- stats::setNames(mapping, gsub(" ", "_", names(mapping)))
    mapping <- c(mapping, mapping_underscore)
  } else {
    mapping <- c(
      "B Cells Naive" = "B cells",
      "B Cells Memory" = "B cells",
      "Plasma Cells" = "B cells",
      "T Cells CD8 Actived" = "CD8 T cells",
      "T Cells CD8 Naive" = "CD8 T cells",
      "T Cells CD8 Memory" = "CD8 T cells",
      "T Cells CD4 Memory" = "CD4 T cells",
      "T Cells CD4 Naive" = "CD4 T cells",
      "T Cells CD4 Follicular" = "CD4 T cells",
      "Treg Cells" = "CD4 T cells",
      "Th1 Cells" = "CD4 T cells",
      "Th17 Cells" = "CD4 T cells",
      "Th2 Cells" = "CD4 T cells",
      "GammaDelta T Cells" = "Other T cells",
      "NK Resting" = "NK cells",
      "NK.Actived" = "NK cells",
      "Monocyte" = "Monocytes",
      "M0 Macrophage" = "Macrophages",
      "M1 Macrophage" = "Macrophages",
      "M2 Macrophage" = "Macrophages",
      "DC Actived" = "DCs",
      "DC Immature" = "DCs",
      "Mast Cells" = "Mast cells",
      "Neutrophil Cells" = "Neutrophils",
      "Eosinophil Cells" = "Eosinophils"
    )
  }

  data.frame(
    cell_type = names(mapping),
    category = unname(mapping),
    stringsAsFactors = FALSE
  )
}

# Return a mapping from xCell fine-grained cell types (IOBR naming) to broad
# categories. xCell column names from IOBR carry the "_xCell" suffix; matching
# is done after stripping that suffix.
get_xcell_category_map <- function() {
  mapping <- c(
    "B-cells" = "B cells",
    "Class-switched memory B-cells" = "B cells",
    "Memory B-cells" = "B cells",
    "naive B-cells" = "B cells",
    "pro B-cells" = "B cells",
    "CD4+ T-cells" = "CD4 T cells",
    "CD4+ memory T-cells" = "CD4 T cells",
    "CD4+ naive T-cells" = "CD4 T cells",
    "CD4+ Tcm" = "CD4 T cells",
    "CD4+ Tem" = "CD4 T cells",
    "Th1 cells" = "CD4 T cells",
    "Th2 cells" = "CD4 T cells",
    "Tregs" = "CD4 T cells",
    "CD8+ T-cells" = "CD8 T cells",
    "CD8+ naive T-cells" = "CD8 T cells",
    "CD8+ Tcm" = "CD8 T cells",
    "CD8+ Tem" = "CD8 T cells",
    "NKT" = "NK cells",
    "NK cells" = "NK cells",
    "Tgd cells" = "Other T cells",
    "Monocytes" = "Monocytes",
    "Macrophages" = "Macrophages",
    "Macrophages M1" = "Macrophages",
    "Macrophages M2" = "Macrophages",
    "DC" = "DCs",
    "cDC" = "DCs",
    "pDC" = "DCs",
    "iDC" = "DCs",
    "aDC" = "DCs",
    "Mast cells" = "Mast cells",
    "Eosinophils" = "Eosinophils",
    "Neutrophils" = "Neutrophils",
    "Basophils" = "Basophils",
    "Fibroblasts" = "Fibroblasts",
    "Endothelial cells" = "Endothelials",
    "ly Endothelial cells" = "Endothelials",
    "mv Endothelial cells" = "Endothelials",
    "Pericytes" = "Other Stromal",
    "Smooth muscle" = "Other Stromal",
    "Myocytes" = "Other Stromal",
    "Skeletal muscle" = "Other Stromal",
    "Mesangial cells" = "Other Stromal",
    "Epithelial cells" = "Epithelial",
    "Keratinocytes" = "Epithelial",
    "Sebocytes" = "Epithelial",
    "Hepatocytes" = "Parenchymal",
    "Neurons" = "Neuronal",
    "Astrocytes" = "Neuronal",
    "Adipocytes" = "Adipose",
    "Preadipocytes" = "Adipose",
    "Osteoblast" = "Osteoblast",
    "Megakaryocytes" = "Megakaryocytes",
    "Erythrocytes" = "Erythrocytes",
    "Platelets" = "Platelets",
    "CLP" = "Progenitors",
    "CMP" = "Progenitors",
    "GMP" = "Progenitors",
    "HSC" = "Progenitors",
    "MEP" = "Progenitors",
    "MPP" = "Progenitors",
    "MSC" = "Progenitors",
    "Chondrocytes" = "Chondrocytes"
  )

  df <- data.frame(
    cell_type = names(mapping),
    category = unname(mapping),
    stringsAsFactors = FALSE
  )
  # IOBR prefixes xCell columns with an underscore in place of spaces.
  # Keep the category labels readable, but make matching keys IOBR-compatible.
  df$cell_type <- gsub(" ", "_", df$cell_type)
  df
}

# Aggregate a wide TME result table by broad cell-type category.
#
# For CIBERSORT fractions, use method = "sum" because fractions within a
# category are additive. For xCell enrichment scores, use method = "mean"
# because scores are not additive.
#
# Returns a list with:
#   wide: aggregated data frame with id_column + one column per category
#   long: long-format data frame with id_column, category, value
aggregate_tme_by_category <- function(tme_df,
                                      id_column = "ID",
                                      category_map = NULL,
                                      method = c("sum", "mean"),
                                      suffix_strip = "_xCell$") {
  method <- match.arg(method)

  if (!id_column %in% colnames(tme_df)) {
    stop("id_column '", id_column, "' not found in tme_df.")
  }
  if (is.null(category_map)) {
    stop("category_map must be provided. Use get_cibersort_category_map() or get_xcell_category_map().")
  }

  category_map <- as.data.frame(category_map, stringsAsFactors = FALSE)
  if (!all(c("cell_type", "category") %in% colnames(category_map))) {
    stop("category_map must contain 'cell_type' and 'category' columns.")
  }

  # Identify numeric cell-type columns, excluding diagnostic/meta columns
  numeric_cols <- setdiff(colnames(tme_df), id_column)
  numeric_cols <- numeric_cols[!grepl("^(P-value|Correlation|RMSE|ImmuneScore|StromaScore|MicroenvironmentScore)(_xCell)?$", numeric_cols)]
  is_num <- vapply(tme_df[numeric_cols], is.numeric, logical(1))
  numeric_cols <- numeric_cols[is_num]

  if (length(numeric_cols) == 0) {
    stop("No numeric cell-type columns found in tme_df.")
  }

  # Strip suffixes such as "_xCell" before matching
  stripped <- sub(suffix_strip, "", numeric_cols)
  match_idx <- match(stripped, category_map$cell_type)

  if (all(is.na(match_idx))) {
    stop("No cell-type columns could be matched to category_map. ",
         "First unmatched names: ", paste(head(stripped, 5), collapse = ", "))
  }

  if (any(is.na(match_idx))) {
    unmatched <- stripped[is.na(match_idx)]
    warning("The following cell types were not matched and will be dropped: ",
            paste(unmatched, collapse = ", "))
  }

  matched_cols <- numeric_cols[!is.na(match_idx)]

  long_df <- tme_df |>
    dplyr::select(dplyr::all_of(c(id_column, matched_cols))) |>
    tidyr::pivot_longer(cols = dplyr::all_of(matched_cols), names_to = "cell_type", values_to = "value")
  long_df$category <- category_map$category[match(sub(suffix_strip, "", long_df$cell_type), category_map$cell_type)]
  long_df <- long_df[!is.na(long_df$category), , drop = FALSE]

  agg_long <- long_df |>
    dplyr::group_by(.data[[id_column]], .data$category) |>
    dplyr::summarise(
      value = if (method == "sum") {
        sum(.data$value, na.rm = TRUE)
      } else {
        mean(.data$value, na.rm = TRUE)
      },
      .groups = "drop"
    )

  agg_wide <- agg_long |>
    tidyr::pivot_wider(names_from = "category", values_from = "value")

  list(
    wide = as.data.frame(agg_wide, stringsAsFactors = FALSE),
    long = as.data.frame(agg_long, stringsAsFactors = FALSE)
  )
}
