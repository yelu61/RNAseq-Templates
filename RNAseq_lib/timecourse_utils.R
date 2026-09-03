# Time-course / temporal pattern helpers for bulk RNA-seq templates.

# Aggregate expression to group means for time-course clustering.
# expr: matrix genes x samples
# group_vec: vector assigning each sample to a time/group
aggregate_expr_by_group <- function(expr, group_vec, fun = mean) {
  group_vec <- factor(group_vec)
  t(apply(expr, 1, function(x) {
    vapply(split(x, group_vec), fun, numeric(1))
  }))
}

# Prepare an ExpressionSet for Mfuzz from a gene x sample matrix.
prepare_mfuzz_eset <- function(expr, na_thres = 0.25, fill_mode = "mean", min_std = 0) {
  if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("Package 'Mfuzz' is required for time-course clustering.")
  }
  # ExpressionSet lives in Biobase (Mfuzz depends on it but does not re-export it
  # in current Bioconductor), so construct it via Biobase rather than Mfuzz::.
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package 'Biobase' is required to build the ExpressionSet for Mfuzz.")
  }
  dat <- Biobase::ExpressionSet(assayData = as.matrix(expr))
  dat <- Mfuzz::filter.NA(dat, thres = na_thres)
  dat <- Mfuzz::fill.NA(dat, mode = fill_mode)
  # Mfuzz::filter.std() plots its ordered-SD diagnostic by default, which opens
  # Rplots.pdf in a headless runner. The pipeline exports explicit diagnostics,
  # so keep this preprocessing step non-visual and device-safe.
  dat <- Mfuzz::filter.std(dat, min.std = min_std, visu = FALSE)
  dat <- Mfuzz::standardise(dat)
  dat
}

# Run Mfuzz soft clustering.
run_mfuzz <- function(eset, n_clusters = 5, seed = 2025, m = NULL) {
  if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("Package 'Mfuzz' is required.")
  }
  set.seed(seed)
  if (is.null(m)) {
    m <- Mfuzz::mestimate(eset)
  }
  Mfuzz::mfuzz(eset, c = n_clusters, m = m)
}

# Extract cluster assignment and membership for all genes.
extract_mfuzz_clusters <- function(mfuzz_result, eset = NULL, min_acore = NULL) {
  clusters <- mfuzz_result$cluster
  membership <- mfuzz_result$membership
  out <- data.frame(
    gene_name = names(clusters),
    cluster = clusters,
    stringsAsFactors = FALSE
  )
  out$max_membership <- apply(membership, 1, max)
  if (!is.null(min_acore)) {
    core <- Mfuzz::acore(eset, mfuzz_result, min.acore = min_acore)
    core_genes <- unique(unlist(lapply(core, function(x) x$NAME)))
    out$core_gene <- out$gene_name %in% core_genes
  }
  out
}

# Plot Mfuzz cluster trends as PDF.
plot_mfuzz_trends_pdf <- function(eset, mfuzz_result, filename,
                                    time_labels = NULL, min_mem = 0,
                                    mfrow = NULL, width = 12, height = 7) {
  if (!requireNamespace("Mfuzz", quietly = TRUE) ||
      !requireNamespace("Biobase", quietly = TRUE)) return(invisible(NULL))
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  n_clusters <- length(mfuzz_result$size)
  if (is.null(mfrow)) {
    cols <- ceiling(sqrt(n_clusters))
    rows <- ceiling(n_clusters / cols)
    mfrow <- c(rows, cols)
  }
  if (is.null(time_labels)) {
    time_labels <- colnames(eset)
  }
  expr <- Biobase::exprs(eset)
  membership <- as.matrix(mfuzz_result$membership)
  assignments <- mfuzz_result$cluster
  cluster_colors <- rep(c(
    "#3B6FB6", "#D0604C", "#4C956C", "#8C6BB1", "#C58B35", "#4F8C8D"
  ), length.out = n_clusters)
  save_pdf_device(filename, width = width, height = height, draw = function() {
    old_par <- graphics::par(
      mfrow = mfrow, family = "sans", mar = c(3.8, 4.2, 2.4, 0.8),
      mgp = c(2.35, 0.7, 0), tcl = -0.25,
      cex.axis = 0.78, cex.lab = 0.88, cex.main = 0.92,
      las = 1
    )
    on.exit(graphics::par(old_par), add = TRUE)
    x <- seq_len(ncol(expr))
    for (cluster_id in seq_len(n_clusters)) {
      gene_ids <- names(assignments)[assignments == cluster_id]
      if (min_mem > 0) {
        gene_ids <- gene_ids[membership[gene_ids, cluster_id] >= min_mem]
      }
      gene_ids <- intersect(gene_ids, rownames(expr))
      if (!length(gene_ids)) {
        graphics::plot.new()
        graphics::title(main = paste("Cluster", cluster_id, "(n = 0)"))
        next
      }
      cluster_expr <- expr[gene_ids, , drop = FALSE]
      y_range <- range(cluster_expr, finite = TRUE)
      if (!all(is.finite(y_range)) || diff(y_range) == 0) y_range <- c(-1, 1)
      pad <- max(diff(y_range) * 0.06, 0.08)
      graphics::plot(
        x, rep(NA_real_, length(x)), type = "n", xaxt = "n",
        xlab = "Time", ylab = "Standardized expression",
        xlim = range(x), ylim = y_range + c(-pad, pad),
        main = sprintf("Cluster %d  (n = %d)", cluster_id, length(gene_ids)),
        bty = "l"
      )
      graphics::axis(1, at = x, labels = time_labels, las = 1)
      graphics::abline(h = 0, col = "#D9D9D9", lwd = 0.6)
      line_color <- grDevices::adjustcolor(cluster_colors[[cluster_id]], alpha.f = 0.18)
      graphics::matlines(x, t(cluster_expr), col = line_color, lty = 1, lwd = 0.45)
      weights <- membership[gene_ids, cluster_id]
      centroid <- colSums(cluster_expr * weights, na.rm = TRUE) / colSums(!is.na(cluster_expr) * weights)
      graphics::lines(x, centroid, col = cluster_colors[[cluster_id]], lwd = 2.2)
      graphics::points(x, centroid, col = cluster_colors[[cluster_id]], pch = 16, cex = 0.55)
    }
  })
  invisible(filename)
}

# Plot a time-course heatmap with cluster-based row split.
plot_timecourse_heatmap_pdf <- function(expr, cluster_df, group_vec, group_levels,
                                         filename, group_colors = NULL,
                                         z_cap = 2, width = mm_to_in(183), height = mm_to_in(247),
                                         row_font_size = 6) {
  if (nrow(expr) == 0 || ncol(expr) == 0) return(invisible(NULL))
  cluster_df <- cluster_df[cluster_df$gene_name %in% rownames(expr), ]
  mat <- expr[cluster_df$gene_name, ]
  z <- t(scale(t(mat)))
  z[z > z_cap] <- z_cap
  z[z < -z_cap] <- -z_cap
  z[is.na(z)] <- 0

  cluster_df$cluster <- factor(paste0("Cluster", cluster_df$cluster))
  row_order <- order(cluster_df$cluster, cluster_df$max_membership, decreasing = c(FALSE, TRUE))
  z <- z[row_order, ]
  row_split <- cluster_df$cluster[row_order]

  top_annotation <- NULL
  if (!is.null(group_vec)) {
    group_vec <- factor(group_vec, levels = group_levels %||% unique(group_vec))
    top_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Group = group_vec,
      col = if (!is.null(group_colors)) list(Group = group_colors) else NULL,
      annotation_name_side = "left"
    )
  }

  col_fun <- circlize::colorRamp2(c(-z_cap, 0, z_cap), c("#2166ac", "white", "#b2182b"))
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ht <- ComplexHeatmap::Heatmap(
    z,
    name = "Z-score",
    col = col_fun,
    top_annotation = top_annotation,
    row_split = row_split,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    cluster_column_slices = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    use_raster = nrow(z) > 500 && .has_working_cairo(),
    row_names_gp = grid::gpar(fontsize = row_font_size)
  )
  save_pdf_device(filename, width = width, height = height, draw = function() {
    ComplexHeatmap::draw(ht)
  })
  invisible(ht)
}

# Run ORA on each Mfuzz cluster.
run_mfuzz_cluster_ora <- function(cluster_df, org_db, universe = NULL,
                                  min_genes = 5, p_cutoff = 0.05, q_cutoff = 0.2,
                                  ontology = "BP") {
  cluster_ids <- sort(unique(cluster_df$cluster))
  results <- list()
  for (cl in cluster_ids) {
    genes <- cluster_df$gene_name[cluster_df$cluster == cl]
    if (length(genes) < min_genes) next
    mapped <- map_symbols_to_entrez(genes, org_db)
    if (nrow(mapped) < min_genes) next
    ego <- clusterProfiler::enrichGO(
      gene = mapped$ENTREZID,
      universe = universe,
      OrgDb = org_db,
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = p_cutoff,
      qvalueCutoff = q_cutoff,
      readable = TRUE
    )
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      results[[paste0("Cluster_", cl)]] <- ego
    }
  }
  results
}

# Time-point vs baseline DESeq2.
# Fits a single DESeq2 model and extracts contrasts for every non-baseline time point vs baseline.
run_timepoint_vs_baseline_deseq2 <- function(count_data, col_data, time_col = "time",
                                              baseline_time = NULL, condition_col = NULL,
                                              subject_col = NULL, design_formula = NULL,
                                              alpha = 0.05, shrink_type = "ashr") {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for time-point DEG analysis.")
  }
  if (!time_col %in% colnames(col_data)) {
    stop("time_col not found in col_data: ", time_col)
  }

  time_levels <- levels(factor(col_data[[time_col]]))
  if (is.null(baseline_time)) {
    baseline_time <- time_levels[1]
    message("Using earliest time level as baseline: ", baseline_time)
  }
  if (!baseline_time %in% time_levels) {
    stop("baseline_time ", baseline_time, " not found in time levels: ", paste(time_levels, collapse = ", "))
  }
  non_baseline <- setdiff(time_levels, baseline_time)
  if (length(non_baseline) == 0) {
    stop("Only one time level found; cannot compare time points vs baseline.")
  }

  col_data[[time_col]] <- factor(col_data[[time_col]], levels = time_levels)
  if (!is.null(condition_col) && condition_col %in% colnames(col_data)) {
    col_data[[condition_col]] <- factor(col_data[[condition_col]])
  }
  if (!is.null(subject_col) && subject_col %in% colnames(col_data)) {
    col_data[[subject_col]] <- factor(col_data[[subject_col]])
  }

  if (is.null(design_formula)) {
    if (!is.null(subject_col) && subject_col %in% colnames(col_data)) {
      design_formula <- stats::as.formula(paste0("~ ", subject_col, " + ", time_col))
      message("Using paired design formula: ", deparse(design_formula))
    } else if (!is.null(condition_col) && condition_col %in% colnames(col_data)) {
      design_formula <- stats::as.formula(paste0("~ ", condition_col, " + ", time_col))
      message("Using design formula: ", deparse(design_formula))
    } else {
      design_formula <- stats::as.formula(paste0("~ ", time_col))
      message("Using design formula: ", deparse(design_formula))
    }
  }

  # Result extraction only supports additive designs: every non-baseline time
  # point vs baseline via Wald contrasts. Interaction designs (condition x
  # time, LRT-style) are not implemented; refuse them explicitly instead of
  # fitting the interaction model and silently extracting main-effect contrasts.
  interaction_terms <- grep(":", attr(stats::terms(design_formula), "term.labels"), value = TRUE)
  if (length(interaction_terms) > 0) {
    stop("Interaction terms are not supported by timepoint-vs-baseline extraction: ",
         paste(interaction_terms, collapse = ", "),
         "\nUse an additive design (e.g. ~ condition + time). Condition x time questions require an LRT workflow that this template does not implement.")
  }

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_data,
    colData = col_data,
    design = design_formula
  )
  dds <- DESeq2::DESeq(dds)

  comparisons <- lapply(non_baseline, function(tp) {
    c(paste0("time_", tp, "_vs_", baseline_time), tp, baseline_time)
  })
  res_list <- extract_deseq2_results(
    dds = dds,
    comparisons = comparisons,
    condition_col = time_col,
    alpha = alpha,
    shrink_type = shrink_type
  )
  res_list
}

# Summarize DEG counts per time-point vs baseline.
summarize_timepoint_deg <- function(res_list, padj_cutoff = 0.05, lfc_cutoff = 0.5,
                                     pvalue_column = "padj", lfc_column = "log2FoldChange") {
  rows <- lapply(names(res_list), function(comp_name) {
    res <- res_list[[comp_name]]
    up <- sum(!is.na(res[[pvalue_column]]) & res[[pvalue_column]] < padj_cutoff & res[[lfc_column]] > lfc_cutoff, na.rm = TRUE)
    down <- sum(!is.na(res[[pvalue_column]]) & res[[pvalue_column]] < padj_cutoff & res[[lfc_column]] < -lfc_cutoff, na.rm = TRUE)
    data.frame(Comparison = comp_name, UP = up, DOWN = down, Total = up + down, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Write all-gene and threshold-specific time-point DEG results.
write_timepoint_deg_results <- function(res_list, outdir,
                                         threshold_grid = NULL,
                                         pvalue_column = "padj", lfc_column = "log2FoldChange") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  write_all_gene_deg_results(res_list, outdir = outdir)
  summary_df <- summarize_timepoint_deg(res_list, pvalue_column = pvalue_column, lfc_column = lfc_column)
  utils::write.csv(summary_df, file.path(outdir, "Timepoint_DEG_summary.csv"), row.names = FALSE)

  if (!is.null(threshold_grid)) {
    deg_by_threshold <- build_deg_threshold_sets(
      res_list,
      threshold_grid,
      pvalue_column = pvalue_column,
      lfc_column = lfc_column
    )
    write_deg_threshold_outputs(
      deg_by_threshold,
      outdir = outdir,
      threshold_grid = threshold_grid,
      pvalue_column = pvalue_column,
      lfc_column = lfc_column
    )
  }
  invisible(summary_df)
}

# Plot time-point DEG summary barplot.
plot_timepoint_deg_summary_pdf <- function(summary_df, filename,
                                            width = 8, height = 6) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(invisible(NULL))
  df <- tidyr::pivot_longer(
    summary_df,
    cols = c("UP", "DOWN"),
    names_to = "Direction",
    values_to = "Count"
  )
  df$Direction <- factor(df$Direction, levels = c("UP", "DOWN"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Comparison, y = Count, fill = Direction)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::scale_fill_manual(values = c("UP" = "#d6604d", "DOWN" = "#4393c3")) +
    ggplot2::labs(x = NULL, y = "Number of DEGs", title = "Time-point vs Baseline DEGs") +
    theme_publication(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  save_pdf_plot(p, filename, width = width, height = height)
  invisible(p)
}

# Write Mfuzz cluster assignments and memberships to CSV.
write_mfuzz_cluster_table <- function(cluster_df, filename) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(cluster_df, filename, row.names = FALSE)
  invisible(filename)
}

# Summarize cluster sizes.
summarize_mfuzz_clusters <- function(cluster_df) {
  cluster_df |>
    dplyr::count(.data$cluster, name = "n_genes") |>
    dplyr::arrange(.data$cluster)
}
