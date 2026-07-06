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
  dat <- Mfuzz::ExpressionSet(assayData = as.matrix(expr))
  dat <- Mfuzz::filter.NA(dat, thres = na_thres)
  dat <- Mfuzz::fill.NA(dat, mode = fill_mode)
  dat <- Mfuzz::filter.std(dat, min.std = min_std)
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
  if (!requireNamespace("Mfuzz", quietly = TRUE)) return(invisible(NULL))
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
  grDevices::pdf(filename, width = width, height = height)
  Mfuzz::mfuzz.plot(
    eset,
    mfuzz_result,
    mfrow = mfrow,
    new.window = FALSE,
    time.labels = time_labels,
    min.mem = min_mem
  )
  grDevices::dev.off()
  invisible(filename)
}

# Plot a time-course heatmap with cluster-based row split.
plot_timecourse_heatmap_pdf <- function(expr, cluster_df, group_vec, group_levels,
                                         filename, group_colors = NULL,
                                         z_cap = 2, width = 9, height = 10,
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
  grDevices::pdf(filename, width = width, height = height)
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
    use_raster = nrow(z) > 500,
    row_names_gp = grid::gpar(fontsize = row_font_size)
  )
  ComplexHeatmap::draw(ht)
  grDevices::dev.off()
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

# Write cluster assignments and memberships to CSV.
write_mfuzz_cluster_table <- function(cluster_df, filename) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(cluster_df, filename, row.names = FALSE)
}

# Summarize cluster sizes.
summarize_mfuzz_clusters <- function(cluster_df) {
  cluster_df |>
    dplyr::count(.data$cluster, name = "n_genes") |>
    dplyr::arrange(.data$cluster)
}
