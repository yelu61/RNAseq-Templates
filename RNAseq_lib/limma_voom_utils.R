# limma-voom differential expression helpers for bulk RNA-seq templates.

# Build a DGEList, filter low counts, and normalize. The explicit raw-count
# fraction is additional to edgeR's group-aware CPM/total-count filter; set it
# to zero to use edgeR filtering alone (the behavior before this option worked).
prepare_dge_for_voom <- function(counts, group = NULL,
                                  min_counts_per_sample = 10,
                                  min_sample_frac = 0.5,
                                  normalize_method = "TMM") {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for limma-voom workflow.")
  }
  if (!is.numeric(min_sample_frac) || length(min_sample_frac) != 1L ||
      !is.finite(min_sample_frac) || min_sample_frac < 0 || min_sample_frac > 1) {
    stop("min_sample_frac must be a finite number in [0, 1].")
  }
  if (!is.numeric(min_counts_per_sample) || length(min_counts_per_sample) != 1L ||
      !is.finite(min_counts_per_sample) || min_counts_per_sample < 0) {
    stop("min_counts_per_sample must be a finite nonnegative number.")
  }
  counts <- as.matrix(counts)
  mode(counts) <- "numeric"

  dge <- edgeR::DGEList(counts = counts, group = group)
  keep <- edgeR::filterByExpr(dge, min.count = min_counts_per_sample)
  min_samples <- ceiling(min_sample_frac * ncol(counts))
  keep <- keep & rowSums(counts >= min_counts_per_sample) >= min_samples
  if (!any(keep)) stop("No genes remain after edgeR and sample-fraction filtering.")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge, method = normalize_method)
  dge
}

# Run voom with an optional design matrix.
run_voom <- function(dge, design = NULL, plot_file = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.")
  }
  if (is.null(plot_file)) {
    return(limma::voom(dge, design = design, plot = FALSE))
  }

  # Run voom exactly once, inside the managed graphics device. The previous
  # implementation plotted once to R's default device and recomputed voom for
  # the canonical PDF, leaving an ignored Rplots.pdf beside the runner.
  v <- NULL
  save_pdf_device(plot_file, width = 5.2, height = 4.6, draw = function() {
    v <<- limma::voom(dge, design = design, plot = TRUE, save.plot = TRUE)
  })
  v
}

# Optional batch correction using limma::removeBatchEffect.
remove_batch_effect_voom <- function(v, batch = NULL, batch2 = NULL, covariates = NULL, design = NULL) {
  if (is.null(batch)) return(v)
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.")
  }
  corrected_E <- limma::removeBatchEffect(
    v$E,
    batch = batch,
    batch2 = batch2,
    covariates = covariates,
    design = design
  )
  v$E <- corrected_E
  v
}

# Fit limma model and extract all requested contrasts.
# comparisons: list of 3-element vectors c(name, numerator, denominator)
run_limma_contrasts <- function(v, design, comparisons) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required.")
  }
  fit <- limma::lmFit(v, design)
  res_list <- list()

  for (comp in comparisons) {
    comp_name <- comp[1]
    treat <- comp[2]
    ctrl <- comp[3]
    missing_terms <- setdiff(c(treat, ctrl), colnames(design))
    if (length(missing_terms) > 0) {
      # A missing comparison is a config error, not a skippable condition:
      # warning + next would silently drop it from the results.
      stop("Contrast terms not found in design matrix for comparison '", comp_name, "': ",
           paste(missing_terms, collapse = ", "),
           "\nAvailable design columns: ", paste(colnames(design), collapse = ", "),
           "\nCheck that COMPARISONS group names match GROUP_LEVELS exactly.")
    }
    # Build the contrast vector programmatically instead of parsing a
    # "treat-ctrl" string, so non-syntactic group names (e.g. "T-cell") work.
    contrast_vec <- stats::setNames(rep(0, ncol(design)), colnames(design))
    contrast_vec[treat] <- 1
    contrast_vec[ctrl] <- -1
    fit2 <- limma::contrasts.fit(fit, contrast_vec)
    fit2 <- limma::eBayes(fit2)
    top <- limma::topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "p")
    top$gene_name <- rownames(top)
    top$Comparison <- comp_name
    # Standardize column names to match DESeq2 output where possible
    top <- top |>
      dplyr::rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val) |>
      dplyr::select(gene_name, log2FoldChange, AveExpr, t, B, pvalue, padj, Comparison, dplyr::everything())
    res_list[[comp_name]] <- top
  }
  res_list
}

# Build a no-intercept design matrix from a grouping vector.
# Batch belongs in the fitted model. removeBatchEffect() is intended for
# visualization and must not replace batch adjustment in differential testing.
make_group_design <- function(group, batch = NULL) {
  group <- factor(group)
  if (is.null(batch)) {
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
  } else {
    batch <- factor(batch)
    if (length(batch) != length(group)) stop("batch length must equal group length.")
    if (nlevels(batch) < 2) stop("batch must contain at least two levels.")
    design <- stats::model.matrix(~ 0 + group + batch)
    colnames(design)[seq_len(nlevels(group))] <- levels(group)
    if (qr(design)$rank < ncol(design)) {
      stop("The group + batch design matrix is not full rank. Group and batch may be confounded.")
    }
  }
  design
}

# Summarize DEG counts from limma result list.
summarize_limma_deg <- function(res_list, padj_cutoff = 0.05, lfc_cutoff = 0.5,
                                   pvalue_column = "padj", lfc_column = "log2FoldChange") {
  rows <- lapply(names(res_list), function(comp_name) {
    res <- res_list[[comp_name]]
    up <- sum(!is.na(res[[pvalue_column]]) & res[[pvalue_column]] < padj_cutoff & res[[lfc_column]] > lfc_cutoff, na.rm = TRUE)
    down <- sum(!is.na(res[[pvalue_column]]) & res[[pvalue_column]] < padj_cutoff & res[[lfc_column]] < -lfc_cutoff, na.rm = TRUE)
    data.frame(Comparison = comp_name, UP = up, DOWN = down, Total = up + down, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Write limma result tables and a summary.
write_limma_results <- function(res_list, outdir = "1-DEG") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  for (comp_name in names(res_list)) {
    utils::write.csv(res_list[[comp_name]],
                     file.path(outdir, paste0("limma_voom_DEG_", comp_name, ".csv")),
                     row.names = FALSE)
  }
  summary_df <- summarize_limma_deg(res_list)
  utils::write.csv(summary_df, file.path(outdir, "limma_voom_DEG_summary.csv"), row.names = FALSE)
  invisible(summary_df)
}

# Make a ranked gene list for GSEA from a limma result table.
ranked_gene_list_limma <- function(res_df, rank_column = "t") {
  if (!rank_column %in% colnames(res_df)) {
    stop("Rank column not found: ", rank_column)
  }
  gene_list <- res_df[[rank_column]]
  names(gene_list) <- res_df$gene_name
  gene_list <- gene_list[!is.na(gene_list)]
  sort(gene_list, decreasing = TRUE)
}
