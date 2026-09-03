# Pathway / gene-set scoring and visualization helpers for bulk RNA-seq templates.
# GSVA-based scoring, group comparison and summary plots for custom gene sets.
# For color conventions and figure sizing see references/VISUALIZATION_STYLE_GUIDE.md.
#
# Color/theme single source of truth: group colors come from make_group_colors(),
# up/down and heatmap diverging colors from RNAseq_PALETTE, all ggplot helpers use
# theme_publication() and save via save_pdf_plot(). Do not hardcode hex values here.

# Align a per-sample vector to a matrix. Named vectors are reordered by sample
# ID; unnamed vectors retain positional semantics for backward compatibility.
.align_sample_vector <- function(x, sample_names, arg = "group") {
  if (is.null(sample_names) || anyNA(sample_names) || any(sample_names == "") ||
      anyDuplicated(sample_names)) {
    stop("Matrix sample names must be present, non-empty, and unique.")
  }
  if (length(x) != length(sample_names)) {
    stop(sprintf("`%s` must have length %d (one value per matrix column).",
                 arg, length(sample_names)))
  }
  if (!is.null(names(x))) {
    x_names <- names(x)
    if (anyNA(x_names) || any(x_names == "") || anyDuplicated(x_names)) {
      stop(sprintf("Named `%s` must have non-empty, unique sample names.", arg))
    }
    missing_samples <- setdiff(sample_names, x_names)
    extra_samples <- setdiff(x_names, sample_names)
    if (length(missing_samples) > 0 || length(extra_samples) > 0) {
      stop(sprintf(
        "Named `%s` does not match matrix samples (missing: %s; extra: %s).",
        arg,
        if (length(missing_samples)) paste(missing_samples, collapse = ", ") else "none",
        if (length(extra_samples)) paste(extra_samples, collapse = ", ") else "none"
      ))
    }
    x <- x[sample_names]
  }
  if (anyNA(x)) stop(sprintf("`%s` contains missing values.", arg))
  x
}

.validate_feature_matrix <- function(mat, label = "expr") {
  mat <- as.matrix(mat)
  if (!is.numeric(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    stop(sprintf("`%s` must be a non-empty numeric feature-by-sample matrix.", label))
  }
  if (is.null(rownames(mat)) || anyNA(rownames(mat)) || any(rownames(mat) == "") ||
      anyDuplicated(rownames(mat))) {
    stop(sprintf("`%s` feature names must be present, non-empty, and unique.", label))
  }
  if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(colnames(mat) == "") ||
      anyDuplicated(colnames(mat))) {
    stop(sprintf("`%s` sample names must be present, non-empty, and unique.", label))
  }
  if (any(!is.finite(mat))) stop(sprintf("`%s` contains non-finite values.", label))
  mat
}

# Score custom gene sets per sample with the GSVA package. `expr` is a numeric
# matrix (features x samples, feature names as rownames); `gene_sets` is a named
# list mapping set name -> feature vector. Returns a gene_set x sample matrix.
score_gene_sets <- function(expr, gene_sets,
                            method = c("gsva", "ssgsea", "zscore", "plage"),
                            kcdf = "Gaussian", min_size = 5, verbose = FALSE) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package 'GSVA' is required for gene-set scoring.")
  }
  method <- match.arg(method)
  expr <- .validate_feature_matrix(expr, "expr")
  if (!is.list(gene_sets) || length(gene_sets) == 0 || is.null(names(gene_sets)) ||
      anyNA(names(gene_sets)) || any(names(gene_sets) == "") || anyDuplicated(names(gene_sets))) {
    stop("`gene_sets` must be a non-empty named list with unique set names.")
  }
  if (length(min_size) != 1 || !is.finite(min_size) || min_size < 1 || min_size != as.integer(min_size)) {
    stop("`min_size` must be one positive integer.")
  }
  gene_sets <- lapply(gene_sets, function(x) unique(as.character(stats::na.omit(x))))
  overlap_n <- vapply(gene_sets, function(x) sum(x %in% rownames(expr)), integer(1))
  audit <- data.frame(
    gene_set = names(gene_sets),
    input_size = vapply(gene_sets, length, integer(1)),
    overlap_size = overlap_n,
    overlap_fraction = overlap_n / pmax(vapply(gene_sets, length, integer(1)), 1L),
    retained = overlap_n >= min_size,
    stringsAsFactors = FALSE
  )
  gene_sets <- gene_sets[audit$retained]
  if (length(gene_sets) == 0) {
    stop(sprintf("No gene set has at least %d features overlapping `expr`.", min_size))
  }
  param <- switch(method,
    gsva   = GSVA::gsvaParam(expr, gene_sets, kcdf = kcdf, minSize = min_size),
    ssgsea = GSVA::ssgseaParam(expr, gene_sets, minSize = min_size),
    zscore = GSVA::zscoreParam(expr, gene_sets, minSize = min_size),
    plage  = GSVA::plageParam(expr, gene_sets, minSize = min_size)
  )
  scores <- GSVA::gsva(param, verbose = verbose)
  attr(scores, "gene_set_audit") <- audit
  scores
}

# Compare gene-set scores between groups for every set. `score_mat` is a
# gene_set x sample matrix; `group` is a per-sample grouping vector aligned to
# colnames(score_mat). `comparisons` is a list of c(control, treatment) pairs, so
# delta = mean(treatment) - mean(control) and Comparison = "treatment_vs_control".
# Per-comparison statistics reuse pairwise_effect_table(); the p-value adjustment
# is applied once across ALL set x comparison tests (global), matching the
# convention of a single BH correction over the whole pathway panel.
pathway_group_comparison <- function(score_mat, group, group_levels = NULL,
                                     comparisons = NULL, method = "t.test",
                                     p_adjust_method = "BH",
                                     exact_permutation = FALSE,
                                     include_effect_size = FALSE,
                                     pair_id = NULL) {
  score_mat <- .validate_feature_matrix(score_mat, "score_mat")
  group <- .align_sample_vector(group, colnames(score_mat), "group")
  if (!is.null(pair_id)) {
    pair_id <- .align_sample_vector(pair_id, colnames(score_mat), "pair_id")
  }
  if (!is.null(group_levels) && (anyNA(group_levels) || anyDuplicated(group_levels))) {
    stop("`group_levels` must contain unique, non-missing values.")
  }
  group <- factor(group, levels = group_levels %||% unique(group))
  if (anyNA(group)) stop("`group` contains values absent from `group_levels`.")
  long <- do.call(rbind, lapply(rownames(score_mat), function(pw) {
    df <- data.frame(Pathway = pw, sample = colnames(score_mat),
                     score = as.numeric(score_mat[pw, ]),
                     condition = group, stringsAsFactors = FALSE)
    if (!is.null(pair_id)) df$pair <- as.character(pair_id)
    df
  }))
  rows <- list()
  for (pw in unique(long$Pathway)) {
    sub <- long[long$Pathway == pw, , drop = FALSE]
    tab <- pairwise_effect_table(sub, value_col = "score", group_col = "condition",
                                 comparisons = comparisons, method = method,
                                 p_adjust_method = "none",
                                 pair_col = if (is.null(pair_id)) NULL else "pair")
    if (nrow(tab) == 0) next
    tab$Pathway <- pw
    rows[[length(rows) + 1]] <- tab
  }
  if (length(rows) == 0) return(data.frame())
  out <- dplyr::bind_rows(rows)
  out$p.adj <- stats::p.adjust(out$p, method = p_adjust_method)
  out$Comparison <- paste0(out$group2, "_vs_", out$group1)
  out$direction <- ifelse(out$delta > 0, "UP",
                          ifelse(out$delta < 0, "DOWN", "NO_CHANGE"))
  keep <- c("Pathway", "Comparison", "group1", "group2", "mean1", "mean2",
            "delta", "p", "p.adj", "direction")

  if (isTRUE(exact_permutation) || isTRUE(include_effect_size)) {
    robust_rows <- lapply(seq_len(nrow(out)), function(i) {
      pw <- out$Pathway[[i]]
      g1 <- out$group1[[i]]
      g2 <- out$group2[[i]]
      x1 <- as.numeric(score_mat[pw, group == g1])
      x2 <- as.numeric(score_mat[pw, group == g2])
      n1 <- length(x1); n2 <- length(x2)
      sd1 <- if (n1 > 1) stats::sd(x1) else NA_real_
      sd2 <- if (n2 > 1) stats::sd(x2) else NA_real_
      pooled_num <- (n1 - 1) * sd1^2 + (n2 - 1) * sd2^2
      pooled_den <- n1 + n2 - 2
      pooled_sd <- if (pooled_den > 0 && is.finite(pooled_num) && pooled_num > 0) {
        sqrt(pooled_num / pooled_den)
      } else {
        NA_real_
      }
      cohen_d <- if (is.finite(pooled_sd) && pooled_sd > 0) (mean(x2) - mean(x1)) / pooled_sd else NA_real_
      correction <- if ((n1 + n2) > 3) 1 - 3 / (4 * (n1 + n2) - 9) else NA_real_
      hedges_g <- cohen_d * correction

      perm_p <- NA_real_
      n_perm <- NA_integer_
      if (isTRUE(exact_permutation)) {
        values <- c(x1, x2)
        assignments <- utils::combn(seq_along(values), n2)
        perm_delta <- apply(assignments, 2, function(idx2) {
          mean(values[idx2]) - mean(values[-idx2])
        })
        observed <- mean(x2) - mean(x1)
        tolerance <- sqrt(.Machine$double.eps)
        perm_p <- mean(abs(perm_delta) + tolerance >= abs(observed))
        n_perm <- ncol(assignments)
      }
      data.frame(n1 = n1, n2 = n2, sd1 = sd1, sd2 = sd2,
                 hedges_g = hedges_g, p.permutation = perm_p,
                 permutation_count = n_perm)
    })
    robust <- dplyr::bind_rows(robust_rows)
    out <- cbind(out, robust)
    if (isTRUE(exact_permutation)) {
      out$p.permutation.adj <- stats::p.adjust(out$p.permutation, method = p_adjust_method)
    }
    keep <- c(keep, "n1", "n2", "sd1", "sd2", "hedges_g")
    if (isTRUE(exact_permutation)) {
      keep <- c(keep, "p.permutation", "p.permutation.adj", "permutation_count")
    }
  }
  out[, keep]
}

# Read and validate a versioned gene-set registry. The registry is intentionally
# richer than a two-column gene list so project reports can preserve biological
# question, source, direction and version metadata alongside the scored genes.
read_gene_set_registry <- function(path,
                                   set_col = "gene_set",
                                   gene_col = "gene",
                                   required_metadata = c("theme", "version", "source_type", "source_id")) {
  if (!file.exists(path)) stop("Gene-set registry not found: ", path)
  registry <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- unique(c(set_col, gene_col, required_metadata))
  missing_cols <- setdiff(required, colnames(registry))
  if (length(missing_cols) > 0) {
    stop("Gene-set registry is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  registry[[set_col]] <- trimws(as.character(registry[[set_col]]))
  registry[[gene_col]] <- trimws(as.character(registry[[gene_col]]))
  registry <- registry[!is.na(registry[[set_col]]) & nzchar(registry[[set_col]]) &
                         !is.na(registry[[gene_col]]) & nzchar(registry[[gene_col]]), , drop = FALSE]
  registry <- registry[!duplicated(registry[, c(set_col, gene_col), drop = FALSE]), , drop = FALSE]
  if (nrow(registry) == 0) stop("Gene-set registry contains no valid set-gene rows.")
  for (column in required_metadata) {
    if (anyNA(registry[[column]]) || any(!nzchar(trimws(as.character(registry[[column]]))))) {
      stop("Gene-set registry metadata column contains missing/blank values: ", column)
    }
  }
  gene_sets <- split(registry[[gene_col]], registry[[set_col]])
  attr(gene_sets, "registry") <- registry
  gene_sets
}

# Summarize registry content and expression overlap without scoring. This is a
# content/provenance audit, not evidence that a gene set is biologically valid.
audit_gene_set_registry <- function(gene_sets, expression_features = NULL, min_size = 5) {
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || anyDuplicated(names(gene_sets))) {
    stop("`gene_sets` must be a uniquely named list.")
  }
  registry <- attr(gene_sets, "registry")
  rows <- lapply(names(gene_sets), function(set_name) {
    genes <- unique(as.character(gene_sets[[set_name]]))
    overlap <- if (is.null(expression_features)) NA_integer_ else sum(genes %in% expression_features)
    metadata <- if (!is.null(registry)) registry[registry$gene_set == set_name, , drop = FALSE] else NULL
    data.frame(
      gene_set = set_name,
      theme = if (!is.null(metadata) && "theme" %in% names(metadata)) metadata$theme[[1]] else NA_character_,
      version = if (!is.null(metadata) && "version" %in% names(metadata)) metadata$version[[1]] else NA_character_,
      source_type = if (!is.null(metadata) && "source_type" %in% names(metadata)) metadata$source_type[[1]] else NA_character_,
      source_id = if (!is.null(metadata) && "source_id" %in% names(metadata)) metadata$source_id[[1]] else NA_character_,
      input_size = length(genes),
      overlap_size = overlap,
      overlap_fraction = if (is.na(overlap)) NA_real_ else overlap / max(length(genes), 1L),
      retained = if (is.na(overlap)) length(genes) >= min_size else overlap >= min_size,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

# Visualize every module x contrast result in one compact audit matrix. Tile
# color encodes the estimated score difference; text keeps the parametric and
# exact-permutation multiplicity results visible without implying that either
# is a separate independent analysis.
plot_pathway_sensitivity_matrix_pdf <- function(
    comp_df, filename,
    pathway_col = "Pathway", comparison_col = "Comparison",
    delta_col = "delta", param_padj_col = "p.adj",
    permutation_padj_col = "p.permutation.adj",
    pathway_levels = NULL,
    title = "Gene-set score sensitivity across comparisons",
    subtitle = "Tile = delta score; labels = parametric BH / exact-permutation BH",
    width = 9.5, height = 8.5) {
  required <- c(pathway_col, comparison_col, delta_col,
                param_padj_col, permutation_padj_col)
  missing_cols <- setdiff(required, colnames(comp_df))
  if (length(missing_cols)) {
    stop("Missing pathway sensitivity columns: ", paste(missing_cols, collapse = ", "))
  }
  df <- comp_df[is.finite(comp_df[[delta_col]]), , drop = FALSE]
  if (!nrow(df)) return(invisible(NULL))
  key <- paste(df[[pathway_col]], df[[comparison_col]], sep = "\r")
  if (anyDuplicated(key)) stop("`comp_df` contains duplicate pathway-by-comparison rows.")

  fmt_p <- function(x) {
    ifelse(is.na(x), "NA",
           ifelse(x < 0.001, formatC(x, format = "e", digits = 1),
                  formatC(x, format = "f", digits = 3)))
  }
  df$AuditLabel <- paste0(
    "BH ", fmt_p(df[[param_padj_col]]), "\nperm ",
    fmt_p(df[[permutation_padj_col]])
  )
  if (is.null(pathway_levels)) pathway_levels <- unique(df[[pathway_col]])
  df[[pathway_col]] <- factor(df[[pathway_col]], levels = rev(pathway_levels))
  df[[comparison_col]] <- factor(df[[comparison_col]], levels = unique(df[[comparison_col]]))
  cap <- max(abs(df[[delta_col]]), na.rm = TRUE)
  if (!is.finite(cap) || cap == 0) cap <- 1

  p <- ggplot2::ggplot(
    df, ggplot2::aes(x = .data[[comparison_col]], y = .data[[pathway_col]],
                     fill = .data[[delta_col]])
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = .data$AuditLabel),
                       size = 2.75, lineheight = 0.9, color = "#202124") +
    ggplot2::scale_fill_gradient2(
      low = RNAseq_PALETTE$down_regulated,
      mid = "#F7F7F7",
      high = RNAseq_PALETTE$up_regulated,
      midpoint = 0, limits = c(-cap, cap), name = "Delta score"
    ) +
    ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_publication(base_size = 8) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 7, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 6.5, color = "#202124"),
      plot.subtitle = ggplot2::element_text(size = 7, color = "#5F6368"),
      legend.position = "right"
    )
  save_pdf_plot(p, filename, width = width, height = height)
  invisible(p)
}

# Registry QC separates provenance/coverage from biological inference. Bar
# length is the fraction of registry genes found in the expression matrix;
# direct labels preserve the input and overlap denominators.
plot_gene_set_registry_qc_pdf <- function(
    audit_df, filename,
    set_col = "gene_set", source_col = "source_type",
    input_col = "input_size", overlap_col = "overlap_size",
    fraction_col = "overlap_fraction",
    title = "Custom gene-set registry coverage",
    subtitle = "Bar = expressed-gene overlap fraction; label = overlap / registry size",
    source_colors = NULL, width = 9, height = 7.5) {
  required <- c(set_col, source_col, input_col, overlap_col)
  missing_cols <- setdiff(required, colnames(audit_df))
  if (length(missing_cols)) {
    stop("Missing registry QC columns: ", paste(missing_cols, collapse = ", "))
  }
  df <- audit_df
  if (!fraction_col %in% colnames(df)) {
    df[[fraction_col]] <- df[[overlap_col]] / pmax(df[[input_col]], 1)
  }
  df <- df[is.finite(df[[fraction_col]]), , drop = FALSE]
  if (!nrow(df)) return(invisible(NULL))
  df$CoverageLabel <- paste0(df[[overlap_col]], " / ", df[[input_col]])
  ord <- order(df[[fraction_col]], df[[overlap_col]], decreasing = FALSE)
  df[[set_col]] <- factor(df[[set_col]], levels = df[[set_col]][ord])
  sources <- unique(as.character(df[[source_col]]))
  if (is.null(source_colors)) {
    palette <- rep(RNAseq_PALETTE$group_three, length.out = length(sources))
    source_colors <- stats::setNames(palette, sources)
  }

  p <- ggplot2::ggplot(
    df, ggplot2::aes(x = .data[[fraction_col]], y = .data[[set_col]],
                     fill = .data[[source_col]])
  ) +
    ggplot2::geom_col(width = 0.72, color = "#333333", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = .data$CoverageLabel),
      hjust = 1.08, size = 3.0, color = "white", fontface = "bold"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), breaks = seq(0, 1, 0.25),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_fill_manual(values = source_colors, name = "Source") +
    ggplot2::labs(x = "Expressed-gene overlap", y = NULL,
                  title = title, subtitle = subtitle) +
    theme_publication(base_size = 8) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6.8, color = "#202124"),
      plot.subtitle = ggplot2::element_text(size = 7, color = "#5F6368"),
      legend.position = "top"
    )
  save_pdf_plot(p, filename, width = width, height = height)
  invisible(p)
}

# Bidirectional bar chart summarising delta scores across pathways and
# comparisons in a single figure. One panel per comparison; bars are filled by
# direction and annotated with significance stars. `comp_df` is the output of
# pathway_group_comparison() (or any table with pathway/comparison/delta/p.adj).
plot_pathway_delta_summary_pdf <- function(comp_df, filename,
                                           pathway_col = "Pathway",
                                           comparison_col = "Comparison",
                                           delta_col = "delta", padj_col = "p.adj",
                                           title = "Pathway activity change across comparisons",
                                           up_color = RNAseq_PALETTE$up_regulated,
                                           down_color = RNAseq_PALETTE$down_regulated,
                                           neutral_color = RNAseq_PALETTE$not_significant,
                                           show_ns = FALSE,
                                           width = 7.2, height = 4.5) {
  df <- comp_df
  if (nrow(df) == 0) return(invisible(NULL))
  required <- c(pathway_col, comparison_col, delta_col, padj_col)
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing pathway summary columns: ", paste(missing_cols, collapse = ", "))
  }
  df <- df[is.finite(df[[delta_col]]), , drop = FALSE]
  if (nrow(df) == 0) return(invisible(NULL))
  df$direction <- ifelse(df[[delta_col]] > 0, "UP",
                         ifelse(df[[delta_col]] < 0, "DOWN", "NO_CHANGE"))
  df$sig <- format_p_stars(df[[padj_col]])
  if (!isTRUE(show_ns)) df$sig[df$sig == "ns"] <- ""
  df[[pathway_col]] <- factor(df[[pathway_col]], levels = unique(df[[pathway_col]]))
  df[[comparison_col]] <- factor(df[[comparison_col]], levels = unique(df[[comparison_col]]))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[delta_col]], y = .data[[pathway_col]],
                                        fill = .data$direction)) +
    ggplot2::geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = .data$sig,
                                    hjust = ifelse(.data[[delta_col]] > 0, -0.2, 1.2)), size = 2.6) +
    ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
    ggplot2::facet_wrap(
      ~ .data[[comparison_col]], nrow = 1,
      labeller = ggplot2::labeller(.default = function(x) gsub("_", " ", x, fixed = TRUE))
    ) +
    ggplot2::scale_fill_manual(values = c(UP = up_color, DOWN = down_color,
                                         NO_CHANGE = neutral_color), drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.18))) +
    ggplot2::scale_y_discrete(labels = function(x) gsub("_", " ", x, fixed = TRUE)) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_publication(base_size = 8) +
    ggplot2::labs(x = "Delta score (treatment - control)", y = NULL,
                  title = title, fill = "Direction") +
    ggplot2::theme(legend.position = "top",
                   plot.margin = ggplot2::margin(8, 18, 8, 12))
  save_pdf_plot(p, filename, width = width, height = height)
  invisible(p)
}

# Heatmap of log2 fold-changes for a set of key genes across comparisons.
# `fc_long` is a long data frame with one row per gene x comparison. `row_group`
# optionally supplies a per-gene grouping (e.g. pathway) used to order, annotate
# and split rows. Unlike plot_expression_heatmap_pdf this does NOT row z-score
# (log2FC values are already comparable) and annotates rows instead of columns.
plot_keygenes_log2fc_heatmap_pdf <- function(fc_long, filename,
                                             gene_col = "Gene",
                                             comparison_col = "Comparison",
                                             lfc_col = "log2FC",
                                             row_group = NULL, row_group_levels = NULL,
                                             lfc_cap = 2, cluster_rows = FALSE,
                                             cluster_cols = FALSE,
                                             title = "Key genes log2 fold-change",
                                             width = 5.8, height = NULL,
                                             row_font_size = 7.5, column_font_size = 9,
                                             na_color = "#E5E5E5") {
  if (nrow(fc_long) == 0) return(invisible(NULL))
  required <- c(gene_col, comparison_col, lfc_col)
  missing_cols <- setdiff(required, colnames(fc_long))
  if (length(missing_cols) > 0) {
    stop("Missing heatmap columns: ", paste(missing_cols, collapse = ", "))
  }
  key <- paste(fc_long[[gene_col]], fc_long[[comparison_col]], sep = "\r")
  if (anyDuplicated(key)) {
    stop("`fc_long` contains duplicate gene-by-comparison rows; summarize or resolve them first.")
  }
  wide <- tidyr::pivot_wider(fc_long, id_cols = dplyr::all_of(gene_col),
                             names_from = dplyr::all_of(comparison_col),
                             values_from = dplyr::all_of(lfc_col))
  m <- as.matrix(as.data.frame(wide[, -1, drop = FALSE]))
  rownames(m) <- wide[[gene_col]]

  # Order and optionally split rows by the supplied per-gene grouping.
  left_annotation <- NULL
  row_split <- NULL
  if (!is.null(row_group)) {
    if (length(row_group) != nrow(fc_long)) {
      stop("`row_group` must have one value per row of `fc_long`.")
    }
    group_map <- unique(data.frame(
      gene = as.character(fc_long[[gene_col]]),
      group = as.character(row_group),
      stringsAsFactors = FALSE
    ))
    if (anyDuplicated(group_map$gene)) {
      stop("Each gene must map to exactly one `row_group`.")
    }
    grp <- stats::setNames(group_map$group, group_map$gene)[rownames(m)]
    if (anyNA(grp)) stop("`row_group` is missing for one or more plotted genes.")
    grp <- factor(grp, levels = row_group_levels %||% unique(grp))
    ord <- order(grp)
    m <- m[ord, , drop = FALSE]
    grp <- grp[ord]
    row_split <- grp
    pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(levels(grp)))
    left_annotation <- ComplexHeatmap::rowAnnotation(
      Group = grp,
      col = list(Group = stats::setNames(pal, levels(grp))),
      show_annotation_name = FALSE,
      simple_anno_size = grid::unit(3, "mm")
    )
  }

  if (is.null(height)) height <- max(3.5, 0.24 * nrow(m) + 1.5)
  col_fun <- circlize::colorRamp2(
    c(-lfc_cap, 0, lfc_cap),
    c(RNAseq_PALETTE$heatmap_down, RNAseq_PALETTE$heatmap_mid, RNAseq_PALETTE$heatmap_up)
  )
  ht <- ComplexHeatmap::Heatmap(
    m, name = "log2FC", col = col_fun, na_col = na_color,
    left_annotation = left_annotation, row_split = row_split,
    cluster_rows = cluster_rows, cluster_columns = cluster_cols,
    show_row_names = TRUE, show_column_names = TRUE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    row_names_gp = grid::gpar(fontsize = row_font_size),
    column_names_gp = grid::gpar(fontsize = column_font_size),
    row_title_side = "left", row_title_rot = 0,
    row_title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
    row_gap = grid::unit(2.5, "mm"),
    column_title = title
  )
  save_pdf_device(filename, width = width, height = height, draw = function() {
    ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  })
  invisible(ht)
}

# Melt selected genes of an expression matrix to long format for grouped plots.
# Returns sample / Gene / value / condition. `group` aligns to colnames(expr).
melt_gene_expression <- function(expr, genes, group) {
  expr <- .validate_feature_matrix(expr, "expr")
  group <- .align_sample_vector(group, colnames(expr), "group")
  genes <- intersect(genes, rownames(expr))
  if (length(genes) == 0) return(data.frame())
  do.call(rbind, lapply(genes, function(g) {
    data.frame(Gene = g, sample = colnames(expr),
               value = as.numeric(expr[g, ]),
               condition = as.character(group), stringsAsFactors = FALSE)
  }))
}

# Per-gene grouped expression plots for a set of key genes. With facet = TRUE a
# single faceted PDF is produced via plot_group_boxplot_pdf(); with facet = FALSE
# one file per gene is written using `filename` as a prefix, and `plot` selects
# the boxplot or violin+box style. Colors always come from `group_colors`.
plot_gene_expression_pdf <- function(expr, genes, group, group_levels, group_colors,
                                     filename, plot = c("boxplot", "violin"),
                                     facet = TRUE, comparisons = NULL,
                                     method = "t.test", p_adjust_method = "BH",
                                     title = NULL, ylab = "Expression",
                                     width = NULL, height = 6,
                                     label_style = c("stars", "full")) {
  plot <- match.arg(plot)
  label_style <- match.arg(label_style)
  long <- melt_gene_expression(expr, genes, group)
  if (nrow(long) == 0) return(invisible(NULL))
  long$condition <- factor(long$condition, levels = group_levels)
  if (isTRUE(facet)) {
    n_genes <- length(unique(long$Gene))
    if (is.null(width)) width <- max(7, ceiling(n_genes / 2) * 2.4)
    p <- plot_group_boxplot_pdf(long, value_col = "value", group_col = "condition",
                                filename = filename, facet_col = "Gene",
                                comparisons = comparisons, method = method,
                                p_adjust_method = p_adjust_method, title = title,
                                ylab = ylab, group_colors = group_colors,
                                width = width, height = height,
                                label_style = label_style)
    return(invisible(p))
  }
  base <- sub("\\.pdf$", "", filename)
  for (g in unique(long$Gene)) {
    sub <- long[long$Gene == g, , drop = FALSE]
    fn <- paste0(base, "_", g, ".pdf")
    if (plot == "violin") {
      plot_group_violin_boxplot_pdf(sub, value_col = "value", group_col = "condition",
                                    filename = fn, comparisons = comparisons,
                                    method = method, p_adjust_method = p_adjust_method,
                                    title = g, ylab = ylab, group_colors = group_colors,
                                    width = width %||% 5.8, height = height,
                                    label_style = label_style)
    } else {
      plot_group_boxplot_pdf(sub, value_col = "value", group_col = "condition",
                             filename = fn, comparisons = comparisons,
                             method = method, p_adjust_method = p_adjust_method,
                             title = g, ylab = ylab, group_colors = group_colors,
                             width = width %||% 5, height = height,
                             label_style = label_style)
    }
  }
  invisible(NULL)
}
