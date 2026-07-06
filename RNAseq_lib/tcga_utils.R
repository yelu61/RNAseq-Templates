# TCGA / public data mining helpers for bulk RNA-seq templates.

# Build a TCGAbiolinks GDC query for a given project.
build_tcga_query <- function(project = "TCGA-STAD",
                               data.category = "Transcriptome Profiling",
                               data.type = "Gene Expression Quantification",
                               workflow.type = "STAR - Counts") {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("Package 'TCGAbiolinks' is required for TCGA data download.")
  }
  TCGAbiolinks::GDCquery(
    project = project,
    data.category = data.category,
    data.type = data.type,
    workflow.type = workflow.type
  )
}

# Extract assay matrices from a GDC prepared SummarizedExperiment object.
extract_tcga_assays <- function(se, assay_names = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required.")
  }
  if (is.null(assay_names)) {
    assay_names <- SummarizedExperiment::assayNames(se)
  }
  assays <- lapply(assay_names, function(n) SummarizedExperiment::assay(se, i = n))
  names(assays) <- assay_names
  assays
}

# Convert ENSEMBL/gene_id matrix to gene-symbol matrix, de-duplicate by highest average expression.
symbolize_and_dedup <- function(mat, id_map,
                                  id_col = "gene_id", symbol_col = "gene_name",
                                  biotype_col = NULL, biotype_filter = NULL) {
  df <- as.data.frame(mat, check.names = FALSE)
  df[[id_col]] <- rownames(df)
  df <- dplyr::inner_join(df, id_map, by = id_col)

  if (!is.null(biotype_col) && biotype_col %in% colnames(df) && !is.null(biotype_filter)) {
    df <- df[df[[biotype_col]] == biotype_filter, , drop = FALSE]
  }

  expr_cols <- setdiff(colnames(df), c(id_col, symbol_col, biotype_col))
  df$mean_expr <- rowMeans(df[, expr_cols, drop = FALSE])
  df <- df[order(df$mean_expr, decreasing = TRUE), ]
  df <- df[!duplicated(df[[symbol_col]]), ]
  rownames(df) <- df[[symbol_col]]
  df[, c(symbol_col, expr_cols), drop = FALSE] |>
    dplyr::select(-dplyr::all_of(symbol_col))
}

# Build an id→symbol map from a SummarizedExperiment object's rowData.
build_id_map_from_se <- function(se, id_col = "gene_id", symbol_col = "gene_name", type_col = "gene_type") {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required.")
  }
  rd <- SummarizedExperiment::rowData(se)
  out_cols <- intersect(c(id_col, symbol_col, type_col), colnames(rd))
  if (length(out_cols) == 0) {
    stop("rowData does not contain expected columns: ", paste(c(id_col, symbol_col, type_col), collapse = ", "),
         "\nAvailable: ", paste(colnames(rd), collapse = ", "))
  }
  df <- as.data.frame(rd[, out_cols, drop = FALSE], stringsAsFactors = FALSE)
  if (symbol_col %in% colnames(df) && id_col %in% colnames(df)) {
    df <- df[!is.na(df[[id_col]]) & df[[id_col]] != "" &
                 !is.na(df[[symbol_col]]) & df[[symbol_col]] != "", , drop = FALSE]
  }
  df
}

# Extract and clean TCGA clinical metadata from a SummarizedExperiment object.
extract_tcga_clinical <- function(se) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required.")
  }
  cd <- SummarizedExperiment::colData(se)
  clinical <- as.data.frame(cd, stringsAsFactors = FALSE)
  # Coerce list columns to character to avoid downstream issues
  for (col in colnames(clinical)) {
    if (is.list(clinical[[col]])) {
      clinical[[col]] <- vapply(clinical[[col]], function(x) paste(unlist(x), collapse = ";"), character(1))
    }
  }
  clinical
}

# Infer Tumor/Normal grouping from TCGA barcode (sample type code in positions 14-15).
infer_tcga_tumor_normal <- function(barcodes, tumor_label = "Tumor", normal_label = "Normal") {
  code <- as.numeric(substr(barcodes, 14, 15))
  ifelse(code < 10, tumor_label, normal_label)
}

# Prepare a clean survival data frame from TCGA clinical metadata.
prepare_tcga_survival <- function(clinical,
                                    barcode_col = "barcode",
                                    vital_col = "vital_status",
                                    death_days_col = "days_to_death",
                                    followup_days_col = "days_to_last_follow_up",
                                    tissue_col = "tissue_type",
                                    tumor_only = TRUE) {
  required <- c(barcode_col, vital_col, death_days_col, followup_days_col)
  missing_cols <- setdiff(required, colnames(clinical))
  if (length(missing_cols) > 0) {
    stop("Missing survival columns: ", paste(missing_cols, collapse = ", "))
  }

  out <- clinical[, required, drop = FALSE]
  if (!is.null(tissue_col) && tissue_col %in% colnames(clinical)) {
    out[[tissue_col]] <- clinical[[tissue_col]]
  }

  # Convert NA strings / blanks to NA
  for (col in c(death_days_col, followup_days_col)) {
    out[[col]] <- suppressWarnings(as.numeric(as.character(out[[col]])))
  }
  out[[vital_col]] <- toupper(as.character(out[[vital_col]]))

  out$days_to_death <- out[[death_days_col]]
  out$days_to_last_follow_up <- out[[followup_days_col]]
  out$time <- pmax(out$days_to_death, out$days_to_last_follow_up, na.rm = TRUE)
  out$status <- ifelse(out[[vital_col]] == "DEAD", 1, 0)
  out <- out[!is.na(out$time) & out$time > 0, ]

  if (tumor_only && !is.null(tissue_col) && tissue_col %in% colnames(out)) {
    out <- out[toupper(out[[tissue_col]]) == "TUMOR", ]
  }
  out
}

# Plot Kaplan-Meier curves for clinical variables.
run_clinical_km <- function(surv_df, clinical_df, var_cols,
                            outdir = ".", time_unit = "month",
                            method = "t.test", show_ns = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE) || !requireNamespace("survminer", quietly = TRUE)) {
    stop("Packages 'survival' and 'survminer' are required for KM plots.")
  }
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  results <- data.frame(variable = character(), type = character(), n = integer(), stringsAsFactors = FALSE)

  for (var in var_cols) {
    if (!var %in% colnames(clinical_df)) {
      warning("Clinical variable not found: ", var)
      next
    }
    if (!"barcode" %in% colnames(clinical_df)) {
      stop("clinical_df must contain a 'barcode' column to match surv_df.")
    }

    tmp <- clinical_df[, c("barcode", var), drop = FALSE]
    tmp <- tmp[!is.na(tmp[[var]]) & tmp[[var]] != "", , drop = FALSE]
    merged <- merge(surv_df, tmp, by = "barcode", all.x = FALSE)
    if (nrow(merged) < 5) {
      warning("Too few samples with non-missing ", var, " for KM plot.")
      next
    }

    if (is.numeric(merged[[var]]) && dplyr::n_distinct(merged[[var]]) > 2) {
      merged[[var]] <- stratify_by_quantile(merged[[var]], n_groups = 2, labels = c("Low", "High"))
      type <- "continuous_split"
    } else {
      merged[[var]] <- factor(merged[[var]])
      type <- "categorical"
    }

    plot_km_by_group_pdf(
      merged,
      group_col = var,
      filename = file.path(outdir, paste0("KM_clinical_", var, ".pdf")),
      title = paste("KM survival by", var),
      time_unit = time_unit
    )

    results <- rbind(results, data.frame(variable = var, type = type, n = nrow(merged), stringsAsFactors = FALSE))
  }

  utils::write.csv(results, file.path(outdir, "clinical_KM_summary.csv"), row.names = FALSE)
  invisible(results)
}

# Plot Kaplan-Meier curve for a continuous variable split by median.
plot_km_by_median_pdf <- function(surv_df, value_col, filename,
                                   title = NULL, ylab = "Survival probability",
                                   xlab = "Time in months", time_unit = "month",
                                   xlim = NULL, ylim = c(0, 1), palette = c("#3c81ba", "#bb3835")) {
  if (!requireNamespace("survival", quietly = TRUE) || !requireNamespace("survminer", quietly = TRUE)) {
    stop("Packages 'survival' and 'survminer' are required for KM plots.")
  }
  if (!all(c("time", "status") %in% colnames(surv_df))) {
    stop("surv_df must contain 'time' and 'status' columns.")
  }
  if (!value_col %in% colnames(surv_df)) {
    stop("Value column not found: ", value_col)
  }

  plot_data <- surv_df[!is.na(surv_df[[value_col]]), ]
  if (nrow(plot_data) == 0) return(invisible(NULL))
  plot_data$Group <- ifelse(plot_data[[value_col]] > median(plot_data[[value_col]], na.rm = TRUE),
                            "High", "Low")

  if (time_unit == "month") {
    plot_data$time <- plot_data$time / 30.4375
  } else if (time_unit == "year") {
    plot_data$time <- plot_data$time / 365.25
  }

  fit <- survival::survfit(survival::Surv(time, status) ~ Group, data = plot_data)
  n_high <- sum(plot_data$Group == "High")
  n_low <- sum(plot_data$Group == "Low")

  p <- survminer::ggsurvplot(
    fit,
    data = plot_data,
    title = title %||% paste(value_col, "survival"),
    xlab = xlab,
    ylab = ylab,
    legend.title = value_col,
    legend.labs = c(paste0("High (N=", n_high, ")"), paste0("Low (N=", n_low, ")")),
    palette = palette,
    pval = TRUE,
    pval.size = 5,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = FALSE,
    surv.scale = "percent",
    ggtheme = ggplot2::theme_classic(base_size = 12),
    xlim = xlim,
    ylim = ylim
  )

  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(filename, width = 7, height = 8)
  print(p, newpage = FALSE)
  grDevices::dev.off()
  invisible(p)
}

# Plot expression of a single gene by group (tumor/normal or user-defined).
plot_tcga_gene_boxplot_pdf <- function(expr_df, gene, group_vec, filename,
                                       title = NULL, ylab = "log2(TPM+1)",
                                       method = "t.test", show_ns = FALSE,
                                       width = 6, height = 6) {
  if (!gene %in% rownames(expr_df)) {
    warning("Gene not found in expression matrix: ", gene)
    return(invisible(NULL))
  }
  plot_data <- data.frame(
    sample = colnames(expr_df),
    expression = as.numeric(expr_df[gene, ]),
    group = factor(group_vec)
  )
  comparisons <- utils::combn(levels(plot_data$group), 2, simplify = FALSE)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = expression, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.55) +
    ggplot2::geom_jitter(width = 0.12, size = 2, shape = 21, color = "#222222", stroke = 0.25, fill = "white", alpha = 0.9) +
    ggplot2::labs(x = NULL, y = ylab, title = title %||% gene) +
    theme_publication(base_size = 12) +
    ggplot2::theme(legend.position = "none")
  if (length(comparisons) > 0) {
    p <- p + ggpubr::stat_compare_means(
      comparisons = comparisons,
      method = method,
      hide.ns = !show_ns,
      tip.length = 0.015,
      bracket.size = 0.4
    )
  }
  save_pdf_plot(p, filename, width = width, height = height)
  p
}
