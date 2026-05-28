# Input/output helpers for bulk RNA-seq templates.

read_count_table <- function(input_file, input_format = c("tsv", "csv", "excel")) {
  input_format <- match.arg(input_format)
  if (input_format == "excel") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for Excel input.")
    }
    openxlsx::read.xlsx(input_file)
  } else if (input_format == "csv") {
    readr::read_csv(input_file, show_col_types = FALSE)
  } else {
    readr::read_tsv(input_file, show_col_types = FALSE)
  }
}

summarize_raw_input <- function(rawcount, gene_name_col, biotype_col = NULL, biotype_filter = NULL) {
  gene_values <- rawcount[[gene_name_col]]
  data.frame(
    step = c("raw_rows", "missing_gene_symbol", "duplicated_gene_symbol", "target_biotype_rows"),
    value = c(
      nrow(rawcount),
      sum(is.na(gene_values) | gene_values == ""),
      sum(duplicated(gene_values[!is.na(gene_values) & gene_values != ""])),
      if (!is.null(biotype_col) && biotype_col %in% colnames(rawcount) && !is.null(biotype_filter)) {
        sum(rawcount[[biotype_col]] == biotype_filter, na.rm = TRUE)
      } else {
        NA_integer_
      }
    )
  )
}

detect_count_columns <- function(rawcount, gene_name_col, count_cols = NULL, annotation_cols = NULL) {
  if (is.null(annotation_cols)) {
    annotation_cols <- c(
      "Gene", "gene_id", "gene_name", "gene_chr", "gene_start", "gene_end",
      "gene_strand", "gene_biotype", "gene_description", "external_gene_name",
      "hgnc_symbol", "chromosome_name", "start_position", "end_position",
      "length", "gene_length", "gene_type"
    )
  }

  if (is.null(count_cols)) {
    count_col_names <- setdiff(colnames(rawcount), annotation_cols)
    count_col_names <- setdiff(count_col_names, gene_name_col)
  } else if (is.numeric(count_cols)) {
    count_col_names <- colnames(rawcount)[count_cols]
  } else {
    count_col_names <- count_cols
  }

  missing_cols <- setdiff(c(gene_name_col, count_col_names), colnames(rawcount))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  count_col_names
}

validate_sample_design <- function(sample_names, groups, group_levels, comparisons, count_col_names) {
  if (length(sample_names) != length(count_col_names)) {
    stop("SAMPLE_NAMES length (", length(sample_names), ") must equal count columns (", length(count_col_names), ").")
  }
  if (length(groups) != length(sample_names)) {
    stop("GROUPS length (", length(groups), ") must equal SAMPLE_NAMES length (", length(sample_names), ").")
  }
  if (!all(groups %in% group_levels)) {
    stop("GROUPS contains values not present in GROUP_LEVELS: ", paste(setdiff(unique(groups), group_levels), collapse = ", "))
  }
  comparison_groups <- unique(unlist(lapply(comparisons, function(x) x[2:3])))
  if (!all(comparison_groups %in% group_levels)) {
    stop("COMPARISONS contains group names not present in GROUP_LEVELS: ", paste(setdiff(comparison_groups, group_levels), collapse = ", "))
  }
  invisible(TRUE)
}

calculate_sample_qc <- function(count_data, col_data = NULL, group_col = "condition") {
  qc <- data.frame(
    sample = colnames(count_data),
    library_size = colSums(count_data),
    detected_genes = colSums(count_data > 0),
    zero_fraction = colMeans(count_data == 0),
    stringsAsFactors = FALSE
  )

  if (!is.null(col_data) && group_col %in% colnames(col_data)) {
    qc$group <- as.character(col_data[qc$sample, group_col])
  }

  log_counts <- log2(as.matrix(count_data) + 1)
  sample_cor <- stats::cor(log_counts, method = "spearman", use = "pairwise.complete.obs")
  diag(sample_cor) <- NA
  qc$median_sample_correlation <- apply(sample_cor, 2, stats::median, na.rm = TRUE)
  qc$library_size_z <- as.numeric(scale(log10(qc$library_size + 1)))
  qc$detected_genes_z <- as.numeric(scale(qc$detected_genes))
  qc$median_correlation_z <- as.numeric(scale(qc$median_sample_correlation))

  attr(qc, "sample_cor") <- sample_cor
  qc
}

flag_sample_qc <- function(sample_qc,
                           min_library_size = NULL,
                           min_detected_genes = NULL,
                           max_zero_fraction = NULL,
                           min_median_correlation_z = -3) {
  qc <- sample_qc
  qc$qc_flag <- FALSE
  reasons <- vector("list", nrow(qc))

  add_reason <- function(mask, reason) {
    qc$qc_flag[mask] <<- TRUE
    reasons[mask] <<- Map(c, reasons[mask], reason)
  }

  if (!is.null(min_library_size)) {
    add_reason(qc$library_size < min_library_size, paste0("library_size<", min_library_size))
  }
  if (!is.null(min_detected_genes)) {
    add_reason(qc$detected_genes < min_detected_genes, paste0("detected_genes<", min_detected_genes))
  }
  if (!is.null(max_zero_fraction)) {
    add_reason(qc$zero_fraction > max_zero_fraction, paste0("zero_fraction>", max_zero_fraction))
  }
  if (!is.null(min_median_correlation_z)) {
    add_reason(qc$median_correlation_z < min_median_correlation_z, paste0("median_correlation_z<", min_median_correlation_z))
  }

  qc$qc_reason <- vapply(reasons, function(x) {
    x <- unique(unlist(x))
    if (length(x) == 0) "" else paste(x, collapse = ";")
  }, character(1))
  qc
}

apply_sample_exclusion <- function(count_data, col_data, sample_exclude = character(0)) {
  sample_exclude <- intersect(sample_exclude, colnames(count_data))
  keep_samples <- setdiff(colnames(count_data), sample_exclude)
  if (length(keep_samples) < 2) {
    stop("Fewer than 2 samples remain after applying SAMPLE_EXCLUDE.")
  }
  count_data <- count_data[, keep_samples, drop = FALSE]
  col_data <- col_data[keep_samples, , drop = FALSE]
  list(count_data = count_data, col_data = col_data, excluded_samples = sample_exclude)
}

build_count_matrix <- function(rawcount, gene_name_col, count_col_names, sample_names, duplicate_report_file = NULL) {
  count_data <- as.data.frame(rawcount[, c(gene_name_col, count_col_names)], check.names = FALSE)
  count_data[, count_col_names] <- lapply(count_data[, count_col_names, drop = FALSE], function(x) as.numeric(as.character(x)))

  if (anyNA(count_data[, count_col_names])) {
    stop("Count columns contain non-numeric or missing values after conversion.")
  }
  if (any(count_data[, count_col_names] < 0)) {
    stop("Count matrix contains negative values.")
  }
  if (any(count_data[, count_col_names] %% 1 != 0)) {
    stop("DESeq2 requires raw integer counts; non-integer values detected.")
  }

  count_data[, count_col_names] <- lapply(count_data[, count_col_names, drop = FALSE], as.integer)
  count_data <- count_data[!is.na(count_data[[gene_name_col]]) & count_data[[gene_name_col]] != "", ]

  expression_cols <- setdiff(names(count_data), gene_name_col)
  count_data$mean_count_for_dedup <- rowMeans(count_data[, expression_cols, drop = FALSE])
  duplicated_symbols <- count_data[[gene_name_col]][duplicated(count_data[[gene_name_col]]) | duplicated(count_data[[gene_name_col]], fromLast = TRUE)]
  if (length(duplicated_symbols) > 0 && !is.null(duplicate_report_file)) {
    duplicate_report <- count_data[count_data[[gene_name_col]] %in% duplicated_symbols, c(gene_name_col, "mean_count_for_dedup", expression_cols), drop = FALSE]
    duplicate_report$kept <- ave(
      duplicate_report$mean_count_for_dedup,
      duplicate_report[[gene_name_col]],
      FUN = function(x) x == max(x, na.rm = TRUE)
    )
    dir.create(dirname(duplicate_report_file), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(duplicate_report, duplicate_report_file, row.names = FALSE)
  }

  index <- order(rowMeans(count_data[, expression_cols, drop = FALSE]), decreasing = TRUE)
  count_data <- count_data[index, ]
  count_data <- count_data[!duplicated(count_data[[gene_name_col]]), ]
  rownames(count_data) <- count_data[[gene_name_col]]
  count_data[[gene_name_col]] <- NULL
  count_data$mean_count_for_dedup <- NULL
  colnames(count_data) <- sample_names
  as.data.frame(count_data, check.names = FALSE)
}

preview_count_matrix <- function(count_data, gene_name_col = "gene_name", n = 6) {
  preview <- as.data.frame(utils::head(count_data, n), check.names = FALSE)
  preview[[gene_name_col]] <- rownames(preview)
  preview <- preview[, c(gene_name_col, setdiff(colnames(preview), gene_name_col)), drop = FALSE]
  rownames(preview) <- NULL
  preview
}

make_col_data <- function(count_data, sample_names, groups, group_levels) {
  group <- factor(groups, levels = group_levels)
  data.frame(
    row.names = colnames(count_data),
    sample = colnames(count_data),
    condition = group
  )
}

filter_low_count_genes <- function(count_data, groups, min_count = 10) {
  group <- factor(groups)
  min_replicates <- min(table(group))
  keep <- rowSums(count_data >= min_count) >= min_replicates
  list(count_data = count_data[keep, , drop = FALSE], keep = keep, min_replicates = min_replicates)
}

calculate_filter_retention <- function(raw_count_data, filtered_count_data) {
  common_samples <- intersect(colnames(raw_count_data), colnames(filtered_count_data))
  raw_library <- colSums(raw_count_data[, common_samples, drop = FALSE])
  filtered_library <- colSums(filtered_count_data[, common_samples, drop = FALSE])
  data.frame(
    sample = common_samples,
    raw_library_size = raw_library,
    filtered_library_size = filtered_library,
    retained_count_fraction = filtered_library / raw_library,
    stringsAsFactors = FALSE
  )
}

write_preprocessing_summary <- function(summary_file, rows) {
  summary_df <- do.call(rbind, rows)
  dir.create(dirname(summary_file), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(summary_df, summary_file, row.names = FALSE)
  summary_df
}
