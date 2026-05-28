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

build_count_matrix <- function(rawcount, gene_name_col, count_col_names, sample_names) {
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
  index <- order(rowMeans(count_data[, expression_cols, drop = FALSE]), decreasing = TRUE)
  count_data <- count_data[index, ]
  count_data <- count_data[!duplicated(count_data[[gene_name_col]]), ]
  rownames(count_data) <- count_data[[gene_name_col]]
  count_data[[gene_name_col]] <- NULL
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
