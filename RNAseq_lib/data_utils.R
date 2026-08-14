# Unified data loading, validation, and gene-identifier conversion helpers for
# bulk RNA-seq templates.
#
# These functions centralize routines that were previously duplicated inline
# across template notebooks: loading expression matrices / metadata, validating
# sample overlap, converting counts to TPM/FPKM, detecting expression scale,
# and converting gene identifiers (Ensembl, MGI, HGNC).

# -----------------------------------------------------------------------------
# Expression matrix loading
# -----------------------------------------------------------------------------

#' Read an expression matrix from CSV/TSV/Excel.
#'
#' The first column is used as row names if it is non-numeric and `gene_column`
#' is NULL. If `gene_column` is provided, that column becomes row names. The
#' remaining columns are coerced to numeric expression values.
#'
#' @param file Path to the expression file.
#' @param gene_column Name of the gene identifier column, or NULL to use the
#'   first non-numeric column.
#' @param sample_ids Optional vector of sample IDs to subset/reorder columns.
#'   Missing samples trigger a warning.
#' @param file_format "auto" (infer from extension), "csv", "tsv", or "excel".
#' @return A numeric matrix with gene identifiers as row names.
read_expression_matrix <- function(file,
                                   gene_column = NULL,
                                   sample_ids = NULL,
                                   file_format = c("auto", "csv", "tsv", "excel")) {
  file_format <- match.arg(file_format)
  if (!file.exists(file)) {
    stop("Expression file not found: ", file)
  }

  if (file_format == "auto") {
    ext <- tolower(tools::file_ext(file))
    file_format <- dplyr::case_when(
      ext == "csv" ~ "csv",
      ext %in% c("tsv", "txt") ~ "tsv",
      ext %in% c("xlsx", "xls") ~ "excel",
      TRUE ~ "tsv"
    )
  }

  if (file_format == "excel") {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for Excel input.")
    }
    df <- openxlsx::read.xlsx(file)
  } else if (file_format == "csv") {
    df <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    df <- utils::read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  }

  df <- as.data.frame(df, check.names = FALSE)
  if (nrow(df) == 0) {
    stop("Expression file is empty: ", file)
  }

  # Determine gene identifier column
  if (!is.null(gene_column)) {
    if (!gene_column %in% colnames(df)) {
      stop("Gene column not found in expression file: ", gene_column)
    }
    gene_ids <- as.character(df[[gene_column]])
    expr_df <- df[, setdiff(colnames(df), gene_column), drop = FALSE]
  } else if (!is.numeric(df[[1]]) && !all(grepl("^[0-9.]+$", as.character(df[[1]])))) {
    gene_ids <- as.character(df[[1]])
    expr_df <- df[, -1, drop = FALSE]
  } else {
    stop("The first column is numeric, so the gene identifier column cannot be inferred safely. ",
         "Set gene_column explicitly (numeric Entrez IDs are supported when named explicitly).")
  }

  # Coerce expression values to numeric
  expr_df[] <- lapply(expr_df, function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  })

  if (any(vapply(expr_df, function(x) anyNA(x), logical(1)))) {
    stop("Expression matrix contains non-numeric values outside the gene column.")
  }

  # Filter empty/NA gene IDs, then deduplicate by mean expression before
  # assigning row names (avoids duplicate-row-name assignment error).
  keep <- !is.na(gene_ids) & gene_ids != ""
  gene_ids <- gene_ids[keep]
  expr_df <- expr_df[keep, , drop = FALSE]

  if (any(duplicated(gene_ids))) {
    mean_expr <- rowMeans(expr_df, na.rm = TRUE)
    ord <- order(mean_expr, decreasing = TRUE)
    expr_df <- expr_df[ord, , drop = FALSE]
    gene_ids <- gene_ids[ord]
    keep_dup <- !duplicated(gene_ids)
    expr_df <- expr_df[keep_dup, , drop = FALSE]
    gene_ids <- gene_ids[keep_dup]
  }

  rownames(expr_df) <- gene_ids

  if (anyDuplicated(colnames(expr_df))) {
    stop("Expression matrix contains duplicated sample column names: ",
         paste(unique(colnames(expr_df)[duplicated(colnames(expr_df))]), collapse = ", "))
  }

  # Subset/reorder samples if requested
  if (!is.null(sample_ids)) {
    missing_samples <- setdiff(sample_ids, colnames(expr_df))
    if (length(missing_samples) > 0) {
      warning("Requested samples missing from expression matrix: ",
              paste(missing_samples, collapse = ", "))
    }
    available_samples <- intersect(sample_ids, colnames(expr_df))
    if (length(available_samples) == 0) {
      stop("None of the requested samples were found in the expression matrix.")
    }
    expr_df <- expr_df[, available_samples, drop = FALSE]
  }

  as.matrix(expr_df)
}

# -----------------------------------------------------------------------------
# Metadata loading
# -----------------------------------------------------------------------------

#' Read sample metadata with validation.
#'
#' Handles duplicated sample columns (renames extras), validates required
# columns, and optionally factorizes group/time columns.
#'
#' @param file Path to metadata file.
#' @param sample_column Name of the sample ID column.
#' @param required_columns Character vector of required column names.
#' @param group_column Optional group column to factorize.
#' @param group_levels Optional levels for the group factor.
#' @param time_column Optional time column to factorize.
#' @param time_levels Optional levels for the time factor.
#' @return A data frame.
read_metadata <- function(file,
                          sample_column = "sample",
                          required_columns = NULL,
                          group_column = NULL,
                          group_levels = NULL,
                          time_column = NULL,
                          time_levels = NULL) {
  if (!file.exists(file)) {
    stop("Metadata file not found: ", file)
  }

  ext <- tolower(tools::file_ext(file))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("Package 'openxlsx' is required for Excel metadata.")
    }
    meta <- openxlsx::read.xlsx(file)
  } else if (ext == "csv") {
    meta <- utils::read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    meta <- utils::read.delim(file, check.names = FALSE, stringsAsFactors = FALSE)
  }

  meta <- as.data.frame(meta, check.names = FALSE)
  if (nrow(meta) == 0) {
    stop("Metadata file is empty: ", file)
  }

  # Handle duplicated sample columns (e.g. colData.csv exported by General notebook)
  if (sum(colnames(meta) == sample_column) > 1) {
    dup_idx <- which(colnames(meta) == sample_column)
    colnames(meta)[dup_idx[-1]] <- paste0(sample_column, "_", seq_len(length(dup_idx) - 1))
  }

  if (!sample_column %in% colnames(meta)) {
    stop("Sample column not found in metadata: ", sample_column)
  }

  meta[[sample_column]] <- as.character(meta[[sample_column]])
  if (any(is.na(meta[[sample_column]]) | meta[[sample_column]] == "")) {
    stop("Metadata sample column contains missing or empty values.")
  }
  if (any(duplicated(meta[[sample_column]]))) {
    stop("Metadata sample column contains duplicate sample IDs.")
  }

  if (!is.null(required_columns)) {
    missing_cols <- setdiff(required_columns, colnames(meta))
    if (length(missing_cols) > 0) {
      stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
    }
  }

  if (!is.null(group_column) && group_column %in% colnames(meta)) {
    if (is.null(group_levels)) group_levels <- unique(as.character(meta[[group_column]]))
    unknown_groups <- setdiff(unique(as.character(meta[[group_column]])), group_levels)
    if (length(unknown_groups) > 0) {
      stop("Metadata contains group values not listed in group_levels: ", paste(unknown_groups, collapse = ", "))
    }
    meta[[group_column]] <- factor(as.character(meta[[group_column]]), levels = group_levels)
  }

  if (!is.null(time_column) && time_column %in% colnames(meta)) {
    if (is.null(time_levels)) {
      time_levels <- unique(as.character(meta[[time_column]]))
    } else {
      unknown_times <- setdiff(unique(as.character(meta[[time_column]])), time_levels)
      if (length(unknown_times) > 0) {
        stop("Metadata contains time values not listed in time_levels: ", paste(unknown_times, collapse = ", "))
      }
      missing_levels <- setdiff(time_levels, unique(as.character(meta[[time_column]])))
      if (length(missing_levels) > 0) {
        warning("Time levels specified but not observed in data: ",
                paste(missing_levels, collapse = ", "))
      }
    }
    meta[[time_column]] <- factor(as.character(meta[[time_column]]), levels = time_levels)
  }

  meta
}

#' Encode metadata traits for WGCNA module--trait correlations.
#'
#' Sample identifiers (including duplicated mirror columns such as `sample_1`)
#' are excluded. Numeric traits are retained as-is; each non-reference level of
#' a categorical trait becomes a clearly labelled 0/1 indicator. All-missing
#' and invariant columns are dropped because their correlations are undefined.
#'
#' @param traits Sample metadata data frame.
#' @param sample_column Name of the primary sample identifier column.
#' @return A numeric matrix with samples in rows and encoded traits in columns.
encode_wgcna_traits <- function(traits, sample_column = "sample") {
  traits <- as.data.frame(traits, check.names = FALSE)
  if (!sample_column %in% colnames(traits)) {
    stop("Sample column not found in metadata: ", sample_column)
  }

  sample_ids <- as.character(traits[[sample_column]])
  candidate_names <- setdiff(colnames(traits), sample_column)
  # read_metadata() deliberately preserves duplicate columns by renaming them.
  # Do not accidentally correlate module eigengenes with a renamed copy of the
  # sample identifier, a failure mode observed in a real template-derived run.
  mirrored_ids <- candidate_names[vapply(candidate_names, function(nm) {
    identical(as.character(traits[[nm]]), sample_ids)
  }, logical(1))]
  candidate_names <- setdiff(candidate_names, mirrored_ids)

  blocks <- list()
  for (nm in candidate_names) {
    x <- traits[[nm]]
    if (is.numeric(x) || is.integer(x) || is.logical(x)) {
      values <- as.numeric(x)
      finite_values <- values[is.finite(values)]
      if (length(finite_values) > 1 && stats::sd(finite_values) > 0) {
        block <- matrix(values, ncol = 1, dimnames = list(NULL, nm))
        blocks[[length(blocks) + 1L]] <- block
      }
      next
    }

    observed <- as.character(x)
    levels_x <- if (is.factor(x)) levels(droplevels(x)) else unique(observed[!is.na(observed)])
    if (length(levels_x) < 2) next
    reference <- levels_x[[1]]
    for (level in levels_x[-1]) {
      values <- ifelse(is.na(observed), NA_real_, as.numeric(observed == level))
      label <- paste0(nm, ": ", level)
      block <- matrix(values, ncol = 1, dimnames = list(NULL, label))
      attr(block, "reference_level") <- reference
      blocks[[length(blocks) + 1L]] <- block
    }
  }

  if (!length(blocks)) {
    return(matrix(numeric(), nrow = nrow(traits), ncol = 0,
                  dimnames = list(rownames(traits), character())))
  }
  out <- do.call(cbind, blocks)
  rownames(out) <- rownames(traits)
  storage.mode(out) <- "double"
  out
}

# -----------------------------------------------------------------------------
# Sample name validation
# -----------------------------------------------------------------------------

#' Validate that expression and metadata sample names match.
#'
#' @param expr_samples Sample names from the expression matrix.
#' @param meta_samples Sample names from metadata.
#' @param context Description used in error/warning messages.
#' @param strict_order If TRUE, require identical order as well as identical set.
#' @return Invisible TRUE.
validate_samples_match <- function(expr_samples,
                                   meta_samples,
                                   context = "expression vs metadata",
                                   strict_order = FALSE) {
  if (anyNA(expr_samples) || any(expr_samples == "") || anyDuplicated(expr_samples)) {
    stop("Expression sample names must be non-empty and unique.")
  }
  if (anyNA(meta_samples) || any(meta_samples == "") || anyDuplicated(meta_samples)) {
    stop("Metadata sample names must be non-empty and unique.")
  }
  common <- intersect(expr_samples, meta_samples)
  if (length(common) == 0) {
    stop("No common samples between ", context, ".\n",
         "Expression samples (first 5): ", paste(head(expr_samples, 5), collapse = ", "), "\n",
         "Metadata samples (first 5): ", paste(head(meta_samples, 5), collapse = ", "))
  }

  only_expr <- setdiff(expr_samples, meta_samples)
  only_meta <- setdiff(meta_samples, expr_samples)
  if (length(only_expr) > 0) {
    warning("Samples in expression but not metadata (", context, "): ",
            paste(head(only_expr, 10), collapse = ", "))
  }
  if (length(only_meta) > 0) {
    warning("Samples in metadata but not expression (", context, "): ",
            paste(head(only_meta, 10), collapse = ", "))
  }

  if (strict_order && !identical(expr_samples, meta_samples)) {
    stop("Sample order mismatch between ", context, ".")
  }

  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Expression scale / unit detection
# -----------------------------------------------------------------------------

#' Heuristically classify the scale of an expression matrix.
#'
#' Returns a list with scale ("raw_counts", "log2_tpm", "tpm", "vst", "unknown"),
#' confidence ("high"/"medium"/"low"), and explanatory notes.
#'
#' @param mat Numeric expression matrix (genes x samples).
#' @param sample_n Number of columns to subsample for speed.
#' @param gene_n Number of rows to subsample for speed.
detect_expression_scale <- function(mat,
                                    sample_n = min(ncol(mat), 500),
                                    gene_n = min(nrow(mat), 2000)) {
  if (!is.matrix(mat) && !is.data.frame(mat)) {
    stop("Input must be a matrix or data.frame.")
  }
  mat <- as.matrix(mat)
  if (!is.numeric(mat)) {
    stop("Expression matrix must be numeric.")
  }

  if (ncol(mat) > sample_n) {
    mat <- mat[, sort(sample.int(ncol(mat), sample_n)), drop = FALSE]
  }
  if (nrow(mat) > gene_n) {
    mat <- mat[sort(sample.int(nrow(mat), gene_n)), , drop = FALSE]
  }

  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return(list(scale = "unknown", confidence = "high", notes = "No finite values found"))
  }

  any_negative <- any(vals < 0, na.rm = TRUE)
  max_val <- max(vals, na.rm = TRUE)
  min_val <- min(vals, na.rm = TRUE)
  mean_val <- mean(vals, na.rm = TRUE)
  median_val <- stats::median(vals, na.rm = TRUE)

  # Fraction of integer values among values >= 1
  countable_vals <- vals[vals >= 1]
  integer_fraction <- if (length(countable_vals) > 0) {
    mean(countable_vals %% 1 < 1e-9, na.rm = TRUE)
  } else {
    0
  }

  # Row sums for TPM-like detection
  mat_non_neg <- pmax(mat, 0)
  col_sums <- colSums(mat_non_neg, na.rm = TRUE)
  median_col_sum <- stats::median(col_sums, na.rm = TRUE)
  col_sum_near_1e6 <- median_col_sum > 5e5 && median_col_sum < 5e6

  # Zero fraction (raw counts tend to have many zeros)
  zero_fraction <- mean(vals == 0, na.rm = TRUE)

  notes <- character(0)
  scale <- "unknown"
  confidence <- "low"

  if (any_negative) {
    if (max_val < 100 && min_val < 0) {
      scale <- "vst"
      confidence <- "medium"
      notes <- c(notes, "Negative values and bounded range suggest VST/rlog or similar transformed data")
    } else {
      scale <- "unknown"
      confidence <- "low"
      notes <- c(notes, "Negative values present; scale cannot be reliably determined")
    }
  } else if (max_val > 500 && integer_fraction > 0.75 && zero_fraction > 0.1) {
    scale <- "raw_counts"
    confidence <- "high"
    notes <- c(notes, "Large integer-dominant values with many zeros suggest raw counts")
  } else if (max_val <= 50 && integer_fraction < 0.25) {
    if (col_sum_near_1e6 && mean_val > 1) {
      scale <- "tpm"
      confidence <- "medium"
      notes <- c(notes, "Values bounded and column sums near 1e6 suggest TPM/FPKM-like data")
    } else {
      scale <- "log2_tpm"
      confidence <- "medium"
      notes <- c(notes, "Small non-integer values suggest log2(TPM+1) or similar log-scale data")
    }
  } else if (max_val > 50 && max_val <= 500 && integer_fraction < 0.5 && zero_fraction < 0.3) {
    scale <- "log2_tpm"
    confidence <- "low"
    notes <- c(notes, "Mixed-range non-integer values; possibly log-scaled data")
  } else if (col_sum_near_1e6 && integer_fraction < 0.5) {
    scale <- "tpm"
    confidence <- "medium"
    notes <- c(notes, "Column sums near 1e6 suggest TPM/FPKM-like data")
  } else if (max_val > 500 && integer_fraction > 0.75) {
    scale <- "raw_counts"
    confidence <- "medium"
    notes <- c(notes, "Large integer-dominant values suggest raw counts")
  } else {
    notes <- c(notes, "Could not confidently determine expression scale")
  }

  list(
    scale = scale,
    confidence = confidence,
    notes = paste(notes, collapse = "; ")
  )
}

# -----------------------------------------------------------------------------
# Count normalization helpers
# -----------------------------------------------------------------------------

#' Validate and align gene lengths to a count matrix.
#'
#' @param counts_mat Numeric count matrix.
#' @param gene_lengths_kb Gene lengths in kilobases. Either a named vector or a
#'   vector whose length equals nrow(counts_mat).
.validate_gene_lengths <- function(counts_mat, gene_lengths_kb) {
  if (!is.numeric(gene_lengths_kb) || any(is.na(gene_lengths_kb)) || any(gene_lengths_kb <= 0)) {
    stop("gene_lengths_kb must be positive numeric values with no missing values.")
  }

  if (!is.null(names(gene_lengths_kb))) {
    aligned <- gene_lengths_kb[match(rownames(counts_mat), names(gene_lengths_kb))]
    if (any(is.na(aligned))) {
      stop("Some row names of the count matrix have no matching gene length entry: ",
           paste(head(rownames(counts_mat)[is.na(aligned)], 10), collapse = ", "))
    }
    return(aligned)
  }

  if (length(gene_lengths_kb) != nrow(counts_mat)) {
    stop("gene_lengths_kb length (", length(gene_lengths_kb),
         ") must equal number of rows in counts_mat (", nrow(counts_mat),
         "), or gene_lengths_kb must be a named vector.")
  }
  gene_lengths_kb
}

#' Convert raw counts to TPM.
#'
#' @param counts_mat Numeric count matrix (genes x samples).
#' @param gene_lengths_kb Gene lengths in kilobases, aligned to row names or as
#'   a vector of equal length.
#' @return Numeric TPM matrix.
counts_to_tpm <- function(counts_mat, gene_lengths_kb) {
  counts_mat <- as.matrix(counts_mat)
  if (!is.numeric(counts_mat) || any(is.na(counts_mat)) || any(counts_mat < 0)) {
    stop("counts_mat must be a numeric matrix with non-negative values.")
  }

  lengths <- .validate_gene_lengths(counts_mat, gene_lengths_kb)
  rpk <- counts_mat / lengths
  rpk_sums <- colSums(rpk, na.rm = TRUE)
  if (any(!is.finite(rpk_sums) | rpk_sums <= 0)) {
    stop("Cannot compute TPM: one or more samples have zero total RPK.")
  }
  tpm <- t(t(rpk) / rpk_sums) * 1e6
  tpm
}

#' Collapse a count/expression matrix to unique gene symbols.
#'
#' Keeps, for each duplicated symbol, the row with the highest total signal
#' (row sum). Rows with missing/empty symbols are dropped. Optionally carries a
#' per-row annotation (e.g. gene length, biotype) through to the retained row so
#' downstream steps (TPM, TME) stay aligned to the kept feature.
#'
#' @param symbols Character vector of gene symbols, one per row of `mat`.
#' @param mat Numeric matrix/data.frame (features x samples) to collapse.
#' @param extra Optional vector (or single-column frame) aligned to the rows of
#'   `mat`; the value from the retained row is returned. `NULL` to skip.
#' @return A list with `gene_name` (unique symbols), `counts` (collapsed matrix),
#'   `n_in`/`n_out` (input/output feature counts), and `extra` when supplied.
collapse_by_symbol <- function(symbols, mat, extra = NULL) {
  mat <- as.matrix(mat)
  keep <- !is.na(symbols) & symbols != ""
  symbols <- symbols[keep]; mat <- mat[keep, , drop = FALSE]
  if (!is.null(extra)) extra <- extra[keep]
  ord <- order(rowSums(mat, na.rm = TRUE), decreasing = TRUE)
  symbols <- symbols[ord]; mat <- mat[ord, , drop = FALSE]
  if (!is.null(extra)) extra <- extra[ord]
  dup <- duplicated(symbols)
  out <- list(gene_name = symbols[!dup], counts = mat[!dup, , drop = FALSE],
              n_in = length(keep), n_out = sum(!dup))
  if (!is.null(extra)) out$extra <- extra[!dup]
  out
}

#' Convert raw counts to FPKM.
#'
#' @param counts_mat Numeric count matrix (genes x samples).
#' @param gene_lengths_kb Gene lengths in kilobases.
#' @return Numeric FPKM matrix.
counts_to_fpkm <- function(counts_mat, gene_lengths_kb) {
  counts_mat <- as.matrix(counts_mat)
  if (!is.numeric(counts_mat) || any(is.na(counts_mat)) || any(counts_mat < 0)) {
    stop("counts_mat must be a numeric matrix with non-negative values.")
  }

  lengths <- .validate_gene_lengths(counts_mat, gene_lengths_kb)
  library_size <- colSums(counts_mat, na.rm = TRUE)
  if (any(library_size == 0)) {
    stop("Cannot compute FPKM: one or more samples have total count of zero.")
  }
  fpkm <- counts_mat / lengths / (library_size / 1e6)
  fpkm
}

#' Extract gene lengths (in kb) from an annotation data frame.
#'
#' @param raw_annot Data frame with annotation columns.
#' @param id_col Column containing gene identifiers that match count row names.
#' @param length_col Optional column with gene lengths (in bp or kb; see
#'   length_unit).
#' @param start_col Column with gene start coordinates.
#' @param end_col Column with gene end coordinates.
#' @param length_unit Unit of length_col if provided: "bp" or "kb".
#' @return Named numeric vector of lengths in kb.
extract_gene_lengths <- function(raw_annot,
                                 id_col = "gene_id",
                                 length_col = NULL,
                                 start_col = "gene_start",
                                 end_col = "gene_end",
                                 length_unit = c("bp", "kb")) {
  raw_annot <- as.data.frame(raw_annot, check.names = FALSE)
  length_unit <- match.arg(length_unit)

  if (!id_col %in% colnames(raw_annot)) {
    stop("ID column not found in annotation data: ", id_col)
  }
  ids <- as.character(raw_annot[[id_col]])

  if (!is.null(length_col) && length_col %in% colnames(raw_annot)) {
    lengths <- as.numeric(as.character(raw_annot[[length_col]]))
    if (length_unit == "bp") lengths <- lengths / 1000
  } else if (all(c(start_col, end_col) %in% colnames(raw_annot))) {
    starts <- as.numeric(as.character(raw_annot[[start_col]]))
    ends <- as.numeric(as.character(raw_annot[[end_col]]))
    lengths <- (ends - starts + 1) / 1000
  } else {
    stop("Gene length information not found. Provide length_col or start_col/end_col.")
  }

  if (any(is.na(lengths)) || any(lengths <= 0)) {
    stop("Some gene lengths are missing or non-positive. Check length/start/end columns.")
  }

  names(lengths) <- ids
  lengths
}

# -----------------------------------------------------------------------------
# Count matrix validation
# -----------------------------------------------------------------------------

#' Validate a count matrix.
#'
#' @param mat Numeric matrix or data.frame.
#' @param require_integer Require integer values.
#' @param require_non_negative Require non-negative values.
#' @param min_samples Minimum number of sample columns.
validate_count_matrix <- function(mat,
                                  require_integer = TRUE,
                                  require_non_negative = TRUE,
                                  min_samples = 2) {
  mat <- as.matrix(mat)
  if (!is.numeric(mat)) {
    stop("Count matrix must be numeric.")
  }
  if (anyNA(mat)) {
    stop("Count matrix contains missing values.")
  }
  if (require_non_negative && any(mat < 0)) {
    stop("Count matrix contains negative values. Counts must be >= 0.")
  }
  if (require_integer && any(mat %% 1 != 0)) {
    stop("Count matrix contains non-integer values. DESeq2/limma-voom require raw integer counts.")
  }
  if (ncol(mat) < min_samples) {
    stop("Count matrix must have at least ", min_samples, " sample columns.")
  }
  if (is.null(rownames(mat)) || any(rownames(mat) == "")) {
    stop("Count matrix must have non-empty row names (gene identifiers).")
  }
  if (anyDuplicated(rownames(mat))) {
    stop("Count matrix gene identifiers must be unique. Use stable gene IDs for TPM conversion or collapse duplicates explicitly.")
  }
  invisible(TRUE)
}

#' Validate a normalized expression matrix against an explicit input contract.
#'
#' Heuristic scale detection is useful for diagnostics, but must not silently
#' choose a transformation for an analysis method.
validate_expression_contract <- function(mat,
                                         expected = c("normalized", "tpm", "log2_tpm", "vst"),
                                         tolerance = 0.05) {
  expected <- match.arg(expected)
  mat <- as.matrix(mat)
  if (!is.numeric(mat) || anyNA(mat) || any(!is.finite(mat))) {
    stop("Expression matrix must contain finite numeric values with no missing values.")
  }
  if (is.null(rownames(mat)) || anyNA(rownames(mat)) || any(rownames(mat) == "") || anyDuplicated(rownames(mat))) {
    stop("Expression matrix must have non-empty, unique gene identifiers as row names.")
  }
  if (is.null(colnames(mat)) || anyNA(colnames(mat)) || any(colnames(mat) == "") || anyDuplicated(colnames(mat))) {
    stop("Expression matrix must have non-empty, unique sample names as column names.")
  }

  scale_info <- detect_expression_scale(mat)
  if (expected == "tpm") {
    if (any(mat < 0)) stop("TPM input cannot contain negative values.")
    col_sums <- colSums(mat)
    if (any(col_sums <= 0)) stop("TPM input contains a sample with zero total expression.")
    relative_error <- abs(col_sums - 1e6) / 1e6
    if (any(relative_error > tolerance)) {
      warning("TPM column sums differ from 1e6 by more than ", tolerance * 100,
              "%. Confirm that the input is TPM rather than FPKM or another normalized unit.")
    }
  } else if (expected == "log2_tpm") {
    if (any(mat < 0)) stop("log2(TPM+1) input cannot contain negative values.")
    if (identical(scale_info$scale, "raw_counts")) stop("Expected log2(TPM+1) but input looks like raw counts.")
  } else if (expected == "vst") {
    if (identical(scale_info$scale, "raw_counts")) stop("Expected VST/rlog expression but input looks like raw counts.")
  } else if (identical(scale_info$scale, "raw_counts")) {
    stop("Expected normalized expression but input looks like raw counts.")
  }
  invisible(scale_info)
}

# -----------------------------------------------------------------------------
# Gene identifier conversion
# -----------------------------------------------------------------------------

#' Detect whether gene identifiers are Ensembl IDs or symbols.
#'
#' @param ids Character vector of gene identifiers.
#' @return "ensembl" or "symbol".
detect_gene_id_type <- function(ids) {
  ids <- as.character(ids)
  ids <- ids[!is.na(ids) & ids != ""]
  if (length(ids) == 0) return("symbol")

  ensembl_like <- grepl("^(ENS|ENSMUS|ENSMUST|ENSG)", ids, ignore.case = TRUE)
  if (mean(ensembl_like) > 0.5) {
    "ensembl"
  } else {
    "symbol"
  }
}

#' Get the species-specific org.*.eg.db package name.
#'
#' @param species "human" or "mouse".
.get_org_db_name <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  switch(species,
         human = "org.Hs.eg.db",
         mouse = "org.Mm.eg.db")
}

#' Convert gene identifiers using AnnotationDbi.
#'
#' @param ids Character vector of IDs.
#' @param from Keytype ("ENSEMBL" or "SYMBOL").
#' @param to Keytype ("SYMBOL" or "ENSEMBL").
#' @param species "human" or "mouse".
#' @return Data frame with columns id, converted, unmapped.
.convert_gene_ids_annotationdbi <- function(ids, from, to, species) {
  pkg <- .get_org_db_name(species)
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for offline gene ID conversion.")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))

  db <- get(pkg)
  lookup_ids <- if (identical(from, "ENSEMBL")) sub("[.][0-9]+$", "", ids) else ids
  mapped <- tryCatch(
    AnnotationDbi::mapIds(
      db,
      keys = lookup_ids,
      column = to,
      keytype = from,
      multiVals = "first"
    ),
    error = function(e) {
      if (grepl("valid keys", conditionMessage(e), ignore.case = TRUE)) {
        stats::setNames(rep(NA_character_, length(ids)), ids)
      } else {
        stop(e)
      }
    }
  )

  data.frame(
    id = ids,
    converted = as.character(mapped),
    unmapped = is.na(mapped),
    stringsAsFactors = FALSE
  )
}

#' Convert MGI symbols to HGNC symbols using babelgene.
#'
#' @param ids Character vector of MGI symbols.
#' @return Data frame with columns id, converted, unmapped.
.convert_mouse_to_human_babelgene <- function(ids) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    stop("Package 'babelgene' is required for offline mouse-to-human conversion.")
  }

  mapping <- babelgene::orthologs(
    genes = ids,
    species = "mouse",
    human = FALSE,
    top = TRUE
  )

  if (is.null(mapping) || nrow(mapping) == 0) {
    return(data.frame(id = ids, converted = NA_character_, unmapped = TRUE, stringsAsFactors = FALSE))
  }

  mapping <- data.frame(
    id = toupper(as.character(mapping$symbol)),
    converted = as.character(mapping$human_symbol),
    stringsAsFactors = FALSE
  )
  mapping <- mapping[mapping$converted != "" & !is.na(mapping$converted), ]

  res <- data.frame(
    id = toupper(ids),
    converted = NA_character_,
    unmapped = TRUE,
    stringsAsFactors = FALSE
  )
  idx <- match(mapping$id, res$id)
  res$converted[idx] <- mapping$converted
  res$unmapped <- is.na(res$converted)
  res
}

#' Convert gene identifiers between Ensembl, symbol, and human orthologs.
#'
#' @param ids Character vector of gene identifiers.
#' @param from "auto", "ensembl", or "symbol".
#' @param to "symbol", "ensembl", or "human_symbol" (mouse MGI → HGNC).
#' @param species "human" or "mouse".
#' @param method "auto", "annotationdbi", or "babelgene".
#' @return Data frame with columns id, converted, unmapped.
convert_gene_ids <- function(ids,
                             from = c("auto", "ensembl", "symbol"),
                             to = c("symbol", "ensembl", "human_symbol"),
                             species = c("human", "mouse"),
                             method = c("auto", "annotationdbi", "babelgene")) {
  from <- match.arg(from)
  to <- match.arg(to)
  species <- match.arg(species)
  method <- match.arg(method)

  if (from == "auto") {
    from <- detect_gene_id_type(ids)
  }

  if (to == "human_symbol") {
    if (species != "mouse") {
      stop("to='human_symbol' requires species='mouse'.")
    }
    if (from != "symbol") {
      stop("to='human_symbol' currently requires from='symbol' (MGI symbols).")
    }
    if (method %in% c("auto", "babelgene")) {
      return(.convert_mouse_to_human_babelgene(ids))
    }
    if (method == "annotationdbi") {
      stop("AnnotationDbi cannot directly convert MGI to HGNC; use method='babelgene'.")
    }
  }

  keytype_from <- if (from == "ensembl") "ENSEMBL" else "SYMBOL"
  keytype_to <- if (to == "ensembl") "ENSEMBL" else "SYMBOL"

  if (method %in% c("auto", "annotationdbi")) {
    pkg <- .get_org_db_name(species)
    if (requireNamespace(pkg, quietly = TRUE)) {
      return(.convert_gene_ids_annotationdbi(ids, keytype_from, keytype_to, species))
    }
    if (method == "annotationdbi") {
      stop("Package '", pkg, "' is not installed.")
    }
  }

  stop("Could not convert gene IDs with method='", method, "'. Install the appropriate annotation package.")
}

#' Convert row names of an expression matrix and deduplicate by mean expression.
#'
#' @param expr Numeric expression matrix.
#' @param species "human" or "mouse".
#' @param target "symbol" or "human_symbol".
#' @param method Conversion backend.
#' @return Expression matrix with converted row names.
convert_expression_rownames <- function(expr,
                                        species = c("human", "mouse"),
                                        target = c("symbol", "human_symbol"),
                                        method = c("auto", "annotationdbi", "babelgene")) {
  species <- match.arg(species)
  target <- match.arg(target)
  method <- match.arg(method)

  ids <- rownames(expr)
  if (is.null(ids) || any(ids == "")) {
    stop("Expression matrix must have non-empty row names for gene ID conversion.")
  }

  mapping <- convert_gene_ids(ids, from = "auto", to = target, species = species, method = method)
  mapped <- mapping$converted
  names(mapped) <- mapping$id

  if (all(is.na(mapped))) {
    stop("No gene IDs could be converted. Check that row names are valid gene identifiers.")
  }

  if (any(is.na(mapped))) {
    warning("Unmapped gene identifiers dropped: ", sum(is.na(mapped)))
  }

  expr_df <- as.data.frame(expr, check.names = FALSE)
  expr_df$.__converted_symbol <- mapped
  expr_df <- expr_df[!is.na(expr_df$.__converted_symbol), , drop = FALSE]

  # Deduplicate on converted symbol by highest mean expression
  mean_expr <- rowMeans(expr_df[, setdiff(colnames(expr_df), ".__converted_symbol"), drop = FALSE], na.rm = TRUE)
  ord <- order(mean_expr, decreasing = TRUE)
  expr_df <- expr_df[ord, , drop = FALSE]
  expr_df <- expr_df[!duplicated(expr_df$.__converted_symbol), , drop = FALSE]
  rownames(expr_df) <- expr_df$.__converted_symbol
  expr_df$.__converted_symbol <- NULL

  as.matrix(expr_df)
}
