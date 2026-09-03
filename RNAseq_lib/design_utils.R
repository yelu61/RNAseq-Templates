# Experimental design helpers for bulk RNA-seq templates.

# Build colData with condition and optional pair ID for paired designs.
make_paired_col_data <- function(sample_names, groups, group_levels, pair_id) {
  pair_id <- align_to_sample_order(pair_id, sample_names, "PAIR_ID")
  groups <- align_to_sample_order(groups, sample_names, "GROUPS")
  if (length(pair_id) != length(sample_names)) {
    stop("pair_id length (", length(pair_id), ") must equal SAMPLE_NAMES length (", length(sample_names), ").")
  }
  if (length(groups) != length(sample_names)) {
    stop("GROUPS length must equal SAMPLE_NAMES length.")
  }

  pair_table <- data.frame(
    pair_id = factor(pair_id),
    condition = factor(groups, levels = group_levels),
    stringsAsFactors = FALSE
  )
  if (any(table(pair_table$pair_id, pair_table$condition) < 1)) {
    warning("Some pair_id groups do not contain all conditions; paired design may be unbalanced.")
  }

  data.frame(
    row.names = sample_names,
    sample = sample_names,
    condition = pair_table$condition,
    pair_id = pair_table$pair_id,
    stringsAsFactors = FALSE
  )
}

# Return a paired design formula for DESeq2 / limma.
build_paired_design_formula <- function(condition_col = "condition", pair_col = "pair_id") {
  stats::as.formula(paste0("~ ", pair_col, " + ", condition_col))
}

# Validate that pair_id covers all samples and has both conditions within each pair.
validate_paired_design <- function(sample_names, groups, pair_id, group_levels = NULL) {
  if (length(pair_id) != length(sample_names)) {
    stop("pair_id length must equal sample_names length.")
  }
  df <- data.frame(
    sample = sample_names,
    condition = groups,
    pair_id = pair_id,
    stringsAsFactors = FALSE
  )
  cond_per_pair <- stats::aggregate(condition ~ pair_id, data = df, FUN = function(x) length(unique(x)))
  if (any(cond_per_pair$condition < 2)) {
    bad_pairs <- cond_per_pair$pair_id[cond_per_pair$condition < 2]
    stop("Paired design requires both conditions in each pair. Problematic pair_id(s): ",
         paste(bad_pairs, collapse = ", "))
  }
  if (!is.null(group_levels) && !all(groups %in% group_levels)) {
    stop("Some groups not in group_levels: ", paste(setdiff(groups, group_levels), collapse = ", "))
  }
  invisible(TRUE)
}

formula_term_labels <- function(design_formula) {
  if (is.null(design_formula)) return(character(0))
  labels <- attr(stats::terms(design_formula), "term.labels")
  if (is.null(labels)) character(0) else labels
}

# Append a batch covariate column to colData so DESIGN_FORMULA can reference it.
# Pairs with validate_batch_design(): that helper warns when the formula expects a
# batch term, this one actually writes it. batch_vector is aligned to
# rownames(colData) by name when named, otherwise positionally to match
# SAMPLE_NAMES / rownames(colData).
add_batch_col <- function(colData, batch_vector, batch_term = "batch") {
  if (is.null(batch_vector)) {
    return(colData)
  }
  batch_vector <- align_to_sample_order(batch_vector, rownames(colData), "BATCH_VECTOR")
  if (nrow(colData) != length(batch_vector)) {
    stop("BATCH_VECTOR length (", length(batch_vector),
         ") does not match number of samples in colData (", nrow(colData), ").")
  }
  colData[[batch_term]] <- factor(batch_vector)
  colData
}

validate_batch_design <- function(batch_vector = NULL, design_formula = ~ condition,
                                  sample_names = NULL,
                                  batch_term = "batch") {
  if (is.null(batch_vector)) {
    return(invisible(TRUE))
  }
  if (length(batch_vector) == 0) {
    stop("BATCH_VECTOR is empty. Use NULL when there is no batch variable.")
  }
  if (any(is.na(batch_vector) | batch_vector == "")) {
    stop("BATCH_VECTOR contains NA or empty values.")
  }
  if (!is.null(sample_names) && length(batch_vector) != length(sample_names)) {
    stop("BATCH_VECTOR length (", length(batch_vector),
         ") must equal SAMPLE_NAMES length (", length(sample_names), ").")
  }
  if (length(unique(batch_vector)) < 2) {
    warning("BATCH_VECTOR has fewer than two unique batches; batch adjustment is not meaningful.")
    return(invisible(TRUE))
  }

  terms <- formula_term_labels(design_formula)
  if (!batch_term %in% terms) {
    warning(
      "BATCH_VECTOR is provided but DESIGN_FORMULA does not include '", batch_term, "'.\n",
      "Batch diagnostics will still be shown, but DEG testing will not adjust for batch.\n",
      "Use DESIGN_FORMULA <- ~ ", batch_term, " + condition if batch is a real design covariate."
    )
  }
  invisible(TRUE)
}
