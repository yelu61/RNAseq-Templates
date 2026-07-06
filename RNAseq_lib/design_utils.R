# Experimental design helpers for bulk RNA-seq templates.

# Build colData with condition and optional pair ID for paired designs.
make_paired_col_data <- function(sample_names, groups, group_levels, pair_id) {
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
