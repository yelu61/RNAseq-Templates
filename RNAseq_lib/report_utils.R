# Report generation helpers for RNAseq-Templates.
# Renders a self-contained HTML report from the CSV/PDF outputs that the
# notebooks already produce, without re-running any analysis.

# Render the shared Quarto/R Markdown analysis report.
#
# @param outdir Project output directory (contains 1-DEG/, 2-GSEA/, 3-Visualization/).
# @param report_file Path of the HTML report to write.
# @param template Path to the .qmd/.Rmd template (default: reports/analysis_report.qmd
#   resolved relative to the repository root discovered via rprojroot, else CWD).
# @param params Named list of parameters passed to the report (title, author, etc.).
# @return Invisibly, the path to the rendered report.
render_analysis_report <- function(outdir = ".",
                                   report_file = file.path(outdir, "RNAseq_report.html"),
                                   template = NULL,
                                   params = list()) {
  if (is.null(template)) {
    repo_root <- tryCatch(rprojroot::find_root(rprojroot::is_git_root), error = function(e) getwd())
    candidate <- file.path(repo_root, "reports", "analysis_report.qmd")
    template <- if (file.exists(candidate)) candidate else file.path("reports", "analysis_report.qmd")
  }
  if (!file.exists(template)) {
    stop("Report template not found: ", template)
  }

  params <- utils::modifyList(
    list(
      title = "RNA-seq Analysis Report",
      author = "",
      outdir = normalizePath(outdir, mustWork = FALSE)
    ),
    params
  )

  # A report is complete only when every analysis domain is accounted for.
  # Generate the machine-readable coverage manifest before rendering so the
  # HTML can display uncovered optional domains instead of silently omitting
  # them. Project-specific workflows may replace rows with explicit
  # not_applicable / omitted_with_reason decisions before calling this helper.
  coverage_file <- file.path(outdir, "report_coverage_manifest.csv")
  if (!file.exists(coverage_file)) {
    write_report_coverage_manifest(outdir, coverage_file)
  }

  template_ext <- tolower(tools::file_ext(template))
  quarto_available <- requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))
  if (quarto_available) {
    quarto::quarto_render(
      input = template,
      output_file = basename(report_file),
      execute_params = params,
      quiet = TRUE
    )
    # quarto renders next to the template; move to the requested location.
    produced <- file.path(dirname(template), basename(report_file))
    if (normalizePath(produced, mustWork = FALSE) != normalizePath(report_file, mustWork = FALSE)) {
      file.copy(produced, report_file, overwrite = TRUE)
      unlink(produced)
    }
  } else if (template_ext %in% c("rmd", "md") && requireNamespace("rmarkdown", quietly = TRUE)) {
    rmarkdown::render(
      input = template,
      output_file = basename(report_file),
      output_dir = dirname(report_file),
      params = params,
      envir = new.env(parent = globalenv()),
      quiet = TRUE
    )
  } else if (template_ext == "qmd") {
    stop("The selected .qmd report requires the Quarto CLI and R 'quarto' package. ",
         "Install Quarto or supply an .Rmd template explicitly.")
  } else {
    stop("HTML report requires the Quarto CLI, or rmarkdown for an .Rmd template.")
  }

  if (!file.exists(report_file)) {
    stop("Report rendering did not produce: ", report_file)
  }
  message("HTML report written to: ", report_file)
  invisible(report_file)
}

# Default analysis-domain contract for a comprehensive bulk RNA-seq report.
# "Required" means the generic General workflow should normally produce the
# domain; optional domains must still appear in the coverage manifest with an
# explicit availability/status instead of disappearing from the report.
default_report_coverage_rules <- function() {
  data.frame(
    domain = c(
      "study_design", "input_filtering", "sample_qc", "differential_expression",
      "threshold_sensitivity", "ora", "gsea", "custom_gene_sets",
      "cross_contrast", "biological_synthesis", "reproducibility", "validation"
    ),
    label = c(
      "Study design and sample inclusion", "Input filtering and retained features",
      "Sample quality control", "Differential expression",
      "Threshold sensitivity", "Over-representation analysis",
      "Rank-based enrichment (GSEA)", "Custom gene-set scoring",
      "Cross-contrast comparison", "Claim-evidence synthesis",
      "Reproducibility", "Independent output validation"
    ),
    required = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
                 FALSE, FALSE, TRUE, FALSE),
    pattern = c(
      "(colData[.]csv|sample.*metadata|0-Config|Analysis_summary[.]txt)",
      "(Filter_retention|Raw_count_distribution|Duplicated_gene_symbols|filter.*summary)",
      "(PCA|Sample_distance|sample.*correlation|Sample_QC_metrics)",
      "(DEG_threshold_summary|DESeq2_all_genes|Volcano_|Heatmap_topDEG)",
      "(DEG_threshold_summary|DEG_lfc_strategy|/(strict|standard|loose)/)",
      "(ORA|CompareCluster)",
      "(GSEA|gsea_)",
      "(GSVA|CustomGeneSets|pathway_gsva|gene_set)",
      "(CompareCluster|overlap|Jaccard|treatment_effect)",
      "(claim_evidence_ledger|ReasoningBrief|report_narrative)",
      "(sessionInfo|analysis_manifest|config[.]R)",
      "(report_validation|figure_manifest|validation_checks)"
    ),
    stringsAsFactors = FALSE
  )
}

# Inventory all durable analysis outputs with project-relative paths.
build_report_inventory <- function(outdir = ".") {
  if (!dir.exists(outdir)) stop("Report output directory not found: ", outdir)
  files <- list.files(outdir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  root <- normalizePath(outdir, mustWork = TRUE)
  rel <- sub(paste0("^", gsub("([][{}()+*^$|\\?.])", "\\\\\\1", root),
                    .Platform$file.sep, "?"), "", normalizePath(files, mustWork = FALSE))
  data.frame(
    path = rel,
    extension = tolower(tools::file_ext(rel)),
    bytes = as.numeric(file.info(files)$size),
    stringsAsFactors = FALSE
  )
}

# Build one row per report domain. Every domain receives an explicit status;
# absence is visible and can later be resolved as not_applicable or
# omitted_with_reason by a project-specific report workflow.
build_report_coverage_manifest <- function(outdir = ".",
                                           rules = default_report_coverage_rules(),
                                           overrides = NULL) {
  inventory <- build_report_inventory(outdir)
  inventory <- inventory[!grepl("(^|/)report_coverage_manifest[.]csv$", inventory$path), , drop = FALSE]
    rows <- lapply(seq_len(nrow(rules)), function(i) {
    hits <- inventory$path[grepl(rules$pattern[[i]], inventory$path, ignore.case = TRUE, perl = TRUE)]
    ext <- unique(inventory$extension[inventory$path %in% hits])
    data.frame(
      domain = rules$domain[[i]],
      label = rules$label[[i]],
      required = rules$required[[i]],
      status = if (length(hits)) "covered" else "not_available",
      report_location = if (length(hits)) "main_body_or_appendix_index" else "unresolved",
      representation = if (!length(hits)) "none" else paste(sort(ext[nzchar(ext)]), collapse = "+"),
      evidence_count = length(hits),
      evidence_files = paste(utils::head(hits, 12), collapse = ";"),
      omission_reason = if (length(hits)) "" else "No matching saved output was found.",
      stringsAsFactors = FALSE
    )
  })
  manifest <- do.call(rbind, rows)
  if (!is.null(overrides) && nrow(overrides) > 0) {
    if (!"domain" %in% colnames(overrides)) stop("Coverage overrides require a `domain` column.")
    for (i in seq_len(nrow(overrides))) {
      idx <- match(overrides$domain[[i]], manifest$domain)
      if (is.na(idx)) stop("Unknown report coverage domain: ", overrides$domain[[i]])
      for (column in setdiff(intersect(colnames(overrides), colnames(manifest)), "domain")) {
        manifest[idx, column] <- overrides[i, column]
      }
    }
  }
  manifest
}

write_report_coverage_manifest <- function(outdir = ".",
                                            path = file.path(outdir, "report_coverage_manifest.csv"),
                                            rules = default_report_coverage_rules(),
                                            overrides = NULL) {
  manifest <- build_report_coverage_manifest(outdir, rules = rules, overrides = overrides)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(manifest, path, row.names = FALSE, na = "")
  invisible(manifest)
}

validate_report_coverage <- function(manifest, require_resolved_optional = FALSE) {
  required_cols <- c("domain", "required", "status", "report_location", "evidence_count", "omission_reason")
  missing_cols <- setdiff(required_cols, colnames(manifest))
  if (length(missing_cols)) stop("Coverage manifest is missing columns: ", paste(missing_cols, collapse = ", "))
  allowed <- c("covered", "not_applicable", "omitted_with_reason", "not_available")
  if (any(!manifest$status %in% allowed)) stop("Coverage manifest contains an unsupported status.")
  unresolved_required <- manifest$required & manifest$status == "not_available"
  if (any(unresolved_required)) {
    stop("Required report domains are not covered: ", paste(manifest$domain[unresolved_required], collapse = ", "))
  }
  needs_reason <- manifest$status %in% c("not_applicable", "omitted_with_reason")
  if (any(needs_reason & !nzchar(trimws(manifest$omission_reason)))) {
    stop("Omitted/not-applicable report domains require a reason.")
  }
  covered_without_location <- manifest$status == "covered" &
    (!nzchar(trimws(manifest$report_location)) | manifest$report_location == "unresolved")
  if (any(covered_without_location)) {
    stop("Covered report domains require a main-body or appendix location.")
  }
  if (isTRUE(require_resolved_optional) && any(manifest$status == "not_available")) {
    stop("Optional report domains remain unresolved: ",
         paste(manifest$domain[manifest$status == "not_available"], collapse = ", "))
  }
  invisible(TRUE)
}

# Validate the narrative handoff used for model-assisted scientific synthesis.
# The ledger keeps generated interpretation subordinate to saved evidence.
validate_claim_evidence_ledger <- function(ledger) {
  if (is.character(ledger) && length(ledger) == 1) {
    if (!file.exists(ledger)) stop("Claim-evidence ledger not found: ", ledger)
    ledger <- utils::read.csv(ledger, stringsAsFactors = FALSE, check.names = FALSE)
  }
  required <- c("claim_id", "claim", "claim_level", "scope", "evidence",
                "source_files", "assumptions", "alternatives", "falsifier", "next_evidence")
  missing_cols <- setdiff(required, colnames(ledger))
  if (length(missing_cols)) stop("Claim-evidence ledger is missing columns: ", paste(missing_cols, collapse = ", "))
  if (anyNA(ledger$claim_id) || any(!nzchar(trimws(ledger$claim_id))) || anyDuplicated(ledger$claim_id)) {
    stop("`claim_id` must be present and unique.")
  }
  allowed_levels <- c("E0_observed", "E1_association", "E2_directional",
                      "E3_causal_model", "E4_mechanism", "E5_translational",
                      "unsupported", "superseded")
  if (any(!ledger$claim_level %in% allowed_levels)) {
    stop("Claim ledger contains unsupported evidence levels.")
  }
  if (any(!nzchar(trimws(ledger$source_files)))) stop("Every claim requires one or more source files.")
  invisible(TRUE)
}

# List PDF figures under a directory, optionally filtered by a regex pattern.
# Returns paths relative to `base` so the report can embed them with <img>.
list_report_figures <- function(dir, pattern = "\\.pdf$", base = dir) {
  if (!dir.exists(dir)) return(character(0))
  files <- list.files(dir, pattern = pattern, full.names = TRUE, recursive = TRUE)
  if (length(files) == 0) return(character(0))
  # Convert to paths relative to base for portable embedding.
  rel <- sub(paste0("^", normalizePath(base, mustWork = FALSE), .Platform$file.sep, "?"), "",
             normalizePath(files, mustWork = FALSE))
  sort(rel)
}

# Read a CSV if it exists, else return NULL (used to conditionally show tables).
read_csv_safe <- function(path, ...) {
  if (!file.exists(path)) return(NULL)
  tryCatch(utils::read.csv(path, stringsAsFactors = FALSE, ...), error = function(e) NULL)
}
