# Report generation helpers for RNAseq-Templates.
# Renders a self-contained HTML report from the CSV/PDF outputs that the
# notebooks already produce, without re-running any analysis.

.report_or <- function(x, fallback) if (is.null(x)) fallback else x

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
  overrides <- NULL
  if (file.exists(coverage_file)) {
    previous <- tryCatch(utils::read.csv(coverage_file, stringsAsFactors = FALSE,
                                         check.names = FALSE), error = function(e) NULL)
    if (!is.null(previous) && all(c("domain", "status") %in% colnames(previous))) {
      decision_columns <- intersect(c("domain", "status", "report_location", "omission_reason"),
                                    colnames(previous))
      overrides <- previous[previous$status %in% c("not_applicable", "omitted_with_reason"),
                            decision_columns, drop = FALSE]
    }
  }
  write_report_coverage_manifest(outdir, coverage_file, overrides = overrides)
  if (!file.exists(file.path(outdir, "report_review_checklist.csv"))) {
    scaffold_report_review(outdir)
  }
  ledger_file <- file.path(outdir, "claim_evidence_ledger.csv")
  if (file.exists(ledger_file)) validate_claim_evidence_ledger(ledger_file, base_dir = outdir)

  template_ext <- tolower(tools::file_ext(template))
  quarto_available <- requireNamespace("quarto", quietly = TRUE) && nzchar(Sys.which("quarto"))
  if (!quarto_available && template_ext == "qmd" && requireNamespace("rmarkdown", quietly = TRUE)) {
    # Out-of-the-box fallback: a sibling .Rmd twin of the .qmd template renders
    # via rmarkdown when the Quarto CLI is not installed (the common case on
    # fresh machines; surfaced by three real projects that silently skipped the
    # report). The twin must carry the same body with rmarkdown YAML.
    sibling <- sub("[.]qmd$", ".Rmd", template, ignore.case = TRUE)
    if (file.exists(sibling)) {
      message("Quarto CLI not found; rendering rmarkdown twin instead: ", sibling)
      template <- sibling
      template_ext <- "rmd"
    }
  }
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
    stop("The selected .qmd report requires the Quarto CLI and R 'quarto' package, ",
         "or an .Rmd twin next to the template plus the 'rmarkdown' package.")
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
      "(sessionInfo|analysis_manifest|run_manifest|config[.]R)",
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
      editable <- c("status", "report_location", "omission_reason")
      for (column in intersect(colnames(overrides), editable)) {
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
validate_claim_evidence_ledger <- function(ledger, base_dir = NULL) {
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
  if (!is.null(base_dir)) {
    missing_sources <- unique(unlist(lapply(ledger$source_files, function(value) {
      files <- trimws(strsplit(value, ";", fixed = TRUE)[[1]])
      files <- files[nzchar(files)]
      files[!file.exists(file.path(base_dir, files))]
    })))
    if (length(missing_sources)) {
      stop("Claim ledger references missing source files: ", paste(missing_sources, collapse = ", "))
    }
  }
  invisible(TRUE)
}

# Independently reconcile the text-summary headline DEG counts against the
# machine-readable threshold summary. This catches stale prose after a rerun.
validate_report_headline_values <- function(outdir = ".") {
  summary_file <- file.path(outdir, "Analysis_summary.txt")
  deg_file <- file.path(outdir, "1-DEG", "DEG_threshold_summary.csv")
  config_file <- file.path(outdir, "0-Config", "analysis_config_used.R")
  if (!file.exists(summary_file) || !file.exists(deg_file)) {
    stop("Headline validation requires Analysis_summary.txt and DEG_threshold_summary.csv.")
  }
  primary <- "standard"
  if (file.exists(config_file)) {
    config_env <- new.env(parent = baseenv())
    tryCatch(sys.source(config_file, envir = config_env), error = function(e) NULL)
    if (exists("DEFAULT_THRESHOLD", envir = config_env, inherits = FALSE)) {
      primary <- get("DEFAULT_THRESHOLD", envir = config_env, inherits = FALSE)
    }
  }
  deg <- utils::read.csv(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("Threshold", "Comparison", "UP", "DOWN", "Total")
  missing <- setdiff(required, colnames(deg))
  if (length(missing)) stop("DEG threshold summary is missing columns: ", paste(missing, collapse = ", "))
  expected <- deg[deg$Threshold == primary, required[-1], drop = FALSE]
  if (!nrow(expected)) stop("No DEG threshold-summary rows found for primary tier: ", primary)

  lines <- readLines(summary_file, warn = FALSE)
  pattern <- "^\\s*-\\s+(.+):\\s+([0-9]+) DEGs \\(Up: ([0-9]+), Down: ([0-9]+)\\)"
  matches <- regexec(pattern, lines, perl = TRUE)
  parsed <- regmatches(lines, matches)
  parsed <- parsed[lengths(parsed) == 5L]
  if (!length(parsed)) stop("No parseable DEG headline counts found in Analysis_summary.txt.")
  observed <- do.call(rbind, lapply(parsed, function(x) {
    data.frame(Comparison = x[[2]], Total = as.integer(x[[3]]),
               UP = as.integer(x[[4]]), DOWN = as.integer(x[[5]]), stringsAsFactors = FALSE)
  }))
  idx <- match(expected$Comparison, observed$Comparison)
  if (anyNA(idx)) stop("Analysis_summary.txt is missing comparisons: ",
                       paste(expected$Comparison[is.na(idx)], collapse = ", "))
  mismatch <- expected$Total != observed$Total[idx] |
    expected$UP != observed$UP[idx] | expected$DOWN != observed$DOWN[idx]
  if (any(mismatch)) stop("Headline DEG counts disagree with DEG_threshold_summary.csv: ",
                          paste(expected$Comparison[mismatch], collapse = ", "))
  invisible(TRUE)
}

# Verify that the rendered HTML did not silently lose report figures. The
# report template converts referenced PDFs to embedded PNG data URIs; any
# conversion marker is therefore a hard asset failure for publication.
validate_report_html_assets <- function(report_file, expect_figures = FALSE) {
  if (!file.exists(report_file)) stop("Rendered HTML report not found: ", report_file)
  html <- paste(readLines(report_file, warn = FALSE), collapse = "\n")
  markers <- c("Figure not found:", "PDF preview unavailable", "PDF preview failed")
  present <- markers[vapply(markers, grepl, logical(1), x = html, fixed = TRUE)]
  if (length(present)) {
    stop("Rendered HTML contains figure-asset failure markers: ", paste(present, collapse = ", "))
  }
  hits <- gregexpr("data:image/(png|jpeg);base64,", html, perl = TRUE)[[1]]
  n_embedded <- if (identical(hits, -1L)) 0L else length(hits)
  if (isTRUE(expect_figures) && n_embedded == 0L) {
    stop("Run contains PDF figures but the HTML has no embedded figure previews.")
  }
  invisible(list(embedded_figures = n_embedded, passed = TRUE))
}

# Human scientific-review gate. Rendering is allowed while items are pending;
# publication requires every item to be approved or explicitly not applicable.
report_review_items <- function() {
  data.frame(
    item_id = c(
      "design_unit", "sample_qc", "contrast_direction", "deg_inference",
      "ora_scope", "gsea_interpretation", "geneset_statistics",
      "claim_evidence", "alternatives_limitations", "figure_integrity"
    ),
    review_item = c(
      "Experimental unit, design formula and covariates are correct",
      "Sample QC flags and exclusions are justified and documented",
      "Every contrast direction and reference group is verified",
      "Primary FDR/effect-size threshold and sensitivity tiers are labelled",
      "ORA direction, input list and background universe are appropriate",
      "GSEA rank metric, NES direction and leading-edge evidence are reviewed",
      "GSVA/TME statistics, multiplicity and gene-set provenance are reviewed",
      "Headline claims trace to saved files and calibrated evidence levels",
      "Alternative explanations, limitations and falsifiers are stated",
      "Primary figures are readable, non-empty and match their source tables"
    ),
    status = "pending",
    reviewer = "",
    reviewed_at = "",
    notes = "",
    stringsAsFactors = FALSE
  )
}

scaffold_report_review <- function(outdir = ".", overwrite = FALSE) {
  if (!dir.exists(outdir)) stop("Report output directory not found: ", outdir)
  path <- file.path(outdir, "report_review_checklist.csv")
  if (file.exists(path) && !isTRUE(overwrite)) return(invisible(path))
  utils::write.csv(report_review_items(), path, row.names = FALSE, na = "")
  invisible(path)
}

validate_report_review <- function(checklist, require_signoff = TRUE) {
  if (is.character(checklist) && length(checklist) == 1L) {
    if (!file.exists(checklist)) stop("Report review checklist not found: ", checklist)
    checklist <- utils::read.csv(checklist, stringsAsFactors = FALSE, check.names = FALSE)
  }
  required <- c("item_id", "review_item", "status", "reviewer", "reviewed_at", "notes")
  missing <- setdiff(required, colnames(checklist))
  if (length(missing)) stop("Report review checklist is missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(checklist$item_id) || any(!nzchar(trimws(checklist$item_id)))) {
    stop("Report review item_id values must be present and unique.")
  }
  allowed <- c("pending", "approved", "not_applicable")
  if (any(!checklist$status %in% allowed)) stop("Report review checklist contains an unsupported status.")
  if (isTRUE(require_signoff)) {
    if (any(checklist$status == "pending")) {
      stop("Report review remains pending: ", paste(checklist$item_id[checklist$status == "pending"], collapse = ", "))
    }
    signed <- checklist$status == "approved"
    if (any(signed & (!nzchar(trimws(checklist$reviewer)) | !nzchar(trimws(checklist$reviewed_at))))) {
      stop("Approved report-review items require reviewer and reviewed_at.")
    }
    na_without_reason <- checklist$status == "not_applicable" & !nzchar(trimws(checklist$notes))
    if (any(na_without_reason)) stop("Not-applicable report-review items require a reason in notes.")
  }
  invisible(TRUE)
}

# Deterministic technical QA plus the optional human sign-off gate. Results are
# always written to report_validation.csv before a failed validation stops.
validate_analysis_report <- function(outdir = ".",
                                     report_file = file.path(outdir, "RNAseq_report.html"),
                                     publish = FALSE) {
  checks <- list()
  add_check <- function(check, passed, detail) {
    checks[[length(checks) + 1L]] <<- data.frame(
      check = check, passed = isTRUE(passed), detail = detail, stringsAsFactors = FALSE
    )
  }

  html_ok <- file.exists(report_file) && is.finite(file.info(report_file)$size) && file.info(report_file)$size > 10000
  add_check("rendered_html", html_ok, if (html_ok) "HTML exists and is non-trivial." else "HTML missing or too small.")

  pdfs <- list.files(outdir, pattern = "[.]pdf$", recursive = TRUE, full.names = TRUE)
  html_asset_error <- tryCatch({
    validate_report_html_assets(report_file, expect_figures = length(pdfs) > 0L); NULL
  }, error = function(e) conditionMessage(e))
  add_check("html_figure_assets", is.null(html_asset_error),
            .report_or(html_asset_error, "No missing-preview markers; expected figures are embedded."))

  coverage_file <- file.path(outdir, "report_coverage_manifest.csv")
  coverage_error <- tryCatch({
    coverage <- utils::read.csv(coverage_file, stringsAsFactors = FALSE, check.names = FALSE)
    validate_report_coverage(coverage, require_resolved_optional = publish)
    NULL
  }, error = function(e) conditionMessage(e))
  add_check("coverage_contract", is.null(coverage_error), .report_or(coverage_error, "Coverage contract passed."))

  ledger_file <- file.path(outdir, "claim_evidence_ledger.csv")
  ledger_error <- if (file.exists(ledger_file)) tryCatch({
    validate_claim_evidence_ledger(ledger_file, base_dir = outdir); NULL
  }, error = function(e) conditionMessage(e)) else NULL
  add_check("claim_ledger", is.null(ledger_error), .report_or(ledger_error, "No invalid claim ledger detected."))

  reproducibility <- file.exists(file.path(outdir, "sessionInfo.txt")) &&
    (file.exists(file.path(outdir, "0-Config", "analysis_config_used.R")) ||
       file.exists(file.path(outdir, "run_manifest.csv")))
  add_check("reproducibility_metadata", reproducibility,
            if (reproducibility) "Session and configuration provenance found." else "Missing session or configuration provenance.")

  headline_error <- tryCatch({
    validate_report_headline_values(outdir); NULL
  }, error = function(e) conditionMessage(e))
  add_check("headline_recomputation", is.null(headline_error),
            .report_or(headline_error, "Headline DEG counts match the saved threshold table."))

  bad_pdfs <- pdfs[is.na(file.info(pdfs)$size) | file.info(pdfs)$size <= 1500]
  add_check("pdf_structure", !length(bad_pdfs),
            if (!length(bad_pdfs)) paste(length(pdfs), "PDF files passed the non-empty check.") else
              paste("Empty/incomplete PDFs:", paste(basename(bad_pdfs), collapse = ", ")))

  review_file <- file.path(outdir, "report_review_checklist.csv")
  review_error <- tryCatch({
    validate_report_review(review_file, require_signoff = publish); NULL
  }, error = function(e) conditionMessage(e))
  add_check("scientific_review", is.null(review_error), .report_or(
    review_error,
    if (publish) "Scientific review signed off." else "Checklist schema is valid; pending items are allowed for draft rendering."
  ))

  result <- do.call(rbind, checks)
  utils::write.csv(result, file.path(outdir, "report_validation.csv"), row.names = FALSE, na = "")
  if (any(!result$passed)) {
    stop(if (publish) "Report is not ready to publish: " else "Report validation failed: ",
         paste(result$check[!result$passed], collapse = ", "))
  }
  invisible(result)
}

# Canonical report sections that accept a project annotation. Keys map to
# report_interpretation/<key>.md files rendered in place by the report
# template's show_interpretation() calls; labels describe the target section.
report_interpretation_sections <- function() {
  c(
    scope = "Scope, samples, and inferential definitions",
    qc = "Sample quality",
    deg = "Differential expression",
    ora = "Over-representation analysis",
    gsea = "Rank-based enrichment (GSEA)",
    custom_genesets = "Custom gene-set scoring",
    cross_contrast = "Cross-contrast comparison",
    synthesis = "Evidence-linked synthesis and hypotheses",
    limitations = "Limitations, uncertainty, and robustness"
  )
}

# Remove whole-line HTML comments (scaffold guidance) from annotation lines;
# the report template applies the same rule before rendering.
strip_interpretation_comments <- function(lines) {
  lines[!grepl("^\\s*<!--.*-->\\s*$", lines)]
}

# Create report_interpretation/ under `outdir` with one starter file per
# report section. Existing files are never touched unless overwrite = TRUE,
# so re-scaffolding cannot destroy written annotations. A starter file
# contains only whole-line HTML comments, which the report template strips;
# it therefore renders nothing until real prose is added. Annotations are
# subordinate to saved evidence and to the claim-evidence ledger.
#
# @return Invisibly, a list with the directory, created keys and skipped keys.
scaffold_report_interpretation <- function(outdir = ".", overwrite = FALSE) {
  if (!dir.exists(outdir)) stop("Report output directory not found: ", outdir)
  dir <- file.path(outdir, "report_interpretation")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  sections <- report_interpretation_sections()
  created <- character(0)
  skipped <- character(0)
  for (key in names(sections)) {
    path <- file.path(dir, paste0(key, ".md"))
    if (file.exists(path) && !isTRUE(overwrite)) {
      skipped <- c(skipped, key)
      next
    }
    writeLines(c(
      sprintf("<!-- Section: %s -->", sections[[key]]),
      "<!-- Write the reviewed, project-specific annotation for this section below. -->",
      "<!-- Whole-line HTML comments (like these) never render; a file containing -->",
      "<!-- only comments is skipped entirely. Inline `r code` is evaluated at render. -->",
      ""
    ), path)
    created <- c(created, key)
  }
  invisible(list(dir = dir, created = created, skipped = skipped))
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

# Report all saved rows, but rank only identifiable terms with finite statistics.
# Counts refer to these files; older nominal-P-filtered tables do not establish
# the complete tested family size.
gsea_report_tables <- function(gsea_dir, contrast, n = 5) {
  files <- c(GO = file.path(gsea_dir, paste0("GSEA_GO_", contrast, ".csv")),
             KEGG = file.path(gsea_dir, paste0("GSEA_KEGG_", contrast, ".csv")))
  summaries <- terms <- list()
  for (db in names(files)) {
    df <- read_csv_safe(files[[db]])
    if (is.null(df)) next
    audit <- audit_gsea_table(df)
    summaries[[db]] <- cbind(Database = db, audit$summary)
    df <- audit$table[audit$table$result_status != "unusable", , drop = FALSE]
    pick <- function(x, direction) {
      if (!nrow(x)) return(NULL)
      x <- utils::head(x[order(x$p.adjust, -abs(x$NES)), , drop = FALSE], n)
      label <- x$Description
      missing <- is.na(label) | !nzchar(trimws(label))
      label[missing] <- x$ID[missing]
      data.frame(Database = db, Direction = direction, ID = x$ID, Term = label,
                 NES = round(x$NES, 2), BH_FDR = signif(x$p.adjust, 3),
                 Significance = ifelse(x$significant, "significant", "not significant"),
                 stringsAsFactors = FALSE)
    }
    terms[[db]] <- rbind(pick(df[df$NES > 0, , drop = FALSE], "Positive NES"),
                         pick(df[df$NES < 0, , drop = FALSE], "Negative NES"))
  }
  list(summary = do.call(rbind, summaries), terms = do.call(rbind, terms))
}
