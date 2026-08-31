# Enrichment helpers for bulk RNA-seq templates.

map_symbols_to_entrez <- function(symbols, org_db) {
  symbols <- unique(symbols[!is.na(symbols) & symbols != ""])
  if (length(symbols) == 0) {
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  }
  mapped <- clusterProfiler::bitr(symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org_db)
  mapped <- mapped[!is.na(mapped$ENTREZID) & !duplicated(mapped$ENTREZID), ]
  mapped
}

make_entrez_ranked_list <- function(gene_list, org_db) {
  entrez_df <- map_symbols_to_entrez(names(gene_list), org_db)
  if (nrow(entrez_df) == 0) {
    return(numeric())
  }
  entrez_df$score <- gene_list[entrez_df$SYMBOL]
  entrez_df <- entrez_df[!is.na(entrez_df$score), ]
  entrez_df <- entrez_df[order(entrez_df$score, decreasing = TRUE), ]
  entrez_df <- entrez_df[!duplicated(entrez_df$ENTREZID), ]
  entrez_list <- entrez_df$score
  names(entrez_list) <- entrez_df$ENTREZID
  sort(entrez_list, decreasing = TRUE)
}

run_go_ora <- function(symbols, org_db, universe, ont = "ALL", p_cutoff = 0.05, q_cutoff = 0.2, min_genes = 5) {
  mapped <- map_symbols_to_entrez(symbols, org_db)
  if (nrow(mapped) < min_genes) {
    message("GO ORA skipped: only ", nrow(mapped), " mapped Entrez IDs (min_genes = ", min_genes, ")")
    return(NULL)
  }
  tryCatch(clusterProfiler::enrichGO(
    gene = mapped$ENTREZID,
    universe = universe,
    OrgDb = org_db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    readable = TRUE
  ), error = function(e) {
    message("GO ORA skipped because enrichGO failed: ", conditionMessage(e))
    NULL
  })
}

run_kegg_ora <- function(symbols, org_db, universe, organism, p_cutoff = 0.05, min_genes = 5) {
  mapped <- map_symbols_to_entrez(symbols, org_db)
  if (nrow(mapped) < min_genes) {
    message("KEGG ORA skipped: only ", nrow(mapped), " mapped Entrez IDs (min_genes = ", min_genes, ")")
    return(NULL)
  }
  tryCatch(clusterProfiler::enrichKEGG(
    gene = mapped$ENTREZID,
    universe = universe,
    organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff
  ), error = function(e) {
    message("KEGG ORA skipped because enrichKEGG failed: ", conditionMessage(e))
    NULL
  })
}

# Keep every returned row for audit; validity and BH significance are separate.
# Missing labels may use the known ID in displays, but missing IDs are never
# reconstructed from row names or guessed pathway names.
audit_gsea_table <- function(gsea_result, fdr_cutoff = NULL) {
  df <- if (is.null(gsea_result)) data.frame() else as.data.frame(gsea_result)
  if (is.null(fdr_cutoff)) {
    saved <- if ("fdr_cutoff" %in% names(df)) unique(df$fdr_cutoff) else numeric()
    fdr_cutoff <- if (length(saved) == 1L) saved else 0.05
  }
  if (length(fdr_cutoff) != 1L || !is.finite(fdr_cutoff) || fdr_cutoff < 0 || fdr_cutoff > 1) {
    stop("fdr_cutoff must be a finite number between 0 and 1.")
  }
  for (nm in c("ID", "Description")) if (!nm %in% names(df)) df[[nm]] <- rep(NA_character_, nrow(df))
  for (nm in c("NES", "pvalue", "p.adjust")) if (!nm %in% names(df)) df[[nm]] <- rep(NA_real_, nrow(df))
  issues <- rep("", nrow(df))
  add_issue <- function(bad, label) {
    issues[bad] <<- ifelse(nzchar(issues[bad]), paste(issues[bad], label, sep = "; "), label)
  }
  id <- as.character(df$ID)
  add_issue(is.na(id) | !nzchar(trimws(id)) | id == "NA", "missing_ID")
  add_issue(!is.na(id) & (duplicated(id) | duplicated(id, fromLast = TRUE)), "duplicate_ID")
  for (nm in c("NES", "pvalue", "p.adjust")) {
    value <- suppressWarnings(as.numeric(as.character(df[[nm]])))
    add_issue(!is.finite(value), paste0("nonfinite_", nm))
    if (nm != "NES") add_issue(is.finite(value) & (value < 0 | value > 1), paste0("out_of_range_", nm))
    df[[nm]] <- value
  }
  if ("backend_record_missing" %in% names(df)) {
    add_issue(!is.na(df$backend_record_missing) & df$backend_record_missing, "backend_record_not_returned")
  }
  valid <- !nzchar(issues)
  df$result_issue <- issues
  df$significant <- valid & df$p.adjust <= fdr_cutoff
  df$result_status <- ifelse(!valid, "unusable", ifelse(df$significant, "significant", "not_significant"))
  df$fdr_cutoff <- rep(fdr_cutoff, nrow(df))
  list(table = df, summary = data.frame(
    table_rows = nrow(df), valid_terms = sum(valid), unusable_rows = sum(!valid),
    significant_terms = sum(df$significant), fdr_cutoff = fdr_cutoff,
    stringsAsFactors = FALSE
  ))
}

significant_gsea_terms <- function(gsea_result, fdr_cutoff = NULL) {
  df <- audit_gsea_table(gsea_result, fdr_cutoff)$table
  df <- df[df$significant, , drop = FALSE]
  df[order(df$p.adjust, -abs(df$NES)), , drop = FALSE]
}

write_gsea_tables <- function(gsea_result, filename) {
  if (is.null(gsea_result)) return(invisible(NULL))
  audit <- audit_gsea_table(gsea_result)
  utils::write.csv(audit$table, filename, row.names = FALSE, na = "")
  quality_file <- sub("[.]csv$", "_quality.csv", filename, ignore.case = TRUE)
  if (identical(quality_file, filename)) quality_file <- paste0(filename, "_quality.csv")
  utils::write.csv(audit$summary, quality_file, row.names = FALSE, na = "")
  invisible(audit)
}

# Both backend generations receive the same input-overlap-filtered gene sets.
# enrichit 0.2.0 otherwise filters raw reference sizes and its nominal-P filter
# loses IDs for NA rows. NULL disables that filter; DOSE/fgsea uses cutoff 1.
run_gsea_with_reference <- function(entrez_list, reference, min_size = 10,
                                    max_size = 500, p_cutoff = 0.05, seed = 123) {
  mapping <- unique(reference@gsid2gene[, c("gsid", "gene")])
  mapping <- mapping[!is.na(mapping$gsid) & !is.na(mapping$gene), , drop = FALSE]
  reference_sizes <- table(mapping$gsid)
  mapping <- mapping[mapping$gene %in% names(entrez_list), , drop = FALSE]
  overlap <- table(mapping$gsid)
  eligible <- names(overlap)[overlap >= min_size & overlap <= max_size]
  if (!length(eligible)) {
    message("GSEA skipped: no gene sets meet the input-overlap size limits.")
    return(NULL)
  }
  reference@gsid2gene <- mapping[mapping$gsid %in% eligible, , drop = FALSE]
  if (!is.null(seed)) set.seed(seed)
  enrichit_backend <- "method" %in% names(formals(clusterProfiler::GSEA))
  result <- clusterProfiler::GSEA(
    geneList = entrez_list, gson = reference, minGSSize = min_size, maxGSSize = max_size,
    pAdjustMethod = "BH", pvalueCutoff = if (enrichit_backend) NULL else 1
  )
  if (is.null(result)) {
    result <- methods::new("gseaResult", result = data.frame(),
      geneList = entrez_list, geneSets = split(reference@gsid2gene$gene, reference@gsid2gene$gsid))
  }
  df <- as.data.frame(result)
  # Older backends may omit uncomputable terms. Their identities are known from
  # the eligible reference, but why no record was returned is not inferred.
  missing_ids <- setdiff(eligible, df$ID)
  if (length(missing_ids)) {
    missing <- data.frame(ID = missing_ids, stringsAsFactors = FALSE)
    for (nm in setdiff(names(df), "ID")) missing[[nm]] <- NA
    missing$backend_record_missing <- TRUE
    df$backend_record_missing <- rep(FALSE, nrow(df))
    df <- dplyr::bind_rows(df, missing)
  } else {
    df$backend_record_missing <- rep(FALSE, nrow(df))
  }
  known <- match(df$ID, reference@gsid2name$gsid)
  if (!"Description" %in% names(df)) df$Description <- rep(NA_character_, nrow(df))
  blank_description <- is.na(df$Description) | !nzchar(trimws(df$Description))
  df$Description[blank_description] <- reference@gsid2name$name[known[blank_description]]
  df$input_overlap <- as.integer(overlap[as.character(df$ID)])
  df$reference_size <- as.integer(reference_sizes[as.character(df$ID)])
  result@result <- audit_gsea_table(df, p_cutoff)$table
  result@params$pvalueCutoff <- if (enrichit_backend) NULL else 1
  result@params$fdr_cutoff <- p_cutoff
  result
}

run_go_gsea <- function(entrez_list, org_db, ont = "ALL", min_size = 10, max_size = 500, p_cutoff = 0.05, seed = 123) {
  if (length(entrez_list) < min_size) {
    message("GO GSEA skipped: ranked list length ", length(entrez_list), " (min_size = ", min_size, ")")
    return(NULL)
  }
  ont <- match.arg(toupper(ont), c("ALL", "BP", "MF", "CC"))
  tryCatch({
    reference <- clusterProfiler::gson_GO(org_db, keytype = "ENTREZID", ont = ont)
    result <- run_gsea_with_reference(entrez_list, reference, min_size, max_size, p_cutoff, seed)
    if (!is.null(result)) {
      result@organism <- reference@species
      result@setType <- ont
      result@keytype <- "ENTREZID"
      result@result$ONTOLOGY <- unname(AnnotationDbi::Ontology(GO.db::GOTERM)[result@result$ID])
    }
    result
  }, error = function(e) {
    message("GO GSEA skipped because analysis failed: ", conditionMessage(e))
    NULL
  })
}

run_kegg_gsea <- function(entrez_list, organism, min_size = 10, max_size = 500, p_cutoff = 0.05, seed = 123) {
  if (length(entrez_list) < min_size) {
    message("KEGG GSEA skipped: ranked list length ", length(entrez_list), " (min_size = ", min_size, ")")
    return(NULL)
  }
  tryCatch({
    reference <- if (inherits(organism, "GSON")) organism else clusterProfiler::gson_KEGG(organism)
    result <- run_gsea_with_reference(entrez_list, reference, min_size, max_size, p_cutoff, seed)
    if (!is.null(result)) {
      result@organism <- reference@species
      result@setType <- "KEGG"
      result@keytype <- reference@keytype
    }
    result
  }, error = function(e) {
    message("KEGG GSEA skipped because analysis failed: ", conditionMessage(e))
    NULL
  })
}

run_threshold_ora <- function(res_list, threshold_grid, org_db, universe, organism, pvalue_column = "padj", lfc_column = "log2FoldChange", outdir = "2-GSEA", plotdir = "3-Visualization") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plotdir, showWarnings = FALSE, recursive = TRUE)
  summary_rows <- list()
  ora_results <- list()

  for (i in seq_len(nrow(threshold_grid))) {
    th <- threshold_grid[i, ]
    th_outdir <- file.path(outdir, th$name)
    th_plotdir <- file.path(plotdir, th$name)
    dir.create(th_outdir, showWarnings = FALSE, recursive = TRUE)
    dir.create(th_plotdir, showWarnings = FALSE, recursive = TRUE)
    ora_results[[th$name]] <- list()

    for (comp_name in names(res_list)) {
      genes <- genes_for_enrichment(res_list[[comp_name]], th$p_cutoff, th$log2fc, pvalue_column = threshold_p_column(th, pvalue_column), lfc_column = lfc_column)
      ora_results[[th$name]][[comp_name]] <- list()

      ego <- run_go_ora(genes$sig, org_db = org_db, universe = universe)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        utils::write.csv(as.data.frame(ego), file.path(th_outdir, paste0("GO_ORA_", comp_name, ".csv")), row.names = FALSE)
        plot_enrich_suite_pdf(ego, file.path(th_plotdir, paste0("GO_ORA_", comp_name)), paste("GO ORA -", comp_name, "-", th$name))
        ora_results[[th$name]][[comp_name]]$go <- ego
        go_terms <- nrow(as.data.frame(ego))
      } else {
        go_terms <- 0
      }

      ego_up <- run_go_ora(genes$up, org_db = org_db, universe = universe)
      ego_down <- run_go_ora(genes$down, org_db = org_db, universe = universe)
      go_up_terms <- if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) nrow(as.data.frame(ego_up)) else 0
      go_down_terms <- if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) nrow(as.data.frame(ego_down)) else 0
      if (go_up_terms > 0) {
        utils::write.csv(as.data.frame(ego_up), file.path(th_outdir, paste0("GO_ORA_UP_", comp_name, ".csv")), row.names = FALSE)
        plot_enrich_suite_pdf(ego_up, file.path(th_plotdir, paste0("GO_ORA_UP_", comp_name)), paste("GO ORA UP -", comp_name, "-", th$name))
        ora_results[[th$name]][[comp_name]]$go_up <- ego_up
      }
      if (go_down_terms > 0) {
        utils::write.csv(as.data.frame(ego_down), file.path(th_outdir, paste0("GO_ORA_DOWN_", comp_name, ".csv")), row.names = FALSE)
        plot_enrich_suite_pdf(ego_down, file.path(th_plotdir, paste0("GO_ORA_DOWN_", comp_name)), paste("GO ORA DOWN -", comp_name, "-", th$name))
        ora_results[[th$name]][[comp_name]]$go_down <- ego_down
      }
      if (go_up_terms > 0 || go_down_terms > 0) {
        plot_enrich_bidirectional_barplot_pdf(
          ego_up, ego_down,
          file.path(th_plotdir, paste0("GO_ORA_bidirectional_", comp_name, ".pdf")),
          paste("GO ORA UP/DOWN -", comp_name, "-", th$name)
        )
      }

      ekegg <- run_kegg_ora(genes$sig, org_db = org_db, universe = universe, organism = organism)
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        utils::write.csv(as.data.frame(ekegg), file.path(th_outdir, paste0("KEGG_ORA_", comp_name, ".csv")), row.names = FALSE)
        plot_enrich_suite_pdf(ekegg, file.path(th_plotdir, paste0("KEGG_ORA_", comp_name)), paste("KEGG ORA -", comp_name, "-", th$name))
        ora_results[[th$name]][[comp_name]]$kegg <- ekegg
        kegg_terms <- nrow(as.data.frame(ekegg))
      } else {
        kegg_terms <- 0
      }

      ekegg_up <- run_kegg_ora(genes$up, org_db = org_db, universe = universe, organism = organism)
      ekegg_down <- run_kegg_ora(genes$down, org_db = org_db, universe = universe, organism = organism)
      kegg_up_terms <- if (!is.null(ekegg_up) && nrow(as.data.frame(ekegg_up)) > 0) nrow(as.data.frame(ekegg_up)) else 0
      kegg_down_terms <- if (!is.null(ekegg_down) && nrow(as.data.frame(ekegg_down)) > 0) nrow(as.data.frame(ekegg_down)) else 0
      if (kegg_up_terms > 0) {
        utils::write.csv(as.data.frame(ekegg_up), file.path(th_outdir, paste0("KEGG_ORA_UP_", comp_name, ".csv")), row.names = FALSE)
        plot_enrich_suite_pdf(ekegg_up, file.path(th_plotdir, paste0("KEGG_ORA_UP_", comp_name)), paste("KEGG ORA UP -", comp_name, "-", th$name))
        ora_results[[th$name]][[comp_name]]$kegg_up <- ekegg_up
      }
      if (kegg_down_terms > 0) {
        utils::write.csv(as.data.frame(ekegg_down), file.path(th_outdir, paste0("KEGG_ORA_DOWN_", comp_name, ".csv")), row.names = FALSE)
        plot_enrich_suite_pdf(ekegg_down, file.path(th_plotdir, paste0("KEGG_ORA_DOWN_", comp_name)), paste("KEGG ORA DOWN -", comp_name, "-", th$name))
        ora_results[[th$name]][[comp_name]]$kegg_down <- ekegg_down
      }
      if (kegg_up_terms > 0 || kegg_down_terms > 0) {
        plot_enrich_bidirectional_barplot_pdf(
          ekegg_up, ekegg_down,
          file.path(th_plotdir, paste0("KEGG_ORA_bidirectional_", comp_name, ".pdf")),
          paste("KEGG ORA UP/DOWN -", comp_name, "-", th$name)
        )
      }

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Threshold = th$name,
        Comparison = comp_name,
        SignificantGenes = length(genes$sig),
        UpGenes = length(genes$up),
        DownGenes = length(genes$down),
        GO_terms = go_terms,
        GO_up_terms = go_up_terms,
        GO_down_terms = go_down_terms,
        KEGG_terms = kegg_terms,
        KEGG_up_terms = kegg_up_terms,
        KEGG_down_terms = kegg_down_terms,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  utils::write.csv(summary_df, file.path(outdir, "ORA_threshold_summary.csv"), row.names = FALSE)
  list(results = ora_results, summary = summary_df)
}

# Convert a clusterProfiler enrichResult/gseaResult (or compareClusterResult)
# to a plain data frame with standard columns.
enrich_result_to_df <- function(enrich_result) {
  if (is.null(enrich_result)) return(data.frame())
  df <- tryCatch(as.data.frame(enrich_result), error = function(e) data.frame())
  if (nrow(df) == 0) return(df)
  required <- c("ID", "Description", "p.adjust")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("enrich_result is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!"Count" %in% colnames(df) && "setSize" %in% colnames(df)) {
    df$Count <- df$setSize
  }
  if (!"Count" %in% colnames(df)) {
    df$Count <- NA_integer_
  }
  if (!"ONTOLOGY" %in% colnames(df)) {
    # For enrichResult objects with a single ontology, fill from the object slot.
    ontology <- tryCatch(slot(enrich_result, "ontology"), error = function(e) NA_character_)
    df$ONTOLOGY <- if (length(ontology) == 1 && !is.na(ontology)) ontology else NA_character_
  }
  df
}

# Build a combined enrichment data frame across multiple comparisons.
# result_map: named list where each element is an enrichResult/gseaResult or data frame.
build_multi_comparison_enrich_df <- function(result_map,
                                             comparison_col = "comparison",
                                             direction_col = NULL,
                                             ontology_filter = NULL) {
  rows <- list()
  for (comp_name in names(result_map)) {
    obj <- result_map[[comp_name]]
    df <- enrich_result_to_df(obj)
    if (nrow(df) == 0) next
    df[[comparison_col]] <- comp_name
    if (!is.null(direction_col)) {
      df[[direction_col]] <- NA_character_
    }
    rows[[length(rows) + 1]] <- df
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) return(out)
  if (!is.null(ontology_filter) && "ONTOLOGY" %in% colnames(out)) {
    # Keep NA-ontology rows: enrichResult objects already carry a real ONTOLOGY,
    # but a plain data.frame (e.g. read back from a saved ORA csv) gets NA here
    # and would otherwise be silently dropped by the filter.
    out <- out[is.na(out$ONTOLOGY) | out$ONTOLOGY %in% ontology_filter, , drop = FALSE]
  }
  out
}
