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

run_go_gsea <- function(entrez_list, org_db, ont = "ALL", min_size = 10, max_size = 500, p_cutoff = 0.05) {
  if (length(entrez_list) < min_size) {
    message("GO GSEA skipped: ranked list length ", length(entrez_list), " (min_size = ", min_size, ")")
    return(NULL)
  }
  tryCatch(clusterProfiler::gseGO(
    geneList = entrez_list,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = ont,
    minGSSize = min_size,
    maxGSSize = max_size,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff
  ), error = function(e) {
    message("GO GSEA skipped because gseGO failed: ", conditionMessage(e))
    NULL
  })
}

run_kegg_gsea <- function(entrez_list, organism, min_size = 10, max_size = 500, p_cutoff = 0.05) {
  if (length(entrez_list) < min_size) {
    message("KEGG GSEA skipped: ranked list length ", length(entrez_list), " (min_size = ", min_size, ")")
    return(NULL)
  }
  tryCatch(clusterProfiler::gseKEGG(
    geneList = entrez_list,
    organism = organism,
    minGSSize = min_size,
    maxGSSize = max_size,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff
  ), error = function(e) {
    message("KEGG GSEA skipped because gseKEGG failed: ", conditionMessage(e))
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
      genes <- genes_for_enrichment(res_list[[comp_name]], th$p_cutoff, th$log2fc, pvalue_column = pvalue_column, lfc_column = lfc_column)
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
