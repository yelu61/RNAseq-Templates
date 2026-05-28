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
  if (nrow(mapped) < min_genes) return(NULL)
  clusterProfiler::enrichGO(
    gene = mapped$ENTREZID,
    universe = universe,
    OrgDb = org_db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff,
    qvalueCutoff = q_cutoff,
    readable = TRUE
  )
}

run_kegg_ora <- function(symbols, org_db, universe, organism, p_cutoff = 0.05, min_genes = 5) {
  mapped <- map_symbols_to_entrez(symbols, org_db)
  if (nrow(mapped) < min_genes) return(NULL)
  clusterProfiler::enrichKEGG(
    gene = mapped$ENTREZID,
    universe = universe,
    organism = organism,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff
  )
}

run_go_gsea <- function(entrez_list, org_db, ont = "ALL", min_size = 10, max_size = 500, p_cutoff = 0.05) {
  if (length(entrez_list) < min_size) return(NULL)
  clusterProfiler::gseGO(
    geneList = entrez_list,
    OrgDb = org_db,
    keyType = "ENTREZID",
    ont = ont,
    minGSSize = min_size,
    maxGSSize = max_size,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff
  )
}

run_kegg_gsea <- function(entrez_list, organism, min_size = 10, max_size = 500, p_cutoff = 0.05) {
  if (length(entrez_list) < min_size) return(NULL)
  clusterProfiler::gseKEGG(
    geneList = entrez_list,
    organism = organism,
    minGSSize = min_size,
    maxGSSize = max_size,
    pAdjustMethod = "BH",
    pvalueCutoff = p_cutoff
  )
}

run_threshold_ora <- function(res_list, threshold_grid, org_db, universe, organism, outdir = "2-GSEA", plotdir = "3-Visualization") {
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
      genes <- genes_for_enrichment(res_list[[comp_name]], th$padj, th$log2fc)
      ora_results[[th$name]][[comp_name]] <- list()

      ego <- run_go_ora(genes$sig, org_db = org_db, universe = universe)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        utils::write.csv(as.data.frame(ego), file.path(th_outdir, paste0("GO_ORA_", comp_name, ".csv")), row.names = FALSE)
        plot_enrich_dotplot(ego, file.path(th_plotdir, paste0("GO_dotplot_", comp_name, ".pdf")), paste("GO ORA -", comp_name, "-", th$name))
        ora_results[[th$name]][[comp_name]]$go <- ego
        go_terms <- nrow(as.data.frame(ego))
      } else {
        go_terms <- 0
      }

      ekegg <- run_kegg_ora(genes$sig, org_db = org_db, universe = universe, organism = organism)
      if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
        utils::write.csv(as.data.frame(ekegg), file.path(th_outdir, paste0("KEGG_ORA_", comp_name, ".csv")), row.names = FALSE)
        plot_enrich_dotplot(ekegg, file.path(th_plotdir, paste0("KEGG_dotplot_", comp_name, ".pdf")), paste("KEGG ORA -", comp_name, "-", th$name))
        ora_results[[th$name]][[comp_name]]$kegg <- ekegg
        kegg_terms <- nrow(as.data.frame(ekegg))
      } else {
        kegg_terms <- 0
      }

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Threshold = th$name,
        Comparison = comp_name,
        SignificantGenes = length(genes$sig),
        GO_terms = go_terms,
        KEGG_terms = kegg_terms,
        stringsAsFactors = FALSE
      )
    }
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  utils::write.csv(summary_df, file.path(outdir, "ORA_threshold_summary.csv"), row.names = FALSE)
  list(results = ora_results, summary = summary_df)
}
