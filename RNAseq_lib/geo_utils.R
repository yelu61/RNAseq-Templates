# GEO public data helpers for bulk RNA-seq templates.

# Download a GEO SeriesMatrix file to a local directory.
# Returns the GEOquery object (GSE).
download_geo_series_matrix <- function(geo_accession, destdir = "./0-Data", getGPL = TRUE) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package 'GEOquery' is required for GEO download.\nRun BiocManager::install('GEOquery')")
  }
  if (!dir.exists(destdir)) dir.create(destdir, showWarnings = FALSE, recursive = TRUE)
  GEOquery::getGEO(GEO = geo_accession, destdir = destdir, getGPL = getGPL)
}

# Parse a GEOquery GSE object into counts, phenotype, and feature annotation.
# Attempts to detect the assay matrix; for SeriesMatrix this is usually the ExpressionSet assay.
parse_geo_series_matrix <- function(gse) {
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Package 'Biobase' is required to parse GEO objects.")
  }
  # gse is typically a list of ExpressionSet objects; take the first
  eset <- if (is.list(gse)) gse[[1]] else gse

  expr <- as.data.frame(Biobase::exprs(eset), check.names = FALSE)
  pdata <- as.data.frame(Biobase::pData(eset), check.names = FALSE)
  feature <- as.data.frame(Biobase::fData(eset), check.names = FALSE)

  list(expr = expr, pdata = pdata, feature = feature)
}

# Collapse probe-level GEO expression to gene symbols by maximum average expression.
prepare_geo_counts <- function(expr, feature, probe_col = "ID", symbol_col = "GeneSymbol",
                                 keep_probes = "max") {
  if (!symbol_col %in% colnames(feature)) {
    # Try common alternatives
    symbol_col <- intersect(c("Gene Symbol", "GeneSymbol", "SYMBOL", "symbol"), colnames(feature))[1]
    if (is.na(symbol_col)) {
      stop("Gene symbol column not found in GEO feature table. Available columns:\n",
           paste(colnames(feature), collapse = ", "))
    }
  }
  if (!probe_col %in% colnames(feature)) {
    stop("Probe ID column not found in GEO feature table: ", probe_col)
  }

  df <- expr
  df[[probe_col]] <- rownames(df)
  df <- dplyr::inner_join(df, feature[, c(probe_col, symbol_col), drop = FALSE], by = probe_col)
  df <- df[!is.na(df[[symbol_col]]) & df[[symbol_col]] != "", , drop = FALSE]

  value_cols <- setdiff(colnames(df), c(probe_col, symbol_col))
  df$mean_expr <- rowMeans(df[, value_cols, drop = FALSE])
  df <- df[order(df$mean_expr, decreasing = TRUE), ]
  df <- df[!duplicated(df[[symbol_col]]), ]
  rownames(df) <- df[[symbol_col]]
  df[, value_cols, drop = FALSE]
}
