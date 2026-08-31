args <- commandArgs(trailingOnly = TRUE)
run <- normalizePath(args[[1]], mustWork = TRUE)
html <- paste(readLines(file.path(run, 'RNAseq_report.html'), warn = FALSE), collapse = '\n')
stopifnot(nchar(html) > 10000, grepl('data:image/png;base64', html, fixed = TRUE))
files <- list.files(file.path(run, '2-GSEA'), pattern = '^GSEA_.*[.]csv$', full.names = TRUE)
files <- files[!grepl('_quality[.]csv$', files)]
stopifnot(length(files) > 0)
for (file in files) {
  tab <- read.csv(file, stringsAsFactors = FALSE)
  q <- read.csv(sub('[.]csv$', '_quality.csv', file))
  stopifnot(!anyNA(tab$ID), !anyDuplicated(tab$ID), all(nzchar(tab$ID)),
            q$table_rows == nrow(tab), q$significant_terms == sum(tab$significant),
            q$valid_terms + q$unusable_rows == nrow(tab),
            all(tab$input_overlap >= 10 & tab$input_overlap <= 500),
            isTRUE(all.equal(tab$p.adjust, p.adjust(tab$pvalue, 'BH'), tolerance = 1e-10)))
  cat(basename(file), ': ', q$table_rows, ' tested, ', q$valid_terms, ' valid, ',
      q$unusable_rows, ' unusable, ', q$significant_terms, ' BH-significant\n', sep = '')
}
cat('HTML QA: nontrivial standalone report with embedded PNG previews.\n')
