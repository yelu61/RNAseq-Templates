# Synthetic layout fixtures, not biological findings. Run from backend root.
source('RNAseq_lib/plot_utils.R')
source('RNAseq_lib/enrichment_utils.R')
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else tempfile('rnaseq-bar-label-qa-')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
terms <- c('Regulation of leukocyte migration and activation in response to inflammatory stimuli',
           'Short bar with a long pathway description that must remain legible',
           'Mitochondrial translation',
           'Antigen processing and presentation of peptide antigen via major histocompatibility complex')
ora <- data.frame(ID=paste0('GO:',1:4), Description=terms, Count=c(40,10,22,30),
                  GeneRatio=c('40/100','10/100','22/100','30/100'),
                  BgRatio=c('40/1000','100/1000','60/1000','80/1000'),
                  FoldEnrichment=c(8,.12,3,6),
                  pvalue=c(1e-20,.025,1e-5,1e-10), p.adjust=c(1e-18,.04,2e-4,1e-8))
down <- ora
# A much smaller negative range exercises unequal direction magnitudes.
down$ID <- paste0('GO:',5:8)
down$FoldEnrichment <- c(.8,.015,.35,.6)
p <- plot_enrich_barplot_pdf(ora,file.path(out,'ordinary-ora.pdf'),'ORA: existing in-bar labels (synthetic layout)')
p <- plot_enrich_bidirectional_barplot_pdf(ora,down,file.path(out,'bidir-both.pdf'),'Bidirectional ORA: unequal ranges (synthetic layout)')
p <- plot_enrich_bidirectional_barplot_pdf(ora,NULL,file.path(out,'bidir-up-only.pdf'),'ORA UP only: short and long bars (synthetic layout)')
p <- plot_enrich_bidirectional_barplot_pdf(NULL,down,file.path(out,'bidir-down-only.pdf'),'ORA DOWN only: short and long bars (synthetic layout)')
gsea <- data.frame(ID=paste0('SET:',1:8), Description=rep(terms,2), NES=c(3,.06,1.4,2.5,-2.8,-.045,-1.6,-2.1),
                   pvalue=rep(c(1e-20,.004,1e-5,1e-10),2),p.adjust=rep(c(1e-18,.04,2e-4,1e-8),2))
p <- plot_gsea_nes_barplot_pdf(gsea,file.path(out,'gsea-both.pdf'),'GSEA: positive and negative NES (synthetic layout)')
p <- plot_gsea_nes_barplot_pdf(gsea[gsea$NES>0,],file.path(out,'gsea-positive-only.pdf'),'GSEA: positive NES only (synthetic layout)')
p <- plot_gsea_nes_barplot_pdf(gsea[gsea$NES<0,],file.path(out,'gsea-negative-only.pdf'),'GSEA: negative NES only (synthetic layout)')
for (f in list.files(out,pattern='[.]pdf$',full.names=TRUE)) {
  status <- system2('pdftoppm',c('-png','-singlefile','-r','120',shQuote(f),shQuote(sub('[.]pdf$','',f))))
  stopifnot(status == 0)
}
writeLines(c('Synthetic layout fixtures only; these are not biological results.',
             'Single-sided terms overlay their bars; two-sided terms sit across zero from their bars.',
             'GSEA FDR stays in a separate outer column on the bar side.',
             'No analysis values or GSEA significance thresholds were changed.',
             'All seven PDFs were rasterized at 120 dpi for visual QA.'),file.path(out,'README.txt'))

cat("Synthetic preview directory:", out, "\n")
