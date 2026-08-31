root <- normalizePath(getwd())
source(file.path(root, 'RNAseq_lib', 'tme_utils.R'))
base <- tempfile('tme-native-cli-', tmpdir='/tmp')
dir.create(base)
writeLines(base, '/tmp/tme-native-cli-last-path.txt')
for (species in c('human', 'mouse')) {
  work <- file.path(base, species)
  dir.create(work)
  sig_path <- file.path(root, 'references', 'CIBERSORT', if (species=='human') 'LM22.txt' else 'cibersort_mouse_22.csv')
  sig <- read_cibersort_signature(sig_path)
  set.seed(37)
  weights <- matrix(runif(ncol(sig)*6, .1, 1), ncol(sig), 6)
  weights <- sweep(weights, 2, colSums(weights), '/')
  tpm <- sig %*% weights
  tpm <- sweep(tpm, 2, colSums(tpm), '/') * 1e6
  colnames(tpm) <- paste0('S', 1:6)
  expr_path <- file.path(work, 'input.csv')
  write.csv(data.frame(gene_name=rownames(tpm), log2(tpm+1), check.names=FALSE), expr_path, row.names=FALSE)
  meta_path <- file.path(work, 'metadata.csv')
  write.csv(data.frame(sample=colnames(tpm), condition=rep(c('Control','Treatment'), each=3)), meta_path, row.names=FALSE)
  settings <- list(INPUT_MODE='expression', EXPR_FILE=expr_path, EXPR_UNIT='log2_tpm',
    META_FILE=meta_path, GENE_COLUMN='gene_name', SAMPLE_COLUMN='sample', GROUP_COLUMN='condition',
    GROUP_LEVELS=c('Control','Treatment'), SPECIES=species, GROUP_COLORS=NULL,
    RUN_ESTIMATE=FALSE, RUN_IOBR=FALSE, RUN_SSGSEA=TRUE, RUN_CIBERSORT=TRUE,
    CIBERSORT_SCRIPT=file.path(root,'references','CIBERSORT','CIBERSORT.R'), CIBERSORT_SIGNATURE=sig_path,
    CIBERSORT_PERM=0, CIBERSORT_QN=FALSE, RUN_CIBERSORT_COMPARISON=TRUE, ORTHOLOG_CACHE=NULL,
    OUTDIR=file.path(work,'run'), GENERATE_HTML_REPORT=FALSE, REPORT_TITLE='Native CIBERSORT validation')
  config <- file.path(work,'config.R')
  writeLines(c('options(parallelly.availableCores.methods="fallback", parallelly.availableCores.fallback=1L)',
    unlist(lapply(names(settings), function(n) paste0(n,' <- ', paste(deparse(settings[[n]]),collapse='\n'))))), config)
  runner <- file.path(work, 'run_analysis.R')
  file.copy(file.path(root, 'templates','TME','run_analysis.R'),runner)
  log <- file.path(work, 'validation.log')
  Sys.setenv(RNASEQ_LIB_DIR=file.path(root,'RNAseq_lib'))
  status <- system2(file.path(R.home('bin'),'Rscript'), c(shQuote(runner),shQuote(config)),stdout=log,stderr=log)
  if(status!=0L) stop('CLI failed: ',log)
  out <- file.path(settings$OUTDIR,'4-TME')
  native <- read.csv(file.path(out,'CIBERSORT_native_results.csv'),check.names=FALSE)
  coverage <- read.csv(file.path(out,'CIBERSORT_native_reference_coverage.csv'))
  scores <- as.matrix(read.csv(file.path(out,'ssGSEA_immune_scores.csv'),row.names=1,check.names=FALSE))
  prepared <- read.csv(file.path(out,'TPM_human_symbols.csv'),row.names=1,check.names=FALSE)
  expected <- run_tme_ssgsea(prepared)$scores
  manifest <- read.csv(file.path(settings$OUTDIR,'run_manifest.csv'))
  inputs <- read.csv(file.path(settings$OUTDIR,'run_inputs.csv'))
  stopifnot(nrow(native)==6,ncol(native)==ncol(sig)+4,all(is.na(native[['P-value']])),
    all(abs(rowSums(native[,colnames(sig)])-1)<1e-7),coverage$matched_genes==nrow(sig),
    all(is.finite(scores)),isTRUE(all.equal(scores,expected,tolerance=1e-10,check.attributes=FALSE)),
    manifest$inputs_complete==(species=='human'))
  cat(species,':',nrow(native),'samples,',ncol(sig),'cell types,',coverage$matched_genes,'reference genes,',nrow(scores),'ssGSEA signatures, complete inputs=',manifest$inputs_complete,'\n')
}
cat('VALIDATION_ROOT=',base,'\n',sep='')
