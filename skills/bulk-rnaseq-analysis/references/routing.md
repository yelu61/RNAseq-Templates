# Backend Routing

## Decision matrix

| Data or primary task | Backend | Reason |
| --- | --- | --- |
| Local raw-count matrix, standard groups | RNAseq-Templates | General DESeq2 workflow and shared plotting/enrichment helpers |
| Local counts with complex contrasts or batch covariates | RNAseq-Templates | limma-voom or design-aware General workflow |
| Local longitudinal/time-series matrix | RNAseq-Templates | TimeCourse template |
| Local immune/stromal deconvolution | RNAseq-Templates | TME template with TPM contract |
| Local co-expression network | RNAseq-Templates | WGCNA template |
| Ordinary GEO expression study | RNAseq-Templates | GEO ingestion plus generic downstream workflows |
| TCGA/TARGET cohort expression or clinical data | TCGA toolkit | Cohort-aware preparation, clinical fields, and task runner |
| TCGA vs matched GTEx | TCGA toolkit | Curated tissue matching and targeted GTEx extraction |
| Pan-cancer expression or survival | TCGA toolkit | Multi-project task support |
| MAF, mutation survival, TMB | TCGA toolkit | Mutation-aware tasks and outputs |
| CNV/GISTIC or methylation | TCGA toolkit | Cancer multi-omic tasks |
| TCGA prognostic model and external validation | TCGA toolkit | Model/export/validation pipeline |
| Local discovery followed by TCGA validation | Both | Signature or gene list is the explicit handoff |

## Routing priority

1. Prefer source-specific handling when the source is TCGA, TARGET, or GTEx.
2. Prefer the TCGA toolkit when the requested analysis is cancer-cohort
   specific even if the input was exported to a generic CSV.
3. Prefer RNAseq-Templates for ordinary local/GEO expression analysis.
4. Use both only when the study contains an explicit discovery and validation
   boundary.
5. Ask a focused question only when the choice changes the statistical model
   or required inputs and cannot be inferred from the data.

## Common mixed workflows

### Local discovery to TCGA validation

1. Run local DEG, WGCNA, or time-course analysis in RNAseq-Templates.
2. Export a versioned gene/signature TSV with `gene`, optional `weight`, and
   optional `direction`.
3. Validate the signature with `validate_sc_signatures`,
   `run_clinical_assoc`, `survival_map`, or a toolkit pipeline.

### TCGA discovery to external cohort

1. Build the cohort and model through the TCGA toolkit.
2. Export coefficients and preprocessing metadata.
3. Run `external_validate` against the named cohort.
4. Do not recompute feature scaling independently in a generic notebook.

### Shared reporting

Keep backend-native result directories intact. Build the final narrative from
their manifests and tables; do not move selected figures into an untracked
folder without provenance.

## Conflict handling

- If the user asks for "TCGA DEG" through a generic notebook, choose the TCGA
  toolkit unless they explicitly need the teaching notebook.
- If the user asks for TCGA and GEO in one validation study, use the TCGA
  toolkit for TCGA/model logic and use RNAseq-Templates only for generic GEO
  preprocessing that the external-validation task does not cover.
- If the user provides TPM but requests DESeq2, stop and request raw counts or
  choose a statistically compatible method with explicit approval.
- If TME input is VST/rlog, request TPM or raw counts plus gene lengths.
