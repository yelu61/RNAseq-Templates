# TCGA Toolkit Backend

## Repository contract

Validate the TCGA root by checking for:

- `tcga_toolkit/VERSION`
- `tcga_toolkit/scripts/run_task.R`
- `tcga_toolkit/specs/README.md`
- `tcga_toolkit/references/task_index.md`

Treat `tcga_toolkit/` as the canonical implementation. Read
`tcga_toolkit/specs/README.md` for configuration schemas and the task index
only for the selected task.

## Task families

| Goal | Tasks |
| --- | --- |
| Inspect and prepare | `audit_data`, `prepare_bulk_rna`, `cohort_qc` |
| Expression and enrichment | `run_deg`, `run_enrichment`, `ssgsea_score` |
| Clinical and survival | `run_clinical_assoc`, `stage_analysis`, `subtype_analysis`, `survival_map` |
| TME and immune | `run_tme`, `immune_phenotype` |
| Cross-cohort | `run_gtex_compare`, `pan_cancer_expression`, `validate_sc_signatures` |
| Mutation | `maf_summary`, `mutation_survival`, `tmb_analysis` |
| Multi-omics | `cnv_summary`, `methylation_diff` |
| Networks and models | `gene_correlation_heatmap`, `wgcna_modules`, `prognostic_model` |
| Validation and treatment | `external_validate`, `drug_response` |
| Orchestration and reports | `pipeline`, `render_report` |

Use `Rscript tcga_toolkit/scripts/list_tasks.R --json` as the
machine-readable source of the current task list.

## Standard workflow

1. Inspect the project and dependencies:

   ```bash
   Rscript tcga_toolkit/scripts/inspect_project.R --project TCGA-XXX
   Rscript tcga_toolkit/scripts/check_deps.R --task <task>
   Rscript tcga_toolkit/scripts/list_runs.R --task <task>
   ```

2. Reuse or copy a JSON configuration from `tcga_toolkit/templates/`.
3. Validate without running:

   ```bash
   Rscript tcga_toolkit/scripts/validate_config.R --config <config.json>
   ```

4. Execute:

   ```bash
   Rscript tcga_toolkit/scripts/run_task.R \
     --config <config.json> \
     --output-root <analysis-project>
   ```

5. Read the native run bundle under `tcga_runs/<task_id>/`, including
   `report.md`, `run_metadata.json`, and results/plots/objects.

Use `--overwrite` only when the user explicitly wants to reuse an existing run
directory.

## Extension rule

If a requested cancer analysis is absent from the task list, follow
`tcga_toolkit/references/extension_guide.md` and add a reusable toolkit task.
Do not implement a one-off copy inside RNAseq-Templates or the Skill.

## Performance and data rules

- Inspect actual clinical columns before choosing group, subtype, or stage
  fields.
- Prefer targeted extraction from `GTEX/gtex_RSEM_gene_tpm.gz`; avoid loading
  the full multi-gigabyte GTEx matrix unless necessary.
- Run cohort QC before prognostic modelling.
- Record the toolkit version from `tcga_toolkit/VERSION` and
  `run_metadata.json`.
- Keep local restricted datasets out of the Skill and out of the
  RNAseq-Templates repository.
