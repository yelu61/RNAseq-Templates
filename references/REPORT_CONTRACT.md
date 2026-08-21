# RNA-seq HTML Report Contract

`RNAseq_report.html` is a technical analysis report, not a directory viewer and
not an unconstrained model-written narrative. Its job is to make every analysis
domain discoverable, connect conclusions to saved evidence, and keep inference
within the limits of the experimental design.

## Coverage rule

Every analysis domain must have one explicit state in
`report_coverage_manifest.csv`:

- `covered`: represented in the main report, a compact summary, or appendix;
- `not_applicable`: incompatible with the data or biological question;
- `omitted_with_reason`: technically possible but intentionally excluded;
- `not_available`: no matching durable output was found and the gap remains
  unresolved.

The main narrative should show one decision-useful view per analysis question.
Redundant threshold, plot-type, and single-term variants belong in the output
index or a collapsible appendix rather than the primary reading path.

## Evidence rule

Every major section follows the same pattern:

1. result headline;
2. precise quantitative evidence;
3. interpretation;
4. assumption or limitation that changes the interpretation;
5. implication, falsifier, or next validation step.

Primary numerical claims must be read from saved tables during rendering or
independently recomputed by a validation script. Do not type headline values
manually into project reports.

## Claim levels

Project-specific narrative is stored in `claim_evidence_ledger.csv` using the
following levels:

| Level | Allowed language |
| --- | --- |
| `E0_observed` | observed pattern in a visualization or object |
| `E1_association` | replicated sample-level difference or association |
| `E2_directional` | directional support from time, dose, or stronger design |
| `E3_causal_model` | causal contribution in the tested perturbation system |
| `E4_mechanism` | convergent functional and orthogonal mechanism evidence |
| `E5_translational` | independently validated clinical/translational relevance |
| `unsupported` | proposed but not supported by current evidence |
| `superseded` | excluded after a contract or provenance failure |

Required ledger columns are:

`claim_id`, `claim`, `claim_level`, `scope`, `evidence`, `source_files`,
`assumptions`, `alternatives`, `falsifier`, and `next_evidence`.

An LLM may help draft or reorganize this ledger, but it must not invent source
files, numerical results, literature support, or evidence levels. A deterministic
validator checks the ledger before rendering.

## Annotation rule

Reviewed, project-specific prose lives in durable inputs that survive
re-rendering, never in the shared template and never by editing the rendered
HTML:

- `report_interpretation.md` — one whole-report narrative, rendered near the
  top;
- `report_interpretation/<section>.md` — per-section annotations rendered in
  place at the end of the matching section. Canonical section keys come from
  `report_interpretation_sections()`; `scaffold_report_interpretation()`
  creates starter files. Whole-line HTML comments are treated as guidance and
  never render, and a file containing only comments is skipped, so untouched
  starters stay invisible.

Annotations are subordinate to the statistical evidence and to the
claim-evidence ledger: they may interpret, contextualize, and flag caveats,
but must not introduce numerical claims that are not read from saved outputs,
and must not raise a conclusion above its ledger evidence level. Editing an
annotation plus re-rendering is the supported modification loop; the rendered
HTML remains a regenerable artifact.

Use `Rscript tools/render_report.R <run_dir> --scaffold` to create durable
annotation files and a review checklist, then re-render after each Markdown
edit without rerunning the statistical pipeline. Do not edit the shared
`reports/analysis_report.Rmd`/`.qmd` or the generated HTML inside a project.

## Statistical rule

- Define the experimental unit, model, contrast, primary threshold and
  multiplicity family before results.
- Report effect estimates and uncertainty alongside P/FDR values.
- Treat additional DEG thresholds as sensitivity analyses.
- A threshold tier using unadjusted (nominal) P values is exploratory
  hypothesis-generation only: it exists to keep ORA/gene-list analyses alive
  in low-DEG projects, must be labelled as nominal, must never serve as the
  primary threshold, and anything derived from it requires confirmation at an
  adjusted-P tier or by rank-based (GSEA/GSVA) evidence before supporting a
  claim above E0/E1.
- Separate database enrichment from custom or model-defined gene sets.
- With very small groups, supplement parametric gene-set comparisons with
  individual samples, standardized effect size, and exact label-permutation
  sensitivity. The coarse permutation resolution must be reported.
- “Not significant” is not equivalence; transcript scores are not functional
  activity; enrichment terms are not diagnoses or mechanisms.

## Gene-set rule

Custom gene sets should use a versioned registry containing at least:

`gene_set`, `theme`, `version`, `source_type`, `source_id`, and `gene`.

Recommended additional columns are `module_role`, `expected_direction`,
`biological_question`, `source_url`, `curation_note`, and `frozen_at`.
Previously executed versions are immutable. A revised set is a new version and
must be compared with the earlier result; a set revised after viewing outcomes
remains post hoc until tested in independent data.

## Validation gate

A shareable report must pass:

1. required analysis-domain coverage;
2. claim-ledger schema and provenance checks when narrative claims are present;
3. independent recomputation of headline values;
4. PDF structure and first-page rendering checks for referenced figures;
5. final HTML checks for required sections, embedded assets, overflow and
   browser errors;
6. a complete result inventory and reproducibility metadata.

The deterministic checks are recorded in `report_validation.csv`. In addition,
headline DEG counts are independently reconciled between
`Analysis_summary.txt` and `DEG_threshold_summary.csv`, and every project-relative
file named by the claim ledger must exist.
`report_review_checklist.csv` requires human sign-off for design and
experimental unit, sample QC/exclusions, contrast direction, primary FDR and
effect-size choices, ORA/GSEA/GSVA interpretation, claim calibration,
alternative explanations, limitations and figure integrity. Draft rendering
allows pending rows; `tools/render_report.R <run_dir> --publish` does not.
