# Best Practices — end-to-end project lifecycle

This is the canonical, ordered path for running a **real** project with
RNAseq-Templates. It ties together steps that are documented in detail
elsewhere; follow the links for the full contract of each stage. For
first-time setup and the interactive notebook path, start with
[GETTING_STARTED.md](GETTING_STARTED.md) instead.

The lifecycle has one guiding rule: **every number and figure a reader sees is
regenerated from a saved run output, never hand-typed or hand-copied.** Each
stage below exists to protect that rule.

```text
configure → run → register → technical report → publish gate
        → narrative report → deck → curate into results/
```

## 0 · Route and configure

Let `bulk-rnaseq-analysis` inspect the data scale, design, and research
question, then pick the workflow ([references/TEMPLATE_SELECTION.md](references/TEMPLATE_SELECTION.md)).
For routine local bulk work, prefer the **General CLI runner** over a notebook.
Edit a project `config.R` — input paths, groups, contrasts, `THRESHOLD_GRID`
(including an optional exploratory nominal-p tier for low-DEG projects), and
any batch/paired design.

## 1 · Reproducible run

Run inside a fresh, versioned native output root so nothing lands in the
template source tree:

```bash
cd /path/to/project/analysis/runs/<run_id>
Rscript run_analysis.R        # copied from templates/General/
```

This produces the complete run bundle (`1-DEG/ 2-GSEA/ 3-Visualization/ 4-TME/
… sessionInfo.txt`). Never overwrite a completed run — a changed question or
config gets a **new `run_id`**. Notebooks are for exploration only and must not
use `notebooks/` as an output root.

## 2 · Register the run

Every General run writes `run_manifest.csv` (input/config checksums, analysis
and backend-code signatures, git revision/dirty state, lifecycle role, parent,
retention class). After one or more runs, rebuild the project registry:

```bash
Rscript tools/build_run_registry.R /path/to/project/analysis/runs
```

`RUN_REGISTRY.csv` marks exact duplicates via `duplicate_of` and recommends
review — it never deletes anything. See
[references/PROJECT_OUTPUT_LAYOUT.md](references/PROJECT_OUTPUT_LAYOUT.md) for
the lifecycle/retention model.

## 3 · Technical report + annotations

With `GENERATE_HTML_REPORT <- TRUE` the runner renders the generic technical
`RNAseq_report.html` (complete, archival, machine-checked against
[references/REPORT_CONTRACT.md](references/REPORT_CONTRACT.md)). Add durable,
per-section project interpretation by setting
`SCAFFOLD_REPORT_INTERPRETATION <- TRUE` once, then editing
`report_interpretation/<section>.md`. Re-render without re-analysis:

```bash
Rscript tools/render_report.R analysis/runs/<run_id> --scaffold
```

## 4 · Scientific publish gate

Rendering a draft is automatic; scientific approval is not. Complete
`report_review_checklist.csv` (design, QC, contrast direction, FDR/effect-size,
enrichment interpretation, evidence calibration, limitations, figure
integrity), then:

```bash
Rscript tools/render_report.R analysis/runs/<run_id> --publish
```

`--publish` fails unless coverage decisions are resolved and every checklist
item is approved or explicitly N/A. Validation lands in `report_validation.csv`.

## 5 · Narrative report (recommended for collaboration)

The technical report is archival; it is not the document you discuss with
collaborators. For that, write the project-specific, conclusion-led narrative
report on top of the shared machinery:

```r
source(file.path(Sys.getenv("RNASEQ_LIB_DIR"), "narrative_utils.R"))
scaffold_narrative_report(project_root)   # reports/narrative_report.Rmd
# write the prose; cite every number via sv/gl/gp/gnes/ora_terms
rmarkdown::render(...)
validate_narrative_report(html, min_figures = 5, required_sections = c(...))
```

Headline numbers are pulled from the run tables at render time, and the QA
gate must pass before curation. Full workflow and writing-style rules:
[references/NARRATIVE_REPORT.md](references/NARRATIVE_REPORT.md).

## 6 · Presentation deck (optional)

For discussion meetings, an [open-slide](https://github.com/1weiho/open-slide)
deck gives a comment-based revision loop. Its PNG figures are a **derived
cache** re-converted from the same run-bundle PDFs — never a second editable
copy. The single source of truth for every figure is the run-bundle PDF under
`3-Visualization/` / `4-TME/`.

## 7 · Curate into `results/`

Keep two layers: the complete native run bundle under `analysis/runs/`, and a
curated `results/{tables,figures,reports}/` holding only reviewed deliverables.
Every curated artifact records its source run and source file (a `SOURCE.txt`
or manifest). Full ownership/duplication rules:
[references/PROJECT_OUTPUT_LAYOUT.md](references/PROJECT_OUTPUT_LAYOUT.md).

## Definition of done

A project reporting cycle is complete when **all** of these hold:

- [ ] The run is a fresh `analysis/runs/<run_id>/` bundle with a `run_manifest.csv`, and `RUN_REGISTRY.csv` is rebuilt.
- [ ] The technical report validates and, if published, the review checklist is signed (`--publish` succeeded).
- [ ] Every headline number in every document is citer-generated or contract-validated — none hand-typed.
- [ ] Every figure traces to a run-bundle PDF; no second editable image copy exists.
- [ ] Exploratory nominal-p results are framed as hypothesis-generating (E0/E1), never headline.
- [ ] Curated `results/` artifacts carry provenance back to their source run.

## Where each topic lives

| Topic | Doc |
|-------|-----|
| First-time setup / notebook path | [GETTING_STARTED.md](GETTING_STARTED.md) |
| Choosing a template | [references/TEMPLATE_SELECTION.md](references/TEMPLATE_SELECTION.md) |
| All parameters | [references/PARAMETER_REFERENCE.md](references/PARAMETER_REFERENCE.md) |
| Two-layer output & run retention | [references/PROJECT_OUTPUT_LAYOUT.md](references/PROJECT_OUTPUT_LAYOUT.md) |
| Technical-report contract & evidence levels | [references/REPORT_CONTRACT.md](references/REPORT_CONTRACT.md) |
| Narrative-report workflow & style | [references/NARRATIVE_REPORT.md](references/NARRATIVE_REPORT.md) |
| Figure style / publication standard | [references/VISUALIZATION_STYLE_GUIDE.md](references/VISUALIZATION_STYLE_GUIDE.md) |
| Errors | [references/TROUBLESHOOTING.md](references/TROUBLESHOOTING.md) |
