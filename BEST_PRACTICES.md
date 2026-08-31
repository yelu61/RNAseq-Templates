# Best Practices — end-to-end project lifecycle

This is the canonical, ordered path for running a **real** project with
RNAseq-Templates. It ties together steps that are documented in detail
elsewhere; follow the links for the full contract of each stage. For
first-time setup and the interactive notebook path, start with
[GETTING_STARTED.md](GETTING_STARTED.md) instead.

**Release boundary:** `v0.11.1` closes the reproduced maintenance defects; the
`v0.11.x` feature scope stays frozen. Read the
[freeze audit](references/FREEZE_AUDIT.md) for execution evidence, migration
requirements and untested optional paths. A stable backend still needs a
reviewed project design and a tested, restorable project environment.

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

Before running, record the biological replicate, species, annotation version,
input scale, sample inclusion, primary contrast and primary FDR/effect-size
threshold in the project analysis plan. Match count columns to metadata by ID
and order; General uses the configured vectors, not an automatically loaded
metadata CSV. Explicitly set `COUNT_COLS` when the input contains other numeric
annotation columns. Keep raw files unchanged and use project-owned configs and
focused scripts; do not adapt the shared backend for each dataset.

General and limma-voom require unnormalized counts, not TPM/FPKM/VST. General's
`PAIR_ID` replaces the design with `~ pair_id + condition`; it does not preserve
an additional batch term. Review a project-specific model before using combined
pairing/batch, interactions, nested designs or other unsupported covariates.
Plot-level t tests are descriptive follow-up and do not replace the fitted DEG
model. These input/model distinctions follow the
[DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

Run only modules that answer the question. Keep optional TME/TF analyses off
until their input, species, reference data and dependencies have been checked.
The default installer does not include `dorothea`/`viper`; TF analysis needs
them separately. A low DEG count is not a reason to replace the predeclared
primary FDR threshold with nominal P values.

For enrichment, check valid gene-set identifiers, finite statistics and the
predeclared adjusted-P threshold. A non-empty CSV or top-term plot does not
establish significance. GSEA tables preserve tested rows and mark invalid
statistics, with a companion `_quality.csv` counting valid, unavailable and
FDR-significant terms. Automatic GSEA figures use significant terms; a manual
ranked selection remains exploratory. Report zero significant results explicitly.
Historical GSEA results affected by the dependency bug need a new run; see the
freeze audit.

Pin one backend commit for the project, and record any local patch. A shared
checkout that receives updates during a study is not a frozen dependency. Save
a restorable project environment (for example, a tested `renv` snapshot plus R,
system-tool and reference-data versions, or an archived container). The
repository does not ship a lockfile; `sessionInfo.txt` records versions but
cannot restore them. See the [renv guide and its limitations](https://rstudio.github.io/renv/articles/renv.html).

## 1 · Reproducible run

Run inside a fresh, versioned native output root so nothing lands in the
template source tree:

After copying `templates/General/config.R` to
`analysis/config/general_primary.R` and editing it, use this shell pattern.
Replace the two absolute paths and the run ID; keep the backend at the recorded
commit. Inputs in the config should be absolute, or relative to the new run
directory (not relative to the config file).

```bash
set -e
RNASEQ_BACKEND="/path/to/pinned/RNAseq-Templates"
RNASEQ_PROJECT="/path/to/project"
RNASEQ_RUN_ID="20260831_general_primary_01"
RNASEQ_RUN_DIR="$RNASEQ_PROJECT/analysis/runs/$RNASEQ_RUN_ID"
RNASEQ_CONFIG="$RNASEQ_PROJECT/analysis/config/general_primary.R"
export RNASEQ_LIB_DIR="$RNASEQ_BACKEND/RNAseq_lib"
test -f "$RNASEQ_CONFIG"
mkdir -p "$RNASEQ_PROJECT/analysis/runs"
mkdir "$RNASEQ_RUN_DIR"  # deliberately fails if this run already exists
cd "$RNASEQ_RUN_DIR"
Rscript "$RNASEQ_BACKEND/templates/General/run_analysis.R" "$RNASEQ_CONFIG" > run.log 2>&1
```

This produces the complete run bundle (`1-DEG/ 2-GSEA/ 3-Visualization/ 4-TME/
… sessionInfo.txt`; `4-TME/` exists only when requested). Never overwrite a
completed run — a changed question or config gets a **new `run_id`**. Runners
reject existing native output directories and provenance files before writing
a config snapshot. Keep the fresh-directory check above too: the guard is not
an operating-system lock against concurrent jobs or later manual edits.

Full notebooks are for exploration under `analysis/notebook_output/`; avoid
duplicating a complete CLI run with a complete notebook run. To inspect a
completed General run, use `RNAseq_General_RunReview.ipynb` with immutable
`RUN_DIR` and a separate `REVIEW_OUTDIR`. `visualize_results.R` also avoids core
refitting, but it changes to its own script directory: use a project-local copy
with explicit source and derivative-output paths, not the central script with
its defaults. Freeze final figures after this review stage.

Check the log, expected tables, non-empty figures and actual model before
accepting completion. An unexpectedly skipped or failed method required by the
analysis plan (including KEGG when required) prevents completion of that scope;
intentionally disabled modules are not required. On cloud-synced storage, verify files on disk; for heavy
runs prefer local scratch and copy the verified complete bundle back.

## 2 · Register the run

Every production runner writes `run_manifest.csv`, `run_inputs.csv` and
`run_runtime.csv` (declared input/reference checksums, R/package versions,
analysis and backend signatures, git revision/dirty state, lifecycle role,
parent and retention class). After one or more runs, rebuild the project registry:

```bash
Rscript "$RNASEQ_BACKEND/tools/build_run_registry.R" "$RNASEQ_PROJECT/analysis/runs"
```

`RUN_REGISTRY.csv` marks candidate duplicates via `duplicate_of` and recommends
review — it never deletes anything. Only complete schema-v2 inventories with
matching scientific inputs, analysis, code and runtime signatures are eligible.
Unknown online references are recorded as unknown and block duplicate grouping;
this includes KEGG in General/limma/TCGA-GEO enrichment runs and uncaptured IOBR
references. Older manifests cannot prove a duplicate. Keep a project manifest
for extra custom-script inputs not known to the runner, and never delete or
slim from `duplicate_of` alone. See
[references/PROJECT_OUTPUT_LAYOUT.md](references/PROJECT_OUTPUT_LAYOUT.md) for
the lifecycle/retention model.

## 3 · Technical report + annotations

For General runs, `GENERATE_HTML_REPORT <- TRUE` renders the generic technical
`RNAseq_report.html` (complete, archival, machine-checked against
[references/REPORT_CONTRACT.md](references/REPORT_CONTRACT.md)). Add durable,
per-section project interpretation by setting
`SCAFFOLD_REPORT_INTERPRETATION <- TRUE` once, then editing
`report_interpretation/<section>.md`. Re-render without re-analysis:

```bash
Rscript "$RNASEQ_BACKEND/tools/render_report.R" "$RNASEQ_RUN_DIR" --primary-threshold=standard --scaffold
```

## 4 · Scientific publish gate

Rendering a draft is automatic; scientific approval is not. Complete
`report_review_checklist.csv` (design, QC, contrast direction, FDR/effect-size,
enrichment interpretation, evidence calibration, limitations, figure
integrity), then:

```bash
Rscript "$RNASEQ_BACKEND/tools/render_report.R" "$RNASEQ_RUN_DIR" --primary-threshold=standard --publish
```

`--publish` fails unless coverage decisions are resolved and every checklist
item is approved or explicitly N/A. Validation lands in `report_validation.csv`.
Pass the actual `DEFAULT_THRESHOLD` to `--primary-threshold`; the tool defaults
to `standard` rather than inferring another primary tier. Here `--publish` is a
local validation gate, not a web deployment. Do not manufacture scientific
approval just to make the gate pass. Once reviewed, archive the report,
annotations, checklist and validation together; later revisions need a separate
version. The generated analysis tables remain unchanged during report editing.

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
- [ ] All scientific inputs and reference files have project-level checksums; the exact backend and a restorable environment are retained.
- [ ] The biological replicate, sample exclusions, actual fitted design, contrast direction, primary threshold and skipped/failed modules are reviewed.
- [ ] Relevant open issues in the freeze audit are repaired and verified, or the affected paths are explicitly excluded.
- [ ] For a General run, the technical report validates and, if published, the review checklist is signed (`--publish` succeeded); topic runs retain a reviewed summary until a topic-specific report contract exists.
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
