# Real-project Output Layout

## Recommended production model

Use `bulk-rnaseq-analysis` to inspect the input scale, sample design, and
research question, then select the backend and create a project-specific
configuration. For routine local bulk RNA-seq, prefer the General command-line
runner as the reproducible execution entry point. Copy a notebook only when an
interactive, auditable exploration is genuinely useful.

Keep two explicit output layers:

```text
analysis/
  config/                         project-specific parameters and contrasts
  scripts/                        copied/adapted runner and focused analyses
  notebooks/                      source notebooks only; no bulk output
  runs/<run_id>/                  complete backend-native run bundle
    0-Config/ 1-DEG/ 2-GSEA/ 3-Visualization/ 4-TME/
    Analysis_summary.txt sessionInfo.txt run_manifest.csv
  runs/RUN_REGISTRY.csv           rebuilt index; duplicate/lifecycle overview
  notebook_output/<analysis>/     exploratory notebook outputs only
results/                          curated, canonical deliverables
  tables/ figures/ reports/
  report_assets/                  derived browser previews; rebuildable cache
```

All six production runners write the listed config snapshot, summary,
`sessionInfo.txt`, and `run_manifest.csv`. The native run bundle preserves completeness and reproducibility. `results/`
contains only reviewed tables, figures, and reports selected for scientific
communication. Every curated artifact records its source run and source file.

## Why `notebooks/` becomes messy

Relative paths such as `./1-DEG` or `./3-Visualization` resolve against the
process working directory, not the notebook file's conceptual role. If a
notebook is executed while the working directory is `notebooks/`, hundreds of
tables and figures are created beside the notebook. Adding a later curated
`results/` layer then appears to duplicate the analysis.

Avoid this by setting one explicit run root before execution, verifying it in
the first cell, and never treating `notebooks/` as an output root. A notebook
copied into a project should be source code; full exploratory output belongs in
`analysis/notebook_output/<analysis>/`, while production CLI bundles belong in
`analysis/runs/<run_id>/`. A run-review notebook reads that CLI bundle and writes
only its derivatives to a separate review directory.

## Ownership and duplication rules

- One canonical owner per artifact. Compatibility names should be links or
  manifest aliases, not byte-for-byte copies.
- Versioned gene-set results (`v1`, `v2`) are legitimate when their registry,
  checksum, and rationale are recorded.
- PNG report previews are derived assets; PDF/SVG files are figure masters.
- Never overwrite a completed run. Use a new `run_id` and record backend
  revision, configuration hash, input checksums, session information, and
  worktree state. Runners reject existing native artifacts before writing; this
  does not lock against concurrent jobs or later manual edits. Require a fresh
  directory and archive the finalized bundle after report/figure review.
- Classify each run as `candidate`, `canonical`, `sensitivity`, `repro_check`
  or `superseded`; record its parent and concise change note. Rebuild
  `RUN_REGISTRY.csv` from per-run manifests after batch execution. A candidate
  duplicate requires complete schema-v2 input/reference inventories and matching
  analysis, backend and runtime signatures. `run_inputs.csv` hashes declared
  auxiliary files; `run_runtime.csv` records R/packages. Unknown references,
  missing checksums and old schemas block grouping. This is not proof of
  identical results: custom inputs still need project-level provenance. Do not
  prune from `duplicate_of` alone. The manifest also records backend revision and dirty
  state. Within one duplicate family, a full canonical run is preferred as the
  retained owner;
  other rows point to it through `duplicate_of`.
- Retention is explicit and non-destructive: canonical runs stay `full`;
  sensitivity runs may become `slim`; superseded or exact reproducibility
  copies may become `metadata_only` only after curated provenance points to a
  retained owner. The template never prunes automatically.
- Separate exhaustive diagnostic plots from manuscript-facing primary figures.
  A figure manifest should record role, source table, dimensions, format, and
  QA status.

## Lesson from `202607XXR_RPE`

The project is scientifically coherent and its curated `results/` layer is
well organized. The apparent disorder comes from historical layering: the
General notebook first generated a complete native bundle under `notebooks/`;
project scripts later read that bundle and built a cleaner `results/` package.
Some root-level PDFs are compatibility copies, while `report_assets/` contains
intentional HTML previews. The problem is therefore output ownership and
naming—not a broken differential-analysis design.

For future projects, route standard execution through the CLI runner, reserve
notebooks for exploration, and declare the two output layers from day one.
