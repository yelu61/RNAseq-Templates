# Changelog

All notable changes to RNAseq-Templates are documented in this file.

## [Unreleased]

The `v0.11.x` line is feature-frozen. Future entries are limited to maintenance
fixes unless a new development line is opened explicitly.

## [0.11.1] - 2026-08-31

### Fixed

- WGCNA eigengenes now use module-color names; hub exports contain one true
  gene/eigengene correlation per selected gene. Module-trait P values use
  pairwise nonmissing sample counts, and empty selections retain table schemas.
- TME computes TPM over the original feature universe before ID/ortholog
  mapping, preserves native mouse symbols for the mouse signature, and shares
  preparation/coverage exports between General and TME entry points.
- Immune ssGSEA uses `ssgseaParam()` rather than the GSVA algorithm. Native
  CIBERSORT handles one sample, reports unavailable permutation P values as NA,
  and compares against IOBR only with matching human reference/preprocessing.
- All six runners reject existing native outputs before writing config or
  provenance. Schema-v2 manifests hash declared auxiliary/reference files,
  record R/package runtime signatures, and exclude incomplete/old inventories
  from duplicate grouping. Uncaptured online references remain explicitly unknown.
- limma `MIN_SAMPLE_FRAC` now enforces the documented raw-count sample fraction
  in addition to edgeR filtering; use zero for the previous edgeR-only behavior.
- GSEA preserves tested term IDs and unavailable statistics, filters sets on
  actual ranked-input overlap on the enrichit backend, and exports quality
  counts. Automatic figures and report interpretation distinguish BH FDR
  significance from descriptive ranking.

- Standardized pathway term labels inside ordinary/bidirectional ORA and GSEA
  NES bar directions, keeping true bar lengths and separate FDR annotations.

### Documentation and validation

- Updated the lifecycle guide, root documentation, template READMEs, parameter
  and output contracts, and migration/freeze assessment. Retained the original
  `v0.11.0` audit separately from the post-fix evidence.
- Added correctness regressions for hub values, species routing, ssGSEA,
  single-sample CIBERSORT, input/runtime provenance, directory protection,
  filtering and GSEA integrity; see the freeze audit for executed results.
- No new analysis module, global dependency upgrade, Git tag or hosted release
  is implied. Historical affected analyses require a new run; environment
  restoration and project scientific review are still project responsibilities.

## [0.11.0] - 2026-08-24

### Changed

- Entered feature-frozen maintenance mode: only bug, security, compatibility,
  database-adaptation, test, and documentation fixes are accepted on `v0.11.x`.
- All six production runners now write the same immutable `run_manifest.csv`
  provenance contract and expose the same lifecycle fields.
- CI smoke coverage now directly executes the General production runner in
  addition to the helper-level General workflow and the five topic runners.
- Aligned runner-first usage, version, parameter, and output-contract
  documentation with the production implementation.

### Added

- **Two-notebook General workflow**: `RNAseq_General.ipynb` remains the complete,
  independently executable pipeline, while the new Jupytext-paired
  `RNAseq_General_RunReview.R/.ipynb` reads an immutable CLI `RUN_DIR`, never
  refits the core differential/enrichment model, and writes selected-gene plots
  or optional custom-GSVA derivatives to a separate `REVIEW_OUTDIR`. A static
  regression test rejects core-analysis calls in the review notebook, and the
  notebook was executed against the bundled General demo.
- **Run manifest and registry (`run_utils.R`, `tools/build_run_registry.R`)**:
  General runs now record input/config checksums, a statistical analysis
  signature, lifecycle role, parent, change note, retention class and size in
  `run_manifest.csv`. The project registry is rebuilt after runs finish, safely
  marks exact duplicates through `duplicate_of`, prefers a full canonical run
  as the family owner, and never prunes automatically.
- **Independent report rendering and scientific publish gate**:
  `tools/render_report.R` re-renders a completed run after project-owned
  Markdown edits without re-analysis. `report_review_checklist.csv` separates
  draft generation from human approval, while `validate_analysis_report()`
  writes `report_validation.csv` and checks HTML, coverage, claims,
  reproducibility, PDF integrity, headline DEG-count reconciliation and review
  status; claim-ledger source paths must resolve. `--publish` requires all coverage
  decisions and checklist items to be resolved and signed.

- **Narrative-report layer machinery (`narrative_utils.R`,
  `reports/narrative_report_scaffold.Rmd`, `references/NARRATIVE_REPORT.md`)**:
  the bespoke, conclusion-led narrative report (the RPE/CF-style Chinese
  document layered above the generic technical `RNAseq_report.html`) was
  previously re-authored from scratch per project with no shared plumbing and
  no quality gate. The repo now ships the machinery that guarantees its
  quality: `make_report_citers(run_dir)` builds inline-citation closures
  (`sv`/`gl`/`gp`/`gnes`/`ora_terms`) that read every headline number from
  saved run outputs at render time (hand-typed values remain contract-banned);
  `show_pdf()` embeds run-bundle PDF figures as the single source of truth;
  `validate_narrative_report()` is an automated QA gate that hard-fails on
  missing figures/sections and warns on failed lookups before a report may be
  curated; `scaffold_narrative_report()` copies the generic skeleton
  (conclusion-led headings, per-section bio-boxes, contract reminders in
  non-rendering comments) into a project. `NARRATIVE_REPORT.md` documents the
  two-layer model, the quality mechanisms, and the standard workflow.
  Backfilled from the reference implementations in 202607XXR_RPE and 202605CF.
- **Per-section project annotations for the HTML report
  (`report_utils.R`, `reports/analysis_report.qmd` + `.Rmd` twin)**:
  projects can now annotate individual report sections via
  `report_interpretation/<section>.md` files (canonical keys:
  `scope`, `qc`, `deg`, `ora`, `gsea`, `custom_genesets`,
  `cross_contrast`, `synthesis`, `limitations`). Annotations render in place
  at the end of the matching section in a distinct callout box, support full
  markdown and inline R, and persist across re-renders — so reviewed
  interpretation is edited as a durable input instead of being re-typed after
  every regeneration. New `report_interpretation_sections()` and
  `scaffold_report_interpretation()` helpers create starter files (never
  overwriting written prose); starter guidance lives in whole-line HTML
  comments that never render, and comment-only files are skipped so untouched
  scaffolds stay invisible. `REPORT_CONTRACT.md` gains an "Annotation rule"
  section keeping annotations subordinate to saved evidence and the
  claim-evidence ledger.
- **Exploratory nominal-P DEG tiers via per-row `p_column` in `THRESHOLD_GRID`
  (`deg_utils.R`, `enrichment_utils.R`, `templates/General/config.R`,
  report template)**: the threshold grid was previously interpreted through a
  single run-wide `DEG_PVALUE_COLUMN`, so a project either called all DEGs on
  adjusted P or all on nominal P. Each grid row can now carry its own
  `p_column` ("padj"/"pvalue", NA/absent = inherit the run-wide default), so a
  run can keep inferential padj tiers *and* append an exploratory nominal-p
  tier (e.g. `exploratory: 0.05, p_column = "pvalue"`) that keeps ORA alive in
  low-DEG projects. DEG marking, the threshold summary (`pvalue_column` +
  `definition` columns), the raw-vs-shrunken strategy summary and
  `run_threshold_ora()` all resolve the column per row via the new
  `threshold_p_column()` helper. Guardrails: config rejects invalid values and
  warns if the default threshold is nominal; the report renders a caution box
  next to the DEG-count table whenever a nominal-p tier is present; the
  report's primary-threshold picker and `REPORT_CONTRACT.md` keep nominal
  tiers out of the headline reading path (confirmation at an adjusted-P tier
  or by rank-based evidence required for claims above E0/E1).
- **`TME_ORTHOLOG_CACHE` runner option (`templates/General/`)**: the General
  runner's TME block now passes `ortholog_cache = TME_ORTHOLOG_CACHE` through to
  `prepare_tme_expression()`, exposing the offline-deterministic ortholog cache
  (see below) to config-driven production runs. `config.R` documents the option
  (default `NULL` = always online, the old behaviour); the runner defaults it
  defensively for configs written before the option existed, and the config
  snapshot records it. Backfilled from 202603AK112, where the cache had only
  existed in project-local notebooks.
- **rmarkdown fallback for the HTML report (`report_utils.R`,
  `reports/analysis_report.Rmd`)**: `render_analysis_report()` previously
  required the Quarto CLI for the default `.qmd` template, so on machines
  without Quarto `GENERATE_HTML_REPORT <- TRUE` silently produced no report —
  observed in all three real-template projects (202507LPJ, 202607XXR_RPE,
  202603AK112). The repo now ships an `.Rmd` twin of
  `reports/analysis_report.qmd` (identical body, rmarkdown YAML); when Quarto
  is unavailable but `rmarkdown` is, the helper automatically renders the twin
  instead of stopping. Regression tests pin the twin body to the `.qmd` (drift
  guard) and exercise the fallback render end-to-end when Quarto is absent.
- **GSVA per-set group-comparison statistics in the runner
  (`templates/General/run_analysis.R`)**: the runner previously wrote
  `GSVA_scores.csv` and boxplots but never *tested* the scores, so "does this
  gene set differ by group?" had no statistical answer — surfaced by the
  202504CF_EV Pyro-EV macrophage runs, where the M1/M2/cGAS-STING sets looked
  different on the boxplot but the group-mean impression did not survive a
  proper test. After writing `GSVA_scores.csv` the runner now calls
  `pathway_group_comparison()` (existing library helper: per-set t-test/wilcoxon
  + BH across `COMPARISONS`) and writes `2-GSEA/GSVA_group_comparison.csv`
  (pathway × comparison, p, BH p.adj, direction). Guarded with `tryCatch` so a
  stats failure never aborts the run.

### Fixed

- **Code-aware duplicate detection and mixed manifest schemas**: General run
  manifests now record the backend git revision, dirty state and a checksum of
  the runner plus shared R modules. `RUN_REGISTRY.csv` includes that code
  signature in the exact-duplicate key, so runs made by different template
  code are not falsely collapsed. Registry rebuilding now fills missing
  columns when old and current manifest schemas coexist.
- **Report asset and coverage-refresh validation**: report QA now rejects HTML
  containing missing/failed PDF-preview markers and requires embedded figure
  previews when a run contains PDFs. Coverage decisions can no longer overwrite
  refreshed evidence counts/files with stale values during re-rendering.
- **Explicit sample IDs in General `colData.csv`**: the runner and complete
  notebook now persist row names as a `sample` column instead of relying on row
  order. The review notebook remains backward-compatible with older runs that
  lack the column.
- **CLI self-location under iCloud paths containing spaces**: on affected macOS
  Rscript installations, `--file=` encodes spaces as `~+~`. All runner,
  re-visualization and report/registry entry points now restore the path before
  normalization; converter-generated runners inherit the same fix.

- **Setup-chunk source leaked into rmarkdown-rendered reports
  (`reports/analysis_report.qmd` + `.Rmd` twin)**: the setup chunk relied on
  `knitr::opts_chunk$set(echo = FALSE)`, which only takes effect for
  *subsequent* chunks, and the `.Rmd` twin's YAML has no global echo setting —
  so under the rmarkdown fallback (the path used on machines without the
  Quarto CLI) the entire helper-code block was dumped at the top of the
  report. The setup chunk now carries an explicit `echo=FALSE` chunk option,
  which both Quarto and rmarkdown honor.
- **`4-TME/` missing from run summaries (`templates/General/run_analysis.R`)**:
  with `RUN_TME <- TRUE` the runner wrote a full `4-TME/` deconvolution bundle
  but neither `Analysis_summary.txt` ("3. Output Files") nor the final
  `Outputs:` log line listed it — surfaced by the 202603AK112 three-model runs.
  Both now include `4-TME/` when `RUN_TME` is on.
- **Mouse TME aborted when any symbol was un-cached and the network was down
  (`tme_utils.R`)**: `convert_mouse_symbols_to_human()` queried the network
  whenever *any* input symbol was absent from `ortholog_cache`, and a failed
  lookup (e.g. Ensembl endpoint HTTP 404) aborted the whole run even when the
  cache already covered the mappable majority. Two hardenings, backfilled from
  the 202605CF Ta1 run (whose `gene_name` column carries Ensembl-ID placeholders
  for ~1,600 features with no MGI symbol): (1) rownames matching `^ENSMUSG` are
  recognised as non-symbol placeholders and skipped rather than sent to the
  `mgi_symbol` filter; (2) a failed online lookup now degrades to the cached
  mapping (with a warning) instead of stopping, and — critically — does *not*
  mark the un-fetched symbols as queried or write them to the cache, so a later
  online run still retries them. Regression tests cover the Ensembl-skip path,
  the cache-only degradation, and cache non-poisoning.
- **Misleading Ensembl-ID warning for mixed mouse symbols (`tme_utils.R`)**:
  a small number of Ensembl placeholders is now handled quietly by the
  mouse-to-human converter; the input warning is reserved for matrices whose
  rownames are predominantly Ensembl IDs. The General runner's default IOBR
  panel now avoids online-only xCell (`estimate/cibersort/epic/mcpcounter`);
  xCell remains an explicit opt-in when its network dependency is verified.

### Added (continued)

- **Offline-deterministic mouse→human ortholog conversion (`tme_utils.R`)**:
  `convert_mouse_symbols_to_human()` gains `fallback_hosts` and `ortholog_cache`,
  and `prepare_tme_expression()` forwards both. The pinned archive host
  (`dec2021.archive.ensembl.org`) avoids dataset drift but is intermittently
  unreachable (and `www.ensembl.org` retired its biomart endpoint with HTTP 404),
  which previously aborted any mouse TME run mid-pipeline. The converter now tries
  each host in turn before failing, and — when `ortholog_cache` (an `.rds` path) is
  set — caches the mapping so repeat runs are offline-deterministic: only genes
  absent from the cache hit the network, and newly mapped genes are written back.
  Genes queried but having no ortholog are recorded so they are not re-queried.
  Backfilled after a real project run failed when the archive host timed out;
  regression tests cover the offline cache path and the unmapped-but-queried case.

- **`collapse_by_symbol()` (`data_utils.R`)**: collapse a count/expression matrix
  to unique gene symbols, keeping the highest-total-signal row per symbol and
  optionally carrying a per-row annotation (e.g. gene length, biotype) through to
  the retained row. Generalizes `deduplicate_expression_by_symbol()` (which does
  not carry annotation) for use when preparing runner-ready count tables that also
  need a gene-length column for TPM/TME. Backfilled from real project use;
  regression tests cover max-count selection, annotation carrying, and NA/empty
  symbol dropping.

### Fixed

- **Blank grid-composite figures in non-interactive runs (`deg_utils.R`,
  `plot_utils.R`)**: `UpSetR::upset()` and `enrichplot::gseaplot2()` both return
  grid objects that are only auto-printed at the top level. Called inside a
  `save_pdf_device()`/`ggsave()` path from `Rscript`, they produced a blank
  single-page PDF (< 1500 bytes), which `.promote_validated_figure()` then
  refused — silently skipping every GSEA running-curve and, for the DEG-overlap
  UpSet, aborting the whole run before `Analysis_summary.txt`. This only
  surfaced on the first multi-comparison design (which actually reaches the
  UpSet branch). `plot_deg_upset_pdf()` now wraps `upset()` in `print()`, and
  `plot_gsea_term_figure_pdf()` renders its re-classed `gglist` via
  `save_pdf_device(..., draw = function() print(p))` instead of `ggsave()`
  (which has no `grid.draw` method for `gglist`). Both now emit valid PDFs.
- **DEG-overlap block no longer fatal (`templates/General/run_analysis.R`)**:
  the optional DEG-overlap UpSet/Jaccard figures are wrapped in `tryCatch`, so a
  single diagnostic-figure failure logs a message and skips instead of killing
  the completed DEG/GSEA/GSVA work and the run summary.
- **GSEA run-to-run jitter (`enrichment_utils.R`)**: `run_go_gsea()` and
  `run_kegg_gsea()` wrap `clusterProfiler::gseGO()`/`gseKEGG()`, which use
  adaptive fgsea permutation with no exposed seed argument (clusterProfiler
  4.20), so NES/pvalue/p.adjust jittered between otherwise-identical runs and
  two runs of the same config could look spuriously different. Both helpers now
  take a `seed` argument (default `123`) and call `set.seed(seed)` immediately
  before the permutation, making GSEA output bit-reproducible across the
  General, Limma_Voom and TCGA_GEO runners without any runner change. Pass
  `seed = NULL` to restore stochastic behaviour. New regression test asserts
  two calls with the default seed return identical result frames.

## [0.10.0] - 2026-08-13

### Added

- **Command-line runners for the remaining five templates**, completing the
  CLI coverage that previously existed only for General. Each ships the same
  four files as `templates/General/` — `config.R` (the only file you edit),
  `run_analysis.R` (headless pipeline), `visualize_results.R` (re-plot from
  saved results, no recompute), and `README.md` — and follows the same
  conventions: config-as-argument bootstrap, 4-step `RNAseq_lib` resolution,
  numbered output layout, `0-Config/analysis_config_used.R` snapshot, targeted
  intermediate `.Rdata` saves, `Analysis_summary.txt` + `sessionInfo.txt` tail,
  and an optional `render_analysis_report()` block (off by default for these
  topics since the shared `.qmd` is General-oriented):
  - `templates/Limma_Voom/` — full voom → contrasts → volcano/heatmap → ORA +
    GSEA → theme-dotheatmap/single-term pipeline. Species-aware OrgDb.
  - `templates/WGCNA/` — direct WGCNA calls (no `wgcna_utils`), soft-threshold
    auto-pick with max-R² fallback and a `SOFT_POWER` override; outputs under
    `5-WGCNA/`.
  - `templates/TME/` — ESTIMATE / IOBR / native-CIBERSORT / ssGSEA under
    `4-TME/` with parent+sub-switch `requireNamespace` degradation.
  - `templates/TimeCourse/` — Mfuzz soft clustering (seeded, `5-TimeCourse/`)
    plus time-point-vs-baseline DESeq2 (`1-DEG_Timepoint/`).
  - `templates/TCGA_GEO/` — GDC / GEO / local-file acquisition, Tumor-vs-Normal
    DESeq2, single-gene + clinical KM, uni/multivariate Cox (`6-Survival/`),
    ORA + GSEA. Local-file mode runs fully offline.
- **Notebook → script converter** `tools/notebook_to_runner.R`: parses a
  template notebook, extracts the parameter cell into `config.R`, flattens the
  remaining code cells into a `run_analysis.R` draft (with the standard
  bootstrap), and writes a `conversion_report.txt` flagging side effects,
  hard-coded OrgDb, display-only calls and OUTDIR routing to review. Sourceable
  functions plus a CLI (`Rscript tools/notebook_to_runner.R in.ipynb out_dir`).
- Unit tests for the converter (`tests/testthat/test-notebook_to_runner.R`) and
  new regression tests for the enrichment/TME fixes below.

### Fixed

- **`prepare_enrich_df()`** dropped `factor level is duplicated` errors when an
  enrichment result carries duplicate `Description` values (KEGG GSEA returns
  the same pathway name under different IDs): all term rows are now preserved,
  missing labels fall back to `ID`, and repeated display labels are
  disambiguated without changing the enrichment table. The same row-preserving
  label contract is used by dot/bar/bidirectional and NES plots.
- **Runner path contract**: all six CLI runners resolve a relative config path
  before loading it and preserve the invocation directory as the explicit run
  root. Relative inputs/outputs therefore no longer resolve accidentally inside
  the template source tree. Converter-generated runners use the same bootstrap.
- **Graphics-device containment**: plotting helpers return invisibly, limma
  `voom()` runs once inside the managed PDF device, and WGCNA base plots use the
  shared transactional export helpers. The smoke test now fails on any leaked
  `Rplots.pdf` instead of hiding the side effect via `.gitignore`.
- **WGCNA trait encoding**: renamed copies of the sample-ID column are excluded
  from module--trait correlations, invariant traits are dropped, and
  categorical contrasts receive readable labels. This prevents malformed or
  legacy metadata from producing one heatmap column per sample.
- **Mfuzz figure styling**: saturated membership-rainbow trend panels are
  replaced by low-opacity cluster trajectories with a weighted centroid,
  journal-scale sans-serif typography, and stable final-size layout.
- **Strict demo exits and offline TME smoke**: topic demos require runner status
  zero. The required TME smoke path uses native ESTIMATE + ssGSEA without IOBR;
  IOBR remains an explicit integration path requiring cached reference bundles
  or network access.
- **`plot_gsea_term_figures_from_df()`** no longer aborts the whole loop when a
  single term produces a degenerate figure on small data — each term is wrapped
  in `tryCatch` and skipped with a message.
- **`compare_native_iobr_cibersort()`** now strips IOBR's `_CIBERSORT` (and other
  method) column suffix before normalizing cell-type names; previously the
  suffix survived normalization as a trailing token so no cell types matched the
  native LM22 names ("No common cell types").
- **TME template ordering**: `group_df_for_plot` is built right after metadata
  loading, before the native-CIBERSORT section that references it (the notebook
  defined it later, so a top-to-bottom `RUN_CIBERSORT=TRUE` run failed).
- **TimeCourse species handling**: the Mfuzz cluster-ORA step no longer
  hard-codes `org.Hs.eg.db`; the runner selects the OrgDb from `SPECIES`.
- **TCGA_GEO median-split KM** is computed once (the notebook duplicated it in
  two sections).

### Changed

- **Skills housekeeping**: the unified `bulk-rnaseq-analysis` skill lives in the
  central `agent-ready-research-skills` collection; the dangling project-level
  `.agents/skills/` + `.claude/skills/` symlinks and the stale
  `tests/__pycache__/` artifact are removed, the empty `skills/` now holds a
  `README.md` pointing at the canonical location, and `GETTING_STARTED.md` no
  longer describes the skill as bundled with this repo.
- **Topic demos now drive the template runners** (`examples/demo_RNAseq_*/
  run_demo.R` call `templates/<Topic>/run_analysis.R` with a demo config and
  assert the numbered-layout outputs), so CI exercises the production runners
  rather than a duplicated inline pipeline.

### Added (earlier in this cycle)

- **Final-size publication figure contract**: shared 89/183 mm dimensions,
  portable sans-serif font resolution, 8 pt typography hierarchy, safe vector
  and 600-dpi raster bundles, and transactional export for ggplot plus
  grid/base devices.
- **Real-project output contract** informed by `202607XXR_RPE`: complete
  backend-native run bundles are separated from curated deliverables;
  notebooks are source artifacts rather than output roots; report previews are
  explicitly derived assets.

- **`pathway_utils.R`** — validated custom gene-set scoring and pathway-focused
  visualization, so projects no longer hand-roll GSVA/statistics/plots:
  - `score_gene_sets()` wraps GSVA (gsva/ssgsea/zscore/plage), requires at least
    five overlapping features by default, and attaches a gene-set overlap audit.
  - `pathway_group_comparison()` aligns named groups by sample ID and computes
    per-set statistics with one global BH pass.
  - `plot_pathway_delta_summary_pdf()` renders a clipping-safe bidirectional
    delta chart across pathways × comparisons.
  - `plot_keygenes_log2fc_heatmap_pdf()` preserves missing log2FC as `NA`,
    rejects duplicate gene-by-comparison rows, and supports readable row groups.
  - `melt_gene_expression()` / `plot_gene_expression_pdf()` produce grouped
    per-gene expression plots with validated sample alignment.
- Group box/violin helpers accept a precomputed `stat_table`, allowing figures
  to reuse the exact multiple-testing correction reported in result tables.
- **`bulk-rnaseq-analysis` Codex/Claude Code Skill**:
  - Adds one user-facing entry point for generic local/GEO matrix workflows
    and the independent TCGA/TARGET/GTEx toolkit.
  - Includes deterministic repository discovery and routing with explicit
    warnings for invalid contracts such as DESeq2 on normalized values or TME
    deconvolution from VST/rlog.
  - Keeps one canonical Skill source under `skills/`, exposed to Codex and
    Claude Code through project-level discovery links.
  - Defines backend routing, scientific guardrails, cross-backend handoffs,
    and a reproducible output/provenance contract without copying analysis
    implementation between repositories.
  - Adds Python unit coverage for generic, TCGA, mixed, warning, and ambiguous
    routing cases; CI runs it before the R workflow.
- **Command-line runner for the General template** under `templates/General/`:
  - `run_analysis.R` + `config.R` run the full General pipeline (all 15 sections) non-interactively via `Rscript`, driven by a single editable config file. `RNAseq_lib` is located via `RNASEQ_LIB_DIR`, the project directory, the repository root, or the parent directory. Intended for reproducible / batch / headless execution alongside the notebook.
  - **Optional TME deconvolution switch** (`RUN_TME`) in `run_analysis.R`: builds TPM from raw counts + gene lengths, then runs native ESTIMATE, IOBR (`estimate`/`cibersort`/`epic`/`xcell`), and ssGSEA immune signatures, writing to `4-TME/`. Sub-switches `RUN_TME_ESTIMATE` / `RUN_TME_IOBR` / `RUN_TME_SSGSEA`; IOBR/estimate absence degrades gracefully. Mouse input is converted to human orthologs via biomaRt.
  - `run_analysis.R` now caches `gseaResult` objects to `2-GSEA/gsea_results.rds` and saves `GROUP_LEVELS`/`COMPARISONS`/`SPECIES`/`colData` into `1-DEG/DEG_results.Rdata` so downstream figures can be regenerated without recompute.
  - **`visualize_results.R`**: targeted re-visualization from saved results with no DESeq2/ORA/GSEA recompute. Three independent sections — key-gene bar/SEM + heatmap (from `DEG_results.Rdata`), single-term gseaplot2 figures (from cached `gsea_results.rds`), and ORA theme dot-heatmaps (from saved ORA csvs).
  - `templates/General/README.md` documenting usage, batch execution, switches, and when to use the notebook instead.
- **Unit tests for `plot_utils.R`** (`tests/testthat/test-plot_utils.R`, 23 test blocks): palette/theme, label/ratio/z-score/SEM helpers, `pairwise_effect_table`, `prepare_enrich_df`, and theme matching.

### Fixed

- Standardized mixed Times/sans typography, oversized titles/canvases, clipped
  annotation headroom, and inconsistent label sizing across QC, DEG, pathway,
  enrichment, GSEA, TME, time-course, TCGA, and survival figures.
- ComplexHeatmap, pheatmap, Mfuzz, UpSet, voom, and KM figures now use the same
  validated temporary-file promotion path as ggplot exports, preventing failed
  devices from replacing canonical figures.

- **Empty/corrupt PDF promotion**: `save_pdf_plot()` now writes to a sibling
  temporary PDF, validates that the device produced non-trivial content, and
  only then promotes the file. A plotting error can no longer leave a
  zero-page or blank canonical PDF behind.
- **GSEA multi-term blank pages**: the overview suite no longer sends multiple
  term IDs through a single-term running-curve helper. Overview plots are
  dotplot/NES/ridgeplot; running curves are one explicitly selected term per
  file via `plot_gsea_term_figures_from_df()`.
- **Enrichment label readability**: ORA dot/bar/bidirectional plots and GSEA
  NES bars use external y-axis labels capped at two lines with ellipsis, rather
  than long text drawn inside bars. Saved ORA CSVs can be replotted without the
  original enrichment object.
- **GSEA output ownership**: the General notebook now writes overview, theme
  maps and selected running curves under separate `3-Visualization/GSEA/`
  subdirectories.
- **General notebook theme-enrichment execution**: restored a code cell that
  had been serialized as one comment containing literal `\\n` separators;
  notebook validation now detects this otherwise syntactically valid failure.
- **Theme dot-heatmap oversizing**: automatic PDF dimensions are capped at a
  readable/safe maximum with an explicit warning instead of failing late when
  `ggsave()` rejects dimensions over 50 inches.
- **Pathway helper integrity**: named sample groups are aligned by sample ID;
  malformed matrices, undersized gene-set overlap, duplicate heatmap keys, and
  ambiguous row-group mappings fail explicitly. Missing log2FC remains `NA`
  instead of being rendered as zero, and zero deltas are labelled `NO_CHANGE`.
- `plot_gsea_nes_barplot_pdf()` no longer fails with `factor level is duplicated` when KEGG GSEA results carry `NA` or duplicated `Description` values (the online KEGG map supplies IDs without names); missing labels now fall back to the term ID and are deduplicated.
- `build_multi_comparison_enrich_df()` no longer drops enrichment rows whose `ONTOLOGY` is `NA` when an `ontology_filter` is set — this previously blanked ORA theme dot-heatmaps built from csv-read results (e.g. in `visualize_results.R`).
- Native-ESTIMATE table parsing in `run_analysis.R` now selects real sample columns explicitly instead of positionally, so the extra `Description.1` column ESTIMATE writes is not mistaken for a sample.
- Removed a duplicated `tme_utils.R` section in `references/FUNCTION_CATALOG.md`.
- `Rplots.pdf` remains ignored as a defensive cleanup rule, while the smoke
  suite now treats creation of one as a graphics-device leak and fails.

### Added

- **Runnable demos for all remaining topic templates**, each with a `regenerate_demo_data.R` + `run_demo.R` pair under `examples/`:
  - `examples/demo_RNAseq_TME/`: raw counts + gene lengths → TPM → ESTIMATE / IOBR (estimate, cibersort, epic) / ssGSEA. Human symbols avoid online ortholog mapping; some IOBR reference datasets may still require network access or a populated local cache.
  - `examples/demo_RNAseq_TimeCourse/`: VST expression across 4 time points with injected temporal patterns → Mfuzz clustering, plus raw-count time-point-vs-baseline DESeq2.
  - `examples/demo_RNAseq_TCGA_GEO/`: local-file mode (no network) with a synthetic TCGA-like cohort (TCGA barcodes, raw counts + TPM + clinical) → Tumor-vs-Normal DESeq2, KM survival, clinical KM, Cox, ORA/GSEA.
- `examples/run_demo_smoke_test.R` now drives every topic demo after the General demo and fails if any one fails, so each push exercises all templates.
- **Unified HTML report**: `reports/analysis_report.qmd` + `RNAseq_lib/report_utils.R` (`render_analysis_report()`). Assembles the saved DEG/ORA/GSEA/QC CSV and PDF outputs into a single self-contained HTML document without re-running the analysis. `RNAseq_General.ipynb` gains a `GENERATE_HTML_REPORT` parameter and a report cell (Section 15). Renders with the quarto CLI when available, else falls back to `rmarkdown`.
- `immune_gene_sets` built-in 28-cell-type immune signature collection (Charoentong et al. 2017) in `RNAseq_lib/tme_utils.R`, used by the TME ssGSEA step.
- CI (`smoke-test.yml`) installs the extra packages the new demos need (Mfuzz, WGCNA, survival, survminer, estimate, corrplot, e1071, babelgene, rprojroot) and rmarkdown for the report.

### Fixed

- **`immune_gene_sets` was undefined** in `RNAseq_TME_Deconvolution_Template.ipynb` — the ssGSEA cell referenced a variable that was never created and is not exported by IOBR, which would error for every user. Now provided as a built-in constant in `tme_utils.R`.
- **Native ESTIMATE failed with "找不到对象 'common_genes'"** in the TME template: `requireNamespace("estimate")` does not resolve the package's lazy-data objects, but `filterCommonGenes()`/`estimateScore()` reference them directly. Both the notebook and the demo now call `utils::data("common_genes"/"SI_geneset", package = "estimate")` first.
- **`plot_km_by_group_pdf()` failed with "object of type 'symbol' is not subsettable"** under R ≥ 4.x: the survfit call captured the formula as a local symbol, which `ggsurvplot()` could not re-evaluate. The formula is now inlined into the call via `eval(substitute(...))`. This broke all clinical-variable KM plots in the TCGA-GEO template.
- TimeCourse time-point-vs-baseline DEG no longer includes a single-level `condition` column in the DESeq2 design (which errored with "design contains variables with all samples having the same value"); the condition covariate is only added when it actually varies.

### Added

- **New `RNAseq_lib/data_utils.R`** with reusable data loading, validation, and gene-conversion helpers:
  - `read_expression_matrix()`, `read_metadata()`, `validate_samples_match()`
  - `detect_expression_scale()` for heuristic classification of raw counts / TPM / log-scale / VST
  - `counts_to_tpm()`, `counts_to_fpkm()`, `extract_gene_lengths()`, `validate_count_matrix()`
  - `convert_gene_ids()`, `convert_expression_rownames()` supporting human/mouse Ensembl↔symbol and MGI→HGNC
- **Unit tests** for `data_utils.R` in `tests/testthat/test-data_utils.R`.
- `validate_expression_contract()` enforces declared TPM/log2(TPM+1)/VST input units and gene/sample-name integrity.
- `tests/validate_notebooks.R` parses every notebook code cell and every shared R script; CI runs it on push and pull requests.
- `RNAseq_TME_Deconvolution_Template.ipynb` now defaults to raw integer counts plus gene lengths and computes TPM internally.
- `RNAseq_limma_voom_Template.ipynb`: new `SPECIES` parameter (default `"human"`) drives species-aware org.db and KEGG organism code selection.
- `RNAseq_TCGA_GEO_Template.ipynb`: new `GDC_COUNTS_ASSAY` and `GDC_TPM_ASSAY` parameters to explicitly select assays instead of relying on hard-coded `assays_list[[4]]`.

### Changed

- `RNAseq_TME_Deconvolution_Template.ipynb`, `RNAseq_TimeCourse_Template.ipynb`, `RNAseq_WGCNA_Template.ipynb`, and `RNAseq_TCGA_GEO_Template.ipynb` now load expression and metadata via `data_utils.R` helpers and perform explicit sample-name validation.
- `RNAseq_limma_voom_Template.ipynb` uses `validate_count_matrix()` after `build_count_matrix()`.
- limma-voom batch adjustment is now fitted as a model covariate; `removeBatchEffect()` is reserved for visualization.
- `RNAseq_General.ipynb` sources `data_utils.R` and calls `validate_count_matrix()` after building the count matrix.

### Fixed

- TME no longer treats VST/rlog as invertible log2(TPM+1); VST is rejected for CIBERSORT/EPIC/ESTIMATE input.
- Restored the truncated TME visualization cell and its CIBERSORT/EPIC/xCell/ESTIMATE outputs.
- Fixed successful gene-ID conversion/deduplication in `convert_expression_rownames()` and rejected unsafe numeric first-column inference.
- GEO SeriesMatrix assays must pass raw-integer-count validation before DESeq2; normalized GEO data are directed to limma.
- TPM conversion now rejects zero-total-RPK samples, and metadata rejects unknown factor levels.

- `plot_tme_heatmap_pdf()` now orders samples by group, then clusters within each group, so heatmaps keep biological replicates together while preserving within-group structure.
- `melt_estimate_scores()` now recognizes IOBR's `_estimate`-suffixed score columns and strips the suffix for consistent plotting.
- `get_cibersort_category_map("human")` now matches both canonical spaced LM22 names and IOBR's underscore-separated column names, fixing empty/wrong broad-category aggregation for IOBR CIBERSORT.
- `plot_estimate_boxplot_pdf()` supports `ncol`, `save_individual`, and `individual_prefix` for a combined 1×4 layout plus per-score single plots.

- **Mouse TME deconvolution support** in `RNAseq_TME_Deconvolution_Template.ipynb` and `RNAseq_lib/tme_utils.R`:
  - New `SPECIES` parameter (`"human"` / `"mouse"`) in the TME notebook.
  - `convert_mouse_symbols_to_human()` uses `biomaRt::getLDS` to map MGI symbols to HGNC symbols.
  - `prepare_tme_expression()` unifies log reversal + optional mouse-to-human conversion.
  - `deduplicate_expression_by_symbol()` keeps the highest-mean duplicate after ortholog conversion.
  - `validate_tme_input()` checks numeric/non-negative values, rownames, and warns on Ensembl IDs.
  - All TME methods (native ESTIMATE, native CIBERSORT, IOBR `estimate`/`cibersort`/`epic`/`xcell`, and ssGSEA immune signatures) now route through the prepared human-symbol matrix.
- **TME visualization improvements**:
  - New `GROUP_COLORS` parameter for user-defined group colors; falls back to `make_group_colors()`.
  - `calc_tme_barplot_size()` and `calc_tme_boxplot_size()` auto-adjust figure size by sample/cell-type counts.
  - New `plot_tme_per_celltype_pdf()` generates one focused violin/box PDF per cell type for CIBERSORT and EPIC.
  - New `plot_tme_heatmap_pdf()` wrapper produces xCell and IOBR ESTIMATE heatmaps with consistent group annotation colors.
  - All boxplots (ESTIMATE, CIBERSORT, EPIC, ssGSEA) now apply `group_colors` consistently.
- Metadata loader now renames duplicated `sample` columns in `colData.csv` to avoid downstream subsetting errors.
- `tests/testthat/test-tme_utils.R`: unit tests for TME helpers.
- `install_dependencies.R`: added `biomaRt` for mouse-to-human ortholog lookup.
- `references/PARAMETER_REFERENCE.md`: documented TME parameters and per-method input requirements.
- `references/FUNCTION_CATALOG.md`: documented new TME helpers.

### Changed

- `notebooks/RNAseq_TME_Deconvolution_Template.ipynb`:
  - Replaced per-method `undo_log_expr()` calls with a single `prepare_tme_expression()` step.
  - Native ESTIMATE, IOBR, native CIBERSORT, and ssGSEA now all use the prepared `expr_tme` / `expr_for_ssgsea` matrix.
  - Visualization block now uses dynamic sizing, per-cell-type outputs, and IOBR ESTIMATE heatmap.

### Fixed

- TME template no longer silently assumes human gene symbols; mouse data is now explicitly converted before deconvolution.
- TME plots no longer ignore user color preferences; `GROUP_COLORS` is propagated throughout.
- `plot_tme_heatmap_pdf()` now orders samples by group (using the factor levels) and clusters within each group, so heatmaps keep biological replicates together while preserving within-group structure.
- `melt_estimate_scores()` now recognizes IOBR's `_estimate`-suffixed score columns and strips the suffix for consistent plotting.
- `get_cibersort_category_map("human")` now matches both canonical spaced LM22 names and IOBR's underscore-separated column names, fixing empty/wrong broad-category aggregation for IOBR CIBERSORT.
- `plot_estimate_boxplot_pdf()` supports `ncol`, `save_individual`, and `individual_prefix` for a combined 1×4 layout plus per-score single plots.

- `examples/run_demo_smoke_test.R`: automated smoke test that runs the General notebook core pipeline on demo data and validates outputs.
- `RNAseq_lib/batch_utils.R`: batch-effect PCA (`plot_pca_by_batch_pdf`), batch variance-explained summary (`summarize_pve_by_batch`), and PVE barplot (`plot_batch_pve_pdf`).
- `RNAseq_lib/design_utils.R`: paired/repeated-measures design helpers (`make_paired_col_data`, `build_paired_design_formula`, `validate_paired_design`).
- `RNAseq_lib/geo_utils.R`: GEO SeriesMatrix download and parsing helpers (`download_geo_series_matrix`, `parse_geo_series_matrix`, `prepare_geo_counts`).
- `RNAseq_lib/io_utils.R`: Excel export helpers (`write_deg_excel`, `write_all_genes_excel`).
- `RNAseq_lib/timecourse_utils.R`: time-point vs baseline DESeq2 helpers (`run_timepoint_vs_baseline_deseq2`, `summarize_timepoint_deg`, `write_timepoint_deg_results`, `plot_timepoint_deg_summary_pdf`).
- `RNAseq_lib/tcga_utils.R`: `build_id_map_from_se()` to derive ENSEMBL→symbol maps from `SummarizedExperiment::rowData`, and `run_clinical_km()` for clinical-variable Kaplan–Meier curves.
- `RNAseq_PALETTE` constant and visualization style guide at `references/VISUALIZATION_STYLE_GUIDE.md`.
- Reference documentation: `references/PARAMETER_REFERENCE.md`, `references/FUNCTION_CATALOG.md`, `references/TEMPLATE_SELECTION.md`, `references/TROUBLESHOOTING.md`.
- `RNAseq_lib/timecourse_utils.R`: restored orphaned code fragment as `write_mfuzz_cluster_table()`.
- `examples/demo_data/demo_counts.tsv`: regenerated with real mouse gene symbols so enrichment steps pass in the smoke test.
- `examples/run_demo_smoke_test.R`: made enrichment assertions data-dependent; now verifies that ORA summary exists rather than requiring a specific threshold-specific CSV.
- `.github/workflows/smoke-test.yml`: GitHub Actions workflow to run the smoke test and unit tests on push/PR.
- `tests/testthat/` and `tests/testthat.R`: initial unit-test skeleton covering deg, io, and enrichment helpers.
- `examples/regenerate_demo_data.R`: helper to regenerate demo count table with real mouse symbols.
- `examples/demo_RNAseq_limma_voom/` and `examples/demo_RNAseq_WGCNA/`: pre-configured runnable demos for the limma-voom and WGCNA templates.

### Changed

- `RNAseq_lib/plot_utils.R`: `plot_deg_summary_pdf()` now tolerates single-threshold summaries (no `Threshold` column).
- `notebooks/RNAseq_WGCNA_Template.ipynb`: force first column to character when `GENE_COLUMN` is NULL, and convert expression matrix to numeric matrix before WGCNA input, preventing failures with numeric-looking gene symbols.

- `RNAseq_lib/deg_utils.R`: fixed `extract_deseq2_results()` referencing `raw_df` before it was defined.
- `notebooks/RNAseq_limma_voom_Template.ipynb`: replaced broken `BATCH_COLUMN` logic with explicit `BATCH_VECTOR` parameter.
- `notebooks/RNAseq_General.ipynb`:
  - Added `BATCH_VECTOR`, `PAIR_ID`, and `EXPORT_EXCEL` parameters.
  - Sources `batch_utils.R` and `design_utils.R`.
  - Supports paired design when `PAIR_ID` is supplied.
  - Generates batch-effect diagnostics when `BATCH_VECTOR` is supplied.
  - Exports `1-DEG/DEG_results.xlsx` when `EXPORT_EXCEL = TRUE`.
- `notebooks/RNAseq_TimeCourse_Template.ipynb`:
  - Added raw-count and metadata parameters for Section 6.
  - Implemented runnable time-point vs baseline DEG with optional paired design.
- `notebooks/RNAseq_TCGA_GEO_Template.ipynb`:
  - Added `DOWNLOAD_FROM_GEO` / `GEO_ACCESSION` parameters.
  - `GENE_ID_MAP_FILE` now defaults to `NULL` and falls back to `build_id_map_from_se()`.
  - Added `CLINICAL_VARS_FOR_KM` parameter and clinical-variable KM output.
- `install_dependencies.R`: added `ashr`, `ggrepel`, `data.table`, `matrixStats`, `msigdbr`, `DOSE`, `BiocParallel`, and `GEOquery`.
- `RNAseq_lib/enrichment_utils.R` and `RNAseq_lib/plot_utils.R`: added friendly `message()` logs when enrichment or plots are skipped due to empty results.

### Fixed

- `extract_deseq2_results()` no longer fails with "object 'raw_df' not found".
- limma-voom batch correction no longer checks nonsensical `names(GROUPS)` / `Sys.getenv()` conditions.

## [0.9.0] - 2026-07-27

Reliability and CI hardening release. The headline: the GitHub Actions smoke
test went from never passing to fully green across all six templates, with the
unit-test suite expanded to every helper module and two real integration bugs
fixed in the General template.

### Fixed

- **Paired-design integration bug (General template)**: the paired branch in
  `templates/General/run_analysis.R` and `RNAseq_General.ipynb` built `colData`
  but never created `group`, so `filter_low_count_genes()` and the
  `save(..., group, ...)` call failed with `object 'group' not found` on any
  paired run. `group <- colData$condition` is now set in the paired branch.
- **Batch-correction loop not closed (General template)**: `BATCH_VECTOR` fed
  PCA/PVE diagnostics but was never written into `colData`, so a batch-aware
  `DESIGN_FORMULA` (e.g. `~ batch + condition`) could not resolve. New
  `add_batch_col()` in `design_utils.R` (paired with the existing
  `validate_batch_design()`) writes the batch covariate into `colData`; both the
  runner and the notebook call it.
- **CI dependency install (the reason CI never passed)**:
  - The hand-written package list in `smoke-test.yml` had drifted from
    `install_dependencies.R` and used wrong namespaces (`bioc::survival` is a
    CRAN package; `any::estimate` is R-Forge-only). The workflow now reuses
    `install_dependencies.R` as the single source of truth.
  - `install_dependencies.R` warned instead of failing on missing packages, so
    the install step reported success and the smoke test was the first thing to
    actually stop. It now `stop()`s with a clear, attributable message.
  - `install_dependencies.R` hard-coded the P3M *jammy* repo while
    `ubuntu-latest` is now *noble* (24.04), pulling a `libMagick++` soname
    (`.so.8`) that noble's apt packages do not provide — this broke
    `SpatialExperiment` → `GSVA` → `IOBR` at load time. The distro codename is
    now detected from `/etc/os-release`, and `libmagick++-dev` (plus other Bioc
    system libraries) is installed explicitly.

### Added

- **Unit tests for the 6 previously-untested modules** (`batch_utils`,
  `geo_utils`, `report_utils`, `survival_utils`, `tcga_utils`,
  `timecourse_utils`), giving all 14 `RNAseq_lib` modules direct coverage.
  Suite: 293 assertions, 0 fail/warn/skip.
- **Integration tests** (`tests/testthat/test-design-integration.R`)
  reproducing the paired-`group` and batch-`colData` failure modes above.
- **limma-voom and WGCNA demos added to the smoke test**; the smoke test now
  exercises all six templates. Their demo input data is committed (with targeted
  `.gitignore` exemptions), matching how `examples/demo_data` is handled.
- **CI diagnosability**: per-demo logs (`examples/demo_logs/`) echoed into the
  main log, and a `failure()` artifact upload of demo outputs, so future CI
  failures are self-describing.
- The bundled CIBERSORT mouse signature
  (`references/CIBERSORT/cibersort_mouse_22.csv`) is now tracked, so the native
  CIBERSORT test runs on CI instead of skipping.

## Earlier releases

- Initial public release with six notebook templates, `RNAseq_lib` helpers, demo data, and Chinese Getting Started guide.
- Added IOBR TME, TCGA/GEO, limma-voom, time-course, and survival analysis modules.
- Added publication-grade enrichment visualization (theme dot-heatmaps and single-term GSEA figures).
