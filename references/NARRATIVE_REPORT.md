# Narrative report — project-specific reporting layer

RNAseq-Templates produces **two** report layers per project, and they are not
interchangeable:

| Layer | File | Origin | Purpose |
|-------|------|--------|---------|
| Technical (generic) | `RNAseq_report.html` | **Template default step** — `GENERATE_HTML_REPORT <- TRUE` in `config.R` → `render_analysis_report()` renders the shared `reports/analysis_report.qmd` / `.Rmd` twin | Complete, archival coverage of every analysis module; it becomes publish-ready only after the scientific-review checklist is signed |
| Narrative (bespoke) | `<PROJECT>_analysis_report.Rmd` → `.html` | **Project-authored** — *not* a template default step | Conclusion-led, plain-language (usually Chinese) interpretation for discussion with collaborators |

The narrative layer is deliberately *not* fully generated: a conclusion-led
reading requires project-specific scientific judgment that a shared template
cannot supply. What the template **does** guarantee is the machinery and the
quality bar, so that each project's prose is the only thing written from
scratch. That machinery lives in `RNAseq_lib/narrative_utils.R` and
`reports/narrative_report_scaffold.Rmd` (backfilled from the 202607XXR_RPE and
202605CF bespoke reports, the reference implementations).

## How quality is guaranteed

Four mechanisms, in increasing order of strength:

1. **No hand-typed numbers.** Every headline value (DEG counts, per-gene
   log2FC/padj, GSEA NES, top ORA terms) is read from the saved run outputs at
   render time through the citers from `make_report_citers(run_dir)`
   (`sv` / `gl` / `gp` / `gnes` / `ora_terms`). Re-rendering after a re-run
   automatically updates every number, so the report can never silently
   disagree with the tables it cites. This is the single most important rule.

2. **Figures have one source of truth.** Narrative reports embed figures only
   from the run-bundle PDFs via `show_pdf("3-Visualization/....pdf", base_dir = run_dir)`.
   The base64 images inside the HTML are ephemeral render-time conversions —
   there is no second, editable copy of any figure. (The open-slide deck's
   PNGs are likewise a derived cache, refreshed from the same PDFs.)

3. **The claim-evidence ledger still applies.** The narrative report inherits
   `references/REPORT_CONTRACT.md`: headline claims must be supported at the
   stated evidence level (E0–E5); nominal-p exploratory-tier results stay out
   of the headline reading path; visual QC separation is not mechanistic
   evidence; "not detected" is not "equivalent". The scaffold's section stubs
   prompt for these framings.

4. **An automated QA gate.** After rendering, run
   `validate_narrative_report(html, min_figures, required_sections)`. It hard-fails on:
   the file missing, fewer than the expected number of embedded figures, any
   missing-figure marker (`图不可用`, `PDF 预览失败`, …), or a required section
   heading absent. It warns on suspicious inline `NA` values (a failed citer
   lookup, e.g. `NES ≈ NA`) so a typo'd gene symbol or contrast name cannot
   slip through silently. A report must pass this gate before it is copied
   into `results/reports/`.

The technical report has a separate human gate in
`report_review_checklist.csv`. Rendering a draft is automatic; scientific
approval is not. `Rscript tools/render_report.R <run_dir> --publish` requires
all design, QC, contrast, multiplicity, enrichment, claim, limitation and
figure-integrity items to be approved (or explicitly not applicable).

## Standard workflow for a new project

```r
# 1. Generate the technical report (template default; config flag).
#    Nothing to do beyond GENERATE_HTML_REPORT <- TRUE.

# 2. Scaffold the narrative report into the project.
source(file.path(Sys.getenv("RNASEQ_LIB_DIR"), "narrative_utils.R"))
scaffold_narrative_report(project_root = "/path/to/PROJECT")

# 3. Edit reports/narrative_report.Rmd: set the title, params$run_id, and
#    write the project-specific prose. Keep the citers; delete N/A sections.

# 4. Render + gate (typically wrapped in scripts/R/09_render_narrative_report.R):
rmarkdown::render("reports/narrative_report.Rmd",
                  params = list(project_root = ".", run_id = "<run_id>"))
validate_narrative_report("reports/narrative_report.html",
                          min_figures = 5,
                          required_sections = c("技术摘要", "差异表达", "综合"))

# 5. Curate the passing HTML into results/reports/ with a SOURCE.txt entry.
```

## Writing-style rules (what made RPE/CF read well)

- **Conclusion-led headings.** `"Combo vs NIR 是结构性零结果"`, not
  `"差异表达结果"`. A reader skimming headings alone should get the story.
- **One bio-box per major section** (`::: {.bio-box}` with a
  `<span class="bio-tag">解读</span>` label): what the module means biologically,
  in plain language, plus its caveat.
- **Lead with the answer.** The first content block is a summary box with the
  3–5 sentence bottom line.
- **Hedge calibrated to evidence.** Distinguish "no signal detected" from
  "no effect exists"; flag exploratory/nominal-p and working-model statements.
- **Keep it Chinese** when that is the collaboration language; keep gene
  symbols, contrast names, and statistical terms in their original form.

## When to revise vs regenerate

The narrative Rmd is a durable, project-owned input. Edit prose and re-render
freely — numbers and figures re-pull from the current run outputs. Only the
prose is manual; everything quantitative is regenerated each render.
