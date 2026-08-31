# v0.11.1 maintenance validation — 2026-08-31

This is post-fix evidence, separate from the [v0.11.0 baseline audit](../2026-08-31/README.md).
The current acceptance decision and limitations are in [FREEZE_AUDIT.md](../../FREEZE_AUDIT.md).

Tests ran in an isolated source copy. A temporary `.git` root provided resource
location only; it had no release revision. [SOURCE_SHA256SUMS.txt](SOURCE_SHA256SUMS.txt)
identifies the tested backend, notebook, report and test sources. A SHA file
identifies bytes, not an independently restored environment. Existing project
and repository demo outputs were not overwritten.

| Evidence | Result |
|---|---|
| [Full unit suite](unit-tests-final.log) | 1,021 pass, 0 fail, 0 skip; 3 UpSetR/ggplot2 deprecation warnings |
| [Notebook/library parsing](notebook-syntax.log) / [Final R parsing](final-r-syntax.log) | 10 notebooks, 17 helper files; 85 final repository R files, including verification scripts |
| [Full smoke](smoke.log) | General helper and six production runners pass; per-demo logs retained alongside this file |
| [Native human/mouse TME](native-tme-cli.log) | Both complete CLI runs pass; 6 samples each, matching 547/511 reference genes, 22/25 cell types, 27/26 true ssGSEA signatures |
| [General report render](report-draft.log) | rmarkdown fallback succeeds with embedded previews |
| [Unsigned publish gate](report-unsigned-publish.log) / [validation rows](report-validation.log) | Expected exit 1: unresolved coverage and scientific review; no approval was manufactured |
| [Result QA](output-qa.log) | 1,513 GO GSEA IDs, overlap limits and whole-family BH values verified; 0 unusable / 0 significant terms |
| [PDF rendering QA](pdf-qa.log) | 20 PDFs render with nonblank first pages; PCA, volcano, DEG heatmap and ORA previews inspected |
| [Bar-label regressions](bar-label-tests.log) | 158 plot assertions pass; 7 synthetic ordinary/one-sided/two-sided ORA/GSEA PDFs visually checked |
| [Analysis environment](analysis-session.log) | R 4.6.1 / macOS aarch64 and loaded package versions |
| [Project-memory audit](project-structure.json) | 0 errors, 15 warnings for absent research-project files; software-repository equivalents accepted |

KEGG REST access failed and was explicitly skipped in these smoke runs. This
is not validation of KEGG network access, public-cohort downloads, IOBR reference
downloads, Quarto, all operating systems, or old DOSE/fgsea. The local GO and
offline KEGG-style regression test the GSEA overlap/quality contract. Synthetic
TME mixtures test routing and calculation, not biological validity in a cohort.

To repeat with the recorded environment, from a fresh disposable backend copy:

```bash
Rscript tests/testthat.R
Rscript tests/validate_notebooks.R
Rscript examples/run_demo_smoke_test.R
Rscript references/validation/2026-08-31-maintenance/validate_tme_native_cli.R
Rscript tools/render_report.R examples/demo_RNAseq_General/RNAseq_General_Output --scaffold
Rscript references/validation/2026-08-31-maintenance/verify_release_outputs.R examples/demo_RNAseq_General/RNAseq_General_Output
python3 references/validation/2026-08-31-maintenance/verify_pdfs.py examples/demo_RNAseq_General/RNAseq_General_Output /tmp/rnaseq-pdf-qa
# This must fail until real coverage/review decisions are supplied:
Rscript tools/render_report.R examples/demo_RNAseq_General/RNAseq_General_Output --publish
```

The PDF checker needs Pillow and the Poppler `pdfinfo` / `pdftoppm` tools.
The demo script intentionally regenerates disposable demo outputs; never use it
inside a real run. Paths printed in retained logs refer to temporary validation
directories and need not remain available. Evidence checksums are recorded in
[SHA256SUMS.txt](SHA256SUMS.txt).

The [bar-label fixture](bar_label_previews.R) can regenerate the seven synthetic
layout examples from the backend root. It preserves actual bar lengths; a long
term on a very short bar uses adjacent same-row whitespace. These are layout
fixtures, not scientific results. Representative retained previews:

- [Bidirectional ORA](bar-label-bidir-both.png)
- [Positive/negative GSEA](bar-label-gsea-both.png)
