# 2026-08-31 validation evidence

Historical audit of `v0.11.0`, commit
`9418686a9e110ae44dd00078ee3a63af4047351e`. Current interpretation and maintenance
status live in [FREEZE_AUDIT.md](../../FREEZE_AUDIT.md).
These logs describe one local environment; they are not a claim of current
remote CI status, all-platform compatibility or scientific sign-off.

| Evidence | Result |
|---|---|
| [Unit suite](audit-unit-tests-git-root.log) | 603 pass, 0 fail, 0 skip, 11 warnings |
| [Notebook/library syntax](audit-syntax.log) | 10 notebooks, 17 library files |
| [Tracked R syntax](audit-tracked-r-syntax.log) | 76 tracked R files |
| [Full smoke workflow](audit-smoke.log) | General helper and six production runners pass; unavailable KEGG outputs are skipped |
| [Draft report](audit-report-draft.log) | General complete-demo report renders with rmarkdown fallback |
| [Unsigned publish attempt](audit-report-publish.log) / [validation rows](audit-report-validation.log) | Expected rejection; no scientific sign-off was supplied |
| [Analysis session](audit-analysis-session.log) | Installed versions used for the General run |
| [Research-project scaffold audit](audit-project-structure.json) | 0 errors, 9 missing research-memory file warnings; software-repository equivalents accepted |
| [WGCNA reproduction](audit-wgcna-repro.R) / [log](audit-wgcna-repro.log) | 60 genes produce 3,600 hub rows rather than 60 |
| [Mouse TME reproduction](audit-tme-repro.R) / [log](audit-tme-repro.log) | Prepared human IDs intersect the selected mouse reference in only 3 rows |
| [Registry reproduction](audit-registry-repro.R) / [log](audit-registry-repro.log) | Changed auxiliary metadata and missing primary-input hashes are incorrectly grouped |
| [GSEA output QC](audit-gsea-qc.R) / [log](audit-gsea-qc.log) | 2,028 all-NA exported rows; 309 finite nominal-P-filtered rows; zero BH<0.05 rows |

Tests ran in a temporary source export with a temporary Git root for resource
discovery. The source commit above, not that temporary repository's revision,
identifies the tested implementation. Existing project outputs were untouched.
Log paths point to that temporary run and need not continue to exist.

The three defect reproductions intentionally assert the **old defect**; they
are archived evidence, not new regression tests for corrected code. After a
fix, add proper expected-correctness tests to `tests/testthat/` rather than
requiring these historical assertions to keep passing. To repeat against the
affected version with the recorded dependencies, from the repository root:

```bash
Rscript references/validation/2026-08-31/audit-wgcna-repro.R
Rscript references/validation/2026-08-31/audit-tme-repro.R .
Rscript references/validation/2026-08-31/audit-registry-repro.R .
Rscript references/validation/2026-08-31/audit-gsea-qc.R /path/to/completed/General/run
```

The GSEA check is read-only and accepts a completed General run; a different
dataset or dependency version need not reproduce these counts. The original
temporary full run was not copied into this repository. Retained evidence
checksums are listed in [SHA256SUMS.txt](SHA256SUMS.txt).
