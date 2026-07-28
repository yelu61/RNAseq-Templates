# Analysis Output Contract

## Required handoff

At completion, report:

1. Backend and repository revision used.
2. Input files, assay scale, species, and identifier type.
3. Sample inclusion/exclusion and final group sizes.
4. Statistical design, contrasts, thresholds, and correction methods.
5. Configuration or notebook path.
6. Result directory and primary tables/figures.
7. Dependency, data-quality, and interpretation warnings.
8. Reproduction command and session/toolkit version.

## Recommended manifest

When the backend does not already create equivalent metadata, write a small
machine-readable manifest beside the outputs:

```json
{
  "backend": "rnaseq-templates",
  "backend_revision": "<git-commit>",
  "analysis_type": ["deg", "enrichment"],
  "species": "human",
  "input_scale": "raw-counts",
  "config": "<relative-path>",
  "results": ["<relative-result-path>"],
  "warnings": []
}
```

Prefer relative paths inside the analysis project. Never include credentials,
tokens, private remote URLs, or unrelated absolute home-directory paths.

## Cross-backend handoffs

Use explicit, versioned files:

- Gene/signature TSV: `gene`, optional `weight`, optional `direction`.
- DEG table: gene identifier, effect estimate, test statistic, raw P value,
  adjusted P value, and comparison label.
- Risk model: feature coefficients plus preprocessing/scaling metadata.
- Expression matrix: genes by samples, with a separate sample metadata table.

Document identifier conversion and dropped features. Never infer that symbols,
Ensembl IDs, transcript IDs, and platform probe IDs are interchangeable.
