# Publication Figure Standard

This is the shared visual contract for notebooks, command-line workflows, and
project-specific re-plotting. It targets figures that remain legible at their
**final printed size**, not only when enlarged in a notebook.

## Final-size geometry

`RNAseq_FIGURE_SPEC` and `publication_dimensions()` define conservative
defaults: 89 mm single column, 183 mm double column, and 247 mm maximum page
height. A journal's current author instructions always take precedence.

| Role | Recommended size | Use |
| --- | --- | --- |
| Single panel | 89 mm wide | simple two-group plot, compact schematic |
| Main multi-panel | 183 mm wide | PCA, volcano, heatmap, pathway comparison |
| Tall supplementary | 183 mm wide, ≤247 mm high | long term lists or many facets |
| Analysis archive | may exceed final size | exhaustive diagnostics; do not submit unchanged |

Do not shrink a 12–24 inch diagnostic canvas into a journal column. Reduce the
number of terms/facets or split it into panels first.

## Typography and line hierarchy

- `theme_publication()` uses 8 pt base text, 6.5–7.5 pt secondary text, and a
  9 pt left-aligned title at final size.
- `resolve_publication_font()` defaults to the portable `sans` alias
  (Helvetica-compatible on standard devices); a verified installed family can
  be requested explicitly. A figure must not mix serif and sans fonts.
- Axis and tick lines are 0.35 pt; decorative frames and heavy grid lines are
  removed. Titles explain the panel; methods and conclusions belong in the
  caption.
- Panel labels (`a`, `b`, …) are added during figure assembly, 9 pt bold, at a
  shared top-left anchor—not placed separately by every plotting helper.

## Labels and annotations

- Prefer direct labels only when they shorten the reader's eye movement.
- Wrap pathway terms to two lines and retain the complete term in the result
  table. Do not reduce text below 6 pt to fit more terms.
- Statistical brackets receive explicit upper-axis expansion. Use compact
  significance stars in dense panels and put exact adjusted P values/effect
  sizes in the table or caption.
- Volcano labels are a small, preselected set; connectors must not cross the
  title or leave the plotting panel. Bar labels require expansion or an inside
  position with adequate contrast.
- For heatmaps, hide unreadable row names rather than rendering illegible text.

## Colour

`RNAseq_PALETTE` is the single source for directional colours. Up/activated is
red-orange (`#d6604d`), down/suppressed is blue (`#4393c3`), and neutral is grey
(`#999999`). Heatmaps use a perceptually balanced blue–white–red diverging
scale centred at zero. Group colours must remain stable across all panels.

Never encode the only distinction with red versus green. Every important
colour distinction also needs position, sign, label, or shape.

## Export contract

- `save_pdf_plot()` and `save_pdf_device()` write to a temporary sibling file,
  reject empty/incomplete output, and only then promote the canonical PDF.
- `save_plot_bundle()` is for selected manuscript-facing ggplots. It exports
  PDF/SVG vector masters and optional 600-dpi TIFF; PNG is a review/report
  preview, not the primary master.
- Fonts are resolved once and applied at export. White is the default
  background. Raster panels are produced at final dimensions, not resized
  screenshots.
- Report previews belong in a derived `report_assets/` cache. The editable PDF
  or SVG under the canonical figure directory remains the source of truth.

## Required QA before delivery

Inspect representative figures at final size and verify: no clipped labels,
no overlap, readable smallest text, consistent font and group colours, correct
legend ordering, visible points/error bars, sensible whitespace, and no blank
pages. Also inspect dense edge cases (long terms, many samples, many facets).

Run the automated tests and notebook validation after any shared style change.
Automated checks detect syntax, dimensions, and corrupt exports; they do not
replace visual inspection.

## Enrichment bar labels

For single-sided ORA or GSEA NES plots, place the wrapped term **on its bar**
near zero. When both directions are present, place each term **across zero from
its bar**: a rightward bar has its term on the left, and a leftward bar has its
term on the right. Choose the layout from the directions actually plotted,
not the function name, and do not duplicate terms on the y axis.
Preserve measured bar lengths. A long label on a very short single-sided bar
may extend into that row's whitespace rather than distort the statistic or
shrink text until unreadable. GSEA FDR annotations stay in a separate outer
column on the bar's side, independently of term alignment.
Check both directions, one-sided results, short bars and long/duplicate labels.
