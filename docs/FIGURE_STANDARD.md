# Academic figure standard

Use traceable solver data and reproducible Matplotlib scripts. Preserve earlier
figures and generate a fresh versioned export directory for every revision.

- Single-column width about 88 mm or two-column width about 179 mm. Use a
  consistent serif/STIX font with mathematical notation, legible at final size.
- Label physical quantities and units; use nondimensional variables when useful.
  Use panel labels, restrained colourblind-friendly colours, thin lines and
  markers that remain distinguishable in monochrome.
- Put detailed captions outside the plotting area. Define symbols, conditions,
  mesh, normalization, iteration/time meaning and any omitted startup interval.
- Preserve relevant scales. Use logarithmic axes deliberately, avoid crowded
  minor-tick labels, and do not imply a physical time history from solver
  iterations. Common contour scales are required for quantitative field comparisons.
- Export vector PDF/SVG and high-resolution PNG. Embed PDF fonts; render PDFs
  for visual QA. Keep full native solver exports as the audit evidence.
- Archive source-data hashes, plotting code and captions. Never invent
  uncertainty/error bars or label sensitivity differences as statistical error.
- Distinguish calculations, numerical convergence, cross-code agreement,
  benchmark comparison and physical validation in captions and prose.

The first revised set is `exports/figures/2026-11-01-academic-v003` in the local
project. `scripts/plot_cylinder_academic.py` reproduces its three figures from
the repository's existing cylinder reports. Presentation changes do not alter
the underlying CFD data.
