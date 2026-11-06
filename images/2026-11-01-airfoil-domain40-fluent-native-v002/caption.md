# Unchanged native Fluent exports: R/c = 40 airfoil domain

These three PNGs are byte-identical exports from the final Fluent saved-case continuation, at iteration 810, for the Re = 100, zero-angle airfoil calculation on 531,018 cells. They retain the solver's original rendering, labels, interpolation and individual display ranges. No image cropping, recolouring or annotation has been applied.

- **residuals.png**: native scaled-residual history through iteration 810, including the retained earlier history. The logarithmic vertical display spans 1e-14 to 1e0; continuity and the two momentum equations use Fluent's own residual definitions.
- **velocity.png**: native velocity-magnitude contour, in m/s. The displayed colour legend runs from 0.00 to 1.09 m/s (rounded labels).
- **pressure.png**: native static gauge-pressure contour, in Pa, with operating pressure 101325 Pa and outlet gauge/reference pressure 0 Pa. The displayed colour legend runs from -0.667 to 0.842 Pa (rounded labels); this is pressure, not a pressure-coefficient plot.

These native exports use their own display scales. The separate paper figures provide the common-range field comparison. The images document this numerical calculation and do not establish physical validation or acceptance of the OpenFOAM comparison.

This v002 caption corrects the earlier v001 operating-pressure description; the three PNG files and computed solution are unchanged. See the [project erratum](../../docs/ERRATA.md).
