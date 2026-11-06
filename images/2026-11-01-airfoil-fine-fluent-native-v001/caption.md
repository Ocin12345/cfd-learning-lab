# Unchanged native Fluent exports: fine R20/c airfoil

These three PNGs are byte-identical native Fluent exports from the saved-case continuation for Re = 100, zero incidence, R/c = 20 and 944,026 triangular cells. The saved source ends at iteration 172; native reconvergence is reported at 173, followed by 100 observation iterations through 273. The original legends, solver interpolation, individual colour ranges, background and ANSYS Student labels are preserved without image edits.

- **residuals.png**: native scaled continuity and two momentum residual histories through iteration 273, retaining the source history. The logarithmic display spans 1e-18 to 1e0, with the horizontal axis ending at 300.
- **velocity.png**: native velocity-magnitude contour in m/s; rounded displayed legend 0.00 to 1.09 m/s.
- **pressure.png**: native static gauge-pressure contour in Pa; rounded displayed legend -0.731 to 0.857 Pa. The operating pressure is 101325 Pa and outlet gauge pressure/reference pressure is 0 Pa. This image is not a pressure-coefficient plot.

The images retain their own native display scales and do not constitute a common-scale solver comparison. All three were visually inspected: legends, airfoil silhouette and residual axes are legible, with no clipped plot content observed. These numerical records do not establish physical validation, mesh independence, aircraft performance or independently authored owner work.
