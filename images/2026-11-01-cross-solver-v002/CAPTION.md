# Figure 4. Matched Fluent–OpenFOAM cylinder calculation at Re = 20

(a) The continuation of OpenFOAM drag convergence after the initial
1,200-iteration cap; the earlier startup interval is omitted here but retained
in the source CSV. The dashed line is the terminal Fluent value, not a Fluent
iteration history. (b) OpenFOAM initial linear-solver residuals for pressure and
the two in-plane velocity components. The vertical dotted line marks the
initial cap and the horizontal dotted line the declared 1e-7 residual target.
These residuals do not share Fluent's scaling and are not overlaid with its
residuals as though they were equivalent error norms.

The same 98,304-cell O-grid was used, with outer radius R/D = 20. OpenFOAM used
a one-cell, unit-span extrusion with empty front/back patches; converted x/y
coordinates agreed with the Fluent source mesh to the stored precision. Both
cases were steady, incompressible and laminar, with D = 1 m, U = 1 m/s,
rho = 1 kg/m3, mu = 0.05 Pa s and unit reference area. The upstream velocity
inlet, downstream zero-pressure outlet, no-slip cylinder and force directions
were matched. Numerical operators differ: coupled second-order Fluent versus
steady SIMPLE and bounded linear-upwind OpenFOAM.

Terminal Cd was 2.0731581 in Fluent and
2.0748298 in OpenFOAM; the relative difference normalized by
Fluent was 0.08064%. OpenFOAM's final
100-iteration Cd range was 0.000037%
of its window mean. Meeting the declared 0.5% cross-code criterion is numerical
agreement for this confined setup, not physical validation or an uncertainty
estimate. Fluent reused a saved converged state; OpenFOAM started from uniform
flow. Their iteration counts and run times are not a performance benchmark.

Codex prepared the cases, analysis, figure and caption. ANSYS Fluent and
OpenFOAM computed their respective solutions. The owner requested the paired
workflow and academic presentation; independent owner analysis is not inferred.
