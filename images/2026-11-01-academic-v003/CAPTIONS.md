# Figure captions

**Figure 1. Cylinder domain-size sensitivity at Re = 20.** (a) Drag coefficient
as a function of the outer-boundary radius. (b) Successive relative difference,
epsilon_R = 100 |C_D(R_i) - C_D(R_{i-1})| / |C_D(R_i)|. The dashed line denotes
the 0.2% criterion declared before the subsequent comparisons. The original
fine-grid coordinates through R/D = 20 were retained exactly; additional outer
rings were generated for larger domains. Each new case reached the 1e-8 residual
target and was extended by 100 iterations. Markers denote computed values;
connecting lines are guides to the eye. Meeting this criterion for drag does
not establish an error bound, convergence of all fields, or physical validation.

**Figure 2. Iterative force stability of the saved Re = 20 cylinder calculation.**
The saved solution at k = 51 was extended by 100 iterations on 98,304 cells with
R/D = 20. (a) Drag increment relative to the saved state, multiplied by 10^7.
(b) Magnitude of the lift coefficient. The forces were sampled every iteration
from Fluent reports with approximately eight significant digits; the staircase
in (a) and zero last-window drag range reflect that reporting precision. Very
small late lift values are numerical, rather than evidence of physical flow
oscillations. These are solver iterations, not physical time steps.

**Figure 3. Scaled residual histories for the Re = 20 cylinder replay.**
The vertical dotted line marks the saved state at k = 51; the horizontal dotted
line marks the original 1e-7 stopping criterion. Residual stopping checks were
disabled during the 100-iteration observation window while residual monitoring
continued. Fluent's residuals are shown using their exported scaling and cannot
be compared numerically to OpenFOAM residuals without reconciling definitions.

**Common conditions.** Steady planar incompressible laminar flow, D = 1 m,
U_infinity = 1 m/s, rho = 1 kg/m^3, mu = 0.05 Pa s, unit span and reference
area 1 m^2; no-slip cylinder; upstream velocity inlet and downstream zero-gauge
pressure outlet. Fluent Student 2026 R1, two solver processes, pressure-based
coupled method and second-order pressure/momentum discretization. No uncertainty
bars have been drawn because a defensible uncertainty estimate is not available.

**Attribution.** Codex prepared the plotting code and captions. ANSYS Fluent
computed the flow and exported the numerical reports. The owner requested the
study and the revision to academic presentation; no independent owner analysis
or research novelty is inferred.
