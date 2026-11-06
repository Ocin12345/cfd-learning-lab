# Saved channel-flow fields

Solvers shown: ANSYS Fluent, OpenFOAM. Uniform inlet 1 m/s; L = 1 m, H = 0.1 m; rho = 1 kg/m3, mu = 0.01 Pa.s; Re_H = 10; stationary no-slip plates and zero-gauge-pressure outlet. Each solver has 240 x 32 = 7,680 cells.

channel-computed-fields.png/svg uses actual native cell-centred values with flat per-cell colours. Velocity uses a common 0..1.6 m/s scale; pressure limits include every shown value and are shared. The vertical display scale is expanded; axis labels retain metres. White streamline curves and arrows use Matplotlib integration through the mapped saved vector field, without changing contour values. Fluent centres are reconstructed from four unique saved face-connected vertices per rectangular cell; OpenFOAM centres come from its saved C field. Both centre sets must map bijectively within 1e-9 m. OpenFOAM kinematic pressure is multiplied by rho = 1 kg/m3.

channel-velocity-profiles.png/svg shows nearest saved columns to x = 0.02, 0.1, 0.3, 0.8 m; exact column positions are printed in the panels and CSV. The dotted parabola is an analytical fully developed reference, not simulated data. Curves connect actual cell-centred samples and do not append wall values. The near-inlet profiles are developing and are not required to match that downstream reference.

Analytical targets passed: true. See checks.json for measured profile and pressure-gradient checks. Residual, mass-balance and iteration-stability evidence is checked separately. This is a bounded educational analytical check, not mesh independence or experimental validation. Fluent/OpenFOAM calculated the fields; Codex reconstructed their saved data, checked the mappings and prepared the figures.

Native Fluent screenshots fluent-native-velocity.png and fluent-native-pressure.png are unchanged exports of the completed channel case. Their own native colour maps and interpolation are preserved; use the common-scale paired figures for comparisons. Source case/data and numerical records remain outside Git on E:.
