# Actual native planar elbow visualization

A uniform.1m/s inlet stream enters from the left and turns upward toward the top outlet. Density1kg/m3, dynamic viscosity.001Pa.s, inlet-width Re20. Horizontal arm x[-.6,.2],y[0,.2]m, vertical arm x[0,.2],y[.2,1]m. Shared3200quads with area.32m2; OpenFOAM one.01m centred extrusion. This is a planar sharp channel bend, with no three-dimensional secondary flow or Dean-vortex claim.

Complete native vertex sets map bijectively to common source quads; saved OpenFOAM C verifies true centres. Raw pressure retains the shared outlet gauge reference0Pa. No area mean is subtracted. Common scales include every native cell value. Sparse-grid voids remain masked; there are no invented field samples in the void.

RK4 uses actual bilinearly interpolated native-centre velocity only when all four surrounding fluid centres exist. Every stage and segment remains within the true union of rectangles. Paths stop at the native-centre envelope, walls or voids. Six seeds lie atx=-.585m, y=.025..175m inside the actual interpolation envelope. Repeated injections reuse these same native-velocity trajectories. Colour identifies particles, not concentration or composition. The fields stay fixed and no trajectory is forced closed.

Common travel-clock duration: 24.0s, replayed in2.4s. Seed injection times(s): [0.0, 6.0, 12.0, 18.0].

Native numerical audits and every missed criterion remain explicit. Pictures and solver agreement do not establish accuracy, benchmark agreement, mesh independence or physical validation. ANSYS Fluent/OpenFOAM computed the fields.

ANSYS Fluent: declared numerical checks passed
OpenFOAM: declared numerical checks passed
