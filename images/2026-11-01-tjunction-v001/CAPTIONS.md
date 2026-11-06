# Actual native T-junction merge visualization

Two equal uniform .1m/s inlet streams enter the horizontal left arm and the upper stem. Density1kg/m3, dynamic viscosity.001Pa.s, inlet-width Re20; no scalar-mixing equation. Horizontal channel x[-.4,1.2],y[-.1,.1]m, upper stem x[-.1,.1],y[.1,.7]m. Shared4400quads with area.44m2; OpenFOAM one.01m centred extrusion.

Complete native vertex sets map bijectively to common source quads; saved OpenFOAM C verifies true centres. Raw pressure retains the shared outlet gauge reference0Pa. No area mean is subtracted. Common scales include every native cell value. Sparse-grid voids remain masked; there are no invented field samples in the void.

RK4 uses actual bilinearly interpolated native-centre velocity only when all four surrounding fluid centres exist. Every stage and segment remains within the true union of rectangles. Paths stop at the native-centre envelope, walls or voids. Seeds lie atx=-.385m ory=.685m inside the actual interpolation envelope. Repeated injections at travel-clock times0/4/8/12s reuse these same native-velocity trajectories. Entry colours are origin labels only, not concentration, composition or mixing. The common travel clock spans16.45s in2.4s playback; the fields stay fixed. No trajectory is forced closed.

Native numerical audits and every missed criterion remain explicit. Pictures and solver agreement do not establish accuracy, benchmark agreement, mesh independence or physical validation. ANSYS Fluent/OpenFOAM computed the fields; Codex prepared code/operations/visualization; no independent user execution or reflection is claimed.

ANSYS Fluent: declared numerical checks passed
OpenFOAM: declared numerical checks passed
