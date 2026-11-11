# Lid-driven cavity: actual native field visualization

Closed 1 m square, lid velocity (1,0) m/s, stationary side/bottom walls, density 1 kg/m3, dynamic viscosity.01 Pa s, Re100. Both solvers compute a steady laminar field on the same 80x80 source quads; OpenFOAM uses one.01 m extruded layer.

Native Fluent face connectivity and native OpenFOAM faces/owner/neighbour reconstruct complete cell vertex sets. A verified bijection maps them to canonical quads. Saved OpenFOAM C independently verifies true polygon centres. Native pressure gauges differ; each pressure field has its own spatial area mean removed. OpenFOAM kinematic pressure is multiplied by rho. Shared scales include all values.

Streamlines and RK4 tracers use only bilinearly interpolated actual saved cell-centre velocity. No boundary values, synthetic vortex, forced closed orbit or solver-time evolution is inserted. Tracers stop at the interpolation envelope or wall, and use a common travel-time clock. GIF replay restarts seeds. Integration travel time differs from playback time. Vorticity is a derived second-order cell-centre gradient, one-sided at the outermost centres.

Every failed numerical criterion from the hash-bound audits remains visible in the figures. A closed-domain zero boundary flux alone does not demonstrate local mass conservation. Solver agreement and attractive flow pictures do not establish benchmark accuracy, mesh independence or physical validation.

Selected final steady snapshots: Fluent iteration 2500, OpenFOAM iteration 2000. Solver iteration counts are independent convergence indices, not a common physical clock.

ANSYS Fluent and OpenFOAM computed the fields.

ANSYS Fluent: declared numerical checks passed
OpenFOAM: declared numerical checks passed
