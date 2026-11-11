# Actual native nozzle visualization

Planar converging-diverging nozzle, length 1 m, full inlet/throat/exit heights 0.20/0.08/0.16 m, throat x = 0.4 m. Both codes solve inviscid ideal-gas Euler flow with gamma = 1.4, R = 287.05 J/(kg K), total inlet pressure 200 kPa, total inlet temperature 300 K and absolute outlet back pressure 140 kPa. Shared 240 x 40 wall-fitted quads. Numerical acceptance is recorded separately.

Native Fluent cells are reconstructed from saved face connectivity. Native OpenFOAM hexahedra are reconstructed from saved faces/owner/neighbour and projected to source quads. Sorted vertex sets give a bijective common-cell map; imported ordering is not assumed. Saved OpenFOAM C, when available, independently checks true polygon centroids. Native pressure is absolute (Fluent operating pressure 0 Pa). Native temperature and velocity give Mach = |U| / sqrt(gamma R T). No one-dimensional solution, shock position, wall velocity or synthetic flow is inserted. Both codes use common unclipped colour scales.

Centreline curves interpolate between the two adjacent actual cell rows at y = 0. An adjacent M > 1 to M < 1 transition and positive pressure jump is only a computed compression/shock candidate. It does not validate shock location or strength. The CSV preserves all samples.

Animation particles follow RK4 integration through velocity interpolated bilinearly in axial and normalized-height coordinates using mapped cell centres. No boundary velocity values are appended. Paths stop at the native cell-centre interpolation envelope or actual wall. Recycled phases and playback timing move particles through each fixed saved snapshot. This is not physical solver-time evolution; OpenFOAM localEuler is a pseudo-time steady-convergence procedure. Exported integration travel times differ from playback and solver pseudo-iterations.

ANSYS Fluent and OpenFOAM computed the fields. No user reflection, independent user execution, mesh/domain independence, shock accuracy, benchmark agreement or physical validation is claimed. Source hashes and executed source accompany the figures; all inputs are rehashed after output.

Actual numerical status: ANSYS Fluent: 1 numerical check missed | OpenFOAM: 3 numerical checks missed.

Actual native scheme readbacks: {
  "fluent": {
    "source": "E:\\CFD-Learning-Lab\\native\\runs\\2026-11-01-nozzle-v001\\fluent-order-control-v004\\final-settings.json",
    "amg_c": "first-order-upwind",
    "flux_type": "Roe-FDS",
    "courant_number": 3.0,
    "short_label": "first order upwind",
    "first_order": true
  },
  "openfoam": {
    "source": "E:\\CFD-Learning-Lab\\native\\runs\\2026-11-01-nozzle-v001\\openfoam-order-control-v002\\system\\fvSchemes",
    "reconstruct": {
      "rho": "upwind",
      "U": "upwind",
      "T": "upwind"
    },
    "flux_scheme": "Kurganov",
    "time_scheme": "localEuler pseudo-time",
    "short_label": "first-order upwind",
    "first_order": true
  }
}.

A selected first-order control introduces additional numerical diffusion and can smear compression/shock structure. It is a robustness diagnostic, not a mesh-independence or shock-accuracy result.
