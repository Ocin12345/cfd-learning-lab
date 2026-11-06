# Finite-wing OpenFOAM inputs: prepared version002

PREPARED ONLY. No version002 native conversion or flow is claimed.
The eight dictionaries are byte-identical to prepared version001 and are
reused unchanged for the one authorized local trailing-cap mesh repair.
Geometry, physics,1200iteration cap, all numerical targets and force references
remain unchanged. The repaired shared mesh hash, native quality and executed
sources must be recorded independently before a flow calculation is accepted.

Preparation helper: scripts/prepare_openfoam_finite_wing_v002.py.
Postprocessor: scripts/analyse_finite_wing_openfoam_v002.py.
Native case scope: native/runs/2026-11-01-finite-wing-v002/openfoam.
Separate root-authorized PrepareOnly and Run phases are required; no old
version001 authorization applies. All original version001 projects are preserved.

This prepared setup is not a converged, mesh-independent or validated prediction.
See the local version002 DESIGN.md and OPENFOAM-METHOD.md for the method,
software primary sources and the unchanged limitations.
