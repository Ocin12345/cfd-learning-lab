# Selected nozzle OpenFOAM inputs

These eight files are exact byte copies from the completed selected native case
`openfoam-order-control-v002`. Numerical outcome: **MISS (retained diagnostic)**. Their presence does
not imply convergence or physical validation. The original higher-order case is
preserved separately. Full methods:../../docs/NOZZLE_METHODS.md.

Fixed cap 6000; maxCo 0.4; localEuler pseudo-time; actual reconstruction
{"rho": "upwind", "U": "upwind", "T": "upwind"}. Native geometry/fields remain external;
use fresh guarded versioned paths for replay. No credentials or installation
files are included.
