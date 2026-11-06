# T-junction methods and evidence

[Visual model](../lessons/11-t-junction/README.md). Both original fixed2000 calculations completed and passed their independent declared native numerical checks. Selected actual fields, unchanged source hashes and raw outlet gauge pressure0 are bound to the paired visualization; no additional iterations or changed targets were needed.
Geometry is x[-0.4,1.2], y[-0.1,0.1] plus stem x[-0.1,0.1], y[0.1,0.7] m.
Uniform 0.01 m cells give 4,400 source quads, 4,641 nodes and area 0.44 m².
Native topology maps to the sparse source; void background cells contain no flow.
Fluent depth is 1 m; OpenFOAM extrusion depth is 0.01 m with empty planes.
Left U=(0.1,0), top U=(0,-0.1) m/s; outlet gauge p=0; other boundaries no-slip.
Density 1 kg/m³, mu=0.001 Pa s, H=0.2 m. No energy, turbulence or scalar equation.
Fluent uses steady pressure-based SIMPLE, least-square gradients, second-order
pressure/momentum. Exact OpenFOAM SIMPLE/Gauss linear/schemes are retained in
[eight dictionaries](../cases/openfoam-tjunction-selected-v001/README.md).

Both caps are 2,000; stability compares 1800/2000. Residual targets are 1e-6,
using each solver's own norm. All-cell oriented mass/volume flux divided by
rho*V or V is scaled by H/U=2 s: max target 1e-5, RMS 1e-6. Global imbalance is
absolute net over the actual sum of two inlet flows, target 1e-4; no closed-wall
zero-flux substitute is used. Velocity relative L2 target is 1e-3, weighting
area (Fluent, denominator U1800) or volume (OpenFOAM, denominator U2000).

| Actual quantity | Fluent | OpenFOAM |
|---|---|---|
| Final residuals | {"continuity": 1.8000808383934613e-14, "x-velocity": 6.614958796638744e-17, "y-velocity": 5.339744320163698e-17} | {"p": 6.9621972978725e-11, "Ux": 2.5162252458652e-16, "Uy": 1.5850248017078e-16} |
| Max scaled divergence | 1.7347234759768175e-14 | 1.16918329401762e-09 |
| RMS scaled divergence | 3.589591977262608e-15 | 1.9324735774421951e-10 |
| Relative mass imbalance | 1.6263032587282564e-16 | 1.1460694589529585e-11 |
| 1800→2000 velocity relative L2 | 2.5904584844117e-16 | 2.8628525205666777e-15 |

Original [Fluent audit](../results/tjunction-v001/fluent/fluent-checks.json)
and [OpenFOAM audit](../results/tjunction-v001/openfoam/checks.json)
retain flags unchanged. Structural/parser/source errors block selection; genuine
numerical misses remain false. Pressure is raw gauge Pa (OpenFOAM rho*p), with
common outlet zero. Figures use actual native fields; seed colours are entry
labels only. Fixed-field RK4 travel time differs from GIF playback and solver time.
Captions/plotter hashes retain the exact interpolation and invocation.

Use fresh guarded versioned E: outputs and existing E: Python -B/D: Fluent;
sequential hidden two-process/rank runs require resource/licence admission.
Native HDF/mesh/NPZ, dependencies, full traces and private launch data remain
external. Original sources/plan and prior projects are preserved; copy hashes,
backups and separate root image QA/artifact receipt remain private. No validation,
mesh independence, owner reflection, Git commit, publication or submission is claimed.
