# planar channel bend methods and evidence

[Visual model](../lessons/12-channel-bend/README.md). Both original fixed2000 calculations completed and passed all declared native numerical checks without extra iterations or target changes. Source-bound native fields retain raw outlet gauge pressure0. The24s tracer clock and6s seed period cover actual10.23-17.18s computed transit times in both saved velocity fields.
Geometry is horizontal x[-0.6,0.2], y[0,0.2] plus vertical x[0,0.2], y[0.2,1.0] m.
Uniform 0.01 m cells give 3,200 source quads, 3,381 nodes and area 0.32 m².
Native topology maps to the 100×80 sparse source; 4,800 void cells contain no flow.
Fluent depth is 1 m; OpenFOAM extrusion depth is 0.01 m with empty planes.
Left inlet U=(0.1,0) m/s; upper outlet gauge p=0; other boundaries no-slip.
Density 1 kg/m³, mu=0.001 Pa s, H=0.2 m. No energy, turbulence or scalar equation.
Fluent uses steady pressure-based SIMPLE, least-square gradients, second-order
pressure/momentum. Exact OpenFOAM SIMPLE/Gauss linear/schemes are retained in
[eight dictionaries](../cases/openfoam-elbow-selected-v001/README.md).

Both caps are 2,000; stability compares 1800/2000. Residual targets are 1e-6,
using each solver's own norm. All-cell oriented mass/volume flux divided by
rho*V or V is scaled by H/U=2 s: max target 1e-5, RMS 1e-6. Global imbalance is
absolute net over the actual single inlet flow, target 1e-4; no closed-wall
zero-flux substitute is used. Velocity relative L2 target is 1e-3, weighting
area (Fluent, denominator U1800) or volume (OpenFOAM, denominator U2000).

| Actual quantity | Fluent | OpenFOAM |
|---|---|---|
| Final residuals | {"continuity": 1.3311030732529432e-14, "x-velocity": 4.46408572826261e-17, "y-velocity": 4.955788215864678e-17} | {"p": 1.0303925290148e-10, "Ux": 1.7891856057138e-12, "Uy": 4.8483386785619e-13} |
| Max scaled divergence | 8.673617379884088e-15 | 1.770552969645501e-10 |
| RMS scaled divergence | 2.1310034219642754e-15 | 2.4637770345784645e-11 |
| Relative mass imbalance | 4.336808689942018e-17 | 7.360377498460968e-13 |
| 1800→2000 velocity relative L2 | 2.0402034717111797e-16 | 8.280548402150795e-12 |

Original [Fluent audit](../results/elbow-v001/fluent/fluent-checks.json)
and [OpenFOAM audit](../results/elbow-v001/openfoam/checks.json)
retain flags unchanged. Structural/parser/source errors block selection; genuine
numerical misses remain false. Pressure is raw gauge Pa (OpenFOAM rho*p), with
common outlet zero, without mean subtraction. Only dimensionless metrics are
compared above; native raw flows scale with depth. Figures use actual native
fields. Fixed-field RK4 travel time differs from GIF playback and solver time.
Captions/plotter hashes retain the exact interpolation and invocation.

Use fresh guarded versioned E: outputs and existing E: Python -B/D: Fluent;
sequential hidden two-process/rank runs require resource/licence admission.
Native HDF/mesh/NPZ, dependencies, full traces and private launch data remain
external. Original sources/plan and prior projects are preserved; copy hashes,
backups and separate root image QA/artifact receipt remain private. No validation,
mesh independence, 3D secondary flow, Dean vortices, owner reflection, Git commit,
publication or submission is claimed.
