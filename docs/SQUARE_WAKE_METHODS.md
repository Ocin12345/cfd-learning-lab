# Confined square-obstacle wake methods and evidence

[Actual visual model](../lessons/14-square-obstacle-wake/README.md).
Select the completed fixed 2,000-iteration Fluent run and corrected OpenFOAM v003 setup. Both own numerical audits pass; their separate predeclared paired drag target also passes. Preserve all earlier preparation and parser evidence; no numerical threshold was relaxed.

Geometry is x[0,2], y[-0.4,0.4] m minus square x[0.45,0.55], y[-0.05,0.05] m.
The 200 by 80 full grid omits the actual 100 body cells, leaving 15,900 quads,
16,200 nodes and area 1.59 m². The native boundary has 80 inlet, 80 outlet,
200 top, 200 bottom and 40 body faces. Fluent depth is 1 m (volume 1.59 m³);
OpenFOAM depth is 0.01 m (volume 0.0159 m³), with empty extrusion planes.

Inlet U=(0.1,0) m/s, outlet gauge p=0, body stationary no-slip and top/bottom
symmetry/slip. Density 1 kg/m³, mu = 0.0005 Pa s and D = 0.1 m give Re_D = 20. Initial U/p
are zero in both codes. This is a finite domain of height 8D and 12.5% blockage,
with body-centre distances 5D from inlet and 15D to outlet. No energy, turbulence,
scalar or motion equation. Exact selected SIMPLE/spatial settings are retained
in [OpenFOAM dictionaries](../cases/openfoam-square-wake-selected-v001/README.md)
and Fluent native settings. Both fixed caps are 2,000; no automatic extension.

| Actual quantity | Fluent | OpenFOAM |
|---|---|---|
| Final own residuals | {"continuity": 1.5902690993284656e-14, "x-velocity": 2.8135862221901946e-17, "y-velocity": 8.95145401919294e-18} | {"p": 1.0201347945396e-10, "Ux": 1.7510340110079e-13, "Uy": 3.3075075570158e-13} |
| Max scaled native divergence | 4.336808689942106e-15 | 8.295840685408618e-12 |
| RMS scaled native divergence | 1.1593612836983962e-15 | 4.801502385066377e-13 |
| Relative global mass imbalance | 8.673617379884035e-17 | 1.7753810574450135e-14 |
| 1800 to 2000 velocity relative L2 | 7.269861682063812e-16 | 1.1882748851932893e-12 |

Final residual targets are 1e-6 using each code's own norm. All native oriented
face fluxes contribute to cell divergence, scaled by D/U = 1 s: max target 1e-5,
volume-weighted RMS 1e-6. Global imbalance uses the actual inlet flow, target 1e-4.
Body/symmetry normal flux is independently checked. Velocity stability compares
1800 and 2000; relative L2 target 1e-3, weighting area in Fluent (denominator 1800)
or volume in OpenFOAM (denominator 2000).

Native force histories retain every actual iteration 1..2000. The exact final
window 1801..2000 has 200 samples. Body-only references are rho = 1, U = 0.1, D = 0.1,
area = D × depth: qA = 0.0005 N for Fluent and 0.000005 N for OpenFOAM. The Cd
peak-to-peak / abs(mean Cd) target is 0.001 and abs(mean Cl) target is 0.001. Zero mean drag
produces an undefined relative range and false target, not an invented ratio.

| Actual force-window quantity | Fluent | OpenFOAM |
|---|---|---|
| Mean Cd, 1801..2000 | 2.8144310326314543 | 2.8187666062523586 |
| Mean Cl, 1801..2000 | -2.0555372047493733e-16 | -7.994678121012351e-14 |
| Cd peak-to-peak / abs(mean Cd) | 1.1045303412689837e-15 | 2.930380715057917e-12 |
| Cd stability target passed | True | True |
| Cl symmetry target passed | True | True |

The separate paired mean-Cd difference abs(meanCd_F-meanCd_OF)/abs(meanCd_F)
is **0.154048%**; its predeclared <= 1% target is **PASS**. A zero Fluent
mean is undefined/false. This comparison does not change either solver's own
numerical audit and agreement would not establish physical validation.

At saved endpoints, Fluent pressure force integrates native SV_P with oriented
body face area. Native SV_WALL_SHEAR already contains integrated force; its
sign/units are established against this case's own force reports before use,
with no second face-area multiplication. OpenFOAM native forces/forceCoeffs
are crosschecked every iteration; saved pressure and wallShearStress boundary
integrals independently crosscheck forces. Its wallShearStress is kinematic
traction, requiring -rho*area integration. Zero internal wall-stress storage
is not claimed as an interior stress solution. Native report-precision force
crosscheck tolerances remain exactly as declared before execution.

[Fluent audit](../results/square-wake-v001/fluent/fluent-checks.json) and
[OpenFOAM audit](../results/square-wake-v001/openfoam/checks.json)
retain flags unchanged. Parser/source/structural errors block selection; genuine
numerical misses stay false. Large meshes/HDF/NPZ, full field tables, tracer
paths and logs remain external; compact force/residual CSV and settings are
copied unchanged. Sources and prior projects are preserved.

Figures use actual saved fields at iteration 2000, common unclipped scales and
raw gauge Pa (rho*p in OpenFOAM). No pressure-mean subtraction. The body hole
is masked; all four interpolation centres must be native fluid cells, and RK4
stages/segments touching the body are rejected. There is no invented void data,
forced path closure or transient vortex-shedding claim. Tracer travel
40 s / injection period 6 s differs from
2.4 s GIF playback and steady solver iterations. All five PNGs, the GIF and the
native ANSYS velocity image require actual root inspection before public record writes.

Hidden sequential runs use two Fluent processes or two OpenFOAM ranks, the
existing Fluent installation on D:, Python on E: and versioned outputs, temporary
files and caches on E:. Live resource, licence and native-quality checks precede
each run.
No mesh/domain independence, benchmark/experimental validation, independent
owner implementation/reflection, novelty, Git commit, publication or submission
is claimed.
