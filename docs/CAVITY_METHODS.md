# Cavity methods and saved results

The [visual model](../lessons/10-lid-driven-cavity/README.md) records one closed
Re100 cavity. Selected native folders: `fluent-extended-v003` and
`openfoam`. Fluent original2000 continuity miss is retained. One separately admitted500-iteration saved-state continuation, without reinitialization or changed targets, produces selected2500. OpenFOAM retains its original2000 cap. Both selected native audits pass all declared checks; this is not physical validation.

The source mesh has 6,561 vertices and 6,400 uniform quadrilaterals. Exact native
topology maps bijectively to the source; native numbering is not assumed.
Fluent uses a planar depth of 1 m. OpenFOAM uses one 0.01 m empty-sided layer.
Density is 1 kg/m³ and dynamic viscosity 0.01 Pa s. Initial U and p are zero;
lid U=(1,0) m/s, other walls stationary, all no-slip and impermeable.
The top corner velocity is discontinuous. There is no energy or turbulence model.

Fluent uses pressure-based steady SIMPLE, least-square cell gradients, second-order
pressure and second-order-upwind momentum. OpenFOAM Foundation14 incompressibleFluid
uses steady SIMPLE consistent, Gauss linear gradients and its exact selected
[eight dictionaries](../cases/openfoam-cavity-selected-v001/README.md).
Native settings/readbacks and executed analyser sources accompany the results.

The original fixed caps were 2,000 steady iterations each. Fluent retains its original 2,000-iteration state and real continuity-residual miss. A separately authorized saved-state continuation adds exactly 500 iterations without reinitialization; final selection is at 2,500.
Iteration indices are not physical time. Original Fluent saves are
500/1000/1500/1800/2000; selected stability endpoints are 2300/2500.
OpenFOAM retains its original 2,000 cap and saves every 100, all twenty states.
Both the original parser-error audit and corrected original2000 numerical audit
are copied distinctly under results/cavity-v001/history, with their flags unchanged.
Final continuity/x/y Fluent scaled residual targets are 1e-6. OpenFOAM final printed
initial residual targets for p/Ux/Uy are 1e-6; these norms differ across codes.

Local conservation uses oriented native fluxes on every face: Fluent signed mass
sum divided by rho*cell volume; OpenFOAM signed volumetric phi divided by native
cell volume. Multiply by L/U=1 s. Max absolute target is 1e-5, RMS 1e-6.
Closed physical-wall net flux divided by rho*U*L*depth (mass) or U*L*depth
(volume) targets 1e-10. No zero-inlet normalization is used, and boundary net zero
alone is not an interior conservation check. Velocity relative L2 target is 1e-3:
Fluent 2300→2500, OpenFOAM 1800→2000. Fluent weights source area and
normalizes by U2300; OpenFOAM
weights volume and normalizes by U2000. Endpoint agreement is a bounded stability check.

| Actual quantity | Fluent | OpenFOAM |
|---|---|---|
| Numerical aggregate | PASS | PASS |
| Final residuals | {"continuity": 2.0132718363977406e-07, "x-velocity": 9.340144415900434e-09, "y-velocity": 8.43748404174094e-09} | {"p": 1.051285218369e-07, "Ux": 1.6623662041456e-07, "Uy": 2.0518826770771e-07} |
| Max scaled local divergence | 3.1414926e-09 | 3.653559e-09 |
| RMS scaled local divergence | 4.5468419e-10 | 7.2578593e-11 |
| Velocity relative L2 (Fluent 2300→2500; OpenFOAM 1800→2000) | 2.502203e-05 | 0.00010984871 |

Fluent misses: none.
OpenFOAM misses: none.
The JSON reports preserve all original booleans. Artifact integrity is verified
separately from these numerical outcomes. Earlier native versions and corrections
remain preserved. No target was relaxed; only the recorded, separately admitted
continuation was performed.

Figures use native U and p; OpenFOAM kinematic p is multiplied by rho. Each pressure
field removes its own area mean. Common scales are unclipped. Vorticity uses
second-order centre differences, one-sided at the envelope. Centreline values
interpolate adjacent native centres; no boundary values are appended. Streamlines
and CFL-limited RK4 tracers use actual saved centre velocity; no synthetic vortex
or forced closed trajectory is inserted. The 48×50 ms GIF replays tracer seeds,
not transient CFD. [Captions](../images/2026-11-01-cavity-v001/CAPTIONS.md) and
[plotter metadata](../results/cavity-v001/visuals/input-hash-manifest.json) retain details.

For reproduction use the existing E: Python -B/environment, D: Fluent and guarded
fresh versioned E: paths. Review saved launch/preparation/analyser sources, exact
settings and actual resource/licence admission. Runs are sequential, hidden and
two-process/rank; no additional run is launched by this record helper. Native
HDF/mesh/NPZ, full trace exports, dependencies and private connection/licence files
stay outside the repository. Source/destination copy hashes and the old public
manifest are retained privately. Numerical checks and pictures do not establish
benchmark/experimental validation or mesh independence.
