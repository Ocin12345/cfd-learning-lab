# Nozzle methods and evidence

The [visual model](../lessons/09-nozzle/README.md) records one bounded paired
compressible-flow model. The selected native states are `fluent-order-control-v004`
and `openfoam-order-control-v002`. One first-order control per code at unchanged mesh, physics, initialization, CFL schedule and original iteration caps; better numerical behaviour, with every genuine remaining miss retained. No additional CFD or iteration extensions.

## Geometry, gas and boundary conditions

The full symmetric nozzle is 1 m long. Its half-height is 0.10 m at the inlet,
0.04 m at x=0.4 m and 0.08 m at the exit. Each segment uses a cosine interpolation
with zero endpoint slope. The 240 axial by 40 transverse cell rows use uniform
eta=y/h(x). Exact polygon centroids and native vertex connectivity are checked;
imported cell ordering is not assumed. Native mesh quality and positive-volume
results are retained in each numerical audit.

The model solves inviscid adiabatic Euler flow, not viscous nozzle or engine
performance. Nominal R=287.05 J/(kg K), γ=1.4 and Cp=1004.675 J/(kg K).
Nominal R287.05J/(kgK), gamma1.4 and Cp1004.675 match; inferred selected native EOS R medians are Fluent287.0503401690 and OpenFOAM287.0501716835, about1ppm installed rounding. Exact runtime gas parity is not claimed; tiny rounding is not assumed to cause instability. No exact runtime gas-constant parity is claimed.
Inlet stagnation pressure/temperature are 200000 Pa/300 K; outlet absolute
static pressure is 140000 Pa. Slip walls are impermeable. Fluent operating
pressure is zero, so its native static pressure is absolute. The uniform initial
state is p=140000 Pa, T=300 K, U=(50,0) m/s.

Fluent is planar with reference depth 1 m. OpenFOAM imports one empty-sided
hexahedral layer of depth 0.01 m. Native OpenFOAM phi is mass flow in kg/s through
that layer; dividing by 0.01 gives kg/s per m depth for comparison. Relative mass
imbalance uses the absolute net physical-boundary flux divided by inlet magnitude.
Both slip walls are included. The reported per-depth inlet values are
37.076488 for Fluent and 37.049748 for OpenFOAM; agreement is not an acceptance
or validation claim.

## Actual selected numerical methods and checks

Fluent uses density-based implicit Roe-FDS and a fixed total of 2000 iterations.
The first 300 iterations use first-order upwind at CFL 1. The remaining 1700
use the selected order at CFL 3. The final recorded amg-c scheme is `first-order-upwind` at CFL 3.0.
The source-bound audit retains initialization and stage settings, every saved
checkpoint, native scaled residuals, per-iteration inlet/outlet mass and an
independent mass report. Its declared targets are continuity/x/y residuals
≤1e-5, energy residual ≤1e-7, relative mass imbalance ≤1e-4 and final 100 actual
mass-flow relative range ≤1e-3. Final residuals are {"continuity": 0.0004769478834009375, "x-velocity": 0.0008498024235996194, "y-velocity": 0.001005544188260254, "energy": 0.0004222352198547506}.
Final 100 mass-flow metrics are retained in the report and CSV. Named misses:
final_residual_targets. No target or original iteration cap was increased.

OpenFOAM Foundation 14 uses shockFluid, inviscid mu=0, `Kurganov`
fluxes, `localEuler` local time stepping and reconstruction
{"rho": "upwind", "U": "upwind", "T": "upwind"}. The fixed budget is 6000 pseudo-iterations
at maxCo=0.4, with saved fields every 500. These indices are not physical time.
Positive finite p/T/rho and an EOS consistency check are distinct from convergence.
Targets are final native and last 100 maximum relative mass imbalance ≤1e-4,
last three saved-state relative L2 changes U/p/T ≤1e-3 and compression-crossing
location range ≤0.01 m when applicable. The saved comparisons are 5000/5500/6000;
endpoint comparisons do not establish every-iteration stability. Diagonal linear
solver residual zeros are excluded from acceptance. Named misses: final_native_relative_mass_imbalance_le_1e_4, last100_phi_history_max_relative_imbalance_le_1e_4, last_three_saved_states_relative_L2_U_p_T_le_1e_3.
All flags and source hashes remain in the selected reports.

The selected first-order treatment adds numerical diffusion and can broaden compression features. The change is an order-only numerical control on the same geometry,
initial/boundary conditions and budgets. It does not establish the cause of the
earlier instability, a physical shock thickness or an accurate shock position.
The original Fluent higher-order state and OpenFOAM vanLeer state remain separate
failed diagnostics, with their original flags and overshoots retained externally.
Two Fluent API-only setup failures and a zero-iteration settings probe also remain
preserved; they are not flow solutions. Human-facing summaries omit raw private
tracebacks. Earlier public projects and all rejected native versions are preserved.

## Figures and analytical context

Mach is derived from saved U and static T using |U|/sqrt(γ R T). Figures share
unclipped scales and show actual selected numerical status and schemes. Native
quad/hex vertex sets map bijectively to the source quads. Pressure is absolute.
The centreline samples use the adjacent native cell rows; any compressive M=1
crossing is a candidate observation, not normal-shock validation. Native
interpolation, RK4 path termination and display details are in the
[figure captions](../images/2026-11-01-nozzle-v001/CAPTIONS.md).
The 48-frame GIF has 50 ms delays. Integration travel time, playback time and
OpenFOAM pseudo-iteration indices are distinct; none is a computed transient
nozzle time history. No analytical shock, synthetic velocity or wall sample is
inserted. Sparse mesh/interpolation and first-order diffusion remain limitations.

The [quasi1D context JSON](../results/nozzle-v001/context/quasi1d-context-v001.json)
is explicitly analytical, not CFD and not a matched 2D validation. It assumes
uniform section flow and isentropic branches separated by an ideal normal shock.
Its estimates are not substituted into the CFD plots or convergence checks.
Primary context: [NASA isentropic relations](https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html)
and [NASA normal-shock relations](https://www.grc.nasa.gov/www/k-12/airplane/normal.html).
OpenFOAM implementation sources are recorded with the exact selected dictionaries
and preparation helpers. No mesh/domain independence, matched benchmark accuracy,
experimental validation, research novelty or owner reflection is established.

## Reproduction and preservation

Use the existing E: Python environment and D: ANSYS installation. Native outputs,
working directories and supported caches/temp remain on E:. Launches require
4 GiB available host RAM, 10 GiB free E: space, no other solver and sequential
hidden two-process/rank sessions. Live Fluent startup establishes licence admission.
Save all fields before closing; stop only verified-idle dedicated CFD-Foundation14.
No solver, WSL, Git or publication action is performed by the record helpers.

The selected executed driver/analyser, source hashes, numerical tables and exact
OpenFOAM dictionaries are retained under results/nozzle-v001 and
cases/openfoam-nozzle-selected-v001. The figure export preserves its executed
plotter and exact invocation in the input manifest. Native HDF/mesh/NPZ and full
trace exports remain under the versioned external E: folders. For replay, create
fresh versioned paths and explicitly review guarded path-only script copies;
never overwrite these records or resume a rejected case automatically. The same
unchanged numerical criteria apply to a replay. No extra mesh, parameter sweep,
cap extension, dependency installation, Git staging/commit/push or submission
was added by the record. Reflection remains unprovided.
