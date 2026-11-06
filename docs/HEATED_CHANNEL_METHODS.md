# Heated-channel methods and evidence

[Actual temperature evolution](../lessons/15-heated-channel/README.md).
Cold fluid warms rapidly from the 350 K walls. Both solvers computed 600 thermal steps to 3 s on their separately frozen flow fields; both final outlet bulk temperatures are 349.9541 K. Their primary numerical checks pass. Auxiliary OpenFOAM gradient checks remain unresolved: the strict stencil comparison fails in 509 of 600 steps and the printed-precision interval check fails in 357 of 600. These diagnostics remain separate from the global energy checks; the production mismatch has not been explained by proven print-precision consistency.

The 1 m by 0.1 m channel uses 96 by 24 uniform quads (2,304 cells).
Fluent is planar with depth 1 m; OpenFOAM uses a 0.01 m extrusion with empty
end planes. The generic teaching fluid has density 1 kg/m³, viscosity
0.01 Pa s, heat capacity 1,000 J/(kg K) and conductivity 10 W/(m K).
Thus nu = alpha = 0.01 m²/s, Pr = 1 and Re_Dh = 20 for Dh = 0.2 m.
Uniform inlet velocity is (1,0) m/s; the outlet gauge pressure is zero;
stationary top and bottom walls are no-slip. This is developing flow, with
no fully developed profile imposed as an accuracy gate.

Each solver first computes exactly 2,000 steady SIMPLE iterations from rest.
Native residual, mass, all-cell face-flux divergence and 1800-to-2000 velocity
checks establish its own flow status. Final residual target is 1e-6 in each
code's own norm. Mass imbalance / actual inlet flow target is 1e-5.
Divergence is scaled by Dh/U = 0.2 s; max target 1e-5 and weighted RMS 1e-6.
Velocity relative L2 uses the final 2000 field in both codes; target 1e-3.
The actual native U, p and every face mass/volume flux are then frozen.
Temperature is patched without reinitializing or recomputing flow. No
buoyancy or time-dependent momentum is simulated.

Initial and inlet temperature is 300 K; walls are 350 K. Outlet temperature
has zero gradient, with Fluent backflow temperature 300 K where exposed.
Inlet diffusion is enabled. Gravity, radiation, viscous heating, pressure
work and kinetic-energy work are disabled or verified inactive.
Both codes use backward Euler, dt = 0.005 s and first-order upwind thermal
advection. The accepted one-step pilot supplies the initial state and solved step 1;
production adds 599 steps, for 600 solved steps and 601 native states to 3 s.
Fluent caps 60 inner iterations per step; OpenFOAM uses one scalar matrix
solve. Terminal own energy/scalar residual target is 1e-8 at every step;
a new-step initial linear residual is recorded separately. Native time,
index, dt, all temperatures and all frozen fields are checked at every state.
Bounds are 300–350 K with 0.01 K tolerance; frozen values must change by no
more than 1e-12 in their respective SI units. No automatic extension or
mesh/time accuracy-order claim is made.

| Actual quantity | Fluent | OpenFOAM |
|---|---|---|
| Maximum own terminal thermal residual | 9.885642979e-09 | 9.940264181e-13 |
| Maximum native global Euler defect | 5.124269272e-06 | 3.838312833e-11 |
| Maximum reconstructed PDE defect | 1.915635995e-06 | 3.843478344e-11 |
| Maximum frozen-value change (respective SI units) | 0 | 0 |
| Final energy above 300 K (J/m) | 4474.186341 | 4474.255695 |
| Final convective outlet bulk temperature (K) | 349.9541064 | 349.9541279 |

Energy above 300 K per unit depth is E = rho cp sum(V(T-300))/depth.
Storage is (E_n-E_previous)/0.005 s. Every physical boundary contributes
signed power, positive inward. Normalize the absolute global Euler defect by
max(sum(abs(Qin_by_boundary)), abs(storage), 5,000 W/m), target 1e-5 at every
solved step. The separate reconstructed PDE-flux diagnostic uses the same
normalization and its original 0.005 target. These diagnostics remain
separate. Neither establishes unobserved all-cell energy-flux conservation.

Fluent total native boundary heat includes sensible enthalpy convection.
The actual initial h300 is measured from native enthalpy, and native
h-h300 = cp(T-300) is verified. Reference-corrected integrated native W is
used once, without adding convection again or multiplying by face area.
Inlet conduction is included. Initial unsolved heat arrays have not evaluated
a thermal balance: their open-boundary reference diagnostic remains false
and is not treated as a valid zero or as a required solved-step target.

OpenFOAM diffusive inward power is k snGrad(T) area; convection is
-rho cp phi (T_upwind-300), using actual native outward volumetric fluxes.
Retained rank-native boundary gradT, positive signed faceProcAddressing,
owner-cell/vertex mappings and outward face vectors prove coverage of all
240 physical faces exactly once. Reconstructed serial extrapolatedCalculated
gradT may reset patch gradients on reload and remains diagnostic only.
The original strict 1e-10 stencil comparison and 14-digit precision-interval
checks remain auxiliary outcomes exactly as recorded; a false comparison
is preserved without changing the primary energy targets.

One Fluent aggregate flag was serialized as the JSON string "True" by
NumPy/default=str. Its raw value/type remains in the public summary. Only
that exact key is normalized after independent proofs of all 601 finite
native T/h/time/frozen arrays and all 600 source-bound step-metric records.
All other flags require literal JSON booleans. Parser/source failures block
selection; genuine numerical misses remain false.

Outlet bulk temperature uses actual outward convective-flux weighting,
not an area average. Cross-code temperature differences are descriptive;
no paired pass threshold was declared. Each solver uses its own frozen flow.
Agreement is not experimental physical validation or mesh/time independence.

OpenFOAM temperatures may be reassembled directly from retained native rank
fields using actual cell addressing. This is postprocessing of completed CFD,
not another solve or a fabricated global native file. Both rank clocks, all
rank T/U/p/face fluxes and physical boundary gradients remain source-hash bound.
Available closed serial states and the inherited pilot are crosschecked. Slow
serial reconstruction status remains separate from actual solver completion;
the record does not claim that all 600 serial reconstruction outputs exist.

The main GIF selects actual native steps 0,10,...,600; the early GIF selects
0..100. Both use 50 ms playback per frame, fixed 300–350 K scales and true
10:1 maps. There are no particles or interpolated time slices. Settled and
quantized duplicate map crops are retained. All 2,769,408 paired native
cell temperatures are independently reread for figures. The repository
contains 600 solved-step history rows per code and final maps only; full
native fields, HDF, meshes, NPZ, rank data, full-time cell exports and raw
private audit/launch provenance remain external.

Selected small OpenFOAM inputs are in
[the case record](../cases/openfoam-heated-channel-selected-v001/README.md).
Exact implementation module hashes are in
[the source index](../results/heated-channel-v001/implementation-source-index.json).
These are guarded implementation records, not a turnkey one-command package.
Their recorded CLIs require generated per-run source/parent manifests,
source-hash authorization documents, accepted pilot archives and local
installation/environment paths. Those are run artifacts to regenerate and
review, not hidden supplied dependencies. Native-data audit reproduction
requires the external native archive; small summaries cannot recreate it.
See each script's argparse interface and preparation/source guards before
attempting a fresh versioned run. The declared TARGETS document is supplied
with this record. Python dependencies include NumPy, h5py, Pillow, Matplotlib
and PyFluent; Fluent 2026 R1 and Foundation OpenFOAM 14 are external installs.
Run sequentially with two hidden processes/ranks, save before close and check
RAM, disk, native quality and the actual licence before each launch.

- [analyse_heated_channel_flow_openfoam_v002.py](../scripts/analyse_heated_channel_flow_openfoam_v002.py)
- [analyse_heated_channel_fluent_pilot_v003.py](../scripts/analyse_heated_channel_fluent_pilot_v003.py)
- [analyse_heated_channel_fluent_production_v001.py](../scripts/analyse_heated_channel_fluent_production_v001.py)
- [analyse_heated_channel_pilot_openfoam_v003.py](../scripts/analyse_heated_channel_pilot_openfoam_v003.py)
- [analyse_heated_channel_production_openfoam_v006.py](../scripts/analyse_heated_channel_production_openfoam_v006.py)
- [analyse_nozzle_fluent_v005.py](../scripts/analyse_nozzle_fluent_v005.py)
- [heated_channel_fluent_common_v001.py](../scripts/heated_channel_fluent_common_v001.py)
- [heated_channel_fluent_common_v002.py](../scripts/heated_channel_fluent_common_v002.py)
- [heated_channel_fluent_common_v003.py](../scripts/heated_channel_fluent_common_v003.py)
- [heated_channel_fluent_pilot_v001.py](../scripts/heated_channel_fluent_pilot_v001.py)
- [heated_channel_fluent_pilot_v002.py](../scripts/heated_channel_fluent_pilot_v002.py)
- [heated_channel_fluent_production_common_v001.py](../scripts/heated_channel_fluent_production_common_v001.py)
- [heated_channel_fluent_production_v001.py](../scripts/heated_channel_fluent_production_v001.py)
- [heated_channel_mesh_v001.py](../scripts/heated_channel_mesh_v001.py)
- [plot_cavity_native_v003.py](../scripts/plot_cavity_native_v003.py)
- [plot_heated_channel_native_scaffold_v001.py](../scripts/plot_heated_channel_native_scaffold_v001.py)
- [plot_heated_channel_native_v005.py](../scripts/plot_heated_channel_native_v005.py)
- [plot_nozzle_native_v004.py](../scripts/plot_nozzle_native_v004.py)
- [prepare_openfoam_heated_channel_pilot_v003.py](../scripts/prepare_openfoam_heated_channel_pilot_v003.py)
- [prepare_openfoam_heated_channel_production_v005.py](../scripts/prepare_openfoam_heated_channel_production_v005.py)
- [prepare_openfoam_heated_channel_v002.py](../scripts/prepare_openfoam_heated_channel_v002.py)

Rejected setup/audit versions remain preserved. The owner directs progression;
Codex selects, implements and operates; ANSYS/OpenFOAM compute. No independent
owner execution, reflection, research novelty, physical validation, mesh/time
independence, Git publication or university submission is claimed.
