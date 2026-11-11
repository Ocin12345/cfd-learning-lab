# Passive T-junction tracer: methods

Both solvers show a dye pulse travelling down the upper branch, turning into the outlet channel and fading as unmarked fluid replaces it.

Paired numerical outcomes: ANSYS Fluent PASS; OpenFOAM PASS. Each code passed its 11 primary numerical checks.

## Model and transport

A dimensionless tracer pulse enters the upper branch of a two-dimensional T-junction; unmarked fluid enters from the left. Advection carries the tracer, while diffusion spreads concentration. Each solver uses its own previously accepted 4,400-quad frozen flow. Channel width H=0.2 m, inlet speed U=0.1 m/s, rho=1 kg/m³ and mu=0.001 Pa s give inlet Re=20. Mean outlet speed is 0.2 m/s and Re=40. Fluent depth is 1 m; OpenFOAM depth is 0.01 m with empty front/back patches.

D=0.001 m²/s and nu=0.001 m²/s give Sc=1 and inlet/outlet Pe=20/40. Fluent uses Gamma=rho D=0.001 kg/(m s); OpenFOAM uses D with volumetric phi. Scalar transport leaves velocity, pressure, density, viscosity, mesh and all native face fluxes fixed. No momentum, pressure, energy, turbulence, reaction or feedback equation is solved.

Initial fluid c=0; the left inlet stays 0. Upper-inlet c=1 at solved steps 1–80, then 0 at 81–800, with dt=0.025 s. Fluent switches after step 80; OpenFOAM's step table falls at 2.0125 s. Actual native 80/81 boundary samples and clocks were checked in each accepted audit. This is 80 nonzero endpoint samples, nominal discrete exposure 2 s, not an exact 2 s continuum-edge claim. Inlet diffusion remains enabled, including back-diffusion after switch-off. Walls have zero diffusion; the qualified outward-flow outlet has convection plus zero diffusion. Unsupported backflow prevents selection.

## Numerical method and checks

Each pulse retains five accepted pilot steps and adds 795 without initialization or replay: 800 steps, 801 native states, 20 s. Separate five-step uniform c = 1 controls are excluded from pulse ancestry. Euler and scalar upwind are first order. Fluent permits 60 inner iterations per step, own new terminal scaled UDS residual<=1e-8, excluding inherited 2000 unsolved rows. Saved M1 must match preceding current c within 1e-12. OpenFOAM uses qualified PBiCGStab/DILU, relaxation = 1, nCorrectors = 0 and one matrix solve, final native residual<=1e-8; no MULES/subcycles. Residual norms differ.

The original OpenFOAM zero-step c/C filename collision and subsequent smoothSolver residual miss remain preserved. Native scalar is now named `tracer`; new controls/pilots passed unchanged thresholds. Fluent's qualified printed UDS report matched total inward transport at printed precision; absent unit headers and depth 1 integrated/per-depth ambiguity remain explicit. Missing production reports are not invented.

Full audits check native clocks, ancestry, mesh quality/identity, source maps, source/copy hashes and frozen fields within 1e-12 absolute SI. Scalar values must be finite and within[-1e-8,1+1e-8]; extrema, undershoot/overshoot counts and roundoff remain unmodified. Courant number near 0.89 is descriptive, not time-accuracy evidence.

Per-depth inventory I=sum(rho V c)/depth has units kg/m of dimensionless mass-weighted tracer, not calibrated dye mass. Every solved step reconstructs upwind advection and orthogonal diffusion from actual c, geometry and phi. With storage S=(I_now−I_previous)/dt and boundary rates B positive inward, require abs(S−sum B)/max(abs S,sum abs B,0.02 kg/(m s))<=1e-5. Retain separate transfers, inlet back-diffusion, raw defect/denominator, cumulative inventory and flux-weighted outlet c. Local defects normalized with H/U=2 s are descriptive. All transfers/defects are RECONSTRUCTED, not saved native matrix flux. OpenFOAM uses 9040 observed active faces; 8800 empty-face placeholders are excluded. No persisted OpenFOAM old-time or incompressible density array is claimed. Initial solved budgets/residuals are unavailable.

## Saved evidence and limitations

Fluent's slow archive check was interrupted after native step 800 was saved. Original partial status, absent summary, 76 available producer budgets (6–81) and 719 missing (82–800) remain preserved. Independent auditing accepted all 801 states/all 800 budgets and 11 checks; the checkpoint is the original autosave. No physical replay or threshold change occurred.

OpenFOAM reached 20 s before its original reconstruction was interrupted (exit 130). Native reconstruction finished in a verified preserved copy (exit 0), without CFD replay. Its independent full audit passed 11 checks across 801 states/800 steps/2,400 clocks; 21,399 source inputs remained unchanged. Original failed/pending flags and the absent execution check remain separate.

The GIF selects 101 actual states 0:8:800, fixed concentration 0–1, equal physical axes and exact T-domain cells. No fabricated particles, interpolated states or clipping are used. Snapshots and inventory/outlet/budget/own-residual histories accompany it. Independent QA verified 7,048,800 full CSV rows, 8,800 final rows, 1,602 history rows, 101 decoded frames and source/output hashes. Data and other artwork remained unchanged. CSV geometry describes canonical plotting polygons; NPZ coordinates/areas retain authenticated native geometry. Mappings and native scalar values were checked independently. Raw stage-time pending flags are preserved alongside the separate completed QA receipts.

First-order smearing, one mesh/time step, illustrative properties and separate frozen flows limit accuracy. Twenty seconds is arrival/early washout. No complete washout, cross-code agreement threshold, mesh/time independence, physical/aircraft validation or novelty is established.

[Targets](../results/tracer-tjunction-v001/PRODUCTION-TARGETS.md), [code](../results/tracer-tjunction-v001/source-index.json) and [dictionaries](../cases/openfoam-tracer-tjunction-selected-v001/README.md) record the method; replay also requires the external native archive and installed environments. All project data/temp stay on E:, Fluent installation on D:. Use fresh paths, hidden sequential two-process/rank sessions, >=4 GiB available RAM/10 GiB E free, saving before closure. Large native files, dependencies and full CSV stay outside Git. The user directs. Fluent/OpenFOAM compute.
