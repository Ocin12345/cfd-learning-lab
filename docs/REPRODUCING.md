# Reproduction and storage

Requires a licensed Fluent installation and Python with
`ansys-fluent-core==0.43.0`. The original runs used Student 2026 R1.

## Local E: setup

Project root: `E:\CFD-Learning-Lab`. The repository is in `repo/cfd-learning-lab`;
native source imports and new runs are in `native`; Python is in `environments`;
project temporary files/caches are in `temp`; large exports are in `exports`.
ANSYS remains installed on D:. The local private copy receipt gives original
source paths and per-file SHA-256 checks. The C: originals were copied, not moved.

In a fresh PowerShell session:

```powershell
. E:\CFD-Learning-Lab\local\Use-Project.ps1
Set-Location E:\CFD-Learning-Lab\repo\cfd-learning-lab
python -c "import ansys.fluent.core as pf; print(pf.__version__)"
```

This selects the Python environment, the D: Fluent executable and a fresh
timestamped `CFD_OUTPUT_ROOT` outside the repository. It sets TEMP, TMP, TMPDIR,
Python bytecode, pip/uv and Matplotlib cache locations for this shell and child
processes. It does not change global machine settings or guarantee that ANSYS
installation-level caches follow those variables. The owner now authorizes
automatic continuation after each actual solver image. Predictions are optional
and must not be invented. Resource/licence checks still precede each launch.
Python bytecode writes are disabled for this project after slow atomic cache
writes were observed on E:; other supported caches remain under `temp`.

## Portable workflow

In another installation, install requirements in your own environment and set
`FLUENT_EXECUTABLE` and `CFD_OUTPUT_ROOT` explicitly. Choose a new output directory
outside the checkout. Original copied workflow scripts now require
`CFD_OUTPUT_ROOT`; they do not default to writing inside the repository.

The original reproduction sequence, for use only when its lesson is authorized:

```powershell
python scripts/cylinder_mesh.py --name ansys-cylinder-refinement/medium --nt 384 --nr 144
python scripts/cylinder_mesh.py --name ansys-cylinder-refinement/fine --nt 512 --nr 192
python scripts/cylinder_background.py
```

For the coarse cylinder use `--nt 256 --nr 96` in a separate run root's medium
directory. Airfoil generation uses `scripts/airfoil_mesh.py` with
`--name ansys-airfoil/alpha0 --alpha 0` and `--name ansys-airfoil/alpha5 --alpha 5`,
followed by `scripts/airfoil_background.py`. These are reference workflows, not
an instruction to run all cases now.

Use one hidden Fluent session at a time, initially with two solver processes;
concurrent launch previously encountered a licence error. Check disk space,
RAM and existing sessions immediately before launching. Save case/data before
closing. Each new investigation uses a new versioned output root.

## What is and is not reproduced

The original meshes and saved case/data pairs are preserved in the external
native imports. The original portfolio is also preserved there as an immutable
import snapshot. The imported baseline scripts were adapted for external output storage and
syntax-checked. Later lesson scripts have separate recorded execution evidence;
consult each lesson for exactly which meshes and calculations were run. They can save
results after an iteration cap; inspect convergence and numerical reports.

Force histories were not configured in those original scripts. Lesson 1 added a
separate saved-case replay and 100 further iterations with force sampling.
The initial licence query returned no feature list, but later hidden two-process
Fluent launches succeeded for the replay and the radius-40 case. Thus live
Student use has now been demonstrated; repeat checks before subsequent launches.

## Matched OpenFOAM and Fluent case

Lesson 4 records the completed paired cylinder pilot. OpenFOAM Foundation 14
runs in the `CFD-Foundation14` WSL distribution; the Windows Fluent installation
remains on D:. Native cases live outside Git on E:. Source the installed
OpenFOAM environment before conversion/checking/running, and direct supported
temporary files to the project's E: temp directory. The recorded dictionaries
are in `cases/openfoam-cylinder-re20-recorded`. The initial and extended control
dictionaries are separate records, not interchangeable startup instructions.

Use `scripts/prepare_openfoam_cylinder.py --help` to prepare a fresh external
case. Convert its mesh with `fluentMeshToFoam -2D 1 input-cylinder.msh`, confirm
patch names/types and pass `checkMesh -allGeometry -allTopology`. Decompose once,
run two MPI ranks, and reconstruct saved fields. Preserve completed states
when continuing from 1,200 to 3,000 iterations; do not restart from zero or
re-decompose over existing processor states. Lesson 4 explains the continuation.

Run `scripts/analyse_cross_solver_cylinder.py --help` for independent numerical
checks on a completed native run. It requires native evidence, including saved
fluxes, settings and force reports; the small Git tables alone are insufficient
to repeat every check. Plotting scripts can regenerate figures from the public
tables into a fresh external directory. Residual definitions and starting
fields differ between the codes, so the pilot is not a speed comparison.


## Airfoil mesh comparison completed on 1 November 2026

The accepted zero-incidence mesh family uses the imported original closed
polygon, not a new analytic outline or a rounded trailing edge. Use
`scripts/airfoil_triangle_mesh.py --help` with `--local-edge-sizing` and
`--scale-edge-limits`, at scales 1.4, 1.0 and 0.55, into three fresh external
mesh directories. At scale 1 the additional edge-scaling option is identical
to the earlier saved medium settings. Actual results were 109,600, 383,562 and
944,026 triangular cells. Retain `source-gmsh.msh`, `airfoil.msh`, `setup.json`
and the executed mesher snapshot.

For each fresh flow run, `scripts/airfoil_run.py` takes the original zero-angle
case/data pair as `--case`, the generated `airfoil.msh` as `--mesh`, the recorded
loaded-settings JSON as `--template-settings`, and a fresh external `--out`.
The mesh-only option exits after Fluent quality checks; omit it for an intended
flow calculation. This comparison used one hidden double-precision Fluent
session at a time, with two processes. The script rejects minimum orthogonal
quality below 0.2, solves the fresh mesh and samples another 100 iterations.
Check its completion status and independent analysis rather than inferring
success from process exit alone.

Use `scripts/airfoil_mesh_metrics.py --help` to measure source-polygon/segment
preservation, patch identities, wall-adjacent normal heights and regional
spacing. Then use `analyse_airfoil.py` for each native saved run and
`analyse_airfoil_sensitivity.py` to aggregate the three checked records.
The aggregate requires matching methods and geometry, increasing cell and wall
face counts, decreasing regional spacings and unchanged source hashes. This
nonuniform sequence does not establish a uniform refinement ratio or GCI.

The final native family is `native/runs/2026-11-01-airfoil-sensitivity-v002`,
with `coarse-solve` and `fine055-solve`; the reused medium is
`native/runs/2026-11-01-lesson03-airfoil-mesh-v002/alpha0-local-quality`.
That medium directory contains a full flow solve despite its historical name;
its original launch-intent deviation remains recorded. Earlier unsolved mesh
trials are retained. See [airfoil status](AIRFOIL_MODEL_STATUS.md) for results.

`plot_airfoil_sensitivity.py` uses the checked tables. `plot_airfoil_fields.py`
uses actual saved Fluent HDF5 face/cell connectivity and fields, without assuming
Gmsh cell order. Native field files remain external. Figure directories must be
fresh. Default common ranges are speed/U_inf 0 to 1.2 and Cp -2 to 2; any revised
common Cp limits must contain all values and be recorded in both solver figures.
Render and inspect vector PDFs as well as the PNGs.


## Matched airfoil OpenFOAM execution

Use `prepare_openfoam_airfoil.py --help` to prepare a fresh external case from
the medium mesh. In Foundation 14, convert with
`fluentMeshToFoam -2D 1 input-airfoil.msh`, then pass
`checkMesh -allGeometry -allTopology`. Decompose once into two MPI ranks, run
`mpirun -np 2 foamRun -parallel` and reconstruct saved fields. Direct supported
temporary files/caches to the E: project. This steady iteration index is not
physical transient time.

The actual uniform-field run saved 2,000, continued to 2,500 with U relaxation
0.7 to 0.9, then to 5,000 at 0.9. Recorded dictionaries preserve each stage and
initial declarations. Preserve saved processor states; do not reinitialize or
re-decompose over them. Native wrappers, checkpoints and raw logs remain on E:.

Use `analyse_openfoam_airfoil.py --case <native-case> --fluent-run
<native-reference> --out <fresh-analysis>`. It requires completed native
execution and evidence. Original targets are: maximum initial p/U residuals
<=1e-7, last-100 relative Cd range <=1e-4, absolute Cl range <=1e-5, maximum
absolute zero-angle Cl <=0.001 and relative boundary mass imbalance <=1e-6.
Pressure includes every corrector. The separate cross-code Cd target is 0.5%.

An explicit `--retain-failed-checks` saves failed numerical/symmetry criteria
with false acceptance and named failures; structural errors or missing evidence
still raise. Field/history plotters require `--diagnostic` for failed records,
visibly mark the failure and retain original thresholds. Default analysis and
plotting still refuse failed comparisons. Use the current airfoil status and
machine-readable record for actual pass/fail outcomes.


## Controlled R40 airfoil extension

See [the domain study](AIRFOIL_DOMAIN_STUDY.md) for exact source preservation, radial annulus parameters, Fluent saved-field continuation, matching OpenFOAM stages, diagnostics and actual reproduction commands. All native states and earlier failed attempts are retained externally on E:.

## Small visual channel model

See lessons/05-channel-flow/README.md for conditions and measured limitations. Generate a fresh external mesh with scripts/channel_mesh.py --out <native-run>/mesh, then channel_fluent.py --mesh <mesh>/channel.msh --out <native-run>/fluent. Source the local E: project environment first. Check actual resources and idle sessions immediately before hidden launch; the small Fluent case used4GiB available RAM/10GiB E: admission. The small OpenFOAM case used a separately documented2.5GiB RAM budget; the paused fine-airfoil6GiB requirement was unchanged.

Prepare OpenFOAM with prepare_openfoam_channel.py --mesh <mesh>/channel.msh --case <native-run>/openfoam. Source Foundation14, convert with fluentMeshToFoam -2D 1, pass full checkMesh, decompose into two ranks, run foamRun -parallel, reconstruct750:800 and use foamPostProcess -func writeCellCentres -latestTime. The exact dictionaries are cases/openfoam-channel-recorded-v001. A native solver iteration index is not physical transient time.

Use plot_channel_flow.py --run <native-run> --out <fresh-export-dir> for paired PNG/SVG fields and profiles; --solver openfoam or fluent provides a clearly single-code view. The analyser scripts read native evidence and never start CFD. Native files are required for full replay; small public CSV/JSON files remain inspectable without ANSYS. Fluent's final XY export was empty in this run; complete residuals were recovered from its saved HDF history, with the original executed script preserved in results/channel-flow-v001/executed-fluent.py. Public scripts include that recovery. All original/native files remain unchanged. No advanced sensitivity campaign was added to this introductory model.


## Small backward-facing-step model

Use fresh versioned directories under E:/CFD-Learning-Lab; do not overwrite the
recorded run. Source local/Use-Project.ps1 and use the E: Python environment.
`scripts/step_mesh.py --out <fresh-native-run>/mesh` generates step.msh and
setup.json. `scripts/step_fluent.py --mesh <mesh>/step.msh --out
<fresh-native-run>/fluent` launches the hidden double-precision two-process
Fluent case, saves its convergence endpoint and another 100 iterations, and
closes. It checks 4 GiB available host RAM and 10 GiB free disk space; also
verify no existing Fluent session and the live Student licence before flow.
The separately budgeted OpenFOAM admission was 2.5 GiB host RAM and 10 GiB disk.
The earlier fine-airfoil 6 GiB limit is unchanged; that continuation is paused.

Prepare the external OpenFOAM directory using
`scripts/prepare_openfoam_step.py --mesh <mesh>/step.msh --case
<fresh-native-run>/openfoam`. In WSL Foundation 14, source the OpenFOAM bashrc,
work in that case and run `fluentMeshToFoam -2D 1 input-step.msh`, then
`checkMesh -allGeometry -allTopology`. After Mesh OK, run `decomposePar`,
`mpirun -np 2 foamRun -parallel`, `reconstructPar -time 1900:2000` and
`foamPostProcess -func writeCellCentres -latestTime`. Use one solver at a time,
with supported temp/cache paths on E:. The exact recorded dictionaries are in
`cases/openfoam-backward-step-recorded-v001`. These steady iteration indices
are not physical transient time.

Each analyser uses `--run <native-run> --out <fresh-analysis-dir>`:
`analyse_step_fluent.py` and `analyse_step_openfoam.py`. They read native
evidence and do not launch solvers. A Fluent exit code 2 here means the retained
additional history diagnostic failed; read the JSON rather than suppressing
the flag. Final residuals must be <=1e-6, relative boundary mass imbalance
<=1e-5, and the 100-iteration saved-field maximum component change/Uin and
section pressure-drop relative change each <=0.001. The additional last-100
residual maximum is also tested at 1e-6; the recorded Fluent result narrowly
fails it. Native execution status files precede analysis; public check JSONs
provide the completed postprocessing outcome.

`scripts/plot_step_flow.py --run <native-run> --out <fresh-export-dir>` creates
paired PNG/SVG fields, CSVs, mapping checks and source hashes. It reconstructs
native cell centres instead of assuming import order and masks the solid step.
The first-cell zero-u crossing is only a near-wall proxy. The exact initial
helper used by the Fluent analyser is archived as
`results/backward-step-v001/executed-initial-field-reader.py`; the refined
figure source is `executed-plotter.py`. Methods, limitations, material/reference
viscosity distinction and residual-export recovery are in the model record.


## Transient cylinder vortex shedding

Use E: for fresh versioned native/output folders and supported caches/temp,
source local/Use-Project.ps1 and use the existing D: ANSYS installation. Check
resources, actual licence and idle sessions immediately before launch. This
model used 4 GiB available host RAM for Fluent, 2.5 GiB for OpenFOAM, and 10 GiB
free E: space. Close the saved, verified-idle dedicated WSL session when needed
for memory; do not stop other work or restart an idle VM just to probe it.

Generate the mesh with scripts/vortex_mesh.py --out <fresh-native>/mesh.
Run scripts/vortex_fluent.py --mesh <mesh>/vortex.msh --out <fresh-native>/fluent
--end-time 2 for the recorded diagnostic pilot. The stronger continuation uses
vortex_fluent_v002.py with --end-time 4 --resume-case
<pilot>/vortex-final.cas.h5 --out <fresh-native>/fluent-corrected-v002.
It preserves all pilot snapshots/history and changes the inner-iteration cap
20 to 60 while retaining the 1e-5 residual target. Always inspect actual checks;
raising a cap alone does not establish convergence.

prepare_openfoam_vortex.py prepares the exact recorded pilot. Convert with
fluentMeshToFoam -2D 1, pass checkMesh -allGeometry -allTopology, then apply the
inlet with --case <case> --apply-inlet; face order is read from native polyMesh.
Foundation 14 places this converted mesh in 0/polyMesh. Run two MPI ranks and
reconstruct saved fields. For the preserved 4 s continuation, a hash-verified
copy becomes openfoam-extended. prepare_openfoam_vortex_v003.py --case <clone>
--prepare-continuation --authorization <root-authorization.json> requires
approved=true, approved_by=root, action=continue_openfoam_vortex_to_4s. It
archives the pilot dictionaries and sets endTime4, outer8, dt0.005, snapshots
every10 steps (0.05 s), forces every step. No numerical target is relaxed.
The portable recorded dictionaries are under cases/openfoam-vortex-*-recorded*;
the source scripts retain the exact E: case guards used for this local run.
For a fresh replay, preserve those scripts and make a new version with the
fresh case paths explicitly allowed; record that path-only change. Exact
executed OpenFOAM wrappers are archived in results/vortex-shedding-v001.

Run analyse_vortex_fluent_v003.py --run <fluent-folder> --out <fresh-analysis>
and analyse_vortex_openfoam_v002.py --run <foam-folder> --out <fresh-analysis>.
These do not start solvers. Exit2 retains failed acceptance flags; exit0 means
the declared bounded checks passed, not physical validation. Native residual
definitions differ. Earlier startup failures remain diagnostic; transient force
oscillation is not required to vanish.

plot_vortex_shedding_v002.py --run <native-root> --fluent-dir
<native-root>/fluent-corrected-v002 --openfoam-dir <native-root>/openfoam-extended
--out <fresh-figure-export> --start-time 2.05 --min-end-time 4
--vorticity-limit 60 --no-mp4 produces the paired corrected-window GIF/PNG/SVG
record. It matches actual physical times without field interpolation or phase
shifting and leaves full numeric values unchanged. Large per-frame CSV exports
remain external to Git. The record distinguishes display saturation, original
under-converged startup, practical outlet differences and untested mesh/time
sensitivity. No additional model or accuracy campaign is launched automatically.


## Finite-span wing

The completed selected source is native/runs/2026-11-01-finite-wing-v003/mesh-export-v002,
not the preserved partial mesh folder. All persistent/native/temp/cache outputs
stay on E:; use the existing E: Python environment and D: Fluent installation.
Preserve recorded runs. A replay needs fresh versioned paths and an explicitly
reviewed path-only copy of the guarded scripts; do not overwrite these cases.

The native selected source has wing.msh, source-gmsh.msh, geometry.npz, setup.json,
profile.json, executed mesher/helper sources and a two-interior-node relocation
witness. Selected source hashes, candidate checks and native quality checks are
separate evidence. The rounded capped shape is not exact NACA 0012.

wing_fluent_v001.py --mesh <selected>/wing.msh --out <fresh-Fluent> is the
recorded hidden two-process 3D steady driver. Its fixed cap includes up to 100
observation iterations after native convergence; it never exceeds 1200. It saves
native cases/data, all-component residuals, native per-iteration Cd/Cl, mass and
settings. analyse_finite_wing_fluent_v002.py --run <saved-Fluent> --out <fresh-audit>
reads native evidence only and retains missed criteria at exit 2.

prepare_openfoam_finite_wing_v003.py --mesh <selected>/wing.msh --case
<native-v003>/openfoam prepares the recorded eight dictionaries. The current
local guard accepts only the explicitly selected mesh-export-v002 source and
version 003 case. The guarded private launcher runs a separate PrepareOnly
conversion/checkMesh phase, then later a separately authorized Run phase.
Both require 4 GiB available RAM, 10 GiB free E: space and no Fluent session.
Authorization action is run_finite_wing_openfoam_v003; approved=true,
approved_by=root and allowed_mode must exactly match PrepareOnly or Run.
fluentMeshToFoam imports full 3D without -2D; full checkMesh, actual patch counts,
3D bounds, positive volumes, nonorthogonality<=75 and skewness<=4 precede flow.

Foundation 14 uses incompressibleFluid/steadyState/SIMPLE, two MPI ranks, fixed
1200 iterations, saved fields 300/600/900/1200 and force coefficients each
iteration. Actual dictionaries remain in cases/openfoam-finite-wing-prepared-v003;
their prepared label is retained, with executed sources and completed audit in
results/finite-wing-v003. analyse_finite_wing_openfoam_v003.py --run <saved-case>
--out <fresh-audit> does not launch CFD. Missed residual/mass/force targets remain
false at exit 2. Native OpenFOAM p is kinematic; Pa=rho*p with rho 1 here.

The final export's executed-plotter.py and captions record the actual invocation,
native tet mappings, interpolation, Cp sampling, downstream gradients and shared
scales. The 48-frame/50 ms GIF depicts tracer integration and camera movement
through saved steady fields, not transient solver evolution. Large native/NPZ
and trace/per-frame data stay outside Git. Numerical checks do not establish
mesh independence, unsteady stability, benchmark accuracy or physical validation.
No further run or accuracy campaign is launched by this record.


## Compressible nozzle

See [full nozzle methods](NOZZLE_METHODS.md) for exact selected native folders,
executed scripts, gas/readbacks, pressure/depth conventions, boundaries, schemes,
unchanged caps/targets, numerical misses and reproduction precautions. Selected
small numerical records are results/nozzle-v001; dictionaries are
cases/openfoam-nozzle-selected-v001. Root final-source selection is private and
explicit; no guessed latest folder is accepted. Large native/NPZ/trace exports
stay external. Use fresh guarded E: paths, the existing E: Python and D: Fluent,
sequential hidden two-process/rank sessions and resource/licence admission.
Do not treat quasi1D context or a fixed-field animation as physical validation.

## Lid-driven cavity

See [cavity methods](CAVITY_METHODS.md) for exact settings, selected source/audit hashes, depth/pressure gauge conventions, local native face-flux divergence and the fixed budgets. Use fresh guarded E: paths and the existing environment. Large native/NPZ/trace files remain external. The record helpers launch no CFD, WSL or Git.

## Paired T-junction - 1 November 2026

The [actual Re20 inlet merge](../lessons/11-t-junction/README.md) records Fluent PASS; OpenFOAM PASS, unchanged 2,000-iteration caps and all numerical misses. Seed colours identify entry only; no scalar mixing or transient evolution is simulated. No validation, reflection or publication claim is added.

## Paired planar channel bend - 1 November 2026

The [actual Re20 channel bend](../lessons/12-channel-bend/README.md) records Fluent PASS; OpenFOAM PASS, unchanged 2,000-iteration caps and all numerical misses. Tracers visualize fixed fields. No scalar mixing, transient evolution or 3D secondary flow is simulated; no Dean-vortex, validation, reflection or publication claim is added.

## Paired solid-slab transient conduction - 1 November 2026

The [actual thermal slab](../lessons/13-solid-conduction/README.md) records Fluent PASS; OpenFOAM PASS, fixed Euler2s/100steps/101savedstates. Its GIF uses actual saved temperatures at0:4:200s. Global energy checks and final-only analytical error targets retain every miss and all-time errors. The initial boundary jump and coarse early time resolution remain limits. No flow, experimental physical validation, mesh/time independence, reflection or publication claim is added.

## Paired confined square-obstacle wake - 2 November 2026

The [actual Re = 20 square wake](../lessons/14-square-obstacle-wake/README.md) records Fluent own checks PASS; OpenFOAM own checks PASS, fixed 2,000 iteration caps. The separate paired last-200 mean Cd target is PASS, actual difference 0.154048% against 1%. Native force histories at iterations 1..2000 and all misses remain explicit. Fixed-field tracers do not simulate transient shedding. Height 8D / 12.5% blockage limits interpretation; no mesh/domain independence, benchmark/physical validation, owner reflection or publication claim is added.

## Heated channel — 2 November 2026

The [actual heated channel](../lessons/15-heated-channel/README.md) uses 2,304 cells, a separately computed 2,000-iteration frozen flow in each solver and 600 actual thermal steps of 0.005 s. Own primary checks: Fluent PASS; OpenFOAM PASS. Two auxiliary OpenFOAM gradient diagnostics remain unresolved. Main/early animations use actual native temperature states. Initial unsolved heat is unavailable; native Euler and reconstructed PDE defects stay separate. Raw audit corrections, auxiliary false diagnostics and all primary misses remain explicit. No paired pass threshold, transient momentum, physical validation, mesh/time independence or owner reflection is claimed.

## Annular Couette spin-up — 2026-11-02

The [actual paired annulus](../lessons/16-annular-couette/README.md) records 400 actual momentum steps to 20 s on 1,536 cells at Re = 20. Motion spreads outward from the spinning inner wall; both solvers settle to very similar circular shear profiles. Own primary checks: fluent PASS; openfoam PASS. The original OpenFOAM 10-outer pressure-coupling miss and explicit 30-outer-only budget amendment remain recorded. Its separate derived-copy preservation diagnostic and raw wrapper postcheck remain false, despite actual CFD completion. Actual native frames retain genuine plateaus; no physical validation, mesh/time independence, torque balance or 3D Taylor-vortex claim is added.

## Offset tandem-square wake — 2026-11-02

The [actual paired wake](../lessons/17-tandem-offset-wake/README.md) records 800 actual momentum steps to 40 s on 9,472 shared cells at Re = 100. The first square produces a low-speed wake that reaches the offset second square; the final native maps and profiles show interacting downstream velocity deficits. Finite-time lift histories vary during the run; sustained periodicity is not established. Own primary checks: ANSYS Fluent PASS; OpenFOAM PASS. Separate native body-force histories, prescribed unsolved initialization and exact derived-copy diagnostics are retained. No periodicity, physical validation or mesh/time independence is claimed.

## Rotating-cylinder startup — 2026-11-02

The [actual paired startup](../lessons/18-magnus-cylinder/README.md) records 500 actual momentum steps to 10 s on 8,192 shared cells at Re = 20, with clockwise wall rotation and nominal spin ratio 1. Both codes show accelerated upper-side flow and a slower asymmetric lower wake at 10 s. Own primary checks: ANSYS Fluent PASS; OpenFOAM PASS. Native forces, actual surface pressure and separate artifact diagnostics are retained. No prescribed lift sign, ideal-flow accuracy, periodicity, physical validation or mesh/time independence is claimed.

## Passive tracer through a T-junction — 2026-11-02

The [actual paired pulse](../lessons/19-tracer-tjunction/README.md) records 800 scalar steps to 20 s on 4,400 shared cells, using each solver's own frozen Re20 inlet flow. Both solvers show a dye pulse travelling down the upper branch, turning into the outlet channel and fading as unmarked fluid replaces it. Own primary checks: ANSYS Fluent PASS; OpenFOAM PASS. The GIF uses 101 actual concentration fields. Actual80/81 inlet samples, reconstructed scalar conservation, distinct native residuals and preserved failed controls accompany the record. No complete washout, physical validation or mesh/time independence is claimed; owner reflection remains unprovided.

## Natural-convection qualification

The [startup record](../results/natural-convection-qualification-v001/README.md) contains executed scripts, actual OpenFOAM dictionaries and native-source hashes. Accepted fields are Fluent qualification v004 and OpenFOAM qualification v002/heated. See the preserved corrections before rerunning; historical scripts include known failed audit/iteration versions. Native data remains on E: outside the repository.

## Natural convection at 20 seconds

Use [model 30's record](../lessons/30-natural-convection/README.md) for the corrected energy formulation, actual accepted dictionaries, executed code and native paths. The original one-second qualification remains historical evidence. Accepted fine results are natural-convection-cavity-48-v001/fluent and openfoam-v002; the latter uses 16 outer correctors. The next study must preserve these outputs and vary the time step in a new directory.

## Natural-convection time-step comparison

See [model 31](../lessons/31-convection-time-step/README.md) for the same-initial-field restart, exact executed scripts, actual dictionaries and source hashes. Accepted output is natural-convection-dt001-v001/fluent-v002 and openfoam; the first Fluent attempt stopped before advancing because of material serialization roundoff. Reproduce into fresh versioned directories.

## Localized floor-heater startup

The [heater qualification](../results/hot-patch-plume-qualification-v001/README.md) lists actual accepted native paths, source hashes, numerical changes and small dictionaries. The accepted Fluent run is fluent-v003 with PISO pressure/momentum relaxation 1.0; earlier attempts remain historical. A floating-pressure flag test is not an accepted continuation. Preserve the one-second references and qualify a longer-run method before producing the developed-plume animation.

The heater [startup animation record](../results/hot-patch-plume-qualification-v001/startup-motion-v001/README.md) includes the executed renderer and instructions. Its --qualification flag reads the completed 1 s references. Omitting that flag requires both 20 s runs to have completed and passed; no partial production is presented as a full animation.

## 2026-11-07: localized floor-heater plume

Model 32 records the unchanged initial fields, stricter Fluent continuity criterion, executed scripts and actual small OpenFOAM dictionaries. Retain the native startup qualification as the restart source and use fresh output paths. See [the model record](../lessons/32-heater-plume/README.md).

## 2026-11-07: water-tank sloshing

Model 33 includes exact executed sources, small accepted OpenFOAM inputs and all-step/saved-state auditors. Use fresh versioned paths, the recorded local qualification archives and the separate native outputs. The code is not a complete public native restart archive. Accepted OpenFOAM source is openfoam-v003, with alpha/p_rgh linear tolerances 1e-11/1e-10; earlier bound misses are retained. See [the model record](../lessons/33-tank-sloshing/README.md).

## 2026-11-07: sloshing mesh sensitivity

Model 34 retains exact fine-grid source versions, small OpenFOAM dictionaries, four short qualification audits, full production audits and matched-strip CSVs. Original native and coarse reference archives remain outside Git. Use fresh versioned outputs and regenerated admission receipts; the archived original pipeline refuses output reuse. See the model reproduction guide for its local archive dependencies. See [the model record](../lessons/34-sloshing-mesh-sensitivity/README.md).

## Study-35 qualification snapshot and publication scope

The study-35 page contains four passing 200-step qualification records, their actual images, executed sources and small inputs. Its longer production pair is unfinished in this snapshot. Source hash receipts refer to preserved E: originals and require the recorded local dependency layout. This publication is a methods/results archive, not a complete native solver restart bundle. Historical statements that a record was prepared locally or had not yet been uploaded describe its preparation stage; consult the current GitHub commit for publication status.
