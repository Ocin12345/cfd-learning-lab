# Earlier aircraft projects and the move to smaller verification cases

This portfolio is part of a longer, AI-assisted aircraft-CFD effort. The following history was checked against accessible project records on 1 November 2026. It is not a claim that every past project or contribution has been recovered.

## Concorde investigation — work in progress

The project explores how Concorde enabled supersonic passenger flight and the engineering and operational tradeoffs that limited wider use. A source-supported introduction was revised and compiled in October 2026. Its scope distinguishes external aerodynamics from propulsion, economics and retirement claims that require other evidence.

Earlier project notes document OpenFOAM wall-layer and mesh-quality investigations, including rejected trials. Those older notes are historical evidence, not a fresh solver audit. The current portfolio does not claim validated viscous performance, mesh/domain independence or a completed research paper for Concorde.

## Ornithopter investigation — mesh qualification in progress

The readiness record dated 16 October 2026 documents a 378,170-cell Gmsh tetrahedral mesh for a simplified fixed-configuration wing-body pilot. All eight wings were retained. The mesh passed several checks, including concavity, but failed two native checks: 6,524 low-determinant cells and 486 low-interpolation-weight faces. It also reported 1,821 severely non-orthogonal faces. A subsequent dual-mesh conversion was rejected.

No flow solver ran on that pilot mesh. This is geometry/meshing progress, not a prediction of ornithopter flight or a flapping-wing simulation.

## Controlled Fluent demonstrations — completed 1 November 2026

The cylinder and low-Reynolds-number airfoil cases in this repository establish a smaller end-to-end workflow: generate a mesh, run the actual solver, inspect residuals and mass balance, export native plots and preserve results. The cylinder additionally has a three-mesh drag comparison. These simpler results do not validate either aircraft project.

## Contribution record

The accessible conversation records show the owner choosing the aircraft topics, requesting continuation, shaping the Concorde paper question and requesting the progression to background Fluent demonstrations. They also show substantial AI implementation, execution and drafting. They do not establish which additional calculations, modelling or physics explanations the owner independently completed outside those conversations. That personal contribution should be described by the owner, with examples or artifacts, before this portfolio is used in an application.

## Evidence inspected

- Ornithopter current-readiness document dated 16 October 2026 and the corresponding chat's latest result.
- Concorde introduction source and the corresponding chat recording its revision and PDF compilation.
- Earlier Concorde mesh-qualification notes, explicitly treated as historical and potentially stale.

The original aircraft files remain in their existing project folders. This repository contains a summary of that history; it does not yet distribute the aircraft geometry, manuscript or their full provenance records.

## Learning Lab extension on E: — 1 November 2026

Sources and exact dependencies were copied and hash-verified without removing
originals. The owner requested automatic continuation after actual solver
images, later adding paired Fluent/OpenFOAM work and academic figures. Codex
implemented and operated the workflows, analysed records and drafted captions.
Fluent computed a saved-cylinder replay, a domain study through 320D and the
repaired zero-incidence airfoil. OpenFOAM Foundation 14 computed the same
98,304-cell, radius-20D cylinder as Fluent; both passed recorded numerical
checks with a 0.0806353% Cd difference. Airfoil mesh sensitivity and its matched
OpenFOAM comparison are unfinished. Owner predictions/reflections are not
invented. None of these smaller cases validates Concorde or the ornithopter.


## Airfoil model extension - 1 November 2026

The owner prioritized running models and showing results over formal lesson
checkpoints, and selected continued airfoil work. Codex implemented and operated
a polygon-preserving nonuniform mesh comparison. Fluent computed 109,600,
383,562 and 944,026-cell Re=100 zero-incidence cases, all passing their recorded
checks. Medium-to-fine Cd changed 0.00435432%, with opposing pressure/viscous
changes and small nonmonotonic lift bias. No uniform convergence order or error
bound is claimed.

The matched OpenFOAM model completed and saved 5,000 steady iterations on the
same 383,562-cell medium mesh. Cd = 0.4086185961, Cl = -0.002897496409; the absolute
Cd difference from Fluent is 1.74852%.
This is a retained diagnostic result. Failed checks: openfoam.zero_angle_lift_diagnostic.
The separate declared 0.5% cross-code target failed.
Numerical convergence status is True; the
zero-angle symmetry diagnostic status is False.

Codex checked saved data and prepared academic figures and source hashes.
The owner's direction is recorded without inventing predictions, reflections,
independent implementation, research novelty or aircraft validation.


## Paired airfoil domain extension - 1 November 2026

The owner requested the next calculation in both solvers. Codex extended only the outer domain, preserving all original medium cells. The initial annulus failed actual Fluent quality before flow; a radial replacement passed both mesh checks. ANSYS computed a saved-field continuation through 810; OpenFOAM computed its recorded cold-start/relaxation sequence. R40 Fluent Cd 0.4146461; OpenFOAM Cd 0.407369962. The cross-code difference is 1.75478%, against the unchanged 0.5% comparison target. See [the full record](AIRFOIL_DOMAIN_STUDY.md). No physical validation or independent owner implementation is claimed.

## Visual starter channel model - 1 November 2026

At the owner's request for basic student models with images first, Codex selected a small developing planar-channel example. Both solvers used the same7,680-cell mesh. OpenFOAM completed800 iterations; Fluent Student2026R1 acquired a live two-process session, converged at222 and saved after50 further iterations through272. Saved-field figures show the velocity development, streamlines, pressure and profiles. Downstream centre velocities are about1.497m/s versus the planar reference1.5m/s. Detailed checks and source hashes accompany the [channel record](../lessons/05-channel-flow/README.md).

The earlier944,026-cell airfoil continuation remains saved at7500 and is on hold after the owner's scope correction; no10000 continuation was launched. The channel does not resolve or validate that airfoil study. Native data remain separate on E:, original projects are preserved, and no publication or university submission occurred. The optional empty Fluent final residual XY export is recorded and recovered from its full native HDF residual history without rerunning the calculation. Owner reflection remains unprovided.


## Small backward-facing step - 1 November 2026

Following the owner's request for more interesting visual models, Codex built
and operated a shared 8,800-cell Re_h=100 step model. Fluent acquired a live
two-process Student session, reached its native convergence marker at 804 and
saved after 100 further iterations through 904. OpenFOAM ran two MPI ranks to
2,000 and saved the 1,900/2,000 endpoints. Sessions ran sequentially and closed
after saving. Both actual fields show reverse flow of about -0.163 m/s behind
the step. Sparse streamline figures and native-derived CSVs are preserved.

Final residual, mass and saved-endpoint checks pass. The additional Fluent
last-100 continuity maximum is 1.02929e-6, narrowly above 1e-6; its aggregate
check flag remains false. OpenFOAM's declared checks pass. No mesh independence
or physical validation is claimed. An exact old field-reader helper was
archived when graphics were improved; numerical outputs and native evidence
were unchanged. See the [model record](../lessons/06-backward-step/README.md).
The owner directed topics and presentation; Codex implemented and analysed;
Fluent/OpenFOAM computed. Owner reflection is unprovided. Original projects and
the paused fine-airfoil state are preserved. No commit or publication occurred.


## Transient cylinder visual model - 1 November 2026

The owner approved the researched vortex-shedding recommendation. Codex built
a shared 21,598-cell DFG-inspired Re=100 channel/cylinder model, checked native
mesh/inlet evidence and operated sequential hidden two-process Fluent and
OpenFOAM calculations. Both 0-2 s pilots failed the unchanged residual target,
while mass and finite fields passed. Their data and false flags were preserved.
Separate 2-4 s continuations used 60 Fluent inner iterations maximum and 8
OpenFOAM outer corrections. Final aggregate outcomes: Fluent Residual target not met (retained diagnostic);
OpenFOAM PASS. See the [full record](../lessons/07-vortex-shedding/README.md).

Figures use actual common timestamps and native fluid connectivity; signed
vorticity reconstruction and colour saturation are explicit. Original projects,
diagnostic pilots and versioned source corrections remain intact. No settled
periodic, grid/time-independent, experimental or aircraft-validation claim is
made. Owner reflection is unprovided; no commit or publication occurred.


## Finite-span 3D wing - 1 November 2026

The owner requested the next visual model. Codex selected an illustrative AR=4,
rounded NACA 0012-like wing at 8 degrees and Re_c=400. A 94,827-cell mesh was rejected
by native Fluent quality at 0.0297565 before flow. A 169,800-cell cap-refined
candidate retained two cells below 0.10. Two interior-node relocations preserved
boundary/CAD/profile geometry and raised the all-tet metric to 0.1068563745.
The selected mesh passed both native quality gates. Earlier versions remain.

Hidden two-process Fluent and two-rank OpenFOAM calculations ran sequentially
with the unchanged 1200-iteration cap. Final aggregate outcomes: Fluent MISS (retained diagnostic);
OpenFOAM PASS. Actual steady velocity/pressure fields support the paired
tracers, surface Cp and downstream figures. Animated tracer/camera motion is
labelled steady-field visualization. See [the record](../lessons/08-finite-wing/README.md).
Numerical misses remain false; no sensitivity or validation campaign, induced-
drag value, independent owner implementation or publication is claimed.
Owner reflection remains unprovided; no Git or university submission occurred.


## Paired compressible nozzle - 1 November 2026

The owner requested the next visual model. Codex implemented one 9,600-cell
inviscid ideal-gas nozzle. The selected fixed-budget states are
fluent-order-control-v004 and openfoam-order-control-v002: Fluent MISS (retained diagnostic); OpenFOAM
MISS (retained diagnostic). First-order diffusion and all missed numerical targets are retained.
Earlier higher-order diagnostics and preflow attempts remain preserved. Actual
fields support the [visual record](../lessons/09-nozzle/README.md); tracers are
fixed-snapshot visualization, not transient CFD. Full methods and source hashes
remain inspectable; no validation, owner reflection, publication or submission
is claimed. No additional CFD or cap extension was launched by the record.

## Paired lid-driven cavity - 1 November 2026

One 6,400-quad Re100 cavity was computed. Fluent retains its original 2,000-iteration state and real continuity-residual miss. A separately authorized saved-state continuation adds exactly 500 iterations without reinitialization; final selection is at 2,500. OpenFOAM remains at 2,000. Fluent PASS; OpenFOAM PASS. [Actual fields](../lessons/10-lid-driven-cavity/README.md), local face-flux checks and preserved numerical misses accompany the record. Tracers are fixed-field visualization. Owner reflection remains unprovided; no validation or publication is claimed.

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

## Heated impinging jet — 6 November 2026

[Model 26](../lessons/26-heated-impinging-jet/README.md) adds temperature transport
to each solver's fixed model-25 flow. Both codes completed 600 steps to 3 s,
with 31 matched animation frames. Residuals, saved-temperature bounds, frozen
flow, mass balance and final-step energy balance pass. Total plate heat transfer
differs by 0.0666%. Side-backflow edge effects, mesh/time sensitivity and physical
validation remain unresolved; thermal steady state is not established.

## 3D wing convergence follow-up — 7 November 2026

[Model 27](../lessons/27-wing-convergence/README.md) preserves the original
failed-residual state and a tighter-pressure-AMG test, then records a separate
300-iteration coupled Fluent restart. The new state passes the original
residual, mass and force-stability criteria. The previous OpenFOAM result was
reused after rechecking all 42 audited source hashes. Updated paired steady-field
visuals retain the unresolved drag/lift differences of 7.89%/10.71%.

## Finite-wing incidence comparison - 7 November 2026

[Model 28](../lessons/28-wing-incidence/README.md) adds paired effective-zero-incidence
runs on the preserved 169,800-cell wing mesh. OpenFOAM passes its checks.
Fluent's separately saved 300-iteration continuation ends at continuity
1.0096634273447894e-5 and retains its residual miss against 1e-5; mass and force
stability pass. Four-panel visuals compare these with model 27's 8-degree cases.
Tracer animation represents steady fields. Force differences, mesh and domain
sensitivity remain unresolved; no physical validation or upload is claimed.

## Finite-wing mesh sensitivity - 7 November 2026

[Model 29](../lessons/29-wing-mesh-sensitivity/README.md) records a shared
224,972-cell refinement at the same 8-degree/Re400 conditions. One interior
node was repaired without changing boundary coordinates or connectivity.
Both solvers saved 1,500 iterations. Fluent MISS; OpenFOAM
PASS. Native checks and misses accompany paired steady-field
visuals. Two grids do not establish mesh independence or physical validation.

## Natural-convection startup qualification - 7 November 2026

A [paired one-second cavity pilot](../results/natural-convection-qualification-v001/README.md) adds a different physical mechanism: circulation driven by hot/cold walls. The accepted Fluent v004 and OpenFOAM v002 results pass their recorded control and saved-state checks; all 50 Fluent time steps meet the strict residual limits. Earlier audit/iteration misses are retained. Full energy-equation matching, sensitivity and longer-time production remain open. This qualification is not counted as another finished model and has not been uploaded.

## Natural-convection development - 7 November 2026

[Model 30](../lessons/30-natural-convection/README.md) adds paired 20-second ideal-gas cavity runs on 2,304 cells. A separate vertical-flow control established the explicit Fluent gravitational-work source needed to correspond to OpenFOAM's energy term. All 1,000 time steps pass the declared criteria in both accepted runs. Earlier startup coupling misses and a recovered native-image export are preserved. The 24-to-48-cell-per-side comparison shows about 7.5% startup velocity change, so mesh/time independence and physical validation remain unclaimed. Figures, code and hashes are prepared locally; no new upload is claimed.

## 2026-11-07: natural-convection time-step sensitivity

[Model 31](../lessons/31-convection-time-step/README.md) halves dt to 0.01 s from the same saved initial fields. Both solvers complete 2,000 steps to 20 s and pass all declared step checks. Field, mass, EOS and heat comparisons are recorded. The t=0 restart rounding diagnostic is preserved. This comparison does not establish time-step independence or physical validation. Results are prepared locally; no new upload is claimed.

## 2026-11-07: localized floor-heater startup qualification

A new 60 x 100 mm closed cavity with a centered 15 mm heater has paired one-second results on 3,840 shared cells. Both accepted runs pass 100 steps and native field checks. Two earlier Fluent convergence misses are preserved; PISO relaxation was changed without weakening criteria. A separate pressure-treatment efficiency diagnostic remains unresolved. No developed-plume or physical-validation claim is made. See [the qualification record](../results/hot-patch-plume-qualification-v001/README.md).

## 2026-11-07: heater startup animation

Added six matched native states through 1 s, paired temperature/speed panels and numerical histories. Native physical times and cell mapping are independently checked. The 20-second production remains separate and unfinished. See [the animation record](../results/hot-patch-plume-qualification-v001/startup-motion-v001/README.md).

## 2026-11-07: localized floor-heater plume

Model 32 adds the paired 20-second localized-heater calculation; all 2,000 steps pass the stated criteria. Native fields, wall fluxes and visual provenance are checked; physical and sensitivity validation remain open. See [the model record](../lessons/32-heater-plume/README.md).

## 2026-11-07: water-tank sloshing

Model 33 records full 1.2 s paired sloshing runs, 2,400 residual-accepted steps each and 121 checked native states. The original OpenFOAM alpha-bound misses are preserved; a pressure-tolerance repeat passes the unchanged bound. Actual-state MP4/GIF, CSV histories, input dictionaries and source hashes are prepared locally. Mesh/time independence and physical validation remain open. See [the model record](../lessons/33-tank-sloshing/README.md).

## 2026-11-07: sloshing mesh sensitivity

Model 34 completes the paired 160 x 80 sloshing repeat: both 1.2 s cases pass 2,400 step checks and 121 native-state checks. Matched 2.5 mm strip comparisons show maximum fine/coarse depth changes 0.168348 mm in Fluent and 0.059781 mm in OpenFOAM. Native-state MP4/GIF and numerical CSVs are prepared locally. Two grids do not establish mesh independence or physical validation. See [the model record](../lessons/34-sloshing-mesh-sensitivity/README.md).

## Portfolio publication preparation — 8 November 2026

The owner requested upload of the complete prepared portfolio. Records through study 34 and the explicitly unfinished study-35 qualification snapshot are included. A study index links all records. Native simulations, installed environments, credentials and licence files remain outside Git. Existing numerical misses and contribution boundaries are retained. Publication success is established by the GitHub commit, not by this preparation note.
