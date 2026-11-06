# Next investigations

The current nozzle model is recorded with actual paired fields and retained
numerical limitations in lessons/09-nozzle/README.md. Only one physical model
was developed in this step; its numerical order controls remain separate records.

1. A compact airfoil angle comparison at stated conditions, retaining prior
   cross-code discrepancies and open sensitivity checks.
2. A prescribed pitching/heaving 2D airfoil before full ornithopter motion.
3. Matching benchmarks and controlled sensitivity only as a separately chosen
   investigation, with unchanged acceptance targets and modest resource budgets.

These are recommendations, not launched tasks. Fine-airfoil continuation remains
on hold. Concorde and ornithopter remain work in progress. Preserve native and
failed records, show actual solver-derived figures, leave owner reflections
unprovided, and use sequential hidden E: runs. No automatic cap escalation,
publication or university submission occurs.

The Re100 cavity is now recorded in [model 10](../lessons/10-lid-driven-cavity/README.md) with actual numerical outcomes. Any benchmark or mesh study remains a separately chosen investigation; no new run is authorized by this note.

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
