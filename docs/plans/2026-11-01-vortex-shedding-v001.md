# Cylinder vortex-shedding implementation plan

> For agentic workers: use executing-plans or subagent-driven-development for this bounded model.

Goal: actual Re=100 cylinder shedding imagery in Fluent and OpenFOAM.
Architecture: shared small triangle mesh; separate hidden sequential two-process transient runs; native field/force postprocessing.
Tech stack: E: Python/PyFluent/Gmsh/Matplotlib, Fluent Student2026R1, OpenFOAMFoundation14.
Spec: E:/CFD-Learning-Lab/local/vortex-shedding-v001/DESIGN.md.

Constraints: preserve prior files; E: artifacts/temp, D: ANSYS; no concurrent solver sessions; actual resource/licence admission; no publication/commit. Practical outlets, transient startup and mesh/time limitations stay explicit. Resource floors4GiB Fluent/2.5GiB OF and10GiB disk. Numerical misses stayfalse.

Review focus: native inlet face ordering and positive flow; pressure-outlet vs DFG traction distinction; correct physical timestamps and shared connectivity; real per-step force sampling; visual startup versus mature periodic shedding.

- [x] Generate/inspect/hash shared mesh and parabolic inlet; check native geometry/quality.
- [x] Prepare/review and run a saved Fluent2s pilot; export actual fields/history.
- [x] Prepare/review and run a saved OpenFOAM2s pilot sequentially; check actual inlet and mesh.
- [x] Continue separate saved states to4s with stronger iterations; preserve under-converged pilots.
- [x] Produce/show actual shared-scale animation/stills and lift/drag traces; retain check outcomes.
- [x] Save methods, sources, limitations and contribution record; verify copies/hashes/links.


Delivered with an evidence-based scope change: preserved diagnostic 2 s pilots,
stronger iterative continuations to 4 s, and actual shared-window animation.
The planned 8 s extension was unnecessary for the bounded visual deliverable
and was not launched. Numerical outcomes are recorded separately; completed
implementation tasks do not mean all numerical criteria or validation passed.
