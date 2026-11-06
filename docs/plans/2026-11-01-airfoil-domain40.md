# Paired airfoil domain-radius comparison plan

Goal: compute a fresh Re=100 zero-incidence case in ANSYS Fluent and OpenFOAM
with outer radius 40c, preserving the accepted radius-20c medium inner mesh.
This controlled extension measures boundary-placement sensitivity; it does not
assume it will explain the previous code discrepancy or validate either solver.

User scope: run models with actual images, no prediction/reflection pauses.
Keep all native data, dependencies, working files and supported caches on E:;
use the existing D: ANSYS installation. Preserve all originals and earlier runs.

- [x] Generate only an outer annulus. Independently verify all original node
  coordinates and triangle connectivities, wall patch and former outer interface;
  no original inner cell is changed. Check total area/topology/patch assignments,
  positive cells, no holes/overlaps and counts below the 950,000 admission ceiling.
- [x] Before launch, verify free E: space >=10 GiB, RAM >=4 GiB and no existing
  Fluent session. Start one hidden Fluent session with two processes. Save native
  quality/convergence/mass/force reports, case/data and actual solver images.
- [x] Convert the exact same mesh in OpenFOAM Foundation14, check all geometry and
  topology, then decompose into two ranks. Use the same R20 cold-start schedule:
  initial U relaxation0.7 to2000, then0.9 with separately recorded continuation.
  Start OpenFOAM only after resource admission; run sequentially if necessary.
- [x] Keep original numerical thresholds: p/U initial residuals <=1e-7 (OpenFOAM
  pressure is maximum across correctors), last100 Cd relative range <=1e-4,
  Cl absolute range <=1e-5, zero-angle |Cl| <=0.001, mass imbalance <=1e-6.
  Retain failed diagnostics explicitly. Fluent observes100 extra iterations.
- [x] Compare each R40 result with its own R20 baseline and compare the two codes
  at R40. Measure pressure/viscous changes where available; do not assert a cause
  for discrepancies without supporting evidence. Preserve full startup histories.
- [x] Export checked actual fields and academic charts with common colour ranges,
  vector PDF/SVG and source hashes. Render/inspect PDFs and show each solver image.
- [x] Save methods, contributions, results and limitations in the repository;
  exclude native solver files, private details and dependencies. Verify manifests,
  relevant source checks and local staged bytes. No GitHub publication this step.

Implementation: new extend_airfoil_domain.py handles annulus generation and
preservation receipts; existing airfoil_run.py and prepare_openfoam_airfoil.py
handle solves/preparation. Independent review covers source preservation and
actual saved native cases. Generalize plotting captions only where required
for the new radius and bind them to checked record values.


Ruling after actual Fluent admission: the initial Gmsh annulus in native v001
failed minimum orthogonal quality (0.0485061 vs0.2), before any flow iterations.
Its rejected native case and mesh are preserved. Use a fresh v002 radial annulus
with96 geometric layers and768 angular segments, split into triangles. All
original inner cells remain unchanged. Actual Fluent admission still governs;
mesh-only Gmsh/cosine measures are not treated as Fluent orthogonal quality.
The exact same v002 mesh will be used by both solvers if it passes.

Execution completed through the saved checked/diagnostic states in AIRFOIL_DOMAIN_STUDY.md; source/hash/visual checks are recorded separately.
