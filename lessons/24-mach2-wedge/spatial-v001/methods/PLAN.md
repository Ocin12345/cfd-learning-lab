# Mach 2 wedge spatial-sensitivity implementation plan

> Execute inline using superpowers:executing-plans. No implementer agents; one fresh whole-improvement review at completion.

Goal: Measure one controlled refinement's shock-angle/pressure-jump response and numerical transition thickness, with candid OpenFOAM adaptive-time limitations.
Architecture: Preserve accepted study24 baseline; double grid counts only, native geometry/source audits, finite sequential qualifications, source-bound measurements/visuals and versioned existing-study record.
Tech stack: existingE:Python/NumPy/HDF5/PyFluent/Matplotlib,existingD:Fluent261,WSLFoundation14.
Spec: SPEC.md

## Global constraints
Use exact physical/domain/boundary/scheme values and unchanged gates in SPEC.md. No overwrites or frozen-source edits. E: storage/temp,hidden2ranks/one solver,RAM>=4GiB/Efree>=10GiB afterimports/immediateprelaunch/afterwaits. No newmodelcount/physicalvalidation/publication claims.

## Review focus
- Refinement preserves actual boundary polygon/corner and is common to both codes; native field ordering must be bijective.
- Steady iteration is not physical time; adaptive-CFL spatial/time confounding must stay explicit.
- Threshold interpolation must be connected to the selected shock and retain missing/ambiguous columns instead of cherry picking.
- Thickness in metres/cell intervals is a numerical diagnostic, not physical shock thickness or observed order.
- Full native convergence/mass/positivity/thermodynamics/outlet/source gates before accepted visuals and a record.

## Tasks
- [ ] 1. Verify all186 old native-source hashes and current2327-file repository baseline. Read original scripts/settings and primary NASA reference. Write independent geometry and thickness-bracket tests and observeRED. Implement minimal helpers; confirmGREEN. Generate versioned89600-cell mesh, compare exact/collinear boundaries and old-node subset. Native Fluent quality-only qualification, save/close and inspect actual mesh.
- [ ] 2. Prepare guarded finite fine Fluent solve reusing the qualified mesh and original schemes/checkpoint schedule; audit all saved field/material/BC/mass/residual/stability evidence. Preserve failures. Admit finite fine OpenFOAM only after acceptance; native topology/ordering/source/resource/time/Courant/field/mass checks. Inspect main endpoint before deciding the bounded exact-checkpoint late-CFL diagnostic. No whole-shortlist launch.
- [ ] 3. Replay unchanged reference angle/pressure metrics plus tested connected10–90%thickness on actual native fields. Actual coarse/fine paired fixed-scale images/profiles/numerical plots; explicit limitations. Stage a versioned study24 record, prior-file/checksum/exclusion verification and ONE fresh technical review. No upload without authorization and remote verification.

## Interfaces and verification
Geometry helper emits wedge.msh,geometry.npz,geometry-proof.json consumed by native quality/flow drivers. Measurement helper consumes actual mesh-aligned p field and upstream/downstream plateau estimates; analytic profiles must recover known widths without smoothing or extrapolation. Source/launch manifests bind exact selected scripts and input files. Numerical audits precede image selection. Baseline same-CFL OpenFOAM comparison is explicitly combined spatial/adaptive-time response; late-CFL branch is a separate diagnostic. Preserve the existing codex checkout; no newworktree/commit/cleanup requested.
