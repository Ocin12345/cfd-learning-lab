# Model 34: sloshing mesh sensitivity implementation plan

> **For agentic workers:** Use superpowers:executing-plans for native autonomous execution. Retain the user-authorized E: checkout and all versioned outputs.

**Goal:** Compute the same 1.2 s sloshing case on a 160x80 grid in Fluent and OpenFOAM, compare with accepted model33 at matching times/physical strips, show real motion and prepare a local record.
**Architecture:** New versioned inputs and sequential hidden two-rank jobs. Short flat controls and tilted qualifications precede full production. Reuse model33's raw native data after source verification; never alter it.
**Tech stack:** Existing E: Python/PyFluent and Foundation14, D: Fluent261, native HDF5/ASCII, Matplotlib and verified FFmpeg.
**Spec:** Model33 README/selected native inputs and local/sloshing-production-v001/accepted-pair-selection-v001.json; user requests continued paired models, real images/video, preservation and autonomous operation.

## Constraints and measurement choices
- Same 200x100mm tank, stationary no-slip walls, 50mm mean depth, initial slope0.025, zero velocity, gravity(0,-9.81), constant air/water properties, laminar incompressible2D and sigma0.
- Fine grid160x80,12800 square1.25mm cells; coarse80x40 remains preserved. OpenFOAM one10mm extruded layer with empty front/back.
- dt0.0005s, production2400steps/1.2s,121 native saved states every20steps. Do not change only the fine time step if a Courant miss occurs.
- Fluent dynamic scaled continuity/Ux/Uy <=1e-4,cap400; OF final-outer initial Ux/Uy/p_rgh <=1e-4,outer8/pressure3, alpha linear1e-11,p_rgh1e-10. Interface Co<=0.25. Alpha within1e-6 of[0,1]; water-area drift<=0.5%.
- Qualify each code with40steps at the production dt, saving every step, first flat then tilted. Flat Umax<=1e-4m/s and maximum fraction change<=0.01. Fluent quiet raw residual limits remain continuity1e-12,momentum1e-10,with automatic convergence and400 cap. The previous forced400 budget was a reporting-equivalence diagnostic, not a general quiet-control requirement; no limits are relaxed.
- Verify analytic initial fractions in every native cell before advancing. Check that four fine child fractions average to each original coarse fraction within1e-10; preserve water area0.01m2.
- Compare column-integrated water depths by averaging adjacent fine columns into the same2.5mm physical strips as the coarse case. The right-most comparison strip is x197.5..200mm, center198.75mm, in both grids. Do not compare different probe locations under the same label.
- Report maximum/RMS depth and trace differences and their histories. Two grids do not establish mesh independence. Time-step independence and physical validation remain open; no aircraft claim.
- Recheck prerequisites/native hashes after memory waits and immediately before solver launch. >=4GiB available RAM,>=10GiB E: free, one solver at a time,two ranks,hidden. Project temps/caches E:,MPI IPC /dev/shm only.
- Native paths: native/runs/sloshing-refinement-v001/qualification/{fluent-flat,openfoam-flat,fluent-tilted,openfoam-tilted}; production/{fluent,openfoam}. No output reuse or deletion.

## Tasks
- [ ] 1. Prepare parameterized fine geometry, Fluent driver, OpenFOAM cases and native auditors. AST/input/nested-initial checks, recorded source hashes, real launch admission checks.
- [ ] 2. Run/audit four short qualifications sequentially; preserve failures. Show actual fine-grid fields before claiming a full result.
- [ ] 3. Run both full1.2s cases only after qualifications pass; audit2400 logged steps and121 native states each.
- [ ] 4. Render matched-time/matched-strip coarse/fine comparisons and decode-check the actual-state video. Review representative frames, numerical histories and all claims; prepare model34 locally with sources/contributions/limitations and no invented reflection.

## Review focus
Stale prerequisites after waiting; wrong inherited3200-cell helper; mismatched probe positions; changes to both mesh and physics/time; early field/alpha/Courant failures hidden by final-state-only checks; confusing numerical agreement with independence or validation.
