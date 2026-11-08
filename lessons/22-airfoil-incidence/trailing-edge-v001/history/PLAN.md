# Local airfoil trailing-edge sensitivity plan

> Use superpowers:executing-plans inline; no implementeragents. One fresh whole-improvement review at the end.

**Goal:** Measure the effect of controlled local TE resolution at accepted5degreeRe100 conditions.
**Architecture:** Preserve actual baseline polygon/outside mesh; conforming red/green edge-midpoint refinement, nativequality gates, finite steady paired runs, source-bound numerical/visual comparison.
**Tech stack:** Existing E:Python/NumPy/HDF5/PyFluent, D:Fluent261, Foundation14 WSL, Matplotlib.
**Spec:** SPEC.md

## Global constraints
Exactpolygon/outerboundary, no movement ofoldnodes or originalprojects. E:tools/data/temp; D:existingANSYS. Hidden2ranks/onesolver; RAM>=4GiB/Efree>=10GiBafterimports/immediatelyprelaunch andsource/processchecksafterwaits. Nativequality>=.2 beforeflow. Unchangedresidual1e-7,mass1e-6,forceCdrelative1e-4/Clabsolute1e-5 checks; no symmetrygateat5degrees. No validation/novelty/independentowner/automaticpublication claims.

## Review focus
- Boundary midpoint subdivisions must stay on original polygon segments and preserve far-field edges exactly.
- Local split/green closure must remain conforming, positive-area and bounded inside0.15chord; unchangedoutsidecells are evidence of controlled locality.
- Nativecell-to-field maps and force area/directions must remain correct after cell renumbering/refinement.
- Steady tracer animations must not be called transientfluid simulations; adjacentcell pressure is not silently relabelled wallpressure.
- Source/resource/licence/process gates and native quality must precede flow, with all failedversions retained.

## Tasks
- [ ] 1. Beforeimplementation, write independent split/conformity/area/boundary tests and observeRED. Implement localrefiner, verifyGREEN. Parse baseline109600mesh, createversionedL1/L2 meshes with proof of originalvertex/farboundary/outsidecell/area preservation andaffectedextent. Native Fluent quality-only2process gate onL1; show actualmesh image. If qualitymiss, diagnose/newversion beforeflow.
- [ ] 2. Onlyafterqualification, finite Fluent5degree localrefinement solves with existingtemplate and after-import/immediateprelaunch guards; inspectactualfields/nativeforce/residual/massstability. MatchedOpenFOAMserially afterFacceptance. Baseline cases reused, not rerun unnecessarily.
- [ ] 3. Source-boundforces/Cp/TE mesh comparisons and actualvisuals; honestlimitations, versioned existing-studyrecord, preservation/privatefile/hashchecks; onefreshtechnicalreview withoneverifiedImportant/Criticalfixpass ifneeded. No newmodelcount orupload.

## Verification
Trianglearea/conformity/edge-tag tests, exactoldcoordinates/farboundary/polygon andoutsidecell comparisons, positive native quality beforeflow, sourcehashes/finitefields/nativecontrols and independent force normalization. Force/stability changes measured atsamemodel/domain. No universalaccuracythreshold inferred from a localmeshresponse. Existingcodexcheckout retained; no newworktree/commit/cleanup requested.
