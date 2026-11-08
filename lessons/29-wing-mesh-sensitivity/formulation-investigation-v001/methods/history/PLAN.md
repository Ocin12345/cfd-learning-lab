# Refined-wing convergence repair implementation plan

> For agentic workers: use superpowers:executing-plans inline. One fresh whole-improvement review only at the end.

Goal: qualify a numerically converged endpoint on the unchanged refined wing, confirm iterative stability and improve existing study29's paired evidence.
Architecture: immutable source/native inventory; tested continuation/settings/monitor helpers; two sequential same-checkpoint100iteration tests; actual accepted visuals; only then finite100iteration confirmation; read-only paired record using accepted OpenFOAM.
Tech stack: existing E:Python/PyFluent/h5py/numpy/matplotlib/tetrahedral visualization code,D:Fluent261, existing Foundation14 native data (no OFsolve).
Spec: SPEC.md. All physical/resource/preservation/numerical constraints there apply unchanged.

Review focus: inner solver control versus outer acceptance; exact1500restart/globaliterationprefix; monitor-output path redirection; force/reference/BC equivalence without accuracy ranking; no earlier-checkpoint selection or inherited miss relabelling.

- [x] Task1 native/source inventory. Actualnative refinedF1500MISS andacceptedO1500, commonTetGeometry, materials/BCs/forcearea/sourceguards verified.2798repo baseline and70explicitpersistentbindings captured. No private/runtime metadata opened/copied.
- [ ] Task2 tests and helpers. RED thenGREEN for strict per-iteration force parsing with only an identical restart echo, full final100force window, native residual coverage/prefix/scaling checks, and temporarily active residual target getter with restoration onfailure. Reuse preserved native readers; no solver.
- [ ] Task3 qualification pair. Fresh exact1500copies; tested loaded-settings/geometry/field verification, output-path guards, only declaredAMGcriterion difference. Freeze source/parent/input manifest. Hidden sequential2rank100baseline then100tight, checkpoint every20, savebeforeclose; closed native audits distinguishstructuralfailure/numericalMISS/PASS. Parentprocess/resource rechecksafterwaits. No OFsolve.
- [ ] Task4 actual pilot visuals/selection then finite confirmation. Sourcebind accepted case/nativefields vs unchangedO, inspect scaledresidual/forces and actual pressureview. Only then new exact1600checkpoint100iterationconfirmation to1700 with ALLnew residuals passing. Preserveallhistory and numericalmisses, no silent cap extension.
- [ ] Task5 existingstudy29record with real fixedscale paired fields/nativeforces/history and optional honestly labelled steady-fieldmedia. Rehash matchedacceptedO, recordsourcebindings/methods/results/limits/contributions, preserve2798baseline/privateexclusions/manifest. No new model/upload/physicalvalidation.
- [ ] Task6 ONEfreshwhole-improvementreview, oneCritical/Importantfixpass ifrequired,Minorledger; finalsource/preservation/exclusion/processchecks andcompletionproof. No repeatedreviews or cleanup of sourceworkspace.

Interfaces: Task2 reads exactnativeglobaliteration andforce data forTask3; parser mustreject missing/conflictingduplicate rows rather than hiding them. Task3 endpoint audit accepts onlydeclared1600; Task4 confirms onlyselected1600 for exactly100more. Task5 derives images/data from explicit selectedcase/audit and untouchedO hashes. No task may infer PASS from solver process exit alone.
