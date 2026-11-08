# Sloshing time-step comparison implementation plan

> For agentic workers: Use superpowers:executing-plans inline, task by task. No implementation subagents. ONE fresh whole-improvement reviewer at the end.

Goal: quantify fixed-fine-grid temporal response in paired Fluent/OpenFOAM with actual native-state media.
Architecture: source-bound exactinitialcase copies, finite80stepqualifications, conditional4720stepcheckpointcontinuations, read-onlymatchedtimefieldcomparison and versionedexistingstudyrecord.
Techstack: existingE:PyFluent/Python/HDF5, Foundation14WSL, Matplotlib/verifiedFFmpeg;D:Fluent261.
Spec: SPEC.md; parent ../PLAN.md. Userstandingautonomyoverridesdesign/predictionapprovalcheckpoints. Preserve existingcodex/initial-learning-lab checkout;no merge/push/newworktree/cleanup.

Globalconstraints: exactsame12800mesh/initialfields/physics/numerics exceptdt.0005->.00025 andoutput/end controls; alloriginalnumericalcriteria; hidden2ranks/ONEsolver,sequentialF/O; RAM>=4GiB/Efree>=10GiB afterimports/immediateprelaunch; freshsource/parent/process/licensegates afterwaits; allE:outputs/temp, D:Fluent; immutablebaseline/frozeninputs/failures; actualvisuals/honestlimitations/contributions.

Reviewfocus: copied0checkpointexactnesswithoutrepatch; phasealpha/pressure/velocitynativeordering anddepth; unchangedmethods/residualnorms everyactualstep includingCourantcounts; continuationactualglobalstep/counter/old-timefluxbinding; matchedphysicalstrip/probe and honestcrossinginterpolation/phaseinterpretation ratherthansilentphasealignment; sourcegatesafterwaits and no frozen-input mutation duringmesh/decomposition.

Task1 — reusable numerical analysis and source-bound inputs:
- [ ] Capture2567repohashes andallacceptedfine native/inputauditbindings; inspectinitialclocks/settings.
- [ ] RED: tests forcomplete-step/residual/Courantcoverage parser usingactualbaseline snippets and missing/duplicate-report rejection; exactinitialfield/settings comparison mustrejectalteredmethods.
- [ ] GREEN: minimal helpers; alltestsrun, retainfailures.
- [ ] CopyexactsourceFluentcase/data0/geometry andOFroot0/constant/system/nativepartition0meshes innewqualificationfolders; no decomposition/buildMesh ifexactpartitioninitialcopiesavailable. Changeonlyrootdt/end/write controls andnativecase dt; documentedheader-onlychangesifrequired.
Task2 — finite paired qualification:
- [ ] Freezeallinputs/driver/helpers/auditors; importedresource/process/licencegatesimmediateprelaunch.
- [ ] HiddenF80steps to.02, saveallstates/recoverybeforeclose; completeclosednativeaudit.
- [ ] OnlyafterFacceptance, O80steps withfreshsource/process/resourcegates; reconstructallsavedstates, auditnativephi/clocks/bothpartitions/residual/alpha/Courantcoverage.
- [ ] Inspectactualpaired.02fields/numericalplot, bindvisualreceipt. No fullcaseautomaticallyqueued.
Task3 — bounded fullinterval afterqualifications:
- [ ] Admitexact .02checkpointcontinuations, testedactualsteps81..4800/4720new, preserveexactsource andmerge80+4720.
- [ ] SequentialFthenO, savesat.01s physicaltimes,all4800stepsand121matchednativefields passoriginalgates. Finite1.2scap, no silentextension.
Task4 — comparison/media/localrecord/review:
- [ ] Testcrossinginterpolation/completecycle/countrejection; computeallmatchedstripdepth metrics/phaseevidencewithoutaccuracygate.
- [ ] Actualfixedscale4caseimages/depth/stepchecks, nativeframes-onlymovie/GIF/allframedecode; inspectrepresentatives.
- [ ] Existingstudy34time-step-v001 record, sourcecopies/baseline/exclusions/manifestchecks,ONEfreshtechnicalreview. OneCritical/Importantfixpass,Minorledger,finalchecks/completionproof. No uploadclaim.

Ruling: target fine-grid temporal response ratherthannewgeometry; previousfineCo<.05 is nottemporalsensitivityevidence. Qualificationusesexactacceptedtiltedinitialstateswithoutnewflatcontrol becausephysics/mesh/initialfieldsunchanged; cost explicit inSpec. Pendingevidencecanrejectthisroutewithoutweakenedlimits. Sourceappend-onlyledgerprogress.md.
