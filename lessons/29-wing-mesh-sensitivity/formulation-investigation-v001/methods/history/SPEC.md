# Existing refined-wing iterative convergence repair, v001

## Question and bounded choice

Can a tighter inner pressure-coupled AMG termination criterion remove the refined Fluent residual floor, compared with the same number of additional iterations using its original settings? This improves existing study 29, with accepted baseline study 27 and matched refined OpenFOAM reused. The separate zero-incidence miss in study 28 remains open.

Fresh native inventory confirms 224,972 common tetrahedral cells and original iteration 1500 continuity 1.5753945420957494e-5. Its last100 range is about1.44e-5…2.89e-5; force coefficients are already stable. Summed absolute cell mass imbalance is6.536745231970203e-7kg/s. Thus simply requesting more original iterations is not assumed to fix convergence. The native Coupled F-cycle uses pressure-coupled termination0.1. Hypothesis: tightening that INNER linear-block control to0.001 may reduce the outer residual floor. It changes no outer acceptance target and success is not predicted.

## Fixed physical and spatial model

Keep the exact saved 1500-iteration refined case/data and native triangulation. No remeshing, repatching, reinitialization, residual renormalization or counter reset. Wing chord1m/span4m/area4m2/AR4,8degrees clockwise pitch, rounded-cap NACA0012-like section (not exact NACA0012). Boxx=-3..9,y=-4..4,z=-6..6m. U=(1,0,0)m/s,rho1kg/m3,mu.0025Pa s,Re400. Steady pressure-based incompressible laminar flow; stationary no-slipwing,inlet/farfield prescribedU,outlet0gaugepressure. Reference density1,speed1,length1,area4,dragx/lifty gives coefficient denominator2N. All native settings/materials/references/BCs and geometrical mappings must match source.

Coupled algorithm,pseudo-timeoff,Courant20,explicit pressure/momentum relaxation0.5,least-squares cell gradient,second-orderpressure/second-orderupwindmomentum remain fixed. Original residual normalization/scaling flags and native scaling history retained. Coupled maxcycle10 and other AMG controls unchanged. Existing OpenFOAM steady consistentSIMPLE/ownnumerics/force definitions remain unchanged and its accepted native1500endpoint is rehashed/reused, not rerun.

## Finite sequential qualification

From two fresh byte-verified copies of EXACT same1500checkpoint run100actualiterations each to1600, one session atatime/two ranks. Baseline retains pressure-coupled AMG termination0.1; candidate changes ONLY that criterion to0.001. Verbosity1 in both is an output-only diagnostic change. Save iteration1500initialfield and checkpoints1520,1540,1560,1580,1600,per-iteration native force reports, native mass/drag/lift reports and loaded/final settings. Redirect ALL monitor output paths into the new case before any iteration. Preserve all original files and refuse output reuse. Enable/restore residual check flags only to verify the unchanged stored1e-5 thresholds; fixed-run stopping remains off, as in source. No source history rescaling.

Qualification PASS per case requires endpoint all4scaledresiduals<=1e-5 AND all final20new iterations<=1e-5, relative global boundary mass imbalance<=1e-4, final100NEW force samples consecutive1501..1600 with range/max(abs(mean),.001)<=.001 forCd/Cl, finite U/p/rho/mu, rho1/mu.0025 to1e-12, peak speed<10m/s, native mesh quality>=.10 with positive volumes, exact sourcegeometry and physical/spatial settings, initialfieldexactness, inherited residual-prefix preservation, correct actualglobaliteration counts, and nativeforce coefficient/report consistency (rtol1e-6/atol1e-10). Additionally terminal summedabsolute cellmassimbalance must not exceed source6.536745231970203e-7kg/s. This strengthens evidence; it never relaxes original targets. Save every MISS; baseline expected MISS does not block the candidate if its structure/physics/finite-field/resource gates pass. Any structural/physics failure blocks the pipeline. Runtime divergenceguard at savedstates:maxscaledres<0.1; this is a stop guard,not acceptance.

Select baseline if it passes; otherwise select tight if it passes; otherwise neither is accepted and no extension follows. Explain the evidence, not an assumed unique cause. Actual accepted endpoint native fields/residual/force plots must be generated, source-bound and inspected before full confirmation admission.

## Finite confirmation, not yet launched

Only selected accepted1600case/data may continue in a NEW directory for100further iterations to1700. No reset/init. Require EVERYnew confirmation iteration1601..1700 all4residuals<=1e-5 plus unchanged mass/force/mesh/settings/field/restart/localclosure checks. This is a finite iterative stability window, not physical-time development or validation. Never silently extend beyond1700. The existing1500iteration MISS stays preserved and explicitly inherited.

## Final evidence and operation

Reuse accepted OpenFOAM and baselinewing fields with source hashes, matched force/reference/BC audit, real shared-scale native pressure/velocity views and numerical histories. Steady-field tracers/camera animation are not simulated transients. No induced-drag decomposition, exactNACA0012, independent full rerun, mesh/time/domain independence, aircraft validation or forced code agreement. Owner execution/reflection not invented. Complete a versioned study29record,2798repo baseline preservation/manifest/privateexclusions and ONEfreshwhole-improvementreview only at the end.

All project code/env/native/data/temp/exports E:, existing D:Fluent261,hidden/background,two ranks,ONEsolver session atatime. Afterimports and immediateprelaunch/afterwaits source/parent/executable/process/RAM>=4GiB/Efree>=10GiB gates. Verify actual licence/version connection before advancing. Save beforeclosing. No simultaneous Fluent/OpenFOAM, desktop takeover, original deletion/overwrites, native/environment/licence/credentials/private app material inGit or unverified upload.
