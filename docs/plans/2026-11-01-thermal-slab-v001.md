# Paired transient solid conduction implementation plan

Spec E:/CFD-Learning-Lab/local/thermal-slab-v001/DESIGN.md. Constant-property heat model, true200s transient and actual native temperature animation.

1. Shared2000quad all-wall geometry with namedslab cell zone; root independently checks rectangle,areas,topology/hashes. Separate canonical analytic cell-average Fourier solution with independent series/mean checks.
2. Freeze guarded2processFluent solid-energy driver/analyser and OpenFOAM14solid/Euler driver/analyser. Source review plus separate resource/nativequality/startup admission. Save0+100states, native wallQ/time/dt/material proofs. No concurrent solvers.
3. Validate actual native temperature/time/solidproperties, residuals at all solvedsteps, bounds/1Dyuniformity, native global Euler energy balance with declared Qref and final analytic error targets. Preserve every miss and early-time error history; no target relaxation.
4. Actual51state0:4:200GIF50ms plus clear temperature profiles, energy/analytic-error figures. Main fields always solver output; reference clearly labelled. Root inspect actualPNG/decodedGIF and verify hashes.
5. Back up876baseline/globaldocs, record small reproducible methods/code/inputs/metrics/visuals with limitations/contributions, independently verify copies/links/AST/exclusions/rootQA and regenerate manifest LAST. Large native snapshots remain separate. No Git or publication.

Euler chosen to match unambiguous startup across codes. Boundary-step singularity and coarse early timestep are explicit limits; final analytical verification is distinct from physical validation and mesh/time independence.
