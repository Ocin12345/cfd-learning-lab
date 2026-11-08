# Methods: local trailing-edge sensitivity at 5° and Re = 100

This extends existing studies 03 and 22. It is a local spatial sensitivity experiment on the accepted 5° mesh, not a new model, a global mesh-independence study, or an aircraft validation.

## Geometry and mesh controls

The accepted baseline has 109,600 triangular cells. The airfoil polygon, unit chord, 5° incidence and outer boundary at a 20-chord radius are reused. All original node coordinates remain bit-identical. New body-boundary nodes are collinear midpoints on existing straight segments; the trailing edge is not rounded or reshaped. Far-boundary coordinates, edges and patch tags are unchanged. Unaffected baseline triangles are preserved exactly.

Triangles near a radius of 0.08 chord from the actual trailing-edge vertex are marked for local red refinement. Conforming green splits close the transition. Quality-aware closure expands red parent marking where a transition child is too poor; it does not move existing nodes. The accepted affected extent is 0.0891698453 chord, inside the predeclared 0.15-chord limit. Level 1 affects 3,378 baseline cells; level 2 affects 3,391. The resulting meshes contain 119,580 and 159,300 cells.

Native Fluent minimum orthogonal quality is 0.325002 for level 1 and 0.225313 for level 2, both above the unchanged 0.20 floor. These are native quality results, not the helper's triangle shape metric. Finer resolution does not imply higher minimum cell quality. A previous transition mesh failed at 0.115009 and was rejected before any flow iteration; its native files and all failed source versions remain preserved.

The emitted setup files inherit some historical fields from the original incidence morph. In particular, `connectivity_unchanged`, `min_area_ratio`, `max_area_ratio` and morph radii describe that earlier morph, not the subsequent local split. Current mesh hashes, native counts and `geometry-proof.json` govern this study. Local connectivity changes; outside-patch triangles and the far boundary do not.

## Physical model and numerical methods

Both codes solve steady planar incompressible laminar flow with density 1 kg/m³, dynamic viscosity 0.01 Pa s, inlet velocity (1, 0) m/s, a stationary no-slip body and zero gauge pressure at the outlet. Chord, span and force-reference area are each unity, giving Re = 100 and reference dynamic pressure 0.5 Pa. Lift is the global positive y force; drag is the global positive x force.

Fluent uses the accepted coupled, second-order numerical template and its own hybrid initialization on each fresh mesh. OpenFOAM Foundation 14 uses a one-layer unit-depth extrusion with empty planes, steady consistent SIMPLE, least-squares gradients, bounded linear-upwind momentum, corrected diffusion, two nonorthogonal pressure correctors and U equation relaxation 0.7. Its initial state is uniform U = (1, 0, 0), kinematic pressure zero. These initializations and residual norms differ between codes.

Fluent runs have a finite 600-iteration initial cap followed by 100 extra observation iterations; normal convergence may stop the initial block earlier. Native force reports are collected at all 101 observation states. OpenFOAM begins with a finite 3,000-iteration cap. Underconverged endpoints are retained as misses. Any continuation is an exact byte-verified copy of the native root and both partition checkpoints; no reinitialization occurs. Only root restart, end and write controls may change, with force-function and physical/numerical dictionaries preserved. Actual selected endpoint and provenance are recorded separately after acceptance.

## Numerical acceptance and native-field checks

The predeclared final residual target is 1e-7, relative signed boundary mass imbalance at most 1e-6, drag relative range at most 1e-4 and lift absolute range at most 1e-5. Fluent uses its last 20 native force samples; OpenFOAM retains its original stricter 100-sample window. The 5° case has no zero-incidence lift-symmetry gate.

Fluent audits inspect saved HDF case/data ranges, coordinate bijection and face connectivity against the emitted mesh, stationary no-slip wall flags, physical settings, finite velocity/pressure/material fields, every native scalar force report, force-reference normalization, terminal residuals and signed mass report. OpenFOAM audits inspect the imported mesh hash, every projected prism and face-owner/neighbour/patch identity, reconstructed native cell ordering, both partition endpoints, complete residual/force iteration history, physical/numerical dictionaries and native endpoint phi boundary closure. The OpenFOAM residual statistic is the maximum printed initial Ux, Uy and pressure residual within the final iteration, including every pressure-correction solve.

A printed last-window force range of zero means stability at the native force report's printed precision. It is not zero discretization error. A completed solver process is not accepted until its independent audit passes.

## Images and pressure samples

Comparison maps read actual native cell fields with no temporal interpolation. All six cases use identical spatial views and shared speed/pressure colour ranges that include every finite saved value. Source hashes and field extrema are retained. Steady maps do not describe transient fluid motion.

Pressure coefficient is Cp = (p - p_ref)/(0.5 rho U²). Fluent stores gauge pressure in Pa. OpenFOAM stores kinematic pressure; rho = 1 makes its numerical Cp conversion the same. Surface plots use Fluent's saved face SV_P and the evaluated OpenFOAM zeroGradient wall pressure, which equals the owner-cell value. The adjacent-cell samples are retained separately. In these Fluent cases the saved wall-face pressure also equals the adjacent-cell pressure exactly. Neither provides an independent extrapolated continuum surface pressure.

## Interpretation and provenance

This controls local resolution while holding the polygon and far domain fixed. It does not establish global mesh or domain independence, a converged sharp-tip point pressure, a unique cause for cross-code differences, experimental validation or aircraft performance. No Richardson/GCI estimate is inferred from these local, nonuniform mesh changes.

The user selected and directed the comparison. Codex prepared and executed the scripts, checks, figures and record. ANSYS Fluent and OpenFOAM computed their respective fields. No owner execution, independent work, research novelty or personal reflection is invented. User reflection remains unprovided. Native solver files, installed environments, licence information, credentials and private application material stay outside the GitHub-ready content. No new upload is claimed.

Sampled wall-pressure minima occur at mesh-dependent face locations. Refinement moves samples closer to the sharp trailing-edge vertex. A changing sampled minimum does not compare pressure at one fixed physical point, prove a mathematical singularity, or establish pointwise nonconvergence. The full traces and actual sample coordinates must be retained alongside extrema.

The geometry-proof flag native_orthogonal_quality_not_yet_measured is a preserved snapshot from mesh creation. Later source-bound native quality qualifications are authoritative; historical proof files are not rewritten. Project temporary files use E: where supported; OpenMPI FIFO metadata uses an isolated project path in /dev/shm because that mechanism is unsupported on E:.

The inherited 0.5% cross-code drag learning goal is a separate diagnostic from each solver's numerical endpoint acceptance. It remains a miss wherever the recorded drag gap exceeds that goal. Passing each solver's residual/mass/stability limits does not establish cross-code agreement or physical accuracy. Drag gaps use OpenFOAM as the denominator, consistently with the parent study.
