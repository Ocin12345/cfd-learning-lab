# Reproducing model 34

This packet records the actual executed cases and their small inputs. Large native cases, installed dependencies, credentials and licence/connection files are excluded. It is not a complete portable native restart archive.

## Original environment and outputs

Project root: E:/CFD-Learning-Lab. Fluent: D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe. Python: environments/fluent/Scripts/python.exe. OpenFOAM: CFD-Foundation14 WSL distribution, /opt/openfoam14/etc/bashrc. Version details are in environment.json.

Accepted fine outputs are native/runs/sloshing-refinement-v001/production/fluent and production/openfoam. Forty-step controls are under qualification/{fluent-flat,openfoam-flat,fluent-tilted,openfoam-tilted}. The coarse references remain in native/runs/sloshing-production-v001/fluent and openfoam-v003, selected by local/sloshing-production-v001/accepted-pair-selection-v001.json. All originals are preserved.

## Execution order and code

The actual source directory is local/sloshing-refinement-v001. The archived code names retain their executed versions. Set PYTHONPATH to that directory plus repo/cfd-learning-lab/scripts, and source local/Use-Project.ps1 so project caches/temp remain on E:. Existing repository helpers analyse_nozzle_fluent_v005.py and wing_fluent_v001.py provide native cell mapping and lightweight memory admission; their installed dependencies are not bundled.

prepare_refined_openfoam_v001.py builds each mode from CFD_FINE_MODE=flat, tilted or production. sloshing_refinement_common_v001.py creates the 160x80 Fluent mesh and initialization. run_refinement_pipeline_v001.py ran flat Fluent, flat OpenFOAM, tilted Fluent, tilted OpenFOAM, full Fluent and full OpenFOAM sequentially. The Fluent driver also rechecks native prerequisite hashes immediately before launch, after resource waits and mesh preparation. The OpenFOAM wrapper rechecks them immediately before mpirun. Required admission is at least 4 GiB RAM and 10 GiB E: free; each solver uses two ranks. MPI IPC alone uses /dev/shm.

The copied prepared-inputs receipt and plan bind the original input/source versions. A rerun needs the separate local qualification and coarse native archives referenced by those hashes. Do not execute into accepted output folders or edit their frozen receipts. Prepare a fresh versioned root, update its task-specific paths and regenerate meshes, qualifications and source-hash receipts before production. The archived driver intentionally rejects reused outputs and changed prerequisites; it is not a one-command standalone replay from this Git packet.

## Audits and figures

The actual auditor is audit_refined_run_v001.py, called with solver, native case, fresh audit path and --mode flat, tilted or production. It checks actual native clocks, finite fields, initial fractions, all declared residual records, alpha bounds, conservation, Courant and mesh diagnostics. The public audits include the native source hashes. Do not relabel failures or clip alpha values to obtain a pass.

render_fine_pair_v001.py requires both passing production audits and source checks. compare_coarse_fine_v001.py computes matched-strip differences against the accepted coarse outputs. render_numerical_checks_v001.py exports the numerical PNG/PDF and CSVs. export_sloshing_video_v001.py and export_fine_sloshing_gif_v001.py encode/decode-check all 121 native-state PNGs. Their original local figure/video paths are retained. The public CSVs allow numerical results to be inspected without either solver; the native plotting scripts additionally require the separate native archives.

The fine animation's last-column centre is 199.375 mm. Coarse/fine measurement instead averages adjacent fine columns into 2.5 mm strips, including the rightmost 197.5–200 mm strip. Do not compare these two differently located traces under one label.

The code folder includes the failed initial mapping-check v001 and corrected v002 as history. The corrected check verifies coarse helper nodes, quads and centres against the analytic grid and hash-bound native mesh because the older coarse audit did not hash that helper file. Solver inputs and results were not modified by this correction.

## Interpretation

These checks establish the recorded numerical completion and two-grid differences at the stated conditions. They do not establish mesh/time independence, damping accuracy, a matched physical benchmark or aircraft validity. No personal reflection or independent owner work is inferred.

The archived PLAN files and prepared-input receipt are immutable preparation snapshots, so their initial unchecked tasks and launch flags remain historical. Actual completion is established by the full production audits and media records.

Courant-report coverage was separately verified for these completed logs:2,400 Fluent reports and2,401 OpenFOAM reports including initialization. The executed auditor requires reports but does not enforce their count; a future auditor should add that guard. See data/courant-report-coverage.json. The executed source is retained unchanged.
