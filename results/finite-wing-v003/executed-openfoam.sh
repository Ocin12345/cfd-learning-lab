#!/usr/bin/env bash
# Fresh sequential 3D finite-wing execution only; root authorization is required.
set -e -o pipefail
case_dir=/mnt/e/CFD-Learning-Lab/native/runs/2026-11-01-finite-wing-v003/openfoam
prep=/mnt/e/CFD-Learning-Lab/repo/cfd-learning-lab/scripts/prepare_openfoam_finite_wing_v003.py
mode=${1:---prepare-only}
case "$mode" in --prepare-only|--run) test "$#" -eq 2 ;; *) exit 2 ;; esac
authorization=$2
test -f "$authorization"
python3 -B - "$authorization" "$mode" <<'PY'
import json,sys
a=json.load(open(sys.argv[1]))
assert a.get('approved') is True and a.get('approved_by')=='root'
assert a.get('action')=='run_finite_wing_openfoam_v003'
assert a.get('allowed_mode')=={'--prepare-only':'PrepareOnly','--run':'Run'}[sys.argv[2]]
PY
cd "$case_dir"
export TMPDIR=/mnt/e/CFD-Learning-Lab/temp/finite-wing-v003
export TEMP="$TMPDIR" TMP="$TMPDIR"
export XDG_CACHE_HOME=/mnt/e/CFD-Learning-Lab/temp/openfoam-cache
export WM_PROJECT_USER_DIR=/mnt/e/CFD-Learning-Lab/environments/openfoam-user
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$WM_PROJECT_USER_DIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=/mnt/e/CFD-Learning-Lab/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
export FOAM_RUN="$WM_PROJECT_USER_DIR/run"
export OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
# Optional environment probes can return nonzero inside substitutions. The trap
# begins only after bashrc, and is not inherited into command substitutions.
stage=environment_ready
status() { printf '{"stage":"%s","processes":2,"mode":"%s","solver_started":false}\n' "$stage" "$mode" > run-status.json; }
trap 'code=$?; printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s}\n" "$stage" "$code" > run-status.json; exit "$code"' ERR
status
if pgrep -x foamRun >/dev/null || pgrep -x simpleFoam >/dev/null || pgrep -x pimpleFoam >/dev/null || pgrep -x fluent >/dev/null; then
    printf 'Another solver is running; refuse concurrent solve.\n' >&2
    exit 3
fi
if test "$mode" = --prepare-only; then
    test ! -e solver.log
    test ! -e mesh-conversion.log
    test ! -d constant/polyMesh
    test ! -d 0/polyMesh
    test ! -d processor0
    test ! -d processor1
    stage=converting_3d; status
    foamVersion > installation-version.txt 2>&1
    fluentMeshToFoam input-wing.msh > mesh-conversion.log 2>&1
    stage=checking_mesh; status
    checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
    stage=validating_native_3d; status
    python3 -B "$prep" --case "$case_dir" --validate-native > native-preflight.stdout.json
    test -s provenance/preflight-validation.json
    stage=ready_for_separate_windows_resource_admission; status
    exit 0
fi
test ! -e solver.log
test ! -d processor0
test ! -d processor1
test -s provenance/preflight-validation.json
grep -q 'Mesh OK.' mesh-check.log
test ! -e provenance/run-openfoam-executed.sh
cp "$0" provenance/run-openfoam-executed.sh
cp "$authorization" provenance/root-authorization.json
stage=decomposing; status
decomposePar > decomposition.log 2>&1
stage=solving_fixed_1200_iterations
printf '{"stage":"%s","processes":2,"solver_started":true,"iteration_cap":1200}\n' "$stage" > run-status.json
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
grep -Eq '^Time = 1200(s)?[[:space:]]*$' solver.log
stage=reconstructing_saved_iterations
printf '{"stage":"%s","processes":2,"solver_started":true}\n' "$stage" > run-status.json
reconstructPar -time '300,600,900,1200' > reconstruction.log 2>&1
for iteration in 300 600 900 1200; do
    for field in U p phi; do test -s "$iteration/$field"; done
done
stage=writing_actual_cell_centres
printf '{"stage":"%s","processes":2,"solver_started":true}\n' "$stage" > run-status.json
foamPostProcess -func writeCellCentres -latestTime > cell-centres.log 2>&1
test -s 1200/C
awk '/^Time =|Solving for|continuity errors|ExecutionTime/' solver.log > solver-convergence.log
printf '{"stage":"finished_execution","processes":2,"final_iteration":1200,"saved_iterations":[300,600,900,1200],"solver_started":true,"numerical_analysis_pending":true,"physical_validation":false}\n' > run-status.json
