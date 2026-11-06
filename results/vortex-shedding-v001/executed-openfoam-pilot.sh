#!/usr/bin/env bash
# Run --prepare-only first. The root agent admits resources before --run-pilot.
# This wrapper never starts Fluent and does not authorize a continuation itself.
set -eE -o pipefail
case_dir=/mnt/e/CFD-Learning-Lab/native/runs/2026-11-01-vortex-shedding-v001/openfoam
prep=/mnt/e/CFD-Learning-Lab/repo/cfd-learning-lab/scripts/prepare_openfoam_vortex.py
mode=${1:---prepare-only}
case "$mode" in
    --prepare-only|--run-pilot) test "$#" -le 1 ;;
    --continue-to-8s) test "$#" -eq 2 ;;
    *) printf 'Usage: %s [--prepare-only | --run-pilot | --continue-to-8s ROOT_AUTHORIZATION_JSON]\n' "$0" >&2; exit 2 ;;
esac
cd "$case_dir"
test -f input-vortex.msh
test -f provenance/prepared-input-sha256.json
stage=environment
status() {
    printf '{"stage":"%s","processes":2,"mode":"%s","initial_pilot_end_time_s":2,"delta_t_s":0.005}\n' "$stage" "$mode" > run-status.json
}
trap 'code=$?; printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s}\n" "$stage" "$code" > run-status.json; exit "$code"' ERR
export TMPDIR=/mnt/e/CFD-Learning-Lab/temp/vortex-shedding-v001
export TEMP="$TMPDIR" TMP="$TMPDIR"
export XDG_CACHE_HOME=/mnt/e/CFD-Learning-Lab/temp/openfoam-cache
export WM_PROJECT_USER_DIR=/mnt/e/CFD-Learning-Lab/environments/openfoam-user
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$WM_PROJECT_USER_DIR"
source /opt/openfoam14/etc/bashrc
# Foundation bashrc assigns its own HOME-based user directory; restore E: paths.
export WM_PROJECT_USER_DIR=/mnt/e/CFD-Learning-Lab/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
export FOAM_RUN="$WM_PROJECT_USER_DIR/run"
export OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
status

prepare_native_mesh() {
    test ! -e solver.log
    test ! -e solver-continuation-to-8s.log
    test ! -d processor0
    test ! -d processor1
    test ! -d constant/polyMesh
    test ! -d 0/polyMesh
    test ! -e mesh-conversion.log
    stage=converting; status
    foamVersion > installation-version.txt 2>&1
    fluentMeshToFoam -2D 1 input-vortex.msh > mesh-conversion.log 2>&1
    stage=checking_mesh; status
    checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
    grep -q 'Mesh OK' mesh-check.log
    stage=applying_native_inlet; status
    python3 -B "$prep" --case "$case_dir" --apply-inlet > inlet-profile-application.json
    test -s provenance/native-inlet-profile.json
    stage=ready_for_root_resource_admission; status
}

if test "$mode" = --prepare-only; then
    prepare_native_mesh
    exit 0
fi

# Linux-side guard complements the root's Windows Fluent/resource admission.
# The root must also ensure no concurrent Windows Fluent solver is active.
if pgrep -x foamRun >/dev/null || pgrep -x simpleFoam >/dev/null || pgrep -x pimpleFoam >/dev/null; then
    printf 'An OpenFOAM solver is already running. Refusing concurrent solve.\n' >&2
    exit 3
fi

if test "$mode" = --run-pilot; then
    test ! -e solver.log
    test ! -e numerical-summary-pilot.json
    if test ! -d constant/polyMesh
    test ! -d 0/polyMesh; then
        prepare_native_mesh
    fi
    test -s provenance/native-inlet-profile.json
    grep -q 'Mesh OK' mesh-check.log
    cmp system/controlDict provenance/initial-system/controlDict
    test ! -e provenance/run-openfoam-pilot.sh
    cp "$0" provenance/run-openfoam-pilot.sh
    stage=decomposing; status
    decomposePar > decomposition.log 2>&1
    stage=solving_pilot; status
    mpirun -np 2 foamRun -parallel > solver.log 2>&1
    log=solver.log
    end=2
    reconstruction_log=reconstruction.log
    summary_log=numerical-summary-pilot.stdout.json
    centres_log=cell-centres.log
else
    test -e numerical-summary-pilot.json
    test -d processor0/2
    test -d processor1/2
    test ! -e solver-continuation-to-8s.log
    test ! -e provenance/run-openfoam-continuation-to-8s.sh
    cp "$0" provenance/run-openfoam-continuation-to-8s.sh
    stage=checking_root_continuation_authorization; status
    python3 -B "$prep" --case "$case_dir" --prepare-continuation --authorization "$2"
    stage=solving_authorized_continuation; status
    mpirun -np 2 foamRun -parallel > solver-continuation-to-8s.log 2>&1
    log=solver-continuation-to-8s.log
    end=8
    reconstruction_log=reconstruction-continuation-to-8s.log
    summary_log=numerical-summary-8s.stdout.json
    centres_log=cell-centres-8s.log
fi
awk '/^Time =|Solving for|continuity errors|Courant Number|ExecutionTime/' "$log" > "${log%.log}-convergence.log"
grep -q '^End' "$log"
stage=reconstructing_saved_snapshots; status
if test "$end" = 2; then
    reconstructPar -time '0.025:2' > "$reconstruction_log" 2>&1
else
    reconstructPar -time '2.025:8' > "$reconstruction_log" 2>&1
fi
for field in U p phi; do
    test -s "$end/$field"
done
stage=writing_cell_centres; status
foamPostProcess -func writeCellCentres -time "$end" > "$centres_log" 2>&1
test -s "$end/C"
stage=checking_final_mass_flux; status
python3 -B "$prep" --case "$case_dir" --summarize "$end" > "$summary_log"
printf '{"stage":"finished_execution","processes":2,"final_time_s":%s,"delta_t_s":0.005,"saved_field_interval_s":0.025,"force_sample_interval_s":0.005,"paired_analysis_pending":true,"physical_validation":false}\n' "$end" > run-status.json
