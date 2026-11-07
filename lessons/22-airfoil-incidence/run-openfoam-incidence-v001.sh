#!/usr/bin/env bash
set -eo pipefail
run=/mnt/e/CFD-Learning-Lab/native/runs/2026-11-06-airfoil-incidence-paired-v001
export TMPDIR=/mnt/e/CFD-Learning-Lab/temp/airfoil-incidence-v001-openfoam
export TMP="$TMPDIR" TEMP="$TMPDIR" XDG_CACHE_HOME=/mnt/e/CFD-Learning-Lab/temp/openfoam-cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=/mnt/e/CFD-Learning-Lab/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
export FOAM_RUN="$WM_PROJECT_USER_DIR/run"
export OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
if pgrep -x 'foamRun|reconstructPar|simpleFoam|pimpleFoam' >/dev/null; then exit 3; fi
for alpha in 0 5; do
    cd "$run/openfoam-alpha$alpha"
    test ! -e solver.log
    test ! -d processor0
    test ! -d 0/polyMesh
    test ! -d constant/polyMesh
    stage=mesh_conversion
    status() { printf '{"stage":"%s","angle_degrees":%s,"ranks":2}\n' "$stage" "$alpha" > run-status.json; }
    trap 'rc=$?; printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s}\n" "$stage" "$rc" > run-status.json; exit "$rc"' ERR
    status
    foamVersion > installation-version.txt 2>&1
    fluentMeshToFoam -2D 1 input-airfoil.msh > mesh-conversion.log 2>&1
    stage=mesh_checks; status
    checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
    grep -q 'Mesh OK' mesh-check.log
    stage=decomposition; status
    decomposePar > decomposition.log 2>&1
    stage=solving; status
    mpirun -np 2 foamRun -parallel > solver.log 2>&1
    stage=checking_native_endpoint; status
    grep -Eq '^Time = 3000(s)?[[:space:]]*$' solver.log
    grep -q '^End' solver.log
    stage=reconstructing; status
    reconstructPar -latestTime > reconstruction.log 2>&1
    for f in U p phi; do test -s "3000/$f"; done
    stage=finished_execution; status
done
sync
