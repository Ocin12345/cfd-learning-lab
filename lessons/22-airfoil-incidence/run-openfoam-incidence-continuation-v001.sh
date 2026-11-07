#!/usr/bin/env bash
set -eo pipefail
run=/mnt/e/CFD-Learning-Lab/native/runs/2026-11-06-airfoil-incidence-paired-v001
export TMPDIR=/mnt/e/CFD-Learning-Lab/temp/airfoil-incidence-v001-openfoam
export TEMP="$TMPDIR" TMP="$TMPDIR" XDG_CACHE_HOME=/mnt/e/CFD-Learning-Lab/temp/openfoam-cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
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
    cd "$run/openfoam-alpha$alpha-v002"
    stage=continuation_preflight
    status() { printf '{"stage":"%s","angle_degrees":%s,"ranks":2,"from":3000,"end":4800}\n' "$stage" "$alpha" > continuation-status.json; }
    trap 'rc=$?; printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s}\n" "$stage" "$rc" > continuation-status.json; exit "$rc"' ERR
    status
    test ! -e solver-continuation.log
    test -d processor0/3000
    test -d processor1/3000
    test ! -d processor0/4800
    test ! -d processor1/4800
    stage=solving; status
    mpirun -np 2 foamRun -parallel > solver-continuation.log 2>&1
    stage=checking_endpoint; status
    grep -Eq '^Time = 3001(s)?[[:space:]]*$' solver-continuation.log
    grep -Eq '^Time = 4800(s)?[[:space:]]*$' solver-continuation.log
    grep -q '^End' solver-continuation.log
    stage=reconstructing; status
    reconstructPar -time 4800 > reconstruction-4800.log 2>&1
    for f in U p phi; do test -s "4800/$f"; done
    stage=finished_execution; status
done
sync
