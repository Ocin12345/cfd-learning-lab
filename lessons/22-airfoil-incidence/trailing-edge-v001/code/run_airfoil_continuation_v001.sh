#!/usr/bin/env bash
set -eo pipefail
root=/mnt/e/CFD-Learning-Lab
level="$1";endpoint="$2"
case "$level" in 1|2) ;; *) exit 2;; esac
[[ "$endpoint" =~ ^[0-9]+$ ]] || exit 2
cd "$root/native/runs/showcase-airfoil-trailing-edge-v001/openfoam-level${level}-v002"
export TMPDIR="$root/temp/airfoil-TE-continuation-v001" TMP="$root/temp" TEMP="$root/temp" XDG_CACHE_HOME="$root/temp/openfoam-cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR="$root/environments/openfoam-user"
export OMPI_MCA_orte_tmpdir_base="/dev/shm/cfd-airfoil-TE-L${level}-v002"
test ! -e "$OMPI_MCA_orte_tmpdir_base";mkdir -p "$OMPI_MCA_orte_tmpdir_base"
test ! -e solver-continuation.log
if pgrep -x 'foamRun|reconstructPar|simpleFoam|pimpleFoam|mpirun' >/dev/null; then exit 3; fi
stage=starting
status() { printf '{"stage":"%s","ranks":2}\n' "$stage" > run-status.json; }
trap 'rc=$?;printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s}\n" "$stage" "$rc" > run-status.json;exit "$rc"' ERR
stage=solving_saved_checkpoint;status
mpirun -np 2 foamRun -parallel > solver-continuation.log 2>&1
grep -Eq "^Time = ${endpoint}(s)?[[:space:]]*$" solver-continuation.log;grep -q '^End' solver-continuation.log
stage=reconstructing;status
reconstructPar -latestTime > reconstruction-continuation.log 2>&1
for f in U p phi;do test -s "$endpoint/$f";done
stage=finished_execution;status
