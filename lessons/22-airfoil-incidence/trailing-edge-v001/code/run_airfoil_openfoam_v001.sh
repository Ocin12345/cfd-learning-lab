#!/usr/bin/env bash
set -eo pipefail
root=/mnt/e/CFD-Learning-Lab
level="$1";mode="$2"
case "$level" in 1|2) ;; *) exit 2;; esac
cd "$root/native/runs/showcase-airfoil-trailing-edge-v001/openfoam-level${level}-v001"
export TMPDIR="$root/temp/airfoil-TE-openfoam-v001" TMP="$root/temp" TEMP="$root/temp" XDG_CACHE_HOME="$root/temp/openfoam-cache"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$TMPDIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR="$root/environments/openfoam-user"
test "$WM_PROJECT_VERSION" = 14
if pgrep -x 'foamRun|reconstructPar|simpleFoam|pimpleFoam|mpirun' >/dev/null; then exit 3; fi
stage=starting
status() { printf '{"stage":"%s","ranks":2}\n' "$stage" > run-status.json; }
trap 'rc=$?;printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s}\n" "$stage" "$rc" > run-status.json;exit "$rc"' ERR
if [ "$mode" = convert ];then
 test ! -e solver.log;test ! -d processor0;test ! -d 0/polyMesh;test ! -d constant/polyMesh
 stage=conversion;status
 foamVersion > installation-version.txt 2>&1
 fluentMeshToFoam -2D 1 input-airfoil.msh > mesh-conversion.log 2>&1
 stage=mesh_checks;status
 checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
 grep -q 'Mesh OK.' mesh-check.log
 decomposePar > decomposition.log 2>&1
 stage=converted_only;status
elif [ "$mode" = solve ];then
 test ! -e solver.log;test -d processor0;test -d processor1
 export OMPI_MCA_orte_tmpdir_base="/dev/shm/cfd-airfoil-TE-L${level}-v001"
 test ! -e "$OMPI_MCA_orte_tmpdir_base";mkdir -p "$OMPI_MCA_orte_tmpdir_base"
 stage=solving;status
 mpirun -np 2 foamRun -parallel > solver.log 2>&1
 grep -Eq '^Time = 3000(s)?[[:space:]]*$' solver.log;grep -q '^End' solver.log
 stage=reconstructing;status
 reconstructPar -latestTime > reconstruction.log 2>&1
 for f in U p phi;do test -s "3000/$f";done
 stage=finished_execution;status
else exit 2;fi
