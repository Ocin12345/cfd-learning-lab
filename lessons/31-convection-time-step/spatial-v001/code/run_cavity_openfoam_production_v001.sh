#!/usr/bin/env bash
set -eo pipefail
cd /mnt/e/CFD-Learning-Lab/native/runs/showcase-convection-spatial-v001/openfoam-production-v001
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-convection-spatial-production-v001
test "$WM_PROJECT_VERSION" = 14
test ! -e solver.log;test ! -e 20;test -s 1/p;test -s processor0/1/p;test -s processor1/1/p
test ! -e "$OMPI_MCA_orte_tmpdir_base"
if pgrep -x 'foamRun|mpirun|simpleFoam|pimpleFoam' >/dev/null;then exit 3;fi
trap 'rc=$?;printf "failed (%s)\n" "$rc" > run-status.txt;exit "$rc"' ERR
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
printf 'continuing1to20s\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
reconstructPar -time '1.5:20' > reconstruct.log 2>&1
foamPostProcess -solver fluid -func wallHeatFlux -latestTime > wall-heat-postprocess.log 2>&1
test -s 20/p;test -s 20/phi;test -s 20/wallHeatFlux
printf 'complete\n' > run-status.txt