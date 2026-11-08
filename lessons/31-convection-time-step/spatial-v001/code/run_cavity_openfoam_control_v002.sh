#!/usr/bin/env bash
set -eo pipefail
cd /mnt/e/CFD-Learning-Lab/native/runs/showcase-convection-spatial-v001/openfoam-control-v002
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
test "$WM_PROJECT_VERSION" = 14
if pgrep -x 'foamRun|mpirun|simpleFoam|pimpleFoam' >/dev/null;then exit 3;fi
trap 'rc=$?;printf "failed (%s)\n" "$rc" > run-status.txt;exit "$rc"' ERR
if [ "$1" = mesh ];then
 test ! -e solver.log;test ! -d processor0;test ! -d 0/polyMesh
 fluentMeshToFoam -2D 0.005 input-cavity.msh > mesh-conversion.log 2>&1
 checkMesh -allTopology -allGeometry > mesh-check.log 2>&1
 grep -q 'Mesh OK.' mesh-check.log
 foamPostProcess -func writeCellCentres -time 0 > preflow-centres.log 2>&1
 test -s 0/C
 decomposePar > decompose.log 2>&1
 printf 'mesh qualified\n' > run-status.txt
elif [ "$1" = solve ];then
 test ! -e solver.log;test -d processor0;test -d processor1
 export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-convection-spatial-control-v002
 test ! -e "$OMPI_MCA_orte_tmpdir_base";mkdir -p "$OMPI_MCA_orte_tmpdir_base"
 printf 'steady isothermal control running\n' > run-status.txt
 mpirun -np 2 foamRun -parallel > solver.log 2>&1
 grep -q '^End' solver.log
 reconstructPar -latestTime > reconstruct.log 2>&1
 foamPostProcess -func writeCellCentres -latestTime > centres.log 2>&1
 test -s 600/C;test -s 600/U
 printf 'complete\n' > run-status.txt
else exit 2;fi