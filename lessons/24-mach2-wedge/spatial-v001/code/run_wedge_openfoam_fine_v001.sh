#!/usr/bin/env bash
set -eo pipefail
cd /mnt/e/CFD-Learning-Lab/native/runs/showcase-mach2-wedge-spatial-v001/openfoam-fine-v001
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
test "$WM_PROJECT_VERSION" = 14
if pgrep -x 'foamRun|mpirun|simpleFoam|pimpleFoam' >/dev/null;then exit 3;fi
stage=starting
trap 'rc=$?;printf "failed %s (%s)\n" "$stage" "$rc" > run-status.txt;exit "$rc"' ERR
if [ "$1" = mesh ];then
 test ! -e solver.log;test ! -d constant/polyMesh;test ! -d processor0
 stage=mesh;printf 'meshing\n' > run-status.txt
 blockMesh > mesh.log 2>&1
 checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
 grep -q 'Mesh OK.' mesh-check.log
 foamPostProcess -func writeCellCentres -time 0 > preflow-centres.log 2>&1
 test -s 0/C
 decomposePar > decompose.log 2>&1
 printf 'mesh qualified, not solved\n' > run-status.txt
elif [ "$1" = solve ];then
 test ! -e solver.log;test -d processor0;test -d processor1
 export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-wedge-spatial-fine-v001
 test ! -e "$OMPI_MCA_orte_tmpdir_base";mkdir -p "$OMPI_MCA_orte_tmpdir_base"
 stage=solve;printf 'solving\n' > run-status.txt
 mpirun -np 2 foamRun -parallel > solver.log 2>&1
 grep -q '^End' solver.log
 stage=reconstruct;printf 'reconstructing\n' > run-status.txt
 reconstructPar -noZero > reconstruct.log 2>&1
 foamPostProcess -func writeCellCentres -time 0.00025 > cell-centres.log 2>&1
 test -s 0.00025/C;test -s 0.008/p;test -s 0.008/phi
 printf 'complete\n' > run-status.txt
else exit 2;fi
