#!/usr/bin/env bash
set -eo pipefail
cd /mnt/e/CFD-Learning-Lab/native/runs/hot-patch-plume-qualification-v001/openfoam
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-plume-qualification-v001
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
trap 'rc=$?; printf "failed (%s)\n" "$rc" > run-status.txt; exit "$rc"' ERR
test ! -e solver.log
if pgrep -x foamRun >/dev/null; then exit 3; fi
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
decomposePar > decompose.log 2>&1
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
reconstructPar > reconstruct.log 2>&1
foamPostProcess -func writeCellCentres -latestTime > centres.log 2>&1
printf 'complete\n' > run-status.txt