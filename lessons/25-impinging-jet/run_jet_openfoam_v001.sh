#!/usr/bin/env bash
set -eo pipefail
cd "$1"
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp
export XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=$root/temp
trap 'rc=$?; printf "failed (%s)\n" "$rc" > run-status.txt; exit "$rc"' ERR
test "$WM_PROJECT_VERSION" = 14
if pgrep -x foamRun >/dev/null; then
    echo 'Another OpenFOAM run is active' >&2
    exit 3
fi
printf 'mesh\n' > run-status.txt
blockMesh > mesh.log 2>&1
checkMesh -allTopology -allGeometry > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
printf 'decomposing\n' > run-status.txt
decomposePar > decompose.log 2>&1
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
printf 'reconstructing\n' > run-status.txt
reconstructPar -time 1000:1250 > reconstruct.log 2>&1
foamPostProcess -func writeCellCentres -time 1250 > cell-centres.log 2>&1
test -s 1250/C
printf 'complete\n' > run-status.txt

