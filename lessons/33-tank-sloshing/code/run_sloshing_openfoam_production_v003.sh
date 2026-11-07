#!/usr/bin/env bash
set -eo pipefail
case_dir=$(realpath "$1")
test "$case_dir" = /mnt/e/CFD-Learning-Lab/native/runs/sloshing-production-v001/openfoam-v003 || { echo 'Unexpected production directory' >&2; exit 3; }
cd "$case_dir"
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-sloshing-production-v003
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
test ! -e solver.log
if pgrep -x foamRun >/dev/null; then echo 'Another OpenFOAM solver is active' >&2; exit 3; fi
trap 'rc=$?; printf "failed (%s)\n" "$rc" > run-status.txt; exit "$rc"' ERR
python3 /mnt/e/CFD-Learning-Lab/local/sloshing-production-v001/qualification_admission_v001.py > qualification-admission.log
printf 'meshing\n' > run-status.txt
blockMesh > block-mesh.log 2>&1
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
decomposePar > decompose.log 2>&1
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
printf 'reconstructing\n' > run-status.txt
reconstructPar -noZero > reconstruct.log 2>&1
printf 'complete\n' > run-status.txt
