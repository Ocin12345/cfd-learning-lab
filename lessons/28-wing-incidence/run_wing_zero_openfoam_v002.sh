#!/usr/bin/env bash
set -eo pipefail
root=/mnt/e/CFD-Learning-Lab
cd "$root/native/runs/finite-wing-zero-openfoam-v002"
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp
export XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
# Linux IPC primitives are not supported by the E: filesystem; only IPC uses RAM.
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-wing-zero-v002
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
trap 'rc=$?; printf "failed (%s)\n" "$rc" > run-status.txt; exit "$rc"' ERR
test "$WM_PROJECT_VERSION" = 14
if pgrep -x foamRun >/dev/null; then echo 'Another OpenFOAM solver is active' >&2; exit 3; fi
test ! -e solver.log
printf 'checking_mesh\n' > run-status.txt
checkMesh -allTopology -allGeometry > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
printf 'decomposing\n' > run-status.txt
decomposePar > decomposition.log 2>&1
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
grep -Eq '^Time = 1200(s)?[[:space:]]*$' solver.log
printf 'reconstructing\n' > run-status.txt
reconstructPar -time '300,600,900,1200' > reconstruction.log 2>&1
for t in 300 600 900 1200; do for f in U p phi; do test -s "$t/$f"; done; done
foamPostProcess -func writeCellCentres -latestTime > cell-centres.log 2>&1
test -s 1200/C
printf 'complete\n' > run-status.txt
