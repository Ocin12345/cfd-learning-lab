#!/usr/bin/env bash
set -eo pipefail
cd "$1"
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp
export XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
# MPI sockets/FIFOs require Linux filesystem semantics unavailable on this E: mount.
# Native outputs and ordinary temporary files remain on E:.
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-heated-jet-v002
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
trap 'rc=$?; printf "failed (%s)\n" "$rc" > run-status.txt; exit "$rc"' ERR
test "$WM_PROJECT_VERSION" = 14
if pgrep -x foamRun >/dev/null; then
    echo 'Another OpenFOAM run is active' >&2
    exit 3
fi
printf 'mesh\n' > run-status.txt
checkMesh -allTopology -allGeometry > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
printf 'decomposing\n' > run-status.txt
decomposePar -time '0,1250' > decompose.log 2>&1
for rank in 0 1; do
    test -s "processor$rank/0/T"
    test -s "processor$rank/1250/phi"
done
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
printf 'reconstructing\n' > run-status.txt
reconstructPar -time '0.1:3' > reconstruct.log 2>&1
test -s 2.995/T
for field in T U p phi; do test -s "3/$field"; done
printf 'complete\n' > run-status.txt

