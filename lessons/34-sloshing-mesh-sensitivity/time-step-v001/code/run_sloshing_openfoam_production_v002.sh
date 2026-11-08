#!/usr/bin/env bash
set -eo pipefail
cd /mnt/e/CFD-Learning-Lab/native/runs/showcase-sloshing-time-step-v001/openfoam-production-v001
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-sloshing-dt-production-v001
test "$WM_PROJECT_VERSION" = 14
test ! -e solver.log;test -s processor0/0.02/alpha.water;test -s processor1/0.02/alpha.water
test ! -e "$OMPI_MCA_orte_tmpdir_base"
if pgrep -x 'foamRun|mpirun|simpleFoam|pimpleFoam' >/dev/null;then exit 3;fi
trap 'rc=$?;printf "failed (%s)\n" "$rc" > run-status.txt;exit "$rc"' ERR
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
/mnt/e/CFD-Learning-Lab/environments/fluent/Scripts/python.exe -B E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/sloshing-time-step-v001/gate_sloshing_full_mpi_v002.py > immediate-pre-mpi-admission.json
printf 'continuing.02to1.2s\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
reconstructPar -time '0.03:1.2' -fields '(alpha.water U p_rgh phi)' > reconstruct.log 2>&1
printf 'complete\n' > run-status.txt
