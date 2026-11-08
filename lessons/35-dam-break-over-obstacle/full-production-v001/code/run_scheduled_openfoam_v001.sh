#!/usr/bin/env bash
set -eo pipefail
root=/mnt/e/CFD-Learning-Lab
case_dir="$root/native/runs/dam-break-over-obstacle-v001/production/openfoam-scheduled-urelax07-v001"
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-dam-break-scheduled-v001
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
cd "$case_dir"
test ! -e solver.log
test ! -e run-status.txt
if pgrep -x foamRun >/dev/null;then exit 3;fi
trap 'rc=$?;printf "failed (%s)\n" "$rc" > run-status.txt;exit "$rc"' ERR
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
python3 "$root/local/dam-break-over-obstacle-v001/verify_scheduled_continuation_v001.py" > prerequisite-admission.log
printf 'solving continuation\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
printf 'reconstructing continuation\n' > run-status.txt
reconstructPar -time '0.250001:0.5' -fields '(alpha.water U p_rgh)' > reconstruct.log 2>&1
printf 'complete\n' > run-status.txt