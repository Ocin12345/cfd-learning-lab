#!/usr/bin/env bash
set -eo pipefail
case_dir=$(realpath "$1")
case "$case_dir" in
 /mnt/e/CFD-Learning-Lab/native/runs/dam-break-over-obstacle-v001/diagnostics/openfoam-pressure-cycle-baseline-v002|/mnt/e/CFD-Learning-Lab/native/runs/dam-break-over-obstacle-v001/diagnostics/openfoam-pressure-cycle-urelax07-v002) ;;
 *) exit 3;;
esac
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-dam-break-pressure-cycle-v002
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
cd "$case_dir"
test ! -e solver.log
test ! -e run-status.txt
if pgrep -x foamRun >/dev/null;then exit 3;fi
trap 'rc=$?;printf "failed (%s)\n" "$rc" > run-status.txt;exit "$rc"' ERR
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
python3 "$root/local/dam-break-over-obstacle-v001/verify_pressure_cycle_recipe_v002.py" > prerequisite-admission.log
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
printf 'reconstructing\n' > run-status.txt
reconstructPar -time '0.255001:0.26' -fields '(alpha.water U p_rgh)' > reconstruct.log 2>&1
printf 'complete\n' > run-status.txt