#!/usr/bin/env bash
set -eo pipefail
case_dir=$(realpath "$1")
case "$case_dir" in
 /mnt/e/CFD-Learning-Lab/native/runs/dam-break-over-obstacle-v001/qualification/openfoam-flat-v002) stage=openfoam-flat ;;
 /mnt/e/CFD-Learning-Lab/native/runs/dam-break-over-obstacle-v001/qualification/openfoam-collapse-v002) stage=openfoam-collapse ;;
 /mnt/e/CFD-Learning-Lab/native/runs/dam-break-over-obstacle-v001/production/openfoam-v002) stage=openfoam-production ;;
 *) echo 'Unexpected case path' >&2;exit 3;;
esac
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-dam-break-v002
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
cd "$case_dir"
test ! -e solver.log
test ! -e run-status.txt
grep -q 'Mesh OK.' mesh-check.log
if pgrep -x foamRun >/dev/null;then echo 'Another OpenFOAM solver is active' >&2;exit 3;fi
trap 'rc=$?;printf "failed (%s)\n" "$rc" > run-status.txt;exit "$rc"' ERR
printf 'decomposing\n' > run-status.txt
decomposePar > decompose.log 2>&1
python3 /mnt/e/CFD-Learning-Lab/local/dam-break-over-obstacle-v001/dam_break_admission_v002.py "$stage" > prerequisite-admission.log
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
printf 'reconstructing\n' > run-status.txt
reconstructPar -noZero -fields '(alpha.water U p_rgh)' > reconstruct.log 2>&1
printf 'complete\n' > run-status.txt