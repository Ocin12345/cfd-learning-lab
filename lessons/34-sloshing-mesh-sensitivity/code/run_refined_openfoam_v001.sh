#!/usr/bin/env bash
set -eo pipefail
case_dir=$(realpath "$1")
case "$case_dir" in
 /mnt/e/CFD-Learning-Lab/native/runs/sloshing-refinement-v001/qualification/openfoam-flat) stage=openfoam-flat ;;
 /mnt/e/CFD-Learning-Lab/native/runs/sloshing-refinement-v001/qualification/openfoam-tilted) stage=openfoam-tilted ;;
 /mnt/e/CFD-Learning-Lab/native/runs/sloshing-refinement-v001/production/openfoam) stage=openfoam-production ;;
 *) echo 'Unexpected refined-case directory' >&2;exit 3 ;;
esac
cd "$case_dir"
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-sloshing-refinement-v001
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
test ! -e solver.log
if pgrep -x foamRun >/dev/null;then echo 'Another OpenFOAM solver is active' >&2;exit 3;fi
trap 'rc=$?;printf "failed (%s)\n" "$rc" > run-status.txt;exit "$rc"' ERR
printf 'meshing\n' > run-status.txt
blockMesh > block-mesh.log 2>&1
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
foamPostProcess -func writeCellCentres -time 0 > centres.log 2>&1
decomposePar > decompose.log 2>&1
python3 /mnt/e/CFD-Learning-Lab/local/sloshing-refinement-v001/sloshing_refinement_admission_v001.py "$stage" > prerequisite-admission.log
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
printf 'reconstructing\n' > run-status.txt
reconstructPar -noZero -fields '(alpha.water U p_rgh)' > reconstruct.log 2>&1
printf 'complete\n' > run-status.txt
