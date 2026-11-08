#!/usr/bin/env bash
set -eo pipefail
root=/mnt/e/CFD-Learning-Lab
source_dir="$root/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001"
case_dir="$root/native/runs/showcase-cylinder-convergence-v001/openfoam-developed-v001"
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-cylinder-developed-v001
test ! -e "$OMPI_MCA_orte_tmpdir_base"
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
cd "$case_dir"
test ! -e solver.log
test ! -e run-status.txt
if pgrep -x foamRun >/dev/null || pgrep -x mpirun >/dev/null; then exit 3; fi
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
foamVersion > installation-version.txt
python3 "$source_dir/verify_openfoam_developed_v001.py" > prerequisite-admission.log
if pgrep -x foamRun >/dev/null || pgrep -x mpirun >/dev/null; then exit 3; fi
trap 'rc=$?;printf "failed (%s)\n" "$rc" > run-status.txt;exit "$rc"' ERR
printf 'solving development interval\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
printf 'reconstructing development interval\n' > run-status.txt
reconstructPar -time '4.500001:8' -fields '(U p phi)' > reconstruct.log 2>&1
printf 'complete\n' > run-status.txt