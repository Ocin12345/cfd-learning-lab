#!/bin/bash
set -eo pipefail
export TMPDIR=/mnt/e/CFD-Learning-Lab/temp TMP=/mnt/e/CFD-Learning-Lab/temp TEMP=/mnt/e/CFD-Learning-Lab/temp
export XDG_CACHE_HOME=/mnt/e/CFD-Learning-Lab/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=/mnt/e/CFD-Learning-Lab/environments/openfoam-user
test "$WM_PROJECT_VERSION" = 14
if [ "${1:-}" = "--probe" ]; then printf 'Foundation version=%s\n' "$WM_PROJECT_VERSION"; command -v foamRun; command -v reconstructPar; exit 0; fi
export TMPDIR=/mnt/e/CFD-Learning-Lab/temp
cd /mnt/e/CFD-Learning-Lab/native/runs/showcase-heated-jet-energy-v001/openfoam-qualification-v001
test ! -e solver.log
checkMesh >mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
export PYTHONPATH='E:\CFD-Learning-Lab\local\showcase-improvements-2026-11-08-v001\heated-jet-energy-v001;E:\CFD-Learning-Lab\repo\cfd-learning-lab\scripts'
/mnt/e/CFD-Learning-Lab/environments/fluent/Scripts/python.exe -B E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/heated-jet-energy-v001/jet_energy_prelaunch_v003.py >immediate-pre-mpi-admission.json
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-heated-jet-energy-qualification-v001
mkdir "$OMPI_MCA_orte_tmpdir_base"
printf 'running20thermalsteps\n' >run-status.txt
mpirun -np 2 foamRun -parallel >solver.log 2>&1
reconstructPar -time '0.005:0.1' >reconstruct.log 2>&1
printf 'complete\n' >run-status.txt
