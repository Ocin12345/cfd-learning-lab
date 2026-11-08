#!/bin/bash
set -euo pipefail
source /opt/openfoam14/etc/bashrc
export TMPDIR=/mnt/e/CFD-Learning-Lab/temp
cd /mnt/e/CFD-Learning-Lab/native/runs/showcase-heated-jet-energy-v001/openfoam-qualification-v001
test ! -e solver.log
checkMesh >mesh-check.log 2>&1
grep -q 'Mesh OK.' mesh-check.log
export PYTHONPATH='E:\CFD-Learning-Lab\local\showcase-improvements-2026-11-08-v001\heated-jet-energy-v001;E:\CFD-Learning-Lab\repo\cfd-learning-lab\scripts'
/mnt/e/CFD-Learning-Lab/environments/fluent/Scripts/python.exe -B E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/heated-jet-energy-v001/jet_energy_prelaunch_v002.py >immediate-pre-mpi-admission.json
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-heated-jet-energy-qualification-v001
mkdir "$OMPI_MCA_orte_tmpdir_base"
printf 'running20thermalsteps\n' >run-status.txt
mpirun -np 2 foamRun -parallel >solver.log 2>&1
reconstructPar -time '0.005:0.1' >reconstruct.log 2>&1
printf 'complete\n' >run-status.txt
