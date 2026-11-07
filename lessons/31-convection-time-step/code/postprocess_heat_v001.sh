#!/usr/bin/env bash
set -eo pipefail
cd /mnt/e/CFD-Learning-Lab/native/runs/natural-convection-dt001-v001/openfoam
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TMP=$root/temp TEMP=$root/temp XDG_CACHE_HOME=$root/temp/wsl-cache OMP_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
test "$(cat run-status.txt)" = complete
test ! -e wall-heat-postprocess.log
foamPostProcess -solver fluid -func wallHeatFlux -latestTime > wall-heat-postprocess.log 2>&1
test -s 20/wallHeatFlux