#!/usr/bin/env bash
# Source-only until root admission. Existing saved-field postprocessing, never a CFD solver.
set -e -o pipefail
root=/mnt/e/CFD-Learning-Lab
export TMPDIR=$root/temp TEMP=$root/temp TMP=$root/temp XDG_CACHE_HOME=$root/temp/openfoam-cache
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin" FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib" FOAM_RUN="$WM_PROJECT_USER_DIR/run"
exec python3 -B "$root/local/prospective-tracer-v001/recover_tracer_reconstruction_v001.py" "$@"
