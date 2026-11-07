#!/usr/bin/env bash
set -eo pipefail
root=/mnt/e/CFD-Learning-Lab
case_dir=$root/native/runs/finite-wing-refined-v001/openfoam
prep=$root/local/finite-wing-mesh-sensitivity-v001/prepare_openfoam_refined_v001.py
mode=${1:---prepare-only}
case "$mode" in --prepare-only|--run) ;; *) exit 2 ;; esac
cd "$case_dir"
export TMPDIR=$root/temp/wing-refined-openfoam TEMP=$root/temp/wing-refined-openfoam TMP=$root/temp/wing-refined-openfoam
export XDG_CACHE_HOME=$root/temp/wsl-cache PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export OMPI_MCA_orte_tmpdir_base=/dev/shm/cfd-wing-refined-v001
mkdir -p "$OMPI_MCA_orte_tmpdir_base"
test "$WM_PROJECT_VERSION" = 14
trap 'rc=$?; printf "failed (%s)\n" "$rc" > run-status.txt; exit "$rc"' ERR
if pgrep -x foamRun >/dev/null; then echo 'Another solver is active' >&2; exit 3; fi
if test "$mode" = --prepare-only; then
    test ! -e mesh-conversion.log
    test ! -e solver.log
    fluentMeshToFoam input-wing.msh > mesh-conversion.log 2>&1
    checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
    python3 -B "$prep" --case "$case_dir" --validate-native > native-preflight.stdout.json
    printf 'mesh_checked\n' > run-status.txt
    exit 0
fi
test ! -e solver.log
test -s provenance/preflight-validation.json
grep -q 'Mesh OK.' mesh-check.log
cp "$0" provenance/executed-wrapper.sh
printf 'decomposing\n' > run-status.txt
decomposePar > decomposition.log 2>&1
printf 'solving\n' > run-status.txt
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -q '^End' solver.log
grep -Eq '^Time = 1500(s)?[[:space:]]*$' solver.log
printf 'reconstructing\n' > run-status.txt
reconstructPar -time '300,600,900,1200,1500' > reconstruction.log 2>&1
for t in 300 600 900 1200 1500; do for f in U p phi; do test -s "$t/$f"; done; done
foamPostProcess -func writeCellCentres -latestTime > cell-centres.log 2>&1
test -s 1500/C
printf 'complete\n' > run-status.txt
