#!/usr/bin/env bash
set -eo pipefail
cd /mnt/e/CFD-Learning-Lab/native/runs/2026-11-06-double-lid-cavity-v001/openfoam
export TMPDIR=/mnt/e/CFD-Learning-Lab/temp/double-lid-cavity-v001-openfoam
export TEMP="$TMPDIR" TMP="$TMPDIR" XDG_CACHE_HOME=/mnt/e/CFD-Learning-Lab/temp/openfoam-cache
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$TMPDIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=/mnt/e/CFD-Learning-Lab/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin"
export FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib"
export FOAM_RUN="$WM_PROJECT_USER_DIR/run"
export OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
test ! -e solver.log
test ! -e processor0
test ! -e processor1
if pgrep -x 'foamRun|reconstructPar|simpleFoam|pimpleFoam' >/dev/null; then exit 3; fi
python3 -B - <<'PY'
from pathlib import Path
import hashlib,json
m=json.loads(Path("input-manifest.json").read_text())
for n,h in m["prepared_hashes"].items():
    assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==h, n
PY
stage=mesh_check
trap 'code=$?; printf "{\"stage\":\"%s\",\"exit_code\":%s}\n" "$stage" "$code" > run-status.json; exit "$code"' ERR
foamVersion > installation-version.txt
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
grep -q 'Mesh OK' mesh-check.log
stage=decomposition
decomposePar > decomposition.log 2>&1
stage=solving
mpirun -np 2 foamRun -parallel > solver.log 2>&1
grep -Eq '^Time = 2000(s)?[[:space:]]*$' solver.log
grep -q '^End' solver.log
stage=reconstructing
reconstructPar -time '1800,2000' > reconstruction.log 2>&1
for iter in 1800 2000; do
  for fld in U p phi; do test -s "$iter/$fld"; done
done
stage=writing_centres
foamPostProcess -func writeCellCentres -time 2000 > cell-centres.log 2>&1
test -s 2000/C
sync
printf '{"stage":"finished_execution","processes":2,"final_iteration":2000,"numerical_analysis_pending":true}\n' > run-status.json
