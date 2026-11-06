#!/usr/bin/env bash
# Fresh amended annular Couette pilot: consume proven native rest inputs; no conversion/registration replay.
set -e -o pipefail
root=/mnt/e/CFD-Learning-Lab
case=$root/native/runs/2026-11-02-annular-couette-v001/openfoam-pilot-v003
prep=$root/repo/cfd-learning-lab/scripts/prepare_openfoam_annular_couette_v003.py
test "$#" -eq 2; test "$1" = --pilot-run; authorization=$2; test -s "$authorization"
python3 -B - "$authorization" "$prep" "$0" "$case" <<'PY'
from pathlib import Path
import hashlib,json,sys
h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
a=json.load(open(sys.argv[1],encoding='utf-8-sig'));c=Path(sys.argv[4]);assert a.get('approved') is True and a.get('approved_by')=='root' and a.get('action')=='run_annular_couette_openfoam_v003' and a.get('allowed_mode')=='PilotRun'
for key,path in [('preparation_sha256',sys.argv[2]),('wrapper_sha256',sys.argv[3]),('prepared_input_manifest_sha256',c/'provenance/prepared-input-sha256.json'),('native_initial_copy_receipt_sha256',c/'provenance/native-initial-copy-receipt.json'),('preflight_receipt_sha256',c/'provenance/preflight-validation.json'),('amendment_sha256',c/'provenance/OPENFOAM-PRESSURE-COUPLING-AMENDMENT-v001.md')]:assert a[key]==h(path)
s=json.loads((c/'run-status.json').read_text());assert s['stage']=='prepared_native_initial_copy_30_outer' and s['solver_started'] is False
for prefix in ('','processor0','processor1'):
 assert sorted(p.name for p in (c/prefix).iterdir() if p.is_dir() and p.name.replace('.','',1).isdigit())==['0']
assert not (c/'solver.log').exists()
PY
export TMPDIR=$root/temp/annular-couette-v001 TEMP=$root/temp/annular-couette-v001 TMP=$root/temp/annular-couette-v001 XDG_CACHE_HOME=$root/temp/openfoam-cache WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$WM_PROJECT_USER_DIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin" FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib" FOAM_RUN="$WM_PROJECT_USER_DIR/run" OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
if pgrep -x foamRun >/dev/null || pgrep -x simpleFoam >/dev/null || pgrep -x pimpleFoam >/dev/null || pgrep -x fluent >/dev/null; then printf 'Another CFD solver active\n' >&2; exit 3; fi
cd "$case"
python3 -B "$prep" --case "$case" --verify-inputs > /dev/null
stage=verifying_native_initial_copy_and_amended_sources; solver_started=false
status() { printf '{"stage":"%s","mode":"PilotRun","solver_started":%s,"mpi_ranks":2,"outer_correctors":30,"production_admitted":false}\n' "$stage" "$solver_started" > run-status.json; }
trap 'code=$?; printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s,\"solver_started\":%s,\"production_admitted\":false}\n" "$stage" "$code" "$solver_started" > run-status.json; exit "$code"' ERR
test ! -e provenance/run-openfoam-pilot-executed.sh; test ! -e provenance/root-pilot-authorization.json
cp "$0" provenance/run-openfoam-pilot-executed.sh; cp "$authorization" provenance/root-pilot-authorization.json
sha256sum -c provenance/installed-source-sha256.txt > installed-source-verification.log
stage=solving_EXACTLY_five_actual_Euler_steps_30_outer; solver_started=true; status
mpirun -np 2 foamRun -parallel > solver.log 2>&1
python3 -B - <<'PY'
from pathlib import Path
import re
s=Path('solver.log').read_text();assert re.search(r'^End\s*$',s,re.M) and 'Finalising parallel run' in s and 'FOAM FATAL' not in s
actual=[float(t) for t in re.findall(r'^Time = ([0-9.eE+-]+)s?\s*$',s,re.M)];assert len(actual)==5 and all(abs(t-i*.05)<1e-12 for i,t in enumerate(actual,1))
for rank in (0,1):
 for i in range(1,6):
  folder=Path(f'processor{rank}')/f'{i*.05:.10g}'
  for n in ('U','p','phi','vorticity'):assert (folder/n).is_file()
  assert (folder/'uniform/time').is_file()
PY
printf '{"stage":"finished_pilot_execution","solver_started":true,"mpi_ranks":2,"actual_steps":5,"physical_time_s":0.25,"deltaT_s":0.05,"native_rank_states":6,"outer_correctors":30,"ancestral_10_outer_numeric_miss_preserved":true,"serial_reconstruction_requested":false,"production_admitted":false,"numerical_analysis_pending":true}\n' > run-status.json
