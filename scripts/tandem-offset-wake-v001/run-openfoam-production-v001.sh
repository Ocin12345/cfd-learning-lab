#!/usr/bin/env bash
# Accepted five-step pilot + exactly795 additional native rank steps; no reinitialization or serial reconstruction.
set -e -o pipefail
root=/mnt/e/CFD-Learning-Lab
case=$root/native/runs/2026-11-02-tandem-offset-wake-v001/openfoam-production-v001
prep=$root/repo/cfd-learning-lab/scripts/prepare_openfoam_tandem_offset_wake_production_v001.py
test "$#" -eq 2; test "$1" = --production-run; authorization=$2; test -s "$authorization"
python3 -B - "$authorization" "$prep" "$0" "$case" <<'PY'
from pathlib import Path
import hashlib,json,sys
h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();a=json.load(open(sys.argv[1],encoding='utf-8-sig'));c=Path(sys.argv[4])
assert a.get('approved') is True and a.get('approved_by')=='root' and a.get('action')=='run_tandem_offset_wake_production_openfoam_v001' and a.get('allowed_mode')=='ProductionRun'
for key,path in [('preparation_sha256',sys.argv[2]),('wrapper_sha256',sys.argv[3]),('prepared_primary_manifest_sha256',c/'provenance/prepared-primary-sha256.json'),('inherited_accepted_pilot_receipt_sha256',c/'provenance/inherited-accepted-pilot.json'),('accepted_pilot_selection_sha256',c/'provenance/accepted-pilot-selection.json')]:assert a[key]==h(path)
p=json.loads((c/'provenance/accepted-pilot-selection.json').read_text())
for key in ('checks','manifest','fields'):assert a['accepted_pilot_'+key+'_sha256']==p['pilot_'+key+'_sha256']
s=json.loads((c/'run-status.json').read_text());assert s['stage']=='prepared_accepted_pilot_continuation_only' and s['solver_started'] is False and s['inherited_actual_steps']==5
for rank in (0,1):
 times=sorted(float(p.name) for p in (c/f'processor{rank}').iterdir() if p.is_dir() and p.name.replace('.','',1).isdigit());assert len(times)==6 and all(abs(t-i*.05)<1e-12 for i,t in enumerate(times))
assert not (c/'solver.log').exists()
PY
export TMPDIR=$root/temp/tandem-offset-wake-v001 TEMP=$root/temp/tandem-offset-wake-v001 TMP=$root/temp/tandem-offset-wake-v001 XDG_CACHE_HOME=$root/temp/openfoam-cache WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$WM_PROJECT_USER_DIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin" FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib" FOAM_RUN="$WM_PROJECT_USER_DIR/run" OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
if pgrep -x foamRun >/dev/null || pgrep -x simpleFoam >/dev/null || pgrep -x pimpleFoam >/dev/null || pgrep -x fluent >/dev/null; then printf 'Another CFD solver active\n' >&2; exit 3; fi
cd "$case"
python3 -B "$prep" --case "$case" --verify-inputs >/dev/null
stage=verifying_accepted_pilot_continuation; solver_started=false
status() { printf '{"stage":"%s","mode":"ProductionRun","solver_started":%s,"mpi_ranks":2,"outer_correctors":30,"production_admitted":true}\n' "$stage" "$solver_started" > run-status.json; }
trap 'code=$?; printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s,\"solver_started\":%s,\"production_admitted\":true}\n" "$stage" "$code" "$solver_started" > run-status.json; exit "$code"' ERR
test ! -e provenance/run-openfoam-production-executed.sh; test ! -e provenance/root-production-authorization.json
cp "$0" provenance/run-openfoam-production-executed.sh; cp "$authorization" provenance/root-production-authorization.json
sha256sum -c provenance/installed-source-sha256.txt > installed-source-verification.log
stage=solving_EXACTLY_795_additional_Euler_steps; solver_started=true; status
mpirun -np 2 foamRun -parallel > solver.log 2>&1
stage=checking_End_795_new_steps_primary_inputs_and_explicit_derived_diagnostics; status
python3 -B "$prep" --case "$case" --finalize-native > native-production-completion.stdout.json
printf '{"stage":"finished_production_execution","solver_started":true,"mpi_ranks":2,"inherited_actual_steps":5,"additional_actual_steps":795,"total_actual_steps":800,"physical_time_s":40,"deltaT_s":0.05,"native_rank_states":801,"outer_correctors":30,"serial_reconstruction_requested":false,"production_admitted":true,"numerical_analysis_pending":true}\n' > run-status.json