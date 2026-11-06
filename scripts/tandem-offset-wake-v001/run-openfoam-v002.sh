#!/usr/bin/env bash
# Offset tandem wake: root-admitted native preflight and separately admitted FIVE-step pilot.
set -e -o pipefail
root=/mnt/e/CFD-Learning-Lab
case=$root/native/runs/2026-11-02-tandem-offset-wake-v001/openfoam-pilot-v002
prep=$root/repo/cfd-learning-lab/scripts/prepare_openfoam_tandem_offset_wake_v002.py
common=$root/repo/cfd-learning-lab/scripts/tandem_offset_openfoam_common_v002.py
test "$#" -eq 2; mode=$1; authorization=$2; test -s "$authorization"
case "$mode" in --prepare-only|--pilot-run) ;; *) exit 2;; esac
python3 -B - "$authorization" "$prep" "$common" "$0" "$mode" "$case" <<'PY'
from pathlib import Path
import hashlib,json,sys
h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();a=json.load(open(sys.argv[1],encoding='utf-8-sig'));phase='PrepareOnly' if sys.argv[5]=='--prepare-only' else 'PilotRun'
assert a.get('approved') is True and a.get('approved_by')=='root' and a.get('action')=='run_tandem_offset_wake_openfoam_v002' and a.get('allowed_mode')==phase
for key,path in [('preparation_sha256',sys.argv[2]),('common_sha256',sys.argv[3]),('wrapper_sha256',sys.argv[4])]:assert a[key]==h(path)
assert a['prepared_input_manifest_sha256']==h(Path(sys.argv[6])/'provenance/prepared-input-sha256.json')
PY
export TMPDIR=$root/temp/tandem-offset-wake-v001 TEMP=$root/temp/tandem-offset-wake-v001 TMP=$root/temp/tandem-offset-wake-v001 XDG_CACHE_HOME=$root/temp/openfoam-cache WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$WM_PROJECT_USER_DIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin" FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib" FOAM_RUN="$WM_PROJECT_USER_DIR/run" OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
if pgrep -x foamRun >/dev/null || pgrep -x fluent >/dev/null || pgrep -x simpleFoam >/dev/null || pgrep -x pimpleFoam >/dev/null; then printf 'Another solver active\n' >&2; exit 3; fi
cd "$case"
python3 -B "$prep" --case "$case" --verify-inputs >/dev/null
stage=verifying_prepared_case; solver_started=false
status() { printf '{"stage":"%s","mode":"%s","solver_started":%s,"mpi_ranks":2,"production_admitted":false}\n' "$stage" "$mode" "$solver_started" > run-status.json; }
trap 'code=$?; printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s,\"solver_started\":%s,\"production_admitted\":false}\n" "$stage" "$code" "$solver_started" > run-status.json; exit "$code"' ERR
if test "$mode" = --prepare-only; then
 test ! -e mesh-check.log; test ! -d 0/polyMesh; test ! -d constant/polyMesh; test ! -e solver.log
 cp "$0" provenance/run-openfoam-preflight-executed.sh; cp "$authorization" provenance/root-prepare-authorization.json
 stage=capturing_actual_installed_solver_force_and_derived_startup_sources; status
 python3 -B - <<'PY'
from pathlib import Path
import json,hashlib
rows=json.loads(Path('provenance/installed-source-list.json').read_text());hashes=[]
for name in rows:
 source=Path('/opt/openfoam14')/name;data=source.read_bytes();target=Path('provenance/installed-sources')/name;target.parent.mkdir(parents=True,exist_ok=True)
 with target.open('xb') as f:f.write(data)
 hashes.append(hashlib.sha256(data).hexdigest()+'  '+str(source))
Path('provenance/installed-source-sha256.txt').write_text('\n'.join(hashes)+'\n')
assert 'executeAtStart_(true)' in Path('provenance/installed-sources/src/OpenFOAM/db/functionObjects/functionObject/functionObject.C').read_text()
assert 'rhoInf' in Path('provenance/installed-sources/src/functionObjects/forces/forcesBase/forcesBase.C').read_text()
assert 'Aref' in Path('provenance/installed-sources/src/functionObjects/forces/forceCoeffs/forceCoeffs.C').read_text()
assert 'devSigma' in Path('provenance/installed-sources/src/functionObjects/field/wallShearStress/wallShearStress.C').read_text()
assert 'writeEntry(os, "value", *this)' in Path('provenance/installed-sources/src/finiteVolume/fields/fvsPatchFields/constraint/symmetry/symmetryFvsPatchField.C').read_text()
PY
 stage=converting_actual_shared_two_hole_mesh; status
 fluentMeshToFoam -2D 0.01 input-tandem-offset-wake.msh > mesh-conversion.log 2>&1
 stage=full_checkMesh; status
 checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
 stage=serial_initial_fields_geometry_and_native_derived_registration_without_steps; status
 foamPostProcess -solver incompressibleFluid -dict system/initialDiagnostics -time 0 > initial-native-registration.log 2>&1
 for n in U p phi C Vc vorticity wallShearStress; do test -s "0/$n"; done
 stage=decomposing_actual_initial_native_state; status
 decomposePar -time 0 > decomposition.log 2>&1
 stage=parallel_force_registration_and_initial_fields_without_steps; status
 mpirun -np 2 foamPostProcess -parallel -solver incompressibleFluid -dict system/forcePreflight -time 0 > parallel-native-registration.log 2>&1
 python3 -B - <<'PY'
from pathlib import Path
import json,re,hashlib
proof={}
for label,log,ranks in [('serial','initial-native-registration.log',1),('parallel','parallel-native-registration.log',2)]:
 s=Path(log).read_text();assert re.search(r'^End\s*$',s,re.M) and 'FOAM FATAL' not in s
 if ranks==2:assert re.search(r'\bnProcs\s*:\s*2\b',s)
 proof[label]=dict(return_code=0,mpi_ranks=ranks,log=log,log_sha256=hashlib.sha256(Path(log).read_bytes()).hexdigest(),written_after_successful_command_return=True)
with Path('provenance/native-registration-process-status.json').open('x') as f:json.dump(proof,f,indent=2)
PY
 python3 -B "$prep" --case "$case" --validate-native > preflight-validation.stdout.json
 printf '{"stage":"native_preflight_passed","solver_started":false,"mpi_ranks":2,"production_admitted":false}\n' > run-status.json
 exit 0
fi
test -s provenance/preflight-validation.json; test ! -e solver.log
python3 -B - "$authorization" <<'PY'
from pathlib import Path
import hashlib,json,sys
h=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();a=json.load(open(sys.argv[1],encoding='utf-8-sig'));r=json.loads(Path('provenance/preflight-validation.json').read_text());assert r['passed'] is True and a['preflight_receipt_sha256']==h('provenance/preflight-validation.json')
assert r['primary_initial_manifest_sha256']==h('provenance/native-initial-primary-sha256.json') and r['derived_initial_manifest_sha256']==h('provenance/native-initial-derived-sha256.json')
for p,digest in json.loads(Path('provenance/native-initial-primary-sha256.json').read_text()).items():assert h(p)==digest
PY
cp "$0" provenance/run-openfoam-pilot-executed.sh; cp "$authorization" provenance/root-pilot-authorization.json
sha256sum -c provenance/installed-source-sha256.txt > installed-source-verification.log
stage=solving_exactly_five_Euler_steps_30_outer; solver_started=true; status
mpirun -np 2 foamRun -parallel > solver.log 2>&1
stage=checking_native_completion_primary_preservation_and_derived_diagnostics; status
python3 -B - "$common" <<'PY'
from pathlib import Path
import re,hashlib,json,sys
sys.path.insert(0,str(Path(sys.argv[1]).parent));import tandem_offset_openfoam_common_v002 as c
h=c.sha;s=Path('solver.log').read_text();times=[float(t) for t in re.findall(r'^Time = ([0-9.eE+-]+)s?\s*$',s,re.M)];assert len(times)==5 and all(abs(t-i*.05)<1e-10 for i,t in enumerate(times,1)) and re.search(r'^End\s*$',s,re.M) and 'Finalising parallel run' in s
for rank in (0,1):
 for i in range(1,6):
  folder=Path(f'processor{rank}')/f'{i*.05:.10g}'
  for n in ('U','p','phi','vorticity','wallShearStress'):assert (folder/n).is_file()
  raw=(folder/'uniform/time').read_text()
  for name,value in [('value',i*.05),('deltaT',.05),('index',i)]:assert abs(float(re.search(r'\b'+name+r'\s+([-+0-9.eE]+)\s*;',raw)[1])-value)<1e-10
for p,digest in json.loads(Path('provenance/native-initial-primary-sha256.json').read_text()).items():assert h(p)==digest
rows=[]
for p,r in json.loads(Path('provenance/native-initial-derived-sha256.json').read_text()).items():
 current=h(p);baseline=Path(r['baseline_path']);assert h(baseline)==r['sha256'];x=c.field(baseline.read_text(),r['cells'],3,r['dimensions']);y=c.field(Path(p).read_text(),r['cells'],3,r['dimensions']);delta=max((abs(a-b) for v,w in zip(x,y) for a,b in zip(v,w)),default=0.);rows.append(dict(path=p,original_sha256=r['sha256'],current_sha256=current,byte_preserved=current==r['sha256'],native_internal_max_abs_difference=delta))
c.json_new(Path('provenance/derived-startup-diagnostic.json'),dict(primary_initial_inputs_unchanged=True,strict_all_derived_copies_preserved=all(r['byte_preserved'] for r in rows),derived_outputs=rows,CFD_replayed=False))
PY
printf '{"stage":"finished_pilot_execution","solver_started":true,"mpi_ranks":2,"actual_steps":5,"physical_time_s":0.25,"deltaT_s":0.05,"native_rank_states":6,"outer_correctors":30,"serial_reconstruction_requested":false,"production_admitted":false,"numerical_analysis_pending":true}\n' > run-status.json
