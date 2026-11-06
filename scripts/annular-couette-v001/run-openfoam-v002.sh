#!/usr/bin/env bash
# Annular Couette: root-admitted mesh/BC preflight and a separately admitted FIVE-step pilot only.
set -e -o pipefail
root=/mnt/e/CFD-Learning-Lab
case=$root/native/runs/2026-11-02-annular-couette-v001/openfoam-pilot-v002
prep=$root/repo/cfd-learning-lab/scripts/prepare_openfoam_annular_couette_v002.py
test "$#" -eq 2; mode=$1; authorization=$2; test -s "$authorization"
case "$mode" in --prepare-only|--pilot-run) ;; *) printf 'Unsupported phase\n' >&2; exit 2;; esac
python3 -B - "$authorization" "$prep" "$0" "$mode" <<'PY'
from pathlib import Path
import hashlib,json,sys
a=json.load(open(sys.argv[1],encoding='utf-8-sig'));phase='PrepareOnly' if sys.argv[4]=='--prepare-only' else 'PilotRun'
assert a.get('approved') is True and a.get('approved_by')=='root' and a.get('action')=='run_annular_couette_openfoam_v002' and a.get('allowed_mode')==phase
for key,path in [('preparation_sha256',sys.argv[2]),('wrapper_sha256',sys.argv[3])]:assert a[key]==hashlib.sha256(Path(path).read_bytes()).hexdigest()
PY
export TMPDIR=$root/temp/annular-couette-v001 TEMP=$root/temp/annular-couette-v001 TMP=$root/temp/annular-couette-v001 XDG_CACHE_HOME=$root/temp/openfoam-cache WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$WM_PROJECT_USER_DIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin" FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib" FOAM_RUN="$WM_PROJECT_USER_DIR/run" OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
if pgrep -x foamRun >/dev/null || pgrep -x simpleFoam >/dev/null || pgrep -x pimpleFoam >/dev/null || pgrep -x fluent >/dev/null; then printf 'Another CFD solver active\n' >&2; exit 3; fi
test -d "$case"; cd "$case"
python3 -B "$prep" --case "$case" --verify-inputs > /dev/null
stage=verifying_prepared_case; solver_started=false
status() { printf '{"stage":"%s","mode":"%s","solver_started":%s,"mpi_ranks":2,"production_admitted":false}\n' "$stage" "$mode" "$solver_started" > run-status.json; }
trap 'code=$?; printf "{\"stage\":\"failed\",\"failed_stage\":\"%s\",\"exit_code\":%s,\"solver_started\":%s,\"production_admitted\":false}\n" "$stage" "$code" "$solver_started" > run-status.json; exit "$code"' ERR
if test "$mode" = --prepare-only; then
 test ! -e mesh-check.log; test ! -d 0/polyMesh; test ! -d constant/polyMesh; test ! -e solver.log
 cp "$0" provenance/run-openfoam-preflight-executed.sh; cp "$authorization" provenance/root-prepare-authorization.json
 stage=capturing_actual_installed_registration_reference_and_solver_sources; status
 python3 -B - <<'PY'
from pathlib import Path
import hashlib,json
p=Path('provenance/official-research-source-manifest.json');rows=json.loads(p.read_text());hashes=[]
for row in rows:
 source=Path('/opt/openfoam14')/row['path'];data=source.read_bytes();target=Path('provenance/installed-sources')/row['path'];target.parent.mkdir(parents=True,exist_ok=True)
 with target.open('xb') as f:f.write(data)
 hashes.append(hashlib.sha256(data).hexdigest()+'  '+str(source))
Path('provenance/installed-source-sha256.txt').write_text('\n'.join(hashes)+'\n')
s=Path('provenance/installed-sources/src/finiteVolume/fields/fvPatchFields/derived/rotatingWallVelocity/rotatingWallVelocityFvPatchVectorField.C').read_text();assert 'updateCoeffs();' in s and 'Up - n*(n & Up)' in s and 'makePatchTypeField' in s
s=Path('provenance/installed-sources/applications/modules/incompressibleFluid/correctPressure.C').read_text();assert 'constrainPressure(p, U, phiHbyA' in s and 'pressureReference.refCell()' in s
omega=Path('provenance/installed-sources/src/OpenFOAM/primitives/functions/Function1/omega1/omega1.C').read_text();assert 'units::radians/dimTime' in omega
s=Path('provenance/installed-sources/src/finiteVolume/cfdTools/general/findRefCell/findRefCell.C').read_text();assert 'sumHasRef != 1' in s and 'refPointName' in s and 'Pstream::master()' in s
PY
 stage=converting_selected_2D_mesh_to_actual_empty_extrusion; status
 fluentMeshToFoam -2D 0.01 input-annular-couette.msh > mesh-conversion.log 2>&1
 stage=full_native_mesh_quality; status
 checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
 stage=actual_serial_BC_registration_initial_rest_and_geometry_without_timesteps; status
 foamPostProcess -solver incompressibleFluid -dict system/preflightDiagnostics -time 0 > initial-native-registration.log 2>&1
 for field in U p phi C Vc vorticity; do test -s "0/$field"; done
 stage=decomposing_actual_initial_native_fields; status
 decomposePar -time 0 > decomposition.log 2>&1
 stage=parallel_closed_pressure_reference_and_rotating_wall_readback_without_timesteps; status
 mpirun -np 2 foamPostProcess -parallel -solver incompressibleFluid -dict system/preflightDiagnostics -time 0 > parallel-native-registration.log 2>&1
 python3 -B - <<'PY'
from pathlib import Path
import re,hashlib,json
proof={}
for label,name,ranks in [('serial','initial-native-registration.log',1),('parallel','parallel-native-registration.log',2)]:
 p=Path(name);text=p.read_text();assert re.search(r'^End\s*$',text,re.M) and 'FOAM FATAL' not in text
 if ranks==2:assert re.search(r'\bnProcs\s*:\s*2\b',text)
 proof[label]=dict(return_code=0,mpi_ranks=ranks,log=name,log_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),written_after_successful_command_return=True)
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
p=Path('provenance/preflight-validation.json');a=json.load(open(sys.argv[1],encoding='utf-8-sig'));assert json.loads(p.read_text())['passed'] is True and a['preflight_receipt_sha256']==hashlib.sha256(p.read_bytes()).hexdigest()
PY
cp "$0" provenance/run-openfoam-pilot-executed.sh; cp "$authorization" provenance/root-pilot-authorization.json
sha256sum -c provenance/installed-source-sha256.txt > installed-source-verification.log
stage=solving_EXACTLY_five_actual_Euler_steps; solver_started=true; status
mpirun -np 2 foamRun -parallel > solver.log 2>&1
python3 -B - <<'PY'
from pathlib import Path
import re
s=Path('solver.log').read_text();assert re.search(r'^End\s*$',s,re.M) and 'Finalising parallel run' in s
actual=[float(t) for t in re.findall(r'^Time = ([0-9.eE+-]+)s?\s*$',s,re.M)];assert len(actual)==5 and all(abs(t-i*.05)<1e-12 for i,t in enumerate(actual,1))
for rank in (0,1):
 for i in range(1,6):
  folder=Path(f'processor{rank}')/f'{i*.05:.10g}'
  for n in ('U','p','phi','vorticity'):assert (folder/n).is_file()
  assert (folder/'uniform/time').is_file()
PY
printf '{"stage":"finished_pilot_execution","solver_started":true,"mpi_ranks":2,"actual_steps":5,"physical_time_s":0.25,"deltaT_s":0.05,"native_rank_states":6,"serial_reconstruction_requested":false,"production_admitted":false,"numerical_analysis_pending":true}\n' > run-status.json
