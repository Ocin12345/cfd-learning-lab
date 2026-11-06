#!/usr/bin/env bash
# SOURCE ONLY until exact root admission. Restart accepted pulse5;795 new two-rank scalar steps.
set -e -o pipefail
root=/mnt/e/CFD-Learning-Lab
private=$root/local/prospective-tracer-v001
prep=$private/prepare_openfoam_tracer_production_v001.py
case_dir=$root/native/runs/2026-11-02-tracer-tjunction-v001/openfoam-production-v001
test "$#" -eq 1
authorization=$1
test -s "$authorization"; test ! -e "$case_dir"
export TMPDIR=$root/temp/prospective-tracer-v001 TEMP=$root/temp/prospective-tracer-v001 TMP=$root/temp/prospective-tracer-v001
export XDG_CACHE_HOME=$root/temp/openfoam-cache WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$TMPDIR" "$XDG_CACHE_HOME" "$WM_PROJECT_USER_DIR"
source /opt/openfoam14/etc/bashrc
export WM_PROJECT_USER_DIR=$root/environments/openfoam-user
export FOAM_USER_APPBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/bin" FOAM_USER_LIBBIN="$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib" FOAM_RUN="$WM_PROJECT_USER_DIR/run"
export OMPI_MCA_orte_tmpdir_base="$TMPDIR"
test "$WM_PROJECT_VERSION" = 14
if pgrep -x foamRun >/dev/null || pgrep -x simpleFoam >/dev/null || pgrep -x pimpleFoam >/dev/null || pgrep -x fluent >/dev/null; then printf 'Another CFD solver active; refuse.\n' >&2; exit 3; fi
python3 -B "$prep" --authorization "$authorization" --verify-sources
stage=preparing_fresh_795_step_pulse_continuation
solver_started=false
status_file=$private/openfoam-production-stage-status-v001.json
test ! -e "$status_file"
status() { printf '{"stage":"%s","processes":2,"solver_started":%s,"pilot_steps":5,"additional_steps_authorized":795,"total_steps":800,"full_signed_scalar_budget_analysis_pending":true,"physical_validation":false}\n' "$stage" "$solver_started" > "$status_file"; }
finish() {
 code=$?
 trap - EXIT
 set +e
 if [ "$code" -ne 0 ]; then printf '{"stage":"failed","failed_stage":"%s","exit_code":%s,"solver_started":%s,"full_signed_scalar_budget_analysis_pending":true,"physical_validation":false}\n' "$stage" "$code" "$solver_started" > "$status_file"; fi
 # All synchronous native commands have returned before sync and any later dedicated VM close.
 sync
 sync_code=$?
 if [ -d "$case_dir" ]; then
  printf '{"native_commands_returned":true,"explicit_sync_exit_code":%s,"wrapper_exit_code":%s,"dedicated_VM_termination_executed":false}\n' "$sync_code" "$code" > "$case_dir/production-sync-completed.json"
 fi
 if [ "$code" -eq 0 ] && [ "$sync_code" -ne 0 ]; then code=$sync_code; fi
 exit "$code"
}
trap finish EXIT
status
python3 -B "$prep" --authorization "$authorization" --prepare > "$private/openfoam-production-preparation-v001.stdout.json"
cd "$case_dir"; status_file=$case_dir/run-status.json
stage=checking_all19_current_installed_source_hashes; status
python3 -B - "$prep" "$authorization" <<'PY'
import importlib.util,pathlib,sys
spec=importlib.util.spec_from_file_location('production_prep',sys.argv[1]);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
a=m.readj(sys.argv[2])
for name,pin in a['installed_source_sha256'].items():
 assert m.sha(pathlib.Path('/opt/openfoam14')/name)==m.sha(pathlib.Path('provenance/installed-sources')/name)==pin,'Current installed source changed: '+name
h=pathlib.Path('provenance/installed-sources/src/functionObjects/solvers/scalarTransport/scalarTransport.H').read_text()
assert m.re.search(r'bool\s+executeAtStart\(\)\s+const\s*\{\s*return false;',h),'No duplicate initial scalar solve permitted'
PY
stage=checking_actual_native_mesh_and_copied_accepted_restart; status
checkMesh -allGeometry -allTopology > production-mesh-check.log 2>&1
python3 -B "$prep" --authorization "$authorization" --verify-baseline > production-baseline-verification.json
# Reuse the accepted saved processor0/1 fields and uniform clocks at0.125. No decompose, initialization or patch.
for rank in 0 1; do
 for field in tracer U p phi; do test -s "processor$rank/0.125/$field"; done
 test -s "processor$rank/0.125/uniform/time"
 for field in U p phi; do test -s "processor$rank/2000/$field"; done
done
stage=solving_exactly795_new_scalar_steps6_through800; solver_started=true; status
mpirun -np 2 foamRun -parallel > production-solver.log 2>&1
stage=reconstructing_exactly795_new_native_states; status
reconstructPar -time '0.15:20' > production-reconstruction.log 2>&1
stage=verifying_all801_native_states_rank_clocks_frozen_fields_and_pilot_bytes; status
python3 -B "$prep" --authorization "$authorization" --verify-completed > production-execution-checks.json
python3 -B - "$prep" "$authorization" <<'PY'
import importlib.util,pathlib,sys
spec=importlib.util.spec_from_file_location('production_prep',sys.argv[1]);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
a=m.readj(sys.argv[2])
for name,pin in a['installed_source_sha256'].items():assert m.sha(pathlib.Path('/opt/openfoam14')/name)==m.sha(pathlib.Path('provenance/installed-sources')/name)==pin
PY
stage=finished_795_step_production_execution; status
printf '{"stage":"finished_795_step_production_execution","processes":2,"solver_started":true,"pilot_steps":5,"additional_steps":795,"total_steps":800,"native_states":801,"deltaT_s":0.025,"physical_time_s":20,"accepted_pilot_all_file_bytes_preserved":true,"actual_step80_top_c":1,"actual_step81_top_c":0,"full_signed_scalar_budget_analysis_pending":true,"root_acceptance_performed":false,"physical_validation":false}\n' > run-status.json
# EXIT trap sync is deliberately before root may close the dedicated CFD distro.
