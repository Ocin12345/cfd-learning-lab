#!/usr/bin/env bash
# SOURCE DRAFT. One separately admitted five-step qualification; no production mode.
set -e -o pipefail
root=/mnt/e/CFD-Learning-Lab
private=$root/local/prospective-tracer-v001
prep=$private/prepare_openfoam_tracer_qualification_v004.py
test "$#" -eq 2
phase=$1; authorization=$2
case "$phase" in
 Control) case_dir=$root/native/runs/2026-11-02-tracer-tjunction-v001/openfoam-control-v003 ;;
 Pulse) case_dir=$root/native/runs/2026-11-02-tracer-tjunction-v001/openfoam-pulse-pilot-v003 ;;
 *) printf 'Only separately authorized Control or Pulse allowed.\n' >&2; exit 2 ;;
esac
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
python3 -B "$prep" --phase "$phase" --authorization "$authorization" --verify-sources
stage=preparing_fresh_five_step_qualification
solver_started=false
status_file=$private/openfoam-${phase,,}-stage-status-v004.json
test ! -e "$status_file"
status() { printf '{"stage":"%s","phase":"%s","processes":2,"solver_started":%s,"actual_steps_authorized":5,"production_admitted":false}\n' "$stage" "$phase" "$solver_started" > "$status_file"; }
finish() {
 code=$?
 trap - EXIT
 if [ "$code" -ne 0 ]; then printf '{"stage":"failed","failed_stage":"%s","exit_code":%s,"phase":"%s","solver_started":%s,"production_admitted":false}\n' "$stage" "$code" "$phase" "$solver_started" > "$status_file"; fi
 # Keep all complete saved states and flush E: before any later dedicated distro close.
 sync
 exit "$code"
}
trap finish EXIT
status
python3 -B "$prep" --phase "$phase" --authorization "$authorization" --prepare > "$private/openfoam-${phase,,}-preparation-v004.stdout.json"
cd "$case_dir"; status_file=$case_dir/run-status.json
stage=capturing_exact_reviewed_current_installed_sources; status
python3 -B - "$prep" "$authorization" <<'PY'
import importlib.util,pathlib,sys
p=pathlib.Path(sys.argv[1]); spec=importlib.util.spec_from_file_location('tracer_prep',p); m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
a=m.readj(sys.argv[2]); rows={}
for name in m.INSTALLED:
 src=pathlib.Path('/opt/openfoam14')/name; dst=pathlib.Path('provenance/installed-sources')/name
 assert m.sha(src)==a['installed_source_sha256'][name], 'Current installed source was not root-reviewed: '+name
 m.copy_new(src,dst); rows[name]=m.sha(dst)
m.json_new(pathlib.Path('provenance/installed-source-sha256.json'),rows)
h=pathlib.Path('provenance/installed-sources/src/functionObjects/solvers/scalarTransport/scalarTransport.H').read_text()
assert m.re.search(r'bool\s+executeAtStart\(\)\s+const\s*\{\s*return false;',h), 'Scalar startup must not solve a duplicate step'
f=pathlib.Path('provenance/installed-sources/src/finiteVolume/fvMatrices/fvMatrix/fvMatrixSolve.C').read_text()
assert 'word(name + "Final")' in f and 'solutionControl::finalIteration' in f
PY
stage=checking_actual_native_mesh; status
checkMesh -allGeometry -allTopology > mesh-check.log 2>&1
stage=writing_native_C_Vc_and_frozen_baseline_without_scalar_solve; status
foamPostProcess -fields '(tracer U p phi)' -dict system/baselineDiagnostics -time 0 > baseline-diagnostics.log 2>&1
for field in tracer U p phi C Vc; do test -s "0/$field"; done
python3 -B "$prep" --phase "$phase" --authorization "$authorization" --verify-baseline > baseline-verification.json
stage=decomposing_physical0_and_frozen_read2000; status
decomposePar -time '0,2000' > decomposition.log 2>&1
for rank in 0 1; do
 for field in tracer U p phi; do test -s "processor$rank/0/$field"; done
 for field in U p phi; do test -s "processor$rank/2000/$field"; done
 for address in cellProcAddressing faceProcAddressing pointProcAddressing; do test -s "processor$rank/0/polyMesh/$address"; done
done
stage=solving_exactly_five_scalar_steps_without_flow_updates; solver_started=true; status
mpirun -np 2 foamRun -parallel > solver.log 2>&1
stage=reconstructing_actual_five_solved_states; status
reconstructPar -time '0.025:0.125' > reconstruction.log 2>&1
stage=verifying_actual_clocks_scalar_fields_frozen_allface_fluxes_and_originals; status
python3 -B "$prep" --phase "$phase" --authorization "$authorization" --verify-completed > qualification-execution-checks.json
python3 -B - "$prep" "$authorization" <<'PY'
import importlib.util,pathlib,sys
spec=importlib.util.spec_from_file_location('tracer_prep',sys.argv[1]);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
a=m.readj(sys.argv[2])
for name,digest in a['installed_source_sha256'].items(): assert m.sha(pathlib.Path('/opt/openfoam14')/name)==m.sha(pathlib.Path('provenance/installed-sources')/name)==digest
PY
stage=finished_five_step_qualification_execution; status
printf '{"stage":"finished_five_step_qualification_execution","phase":"%s","processes":2,"solver_started":true,"actual_steps":5,"native_states":6,"deltaT_s":0.025,"physical_time_s":0.125,"full_signed_scalar_budget_analysis_pending":true,"production_admitted":false,"physical_validation":false}\n' "$phase" > run-status.json
# EXIT trap sync runs only after all native commands and saved readbacks are complete.