import os,json,sys,time,traceback
import wing_import_paths_v001
import ansys.fluent.core as pf
from wing_repair_common_v001 import *
from wing_simplec_prelaunch_v001 import wait_admission,verify
from wing_repair_native_v003 import state,state_keys,raw_residual,same_native
from wing_repair_helpers_v001 import original_target,same_residual_prefix
from wing_simplec_capability_helpers_v001 import validate_state,diff_inventory,assert_zero_advance_source
from wing_fluent_v001 import mesh_quality,hdf_terminal
from sloshing_fixed_settings_v001 import compare_fixed
OUT=N/'fluent-simplec-capability-v001'

def iteration_query(st,label):
 expressions=st.setup.named_expressions;before=expressions.get_state();name='codex-capability-current-iteration';assert name not in expressions.get_object_names()
 try:
  expressions.create(name=name);expressions[name].definition='Iteration';value=expressions[name].get_value();dump(OUT/(label+'-live-iteration-raw.json'),dict(expression='Iteration',query_value=value));number=float(value);assert number==1500
 finally:
  if name in expressions.get_object_names():del expressions[name]
 compare_fixed(expressions.get_state(),before);return int(number)

def main():
 assert_zero_advance_source(Path(__file__).read_text(encoding='utf-8-sig'));assert not (OUT/'solver-transcript.txt').exists();gate=wait_admission(lambda **v:mark(OUT/'status.json','waiting_admission',**v));dump(OUT/'launch-admission.json',gate);session=None;loaded=False
 try:
  pf.config.examples_path=str(R/'temp/fluent-examples');(R/'temp/fluent-examples').mkdir(exist_ok=True)
  session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=3,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True);dump(OUT/'solver-version.json',dict(version=str(session.get_fluent_version()),processor_count=2,connected_before_any_advance=True));session.transcript.start(file_name=str(OUT/'solver-transcript.txt'));st=session.settings
  mark(OUT/'status.json','reading_exact_original1500');st.file.read_case_data(file_name=str(OUT/'source.cas.h5'));loaded=True
  reports=st.solution.monitor.report_files
  for name in reports.get_object_names():reports[name].file_name=str(OUT/(name+'-history.out'))
  actual=state(st);dump(OUT/'loaded-source-settings.json',actual);base=json.loads((OUT/'expected-source-settings.json').read_text());compare_fixed(actual,{k:base[k] for k in state_keys()});assert not st.setup.models.energy.enabled()
  targets={k:original_target(st.solution.monitor.residual.equations[k]) for k in EQ};dump(OUT/'loaded-targets.json',targets);assert all(v==1e-5 for v in targets.values());compare_fixed(state(st),actual)
  live_before=iteration_query(st,'loaded');st.file.write_case_data(file_name=str(OUT/'loaded-restart.cas.h5'));assert hdf_terminal(OUT/'loaded-restart.dat.h5')[0]==1500;dump(OUT/'loaded-field-proof.json',same_native(OUT/'loaded-restart.cas.h5',OUT/'loaded-restart.dat.h5',OUT/'source.cas.h5',OUT/'source.dat.h5'))
  st.mesh.check();st.mesh.quality()
  for attempt in range(10):
   try:q=mesh_quality((OUT/'solver-transcript.txt').read_text(errors='replace'));break
   except ValueError:
    if attempt==9:raise
    time.sleep(.1)
  assert q['minimum_orthogonal_quality']>=.10 and q['minimum_cell_volume_m3']>0;dump(OUT/'geometry-quality.json',q)
  methods=st.solution.methods;allowed=methods.p_v_coupling.flow_scheme.allowed_values();dump(OUT/'allowed-flow-schemes.json',allowed);assert 'SIMPLEC' in allowed;methods.p_v_coupling.flow_scheme='SIMPLEC';dump(OUT/'after-switch-defaults-settings.json',state(st))
  methods.pseudo_time_method.formulation.segregated_solver='off';methods.p_v_coupling.skewness_correction_itr_count=1
  ur=st.solution.controls.under_relaxation;names=ur.get_object_names();dump(OUT/'active-under-relaxation-names.json',names);assert {'pressure','mom'}.issubset(names);ur['pressure'].set_state(.7);ur['mom'].set_state(.7)
  mg=st.solution.controls.advanced.multi_grid;names=mg.mg_controls.get_object_names();dump(OUT/'active-AMG-names.json',names);assert 'pressure' in names;mg.mg_controls['pressure'].termination_criteria=.001;mg.amg_controls.options.verbosity=1;st.solution.controls.advanced.expert.save_cell_residuals=True
  candidate=state(st);dump(OUT/'candidate-settings.json',candidate);dump(OUT/'native-state-differences.json',diff_inventory(actual,candidate));validate_state(actual,candidate)
  targets_after={k:original_target(st.solution.monitor.residual.equations[k]) for k in EQ};dump(OUT/'candidate-targets.json',targets_after);assert targets_after==targets;compare_fixed(state(st),candidate)
  live_after=iteration_query(st,'candidate');assert live_before==live_after==1500;st.file.write_case_data(file_name=str(OUT/'candidate.cas.h5'));dump(OUT/'candidate-field-proof.json',same_native(OUT/'candidate.cas.h5',OUT/'candidate.dat.h5',OUT/'source.cas.h5',OUT/'source.dat.h5'))
  pi,pr=raw_residual(OUT/'source.dat.h5');ci,cr=raw_residual(OUT/'candidate.dat.h5');same_residual_prefix(pi,pr,ci,cr,1500);assert list(ci)==list(pi) and ci[-1]==1500
  dump(OUT/'final-settings.json',state(st));compare_fixed(state(st),candidate);dump(OUT/'zero-advance-proof.json',dict(live_loaded_iteration=live_before,live_candidate_iteration=live_after,native_current_iteration=int(ci[-1]),new_iterations=0,raw_history_exact=True,initial_fields_exact=True,source_contract_pass=True));mark(OUT/'status.json','configured_zero_advance',new_iterations=0,iteration=1500)
 except Exception:
  if loaded:
   try:st.file.write_case_data(file_name=str(OUT/'recovery.cas.h5'))
   except Exception:pass
  mark(OUT/'status.json','failed',error=traceback.format_exc());raise
 finally:
  if session:session.exit()
 verify();mark(OUT/'status.json','closed',iteration=1500,new_iterations=0)
if __name__=='__main__':main()
