import os,json,time,traceback
import wing_import_paths_v001
import ansys.fluent.core as pf
from wing_repair_common_v001 import *
from wing_simplec_qualification_prelaunch_v001 import wait_admission,verify
from wing_repair_native_v003 import state,state_keys,raw_residual,same_native
from wing_repair_helpers_v001 import original_target,same_residual_prefix
from wing_simplec_capability_helpers_v002 import validate_state,diff_inventory
from wing_simplec_qualification_helpers_v001 import checkpoints,scale_guard,query_iteration
from wing_runtime_residual_export_v001 import export_runtime_residuals
from wing_fluent_v001 import mesh_quality,hdf_terminal
from sloshing_fixed_settings_v001 import compare_fixed
from audit_wing_continuation_v001 import native
OUT=N/'fluent-simplec-qualification-v001'

def main():
 assert not (OUT/'solver-transcript.txt').exists();gate=wait_admission(lambda **v:mark(OUT/'status.json','waiting_admission',**v));dump(OUT/'launch-admission.json',gate);session=None;loaded=False
 try:
  pf.config.examples_path=str(R/'temp/fluent-examples');(R/'temp/fluent-examples').mkdir(exist_ok=True)
  session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=3,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True);dump(OUT/'solver-version.json',dict(version=str(session.get_fluent_version()),processor_count=2,connected_before_any_advance=True));session.transcript.start(file_name=str(OUT/'solver-transcript.txt'));st=session.settings
  mark(OUT/'status.json','reading_exact_original1500');st.file.read_case_data(file_name=str(OUT/'source.cas.h5'));loaded=True
  reports=st.solution.monitor.report_files
  for name in reports.get_object_names():reports[name].file_name=str(OUT/(name+'-history.out'))
  assert reports['forces'].active() and reports['forces'].frequency()==1
  actual=state(st);dump(OUT/'loaded-source-settings.json',actual);base=json.loads((OUT/'expected-source-settings.json').read_text());compare_fixed(actual,{k:base[k] for k in state_keys()});assert not st.setup.models.energy.enabled()
  targets={k:original_target(st.solution.monitor.residual.equations[k]) for k in EQ};dump(OUT/'loaded-targets.json',targets);assert all(v==1e-5 for v in targets.values());compare_fixed(state(st),actual)
  query_iteration(st,OUT,'loaded',1500);st.file.write_case_data(file_name=str(OUT/'loaded-restart.cas.h5'));assert hdf_terminal(OUT/'loaded-restart.dat.h5')[0]==1500;dump(OUT/'loaded-field-proof.json',same_native(OUT/'loaded-restart.cas.h5',OUT/'loaded-restart.dat.h5',OUT/'source.cas.h5',OUT/'source.dat.h5'))
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
  mg=st.solution.controls.advanced.multi_grid;names=mg.mg_controls.get_object_names();dump(OUT/'active-AMG-names.json',names);assert set(names)=={'pressure','x-mom','y-mom','z-mom'};mg.mg_controls['pressure'].termination_criteria=.001;mg.amg_controls.options.verbosity=1;st.solution.controls.advanced.expert.save_cell_residuals=True
  candidate=state(st);expected=json.loads((W/'accepted-native-SIMPLEC-settings-v001.json').read_text());dump(OUT/'candidate-settings.json',candidate);dump(OUT/'native-state-differences.json',diff_inventory(actual,candidate));validate_state(actual,candidate);compare_fixed(candidate,expected)
  ta={k:original_target(st.solution.monitor.residual.equations[k]) for k in EQ};dump(OUT/'candidate-targets.json',ta);assert ta==targets;compare_fixed(state(st),expected);query_iteration(st,OUT,'configured',1500)
  st.file.write_case_data(file_name=str(OUT/'configured-start.cas.h5'));dump(OUT/'configured-field-proof.json',same_native(OUT/'configured-start.cas.h5',OUT/'configured-start.dat.h5',OUT/'source.cas.h5',OUT/'source.dat.h5'))
  pi,pr=raw_residual(OUT/'source.dat.h5');ci,cr=raw_residual(OUT/'configured-start.dat.h5');same_residual_prefix(pi,pr,ci,cr,1500);assert list(ci)==list(pi) and ci[-1]==1500;scale=float(pr['continuity'][-1,1]);assert scale==2.2344914086503613
  dump(OUT/'declared-test.json',dict(parent=1500,new_iterations=100,endpoint=1600,checkpoints=checkpoints(),continuity_scale=scale,residual_limit=1e-5,last20_required=True,force_window=100,globalmass_native_checks_at_five_saved_states_only=True))
  verify()
  for end in checkpoints():
   mark(OUT/'status.json','iterating',target_iteration=end);st.solution.run_calculation.iterate(iter_count=20);stem=OUT/f'iter-{end}';st.file.write_case_data(file_name=str(stem)+'.cas.h5');it,res=hdf_terminal(Path(str(stem)+'.dat.h5'));assert it==end and max(res.values())<.1;query_iteration(st,OUT,f'actual{end}',end)
   ni,nr=raw_residual(Path(str(stem)+'.dat.h5'));same_residual_prefix(pi,pr,ni,nr,1500);scale_guard(ni,nr['continuity'],scale,end)
   snapshot=state(st);compare_fixed(snapshot,expected);dump(OUT/f'settings-{end}.json',snapshot)
   r,_,_,_=native(Path(str(stem)+'.cas.h5'),Path(str(stem)+'.dat.h5'));assert r['iteration']==end and r['relative_mass_imbalance']<=1e-4;dump(OUT/f'native-guard-{end}.json',r)
   st.results.report.fluxes.mass_flow(zones=['inlet','outlet','farfield','wing'],write_to_file=True,file_name=str(OUT/'force-reports'/f'mass-{end}.txt'))
   for label,direction in [('drag',[1,0,0]),('lift',[0,1,0])]:st.results.report.forces(option='forces',wall_zones=['wing'],direction_vector=direction,write_to_file=True,file_name=str(OUT/'force-reports'/f'{label}-{end}.txt'))
   capture=OUT/'runtime'/str(end);capture.mkdir();export_runtime_residuals(session,capture);mark(OUT/'status.json','saved_checkpoint',iteration=end,residuals=res)
  dump(OUT/'final-settings.json',state(st));compare_fixed(state(st),expected);st.file.write_case_data(file_name=str(OUT/'wing-final.cas.h5'))
  try:
   pic=st.results.graphics;pic.picture.use_window_resolution=False;pic.picture.x_resolution=1600;pic.picture.y_resolution=1000;pic.views.camera.projection='orthographic';pic.contour.create(name='wing-simplec-pressure');c=pic.contour['wing-simplec-pressure'];c.field='pressure';c.surfaces_list=['wing'];c.display();pic.views.camera.position=[4.,3.,7.];pic.views.camera.target=[.5,0,0];pic.views.camera.field=[6.8,4.25];pic.picture.save_picture(file_name=str(OUT/'wing-pressure-native.png'))
  except Exception:dump(OUT/'optional-picture-error.json',dict(error=traceback.format_exc()))
 except Exception:
  if loaded:
   try:st.file.write_case_data(file_name=str(OUT/'recovery.cas.h5'))
   except Exception:pass
  mark(OUT/'status.json','failed',error=traceback.format_exc());raise
 finally:
  if session:session.exit()
 verify();mark(OUT/'status.json','closed',iteration=1600,new_iterations=100)
if __name__=='__main__':main()
