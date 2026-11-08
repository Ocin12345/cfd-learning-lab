import os,json,sys,traceback,copy,time
import ansys.fluent.core as pf
from wing_repair_common_v001 import *
from wing_repair_prelaunch_v003 import wait_admission,verify
from wing_repair_native_v003 import state,expected_settings,same_native,state_keys,raw_residual
from wing_repair_helpers_v001 import original_target
from wing_fluent_v001 import mesh_quality,hdf_terminal
from sloshing_fixed_settings_v001 import compare_fixed
LABEL=sys.argv[1];assert LABEL=='low-relax25';OUT=N/('fluent-'+LABEL+'-v001');CRITERION=.001

def main():
 assert not (OUT/'solver-transcript.txt').exists();admission=wait_admission(lambda **v:mark(OUT/'status.json','waiting_admission',**v));dump(OUT/'launch-admission.json',admission);session=None;loaded=False
 try:
  pf.config.examples_path=str(R/'temp/fluent-examples');(R/'temp/fluent-examples').mkdir(exist_ok=True)
  session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=3,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True);dump(OUT/'solver-version.json',dict(version=str(session.get_fluent_version()),processor_count=2,connected_before_advance=True));session.transcript.start(file_name=str(OUT/'solver-transcript.txt'))
  st=session.settings;st.file.read_case_data(file_name=str(OUT/'source.cas.h5'));loaded=True
  reports=st.solution.monitor.report_files
  for name in reports.get_object_names():reports[name].file_name=str(OUT/(name+'-history.out'))
  assert reports['forces'].active() and reports['forces'].frequency()==1
  actual=state(st);dump(OUT/'loaded-source-settings.json',actual);base=json.loads((OUT/'expected-source-settings.json').read_text());compare_fixed(actual,{k:base[k] for k in state_keys()});assert not st.setup.models.energy.enabled()
  for n in EQ:assert original_target(st.solution.monitor.residual.equations[n])==1e-5
  compare_fixed(state(st),actual)
  st.file.write_case_data(file_name=str(OUT/'loaded-restart.cas.h5'));dump(OUT/'initial-field-check.json',same_native(OUT/'loaded-restart.cas.h5',OUT/'loaded-restart.dat.h5',OUT/'source.cas.h5',OUT/'source.dat.h5'));assert hdf_terminal(OUT/'loaded-restart.dat.h5')[0]==1500
  st.mesh.check();st.mesh.quality()
  for attempt in range(10):
   try:quality=mesh_quality((OUT/'solver-transcript.txt').read_text(errors='replace'));break
   except ValueError:
    if attempt==9:raise
    time.sleep(.1)
  dump(OUT/'geometry-quality.json',quality)
  mg=st.solution.controls.advanced.multi_grid;assert mg.mg_controls.get_object_names()==['pressure-coupled'] and mg.mg_controls['pressure-coupled'].termination_criteria()==.1
  mg.mg_controls['pressure-coupled'].termination_criteria=.001
  st.solution.controls.p_v_controls.flow_courant_number=10.
  st.solution.controls.p_v_controls.explicit_pressure_under_relaxation=.25
  st.solution.controls.p_v_controls.explicit_momentum_under_relaxation=.25
  mg.amg_controls.options.verbosity=1
  compare_fixed(state(st),expected_settings(CRITERION,10.,.25));dump(OUT/'candidate-settings.json',state(st));dump(OUT/'declared-test.json',dict(parent=1500,additional_iterations=100,terminal=1600,checkpoints=[1520,1540,1560,1580,1600],flow_courant=10.,coupled_pressure_and_momentum_explicit_relaxation=.25,pressure_coupled_AMG_termination=CRITERION,original_outer_residual_limit=1e-5,force_window=100,last20_residual_limit=1e-5))
  verify()
  for end in [1520,1540,1560,1580,1600]:
   mark(OUT/'status.json','iterating',target_iteration=end);st.solution.run_calculation.iterate(iter_count=20);stem=OUT/('iter-'+str(end));st.file.write_case_data(file_name=str(stem)+'.cas.h5');it,res=hdf_terminal(Path(str(stem)+'.dat.h5'));assert it==end and max(res.values())<.1
   st.results.report.fluxes.mass_flow(zones=['inlet','outlet','farfield','wing'],write_to_file=True,file_name=str(OUT/'force-reports'/f'mass-{end}.txt'))
   for label,direction in [('drag',[1,0,0]),('lift',[0,1,0])]:st.results.report.forces(option='forces',wall_zones=['wing'],direction_vector=direction,write_to_file=True,file_name=str(OUT/'force-reports'/f'{label}-{end}.txt'))
   mark(OUT/'status.json','saved_checkpoint',iteration=it,residuals=res)
  dump(OUT/'final-settings.json',state(st));st.file.write_case_data(file_name=str(OUT/'wing-final.cas.h5'))
  try:
   pic=st.results.graphics;pic.picture.use_window_resolution=False;pic.picture.x_resolution=1600;pic.picture.y_resolution=1000;pic.views.camera.projection='orthographic';pic.contour.create(name='wing-repair-pressure');c=pic.contour['wing-repair-pressure'];c.field='pressure';c.surfaces_list=['wing'];c.display();pic.views.camera.position=[4.,3.,7.];pic.views.camera.target=[.5,0,0];pic.views.camera.field=[6.8,4.25];pic.picture.save_picture(file_name=str(OUT/'wing-pressure-native.png'))
  except Exception:dump(OUT/'optional-picture-error.json',dict(error=traceback.format_exc()))
 except Exception:
  if loaded:
   try:st.file.write_case_data(file_name=str(OUT/'recovery.cas.h5'))
   except Exception:pass
  mark(OUT/'status.json','failed',error=traceback.format_exc());raise
 finally:
  if session:session.exit()
 verify();mark(OUT/'status.json','closed',iteration=1600)
if __name__=='__main__':main()
