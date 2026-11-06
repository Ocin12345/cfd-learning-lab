"""Separate root-GO production continuation: accepted native pilot5 ->400, no init.
Unexecuted preparation. Exactly395additional physical steps; no automatic extension.
"""
from pathlib import Path
import argparse,datetime,json,math,os,shutil,tempfile,traceback
import numpy as np
import annular_couette_fluent_common_v003 as c
import annular_couette_fluent_pilot_v004 as pilot
from analyse_annular_couette_fluent_pilot_v005 import settings_gate
ROOT=Path('E:/CFD-Learning-Lab');LOCAL=ROOT/'local/annular-couette-v001';RUN=ROOT/'native/runs/2026-11-02-annular-couette-v001';MESH=RUN/'mesh/annular-couette.msh'
ACTION='annular_couette_fluent_production_v001'
def bound(path,digest):
 path=Path(path).resolve();c.require(path.is_relative_to(ROOT) and path.is_file() and c.sha(path)==digest,'Changed/missing production input '+str(path));return path
def verify_accepted_pilot(auth):
 checks_path=bound(auth['accepted_pilot_checks'],auth['accepted_pilot_checks_sha256']);manifest_path=bound(auth['accepted_pilot_manifest'],auth['accepted_pilot_manifest_sha256'])
 checks=c.read(checks_path);manifest=c.read(manifest_path)
 c.require(checks['all_primary_checks_passed'] is True and all(v is True for v in checks['checks'].values()) and manifest['all_primary_checks_passed'] is True,'Production requires the actually accepted strict pilot audit')
 c.require(checks_path.parent==manifest_path.parent,'Different accepted audit directories')
 for row in manifest['inputs']:bound(row['path'],row['sha256'])
 for name,digest in manifest['outputs'].items():c.require(Path(name).name==name,'Unsafe pilot audit output name');bound(manifest_path.parent/name,digest)
 paths=checks['native_state_paths'];c.require([r['index'] for r in paths]==list(range(6)) and all(abs(r['time_s']-k*.05)<1e-12 for k,r in enumerate(paths)),'Accepted pilot must have all actual native0..5 states')
 parent=Path(auth['accepted_pilot_run']).resolve();c.require(parent.parent==RUN.resolve(),'Accepted pilot is outside fixed model')
 rows=c.read(parent/'states.json');c.require([r['step'] for r in rows]==list(range(6)),'Pilot driver state history differs')
 for row,path in zip(rows,paths):
  case=parent/row['case'];data=parent/row['data'];c.require(case.resolve()==Path(path['case']).resolve() and data.resolve()==Path(path['data']).resolve(),'Accepted pilot native path mismatch')
  bound(case,row['case_sha256']);bound(data,row['data_sha256'])
  row['case']=str(case.resolve());row['data']=str(data.resolve());row['fields']=str((parent/'states'/row['fields']).resolve())
 return parent,rows,checks_path,manifest_path

def preflight(authorization,out):
 auth=c.read(authorization);out=Path(out).resolve()
 c.require(auth.get('approved') is True and auth.get('root')=='root' and auth.get('action')==ACTION and auth.get('phase')=='production-395-additional-steps','No separate root production GO')
 c.require(Path(auth['out']).resolve()==out and out.parent==RUN.resolve() and not out.exists(),'Fresh direct production output required')
 c.require(auth['solver_processes']==2 and auth['start_step']==5 and auth['additional_steps']==395 and auth['final_step']==400 and auth['dt_s']==.05 and auth['end_time_s']==20.,'Changed fixed continuation budget')
 parent,rows,checks_path,manifest_path=verify_accepted_pilot(auth)
 source={Path(__file__).resolve():auth['source_sha256'],Path(c.__file__).resolve():auth['common_sha256'],Path(pilot.__file__).resolve():'7de9f8781f861f5520de5758f2385b2d3052ef44b9486a1a2bb6f300e8efa0b2',Path(settings_gate.__code__.co_filename).resolve():'5976db1c39965c3aa32d0d62bea41dbe80652b55fb1268e2c1b531fbc8dace29',MESH:c.MESH_SHA,MESH.parent/'geometry.npz':c.GEOMETRY_SHA,MESH.parent/'setup.json':c.SETUP_SHA,LOCAL/'TARGETS.md':c.TARGETS_SHA,checks_path:auth['accepted_pilot_checks_sha256'],manifest_path:auth['accepted_pilot_manifest_sha256']}
 for path,digest in source.items():bound(path,digest)
 c.require(pilot.EXECUTABLE.is_file(),'Installed Fluent unavailable');ram=pilot.available_memory();disk=shutil.disk_usage('E:/').free;c.require(ram>=4*1024**3 and disk>=10*1024**3,'Unchanged production4GiB/10GiB admission failed')
 return auth,parent,rows,dict(source_sha256=c.sha(__file__),root_authorization_sha256=c.sha(authorization),inputs=[dict(path=str(p),sha256=d) for p,d in source.items()],available_ram_bytes=ram,free_E_bytes=disk,processor_count=2,snapshot_before_solver_launch=True,accepted_pilot_run=str(parent),accepted_pilot_final_case=rows[-1]['case'],accepted_pilot_final_data=rows[-1]['data'],accepted_pilot_checks_sha256=auth['accepted_pilot_checks_sha256'],accepted_pilot_manifest_sha256=auth['accepted_pilot_manifest_sha256'],no_reinitialization=True,additional_steps=395,final_step=400)
def snapshot_settings(st):
 return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),dynamic_mesh=st.setup.dynamic_mesh.get_state(),material=st.setup.materials.fluid['air'].get_state(),cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),walls={n:st.setup.boundary_conditions.wall[n].get_state() for n in ['inner','outer']},references=st.setup.reference_values.get_state(),methods=st.solution.methods.get_state(),run_calculation=st.solution.run_calculation.get_state(),residual=st.solution.monitor.residual.get_state())
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--prepare-check',action='store_true');a=p.parse_args();a.authorization=a.authorization.resolve()
 auth,parent,rows,proof=preflight(a.authorization,a.out)
 if a.prepare_check:print(json.dumps(dict(preflight_passed=True,solver_launched=False,**proof),indent=2));return 0
 out=a.out.resolve();out.mkdir();(out/'states').mkdir();(out/'temp').mkdir();(out/'examples').mkdir()
 for name in ['TEMP','TMP','TMPDIR']:os.environ[name]=str(out/'temp')
 tempfile.tempdir=str(out/'temp');os.chdir(out);os.environ['PYTHONDONTWRITEBYTECODE']='1'
 for source,target in [(Path(__file__),'executed-production.py'),(Path(c.__file__),'executed-common.py'),(a.authorization,'authorization.json'),(LOCAL/'TARGETS.md','TARGETS.md'),(MESH.parent/'setup.json','mesh-setup.json')]:shutil.copyfile(source,out/target)
 c.write_json(out/'code-provenance.json',proof);c.write_json(out/'declared-targets.json',c.TARGETS);c.write_json(out/'states.json',rows)
 session=None;events=[]
 def update(stage,**values):
  events.append(dict(stage=stage,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**values));c.write_json(out/'status.json',events[-1]);c.write_json(out/'events.json',events)
 try:
  import ansys.fluent.core as pf
  pf.config.examples_path=str(out/'examples');update('launching_one_hidden_continuation_session',processor_count=2)
  session=pf.launch_fluent(fluent_path=str(pilot.EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
  session.transcript.start(file_name=str(out/'solver-transcript.txt'));st=session.settings;st.file.read_case_data(file_name=proof['accepted_pilot_final_case'])
  g=np.load(MESH.parent/'geometry.npz');setup=c.read(MESH.parent/'setup.json');loaded=snapshot_settings(st);settings_gate(loaded,setup)
  accepted=c.read(parent/'final-settings.json')
  for key in ['general','energy','viscous','dynamic_mesh','material','cell_zones','walls','references','methods']:c.require(loaded[key]==accepted[key],'Native resumed physics/methods changed: '+key)
  run=st.solution.run_calculation;c.require(abs(float(run.transient_controls.flow_time())-.25)<1e-12 and float(run.parameters.time_step_size())==.05 and int(run.parameters.max_iter_per_time_step())==100,'Resumed actual0.25s/dt/innercap differs')
  baseline=c.read_state(Path(proof['accepted_pilot_final_case']),Path(proof['accepted_pilot_final_data']),g);c.require(c.state_gates(baseline,5)['all_checks_passed'],'Accepted pilot native fields no longer pass')
  # Save a distinct resume proof at0.25s before any further physical step.
  resume_case=out/'resume-check.cas.h5';resume_data=out/'resume-check.dat.h5';st.file.write_case_data(file_name=str(resume_case));resumed=c.read_state(resume_case,resume_data,g)
  c.require(c.state_gates(resumed,5)['all_checks_passed'],'Actual resumed native state no longer passes pilot gates')
  c.require(np.array_equal(resumed['mapping'],baseline['mapping']) and np.array_equal(resumed['face_nodes'],baseline['face_nodes']) and np.array_equal(resumed['native_cell_centres'],baseline['native_cell_centres']),'Native mapping changed during reload; refuse before additional physical step')
  errors={key:float(np.max(abs(resumed['common'][key]-baseline['common'][key]))) for key in ['u','v','p','rho','mu']};errors['all_native_face_mass_flux']=float(np.max(abs(resumed['flux']-baseline['flux'])))
  c.require(all(v<=1e-12 for v in errors.values()),'Actual saved pilot fields changed during reload')
  c.write_json(out/'resume-preservation-proof.json',dict(actual_resumed_step=5,time_s=.25,no_physical_step_executed=True,no_initialization_called=True,all_fields_and_native_mapping_retained=True,maximum_absolute_field_changes_SI=errors,case=resume_case.name,data=resume_data.name,case_sha256=c.sha(resume_case),data_sha256=c.sha(resume_data),diagnostic_duplicate_excluded_from401_state_history=True))
  previous_iteration=int(rows[-1]['native_terminal_iteration']);rr,last=c.terminal_residuals(Path(proof['accepted_pilot_final_data']),parent/'solver-transcript.txt');c.require(last==previous_iteration and rr==rows[-1]['terminal_scaled_residuals'],'Accepted pilot final residuals differ')
  # Recording retention alone: no numerical criterion/method/physics change.
  st.solution.monitor.residual.options.n_save=50000
  c.write_json(out/'residual-retention-control.json',dict(prior=loaded['residual']['options']['n_save'],current=50000,reason='Retain every native inner residual up to400steps*100inners; no solver criterion changed'))
  prop=st.setup.boundary_conditions.wall['inner'].momentum.rotation_speed.value;actual_units=prop.units();actual_rad=float(prop.as_quantity().to('radian s^-1').value);c.require(prop.is_active() and abs(actual_rad-1)<1e-12,'Resumed actual angular speed changed')
  physical=c.read(parent/'native-physical-readback.json');physical['rotation']['rotation_speed']=dict(active=True,native_value=float(prop()),native_units=actual_units,units_quantity=prop.get_attr('units-quantity'),rad_s=actual_rad);physical['loaded_from_accepted_pilot_without_configuration_setters']=True
  c.write_json(out/'native-physical-readback.json',physical);c.write_json(out/'loaded-settings.json',loaded);shutil.copyfile(parent/'geometry-quality.json',out/'geometry-quality.json')
  for step in range(6,401):
   c.require(shutil.disk_usage('E:/').free>=10*1024**3,'Production E floor crossed');update('advancing_one_actual_step',requested_step=step,additional_step=step-5,total_budget=395)
   run.dual_time_iterate(time_step_count=1,max_iter_per_step=100)
   stem=f'state-{step:04d}';case=out/'states'/(stem+'.cas.h5');data=out/'states'/(stem+'.dat.h5');st.file.write_case_data(file_name=str(case));s=c.read_state(case,data,g);gates=c.state_gates(s,step)
   c.require(np.array_equal(s['mapping'],baseline['mapping']) and np.array_equal(s['face_nodes'],baseline['face_nodes']) and np.array_equal(s['native_cell_centres'],baseline['native_cell_centres']),'Native fixed mesh changed')
   c.save_arrays(out/'states'/(stem+'-fields.npz'),s);residual,index=c.terminal_residuals(data,out/'solver-transcript.txt');c.require(1<=index-previous_iteration<=100,'Production native inner cap exceeded');previous_iteration=index
   gates['checks']['own_terminal_residuals']=bool(all(v<=1e-6 for v in residual.values()));gates['all_checks_passed']=bool(all(gates['checks'].values()))
   rows.append(dict(step=step,time_s=step*.05,case=str(case.relative_to(out)),data=str(data.relative_to(out)),case_sha256=c.sha(case),data_sha256=c.sha(data),fields=stem+'-fields.npz',native_gates=gates,native_terminal_iteration=index,terminal_scaled_residuals=residual));c.write_json(out/'states.json',rows);update('actual_state_saved',step=step,time_s=step*.05,all_state_gates_passed=gates['all_checks_passed']);c.require(gates['all_checks_passed'],'Production actual gate failed; saved diagnostic preserved')
  c.write_json(out/'final-settings.json',snapshot_settings(st));st.solution.monitor.residual.write(file_name=str(out/'residuals.xy'))
  graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1200;graphics.picture.y_resolution=1200;graphics.picture.invert_background=False;graphics.contour.create(name='annular-speed');contour=graphics.contour['annular-speed'];contour.field='velocity-magnitude';contour.surfaces_list=[];contour.display();graphics.views.camera.projection='orthographic';graphics.views.camera.position=[0.,0.,10.];graphics.views.camera.target=[0.,0.,0.];graphics.views.camera.field=[.225,.225];graphics.picture.save_picture(file_name=str(out/'native-velocity.png'))
  # Reprove accepted pilot originals at completion; never overwrite/reinitialize them.
  verify_accepted_pilot(auth);update('production_completed',actual_additional_steps=395,actual_total_steps=400,actual_end_time_s=20.,native_states=401,numerical_analysis_pending=True);return 0
 except Exception as error:
  (out/'failure.txt').write_text(traceback.format_exc(),encoding='utf-8');update('failed_preserved',error=str(error),native_states=len(rows));raise
 finally:
  if session is not None:
   try:session.transcript.stop()
   finally:session.exit()
if __name__=='__main__':raise SystemExit(main())
