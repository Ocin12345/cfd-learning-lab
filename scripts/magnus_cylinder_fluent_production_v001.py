"""Separately authorized native Magnus continuation5->500; exactly495steps, no init."""
from pathlib import Path
import argparse,datetime,hashlib,importlib.util,importlib.metadata,json,os,shutil,sys,tempfile,traceback
ROOT=Path('E:/CFD-Learning-Lab');LOCAL=ROOT/'local/magnus-cylinder-v001';RUN=ROOT/'native/runs/2026-11-02-magnus-cylinder-v001';ACTION='magnus_cylinder_fluent_production_v001'
COMMON_SHA='1a61e08a2f3624e6ccd36b3c5f244c67d3cbd89eee4865c2a67d030378905820';PILOT_SHA='6806e4ff269bfa95452d7af87d387394bbf70c42ae13f0d3504f3dfd66f644cc';AUDIT_SHA='41d5e426609689b2918e93fe7eda25ed08b27ceebd6643cc53fe0a1b02efa48d'
PRESERVED_SETTINGS=('general','energy','viscous','dynamic_mesh','material','cell_zones','walls','inlets','outlet','symmetries','references','methods','controls','rotation')
def require(ok,message):
 if not ok:raise ValueError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
def checked(path,digest):
 require(sha(path)==digest,'Changed source before import '+str(path));spec=importlib.util.spec_from_file_location(path.stem+'_production_checked',path);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
def validate_budget(auth):
 require(auth.get('approved') is True and auth.get('native_pilot_image_reviewed') is True,'Separate root production/image acceptance required')
 for key,value in dict(tandem_pair_and_local_record_complete=True,root='root',action=ACTION,phase='production-495-additional-steps',solver_processes=2,start_step=5,additional_steps=495,final_step=500,dt_s=.02,end_time_s=10.).items():require(auth.get(key)==value,'Changed production budget '+key)
 return True
def validate_pilot_acceptance(checks,manifest):
 require(isinstance(checks,dict) and isinstance(manifest,dict),'Future actual pilot checks/manifest required; source preparation is not acceptance')
 require(checks.get('all_primary_checks_passed') is True and all(v is True for v in checks['checks'].values()) and manifest.get('all_primary_checks_passed') is True,'Actual strict pilot acceptance absent')
 require(checks.get('phase')=='pilot-5-steps' and checks.get('native_velocity_png_available') is True,'Wrong pilot phase/native image');require(manifest.get('common_sha256')==COMMON_SHA and manifest.get('analyser_sha256')==AUDIT_SHA,'Pilot audit source differs');return True
def preserve_settings(loaded,accepted):
 require(all(k in loaded and k in accepted for k in PRESERVED_SETTINGS),'Incomplete actual physics/rotation snapshots')
 for key in PRESERVED_SETTINGS:require(loaded[key]==accepted[key],'Resumed actual physics/methods/rotation changed '+key)
 return True
def normalize_settings_path(parent,row):
 path=(Path(parent)/row['settings']).resolve();require(path.is_relative_to(Path(parent).resolve()),'Pilot settings outside authoritative original parent');row['settings']=str(path);return path
def resume_differences(before,after):
 import numpy as np
 require(set(before)==set(after)=={'u','v','p','rho','mu','all_native_face_mass_flux'},'Incomplete native resume fields');errors={}
 for key in before:
  a=np.asarray(before[key]);b=np.asarray(after[key]);require(a.shape==b.shape and np.isfinite(a).all() and np.isfinite(b).all(),'Invalid native resume field '+key);errors[key]=float(np.max(abs(a-b)))
 require(all(v<=1e-12 for v in errors.values()),'Native primary state changed during resume');return errors
def bound(c,path,digest):
 path=Path(path).resolve();c.require(path.is_relative_to(ROOT) and path.is_file() and c.sha(path)==digest,'Changed/missing accepted input '+str(path));return path
def verify_accepted_pilot(c,auth):
 checks_path=bound(c,auth['accepted_pilot_checks'],auth['accepted_pilot_checks_sha256']);manifest_path=bound(c,auth['accepted_pilot_manifest'],auth['accepted_pilot_manifest_sha256']);checks=c.read(checks_path);manifest=c.read(manifest_path);validate_pilot_acceptance(checks,manifest);c.require(checks_path.parent==manifest_path.parent,'Different pilot audit directories');inputs={}
 for row in manifest['inputs']:p=bound(c,row['path'],row['sha256']);inputs[p]=row['sha256']
 for name,digest in manifest['outputs'].items():c.require(Path(name).name==name,'Unsafe pilot output name');p=bound(c,manifest_path.parent/name,digest);inputs[p]=digest
 fields=bound(c,auth['accepted_pilot_fields'],auth['accepted_pilot_fields_sha256']);c.require(fields==(checks_path.parent/checks['audited_fields_npz']).resolve(),'Explicit pilot fields path differs from accepted audit');c.require(c.sha(fields)==checks['audited_fields_npz_sha256'],'Accepted audited fields differ');png=bound(c,checks['native_velocity_png'],auth['accepted_pilot_native_png_sha256'])
 parent=Path(auth['accepted_pilot_run']).resolve();c.require(parent.parent==RUN.resolve(),'Wrong accepted model parent');status=c.read(parent/'status.json');c.require(status['stage']=='pilot_completed' and status['actual_steps']==5 and status['actual_end_time_s']==.1,'Actual pilot incomplete');proof=c.read(parent/'code-provenance.json');bound(c,parent/'authorization.json',proof['root_authorization_sha256'])
 paths=checks['native_state_paths'];rows=c.read(parent/'states.json');c.require([r['index'] for r in paths]==list(range(6)) and [r['step'] for r in rows]==list(range(6)),'Actual pilot0..5 state coverage absent')
 for row,path in zip(rows,paths):
  case=bound(c,parent/row['case'],row['case_sha256']);data=bound(c,parent/row['data'],row['data_sha256']);c.require(case==Path(path['case']).resolve() and data==Path(path['data']).resolve() and abs(path['time_s']-row['step']*.02)<1e-12,'Pilot native path/clock mismatch');row['case']=str(case);row['data']=str(data);settings=normalize_settings_path(parent,row);bound(c,settings,row['settings_sha256']);c.settings_gate(c.read(settings))
 for p in [checks_path,manifest_path,fields,png]:inputs[p]=c.sha(p)
 return parent,rows,inputs
def preflight(authorization,out):
 auth=json.loads(Path(authorization).read_text(encoding='utf-8-sig'));validate_budget(auth);out=Path(out).resolve();require(auth.get('source_sha256')==sha(__file__) and auth.get('common_sha256')==COMMON_SHA,'Unreviewed production/common source');require(Path(auth['out']).resolve()==out and out.parent==RUN.resolve() and not out.exists(),'Fresh authorized native output required')
 c=checked(Path(__file__).with_name('magnus_cylinder_fluent_common_v001.py'),COMMON_SHA);pilot=checked(Path(__file__).with_name('magnus_cylinder_fluent_pilot_v001.py'),PILOT_SHA)
 inputs={Path(__file__).resolve():auth['source_sha256'],Path(c.__file__).resolve():COMMON_SHA,Path(pilot.__file__).resolve():PILOT_SHA,c.MESH:c.MESH_SHA,c.MESH.parent/'geometry.npz':c.GEOMETRY_SHA,c.MESH.parent/'setup.json':c.SETUP_SHA,LOCAL/'TARGETS.md':c.TARGETS_SHA,LOCAL/'EXECUTION-PLAN-v001.md':c.PLAN_SHA}
 for p,digest in inputs.items():bound(c,p,digest)
 parent,rows,accepted=verify_accepted_pilot(c,auth);inputs.update(accepted);c.require(pilot.EXECUTABLE.is_file(),'Installed Fluent absent');ram=pilot.available_memory();disk=shutil.disk_usage('E:/').free;c.require(ram>=4*1024**3 and disk>=10*1024**3,'Production requires4GiB RAM/10GiBE');idle=pilot.idle_check();versions={name:importlib.metadata.version(name) for name in ['ansys-fluent-core','numpy','h5py']};c.require(versions==c.read(parent/'code-provenance.json')['dependency_versions'],'Actual runtime dependencies changed from accepted pilot')
 return c,pilot,auth,parent,rows,dict(source_sha256=c.sha(__file__),root_authorization_sha256=c.sha(authorization),inputs=[dict(path=str(p),sha256=h) for p,h in inputs.items()],available_ram_bytes=ram,free_E_bytes=disk,processor_count=2,idle_check=idle,snapshot_before_solver_launch=True,python_version=sys.version,dependency_versions=versions,accepted_pilot_run=str(parent),accepted_pilot_final_case=rows[-1]['case'],accepted_pilot_final_data=rows[-1]['data'],no_reinitialization=True,additional_steps=495,final_step=500)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--prepare-check',action='store_true');a=p.parse_args();a.authorization=a.authorization.resolve();c,pilot,auth,parent,rows,proof=preflight(a.authorization,a.out)
 if a.prepare_check:print(json.dumps(dict(preflight_passed=True,solver_launched=False,**proof),indent=2));return 0
 import numpy as np
 out=a.out.resolve();out.mkdir();[(out/n).mkdir() for n in ['states','force-reports','settings','temp','examples']]
 for key in ['TEMP','TMP','TMPDIR']:os.environ[key]=str(out/'temp')
 tempfile.tempdir=str(out/'temp');os.environ['PYTHONDONTWRITEBYTECODE']='1';os.chdir(out)
 for src,name in [(Path(__file__),'executed-production.py'),(Path(c.__file__),'executed-common.py'),(a.authorization,'authorization.json'),(LOCAL/'TARGETS.md','TARGETS.md'),(LOCAL/'EXECUTION-PLAN-v001.md','EXECUTION-PLAN-v001.md'),(c.MESH.parent/'setup.json','mesh-setup.json'),(parent/'geometry-quality.json','geometry-quality.json'),(parent/'initialized-settings.json','initialized-settings.json')]:shutil.copyfile(src,out/name)
 c.write_json(out/'code-provenance.json',proof);c.write_json(out/'states.json',rows);session=None;events=[]
 def update(stage,**values):
  events.append(dict(stage=stage,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**values));c.write_json(out/'status.json',events[-1]);c.write_json(out/'events.json',events);print(json.dumps(events[-1]),flush=True)
 try:
  import ansys.fluent.core as pf
  pf.config.examples_path=str(out/'examples');update('launching_one_hidden_continuation_session',processor_count=2)
  session=pf.launch_fluent(fluent_path=str(pilot.EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True);update('native_session_connected',version=str(session.get_fluent_version()),licence_status='Actual Fluent continuation startup succeeded; native launcher/transcript retained')
  session.transcript.start(file_name=str(out/'solver-transcript.txt'));st=session.settings;st.file.read_case_data(file_name=proof['accepted_pilot_final_case']);loaded=pilot.snapshot(st);c.settings_gate(loaded);accepted=c.read(parent/'final-settings.json')
  preserve_settings(loaded,accepted)
  run=st.solution.run_calculation;c.require(abs(float(run.transient_controls.flow_time())-.1)<1e-12 and float(run.parameters.time_step_size())==.02 and int(run.parameters.max_iter_per_time_step())==100,'Resumed native clock/dt/inner cap differs');g=np.load(c.MESH.parent/'geometry.npz',allow_pickle=False);baseline=c.read_state(Path(proof['accepted_pilot_final_case']),Path(proof['accepted_pilot_final_data']),g);cache={k:baseline[k] for k in ['mapping','node_mapping','face_mapping','owner','neighbour','flux_to_source_sign','face_nodes','face_mid','face_area_vector','native_cell_areas','native_cell_centres','zones','mesh_fingerprint','volume_definition']};c.require(c.state_gates(baseline,5)['all_checks_passed'],'Accepted pilot fields changed')
  case=out/'resume-check.cas.h5';data=out/'resume-check.dat.h5';st.file.write_case_data(file_name=str(case));resumed=c.read_state(case,data,g,cache);c.require(c.state_gates(resumed,5)['all_checks_passed'],'Actual resumed fields fail');errors=resume_differences(dict(**baseline['native'],all_native_face_mass_flux=baseline['flux']),dict(**resumed['native'],all_native_face_mass_flux=resumed['flux']));c.require(c.body_forces(resumed)==c.body_forces(baseline),'Resumed native body pressure/shear forces changed')
  c.write_json(out/'resume-preservation-proof.json',dict(actual_resumed_step=5,time_s=.1,no_physical_step_executed=True,no_initialization_called=True,all_primary_fields_native_mesh_and_body_forces_retained=True,maximum_absolute_field_changes_SI=errors,case=case.name,data=data.name,case_sha256=c.sha(case),data_sha256=c.sha(data),diagnostic_duplicate_excluded_from501_state_history=True))
  del resumed,baseline
  previous_iteration=int(rows[-1]['native_terminal_iteration']);rr,last=c.terminal_residuals(Path(proof['accepted_pilot_final_data']),parent/'solver-transcript.txt');c.require(last==previous_iteration and rr==rows[-1]['terminal_scaled_residuals'],'Pilot final residual proof changed')
  st.solution.monitor.residual.options.n_save=100000;c.write_json(out/'residual-retention-control.json',dict(prior=loaded['residual']['options']['n_save'],current=100000,reason='Retain up to500steps*100inner samples; physics/methods/criteria unchanged'));c.write_json(out/'loaded-settings.json',loaded)
  for step in range(6,501):
   c.settings_gate(pilot.snapshot(st))
   c.require(shutil.disk_usage('E:/').free>=10*1024**3,'E free space floor crossed');update('advancing_one_actual_step',requested_step=step,additional_step=step-5,total_additional_budget=495);run.dual_time_iterate(time_step_count=1,max_iter_per_step=100)
   case=out/'states'/f'state-{step:04d}.cas.h5';data=case.with_name(f'state-{step:04d}.dat.h5');st.file.write_case_data(file_name=str(case));s=c.read_state(case,data,g,cache);gate=c.state_gates(s,step);actual_settings=pilot.snapshot(st);c.settings_gate(actual_settings);settings_path=out/'settings'/f'state-{step:04d}.json';c.write_json(settings_path,actual_settings);row=dict(settings=str(settings_path.relative_to(out)),settings_sha256=c.sha(settings_path),step=step,time_s=s['native_time']['flow-time'],case=str(case.relative_to(out)),data=str(data.relative_to(out)),case_sha256=c.sha(case),data_sha256=c.sha(data),native_gates=gate);rows.append(row);c.write_json(out/'states.json',rows)
   rr,index=c.terminal_residuals(data,out/'solver-transcript.txt');c.require(1<=index-previous_iteration<=100,'Observed native inner cap exceeded');previous_iteration=index;gate['checks']['own_terminal_residuals']=bool(all(v<=1e-6 for v in rr.values()));forces=c.body_forces(s);row.update(native_terminal_iteration=index,terminal_scaled_residuals=rr,body_forces=forces)
   if step==500:row['native_force_reports']=pilot.export_force_reports(st,out,step,forces,c)
   gate['all_checks_passed']=bool(all(gate['checks'].values()));c.write_json(out/'states.json',rows);update('actual_state_saved',step=step,time_s=s['native_time']['flow-time'],all_state_gates_passed=gate['all_checks_passed']);c.require(gate['all_checks_passed'],'Production native gate failed; saved fields retained')
  c.write_json(out/'final-settings.json',pilot.snapshot(st));st.solution.monitor.residual.write(file_name=str(out/'residuals.xy'))
  try:pilot.native_picture(st,out/'native-velocity.png')
  except Exception:(out/'graphics-export-error.txt').write_text(traceback.format_exc(),encoding='utf-8')
  verify_accepted_pilot(c,auth);update('production_completed',actual_additional_steps=495,actual_total_steps=500,actual_end_time_s=10.,native_states=501,analysis_pending=True);return 0
 except Exception as error:(out/'failure.txt').write_text(traceback.format_exc(),encoding='utf-8');update('failed_preserved',error=str(error),native_states=len(rows));raise
 finally:
  if session is not None:
   try:session.transcript.stop()
   finally:session.exit()
if __name__=='__main__':raise SystemExit(main())
