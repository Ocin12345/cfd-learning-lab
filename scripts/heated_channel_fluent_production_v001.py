"""Accepted heated-channel pilot continuation: 599 additional thermal steps to total600.
No flow initialization, flow iteration, T patch, or first-step replay exists.
"""
from pathlib import Path
import argparse,ctypes,json,math,os,re,shutil,subprocess,tempfile,time,traceback,hashlib
from datetime import datetime
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
TARGETS_FILE=Path('E:/CFD-Learning-Lab/local/heated-channel-v001/TARGETS.md')
TARGETS_SHA='0a163061a9836c096d8b30b16c81410285e29bf3507f10d231619c43b64cab02'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')

def available_memory():
    class Status(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(name,ctypes.c_ulonglong) for name in
            ('total_phys','available_phys','total_page','available_page','total_virtual','available_virtual','available_extended')]
    status=Status();status.length=ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):raise ctypes.WinError()
    return status.available_phys

def require_idle():
    query="@(Get-Process fluent*,cortex* -ErrorAction SilentlyContinue) | ForEach-Object { $_.ProcessName + ':' + $_.Id }"
    found=subprocess.run(['powershell.exe','-NoProfile','-NonInteractive','-Command',query],capture_output=True,text=True,check=True).stdout.strip()
    if found:raise RuntimeError('Existing Fluent session must be saved and closed first: '+found)
    result=subprocess.run(['wsl.exe','--list','--running','--quiet'],capture_output=True,check=True)
    raw=result.stdout
    running=raw.decode('utf-16-le' if b'\x00' in raw else 'utf-8',errors='replace').replace('\x00','').strip()
    if running:raise RuntimeError('Require no concurrently running WSL distribution: '+running)
    return dict(existing_fluent_processes=[],running_wsl_distributions=[])

def quality_metrics(text):
    patterns=dict(minimum_orthogonal_quality=r'Minimum Orthogonal Quality\s*=',
        minimum_cell_volume_m3=r'minimum volume\s*\([^)]*\)\s*:',
        maximum_cell_volume_m3=r'maximum volume\s*\([^)]*\)\s*:',
        total_cell_volume_m3=r'total volume\s*\([^)]*\)\s*:',maximum_aspect_ratio=r'Maximum Aspect Ratio\s*=')
    result={}
    for name,pattern in patterns.items():
        matches=re.findall(pattern+r'\s*([-+0-9.eE]+)',text,re.I)
        if matches:result[name]=float(matches[-1])
    if not all(map(math.isfinite,result.values())):raise ValueError('Nonfinite native quality metric')
    if result.get('minimum_orthogonal_quality',-1)<.15 or result.get('minimum_cell_volume_m3',0)<=0:
        raise ValueError('Native positive-volume/quality .15 gate failed or missing')
    if abs(result.get('total_cell_volume_m3',0)-.1)>1e-8:raise ValueError('Native heated-channel total volume must equal0.1 m3 at reference depth1')
    return result

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--mesh',type=Path,required=True);ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--authorization',type=Path,required=True);a=ap.parse_args()
 mesh=a.mesh.resolve();parent=a.parent.resolve();out=a.out.resolve();repo=Path(__file__).resolve().parents[1]
 if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo) or out==parent:raise ValueError('Fresh external E: production output required')
 manifest_path=Path('E:/CFD-Learning-Lab/local/heated-channel-v001/production-parent-inputs-v001.json');manifest=json.loads(manifest_path.read_text(encoding='utf-8'));assert Path(manifest['parent']).resolve()==parent
 for name,digest in manifest['files'].items():assert sha(parent/name)==digest,'Parent changed: '+name
 audit_path=Path(manifest['selected_pilot_audit']);assert sha(audit_path)==manifest['selected_pilot_audit_sha256'];audit=json.loads(audit_path.read_text(encoding='utf-8'));assert audit['all_checks_passed'] and audit['errors']=={}
 auth=json.loads(a.authorization.read_text(encoding='utf-8-sig'));required=dict(action='continue_heated_channel_fluent_production_v001',root='root',approved=True,allowed_mode='Run',source_sha256=sha(__file__),mesh_sha256=sha(mesh),targets_sha256=TARGETS_SHA,parent_manifest_sha256=sha(manifest_path),parent_case_sha256=sha(parent/'thermal-step000001.cas.h5'),parent_data_sha256=sha(parent/'thermal-step000001.dat.h5'),selected_pilot_audit_sha256=sha(audit_path),total_thermal_steps=600,additional_thermal_steps=599,additional_flow_iterations=0,monitor_history_retention_n_save=40000)
 if any(auth.get(k)!=v for k,v in required.items()) or Path(auth.get('out','')).resolve()!=out or Path(auth.get('parent','')).resolve()!=parent:raise ValueError('Root SHA-bound accepted-pilot production authorization required')
 helpers={n:sha(Path(__file__).with_name(n)) for n in ('heated_channel_fluent_production_common_v001.py','analyse_nozzle_fluent_v005.py')}
 if auth.get('helper_sha256')!=helpers or sha(TARGETS_FILE)!=TARGETS_SHA:raise ValueError('Reviewed helper/frozen TARGETS changed')
 from heated_channel_fluent_production_common_v001 import TARGETS,WALLS,BOUNDARIES,source_files,geometry_check,snapshot,frozen_compare,flow_gates,thermal_metrics,native_time,thermal_storage_inventory
 import numpy as np
 geometry_path=mesh.parent/'geometry.npz';setup_path=mesh.parent/'setup.json';g=dict(np.load(geometry_path,allow_pickle=False));geometry_check(g);setup=json.loads(setup_path.read_text(encoding='utf-8-sig'));assert setup['mesh_sha256']==sha(mesh) and setup['geometry_npz_sha256']==sha(geometry_path)
 if not EXECUTABLE.is_file():raise FileNotFoundError(EXECUTABLE)
 memory=available_memory();disk=shutil.disk_usage(out.anchor).free
 if memory<4*1024**3 or disk<10*1024**3:raise RuntimeError('Require4GiBavailableRAM/10GiBEfree')
 idle=require_idle();out.mkdir(parents=True)
 for name in ('sources','examples','states','heat-reports','mass-reports','step-metrics'):(out/name).mkdir()
 temporary=Path('E:/CFD-Learning-Lab/temp/heated-channel-v001-fluent-production-v001');temporary.mkdir(parents=True,exist_ok=True)
 for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
 os.environ['PYTHONDONTWRITEBYTECODE']='1';tempfile.tempdir=str(temporary);os.chdir(out)
 (out/'executed-script.py').write_bytes(Path(__file__).read_bytes())
 for name,path in source_files().items():shutil.copyfile(path,out/'sources'/name)
 for src,dest in ((TARGETS_FILE,out/'TARGETS.md'),(setup_path,out/'mesh-setup.json'),(manifest_path,out/'parent-native-inputs.json'),(audit_path,out/'accepted-pilot-audit.json')):shutil.copyfile(src,dest)
 write_json(out/'authorization.json',auth);write_json(out/'declared-targets.json',TARGETS)
 inherited={}
 def inherit(src,dest):
  shutil.copyfile(src,dest);assert sha(src)==sha(dest);inherited[str(dest.relative_to(out)).replace('\\','/')]=dict(source=str(src),sha256=sha(src))
 for x in sorted(parent.iterdir()):
  if x.is_file() and (x.name.startswith('flow-') or x.name in ('geometry-quality.json','frozen-flow-source-hashes.json','native-flow-geometric-normal-proof.json')):inherit(x,out/x.name)
 inherit(parent/'thermal-initial.cas.h5',out/'thermal-initial.cas.h5');inherit(parent/'thermal-initial.dat.h5',out/'states/step-000000.dat.h5');inherit(parent/'thermal-step000001.dat.h5',out/'states/step-000001.dat.h5');inherit(parent/'thermal-step000001.cas.h5',out/'inherited-step000001.cas.h5');inherit(parent/'solver-transcript.txt',out/'pilot-solver-transcript.txt');inherit(parent/'parent-solver-transcript.txt',out/'flow-solver-transcript.txt')
 for index,label in ((0,'initial'),(1,'step000001')):
  inherit(parent/'reports'/(label+'-total-heat.txt'),out/'heat-reports'/f'step-{index:06d}.txt');inherit(parent/'reports'/(label+'-mass.txt'),out/'mass-reports'/f'step-{index:06d}.txt')
 write_json(out/'inherited-native-files.json',inherited)
 write_json(out/'code-provenance.json',dict(source_sha256=sha(__file__),helper_sha256=helpers,targets_sha256=TARGETS_SHA,mesh=str(mesh),mesh_sha256=sha(mesh),geometry_npz_sha256=sha(geometry_path),mesh_setup_sha256=sha(setup_path),snapshot_before_solver_launch=True,available_ram_before_launch_bytes=memory,free_disk_before_launch_bytes=disk,idle_check=idle,temporary_directory=str(temporary),continuation_parent=str(parent),parent_manifest_sha256=sha(manifest_path),selected_pilot_audit_sha256=sha(audit_path),parent_case_sha256=required['parent_case_sha256'],parent_data_sha256=required['parent_data_sha256'],additional_flow_iterations=0,total_thermal_steps=600,additional_thermal_steps=599,monitor_history_retention_n_save=40000))
 flow_reference=snapshot(parent/'flow-final.cas.h5',parent/'flow-final.dat.h5',g);baseline=snapshot(parent/'flow-iter1800.cas.h5',parent/'flow-iter1800.dat.h5',g);flow=flow_gates(flow_reference,baseline,g,parent/'flow-final.dat.h5');assert flow['all_checks_passed'];initial=snapshot(out/'thermal-initial.cas.h5',out/'states/step-000000.dat.h5',g,True);previous=snapshot(out/'inherited-step000001.cas.h5',out/'states/step-000001.dat.h5',g,True)
 assert abs(initial['common']['h'].mean()-1850)<1e-6 and np.max(abs(initial['common']['T']-300))<=1e-10 and frozen_compare(flow_reference,initial)['passed'] and frozen_compare(flow_reference,previous)['passed']
 # Recompute and retain the already solved seed step, never re-execute it.
 first=thermal_metrics(initial,previous,g,out/'states/step-000000.dat.h5',out/'states/step-000001.dat.h5',out/'heat-reports/step-000001.txt',out/'mass-reports/step-000001.txt',initial_total_report=out/'heat-reports/step-000000.txt',expected_step=1,transcript_path=out/'pilot-solver-transcript.txt');assert first['all_checks_passed'];write_json(out/'step-metrics/step-000001.json',first)
 parent_settings=json.loads((parent/'thermal-step000001-settings.json').read_text(encoding='utf-8'));session=None;loaded=False;final_saved=False;actual_step=1;events=[];errors=[];work_proof={};all_numeric=first['all_checks_passed'];max_frozen={k:0. for k in ('u','v','p','all_face_mass_flux')}
 def update(stage,**values):
  row=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**values);events.append(row);write_json(out/'history.json',events);write_json(out/'status.json',row);print(json.dumps(row,default=str),flush=True)
 def optional(name,fn):
  try:return fn()
  except Exception:errors.append(dict(export=name,error=traceback.format_exc()));write_json(out/'export-errors.json',errors);return None
 try:
  update('importing_pyfluent');import ansys.fluent.core as pf
  pf.config.examples_path=str(out/'examples');update('launching_single_hidden_production_session',processor_count=2)
  session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
  session.transcript.start(file_name=str(out/'solver-transcript.txt'));update('solver_connected',version=str(session.get_fluent_version()))
  st=session.settings;st.file.read_case_data(file_name=str(parent/'thermal-step000001.cas.h5'));loaded=True
  material=st.setup.materials.fluid['air'];inlet=st.setup.boundary_conditions.velocity_inlet['inlet'];outlet=st.setup.boundary_conditions.pressure_outlet['outlet'];methods=st.solution.methods;residual=st.solution.monitor.residual;init=st.solution.initialization;run=st.solution.run_calculation
  energy=st.setup.models.energy
  for name in ('pressure_work','kinetic_energy','viscous_dissipation'):
   active=bool(getattr(energy,name).is_active());actual=energy.get_state();work_proof[name]=dict(active=active,readback=actual.get(name),verified_off_or_inactive=actual.get(name) is False or not active)
  active=bool(st.setup.models.viscous.options.viscous_heating.is_active());value=st.setup.models.viscous.get_state().get('options',{}).get('viscous_heating');work_proof['viscous_heating']=dict(active=active,readback=value,verified_off_or_inactive=value is False or not active)
  active=bool(st.setup.models.radiation.model.is_active());value=st.setup.models.get_state().get('radiation',{}).get('model');work_proof['radiation_model']=dict(active=active,readback=value,verified_off_or_inactive=value=='none' or not active)
  work_proof['inlet_diffusion']=dict(active=bool(energy.inlet_diffusion.is_active()),readback=energy.get_state().get('inlet_diffusion'))
  def state():
   return dict(general=st.setup.general.get_state(),models=st.setup.models.get_state(),material=material.get_state(),cell_zone_conditions=st.setup.cell_zone_conditions.get_state(),inlet=inlet.get_state(),outlet=outlet.get_state(),walls={n:st.setup.boundary_conditions.wall[n].get_state() for n in WALLS},references=st.setup.reference_values.get_state(),methods=methods.get_state(),controls=st.solution.controls.get_state(),residual=residual.get_state(),run_calculation=run.get_state(),initialization=init.get_state(),work_control_proof=work_proof)
  loaded_settings=state();write_json(out/'production-loaded-settings.json',loaded_settings)
  keys=('models','material','cell_zone_conditions','inlet','outlet','walls','references','methods','controls','residual');equal={k:loaded_settings[k]==parent_settings[k] for k in keys};equal['solver']=loaded_settings['general']['solver']==parent_settings['general']['solver'];equal['operating_conditions']=loaded_settings['general']['operating_conditions']==parent_settings['general']['operating_conditions'];write_json(out/'parent-physical-settings-comparison.json',dict(checks=equal,all_checks_passed=all(equal.values())));assert all(equal.values())
  assert loaded_settings['controls']['equations']=={'flow':False,'temperature':True} and loaded_settings['run_calculation']['parameters']['time_step_size']==.005 and loaded_settings['run_calculation']['parameters']['max_iter_per_time_step']==60
  residual.options.n_save=40000
  retained=state();write_json(out/'production-start-settings.json',retained);before_residual=loaded_settings['residual'];after_residual=retained['residual'];assert after_residual['options']['n_save']==40000;expected=json.loads(json.dumps(before_residual));expected['options']['n_save']=40000;assert after_residual==expected;write_json(out/'monitor-retention-only-change.json',dict(before_n_save=before_residual['options']['n_save'],after_n_save=40000,all_other_residual_settings_identical=True,maximum_possible_total_history_rows=38000,reason='Retention only; no convergence criterion or solve changes'))
  st.file.write_case_data(file_name=str(out/'production-start.cas.h5'));live=snapshot(out/'production-start.cas.h5',out/'production-start.dat.h5',g,True);assert frozen_compare(flow_reference,live)['passed'] and np.array_equal(live['common']['T'],previous['common']['T']) and native_time(out/'production-start.dat.h5')['time-step']==1
  def reports(index):
   st.results.report.fluxes.heat_transfer(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'heat-reports'/f'step-{index:06d}.txt'))
   st.results.report.fluxes.mass_flow(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'mass-reports'/f'step-{index:06d}.txt'))
  for step in range(2,601):
   run.dual_time_iterate(time_step_count=1,max_iter_per_step=60);actual_step=step;reports(step);data=out/'states'/f'step-{step:06d}.dat.h5';st.file.write_data(file_name=str(data))
   current=snapshot(out/'production-start.cas.h5',data,g,True);frozen=frozen_compare(flow_reference,current)
   for key,value in frozen['maximum_absolute_changes_SI'].items():max_frozen[key]=max(max_frozen[key],value)
   write_json(out/'last-frozen-flow-check.json',dict(actual_step=step,**frozen));assert frozen['passed'],'Native frozen flow changed; preserve and stop'
   metrics=thermal_metrics(initial,current,g,out/'states/step-000000.dat.h5',data,out/'heat-reports'/f'step-{step:06d}.txt',out/'mass-reports'/f'step-{step:06d}.txt',initial_total_report=out/'heat-reports/step-000000.txt',previous=previous,previous_data=out/'states'/f'step-{step-1:06d}.dat.h5',expected_step=step,transcript_path=out/'solver-transcript.txt');write_json(out/'step-metrics'/f'step-{step:06d}.json',metrics);all_numeric=all_numeric and metrics['all_checks_passed'];previous=current
   if step%10==0 or not metrics['all_checks_passed']:update('thermal_state_saved',actual_thermal_step=step,actual_time_s=metrics['metrics']['native_time_step1']['flow-time'],step_checks_passed=metrics['all_checks_passed'],all_step_checks_so_far_passed=all_numeric,energy_final_scaled=metrics['metrics']['native_energy_final_scaled_residual'],native_Euler_relative_defect=metrics['metrics']['native_Euler_relative_defect'])
  st.file.write_case_data(file_name=str(out/'thermal-final.cas.h5'));final_saved=True;write_json(out/'thermal-final-settings.json',state());write_json(out/'thermal-final-native-storage-inventory.json',thermal_storage_inventory(out/'thermal-final.dat.h5',previous['zones']))
  final=snapshot(out/'thermal-final.cas.h5',out/'thermal-final.dat.h5',g,True);assert np.array_equal(final['common']['T'],previous['common']['T']) and frozen_compare(flow_reference,final)['passed']
  def picture():
   graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1600;graphics.picture.y_resolution=1000;graphics.views.camera.projection='orthographic';graphics.contour.create(name='temperature');contour=graphics.contour['temperature'];contour.field='temperature';contour.surfaces_list=[];contour.display();graphics.views.camera.position=[.5,.05,1.];graphics.views.camera.target=[.5,.05,0.];graphics.views.camera.field=[1.15,.71875];graphics.picture.save_picture(file_name=str(out/'temperature.png'))
  optional('native_temperature_png',picture);optional('native_energy_residual_xy',lambda:residual.write(file_name=str(out/'residuals.xy')));optional('transcript_stop',session.transcript.stop)
  unchanged=all(sha(parent/n)==h for n,h in manifest['files'].items());write_json(out/'parent-native-inputs-unchanged.json',dict(all_unchanged=unchanged));assert unchanged
  summary=dict(actual_flow_iterations=2000,additional_flow_iterations=0,actual_thermal_steps=actual_step,additional_thermal_steps_executed=actual_step-1,native_saved_states=601,actual_time_s=native_time(out/'states/step-000600.dat.h5')['flow-time'],all_numeric_step_checks_passed=all_numeric,maximum_frozen_field_changes_SI=max_frozen,original_parent_files_unchanged=unchanged,export_errors=errors,physical_validation=False)
  write_json(out/'summary.json',summary);update('completed_with_export_errors' if errors else 'completed',**summary)
 except Exception:
  failure=dict(error=traceback.format_exc(),actual_flow_iterations=2000,additional_flow_iterations=0,actual_thermal_steps=actual_step,additional_thermal_steps_executed=actual_step-1,export_errors=errors);write_json(out/'summary.json',dict(stage='failed',**failure));update('failed',**failure);raise
 finally:
  if session is not None:
   if loaded and not final_saved:
    try:session.settings.file.write_case_data(file_name=str(out/'recovery.cas.h5'))
    except Exception:write_json(out/'recovery-error.json',dict(error=traceback.format_exc()))
   try:session.exit()
   except Exception:write_json(out/'session-exit-error.json',dict(error=traceback.format_exc()))
if __name__=='__main__':main()
