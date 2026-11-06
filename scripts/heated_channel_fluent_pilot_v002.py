"""Root-authorized continuation of saved converged flow: temperature patch and ONE thermal step only.
No flow initialization, flow iteration, or production loop exists.
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
 if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo) or out==parent:raise ValueError('Fresh external E: continuation output required')
 parent_manifest_path=Path('E:/CFD-Learning-Lab/local/heated-channel-v001/parent-native-inputs-v002.json')
 parent_manifest=json.loads(parent_manifest_path.read_text(encoding='utf-8'));assert Path(parent_manifest['parent']).resolve()==parent
 for name,digest in parent_manifest['files'].items():assert sha(parent/name)==digest,'Parent changed: '+name
 auth=json.loads(a.authorization.read_text(encoding='utf-8-sig'));required=dict(action='continue_heated_channel_fluent_pilot_v002',root='root',approved=True,allowed_mode='Run',source_sha256=sha(__file__),mesh_sha256=sha(mesh),targets_sha256=TARGETS_SHA,pilot_thermal_steps=1,additional_flow_iterations=0,parent_manifest_sha256=sha(parent_manifest_path),parent_case_sha256=sha(parent/'thermal-before-patch.cas.h5'),parent_data_sha256=sha(parent/'thermal-before-patch.dat.h5'))
 if any(auth.get(k)!=v for k,v in required.items()) or Path(auth.get('out','')).resolve()!=out or Path(auth.get('parent','')).resolve()!=parent:raise ValueError('Root SHA-bound saved-parent one-step continuation authorization required')
 helpers={n:sha(Path(__file__).with_name(n)) for n in ('heated_channel_fluent_common_v002.py','analyse_nozzle_fluent_v005.py')}
 if auth.get('helper_sha256')!=helpers or sha(TARGETS_FILE)!=TARGETS_SHA:raise ValueError('Reviewed helper/frozen TARGETS changed')
 from heated_channel_fluent_common_v002 import TARGETS,WALLS,BOUNDARIES,source_files,geometry_check,snapshot,frozen_compare,flow_gates,thermal_metrics,native_time,thermal_storage_inventory
 import numpy as np
 geometry_path=mesh.parent/'geometry.npz';setup_path=mesh.parent/'setup.json';g=dict(np.load(geometry_path,allow_pickle=False));geometry_check(g);setup=json.loads(setup_path.read_text(encoding='utf-8-sig'))
 assert setup['mesh_sha256']==sha(mesh) and setup['geometry_npz_sha256']==sha(geometry_path) and setup['cells']==2304
 assert json.loads((parent/'status.json').read_text(encoding='utf-8'))['actual_thermal_steps']==0
 if not EXECUTABLE.is_file():raise FileNotFoundError(EXECUTABLE)
 memory=available_memory();disk=shutil.disk_usage(out.anchor).free
 if memory<4*1024**3 or disk<10*1024**3:raise RuntimeError('Require4GiBavailableRAM/10GiBEfree')
 idle=require_idle();out.mkdir(parents=True)
 for name in ('sources','examples','reports'):(out/name).mkdir()
 temporary=Path('E:/CFD-Learning-Lab/temp/heated-channel-v001-fluent-pilot-v002');temporary.mkdir(parents=True,exist_ok=True)
 for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
 os.environ['PYTHONDONTWRITEBYTECODE']='1';tempfile.tempdir=str(temporary);os.chdir(out)
 (out/'executed-script.py').write_bytes(Path(__file__).read_bytes())
 for name,path in source_files().items():shutil.copyfile(path,out/'sources'/name)
 for src,dest in ((TARGETS_FILE,out/'TARGETS.md'),(setup_path,out/'mesh-setup.json'),(parent_manifest_path,out/'parent-native-inputs.json')):shutil.copyfile(src,dest)
 write_json(out/'authorization.json',auth);write_json(out/'declared-targets.json',TARGETS)
 # Exact byte copies retain original flow evidence and its provenance; no flow is recomputed.
 inherited={}
 for p in sorted(parent.iterdir()):
  if p.is_file() and (p.name.startswith('flow-') or p.name.startswith('thermal-before-patch') or p.name in ('geometry-quality.json','frozen-before-patch.json','frozen-flow-source-hashes.json','native-flow-geometric-normal-proof.json')):
   shutil.copyfile(p,out/p.name);inherited[p.name]=dict(source=str(p),sha256=sha(p));assert sha(out/p.name)==sha(p)
 for name in ('before-patch-total-heat.txt','before-patch-mass.txt','before-patch-sensible-command.json'):shutil.copyfile(parent/'reports'/name,out/'reports'/name);inherited['reports/'+name]=dict(source=str(parent/'reports'/name),sha256=sha(parent/'reports'/name))
 shutil.copyfile(parent/'solver-transcript.txt',out/'parent-solver-transcript.txt');inherited['parent-solver-transcript.txt']=dict(source=str(parent/'solver-transcript.txt'),sha256=sha(parent/'solver-transcript.txt'))
 write_json(out/'inherited-native-files.json',inherited)
 write_json(out/'code-provenance.json',dict(source_sha256=sha(__file__),helper_sha256=helpers,targets_sha256=TARGETS_SHA,mesh=str(mesh),mesh_sha256=sha(mesh),geometry_npz_sha256=sha(geometry_path),mesh_setup_sha256=sha(setup_path),snapshot_before_solver_launch=True,available_ram_before_launch_bytes=memory,free_disk_before_launch_bytes=disk,idle_check=idle,temporary_directory=str(temporary),continuation_parent=str(parent),parent_manifest_sha256=sha(parent_manifest_path),additional_flow_iterations=0,parent_case_sha256=required['parent_case_sha256'],parent_data_sha256=required['parent_data_sha256']))
 flow_reference=snapshot(parent/'flow-final.cas.h5',parent/'flow-final.dat.h5',g);baseline=snapshot(parent/'flow-iter1800.cas.h5',parent/'flow-iter1800.dat.h5',g);flow=flow_gates(flow_reference,baseline,g,parent/'flow-final.dat.h5');assert flow==json.loads((parent/'flow-gates.json').read_text(encoding='utf-8')) and flow['all_checks_passed']
 parent_settings=json.loads((parent/'thermal-before-patch-settings.json').read_text(encoding='utf-8'))
 session=None;loaded=False;pilot_saved=False;actual_thermal=0;events=[];errors=[];work_proof={}
 def update(stage,**values):
  row=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**values);events.append(row);write_json(out/'history.json',events);write_json(out/'status.json',row);print(json.dumps(row,default=str),flush=True)
 def optional(name,fn):
  try:return fn()
  except Exception:errors.append(dict(export=name,error=traceback.format_exc()));write_json(out/'export-errors.json',errors);return None
 try:
  update('importing_pyfluent');import ansys.fluent.core as pf
  pf.config.examples_path=str(out/'examples');update('launching_single_hidden_continuation_session',processor_count=2)
  session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
  session.transcript.start(file_name=str(out/'solver-transcript.txt'));update('solver_connected',version=str(session.get_fluent_version()))
  st=session.settings;st.file.read_case_data(file_name=str(parent/'thermal-before-patch.cas.h5'));loaded=True
  material=st.setup.materials.fluid['air'];inlet=st.setup.boundary_conditions.velocity_inlet['inlet'];outlet=st.setup.boundary_conditions.pressure_outlet['outlet'];methods=st.solution.methods;residual=st.solution.monitor.residual;init=st.solution.initialization;run=st.solution.run_calculation
  energy=st.setup.models.energy
  for name in ('pressure_work','kinetic_energy','viscous_dissipation'):
   active=bool(getattr(energy,name).is_active());actual=energy.get_state();work_proof[name]=dict(active=active,readback=actual.get(name),verified_off_or_inactive=actual.get(name) is False or not active)
  active=bool(st.setup.models.viscous.options.viscous_heating.is_active());value=st.setup.models.viscous.get_state().get('options',{}).get('viscous_heating');work_proof['viscous_heating']=dict(active=active,readback=value,verified_off_or_inactive=value is False or not active)
  active=bool(st.setup.models.radiation.model.is_active());value=st.setup.models.get_state().get('radiation',{}).get('model');work_proof['radiation_model']=dict(active=active,readback=value,verified_off_or_inactive=value=='none' or not active)
  work_proof['inlet_diffusion']=dict(active=bool(energy.inlet_diffusion.is_active()),readback=energy.get_state().get('inlet_diffusion'))
  def state():
   return dict(general=st.setup.general.get_state(),models=st.setup.models.get_state(),material=material.get_state(),cell_zone_conditions=st.setup.cell_zone_conditions.get_state(),inlet=inlet.get_state(),outlet=outlet.get_state(),walls={n:st.setup.boundary_conditions.wall[n].get_state() for n in WALLS},references=st.setup.reference_values.get_state(),methods=methods.get_state(),controls=st.solution.controls.get_state(),residual=residual.get_state(),run_calculation=run.get_state(),initialization=init.get_state(),work_control_proof=work_proof)
  loaded_settings=state();write_json(out/'continuation-loaded-settings.json',loaded_settings)
  equal={key:loaded_settings[key]==parent_settings[key] for key in ('models','material','cell_zone_conditions','inlet','outlet','walls','references','methods','controls','residual')}
  equal['solver']=loaded_settings['general']['solver']==parent_settings['general']['solver'];equal['operating_conditions']=loaded_settings['general']['operating_conditions']==parent_settings['general']['operating_conditions']
  write_json(out/'parent-physical-settings-comparison.json',dict(checks=equal,all_checks_passed=all(equal.values())));assert all(equal.values()),'Loaded parent physical settings changed'
  assert loaded_settings['controls']['equations']=={'flow':False,'temperature':True} and all(v['verified_off_or_inactive'] for k,v in work_proof.items() if k!='inlet_diffusion') and work_proof['inlet_diffusion']['active'] and work_proof['inlet_diffusion']['readback'] is True
  st.file.write_case_data(file_name=str(out/'continuation-loaded.cas.h5'));loaded_native=snapshot(out/'continuation-loaded.cas.h5',out/'continuation-loaded.dat.h5',g,thermal=True)
  frozen_loaded=frozen_compare(flow_reference,loaded_native);write_json(out/'frozen-loaded-parent.json',frozen_loaded);assert frozen_loaded['passed'];assert native_time(out/'continuation-loaded.dat.h5')=={'flow-time':0.,'time-step':0.,'physical-time-step':.005}
  write_json(out/'prepatch-uninitialized-storage-inventory.json',thermal_storage_inventory(out/'continuation-loaded.dat.h5',loaded_native['zones']))
  def thermal_reports(label):
   st.results.report.fluxes.heat_transfer(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'reports'/(label+'-total-heat.txt')))
   sensible=st.results.report.fluxes.heat_transfer_sensible;active=bool(sensible.is_active());write_json(out/'reports'/(label+'-sensible-command.json'),dict(active=active))
   if active:sensible(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'reports'/(label+'-sensible-heat.txt')))
   st.results.report.fluxes.mass_flow(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'reports'/(label+'-mass.txt')))
  patch=init.patch.calculate_patch;variables=patch.variable.allowed_values();zones=patch.cell_zones.allowed_values()
  discovery=dict(variable_allowed_values=variables,cell_zone_allowed_values=zones,empty_variable_enumeration_is_unavailable_metadata=True,source_backed_temperature_command=True,required_native_post_patch_T_K=300.,flow_reinitialized=False)
  write_json(out/'temperature-patch-discovery.json',discovery)
  assert (not variables or 'temperature' in variables) and 'fluid' in zones and set(loaded_settings['cell_zone_conditions']['fluid'])=={'fluid'} and loaded_settings['models']['energy']['enabled'] is True
  patch(variable='temperature',value=300.,cell_zones=['fluid'])
  write_json(out/'thermal-initial-settings.json',state());thermal_reports('initial');st.file.write_case_data(file_name=str(out/'thermal-initial.cas.h5'))
  initial=snapshot(out/'thermal-initial.cas.h5',out/'thermal-initial.dat.h5',g,thermal=True);frozen_initial=frozen_compare(flow_reference,initial);write_json(out/'frozen-initial.json',frozen_initial)
  inventory=thermal_storage_inventory(out/'thermal-initial.dat.h5',initial['zones']);write_json(out/'thermal-initial-native-storage-inventory.json',inventory)
  h0=initial['common']['h'];refT=float(loaded_settings['controls']['advanced']['enthalpy_reference_temperature']);proof=dict(native_T300_maximum_absolute_error_K=float(abs(initial['common']['T']-300.).max()),native_H_minimum=float(h0.min()),native_H_maximum=float(h0.max()),native_H_mean=float(h0.mean()),native_H_uniform_finite=bool(np.isfinite(h0).all() and np.ptp(h0)<=1e-6),native_enthalpy_reference_temperature_K=refT,constant_cp_reference_diagnostic_J_kg=1000.*(300.-refT),constant_cp_reference_diagnostic_error_J_kg=float(abs(h0-1000.*(300.-refT)).max()),native_U_V_p_and_all_face_flux_frozen=frozen_initial['passed'],initial_heat_unsolved=True)
  write_json(out/'temperature-only-patch-native-proof.json',proof)
  assert frozen_initial['passed'] and proof['native_T300_maximum_absolute_error_K']<=1e-10 and proof['native_H_uniform_finite'],'Actual native temperature-only patch proof failed'
  update('advancing_EXACTLY_ONE_thermal_step',physical_step_cap=1,dt_s=.005,max_inner_iterations=60,additional_flow_iterations=0)
  run.dual_time_iterate(time_step_count=1,max_iter_per_step=60);actual_thermal=1
  thermal_reports('step000001');write_json(out/'thermal-step000001-settings.json',state());st.file.write_case_data(file_name=str(out/'thermal-step000001.cas.h5'));pilot_saved=True
  def picture():
   graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1600;graphics.picture.y_resolution=1000
   graphics.views.camera.projection='orthographic';graphics.contour.create(name='temperature');contour=graphics.contour['temperature'];contour.field='temperature';contour.surfaces_list=[];contour.display()
   graphics.views.camera.position=[.5,.05,1.];graphics.views.camera.target=[.5,.05,0.];graphics.views.camera.field=[1.15,.71875];graphics.picture.save_picture(file_name=str(out/'temperature.png'))
  optional('native_temperature_png',picture)
  final=snapshot(out/'thermal-step000001.cas.h5',out/'thermal-step000001.dat.h5',g,thermal=True);frozen_final=frozen_compare(flow_reference,final);write_json(out/'frozen-step000001.json',frozen_final)
  write_json(out/'thermal-step000001-native-storage-inventory.json',thermal_storage_inventory(out/'thermal-step000001.dat.h5',final['zones']))
  thermal=thermal_metrics(initial,final,g,out/'thermal-initial.dat.h5',out/'thermal-step000001.dat.h5',out/'reports/step000001-total-heat.txt',out/'reports/step000001-mass.txt',initial_total_report=out/'reports/initial-total-heat.txt');write_json(out/'pilot-thermal-gates.json',thermal)
  optional('native_energy_residual_xy',lambda:residual.write(file_name=str(out/'residuals.xy')));optional('transcript_stop',session.transcript.stop)
  unchanged=all(sha(parent/n)==digest for n,digest in parent_manifest['files'].items());write_json(out/'parent-native-inputs-unchanged.json',dict(all_unchanged=unchanged,parent_manifest_sha256=sha(parent_manifest_path)))
  hashes=json.loads((out/'frozen-flow-source-hashes.json').read_text(encoding='utf-8'));unchanged_copies=all(sha(out/n)==digest for n,digest in hashes.items());write_json(out/'frozen-flow-original-files-check.json',dict(all_unchanged=unchanged_copies,sha256=hashes,original_sources_in_parent=str(parent)))
  summary=dict(kind='heated-channel-one-step-continuation-pilot',actual_flow_iterations=2000,additional_flow_iterations=0,actual_thermal_steps=actual_thermal,actual_time_s=native_time(out/'thermal-step000001.dat.h5')['flow-time'],flow_native_saved=True,flow_native_reused_from_parent=str(parent),pilot_native_saved=pilot_saved,flow_gates_passed=flow['all_checks_passed'],frozen_flow_before_patch=frozen_loaded,frozen_flow_initial=frozen_initial,frozen_flow_step1=frozen_final,pilot_thermal_checks=thermal,all_pilot_checks_passed=flow['all_checks_passed'] and frozen_loaded['passed'] and frozen_initial['passed'] and frozen_final['passed'] and thermal['all_checks_passed'] and unchanged and unchanged_copies,export_errors=errors,production_authorized=False,original_flow_files_unchanged=unchanged and unchanged_copies,physical_validation=False)
  write_json(out/'summary.json',summary);update('pilot_completed',**summary)
 except Exception:
  failure=dict(error=traceback.format_exc(),actual_flow_iterations=2000,additional_flow_iterations=0,actual_thermal_steps=actual_thermal,flow_native_saved=True,pilot_native_saved=pilot_saved,export_errors=errors,production_authorized=False)
  write_json(out/'summary.json',dict(stage='failed',**failure));update('failed',**failure);raise
 finally:
  if session is not None:
   if loaded and not pilot_saved:
    try:session.settings.file.write_case_data(file_name=str(out/'recovery.cas.h5'))
    except Exception:write_json(out/'recovery-error.json',dict(error=traceback.format_exc()))
   try:session.exit()
   except Exception:write_json(out/'session-exit-error.json',dict(error=traceback.format_exc()))
if __name__=='__main__':main()
