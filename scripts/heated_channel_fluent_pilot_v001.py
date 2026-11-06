"""Authorized heated-channel pilot: fixed2000 steadyflow then EXACTLY ONE energy step.

No solver starts on import. Production600steps is NOT implemented or authorized.
Native heat/mass reports are evaluated before saving each thermal state.
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
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--mesh',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);ap.add_argument('--authorization',type=Path,required=True);a=ap.parse_args()
 mesh=a.mesh.resolve();out=a.out.resolve();repo=Path(__file__).resolve().parents[1]
 if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo):raise ValueError('Fresh external E: nativepilot output required')
 auth=json.loads(a.authorization.read_text(encoding='utf-8-sig'));required=dict(action='run_heated_channel_fluent_pilot_v001',root='root',approved=True,allowed_mode='Run',source_sha256=sha(__file__),mesh_sha256=sha(mesh),targets_sha256=TARGETS_SHA,pilot_thermal_steps=1)
 if any(auth.get(k)!=v for k,v in required.items()) or Path(auth.get('out','')).resolve()!=out:raise ValueError('Root SHA-bound one-step pilot authorization required')
 helpers={n:sha(Path(__file__).with_name(n)) for n in ('heated_channel_fluent_common_v001.py','analyse_nozzle_fluent_v005.py')}
 if auth.get('helper_sha256')!=helpers or sha(TARGETS_FILE)!=TARGETS_SHA:raise ValueError('Reviewed helper/frozen TARGETS changed')
 from heated_channel_fluent_common_v001 import (TARGETS,WALLS,BOUNDARIES,EQUATIONS,source_files,geometry_check,snapshot,frozen_compare,flow_gates,terminal_residuals,thermal_metrics,native_time)
 import numpy as np
 geometry_path=mesh.parent/'geometry.npz';setup_path=mesh.parent/'setup.json';g=dict(np.load(geometry_path,allow_pickle=False));geometry_check(g);setup=json.loads(setup_path.read_text(encoding='utf-8-sig'))
 assert setup['mesh_sha256']==sha(mesh) and setup['geometry_npz_sha256']==sha(geometry_path) and setup['cells']==2304
 if not EXECUTABLE.is_file():raise FileNotFoundError(EXECUTABLE)
 memory=available_memory();disk=shutil.disk_usage(out.anchor).free
 if memory<4*1024**3 or disk<10*1024**3:raise RuntimeError('Require4GiBavailableRAM/10GiBEfree')
 idle=require_idle();out.mkdir(parents=True)
 for name in ('sources','examples','reports'): (out/name).mkdir()
 temporary=Path('E:/CFD-Learning-Lab/temp/heated-channel-v001-fluent-pilot');temporary.mkdir(parents=True,exist_ok=True)
 for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
 os.environ['PYTHONDONTWRITEBYTECODE']='1';tempfile.tempdir=str(temporary);os.chdir(out)
 (out/'executed-script.py').write_bytes(Path(__file__).read_bytes())
 for name,path in source_files().items():shutil.copyfile(path,out/'sources'/name)
 shutil.copyfile(TARGETS_FILE,out/'TARGETS.md');shutil.copyfile(setup_path,out/'mesh-setup.json');write_json(out/'authorization.json',auth);write_json(out/'declared-targets.json',TARGETS)
 write_json(out/'code-provenance.json',dict(source_sha256=sha(__file__),helper_sha256=helpers,targets_sha256=TARGETS_SHA,mesh=str(mesh),mesh_sha256=sha(mesh),geometry_npz_sha256=sha(geometry_path),mesh_setup_sha256=sha(setup_path),snapshot_before_solver_launch=True,available_ram_before_launch_bytes=memory,free_disk_before_launch_bytes=disk,idle_check=idle,temporary_directory=str(temporary)))
 session=None;loaded=False;flow_saved=False;pilot_saved=False;actual_flow=0;actual_thermal=0;events=[];errors=[];work_proof={};flow_reference=None
 def update(stage,**values):
  row=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**values);events.append(row);write_json(out/'history.json',events);write_json(out/'status.json',row);print(json.dumps(row,default=str),flush=True)
 def optional(name,fn):
  try:return fn()
  except Exception:errors.append(dict(export=name,error=traceback.format_exc()));write_json(out/'export-errors.json',errors);return None
 try:
  update('importing_pyfluent');import ansys.fluent.core as pf
  pf.config.examples_path=str(out/'examples');update('launching_single_hidden_session',processor_count=2)
  session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
  session.transcript.start(file_name=str(out/'solver-transcript.txt'));update('solver_connected',version=str(session.get_fluent_version()))
  st=session.settings;st.file.read_mesh(file_name=str(mesh));loaded=True
  st.setup.general.solver.type='pressure-based';st.setup.general.solver.two_dim_space='planar';st.setup.general.solver.time='steady';st.setup.models.viscous.model='laminar';st.setup.models.energy.enabled=False
  st.setup.general.operating_conditions.operating_pressure=0.;st.setup.general.operating_conditions.gravity.enable=False
  material=st.setup.materials.fluid['air'];material.density.option='value';material.density.value=1.;material.viscosity.option='value';material.viscosity.value=.01;st.setup.cell_zone_conditions.fluid['fluid'].general.material='air'
  inlet=st.setup.boundary_conditions.velocity_inlet['inlet'];inlet.momentum.velocity_specification_method='Components';inlet.momentum.velocity_components=[dict(value=1.),dict(value=0.)]
  outlet=st.setup.boundary_conditions.pressure_outlet['outlet'];outlet.momentum.gauge_pressure.value=0.
  for name in WALLS:
   wall=st.setup.boundary_conditions.wall[name].momentum;wall.wall_motion='Stationary Wall';wall.shear_condition='No Slip'
  for key,value in dict(density=1.,viscosity=.01,velocity=1.,length=.2,area=.1,depth=1.).items():setattr(st.setup.reference_values,key,value)
  methods=st.solution.methods;methods.p_v_coupling.flow_scheme='SIMPLE';methods.spatial_discretization.gradient_scheme='least-square-cell-based';methods.spatial_discretization.discretization_scheme={'mom':'second-order-upwind','pressure':'second-order'}
  residual=st.solution.monitor.residual;residual.options.criterion_type='absolute';residual.options.n_save=3000
  for name in EQUATIONS:residual.equations[name].absolute_criteria=1e-6;residual.equations[name].check_convergence=False
  def state():
   return dict(general=st.setup.general.get_state(),models=st.setup.models.get_state(),material=material.get_state(),cell_zone_conditions=st.setup.cell_zone_conditions.get_state(),inlet=inlet.get_state(),outlet=outlet.get_state(),walls={n:st.setup.boundary_conditions.wall[n].get_state() for n in WALLS},references=st.setup.reference_values.get_state(),methods=methods.get_state(),controls=st.solution.controls.get_state(),residual=residual.get_state(),run_calculation=st.solution.run_calculation.get_state(),initialization=st.solution.initialization.get_state(),work_control_proof=work_proof)
  write_json(out/'flow-loaded-settings.json',state());check=st.mesh.check();quality=st.mesh.quality()
  for attempt in range(10):
   try:q=quality_metrics(str(check)+'\n'+str(quality)+'\n'+(out/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace'));break
   except ValueError:
    if attempt==9:raise
    time.sleep(.1)
  write_json(out/'geometry-quality.json',dict(**q,checked_before_flow=True))
  init=st.solution.initialization;init.initialization_type='standard'
  for key,value in (('x-velocity',0.),('y-velocity',0.),('pressure',0.)):init.defaults[key]=value
  init.standard_initialize();st.file.write_case_data(file_name=str(out/'flow-initial.cas.h5'));write_json(out/'flow-initialized-settings.json',state())
  for endpoint in (500,1000,1500,1800,2000):
   update('iterating_fixed_flow',current_actual_iteration=actual_flow,requested_additional_iterations=endpoint-actual_flow,total_cap=2000)
   st.solution.run_calculation.iterate(iter_count=endpoint-actual_flow);stem=f'flow-iter{endpoint:04d}';st.file.write_case_data(file_name=str(out/(stem+'.cas.h5')))
   r,index=terminal_residuals(out/(stem+'.dat.h5'),EQUATIONS);assert all(v[-1]==endpoint for v in index.values());actual_flow=endpoint;update('flow_checkpoint_saved',actual_iteration=actual_flow,final_scaled_residuals=r)
  st.file.write_case_data(file_name=str(out/'flow-final.cas.h5'));flow_saved=True;write_json(out/'flow-final-settings.json',state())
  st.results.report.fluxes.mass_flow(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'reports/flow-mass.txt'))
  flow_reference=snapshot(out/'flow-final.cas.h5',out/'flow-final.dat.h5',g);baseline=snapshot(out/'flow-iter1800.cas.h5',out/'flow-iter1800.dat.h5',g);flow=flow_gates(flow_reference,baseline,g,out/'flow-final.dat.h5');write_json(out/'flow-gates.json',flow)
  if not flow['all_checks_passed']:raise RuntimeError('Actual steady-flow numerical gates missed; see flow-gates.json. Thermal step not executed.')
  write_json(out/'native-flow-geometric-normal-proof.json',flow_reference['face_orientation_proof'])
  frozen_file_hashes={n:sha(out/n) for n in ('flow-final.cas.h5','flow-final.dat.h5','flow-iter1800.cas.h5','flow-iter1800.dat.h5')}
  write_json(out/'frozen-flow-source-hashes.json',frozen_file_hashes)
  # Enable energy without another initialization; only native temperature is patched.
  st.setup.models.energy.enabled=True
  radiation=st.setup.models.radiation.model;active=bool(radiation.is_active())
  if active:radiation.set_state('none')
  radiation_parent=st.setup.models.get_state().get('radiation',{});work_proof['radiation_model']=dict(active=active,readback=radiation_parent.get('model'),verified_off_or_inactive=(radiation_parent.get('model')=='none' or not active))
  for key,value in (('density',1.),('viscosity',.01),('specific_heat',1000.),('thermal_conductivity',10.)):
   prop=getattr(material,key);prop.option='value';prop.value=value
  energy=st.setup.models.energy
  for name in ('pressure_work','kinetic_energy','viscous_dissipation'):
   leaf=getattr(energy,name);active=bool(leaf.is_active())
   if active:leaf.set_state(False)
   actual=energy.get_state();work_proof[name]=dict(active=active,readback=actual.get(name),verified_off_or_inactive=(actual.get(name) is False or not active))
  heating=st.setup.models.viscous.options.viscous_heating;active=bool(heating.is_active())
  if active:heating.set_state(False)
  work_proof['viscous_heating']=dict(active=active,readback=st.setup.models.viscous.get_state().get('options',{}).get('viscous_heating'),verified_off_or_inactive=not active or st.setup.models.viscous.get_state().get('options',{}).get('viscous_heating') is False)
  diffusion=energy.inlet_diffusion;work_proof['inlet_diffusion']=dict(active=bool(diffusion.is_active()),before_parent_state=energy.get_state())
  write_json(out/'energy-control-discovery.json',work_proof)
  if not diffusion.is_active():raise RuntimeError('Fluent inlet diffusion control inactive; cannot claim matched inlet conduction. Pilot stops for root review.')
  diffusion.set_state(True);assert energy.get_state().get('inlet_diffusion') is True
  inlet.thermal.temperature.value=300.;outlet.thermal.backflow_total_temperature.value=300.
  for name in WALLS:thermal=st.setup.boundary_conditions.wall[name].thermal;thermal.thermal_condition='Temperature';thermal.temperature.value=350.
  st.setup.general.solver.time='transient';methods.transient_formulation='unsteady-1st-order';methods.spatial_discretization.discretization_scheme={'temperature':'first-order-upwind'}
  eq=st.solution.controls.get_state().get('equations',{});write_json(out/'discovered-thermal-equation-keys.json',eq);assert any('temperature' in str(n).lower() or 'energy' in str(n).lower() for n in eq)
  for name in eq:st.solution.controls.equations[name]=('temperature' in str(name).lower() or 'energy' in str(name).lower())
  active_eq=st.solution.controls.get_state().get('equations',{});assert all(('temperature' in str(n).lower() or 'energy' in str(n).lower()) for n,v in active_eq.items() if v)
  for name in EQUATIONS:
   if name in residual.get_state().get('equations',{}):residual.equations[name].check_convergence=False
  residual.equations['energy'].check_convergence=True;residual.equations['energy'].absolute_criteria=1e-8
  run=st.solution.run_calculation;run.transient_controls.type='Fixed';run.transient_controls.method='User-Specified';run.parameters.time_step_size=.005;run.parameters.max_iter_per_time_step=60
  def thermal_reports(label):
   st.results.report.fluxes.heat_transfer(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'reports'/(label+'-total-heat.txt')))
   sensible=st.results.report.fluxes.heat_transfer_sensible;write_json(out/'reports'/(label+'-sensible-command.json'),dict(active=bool(sensible.is_active())))
   if sensible.is_active():sensible(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'reports'/(label+'-sensible-heat.txt')))
   st.results.report.fluxes.mass_flow(zones=list(BOUNDARIES),write_to_file=True,file_name=str(out/'reports'/(label+'-mass.txt')))
  write_json(out/'thermal-before-patch-settings.json',state());thermal_reports('before-patch');st.file.write_case_data(file_name=str(out/'thermal-before-patch.cas.h5'))
  before_patch=snapshot(out/'thermal-before-patch.cas.h5',out/'thermal-before-patch.dat.h5',g,thermal=True);frozen_before=frozen_compare(flow_reference,before_patch);write_json(out/'frozen-before-patch.json',frozen_before)
  if not frozen_before['passed']:raise RuntimeError('Flow changed during thermal enable/freeze before patch')
  patch=init.patch.calculate_patch;variables=patch.variable.allowed_values();zones=patch.cell_zones.allowed_values();write_json(out/'temperature-patch-discovery.json',dict(variable_allowed_values=variables,cell_zone_allowed_values=zones))
  assert variables and 'temperature' in variables and zones and 'fluid' in zones
  patch(variable='temperature',value=300.,cell_zones=['fluid'])
  write_json(out/'thermal-initial-settings.json',state());thermal_reports('initial');st.file.write_case_data(file_name=str(out/'thermal-initial.cas.h5'))
  initial=snapshot(out/'thermal-initial.cas.h5',out/'thermal-initial.dat.h5',g,thermal=True);frozen_initial=frozen_compare(flow_reference,initial);write_json(out/'frozen-initial.json',frozen_initial)
  if not frozen_initial['passed'] or np.max(abs(initial['common']['T']-300))>1e-10:raise RuntimeError('Temperature-only patch failed frozen-flow/uniform300K proof')
  update('advancing_EXACTLY_ONE_thermal_step',physical_step_cap=1,dt_s=.005,max_inner_iterations=60)
  run.dual_time_iterate(time_step_count=1,max_iter_per_step=60);actual_thermal=1
  thermal_reports('step000001');write_json(out/'thermal-step000001-settings.json',state());st.file.write_case_data(file_name=str(out/'thermal-step000001.cas.h5'));pilot_saved=True
  def picture():
   graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1600;graphics.picture.y_resolution=1000
   graphics.views.camera.projection='orthographic';graphics.contour.create(name='temperature');contour=graphics.contour['temperature'];contour.field='temperature';contour.surfaces_list=[];contour.display()
   graphics.views.camera.position=[.5,.05,1.];graphics.views.camera.target=[.5,.05,0.];graphics.views.camera.field=[1.15,.71875]
   graphics.picture.save_picture(file_name=str(out/'temperature.png'))
  optional('native_temperature_png',picture)
  final=snapshot(out/'thermal-step000001.cas.h5',out/'thermal-step000001.dat.h5',g,thermal=True);frozen_final=frozen_compare(flow_reference,final);write_json(out/'frozen-step000001.json',frozen_final)
  thermal=thermal_metrics(initial,final,g,out/'thermal-initial.dat.h5',out/'thermal-step000001.dat.h5',out/'reports/step000001-total-heat.txt',out/'reports/step000001-mass.txt',initial_total_report=out/'reports/initial-total-heat.txt');write_json(out/'pilot-thermal-gates.json',thermal)
  r,index=terminal_residuals(out/'thermal-step000001.dat.h5',('energy',));write_json(out/'pilot-native-field-inventory.json',dict(initial_cell_fields=initial['cell_field_names'],final_cell_fields=final['cell_field_names'],initial_face_fields=initial['face_field_names'],final_face_fields=final['face_field_names'],initial_geometric_normal_proof=initial['face_orientation_proof'],final_geometric_normal_proof=final['face_orientation_proof']))
  optional('native_energy_residual_xy',lambda:residual.write(file_name=str(out/'residuals.xy')));optional('transcript_stop',session.transcript.stop)
  original_flow_files_unchanged=all(sha(out/n)==digest for n,digest in frozen_file_hashes.items());write_json(out/'frozen-flow-original-files-check.json',dict(all_unchanged=original_flow_files_unchanged,sha256=frozen_file_hashes))
  summary=dict(kind='heated-channel-one-step-pilot',actual_flow_iterations=actual_flow,actual_thermal_steps=actual_thermal,actual_time_s=native_time(out/'thermal-step000001.dat.h5')['flow-time'],flow_native_saved=flow_saved,pilot_native_saved=pilot_saved,flow_gates_passed=flow['all_checks_passed'],frozen_flow_before_patch=frozen_before,frozen_flow_initial=frozen_initial,frozen_flow_step1=frozen_final,pilot_thermal_checks=thermal,all_pilot_checks_passed=flow['all_checks_passed'] and frozen_before['passed'] and frozen_initial['passed'] and frozen_final['passed'] and thermal['all_checks_passed'],export_errors=errors,production_authorized=False,original_flow_files_unchanged=original_flow_files_unchanged,physical_validation=False)
  write_json(out/'summary.json',summary);update('pilot_completed',**summary)
 except Exception:
  failure=dict(error=traceback.format_exc(),actual_flow_iterations=actual_flow,actual_thermal_steps=actual_thermal,flow_native_saved=flow_saved,pilot_native_saved=pilot_saved,export_errors=errors,production_authorized=False)
  write_json(out/'summary.json',dict(stage='failed',**failure));update('failed',**failure);raise
 finally:
  if session is not None:
   if loaded and not pilot_saved:
    try:session.settings.file.write_case_data(file_name=str(out/'recovery.cas.h5'))
    except Exception:write_json(out/'recovery-error.json',dict(error=traceback.format_exc()))
   try:session.exit()
   except Exception:write_json(out/'session-exit-error.json',dict(error=traceback.format_exc()))
if __name__=='__main__':main()
