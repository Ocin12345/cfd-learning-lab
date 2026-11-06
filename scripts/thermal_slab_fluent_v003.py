"""Authorized solid-only Fluent Euler slab,100 actual2s steps; no solver starts on import.

The first solved step must establish native solid/material/time/inward heat W and
Euler storage balance before any later step. Failures retain native recovery data.
"""
from pathlib import Path
import argparse,ctypes,json,os,shutil,subprocess,tempfile,time,traceback,re,math
from datetime import datetime
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')

def available_memory():
 class Status(ctypes.Structure):
  _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ('total_phys','available_phys','total_page','available_page','total_virtual','available_virtual','available_extended')]
 s=Status();s.length=ctypes.sizeof(s)
 if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(s)):raise ctypes.WinError()
 return s.available_phys

def require_idle():
 query="@(Get-Process fluent*,cortex* -ErrorAction SilentlyContinue) | ForEach-Object { $_.ProcessName + ':' + $_.Id }"
 p=subprocess.run(['powershell.exe','-NoProfile','-NonInteractive','-Command',query],capture_output=True,text=True,check=True)
 if p.stdout.strip():raise RuntimeError('Existing Fluent session: '+p.stdout.strip())
 p=subprocess.run(['wsl.exe','--list','--running','--quiet'],capture_output=True,check=True)
 running=p.stdout.decode('utf-16-le' if b'\x00' in p.stdout else 'utf-8',errors='replace').replace('\x00','').strip()
 if running:raise RuntimeError('WSL must be stopped: '+running)
 return dict(existing_fluent_processes=[],running_wsl_distributions=[])

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--mesh',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--authorization',type=Path,required=True)
 a=parser.parse_args();mesh=a.mesh.resolve();out=a.out.resolve();repo=Path(__file__).resolve().parents[1]
 if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo):raise ValueError('Fresh external E: native output required')
 if not mesh.is_file() or mesh.drive.upper()!='E:':raise ValueError('Approved E: mesh required')
 # Guard before importing solver or helper modules.
 import hashlib
 sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
 auth=json.loads(a.authorization.read_text(encoding='utf-8-sig'))
 required=dict(action='run_thermal_slab_fluent_v003',root='root',approved=True,allowed_mode='Run',source_sha256=sha(__file__),mesh_sha256=sha(mesh))
 if any(auth.get(k)!=v for k,v in required.items()) or Path(auth.get('out','')).resolve()!=out:raise ValueError('Root reviewed SHA-bound thermal Run authorization required')
 from thermal_slab_fluent_common_v003 import (write_json,source_files,verify_sources,geometry_check,physical_check,cell_zone_snapshot,topology,native_time,temperatures,state_metrics,step_metrics,TARGETS,WALLS)
 verify_sources(auth.get('helper_sha256',{}))
 import numpy as np
 setup_path=mesh.parent/'setup.json';geometry_path=mesh.parent/'geometry.npz';setup=json.loads(setup_path.read_text(encoding='utf-8-sig'))
 assert setup['mesh_sha256']==sha(mesh) and setup['geometry_npz_sha256']==sha(geometry_path) and setup['cells']==2000 and setup['cell_zone_name']=='slab'
 geometry=dict(np.load(geometry_path,allow_pickle=False));geometry_check(geometry)
 if not EXECUTABLE.is_file():raise FileNotFoundError(EXECUTABLE)
 memory=available_memory();disk=shutil.disk_usage(out.anchor).free
 if memory<4*1024**3 or disk<10*1024**3:raise RuntimeError('Require4GiBavailableRAM/10GiBEfree')
 idle=require_idle();out.mkdir(parents=True);(out/'snapshots').mkdir();(out/'heat-reports').mkdir();(out/'examples').mkdir();(out/'sources').mkdir()
 temporary=Path('E:/CFD-Learning-Lab/temp/thermal-slab-v001-fluent');temporary.mkdir(parents=True,exist_ok=True)
 for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
 os.environ['PYTHONDONTWRITEBYTECODE']='1';tempfile.tempdir=str(temporary);os.chdir(out)
 (out/'executed-script.py').write_bytes(Path(__file__).read_bytes())
 for name,p in source_files().items():shutil.copyfile(p,out/'sources'/name)
 shutil.copyfile(setup_path,out/'mesh-setup.json');write_json(out/'authorization.json',auth);write_json(out/'declared-targets.json',TARGETS)
 write_json(out/'code-provenance.json',dict(source_sha256=sha(__file__),helper_sha256=auth['helper_sha256'],mesh=str(mesh),mesh_sha256=sha(mesh),mesh_setup_sha256=sha(setup_path),geometry_npz_sha256=sha(geometry_path),snapshot_before_solver_launch=True,available_ram_before_launch_bytes=memory,free_disk_before_launch_bytes=disk,idle_check=idle,temporary_directory=str(temporary)))
 write_json(out/'requested-settings.json',dict(material='generic constant solid',density_kg_m3=1000,specific_heat_J_kg_K=1000,thermal_conductivity_W_m_K=10,solver='2D planar solid energy only',time_integration='first-order implicit Euler',time_step_s=2,time_steps=100,max_inner_iterations=20,energy_residual_target=1e-10,initial_temperature_K=300,hot_temperature_K=400,cold_temperature_K=300,insulation_heat_flux_W_m2=0,reference_depth_m=1,processor_count=2,ui_mode='hidden_gui',first_step_gate_required=True))
 session=None;loaded=False;saved=False;actual=0;events=[];records=[];errors=[]
 def update(stage,**values):
  row=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**values);events.append(row);write_json(out/'history.json',events);write_json(out/'status.json',row);print(json.dumps(row,default=str),flush=True)
 def optional(name,fn):
  try:fn()
  except Exception:errors.append(dict(export=name,error=traceback.format_exc()));write_json(out/'export-errors.json',errors)
 try:
  update('importing_pyfluent');import ansys.fluent.core as pf
  pf.config.examples_path=str(out/'examples')
  update('launching_single_session',processor_count=2,ui_mode='hidden_gui')
  session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
  update('solver_connected',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(out/'solver-transcript.txt'))
  st=session.settings;st.file.read_mesh(file_name=str(mesh));loaded=True
  st.setup.general.solver.type='pressure-based';st.setup.general.solver.two_dim_space='planar';st.setup.general.solver.time='transient'
  # Mesh import zone is temporary fluid only; solved cell collection must be solid.
  st.setup.cell_zone_conditions.set_zone_type(zone_list=['slab'],new_type='solid')
  converted=st.setup.cell_zone_conditions.get_state();write_json(out/'converted-cell-zones.json',converted)
  cell_zone_snapshot(converted)
  st.setup.models.energy.enabled=True
  st.setup.materials.solid.create(name='slab-solid');material=st.setup.materials.solid['slab-solid']
  for key,value in (('density',1000.),('specific_heat',1000.),('thermal_conductivity',10.)):
   prop=getattr(material,key);prop.option='value';prop.value=value
  st.setup.cell_zone_conditions.solid['slab'].general.material='slab-solid'
  st.setup.reference_values.depth=1.
  for name in WALLS:
   thermal=st.setup.boundary_conditions.wall[name].thermal
   thermal.thermal_condition='Temperature' if name in ('hot','cold') else 'Heat Flux'
   if name in ('hot','cold'):thermal.temperature.value=400. if name=='hot' else 300.
   else:thermal.heat_flux.value=0.
  methods=st.solution.methods;methods.transient_formulation='unsteady-1st-order'
  methods.spatial_discretization.gradient_scheme='least-square-cell-based'
  eq=st.solution.controls.get_state().get('equations',{});write_json(out/'discovered-equation-controls.json',eq)
  assert isinstance(eq,dict) and any(any(t in str(k).lower() for t in ('energy','temperature')) for k in eq),'No native thermal equation exposed'
  # Discover names from the actual dictionary; never query an inactive flow leaf.
  for key in eq:st.solution.controls.equations[key]=any(t in str(key).lower() for t in ('energy','temperature'))
  residual=st.solution.monitor.residual;residual.options.criterion_type='absolute';residual.options.n_save=2500
  write_json(out/'discovered-residual-state.json',residual.get_state());assert 'energy' in residual.get_state().get('equations',{})
  residual.equations['energy'].check_convergence=True;residual.equations['energy'].absolute_criteria=1e-10
  run=st.solution.run_calculation;run.transient_controls.type='Fixed';run.transient_controls.method='User-Specified'
  run.parameters.time_step_size=2.;run.parameters.max_iter_per_time_step=20;run.options.solid_time_step_size.enable_solid_time_step=False
  def state():
   return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),material=material.get_state(),**cell_zone_snapshot(st.setup.cell_zone_conditions.get_state()),walls={name:st.setup.boundary_conditions.wall[name].get_state() for name in WALLS},references=st.setup.reference_values.get_state(),methods=methods.get_state(),controls=st.solution.controls.get_state(),residual=residual.get_state(),run_calculation=run.get_state(),initialization=st.solution.initialization.get_state())
  native=state();write_json(out/'loaded-settings.json',native);physical_check(native)
  check=st.mesh.check();quality=st.mesh.quality()
  def quality_metrics(text):
   patterns=dict(minimum_orthogonal_quality=r'Minimum Orthogonal Quality\s*=',minimum_cell_volume_m3=r'minimum volume\s*\([^)]*\)\s*:',maximum_cell_volume_m3=r'maximum volume\s*\([^)]*\)\s*:',total_cell_volume_m3=r'total volume\s*\([^)]*\)\s*:',maximum_aspect_ratio=r'Maximum Aspect Ratio\s*=')
   result={}
   for key,pat in patterns.items():
    values=re.findall(pat+r'\s*([-+0-9.eE]+)',text,re.I)
    if values:result[key]=float(values[-1])
   assert all(map(math.isfinite,result.values())) and result.get('minimum_orthogonal_quality',-1)>=.15 and result.get('minimum_cell_volume_m3',0)>0 and abs(result.get('total_cell_volume_m3',0)-.002)<1e-10
   return result
  for attempt in range(10):
   try:q=quality_metrics(str(check)+'\n'+str(quality)+'\n'+(out/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace'));break
   except AssertionError:
    if attempt==9:raise
    time.sleep(.1)
  write_json(out/'geometry-quality.json',dict(**q,native_numeric_gate_pass=True,checked_before_energy_advance=True))
  init=st.solution.initialization;init.initialization_type='standard';init.defaults['temperature']=300.;init.standard_initialize()
  write_json(out/'initialized-settings.json',state());st.file.write_case(file_name=str(out/'slab.cas.h5'))
  initial=out/'snapshots/step-000000.dat.h5';st.file.write_data(file_name=str(initial))
  mapping,zones=topology(out/'slab.cas.h5',geometry);T0=temperatures(initial,mapping);time0=native_time(initial)
  assert time0['time-step']==0 and abs(time0['flow-time'])<1e-12 and abs(time0['physical-time-step']-2)<1e-12 and np.max(abs(T0-300))<1e-10
  m0=state_metrics(T0,geometry,0);records.append(dict(step=0,time_s=0.,native_dt_s=time0['physical-time-step'],native_time_index=0,native_iteration=0,initial_unsolved=True,**m0));write_json(out/'time-history.json',records)
  previous_energy=m0['energy_J'];previous_iteration=0
  for step in range(1,101):
   update('advancing_one_time_step',current_actual_step=actual,next_step=step,cap=100)
   run.dual_time_iterate(time_step_count=1,max_iter_per_step=20)
   data=out/f'snapshots/step-{step:06d}.dat.h5';st.file.write_data(file_name=str(data));actual=step
   report=out/f'heat-reports/step-{step:06d}.txt'
   st.results.report.fluxes.heat_transfer(zones=list(WALLS),write_to_file=True,file_name=str(report))
   T,row=step_metrics(data,report,mapping,zones,geometry,step,previous_energy,previous_iteration)
   actual_live_time=float(run.transient_controls.flow_time());assert abs(actual_live_time-row['time_s'])<1e-9
   records.append(row);write_json(out/'time-history.json',records);previous_energy=row['energy_J'];previous_iteration=row['native_iteration']
   if step==1:
    actual_state=state();write_json(out/'first-step-settings.json',actual_state);physical_check(actual_state)
    gate=dict(passed=True,actual_step=1,actual_time_s=row['time_s'],solid_zone='slab',fluid_zones_absent=True,constant_material_verified=True,native_temperature_time_verified=True,all240_native_boundary_heat_powers_verified=True,report_HDF_integrated_W_inward_sign_verified=True,global_Euler_energy_balance_verified=True,metrics=row)
    write_json(out/'first-step-gate.json',gate);update('first_step_gate_passed',**gate)
   update('time_step_saved_and_verified',**row)
  st.file.write_case_data(file_name=str(out/'slab-final.cas.h5'));saved=True
  final_state=state();write_json(out/'final-settings.json',final_state);physical_check(final_state)
  optional('native_energy_residual_xy',lambda:residual.write(file_name=str(out/'residuals.xy')))
  def picture():
   graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1600;graphics.picture.y_resolution=1000
   graphics.views.camera.projection='orthographic';graphics.contour.create(name='temperature');contour=graphics.contour['temperature'];contour.field='temperature';contour.surfaces_list=[];contour.display()
   graphics.views.camera.position=[.05,.01,1.];graphics.views.camera.target=[.05,.01,0.];graphics.views.camera.field=[.12,.075]
   graphics.picture.save_picture(file_name=str(out/'temperature.png'))
  optional('native_final_temperature_png',picture);optional('transcript_stop',session.transcript.stop)
  summary=dict(kind='transient-solid-conduction-slab',actual_steps=actual,actual_time_s=records[-1]['time_s'],actual_inner_iterations=records[-1]['native_iteration'],native_states=101,native_case='slab.cas.h5',native_final_case='slab-final.cas.h5',native_final_data='slab-final.dat.h5',snapshot_pattern='snapshots/step-{step:06d}.dat.h5',native_final_saved=True,first_step_gate_passed=True,final_metrics=records[-1],export_errors=errors,targets=TARGETS,physical_validation=False,numerical_checks='Pending independent read-only audit of all101 saved native states')
  write_json(out/'summary.json',summary);update('completed_with_export_errors' if errors else 'completed',**summary)
 except Exception:
  failure=dict(error=traceback.format_exc(),native_final_saved=saved,last_observed_actual_step=actual,export_errors=errors)
  write_json(out/'summary.json',dict(stage='failed',**failure));update('failed',**failure);raise
 finally:
  if session is not None:
   if loaded and not saved:
    try:session.settings.file.write_case_data(file_name=str(out/'recovery.cas.h5'))
    except Exception:write_json(out/'recovery-error.json',dict(error=traceback.format_exc()))
   try:session.exit()
   except Exception:write_json(out/'session-exit-error.json',dict(error=traceback.format_exc()))
if __name__=='__main__':main()
