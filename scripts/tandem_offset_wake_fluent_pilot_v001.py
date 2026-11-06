"""Root-authorized hidden Fluent tandem PILOT: exactly5steps; no production entry."""
from pathlib import Path
import argparse,ctypes,datetime,hashlib,importlib.util,json,math,os,re,shutil,subprocess,sys,tempfile,time,traceback
ROOT=Path('E:/CFD-Learning-Lab');LOCAL=ROOT/'local/tandem-offset-wake-v001';RUN=ROOT/'native/runs/2026-11-02-tandem-offset-wake-v001'
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe');ACTION='tandem_offset_wake_fluent_pilot_v001'
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
def available_memory():
 class MS(ctypes.Structure):_fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ['total','available','total_page','available_page','total_virtual','available_virtual','available_extended']]
 m=MS();m.length=ctypes.sizeof(m)
 if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):raise ctypes.WinError()
 return m.available
def idle_check():
 query="@(Get-Process fluent*,cortex* -ErrorAction SilentlyContinue) | ForEach-Object { $_.ProcessName + ':' + $_.Id }";found=subprocess.run(['powershell.exe','-NoProfile','-NonInteractive','-Command',query],capture_output=True,text=True,check=True).stdout.strip()
 if found:raise ValueError('Existing Fluent processes '+found)
 result=subprocess.run(['wsl.exe','--list','--running','--quiet'],capture_output=True,check=True);running=result.stdout.decode('utf-16-le' if b'\x00' in result.stdout else 'utf-8',errors='replace').replace('\x00','').strip()
 if running:raise ValueError('Require root-stopped CFD WSL before Fluent launch: '+running)
 return dict(existing_Fluent=[],running_WSL=[])
def preflight(authorization,out):
 auth=json.loads(Path(authorization).read_text(encoding='utf-8-sig'));out=Path(out).resolve()
 if any(auth.get(k)!=v for k,v in dict(approved=True,root='root',action=ACTION,phase='pilot-5-steps',solver_processes=2,pilot_steps=5,dt_s=.05).items()):raise ValueError('Wrong root pilot authorization or fixed budget')
 if Path(auth.get('out','')).resolve()!=out or out.parent!=RUN.resolve() or out.exists():raise ValueError('Fresh authorized model output required')
 if auth.get('source_sha256')!=sha(__file__):raise ValueError('Runner source changed')
 helper=Path(__file__).with_name('tandem_offset_wake_fluent_common_v001.py')
 if auth.get('common_sha256')!=sha(helper):raise ValueError('Common source changed before import')
 spec=importlib.util.spec_from_file_location('tandem_pilot_common_checked',helper);c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
 sources={Path(__file__).resolve():auth['source_sha256'],helper.resolve():auth['common_sha256'],c.MESH:c.MESH_SHA,c.MESH.parent/'geometry.npz':c.GEOMETRY_SHA,c.MESH.parent/'setup.json':c.SETUP_SHA,LOCAL/'TARGETS.md':c.TARGETS_SHA,LOCAL/'EXECUTION-PLAN-v001.md':c.PLAN_SHA}
 for key,digest in [('mesh_sha256',c.MESH_SHA),('geometry_sha256',c.GEOMETRY_SHA),('setup_sha256',c.SETUP_SHA),('targets_sha256',c.TARGETS_SHA)]:c.require(auth.get(key)==digest,'Authorization source binding differs '+key)
 for path,digest in sources.items():c.require(path.is_file() and c.sha(path)==digest,'Changed source '+str(path))
 import numpy as np
 with np.load(c.MESH.parent/'geometry.npz',allow_pickle=False) as g:c.geometry_check(g)
 c.require(EXECUTABLE.is_file(),'Installed Fluent missing');ram=available_memory();disk=shutil.disk_usage('E:/').free;c.require(ram>=4*1024**3 and disk>=10*1024**3,'Require4GiB availableRAM and10GiBE');idle=idle_check()
 return c,auth,dict(source_sha256=c.sha(__file__),root_authorization_sha256=c.sha(authorization),inputs=[dict(path=str(p),sha256=h) for p,h in sources.items()],available_ram_bytes=ram,free_E_bytes=disk,processor_count=2,idle_check=idle,snapshot_before_solver_launch=True)
def quality(text,c):
 result={}
 for key,pattern in [('minimum_orthogonal_quality',r'Minimum Orthogonal Quality\s*='),('minimum_native_cell_volume',r'minimum volume\s*\([^)]*\)\s*:'),('total_native_cell_volume',r'total volume\s*\([^)]*\)\s*:')]:
  values=re.findall(pattern+r'\s*([-+0-9.eE]+)',text,re.I);c.require(values,'Missing native quality '+key);result[key]=float(values[-1])
 c.require(all(map(math.isfinite,result.values())) and result['minimum_orthogonal_quality']>=.15 and result['minimum_native_cell_volume']>0 and abs(result['total_native_cell_volume']-1.48)<1e-8,'Native quality/positive volume/total area gate failed');return dict(**result,checked_before_flow=True)
def snapshot(st):
 return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),dynamic_mesh=st.setup.dynamic_mesh.get_state(),material=st.setup.materials.fluid['air'].get_state(),cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),walls={n:st.setup.boundary_conditions.wall[n].get_state() for n in ['body1','body2']},inlets={'inlet':st.setup.boundary_conditions.velocity_inlet['inlet'].get_state()},outlet=st.setup.boundary_conditions.pressure_outlet['outlet'].get_state(),symmetries=st.setup.boundary_conditions.get_state().get('symmetry',{}),references=st.setup.reference_values.get_state(),methods=st.solution.methods.get_state(),controls=st.solution.controls.get_state(),run_calculation=st.solution.run_calculation.get_state(),residual=st.solution.monitor.residual.get_state(),initialization=st.solution.initialization.get_state())
def export_force_reports(st,out,step,forces,c):
 reports={}
 for body in c.BODIES:
  for label,direction in [('x',[1,0,0]),('y',[0,1,0])]:
   path=out/'force-reports'/f'step{step:04d}-{body}-{label}.txt';st.results.report.forces(option='forces',wall_zones=[body],direction_vector=direction,write_to_file=True,file_name=str(path));reports[body+'_'+label]=dict(path=str(path.relative_to(out)),sha256=c.sha(path),report=c.force_report_agreement(forces,body,direction,path))
 return reports
def native_picture(st,path):
 graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1600;graphics.picture.y_resolution=750;graphics.picture.invert_background=False;graphics.contour.create(name='tandem-speed');contour=graphics.contour['tandem-speed'];contour.field='velocity-magnitude';contour.surfaces_list=[];contour.display();graphics.views.camera.projection='orthographic';graphics.views.camera.position=[1.,0.,10.];graphics.views.camera.target=[1.,0.,0.];graphics.views.camera.field=[2.15,1.0078125];graphics.picture.save_picture(file_name=str(path))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--prepare-check',action='store_true');a=p.parse_args();a.authorization=a.authorization.resolve();c,auth,proof=preflight(a.authorization,a.out)
 if a.prepare_check:print(json.dumps(dict(preflight_passed=True,solver_launched=False,**proof),indent=2));return 0
 import numpy as np
 out=a.out.resolve();out.mkdir();[(out/n).mkdir() for n in ['states','force-reports','temp','examples']]
 for key in ['TEMP','TMP','TMPDIR']:os.environ[key]=str(out/'temp')
 tempfile.tempdir=str(out/'temp');os.environ['PYTHONDONTWRITEBYTECODE']='1';os.chdir(out)
 for src,name in [(Path(__file__),'executed-pilot.py'),(Path(c.__file__),'executed-common.py'),(a.authorization,'authorization.json'),(LOCAL/'TARGETS.md','TARGETS.md'),(c.MESH.parent/'setup.json','mesh-setup.json')]:shutil.copyfile(src,out/name)
 c.write_json(out/'code-provenance.json',proof);events=[];rows=[];session=None;previous_iteration=0
 def update(stage,**values):
  events.append(dict(stage=stage,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**values));c.write_json(out/'status.json',events[-1]);c.write_json(out/'events.json',events);print(json.dumps(events[-1]),flush=True)
 try:
  import ansys.fluent.core as pf
  pf.config.examples_path=str(out/'examples');update('launching_one_hidden_session',processor_count=2)
  session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True);update('native_session_connected',version=str(session.get_fluent_version()),licence_status='Actual Fluent startup succeeded; launcher/transcript retained')
  session.transcript.start(file_name=str(out/'solver-transcript.txt'));st=session.settings;st.file.read_mesh(file_name=str(c.MESH));st.setup.general.solver.type='pressure-based';st.setup.general.solver.two_dim_space='planar';st.setup.general.solver.time='transient'
  st.setup.models.energy.enabled=False;st.setup.models.viscous.model='laminar';st.setup.general.operating_conditions.gravity.enable=False;st.setup.general.operating_conditions.operating_pressure=0.;st.setup.dynamic_mesh.enabled=False
  air=st.setup.materials.fluid['air'];air.density.value=1.;air.viscosity.value=.0001;zone=st.setup.cell_zone_conditions.fluid['fluid'];zone.general.material='air';zone.reference_frame.frame_motion=False;zone.mesh_motion.enable=False
  inlet=st.setup.boundary_conditions.velocity_inlet['inlet'].momentum;inlet.velocity_specification_method='Components';inlet.velocity_components=[dict(value=.1),dict(value=0.)];st.setup.boundary_conditions.pressure_outlet['outlet'].momentum.gauge_pressure.value=0.
  for name in c.BODIES:wall=st.setup.boundary_conditions.wall[name].momentum;wall.wall_motion='Stationary Wall';wall.shear_condition='No Slip'
  for key,value in dict(density=1.,viscosity=.0001,velocity=.1,length=.1,area=.1,depth=1.,pressure=0.).items():setattr(st.setup.reference_values,key,value)
  methods=st.solution.methods;methods.p_v_coupling.flow_scheme='SIMPLE';methods.spatial_discretization.gradient_scheme='least-square-cell-based';methods.spatial_discretization.discretization_scheme={'mom':'first-order-upwind','pressure':'standard'};methods.transient_formulation='unsteady-1st-order'
  run=st.solution.run_calculation;run.transient_controls.type='Fixed';run.transient_controls.method='User-Specified';run.transient_controls.duration_specification_method='incremental-time-steps';run.parameters.time_step_size=.05;run.parameters.max_iter_per_time_step=100
  residual=st.solution.monitor.residual;residual.options.criterion_type='absolute';residual.options.n_save=1000
  for name in c.EQUATIONS:residual.equations[name].absolute_criteria=1e-6;residual.equations[name].check_convergence=True
  check=st.mesh.check();q=st.mesh.quality()
  for attempt in range(20):
   try:qm=quality(str(check)+'\n'+str(q)+'\n'+(out/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace'),c);break
   except ValueError:
    if attempt==19:raise
    time.sleep(.1)
  c.write_json(out/'geometry-quality.json',qm);loaded=snapshot(st);c.settings_gate(loaded);c.write_json(out/'loaded-settings.json',loaded)
  init=st.solution.initialization;init.initialization_type='standard'
  for name in ['x-velocity','y-velocity','pressure']:init.defaults[name]=0.
  init.standard_initialize();c.write_json(out/'initialized-settings.json',snapshot(st));g=np.load(c.MESH.parent/'geometry.npz',allow_pickle=False);cache=None
  for step in range(6):
   if step:update('advancing_one_actual_step',requested_step=step,max_inner_iterations=100);run.dual_time_iterate(time_step_count=1,max_iter_per_step=100)
   case=out/'states'/f'state-{step:04d}.cas.h5';data=case.with_name(f'state-{step:04d}.dat.h5');st.file.write_case_data(file_name=str(case));c.require(case.is_file() and data.is_file(),'Native case/data absent')
   if cache is None:cache=c.map_native_mesh(case,g)
   s=c.read_state(case,data,g,cache);gate=c.state_gates(s,step);row=dict(step=step,time_s=s['native_time']['flow-time'],case=str(case.relative_to(out)),data=str(data.relative_to(out)),case_sha256=c.sha(case),data_sha256=c.sha(data),native_gates=gate,initial_forces_available=False if step==0 else None)
   rows.append(row);c.write_json(out/'states.json',rows)
   if step:
    rr,index=c.terminal_residuals(data,out/'solver-transcript.txt');c.require(1<=index-previous_iteration<=100,'Actual inner cap exceeded');previous_iteration=index;gate['checks']['own_terminal_residuals']=bool(all(v<=1e-6 for v in rr.values()));forces=c.body_forces(s);row.update(native_terminal_iteration=index,terminal_scaled_residuals=rr,body_forces=forces)
    if step==5:row['native_force_reports']=export_force_reports(st,out,step,forces,c)
   gate['all_checks_passed']=bool(all(gate['checks'].values()));c.write_json(out/'states.json',rows);update('actual_state_saved',step=step,time_s=s['native_time']['flow-time'],all_state_gates_passed=gate['all_checks_passed']);c.require(gate['all_checks_passed'],'Pilot numerical gate failed; actual native state preserved')
  c.write_json(out/'final-settings.json',snapshot(st));residual.write(file_name=str(out/'residuals.xy'))
  try:native_picture(st,out/'native-velocity.png')
  except Exception:(out/'graphics-export-error.txt').write_text(traceback.format_exc(),encoding='utf-8')
  update('pilot_completed',actual_steps=5,actual_end_time_s=.25,native_states=6,analysis_pending=True,production_authorized=False);return 0
 except Exception as error:(out/'failure.txt').write_text(traceback.format_exc(),encoding='utf-8');update('failed_preserved',error=str(error),native_states=len(rows));raise
 finally:
  if session is not None:
   try:session.transcript.stop()
   finally:session.exit()
if __name__=='__main__':raise SystemExit(main())
