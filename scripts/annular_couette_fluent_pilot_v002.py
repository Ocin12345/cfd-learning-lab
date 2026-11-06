"""Root-authorized single-session Fluent annular spin-up PILOT ONLY:5steps,0.25s.
No launch occurs on import or --prepare-check. Production has no entry point here.
"""
from pathlib import Path
import argparse,ctypes,datetime,hashlib,json,math,os,re,shutil,tempfile,time,traceback
import numpy as np
import annular_couette_fluent_common_v002 as c
ROOT=Path('E:/CFD-Learning-Lab');LOCAL=ROOT/'local/annular-couette-v001';RUN=ROOT/'native/runs/2026-11-02-annular-couette-v001';MESH=RUN/'mesh/annular-couette.msh'
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
ACTION='annular_couette_fluent_pilot_v002'
def available_memory():
 class MS(ctypes.Structure):
  _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ['total','available','total_page','available_page','total_virtual','available_virtual','available_extended']]
 m=MS();m.length=ctypes.sizeof(m);c.require(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)),'Cannot query actual RAM');return m.available
def preflight(authorization,out):
 auth=c.read(authorization);out=Path(out).resolve()
 c.require(auth.get('approved') is True and auth.get('root')=='root' and auth.get('action')==ACTION and auth.get('phase')=='pilot-5-steps','No root five-step pilot GO')
 c.require(Path(auth['out']).resolve()==out and out.parent==RUN.resolve() and not out.exists(),'Fresh authorized direct native output required')
 c.require(auth['solver_processes']==2 and auth['pilot_steps']==5 and auth['dt_s']==.05,'Authorization changes pilot cap/ranks/dt')
 sources={Path(__file__).resolve():auth['source_sha256'],Path(c.__file__).resolve():auth['common_sha256'],MESH:c.MESH_SHA,MESH.parent/'geometry.npz':c.GEOMETRY_SHA,MESH.parent/'setup.json':c.SETUP_SHA,LOCAL/'TARGETS.md':c.TARGETS_SHA}
 c.require(auth['mesh_sha256']==c.MESH_SHA and auth['geometry_sha256']==c.GEOMETRY_SHA and auth['targets_sha256']==c.TARGETS_SHA,'Authorization binds different source/targets')
 for path,digest in sources.items():c.require(path.is_file() and c.sha(path)==digest,'Changed pilot source '+str(path))
 c.require(EXECUTABLE.is_file(),'Installed Fluent missing');ram=available_memory();disk=shutil.disk_usage('E:/').free
 c.require(ram>=4*1024**3 and disk>=10*1024**3,'Pilot requires unchanged4GiB RAM/10GiB E admission')
 return dict(root_authorization_sha256=c.sha(authorization),inputs=[dict(path=str(p),sha256=h) for p,h in sources.items()],available_ram_bytes=ram,free_E_bytes=disk,processor_count=2,source_sha256=c.sha(__file__),mesh_sha256=c.MESH_SHA,geometry_npz_sha256=c.GEOMETRY_SHA,snapshot_before_solver_launch=True)
def quality(text):
 a=re.findall(r'Minimum Orthogonal Quality\s*=\s*([-+0-9.eE]+)',text,re.I);b=re.findall(r'minimum volume\s*\([^)]*\)\s*:\s*([-+0-9.eE]+)',text,re.I)
 c.require(a and b,'Missing native quality numbers');orth=float(a[-1]);volume=float(b[-1]);c.require(math.isfinite(orth) and orth>=.15 and math.isfinite(volume) and volume>0,'Native mesh quality failed')
 return dict(minimum_orthogonal_quality=orth,minimum_native_cell_volume=volume,native_numeric_gate_pass=True,checked_before_flow=True)
def active_read(setting):
 active=setting.is_active();return dict(active=active,state=setting.get_state() if active else None)
def configure_rotation(inner):
 c.require('Moving Wall' in inner.wall_motion.allowed_values() and 'Absolute' in inner.relative.allowed_values() and 'Rotational' in inner.velocity_spec.allowed_values(),'Rotational wall API unavailable')
 inner.wall_motion='Moving Wall';inner.relative='Absolute';inner.shear_condition='No Slip';inner.velocity_spec='Rotational';inner.rotation_speed.option='value'
 prop=inner.rotation_speed.value;c.require(prop.is_active(),'Rotational speed inactive');unit=prop.units();quantity=prop.get_attr('units-quantity')
 c.require(unit is not None,'Missing live angular-speed units');prop.set_state((1.,'rad/s'))
 q=prop.as_quantity();c.require(q is not None and abs(float(q.to('rad/s').value)-1)<1e-12,'Actual angular speed is not1rad/s')
 origin=inner.rotation_axis_origin;c.require(origin.is_active(),'Rotation origin inactive');old=origin.get_state();c.require(len(old) in [2,3],'Unexpected2D axis origin length');origin.set_state([0.]*len(old))
 axis=inner.rotation_axis_direction;axis_active=axis.is_active()
 if axis_active:
  oldaxis=axis.get_state();c.require(len(oldaxis)==3,'Active rotation axis must expose3components');axis.set_state([0.,0.,1.])
 return dict(wall_motion=inner.wall_motion(),relative=inner.relative(),velocity_spec=inner.velocity_spec(),shear_condition=inner.shear_condition(),
  rotation_speed=dict(active=prop.is_active(),native_value=float(prop()),native_units=unit,units_quantity=quantity,rad_s=float(prop.as_quantity().to('rad/s').value)),
  origin=active_read(origin),axis=active_read(axis),axis_if_inactive='2D planar implicit z sign must be proved by actual saved wall U,V',
  allowed=dict(wall_motion=inner.wall_motion.allowed_values(),relative=inner.relative.allowed_values(),velocity_spec=inner.velocity_spec.allowed_values()))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--prepare-check',action='store_true');args=p.parse_args()
 provenance=preflight(args.authorization,args.out)
 if args.prepare_check:print(json.dumps(dict(preflight_passed=True,solver_launched=False,**provenance),indent=2));return 0
 out=args.out.resolve();out.mkdir();(out/'temp').mkdir();(out/'examples').mkdir();(out/'states').mkdir()
 for n in ['TEMP','TMP','TMPDIR']:os.environ[n]=str(out/'temp')
 tempfile.tempdir=str(out/'temp');os.environ['PYTHONDONTWRITEBYTECODE']='1';os.chdir(out)
 shutil.copyfile(__file__,out/'executed-pilot.py');shutil.copyfile(c.__file__,out/'executed-common.py');shutil.copyfile(args.authorization,out/'authorization.json')
 shutil.copyfile(LOCAL/'TARGETS.md',out/'TARGETS.md');shutil.copyfile(MESH.parent/'setup.json',out/'mesh-setup.json')
 c.write_json(out/'code-provenance.json',provenance);c.write_json(out/'declared-targets.json',c.TARGETS)
 session=None;rows=[];events=[];previous_iteration=0
 def update(stage,**values):
  events.append(dict(stage=stage,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**values));c.write_json(out/'status.json',events[-1]);c.write_json(out/'events.json',events)
 try:
  import ansys.fluent.core as pf
  pf.config.examples_path=str(out/'examples');update('launching_one_hidden_session',processor_count=2)
  session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
  update('native_session_connected',version=str(session.get_fluent_version()),licence_status='Actual Fluent startup succeeded; native launcher/transcript retained; no licence file copied')
  session.transcript.start(file_name=str(out/'solver-transcript.txt'));st=session.settings;st.file.read_mesh(file_name=str(MESH))
  st.setup.general.solver.type='pressure-based';st.setup.general.solver.two_dim_space='planar';st.setup.general.solver.time='transient'
  st.setup.models.energy.enabled=False;st.setup.models.viscous.model='laminar';st.setup.general.operating_conditions.gravity.enable=False
  st.setup.dynamic_mesh.enabled=False
  op=st.setup.general.operating_conditions;op.operating_pressure=0.;g=np.load(MESH.parent/'geometry.npz');c.geometry_check(g)
  ref=g['centres'][0].tolist();op.reference_pressure_location=ref
  air=st.setup.materials.fluid['air'];air.density.value=1.;air.viscosity.value=.000125
  zone=st.setup.cell_zone_conditions.fluid['fluid'];zone.general.material='air';zone.reference_frame.frame_motion=False;zone.mesh_motion.enable=False
  inner=st.setup.boundary_conditions.wall['inner'].momentum;rotation=configure_rotation(inner)
  outer=st.setup.boundary_conditions.wall['outer'].momentum;outer.wall_motion='Stationary Wall';outer.shear_condition='No Slip'
  for key,value in dict(density=1.,viscosity=.000125,velocity=.05,length=.05,area=1.,depth=1.,pressure=0.).items():setattr(st.setup.reference_values,key,value)
  methods=st.solution.methods;methods.p_v_coupling.flow_scheme='SIMPLE';methods.spatial_discretization.gradient_scheme='least-square-cell-based'
  methods.spatial_discretization.discretization_scheme={'mom':'first-order-upwind','pressure':'standard'};methods.transient_formulation='unsteady-1st-order'
  run=st.solution.run_calculation;run.transient_controls.type='Fixed';run.transient_controls.method='User-Specified';run.transient_controls.duration_specification_method='incremental-time-steps';run.parameters.time_step_size=.05;run.parameters.max_iter_per_time_step=100
  residual=st.solution.monitor.residual;residual.options.criterion_type='absolute';residual.options.n_save=1000
  for name in c.EQUATIONS:residual.equations[name].absolute_criteria=1e-6;residual.equations[name].check_convergence=True
  check=st.mesh.check();q=st.mesh.quality()
  for attempt in range(20):
   try:native_quality=quality(str(check)+'\n'+str(q)+'\n'+(out/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace'));break
   except ValueError:
    if attempt==19:raise
    time.sleep(.1)
  c.write_json(out/'geometry-quality.json',native_quality)
  def state():return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),dynamic_mesh=st.setup.dynamic_mesh.get_state(),material=air.get_state(),cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),walls={n:st.setup.boundary_conditions.wall[n].get_state() for n in ['inner','outer']},references=st.setup.reference_values.get_state(),methods=methods.get_state(),run_calculation=run.get_state(),residual=residual.get_state())
  readback=dict(rotation=rotation,pressure_reference_location=op.reference_pressure_location(),pressure_reference_source_cell=0,pressure_reference_used_native_command=op.used_ref_pressure_location(),
   material=air.get_state(),reference_depth_m=st.setup.reference_values.depth(),frame_motion=zone.reference_frame.frame_motion(),mesh_motion=zone.mesh_motion.enable(),
   dynamic_mesh=st.setup.dynamic_mesh.enabled(),gravity=op.gravity.enable(),energy=st.setup.models.energy.enabled(),planar=st.setup.general.solver.two_dim_space())
  c.require(np.allclose(readback['pressure_reference_location'],ref,rtol=0,atol=1e-14) and not any(readback[n] for n in ['frame_motion','mesh_motion','dynamic_mesh','gravity','energy']) and readback['reference_depth_m']==1.,'Actual inertial/pressure/material controls differ')
  c.write_json(out/'native-physical-readback.json',readback);c.write_json(out/'loaded-settings.json',state())
  init=st.solution.initialization;init.initialization_type='standard'
  for name in ['x-velocity','y-velocity','pressure']:init.defaults[name]=0.
  init.standard_initialize();c.write_json(out/'initialized-settings.json',state())
  for step in range(6):
   if step:
    update('advancing_one_actual_step',requested_step=step,max_inner_iterations=100)
    run.dual_time_iterate(time_step_count=1,max_iter_per_step=100)
   stem=f'state-{step:04d}';case=out/'states'/(stem+'.cas.h5');data=out/'states'/(stem+'.dat.h5');st.file.write_case_data(file_name=str(case))
   c.require(case.is_file() and data.is_file(),'Actual native case/data save failed');snap=c.read_state(case,data,g);gates=c.state_gates(snap,step)
   c.save_arrays(out/'states'/(stem+'-fields.npz'),snap)
   row=dict(step=step,time_s=step*.05,case=str(case.relative_to(out)),data=str(data.relative_to(out)),case_sha256=c.sha(case),data_sha256=c.sha(data),fields=stem+'-fields.npz',native_gates=gates)
   if step:
    rr,iteration=c.terminal_residuals(data,out/'solver-transcript.txt');c.require(1<=iteration-previous_iteration<=100,'Native step inner count outside1..100');previous_iteration=iteration
    row.update(native_terminal_iteration=iteration,terminal_scaled_residuals=rr);row['native_gates']['checks']['own_terminal_residuals']=all(v<=1e-6 for v in rr.values());row['native_gates']['all_checks_passed']=all(row['native_gates']['checks'].values())
   rows.append(row);c.write_json(out/'states.json',rows);update('actual_state_saved',step=step,time_s=step*.05,all_state_gates_passed=row['native_gates']['all_checks_passed'])
   c.require(row['native_gates']['all_checks_passed'],'Pilot native numerical/structural gate failed; saved evidence preserved')
  c.write_json(out/'final-settings.json',state());residual.write(file_name=str(out/'residuals.xy'))
  graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1200;graphics.picture.y_resolution=1200;graphics.picture.invert_background=False
  graphics.contour.create(name='annular-speed');contour=graphics.contour['annular-speed'];contour.field='velocity-magnitude';contour.surfaces_list=[];contour.display()
  graphics.views.camera.projection='orthographic';graphics.views.camera.position=[0.,0.,10.];graphics.views.camera.target=[0.,0.,0.];graphics.views.camera.field=[.225,.225]
  graphics.picture.save_picture(file_name=str(out/'native-velocity.png'))
  c.require((out/'native-velocity.png').is_file(),'Native velocity PNG missing');update('pilot_completed',actual_steps=5,actual_end_time_s=.25,native_states=6,production_executed=False,numerical_analysis_pending=True)
  return 0
 except Exception as error:
  (out/'failure.txt').write_text(traceback.format_exc(),encoding='utf-8');update('failed_preserved',error=str(error),native_states=len(rows));raise
 finally:
  if session is not None:
   try:session.transcript.stop()
   finally:session.exit()
if __name__=='__main__':raise SystemExit(main())
