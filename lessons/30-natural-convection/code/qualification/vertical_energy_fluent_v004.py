"""Compare computed energy slopes with and without gravity in an adiabatic slip duct."""
from pathlib import Path
import json,os,shutil,traceback
import numpy as np
import ansys.fluent.core as pf
from cavity_buoyancy_common_v003 import *
from wing_fluent_v001 import available_memory,mesh_quality
from vertical_energy_reader_v002 import fluent_native
OUT=ROOT/'native/runs/vertical-energy-control-fluent-v003'
MESH=ROOT/'native/meshes/vertical-energy-control-v001/cavity.msh'

def main():
 assert not OUT.exists();assert available_memory()>4*2**30 and shutil.disk_usage(ROOT).free>10*2**30
 OUT.mkdir(parents=True);shutil.copyfile(__file__,OUT/'executed-driver.py');s=None
 def mark(stage,**kw):write(OUT/'status.json',dict(stage=stage,**kw));print(stage,flush=True)
 try:
  mark('launching');s=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True)
  s.transcript.start(file_name=str(OUT/'solver-transcript.txt'));st=s.settings;st.file.read_mesh(file_name=str(MESH));mark('licensed_connected')
  st.setup.general.solver.type='pressure-based';st.setup.general.solver.time='steady';st.setup.models.viscous.model='laminar';st.setup.models.energy.enabled=True
  gas=st.setup.materials.fluid['air'];gas.density.option='ideal-gas';gas.molecular_weight.value=MW
  for name,value in [('specific_heat',CP),('viscosity',MU),('thermal_conductivity',K)]:q=getattr(gas,name);q.option='constant';q.value=value
  st.setup.cell_zone_conditions.fluid['fluid'].general.material='air';op=st.setup.general.operating_conditions;op.operating_pressure=P0;op.gravity.enable=True;op.gravity.components=[0.,-G];op.operating_density.enable=True;op.operating_density.value=0.
  for name,value in [('pressure_work',True),('kinetic_energy',True),('viscous_dissipation',False)]:
   q=getattr(st.setup.models.energy,name)
   if q.is_active():q.set_state(value)
  bc=st.setup.boundary_conditions
  for zones,kind in [(['bottom'],'velocity-inlet'),(['top'],'pressure-outlet'),(['left','right'],'symmetry')]:bc.set_zone_type(zone_list=zones,new_type=kind)
  inlet=bc.velocity_inlet['bottom'];inlet.momentum.velocity_specification_method='Components';inlet.momentum.velocity_components=[{'value':0.},{'value':1.}];inlet.thermal.temperature.value=T0
  sol=st.solution;sol.methods.p_v_coupling.flow_scheme='SIMPLE';sol.methods.spatial_discretization.gradient_scheme='least-square-cell-based';sol.methods.spatial_discretization.discretization_scheme={'pressure':'body-force-weighted','mom':'second-order-upwind','temperature':'second-order-upwind'}
  res=sol.monitor.residual;res.options.n_save=3000
  for name in res.equations.get_object_names():res.equations[name].absolute_criteria=1e-10 if name=='energy' else 1e-8;res.equations[name].check_convergence=False
  st.mesh.check();st.mesh.quality();write(OUT/'mesh-quality.json',mesh_quality((OUT/'solver-transcript.txt').read_text(errors='replace')))
  for label,g,iters in [('with-gravity',9.81,1000)]:
   op.gravity.enable=bool(g)
   if g:op.gravity.components=[0.,-g]
   init=sol.initialization;init.initialization_type='standard'
   for name,value in [('x-velocity',0.),('y-velocity',1.),('pressure',0.),('temperature',T0)]:init.defaults[name]=value
   init.standard_initialize();mark('solving',label=label);sol.run_calculation.iterate(iter_count=iters)
   st.file.write_case_data(file_name=str(OUT/(label+'.cas.h5')))
   st.results.report.fluxes.mass_flow(zones=['top','bottom'],write_to_file=True,file_name=str(OUT/(label+'-mass.txt')))
   q=fluent_native(OUT/(label+'.cas.h5'),OUT/(label+'.dat.h5'),P0);xy=q['nodes'][q['cells'],:2].mean(axis=1);H=CP*q['T']+.5*(q['U']**2).sum(axis=1)
   slope=float(np.polyfit(xy[:,1],H,1)[0]);write(OUT/(label+'-energy-check.json'),dict(g=g,fitted_dH_dy_J_kg_m=slope,expected_if_gravitational_work_included=-g,expected_if_no_gravitational_work=0.,minimum_T_K=float(q['T'].min()),maximum_T_K=float(q['T'].max()),u_min_m_s=float(q['U'][:,1].min()),u_max_m_s=float(q['U'][:,1].max()),native_residuals=q['solve_evidence']))
   gr=st.results.graphics;gr.picture.use_window_resolution=False;gr.picture.x_resolution=800;gr.picture.y_resolution=1200
   if 'temperature' not in gr.contour.get_object_names():gr.contour.create(name='temperature')
   c=gr.contour['temperature'];c.field='temperature';c.surfaces_list=[];c.options.filled=True;c.display();gr.views.camera.projection='orthographic';gr.views.camera.position=[.5,5,30];gr.views.camera.target=[.5,5,0];gr.views.camera.field=[8.,12.];gr.picture.save_picture(file_name=str(OUT/(label+'-temperature.png')))
   write(OUT/(label+'-settings.json'),dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),material=gas.get_state(),boundaries=bc.get_state(),residual=res.get_state()))
 except Exception:
  mark('failed',error=traceback.format_exc());raise
 finally:
  if s:s.exit()
 mark('closed')
if __name__=='__main__':main()
