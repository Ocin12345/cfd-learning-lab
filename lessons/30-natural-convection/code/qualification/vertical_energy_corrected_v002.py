from pathlib import Path
import json,os,shutil,traceback
import numpy as np
import ansys.fluent.core as pf
from cavity_buoyancy_common_v003 import *
from wing_fluent_v001 import available_memory
from vertical_energy_reader_v002 import fluent_native
from fluent_gravity_work_v001 import install_gravity_work
PARENT=ROOT/'native/runs/vertical-energy-control-fluent-v003';OUT=ROOT/'native/runs/vertical-energy-control-fluent-v005'
assert not OUT.exists();assert available_memory()>4*2**30 and shutil.disk_usage(ROOT).free>10*2**30
OUT.mkdir(parents=True);shutil.copyfile(__file__,OUT/'executed-driver.py');s=None
try:
 write(OUT/'status.json',dict(stage='launching'))
 s=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True)
 s.transcript.start(file_name=str(OUT/'solver-transcript.txt'));st=s.settings;st.file.read_case_data(file_name=str(PARENT/'with-gravity.cas.h5'))
 write(OUT/'gravity-work-source.json',install_gravity_work(st))
 res=st.solution.monitor.residual;res.options.criterion_type='absolute'
 for name in res.equations.get_object_names():res.equations[name].check_convergence=True;res.equations[name].absolute_criteria=1e-10 if name=='energy' else 1e-8;res.equations[name].check_convergence=False
 write(OUT/'status.json',dict(stage='solving'));st.solution.run_calculation.iterate(iter_count=800)
 st.file.write_case_data(file_name=str(OUT/'corrected.cas.h5'));st.results.report.fluxes.mass_flow(zones=['top','bottom'],write_to_file=True,file_name=str(OUT/'mass.txt'))
 q=fluent_native(OUT/'corrected.cas.h5',OUT/'corrected.dat.h5',P0);xy=q['nodes'][q['cells'],:2].mean(axis=1);H=CP*q['T']+.5*(q['U']**2).sum(axis=1);slope=float(np.polyfit(xy[:,1],H,1)[0])
 write(OUT/'energy-check.json',dict(fitted_dH_dy_J_kg_m=slope,expected=-G,slope_relative_error=abs(slope/G+1),minimum_T_K=float(q['T'].min()),maximum_T_K=float(q['T'].max()),native_residuals=q['solve_evidence'],gravity_work_source=st.setup.cell_zone_conditions.fluid['fluid'].sources.get_state()))
 gr=st.results.graphics;gr.picture.use_window_resolution=False;gr.picture.x_resolution=800;gr.picture.y_resolution=1200
 if 'temperature' not in gr.contour.get_object_names():gr.contour.create(name='temperature')
 c=gr.contour['temperature'];c.field='temperature';c.surfaces_list=[];c.options.filled=True;c.display();gr.views.camera.projection='orthographic';gr.views.camera.position=[.5,5,30];gr.views.camera.target=[.5,5,0];gr.views.camera.field=[8.,12.];gr.picture.save_picture(file_name=str(OUT/'temperature-native.png'))
 assert abs(slope/G+1)<.02,'Energy slope differs by more than 2% from source prediction'
except Exception:
 write(OUT/'status.json',dict(stage='failed',error=traceback.format_exc()));raise
finally:
 if s:s.exit()
write(OUT/'status.json',dict(stage='closed'))
