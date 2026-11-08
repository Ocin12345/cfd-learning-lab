from pathlib import Path
import os,json,shutil,csv,traceback
import numpy as np,h5py
import ansys.fluent.core as pf
from sloshing_dt_common_v001 import *
from sloshing_prelaunch_v001 import prelaunch,verify
from sloshing_fixed_settings_v001 import compare_fixed
from analyse_nozzle_fluent_v005 import common_mapping
O=N/'fluent-qualification-v001'
def main():
 assert not (O/'solver-transcript.txt').exists()
 resources=prelaunch('qualification-inputs-v001.json');write_json(O/'launch-admission.json',resources);session=None;loaded=False
 try:
  session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(O),start_timeout=180,cleanup_on_exit=True)
  write_json(O/'status.json',dict(stage='licensed_connected_before_flow',version=str(session.get_fluent_version())))
  session.transcript.start(file_name=str(O/'solver-transcript.txt'));st=session.settings
  st.file.read_case(file_name=str(O/'tank.cas.h5'));st.file.read_data(file_name=str(O/'tank-0000.dat.h5'));loaded=True
  assert abs(float(st.solution.run_calculation.transient_controls.flow_time()))<1e-14
  base=json.loads((B/'fluent/setup.json').read_text());mp=st.setup.models.multiphase
  fixed=dict(model=mp.get_state(),operating_conditions=st.setup.general.operating_conditions.get_state(),materials={k:st.setup.materials.fluid[k].get_state() for k in MATERIALS},numerical_methods=st.solution.methods.get_state(),residual_options=st.solution.monitor.residual.options.get_state(),residual_equations=st.solution.monitor.residual.equations.get_state())
  expected={k:base[k] for k in ['model','operating_conditions','materials','numerical_methods']};expected.update(residual_options=base['residual_reporting']['options'],residual_equations=base['residual_reporting']['equations'])
  write_json(O/'loaded-fixed-settings.json',fixed);write_json(O/'fixed-settings-comparison.json',compare_fixed(fixed,expected))
  assert st.setup.models.viscous.model()=='laminar'
  st.mesh.check();st.mesh.quality();run=st.solution.run_calculation
  assert run.parameters.max_iter_per_time_step()==400
  run.parameters.time_step_size=DT
  base.update(case_mode='qualification',total_steps=80,frame_interval=1,time_step_s=DT,run_calculation=run.get_state())
  write_json(O/'setup.json',base);shutil.copyfile(__file__,O/'executed-driver.py');st.file.write_case(file_name=str(O/'active.cas.h5'))
  with (O/'frames.csv').open('x',newline='') as f:
   w=csv.writer(f);w.writerow(['step','time_s','file','water_area_m2','alpha_min','alpha_max'])
   for step in range(81):
    if step:run.dual_time_iterate(time_step_count=1,max_iter_per_step=400)
    t=float(run.transient_controls.flow_time());assert abs(t-step*DT)<1e-10
    p=O/f'state-{step:04d}.dat.h5';st.file.write_data(file_name=str(p))
    with h5py.File(p) as g:a=cell_values(alpha_group(g));uv=[cell_values(g['results/1/phase-1/cells/'+k]) for k in ['SV_U','SV_V','SV_P']]
    assert all(np.isfinite(x).all() for x in uv);area=float(a.sum()*DX**2)
    w.writerow([step,t,p.name,area,float(a.min()),float(a.max())]);f.flush()
    if not step:
     with np.load(O/'geometry.npz') as g:
      m=common_mapping(O/'active.cas.h5',g);xy=g['centres'];expect=np.clip((.05+.025*(xy[:,0]-.1)-(xy[:,1]-DX/2))/DX,0,1)
     aligned=np.empty_like(a);aligned[m]=a;err=float(abs(aligned-expect).max());assert err<1e-10
     write_json(O/'initial-field-check.json',dict(passed=True,maximum_cell_fraction_difference=err,water_area_m2=area))
    assert a.min()>=-1e-6 and a.max()<=1+1e-6 and abs(area/.01-1)<=.005
    write_json(O/'status.json',dict(stage='native_state_saved',step=step,time_s=t))
  st.file.write_case_data(file_name=str(O/'tank-final.cas.h5'))
 except Exception:
  if loaded:
   try:st.file.write_case_data(file_name=str(O/'recovery.cas.h5'))
   except Exception:pass
  write_json(O/'status.json',dict(stage='failed',error=traceback.format_exc()));raise
 finally:
  if session:session.exit()
 verify('qualification-inputs-v001.json');write_json(O/'status.json',dict(stage='closed',final_time_s=.02))
if __name__=='__main__':main()
