import json,os,shutil,csv,traceback
import numpy as np,h5py
import ansys.fluent.core as pf
from jet_energy_common_v001 import *
from jet_energy_full_prelaunch_v001 import wait_admission,verify
from sloshing_fixed_settings_v001 import compare_fixed
from analyse_heated_jet_v002 import values,N as NC,NF
O=N/'fluent-production-v001'
def main():
 assert not (O/'solver-transcript.txt').exists();admission=wait_admission(lambda **v:mark(O/'status.json','waiting_admission',**v));dump(O/'launch-admission.json',admission)
 session=None;loaded=False
 try:
  session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(O),start_timeout=180,cleanup_on_exit=True)
  mark(O/'status.json','licensed_connected_before_any_advance',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(O/'solver-transcript.txt'));st=session.settings
  st.file.read_case(file_name=str(O/'source.cas.h5'));st.file.read_data(file_name=str(O/'source-0020.dat.h5'));loaded=True
  run=st.solution.run_calculation;assert abs(float(run.transient_controls.flow_time())-.1)<1e-10
  actual=dict(energy=st.setup.models.energy.get_state(),material=st.setup.materials.fluid['air'].get_state(),equations=st.solution.controls.get_state()['equations'],methods=st.solution.methods.get_state(),plate=st.setup.boundary_conditions.wall['plate'].thermal.get_state(),ceiling=st.setup.boundary_conditions.wall['ceiling'].thermal.get_state())
  base=json.loads((O/'expected-setup.json').read_text());expected={k:base[k] for k in actual};dump(O/'loaded-fixed-settings.json',actual);dump(O/'fixed-settings-comparison.json',compare_fixed(actual,expected))
  assert run.parameters.time_step_size()==.005 and run.parameters.max_iter_per_time_step()==30
  residual=st.solution.monitor.residual;assert residual.equations.get_object_names()==['energy'];assert residual.equations['energy'].check_convergence() and residual.equations['energy'].absolute_criteria()==1e-8
  for name in ['left','right']:assert st.setup.boundary_conditions.pressure_outlet[name].thermal.backflow_total_temperature.value()==300
  assert st.setup.boundary_conditions.velocity_inlet['inlet'].thermal.temperature.value()==300
  st.mesh.check();st.mesh.quality();dump(O/'run-settings.json',run.get_state());st.file.write_case(file_name=str(O/'active.cas.h5'))
  with (O/'frames.csv').open('x',newline='') as f:
   w=csv.writer(f);w.writerow(['step','time_s','file'])
   for step in range(20,601):
    if step>20:run.dual_time_iterate(time_step_count=1,max_iter_per_step=30)
    t=float(run.transient_controls.flow_time());assert abs(t-step*.005)<1e-10
    p=O/f'heat-{step:04d}.dat.h5';st.file.write_data(file_name=str(p));T=values(p,'SV_T');assert np.isfinite(T).all() and T.min()>=299.99 and T.max()<=320.01
    if step==20:
     for k,kind,n in [('SV_T','cells',NC),('SV_H','cells',NC),('SV_U','cells',NC),('SV_V','cells',NC),('SV_P','cells',NC),('SV_FLUX','faces',NF)]:assert np.array_equal(values(p,k,kind,n),values(O/'source-0020.dat.h5',k,kind,n)),k
    else:
     st.results.report.fluxes.heat_transfer(zones=['inlet','left','right','plate','ceiling'],write_to_file=True,file_name=str(O/f'heat-report-{step:04d}.txt'))
     st.results.report.fluxes.mass_flow(zones=['inlet','left','right','plate','ceiling'],write_to_file=True,file_name=str(O/f'mass-report-{step:04d}.txt'))
    w.writerow([step,t,p.name]);f.flush();mark(O/'status.json','native_step_saved',step=step,time_s=t)
  st.file.write_case_data(file_name=str(O/'production-final.cas.h5'))
 except Exception:
  if loaded:
   try:st.file.write_case_data(file_name=str(O/'recovery.cas.h5'))
   except Exception:pass
  mark(O/'status.json','failed',error=traceback.format_exc());raise
 finally:
  if session:session.exit()
 verify();mark(O/'status.json','closed',step=600,time_s=3.)
if __name__=='__main__':main()
