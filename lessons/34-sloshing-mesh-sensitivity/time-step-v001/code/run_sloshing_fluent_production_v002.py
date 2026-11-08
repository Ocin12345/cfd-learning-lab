from pathlib import Path
from datetime import datetime
import os,json,hashlib,shutil,time,traceback
import ansys.fluent.core as pf
from sloshing_dt_common_v001 import *
from sloshing_prelaunch_v001 import prelaunch,verify
from sloshing_fixed_settings_v001 import compare_fixed
from sloshing_wait_admission_v001 import wait_admission
from sloshing_fluent_steps_v001 import parse_steps
from sloshing_full_audits_v002 import guardF,auditF,fclock
from wing_fluent_v001 import available_memory
P=N/'fluent-qualification-v001';O=N/'fluent-production-v001'
def mark(stage,**v):
 p=S/'fluent-production-v001-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def main():
 assert not (S/'fluent-production-v001-status.json').exists() and not (O/'solver-transcript.txt').exists()
 admission=wait_admission('full-interval-inputs-v002.json',lambda **v:mark('waiting_fresh_prelaunch_gate',**v));mark('launching_exact_checkpoint_continuation',admission=admission);session=None;loaded=False
 try:
  session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(O),start_timeout=180,cleanup_on_exit=True)
  mark('licensed_connected_before_new_flow',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(O/'solver-transcript.txt'));st=session.settings
  st.file.read_case_data(file_name=str(O/'restart-source.cas.h5'));loaded=True;run=st.solution.run_calculation
  assert abs(float(run.transient_controls.flow_time())-.02)<1e-10 and abs(float(run.parameters.time_step_size())-DT)<1e-14 and run.parameters.max_iter_per_time_step()==400
  fixed=dict(model=st.setup.models.multiphase.get_state(),operating_conditions=st.setup.general.operating_conditions.get_state(),materials={k:st.setup.materials.fluid[k].get_state() for k in MATERIALS},numerical_methods=st.solution.methods.get_state(),residual_options=st.solution.monitor.residual.options.get_state(),residual_equations=st.solution.monitor.residual.equations.get_state())
  write_json(O/'loaded-fixed-settings.json',fixed);write_json(O/'fixed-settings-comparison.json',compare_fixed(fixed,json.loads((P/'loaded-fixed-settings.json').read_text())))
  shutil.copyfile(__file__,O/'executed-driver.py');st.mesh.check();st.mesh.quality();st.file.write_case(file_name=str(O/'active.cas.h5'));frames=[];previous=80
  write_json(O/'declared-settings.json',dict(start_step=80,last_step=4800,new_steps=4720,dt_s=DT,initialization='exactaccepted.02scheckpoint,no field reinitialization',frame_steps=[80,81]+list(range(120,4801,40))))
  for step in [80,81]+list(range(120,4801,40)):
   if step>previous:run.dual_time_iterate(time_step_count=step-previous,max_iter_per_step=400)
   p=O/f'state-{step:04d}.dat.h5';st.file.write_data(file_name=str(p));frames.append(guardF(O/'active.cas.h5',p,O/'geometry.npz',step));write_json(O/'frames.json',frames)
   if step>80:
    rows=parse_steps((O/'solver-transcript.txt').read_text(errors='replace'),80,step,DT);write_json(O/'all-step-checks.json',rows);assert all(r['passed'] for r in rows),'Native state saved before residual miss stop'
   previous=step;mark('native_checkpoint_pass',step=step,time_s=step*DT,states=len(frames))
  st.file.write_case_data(file_name=str(O/'continuation-final.cas.h5'))
 except Exception:
  if loaded:
   try:st.file.write_case_data(file_name=str(O/'recovery.cas.h5'))
   except Exception:pass
  raise
 finally:
  if session:session.exit()
 verify('full-interval-inputs-v002.json');mark('auditing_closed_full_interval');a=auditF(O);assert a['numerical_pass'];mark('fluent_full_interval_pass',new_steps=4720,combined_steps=4800,OpenFOAM_full_not_yet_run=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise
