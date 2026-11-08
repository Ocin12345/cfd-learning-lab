from pathlib import Path
from datetime import datetime
import json,os,hashlib,shutil,time,traceback
import numpy as np
import ansys.fluent.core as pf
from cavity_spatial_common_v001 import *
from cavity_prelaunch_v001 import prelaunch,verify
from cavity_continued_steps_v001 import audit_continued_steps
from cavity_residual_settings_v001 import configure_residuals
from cavity_native_heat_v001 import native_fluent_time,fluent_wall_flux
from cavity_native_fields_v001 import fluent,metrics
from run_cavity_fluent_qualification_v004 import saved_stats
from cavity_restart_settings_v001 import compare_restart_settings
def check_settings(actual):
 return compare_restart_settings(actual,json.loads((P/'final-settings.json').read_text()))
from wing_fluent_v001 import available_memory
S=Path(__file__).parent;NATIVE=ROOT/'native/runs/showcase-convection-spatial-v001';P=NATIVE/'fluent-qualification-v004';OUT=NATIVE/'fluent-production-v002'
def mark(stage,**v):
 p=S/'fluent-production-v002-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def state(st):
 return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),material=st.setup.materials.fluid['air'].get_state(),boundaries=st.setup.boundary_conditions.get_state(),references=st.setup.reference_values.get_state(),methods=st.solution.methods.get_state(),controls=st.solution.controls.get_state(),residual=st.solution.monitor.residual.get_state())
def guard(path,step,m0):
 clock=native_fluent_time(path);assert clock['step']==step and abs(clock['time_s']-step*.01)<1e-8 and abs(clock['dt_s']-.01)<1e-12
 row=saved_stats(path);row.update(step=step,time_s=clock['time_s'],file=path.name,mass_relative_drift=abs(row['mass_per_depth_kg_m']/m0-1))
 assert COLD-1e-4<=row['min_T_K']<=row['max_T_K']<=HOT+1e-4 and row['eos_relative_error']<=1e-6 and row['mass_relative_drift']<=1e-5,row
 return row

def audit():
 rows=audit_continued_steps((OUT/'solver-transcript.txt').read_text(errors='replace'),100,2000,.01);prefix=json.loads((P/'time-step-checks.json').read_text());merged=prefix+rows
 assert [r['step'] for r in merged]==list(range(1,2001)) and all(r['passed'] for r in merged)
 check_settings(json.loads((OUT/'loaded-settings.json').read_text()));frames=json.loads((OUT/'frames.json').read_text());assert [r['step'] for r in frames]==[100,101]+list(range(150,2001,50))
 m0=saved_stats(P/'heated-0000.dat.h5')['mass_per_depth_kg_m'];results=[]
 for frame in frames:
  results.append(guard(OUT/frame['file'],frame['step'],m0));q=fluent(OUT/'continuation-final.cas.h5',OUT/frame['file'],96);assert all(np.isfinite(v).all() for v in q.values()) and all((q[k]>0).all() for k in ['p','rho','T'])
 a=fluent(P/'qualification-final.cas.h5',P/'qualification-final.dat.h5',96);b=fluent(OUT/'continuation-final.cas.h5',OUT/'heated-0100.dat.h5',96);difference={k:float(abs(a[k]-b[k]).max()) for k in a}
 assert all(np.allclose(a[k],b[k],rtol=1e-14,atol=1e-14) for k in a)
 walls,_=fluent_wall_flux(OUT/'continuation-final.cas.h5',OUT/'continuation-final.dat.h5',OUT/'wall-heat-2000.txt',WALLS);assert all(v['face_count']==96 and abs(v['wall_length_m']-L)<1e-10 for v in walls.values())
 assert sum(v['absolute_mass_flux_per_depth_kg_ms'] for v in walls.values())<=1e-12
 files=[p for p in OUT.iterdir() if p.is_file() and (p.suffix=='.h5' or p.name in ['restart-copy-provenance.json','executed-driver.py','loaded-settings.json','settings-comparison.json','declared-settings.json','frames.json','all-step-checks.json','solver-transcript.txt'] or p.name.startswith(('wall-heat-','wall-mass-')))]
 a=dict(numerical_pass=True,new_steps=1900,merged_steps=2000,continued_saved_states=len(frames),selected_common_states=41,exact_saved_checkpoint_source=str(P),restart_field_max_absolute_differences=difference,settings_match_baseline=True,saved_native_checks=results,native_wall_heat_at20s=walls,Nu_hot_at20s=walls['left']['heat_per_depth_W_m']/(K*(HOT-COLD)),no_independent_full_rerun=True,physical_validation=False,energy_scope='Native residual/EOS/source and instantaneous face heat, not full transient energy closure',source_sha256={str(p):sha(p) for p in files})
 (OUT/'production-audit-v001.json').write_text(json.dumps(a,indent=2)+'\n');(S/'fluent-combined-2000-step-checks-v001.json').write_text(json.dumps(merged,indent=2)+'\n');return a

def main():
 assert not (S/'fluent-production-v002-status.json').exists() and not (OUT/'solver-transcript.txt').exists()
 while available_memory()<4*2**30 or shutil.disk_usage(ROOT).free<10*2**30:mark('waiting_resources');time.sleep(30)
 verify('fluent-production-inputs-v002.json');admission=prelaunch('fluent-production-inputs-v002.json');mark('launching_exact1s_continuation',resources=admission)
 session=None;loaded=False
 try:
  assert available_memory()>=4*2**30 and shutil.disk_usage(ROOT).free>=10*2**30
  session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True)
  mark('licensed_connected',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(OUT/'solver-transcript.txt'));st=session.settings;st.file.read_case_data(file_name=str(OUT/'restart-source.cas.h5'));loaded=True
  assert abs(float(st.solution.run_calculation.transient_controls.flow_time())-1)<1e-9
  write(OUT/'settings-comparison.json',check_settings(state(st)));write(OUT/'loaded-settings.json',state(st));write(OUT/'executed-driver.json',dict(source=str(Path(__file__)),sha256=sha(Path(__file__))))
  configure_residuals(st.solution.monitor.residual,True);run=st.solution.run_calculation;run.parameters.time_step_size=.01;run.parameters.max_iter_per_time_step=200
  shutil.copyfile(__file__,OUT/'executed-driver.py');write(OUT/'declared-settings.json',dict(specification(),model='96x96 side-heated cavity continuation1..20s',heated_steps=2000,new_steps=1900,dt_s=.01,qualification_only=False,prefix_step=100,last_step=2000))
  m0=saved_stats(P/'heated-0000.dat.h5')['mass_per_depth_kg_m'];frames=[];previous=100
  for step in [100,101]+list(range(150,2001,50)):
   if step>previous:run.dual_time_iterate(time_step_count=step-previous,max_iter_per_step=200)
   st.results.report.fluxes.mass_flow(zones=WALLS,write_to_file=True,file_name=str(OUT/f'wall-mass-{step:04d}.txt'));st.results.report.fluxes.heat_transfer(zones=WALLS,write_to_file=True,file_name=str(OUT/f'wall-heat-{step:04d}.txt'))
   p=OUT/f'heated-{step:04d}.dat.h5';st.file.write_data(file_name=str(p));frames.append(guard(p,step,m0));write(OUT/'frames.json',frames)
   if step>100:
    checks=audit_continued_steps((OUT/'solver-transcript.txt').read_text(errors='replace'),100,step,.01);write(OUT/'all-step-checks.json',checks);assert all(r['passed'] for r in checks),'Actual step missed unchanged native limit; state saved before stopping'
   previous=step;mark('native_checkpoint_pass',step=step,time_s=step*.01,saved_states=len(frames))
  st.file.write_case_data(file_name=str(OUT/'continuation-final.cas.h5'))
 except Exception:
  if loaded:
   try:st.file.write_case_data(file_name=str(OUT/'recovery.cas.h5'))
   except Exception:pass
  raise
 finally:
  if session:session.exit()
 verify('fluent-production-inputs-v002.json');mark('auditing_closed_native_interval');a=audit();assert a['numerical_pass'];mark('fluent_full_interval_pass',audit=str(OUT/'production-audit-v001.json'),new_steps=1900,combined_steps=2000,OpenFOAM_full_not_yet_run=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise