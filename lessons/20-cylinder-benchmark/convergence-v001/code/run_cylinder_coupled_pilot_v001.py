from pathlib import Path
import datetime,json,os,shutil,subprocess,sys,time,traceback
from cylinder_guard_v001 import verify_inputs
from audit_cylinder_window_v002 import audit_window
from cylinder_fluent_coupled_v001 import available_memory
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-cylinder-convergence-v001';OLD=R/'native/runs/2026-11-01-vortex-shedding-v001';FROZEN=S/'coupled-controller-inputs-v001.json'

def mark(state,**kw):
 p=S/'cylinder-coupled-pilot-v001-status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),**kw),indent=2)+'\n');t.replace(p)
def idle():
 command="Get-CimInstance Win32_Process | Where-Object { $_.Name -in @('fluent.exe','cortex.exe') } | Select-Object ProcessId,Name,CreationDate,CommandLine | ConvertTo-Json -Compress"
 x=subprocess.run(['powershell','-NoProfile','-Command',command],capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW);assert x.returncode==0
 if x.stdout.strip():raise RuntimeError('An existing Fluent/cortex process is present; no duplicate launched')
 x=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','pgrep','-af','[f]oamRun|[s]impleFoam|[p]impleFoam|[m]pirun'],capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if x.returncode not in [0,1]:raise RuntimeError('Cannot establish OpenFOAM process state: '+repr(x.stderr))
 if x.returncode==0:raise RuntimeError('An existing OpenFOAM/MPI process is present; no concurrent solver launched')
def main():
 assert not (S/'cylinder-coupled-pilot-v001-status.json').exists();results={}
 for label,termination in [('coupled',None)]:
  out=N/'fluent-coupled-pilot-v001';verify_inputs(FROZEN,out)
  while available_memory()<4*1024**3 or shutil.disk_usage('E:/').free<10*1024**3:
   mark('waiting_resources',variant=label);time.sleep(30)
  idle();verify_inputs(FROZEN,out)
  cmd=[sys.executable,'-B',str(S/'cylinder_fluent_coupled_v001.py'),'--mesh',str(OLD/'mesh/vortex.msh'),'--out',str(out),'--resume-case',str(OLD/'fluent-corrected-v003/vortex-final.cas.h5'),'--end-time','4.05']
  if termination is not None:cmd+=['--pressure-termination',str(termination)]
  mark('running_diagnostic',variant=label,output=str(out),controller_pid=os.getpid())
  with (S/('cylinder-'+label+'-stdout.log')).open('w') as stdout,(S/('cylinder-'+label+'-stderr.log')).open('w') as stderr:
   x=subprocess.run(cmd,stdout=stdout,stderr=stderr,cwd=S,creationflags=subprocess.CREATE_NO_WINDOW)
  if x.returncode:raise RuntimeError('Native diagnostic failed: '+label+'; see preserved stderr')
  verify_inputs(FROZEN);a=audit_window(out);results[label]=dict(audit=str(out/'window-audit.json'),passed=a['numerical_window_pass'],miss_steps=a['residual_miss_steps'],pressure_termination=a['pressure_termination'])
  mark('variant_audited',variant=label,**results[label])
 baseline=N/'fluent-pressure-baseline-v001';selected=N/'fluent-coupled-pilot-v001'
 b=json.loads((baseline/'loaded-settings.json').read_text());c=json.loads((selected/'loaded-settings.json').read_text())
 for key in ['general','viscous','material','cell_zones','inlet','outlet','walls','references']:
  assert b[key]==c[key],'Unexpected physical input change: '+key
 for key in ['spatial_discretization','transient_formulation']:
  assert b['methods'][key]==c['methods'][key]
 assert c['methods']['p_v_coupling']['flow_scheme']=='Coupled'
 assert b['residual']==c['residual'],'Residual criteria/normalization changed'
 controls=json.loads((selected/'pressure-controls.json').read_text());assert controls['flow_courant_number']==20 and controls['explicit_pressure_under_relaxation']==controls['explicit_momentum_under_relaxation']==.5
 (S/'cylinder-coupled-comparison-v001.json').write_text(json.dumps(dict(results=results,physical_inputs_spatial_time_schemes_and_residual_settings_unchanged=True,algorithm_bundle=controls,full_extension_launched=False),indent=2)+'\n')
 mark('coupled_pilot_complete',results=results,full_extension_launched=False)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise