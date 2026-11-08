from pathlib import Path
import datetime,json,os,shutil,subprocess,sys,time,traceback
from cylinder_guard_v001 import verify_inputs
from audit_cylinder_window_v001 import audit_window
from cylinder_fluent_pressure_v002 import available_memory
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-cylinder-convergence-v001';OLD=R/'native/runs/2026-11-01-vortex-shedding-v001';FROZEN=S/'diagnostic-controller-inputs-v001.json'

def mark(state,**kw):
 p=S/'cylinder-diagnostic-pair-v001-status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),**kw),indent=2)+'\n');t.replace(p)
def idle():
 command="Get-CimInstance Win32_Process | Where-Object { $_.Name -in @('fluent.exe','cortex.exe') } | Select-Object ProcessId,Name,CreationDate,CommandLine | ConvertTo-Json -Compress"
 x=subprocess.run(['powershell','-NoProfile','-Command',command],capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW);assert x.returncode==0
 if x.stdout.strip():raise RuntimeError('An existing Fluent/cortex process is present; no duplicate launched')
 x=subprocess.run(['wsl','-d','CFD-Foundation14','--','pgrep','-af','[f]oamRun|[s]impleFoam|[p]impleFoam|[m]pirun'],capture_output=True,text=True,creationflags=subprocess.CREATE_NO_WINDOW)
 if x.returncode not in [0,1]:raise RuntimeError('Cannot establish OpenFOAM process state: '+x.stderr)
 if x.returncode==0:raise RuntimeError('An existing OpenFOAM/MPI process is present; no concurrent solver launched')
def main():
 assert not (S/'cylinder-diagnostic-pair-v001-status.json').exists();results={}
 for label,termination in [('baseline',None),('tight',.001)]:
  out=N/('fluent-pressure-'+label+'-v001');verify_inputs(FROZEN,out)
  while available_memory()<4*1024**3 or shutil.disk_usage('E:/').free<10*1024**3:
   mark('waiting_resources',variant=label);time.sleep(30)
  idle();verify_inputs(FROZEN,out)
  cmd=[sys.executable,'-B',str(S/'cylinder_fluent_pressure_v002.py'),'--mesh',str(OLD/'mesh/vortex.msh'),'--out',str(out),'--resume-case',str(OLD/'fluent-corrected-v003/vortex-final.cas.h5'),'--end-time','4.05']
  if termination is not None:cmd+=['--pressure-termination',str(termination)]
  mark('running_diagnostic',variant=label,output=str(out),controller_pid=os.getpid())
  with (S/('cylinder-'+label+'-stdout.log')).open('w') as stdout,(S/('cylinder-'+label+'-stderr.log')).open('w') as stderr:
   x=subprocess.run(cmd,stdout=stdout,stderr=stderr,cwd=S,creationflags=subprocess.CREATE_NO_WINDOW)
  if x.returncode:raise RuntimeError('Native diagnostic failed: '+label+'; see preserved stderr')
  verify_inputs(FROZEN);a=audit_window(out);results[label]=dict(audit=str(out/'window-audit.json'),passed=a['numerical_window_pass'],miss_steps=a['residual_miss_steps'],pressure_termination=a['pressure_termination'])
  mark('variant_audited',variant=label,**results[label])
 # Control comparison: exact native controls differ only in pressure AMG termination.
 b=json.loads((N/'fluent-pressure-baseline-v001/pressure-controls.json').read_text());t=json.loads((N/'fluent-pressure-tight-v001/pressure-controls.json').read_text());expected=json.loads(json.dumps(b['after']));expected['advanced']['multi_grid']['mg_controls']['pressure']['termination_criteria']=.001;assert expected==t['after']
 (S/'cylinder-pressure-comparison-v001.json').write_text(json.dumps(dict(results=results,only_pressure_termination_changes=True,full_extension_launched=False),indent=2)+'\n')
 mark('diagnostic_pair_complete',results=results,full_extension_launched=False)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise