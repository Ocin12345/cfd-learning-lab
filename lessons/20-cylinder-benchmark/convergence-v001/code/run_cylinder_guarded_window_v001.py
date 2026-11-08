from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,sys,time,traceback
from cylinder_guard_v001 import verify_inputs
from cylinder_fluent_coupled_guarded_v001 import available_memory
from run_cylinder_pressure_pair_v003 import idle
from audit_cylinder_guarded_window_v001 import audit_window
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-cylinder-convergence-v001';FROZEN=S/'guarded-window-controller-inputs-v001.json'
def mark(state,**kw):
 p=S/'cylinder-guarded-window-v001-status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),controller_pid=os.getpid(),**kw),indent=2)+'\n');t.replace(p)
def main():
 assert not (S/'cylinder-guarded-window-v001-status.json').exists();parent=N/'fluent-coupled-pilot-v001';out=N/'fluent-coupled-window-v001';verify_inputs(FROZEN,out)
 a=json.loads((parent/'window-audit.json').read_text());assert a['numerical_window_pass'] and a['steps']==10
 while available_memory()<4*1024**3 or shutil.disk_usage('E:/').free<10*1024**3:
  mark('waiting_resources');time.sleep(30)
 idle();verify_inputs(FROZEN,out)
 cmd=[sys.executable,'-B',str(S/'cylinder_fluent_coupled_guarded_v001.py'),'--mesh',str(R/'native/runs/2026-11-01-vortex-shedding-v001/mesh/vortex.msh'),'--out',str(out),'--resume-case',str(parent/'vortex-final.cas.h5'),'--end-time','4.55']
 mark('running_guarded_window',output=str(out),start_time_s=4.05,end_time_s=4.55,steps=100)
 with (S/'cylinder-window-stdout.log').open('x') as stdout,(S/'cylinder-window-stderr.log').open('x') as stderr:
  result=subprocess.run(cmd,cwd=S,stdout=stdout,stderr=stderr,creationflags=subprocess.CREATE_NO_WINDOW)
 if result.returncode:raise RuntimeError('Guarded native window stopped; preserve recovery and inspect stderr/native evidence')
 verify_inputs(FROZEN);a=audit_window(out,4.05,100);assert a['numerical_window_pass'] and a['all_runtime_guards_pass']
 b=json.loads((parent/'loaded-settings.json').read_text());c=json.loads((out/'loaded-settings.json').read_text())
 for k in ['general','viscous','material','cell_zones','inlet','outlet','walls','references','methods','residual']:assert b[k]==c[k],'Unexpected input change: '+k
 mark('guarded_half_second_window_pass',audit=str(out/'window-audit.json'),steps=100,saved_states=101,next_extension_launched=False,physical_validation=False)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise