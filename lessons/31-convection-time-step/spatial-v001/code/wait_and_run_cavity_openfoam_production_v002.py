from pathlib import Path
from datetime import datetime
import os,json,hashlib,subprocess,shutil,time,traceback
from cavity_prelaunch_v001 import prelaunch,verify
from wing_fluent_v001 import available_memory
from audit_cavity_openfoam_production_v001 import audit
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-convection-spatial-v001';O=N/'openfoam-production-v001';F=N/'fluent-production-v002'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(stage,**v):
 p=S/'openfoam-production-v002-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def parent():
 p=F/'production-audit-v001.json';a=json.loads(p.read_text());assert a['numerical_pass'] and a['new_steps']==1900 and a['merged_steps']==2000
 for n,h in a['source_sha256'].items():assert sha(Path(n))==h,n
 return dict(path=str(p),sha256=sha(p))
def main():
 assert not (S/'openfoam-production-v002-status.json').exists()
 while True:
  p=S/'fluent-production-v002-status.json';state=json.loads(p.read_text()) if p.exists() else {}
  if state.get('stage')=='fluent_full_interval_pass':break
  if state.get('stage')=='failed':raise RuntimeError('Fluent failed; sequential OpenFOAM continuation withheld')
  mark('waiting_for_accepted_full_Fluent',observed_Fluent_stage=state.get('stage'),observed_Fluent_step=state.get('step'));time.sleep(30)
 while available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources_after_Fluent');time.sleep(30)
 bound=parent();admission=prelaunch('openfoam-production-inputs-v002.json');parent();verify('openfoam-production-inputs-v002.json')
 assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30
 mark('running_exact1_to20s_continuation',resources=admission,accepted_Fluent=bound)
 with (S/'openfoam-production-native-v002.log').open('xb') as f:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/convection-spatial-v001/run_cavity_openfoam_production_v001.sh'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode
 args=[str(R/'environments/fluent/Scripts/python.exe'),'-B',str(R/'local/natural-convection-production-v001/audit_openfoam_steps_v001.py'),str(O),'16','1900',str(O/'all-step-checks.json')];p=subprocess.run(args,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW);assert p.returncode==0,p.stdout[-500:]
 verify('openfoam-production-inputs-v002.json');parent();mark('auditing_closed_native_interval');a=audit(O);assert a['numerical_pass'];mark('openfoam_full_interval_pass',audit=str(O/'production-audit-v001.json'),new_steps=1900,combined_steps=2000,paired_media_not_yet_exported=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise