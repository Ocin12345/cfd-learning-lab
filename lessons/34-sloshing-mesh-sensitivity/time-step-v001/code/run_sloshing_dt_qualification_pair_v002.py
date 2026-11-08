from pathlib import Path
import json,hashlib,subprocess,os,time,datetime,traceback,shutil
from wing_fluent_v001 import available_memory
from sloshing_prelaunch_v001 import prelaunch,verify
from audit_sloshing_dt_qualification_v002 import audit
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-sloshing-time-step-v001';PY=R/'environments/fluent/Scripts/python.exe'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def mark(stage,**v):
 p=S/'qualification-pair-v002-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,observed_at=datetime.datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def wait():
 while available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources');time.sleep(30)
def main():
 assert not (S/'qualification-pair-v002-status.json').exists();wait();mark('Fluent_qualification_starting');prelaunch('qualification-inputs-v002.json')
 with (S/'fluent-qualification-native-v002.log').open('xb') as f:p=subprocess.run([str(PY),'-B',str(S/'run_sloshing_dt_fluent_qualification_v002.py')],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify('qualification-inputs-v002.json')
 a=audit('fluent',N/'fluent-qualification-v001',N/'fluent-qualification-v001/qualification-audit-v001.json');assert a['numerical_pass'];mark('Fluent_80step_qualification_pass')
 wait();prelaunch('qualification-inputs-v002.json')
 for n,h in a['source_sha256'].items():assert sha(Path(n))==h,n
 verify('qualification-inputs-v002.json');assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30;mark('OpenFOAM_qualification_starting')
 with (S/'openfoam-qualification-native-v002.log').open('xb') as f:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/sloshing-time-step-v001/run_sloshing_dt_openfoam_qualification_v002.sh'],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify('qualification-inputs-v002.json')
 a=audit('openfoam',N/'openfoam-qualification-v001',N/'openfoam-qualification-v001/qualification-audit-v001.json');assert a['numerical_pass'];mark('paired_80step_qualifications_pass',full_interval_not_launched=True,actual_visual_inspection_pending=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise
