from pathlib import Path
from datetime import datetime
import os,json,hashlib,subprocess,shutil,time,traceback
from cavity_prelaunch_v001 import verify,prelaunch
from wing_fluent_v001 import available_memory
from audit_cavity_openfoam_qualification_v001 import audit
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-convection-spatial-v001';O=N/'openfoam-qualification-v001'
def mark(stage,**v):
 p=S/'openfoam-heated-qualification-v001-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def parent():
 for p,key in [(N/'fluent-qualification-v004/qualification-audit-v002.json','numerical_pass'),(N/'openfoam-control-v004/control-audit-v001.json','control_pass')]:
  a=json.loads(p.read_text());assert a[key]
  for n,h in a['source_sha256'].items():assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==h,n
def main():
 assert not (S/'openfoam-heated-qualification-v001-status.json').exists()
 while available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources');time.sleep(30)
 parent();admission=prelaunch('openfoam-heated-qualification-inputs-v001.json');mark('running_finite_100_heated_steps',resources=admission)
 with (S/'openfoam-heated-native-v001.log').open('xb') as f:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/convection-spatial-v001/run_cavity_openfoam_heated_v001.sh'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode
 args=[str(R/'environments/fluent/Scripts/python.exe'),'-B',str(R/'local/natural-convection-production-v001/audit_openfoam_steps_v001.py'),str(O),'16','100',str(O/'all-step-checks.json')];p=subprocess.run(args,capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW);assert p.returncode==0,p.stdout[-500:]
 verify('openfoam-heated-qualification-inputs-v001.json');parent();mark('auditing_all_native_states');a=audit(O);assert a['numerical_pass']
 mark('openfoam_1s_qualification_pass',audit=str(O/'qualification-audit-v001.json'),full20s_not_admitted=True,paired_visual_inspection_pending=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise