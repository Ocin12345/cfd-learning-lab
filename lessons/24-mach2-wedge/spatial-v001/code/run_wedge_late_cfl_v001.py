from pathlib import Path
from datetime import datetime
import os,json,traceback,subprocess,hashlib,shutil,time
from wing_fluent_v001 import available_memory
from wedge_prelaunch_v001 import prelaunch,verify
from audit_wedge_late_cfl_v001 import audit
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-mach2-wedge-spatial-v001';O=N/'openfoam-late-cfl-v001'
def mark(stage,**v):
 p=S/'openfoam-late-cfl-v001-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def verify_parent():
 p=N/'openfoam-fine-v001/endpoint-audit.json';a=json.loads(p.read_text());assert a['numerical_endpoint_pass']
 for n,h in a['source_sha256'].items():assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==h,n
 return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 assert not (S/'openfoam-late-cfl-v001-status.json').exists()
 while available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources');time.sleep(30)
 parent=verify_parent();r=prelaunch('late-cfl-inputs-v001.json');verify_parent();verify('late-cfl-inputs-v001.json')
 assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30
 mark('running_exact_7_to_8ms_continuation',parent_audit_sha256=parent,resources=r)
 with (S/'late-cfl-native-v001.log').open('xb') as f:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/mach2-wedge-spatial-v001/run_wedge_late_cfl_v001.sh'],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode
 verify('late-cfl-inputs-v001.json');verify_parent();mark('auditing_native_checkpoint_branch');a=audit(O)
 mark('late_cfl_endpoint_pass' if a['numerical_endpoint_pass'] else 'late_cfl_endpoint_numerical_miss',audit=str(O/'endpoint-audit.json'),failed_checks=a['failed_checks'],independent_full_rerun=False)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise