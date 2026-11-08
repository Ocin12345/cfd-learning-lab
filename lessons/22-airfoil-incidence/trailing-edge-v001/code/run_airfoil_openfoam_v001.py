from pathlib import Path
import argparse,datetime,hashlib,json,os,shutil,subprocess,time,traceback
import numpy,h5py
from airfoil_prelaunch_v003 import prelaunch,verify_inputs
from wing_fluent_v001 import available_memory
from audit_refined_openfoam_v001 import audit
import plot_openfoam_airfoil_fields_v002 as top
S=Path(__file__).resolve().parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-airfoil-trailing-edge-v001'
def mark(level,state,**values):
 p=S/f'openfoam-level{level}-v001-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),pid=os.getpid(),**values),indent=2)+'\n');q.replace(p)
def verify():
 a=json.loads((S/'openfoam-flow-inputs-v001.json').read_text())
 for p,h in a.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
 verify_inputs()
 for level in (1,2):
  r=json.loads((S/f'analysis-fluent-level{level}-v001/checks.json').read_text());assert r['numerical_checks_passed'] and not r['failed_checks']
  for p,h in r['source_sha256'].items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
 return len(a)
def gate(level):
 while available_memory()<4*2**30 or shutil.disk_usage('E:/').free<10*2**30:
  mark(level,'waiting_resources');time.sleep(30)
 verify();return prelaunch()
def wsl(level,mode):
 command=['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/airfoil-trailing-edge-v001/run_airfoil_openfoam_v001.sh',str(level),mode]
 with (S/f'openfoam-level{level}-{mode}-v001.log').open('xb') as f:r=subprocess.run(command,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 if r.returncode:raise RuntimeError('Native OpenFOAM stage failed: '+mode)
def main(level):
 assert not (S/f'openfoam-level{level}-v001-status.json').exists()
 mark(level,'admission');resources=gate(level)
 case=N/f'openfoam-level{level}-v001';assert json.loads((case/'run-status.json').read_text())['stage']=='prepared_only'
 mark(level,'converting_mesh',resources=resources);wsl(level,'convert')
 mark(level,'checking_exact_native_geometry');geo=top.load_case_geometry(case)
 files=list(geo['inputs'])+[p for part in ('processor0','processor1') for p in (case/part/'constant/polyMesh').glob('*') if p.is_file()]
 native={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
 resources=gate(level)
 for p,h in native.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
 assert available_memory()>=4*2**30 and shutil.disk_usage('E:/').free>=10*2**30
 (S/f'openfoam-level{level}-launch-admission-v001.json').write_text(json.dumps(dict(resources=resources,source_input_count=verify(),converted_native_mesh=geo['metrics'],native_mesh_sha256=native),indent=2)+'\n')
 mark(level,'solving',resources=resources);wsl(level,'solve');verify()
 for p,h in native.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
 mark(level,'auditing');a=audit(case,S/f'analysis-openfoam-level{level}-v001')
 mark(level,'numerical_pass' if a['numerical_checks_passed'] else 'numerical_miss',audit=str(S/f'analysis-openfoam-level{level}-v001/checks.json'),failed_checks=a['failed_checks'])
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--level',type=int,choices=(1,2),required=True);a=p.parse_args()
 try:main(a.level)
 except Exception:mark(a.level,'failed',error=traceback.format_exc());raise
