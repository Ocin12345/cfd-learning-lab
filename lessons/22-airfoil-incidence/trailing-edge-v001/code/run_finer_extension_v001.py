from pathlib import Path
import argparse,datetime,hashlib,json,os,shutil,subprocess,time,traceback
import numpy,h5py
from airfoil_prelaunch_v003 import prelaunch,verify_inputs
from wing_fluent_v001 import available_memory
from audit_refined_openfoam_v001 import audit
import analyse_openfoam_airfoil_v002 as av
S=Path(__file__).resolve().parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-airfoil-trailing-edge-v001'
def mark(level,state,**values):
 p=S/f'openfoam-level{level}-v003-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),pid=os.getpid(),**values),indent=2)+'\n');q.replace(p)
def verify(level):
 a=json.loads((S/f'openfoam-finer-extension-inputs-v001.json').read_text())
 for p,h in a.items():assert av.sha256(Path(p))==h,p
 verify_inputs();case=N/f'openfoam-level{level}-v003';r=json.loads((case/'continuation-inputs-v003.json').read_text());source=Path(r['source_case'])
 for n,h in r['source_manifest'].items():
  assert av.sha256(source/n)==h,n
  if n not in ('system/controlDict','run-status.json'):assert av.sha256(case/n)==h,n
 assert av.sha256(case/'system/controlDict')==r['revised_control_sha256']
 assert av.control_invariants((source/'system/controlDict').read_text())==av.control_invariants((case/'system/controlDict').read_text())
 assert r['source_start_iteration']==24000 and r['end_iteration']==32000
 return r

def main(level,wait_for):
 assert not (S/f'openfoam-level{level}-v003-status.json').exists()
 verify(level)
 if wait_for:
  p=S/wait_for;assert p.parent==S
  while True:
   state=json.loads(p.read_text())['state'] if p.exists() else None
   if state=='failed':raise RuntimeError('Prior sequential job failed structurally; continuation withheld')
   if state in ('numerical_pass','numerical_miss'):break
   mark(level,'waiting_for_prior_solver',wait_for=str(p));time.sleep(30)
 while available_memory()<4*2**30 or shutil.disk_usage('E:/').free<10*2**30:
  mark(level,'waiting_resources');time.sleep(30)
 r=verify(level);resources=prelaunch();r=verify(level)
 assert available_memory()>=4*2**30 and shutil.disk_usage('E:/').free>=10*2**30
 (S/f'finer-extension-admission-v001.json').write_text(json.dumps(dict(resources=resources,exact_saved_checkpoint=True,from_iteration=24000,end_iteration=r['end_iteration'],no_reinitialization=True),indent=2)+'\n')
 mark(level,'solving',endpoint=r['end_iteration'],resources=resources)
 command=['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/airfoil-trailing-edge-v001/run_finer_extension_v001.sh',str(level),str(r['end_iteration'])]
 with (S/f'finer-extension-wrapper-v001.log').open('xb') as f:p=subprocess.run(command,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 if p.returncode:raise RuntimeError('Continuation native solver failed; inspect preserved evidence')
 verify(level);mark(level,'auditing');a=audit(N/f'openfoam-level{level}-v003',S/f'analysis-openfoam-level{level}-v003',r['end_iteration'],['solver-continuation.log','solver-continuation-v003.log'])
 mark(level,'numerical_pass' if a['numerical_checks_passed'] else 'numerical_miss',endpoint=r['end_iteration'],failed_checks=a['failed_checks'],audit=str(S/f'analysis-openfoam-level{level}-v003/checks.json'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--level',type=int,choices=(2,),required=True);p.add_argument('--wait-for');a=p.parse_args()
 try:main(a.level,a.wait_for)
 except Exception:mark(a.level,'failed',error=traceback.format_exc());raise
