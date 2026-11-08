from pathlib import Path
import datetime,hashlib,json,os,shutil,subprocess,time,traceback
from cylinder_fluent_coupled_guarded_v002 import available_memory
from run_cylinder_pressure_pair_v003 import idle
from verify_openfoam_qualification_v001 import verify
from audit_cylinder_openfoam_v002 import audit_openfoam
S=Path(__file__).resolve().parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-cylinder-convergence-v001'
def mark(state,**kw):
 p=S/'scheduled-openfoam-qualification-v001-status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),controller_pid=os.getpid(),**kw),indent=2)+'\n');t.replace(p)
def fluent_acceptance():
 a=json.loads((N/'fluent-developed-v001/window-audit.json').read_text());assert a['numerical_window_pass'] and a['all_runtime_guards_pass'] and a['steps']==690 and a['saved_states']==691 and abs(a['end_time_s']-8)<1e-8
 for n,h in a['native_source_sha256'].items():assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==h,n
 return dict(audit=str(N/'fluent-developed-v001/window-audit.json'),sha256=hashlib.sha256((N/'fluent-developed-v001/window-audit.json').read_bytes()).hexdigest(),verified_native_files=len(a['native_source_sha256']))
def main():
 assert not (S/'scheduled-openfoam-qualification-v001-status.json').exists();verify();mark('waiting_for_accepted_fluent_interval')
 while True:
  p=S/'cylinder-developed-fluent-v001-status.json'
  if p.exists():
   state=json.loads(p.read_text())['state']
   if state=='failed':raise RuntimeError('Fluent interval failed; OpenFOAM not launched')
   if state=='fluent_development_interval_pass':break
  time.sleep(30)
 mark('verifying_completed_fluent_native_sources');accepted=fluent_acceptance()
 while available_memory()<4*1024**3 or shutil.disk_usage('E:/').free<10*1024**3:
  mark('waiting_resources_after_fluent_acceptance');time.sleep(30)
 idle();verify()
 if available_memory()<4*1024**3 or shutil.disk_usage('E:/').free<10*1024**3:raise RuntimeError('Immediately-prelaunch resource gate failed')
 receipt=dict(accepted_fluent=accepted,available_RAM_bytes=available_memory(),free_E_bytes=shutil.disk_usage('E:/').free,inputs_verified=verify(),two_ranks=True,not_concurrent=True)
 (S/'openfoam-qualification-launch-admission-v001.json').write_text(json.dumps(receipt,indent=2)+'\n');mark('running_openfoam_qualification')
 command=['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001/run_cylinder_openfoam_qualification_v001.sh']
 with (S/'openfoam-qualification-wrapper.log').open('x') as f:r=subprocess.run(command,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 if r.returncode:raise RuntimeError('OpenFOAM qualification command failed; see native/wrapper evidence')
 verify();a=audit_openfoam(N/'openfoam-qualification-v001',4.,100)
 mark('openfoam_qualification_pass' if a['numerical_window_pass'] else 'openfoam_qualification_numerical_miss',audit=str(N/'openfoam-qualification-v001/window-audit.json'),steps=100,states=101,full_extension_launched=False)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise