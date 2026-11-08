from pathlib import Path
from datetime import datetime
import hashlib,json,os,shutil,tempfile,time,traceback
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).resolve().parent;M=R/'native/meshes/showcase-mach2-wedge-spatial-v001/fine-v001';O=R/'native/runs/showcase-mach2-wedge-spatial-v001/fluent-quality-v001'
for n in ('TEMP','TMP','TMPDIR'):os.environ[n]=str(R/'temp')
tempfile.tempdir=str(R/'temp')
from wedge_prelaunch_v001 import prelaunch,verify
from wing_fluent_v001 import available_memory
import ansys.fluent.core as pf
from nozzle_fluent_v004 import quality_metrics
def mark(stage,**values):
 p=S/'fluent-quality-v001-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**values),indent=2)+'\n');q.replace(p)
def main():
 assert not O.exists() and not (S/'fluent-quality-v001-status.json').exists();verify('mesh-qualification-inputs-v001.json')
 O.mkdir(parents=True);session=None;saved=False
 try:
  while available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources');time.sleep(30)
  receipt=prelaunch('mesh-qualification-inputs-v001.json');(O/'launch-admission.json').write_text(json.dumps(receipt,indent=2)+'\n')
  mark('launching_hidden_two_rank_session');pf.config.examples_path=str(O/'examples')
  session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(O),start_timeout=180,cleanup_on_exit=True)
  mark('live_startup_succeeded',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(O/'solver-transcript.txt'))
  session.settings.file.read_mesh(file_name=str(M/'wedge.msh'));session.settings.mesh.check();session.settings.mesh.quality()
  for attempt in range(30):
   try:q=quality_metrics((O/'solver-transcript.txt').read_text(errors='replace'));break
   except ValueError:
    if attempt==29:raise
    time.sleep(.2)
  (O/'native-mesh-quality.json').write_text(json.dumps(q,indent=2)+'\n')
  accepted=q['minimum_orthogonal_quality']>=.2
  session.settings.file.write_case(file_name=str(O/('quality-accepted.cas.h5' if accepted else 'rejected-mesh.cas.h5')));saved=True
  if not accepted:raise RuntimeError('Native quality below unchanged.2 floor')
  verify('mesh-qualification-inputs-v001.json');mark('quality_only_pass',quality=q,flow_iterations_executed=0,native_case_saved=True)
 finally:
  if session is not None:session.exit()
 mark('closed_quality_only_pass',quality=q,flow_iterations_executed=0,native_case_saved=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise
