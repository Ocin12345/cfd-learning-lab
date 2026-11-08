"""Admission for the corrected smaller-step pair; old runs remain immutable."""
from pathlib import Path
import json,sys
from dam_break_admission_v002 import R,P,sha,project_path,wait_resources,verify as verify_original
def verify(stage):
 proof=verify_original(stage)
 receipt=P/'solver-inputs-v004.json';base=json.loads(receipt.read_text())
 for name,h in base['source_sha256'].items():assert sha(project_path(name))==h,'Corrected source/input changed: '+name
 proof[str(receipt)]=sha(receipt)
 for start,steps in [(1900,200),(3200,400)]:
  path=R/'native/runs/dam-break-over-obstacle-v001/diagnostics'/f'coupled-half-dt-window{start:04d}-v001/window-audit.json'
  a=json.loads(path.read_text());assert a['window_checks_pass'] and a['steps']==steps and a['dt_s']==.00005 and a['native_clock_checks_pass'] and a['maximum_relative_mixture_mass_closure']<=.005
  for name,h in a['native_source_sha256'].items():assert sha(project_path(name))==h,'Corrected diagnostic changed: '+name
  normalized={str(project_path(n)):h for n,h in a['native_source_sha256'].items()}
  assert normalized[str((path.parent/'executed-driver.py').resolve())]==sha(P/'run_coupled_half_dt_windows_v001.py')
  proof[str(path)]=sha(path)
 needed=[] if stage.endswith('-flat') else ['fluent-flat','openfoam-flat']
 if stage.endswith('-production'):needed+=['fluent-collapse','openfoam-collapse']
 for key in needed:
  path=P/(key+'-dt005ms-audit.json');a=json.loads(path.read_text())
  assert a['numerical_checks_pass'] and a['steps']==400 and a['dt_s']==.00005 and a['saved_states']==21 and a['mode']==key.split('-',1)[1]
  for name,h in a['native_source_sha256'].items():assert sha(project_path(name))==h,'Corrected qualification changed: '+name
  proof[str(path)]=sha(path)
 return proof
if __name__=='__main__':print(json.dumps(verify(sys.argv[1])))