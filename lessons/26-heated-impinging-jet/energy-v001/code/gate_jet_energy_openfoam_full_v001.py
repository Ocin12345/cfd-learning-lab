import sys,json,time
sys.path.insert(0,'E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts')
from jet_energy_common_v001 import *
from jet_energy_full_prelaunch_v001 import prelaunch
from jet_energy_parent_gate_v001 import wait_checked

def check_parent():
 selected=json.loads((S/'accepted-fluent-full-parent-v001.json').read_text());p=N/'fluent-production-v001/production-audit-v001.json'
 assert selected['audit']==str(p) and sha(p)==selected['audit_sha256']
 a=json.loads(p.read_text());assert a['numerical_pass'] and a['actual_steps']==600 and a['unique_states']==601
 for f,h in a['source_sha256'].items():
  if sha(f)!=h:raise RuntimeError('Accepted Fluent parent changed: '+f)
 return True
if __name__=='__main__':print(json.dumps(wait_checked(check_parent,prelaunch,time.sleep)))
