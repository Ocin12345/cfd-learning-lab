import wing_import_paths_v001
import subprocess,traceback
from wing_repair_common_v001 import *
from wing_simplec_prelaunch_v001 import wait_admission,verify
from audit_simplec_capability_v001 import audit_capability
STATUS=W/'simplec-capability-v001-status.json'
def main():
 assert not STATUS.exists();mark(STATUS,'waiting_zero_advance_admission');wait_admission(lambda **v:mark(STATUS,'waiting_zero_advance_admission',**v))
 with (W/'simplec-capability-native-v001.log').open('xb') as f:p=subprocess.run([str(PY),'-B',str(W/'run_wing_simplec_capability_v001.py')],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();mark(STATUS,'auditing_closed_zero_advance');r=audit_capability();mark(STATUS,'simplec_zero_advance_capability_pass',new_iterations=0,capability_pass=True,qualification_not_launched=True,audit=str(N/'fluent-simplec-capability-v001/capability-audit-v001.json'),native_schema_review_pending=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
