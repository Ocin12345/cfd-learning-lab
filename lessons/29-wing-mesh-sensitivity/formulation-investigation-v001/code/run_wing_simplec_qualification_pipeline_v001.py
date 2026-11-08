import wing_import_paths_v001
import subprocess,traceback
from wing_repair_common_v001 import *
from wing_simplec_qualification_prelaunch_v001 import wait_admission,verify
from audit_simplec_qualification_v001 import audit
STATUS=W/'simplec-qualification-v001-status.json'
def main():
 assert not STATUS.exists();mark(STATUS,'waiting_finite100_admission');wait_admission(lambda **v:mark(STATUS,'waiting_finite100_admission',**v))
 with (W/'simplec-qualification-native-v001.log').open('xb') as f:p=subprocess.run([str(PY),'-B',str(W/'run_wing_simplec_qualification_v001.py')],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();mark(STATUS,'auditing_closed_finite100');r=audit();mark(STATUS,'simplec_qualification_pass' if r['numerical_pass'] else 'simplec_qualification_numerical_miss',new_iterations=100,endpoint=1600,numerical_pass=r['numerical_pass'],audit=str(N/'fluent-simplec-qualification-v001/qualification-audit-v001.json'),actual_visual_inspection_pending=True,confirmation_not_launched=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
