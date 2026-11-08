import wing_import_paths_v001
import subprocess,traceback
from wing_repair_common_v001 import *
from wing_residual_window_prelaunch_v001 import wait_admission,verify
from audit_retained_window_v001 import audit_window
STATUS=W/'retained-window-v001-status.json'
def main():
 assert not STATUS.exists();mark(STATUS,'waiting_retained_window_admission');wait_admission(lambda **v:mark(STATUS,'waiting_retained_window_admission',**v))
 with (W/'retained-window-native-v001.log').open('xb') as f:p=subprocess.run([str(PY),'-B',str(W/'run_wing_retained_window_v001.py'),'retained-residual-window'],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();mark(STATUS,'auditing_closed_eight_states');r=audit_window();mark(STATUS,'retained_window_capture_audited',diagnostic_capture_pass=r['diagnostic_capture_pass'],original_qualification_pass=False,descriptive_persistence_all_eight=r['descriptive_persistence_all_eight'],audit=str(N/'fluent-retained-residual-window-v001/retained-window-audit-v001.json'),actual_visuals_pending=True,no_repair_algorithm_launched=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
