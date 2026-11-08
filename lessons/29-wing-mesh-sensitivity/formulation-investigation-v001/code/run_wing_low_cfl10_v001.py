import subprocess,traceback
from wing_repair_common_v001 import *
from wing_repair_prelaunch_v002 import wait_admission,verify
from wing_repair_native_v002 import audit
STATUS=W/'low-cfl10-qualification-v001-status.json'
def main():
 assert not STATUS.exists();mark(STATUS,'waiting_low_cfl10_admission');wait_admission(lambda **v:mark(STATUS,'waiting_low_cfl10_admission',**v))
 with (W/'low-cfl10-native-v001.log').open('xb') as f:p=subprocess.run([str(PY),'-B',str(W/'run_wing_cfl_qualification_v002.py'),'low-cfl10'],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();a=audit('low-cfl10');mark(STATUS,'low_cfl10_qualification_audited',numerical_pass=a['numerical_pass'],audit=str(N/'fluent-low-cfl10-v001/qualification-audit-v001.json'),actual_visual_inspection_pending=True,confirmation_not_launched=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
