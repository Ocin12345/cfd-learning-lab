import subprocess,traceback
from wing_repair_common_v001 import *
from wing_repair_prelaunch_v003 import wait_admission,verify
from wing_repair_native_v003 import audit
STATUS=W/'low-relax25-qualification-v001-status.json'
def main():
 assert not STATUS.exists();mark(STATUS,'waiting_low_relax25_admission');wait_admission(lambda **v:mark(STATUS,'waiting_low_relax25_admission',**v))
 with (W/'low-relax25-native-v001.log').open('xb') as f:p=subprocess.run([str(PY),'-B',str(W/'run_wing_relax_qualification_v001.py'),'low-relax25'],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();a=audit('low-relax25');mark(STATUS,'low_relax25_qualification_audited',numerical_pass=a['numerical_pass'],audit=str(N/'fluent-low-relax25-v001/qualification-audit-v001.json'),actual_visual_inspection_pending=True,confirmation_not_launched=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
