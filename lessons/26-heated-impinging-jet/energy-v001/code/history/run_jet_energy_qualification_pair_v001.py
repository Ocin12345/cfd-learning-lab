import subprocess,traceback,time
from jet_energy_common_v001 import *
from jet_energy_prelaunch_v001 import wait_admission,verify
from audit_jet_energy_qualification_v001 import audit_fluent,audit_openfoam
STATUS=S/'qualification-pair-v001-status.json'
def main():
 assert not STATUS.exists();mark(STATUS,'waiting_fluent_admission');wait_admission(lambda **v:mark(STATUS,'waiting_fluent_admission',**v))
 with (S/'fluent-qualification-native-v001.log').open('xb') as f:p=subprocess.run([str(PY),'-B',str(S/'run_jet_energy_fluent_qualification_v001.py')],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();a=audit_fluent();assert a['numerical_pass'];mark(STATUS,'Fluent_20step_qualification_pass')
 wait_admission(lambda **v:mark(STATUS,'waiting_openfoam_admission',**v))
 for p,h in a['source_sha256'].items():assert sha(p)==h,p
 mark(STATUS,'OpenFOAM_qualification_starting')
 with (S/'openfoam-qualification-native-v001.log').open('xb') as f:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/heated-jet-energy-v001/run_jet_energy_openfoam_qualification_v001.sh'],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();a=audit_openfoam();assert a['numerical_pass'];mark(STATUS,'paired_20step_qualifications_pass',actual_visual_inspection_pending=True,full_not_launched=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
