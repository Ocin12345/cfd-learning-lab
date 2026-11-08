import subprocess,traceback,json
from jet_energy_common_v001 import *
from jet_energy_prelaunch_v003 import wait_admission,verify
from audit_jet_energy_qualification_v002 import audit_openfoam
STATUS=S/'openfoam-qualification-v003-status.json'
def main():
 assert not STATUS.exists();f=N/'fluent-qualification-v001';a=json.loads((f/'qualification-audit-v002.json').read_text());assert a['numerical_pass'] and json.loads((f/'status.json').read_text())['stage']=='closed'
 for p,h in a['source_sha256'].items():assert sha(p)==h,p
 mark(STATUS,'waiting_OpenFOAM_admission_Fluent_already_accepted');wait_admission(lambda **v:mark(STATUS,'waiting_OpenFOAM_admission',**v))
 for p,h in a['source_sha256'].items():assert sha(p)==h,p
 mark(STATUS,'OpenFOAM_20step_qualification_starting')
 with (S/'openfoam-qualification-native-v003.log').open('xb') as f:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/heated-jet-energy-v001/run_jet_energy_openfoam_qualification_v003.sh'],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();a=audit_openfoam();assert a['numerical_pass'];mark(STATUS,'paired_20step_qualifications_pass',Fluent_reused=True,actual_visual_inspection_pending=True,full_not_launched=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
