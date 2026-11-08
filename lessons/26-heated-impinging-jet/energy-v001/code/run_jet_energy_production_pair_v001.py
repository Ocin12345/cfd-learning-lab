import subprocess,traceback,json
from jet_energy_common_v001 import *
from jet_energy_full_prelaunch_v001 import wait_admission,verify
from audit_jet_energy_full_v001 import audit_fluent,audit_openfoam
STATUS=S/'production-pair-v001-status.json'
def main():
 assert not STATUS.exists();mark(STATUS,'waiting_Fluent_full_admission');wait_admission(lambda **v:mark(STATUS,'waiting_Fluent_full_admission',**v))
 with (S/'fluent-production-native-v001.log').open('xb') as f:p=subprocess.run([str(PY),'-B',str(S/'run_jet_energy_fluent_production_v001.py')],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();a=audit_fluent();assert a['numerical_pass'];mark(STATUS,'Fluent_full_600step_pass');dump(S/'accepted-fluent-full-parent-v001.json',dict(audit=str(N/'fluent-production-v001/production-audit-v001.json'),audit_sha256=sha(N/'fluent-production-v001/production-audit-v001.json')))
 wait_admission(lambda **v:mark(STATUS,'waiting_OpenFOAM_full_admission',**v))
 for p,h in a['source_sha256'].items():assert sha(p)==h,p
 mark(STATUS,'OpenFOAM_580step_continuation_starting')
 with (S/'openfoam-production-native-v001.log').open('xb') as f:p=subprocess.run(['wsl.exe','-d','CFD-Foundation14','-e','bash','/mnt/e/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/heated-jet-energy-v001/run_jet_energy_openfoam_production_v001.sh'],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,p.returncode;verify();a=audit_openfoam();assert a['numerical_pass'];mark(STATUS,'paired_full_600step_audits_pass',actual_media_inspection_pending=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
