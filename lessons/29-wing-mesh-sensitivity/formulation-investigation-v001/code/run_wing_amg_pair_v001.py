import subprocess,traceback
from wing_repair_common_v001 import *
from wing_repair_prelaunch_v001 import wait_admission,verify
from wing_repair_native_v001 import audit
STATUS=W/'qualification-pair-v001-status.json'
def main():
 assert not STATUS.exists();results={}
 for label in ['baseline','tight-amg']:
  mark(STATUS,'waiting_'+label+'_admission');wait_admission(lambda **v:mark(STATUS,'waiting_'+label+'_admission',**v))
  with (W/(label+'-native-v001.log')).open('xb') as f:p=subprocess.run([str(PY),'-B',str(W/'run_wing_amg_qualification_v001.py'),label],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
  assert p.returncode==0,(label,p.returncode);verify();a=audit(label);results[label]=dict(numerical_pass=a['numerical_pass'],audit=str(N/('fluent-'+label+'-v001')/'qualification-audit-v001.json'));mark(STATUS,label+'_audited',result=results[label])
 selected='baseline' if results['baseline']['numerical_pass'] else 'tight-amg' if results['tight-amg']['numerical_pass'] else None;mark(STATUS,'qualification_pair_audited',results=results,selected=selected,actual_visual_inspection_pending=True,confirmation_not_launched=True)
if __name__=='__main__':
 try:main()
 except Exception:mark(STATUS,'failed',error=traceback.format_exc());raise
