from pathlib import Path
from datetime import datetime
import os,json,shutil,time,traceback
import run_cavity_fluent_qualification_v004 as driver
from wing_fluent_v001 import available_memory
from cavity_prelaunch_v001 import verify
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab')
def mark(stage,**v):
 p=S/'fluent-qualification-v004-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)
def main():
 assert not (S/'fluent-qualification-v004-status.json').exists()
 while available_memory()<4*2**30 or shutil.disk_usage(R).free<10*2**30:mark('waiting_resources');time.sleep(30)
 verify('fluent-qualification-inputs-v004.json');mark('running_saved_1200_plus_600_control_and_1s_qualification')
 driver.main();verify('fluent-qualification-inputs-v004.json')
 a=json.loads((driver.OUT/'qualification-audit-v002.json').read_text());assert a['numerical_pass']
 mark('fluent_1s_qualification_pass',audit=str(driver.OUT/'qualification-audit-v002.json'),full_20s_not_admitted=True,openfoam_not_queued=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise