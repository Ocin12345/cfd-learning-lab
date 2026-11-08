from pathlib import Path
from datetime import datetime
import json,hashlib,subprocess,os,time,traceback
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');PY=R/'environments/fluent/Scripts/python.exe';FIG=R/'exports/figures/2026-11-08-sloshing-dt-full-v001';VIDEO=R/'exports/videos/2026-11-08-sloshing-time-step-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
 for p,h in json.loads((S/'media-inputs-v001.json').read_text()).items():assert sha(Path(p))==h,p

def mark(stage,**v):
 p=S/'full-visual-pipeline-v001-status.json';q=p.with_suffix('.tmp');q.write_text(json.dumps(dict(stage=stage,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**v),indent=2)+'\n');q.replace(p)

def run(script,args,label):
 with (S/(label+'-v001.log')).open('xb') as f:p=subprocess.run([str(PY),'-B',str(S/script),*map(str,args)],stdout=f,stderr=subprocess.STDOUT,creationflags=subprocess.CREATE_NO_WINDOW)
 assert p.returncode==0,(label,p.returncode)

def main():
 assert not (S/'full-visual-pipeline-v001-status.json').exists();verify()
 while True:
  statuses={code:json.loads((S/(code+'-production-v001-status.json')).read_text()) for code in ['fluent','openfoam']}
  if any(a.get('stage')=='failed' for a in statuses.values()):raise RuntimeError('A full solver controller failed; comparison export withheld')
  if all(statuses[k].get('stage')==k+'_full_interval_pass' for k in statuses):break
  mark('waiting_for_both_accepted_full_audits',solver_stages={k:v.get('stage') for k,v in statuses.items()});time.sleep(30)
 verify()
 from sloshing_wait_admission_v001 import wait_admission
 wait_admission('full-interval-inputs-v002.json',lambda **v:mark('waiting_media_resources',**v));verify()
 for code in ['fluent','openfoam']:
  p=R/'native/runs/showcase-sloshing-time-step-v001'/(code+'-production-v001')/'production-audit-v001.json';a=json.loads(p.read_text());assert a['numerical_pass'] and a['merged_steps']==4800
  for n,h in a['source_sha256'].items():assert sha(Path(n))==h,n
 mark('rendering_actual_four_case_comparison');run('render_sloshing_full_comparison_v001.py',[],'full-native-render');verify()
 mark('encoding_and_verifying_MP4');run('export_sloshing_dt_video_v001.py',['--frames-dir',FIG/'frames','--out',VIDEO],'full-video-export');verify()
 mark('encoding_and_verifying_GIF');run('export_sloshing_dt_gif_v001.py',['--frames-dir',FIG/'frames','--out',VIDEO/'gif-v001'],'full-gif-export');verify()
 mark('visuals_ready_for_actual_inspection',figures=str(FIG),video=str(VIDEO/'sloshing-time-step-comparison.mp4'),gif=str(VIDEO/'gif-v001/sloshing-time-step-comparison.gif'),whole_improvement_complete=False,fresh_review_not_yet_started=True)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise
