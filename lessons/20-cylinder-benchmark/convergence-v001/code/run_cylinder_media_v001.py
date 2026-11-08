from pathlib import Path
import datetime,hashlib,json,os,subprocess,sys,traceback
S=Path(__file__).resolve().parent;R=Path('E:/CFD-Learning-Lab');FIG=R/'exports/figures/2026-11-08-cylinder-developed-paired-v001';VID=R/'exports/videos/2026-11-08-cylinder-developed-v001'
def verify():
 a=json.loads((S/'cylinder-media-inputs-v001.json').read_text())
 for n,h in a.items():assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==h,n
 return len(a)
def mark(state,**kw):
 p=S/'cylinder-media-v001-status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),controller_pid=os.getpid(),**kw),indent=2)+'\n');t.replace(p)
def run(stage,script,args):
 verify();mark(stage)
 with (S/(stage+'-stdout.log')).open('x') as stdout,(S/(stage+'-stderr.log')).open('x') as stderr:r=subprocess.run([sys.executable,'-B',str(S/script),*map(str,args)],cwd=S,stdout=stdout,stderr=stderr,creationflags=subprocess.CREATE_NO_WINDOW)
 assert r.returncode==0,(stage,r.returncode)
def main():
 assert not (S/'cylinder-media-v001-status.json').exists() and not FIG.exists() and not VID.exists();verify()
 run('rendering_full_cylinder','render_cylinder_developed_v001.py',['--out',FIG])
 run('encoding_cylinder_mp4','export_cylinder_video_v001.py',['--frames-dir',FIG,'--out',VID])
 run('encoding_cylinder_gif','export_cylinder_gif_v001.py',['--frames-dir',FIG,'--out',VID/'gif-preview-v001'])
 mark('media_generated_awaiting_visual_record_and_technical_review',figures=str(FIG),videos=str(VID),frames=161,native_states_interpolated=False,physical_validation=False)
if __name__=='__main__':
 try:main()
 except Exception:mark('failed',error=traceback.format_exc());raise