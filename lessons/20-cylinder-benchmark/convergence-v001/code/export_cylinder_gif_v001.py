"""Encode/verify an aspect-preserving GIF from all actual audited PNG states."""
from pathlib import Path
import argparse,hashlib,json,shutil,subprocess
import numpy as np
from PIL import Image
R=Path('E:/CFD-Learning-Lab')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--frames-dir',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
 FIG=args.frames_dir.resolve();OUT=args.out.resolve()
 assert FIG.is_relative_to((R/'exports/figures').resolve()) and OUT.is_relative_to((R/'exports/videos').resolve()) and not OUT.exists()
 from wing_fluent_v001 import available_memory
 assert available_memory()/2**30>=4 and shutil.disk_usage(R).free>=10*2**30
 audit=FIG/'native-field-audit.json';a=json.loads(audit.read_text());assert a['all_saved_field_checks_pass'] and not a.get('preview',False)
 n=a['rendering']['frames'];fps=1000/a['rendering']['frame_duration_ms'];assert n==161 and abs(fps-20)<1e-10 and a['plot']['times_s'][0]==4 and a['plot']['times_s'][-1]==8
 frames=[FIG/f'frame-{i:03d}.png' for i in range(n)]
 for p in frames:assert a['rendered_png_sha256'][str(p)]==sha(p),p
 with Image.open(frames[0]) as im:width,height=im.size
 out_width=1080;out_height=round(height*out_width/width)
 assert abs(out_height/out_width-height/width)<1e-10
 install=json.loads((R/'tools/video-export-v001/installation-receipt.json').read_text());exe=Path(install['executable']);assert sha(exe)==install['executable_sha256']
 OUT.mkdir(parents=True);palette=OUT/'palette.png';gif=OUT/'cylinder-developed.gif'
 scale=f'scale={out_width}:{out_height}:flags=lanczos'
 first=[str(exe),'-hide_banner','-nostdin','-n','-threads','2','-filter_threads','1','-framerate',str(fps),'-start_number','0','-i',str(FIG/'frame-%03d.png'),'-vf',scale+',palettegen=stats_mode=full','-frames:v','1','-update','1',str(palette)]
 second=[str(exe),'-hide_banner','-nostdin','-n','-threads','2','-filter_complex_threads','1','-framerate',str(fps),'-start_number','0','-i',str(FIG/'frame-%03d.png'),'-i',str(palette),'-filter_complex','[0:v]'+scale+'[s];[s][1:v]paletteuse=dither=sierra2_4a','-frames:v',str(n),'-loop','0',str(gif)]
 for name,command in [('palette.log',first),('gif.log',second)]:
  with (OUT/name).open('wb') as log:subprocess.run(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
 errors=[];durations=[]
 with Image.open(gif) as im:
  assert im.n_frames==n and im.size==(out_width,out_height)
  for index in range(n):
   im.seek(index);decoded=np.asarray(im.convert('RGB'));durations.append(im.info['duration'])
   with Image.open(frames[index]) as original:ref=np.asarray(original.convert('RGB').resize((out_width,out_height),Image.Resampling.LANCZOS))
   errors.append(float(np.sqrt(np.mean((decoded.astype(np.float32)-ref.astype(np.float32))**2))))
   if index in [0,40,80,120,160]:Image.fromarray(decoded).save(OUT/f'decoded-{index:03d}.png')
 assert max(errors)<8 and abs(sum(durations)/1000-n/fps)<.04
 record=dict(frames=n,decoded_frames=n,width=out_width,height=out_height,duration_s=sum(durations)/1000,frame_durations_ms=durations,maximum_rgb_rmse_0_255=max(errors),source_native_field_audit_sha256=sha(audit),source_png_sha256={str(p):sha(p) for p in frames},gif_sha256=sha(gif),gif_bytes=gif.stat().st_size,interpolated_motion=False,fixed_global_palette=True,physical_times_s=a['plot']['times_s'],source_sha256=sha(Path(__file__)))
 (OUT/'gif-export-record.json').write_text(json.dumps(record,indent=2)+'\n');shutil.copyfile(__file__,OUT/'executed-gif-exporter.py');print(json.dumps({k:record[k] for k in ['frames','decoded_frames','width','height','duration_s','maximum_rgb_rmse_0_255','gif_bytes']}))
if __name__=='__main__':main()