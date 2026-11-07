"""Encode audited CFD PNG states as MP4; no generated intermediate flow states."""
from pathlib import Path
import argparse,hashlib,json,re,shutil,subprocess
import numpy as np
from PIL import Image
R=Path('E:/CFD-Learning-Lab')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--frames-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True);args=p.parse_args()
 assert args.frames_dir.is_relative_to(R/'exports/figures') and args.out.is_relative_to(R/'exports/videos') and not args.out.exists()
 audit=args.frames_dir/'native-field-audit.json';a=json.loads(audit.read_text());assert a['all_saved_field_checks_pass']
 n=a['rendering']['frames'];ms=a['rendering']['frame_duration_ms'];fps=1000/ms
 frames=[args.frames_dir/f'frame-{i:03d}.png' for i in range(n)];assert all(p.is_file() for p in frames)
 with Image.open(frames[0]) as im:w,h=im.size
 assert w%2==0 and h%2==0
 for frame in frames:
  with Image.open(frame) as im:assert im.size==(w,h)
 assert shutil.disk_usage(R).free>2*2**30
 install=json.loads((R/'tools/video-export-v001/installation-receipt.json').read_text());exe=Path(install['executable']);assert sha(exe)==install['executable_sha256']
 args.out.mkdir(parents=True);movie=args.out/'sloshing-development.mp4'
 command=[str(exe),'-hide_banner','-nostdin','-n','-threads','2','-framerate',str(fps),'-start_number','0','-i',str(args.frames_dir/'frame-%03d.png'),'-frames:v',str(n),'-an','-c:v','libx264','-threads','2','-preset','medium','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',str(movie)]
 with (args.out/'encode.log').open('wb') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
 decoded=0;errors=[];frame_bytes=w*h*3
 with (args.out/'decode.log').open('wb') as log:
  proc=subprocess.Popen([str(exe),'-hide_banner','-nostdin','-threads','2','-i',str(movie),'-an','-f','rawvideo','-pix_fmt','rgb24','-fps_mode','passthrough','pipe:1'],stdout=subprocess.PIPE,stderr=log,stdin=subprocess.DEVNULL,creationflags=subprocess.CREATE_NO_WINDOW)
  try:
   while True:
    buf=proc.stdout.read(frame_bytes)
    if not buf:break
    assert len(buf)==frame_bytes and decoded<n,'Incomplete or unexpected video frame'
    v=np.frombuffer(buf,dtype=np.uint8).reshape(h,w,3)
    with Image.open(frames[decoded]) as im:original=np.asarray(im.convert('RGB'))
    error=float(np.sqrt(np.mean((v.astype(np.float32)-original.astype(np.float32))**2)));errors.append(error)
    if decoded in {0,n//2,n-1}:Image.fromarray(v).save(args.out/f'decoded-{decoded:03d}.png')
    decoded+=1
  finally:
   proc.stdout.close();rc=proc.wait()
 assert rc==0 and decoded==n,(rc,decoded,n)
 assert max(errors)<5,'Unexpectedly large RGB encoding error'
 log=(args.out/'decode.log').read_text(errors='replace');duration=re.search(r'Duration: (\d+):(\d+):([\d.]+)',log);assert duration
 seconds=int(duration[1])*3600+int(duration[2])*60+float(duration[3]);assert abs(seconds-n/fps)<.02
 assert re.search(r'Video: h264',log) and re.search(r'\b'+re.escape(f'{fps:g}')+r' fps',log)
 record=dict(input_audit=str(audit),input_audit_sha256=sha(audit),frames=n,decoded_frames=decoded,width=w,height=h,fps=fps,video_duration_s=seconds,physical_times_s=a['plot']['times_s'],interpolated_motion=False,encoder=install,encode_command=command,frame_rgb_rmse_0_255=errors,maximum_rgb_rmse_0_255=max(errors),source_png_sha256={str(p):sha(p) for p in frames},video_sha256=sha(movie),video_bytes=movie.stat().st_size,physical_validation=False)
 (args.out/'video-export-record.json').write_text(json.dumps(record,indent=2)+'\n');shutil.copyfile(__file__,args.out/'executed-video-exporter.py')
 print(json.dumps({k:record[k] for k in ['frames','decoded_frames','width','height','fps','video_duration_s','maximum_rgb_rmse_0_255','video_bytes']},indent=2))
if __name__=='__main__':main()
