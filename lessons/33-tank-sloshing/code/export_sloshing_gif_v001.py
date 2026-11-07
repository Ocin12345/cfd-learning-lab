"""Encode a fixed-palette GIF preview from all 121 actual native-state PNGs."""
from pathlib import Path
import hashlib,json,shutil,subprocess
import numpy as np
from PIL import Image
R=Path('E:/CFD-Learning-Lab');FIG=R/'exports/figures/2026-11-07-sloshing-1p2s-paired-v003'
OUT=R/'exports/videos/2026-11-07-sloshing-1p2s-v003/gif-preview-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert not OUT.exists()
    audit=FIG/'native-field-audit.json';a=json.loads(audit.read_text());assert a['all_saved_field_checks_pass'] and a['rendering']['frames']==121
    frames=[FIG/f'frame-{i:03d}.png' for i in range(121)];assert all(p.is_file() for p in frames)
    install=json.loads((R/'tools/video-export-v001/installation-receipt.json').read_text());exe=Path(install['executable']);assert sha(exe)==install['executable_sha256']
    OUT.mkdir(parents=True);palette=OUT/'palette.png';gif=OUT/'sloshing-development.gif'
    first=[str(exe),'-hide_banner','-nostdin','-n','-threads','2','-filter_threads','1','-framerate','12','-start_number','0','-i',str(FIG/'frame-%03d.png'),'-vf','scale=1080:1080:flags=lanczos,palettegen=stats_mode=full','-frames:v','1','-update','1',str(palette)]
    second=[str(exe),'-hide_banner','-nostdin','-n','-threads','2','-filter_complex_threads','1','-framerate','12','-start_number','0','-i',str(FIG/'frame-%03d.png'),'-i',str(palette),'-filter_complex','[0:v]scale=1080:1080:flags=lanczos[s];[s][1:v]paletteuse=dither=sierra2_4a','-frames:v','121','-loop','0',str(gif)]
    for name,command in [('palette.log',first),('gif.log',second)]:
        with (OUT/name).open('wb') as log:subprocess.run(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,check=True,creationflags=subprocess.CREATE_NO_WINDOW)
    errors=[];durations=[]
    with Image.open(gif) as im:
        assert im.n_frames==121 and im.size==(1080,1080)
        for index in range(121):
            im.seek(index);decoded=np.asarray(im.convert('RGB'));durations.append(im.info['duration'])
            with Image.open(frames[index]) as original:ref=np.asarray(original.convert('RGB').resize((1080,1080),Image.Resampling.LANCZOS))
            errors.append(float(np.sqrt(np.mean((decoded.astype(np.float32)-ref.astype(np.float32))**2))))
            if index in [0,60,120]:Image.fromarray(decoded).save(OUT/f'decoded-{index:03d}.png')
    assert max(errors)<8 and abs(sum(durations)/1000-121/12)<.04
    record=dict(frames=121,decoded_frames=121,width=1080,height=1080,duration_s=sum(durations)/1000,frame_durations_ms=durations,maximum_rgb_rmse_0_255=max(errors),source_native_field_audit_sha256=sha(audit),source_png_sha256={str(p):sha(p) for p in frames},gif_sha256=sha(gif),gif_bytes=gif.stat().st_size,interpolated_motion=False,fixed_global_palette=True,physical_times_s=a['plot']['times_s'],source_sha256=sha(Path(__file__)))
    (OUT/'gif-export-record.json').write_text(json.dumps(record,indent=2)+'\n');shutil.copyfile(__file__,OUT/'executed-gif-exporter.py')
    print(json.dumps({k:record[k] for k in ['frames','decoded_frames','width','height','duration_s','maximum_rgb_rmse_0_255','gif_bytes']}))
if __name__=='__main__':main()
