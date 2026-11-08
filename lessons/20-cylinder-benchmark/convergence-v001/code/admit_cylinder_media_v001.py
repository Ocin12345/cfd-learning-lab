from pathlib import Path
import hashlib,json,subprocess,sys
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');R=Path('E:/CFD-Learning-Lab');P=R/'exports/figures/2026-11-08-cylinder-developed-preview-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
results=[]
for kind in ['video','gif']:
 out=R/'exports/videos'/('cylinder-negative-preview-'+kind+'-v001');assert not out.exists()
 p=subprocess.run([sys.executable,'-B',str(S/('export_cylinder_'+kind+'_v001.py')),'--frames-dir',str(P),'--out',str(out)],capture_output=True,creationflags=subprocess.CREATE_NO_WINDOW)
 (S/('preview-rejection-'+kind+'-v001.log')).write_bytes(p.stdout+p.stderr);assert p.returncode!=0 and not out.exists();results.append(dict(export=kind,rejected_preview=True,no_output_created=True))
(S/'preview-export-rejection-v001.json').write_text(json.dumps(results,indent=2)+'\n')
(S/'developed-preview-visual-review-v001.json').write_text(json.dumps(dict(reviewed_by='Codex',visually_inspected_frames=[1,2],physical_times_s=[6,8],shown_in_chat_frame=2,fixed_speed_and_vorticity_scales=True,readable_axes_and_captions=True,raw_velocity_derived_vorticity_saturation_disclosed=True,physical_validation=False,source_sha256={str(P/n):sha(P/n) for n in ['frame-001.png','frame-002.png','native-field-audit.json']}),indent=2)+'\n')
print('Both encoders reject the actual partial preview before creating an output')