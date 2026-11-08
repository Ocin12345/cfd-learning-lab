from pathlib import Path
import hashlib,json
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');R=Path('E:/CFD-Learning-Lab');paths=[S/n for n in ['render_cylinder_developed_v001.py','export_cylinder_video_v001.py','export_cylinder_gif_v001.py','run_cylinder_media_v001.py','developed-preview-visual-review-v001.json','preview-export-rejection-v001.json']]
for solver in ['fluent','openfoam']:
 paths +=[S/('analysis-'+solver+'-v001')/n for n in ['combined-'+solver+'-audit.json',solver+'-snapshots.json']]
paths +=[R/'repo/cfd-learning-lab/scripts'/n for n in ['plot_vortex_shedding_v002.py','wing_fluent_v001.py']]
(S/'cylinder-media-inputs-v001.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2)+'\n');print('Media inputs frozen')