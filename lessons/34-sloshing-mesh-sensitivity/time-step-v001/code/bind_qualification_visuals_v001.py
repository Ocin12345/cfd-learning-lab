from pathlib import Path
import json,hashlib,datetime
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');O=R/'exports/figures/2026-11-08-sloshing-dt-qualification-v001';N=R/'native/runs/showcase-sloshing-time-step-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=json.loads((O/'qualification-render-audit.json').read_text())
for p,h in a['image_sha256'].items():assert sha(Path(p))==h
for p,h in a['selected_audits_sha256'].items():assert sha(Path(p))==h
receipt={'passed':True,'observed_at':datetime.datetime.now().astimezone().isoformat(),'inspected_actual_images_sha256':a['image_sha256'],'selected_audits_sha256':a['selected_audits_sha256'],'observations':'Personally inspected both actual PNGs. Tilted native water fraction and local mixture speeds shown without smoothing; four endpoint depth profiles and signed fine-column differences visible. All native residual histories below original limits, norms labelled separately. Qualification only0.02s, notfull1.2s or validation. Local near-interface velocity features are numerical fields, not proof of physically resolved interface currents.','full_interval_complete':False}
(S/'paired-qualification-visual-review-v001.json').write_text(json.dumps(receipt,indent=2)+'\n')
