from pathlib import Path
import hashlib,json
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');a=json.loads((S/'cylinder-coupled-window-inputs-v001.json').read_text())
for n in ['run_cylinder_guarded_window_v001.py','audit_cylinder_guarded_window_v001.py','run_cylinder_pressure_pair_v003.py','audit_cylinder_window_v002.py','cylinder_fluent_pressure_v002.py','cylinder-coupled-window-inputs-v001.json','cylinder_guard_v001.py']:
 p=S/n;a[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
(S/'guarded-window-controller-inputs-v001.json').write_text(json.dumps(a,indent=2)+'\n');print('Frozen controller dependencies: '+str(len(a)))