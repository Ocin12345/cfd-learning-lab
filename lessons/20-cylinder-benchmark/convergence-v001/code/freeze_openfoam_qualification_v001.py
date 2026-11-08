from pathlib import Path
import hashlib,json
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');R=Path('E:/CFD-Learning-Lab');new=R/'native/runs/showcase-cylinder-convergence-v001/openfoam-qualification-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=json.loads((S/'developed-controller-inputs-v001.json').read_text())
for n in ['verify_openfoam_qualification_v001.py','audit_cylinder_openfoam_v002.py','wait_and_run_openfoam_qualification_v001.py','run_cylinder_openfoam_qualification_v001.sh','openfoam-qualification-copy-proof-v001.json','test_openfoam_coverage_v001.py','openfoam-coverage-green-v001.log']:
 p=S/n;compile(p.read_text(),p.name,'exec') if p.suffix=='.py' else None;a[str(p)]=sha(p)
for p in new.rglob('*'):
 if p.is_file():a[str(p)]=sha(p)
copy=json.loads((S/'openfoam-qualification-copy-proof-v001.json').read_text())
for r in copy['initial_copy_files_verified'].values():a[r['source']]=r['sha256']
for n in ['analyse_vortex_openfoam_v002.py','plot_vortex_shedding_v002.py']:
 p=R/'repo/cfd-learning-lab/scripts'/n;a[str(p)]=sha(p)
(S/'openfoam-qualification-inputs-v001.json').write_text(json.dumps(a,indent=2)+'\n')
assert b'\r' not in (S/'run_cylinder_openfoam_qualification_v001.sh').read_bytes()
print('Frozen qualification sources and inputs: '+str(len(a)))