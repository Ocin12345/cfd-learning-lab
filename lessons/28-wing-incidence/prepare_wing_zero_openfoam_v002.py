"""Create a fresh zero-incidence OpenFOAM case by copying the verified native mesh."""
from pathlib import Path
import hashlib,json,shutil
from wing_zero_common_v001 import U,LIFT
R=Path('E:/CFD-Learning-Lab');P=R/'native/runs/2026-11-01-finite-wing-v003/openfoam'
O=R/'native/runs/finite-wing-zero-openfoam-v002'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def vec(v):return '('+' '.join(f'{a:.17g}' for a in v)+')'
def main():
    assert not O.exists() and shutil.disk_usage(R).free>10*2**30
    O.mkdir();hashes={}
    for n in ['constant/polyMesh/'+s for s in ['points','faces','owner','neighbour','boundary']]+[
        'constant/physicalProperties','constant/momentumTransport','system/fvSchemes','system/fvSolution','system/decomposeParDict','0/p']:
        p=P/(n.replace('constant/polyMesh/','0/polyMesh/') if n.startswith('constant/polyMesh/') else n);q=O/n;q.parent.mkdir(parents=True,exist_ok=True);h=sha(p);shutil.copyfile(p,q);assert sha(p)==sha(q)==h;hashes[n]=h
    text=(P/'0/U').read_text();assert text.count('(1 0 0)')==3
    (O/'0/U').write_text(text.replace('(1 0 0)',vec(U)),encoding='ascii',newline='\n')
    text=(P/'system/controlDict').read_text();assert text.count('liftDir (0 1 0);')==1 and text.count('dragDir (1 0 0);')==1
    text=text.replace('liftDir (0 1 0);','liftDir '+vec(LIFT)+';').replace('dragDir (1 0 0);','dragDir '+vec(U)+';')
    (O/'system/controlDict').write_text(text,encoding='ascii',newline='\n')
    (O/'copied-source-hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
    (O/'declared-settings.json').write_text(json.dumps(dict(effective_incidence_deg=0,body_clockwise_pitch_deg=8,U=U,liftDir=LIFT,
        maximum_iterations=1200,source_mesh_unchanged=True,residual_target=1e-5,mass_target=1e-4,
        final100_force_relative_range_target=.001,force_reference_floor=.001,physical_validation=False),indent=2)+'\n')
    print('Prepared',O)
if __name__=='__main__':main()
