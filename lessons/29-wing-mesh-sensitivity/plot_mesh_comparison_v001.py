from pathlib import Path
import ast,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
R=Path('E:/CFD-Learning-Lab');L=R/'local/finite-wing-mesh-sensitivity-v001'
for p in L.glob('*.py'):ast.parse(p.read_text(encoding='utf-8-sig'))
out=R/'exports/figures/2026-11-07-wing-refined-mesh-v001';assert not out.exists();out.mkdir()
paths=[R/'native/runs/2026-11-01-finite-wing-v003/mesh-export-v002',R/'native/meshes/finite-wing-refined-v002']
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
hashes={}
for j,p in enumerate(paths):
    with (p/'geometry.npz').open('rb') as f:hashes[str(p/'geometry.npz')]=hashlib.file_digest(f,'sha256').hexdigest()
    g=np.load(p/'geometry.npz');xyz=g['points_xyz_m'];tri=xyz[g['wing_face_nodes']];normal=g['wing_normals'];cells=len(g['tets'])
    top=tri[normal[:,1]>.15][:,:,[0,2]]
    for i,limits in enumerate([((-.03,1.03),(-2.05,2.05)),((.86,1.01),(1.7,2.03))]):
        ax=axs[i,j];ax.add_collection(PolyCollection(top,facecolors='#edf3f7',edgecolors='#426b88',linewidths=.23))
        ax.set(xlim=limits[0],ylim=limits[1],xlabel='x (m)',ylabel='span z (m)',title=('Baseline' if j==0 else 'Refined')+f' | {cells:,} cells'+(' | tip/trailing-edge detail' if i else ''))
        ax.set_aspect('equal' if i else 'auto');ax.spines[['top','right']].set_visible(False)
fig.suptitle('Actual upper wing surface meshes | same CAD, 8 deg incidence\nNative tetrahedral boundary triangles projected onto the x-z plane',fontsize=13)
fig.savefig(out/'wing-mesh-comparison.png',dpi=170);fig.savefig(out/'wing-mesh-comparison.svg');plt.close(fig)
(out/'native-source-hashes-private.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(out/'wing-mesh-comparison.png')
