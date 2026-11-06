"""Inspect actual Gmsh triangles in two cropped views; no synthetic mesh."""
from pathlib import Path
import argparse,hashlib,json
import gmsh
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--mesh',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
a=p.parse_args()
if a.output_dir.exists():raise FileExistsError('Use a fresh versioned figure directory')
a.output_dir.mkdir(parents=True)
gmsh.initialize();gmsh.open(str(a.mesh))
tags,xyz,_=gmsh.model.mesh.getNodes()
xy=np.asarray(xyz).reshape(-1,3)[:,:2]
tags=np.asarray(tags)
order=np.argsort(tags)
_,element_nodes=gmsh.model.mesh.getElementsByType(2)
tri=xy[order[np.searchsorted(tags[order],np.asarray(element_nodes))]].reshape(-1,3,2)
gmsh.finalize()
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':9.5,
    'axes.linewidth':.65,'xtick.direction':'in','ytick.direction':'in',
    'xtick.top':True,'ytick.right':True,'pdf.fonttype':42,'svg.fonttype':'path'})
fig,axes=plt.subplots(1,2,figsize=(7.05,3),layout='constrained')
for ax,bounds,label in zip(axes,[(-.03,1.03,-.13,.13),(.990,1.002,-.004,.004)],['(a)','(b)']):
    x0,x1,y0,y1=bounds
    lo=tri.min(axis=1);hi=tri.max(axis=1)
    visible=tri[(hi[:,0]>=x0)&(lo[:,0]<=x1)&(hi[:,1]>=y0)&(lo[:,1]<=y1)]
    segments=np.stack([visible,visible[:,[1,2,0],:]],axis=2).reshape(-1,2,2)
    ax.add_collection(LineCollection(segments,colors='#34495e',linewidths=.23,rasterized=False))
    ax.set(xlim=(x0,x1),ylim=(y0,y1),xlabel=r'$x/c$',ylabel=r'$y/c$')
    ax.set_aspect('equal',adjustable='box')
    ax.text(.02,1.02,label,transform=ax.transAxes,ha='left',va='bottom')
base=a.output_dir/'fig05-airfoil-mesh'
fig.savefig(base.with_suffix('.pdf'),bbox_inches='tight')
fig.savefig(base.with_suffix('.svg'),bbox_inches='tight')
fig.savefig(base.with_suffix('.png'),dpi=600,bbox_inches='tight')
plt.close(fig)
caption=('Figure 5. Actual unstructured triangular mesh for the preserved, closed NACA 0012 polygon at zero incidence. '
 '(a) Near-body view; (b) trailing-edge detail. Coordinates are normalized by the 1 m chord. '
 f'The complete domain contains {len(tri):,} triangles and extends to radius 20 m; only cropped views are shown. '
 'Locally prescribed edge sizes resolve the original polygon vertices. Mesh quality and flow accuracy are distinct checks; '
 'the image alone establishes neither convergence nor mesh independence. Gmsh generated the mesh; Codex prepared this visualization.\n')
(a.output_dir/'CAPTION.md').write_text(caption,encoding='utf-8')
(a.output_dir/'source-manifest.json').write_text(json.dumps({a.mesh.name:hashlib.sha256(a.mesh.read_bytes()).hexdigest()},indent=2)+'\n',encoding='utf-8')
print(caption)
