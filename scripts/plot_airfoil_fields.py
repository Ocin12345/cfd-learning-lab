"""Academic field figure read directly from saved triangular Fluent HDF5 data.

Reconstruct triangles using saved face/cell IDs; never assume Gmsh cell order.
Only the supported single-zone, 2D, triangular Re100 teaching case is accepted.
"""
from pathlib import Path
import argparse,hashlib,json
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--run',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--cp-min',type=float,default=-2.)
p.add_argument('--cp-max',type=float,default=2.)
a=p.parse_args()
assert np.isfinite([a.cp_min,a.cp_max]).all() and a.cp_min<0<a.cp_max
if a.output_dir.exists():raise FileExistsError('Choose a fresh figure output directory')
settings=json.loads((a.run/'loaded-settings.json').read_text())
geometry=json.loads((a.run/'geometry-setup.json').read_text())
assert geometry['angle_degrees']==0 and geometry['reynolds_number']==100
radius=geometry['outer_radius_m']/geometry['chord_m']
refs=settings['references'];rho=settings['material']['density']['value'];uin=refs['velocity']
assert rho==uin==refs['length']==1 and refs['pressure']==0
case=a.run/'airfoil-Re100-final.cas.h5';data=a.run/'airfoil-Re100-final.dat.h5'
with h5py.File(case,'r') as f:
    mesh=f['meshes/1'];xy=mesh['nodes/coords/1'][:]
    assert xy.shape[1]==2 and int(mesh['nodes/coords/1'].attrs['minId'][0])==1
    assert np.all(mesh['faces/nodes/1/nnodes'][:]==2)
    edges=mesh['faces/nodes/1/nodes'][:].reshape(-1,2).astype(np.int64)-1
    c0=mesh['faces/c0/1'][:].astype(np.int64)-1
    c1field=mesh['faces/c1/1'];c1=c1field[:].astype(np.int64)-1
    begin=int(c1field.attrs['minId'][0])-1
    count=int(mesh['cells/zoneTopology/maxId'][0])
    assert int(mesh['cells/zoneTopology/minId'][0])==1
    pairs=np.concatenate([np.column_stack([np.repeat(c0,2),edges.ravel()]),
        np.column_stack([np.repeat(c1,2),edges[begin:begin+len(c1)].ravel()])])
    unique=np.unique(pairs,axis=0)
    ids,counts=np.unique(unique[:,0],return_counts=True)
    assert np.array_equal(ids,np.arange(count)) and np.all(counts==3)
    triangles=unique[:,1].reshape(-1,3)
with h5py.File(data,'r') as f:
    fields=f['results/1/phase-1/cells']
    def field(name):
        dataset=fields[name+'/1']
        assert int(dataset.attrs['minId'][0])==1 and len(dataset)==count
        return dataset[:]
    assert np.allclose(field('SV_DENSITY'),rho) and np.allclose(field('SV_MU_LAM'),.01)
    speed=np.hypot(field('SV_U'),field('SV_V'))/uin
    cp=(field('SV_P')-refs['pressure'])/(.5*rho*uin**2)
assert np.isfinite(speed).all() and np.isfinite(cp).all()
coords=xy[triangles]
ab=coords[:,1]-coords[:,0];ac=coords[:,2]-coords[:,0]
cross=ab[:,0]*ac[:,1]-ab[:,1]*ac[:,0]
assert np.all(np.abs(cross)>0)
triangles[cross<0]=triangles[cross<0][:,[0,2,1]]
assert speed.min()>=0 and speed.max()<=1.2 and cp.min()>=a.cp_min and cp.max()<=a.cp_max
triangulation=mtri.Triangulation(xy[:,0],xy[:,1],triangles)
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':9.5,
 'axes.linewidth':.65,'xtick.direction':'in','ytick.direction':'in',
 'xtick.top':True,'ytick.right':True,'pdf.fonttype':42,'svg.fonttype':'path'})
fig,axes=plt.subplots(2,1,figsize=(7.05,4.6),layout='constrained')
for ax,values,limits,cmap,label,panel in zip(axes,[speed,cp],[(0,1.2),(a.cp_min,a.cp_max)],
 ['viridis','RdBu_r'],[r'$|\mathbf{u}|/U_\infty$',r'$C_p$'],['(a)','(b)']):
    norm=(matplotlib.colors.TwoSlopeNorm(vcenter=0,vmin=limits[0],vmax=limits[1])
          if panel=='(b)' else matplotlib.colors.Normalize(*limits))
    art=ax.tripcolor(triangulation,facecolors=values,shading='flat',
        norm=norm,cmap=cmap,rasterized=True)
    ax.set(xlim=(-.3,2.5),ylim=(-.55,.55),xlabel=r'$x/c$',ylabel=r'$y/c$',aspect='equal')
    ax.text(.01,1.025,panel,transform=ax.transAxes,ha='left',va='bottom')
    cb=fig.colorbar(art,ax=ax,pad=.02,fraction=.035);cb.set_label(label)
a.output_dir.mkdir(parents=True)
base=a.output_dir/'fig06-airfoil-computed-fields'
fig.canvas.draw()
fig.set_layout_engine('none')
for suffix in ['pdf','svg','png']:fig.savefig(base.with_suffix('.'+suffix),dpi=600)
(a.output_dir/'executed-plotter.py').write_bytes(Path(__file__).read_bytes())
plt.close(fig)
caption=('Figure 6. Saved Fluent solution for the closed NACA 0012 polygon at Re=100, zero incidence, '
 f'{count:,} triangular cells and outer radius {radius:g}c. (a) Speed normalized by U_inf=1 m/s; '
 '(b) Cp=(p-p_ref)/(rho U_inf^2/2), with rho=1 kg/m3 and p_ref=0 Pa gauge. '
 'Each filled triangle uses its actual saved cell-centred value; no spatial smoothing or remeshing is applied. '
 f'The solid interior remains blank. Fixed speed range0..1.2 and Cp range[{a.cp_min:g},{a.cp_max:g}] include all saved values; zero Cp is the divergent colour midpoint. '
 'Only cropped near-body/wake views are shown. Field layers are rasterized at 600 dpi within PDF/SVG vector axes. '
 'Fluent computed the fields; Codex reconstructed the saved face/cell connectivity and prepared this figure. '
 'Numerical checks do not establish mesh/domain independence or physical validation.\n')
(a.output_dir/'CAPTION.md').write_text(caption,encoding='utf-8')
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [case,data,a.run/'loaded-settings.json',a.run/'geometry-setup.json']}
(a.output_dir/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(cells=count,speed_range=[float(speed.min()),float(speed.max())],
 Cp_range=[float(cp.min()),float(cp.max())],connectivity_verified=True),indent=2))
