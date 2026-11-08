"""Replay selected native fields, all-column widths, fixed-scale figures. No solvers."""
from pathlib import Path
import hashlib,json,sys,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LogNorm
from wedge_width_columns_v001 import measure_columns
from wedge_physical_gradient_v001 import gradient_magnitude
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-mach2-wedge-spatial-v001'
LESSON=R/'repo/cfd-learning-lab/lessons/24-mach2-wedge'
sys.path.append(str(LESSON))
import plot_wedge_fluent_v001 as coarse
import plot_wedge_fluent_fine_v001 as fine
from analyse_plot_wedge_v002 import read_openfoam as read_coarse
from wedge_openfoam_fine_io_v001 import read_openfoam as read_fine

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(a):
 for n,h in a.items():assert sha(n)==h,n

def collect():
 verify(json.loads((S/'baseline-native-hashes-v001.json').read_text()))
 for label in ('fluent','openfoam'):
  a=json.loads((N/(label+'-fine-v001')/'endpoint-audit.json').read_text());assert a['numerical_endpoint_pass'] and not a['failed_checks'];verify(a['source_sha256'])
 rows=[]
 for code,resolution,run in [('Fluent','Coarse',R/'native/runs/mach2-wedge-fluent-v001'),('Fluent','Fine',N/'fluent-fine-v001'),('OpenFOAM','Coarse',R/'native/runs/mach2-wedge-openfoam-v001'),('OpenFOAM','Fine',N/'openfoam-fine-v001')]:
  mod=coarse if resolution=='Coarse' else fine
  gp=(run/'geometry.npz') if code=='Fluent' else (R/'native/runs/mach2-wedge-fluent-v001/geometry.npz' if resolution=='Coarse' else R/'native/meshes/showcase-mach2-wedge-spatial-v001/fine-v001/geometry.npz')
  with np.load(gp) as z:g={n:z[n] for n in z.files}
  if code=='Fluent':data,previous,native=mod.read_fluent(run,R/'repo/cfd-learning-lab/scripts',g)
  else:
   times,frames,masses,native=(read_coarse if resolution=='Coarse' else read_fine)(run,g);data,previous=frames[-1],frames[-5]
  metrics=mod.measure(data,previous,g);width=measure_columns(data['p'],g,metrics)
  xy=g['centres'].reshape(*data['p'].shape,2)
  data['gradient']=gradient_magnitude(data['rho'],xy[...,0],xy[...,1])
  rows.append(dict(code=code,resolution=resolution,cells=data['p'].size,run=str(run),geo=g,data=data,metrics=metrics,width=width,native=native,geometry_file=str(gp),geometry_sha256=sha(gp)))
 return rows

def save_fig(fig,out,name):
 fig.savefig(out/(name+'.png'),dpi=175);fig.savefig(out/(name+'.pdf'));plt.close(fig)

def render(rows,out):
 out.mkdir(parents=True,exist_ok=False)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'#193245','axes.labelcolor':'#193245','axes.spines.top':False,'axes.spines.right':False})
 # Four accepted cases, common scales, direct native-cell colour mapping.
 fig,axes=plt.subplots(2,4,figsize=(15,8),layout='constrained',sharex=True,sharey=True)
 for j,row in enumerate(rows):
  g=row['geo'];x=g['x_nodes'];y=g['y_nodes'];X=np.broadcast_to(x,y.shape)
  for i,(key,cmap,norm) in enumerate([('pressure','cividis',Normalize(.98,1.74)),('mach','viridis',Normalize(1.6,2.03))]):
   ax=axes[i,j];im=ax.pcolormesh(X,y,row['data'][key],cmap=cmap,norm=norm,shading='flat',rasterized=True)
   ax.fill_between(x,-.025,y[0],color='#233744');ax.set(xlim=(0,1),ylim=(-.02,.82),aspect='equal',xlabel='x (m)')
   if j==0:ax.set_ylabel('y (m)')
   if i==0:ax.set_title(f"{row['code']} | {row['resolution'].lower()}\n{row['cells']:,} cells",fontsize=10)
   if j==3:fig.colorbar(im,ax=axes[i,:].tolist(),shrink=.75,label='p / freestream p' if i==0 else 'Mach number')
 fig.suptitle('Mach 2 over a 10° wedge: one refinement in each solver',fontsize=17,weight='bold')
 fig.supxlabel('Fluent: saved steady endpoints. OpenFOAM: 8 ms; identical maxCo = 0.3 gives a smaller fine-grid time step.',fontsize=10)
 save_fig(fig,out,'wedge-coarse-fine-native-fields')
 fig,axes=plt.subplots(2,2,figsize=(10.8,8.5),layout='constrained',sharex=True,sharey=True)
 for k,row in enumerate(rows):
  ax=axes[k//2,k%2];g=row['geo'];x=g['x_nodes'];y=g['y_nodes']
  im=ax.pcolormesh(np.broadcast_to(x,y.shape),y,row['data']['gradient'],cmap='gray_r',norm=LogNorm(.3,250),shading='flat',rasterized=True)
  ax.set(xlim=(.4,.8),ylim=(.10,.51),aspect='equal',xlabel='x (m)',ylabel='y (m)',title=f"{row['code']} | {row['resolution'].lower()}, {row['cells']:,} cells")
 fig.colorbar(im,ax=axes.ravel().tolist(),label=r'Numerical density-gradient magnitude (kg m$^{-4}$)',shrink=.8)
 fig.suptitle('The captured shock occupies fewer metres on the finer grid',fontsize=16,weight='bold')
 fig.supxlabel('Physical-coordinate Jacobian gradient of native cell-centred density; fixed colour scale and view. Numerical transition only.',fontsize=9)
 save_fig(fig,out,'wedge-coarse-fine-shock-detail')
 fig,axes=plt.subplots(1,2,figsize=(12,5.5),layout='constrained')
 for row in rows:
  color='#0b667a' if row['code']=='Fluent' else '#ae5729';style='--' if row['resolution']=='Coarse' else '-';label=row['code']+' '+row['resolution'].lower()
  w=row['width'];axes[0].plot([r['x_mean_m'] for r in w['columns']],[1000*r['width_m'] for r in w['columns']],style,c=color,label=label,lw=1.6)
  p=row['data']['p'];ny,nx=p.shape;xy=row['geo']['centres'].reshape(ny,nx,2);j=np.argmin(abs(xy[...,0].mean(axis=0)-.6));m=row['metrics']
  d=(m['shock_fit_slope']*xy[:,j,0]+m['shock_fit_intercept_m']-xy[:,j,1])/np.sqrt(1+m['shock_fit_slope']**2)
  q=(p[:,j]-m['pressure_upstream_pa'])/(m['pressure_downstream_pa']-m['pressure_upstream_pa']);order=np.argsort(d)
  axes[1].plot(1000*d[order],q[order],style,c=color,label=label,lw=1.4,marker='.',ms=3)
 axes[0].set(xlabel='Native column mean x (m)',ylabel='Connected 10–90% projected width (mm)',title='All prescribed columns; no dropped profiles')
 axes[1].set(xlim=(-20,20),ylim=(-.08,1.08),xlabel='Signed distance from each fitted shock (mm)',ylabel='Observed pressure jump fraction',title='Native column nearest x = 0.60 m')
 for q in (.1,.9):axes[1].axhline(q,color='#9ca8ae',lw=.7)
 for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=9)
 fig.suptitle('Measuring numerical shock sharpness',fontsize=17,weight='bold')
 fig.supxlabel('Linear interpolation between native centres. Profile distance is relative to each case’s own fitted line. Not molecular shock thickness.',fontsize=9)
 save_fig(fig,out,'wedge-shock-width-and-profiles')

if __name__=='__main__':
 out=R/'exports/figures/2026-11-08-wedge-spatial-paired-v001'
 rows=collect()
 # All-column data and native bindings are saved before figures are inspected.
 payload=[]
 for r in rows:payload.append({k:v for k,v in r.items() if k not in ('geo','data')})
 selection=S/'wedge-selected-measurements-v001.json';assert not selection.exists();selection.write_text(json.dumps(payload,indent=2)+'\n')
 render(rows,out)
 (out/'selected-source-binding.json').write_text(json.dumps(dict(measurement_record=str(selection),sha256=sha(selection),script_sha256=sha(__file__),plots={p.name:sha(p) for p in out.glob('*') if p.is_file()},actual_inspection_pending=True),indent=2)+'\n')
 print(json.dumps([dict(code=r['code'],resolution=r['resolution'],median_width_mm=1000*r['width']['median_width_m'],all_columns=r['width']['valid_columns'],angle_deg=r['metrics']['shock_angle_deg'],pressure_ratio=r['metrics']['pressure_ratio']) for r in rows],indent=2))