from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from analyse_nozzle_openfoam_v002 import field,FIELDS
from wedge_foam_mapping_v002 import native_mapping
R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-mach2-wedge-spatial-v001';O=N/'openfoam-late-cfl-v001';P=N/'openfoam-fine-v001';S=Path(__file__).parent

def main():
 a=json.loads((O/'endpoint-audit.json').read_text());assert a['numerical_endpoint_pass']
 for n,h in a['source_sha256'].items():assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==h,n
 out=R/'exports/figures/2026-11-08-wedge-late-cfl-v002';out.mkdir(parents=True,exist_ok=False)
 with np.load(R/'native/meshes/showcase-mach2-wedge-spatial-v001/fine-v001/geometry.npz') as g:
  C=field(O/'0/C',89600,*FIELDS['C']);order,metrics=native_mapping(C,g['centres'],.01)
  pm=field(P/'0.008/p',89600,*FIELDS['p'])[order].reshape(224,400);ph=field(O/'0.008/p',89600,*FIELDS['p'])[order].reshape(224,400);delta=(ph-pm)/101325*100
  plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
  fig,axes=plt.subplots(2,2,figsize=(11.5,9),layout='constrained');x=g['x_nodes'];y=g['y_nodes'];X=np.broadcast_to(x,y.shape)
  for ax,p,title in [(axes[0,0],pm,'Main fine run | maxCo 0.30'),(axes[0,1],ph,'Saved 7 ms checkpoint branch | maxCo 0.15')]:
   im=ax.pcolormesh(X,y,p/101325,cmap='cividis',norm=Normalize(.98,1.74),shading='flat',rasterized=True);ax.fill_between(x,-.02,y[0],color='#233744');ax.set(xlim=(0,1),ylim=(-.02,.82),aspect='equal',xlabel='x (m)',ylabel='y (m)',title=title)
  fig.colorbar(im,ax=axes[0,:].tolist(),label='p / freestream p',shrink=.8)
  ax=axes[1,0];im=ax.pcolormesh(X,y,delta,cmap='RdBu_r',norm=Normalize(-8,8),shading='flat',rasterized=True);ax.fill_between(x,-.02,y[0],color='#233744');ax.set(xlim=(0,1),ylim=(-.02,.82),aspect='equal',xlabel='x (m)',ylabel='y (m)',title='Endpoint pressure difference: branch − main');fig.colorbar(im,ax=ax,label='Difference / freestream pressure (%)',extend='both',shrink=.8)
  records=a['courant']['records'];ax=axes[1,1];ax.plot([1000*r['pre_advance_time_s'] for r in records],[r['native_raw_courant'] for r in records],c='#a85124',lw=1.1,label='All native pre-advance reports')
  ax.scatter([7,7],a['courant']['raw_startup_courant'],s=32,c='#243e50',label='Retained inherited startup (two reports)');ax.axhline(.153,color='#777',ls='--',lw=.8,label='Established reduced-CFL bound');ax.set(xlim=(6.98,8.01),ylim=(0,.32),xlabel='Actual pre-advance state time (ms)',ylabel='Maximum native Courant number',title='Inherited deltaT is explicit, not discarded');ax.legend(fontsize=8);ax.grid(alpha=.2)
  fig.suptitle('A bounded late-time step-sensitivity check',fontsize=17,weight='bold')
  fig.supxlabel('Same 89,600-cell mesh; exact saved 7 ms state, five native states through 8 ms. This does not test the full trajectory or establish time independence.',fontsize=9)
  fig.savefig(out/'wedge-late-cfl-native-evidence.png',dpi=175);fig.savefig(out/'wedge-late-cfl-native-evidence.pdf');plt.close(fig)
 record=dict(source_audit=str(O/'endpoint-audit.json'),audit_sha256=hashlib.sha256((O/'endpoint-audit.json').read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),signed_pressure_difference_percent_min=float(delta.min()),signed_pressure_difference_percent_max=float(delta.max()),difference_colour_limits_percent=[-8,8],clipped_cell_count=int((abs(delta)>8).sum()),native_data_sources_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [P/'0.008/p',O/'0.008/p',O/'0/C']},actual_inspection_pending=True)
 (out/'source-binding.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
if __name__=='__main__':main()