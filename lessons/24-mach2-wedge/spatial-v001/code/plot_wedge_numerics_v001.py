from pathlib import Path
import json,hashlib,re,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyse_nozzle_fluent_v005 import residual_rows
from wedge_endpoint_v001 import load
from plot_wedge_fluent_fine_v001 import measure
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent;N=R/'native/runs/showcase-mach2-wedge-spatial-v001';F=N/'fluent-fine-v001';O=N/'openfoam-fine-v001'
def main():
 out=R/'exports/figures/2026-11-08-wedge-spatial-numerics-v001';out.mkdir(parents=True,exist_ok=False)
 checkpoints=json.loads((F/'checkpoints.json').read_text());history=residual_rows(F/checkpoints[-1]['data']);windows=[]
 with np.load(F/'geometry.npz') as geo:
  previous=None
  for cp in checkpoints:
   current=load(F/cp['case'],F/cp['data'],geo)
   if previous is not None:
    m=measure(current,previous,geo)
    windows.append(dict(iteration=cp['iteration'],**m['relative_field_change_l2']))
   previous=current
 log=(O/'solver.log').read_text();number=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
 times=np.array([float(t) for t in re.findall(r'^Time = ('+number+r')s?\s*$',log,re.M)]);co=np.array([float(t) for t in re.findall(r'^Courant Number mean: '+number+r' max: ('+number+r')',log,re.M)])
 assert len(co)==len(times)+1
 # Constructor at t0 plus loop report at t0; each subsequent report is native pre-advance time.
 co_times=np.r_[0.,0.,times[:-1]];assert len(co_times)==len(co)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(2,2,figsize=(12.2,8),layout='constrained')
 for key in ('continuity','x-velocity','y-velocity','energy'):
  axes[0,0].semilogy([r['iteration'] for r in history],[r[key] for r in history],label=key,lw=1)
 axes[0,0].axhline(1e-5,c='#555',ls='--',lw=.8,label='flow limit 1e-5');axes[0,0].axhline(1e-6,c='#888',ls=':',lw=.8,label='energy limit 1e-6');axes[0,0].axvline(120,c='#999',lw=.7)
 axes[0,0].set(xlabel='Fluent steady iteration',ylabel='Native scaled residual',title='Fine Fluent: native residual history');axes[0,0].legend(fontsize=8,ncol=2)
 for key in ('p','mach'):axes[0,1].semilogy([r['iteration'] for r in windows],[r[key] for r in windows],'.-',label=key)
 axes[0,1].axhline(.002,c='#555',ls='--',label='original limit 0.002');axes[0,1].set(xlabel='Ending steady iteration',ylabel='Volume-weighted relative L2 change',title='Saved 200-iteration windows');axes[0,1].legend(fontsize=8)
 axes[1,0].plot(1000*co_times,co,c='#a7512a',lw=.8);axes[1,0].axhline(.306,c='#555',ls='--',label='original 1.02 × 0.3 bound');axes[1,0].set(xlabel='Native pre-advance state time (ms)',ylabel='Maximum native Courant number',title='Fine OpenFOAM: every printed Courant report');axes[1,0].legend(fontsize=8)
 selected=json.loads((S/'wedge-selected-measurements-v001.json').read_text());labels=[r['code']+'\n'+r['resolution'].lower() for r in selected];mass=[r['native']['mass_imbalance'] for r in selected]
 axes[1,1].bar(labels,mass,color=['#4a8494','#0b667a','#cd9777','#ae5729']);axes[1,1].set_yscale('log');axes[1,1].axhline(.001,c='#555',ls='--',label='original limit 0.001');axes[1,1].set(ylabel='Absolute boundary imbalance / inlet',title='Endpoint native boundary mass balance');axes[1,1].legend(fontsize=8)
 for ax in axes.ravel():ax.grid(alpha=.2)
 fig.suptitle('Numerical evidence for the accepted comparison',fontsize=17,weight='bold')
 fig.supxlabel('OpenFOAM diagonal-solve zeros are not nonlinear convergence. Late 7–8 ms pressure change = 0.0019297; criterion = 0.002.',fontsize=9)
 fig.savefig(out/'wedge-numerical-checks.png',dpi=175);fig.savefig(out/'wedge-numerical-checks.pdf');plt.close(fig)
 data=dict(fluent_native_residual_rows=history,fluent_saved_window_changes=windows,openfoam_native_pre_advance_courant_records=[dict(time_s=float(t),maximum_courant=float(c)) for t,c in zip(co_times,co)],endpoint_mass_balance=mass,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [F/checkpoints[-1]['data'],F/'checkpoints.json',O/'solver.log',S/'wedge-selected-measurements-v001.json',Path(__file__)]})
 (out/'numerical-data.json').write_text(json.dumps(data,indent=2)+'\n');print('Rendered native numerical evidence',len(history),'Fluent history rows',len(co),'native Courant records')
if __name__=='__main__':main()