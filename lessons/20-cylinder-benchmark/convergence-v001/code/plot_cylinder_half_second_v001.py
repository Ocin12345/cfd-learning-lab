from pathlib import Path
import csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from plot_vortex_shedding_v002 import fluent_mesh,native_field
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');N=Path('E:/CFD-Learning-Lab/native/runs/showcase-cylinder-convergence-v001');O=Path('E:/CFD-Learning-Lab/exports/figures/2026-11-08-cylinder-half-second-v001');assert not O.exists();O.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
variants=[('Coupled','fluent-coupled-window-v001','#167d83')];audits={};sources={}
for label,name,col in variants:
 a=json.loads((N/name/'window-audit.json').read_text());assert a['steps']==100 and a['saved_states']==101
 for p,h in a['native_source_sha256'].items():assert sha(Path(p))==h;sources[p]=h
 audits[label]=a
assert audits['Coupled']['numerical_window_pass'] and audits['Coupled']['all_runtime_guards_pass']
case=N/'fluent-coupled-window-v001/vortex-final.cas.h5';data=N/'fluent-coupled-window-v001/vortex-final.dat.h5';mesh=fluent_mesh(case);field=native_field('fluent',dict(path=data),mesh)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'normal','savefig.facecolor':'white'})
f=plt.figure(figsize=(12.8,7.6),layout='constrained');gs=f.add_gridspec(2,2,height_ratios=[.8,1.2]);ax=f.add_subplot(gs[0,:]);vmax=float(np.ceil(field['speed'].max()*10)/10);im=ax.tripcolor(mesh['triangulation'],facecolors=field['speed'],cmap='viridis',vmin=0,vmax=vmax,rasterized=True);ax.set(xlim=(0,2.2),ylim=(0,.41),xlabel='x (m)',ylabel='y (m)',title='Actual Fluent speed at t = 4.550 s · Coupled continuation');ax.set_aspect('equal');f.colorbar(im,ax=ax,label='Speed (m/s)',fraction=.025,pad=.018)
a1=f.add_subplot(gs[1,0]);a2=f.add_subplot(gs[1,1])
for label,name,col in variants:
 a=audits[label];times=[r['time_s'] for r in a['step_records']];res=[r['continuity'] for r in a['step_records']];a1.plot(times,res,'o-',ms=3,lw=1.4,color=col,label=label)
 rows=list(csv.DictReader((N/name/'force-history.csv').open()));a2.plot([float(r['time_s']) for r in rows],[float(r['Cl']) for r in rows],'.-',ms=3,lw=1.2,color=col,label=label)
a1.axhline(1e-5,color='#272727',ls='--',lw=1,label='Unchanged target');a1.set(xlabel='Physical time (s)',ylabel='Terminal scaled continuity residual',title='Every newly computed time step',ylim=(0,1.5e-5));a1.ticklabel_format(axis='y',style='sci',scilimits=(0,0));a1.legend(fontsize=8,loc='lower left');a1.grid(alpha=.16)
a2.set(xlabel='Physical time (s)',ylabel='Lift coefficient, $C_L$',title='Actual native force samples');a2.ticklabel_format(axis='y',useOffset=False);a2.grid(alpha=.16);a2.legend(fontsize=8)
f.suptitle('Cylinder wake: the longer convergence check',fontsize=17,x=.5);f.supxlabel('21,598 cells · Re = 100 · dt = 0.005 s · continues the accepted 4.05 s checkpoint · 120-inner cap\nCoupled passes all 100 new steps. Inherited startup misses, developed-cycle evidence and benchmark matching remain open.',fontsize=9)
f.savefig(O/'cylinder-window.png',dpi=160);f.savefig(O/'cylinder-window.pdf');plt.close(f)
record=dict(actual_native_fields=True,solver='ANSYS Fluent',physical_time_s=4.55,source_sha256=sources,plot_source_sha256=sha(Path(__file__)),speed_max_m_s=float(field['speed'].max()),plot_speed_scale_max_m_s=vmax,convergence_result={k:dict(passed=a['numerical_window_pass'],miss_steps=a['residual_miss_steps'],mass_imbalance=a['maximum_relative_mass_imbalance']) for k,a in audits.items()},physical_validation=False,benchmark_pass=False,figure_sha256={p.name:sha(p) for p in O.iterdir() if p.is_file()})
(O/'plot-receipt.json').write_text(json.dumps(record,indent=2)+'\n');print(str(O/'cylinder-window.png'))