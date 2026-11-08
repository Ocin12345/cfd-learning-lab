from pathlib import Path
import csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from plot_vortex_shedding_v002 import fluent_mesh,openfoam_mesh,native_field
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-cylinder-convergence-v001';A=S/'analysis-fluent-v001';O=R/'exports/figures/2026-11-08-cylinder-paired-qualification-v001';assert not O.exists();O.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
fa=json.loads((A/'combined-fluent-audit.json').read_text());oa=json.loads((N/'openfoam-qualification-v001/window-audit.json').read_text());assert fa['all_new_steps_pass'] and oa['numerical_window_pass'] and oa['steps']==100
for a,key in [(fa,'source_sha256'),(oa,'native_source_sha256')]:
 for p,h in a[key].items():assert sha(Path(p))==h,p
snap=next(r for r in json.loads((A/'fluent-snapshots.json').read_text()) if abs(r['time_s']-4.5)<1e-9);fm=fluent_mesh(Path(snap['case']));om=openfoam_mesh(N/'openfoam-qualification-v001');ff=native_field('fluent',dict(path=Path(snap['data'])),fm);of=native_field('openfoam',dict(path=N/'openfoam-qualification-v001/4.5'),om);vmax=float(np.ceil(max(ff['speed'].max(),of['speed'].max())*10)/10)
m=json.loads((A/'fluent-cycle-metrics.json').read_text());assert m['last_six_cycles_stable_under_declared_criterion']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
f=plt.figure(figsize=(12.8,9.0),layout='constrained');gs=f.add_gridspec(3,2,height_ratios=[.85,.85,1.1])
for k,label,mesh,field in [(0,'ANSYS Fluent',fm,ff),(1,'OpenFOAM Foundation 14',om,of)]:
 ax=f.add_subplot(gs[k,:]);im=ax.tripcolor(mesh['triangulation'],facecolors=field['speed'],cmap='viridis',vmin=0,vmax=vmax,rasterized=True);ax.set(xlim=(0,2.2),ylim=(0,.41),xlabel='x (m)',ylabel='y (m)',title=label+' · actual speed at t = 4.500 s');ax.set_aspect('equal');f.colorbar(im,ax=ax,label='Speed (m/s)',fraction=.025,pad=.018)
a1=f.add_subplot(gs[2,0]);a2=f.add_subplot(gs[2,1])
for label,path,color in [('ANSYS Fluent',A/'fluent-force-history.csv','#1765a0'),('OpenFOAM',N/'openfoam-qualification-v001/force-history.csv','#b26729')]:
 rows=[r for r in csv.DictReader(path.open()) if 4+1e-8<float(r['time_s'])<=4.5+1e-8];assert len(rows)==100;a1.plot([float(r['time_s']) for r in rows],[float(r['Cl']) for r in rows],'.-',color=color,lw=1.1,ms=2.5,label=label)
a1.set(xlabel='Physical time (s)',ylabel='Lift coefficient, $C_L$',title='Native forces at their own physical times');a1.grid(alpha=.16);a1.legend(fontsize=9)
cycles=m['complete_cycles'];tx=[.5*(r['start_time_s']+r['end_time_s']) for r in cycles];amps=[r['sampled_Cl_peak_to_peak'] for r in cycles];a2.plot(tx,amps,'o-',ms=4,color='#167d83');a2.axvspan(tx[-6],tx[-1],color='#167d83',alpha=.08);a2.set(xlabel='Cycle midpoint (s)',ylabel='Sampled lift peak-to-peak',title='Fluent cycle amplitude over 4–8 s');a2.grid(alpha=.16);a2.ticklabel_format(axis='y',useOffset=False)
a2.text(.03,.96,'Last six-cycle relative range: %.4f%%'%(100*m['last_six_relative_ranges']['sampled_Cl_peak_to_peak']),transform=a2.transAxes,va='top',fontsize=9)
f.suptitle('Cylinder wake: stronger convergence evidence',fontsize=17);f.supxlabel('21,598 cells · Re = 100 · dt = 0.005 s · shared speed scale · no phase alignment\nFluent: 800 new steps to 8 s pass. OpenFOAM: 100-step qualification passes; its 8 s extension is pending.\nThe declared cycle-stability check is not a matched benchmark pass, discretization independence or physical validation.',fontsize=9)
f.savefig(O/'paired-cylinder-check.png',dpi=160);f.savefig(O/'paired-cylinder-check.pdf');plt.close(f)
(O/'plot-receipt.json').write_text(json.dumps(dict(actual_native_fields=True,time_s=4.5,speed_scale_max_m_s=vmax,fluent_full_steps_passed=800,openfoam_qualification_steps_passed=100,fluent_cycle_metrics=m,source_audits={str(p):sha(p) for p in [A/'combined-fluent-audit.json',A/'fluent-cycle-metrics.json',N/'openfoam-qualification-v001/window-audit.json']},source_code_sha256=sha(Path(__file__)),figure_sha256={p.name:sha(p) for p in O.iterdir() if p.is_file()},physical_validation=False),indent=2)+'\n');print(str(O/'paired-cylinder-check.png'))