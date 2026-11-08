from pathlib import Path
import csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
S=Path(__file__).parent;N=Path('E:/CFD-Learning-Lab/native/runs/showcase-cylinder-convergence-v001');O=Path('E:/CFD-Learning-Lab/exports/figures/2026-11-08-cylinder-numerical-record-v002');assert not O.exists();O.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
D={};sources={};summary={}
for solver,label,color in [('fluent','ANSYS Fluent','#1765a0'),('openfoam','OpenFOAM','#b26729')]:
 a=S/('analysis-'+solver+'-v001');forces=list(csv.DictReader((a/(solver+'-force-history.csv')).open()));metrics=json.loads((a/(solver+'-cycle-metrics.json')).read_text());p=list(csv.DictReader((S/'analysis-pressure-v001'/(solver+'-wall-pressure-estimates.csv')).open()));assert len(forces)==len(p)==801
 if solver=='fluent':res=list(csv.DictReader((a/'fluent-step-residuals.csv').open()));rt=[float(r['time_s']) for r in res];rv=[max(float(r[k]) for k in ['continuity','x-velocity','y-velocity']) for r in res];mt=[float(r['time_s']) for r in forces];mv=[float(r['relative_mass_imbalance']) for r in forces]
 else:
  res=json.loads((a/'openfoam-step-records.json').read_text());rt=[r['time_s'] for r in res];rv=[max(r['variables'][k]['terminal_initial'] for k in ['Ux','Uy','p']) for r in res];mass={}
  for folder in ['openfoam-qualification-v001','openfoam-developed-v001']:
   b=json.loads((N/folder/'window-audit.json').read_text())
   for r in b['fields_and_mass']:mass.setdefault(round(r['time_s'],8),r['mass']['relative_all_physical_boundary_imbalance'])
  mt=sorted(mass);mv=[mass[t] for t in mt]
 D[solver]=dict(label=label,color=color,t=[float(r['time_s']) for r in forces],Cd=[float(r['Cd']) for r in forces],Cl=[float(r['Cl']) for r in forces],pressure=[float(r['difference_Pa']) for r in p],metrics=metrics,residual_times=rt,residual=rv,mass_times=mt,mass=mv)
 for f in [a/(solver+'-force-history.csv'),a/(solver+'-cycle-metrics.json'),S/'analysis-pressure-v001'/(solver+'-wall-pressure-estimates.csv')]:sources[str(f)]=sha(f)
 summary[solver]=dict(mean_period_s=metrics['last_six_mean_period_s'],St=.1/metrics['last_six_mean_period_s'],time_mean_Cd=metrics['last_six_mean_Cd'],mean_sampled_Cl_peak_to_peak=float(np.mean([c['sampled_Cl_peak_to_peak'] for c in metrics['last_six_cycles']])),last_six_stable=metrics['last_six_cycles_stable_under_declared_criterion'],relative_cycle_ranges=metrics['last_six_relative_ranges'],maximum_relative_boundary_mass_imbalance=max(mv),new_converged_steps=800,native_states=801)
f,axes=plt.subplots(3,1,figsize=(12.8,8.4),layout='constrained')
for d in D.values():
 for ax,key in zip(axes,['Cd','Cl','pressure']):ax.plot(d['t'],d[key],color=d['color'],lw=1.1,label=d['label'])
for ax,label in zip(axes,['Drag coefficient, $C_D$','Lift coefficient, $C_L$','Estimated front–rear pressure (Pa)']):ax.set(xlabel='Physical time (s)',ylabel=label,xlim=(4,8));ax.grid(alpha=.16);ax.legend(fontsize=9,loc='lower right',bbox_to_anchor=(1,1.02),ncol=2,frameon=False) if ax is axes[0] else None;ax.ticklabel_format(axis='y',useOffset=False)
f.suptitle('Cylinder wake: native forces and wall-pressure estimates',fontsize=16);f.supxlabel('Both codes at their own native times; no phase alignment. Pressure uses matching adjacent wall-face angular stencils.\nSpatial wall-pressure interpolation and cycle-boundary temporal interpolation are disclosed in the method record.\nConverged continuation from inherited startup; benchmark matching and physical validation remain open.',fontsize=9)
f.savefig(O/'forces-and-pressure.png',dpi=140);f.savefig(O/'forces-and-pressure.pdf');plt.close(f)
f,axes=plt.subplots(2,2,figsize=(12.8,8),layout='constrained')
for j,solver in enumerate(['fluent','openfoam']):
 d=D[solver];ax=axes[0,j];ax.semilogy(d['residual_times'],d['residual'],color=d['color'],lw=.7);ax.axhline(1e-5,color='#272727',ls='--',lw=1);ax.set(xlabel='Physical time (s)',ylabel='Maximum declared terminal residual',title=d['label']+' | all 800 new steps');ax.grid(alpha=.15)
axes[0,0].text(.02,.04,'Scaled continuity / momentum',transform=axes[0,0].transAxes,fontsize=9);axes[0,1].text(.02,.87,'Last initial linear p / Ux / Uy residuals',transform=axes[0,1].transAxes,fontsize=9)
for d in D.values():
 axes[1,0].semilogy(d['mass_times'],np.maximum(d['mass'],1e-16),color=d['color'],lw=.8,label=d['label']);cycles=d['metrics']['complete_cycles'];axes[1,1].plot([.5*(c['start_time_s']+c['end_time_s']) for c in cycles],[c['sampled_Cl_peak_to_peak'] for c in cycles],'o-',ms=3,lw=1,color=d['color'],label=d['label'])
axes[1,0].axhline(1e-4,color='#272727',ls='--',lw=1,label='Unchanged limit');axes[1,0].set(xlabel='Physical time (s)',ylabel='Relative boundary mass imbalance',title='Signed native flux checks at every state');axes[1,0].legend(fontsize=9);axes[1,0].grid(alpha=.15)
axes[1,1].set(xlabel='Cycle midpoint (s)',ylabel='Sampled lift peak-to-peak',title='Cycle-by-cycle amplitude development');axes[1,1].legend(fontsize=9);axes[1,1].grid(alpha=.15);axes[1,1].ticklabel_format(axis='y',useOffset=False)
f.suptitle('Numerical checks and cycle development',fontsize=16);f.supxlabel('Residual definitions differ between codes; magnitudes are not an accuracy ranking. Thresholds are unchanged.\nLast six cycles meet the declared 0.5% range/mean criterion for period, lift amplitude and time-mean drag.\nCycle stability does not establish mesh/time independence, a matched benchmark pass or physical validation.',fontsize=9)
f.savefig(O/'numerical-checks.png',dpi=140);f.savefig(O/'numerical-checks.pdf');plt.close(f)
pa=json.loads((S/'analysis-pressure-v001/wall-pressure-audit.json').read_text())
for solver in summary:summary[solver]['last_six_mean_wall_pressure_difference_Pa']=pa['summary'][solver]['last_six_mean_difference_Pa']
(O/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(O/'plot-receipt.json').write_text(json.dumps(dict(source_sha256=sources,code_sha256=sha(Path(__file__)),figure_sha256={p.name:sha(p) for p in O.iterdir() if p.suffix in ['.png','.pdf']},zero_mass_points_plotted_at_floor_1e_minus16=sum(sum(float(x)==0 for x in d['mass']) for d in D.values()),physical_validation=False),indent=2)+'\n');print(json.dumps(summary,indent=2))