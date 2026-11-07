"""Compare actual fields at matched times; do not infer independence from two steps."""
from pathlib import Path
import json,re,hashlib,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from cavity_native_fields_v001 import fluent,foam,metrics
from cavity_production_common_v004 import ROOT,K,HOT,COLD,DEPTH
LOCAL=ROOT/'local/natural-convection-time-step-v001'
OUT=ROOT/'exports/figures/2026-11-07-natural-convection-dt001-paired-v001'
BASE=ROOT/'native/runs/natural-convection-cavity-48-v001'
NEW=ROOT/'native/runs/natural-convection-dt001-v001'
F=NEW/'fluent-v002';O=NEW/'openfoam'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 assert json.loads((F/'status.json').read_text())['stage']=='closed'
 assert (O/'run-status.txt').read_text().strip()=='complete'
 fs=json.loads((F/'time-step-checks.json').read_text());os=json.loads((LOCAL/'openfoam-steps.json').read_text())
 assert len(fs)==2000 and [v['step'] for v in fs]==list(range(1,2001)) and all(v['passed'] for v in fs)
 assert os['all_passed'] and os['steps_found']==2000
 assert np.max(abs(np.array([x['time_s'] for x in fs])-np.arange(1,2001)*.01))<1e-8
 assert np.max(abs(np.array([x['time_s'] for x in os['steps']])-np.arange(1,2001)*.01))<1e-8
 for n,h in json.loads((F/'source-hashes.json').read_text()).items():assert sha(BASE/'fluent'/n)==h
 copy=json.loads((LOCAL/'copied-openfoam-inputs.json').read_text())
 for item in copy['files']:
  assert sha(BASE/'openfoam-v002'/item['path'])==item['original_sha256']
  if item['path']!='system/controlDict':assert sha(O/item['path'])==item['original_sha256']
 expected=(BASE/'openfoam-v002/system/controlDict').read_text().replace('deltaT 0.02;','deltaT 0.01;').replace('writeInterval 25;','writeInterval 50;')
 assert (O/'system/controlDict').read_text()==expected
 checks={};rows=[];mass=[];source=[];fields={};initial_comparison={}
 for name in ['ANSYS Fluent','OpenFOAM']:
  fields[name]={};frames=json.loads((F/'frames.json').read_text()) if name=='ANSYS Fluent' else [dict(time_s=i*.5) for i in range(41)]
  samples=[]
  for frame in frames:
   t=frame['time_s']
   if name=='ANSYS Fluent':
    p=F/frame['file'];q=fluent(F/'qualification-final.cas.h5',p,48);source.append(p)
   else:
    q=foam(O,t,48);source.extend(O/f'{t:g}'/n for n in ['U','p','T']+(['rho'] if t else []))
   fields[name][round(t,8)]=q;m=dict(time_s=t,**metrics(q,48));samples.append(m)
  m0=samples[0]['mass_per_depth_kg_m']
  for m in samples:m['mass_relative_change']=m['mass_per_depth_kg_m']/m0-1;mass.append(dict(solver=name,**m))
  checks[name]=dict(steps=2000,saved_states=len(samples),temperature_bounded=all(COLD-1e-4<=m['Tmin_K']<=m['Tmax_K']<=HOT+1e-4 for m in samples),max_mass_relative_drift=max(abs(m['mass_relative_change']) for m in samples),max_eos_error=max(m['eos_relative_error'] for m in samples))
  assert checks[name]['temperature_bounded'] and checks[name]['max_mass_relative_drift']<=1e-5 and checks[name]['max_eos_error']<=1e-6
  for i in range(41):
   t=i*.5;q=fields[name][t]
   if name=='ANSYS Fluent':
    p=BASE/'fluent'/f'heated-{round(t/.02):04d}.dat.h5';b=fluent(BASE/'fluent/qualification-final.cas.h5',p,48);source.append(p)
   else:
    b=foam(BASE/'openfoam-v002',t,48);source.extend(BASE/'openfoam-v002'/f'{t:g}'/n for n in ['U','p','T']+(['rho'] if t else []))
   if t==0:
    initial_comparison[name]={key:float(np.max(abs(q[key]-b[key]))) for key in q}
    for key in q:assert np.allclose(q[key],b[key],rtol=1e-14,atol=1e-14),(name,key,initial_comparison[name])
   rows.append(dict(solver=name,time_s=t,velocity_relative_L2=float(np.linalg.norm(q['U']-b['U'])/np.linalg.norm(q['U'])),max_temperature_difference_K=float(np.max(abs(q['T']-b['T']))),baseline_max_speed_mm_s=float(np.linalg.norm(b['U'],axis=1).max()*1000),refined_max_speed_mm_s=float(np.linalg.norm(q['U'],axis=1).max()*1000)))
 heat={}
 for name,folder,depth in [('ANSYS Fluent',F,1.),('OpenFOAM',O,DEPTH)]:
  log=(folder/('wall-heat-final.txt' if name=='ANSYS Fluent' else 'wall-heat-postprocess.log')).read_text();wall={}
  for w in ['left','right','top','bottom']:
   pattern=r'^\s*'+w+r'\s+([-+\d.eE]+)\s*$' if name=='ANSYS Fluent' else r'for patch '+w+r' = [^,]+, [^,]+, ([^,]+),'
   wall[w]=float(re.search(pattern,log,re.M)[1])/depth
  heat[name]=dict(wall_heat_per_depth_W_m=wall,hot_wall_Nu=wall['left']/(K*(HOT-COLD)),net_heat_W_m=sum(wall.values()))
 bh=json.loads((ROOT/'exports/figures/2026-11-07-natural-convection-20s-paired-v001/final-audit.json').read_text())['heat_at_20s']
 for name in heat:
  heat[name]['baseline_hot_wall_Nu']=bh[name]['hot_wall_Nu'];heat[name]['relative_Nu_change']=(heat[name]['hot_wall_Nu']/bh[name]['hot_wall_Nu']-1)
 co=[float(x) for x in re.findall(r'Courant Number mean: [^ ]+ max: ([-+\d.eE]+)',(O/'solver.log').read_text())]
 summary=dict(checks=checks,initial_field_max_absolute_differences=initial_comparison,input_preservation_verified=True,baseline_dt_s=.02,refined_dt_s=.01,physical_duration_s=20,cells=2304,accepted_steps={'ANSYS Fluent':2000,'OpenFOAM':2000},heat_at_20s=heat,max_OpenFOAM_Courant=max(co),final_comparisons=[v for v in rows if v['time_s']==20],max_sampled_velocity_change={name:max((v for v in rows if v['solver']==name and v['time_s']>0),key=lambda v:v['velocity_relative_L2']) for name in fields},relative_velocity_difference_denominator='L2 norm of dt=0.01 velocity at same time; all matched native cells',physical_validation=False,mesh_independence=False,time_step_independence=False,limitation='Two time steps quantify sensitivity but do not establish asymptotic convergence. Mesh sensitivity remains unresolved; no transient global energy audit or benchmark validation claimed.')
 source.extend([F/'qualification-final.cas.h5',F/'qualification-final.dat.h5',BASE/'fluent/qualification-final.cas.h5',F/'heated-settings.json',F/'restart-settings-check.json',F/'solver-transcript.txt',F/'wall-heat-final.txt',F/'time-step-checks.json',O/'solver.log',O/'wall-heat-postprocess.log',O/'system/controlDict',O/'system/fvSolution',O/'constant/physicalProperties'])
 summary['native_source_sha256']={str(p):sha(p) for p in source}
 (OUT/'time-step-audit.json').write_text(json.dumps(summary,indent=2)+'\n');(OUT/'time-step-history.json').write_text(json.dumps(rows,indent=2)+'\n');(OUT/'saved-state-metrics.json').write_text(json.dumps(mass,indent=2)+'\n')
 with (OUT/'time-step-history.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 fig,axs=plt.subplots(3,2,figsize=(12,12),layout='constrained');colors={'ANSYS Fluent':'#1565a8','OpenFOAM':'#b65718'}
 for j,name in enumerate(fields):
  rr=[v for v in rows if v['solver']==name];tt=[v['time_s'] for v in rr]
  axs[0,j].plot(tt,[v['baseline_max_speed_mm_s'] for v in rr],color=colors[name],lw=2,label='dt = 0.02 s')
  axs[0,j].plot(tt,[v['refined_max_speed_mm_s'] for v in rr],color='black',lw=1.2,ls='--',label='dt = 0.01 s')
  axs[0,j].set(title=name,ylabel='Maximum speed (mm/s)')
  axs[1,0].plot(tt,[v['velocity_relative_L2']*100 for v in rr],label=name,color=colors[name]);axs[1,1].plot(tt,[v['max_temperature_difference_K']*1000 for v in rr],label=name,color=colors[name])
 axs[1,0].set(title='Sensitivity at matched physical times',ylabel='Velocity L2 difference (%)');axs[1,1].set(title='Largest cellwise temperature change',ylabel='Maximum |T(0.01) - T(0.02)| (mK)')
 for eq,i,lim in [('continuity',0,1e-6),('x velocity',1,1e-6),('y velocity',2,1e-6),('energy',3,1e-9)]:axs[2,0].semilogy([v['time_s'] for v in fs],[max(v['residuals'][i]/lim,1e-12) for v in fs],lw=.8,label=eq)
 for eq,lim in os['targets'].items():axs[2,1].semilogy([v['time_s'] for v in os['steps']],[max(v['residuals'][eq]/lim,1e-12) for v in os['steps']],lw=.8,label=eq)
 axs[2,0].set(title='Fluent: each smaller time-step endpoint',ylabel='Native residual / declared limit');axs[2,1].set(title='OpenFOAM: final outer corrector',ylabel='Native initial residual / declared limit')
 for ax in axs[2]:ax.axhline(1,color='black',ls='--',lw=1)
 for ax in axs.ravel():ax.set_xlabel('Physical time (s)');ax.grid(alpha=.2);ax.legend(fontsize=8)
 fig.suptitle('Time-step sensitivity | same 2,304 cells and initial fields\nResidual definitions differ between solvers',fontsize=15);fig.savefig(OUT/'time-step-comparison.png',dpi=160);fig.savefig(OUT/'time-step-comparison.pdf');plt.close(fig)
 with Image.open(OUT/'convection-development.gif') as im:
  assert im.n_frames==41
  for i in range(41):im.seek(i);im.load();assert im.info['duration']==250
 print(json.dumps({k:v for k,v in summary.items() if k!='native_source_sha256'},indent=2))
if __name__=='__main__':main()
