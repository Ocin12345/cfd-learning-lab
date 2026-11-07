"""Audit accepted native cavity fields, coupling, heat reports and the same-time two-grid comparison."""
from pathlib import Path
import json,re,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from cavity_native_fields_v001 import *
from cavity_production_common_v004 import ROOT,CP,K,HOT,COLD,DEPTH
LOCAL=ROOT/'local/natural-convection-production-v001';OUT=ROOT/'exports/figures/2026-11-07-natural-convection-20s-paired-v001'
F=ROOT/'native/runs/natural-convection-cavity-48-v001/fluent';O=ROOT/'native/runs/natural-convection-cavity-48-v001/openfoam-v002'

def main():
 fsteps=json.loads((F/'time-step-checks.json').read_text());osteps=json.loads((LOCAL/'openfoam-time-step-checks-v002.json').read_text());assert len(fsteps)==1000 and all(x['passed'] for x in fsteps) and osteps['all_passed']
 assert (O/'run-status.txt').read_text().strip()=='complete'
 sources=[F/'qualification-final.cas.h5',F/'qualification-final.dat.h5',F/'time-step-checks.json',F/'heated-settings.json',F/'heated-gravity-work-source.json',F/'wall-heat-final.txt',F/'solver-transcript.txt',O/'solver.log',O/'wall-heat-postprocess.log',O/'system/fvSolution',O/'system/controlDict',O/'constant/physicalProperties']
 rows={};fields={};checks={}
 for name,frames in [('ANSYS Fluent',json.loads((F/'frames.json').read_text())),('OpenFOAM',[dict(time_s=i*.5) for i in range(41)])]:
  values=[];fields[name]={}
  for frame in frames:
   t=frame['time_s']
   if name=='ANSYS Fluent':
    data=F/frame['file'];q=fluent(F/'qualification-final.cas.h5',data,48);sources.append(data)
   else:
    q=foam(O,t,48);sources.extend(O/f'{t:g}'/n for n in ['U','p','T']+(['rho'] if t else []))
   v=metrics(q,48);v['time_s']=t;values.append(v);fields[name][t]=q
  initial=values[0]['mass_per_depth_kg_m']
  for v in values:v['mass_drift']=v['mass_per_depth_kg_m']/initial-1
  rows[name]=values
  checks[name]=dict(saved_states=len(values),temperature_bounded=all(COLD-1e-4<=v['Tmin_K']<=v['Tmax_K']<=HOT+1e-4 for v in values),mass=all(abs(v['mass_drift'])<=1e-5 for v in values),eos=all(v['eos_relative_error']<=1e-6 for v in values),max_mass_drift=max(abs(v['mass_drift']) for v in values),last_half_second_relative_velocity_change=float(np.linalg.norm(fields[name][20.]['U']-fields[name][19.5]['U'])/np.linalg.norm(fields[name][20.]['U'])))
  assert all(checks[name][k] for k in ['temperature_bounded','mass','eos'])
 cf=ROOT/'native/runs/natural-convection-qualification-fluent-v005';co=ROOT/'native/runs/natural-convection-qualification-openfoam-v003/heated'
 assert json.loads((LOCAL/'openfoam-coarse-time-step-checks-v001.json').read_text())['all_passed']
 coarse={'ANSYS Fluent':fluent(cf/'qualification-final.cas.h5',cf/'heated-0050.dat.h5',24),'OpenFOAM':foam(co,1.,24)}
 mesh={name:mesh_comparison(coarse[name],fields[name][1.]) for name in coarse}
 sources.extend([cf/'qualification-final.cas.h5',cf/'heated-0050.dat.h5',co/'1/U',co/'1/T',co/'1/p',co/'1/rho',co/'solver.log'])
 heat={};raw=(F/'wall-heat-final.txt').read_text()
 for name,log,depth in [('ANSYS Fluent',raw,1.),('OpenFOAM',(O/'wall-heat-postprocess.log').read_text(),DEPTH)]:
  wall={}
  for w in ['top','bottom','left','right']:
   if name=='ANSYS Fluent':v=float(re.search(r'^\s*'+w+r'\s+([-+\d.eE]+)\s*$',log,re.M)[1])
   else:v=float(re.search(r'for patch '+w+r' = [^,]+, [^,]+, ([^,]+),',log)[1])
   wall[w]=v/depth
  heat[name]=dict(wall_heat_per_depth_W_m=wall,hot_wall_Nu=wall['left']/(K*(HOT-COLD)),net_wall_heat_per_depth_W_m=sum(wall.values()),net_wall_heat_fraction_of_hot=abs(sum(wall.values()))/abs(wall['left']),heat_balance_is_transient_not_steady_proof=True)
 raw=(O/'solver.log').read_text();Co=[float(m[1]) for m in re.finditer(r'Courant Number mean: [^ ]+ max: ([-+\d.eE]+)',raw)]
 summary=dict(checks=checks,accepted_steps={'ANSYS Fluent':1000,'OpenFOAM':1000},openfoam_outer_correctors=16,mesh_comparison_at_1s=mesh,heat_at_20s=heat,max_openfoam_Courant=max(Co),physical_validation=False,mesh_independence=False,time_step_independence=False,gravity_work_alignment_verified=True,fluent_native_export_recovered_without_advancing_solution='E:/CFD-Learning-Lab/exports/figures/2026-11-07-natural-convection-20s-native-v001/export-record.json')
 with Image.open(OUT/'convection-development.gif') as im:
  assert im.n_frames==41
  for i in range(41):im.seek(i);im.load();assert im.info['duration']==250
 sources.extend([ROOT/'native/runs/vertical-energy-control-fluent-v005/energy-check.json',ROOT/'exports/figures/2026-11-07-gravity-work-control-v001/evidence.json'])
 summary['native_source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
 (OUT/'final-audit.json').write_text(json.dumps(summary,indent=2)+'\n')
 (OUT/'saved-state-metrics.json').write_text(json.dumps(rows,indent=2)+'\n')
 fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
 for name,v in rows.items():
  axs[0,0].plot([x['time_s'] for x in v],[x['Umax_m_s']*1000 for x in v],label=name)
  axs[0,1].plot([x['time_s'] for x in v],[x['mass_drift']*1e6 for x in v],label=name)
 for eq,i,lim in [('continuity',0,1e-6),('x-velocity',1,1e-6),('y-velocity',2,1e-6),('energy',3,1e-9)]:axs[1,0].semilogy([x['time_s'] for x in fsteps],[x['residuals'][i]/lim for x in fsteps],label=eq)
 for eq,lim in osteps['targets'].items():axs[1,1].semilogy([x['time_s'] for x in osteps['steps']],[max(x['residuals'][eq]/lim,1e-12) for x in osteps['steps']],label=eq)
 axs[0,0].set(title='Circulation develops and relaxes',ylabel='Maximum speed (mm/s)');axs[0,1].set(title='Closed-cavity mass',ylabel='Mass change from heating start (ppm)')
 axs[1,0].set(title='Fluent: each time-step endpoint',ylabel='Reported residual / declared limit');axs[1,1].set(title='OpenFOAM: final outer corrector',ylabel='Initial residual / declared limit')
 for ax in axs[1]:ax.axhline(1,ls='--',color='black',lw=1)
 for ax in axs.ravel():ax.set_xlabel('Physical time (s)');ax.grid(alpha=.2);ax.legend(fontsize=8)
 fig.suptitle('Native numerical checks | solver residual norms are not interchangeable',fontsize=14);fig.savefig(OUT/'numerical-checks.png',dpi=170);fig.savefig(OUT/'numerical-checks.pdf');plt.close(fig)
 print(json.dumps({k:v for k,v in summary.items() if k!='native_source_sha256'},indent=2))
if __name__=='__main__':main()
