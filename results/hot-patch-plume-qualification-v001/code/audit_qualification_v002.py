from pathlib import Path
import json,re,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plume_common_v002 import *
from plume_native_fields_v001 import fluent,foam,metrics
F=ROOT/'native/runs/hot-patch-plume-qualification-v001/fluent-v003';O=F.parent/'openfoam';LOCAL=ROOT/'local/hot-patch-plume-v001';OUT=ROOT/'exports/figures/2026-11-07-hot-patch-plume-qualification-v001'
def main():
 assert json.loads((F/'status.json').read_text())['stage']=='closed'
 assert (O/'run-status.txt').read_text().strip()=='complete'
 fs=json.loads((F/'time-step-checks.json').read_text());os=json.loads((LOCAL/'openfoam-steps.json').read_text());assert len(fs)==100 and all(x['passed'] for x in fs) and os['all_passed']
 assert all(json.loads((F.parent/'fluent/control-checks.json').read_text())['checks'].values())
 assert all(json.loads((F.parent/'openfoam-control/control-checks.json').read_text())['checks'].values())
 mesh=O/'constant/polyMesh';mesh=mesh if mesh.exists() else O/'0/polyMesh'
 b=(mesh/'boundary').read_text();assert int(re.search(r'heater\s*\{[^{}]*nFaces\s+(\d+)',b,re.S)[1])==12
 assert not OUT.exists();OUT.mkdir(parents=True);rows={};fields={};sources=[]
 for name,frames in [('ANSYS Fluent',json.loads((F/'frames.json').read_text())),('OpenFOAM',[dict(time_s=i*.2) for i in range(6)])]:
  rows[name]=[];fields[name]={}
  for frame in frames:
   t=frame['time_s']
   if name=='ANSYS Fluent':p=F/frame['file'];q=fluent(F/'qualification-final.cas.h5',p);sources.append(p)
   else:q=foam(O,t);sources.extend(O/f'{t:g}'/n for n in ['U','T','p']+(['rho'] if t else []))
   m=dict(time_s=t,**metrics(q));rows[name].append(m);fields[name][round(t,8)]=q
  m0=rows[name][0]['mass_per_depth_kg_m']
  for m in rows[name]:
   m['relative_mass_drift']=abs(m['mass_per_depth_kg_m']/m0-1)
   assert m['relative_mass_drift']<=1e-5 and m['eos_relative_error']<=1e-6 and COLD-1e-4<=m['Tmin_K']<=m['Tmax_K']<=HOT+1e-4
 f=fields['ANSYS Fluent'][1.];o=fields['OpenFOAM'][1.]
 a=dict(specification=specification(),accepted_steps={'ANSYS Fluent':100,'OpenFOAM':100},quiet_controls_pass=True,fluent_pressure_momentum_relaxation=[1.,1.],prior_attempts=[dict(folder="fluent",missed_steps=[1],inner_cap=200),dict(folder="fluent-v002",missed_steps=[4,5,9],inner_cap=800)],all_saved_field_checks_pass=True,heater_faces=12,heater_area_per_depth_m=12*DX,geometric_mapping='native cell centres and spans checked bijectively against 48x80 shared grid',final_comparison=dict(max_T_difference_K=float(np.max(abs(f['T']-o['T']))),velocity_relative_L2=float(np.linalg.norm(f['U']-o['U'])/np.linalg.norm(f['U']))),max_mass_drift={n:max(v['relative_mass_drift'] for v in r) for n,r in rows.items()},records=rows,qualification_only=True,physical_validation=False,mesh_independence=False,time_step_independence=False)
 sources.extend([F/'qualification-final.cas.h5',F/'qualification-final.dat.h5',F/'time-step-checks.json',F/'heated-settings.json',F/'heated-gravity-work-source.json',F/'source-hashes.json',F/'restart-settings-check.json',O/'solver.log',O/'system/fvSolution',O/'system/controlDict',O/'constant/physicalProperties',mesh/'boundary'])
 a['native_source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};write(OUT/'qualification-audit.json',a)
 fig,axs=plt.subplots(1,2,figsize=(9,7.8));fig.subplots_adjust(left=.09,right=.85,top=.84,bottom=.14,wspace=.34)
 x=(np.arange(NX)+.5)*DX*1000;y=(np.arange(NY)+.5)*DY*1000
 for ax,name in zip(axs,['ANSYS Fluent','OpenFOAM']):
  q=fields[name][1.];T=q['T'].reshape(NY,NX);U=q['U'].reshape(NY,NX,2)
  im=ax.pcolormesh(np.linspace(0,W*1000,NX+1),np.linspace(0,H*1000,NY+1),T-T0,cmap='inferno',vmin=0,vmax=1,shading='flat',rasterized=True)
  ax.quiver(x[2::5],y[2::5],U[2::5,2::5,0],U[2::5,2::5,1],color='#71ccff',scale=.03,width=.004)
  ax.plot([HEATER_LEFT*1000,HEATER_RIGHT*1000],[0,0],color='#ff392e',lw=5,clip_on=False)
  ax.set(title=name,xlabel='x (mm)',ylabel='y (mm)',xlim=(0,W*1000),ylim=(0,H*1000));ax.set_aspect('equal')
 fig.colorbar(im,cax=fig.add_axes([.88,.18,.022,.6]),label='Temperature above 300 K (K)')
 fig.suptitle('A localized heater starts the flow',fontsize=20,fontweight='bold');fig.text(.09,.89,'Actual native fields at t = 1.0 s | 15 mm heater | Ra(H) = 99,238',fontsize=11)
 fig.text(.09,.04,'Initial qualification only; the developed plume has not yet been simulated.\n3,840 cells | dt = 0.01 s | matched arrow scale | no physical validation',fontsize=10)
 fig.savefig(OUT/'paired-startup.png',dpi=160);fig.savefig(OUT/'paired-startup.pdf');plt.close(fig)
 print(json.dumps({k:v for k,v in a.items() if k not in ['records','native_source_sha256']},indent=2))
if __name__=='__main__':main()
