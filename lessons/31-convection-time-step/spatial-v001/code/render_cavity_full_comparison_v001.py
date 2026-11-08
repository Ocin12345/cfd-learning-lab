"""Select audited native segments; render actual states, no solver or temporal synthesis."""
from pathlib import Path
import csv,json,hashlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image
from cavity_spatial_common_v001 import K,HOT,COLD
from cavity_native_fields_v001 import fluent,foam,metrics
from cavity_native_heat_v001 import native_fluent_time,fluent_wall_flux
from cavity_openfoam_heat_v001 import native_foam_time,foam_walls
from cavity_restriction_v001 import compare_fields
from cavity_wall_gradient_v001 import wall_gradient_heat
from wing_fluent_v001 import available_memory
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent;N=R/'native/runs/showcase-convection-spatial-v001';B=R/'native/runs/natural-convection-dt001-v001';FQ=N/'fluent-qualification-v004';OQ=N/'openfoam-qualification-v001';F=N/'fluent-production-v002';O=N/'openfoam-production-v001'
TIMES=np.arange(41)*.5

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(a):
 for p,h in a.items():assert sha(p)==h,p

def collect():
 verify(json.loads((S/'baseline-native-source-hashes-v001.json').read_text()));audits={}
 for name,p in [('Fluent',F/'production-audit-v001.json'),('OpenFOAM',O/'production-audit-v001.json')]:
  a=json.loads(p.read_text());assert a['numerical_pass'] and a['new_steps']==1900 and a['merged_steps']==2000;verify(a['source_sha256']);audits[name]=a
 selected={};bindings={};histories=[];responses=[];heat_validation=[]
 for code in ['Fluent','OpenFOAM']:
  for grid,n in [('Coarse',48),('Fine',96)]:
   key=code+' '+grid;frames=[];m0=None
   for t in TIMES:
    step=int(round(t/.01));sources=[]
    if code=='Fluent':
     case=(B/'fluent-v002') if grid=='Coarse' else (FQ if t<=1 else F)
     mesh=case/('continuation-final.cas.h5' if case==F else 'qualification-final.cas.h5');data=case/f'heated-{step:04d}.dat.h5';clock=native_fluent_time(data);assert abs(clock['time_s']-t)<1e-8 and clock['step']==step and abs(clock['dt_s']-.01)<1e-12
     q=fluent(mesh,data,n);sources=[mesh,data]
    else:
     case=(B/'openfoam') if grid=='Coarse' else (OQ if t<=1 else O)
     folder=case/f'{t:g}';q=foam(case,float(t),n);sources=[folder/k for k in ['U','T','p']+(['rho'] if t else [])]
     if t:clock=native_foam_time(folder/'uniform/time',float(t));assert clock['step']==step;sources.append(folder/'uniform/time')
    assert all(np.isfinite(v).all() for v in q.values()) and all((q[k]>0).all() for k in ['p','rho','T'])
    m=metrics(q,n)
    if m0 is None:m0=m['mass_per_depth_kg_m']
    m['mass_relative_drift']=abs(m['mass_per_depth_kg_m']/m0-1);assert m['mass_relative_drift']<=1e-5 and m['eos_relative_error']<=1e-6 and COLD-1e-4<=m['Tmin_K']<=m['Tmax_K']<=HOT+1e-4
    heat=wall_gradient_heat(q['T'].reshape(n,n),K) if t>0 else None
    histories.append(dict(code=code,grid=grid,time_s=float(t),**m,Q_hot_gradient_reconstruction_W_m=heat['left'] if heat else None,Q_cold_gradient_reconstruction_W_m=heat['right'] if heat else None,Nu_gradient_reconstruction=heat['left']/(K*(HOT-COLD)) if heat else None))
    frames.append(q);bindings[key+f'@{t:g}s']=dict(case=str(case),native_sources=[str(p) for p in sources],source_sha256={str(p):sha(p) for p in sources},initial_density_method='native rho except prepared OpenFOAM0s uses native p/(R*T)' if code=='OpenFOAM' and t==0 else 'native rho')
   selected[key]=frames
   for t in ([1.,20.] if grid=='Fine' else [20.]):
    q=frames[int(t*2)];derived=wall_gradient_heat(q['T'].reshape(n,n),K)
    if grid=='Fine':
     a=audits[code] if t==20 else json.loads(((FQ/'qualification-audit-v002.json') if code=='Fluent' else (OQ/'qualification-audit-v001.json')).read_text());walls=a['native_wall_heat_at20s'] if t==20 else a['native_wall_heat']
    elif code=='Fluent':walls,_=fluent_wall_flux(B/'fluent-v002/qualification-final.cas.h5',B/'fluent-v002/qualification-final.dat.h5',B/'fluent-v002/wall-heat-final.txt',['left','right','top','bottom'])
    else:walls,_=foam_walls(B/'openfoam',['left','right','top','bottom'],20.)
    for name in ['left','right']:
     exact=walls[name]['heat_per_depth_W_m'];error=abs(derived[name]-exact);assert error<=1e-10+1e-6*abs(exact),(key,t,name,derived[name],exact)
     heat_validation.append(dict(code=code,grid=grid,time_s=t,wall=name,native_face_heat_W_m=exact,reconstructed_heat_W_m=derived[name],absolute_difference_W_m=error,relative_difference=error/abs(exact)))
 for code in ['Fluent','OpenFOAM']:
  for i,t in enumerate(TIMES):responses.append(dict(code=code,time_s=float(t),**compare_fields(selected[code+' Coarse'][i],selected[code+' Fine'][i],48,relative=t>0)))
 peak=max(np.linalg.norm(q['U'],axis=1).max() for f in selected.values() for q in f);speed_limit_mm_s=float(np.ceil(peak*10000)/10)
 return selected,dict(times_s=TIMES.tolist(),source_bindings=bindings,history=histories,mesh_responses=responses,native_heat_reconstruction_checks=heat_validation,global_native_peak_speed_m_s=float(peak),speed_colour_limit_mm_s=speed_limit_mm_s,physical_validation=False,full_transient_energy_closure_claimed=False)

def save(fig,out,name):
 for ext in ('png','pdf'):fig.savefig(out/(name+'.'+ext),dpi=150)
 plt.close(fig)

def maps(ax,q,n,field,speed_limit,stream=False):
 edges=np.linspace(0,50,n+1);array=q['T'].reshape(n,n) if field=='T' else np.linalg.norm(q['U'],axis=1).reshape(n,n)*1000
 im=ax.pcolormesh(edges,edges,array,cmap='coolwarm' if field=='T' else 'viridis',norm=Normalize(299.6,300.4) if field=='T' else Normalize(0,speed_limit),shading='flat',rasterized=True)
 ax.set(aspect='equal',xlim=(0,50),ylim=(0,50),xlabel='x (mm)',ylabel='y (mm)')
 if stream and np.linalg.norm(q['U'],axis=1).max()*1000>.005*speed_limit:
  centres=(np.arange(n)+.5)*50/n;U=q['U'].reshape(n,n,2);seed=np.array([(x,y) for y in np.linspace(7,43,5) for x in np.linspace(7,43,5)])
  ax.streamplot(centres,centres,U[...,0],U[...,1],color='white',linewidth=.55,arrowsize=.65,density=.9,start_points=seed,broken_streamlines=True)
 return im

def render(selected,a,out):
 out.mkdir(parents=True,exist_ok=False);framesdir=out/'frames';framesdir.mkdir();plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':'#193245','axes.labelcolor':'#193245','axes.spines.top':False,'axes.spines.right':False})
 speed=a['speed_colour_limit_mm_s']
 for index,t in enumerate(TIMES):
  fig,axes=plt.subplots(2,2,figsize=(12.8,9.6),layout='constrained')
  for j,code in enumerate(['Fluent','OpenFOAM']):
   q=selected[code+' Fine'][index]
   for i,field in enumerate(['T','U']):
    im=maps(axes[i,j],q,96,field,speed,stream=field=='U');axes[i,j].set_title(('ANSYS Fluent' if code=='Fluent' else code)+(' | temperature' if field=='T' else ' | speed + instantaneous curves'),fontsize=12)
    if j==1:fig.colorbar(im,ax=axes[i,:].tolist(),label='Temperature (K)' if field=='T' else 'Speed (mm/s)',shrink=.85)
  fig.suptitle(f'Heat drives circulation | t = {t:05.2f} s',fontsize=17,weight='bold');fig.supxlabel('9,216 native cells | dt = 0.01 s | fixed scales | white curves are instantaneous streamlines, not particles',fontsize=10)
  fig.savefig(framesdir/f'frame-{index:03d}.png',dpi=100);plt.close(fig)
 # Both grids and codes at the same physical endpoint.
 fig,axes=plt.subplots(2,4,figsize=(16,9),layout='constrained',sharex=True,sharey=True)
 for j,key in enumerate(['Fluent Coarse','Fluent Fine','OpenFOAM Coarse','OpenFOAM Fine']):
  n=48 if key.endswith('Coarse') else 96;q=selected[key][-1]
  for i,f in enumerate(['T','U']):
   im=maps(axes[i,j],q,n,f,speed,stream=f=='U')
   if i==0:axes[i,j].set_title(f'{key}\n{n*n:,} native cells',fontsize=11)
   if j==3:fig.colorbar(im,ax=axes[i,:].tolist(),label='Temperature (K)' if f=='T' else 'Speed (mm/s)',shrink=.75)
 fig.suptitle('Mesh response at 20 s | identical dt = 0.01 s',fontsize=18,weight='bold');fig.supxlabel('Actual saved cell fields. Both meshes share the same square and physical conditions; two grids do not establish mesh independence.',fontsize=11);save(fig,out,'coarse-fine-at20s')
 # Actual field/heat histories and mesh response.
 fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
 for key in ['Fluent Coarse','Fluent Fine','OpenFOAM Coarse','OpenFOAM Fine']:
  code,grid=key.split();rows=[r for r in a['history'] if r['code']==code and r['grid']==grid];color='#16677e' if code=='Fluent' else '#b65c2b';style='--' if grid=='Coarse' else '-'
  axes[0,0].plot(TIMES,[r['Umax_m_s']*1000 for r in rows],style,c=color,label=key)
  axes[0,1].plot(TIMES[1:],[r['Nu_gradient_reconstruction'] for r in rows[1:]],style,c=color,label=key)
 for code,color in [('Fluent','#16677e'),('OpenFOAM','#b65c2b')]:
  rows=[r for r in a['mesh_responses'] if r['code']==code];axes[1,0].plot(TIMES[1:],[r['velocity_relative_L2']*100 for r in rows[1:]],c=color,label=code);axes[1,1].plot(TIMES,[r['max_temperature_difference_K']*1000 for r in rows],c=color,label=code)
 axes[0,0].set(title='Native maximum speed',ylabel='Speed (mm/s)');axes[0,1].set(title='Wall-gradient heat reconstruction',ylabel='Hot-wall Nu (dimensionless)');axes[1,0].set(title='Fine field restricted onto matching coarse cells',ylabel='Velocity L2 response (%)');axes[1,1].set(title='Largest same-cell temperature response',ylabel='Maximum |T fine − T coarse| (mK)')
 for ax in axes.ravel():ax.set_xlabel('Physical time (s)');ax.grid(alpha=.2);ax.legend(fontsize=9)
 fig.suptitle('Two grids at matched physical times',fontsize=18,weight='bold');fig.supxlabel('Heat curves reconstruct gradients from native T; checked against native face heat at endpoints. No t=0 heat, order or accuracy claim.',fontsize=10);save(fig,out,'mesh-response-histories')
 # Native residual evidence and closed mass from the selected fine states.
 fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained');f=json.loads((S/'fluent-combined-2000-step-checks-v001.json').read_text());o=json.loads((S/'openfoam-combined-2000-step-checks-v001.json').read_text());assert len(f)==len(o)==2000
 for i,name in enumerate(['continuity','Ux','Uy','energy']):axes[0,0].semilogy([r['time_s'] for r in f],[max(r['residuals'][i]/([1e-6]*3+[1e-9])[i],1e-12) for r in f],label=name,lw=.8)
 for name,lim in [('Ux',1e-6),('Uy',1e-6),('h',1e-8),('p_rgh',1e-6)]:axes[0,1].semilogy([r['time_s'] for r in o],[max(r['residuals'][name]/lim,1e-12) for r in o],label=name,lw=.8)
 for ax,title in zip(axes[0],['Fluent | all 2,000 terminal scaled residuals','OpenFOAM | all 2,000 final-outer initial residuals']):ax.axhline(1,c='black',ls='--',lw=1);ax.set(title=title,ylabel='Native residual / declared limit');ax.legend(fontsize=8,ncol=2)
 for code,color in [('Fluent','#16677e'),('OpenFOAM','#b65c2b')]:
  rr=[r for r in a['history'] if r['code']==code and r['grid']=='Fine'];axes[1,0].plot(TIMES,[r['mass_relative_drift']*1e6 for r in rr],c=color,label=code)
  checks=[r for r in a['native_heat_reconstruction_checks'] if r['code']==code and r['grid']=='Fine' and r['wall']=='left'];axes[1,1].plot([r['time_s'] for r in checks],[r['native_face_heat_W_m'] for r in checks],'o-',c=color,label=code)
 axes[1,0].set(title='Native closed mass at matched saved states',ylabel='Absolute relative drift (ppm)');axes[1,1].set(title='Independent native face heat endpoints',ylabel='Hot-wall heat per depth (W/m)')
 for ax in axes.ravel():ax.set_xlabel('Physical time (s)');ax.grid(alpha=.2)
 for ax in axes[1]:ax.legend()
 fig.suptitle('Numerical checks for the accepted finite trajectories',fontsize=18,weight='bold');fig.supxlabel('Residual norms differ. Endpoint convergence and small mass drift do not establish experimental accuracy or full transient energy closure.',fontsize=10);save(fig,out,'full-numerical-checks')
 # Machine-readable source-bound outputs; raw native solver data remain external.
 with (out/'field-and-heat-history.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(a['history'][0]));w.writeheader();w.writerows(a['history'])
 with (out/'mesh-response.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(a['mesh_responses'][0]));w.writeheader();w.writerows(a['mesh_responses'])
 (out/'selected-native-analysis.json').write_text(json.dumps(a,indent=2)+'\n')
 images=[framesdir/f'frame-{i:03d}.png' for i in range(41)]+list(out.glob('*.png'))
 for p in images:
  with Image.open(p) as im:im.load()
 receipt=dict(all_saved_field_checks_pass=True,preview=False,rendering=dict(frames=41,frame_duration_ms=250),plot=dict(times_s=TIMES.tolist(),temperature_range_K=[299.6,300.4],speed_range_mm_s=[0,speed]),rendered_png_sha256={str(p):sha(p) for p in images},analysis_sha256=sha(out/'selected-native-analysis.json'),selected_full_audit_sha256={str(p):sha(p) for p in [F/'production-audit-v001.json',O/'production-audit-v001.json']},script_sha256=sha(Path(__file__)),interpolated_motion=False,actual_visual_inspection_pending=True)
 (out/'native-field-audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(frames=41,final_mesh_responses=[x for x in a['mesh_responses'] if x['time_s']==20],max_heat_reconstruction_relative_difference=max(x['relative_difference'] for x in a['native_heat_reconstruction_checks'])),indent=2))
if __name__=='__main__':
 assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30
 q,a=collect();render(q,a,R/'exports/figures/2026-11-08-convection-spatial-full-v001')