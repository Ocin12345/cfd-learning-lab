from pathlib import Path
import json,hashlib,csv,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from sloshing_dt_common_v001 import *
from sloshing_full_audits_v002 import ffields,fclock
from audit_sloshing_dt_qualification_v003 import field
from sloshing_depth_compare_v001 import depth_comparison
from sloshing_courant_v001 import courant_records
from sloshing_cycles_v001 import cycle_evidence
T=np.arange(121,dtype=float)*.01
FIG=R/'exports/figures/2026-11-08-sloshing-dt-full-v001'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def verify(m):
 for name,h in m.items():assert sha(name)==h,name

def load_audits():
 a={};selection={}
 for code in ['fluent','openfoam']:
  for kind in ['base','new']:
   p=(R/'local/sloshing-refinement-v001'/f'{code}-production-audit.json') if kind=='base' else N/(code+'-production-v001')/'production-audit-v001.json'
   q=json.loads(p.read_text());assert q.get('numerical_pass',q.get('numerical_checks_pass',False));verify(q['source_sha256'] if 'source_sha256' in q else q['native_source_sha256'])
   assert q.get('merged_steps',q.get('steps'))==(4800 if kind=='new' else 2400)
   a[(code,kind)]=q;selection[code+' '+kind]=dict(native=str(N/(code+'-production-v001') if kind=='new' else B/code),dt_s=.00025 if kind=='new' else .0005,audit=str(p),audit_sha256=sha(p))
  p=N/(code+'-qualification-v001')/'qualification-audit-v001.json';q=json.loads(p.read_text());assert q['numerical_pass'] and q['steps']==80;verify(q['source_sha256']);a[(code,'prefix')]=q
  selection[code+' new']['prefix_native']=str(p.parent);selection[code+' new']['prefix_audit']=str(p);selection[code+' new']['prefix_audit_sha256']=sha(p)
 verify(json.loads((S/'baseline-native-source-hashes-v001.json').read_text()))
 return a,selection

def oclock(path,step,dt):
 import re
 raw=path.read_text()
 for key,val in [('value',step*dt),('index',step),('deltaT',dt)]:
  m=re.search(r'^\s*'+key+r'\s+([-+0-9.eE]+);',raw,re.M);assert m and abs(float(m[1])-val)<(1e-14 if key=='deltaT' else 1e-10)

def collect():
 audits,selection=load_audits();frames={};bindings={};hist=[];depths={}
 for code in ['fluent','openfoam']:
  for kind in ['base','new']:
   key=code+' '+kind;qs=[];ds=[];dt=.0005 if kind=='base' else .00025
   for i,t in enumerate(T):
    step=round(t/dt);case=B/code if kind=='base' else N/(code+('-qualification-v001' if i<=2 else '-production-v001'))
    if code=='fluent':
     data=case/(f'tank-{step:04d}.dat.h5' if kind=='base' else f'state-{step:04d}.dat.h5');mesh=case/('tank.cas.h5' if kind=='base' else 'active.cas.h5');q=ffields(mesh,data,case/'geometry.npz');c=fclock(data);assert c['step']==step and abs(c['time_s']-t)<1e-10 and abs(c['dt_s']-dt)<1e-14;sources=[mesh,data,case/'geometry.npz']
    else:
     centre=field(case/'0/C',3);ij=np.rint(centre[:,:2]/DX-.5).astype(int);ids=ij[:,0]+NX*ij[:,1];assert len(np.unique(ids))==NX*NY and np.max(abs(centre[:,:2]-(ij+.5)*DX))<1e-10 and np.max(abs(centre[:,2]-.005))<1e-10
     folder=case/f'{t:.8g}';q={}
     for k,name,c in [('a','alpha.water',1),('p','p_rgh',1),('U','U',3)]:
      raw=field(folder/name,c);x=np.empty_like(raw);x[ids]=raw;q[k]=x
     q['a']=q['a'].ravel();q['p']=q['p'].ravel();q['u']=q['U'][:,0];q['v']=q['U'][:,1];del q['U'];sources=[folder/name for name in ['alpha.water','U','p_rgh']]+[case/'0/C']
     if step:oclock(folder/'uniform/time',step,dt);sources.append(folder/'uniform/time')
    assert all(np.isfinite(x).all() for x in q.values());a=q['a'];area=float(a.sum()*DX*DX);assert a.min()>=-1e-6 and a.max()<=1+1e-6 and abs(area/.01-1)<=.005
    depth=a.reshape(NY,NX).sum(axis=0)*DX;ds.append(depth);speed=np.hypot(q['u'],q['v']);water=q['a']>=.5
    hist.append(dict(code=code,kind=kind,dt_s=dt,time_s=float(t),water_area_m2=area,relative_water_area_drift=abs(area/.01-1),alpha_min=float(a.min()),alpha_max=float(a.max()),max_all_mixture_speed_m_s=float(speed.max()),max_water_dominant_mixture_speed_m_s=float(speed[water].max()),right_native_column_depth_m=float(depth[-1]),right_matched_strip_depth_m=float(depth[-2:].mean()),right_native_column_centre_m=.199375,right_strip_centre_m=.19875,colour_saturation_count=int(((a<0)|(a>1)).sum())))
    bindings[key+f'@{t:.2f}s']=dict(case=str(case),global_step=int(step),native_time_s=float(t),source_sha256={str(p):sha(p) for p in sources})
    qs.append({k:v.reshape(NY,NX) for k,v in q.items() if k!='p'})
   frames[key]=qs;depths[key]=np.array(ds)
 comparisons={}
 for code in ['fluent','openfoam']:
  assert np.array_equal(depths[code+' base'][0],depths[code+' new'][0]);comparisons[code]=depth_comparison(T,depths[code+' base'],depths[code+' new'],.05)
 code_response=depths['fluent new']-depths['openfoam new']
 summary=dict(times_s=T.tolist(),source_bindings=bindings,selected_cases=selection,history=hist,temporal_comparisons=comparisons,small_dt_cross_code_depth_difference=dict(maximum_absolute_m=float(abs(code_response).max()),sample_RMS_m=float(np.sqrt(np.mean(code_response**2))),endpoint_RMS_m=float(np.sqrt(np.mean(code_response[-1]**2)))),global_water_dominant_peak_speed_m_s=max(r['max_water_dominant_mixture_speed_m_s'] for r in hist),global_all_mixture_peak_speed_m_s=max(r['max_all_mixture_speed_m_s'] for r in hist),physical_validation=False,time_independence=False,phase_alignment_applied=False,full_trajectory_workflow='new dt uses80-step qualification plus exact4720-step native checkpoint continuation')
 return frames,depths,summary,audits

def save(fig,name):
 for ext in ['png','pdf']:fig.savefig(FIG/(name+'.'+ext),dpi=140)
 plt.close(fig)

def render(frames,depths,a,audits):
 assert not FIG.exists();FIG.mkdir(parents=True);out=FIG/'frames';out.mkdir();plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'#193245','axes.labelcolor':'#193245','axes.spines.top':False,'axes.spines.right':False})
 x=(np.arange(NX)+.5)*DX*1000;xe=np.arange(NX+1)*DX*1000;ye=np.arange(NY+1)*DX*1000;cc=np.meshgrid(x,(np.arange(NY)+.5)*DX*1000)
 waterpeak=a['global_water_dominant_peak_speed_m_s'];colors={'fluent':'#16677e','openfoam':'#b65c2b'};values=[np.array(v) for z in a['temporal_comparisons'].values() for v in [z['right_strip_base_depth_m'],z['right_strip_new_depth_m']]];lo=min(v.min() for v in values)*1000;hi=max(v.max() for v in values)*1000;pad=(hi-lo)*.1
 delta=max(max(abs(np.array(v['right_strip_difference_m']))) for v in a['temporal_comparisons'].values())*1e6*1.15
 for i,t in enumerate(T):
  fig=plt.figure(figsize=(14,10),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1,.9]);tankaxs=[]
  for row,kind in enumerate(['base','new']):
   for col,code in enumerate(['fluent','openfoam']):
    ax=fig.add_subplot(gs[row,col]);tankaxs.append(ax);q=frames[code+' '+kind][i];im=ax.pcolormesh(xe,ye,q['a'],cmap='Blues',vmin=0,vmax=1,shading='flat',rasterized=True)
    ax.plot(x,depths[code+' '+kind][i]*1000,c='#dc8a2f',lw=.7)
    sl=np.s_[5::12,6::16];u=q['u'][sl]*1000;v=q['v'][sl]*1000;mask=(q['a'][sl]<.5)|(np.hypot(u,v)<.005*waterpeak*1000)
    ax.quiver(cc[0][sl],cc[1][sl],np.ma.array(u,mask=mask),np.ma.array(v,mask=mask),color='white',angles='xy',scale_units='xy',scale=waterpeak*1000/8,width=.003)
    ax.set(aspect='equal',xlim=(0,200),ylim=(0,100),xlabel='x (mm)',ylabel='y (mm)',title=('ANSYS Fluent' if code=='fluent' else 'OpenFOAM')+(' | dt = 0.50 ms' if kind=='base' else ' | dt = 0.25 ms'))
  fig.colorbar(im,ax=tankaxs,label='Water volume fraction',shrink=.78)
  p=fig.add_subplot(gs[2,0]);d=fig.add_subplot(gs[2,1])
  for code in ['fluent','openfoam']:
   z=a['temporal_comparisons'][code]
   for kind,ls in [('base','--'),('new','-')]:
    val=np.array(z['right_strip_base_depth_m' if kind=='base' else 'right_strip_new_depth_m'])*1000;p.plot(T[:i+1],val[:i+1],ls,c=colors[code],label=code.title()+(' dt' if kind=='base' else ' dt/2'));p.plot(t,val[i],'o',ms=3,c=colors[code])
   val=np.array(z['right_strip_difference_m'])*1e6;d.plot(T[:i+1],val[:i+1],c=colors[code],label=code.title());d.plot(t,val[i],'o',ms=3,c=colors[code])
  p.set(xlim=(0,1.2),ylim=(lo-pad,hi+pad),xlabel='Physical time (s)',ylabel='Water depth (mm)',title='Right physical strip: x = 197.5–200 mm');d.set(xlim=(0,1.2),ylim=(-delta,delta),xlabel='Physical time (s)',ylabel='dt/2 minus dt depth (µm)',title='Same-time step-size response, no phase alignment')
  for ax in [p,d]:ax.grid(alpha=.2);ax.legend(fontsize=8,ncol=2);ax.axvline(t,c='gray',lw=.6,ls=':')
  fig.suptitle(f'Sloshing time-step comparison | t = {t:.3f} s',fontsize=17,weight='bold');fig.supxlabel('Same 12,800 cells | orange curve: column-integrated water depth | white arrows: instantaneous mixture velocity in α ≥ 0.5 cells',fontsize=10)
  fig.savefig(out/f'frame-{i:03d}.png',dpi=100);plt.close(fig)
 # Native endpoint cell-speed maps, with an explicit water-dominant display mask.
 fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained');limit=np.ceil(waterpeak*10000)/10
 for i,kind in enumerate(['base','new']):
  for j,code in enumerate(['fluent','openfoam']):
   q=frames[code+' '+kind][-1];s=np.ma.array(np.hypot(q['u'],q['v'])*1000,mask=q['a']<.5);im=axs[i,j].pcolormesh(xe,ye,s,cmap='viridis',vmin=0,vmax=limit,rasterized=True);axs[i,j].plot(x,depths[code+' '+kind][-1]*1000,c='#dc8a2f',lw=.7);axs[i,j].set(aspect='equal',xlabel='x (mm)',ylabel='y (mm)',title=code.title()+(' | dt = 0.50 ms' if kind=='base' else ' | dt = 0.25 ms'))
 fig.colorbar(im,ax=axs.ravel().tolist(),label='Mixture speed where α ≥ 0.5 (mm/s)',shrink=.8);fig.suptitle('Actual native fields at 1.2 s | identical grid and fixed scale',fontsize=17,weight='bold');fig.supxlabel('Air-dominant cells masked for this display only. Orange curves integrate column water volume; they are not traced interface points.',fontsize=10);save(fig,'native-endpoint-speed')
 # Quantitative depth/phase evidence.
 fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained')
 for code in ['fluent','openfoam']:
  z=a['temporal_comparisons'][code];colour=colors[code]
  for kind,ls in [('base','--'),('new','-')]:axs[0,0].plot(T,np.array(z['right_strip_base_depth_m' if kind=='base' else 'right_strip_new_depth_m'])*1000,ls,c=colour,label=code.title()+(' dt' if kind=='base' else ' dt/2'))
  axs[0,1].plot(T,np.array(z['right_strip_difference_m'])*1e6,c=colour,label=code.title());deltafield=depths[code+' new']-depths[code+' base'];axs[1,0].plot(T,np.max(abs(deltafield),axis=1)*1000,c=colour,label=code.title());axs[1,1].plot(x,deltafield[-1]*1000,c=colour,label=code.title())
 axs[0,0].set(title='Matched right-strip water depth',ylabel='Depth (mm)',xlabel='Physical time (s)');axs[0,1].set(title='Right-strip step-size response',ylabel='dt/2 minus dt (µm)',xlabel='Physical time (s)');axs[1,0].set(title='Largest native-column response at each saved time',ylabel='Maximum absolute depth difference (mm)',xlabel='Physical time (s)');axs[1,1].set(title='Signed endpoint depth response at 1.2 s',ylabel='dt/2 minus dt (mm)',xlabel='Native fine-column centre x (mm)')
 for ax in axs.ravel():ax.grid(alpha=.2);ax.legend(fontsize=9)
 fig.suptitle('Temporal response at matched physical positions and times',fontsize=17,weight='bold');fig.supxlabel('No temporal alignment. Crossing estimates use linearly bracketed 0.01 s samples; two steps do not establish accuracy or time independence.',fontsize=10);save(fig,'depth-and-phase-evidence')
 # Original residual norms, water-area and native pre-advance Courant records.
 fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained')
 for code,j in [('fluent',0),('openfoam',1)]:
  rows=audits[(code,'prefix')]['step_records']+audits[(code,'new')]['step_records'];assert len(rows)==4800 and [r['step'] for r in rows]==list(range(1,4801))
  for k in rows[0]['residuals']:axs[0,j].semilogy([r['time_s'] for r in rows],[max(r['residuals'][k],1e-20) for r in rows],label=k,lw=.7)
  axs[0,j].axhline(1e-4,c='black',ls='--',lw=1);axs[0,j].set(title=code.title()+(' terminal scaled' if code=='fluent' else ' final-outer INITIAL')+' residuals',ylabel='Native residual',xlabel='Physical time (s)');axs[0,j].legend(fontsize=8,ncol=2)
  h=[r for r in a['history'] if r['code']==code and r['kind']=='new'];axs[1,0].plot(T,[r['relative_water_area_drift']*1e6 for r in h],c=colors[code],label=code.title())
  co=audits[(code,'prefix')]['courant_records']+audits[(code,'new')]['courant_records'];axs[1,1].plot([r['actual_pre_advance_time_s'] for r in co],[r['maximum_native_interface_courant'] for r in co],c=colors[code],label=code.title(),lw=.8)
 axs[1,0].set(title='Native water-area conservation (limit 5,000 ppm)',ylabel='Absolute relative drift (ppm)',xlabel='Physical time (s)');axs[1,1].set(title='Native VOF/interface Courant, including restart reports',ylabel='Maximum Courant number',xlabel='Actual pre-advance physical time (s)');axs[1,1].axhline(.25,c='black',ls='--',lw=1)
 for ax in axs.ravel():ax.grid(alpha=.2)
 for ax in axs[1]:ax.legend(fontsize=9)
 fig.suptitle('Accepted 4,800-step trajectories | native numerical evidence',fontsize=17,weight='bold');fig.supxlabel('Residual norms differ. Saved fraction bounds, clocks and wall impermeability pass; numerical convergence is not physical validation.',fontsize=10);save(fig,'full-numerical-evidence')
 with (FIG/'field-history.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(a['history'][0]));w.writeheader();w.writerows(a['history'])
 with (FIG/'native-column-depths.csv').open('w',newline='') as f:
  w=csv.writer(f);w.writerow(['case','time_s','column_centre_m','depth_m'])
  for key,d in depths.items():
   for i,t in enumerate(T):
    for j in range(NX):w.writerow([key,float(t),float((j+.5)*DX),float(d[i,j])])
 (FIG/'selected-native-analysis.json').write_text(json.dumps(a,indent=2,allow_nan=False)+'\n');images=list(out.glob('*.png'))+list(FIG.glob('*.png'));assert len(list(out.glob('*.png')))==121
 for p in images:
  with Image.open(p) as im:im.load()
 receipt=dict(all_saved_field_checks_pass=True,preview=False,rendering=dict(frames=121,frame_duration_ms=1000/12),plot=dict(times_s=T.tolist(),alpha_range=[0,1],arrow_common_water_peak_m_s=waterpeak),rendered_png_sha256={str(p):sha(p) for p in images},analysis_sha256=sha(FIG/'selected-native-analysis.json'),selection_sha256=sha(S/'accepted-full-pair-selection-v001.json'),script_sha256=sha(Path(__file__)),actual_visual_inspection_pending=True,interpolated_motion=False)
 (FIG/'native-field-audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({code:{k:v for k,v in z.items() if not isinstance(v,(list,dict))} for code,z in a['temporal_comparisons'].items()},indent=2))

def main():
 from wing_fluent_v001 import available_memory
 assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30
 frames,depths,a,audits=collect();p=S/'accepted-full-pair-selection-v001.json';assert not p.exists();p.write_text(json.dumps(dict(cases=a['selected_cases'],matched_times_s=T.tolist(),native_analysis_method_sha256=sha(Path(__file__)),both_full_audits_pass=True,physical_validation=False),indent=2)+'\n');render(frames,depths,a,audits)
if __name__=='__main__':main()
