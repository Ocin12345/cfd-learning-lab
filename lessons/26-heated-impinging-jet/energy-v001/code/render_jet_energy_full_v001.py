"""Actual accepted native temperatures/wall heat and complete discrete energy plots."""
import json,csv,re,shutil
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from jet_energy_common_v001 import *
from jet_energy_full_prelaunch_v001 import wait_admission
from audit_jet_energy_full_v001 import locate
from analyse_heated_jet_v002 import values,partial,common_mapping,chunk_array,field,counted,TDIM,NF
from jet_common_v001 import NX,NY,DX,HALF,W,RHO,DEPTH
FIG=R/'exports/figures/2026-11-08-heated-jet-energy-full-v001'
COLORS={'fluent':'#1768a6','openfoam':'#cf692d'};LABELS={'fluent':'ANSYS Fluent','openfoam':'OpenFOAM'}
def csvout(name,cols,rows):
 with (FIG/name).open('x',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(cols);w.writerows(rows)
def main():
 assert not FIG.exists();wait_admission();audits={};selection={}
 for code in ['fluent','openfoam']:
  p=N/(code+'-production-v001')/'production-audit-v001.json';a=json.loads(p.read_text());assert a['numerical_pass'] and a['actual_steps']==600 and a['unique_states']==601
  for src,h in a['source_sha256'].items():assert sha(src)==h,src
  audits[code]=a;selection[code]=dict(case=str(N/(code+'-production-v001')),audit=str(p),audit_sha256=sha(p),prefix_case=str(N/(code+'-qualification-v001')),continuation_not_independent_full_rerun=True)
 dump(S/'accepted-full-pair-selection-v001.json',selection);FIG.mkdir(parents=True);(FIG/'frames').mkdir()
 with np.load(PARENT/'geometry.npz') as geo:fm=common_mapping(N/'fluent-production-v001/active.cas.h5',geo)
 of=N/'openfoam-production-v001';c=field(of/'1250/C',NX*NY,3,[0,1,0,0,0,0,0]);om=np.rint(c[:,1]/DX-.5).astype(int)*NX+np.rint((c[:,0]+HALF)/DX-.5).astype(int);assert len(np.unique(om))==NX*NY
 def align(v,m):
  a=np.empty_like(v);a[m]=v;return a.reshape(NY,NX)
 U={};V={}
 fdata=locate('fluent',600);U['fluent']=align(values(fdata,'SV_U'),fm);V['fluent']=align(values(fdata,'SV_V'),fm)
 ov=field(locate('openfoam',600).parent/'U',NX*NY,3,[0,1,-1,0,0,0,0]);U['openfoam']=align(ov[:,0],om);V['openfoam']=align(ov[:,1],om)
 temps={code:{} for code in audits};native_binding={code:[] for code in audits};steps=list(range(0,601,10));times=[s*.005 for s in steps]
 for code,a in audits.items():
  for s in steps:
   p=locate(code,s);assert a['source_sha256'][str(p)]==sha(p)
   v=values(p,'SV_T') if code=='fluent' else field(p,NX*NY,1,TDIM);assert np.isfinite(v).all() and v.min()>=299.99 and v.max()<=320.01
   temps[code][s]=align(v,fm if code=='fluent' else om);native_binding[code].append(dict(step=s,time_s=s*.005,path=str(p),sha256=sha(p)))
 # Native plate coordinate/face mappings, unchanged geometry.
 with h5py.File(N/'fluent-production-v001/active.cas.h5') as h:
  b=h['meshes/1/faces'];z=b['zoneTopology'];zones={n:slice(int(lo)-1,int(hi)) for n,lo,hi in zip(z['name'][0].decode().split(';'),z['minId'][...],z['maxId'][...])};pts=chunk_array(h['meshes/1/nodes/coords'],(NX+1)*(NY+1),(2,));fnode=b['nodes/1/nodes'][...].reshape(-1,2).astype(int)-1;xf=pts[fnode].mean(axis=1)[zones['plate'],0];fo=np.argsort(xf)
 n,b=counted(of/'constant/polyMesh/boundary');pb={k:dict(count=int(re.search(r'nFaces\s+(\d+)',v)[1]),start=int(re.search(r'startFace\s+(\d+)',v)[1])) for k,v in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',b)};_,b=counted(of/'constant/polyMesh/owner');owner=np.fromstring(b,sep=' ',dtype=int);pa=pb['plate'];op=owner[pa['start']:pa['start']+pa['count']];oo=np.argsort(c[op,0]);assert np.allclose(xf[fo],c[op,0][oo],rtol=0,atol=1e-12)
 with h5py.File(locate('fluent',600)) as h:fq=partial(h['results/1/phase-1/faces/SV_HEAT_FLUX'],NF)[zones['plate']][fo]/DX
 ot=field(locate('openfoam',600),NX*NY,1,TDIM);oq=(.0262*(320-ot[op])/(DX/2))[oo]
 for code,q in [('fluent',fq),('openfoam',oq)]:assert np.isclose(q.sum()*DX,audits[code]['patch_input_W_per_m']['plate'][-1],rtol=1e-10,atol=1e-10)
 csvout('plate-heat-flux-3s.csv',['x_over_W','Fluent_native_wall_q_W_m2','OpenFOAM_scheme_wall_q_W_m2'],zip(xf[fo]/W,fq,oq))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#203648','text.color':'#203648','xtick.color':'#203648','ytick.color':'#203648'})
 x=(np.arange(NX)+.5)*DX-HALF;y=(np.arange(NY)+.5)*DX
 fig=plt.figure(figsize=(12,7),dpi=100,facecolor='white');gs=fig.add_gridspec(2,2,width_ratios=[1,.025],left=.075,right=.93,top=.83,bottom=.15,hspace=.36,wspace=.04);ims=[]
 for i,code in enumerate(['fluent','openfoam']):
  ax=fig.add_subplot(gs[i,0]);im=ax.imshow(temps[code][0],origin='lower',extent=[-6,6,0,3],cmap='inferno',vmin=300,vmax=320,interpolation='nearest',aspect='equal');ims.append(im)
  ax.streamplot(x/W,y/W,U[code],V[code],color='#97a7b1',density=(1.,.65),linewidth=.55,arrowsize=.6);ax.set(title=LABELS[code],ylabel='$y/W$',xlim=(-6,6),ylim=(0,3));ax.plot([-6,6],[0,0],color='#ee783d',lw=2);ax.plot([-.5,.5],[3,3],color='#1eb9c0',lw=2)
  if i:ax.set_xlabel('$x/W$')
  else:ax.tick_params(labelbottom=False)
 fig.colorbar(ims[0],cax=fig.add_subplot(gs[:,1])).set_label('Air temperature (K)');fig.text(.075,.945,'Cold jet over a prescribed hot plate',size=19,weight='bold');fig.text(.075,.895,'Re = 100  |  300 K inlet  |  320 K plate  |  temperature evolves in frozen airflow',size=11)
 clock=fig.text(.075,.065,'Actual saved fields  |  t = 0.00 s',size=12);fig.text(.075,.023,'Grey: streamlines of each code’s own fixed velocity field. Plate temperature is prescribed; it does not cool.',size=9)
 for index,s in enumerate(steps):
  for im,code in zip(ims,['fluent','openfoam']):im.set_data(temps[code][s])
  clock.set_text(f'Actual saved fields  |  t = {s*.005:.2f} s');fig.savefig(FIG/'frames'/f'frame-{index:03d}.png',dpi=100)
  if s==600:fig.savefig(FIG/'temperature-comparison-3s.png',dpi=170);fig.savefig(FIG/'temperature-comparison-3s.pdf')
 plt.close(fig)
 # Common-time wall heat and actual endpoint profile, no accuracy ranking.
 fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
 for code in audits:
  a=audits[code];t=np.array(a['energy_ledger']['time_s']);axes[0].plot(t,a['patch_input_W_per_m']['plate'],color=COLORS[code],label=LABELS[code])
 axes[1].plot(xf[fo]/W,fq,color=COLORS['fluent'],label='ANSYS Fluent: native wall heat');axes[1].plot(xf[fo]/W,oq,color=COLORS['openfoam'],ls='--',label='OpenFOAM: scheme reconstruction');axes[0].set(title='Heat supplied by the fixed-temperature plate',xlabel='Time (s)',ylabel='Heat rate into air (W m$^{-1}$)',xlim=(0,3));axes[1].set(title='Local plate heat flux at t = 3 s',xlabel='$x/W$',ylabel='Heat flux into air (W m$^{-2}$)',xlim=(-6,6))
 for ax in axes:ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8)
 fig.suptitle('Thermal evolution in separately frozen flow fields; prescribed plate remains at 320 K',size=11);fig.savefig(FIG/'wall-heat-transfer.png',dpi=160);fig.savefig(FIG/'wall-heat-transfer.pdf');plt.close(fig)
 fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
 for col,code in enumerate(['fluent','openfoam']):
  a=audits[code];l=a['energy_ledger'];t=np.array(l['time_s']);plate=np.cumsum(.005*np.array(a['patch_input_W_per_m']['plate']));other=np.cumsum(.005*sum(np.array(v) for k,v in a['patch_input_W_per_m'].items() if k!='plate'));ax=axes[0,col]
  for v,label,color,ls in [(plate,'Plate input','#a94f19','-'),(other,'Other boundaries: signed input','#68958a','-'),(l['cumulative_net_input_J_per_m'],'Net boundary input',COLORS[code],'-'),(l['stored_energy_change_J_per_m'],'Stored thermal energy','#2b3642','--')]:ax.plot(t,v,label=label,color=color,ls=ls)
  ax.set(title=LABELS[code],ylabel='Energy relative to 300 K (J m$^{-1}$)');ax.legend(frameon=False,fontsize=8)
  ax=axes[1,col];ax.semilogy(t,np.maximum(l['step_relative_error'],1e-18),label='Each-step mismatch',color=COLORS[code]);ax.semilogy(t,np.maximum(l['cumulative_relative_error'],1e-18),label='Cumulative prefix mismatch',color='#35424f',ls='--');ax.axhline(.005,color='#777777',ls=':',label='Unchanged 0.5% criterion');ax.set(ylabel='Normalized energy mismatch');ax.legend(frameon=False,fontsize=8)
  for ax in axes[:,col]:ax.set(xlabel='Physical time (s)',xlim=(0,3));ax.grid(alpha=.15)
 fig.suptitle('Independent boundary inputs and thermal storage over all 600 actual steps\nDiscrete backward-Euler accounting; OpenFOAM boundary heat is reconstructed from its native fields',size=11);fig.savefig(FIG/'complete-energy-ledger.png',dpi=160);fig.savefig(FIG/'complete-energy-ledger.pdf');plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
 for ax,code in zip(axes,['fluent','openfoam']):
  a=audits[code];ax.semilogy([r['time_s'] for r in a['step_residuals']],[r['residual'] for r in a['step_residuals']],color=COLORS[code]);ax.axhline(a['thermal_residual_limit'],color='#777777',ls=':',label='Original criterion');ax.set(title=LABELS[code],xlabel='Physical time (s)',ylabel='Native terminal thermal residual',xlim=(0,3));ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8)
 fig.suptitle('All600 steps checked: Fluent scaled energy residual / OpenFOAM final linear T residual\nThe two residual norms differ; their magnitudes do not rank physical accuracy.',size=11);fig.savefig(FIG/'thermal-residuals.png',dpi=160);fig.savefig(FIG/'thermal-residuals.pdf');plt.close(fig)
 for code,a in audits.items():
  l=a['energy_ledger'];csvout(code+'-complete-energy-history.csv',['step','time_s','plate_input_W_m','net_input_W_m','storage_rate_W_m','cumulative_net_input_J_m','stored_energy_change_J_m','step_error_relative','prefix_error_relative','terminal_residual'],[(i+1,l['time_s'][i],a['patch_input_W_per_m']['plate'][i],l['net_input_W_per_m'][i],l['storage_rate_W_per_m'][i],l['cumulative_net_input_J_per_m'][i],l['stored_energy_change_J_per_m'][i],l['step_relative_error'][i],l['cumulative_relative_error'][i],a['step_residuals'][i]['residual']) for i in range(600)])
 pngs={str(p):sha(p) for p in (FIG/'frames').glob('*.png')};assert len(pngs)==61
 for p in list((FIG/'frames').glob('*.png'))+list(FIG.glob('*.png')):
  with Image.open(p) as im:im.load();assert im.width>=1000 and im.height>=600
 summary={code:dict(total_plate_input_J_per_m=float(.005*sum(a['patch_input_W_per_m']['plate'])),net_boundary_input_J_per_m=a['energy_ledger']['cumulative_net_input_J_per_m'][-1],stored_energy_change_J_per_m=a['energy_ledger']['stored_energy_change_J_per_m'][-1],max_step_relative=max(a['energy_ledger']['step_relative_error']),max_prefix_relative=max(a['energy_ledger']['cumulative_relative_error']),max_original_T_replay_difference_K=max(a['original_state_T_max_abs_replay_difference_K'].values()),final_plate_heat_W_per_m=a['patch_input_W_per_m']['plate'][-1],max_mass_relative=a['max_mass_relative'],max_terminal_residual=a['max_terminal_residual']) for code,a in audits.items()};dump(FIG/'energy-summary.json',summary)
 dump(FIG/'native-field-audit.json',dict(all_saved_field_checks_pass=True,preview=False,selection=selection,selection_sha256=sha(S/'accepted-full-pair-selection-v001.json'),native_frame_binding=native_binding,rendered_png_sha256=pngs,static_png_sha256={str(p):sha(p) for p in FIG.glob('*.png')},rendering=dict(frames=61,frame_duration_ms=100),plot=dict(times_s=times,temperature_limits_K=[300,320],streamlines='Each code own frozen native U; static postprocessing, not transient tracers'),energy_summary=summary,interpolated_transient_frames=False,physical_validation=False))
 print(json.dumps(dict(figures=str(FIG),actual_frames=61,summary=summary)))
if __name__=='__main__':main()
