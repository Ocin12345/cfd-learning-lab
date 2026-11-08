"""Render only source-bound accepted Fluent/OpenFOAM production fields."""
from pathlib import Path
import argparse,csv,hashlib,json,re,shutil
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from analyse_nozzle_fluent_v005 import common_mapping
from dam_break_common_v001 import *
R=Path('E:/CFD-Learning-Lab');P=R/'local/dam-break-over-obstacle-v001'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def water_group(h):
 groups=[p['cells/SV_VOF'] for p in h['results/1'].values() if 'cells/SV_VOF' in p]
 if len(groups)==1:return groups[0]
 groups=[v for v in groups if '/phase-3/' in v.name];assert len(groups)==1;return groups[0]
def main():
 ap=argparse.ArgumentParser()
 for name in ['fluent-case','openfoam-case','fluent-audit','openfoam-audit','out']:ap.add_argument('--'+name,type=Path,required=True)
 ap.add_argument('--check-only',action='store_true');args=ap.parse_args()
 F=args.fluent_case.resolve();O=args.openfoam_case.resolve();OUT=args.out.resolve()
 assert F.is_relative_to((R/'native/runs').resolve()) and O.is_relative_to((R/'native/runs').resolve())
 assert OUT.is_relative_to((R/'exports/figures').resolve()) and not OUT.exists()
 from wing_fluent_v001 import available_memory
 assert available_memory()/2**30>=4 and shutil.disk_usage(R).free>=10*2**30
 audit_paths=[args.fluent_audit.resolve(),args.openfoam_audit.resolve()]
 selection_path=P/'accepted-full-pair-selection-v001.json';selection=json.loads(selection_path.read_text());assert selection['full_pair_numerical_checks_pass']
 assert Path(selection['fluent_case']).resolve()==F and Path(selection['openfoam_case']).resolve()==O,'Selected case does not match the accepted pair manifest'
 assert Path(selection['fluent_audit']).resolve()==audit_paths[0] and Path(selection['openfoam_audit']).resolve()==audit_paths[1]
 for p in audit_paths:assert selection['audit_sha256'][str(p)]==sha(p)
 audits=[json.loads(p.read_text()) for p in audit_paths]
 for a in audits:
  assert a['numerical_checks_pass'] and a['production_complete'] and a['mode']=='production' and a['cells']==NCELLS and a['saved_states']==101, 'Full paired production audits must pass before rendering'
  assert len(a['records'])==101 and len(a['step_records'])==a['steps']
  assert abs(a['final_time_s']-.5)<1e-10
  for name,h in a['native_source_sha256'].items():assert sha(Path(name))==h,name
 times=np.array([r['time_s'] for r in audits[0]['records']])
 assert np.allclose(times,np.arange(101)*.005,atol=1e-10,rtol=0)
 assert np.allclose(times,[r['time_s'] for r in audits[1]['records']],atol=1e-10,rtol=0)
 proofs=[{str(Path(k).resolve()):v for k,v in a['native_source_sha256'].items()} for a in audits];sources=[]
 def proven(p,j):
  assert proofs[j].get(str(p.resolve()))==sha(p),'Case is not bound to its selected audit: '+str(p)
  sources.append(p)
 proven(F/'dam.cas.h5',0);proven(O/'0/C',1);proven(O/'0/Vc',1);proven(O/'solver.log',1)
 if args.check_only:
  print(json.dumps(dict(full_pair_checks_pass=True,selected_case_binding_pass=True,rendered_frames=0)));return
 g=geometry();fm=common_mapping(F/'dam.cas.h5',g);om=g['openfoam_to_common']
 c=native_field(O/'0/C',3);vol=native_field(O/'0/Vc',1).ravel()
 assert np.max(abs(c[:,:2]-g['centres'][om]))<1e-12 and np.max(abs(vol/DEPTH-g['areas'][om]))<1e-15
 data=np.empty((2,101,4,NCELLS))
 for index,t in enumerate(times):
  fstep=audits[0]['records'][index]['step'];p=F/f'dam-{fstep:04d}.dat.h5';proven(p,0)
  with h5py.File(p) as h:
   for k,group in enumerate([water_group(h),h['results/1/phase-1/cells/SV_U'],h['results/1/phase-1/cells/SV_V'],h['results/1/phase-1/cells/SV_P']]):
    raw=cell_values(group);data[0,index,k,fm]=raw
  folder=O/f'{t:.12g}'
  for name in ['alpha.water','U','p_rgh']:proven(folder/name,1)
  a=native_field(folder/'alpha.water',1).ravel();v=native_field(folder/'U',3);pr=native_field(folder/'p_rgh',1).ravel()
  for k,raw in enumerate([a,v[:,0],v[:,1],pr-(1+999*a)*G*c[:,1]]):data[1,index,k,om]=raw
 assert np.isfinite(data).all()
 alpha=data[:,:,0];u=data[:,:,1];v=data[:,:,2];pressure=data[:,:,3]
 assert alpha.min()>=-1e-6 and alpha.max()<=1+1e-6
 assert np.max(abs(alpha[0,0]-alpha[1,0]))<1e-10 and np.max(abs(pressure[0,0]-pressure[1,0]))<1e-8
 speed=np.hypot(u,v);vmax=float(np.max(np.where(alpha>=.5,speed,0)));assert vmax>0
 ij=g['cell_ij'];selection=(ij[:,0]%8==0)&(ij[:,1]%8==0)
 def grid(values):
  a=np.ma.masked_all((NY,NX));a[ij[:,1],ij[:,0]]=values;return a
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
 labels=['ANSYS Fluent','OpenFOAM 14'];colors=['#176b9c','#ce652b'];styles=['-','--']
 fig,axes=plt.subplots(2,2,figsize=(12,11));fig.subplots_adjust(left=.07,right=.86,bottom=.115,top=.83,hspace=.30,wspace=.27)
 images=[[],[]];arrows=[]
 def decorate(ax):
  ax.add_patch(Rectangle((OB_LEFT*1000,0),(OB_RIGHT-OB_LEFT)*1000,OB_HEIGHT*1000,facecolor='#8d969e',edgecolor='#454b50',zorder=4))
  ax.set(aspect='equal',xlim=[0,584],ylim=[0,584],xlabel='x (mm)',ylabel='y (mm)',xticks=[0,146,292,438,584],yticks=[0,146,292,438,584])
 for j in range(2):
  for row in range(2):
   values=alpha[j,0] if row==0 else np.ma.masked_where(alpha[j,0]<.5,speed[j,0])
   im=axes[row,j].pcolormesh(g['x_edges']*1000,g['y_edges']*1000,grid(values),shading='flat',cmap='Blues' if row==0 else 'viridis',vmin=0,vmax=1 if row==0 else vmax,rasterized=True)
   images[row].append(im);decorate(axes[row,j]);axes[row,j].set_title(labels[j]+(' | water fraction' if row==0 else ' | mixture speed'))
  mask=(alpha[j,0,selection]<.5)|(speed[j,0,selection]<vmax*.005)
  arrows.append(axes[1,j].quiver(g['centres'][selection,0]*1000,g['centres'][selection,1]*1000,np.ma.masked_where(mask,u[j,0,selection]),np.ma.masked_where(mask,v[j,0,selection]),angles='xy',scale_units='xy',scale=vmax/35,color='white',width=.004,minlength=0))
 for row,label in [(0,'Water volume fraction'),(1,'Mixture speed (m/s), alpha >= 0.5')]:
  pos=axes[row,1].get_position();cb=fig.add_axes([.90,pos.y0,.015,pos.height]);fig.colorbar(images[row][0],cax=cb).set_label(label)
 fig.suptitle('Dam break over an obstacle',x=.07,ha='left',y=.98,fontsize=22,fontweight='bold')
 clock=fig.text(.07,.935,'',fontsize=14)
 dts=[a['dt_s']*1000 for a in audits]
 fig.text(.07,.895,f'9,072 cells | Fluent dt = {dts[0]:g} ms | OpenFOAM dt = {dts[1]:g} ms | static open atmosphere',fontsize=10)
 fig.text(.07,.065,'Fixed scales across codes and times. White arrows: instantaneous velocity in cells with alpha >=0.5, not particle tracks.',fontsize=9)
 fig.text(.07,.040,'101 actual native states, 0.005 s apart, played at 12 fps; no generated intermediate motion.',fontsize=9)
 fig.text(.07,.018,'Declared numerical checks pass. Mesh/time independence and experimental validation remain open.',fontsize=9)
 OUT.mkdir(parents=True)
 for index,t in enumerate(times):
  clock.set_text(f'Physical time t = {t:.3f} s')
  for j in range(2):
   images[0][j].set_array(grid(alpha[j,index]).ravel())
   images[1][j].set_array(grid(np.ma.masked_where(alpha[j,index]<.5,speed[j,index])).ravel())
   mask=(alpha[j,index,selection]<.5)|(speed[j,index,selection]<vmax*.005)
   arrows[j].set_UVC(np.ma.masked_where(mask,u[j,index,selection]),np.ma.masked_where(mask,v[j,index,selection]))
  fig.savefig(OUT/f'frame-{index:03d}.png',dpi=120)
  if index in {0,25,50,75,100}:fig.savefig(OUT/f'state-{index:03d}.pdf')
  print(json.dumps(dict(frame=index,physical_time_s=float(t))),flush=True)
 plt.close(fig)
 # Physical gauge pressure reconstructed from each solver's own native state.
 pbound=max(1e-6,float(np.max(np.where(alpha>=.5,abs(pressure),0))))/1000
 fig,ax=plt.subplots(1,2,figsize=(11,6));fig.subplots_adjust(left=.075,right=.86,bottom=.20,top=.80,wspace=.28)
 for j in range(2):
  image=ax[j].pcolormesh(g['x_edges']*1000,g['y_edges']*1000,grid(np.ma.masked_where(alpha[j,-1]<.5,pressure[j,-1]/1000)),cmap='RdBu_r',vmin=-pbound,vmax=pbound,shading='flat',rasterized=True)
  decorate(ax[j]);ax[j].set_title(labels[j])
 cb=fig.add_axes([.90,.25,.015,.45]);fig.colorbar(image,cax=cb).set_label('Physical gauge pressure (kPa), alpha >= 0.5')
 fig.suptitle('Dam break | pressure at t = 0.500 s',x=.075,ha='left',y=.96,fontsize=19,fontweight='bold')
 fig.text(.075,.10,'Common symmetric colour scale uses the maximum absolute water-dominated pressure over all 101 states in both codes.',fontsize=9)
 fig.text(.075,.045,'OpenFOAM p = p_rgh - rho*g*y. Instantaneous cell pressures do not establish resolved impact loads.',fontsize=9)
 fig.savefig(OUT/'pressure-final.png',dpi=140);fig.savefig(OUT/'pressure-final.pdf');plt.close(fig)
 # Numerical histories from actual native records, not identical residual definitions.
 fig,ax=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(left=.10,right=.96,bottom=.16,top=.82,hspace=.45,wspace=.30)
 for j,a in enumerate(audits):
  t=np.array([r['time_s'] for r in a['step_records']]);ratio=np.array([max(r['residuals'].values())/1e-4 for r in a['step_records']]);assert ratio.max()<=1
  closure=np.array([abs(r['relative_mixture_mass_closure']) for r in a['records']])*100
  bounds=np.array([max(0,-r['alpha_min'],r['alpha_max']-1) for r in a['records']])
  ax[0,0].semilogy(t,np.maximum(ratio,1e-16),styles[j],color=colors[j],lw=1,label=labels[j])
  ax[0,1].semilogy(times,np.maximum(closure,1e-16),styles[j],color=colors[j],lw=1.5)
  ax[1,0].semilogy(times,np.maximum(bounds,1e-16),styles[j],color=colors[j],lw=1.5)
  if j==0:
   p=F/'solver-transcript.txt';proven(p,0);ct=p.read_text(errors='replace');co=[float(x) for x in re.findall(r'Global Courant Number \[Explicit VOF Criteria\]\s*:\s*('+NUMBER+')',ct)];tc=np.arange(1,a['steps']+1)*a['dt_s']
  else:
   co=[r['value'] for r in a['courant_records']];tc=np.array([r['time_s'] for r in a['courant_records']]);assert np.all(np.diff(tc)>=-1e-12)
  assert len(co)==len(tc)==a['reported_interface_courant_count'] and max(co)<=.25
  ax[1,1].plot(tc,co,styles[j],color=colors[j],lw=1)
  with (OUT/(['fluent','openfoam'][j]+'-step-residuals.csv')).open('w',newline='') as stream:
   keys=list(a['step_records'][0]['residuals']);w=csv.writer(stream);w.writerow(['step','time_s',*keys,'passed'])
   for row in a['step_records']:w.writerow([row['step'],row['time_s'],*[row['residuals'][k] for k in keys],row['passed']])
  with (OUT/(['fluent','openfoam'][j]+'-courant.csv')).open('w',newline='') as stream:
   w=csv.writer(stream);w.writerow(['time_s','native_reported_interface_courant']);w.writerows(zip(tc,co))
 for axis,limit in [(ax[0,0],1),(ax[0,1],.5),(ax[1,0],1e-6),(ax[1,1],.25)]:axis.axhline(limit,color='black',ls=':',lw=1)
 for axis in ax.ravel():axis.set_xlabel('Physical time (s)');axis.set_xlim(0,.5);axis.grid(alpha=.18)
 ax[0,0].set(title='Worst terminal equation residual / its limit',ylabel='Residual / code-specific limit')
 ax[0,1].set(title='Open-atmosphere mixture-mass closure',ylabel='Absolute closure (% of initial inventory)')
 ax[1,0].set(title='Saved departure of alpha from [0,1]',ylabel='Maximum bound exceedance')
 ax[1,1].set(title='Native VOF/interface Courant reports',ylabel='Reported Courant number',ylim=[0,.29])
 fig.suptitle('Dam break | complete numerical checks',x=.10,ha='left',y=.97,fontsize=20,fontweight='bold')
 handles,names=ax[0,0].get_legend_handles_labels();fig.legend(handles,names,loc='center',bbox_to_anchor=(.5,.88),ncol=2,frameon=False)
 fig.text(.10,.088,'Dotted lines: declared limits. Exact zeros use a 1e-16 plotting floor. Residual normalizations differ between solvers.',fontsize=9)
 fig.text(.10,.045,'Mass closure includes signed, depth-corrected native atmosphere flux at every step. Passing numerics is not physical validation.',fontsize=9)
 fig.savefig(OUT/'numerical-checks.png',dpi=140);fig.savefig(OUT/'numerical-checks.pdf');plt.close(fig)
 keys=sorted(set(k for a in audits for r in a['records'] for k in r))
 with (OUT/'saved-state-checks.csv').open('w',newline='') as stream:
  w=csv.DictWriter(stream,fieldnames=['solver']+keys);w.writeheader()
  for label,a in zip(labels,audits):
   for row in a['records']:w.writerow(dict(solver=label,**row))
 mismatch=np.sum(abs(alpha[0]-alpha[1])*g['areas'],axis=1)/WATER_AREA
 record=dict(accepted_pair_selection_sha256=sha(selection_path),openfoam_continued_from_verified_prefix=True,rendered_png_sha256={str(p):sha(p) for p in sorted(OUT.glob('*.png'))},all_saved_field_checks_pass=True,rendering=dict(frames=101,frame_duration_ms=1000/12),plot=dict(times_s=times.tolist(),global_water_dominated_speed_scale_max_m_s=vmax,global_water_dominated_pressure_bound_kPa=pbound),cross_code_alpha_area_L1_over_initial_water_area=mismatch.tolist(),difference_is_not_accuracy_error=True,numerical_audit_sha256={str(p):sha(p) for p in audit_paths},native_source_sha256={str(p):sha(p) for p in sources},source_sha256={str(Path(__file__)):sha(Path(__file__))},physical_validation=False,mesh_independence=False,time_step_independence=False,interpolated_motion=False)
 (OUT/'native-field-audit.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n');shutil.copyfile(__file__,OUT/'executed-renderer.py');print(OUT)
if __name__=='__main__':main()