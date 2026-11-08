"""Read-only empirical cross-check of coupled matrix, residual history and saved fluxes."""
import wing_import_paths_v001
from wing_repair_common_v001 import *
from wing_repair_native_v003 import raw_residual
from wing_amg_transcript_v001 import parse_blocks
from plot_finite_wing_native_v002 import read_fluent_geometry,indexed
import h5py,numpy as np,csv,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
OUT=W/'formulation-investigation-v001';assert not OUT.exists();OUT.mkdir()
FIG=R/'exports/figures/2026-11-08-wing-formulation-investigation-v001';assert not FIG.exists();FIG.mkdir(parents=True)
labels=['baseline','tight-amg','low-cfl10','low-relax25'];names=['Original controls','Tight AMG','Courant 10','Relaxation 0.25'];colors=['#58758e','#d4833e','#237b72','#8e609c']
inputs=json.loads((W/'qualification-inputs-v004.json').read_text());inputs[str(W/'qualification-inputs-v004.json')]=sha(W/'qualification-inputs-v004.json')
for label in labels:
 p=N/('fluent-'+label+'-v001')/'qualification-audit-v001.json';a=json.loads(p.read_text());assert a['structural_pass'] and not a['numerical_pass'];inputs[str(p)]=sha(p);inputs.update(a['source_sha256'])
for p,h in inputs.items():assert sha(p)==h,p
geom=read_fluent_geometry(F/'wing-final.cas.h5');n=len(geom['tets']);nf=len(geom['faces']);xyz=geom['points'][geom['tets']];vol=np.abs(np.linalg.det(np.stack([xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0],xyz[:,3]-xyz[:,0]],axis=-1)))/6;assert (vol>0).all();del xyz

def fields(path):
 with h5py.File(path) as h:
  c=h['results/1/phase-1/cells'];u=np.column_stack([indexed(c[k],n) for k in ['SV_U','SV_V','SV_W']]);p=indexed(c['SV_P'],n);flux=indexed(h['results/1/phase-1/faces/SV_FLUX'],nf)
  assert np.isfinite(u).all() and np.isfinite(p).all() and np.isfinite(flux).all()
  residual_fields=[k for k in c if 'RES' in k.upper()]
 return u,p,flux,residual_fields
u0,p0,phi0,res_fields0=fields(F/'wing-final.dat.h5');iv,raw0=raw_residual(F/'wing-final.dat.h5');den0=float(raw0['continuity'][-1,1]);co=geom['owner'];cn=geom['neighbour'];internal=cn>=0
reports={};rows=[];checkpoints=[]
for label in labels:
 case=N/('fluent-'+label+'-v001');a=json.loads((case/'qualification-audit-v001.json').read_text());blocks=parse_blocks((case/'solver-transcript.txt').read_text(errors='strict'),1501,1600);its,raw=raw_residual(case/'wing-final.dat.h5');lookup={int(i):j for j,i in enumerate(its)}
 assert np.array_equal(its[:1500],iv) and all(np.array_equal(raw[k][:1500],raw0[k]) for k in EQ)
 for block in blocks:
  j=lookup[block['iteration']];nums=np.array([raw[k][j,0] for k in EQ]);dens=np.array([raw[k][j,1] for k in EQ]);assert raw['continuity'][j,1]==den0
  assert np.allclose(block['initial'],nums,rtol=6e-5,atol=1e-12),(label,block['iteration'],'AMGinitial vs HDF raw')
  assert np.allclose(block['outer_printed'],nums/dens,rtol=6e-5,atol=1e-14),(label,block['iteration'],'printed vs HDF scaled')
  row=dict(pilot=label,iteration=block['iteration'],outer_scaled=float(nums[0]/dens[0]),scale_denominator=float(dens[0]),raw_initial_continuity=float(nums[0]),AMG_printed_initial_continuity=block['initial'][0],AMG_printed_final_continuity=block['terminal'][0],inner_cycles=block['cycles'],inner_reduction=block['terminal'][0]/block['initial'][0]);rows.append(row)
 for end in [1520,1540,1560,1580,1600]:
  u,p,phi,resfields=fields(case/f'iter-{end}.dat.h5');balance=np.bincount(co,weights=phi,minlength=n)-np.bincount(cn[internal],weights=phi[internal],minlength=n);abs_mass=float(np.abs(balance).sum());block=blocks[end-1501]
  assert np.isclose(abs_mass,block['terminal'][0],rtol=6e-5,atol=1e-12),(label,end,'AMGfinal vs actual saved flux')
  delta=u-u0;pressure=p-p0
  checkpoints.append(dict(pilot=label,iteration=end,AMG_printed_final_continuity=block['terminal'][0],native_cell_mass_abs_sum_kg_s=abs_mass,absolute_match_error_kg_s=abs(abs_mass-block['terminal'][0]),velocity_volume_L2_relative_percent=100*float(np.sqrt(np.sum(vol*np.sum(delta*delta,axis=1))/np.sum(vol*np.sum(u0*u0,axis=1)))),maximum_velocity_vector_change_m_s=float(np.max(np.linalg.norm(delta,axis=1))),pressure_volume_RMS_change_Pa=float(np.sqrt(np.sum(vol*pressure**2)/np.sum(vol))),pressure_RMS_relative_to_qinf_percent=100*float(np.sqrt(np.sum(vol*pressure**2)/np.sum(vol)))/.5,maximum_pressure_change_Pa=float(np.max(np.abs(pressure))),saved_per_cell_residual_fields=resfields))
 ratios=[blocks[j]['initial'][0]/blocks[j-1]['terminal'][0] for j in range(1,len(blocks))]
 reports[label]=dict(AMG_to_native_crosschecks_passed=True,actual_outer_rows=100,saved_native_flux_crosschecks=5,AMG_initial_matches_native_unscaled_numerator=True,AMG_final_matches_actual_saved_flux_sum=True,continuity_scale_unchanged=den0,next_outer_initial_to_previous_inner_final_ratio=dict(min=float(min(ratios)),median=float(np.median(ratios)),max=float(max(ratios))),inner_cycle_range=[min(b['cycles'] for b in blocks),max(b['cycles'] for b in blocks)],final_continuity=a['native']['final_residuals']['continuity'],numerical_pass=False)
for fn,data in [('actual-AMG-stage-history.csv',rows),('actual-native-checkpoint-changes.csv',checkpoints)]:
 with (OUT/fn).open('x',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for label,name,c in zip(labels,names,colors):
 r=[v for v in rows if v['pilot']==label];v=[v for v in checkpoints if v['pilot']==label];axs[0,0].semilogy([x['iteration'] for x in r],[x['raw_initial_continuity'] for x in r],label=name,color=c,lw=1.4);axs[0,1].semilogy([x['iteration'] for x in r],[x['AMG_printed_final_continuity'] for x in r],label=name,color=c,lw=1.3);axs[0,1].scatter([x['iteration'] for x in v],[x['native_cell_mass_abs_sum_kg_s'] for x in v],color=c,s=29,marker='x',zorder=5)
 axs[1,0].plot([x['iteration'] for x in v],[x['velocity_volume_L2_relative_percent'] for x in v],color=c,marker='o',ms=3)
 axs[1,1].plot([x['iteration'] for x in v],[x['pressure_volume_RMS_change_Pa'] for x in v],color=c,marker='o',ms=3)
axs[0,0].axhline(den0*1e-5,color='#555',ls='--',lw=1);axs[0,0].set(title='(a) Start of coupled linear solve',ylabel='Unscaled continuity residual (kg/s)');axs[0,0].legend(frameon=False,fontsize=8,ncols=2,loc='upper center',bbox_to_anchor=(.5,1.00))
axs[0,1].set(title='(b) End of linear solve; crosses = saved flux',ylabel='Absolute cell mass sum (kg/s)')
axs[1,0].set(title='(c) Velocity change from original checkpoint',ylabel='Volume-weighted relative L2 change (%)')
axs[1,1].set(title='(d) Pressure change from original checkpoint',ylabel='Volume-weighted RMS change (Pa)')
for ax in axs.flat:ax.set_xlabel('Steady solver iteration');ax.grid(alpha=.15)
fig.suptitle('Wing formulation investigation: the inner correction works, the outer imbalance returns',fontsize=14)
fig.savefig(FIG/'solver-formulation-evidence.png',dpi=170);fig.savefig(FIG/'solver-formulation-evidence.pdf');plt.close(fig)
for p,h in inputs.items():assert sha(p)==h,p
result=dict(read_only=True,new_solver_iterations=0,new_solver_sessions=0,preserved_frozen_bindings=len(inputs),current_continuity_scale=den0,all_400_initial_and_scaled_crosschecks_passed=True,all_20_postsolve_native_mass_crosschecks_passed=True,reports=reports,checkpoint_changes=checkpoints,source_cell_residuals_retained=False,outer_cell_residual_spatial_localization_unavailable=True,root_cause_not_established=True,original_residual_acceptance_unchanged=True,all_pilots_still_numerical_miss=True,figure=str(FIG/'solver-formulation-evidence.png'),figure_sha256=sha(FIG/'solver-formulation-evidence.png'),no_new_upload=True)
dump(OUT/'investigation.json',result);dump(OUT/'persistent-source-bindings.json',inputs)
print(json.dumps(dict(reports=reports,largest_velocity_L2_change_percent=max(r['velocity_volume_L2_relative_percent'] for r in checkpoints),largest_pressure_RMS_change_Pa=max(r['pressure_volume_RMS_change_Pa'] for r in checkpoints),figure=result['figure']),indent=2))
