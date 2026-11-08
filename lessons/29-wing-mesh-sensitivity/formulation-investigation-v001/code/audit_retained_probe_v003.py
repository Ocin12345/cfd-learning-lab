import wing_import_paths_v001
from wing_repair_common_v001 import *
from wing_repair_native_v003 import raw_residual,same_native,state_keys
from wing_retained_residual_helpers_v001 import expected_probe_settings,classify_stage
from wing_retained_residual_prelaunch_v002 import verify
from wing_repair_helpers_v001 import same_residual_prefix,force_window
from wing_fluent_v001 import force_rows
from audit_wing_continuation_v001 import native
from plot_finite_wing_native_v002 import read_fluent,indexed,TetGeometry
from sloshing_fixed_settings_v001 import compare_fixed
from wing_amg_transcript_v001 import parse_blocks
import h5py,numpy as np,csv,matplotlib
from wing_cell_identity_v001 import centroid_bijection
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
O=N/'fluent-retained-residual-probe-v002';P=N/'fluent-retained-residual-probe-v001';OUT=W/'retained-probe-analysis-v002';FIG=R/'exports/figures/2026-11-08-wing-retained-residual-probe-v002';assert not OUT.exists() and not FIG.exists();OUT.mkdir();FIG.mkdir(parents=True)
assert verify()==278;assert json.loads((O/'status.json').read_text())['stage']=='closed'
for x,y in [('wing-final.cas.h5','source.cas.h5'),('wing-final.dat.h5','source.dat.h5'),('final-settings.json','expected-source-settings.json')]:assert sha(P/x)==sha(O/y)
compare_fixed(json.loads((O/'loaded-source-settings.json').read_text()),{k:v for k,v in json.loads((P/'final-settings.json').read_text()).items() if k in state_keys()});compare_fixed(json.loads((O/'candidate-settings.json').read_text()),expected_probe_settings());compare_fixed(json.loads((O/'final-settings.json').read_text()),expected_probe_settings())
same_native(O/'loaded-restart.cas.h5',O/'loaded-restart.dat.h5',O/'source.cas.h5',O/'source.dat.h5')
pi,pr=raw_residual(O/'source.dat.h5');it,raw=raw_residual(O/'wing-final.dat.h5');same_residual_prefix(pi,pr,it,raw,1501);assert int(it[-1])==1502 and int(pi[-1])==1501
n_result,_,_,g=native(O/'wing-final.cas.h5',O/'wing-final.dat.h5');assert n_result['iteration']==1502
baseline=json.loads((N/'fluent-baseline-v001/qualification-audit-v001.json').read_text());expected=baseline['new_residual_history'][1]
for k in EQ:assert np.isclose(raw[k][-1,0]/raw[k][-1,1],expected[k],rtol=1e-8,atol=1e-12),k
parent=force_rows((P/'forces-history.out').read_text())[-1];forces=force_window((O/'forces-history.out').read_text(),1502,1502,parent)
for k in ['Cd','Cl']:assert np.isclose(forces[0][k],baseline['new_force_history'][1][k],rtol=0,atol=1e-10),k
native_field=read_fluent(O/'wing-final.cas.h5',O/'wing-final.dat.h5');n=len(g['tets']);centres=g['points'][g['tets']].mean(axis=1)
with np.load(O/'runtime-residual-arrays.npz') as a:
 api_centres=a['SV_CENTROID'].reshape(n,3);index,distance=centroid_bijection(centres,api_centres);assert np.max(distance)<2e-10 and len(np.unique(index))==n
 api_u=np.column_stack([a[k] for k in ['SV_U','SV_V','SV_W']]);assert np.allclose(api_u,native_field['U'][index],rtol=1e-10,atol=1e-12);assert np.allclose(a['SV_P'],native_field['p'][index],rtol=1e-10,atol=1e-12)
 values={};sums={}
 for var,k in zip(['SV_P_R','SV_U_R','SV_V_R','SV_W_R'],EQ):
  v=np.asarray(a[var]);assert v.shape==(n,) and np.isfinite(v).all();values[var]=np.empty(n);values[var][index]=v;sums[var]=float(np.abs(v).sum());assert np.isclose(sums[var],raw[k][-1,0],rtol=1e-6,atol=1e-12),(var,sums[var],raw[k][-1,0])
res=values['SV_P_R'];stage=classify_stage(res,float(raw['continuity'][-1,0]),n_result['cell_mass_imbalance_abs_sum_kg_s']);assert stage=='assembled_initial'
block=parse_blocks((O/'solver-transcript.txt').read_text(),1502,1502)[0];assert np.isclose(block['initial'][0],sums['SV_P_R'],rtol=6e-5,atol=1e-12);assert np.isclose(block['terminal'][0],n_result['cell_mass_imbalance_abs_sum_kg_s'],rtol=6e-5,atol=1e-12)
geometry=TetGeometry(M/'geometry.npz');mapping,wing_ids=geometry.map_native(native_field);wall=np.zeros(n,dtype=bool);wall[g['owner'][wing_ids]]=True;outer=np.zeros(n,dtype=bool);bf=np.flatnonzero(g['neighbour']<0);outer[g['owner'][np.setdiff1d(bf,wing_ids)]]=True;assert not np.any(wall&outer);inside=~(wall|outer)
mag=np.abs(res);order=np.argsort(mag)[::-1];total=float(mag.sum());stats={}
for f in [.01,.05,.10]:
 selected=order[:int(np.ceil(n*f))];stats[str(f)]=dict(cells=len(selected),fraction_absolute_residual=float(mag[selected].sum()/total),centre_bounds_m=[centres[selected].min(axis=0).tolist(),centres[selected].max(axis=0).tolist()],wing_adjacent_cells=int(wall[selected].sum()),outer_boundary_adjacent_cells=int(outer[selected].sum()))
regions={k:dict(cells=int(mask.sum()),fraction_absolute_residual=float(mag[mask].sum()/total)) for k,mask in [('Wing-adjacent cells',wall),('Outer-boundary-adjacent cells',outer),('Other cells',inside)]}
with (OUT/'retained-continuity-cells.csv').open('x',newline='',encoding='utf-8') as f:
 w=csv.writer(f);w.writerow(['native_cell_0based','x_m','y_m','z_m','retained_SV_P_R_kg_s','wing_adjacent','outer_boundary_adjacent']);w.writerows((j,*centres[j],res[j],int(wall[j]),int(outer[j])) for j in range(n))
np.savez_compressed(OUT/'mapped-retained-residuals.npz',native_centres_m=centres,native_cell_0based=np.arange(n),**values)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axs=plt.subplots(2,2,figsize=(12,8.2),layout='constrained');sel=order[:int(np.ceil(.01*n))];norm=LogNorm(vmin=float(mag[sel].min()),vmax=float(mag[sel].max()));wingpoints=geometry.wing_vertices.reshape(-1,3)
for ax,dim,ylab,title in [(axs[0,0],2,'z (m)','(a) Top 1% of cell residuals: plan view'),(axs[0,1],1,'y (m)','(b) Same cells: side view')]:
 p=np.unique(wingpoints[:,[0,dim]],axis=0);ax.scatter(p[:,0],p[:,1],color='#777',s=.4,alpha=.25,label='Wing surface projection',zorder=1);sc=ax.scatter(centres[sel,0],centres[sel,dim],c=mag[sel],s=8,cmap='magma',norm=norm);ax.set(xlabel='x (m)',ylabel=ylab,title=title);ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.15)
fig.colorbar(sc,ax=list(axs[0]),label='Absolute retained cell continuity residual (kg/s)',shrink=.85,pad=.02)
cumulative=np.cumsum(mag[order])/total;axs[1,0].semilogx(np.arange(1,n+1)/n*100,cumulative*100,color='#237b72',lw=1.6)
for f in [.01,.05,.10]:
 k=int(np.ceil(n*f))-1;axs[1,0].scatter([(k+1)/n*100],[cumulative[k]*100],s=23,color='#237b72')
axs[1,0].set(xlabel='Cells ranked by absolute residual (%)',ylabel='Cumulative absolute residual (%)',title='(c) Residual concentration',ylim=(0,101));axs[1,0].grid(alpha=.15)
vals=[100*v['fraction_absolute_residual'] for v in regions.values()];axs[1,1].barh(list(regions),vals,color=['#d4833e','#58758e','#237b72']);axs[1,1].set(xlabel='Share of absolute continuity residual (%)',title='(d) Actual boundary adjacency',xlim=(0,100))
for j,v in enumerate(vals):axs[1,1].text(min(v+1,92),j,f'{v:.2f}%',va='center')
fig.suptitle('Wing diagnostic at iteration 1502: retained assembled continuity residual',fontsize=14);fig.savefig(FIG/'retained-residual-map.png',dpi=170);fig.savefig(FIG/'retained-residual-map.pdf');plt.close(fig)
paths=[O/x for x in ['source.cas.h5','source.dat.h5','loaded-restart.cas.h5','loaded-restart.dat.h5','wing-final.cas.h5','wing-final.dat.h5','solver-transcript.txt','forces-history.out','loaded-source-settings.json','candidate-settings.json','final-settings.json','runtime-residual-arrays.npz','runtime-residual-catalogue.json','initial-field-check.json','geometry-quality.json','solver-version.json']];paths+=list((O/'force-reports').glob('*.txt'))
report=dict(diagnostic_capture_pass=True,iteration=1502,actual_new_steps_total=2,trajectory='one step original1500->1501; one exact checkpoint continuation1501->1502',first_capture_missing_runtime_export_preserved=True,source_history_and_settings_exact=True,baseline_iteration1502_residual_force_replay_pass=True,runtime_to_native_centroid_bijection=True,maximum_centroid_match_distance_m=float(distance.max()),runtime_U_p_match_native=True,all_four_retained_residual_sums_match_raw_HDF=True,stage=stage,residual_absolute_sums=sums,cell_statistics=stats,regions=regions,native=n_result,forces=forces[0],original_residual_acceptance_pass=max(n_result['final_residuals'].values())<=1e-5,numerical_convergence_not_claimed=True,root_cause_not_yet_proved=True,nine_iteration_window_not_admitted=True,figure=str(FIG/'retained-residual-map.png'),figure_sha256=sha(FIG/'retained-residual-map.png'),source_sha256={str(p):sha(p) for p in paths})
assert verify()==278;dump(O/'retained-probe-audit-v001.json',report);dump(OUT/'diagnostic-summary.json',report)
print(json.dumps(dict(diagnostic_capture_pass=True,stage=stage,top_cells=stats,regions=regions,original_numerical_pass=report['original_residual_acceptance_pass'],figure=report['figure']),indent=2))
