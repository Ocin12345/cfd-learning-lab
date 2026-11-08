"""Source-bound actual nine-map diagnostic figures; never launches a solver."""
import wing_import_paths_v001
from wing_repair_common_v001 import *
from wing_residual_window_prelaunch_v001 import wait_admission,verify
from wing_cell_identity_v001 import centroid_bijection
from plot_finite_wing_native_v002 import read_fluent_geometry,TetGeometry,read_fluent
import numpy as np,csv,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
O=N/'fluent-retained-residual-window-v001';P=N/'fluent-retained-residual-probe-v002';OUT=W/'retained-window-analysis-v001';FIG=R/'exports/figures/2026-11-08-wing-residual-persistence-v001';STATUS=W/'retained-window-visual-v001-status.json'

def main():
 assert not OUT.exists() and not FIG.exists();pipeline=json.loads((W/'retained-window-v001-status.json').read_text());assert pipeline['stage']=='retained_window_capture_audited' and pipeline['diagnostic_capture_pass']
 wait_admission(lambda **v:mark(STATUS,'waiting_postprocessing_admission',**v));a=json.loads((O/'retained-window-audit-v001.json').read_text());pa=json.loads((P/'retained-probe-audit-v001.json').read_text());assert a['diagnostic_capture_pass'] and pa['diagnostic_capture_pass']
 sources={str(O/'retained-window-audit-v001.json'):sha(O/'retained-window-audit-v001.json'),str(P/'retained-probe-audit-v001.json'):sha(P/'retained-probe-audit-v001.json'),**a['source_sha256'],**pa['source_sha256']}
 for p,h in sources.items():assert sha(p)==h,p
 g=read_fluent_geometry(O/'source.cas.h5');centres=g['points'][g['tets']].mean(axis=1);n=len(centres);geometry=TetGeometry(M/'geometry.npz');_,wing_ids=geometry.map_native(read_fluent(O/'source.cas.h5',O/'source.dat.h5'));wing=np.zeros(n,bool);wing[g['owner'][wing_ids]]=True;wing_max=float(geometry.wing_vertices[...,0].max());frames=[];selection=[]
 for end in [1502,*range(1503,1511)]:
  case=P if end==1502 else O;capture=P if end==1502 else O/'runtime'/str(end);path=capture/'runtime-residual-arrays.npz'
  with np.load(path) as arrays:
   ids,dist=centroid_bijection(centres,arrays['SV_CENTROID'].reshape(n,3));res=np.empty(n);res[ids]=arrays['SV_P_R']
  mag=np.abs(res);order=np.argsort(mag)[::-1];top=order[:int(np.ceil(.01*n))];total=float(mag.sum());step=pa['native'] if end==1502 else a['steps'][end-1503]['native']
  assert np.isclose(total,step['final_residuals']['continuity']*step['native_scaling_denominators']['continuity'],rtol=1e-6,atol=1e-12)
  row=dict(iteration=end,top1_fraction=float(mag[top].sum()/total),top1_min_x_m=float(centres[top,0].min()),top1_max_x_m=float(centres[top,0].max()),wing_adjacent_fraction=float(mag[wing].sum()/total),downstream_fraction=float(mag[centres[:,0]>wing_max].sum()/total),scaled_continuity=step['final_residuals']['continuity'],corrected_cell_mass_sum_kg_s=step['cell_mass_imbalance_abs_sum_kg_s'],actual_captured_cells=n)
  row['descriptive_persistence_pass']=row['top1_fraction']>=.9 and row['top1_min_x_m']>wing_max and row['wing_adjacent_fraction']<=.001
  frames.append(dict(row=row,res=res,top=top))
  native_case=case/'wing-final.cas.h5' if end==1502 else case/f'iter-{end}.cas.h5';native_data=case/'wing-final.dat.h5' if end==1502 else case/f'iter-{end}.dat.h5'
  selection.append(dict(iteration=end,case=str(native_case),data=str(native_data),runtime_arrays=str(path),runtime_sha256=sha(path),audit=str(P/'retained-probe-audit-v001.json' if end==1502 else O/'retained-window-audit-v001.json')))
 OUT.mkdir();FIG.mkdir(parents=True);dump(W/'accepted-residual-diagnostic-selection-v001.json',dict(diagnostic_capture_only=True,original_qualification_pass=False,frames=selection,source_sha256=sources,geometry_sha256=sha(M/'geometry.npz'),nine_actual_maps=True,first1501_map_unavailable=True))
 rows=[v['row'] for v in frames]
 with (OUT/'nine-actual-map-metrics.csv').open('x',newline='',encoding='utf-8') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 vmin=min(float(np.abs(v['res'][v['top']]).min()) for v in frames);vmax=max(float(np.abs(v['res']).max()) for v in frames);norm=LogNorm(vmin=vmin,vmax=vmax);points=np.unique(geometry.wing_vertices.reshape(-1,3)[:,[0,2]],axis=0)
 fig,axs=plt.subplots(1,3,figsize=(13.8,5.3),layout='constrained')
 for ax,j in zip(axs,[0,4,8]):
  v=frames[j];sel=v['top'];ax.scatter(points[:,0],points[:,1],s=.35,color='#777',alpha=.3);sc=ax.scatter(centres[sel,0],centres[sel,2],c=np.abs(v['res'][sel]),norm=norm,cmap='magma',s=6);ax.set(xlim=(-.3,9.2),ylim=(-6.1,6.1),xlabel='x (m)',ylabel='z (m)',title=f"Iteration {v['row']['iteration']}\nTop1% share: {100*v['row']['top1_fraction']:.2f}%");ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.15)
 fig.colorbar(sc,ax=list(axs),shrink=.85,label='Absolute retained continuity residual per cell (kg/s)',pad=.015)
 fig.suptitle('Actual residual maps: identical solver controls and colour scale',fontsize=14);fig.savefig(FIG/'three-actual-residual-maps.png',dpi=170);fig.savefig(FIG/'three-actual-residual-maps.pdf');plt.close(fig)
 fig,axs=plt.subplots(2,2,figsize=(11.7,7.8),layout='constrained');x=[r['iteration'] for r in rows]
 axs[0,0].plot(x,[100*r['top1_fraction'] for r in rows],color='#237b72',marker='o',ms=4);axs[0,0].axhline(90,color='#777',ls='--',lw=1);axs[0,0].set(title='(a) Concentration in top1% of cells',ylabel='Share of absolute residual (%)',ylim=(89,100))
 axs[0,1].plot(x,[r['top1_min_x_m'] for r in rows],color='#d4833e',marker='o',ms=4);axs[0,1].axhline(wing_max,color='#777',ls='--',lw=1,label='Wing maximum x');axs[0,1].set(title='(b) Nearest top1% cell in streamwise direction',ylabel='Minimum x (m)');axs[0,1].legend(frameon=False,fontsize=8)
 axs[1,0].semilogy(x,[r['scaled_continuity'] for r in rows],color='#58758e',marker='o',ms=4);axs[1,0].axhline(1e-5,color='#777',ls='--',lw=1,label='Original convergence target');axs[1,0].set(title='(c) Original convergence target remains separate',ylabel='Native scaled continuity residual');axs[1,0].legend(frameon=False,fontsize=8)
 axs[1,1].semilogy(x,[100*r['wing_adjacent_fraction'] for r in rows],color='#8e609c',marker='o',ms=4);axs[1,1].axhline(.1,color='#777',ls='--',lw=1);axs[1,1].set(title='(d) Cells directly touching the wing',ylabel='Share of absolute residual (%)')
 for ax in axs.flat:ax.set_xlabel('Steady solver iteration');ax.grid(alpha=.15)
 fig.suptitle('Finite diagnostic persistence: nine captured iterations, not physical time',fontsize=14);fig.savefig(FIG/'persistence-and-convergence.png',dpi=170);fig.savefig(FIG/'persistence-and-convergence.pdf');plt.close(fig)
 summary=dict(diagnostic_capture_pass=True,nine_observed_iterations=[r['iteration'] for r in rows],eight_new_steps=8,total_new_steps_from1500=10,descriptive_pattern_persists_all_nine=all(r['descriptive_persistence_pass'] for r in rows),top1_fraction_range=[min(r['top1_fraction'] for r in rows),max(r['top1_fraction'] for r in rows)],wing_adjacent_fraction_range=[min(r['wing_adjacent_fraction'] for r in rows),max(r['wing_adjacent_fraction'] for r in rows)],top1_minimum_x_range_m=[min(r['top1_min_x_m'] for r in rows),max(r['top1_min_x_m'] for r in rows)],scaled_continuity_range=[min(r['scaled_continuity'] for r in rows),max(r['scaled_continuity'] for r in rows)],original_convergence_all_nine=all(r['scaled_continuity']<=1e-5 for r in rows),one_snapshot_does_not_prove_root_cause=True,nine_observations_not_full_cycle_or_validation=True,all_numerical_limits_unchanged=True,figures=[str(FIG/f) for f in ['three-actual-residual-maps.png','persistence-and-convergence.png']],figures_sha256={str(FIG/f):sha(FIG/f) for f in ['three-actual-residual-maps.png','persistence-and-convergence.png']})
 for p,h in sources.items():assert sha(p)==h,p
 verify();dump(OUT/'summary.json',summary);mark(STATUS,'visuals_ready_for_actual_inspection',summary=str(OUT/'summary.json'),original_qualification_pass=False);print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
