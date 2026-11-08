import wing_import_paths_v001
from wing_repair_common_v001 import *
import numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm,Normalize
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plot_finite_wing_native_v002 import TetGeometry,read_fluent,read_foam,display,decorate_3d
from wing_simplec_qualification_prelaunch_v001 import verify
verify();O=N/'fluent-simplec-qualification-v001';ap=O/'qualification-audit-v001.json';a=json.loads(ap.read_text());assert a['structural_pass'] and not a['numerical_pass']
for p,h in a['source_sha256'].items():assert sha(p)==h,p
out=R/'exports/figures/2026-11-08-wing-SIMPLEC-qualification-v002';assert not out.exists();out.mkdir(parents=True)
geom=TetGeometry(M/'geometry.npz');fields={k:geom.match(v) for k,v in {'Fluent':read_fluent(O/'wing-final.cas.h5',O/'wing-final.dat.h5'),'OpenFOAM':read_foam(OFO,'1500','wing')}.items()}
before={str(p):sha(p) for p in [ap,M/'geometry.npz',W/'native-evidence-inventory-v001.json',N/'fluent-baseline-v001/qualification-audit-v001.json']+[p for v in fields.values() for p in v['inputs']]}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
limit=float(np.ceil(max(np.max(np.abs(2*v['wing_p'])) for v in fields.values())*10)/10);norm=TwoSlopeNorm(vmin=-limit,vcenter=0,vmax=limit);cmap=plt.get_cmap('RdBu_r');fig=plt.figure(figsize=(12,6.2))
for i,(name,v) in enumerate(fields.items()):
 ax=fig.add_axes([.015+i*.49,.17,.48,.70],projection='3d');ax.add_collection3d(Poly3DCollection(display(geom.wing_vertices),facecolors=cmap(norm(2*v['wing_p'])),edgecolor='none'));decorate_3d(ax,'ANSYS Fluent: SIMPLEC 1600 — continuity MISS' if name=='Fluent' else 'OpenFOAM: accepted original 1500 — reused');ax.set(xlim=(-.1,1.25),ylim=(-2.3,2.3),zlim=(-.35,.35));ax.set_box_aspect((1.35,4.6,.7))
cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.23,.14,.54,.025]),orientation='horizontal');cb.set_label(r'$C_p=(p-0)/(\rho U_\infty^2/2)$; identical colour scale')
fig.suptitle('Same refined wing: actual saved pressure',fontsize=16,y=.98);fig.text(.5,.025,'Rounded-cap NACA 0012-like wing · 8° · Re = 400 · no accuracy ranking or physical validation',ha='center',fontsize=10);fig.savefig(out/'paired-wing-pressure.png',dpi=170);fig.savefig(out/'paired-wing-pressure.pdf');plt.close(fig)
# Actual native-cell cross-flow magnitude on a shared tetrahedral cut, no velocity interpolation.
xyz=geom.points[geom.tets];chosen=np.flatnonzero((xyz[:,:,0].min(1)<=1.5)&(xyz[:,:,0].max(1)>=1.5));tri=[];owners=[]
for cell in chosen:
 vertices=[]
 for j in np.flatnonzero(np.abs(xyz[cell,:,0]-1.5)<=1e-12):vertices.append(xyz[cell,j])
 for j,k in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]:
  xj,xk=xyz[cell,j,0],xyz[cell,k,0]
  if (xj-1.5)*(xk-1.5)<0:vertices.append(xyz[cell,j]+(1.5-xj)/(xk-xj)*(xyz[cell,k]-xyz[cell,j]))
 if vertices:
  _,ids=np.unique(np.rint(np.asarray(vertices)/1e-10).astype(np.int64),axis=0,return_index=True);v=np.asarray(vertices)[ids]
  if len(v) not in (3,4):continue
  c=v.mean(0);v=v[np.argsort(np.arctan2(v[:,1]-c[1],v[:,2]-c[2]))]
  for j in range(1,len(v)-1):tri.append(v[[0,j,j+1]][:,[2,1]]);owners.append(cell)
tri=np.asarray(tri);owners=np.asarray(owners);assert len(tri)>10
values={k:np.hypot(v['U'][owners,1],v['U'][owners,2]) for k,v in fields.items()};vmax=max(float(x.max()) for x in values.values());norm=Normalize(0,vmax);fig,axs=plt.subplots(2,1,figsize=(11,6.5));fig.subplots_adjust(top=.87,bottom=.25,hspace=.4)
for ax,(name,v) in zip(axs,values.items()):
 ax.add_collection(PolyCollection(tri,array=v,cmap='viridis',norm=norm,edgecolor='none',rasterized=True));ax.set(xlim=(-2.7,2.7),ylim=(-1,1),aspect='equal',xlabel='Spanwise z (m)',ylabel='Vertical y (m)',title=name+(' SIMPLEC 1600: continuity MISS' if name=='Fluent' else ' original accepted 1500: reused'))
cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap='viridis'),cax=fig.add_axes([.26,.13,.48,.025]),orientation='horizontal');cb.set_label(r'Cross-flow speed $\sqrt{U_y^2+U_z^2}$ (m/s); shared 0–'+f'{vmax:.4g}')
fig.suptitle('Actual downstream cross-flow: x = 1.5 m',fontsize=16,y=.98);fig.text(.5,.025,'Native cell values on common tetrahedral intersections · steady saved fields, not a transient simulation',ha='center',fontsize=9);fig.savefig(out/'paired-downstream-crossflow.png',dpi=170);fig.savefig(out/'paired-downstream-crossflow.pdf');plt.close(fig)
b=json.loads((N/'fluent-baseline-v001/qualification-audit-v001.json').read_text());fig,axs=plt.subplots(2,2,figsize=(12,8.4));fig.subplots_adjust(top=.88,bottom=.11,hspace=.44,wspace=.3)
for data,label,color in [(b,'Coupled: original controls','#3b6a93'),(a,'SIMPLEC: declared bundle','#ad512e')]:
 h=data['new_residual_history'];f=data['new_force_history'];it=[r['iteration'] for r in h];axs[0,0].semilogy(it,[r['continuity'] for r in h],label=label,color=color,lw=1.6)
 for ax,k in [(axs[0,1],'Cd'),(axs[1,0],'Cl')]:ax.plot([r['iteration'] for r in f],[(r[k]/json.loads((W/'native-evidence-inventory-v001.json').read_text())['Fluent']['forces'][k]-1)*1e6 for r in f],color=color,lw=1.6,label=label);ax.ticklabel_format(axis='y',style='plain',useOffset=False)
axs[0,0].axhline(1e-5,color='#222',ls='--',lw=1,label='Unchanged target: 1e−5');axs[0,0].set(ylabel='Scaled continuity',title='Same retained continuity denominator');axs[0,0].legend(fontsize=8)
axs[0,1].set(ylabel='Drag change from original 1500 (ppm)',title='100 actual new drag samples');axs[1,0].set(ylabel='Lift change from original 1500 (ppm)',title='100 actual new lift samples')
s=a['states'];axs[1,1].semilogy([r['iteration'] for r in s],[r['native']['cell_mass_imbalance_abs_sum_kg_s'] for r in s],'o-',color='#ad512e',label='Five saved native states');axs[1,1].axhline(6.536745231970203e-7,color='#222',ls='--',label='Original endpoint ceiling');axs[1,1].set(ylabel='Σ |cell mass imbalance| (kg/s)',title='Corrected native face-flux evidence');axs[1,1].legend(fontsize=8)
for ax in axs.flat:ax.set_xlabel('Actual steady solver iteration');ax.grid(alpha=.2)
fig.suptitle('SIMPLEC qualification: original convergence target still missed',fontsize=16,y=.97);fig.text(.5,.025,'100 iterations from exact original 1500 · stable forces and small mass error do not establish residual convergence',ha='center',fontsize=10);fig.savefig(out/'qualification-numerics.png',dpi=170);fig.savefig(out/'qualification-numerics.pdf');plt.close(fig)
np.savez_compressed(out/'native-fields.npz',wing_vertices_m=geom.wing_vertices,fluent_pressure_Pa=fields['Fluent']['wing_p'],openfoam_pressure_Pa=fields['OpenFOAM']['wing_p'],cut_triangles_zy_m=tri,cut_cell_ids=owners,fluent_crossflow_speed_m_s=values['Fluent'],openfoam_crossflow_speed_m_s=values['OpenFOAM'])
for p,h in before.items():assert sha(p)==h,p
verify();dump(W/'simplec-qualification-visual-binding-v002.json',dict(actual_selected_case=str(O),native_inputs_sha256=before,figures_sha256={str(out/n):sha(out/n) for n in ['paired-wing-pressure.png','paired-downstream-crossflow.png','qualification-numerics.png','native-fields.npz']},shared_Cp_limits=[-limit,limit],shared_crossflow_limits_m_s=[0,vmax],pressure_methods={k:v['pressure_method'] for k,v in fields.items()},crossflow_method='piecewise constant native cell velocities on actual common tetrahedral cut',fluent_numerical_pass=False,openfoam_reused_not_rerun=True,no_transient_or_video_claim=True,actual_inspection_pending=True));print(out)
