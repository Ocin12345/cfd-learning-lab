from pathlib import Path
import json,hashlib,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from cavity_native_fields_v001 import fluent,foam
from cavity_restriction_v001 import compare_fields
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent;N=R/'native/runs/showcase-convection-spatial-v001';F=N/'fluent-qualification-v004';O=N/'openfoam-qualification-v001';B=R/'native/runs/natural-convection-dt001-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 fa=json.loads((F/'qualification-audit-v002.json').read_text());oa=json.loads((O/'qualification-audit-v001.json').read_text());assert fa['numerical_pass'] and oa['numerical_pass']
 for a in (fa,oa):
  for p,h in a['source_sha256'].items():assert sha(Path(p))==h,p
 out=R/'exports/figures/2026-11-08-convection-spatial-qualification-paired-v001';out.mkdir(parents=True,exist_ok=False)
 ff=fluent(F/'qualification-final.cas.h5',F/'heated-0100.dat.h5',96);of=foam(O,1.,96)
 comparisons=dict(Fluent=compare_fields(fluent(B/'fluent-v002/qualification-final.cas.h5',B/'fluent-v002/heated-0100.dat.h5',48),ff,48),OpenFOAM=compare_fields(foam(B/'openfoam',1.,48),of,48))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 fig,axes=plt.subplots(2,2,figsize=(11,10),layout='constrained');edges=np.linspace(0,50,97);maxspeed=max(np.linalg.norm(q['U'],axis=1).max() for q in (ff,of))*1000
 for j,(name,q) in enumerate([('ANSYS Fluent',ff),('OpenFOAM',of)]):
  for i,(a,cmap,norm,label) in enumerate([(q['T'].reshape(96,96),'coolwarm',Normalize(299.6,300.4),'Temperature (K)'),(np.linalg.norm(q['U'],axis=1).reshape(96,96)*1000,'viridis',Normalize(0,maxspeed),'Speed (mm/s)')]):
   ax=axes[i,j];im=ax.pcolormesh(edges,edges,a,shading='flat',cmap=cmap,norm=norm,rasterized=True);ax.set(aspect='equal',xlabel='x (mm)',ylabel='y (mm)',title=name if i==0 else 'Actual native velocity magnitude')
   if j==1:fig.colorbar(im,ax=axes[i,:].tolist(),label=label,shrink=.8)
 fig.suptitle('Refined cavity: both solvers at 1.00 s',fontsize=17,weight='bold');fig.supxlabel('9,216 cells, dt = 0.01 s. Both finite 100-step qualifications passed; full 20 s comparison remains pending.',fontsize=10)
 for ext in ('png','pdf'):fig.savefig(out/('paired-qualification-fields.'+ext),dpi=160)
 plt.close(fig)
 fs=json.loads((F/'time-step-checks.json').read_text());os=json.loads((O/'all-step-checks.json').read_text());fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for i,label in enumerate(['continuity','Ux','Uy','energy']):axes[0,0].semilogy([r['time_s'] for r in fs],[max(r['residuals'][i]/([1e-6]*3+[1e-9])[i],1e-12) for r in fs],label=label)
 for key,lim in os['targets'].items():axes[0,1].semilogy([r['time_s'] for r in os['steps']],[max(r['residuals'][key]/lim,1e-12) for r in os['steps']],label=key)
 for ax,title in zip(axes[0],['Fluent: terminal native scaled residuals','OpenFOAM: last-outer initial residuals']):ax.axhline(1,c='black',ls='--',lw=1);ax.set(title=title,xlabel='Physical time (s)',ylabel='Native residual / its declared limit');ax.legend(fontsize=8);ax.grid(alpha=.2)
 for name,a,color in [('Fluent',fa,'#17677f'),('OpenFOAM',oa,'#b45a2b')]:
  rows=a['saved_native_states'] if name=='Fluent' else a['native_state_checks'];axes[1,0].plot([r['clock']['time_s'] for r in rows],[r['metrics']['relative_mass_drift']*1e6 for r in rows],'.-',label=name,c=color)
 axes[1,0].set(xlabel='Physical time (s)',ylabel='Absolute relative mass drift (ppm)',title='Closed mass: actual saved native states');axes[1,0].legend();axes[1,0].grid(alpha=.2)
 axes[1,1].bar(['Fluent','OpenFOAM'],[fa['Nu_hot_at_1s'],oa['Nu_hot_at_1s']],color=['#17677f','#b45a2b']);axes[1,1].set(ylabel='Hot-wall Nusselt number',title='Native face heat integration at 1 s');axes[1,1].grid(axis='y',alpha=.2)
 fig.suptitle('Qualification evidence | unchanged numerical limits',fontsize=17,weight='bold');fig.supxlabel('Native residual norms differ. Instantaneous heat transfer is not a full transient energy audit or physical validation.',fontsize=10)
 for ext in ('png','pdf'):fig.savefig(out/('paired-qualification-checks.'+ext),dpi=160)
 plt.close(fig)
 payload=dict(Fluent_audit=str(F/'qualification-audit-v002.json'),OpenFOAM_audit=str(O/'qualification-audit-v001.json'),original_metadata_guard_failure_retained=True,metadata_correction=str(S/'decomposition-header-correction-v002.json'),both_numerical_qualifications_passed=True,time_s=1.,dt_s=.01,cells=9216,coarse_fine_comparisons=comparisons,Nu_hot_at_1s=dict(Fluent=fa['Nu_hot_at_1s'],OpenFOAM=oa['Nu_hot_at_1s']),source_sha256={str(p):sha(p) for p in [F/'qualification-audit-v002.json',O/'qualification-audit-v001.json',S/'decomposition-header-correction-v002.json',Path(__file__)]},full20s_complete=False,physical_validation=False)
 (out/'selection-and-measurements.json').write_text(json.dumps(payload,indent=2)+'\n');(S/'accepted-qualification-pair-v001.json').write_text(json.dumps(payload,indent=2)+'\n');print(json.dumps(comparisons,indent=2))
if __name__=='__main__':main()