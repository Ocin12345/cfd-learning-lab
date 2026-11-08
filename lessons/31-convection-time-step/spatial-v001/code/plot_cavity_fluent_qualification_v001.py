from pathlib import Path
import json,hashlib,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from cavity_native_fields_v001 import fluent
from cavity_restriction_v001 import compare_fields
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent;BASE=R/'native/runs/natural-convection-dt001-v001/fluent-v002'
def main(run,out):
 audit=run/'qualification-audit-v002.json';a=json.loads(audit.read_text());assert a['numerical_pass'] and a['time_steps']==100
 for n,h in a['source_sha256'].items():assert hashlib.sha256(Path(n).read_bytes()).hexdigest()==h,n
 coarse=fluent(BASE/'qualification-final.cas.h5',BASE/'heated-0100.dat.h5',48);fine=fluent(run/'qualification-final.cas.h5',run/'heated-0100.dat.h5',96);comparison=compare_fields(coarse,fine,48,relative=True)
 out.mkdir(parents=True,exist_ok=False);fig,axes=plt.subplots(2,2,figsize=(11,10),layout='constrained');plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
 limit=max(float(np.linalg.norm(q['U'],axis=1).max()) for q in (coarse,fine))*1000
 for j,(q,n,name) in enumerate([(coarse,48,'Coarse | 2,304 cells'),(fine,96,'Refined | 9,216 cells')]):
  edges=np.linspace(0,50,n+1)
  for i,(field,cmap,norm,label) in enumerate([(q['T'].reshape(n,n),'coolwarm',Normalize(299.6,300.4),'Temperature (K)'),(np.linalg.norm(q['U'],axis=1).reshape(n,n)*1000,'viridis',Normalize(0,limit),'Speed (mm/s)')]):
   ax=axes[i,j];im=ax.pcolormesh(edges,edges,field,shading='flat',cmap=cmap,norm=norm,rasterized=True);ax.set(aspect='equal',xlim=(0,50),ylim=(0,50),xlabel='x (mm)',ylabel='y (mm)',title=name if i==0 else 'Actual saved velocity magnitude')
   if j==1:fig.colorbar(im,ax=axes[i,:].tolist(),label=label,shrink=.8)
 fig.suptitle('One second after heating begins | ANSYS Fluent',fontsize=17,weight='bold')
 fig.supxlabel('Native cells, fixed dt = 0.01 s. Finite 100-step qualification; the 20 s comparison has not yet been admitted.',fontsize=10)
 fig.savefig(out/'fluent-cavity-mesh-qualification.png',dpi=160);fig.savefig(out/'fluent-cavity-mesh-qualification.pdf');plt.close(fig)
 record=dict(native_audit=str(audit),audit_sha256=hashlib.sha256(audit.read_bytes()).hexdigest(),fine_case=str(run),coarse_case=str(BASE),same_physical_time_s=1.,same_dt_s=.01,comparison=comparison,fine_initial_comparison=a['initial_fields'],fine_steady_control=json.loads((run/'control-checks.json').read_text())['control'],source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [BASE/'qualification-final.cas.h5',BASE/'heated-0100.dat.h5',run/'qualification-final.cas.h5',run/'heated-0100.dat.h5',Path(__file__)]},actual_inspection_pending=True,full20s_complete=False,physical_validation=False)
 (out/'source-binding.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(comparison,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.run,a.out)