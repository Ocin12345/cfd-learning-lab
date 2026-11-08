from pathlib import Path
import argparse,csv,hashlib,json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent
C=['#637383','#087e8b','#d17b24'];L=['Baseline','Local 1','Local 2']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load_csv(p):
 with p.open(newline='') as f:return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def main(selection,out):
 assert not out.exists();chosen=json.loads(selection.read_text());base=R/'native/runs/2026-11-06-airfoil-incidence-paired-v001';inputs={str(selection):sha(selection)};data={};summary=[]
 import analyse_openfoam_airfoil_v002 as av
 old=R/'exports/figures/2026-11-06-airfoil-incidence-paired-v002/native-input-hashes-private.json';oldhash=json.loads(old.read_text())
 for n,h in oldhash.items():assert sha(Path(n))==h
 inputs[str(old)]=sha(old);inputs.update(oldhash)
 for level in (0,1,2):
  f=base/'fluent-alpha5' if level==0 else Path(chosen['fluent'][str(level)]['case'])
  forces=load_csv(f/'force-history.csv');text=(f/'solver-transcript.txt').read_text();raw=re.findall(r'^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)',text,re.M);res={int(r[0]):max(map(float,r[1:])) for r in raw};data['fluent',level]=(res,forces)
  for p in [f/'force-history.csv',f/'solver-transcript.txt']:inputs[str(p)]=sha(p)
  if level==0:
   o=base/'openfoam-alpha5-v002';resrows=av.residual_history([o/'solver.log',o/'solver-continuation.log']);of,_=av.force_history(o/'postProcessing/forceCoefficients')
  else:
   check=Path(chosen['openfoam'][str(level)]['audit']);a=json.loads(check.read_text());assert a['numerical_checks_passed']
   for n,h in a['source_sha256'].items():assert sha(Path(n))==h
   inputs[str(check)]=sha(check);inputs.update(a['source_sha256']);resrows=load_csv(check.parent/'residuals.csv');of=load_csv(check.parent/'forces.csv')
   for p in [check.parent/'residuals.csv',check.parent/'forces.csv']:inputs[str(p)]=sha(p)
  data['openfoam',level]=({int(r['iteration']):max(r[k] for k in ('Ux','Uy','p')) for r in resrows},of)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
 fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained');fig.suptitle('Airfoil numerical evidence | accepted steady endpoints\n5° · Re = 100 · local trailing-edge sensitivity',fontsize=16)
 for j,code in enumerate(('fluent','openfoam')):
  for level in (0,1,2):
   res,forces=data[code,level];x=np.array(sorted(res));y=np.array([res[k] for k in x]);ax[0,j].semilogy(x,y,color=C[level],lw=1.2,label=L[level])
   tail=forces[-100:];cd=np.array([r['drag_coefficient' if code=='fluent' else 'Cd'] for r in tail]);xx=np.arange(1-len(tail),1);change=100*(cd-cd[-1])/abs(cd[-1]);ax[1,j].plot(xx,change,color=C[level],lw=1.3,label=L[level])
  ax[0,j].axhline(1e-7,color='#263642',ls='--',lw=.9,label='Acceptance target');ax[0,j].set(xlabel='Steady iteration',ylabel='Largest terminal equation residual',title=('Fluent · scaled residuals' if code=='fluent' else 'OpenFOAM · maximum initial residual per iteration'));ax[0,j].grid(alpha=.15)
  ax[1,j].set(xlabel='Iteration relative to accepted endpoint',ylabel='Cd change from final value (%)',title='Final 100 force samples');ax[1,j].grid(alpha=.15)
 ax[0,1].legend(frameon=False,fontsize=9);fig.supxlabel('Residual definitions differ between codes. OpenFOAM uses the maximum initial residual over all pressure corrections.\nFlat Fluent force windows reflect printed report precision. Numerical convergence does not establish discretization accuracy.',fontsize=9)
 out.mkdir(parents=True);fig.savefig(out/'airfoil-numerical-checks.png',dpi=180);fig.savefig(out/'airfoil-numerical-checks.pdf');plt.close(fig)
 inputs[str(Path(__file__))]=sha(Path(__file__));(out/'numerical-plot-audit.json').write_text(json.dumps(dict(source_sha256=inputs,source_bound=True,residual_norms_equivalent=False,physical_validation=False),indent=2)+'\n')
 print(str(out/'airfoil-numerical-checks.png'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--selection',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();main(a.selection,a.out)
