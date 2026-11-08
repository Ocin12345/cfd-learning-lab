from pathlib import Path
from datetime import datetime
import csv,hashlib,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
R=Path('E:/CFD-Learning-Lab');S=Path(__file__).parent;OLD=R/'exports/figures/2026-11-08-airfoil-TE-paired-v002';O=R/'exports/figures/2026-11-08-airfoil-TE-paired-v003'
assert not O.exists()
source=OLD/'wall-pressure-samples.csv';digest=hashlib.sha256(source.read_bytes()).hexdigest()
with source.open(newline='') as f:rows=[{k:v if k in ('solver','side') else float(v) for k,v in r.items()} for r in csv.DictReader(f)]
colors=['#637383','#087e8b','#d17b24'];labels=['Baseline','Local level 1','Local level 2']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
fig,axes=plt.subplots(2,2,figsize=(12,9.8));fig.subplots_adjust(left=.08,right=.98,bottom=.14,top=.82,hspace=.34,wspace=.18)
fig.suptitle('Pressure near the sharp trailing edge\nActual finite-volume wall samples | local refinement sensitivity',y=.985,fontsize=15)
for j,code in enumerate(('ANSYS Fluent','OpenFOAM')):
 for level in (0,1,2):
  for side in ('upper','lower'):
   selected=sorted([r for r in rows if r['solver']==code and int(r['level'])==level and r['side']==side],key=lambda r:r['x_over_c']);assert selected
   x=[r['x_over_c'] for r in selected];cp=[r['wall_Cp'] for r in selected]
   for i,ax in enumerate(axes[:,j]):ax.plot(x,cp,color=colors[level],ls='-' if side=='upper' else '--',lw=1.2,marker='.' if i else None,ms=2)
 for i,ax in enumerate(axes[:,j]):
  ax.set(xlim=(0,1) if i==0 else (.985,1.0002),xlabel='x / chord in airfoil frame',ylabel=r'Wall $C_p$',title=code+(' · full chord' if i==0 else ' · final 1.5% chord'));ax.grid(alpha=.17)
for i in (0,1):
 values=[r['wall_Cp'] for r in rows if i==0 or .985<=r['x_over_c']<=1.0002];lo,hi=min(values),max(values);m=max(.02,.04*(hi-lo))
 for ax in axes[i]:ax.set_ylim(hi+m,lo-m)
handles=[Line2D([],[],color=colors[i],label=labels[i]) for i in (0,1,2)]+[Line2D([],[],color='#354451',ls='-',label='Upper side'),Line2D([],[],color='#354451',ls='--',label='Lower side')]
fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.91),ncols=5,frameon=False,fontsize=9)
fig.text(.5,.035,'Fluent: native face pressure. OpenFOAM: zeroGradient wall value = owner-cell pressure.\nSample coordinates move with refinement; these are not exact continuum point pressures.',ha='center',va='center',fontsize=9)
O.mkdir()
copies={}
for p in OLD.iterdir():
 if p.name in ('airfoil-wall-pressure-sensitivity.png','airfoil-wall-pressure-sensitivity.pdf'):continue
 assert p.is_file();shutil.copy2(p,O/p.name);copies[p.name]=hashlib.sha256(p.read_bytes()).hexdigest();assert hashlib.sha256((O/p.name).read_bytes()).hexdigest()==copies[p.name]
fig.savefig(O/'airfoil-wall-pressure-sensitivity.png',dpi=180);fig.savefig(O/'airfoil-wall-pressure-sensitivity.pdf');plt.close(fig)
assert hashlib.sha256(source.read_bytes()).hexdigest()==digest
record=dict(created_at=datetime.now().astimezone().isoformat(),reason='Original v002 wall-pressure legend overlapped its footer; preserve that version and repair only layout.',data_changed=False,source_csv_sha256=digest,unchanged_copies_sha256=copies,source_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(O/'presentation-correction-v001.json').write_text(json.dumps(record,indent=2)+'\n')
print(str(O/'airfoil-wall-pressure-sensitivity.png'))
