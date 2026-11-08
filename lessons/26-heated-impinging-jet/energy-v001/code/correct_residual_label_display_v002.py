from pathlib import Path
import json,shutil,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab');sys.path.insert(0,str(S))
from jet_energy_common_v001 import sha,dump
A=R/'exports/figures/2026-11-08-heated-jet-energy-full-v001';B=R/'exports/figures/2026-11-08-heated-jet-energy-full-v002';assert not B.exists();B.mkdir()
for p in A.iterdir():
 if p.is_file() and p.name not in ['thermal-residuals.png','thermal-residuals.pdf','native-field-audit.json']:shutil.copyfile(p,B/p.name)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':'#203648','text.color':'#203648','xtick.color':'#203648','ytick.color':'#203648'})
fig,axes=plt.subplots(1,2,figsize=(11,4.2),layout='constrained')
for ax,code,label,color in zip(axes,['fluent','openfoam'],['ANSYS Fluent','OpenFOAM'],['#1768a6','#cf692d']):
 p=R/'native/runs/showcase-heated-jet-energy-v001'/(code+'-production-v001')/'production-audit-v001.json';a=json.loads(p.read_text());assert a['numerical_pass'];ax.semilogy([r['time_s'] for r in a['step_residuals']],[r['residual'] for r in a['step_residuals']],color=color);ax.axhline(a['thermal_residual_limit'],color='#777777',ls=':',label='Original criterion');ax.set(title=label,xlabel='Physical time (s)',ylabel='Native terminal thermal residual',xlim=(0,3));ax.grid(alpha=.15);ax.legend(frameon=False,fontsize=8)
fig.suptitle('All 600 steps checked: Fluent scaled energy residual / OpenFOAM final linear T residual\nThe two residual norms differ; their magnitudes do not rank physical accuracy.',size=11);fig.savefig(B/'thermal-residuals.png',dpi=160);fig.savefig(B/'thermal-residuals.pdf');plt.close(fig)
a=json.loads((A/'native-field-audit.json').read_text());a['display_parent']=dict(path=str(A/'native-field-audit.json'),sha256=sha(A/'native-field-audit.json'));a['static_png_sha256']={str(p):sha(p) for p in B.glob('*.png')};a['display_changes']=['Insert a missing space in All 600 title; no numerical data, frame or video change'];dump(B/'native-field-audit.json',a)
unchanged=[p.name for p in B.iterdir() if p.is_file() and p.name not in ['thermal-residuals.png','thermal-residuals.pdf','native-field-audit.json']]
assert all(sha(B/n)==sha(A/n) for n in unchanged)
with Image.open(B/'thermal-residuals.png') as im:im.load()
dump(B/'display-derivation-v002.json',dict(parent=str(A),source_sha256=sha(Path(__file__)),unchanged_files={n:sha(B/n) for n in unchanged},data_and_all61_native_frames_unchanged=True,video_still_from_v001=True,changed_files=['thermal-residuals.png','thermal-residuals.pdf','native-field-audit.json']))
print('Display-only corrected residual label; all data and native video frames unchanged')
