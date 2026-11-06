"""Plot measured airfoil force sensitivity from checked source tables."""
from pathlib import Path
import argparse,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,ScalarFormatter
p=argparse.ArgumentParser()
p.add_argument('--data-dir',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
a=p.parse_args()
data=json.loads((a.data_dir/'comparison.json').read_text())
rows=data['records'];n=[r['cells'] for r in rows]
assert len(rows)==3 and not data['uniform_refinement']
if a.output_dir.exists():raise FileExistsError('Use a fresh versioned figure directory')
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':9.5,
 'axes.labelsize':10,'axes.linewidth':.65,'xtick.direction':'in','ytick.direction':'in',
 'xtick.top':True,'ytick.right':True,'lines.linewidth':1.05,'legend.fontsize':8.5,
 'legend.frameon':False,'pdf.fonttype':42,'svg.fonttype':'path'})
fig,axes=plt.subplots(1,2,figsize=(7.05,3.0),layout='constrained')
axes[0].plot(n,[r['Cd'] for r in rows],'o-',color='#0072B2',mfc='white',ms=4,mew=.8)
axes[0].set_ylabel(r'Drag coefficient, $C_D$')
axes[0].ticklabel_format(axis='y',style='plain',useOffset=False)
for key,label,color,marker,linestyle in [('pressure_Cd','Pressure','#0072B2','s','-'),
 ('viscous_Cd','Viscous','#D55E00','^','--'),('Cd','Total','.2','o',':')]:
 axes[1].plot(n,[r[key]-rows[-1][key] for r in rows],label=label,color=color,
  marker=marker,ls=linestyle,ms=4,mfc='white',mew=.8)
axes[1].axhline(0,color='.6',lw=.6)
axes[1].set_ylabel(r'$C_{D,i}-C_{D,i,\mathrm{fine}}$')
formatter=ScalarFormatter(useMathText=True);formatter.set_powerlimits((-3,3));formatter.set_useOffset(False)
axes[1].yaxis.set_major_formatter(formatter)
axes[1].yaxis.set_offset_position('right')
axes[1].legend(loc='best')
for ax,panel in zip(axes,['(a)','(b)']):
 ax.set_xscale('log');ax.set_xlim(9e4,1.1e6)
 ax.set_xticks([1e5,3e5,1e6],labels=[r'$10^5$',r'$3\times10^5$',r'$10^6$'])
 ax.xaxis.set_minor_formatter(NullFormatter())
 ax.set_xlabel(r'Number of fluid cells, $N$')
 ax.text(.03,1.04,panel,transform=ax.transAxes,va='bottom')
a.output_dir.mkdir(parents=True)
base=a.output_dir/'fig07-airfoil-mesh-sensitivity'
for extension in ['pdf','svg','png']:fig.savefig(base.with_suffix('.'+extension),dpi=600,bbox_inches='tight')
plt.close(fig)
caption=('Figure 7. Fluent airfoil force sensitivity at Re=100 and zero incidence on the same closed NACA0012 polygon '
 'and radius20c domain. The nonuniform triangular meshes have '+', '.join(f"{r['cells']:,}" for r in rows)+' cells. '
 '(a) Total drag coefficient. (b) Pressure, viscous and total drag-coefficient differences relative to the finest tested mesh; '
 'opposing changes in pressure and viscous contributions can cancel. Lines guide the eye and imply no fitted convergence order. '
 f"The medium-to-fine total Cd difference is {data['medium_fine_Cd_difference_percent']:.6g}% using the fine result as denominator. "
 'All three calculations passed their residual, force-stability, mass-balance and quality checks. '
 'Actual wall and wake spacings were independently measured. Fixed original polygon vertices and Gmsh subdivision quantization '
 'prevent a uniform-refinement interpretation. No statistical error bars, GCI, universal mesh independence, '
 'domain independence or physical validation are established. Fluent computed the forces; Codex prepared the checks and figure.\n')
(a.output_dir/'CAPTION.md').write_text(caption,encoding='utf-8')
manifest={name:hashlib.sha256((a.data_dir/name).read_bytes()).hexdigest() for name in ['comparison.json','comparison.csv']}
(a.output_dir/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
print(caption)
