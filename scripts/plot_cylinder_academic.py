"""Publication-style figures from existing, inspectable Fluent records.

This changes presentation only. It does not rerun CFD, infer uncertainty, or
overwrite the earlier figures. Captions accompany each standalone vector figure.
"""
from pathlib import Path
import argparse,csv,hashlib,json,re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter,FixedLocator,FuncFormatter,NullFormatter
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument('--output-dir',type=Path,required=True)
args=parser.parse_args()
out=args.output_dir.resolve()
if out.exists():raise FileExistsError('Use a new versioned figure directory')
out.mkdir(parents=True)
repo=Path(__file__).resolve().parents[1]

plt.rcParams.update({
    'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':9.5,
    'axes.labelsize':10,'axes.titlesize':10,'axes.linewidth':.65,
    'xtick.labelsize':9,'ytick.labelsize':9,'xtick.direction':'in','ytick.direction':'in',
    'xtick.top':True,'ytick.right':True,'xtick.major.size':3,'ytick.major.size':3,
    'xtick.major.width':.65,'ytick.major.width':.65,'axes.grid':False,
    'lines.linewidth':1.05,'lines.markersize':4.2,'legend.fontsize':8.5,
    'legend.frameon':False,'figure.facecolor':'white','savefig.facecolor':'white',
    'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'path',
})
blue='#0072B2'; orange='#D55E00'; green='#009E73'
sources={}
def read_source(relative):
    p=repo/relative
    sources[relative]=hashlib.sha256(p.read_bytes()).hexdigest()
    return p.read_text()
def panels(axes):
    for tag,ax in zip('abcdefghijklmnopqrstuvwxyz',axes):
        ax.text(.035,.955,f'({tag})',transform=ax.transAxes,va='top',ha='left',fontsize=10)
def save(fig,name):
    for suffix in ['pdf','svg','png']:
        fig.savefig(out/f'{name}.{suffix}',dpi=600)
    plt.close(fig)

comparison=json.loads(read_source('results/lesson02-domain-comparison.json'))
rows=comparison['rows']
radius=np.array([r['outer_radius_m'] for r in rows])
cd=np.array([r['drag_coefficient'] for r in rows])
changes=np.array([r['change_from_preceding_percent'] for r in rows[1:]])
fig,axes=plt.subplots(1,2,figsize=(7.05,3.0),layout='constrained')
axes[0].plot(radius,cd,'o-',color=blue,markerfacecolor='white',markeredgewidth=.9)
axes[0].set_xscale('log',base=2)
axes[0].set_xticks(radius,labels=[str(int(r)) for r in radius])
axes[0].set(xlabel=r'Outer-boundary radius, $R_{\infty}/D$',ylabel=r'Drag coefficient, $C_D$',ylim=(1.995,2.082))
axes[0].yaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:.2f}'))
axes[1].plot(radius[1:],changes,'s-',color=orange,markerfacecolor='white',markeredgewidth=.9,label='Successive-domain change')
axes[1].axhline(.2,color='.25',ls='--',lw=.85,label='Declared criterion: 0.2%')
axes[1].set_xscale('log',base=2)
axes[1].set_yscale('log')
axes[1].set_xticks(radius[1:],labels=[str(int(r)) for r in radius[1:]])
axes[1].set_yticks([.1,.2,.5,1,2],labels=['0.1','0.2','0.5','1','2'])
axes[1].yaxis.set_minor_formatter(NullFormatter())
axes[1].set(xlabel=r'Outer-boundary radius, $R_{\infty}/D$',ylabel=r'Successive change, $\varepsilon_R$ (%)',ylim=(.1,3.1))
axes[1].legend(loc='upper right',bbox_to_anchor=(1.01,.99),handlelength=2)
panels(axes)
save(fig,'fig01-cylinder-domain-sensitivity')

history=list(csv.DictReader(read_source('results/lesson01-cylinder-replay/force-history.csv').splitlines()))
k=np.array([int(r['additional_iteration']) for r in history])
drag=np.array([float(r['drag_coefficient']) for r in history])
lift=np.array([float(r['lift_coefficient']) for r in history])
fig,axes=plt.subplots(1,2,figsize=(7.05,3.0),layout='constrained')
axes[0].step(k,(drag-drag[0])*1e7,where='mid',color=blue)
axes[0].set(xlabel=r'Additional solver iterations, $\Delta k$',ylabel=r'Drag increment, $10^7\Delta C_D$',xlim=(-2,102),ylim=(-.3,8.2))
axes[1].plot(k,np.maximum(np.abs(lift),1e-18),color=orange)
axes[1].set_yscale('log')
axes[1].set_yticks([1e-6,1e-9,1e-12,1e-15])
axes[1].set(xlabel=r'Additional solver iterations, $\Delta k$',ylabel=r'Absolute lift coefficient, $|C_L|$',xlim=(-2,102),ylim=(1e-16,2e-6))
panels(axes)
save(fig,'fig02-cylinder-force-stability')

residual_text=read_source('results/lesson01-cylinder-replay/residuals.xy')
blocks=re.findall(r'\(\(xy/key/label "([^"]+)"\)\s*([^)]*)\)',residual_text)
assert len(blocks)==3,'Unexpected residual export format'
fig,ax=plt.subplots(figsize=(3.46,2.9),layout='constrained')
styles=[(blue,'-',r'Continuity'),(orange,'--',r'$x$ momentum'),(green,'-.',r'$y$ momentum')]
for (name,block),(colour,style,label) in zip(blocks,styles):
    values=np.array([[float(v) for v in line.split()] for line in block.strip().splitlines()])
    assert values.shape[1]==2 and len(values)>=151
    ax.plot(values[:,0],values[:,1],color=colour,ls=style,label=label)
ax.axhline(1e-7,color='.4',ls=':',lw=.8)
ax.axvline(51,color='.55',ls=':',lw=.8)
ax.set_yscale('log')
ax.set(xlabel=r'Solver iteration, $k$',ylabel='Scaled residual',xlim=(0,155),ylim=(1e-18,2))
ax.set_yticks([1,1e-3,1e-6,1e-9,1e-12,1e-15,1e-18])
ax.legend(loc='upper right',handlelength=2.6)
save(fig,'fig03-cylinder-scaled-residuals')

captions=r'''# Figure captions

**Figure 1. Cylinder domain-size sensitivity at Re = 20.** (a) Drag coefficient
as a function of the outer-boundary radius. (b) Successive relative difference,
epsilon_R = 100 |C_D(R_i) - C_D(R_{i-1})| / |C_D(R_i)|. The dashed line denotes
the 0.2% criterion declared before the subsequent comparisons. The original
fine-grid coordinates through R/D = 20 were retained exactly; additional outer
rings were generated for larger domains. Each new case reached the 1e-8 residual
target and was extended by 100 iterations. Markers denote computed values;
connecting lines are guides to the eye. Meeting this criterion for drag does
not establish an error bound, convergence of all fields, or physical validation.

**Figure 2. Iterative force stability of the saved Re = 20 cylinder calculation.**
The saved solution at k = 51 was extended by 100 iterations on 98,304 cells with
R/D = 20. (a) Drag increment relative to the saved state, multiplied by 10^7.
(b) Magnitude of the lift coefficient. The forces were sampled every iteration
from Fluent reports with approximately eight significant digits; the staircase
in (a) and zero last-window drag range reflect that reporting precision. Very
small late lift values are numerical, rather than evidence of physical flow
oscillations. These are solver iterations, not physical time steps.

**Figure 3. Scaled residual histories for the Re = 20 cylinder replay.**
The vertical dotted line marks the saved state at k = 51; the horizontal dotted
line marks the original 1e-7 stopping criterion. Residual stopping checks were
disabled during the 100-iteration observation window while residual monitoring
continued. Fluent's residuals are shown using their exported scaling and cannot
be compared numerically to OpenFOAM residuals without reconciling definitions.

**Common conditions.** Steady planar incompressible laminar flow, D = 1 m,
U_infinity = 1 m/s, rho = 1 kg/m^3, mu = 0.05 Pa s, unit span and reference
area 1 m^2; no-slip cylinder; upstream velocity inlet and downstream zero-gauge
pressure outlet. Fluent Student 2026 R1, two solver processes, pressure-based
coupled method and second-order pressure/momentum discretization. No uncertainty
bars have been drawn because a defensible uncertainty estimate is not available.

**Attribution.** Codex prepared the plotting code and captions. ANSYS Fluent
computed the flow and exported the numerical reports. The owner requested the
study and the revision to academic presentation; no independent owner analysis
or research novelty is inferred.
'''
(out/'CAPTIONS.md').write_text(captions,encoding='utf-8')
(out/'source-manifest.json').write_text(json.dumps(dict(sources=sources,matplotlib=matplotlib.__version__,presentation_only=True),indent=2))
print(json.dumps(dict(output_directory=str(out),figures=3,formats=['pdf','svg','png'],data_unchanged=True),indent=2))
