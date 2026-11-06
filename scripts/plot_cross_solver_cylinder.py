"""Academic comparison figure from accepted, matched-cylinder solver records."""
from pathlib import Path
import argparse,csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

p=argparse.ArgumentParser()
p.add_argument('--data-dir',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
a=p.parse_args()
if a.output_dir.exists():raise FileExistsError('Use a fresh versioned output directory')
a.output_dir.mkdir(parents=True)
data=json.loads((a.data_dir/'comparison.json').read_text())
assert data['numerical_checks_passed'] and data['cross_code_target_pass']
forces=list(csv.DictReader((a.data_dir/'openfoam-force-history.csv').open()))
residuals=list(csv.DictReader((a.data_dir/'openfoam-initial-residuals.csv').open()))

plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':9.5,
    'axes.labelsize':10,'axes.linewidth':.65,'xtick.labelsize':9,'ytick.labelsize':9,
    'xtick.direction':'in','ytick.direction':'in','xtick.top':True,'ytick.right':True,
    'xtick.major.width':.65,'ytick.major.width':.65,'axes.grid':False,
    'lines.linewidth':1.05,'legend.fontsize':8.5,'legend.frameon':False,
    'pdf.fonttype':42,'svg.fonttype':'path','savefig.facecolor':'white'})
fig,axes=plt.subplots(1,2,figsize=(7.05,3.0),layout='constrained')
late=[r for r in forces if int(r['iteration'])>=1200]
k=[int(r['iteration']) for r in late]
cd=[float(r['Cd']) for r in late]
axes[0].plot(k,cd,color='#0072B2',label='OpenFOAM')
axes[0].plot(k[-1],cd[-1],'o',ms=4,mfc='white',mec='#0072B2',mew=.8)
axes[0].axhline(data['fluent']['Cd'],color='.2',ls='--',lw=.9,label='Fluent terminal value')
axes[0].set(xlabel=r'OpenFOAM solver iteration, $k$',ylabel=r'Drag coefficient, $C_D$',
            xlim=(1200,3050),ylim=(2.0725,2.0797))
axes[0].set_xticks([1200,1800,2400,3000])
axes[0].legend(loc='upper right',handlelength=2.7)
for name,color,style in [('p','#0072B2','-'),('Ux','#D55E00','--'),('Uy','#009E73','-.')]:
    axes[1].plot([int(r['iteration']) for r in residuals],[float(r[name]) for r in residuals],
                 color=color,ls=style,label={'p':r'$p$','Ux':r'$U_x$','Uy':r'$U_y$'}[name])
axes[1].axhline(1e-7,color='.45',ls=':',lw=.8)
axes[1].axvline(1200,color='.55',ls=':',lw=.8)
axes[1].set_yscale('log')
axes[1].set(xlabel=r'OpenFOAM solver iteration, $k$',ylabel='Initial linear-solver residual',xlim=(0,3050),ylim=(1e-11,2))
axes[1].set_yticks([1,1e-3,1e-6,1e-9])
axes[1].set_xticks([0,1000,2000,3000])
axes[1].legend(loc='upper right',handlelength=2.7)
for label,ax in zip(['(a)','(b)'],axes):
    ax.text(.035,.955,label,transform=ax.transAxes,va='top',ha='left',fontsize=10)
name='fig04-fluent-openfoam-cylinder'
for extension in ['pdf','svg','png']:
    fig.savefig(a.output_dir/f'{name}.{extension}',dpi=600)
plt.close(fig)
caption=f'''# Figure 4. Matched Fluent–OpenFOAM cylinder calculation at Re = 20

(a) The continuation of OpenFOAM drag convergence after the initial
1,200-iteration cap; the earlier startup interval is omitted here but retained
in the source CSV. The dashed line is the terminal Fluent value, not a Fluent
iteration history. (b) OpenFOAM initial linear-solver residuals for pressure and
the two in-plane velocity components. The vertical dotted line marks the
initial cap and the horizontal dotted line the declared 1e-7 residual target.
These residuals do not share Fluent's scaling and are not overlaid with its
residuals as though they were equivalent error norms.

The same 98,304-cell O-grid was used, with outer radius R/D = 20. OpenFOAM used
a one-cell, unit-span extrusion with empty front/back patches; converted x/y
coordinates agreed with the Fluent source mesh to the stored precision. Both
cases were steady, incompressible and laminar, with D = 1 m, U = 1 m/s,
rho = 1 kg/m3, mu = 0.05 Pa s and unit reference area. The upstream velocity
inlet, downstream zero-pressure outlet, no-slip cylinder and force directions
were matched. Numerical operators differ: coupled second-order Fluent versus
steady SIMPLE and bounded linear-upwind OpenFOAM.

Terminal Cd was {data['fluent']['Cd']:.7f} in Fluent and
{data['openfoam']['Cd']:.7f} in OpenFOAM; the relative difference normalized by
Fluent was {data['cross_code_Cd_difference_percent']:.5f}%. OpenFOAM's final
100-iteration Cd range was {data['openfoam']['last100_Cd_relative_range']*100:.6f}%
of its window mean. Meeting the declared 0.5% cross-code criterion is numerical
agreement for this confined setup, not physical validation or an uncertainty
estimate. Fluent reused a saved converged state; OpenFOAM started from uniform
flow. Their iteration counts and run times are not a performance benchmark.

Codex prepared the cases, analysis, figure and caption. ANSYS Fluent and
OpenFOAM computed their respective solutions. The owner requested the paired
workflow and academic presentation; independent owner analysis is not inferred.
'''
(a.output_dir/'CAPTION.md').write_text(caption,encoding='utf-8')
manifest={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in a.data_dir.iterdir() if f.is_file()}
(a.output_dir/'source-manifest.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(dict(figure=name,formats=['pdf','svg','png'],difference_percent=data['cross_code_Cd_difference_percent']),indent=2))
