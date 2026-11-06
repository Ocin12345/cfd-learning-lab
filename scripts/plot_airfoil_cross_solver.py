"""Plot saved matched-airfoil histories, including explicitly labelled failed diagnostics."""
from pathlib import Path
import argparse,csv,hashlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=argparse.ArgumentParser()
p.add_argument('--data-dir',type=Path,required=True)
p.add_argument('--output-dir',type=Path,required=True)
p.add_argument('--diagnostic',action='store_true',help='Explicitly plot a retained failed numerical record')
a=p.parse_args()
record=json.loads((a.data_dir/'comparison.json').read_text())
accepted=record['numerical_checks_passed']
if not accepted and not (a.diagnostic and record.get('failed_checks')):
 raise ValueError('Failed checks require explicit --diagnostic and recorded failed_checks')
with (a.data_dir/'openfoam-force-history.csv').open(newline='') as f:forces=list(csv.DictReader(f))
with (a.data_dir/'openfoam-initial-residuals.csv').open(newline='') as f:residuals=list(csv.DictReader(f))
targets=json.loads((a.data_dir/'declared-comparison.json').read_text())
fluent_cd=record['fluent']['Cd'];target=record['cross_code_learning_target_percent']
plt.rcParams.update({'font.family':'STIXGeneral','mathtext.fontset':'stix','font.size':9.5,
 'axes.labelsize':10,'axes.linewidth':.65,'xtick.direction':'in','ytick.direction':'in',
 'xtick.top':True,'ytick.right':True,'lines.linewidth':1.05,'legend.fontsize':8.5,
 'legend.frameon':False,'pdf.fonttype':42,'svg.fonttype':'path'})
fig=plt.figure(figsize=(7.05,5.4),layout='constrained')
gs=fig.add_gridspec(2,2)
axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
k=[int(r['iteration']) for r in forces]
difference=[100*(float(r['Cd'])-fluent_cd)/abs(fluent_cd) for r in forces]
axes[0].plot(k,difference,color='#0072B2')
axes[0].set_yscale('symlog',linthresh=.1)
axes[0].axhline(0,color='.35',ls='--',lw=.7)
for sign in [-1,1]:axes[0].axhline(sign*target,color='.55',ls=':',lw=.7)
axes[0].plot(k[-1],difference[-1],'o',color='#0072B2',mfc='white',ms=4,mew=.8)
axes[0].set_ylabel(r'Drag difference from Fluent [%]')
axes[1].plot(k,[float(r['Cl']) for r in forces],color='#D55E00')
axes[1].set_yscale('symlog',linthresh=1e-3)
axes[1].axhline(0,color='.35',ls='--',lw=.7)
for sign in [-1,1]:axes[1].axhline(sign*targets['zero_angle_abs_lift_diagnostic_target'],color='.55',ls=':',lw=.7)
axes[1].set_ylabel(r'OpenFOAM lift coefficient, $C_L$')
for name,label,color,ls in [('p',r'$p$','#0072B2','-'),('Ux',r'$U_x$','#D55E00','--'),('Uy',r'$U_y$','#009E73','-.')]:
 axes[2].semilogy([int(r['iteration']) for r in residuals],[float(r[name]) for r in residuals],label=label,color=color,ls=ls)
axes[2].axhline(targets['initial_p_U_residual_learning_target'],color='.5',ls=':',lw=.7)
axes[2].set_ylabel('Initial linear-solver residual');axes[2].legend(loc='upper right')
for ax,panel in zip(axes,['(a)','(b)','(c)']):
 ax.set_xlabel(r'OpenFOAM steady solver iteration, $k$');ax.set_xlim(0,k[-1]*1.015)
 ax.text(.02,1.04,panel,transform=ax.transAxes,va='bottom')
 if record.get('method_continuation_stage02'):ax.axvline(2000,color='.55',ls=':',lw=.7)
if not accepted:fig.suptitle('DIAGNOSTIC - checks failed',fontsize=11,color='#9c281f')
fig.canvas.draw();fig.set_layout_engine('none')
if a.output_dir.exists():raise FileExistsError('Use a fresh figure directory')
a.output_dir.mkdir(parents=True)
base=a.output_dir/'fig08-airfoil-cross-solver'
for extension in ['pdf','svg','png']:fig.savefig(base.with_suffix('.'+extension),dpi=600)
plt.close(fig)
status=('Both saved numerical solutions passed their checks.' if accepted else
 'DIAGNOSTIC: OpenFOAM failed checks: '+', '.join(record['failed_checks'])+'. This figure is not an accepted comparison.')
caption=(f"Figure 8. Matched Fluent/OpenFOAM airfoil model at Re=100 and zero incidence on the same {record['cells']:,}-cell "
 f"triangular mesh and radius {targets['outer_radius_m']/targets['chord_m']:g}c domain, with unit density, speed, chord, reference area and span. "
 '(a) Signed difference 100*(Cd_OF-Cd_Fluent)/abs(Cd_Fluent), on a symmetric-logarithmic axis '
 'linear within 0.1 percentage points. Zero is the terminal Fluent reference; dotted lines show '
 f'the declared +/-{target:g}% code-comparison criterion. '
 '(b) OpenFOAM Cl, with a symmetric-logarithmic axis linear within +/-0.001; dotted lines mark the '
 'declared zero-angle lift diagnostic. Full startup histories are included in both panels. '
 '(c) OpenFOAM initial linear-solver residuals, with pressure the maximum across all correction solves '
 'within each iteration. The horizontal dotted line is the unchanged residual criterion. '
 'The vertical line marks U relaxation 0.7 to 0.9 after iteration 2000; discretization and acceptance '
 'criteria remain unchanged. '
 f"Final absolute drag difference is {record['cross_code_Cd_difference_percent']:.6g}%. "
 'OpenFOAM starts from uniform fields. Fluent uses its own converged solution followed by a further '
 '100-iteration observation window. Residual norms and initialization differ, so this is neither a speed '
 'nor residual-magnitude ranking. '+status+' Code agreement is not physical validation.\n')
(a.output_dir/'CAPTION.md').write_text(caption,encoding='utf-8')
(a.output_dir/'source-manifest.json').write_text(json.dumps(dict(
 numerical_checks_passed=accepted,diagnostic=not accepted,failed_checks=record.get('failed_checks',[]),
 sources={name:hashlib.sha256((a.data_dir/name).read_bytes()).hexdigest() for name in
 ['comparison.json','declared-comparison.json','openfoam-force-history.csv','openfoam-initial-residuals.csv']}),indent=2)+'\n',encoding='utf-8')
(a.output_dir/'executed-plotter.py').write_bytes(Path(__file__).read_bytes())
print(caption)
