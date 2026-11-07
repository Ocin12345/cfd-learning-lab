"""Academic plots and CSVs from fully accepted native numerical audits."""
from pathlib import Path
import argparse,csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('out',type=Path);args=ap.parse_args();assert not args.out.exists()
    paths=[P/'fluent-production-audit.json',P/'openfoam-production-audit.json'];audits=[json.loads(p.read_text()) for p in paths]
    for a in audits:
        assert a['production_complete'] and a['numerical_checks_pass'] and a['steps']==2400 and a['saved_states']==121
        for name,h in a['native_source_sha256'].items():assert sha(Path(name))==h,name
    labels=['ANSYS Fluent','OpenFOAM 14'];colors=['#176b9c','#ce652b'];styles=['-','--'];summary={}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(left=.09,right=.96,bottom=.16,top=.83,hspace=.48,wspace=.30)
    for j,a in enumerate(audits):
        steps=a['step_records'];records=a['records'];t=np.array([r['time_s'] for r in steps]);tf=np.array([r['time_s'] for r in records]);targets=a['extra']['residual_targets'];ratio=np.array([max(r['residuals'][k]/v for k,v in targets.items()) for r in steps]);assert ratio.max()<=1
        area=np.array([abs(r['water_area_m2']/.01-1) for r in records]);saved=np.array([max(0.,-r['alpha_min'],r['alpha_max']-1) for r in records])
        ax[0,0].semilogy(t,np.maximum(ratio,1e-16),styles[j],color=colors[j],lw=1.1,label=labels[j]);ax[0,1].semilogy(tf,np.maximum(area,1e-16),styles[j],color=colors[j],lw=1.5)
        ax[1,0].semilogy(tf,np.maximum(saved,1e-16),styles[j],color=colors[j],lw=1.5)
        if j==1:
            terminal=np.array([max(0.,-r['alpha_min'],r['alpha_max']-1) for r in steps]);ax[1,0].semilogy(t,np.maximum(terminal,1e-16),':',color='#7544a0',lw=.9,label='OpenFOAM terminal alpha (all steps)')
        summary[labels[j]]=dict(steps=len(steps),saved_states=len(records),all_numerical_checks_pass=True,maximum_residual_limit_ratio=float(ratio.max()),maximum_saved_relative_water_area_drift=float(area.max()),maximum_saved_alpha_bound_exceedance=float(saved.max()),peak_interface_courant=a['peak_interface_courant'])
    for a,limit in [(ax[0,0],1),(ax[0,1],.005),(ax[1,0],1e-6)]:a.axhline(limit,color='black',ls=':',lw=1)
    ax[0,0].set(title='Worst terminal equation residual / required limit',xlabel='Physical time (s)',ylabel='Residual / its code-specific limit',xlim=[0,1.2]);ax[0,1].set(title='Saved water-area conservation',xlabel='Physical time (s)',ylabel='Absolute relative water-area drift',xlim=[0,1.2]);ax[1,0].set(title='Departure of alpha from [0,1]',xlabel='Physical time (s)',ylabel='Maximum bound exceedance',xlim=[0,1.2]);ax[1,0].legend(loc='upper right',fontsize=8)
    peaks=[a['peak_interface_courant'] for a in audits];bars=ax[1,1].bar(labels,peaks,color=colors,width=.55);ax[1,1].axhline(.25,color='black',ls=':',lw=1);ax[1,1].set(title='Peak reported VOF/interface Courant number',ylabel='Courant number',ylim=[0,.29])
    for bar,value in zip(bars,peaks):ax[1,1].text(bar.get_x()+bar.get_width()/2,value+.006,f'{value:.5f}',ha='center',fontsize=10)
    for a in ax.ravel():a.grid(alpha=.18);a.set_axisbelow(True)
    fig.suptitle('Sloshing | complete fine-grid numerical checks',x=.09,ha='left',y=.97,fontsize=20,fontweight='bold');fig.text(.09,.922,'12,800 cells | dt = 0.5 ms | 2,400 steps and 121 saved native states per solver',fontsize=12)
    handles,ls=ax[0,0].get_legend_handles_labels();fig.legend(handles,ls,loc='center',bbox_to_anchor=(.5,.875),ncol=2,frameon=False)
    fig.text(.09,.092,'Black dotted lines are acceptance limits. Exact zeros are shown at 1e-16 on logarithmic axes.',fontsize=9)
    fig.text(.09,.063,'Residual normalizations differ between codes. Passing these checks does not establish mesh/time independence or physical validation.',fontsize=9)
    args.out.mkdir(parents=True);fig.savefig(args.out/'numerical-checks.png',dpi=150);fig.savefig(args.out/'numerical-checks.pdf');plt.close(fig)
    for label,a in zip(['fluent','openfoam'],audits):
        rows=a['step_records'];keys=list(rows[0]['residuals'])
        with (args.out/(label+'-step-residuals.csv')).open('w',newline='') as stream:
            w=csv.writer(stream);w.writerow(['step','time_s',*keys,'passed'])
            for r in rows:w.writerow([r['step'],r['time_s'],*[r['residuals'][k] for k in keys],r['passed']])
    with (args.out/'saved-state-checks.csv').open('w',newline='') as stream:
        keys=list(audits[0]['records'][0]);w=csv.writer(stream);w.writerow(['solver',*keys])
        for label,a in zip(labels,audits):
            for r in a['records']:assert list(r)==keys;w.writerow([label,*r.values()])
    record=dict(summary=summary,numerical_checks_pass=True,physical_validation=False,audit_sha256={str(p):sha(p) for p in paths},code_sha256=sha(Path(__file__)))
    (args.out/'numerical-summary.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
