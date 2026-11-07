"""Publish inspectable numerical histories for the accepted sloshing pair and retained misses."""
from pathlib import Path
import csv,hashlib,json,re,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('E:/CFD-Learning-Lab');P=R/'local/sloshing-production-v001';N=R/'native/runs/sloshing-production-v001'
OUT=R/'exports/figures/2026-11-07-sloshing-1p2s-paired-v003'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    targets=[OUT/'numerical-checks.png',OUT/'numerical-checks.pdf',OUT/'fluent-step-residuals.csv',OUT/'openfoam-step-residuals.csv',OUT/'saved-state-checks.csv',OUT/'column-water-depth.csv']
    assert not any(p.exists() for p in targets)
    fa_path=P/'fluent-production-audit-v001.json';oa_path=P/'openfoam-production-audit-v003.json'
    f=json.loads(fa_path.read_text());o=json.loads(oa_path.read_text())
    assert f['production_checks_pass'] and o['production_checks_pass']
    logpath=N/'fluent/solver-transcript.txt';assert sha(logpath)==f['native_source_sha256'][str(logpath)]
    text=logpath.read_text(errors='replace');num=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
    stops=list(re.finditer(r'Flow time = ('+num+r')s, time step = (\d+)',text))
    rowpattern=re.compile(r'^\s*(\d+)\s+('+num+r')\s+('+num+r')\s+('+num+r')\s+\d+:\d+:\d+\s+\d+\s*$',re.M)
    assert len(stops)==2400
    rows=[];start=0;last_iteration=0
    for step,stop in enumerate(stops,1):
        assert int(stop[2])==step and abs(float(stop[1])-step*.0005)<1e-10
        candidates=list(rowpattern.finditer(text,start,stop.start()));assert candidates
        last=candidates[-1];iteration=int(last[1]);values=[float(last[k]) for k in [2,3,4]]
        assert max(values)<=1e-4 and iteration>last_iteration
        rows.append([step,float(stop[1]),*values,iteration-last_iteration]);last_iteration=iteration;start=stop.end()
    with targets[2].open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['step','time_s','scaled_continuity','scaled_Ux','scaled_Uy','inner_iterations']);writer.writerows(rows)
    with targets[3].open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['step','time_s','last_outer_initial_Ux','last_outer_initial_Uy','last_outer_initial_p_rgh','alpha_min_final_outer','alpha_max_final_outer'])
        for step,row in enumerate(o['steps'],1):writer.writerow([step,row['time_s'],*[row['residuals'][k] for k in ['Ux','Uy','p_rgh']],row['alpha_min_final_outer'],row['alpha_max_final_outer']])
    with targets[4].open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['solver','step','time_s','water_area_m2','relative_area_drift','alpha_min','alpha_max','maximum_mixture_speed_m_s'])
        for label,a in [('Fluent',f),('OpenFOAM',o)]:
            for row in a['records']:writer.writerow([label,row['step'],row['time_s'],row['water_area_m2'],row['water_area_m2']/.01-1,row['alpha_min'],row['alpha_max'],row['max_speed_m_s']])
    depth_path=OUT/'column-depth-histories.json';depth=json.loads(depth_path.read_text())
    with targets[5].open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['solver','time_s','x_mm','column_depth_minus_50mm'])
        for j,label in enumerate(depth['solver_labels']):
            for i,t in enumerate(depth['times_s']):
                for k,x in enumerate(depth['x_mm']):writer.writerow([label,t,x,depth['depth_minus_50mm'][j][i][k]])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(12,8.5));fig.subplots_adjust(left=.09,right=.98,bottom=.13,top=.83,wspace=.28,hspace=.48)
    data=np.array(rows)
    for index,label in enumerate(['Continuity','x momentum','y momentum']):axes[0,0].semilogy(data[:,1],data[:,index+2],lw=.9,label=label)
    for name in ['Ux','Uy','p_rgh']:axes[0,1].semilogy([r['time_s'] for r in o['steps']],[r['residuals'][name] for r in o['steps']],lw=.9,label=name)
    axes[0,0].set(title='Fluent: terminal scaled residuals',ylabel='Scaled residual')
    axes[0,1].set(title='OpenFOAM: last-outer initial residuals',ylabel='Native normalized residual')
    for ax in axes[0]:ax.axhline(1e-4,color='#555555',ls='--',lw=1,label='Declared limit');ax.legend(fontsize=8,framealpha=.95,loc='lower left')
    for a,color,label in [(f,'#087eb5','Fluent'),(o,'#cc6423','OpenFOAM v003')]:axes[1,0].plot([r['time_s'] for r in a['records']],[(r['water_area_m2']/.01-1)*1e6 for r in a['records']],color=color,lw=1.5,label=label)
    axes[1,0].set(title='Saved-state water conservation',ylabel='Relative water-area change (ppm)');axes[1,0].legend(fontsize=8,framealpha=.95)
    prior_paths=[P/'openfoam-production-audit-v001-complete.json',P/'openfoam-production-audit-v002.json',oa_path]
    for path,label,color in zip(prior_paths,['Original v001: MISS','Tighter alpha v002: MISS','Tighter pressure v003: PASS'],['#888888','#cc6423','#087eb5']):
        a=json.loads(path.read_text());excess=np.array([max(0,r['alpha_max']-1) for r in a['records']]);mask=excess>0
        axes[1,1].semilogy(np.array([r['time_s'] for r in a['records']])[mask],excess[mask],color=color,lw=1.3,label=label)
    axes[1,1].axhline(1e-6,color='#555555',ls='--',lw=1,label='Unchanged bound tolerance');axes[1,1].set(title='OpenFOAM saved alpha overshoot',ylabel='Positive excess above alpha = 1');axes[1,1].legend(fontsize=8,framealpha=.95,loc='lower left')
    for ax in axes.flat:ax.set(xlabel='Physical time (s)',xlim=[0,1.2]);ax.grid(alpha=.2)
    fig.suptitle('Sloshing: numerical checks and retained misses',x=.09,ha='left',y=.97,fontsize=20,fontweight='bold')
    fig.text(.09,.90,'2,400 time steps in each full run; 121 saved states | residual definitions differ between solvers',fontsize=11)
    fig.text(.09,.052,'Residual values are not a cross-code accuracy ranking. Zero alpha excess is omitted from the logarithmic panel.',fontsize=9)
    fig.text(.09,.023,'Original attempts remain preserved. Numerical acceptance does not establish mesh/time independence or physical validation.',fontsize=9)
    fig.savefig(targets[0],dpi=170);fig.savefig(targets[1]);plt.close(fig)
    summary=dict(fluent_steps=2400,openfoam_steps=2400,saved_states_each=121,fluent_max_inner_iterations=max(r[-1] for r in rows),fluent_max_relative_water_area_drift=max(abs(r['water_area_m2']/.01-1) for r in f['records']),openfoam_max_relative_water_area_drift=o['maximum_relative_water_area_drift'],source_sha256={str(p):sha(p) for p in [fa_path,oa_path,logpath,depth_path,*prior_paths,Path(__file__)]},artifact_sha256={p.name:sha(p) for p in targets})
    (OUT/'numerical-figure-record.json').write_text(json.dumps(summary,indent=2)+'\n');shutil.copyfile(__file__,OUT/'executed-numerical-plotter.py')
    print(json.dumps({k:v for k,v in summary.items() if k not in ['source_sha256','artifact_sha256']}))
if __name__=='__main__':main()
