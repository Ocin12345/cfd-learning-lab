"""Plot the final Fluent state only after its complete native numerical audit passes."""
from pathlib import Path
import argparse,csv,hashlib,json,re
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ansys.fluent.core.offline import _lispy
from analyse_nozzle_fluent_v005 import common_mapping
from sloshing_refinement_common_v001 import cell_values,alpha_group
from audit_refined_run_v001 import NUM
R=Path('E:/CFD-Learning-Lab');RUN=R/'native/runs/sloshing-refinement-v001/production/fluent'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('step',type=int);ap.add_argument('out',type=Path);args=ap.parse_args();step=args.step;assert step>0 and step%20==0 and step==2400 and not args.out.exists()
    audit_path=Path(__file__).resolve().parent/'fluent-production-audit.json';audit=json.loads(audit_path.read_text());assert audit['production_complete'] and audit['numerical_checks_pass'] and audit['steps']==2400 and audit['saved_states']==121
    for name,h in audit['native_source_sha256'].items():assert sha(Path(name))==h,name
    from wing_fluent_v001 import available_memory
    assert available_memory()/2**30>=4
    setup=json.loads((RUN/'setup.json').read_text());assert setup['cells']==[160,80] and setup['case_mode']=='production' and setup['time_step_s']==.0005
    with (RUN/'frames.csv').open(newline='') as f:ledger=list(csv.DictReader(f))
    row=next(r for r in ledger if int(r['step'])==step);assert abs(float(row['time_s'])-step*.0005)<1e-10
    with np.load(RUN/'geometry.npz') as g:fm=common_mapping(RUN/'tank.cas.h5',g)
    paths=[RUN/'tank.cas.h5',RUN/'geometry.npz',RUN/'setup.json'];data=[]
    for selected in [0,step]:
        path=RUN/f'tank-{selected:04d}.dat.h5';paths.append(path)
        with h5py.File(path) as f:
            clock={v[0]:v[1] for v in _lispy.parse(f['settings/Data Variables'][0].decode())[1]};assert int(clock['time-step'])==selected and abs(float(clock['flow-time'])-selected*.0005)<1e-10
            assert abs(float(clock['physical-time-step'])-.0005)<1e-14
            values=[]
            for group in [alpha_group(f),f['results/1/phase-1/cells/SV_U'],f['results/1/phase-1/cells/SV_V']]:
                raw=cell_values(group);aligned=np.empty(12800);aligned[fm]=raw;values.append(aligned.reshape(80,160))
            assert np.isfinite(cell_values(f['results/1/phase-1/cells/SV_P'])).all()
            a=values[0];assert a.min()>=-1e-6 and a.max()<=1+1e-6 and abs(a.sum()*.00125**2/.01-1)<=.005
            data.append(values)
    text=(RUN/'solver-transcript.txt').read_text(errors='replace');bounds=list(re.finditer(r'Flow time = ('+NUM+r')s, time step = (\d+)',text));assert len(bounds)>=step;end=bounds[step-1];assert int(end[2])==step and abs(float(end[1])-step*.0005)<1e-10
    prefix=text[:end.end()];pattern=re.compile(r'^\s*(\d+)\s+('+NUM+r')\s+('+NUM+r')\s+('+NUM+r')\s+\d+:\d+:\d+\s+\d+\s*$',re.M);hist=[];start=0;previous=0
    for k,b in enumerate(bounds[:step],1):
        assert int(b[2])==k;matches=list(pattern.finditer(prefix,start,b.start()));assert matches;last=matches[-1];values=[float(last[i]) for i in [2,3,4]];assert np.isfinite(values).all() and max(values)<=1e-4
        iteration=int(last[1]);assert 0<iteration-previous<=400;hist.append(values);start=b.end();previous=iteration
    initial,last=data;x=(np.arange(160)+.5)*1.25;speed=np.hypot(last[1],last[2]);depth=[v[0].sum(axis=0)*1.25-50 for v in data];t=np.arange(1,step+1)*.0005
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for j,values in enumerate([last[0],np.ma.masked_where(last[0]<.5,speed)]):
        im=ax[0,j].imshow(values,origin='lower',extent=[0,200,0,100],interpolation='nearest',vmin=0,vmax=1 if j==0 else None,cmap='Blues' if j==0 else 'viridis',aspect='equal');fig.colorbar(im,ax=ax[0,j],label='Water volume fraction' if j==0 else 'Mixture speed (m/s)');ax[0,j].set(xlabel='x (mm)',ylabel='y (mm)',title='Water fraction' if j==0 else 'Mixture speed | water fraction >= 0.5')
    ax[1,0].plot(x,depth[0],ls='--',color='#777777',label='Initial');ax[1,0].plot(x,depth[1],color='#176b9c',label=f'{step*.0005:.2f} s');ax[1,0].set(title='Column-integrated water depth',xlabel='Column centre x (mm)',ylabel='Depth minus 50 mm (mm)',xlim=[0,200]);ax[1,0].legend()
    for i,label in enumerate(['continuity','x-velocity','y-velocity']):ax[1,1].semilogy(t,np.maximum(np.asarray(hist)[:,i]/1e-4,1e-16),label=label)
    ax[1,1].axhline(1,color='black',ls=':',lw=1);ax[1,1].set(title='Terminal residual divided by required limit',xlabel='Physical time (s)',ylabel='Residual / 1e-4',xlim=[0,step*.0005]);ax[1,1].legend(fontsize=9)
    for a in ax[1]:a.grid(alpha=.18)
    fig.suptitle(f'ANSYS Fluent | final numerical state at {step*.0005:.2f} s | 12,800 cells',fontsize=18)
    fig.get_layout_engine().set(rect=(0,.07,1,.87));fig.text(.035,.025,'2,400 steps and 121 native states pass the full numerical audit. Mesh/time independence and physical validation remain open.',fontsize=10)
    args.out.mkdir(parents=True);fig.savefig(args.out/'fluent-final.png',dpi=140);fig.savefig(args.out/'fluent-final.pdf');plt.close(fig)
    (args.out/'solver-prefix.txt').write_text(prefix,encoding='utf-8')
    record=dict(actual_step=step,actual_time_s=step*.0005,planned_final_time_s=1.2,selected_state_finite_bounds_and_area_pass=True,logged_steps_meeting_residual_limits=step,water_area_m2=float(last[0].sum()*.00125**2),full_run_accepted=True,numerical_audit_sha256=sha(audit_path),physical_validation=False,native_source_sha256={str(p):sha(p) for p in paths},solver_prefix_sha256=sha(args.out/'solver-prefix.txt'),code_sha256=sha(Path(__file__)))
    (args.out/'final-receipt.json').write_text(json.dumps(record,indent=2)+'\n');print(args.out/'fluent-final.png')
if __name__=='__main__':main()
