"""Plot actual accepted fine-grid Fluent qualification fields and residuals."""
from pathlib import Path
import argparse,hashlib,json
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyse_nozzle_fluent_v005 import common_mapping
from sloshing_refinement_common_v001 import cell_values,alpha_group
P=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['flat','tilted']);ap.add_argument('out',type=Path);args=ap.parse_args();assert not args.out.exists()
    audit=P/('fluent-'+args.mode+'-audit.json');a=json.loads(audit.read_text());assert a['numerical_checks_pass'] and a['cells']==12800 and a['case_mode']==args.mode
    for name,h in a['native_source_sha256'].items():assert sha(Path(name))==h,name
    run=Path('E:/CFD-Learning-Lab/native/runs/sloshing-refinement-v001/qualification')/('fluent-'+args.mode)
    with np.load(run/'geometry.npz') as g:fm=common_mapping(run/'tank.cas.h5',g)
    fields=[]
    for step in [0,40]:
        with h5py.File(run/f'tank-{step:04d}.dat.h5') as f:
            data=[]
            for group in [alpha_group(f),f['results/1/phase-1/cells/SV_U'],f['results/1/phase-1/cells/SV_V']]:
                raw=cell_values(group);aligned=np.empty(12800);aligned[fm]=raw;data.append(aligned.reshape(80,160))
            fields.append(data)
    initial,last=fields;x=(np.arange(160)+.5)*1.25;speed=np.hypot(last[1],last[2]);depth=[v[0].sum(axis=0)*1.25 for v in fields]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    im=ax[0,0].imshow(last[0],origin='lower',extent=[0,200,0,100],interpolation='nearest',vmin=0,vmax=1,cmap='Blues',aspect='equal');fig.colorbar(im,ax=ax[0,0],label='Water volume fraction');ax[0,0].set(title='Actual water fraction | 0.020 s',xlabel='x (mm)',ylabel='y (mm)')
    im=ax[0,1].imshow(np.ma.masked_where(last[0]<.5,speed),origin='lower',extent=[0,200,0,100],interpolation='nearest',vmin=0,vmax=1e-4 if args.mode=='flat' else None,cmap='viridis',aspect='equal');fig.colorbar(im,ax=ax[0,1],label='Mixture speed (m/s)');ax[0,1].set(title='Actual mixture speed | water fraction >= 0.5',xlabel='x (mm)',ylabel='y (mm)')
    ax[1,0].plot(x,depth[0]-50,color='#777777',ls='--',label='Initial');ax[1,0].plot(x,depth[1]-50,color='#176b9c',label='0.020 s');ax[1,0].set(title='Column-integrated water depth',xlabel='Column centre x (mm)',ylabel='Depth minus 50 mm (mm)',xlim=[0,200]);ax[1,0].legend()
    hist=a['step_records'];limits=a['extra']['residual_targets']
    for key,target in limits.items():ax[1,1].semilogy([r['time_s'] for r in hist],[max(r['residuals'][key]/target,1e-16) for r in hist],label=key)
    ax[1,1].axhline(1,color='black',ls=':',lw=1);ax[1,1].set(title='Terminal residual divided by its required limit',xlabel='Physical time (s)',ylabel='Residual / acceptance limit',xticks=[0,.005,.01,.015,.02],xlim=[0,.02]);ax[1,1].legend(fontsize=9)
    for v in ax[1]:v.grid(alpha=.2)
    fig.suptitle('ANSYS Fluent | '+args.mode+'-water qualification | 160 x 80 cells',fontsize=18)
    fig.get_layout_engine().set(rect=(0,.06,1,.96));fig.text(.035,.02,'40 actual steps; dt = 0.5 ms. This short numerical qualification does not establish full-run accuracy or physical validation.',fontsize=10)
    args.out.mkdir(parents=True);fig.savefig(args.out/'fluent-qualification.png',dpi=140);fig.savefig(args.out/'fluent-qualification.pdf');plt.close(fig)
    (args.out/'plot-receipt.json').write_text(json.dumps(dict(audit_sha256=sha(audit),numerical_checks_pass=True,case_mode=args.mode,actual_steps=[0,40],actual_time_s=[0,.02],physical_validation=False,code_sha256=sha(Path(__file__))),indent=2)+'\n')
    print(args.out/'fluent-qualification.png')
if __name__=='__main__':main()
