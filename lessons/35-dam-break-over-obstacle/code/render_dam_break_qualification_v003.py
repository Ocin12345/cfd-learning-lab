"""Render actual accepted short qualification fields on the nonuniform mesh."""
from pathlib import Path
import argparse,hashlib,json
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter
from analyse_nozzle_fluent_v005 import common_mapping
from dam_break_common_v001 import *
R=Path('E:/CFD-Learning-Lab');P=R/'local/dam-break-over-obstacle-v001';N=R/'native/runs/dam-break-over-obstacle-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['flat','collapse']);ap.add_argument('out',type=Path);args=ap.parse_args();assert not args.out.exists() and args.out.resolve().is_relative_to((R/'exports/figures').resolve())
    from wing_fluent_v001 import available_memory
    assert available_memory()/2**30>=4,'Insufficient analysis RAM; leave solver running and retry later.'
    apaths=[P/('fluent-'+args.mode+'-audit.json'),P/('openfoam-'+args.mode+'-audit.json')];audits=[json.loads(p.read_text()) for p in apaths]
    for a in audits:
        assert a['numerical_checks_pass'] and a['steps']==200 and a['saved_states']==21 and a['mode']==args.mode
        for name,h in a['native_source_sha256'].items():assert sha(Path(name))==h,name
    f=N/('qualification/fluent-flat-v003' if args.mode=='flat' else 'qualification/fluent-collapse-v001');o=N/('qualification/openfoam-flat-v002' if args.mode=='flat' else 'qualification/openfoam-collapse-v002')
    g=geometry();fm=common_mapping(f/'dam.cas.h5',g);om=g['openfoam_to_common'];data=[];sources=[]
    def proven(path,j):assert audits[j]['native_source_sha256'][str(path)]==sha(path);sources.append(path)
    p=f/'dam-0200.dat.h5';proven(p,0)
    with h5py.File(p) as h:
        groups=[phase['cells/SV_VOF'] for phase in h['results/1'].values() if 'cells/SV_VOF' in phase];water=groups[0] if len(groups)==1 else next(v for v in groups if '/phase-3/' in v.name)
        aligned=[]
        for group in [water,h['results/1/phase-1/cells/SV_U'],h['results/1/phase-1/cells/SV_V']]:
            raw=cell_values(group);v=np.empty(NCELLS);v[fm]=raw;aligned.append(v)
        data.append(aligned)
    p=o/'0.02/alpha.water';proven(p,1);alpha=native_field(p,1).ravel();p=o/'0.02/U';proven(p,1);velocity=native_field(p,3);aligned=[]
    for raw in [alpha,velocity[:,0],velocity[:,1]]:
        v=np.empty(NCELLS);v[om]=raw;aligned.append(v)
    data.append(aligned);assert np.isfinite(np.asarray(data)).all()
    maximum=max(float(np.max(np.where(d[0]>=.5,np.hypot(d[1],d[2]),0))) for d in data)
    vmax=max(1e-4,maximum) if args.mode=='flat' else maximum;assert vmax>0
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42})
    fig,ax=plt.subplots(2,2,figsize=(11,11));fig.subplots_adjust(left=.075,right=.86,bottom=.12,top=.83,hspace=.35,wspace=.27);labels=['ANSYS Fluent','OpenFOAM 14'];images=[]
    for j,d in enumerate(data):
        speed=np.hypot(d[1],d[2])
        for row,values in enumerate([d[0],np.ma.masked_where(d[0]<.5,speed)]):
            field=np.ma.masked_all((NY,NX));ij=g['cell_ij'];field[ij[:,1],ij[:,0]]=values
            im=ax[row,j].pcolormesh(g['x_edges']*1000,g['y_edges']*1000,field,shading='flat',cmap='Blues' if row==0 else 'viridis',vmin=0,vmax=1 if row==0 else vmax)
            if j==0:images.append(im)
            ax[row,j].add_patch(Rectangle((OB_LEFT*1000,0),(OB_RIGHT-OB_LEFT)*1000,OB_HEIGHT*1000,facecolor='#8d969e',edgecolor='#454b50',zorder=4))
            ax[row,j].set(aspect='equal',xlim=[0,584],ylim=[0,584],xlabel='x (mm)',ylabel='y (mm)',xticks=[0,146,292,438,584],yticks=[0,146,292,438,584],title=labels[j]+(' | water fraction' if row==0 else ' | mixture speed'))
        selection=(g['cell_ij'][:,0]%8==0)&(g['cell_ij'][:,1]%8==0)&(d[0]>=.5)&(speed>=vmax*.005)
        if selection.any():ax[1,j].quiver(g['centres'][selection,0]*1000,g['centres'][selection,1]*1000,d[1][selection],d[2][selection],angles='xy',scale_units='xy',scale=vmax/35,color='white',width=.004,minlength=0)
    for row,label in [(0,'Water volume fraction'),(1,'Mixture speed (m/s), alpha >= 0.5')]:
        pos=ax[row,1].get_position();cb=fig.add_axes([.90,pos.y0,.015,pos.height]);colorbar=fig.colorbar(images[row],cax=cb);colorbar.set_label(label)
        if row==1 and vmax<.001:
            colorbar.formatter=ScalarFormatter(useMathText=True);colorbar.formatter.set_powerlimits((0,0));colorbar.update_ticks()
    label='quiet-water control' if args.mode=='flat' else 'released-column qualification'
    fig.suptitle('Dam break | '+label,x=.075,ha='left',y=.98,fontsize=19,fontweight='bold')
    fig.text(.075,.936,'Actual native t = 0.020 s | 200 steps | 9,072 cells | dt = 0.1 ms',fontsize=12)
    fig.text(.075,.895,'Both short-run numerical audits passed. This 0.02 s qualification is separate from the planned 0.5 s case.',fontsize=10)
    fig.text(.075,.073,'Grey: solid obstacle. Velocity shown only in cells with alpha >=0.5; arrows below 0.5% of the common scale are hidden.',fontsize=9)
    note='Quiet-control speed scale is fixed at its 1e-4 m/s acceptance limit.' if args.mode=='flat' else 'Speed scale is shared between these two actual snapshots; arrows are instantaneous velocity, not particle paths.'
    fig.text(.075,.046,note,fontsize=9);fig.text(.075,.022,'Source-bound native fields and clocks checked. These qualification results do not establish physical validation.',fontsize=9)
    args.out.mkdir(parents=True);fig.savefig(args.out/'paired-qualification.png',dpi=140);fig.savefig(args.out/'paired-qualification.pdf');plt.close(fig)
    (args.out/'plot-receipt.json').write_text(json.dumps(dict(mode=args.mode,actual_native_time_s=.02,actual_step=200,audits_sha256={str(p):sha(p) for p in apaths},native_source_sha256={str(p):sha(p) for p in sources},common_speed_scale_max_m_s=vmax,full_production_complete=False,physical_validation=False,code_sha256=sha(Path(__file__))),indent=2)+'\n');print(args.out/'paired-qualification.png')
if __name__=='__main__':main()