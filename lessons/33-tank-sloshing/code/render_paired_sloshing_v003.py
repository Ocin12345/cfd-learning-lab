"""Render synchronized actual native sloshing states after both production audits pass."""
from pathlib import Path
import argparse,hashlib,json,shutil
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from analyse_nozzle_fluent_v005 import common_mapping
from audit_sloshing_openfoam_v001 import field
R=Path('E:/CFD-Learning-Lab');P=R/'local/sloshing-production-v001';N=R/'native/runs/sloshing-production-v001'
F=N/'fluent'
STEPS=list(range(0,2401,20));TIMES=np.array(STEPS)*.0005
FPS=12

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def chunk(group):
    a=np.full(3200,np.nan)
    for ds in group.values():
        lo=int(np.asarray(ds.attrs['minId']).flat[0])-1;hi=int(np.asarray(ds.attrs['maxId']).flat[0])
        assert np.isnan(a[lo:hi]).all();a[lo:hi]=np.asarray(ds).ravel()
    assert np.isfinite(a).all()
    return a

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--openfoam-case',type=Path,required=True)
    parser.add_argument('--openfoam-audit',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--check-only',action='store_true')
    args=parser.parse_args()
    O=args.openfoam_case.resolve();OUT=args.out.resolve()
    assert O.is_relative_to(N.resolve()) and OUT.is_relative_to((R/'exports/figures').resolve())
    assert args.openfoam_audit.resolve().is_relative_to(P.resolve())
    assert not OUT.exists()
    audit_paths=[P/'fluent-production-audit-v001.json',args.openfoam_audit.resolve()]
    audits=[json.loads(p.read_text()) for p in audit_paths]
    for a in audits:
        assert a['production_checks_pass'] and a['native_field_checks_pass'] and a['production_complete']
        for name,expected in a['native_source_sha256'].items():assert sha(Path(name))==expected,name
        assert len(a['records'])==121
        assert np.allclose([r['time_s'] for r in a['records']],TIMES,atol=1e-10,rtol=0)
    proofs=[{str(Path(k).resolve()):v for k,v in a['native_source_sha256'].items()} for a in audits]
    def proven(path,solver):
        assert str(path.resolve()) in proofs[solver], 'Selected field has no matching solver-audit provenance: '+str(path)
        assert sha(path)==proofs[solver][str(path.resolve())]
    proven(F/'tank.cas.h5',0);proven(O/'0/C',1);proven(O/'solver.log',1)
    if args.check_only:
        print(json.dumps(dict(selected_case_audit_binding_pass=True,full_numerical_audits_pass=True,rendered_frames=0)))
        return
    with np.load(F/'geometry.npz') as g:
        fm=common_mapping(F/'tank.cas.h5',g);xy=g['centres']
    oc=field(O/'0/C',3);ij=np.rint(oc[:,:2]/.0025-.5).astype(int);om=ij[:,0]+80*ij[:,1]
    assert len(np.unique(om))==3200 and np.all((om>=0)&(om<3200))
    assert np.max(abs(oc[:,:2]-(ij+.5)*.0025))<1e-10
    arrays=np.empty((2,121,3,40,80));sources=[]
    for index,step in enumerate(STEPS):
        p=F/f'tank-{step:04d}.dat.h5'
        proven(p,0)
        with h5py.File(p) as f:
            for k,name in enumerate(['results/1/phase-3/cells/SV_VOF','results/1/phase-1/cells/SV_U','results/1/phase-1/cells/SV_V']):
                a=np.empty(3200);a[fm]=chunk(f[name]);arrays[0,index,k]=a.reshape(40,80)
        sources.append(p)
        folder=O/f'{step*.0005:.8g}'
        proven(folder/'alpha.water',1);proven(folder/'U',1)
        alpha=field(folder/'alpha.water',1).ravel();velocity=field(folder/'U',3)
        for k,raw in enumerate([alpha,velocity[:,0],velocity[:,1]]):
            a=np.empty(3200);a[om]=raw;arrays[1,index,k]=a.reshape(40,80)
        sources.extend([folder/'alpha.water',folder/'U'])
    assert np.isfinite(arrays).all()
    alpha=arrays[:,:,0];u=arrays[:,:,1];v=arrays[:,:,2]
    assert alpha.min()>=-1e-6 and alpha.max()<=1+1e-6
    assert np.max(abs(alpha[0,0]-alpha[1,0]))<1e-10
    speed=np.hypot(u,v);depth=alpha.sum(axis=2)*2.5-50
    volumes=alpha.sum(axis=(2,3))*.0025**2
    assert np.max(abs(volumes/.01-1))<=.005
    maximum_speed=float(np.max(np.where(alpha>=.5,speed,0)))
    assert maximum_speed>0
    bound=max(3.,float(np.max(abs(depth)))*1.1)
    x=xy[:80,0]*1000;y=xy[::80,1]*1000;xx,yy=np.meshgrid(x,y)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(12,12))
    grid=fig.add_gridspec(3,2,left=.08,right=.89,bottom=.10,top=.86,hspace=.49,wspace=.30)
    axes=[[fig.add_subplot(grid[i,j]) for j in range(2)] for i in range(3)]
    water=LinearSegmentedColormap.from_list('water',['#fafafa','#087eb5'])
    velocity_cmap=plt.get_cmap('viridis').with_extremes(bad='#fafafa')
    images=[[],[]];arrows=[]
    labels=['ANSYS Fluent','OpenFOAM 14']
    for solver in range(2):
        for row in [0,1]:
            ax=axes[row][solver]
            values=alpha[solver,0] if row==0 else np.ma.masked_where(alpha[solver,0]<.5,speed[solver,0])
            image=ax.imshow(values,origin='lower',extent=[0,200,0,100],vmin=0,vmax=1 if row==0 else maximum_speed,cmap=water if row==0 else velocity_cmap,interpolation='nearest',aspect='equal')
            images[row].append(image)
            ax.set(xlabel='x (mm)',ylabel='y (mm)',xticks=[0,50,100,150,200],yticks=[0,25,50,75,100],title=labels[solver]+(' | water fraction' if row==0 else ' | mixture speed'))
        selection=(slice(2,None,6),slice(3,None,8))
        mask=(alpha[solver,0][selection]<.5)|(speed[solver,0][selection]<maximum_speed*.005)
        q=axes[1][solver].quiver(xx[selection],yy[selection],np.ma.masked_where(mask,u[solver,0][selection]),np.ma.masked_where(mask,v[solver,0][selection]),color='white',angles='xy',scale_units='xy',scale=maximum_speed/18,width=.004,minlength=0)
        arrows.append(q)
    for row,label in [(0,'Water volume fraction'),(1,'Mixture speed (m/s)')]:
        pos=axes[row][1].get_position();cax=fig.add_axes([.92,pos.y0,.015,pos.height]);fig.colorbar(images[row][0],cax=cax).set_label(label)
    colors=['#087eb5','#cc6423'];profiles=[];histories=[]
    for solver in range(2):
        profiles.append(axes[2][0].plot(x,depth[solver,0],color=colors[solver],lw=1.7,label=labels[solver])[0])
        histories.append(axes[2][1].plot(TIMES[:1],depth[solver,:1,-1],color=colors[solver],lw=1.7,label=labels[solver])[0])
    reference=axes[2][0].plot(x,depth[0,0],color='#777777',ls='--',lw=.8,label='Initial profile')[0]
    axes[2][0].set(xlabel='x (mm)',ylabel='Column depth minus 50 mm (mm)',title='Column-integrated water depth',xlim=[0,200],ylim=[-bound,bound])
    axes[2][1].set(xlabel='Physical time (s)',ylabel='Column depth minus 50 mm (mm)',title='Right-most column: x = 198.75 mm',xlim=[0,1.2],ylim=[-bound,bound])
    for ax in axes[2]:ax.grid(alpha=.2)
    fig.legend([profiles[0],profiles[1],reference],labels+['Initial spatial profile (left)'],bbox_to_anchor=(.485,.329),loc='center',ncol=3,frameon=False,fontsize=9)
    cursor=axes[2][1].axvline(0,color='#555555',lw=.8,ls=':')
    fig.suptitle('Small-amplitude water-tank sloshing',x=.08,ha='left',y=.98,fontsize=21,fontweight='bold')
    clock=fig.text(.08,.93,'',fontsize=14)
    fig.text(.08,.895,'200 x 100 mm tank | 3,200 cells | dt = 0.5 ms | same initial fractions and material properties',fontsize=10)
    fig.text(.08,.047,'Mixture speed uses water fraction >=0.5; velocity arrows below 0.5% of the common speed scale are hidden.',fontsize=9)
    fig.text(.08,.022,'121 native states, 0.01 s apart, shown at 12 fps. Both numerical audits passed; mesh/time independence and physical validation remain open.',fontsize=9)
    OUT.mkdir(parents=True)
    selected={0,30,60,90,120}
    for index,t in enumerate(TIMES):
        clock.set_text(f'Physical time t = {t:.2f} s | native step {STEPS[index]:04d} / 2400')
        for solver in range(2):
            images[0][solver].set_data(alpha[solver,index])
            images[1][solver].set_data(np.ma.masked_where(alpha[solver,index]<.5,speed[solver,index]))
            mask=(alpha[solver,index][selection]<.5)|(speed[solver,index][selection]<maximum_speed*.005)
            arrows[solver].set_UVC(np.ma.masked_where(mask,u[solver,index][selection]),np.ma.masked_where(mask,v[solver,index][selection]))
            profiles[solver].set_ydata(depth[solver,index])
            histories[solver].set_data(TIMES[:index+1],depth[solver,:index+1,-1])
        cursor.set_xdata([t,t])
        fig.savefig(OUT/f'frame-{index:03d}.png',dpi=120)
        if index in selected:fig.savefig(OUT/f'state-{index:03d}.pdf')
        print(json.dumps(dict(rendered_frame=index,physical_time_s=float(t))),flush=True)
    plt.close(fig)
    max_hdiff=float(np.max(abs(depth[0]-depth[1])))
    record=dict(all_saved_field_checks_pass=True,rendering=dict(frames=121,frame_duration_ms=1000/FPS),plot=dict(times_s=TIMES.tolist(),global_water_dominated_speed_scale_max_m_s=maximum_speed,column_depth_range_mm=[-bound,bound]),maximum_column_depth_difference_mm=max_hdiff,final_column_depth_rms_difference_mm=float(np.sqrt(np.mean((depth[0,-1]-depth[1,-1])**2))),maximum_relative_water_area_drift=float(np.max(abs(volumes/.01-1))),numerical_audit_sha256={str(p):sha(p) for p in audit_paths},physical_validation=False,mesh_independence=False,time_step_independence=False,interpolated_motion=False,native_source_sha256={str(p):sha(p) for p in sources},source_sha256={str(Path(__file__)):sha(Path(__file__))})
    (OUT/'native-field-audit.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n')
    (OUT/'column-depth-histories.json').write_text(json.dumps(dict(times_s=TIMES.tolist(),x_mm=x.tolist(),depth_minus_50mm=depth.tolist(),solver_labels=labels),allow_nan=False)+'\n')
    shutil.copyfile(__file__,OUT/'executed-renderer.py')
    print(json.dumps(dict(output=str(OUT),frames=121,maximum_column_depth_difference_mm=max_hdiff)))
if __name__=='__main__':main()
