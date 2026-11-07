"""Four actual wing fields, shown in each case's wind-aligned frame.
The animation integrates tracers through saved steady fields, not transient CFD.
"""
import gzip,hashlib,json,sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
from PIL import Image
from plot_finite_wing_native_v002 import TetGeometry,read_fluent,read_foam,display,decorate_3d,sha
from wing_zero_common_v001 import U,LIFT
R=Path('E:/CFD-Learning-Lab')
ORIGIN=np.array([.25,0,0.]);B0=np.column_stack([U,LIFT,[0,0,1]])
def wind(points,basis):return ORIGIN+(np.asarray(points)-ORIGIN)@basis

def main(out,zero_fluent):
    assert not out.exists()
    geo_path=R/'native/runs/2026-11-01-finite-wing-v003/mesh-export-v002/geometry.npz'
    geo=TetGeometry(geo_path)
    specs=[('Fluent 0 deg','ANSYS Fluent  |  0°',zero_fluent,0,'fluent'),
           ('Fluent 8 deg','ANSYS Fluent  |  8°',R/'native/runs/finite-wing-coupled-fluent-v001',8,'fluent'),
           ('OpenFOAM 0 deg','OpenFOAM  |  0°',R/'native/runs/finite-wing-zero-openfoam-v002',0,'openfoam'),
           ('OpenFOAM 8 deg','OpenFOAM  |  8°',R/'native/runs/2026-11-01-finite-wing-v003/openfoam',8,'openfoam')]
    entries=[];hashes={str(geo_path):sha(geo_path)}
    seeds=np.array([[-.25,y,z] for y in [-.15,.15] for z in [-2.2,-2.06,-1.94,-1.8,0,1.8,1.94,2.06,2.2]])
    for key,label,path,angle,solver in specs:
        nat=read_fluent(path/'wing-final.cas.h5',path/'wing-final.dat.h5') if solver=='fluent' else read_foam(path,'1200','wing')
        if solver=='fluent' and angle==0:
            assert max(nat['solve_evidence']['native_terminal_scaled_residuals'].values())>1e-5
            label += ' (residual miss)'
        field=geo.match(nat);geo.reconstruct(field);basis=B0 if angle==0 else np.eye(3)
        physical=ORIGIN+(seeds-ORIGIN)@basis.T
        traces=[geo.trace(p,field['nodal_U'],5.) for p in physical]
        for p in field['inputs']:hashes[str(p)]=sha(p)
        entries.append(dict(key=key,label=label,angle=angle,basis=basis,field=field,traces=traces,
            wing=display(wind(geo.wing_vertices,basis)),solve_evidence=nat['solve_evidence']))
    out.mkdir(parents=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.7})
    maxspeed=max(float(np.linalg.norm(e['field']['U'],axis=1).max()) for e in entries)
    norm=Normalize(0,maxspeed);cmap=plt.get_cmap('viridis')
    fig=plt.figure(figsize=(12.8,10.0),dpi=100,facecolor='white')
    for i,e in enumerate(entries):
        ax=fig.add_axes([.025+(i%2)*.50,.535-(i//2)*.425,.435,.38],projection='3d')
        ax.add_collection3d(Poly3DCollection(e['wing'],facecolors='#b8c3cc',edgecolor='none',linewidth=0))
        for tr in e['traces']:
            if len(tr['points'])<2:continue
            p=display(wind(tr['points'],e['basis']))
            line=Line3DCollection(np.stack([p[:-1],p[1:]],axis=1),cmap=cmap,norm=norm,linewidths=.8,alpha=.65)
            line.set_array((tr['speeds'][:-1]+tr['speeds'][1:])/2);ax.add_collection3d(line)
        e['dot']=ax.scatter([],[],[],s=13,c=[],cmap=cmap,norm=norm,depthshade=False,edgecolor='none')
        decorate_3d(ax,e['label']);ax.set(xlabel='$x_{flow}$ (m)',ylabel='$z$ (m)',zlabel='$y_{lift}$ (m)');e['ax']=ax
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.24,.070,.52,.018]),orientation='horizontal')
    cb.set_label('Saved speed (m/s)  |  shared scale',size=10)
    fig.text(.5,.976,'A finite wing at 0° and 8° incidence',ha='center',size=18,weight='bold')
    fig.text(.5,.945,'Re = 400  |  same mesh  |  axes aligned with the freestream in each panel',ha='center',size=10,color='#4b6578')
    fig.text(.5,.014,'Tracers through actual steady solutions; camera motion is visualization, not transient CFD.',ha='center',size=9,color='#4b6578')
    images=[]
    for frame in range(48):
        phase=frame/48
        for e in entries:
            good=[t for t in e['traces'] if len(t['points'])>=2 and t['times'][-1]>0]
            positions=[];speeds=[]
            for j,tr in enumerate(good):
                travel=((phase+j/max(1,len(good)))%1)*tr['times'][-1]
                positions.append([np.interp(travel,tr['times'],tr['points'][:,k]) for k in range(3)])
                speeds.append(float(np.interp(travel,tr['times'],tr['speeds'])))
            if positions:
                p=display(wind(np.array(positions),e['basis']));e['dot']._offsets3d=(p[:,0],p[:,1],p[:,2]);e['dot'].set_array(np.array(speeds))
            e['ax'].view_init(elev=24+2*np.sin(2*np.pi*phase),azim=-55+12*np.sin(2*np.pi*phase))
        fig.canvas.draw();images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()).quantize(colors=256))
        if frame==0:fig.savefig(out/'wing-incidence-tracers.png',dpi=145)
    images[0].save(out/'wing-incidence-paired.gif',save_all=True,append_images=images[1:],duration=70,loop=0,disposal=2,optimize=False);plt.close(fig)
    maximum=max(float(abs(2*e['field']['wing_p']).max()) for e in entries)
    pn=TwoSlopeNorm(vmin=-maximum,vcenter=0,vmax=maximum);pcmap=plt.get_cmap('RdBu_r')
    fig=plt.figure(figsize=(12.8,9),facecolor='white')
    for i,e in enumerate(entries):
        ax=fig.add_axes([.03+(i%2)*.50,.53-(i//2)*.42,.44,.37],projection='3d')
        ax.add_collection3d(Poly3DCollection(e['wing'],facecolors=pcmap(pn(2*e['field']['wing_p'])),edgecolor='none',linewidth=0))
        decorate_3d(ax,e['label']);ax.set(xlim=(-.1,1.25),ylim=(-2.3,2.3),zlim=(-.35,.35),xlabel='$x_{flow}$ (m)',ylabel='$z$ (m)',zlabel='$y_{lift}$ (m)');ax.set_box_aspect((1.35,4.6,.7))
        method='Owner-cell pressure proxy' if 'proxy' in e['field']['pressure_method'] else 'Native face pressure'
        ax.text2D(.5,0,method,transform=ax.transAxes,ha='center',fontsize=8)
    cb=fig.colorbar(matplotlib.cm.ScalarMappable(norm=pn,cmap=pcmap),cax=fig.add_axes([.23,.055,.54,.021]),orientation='horizontal');cb.set_label('Pressure coefficient $C_p$  |  shared scale')
    fig.suptitle('Wing surface pressure at 0° and 8° incidence',fontsize=17,y=.97)
    fig.text(.5,.925,'Flow-aligned display coordinates; the computational mesh is unchanged.',ha='center',fontsize=10,color='#4b6578')
    fig.savefig(out/'wing-incidence-pressure.png',dpi=150);plt.close(fig)
    with Image.open(out/'wing-incidence-paired.gif') as im:
        count=im.n_frames;unique=set();delays=[]
        for i in range(count):im.seek(i);unique.add(hashlib.sha256(im.convert('RGB').tobytes()).hexdigest());delays.append(im.info['duration'])
        assert count==48 and len(unique)==48 and delays==[70]*48
        im.seek(24);im.convert('RGB').save(out/'gif-midpoint.png')
    info=dict(cases=[dict(key=e['key'],angle=e['angle'],basis_columns=e['basis'].tolist(),pressure_method=e['field']['pressure_method'],native_solve_evidence=e['solve_evidence']) for e in entries],
        native_geometry_matches=True,physical_validation=False,zero_fluent_residual_check=False,transient_simulation=False,frames=48,unique_frames=48,
        frame_delay_ms=70,speed_limits=[0,maxspeed],pressure_Cp_limits=[-maximum,maximum],
        transform_origin=ORIGIN.tolist(),display_formula='p_wind = origin + (p_global-origin) @ [drag, lift, span]; display axes permuted (x_flow,z,y_lift)',
        tracer_method='Same 18 seeds in wind coordinates; RK4 through native volume-weighted nodal velocities in original fluid tetrahedra; boundary crossings stop the path')
    (out/'visual-checks.json').write_text(json.dumps(info,indent=2)+'\n')
    with gzip.open(out/'actual-traces.json.gz','wt',encoding='utf-8') as f:
        json.dump({e['key']:[{k:(v.tolist() if isinstance(v,np.ndarray) else v) for k,v in tr.items()} for tr in e['traces']] for e in entries},f)
    assert all(sha(Path(p))==h for p,h in hashes.items())
    (out/'native-source-hashes-private.json').write_text(json.dumps(hashes,indent=2)+'\n')
    print(json.dumps(info,indent=2))
if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
