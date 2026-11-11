"""Audit the steady endpoints and render a matched formula-car flap comparison."""
import argparse
import csv
import json
import runpy
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

WORK = Path(__file__).resolve().parent
native_audit = runpy.run_path(str(WORK/'audit-flap-recovery-v002.py'))
ANGLES = (8,12,18)
COLORS = ('#16768d','#b65732','#6a508d')


def write_csv(path, headings, rows):
    with path.open('x',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream)
        writer.writerow(headings)
        writer.writerows(rows)


def plot_fields(output, results):
    fig=plt.figure(figsize=(14,8.5),layout='constrained')
    fig.suptitle('Front-flap angle on an open-wheel formula car',fontsize=19,fontweight='bold')
    grid=fig.add_gridspec(3,len(results),height_ratios=(1.3,1,1))
    pressure_axes, full_axes, wing_axes = [],[],[]
    normal=TwoSlopeNorm(vmin=-1.8,vcenter=0,vmax=1.2)
    for column,(angle,(_,geometry,centres,velocity,*_)) in enumerate(results):
        ax=fig.add_subplot(grid[0,column],projection='3d')
        pressure_axes.append(ax)
        for polygons,_,_,_,cp in geometry.values():
            reflected=[poly*np.array([1.,-1.,1.]) for poly in polygons]
            surface=Poly3DCollection(polygons+reflected,cmap='coolwarm',norm=normal,
                                     edgecolors='none',linewidth=0,antialiased=False)
            surface.set_array(np.concatenate((cp,cp)))
            ax.add_collection3d(surface)
        ax.set(xlim=(-.08,4.65),ylim=(-1.12,1.12),zlim=(0,1.24),title=f'{angle}° front flap')
        ax.set_box_aspect((4.73,2.24,1.24),zoom=1.28)
        ax.view_init(elev=25,azim=143)
        ax.set_axis_off()
        ax.set_anchor('C')
        x,y,z=centres.T
        speed=np.linalg.norm(velocity,axis=1)/30
        ax=fig.add_subplot(grid[1,column])
        full_axes.append(ax)
        view=(abs(y)<.052)&(x>-.5)&(x<4.95)&(z<1.45)
        flow=ax.scatter(x[view],z[view],c=speed[view],s=6,cmap='viridis',vmin=0,vmax=1.6,
                        edgecolors='none',rasterized=True)
        ax.axvline(3.2,color='.4',ls=':',lw=.7)
        ax.set(xlim=(-.5,4.95),ylim=(0,1.45),title='Near-centreplane flow',xlabel='Streamwise position (m)')
        ax=fig.add_subplot(grid[2,column])
        wing_axes.append(ax)
        view=(abs(y)<.023)&(x>-.13)&(x<1.04)&(z<.57)
        ax.scatter(x[view],z[view],c=speed[view],s=7,cmap='viridis',vmin=0,vmax=1.6,
                   edgecolors='none',rasterized=True)
        for part in ('frontMain','frontFlap'):
            _,xyz,*_=geometry[part]
            mask=abs(xyz[:,1])<.024
            ax.scatter(xyz[mask,0],xyz[mask,2],s=5,c='#273346')
        ax.set(xlim=(-.13,1.04),ylim=(0,.57),title='Front-wing detail',xlabel='Streamwise position (m)')
    for axes in (full_axes,wing_axes):
        axes[0].set_ylabel('Height above road (m)')
    fig.colorbar(surface,ax=pressure_axes,label='Surface pressure coefficient (mirrored)',shrink=.85,extend='both',panchor=False)
    fig.colorbar(flow,ax=full_axes+wing_axes,label='Speed / inlet speed',shrink=.82,extend='max')
    for ext in ('png','pdf'):
        fig.savefig(output/f'flap-angle-flow.{ext}',dpi=180)
    plt.close(fig)


def plot_forces(output, results):
    fig, axes=plt.subplots(2,2,figsize=(11,7.8),layout='constrained')
    fig.suptitle('Drag and downforce by front-flap angle',fontsize=18,fontweight='bold')
    groups=('Total','Front wing','Chassis','Rear wing','Wheels')
    names=('frontWingForce','chassisForce','rearWingForce','wheelsCombined')
    positions=np.arange(len(groups))
    for col,(key,label,component,coordinate,sign) in enumerate((
        ('last_50_mean_drag_n','Drag','mean_drag_n',0,1),
        ('last_50_mean_downforce_n','Downforce','mean_downforce_n',2,-1),
    )):
        bars,history=axes[:,col]
        width=.78/len(results)
        for idx,(angle,(case,_,_,_,_,force,*_)) in enumerate(results):
            loads=[case[key],*(case['component_forces'][name][component] for name in names)]
            bars.bar(positions+(idx-(len(results)-1)/2)*width,loads,width,
                     label=f'{angle}°',color=COLORS[idx])
            history.plot(np.arange(-49,1),sign*force[-50:,coordinate],label=f'{angle}°',color=COLORS[idx],lw=1.5)
        bars.axhline(0,color='.4',lw=.7)
        bars.set(xticks=positions,xticklabels=groups,ylabel=f'{label} (N)',title=f'{label} by component')
        history.set(xlabel='Steady iterations relative to endpoint',ylabel=f'{label} (N)',title=f'Final 50 {label.lower()} samples')
        for ax in (bars,history):
            ax.grid(axis='y',alpha=.19)
            ax.legend(frameon=False,ncol=len(results),fontsize=8)
    for ext in ('png','pdf'):
        fig.savefig(output/f'flap-angle-loads.{ext}',dpi=180)
    plt.close(fig)


def finish(results, output):
    output.mkdir()
    summary={'pass':all(data[0]['pass'] for _,data in results),'cases':{str(angle):data[0] for angle,data in results}}
    baseline=results[0][1][0]
    summary['method']=baseline['restart_stages'][-1]['phase'] if baseline['restart_stages'] else 'original-steady'
    summary['differing_pre_mesh_input_files']=['constant/geometry/frontFlap.stl']
    summary['changes_vs_8_degrees']={str(angle):{
        'drag_n': data[0]['last_50_mean_drag_n']-baseline['last_50_mean_drag_n'],
        'downforce_n': data[0]['last_50_mean_downforce_n']-baseline['last_50_mean_downforce_n'],
    } for angle,data in results[1:]}
    (output/'checks.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    if not summary['pass']:
        return False
    def phase(iteration):
        if iteration<=600:
            return 'original'
        if summary['method']=='blended':
            return 'blendedSteady'
        return 'localEuler' if iteration<=1500 else 'steadyVerification'

    write_csv(output/'forces.csv',('angle_deg','iteration','phase','drag_n','downforce_n'),(
        (angle,i,phase(i),x,-z)
        for angle,(_,_,_,_,_,force,_,_) in results for i,(x,_,z) in enumerate(force)))
    write_csv(output/'initial-residuals.csv',('angle_deg','iteration','phase',*native_audit['EQS'],'p_second_correction'),(
        (angle,i,phase(i),*row,p[i-601] if phase(i)=='localEuler' else '')
        for angle,(_,_,_,_,_,_,res,p) in results for i,row in enumerate(res,1)))
    write_csv(output/'component-loads.csv',('angle_deg','component','mean_drag_n','mean_downforce_n'),(
        (angle,name,q['mean_drag_n'],q['mean_downforce_n'])
        for angle,(case,*_) in results for name,q in case['component_forces'].items()))
    plot_fields(output,results)
    plot_forces(output,results)
    return True


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('native','reader','output'):
        parser.add_argument(name,type=Path)
    parser.add_argument('--blended',action='store_true')
    args=parser.parse_args()
    originals=[args.native/f'flap{angle:02}-symmetric-v001' for angle in ANGLES]
    manifests=[dict(line.split('  ',1)[::-1] for line in (folder/'pre-mesh-sha256.txt').read_text().splitlines())
               for folder in originals]
    baseline=manifests[0]
    for current in manifests[1:]:
        assert current.keys()==baseline.keys()
        assert {key for key in baseline if baseline[key]!=current[key]}=={'constant/geometry/frontFlap.stl'}
    results=[]
    for angle in ANGLES:
        prefix=f'flap{angle:02}'
        original=args.native/f'{prefix}-symmetric-v001'
        stages=() if angle==8 else (args.native/f'{prefix}-localEuler-v001',args.native/f'{prefix}-steady-verification-v001')
        cap=1800 if stages else 600
        if args.blended:
            stages=(args.native/f'{prefix}-blended-v001',)
            cap=1200
        data=native_audit['audit'](original,args.reader,cap,stages)
        q=data[0]
        print(angle,q['gates'],q['last_50_mean_drag_n'],q['last_50_mean_downforce_n'],flush=True)
        results.append((angle,data))
    raise SystemExit(0 if finish(results,args.output) else 2)


if __name__=='__main__':
    main()
