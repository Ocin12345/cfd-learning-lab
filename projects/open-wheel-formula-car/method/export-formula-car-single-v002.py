"""Render and export the audited 8-degree symmetric formula-car model."""
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


def export(case, reader, output):
    audit = runpy.run_path(str(WORK / 'audit-formula-car-symmetric-v002.py'))
    q, geometry, centres, velocity, _, force, residual = audit['audit'](case, reader, 600)
    assert q['pass'] and all(q['gates'].values())
    output.mkdir()
    (output / 'checks.json').write_text(json.dumps({'pass': True, 'flap_angle_deg': 8, 'case': q}, indent=2) + '\n',
                                        encoding='utf-8')
    with (output / 'forces.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('iteration','full_equivalent_drag_n','full_equivalent_downforce_n',
                         'drag_coefficient','downforce_coefficient'))
        ref = q['reference_force_n']
        writer.writerows((i,x,-z,x/ref,-z/ref) for i,(x,_,z) in enumerate(force))
    with (output / 'initial-residuals.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('iteration',*audit['EQS']))
        writer.writerows((i,*row) for i,row in enumerate(residual,1))
    with (output / 'component-loads.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(('component','mean_drag_n','mean_downforce_n'))
        writer.writerows((name,row['mean_drag_n'],row['mean_downforce_n'])
                         for name,row in q['component_forces'].items())

    fig = plt.figure(figsize=(12,8.2),layout='constrained')
    fig.suptitle('Ground effect on an open-wheel formula car',fontsize=19,fontweight='bold')
    grid = fig.add_gridspec(3,2,height_ratios=(1.65,1,1))
    panels = []
    cp_norm = TwoSlopeNorm(vmin=-1.8,vcenter=0,vmax=1.2)
    mirrored = []
    for polygons,_,_,_,cp in geometry.values():
        reflected = [poly*np.array([1.,-1.,1.]) for poly in polygons]
        mirrored.append((polygons+reflected,np.concatenate([cp,cp])))
    for column,(view,azimuth) in enumerate((('Front three-quarter',143),('Rear three-quarter',34))):
        body = fig.add_subplot(grid[0,column],projection='3d')
        panels.append(body)
        for polygons,cp in mirrored:
            surface = Poly3DCollection(polygons,cmap='coolwarm',norm=cp_norm,
                                       edgecolors='none',linewidth=0,antialiased=False)
            surface.set_array(cp)
            body.add_collection3d(surface)
        body.set(xlim=(-.08,4.65),ylim=(-1.12,1.12),zlim=(0,1.24),title=view)
        body.set_box_aspect((4.73,2.24,1.24),zoom=1.32)
        body.view_init(elev=25,azim=azimuth)
        body.set_axis_off()
        body.set_anchor('C')
    fig.colorbar(surface,ax=panels,label='Surface pressure coefficient (mirrored)',
                 shrink=.8,extend='both',panchor=False)
    colours = ('#136f86','#a74d30')
    speed = np.linalg.norm(velocity,axis=1)/30
    flow = fig.add_subplot(grid[1,:])
    x,y,z = centres.T
    mask = (abs(y)<.052)&(x>-.5)&(x<4.95)&(z<1.45)
    dots = flow.scatter(x[mask],z[mask],c=speed[mask],s=7,cmap='viridis',vmin=0,vmax=1.5,
                        edgecolors='none',rasterized=True)
    flow.axvline(3.2,color='.4',ls=':',lw=.7)
    flow.set(xlim=(-.5,4.95),ylim=(0,1.45),xlabel='Streamwise position (m)',
             ylabel='Height (m)',title='Near-centreplane cell samples · 8° flap · 30 m/s · zero yaw')
    fig.colorbar(dots,ax=flow,label='Speed / inlet speed',extend='max')
    loads = fig.add_subplot(grid[2,0])
    iterations = np.arange(len(force))
    for label,values,colour in zip(('Drag','Downforce'),(force[:,0],-force[:,2]),colours):
        loads.plot(iterations[20:],values[20:],label=label,color=colour,lw=1.4)
    loads.set(xlabel='Steady iteration',ylabel='Full-car equivalent (N)',title='Vehicle force history',xlim=(20,600))
    loads.legend(frameon=False,fontsize=8)
    loads.grid(alpha=.17)
    errors = fig.add_subplot(grid[2,1])
    for label,values in zip(audit['EQS'],residual.T):
        errors.semilogy(iterations[1:],values,label=label,lw=1)
    errors.axhline(1e-3,color='.35',ls=':',lw=1)
    errors.set(xlabel='Steady iteration',ylabel='Initial residual',title='Pressure, momentum and turbulence')
    errors.legend(frameon=False,fontsize=8,ncol=3)
    errors.grid(alpha=.17)
    fig.savefig(output/'formula-car-flow.png',dpi=180)
    fig.savefig(output/'formula-car-flow.pdf')
    plt.close(fig)

    fig, axes = plt.subplots(2,1,figsize=(9,7),layout='constrained')
    fig.suptitle('Front-wing flow and vehicle load distribution',fontsize=17,fontweight='bold')
    wing, parts = axes
    region = (abs(y)<.023)&(x>-.13)&(x<1.04)&(z<.57)
    dots = wing.scatter(x[region],z[region],c=speed[region],s=9,cmap='viridis',vmin=0,vmax=1.6,
                        edgecolors='none',rasterized=True)
    for name in ('frontMain','frontFlap'):
        _, xyz, *_ = geometry[name]
        mid = abs(xyz[:,1])<.024
        wing.scatter(xyz[mid,0],xyz[mid,2],s=6,c='#273346')
    for label,target,position in (
        ('Main element',(.255,.254),(.17,.425)), ('Flap',(.65,.331),(.57,.48))
    ):
        wing.annotate(label,xy=target,xytext=position,ha='center',fontsize=9,
                      bbox={'facecolor':'white','edgecolor':'none','pad':2},
                      arrowprops={'arrowstyle':'-','color':'#273346','lw':.8})
    wing.set(xlabel='Streamwise position (m)',ylabel='Height above road (m)',
             xlim=(-.13,1.04),ylim=(0,.57),title='8° front flap · near-centreplane cell samples')
    fig.colorbar(dots,ax=wing,label='Speed / inlet speed',extend='max')
    groups = ('frontWingForce','chassisForce','rearWingForce','wheelsCombined')
    positions = np.arange(len(groups))
    for i,(key,label,colour) in enumerate(zip(('mean_drag_n','mean_downforce_n'),('Drag','Downforce'),colours)):
        bars = parts.bar(positions+(i-.5)*.36,[q['component_forces'][name][key] for name in groups],
                         .36,label=label,color=colour)
        parts.bar_label(bars,fmt='%.0f',padding=3,fontsize=8)
    parts.axhline(0,color='.4',lw=.7)
    parts.set(xticks=positions,xticklabels=('Front wing','Chassis','Rear wing','Wheels'),
              ylabel='Full-car equivalent (N)',title='Mean over the final 50 iterations')
    parts.margins(y=.18)
    parts.grid(axis='y',alpha=.17)
    parts.legend(frameon=False,fontsize=9)
    fig.savefig(output/'front-wing-detail.png',dpi=180)
    fig.savefig(output/'front-wing-detail.pdf')
    plt.close(fig)
    print({key:q[key] for key in ('cells','pass','last_50_mean_drag_n','last_50_mean_downforce_n','component_forces')})


if __name__=='__main__':
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('case',type=Path)
    cli.add_argument('reader',type=Path)
    cli.add_argument('output',type=Path)
    args = cli.parse_args()
    export(args.case,args.reader,args.output)
