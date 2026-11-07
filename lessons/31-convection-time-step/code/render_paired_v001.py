"""Read actual native cell fields, verify geometric correspondence, render qualification only."""
from pathlib import Path
import json,re,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from cavity_production_common_v004 import *
DT=.01; HEATED_STEPS=2000
import plot_nozzle_native_v004 as native
F=ROOT/'native/runs/natural-convection-dt001-v001/fluent-v002'
O=ROOT/'native/runs/natural-convection-dt001-v001/openfoam'
OUT=ROOT/'exports/figures/2026-11-07-natural-convection-dt001-paired-v001'

def reorder(q):
    xy=q['nodes'][q['cells'],:2].mean(axis=1)
    ij=np.rint(xy/(L/N)-.5).astype(int)
    assert np.all((ij>=0)&(ij<N))
    assert np.max(abs(xy-(ij+.5)*(L/N)))<1e-9
    slots=ij[:,1]*N+ij[:,0]
    assert np.array_equal(np.sort(slots),np.arange(N*N))
    # Every native cell must span exactly one uniform source quad in x and y.
    verts=q['nodes'][q['cells'],:2]
    assert np.max(abs(np.ptp(verts,axis=1)-L/N))<1e-9
    idx=np.argsort(slots)
    for key in ['U','T','p','rho']:q[key]=q[key][idx]
    q['mapping']=f'native vertex geometry mapped bijectively to {N}x{N} uniform cells'
    return q

def ofields(t):
    q=native.foam_mesh(O);count=N*N;p=O/f'{t:g}'
    for key,w,dim in [('U',3,[0,1,-1,0,0,0,0]),('T',1,[0,0,0,1,0,0,0]),('p',1,[1,-1,-2,0,0,0,0])]:q[key]=native.foam_field(p/key,count,w,dim)
    q['U']=q['U'][:,:2]
    if t==0:
        q['rho']=q['p']/(R*q['T']);q['density_source']='EOS from saved start p,T (no start rho file)'
    else:q['rho']=native.foam_field(p/'rho',count,1,[1,-3,0,0,0,0,0]);q['density_source']='native rho'
    assert all(np.isfinite(q[k]).all() for k in ['U','T','p','rho'])
    return reorder(q)

def stats(q):
    return dict(max_speed_m_s=float(np.linalg.norm(q['U'],axis=1).max()),min_T_K=float(q['T'].min()),max_T_K=float(q['T'].max()),mass_per_depth_kg_m=float(q['rho'].sum()*(L/N)**2),eos_relative_error=float(np.max(abs(q['rho']*R*q['T']/q['p']-1))))

def main():
    assert (O/'run-status.txt').read_text().strip()=='complete'
    assert not OUT.exists();OUT.mkdir(parents=True)
    times=[i*.5 for i in range(41)];data={};records=[];initial={}
    for t in times:
        f=reorder(native.fluent_native(F/'qualification-final.cas.h5',F/f'heated-{round(t/DT):04d}.dat.h5',P0));o=ofields(t);data[t]=[f,o]
        for name,q in zip(['ANSYS Fluent','OpenFOAM'],[f,o]):
            s=stats(q)
            if t==0:initial[name]=s['mass_per_depth_kg_m']
            s['relative_mass_drift_from_heating_start']=abs(s['mass_per_depth_kg_m']/initial[name]-1)
            records.append(dict(solver=name,time_s=t,**s))
    audit=dict(specification=dict(specification(),dt_s=DT,heated_steps=HEATED_STEPS),records=records,mapping='Every snapshot mapped by native vertices; unique grid slots checked',physical_energy_terms_alignment_documented=True,energy_limitation='Fluent explicit rho g dot U source verified by the vertical adiabatic control; OpenFOAM includes this source natively. This does not validate physical predictions or identical discretizations.',physical_validation=False,fluent_final_residuals=data[20.][0]['solve_evidence'],openfoam_final_density_source=data[20.][1]['density_source'])
    audit['heated_checks']={name:dict(finite=True,temperature_bounded=all(COLD-1e-4<=s['min_T_K']<=s['max_T_K']<=HOT+1e-4 for s in records if s['solver']==name),mass=all(s['relative_mass_drift_from_heating_start']<=1e-5 for s in records if s['solver']==name),eos=all(s['eos_relative_error']<=1e-6 for s in records if s['solver']==name)) for name in initial}
    raw=(O/'solver.log').read_text();last=re.split(r'^Time = 20s?\s*$',raw,flags=re.M)[-1];audit['openfoam_final_step_log']=last
    end=last.split('PIMPLE: Iteration 16')[-1]
    audit['openfoam_last_outer_initial_residuals']={m[0]:float(m[1]) for m in re.findall(r'Solving for (\w+), Initial residual = ([-+\d.eE]+)',end)}
    fu,ou=[x['U'] for x in data[20.]]
    audit['final_comparison']=dict(temperature_max_difference_K=float(np.max(abs(data[20.][0]['T']-data[20.][1]['T']))),velocity_relative_L2_difference=float(np.linalg.norm(fu-ou)/np.linalg.norm(fu)))
    audit['fluent_strict_final_residual_checks']={k:v <= (1e-9 if k=='energy' else 1e-6) for k,v in audit['fluent_final_residuals']['terminal_scaled_residuals'].items()}
    trans=(F/'solver-transcript.txt').read_text(errors='replace')
    step_checks=[]
    for chunk in re.split(r'Updating solution at time level N',trans)[1:]:
        end_time=re.search(r'Flow time = ([0-9.eE+-]+)s, time step = (\d+)',chunk)
        if not end_time:continue
        rows=re.findall(r'^\s*(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+\d+:',chunk,re.M)
        assert rows
        lastrow=rows[-1];r=[float(v) for v in lastrow[1:]]
        step_checks.append(dict(time_s=float(end_time[1]),step=int(end_time[2]),residuals=dict(zip(['continuity','x-velocity','y-velocity','energy'],r)),strict_limits_met=all(v<=lim for v,lim in zip(r,[1e-6,1e-6,1e-6,1e-9]))))
    assert [r['step'] for r in step_checks]==list(range(1,HEATED_STEPS+1))
    audit['fluent_all_step_checks']=step_checks
    flux=(O/'20/phi').read_text();boundary={}
    for wall in WALLS:
        b=re.search(r'\b'+wall+r'\s*\{([^{}]+)\}',flux,re.S)[1]
        m=re.search(r'value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(([^()]*)\)',b,re.S)
        if m:
            v=np.fromstring(m[2],sep=' ');assert len(v)==int(m[1])==N
        else:v=np.full(N,float(re.search(r'value\s+uniform\s+([^;]+)',b)[1]))
        boundary[wall]=dict(net_kg_s=float(v.sum()),sum_absolute_kg_s=float(abs(v).sum()))
    audit['openfoam_final_wall_mass_flux']=boundary
    audit['strict_fluent_all_steps_pass']=all(v['strict_limits_met'] for v in step_checks)
    write(OUT/'numerical-checks.json',audit)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':14,'axes.labelsize':11,'savefig.facecolor':'#f7f9fc'})
    vmax=max(r['max_speed_m_s']*1000 for r in records)*1.03
    x=(np.arange(N)+.5)*L/N*1000;paths=[]
    for i,t in enumerate(times):
        fig,axs=plt.subplots(2,2,figsize=(12.8,10.2),facecolor='#f7f9fc');fig.subplots_adjust(left=.07,right=.88,bottom=.11,top=.87,hspace=.28,wspace=.24)
        for col,(name,q) in enumerate(zip(['ANSYS Fluent','OpenFOAM'],data[t])):
            U=q['U'].reshape(N,N,2);T=q['T'].reshape(N,N);speed=np.linalg.norm(U,axis=2)*1000
            a=axs[0,col];im=a.pcolormesh(np.linspace(0,50,N+1),np.linspace(0,50,N+1),T-T0,cmap='coolwarm',vmin=-.4,vmax=.4,shading='flat',rasterized=True)
            if t>0:a.streamplot(x,x,U[:,:,0],U[:,:,1],color='#24384b',density=.75,linewidth=.65,arrowsize=.8)
            a.set_title(name+' | temperature')
            a=axs[1,col];vel=a.pcolormesh(np.linspace(0,50,N+1),np.linspace(0,50,N+1),speed,cmap='viridis',vmin=0,vmax=vmax,shading='flat',rasterized=True)
            a.quiver(x[2::5],x[2::5],U[2::5,2::5,0]*1000,U[2::5,2::5,1]*1000,color='white',scale=vmax*12,scale_units='width',width=.004)
            a.set_title('Speed and velocity direction')
            for a in axs[:,col]:a.set_aspect('equal');a.set_xlim(0,50);a.set_ylim(0,50);a.set_xlabel('x (mm)');a.set_ylabel('y (mm)')
        fig.colorbar(im,cax=fig.add_axes([.905,.55,.018,.28]),label='Temperature relative to 300 K (K)')
        fig.colorbar(vel,cax=fig.add_axes([.905,.15,.018,.28]),label='Speed (mm/s)')
        fig.suptitle('Heat drives a circulation',x=.07,y=.98,ha='left',fontsize=24,fontweight='bold')
        fig.text(.07,.925,f'Warm left wall / cool right wall | 50 mm cavity | Ra = {RA:,.0f} | t = {t:.2f} s',fontsize=12)
        fig.text(.07,.035,'Actual transient solver fields | 48 x 48 cells | dt = 0.01 s\nSmaller time-step run; see comparison with dt = 0.02 s. Physical validation remains open.',fontsize=10,color='#46566b')
        path=OUT/f'frame-{i:02d}.png';fig.savefig(path,dpi=120);paths.append(path)
        if t==20:fig.savefig(OUT/'paired-final.png',dpi=160);fig.savefig(OUT/'paired-final.pdf')
        plt.close(fig)
    images=[Image.open(p).convert('RGB') for p in paths];images[0].save(OUT/'convection-development.gif',save_all=True,append_images=images[1:],duration=250,loop=0)
    for im in images:im.close()
    write(OUT/'visualization-record.json',dict(source_times_s=times,frame_duration_ms=250,interpolated_states=False,arrows='native velocity with common scale',streamlines='instantaneous native velocity, not particle paths',T_limits_K=[COLD,HOT],speed_limits_mm_s=[0,vmax],sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))
    print(json.dumps(dict(checks=audit['heated_checks'],comparison=audit['final_comparison'],fluent_residuals=audit['fluent_final_residuals'],openfoam_initial_residuals=audit['openfoam_last_outer_initial_residuals']),indent=2))
if __name__=='__main__':main()
