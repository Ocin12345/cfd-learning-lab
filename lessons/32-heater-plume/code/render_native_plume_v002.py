"""Audit and render saved heater-plume fields; never launch or alter a solver."""
from pathlib import Path
import argparse, hashlib, json, re, shutil
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from plume_production_common_v001 import *
from plume_native_fields_v002 import fluent, foam, metrics

NAMES=['ANSYS Fluent','OpenFOAM']
COLORS=['#1261a0','#c85322']

def native_fluent_time(path):
    with h5py.File(path) as h: raw=h['settings/Data Variables'][()][0].decode()
    def val(key):
        a=re.findall(r'^\('+re.escape(key)+r'\s+([-+0-9.eE]+)\)\s*$',raw,re.M)
        assert len(a)==1,(path,key,a)
        return float(a[0])
    return dict(time_s=val('flow-time'),step=val('time-step'),dt_s=val('physical-time-step'))

def native_foam_time(path,t):
    raw=path.read_text()
    def val(key):
        a=re.findall(r'^\s*'+key+r'\s+([-+0-9.eE]+);',raw,re.M)
        assert len(a)==1,(path,key)
        return float(a[0])
    d=dict(time_s=val('value'),step=val('index'),dt_s=val('deltaT'))
    assert abs(d['time_s']-t)<1e-8 and abs(d['step']*DT-t)<1e-8 and abs(d['dt_s']-DT)<1e-12,d
    return d

def load(args):
    qualification=args.qualification
    local=ROOT/'local/hot-patch-plume-v001' if qualification else ROOT/'local/hot-patch-plume-production-v001'
    base=ROOT/'native/runs'/('hot-patch-plume-qualification-v001' if qualification else 'hot-patch-plume-20s-v001')
    F=base/('fluent-v003' if qualification else 'fluent'); O=base/'openfoam'
    duration=1. if qualification else 20.; interval=.2 if qualification else .5
    steps=round(duration/DT); limits=np.array([1e-6 if qualification else 1e-8,1e-6,1e-6,1e-9])
    assert json.loads((F/'status.json').read_text())['stage']=='closed'
    assert (O/'run-status.txt').read_text().strip()=='complete'
    fs=json.loads((F/'time-step-checks.json').read_text()); os=json.loads((local/'openfoam-steps.json').read_text())
    assert len(fs)==steps and [r['step'] for r in fs]==list(range(1,steps+1))
    assert all(r['passed'] and np.all(np.array(r['residuals'])<=limits) for r in fs)
    assert os['all_passed'] and len(os['steps'])==steps
    assert all(r['passed'] and all(v<=os['targets'][k] for k,v in r['residuals'].items()) for r in os['steps'])
    assert np.allclose([r['time_s'] for r in fs],np.arange(1,steps+1)*DT,atol=1e-8,rtol=0)
    assert np.allclose([r['time_s'] for r in os['steps']],np.arange(1,steps+1)*DT,atol=1e-8,rtol=0)
    assert re.search(r'^End\s*$',(O/'solver.log').read_text(),re.M)
    times=np.round(np.arange(round(duration/interval)+1)*interval,8)
    rows={n:[] for n in NAMES}; fields={n:{} for n in NAMES}; sources=[]
    for n in NAMES:
        frames=json.loads((F/'frames.json').read_text()) if n==NAMES[0] else [dict(time_s=float(t)) for t in times]
        for frame in frames:
            t=round(frame['time_s'],8)
            if n==NAMES[0]:
                p=F/frame['file']; meta=native_fluent_time(p)
                assert abs(meta['time_s']-t)<1e-8 and abs(meta['step']*DT-t)<1e-8 and abs(meta['dt_s']-DT)<1e-12,meta
                q=fluent(F/'qualification-final.cas.h5',p); sources.append(p)
            else:
                p=O/f'{t:g}'; q=foam(O,t)
                meta=native_foam_time(p/'uniform/time',t) if t else dict(time_s=0,initial_directory=True)
                sources.extend(p/k for k in ['U','T','p']+(['rho','uniform/time'] if t else []))
            row=dict(time_s=t,native_time=meta,**metrics(q)); rows[n].append(row); fields[n][t]=q
        m0=rows[n][0]['mass_per_depth_kg_m']
        for row in rows[n]:
            row['mass_drift_signed']=row['mass_per_depth_kg_m']/m0-1
            assert abs(row['mass_drift_signed'])<=1e-5 and row['eos_relative_error']<=1e-6
            assert COLD-1e-4<=row['Tmin_K']<=row['Tmax_K']<=HOT+1e-4
        assert set(times).issubset(fields[n]),n
    comparisons=[]
    for t in times:
        f,o=[fields[n][float(t)] for n in NAMES]
        comparisons.append(dict(time_s=float(t),max_T_difference_K=float(np.max(abs(f['T']-o['T']))),velocity_relative_L2=float(np.linalg.norm(f['U']-o['U'])/max(np.linalg.norm(f['U']),1e-30))))
    mesh=O/'constant/polyMesh'
    if not mesh.exists():mesh=O/'0/polyMesh'
    b=(mesh/'boundary').read_text(); assert int(re.search(r'heater\s*\{[^{}]*nFaces\s+(\d+)',b,re.S)[1])==12
    sources.extend([F/'qualification-final.cas.h5',F/'time-step-checks.json',F/'heated-settings.json',F/'heated-gravity-work-source.json',F/'source-hashes.json',O/'solver.log',O/'system/controlDict',O/'system/fvSolution',O/'constant/physicalProperties',mesh/'boundary',local/'openfoam-steps.json',Path(__file__)])
    spec=specification(); spec.update(heated_steps=steps,qualification_only=qualification,fluent_max_inner=800 if qualification else INNER)
    audit=dict(specification=spec,all_saved_field_checks_pass=True,accepted_steps={n:steps for n in NAMES},records=rows,comparisons=comparisons,physical_validation=False,mesh_independence=False,time_step_independence=False,native_source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},plot=dict(times_s=times.tolist(),interpolated_frames=False,instantaneous_streamlines_not_particle_paths=True,openfoam_initial_density='derived from p/(R*T); t=0 EOS is not an independent check',geometry='native cell centres and spans bijectively checked against shared 48x80 grid',wall_mass_balance='closed-domain density integral drift; no transient energy balance asserted'))
    return fields,rows,times,fs,os,limits,audit

def render_fields(fields,times,out,qualification):
    x=(np.arange(NX)+.5)*DX*1000; y=(np.arange(NY)+.5)*DY*1000
    xe=np.linspace(0,W*1000,NX+1); ye=np.linspace(0,H*1000,NY+1)
    vmax=max(float(np.linalg.norm(fields[n][float(t)]['U'],axis=1).max()) for n in NAMES for t in times)*1000
    assert vmax>0
    images=[]
    for i,t in enumerate(times):
        fig,axs=plt.subplots(2,2,figsize=(9.8,12.5))
        fig.subplots_adjust(left=.08,right=.83,bottom=.15,top=.87,hspace=.26,wspace=.28)
        for col,n in enumerate(NAMES):
            q=fields[n][float(t)]; U=q['U'].reshape(NY,NX,2)*1000; speed=np.linalg.norm(U,axis=2)
            axs[0,col].pcolormesh(xe,ye,q['T'].reshape(NY,NX)-T0,cmap='inferno',vmin=0,vmax=1,shading='flat',rasterized=True)
            axs[1,col].pcolormesh(xe,ye,speed,cmap='viridis',vmin=0,vmax=vmax,shading='flat',rasterized=True)
            if t>0:
                u=np.ma.masked_where(speed<.005*vmax,U[:,:,0]); v=np.ma.masked_where(speed<.005*vmax,U[:,:,1])
                axs[1,col].streamplot(x,y,u,v,color='white',density=.9,linewidth=.65,arrowsize=.7)
            for row in range(2):
                ax=axs[row,col];ax.set(xlim=(0,W*1000),ylim=(0,H*1000),xlabel='x (mm)',ylabel='y (mm)');ax.set_aspect('equal')
                ax.plot([HEATER_LEFT*1000,HEATER_RIGHT*1000],[0,0],color='#ff4538',lw=4,clip_on=False)
                ax.plot([0,W*1000],[H*1000,H*1000],color='#32b9e8',lw=3,clip_on=False)
            axs[0,col].set_title(n,fontsize=13,fontweight='bold',pad=12)
        fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='inferno'),cax=fig.add_axes([.87,.55,.018,.29]),label='Temperature above 300 K (K)')
        fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,vmax),cmap='viridis'),cax=fig.add_axes([.87,.16,.018,.29]),label='Speed (mm/s)')
        fig.suptitle('Warm air above a small floor heater',x=.08,ha='left',y=.97,fontsize=20,fontweight='bold')
        phase='One-second startup qualification' if qualification else 'Twenty-second plume development'
        fig.text(.08,.922,f'{phase}  |  t = {t:.2f} s',fontsize=12)
        fig.text(.08,.898,'301 K heater (red)  |  300 K ceiling (blue)  |  Ra(H) = 99,238',fontsize=10)
        fig.text(.08,.025,'Actual cell-centred solver fields; fixed scales across both solvers and all frames.\nWhite lines: instantaneous streamlines, hidden below 0.5% of the global maximum speed.\n3,840 cells | dt = 0.01 s | educational calculation; physical validation remains open.',fontsize=9,linespacing=1.5)
        p=out/f'frame-{i:03d}.png';fig.savefig(p,dpi=120)
        if i==len(times)-1:
            fig.savefig(out/'paired-final.png',dpi=180);fig.savefig(out/'paired-final.pdf')
        plt.close(fig)
        with Image.open(p) as im:images.append(im.convert('RGB').copy())
    ms=200 if qualification else 250
    images[0].save(out/'plume-development.gif',save_all=True,append_images=images[1:],duration=ms,loop=0,optimize=False)
    with Image.open(out/'plume-development.gif') as gif:assert gif.n_frames==len(times)
    return dict(speed_scale_mm_s=[0,vmax],temperature_scale_K=[0,1],frames=len(times),frame_duration_ms=ms,physical_interval_s=float(times[1]-times[0]),playback='discrete saved fields, no interpolated motion')

def render_checks(rows,fs,os,limits,out):
    fig,axs=plt.subplots(2,2,figsize=(12,8));fig.subplots_adjust(left=.09,right=.97,bottom=.14,top=.85,wspace=.28,hspace=.4)
    for n,color in zip(NAMES,COLORS):
        r=rows[n];t=[v['time_s'] for v in r]
        axs[0,0].plot(t,[v['Umax_m_s']*1000 for v in r],label=n,color=color)
        axs[0,1].plot(t,[v['mass_drift_signed']*1e6 for v in r],label=n,color=color)
    axs[0,0].set(ylabel='Maximum speed (mm/s)',title='Flow startup');axs[0,0].legend(frameon=False)
    axs[0,1].set(ylabel='Mass change from heating start (ppm)',title='Closed-domain mass check')
    axs[0,1].axhline(10,color='grey',ls='--',lw=.8);axs[0,1].axhline(-10,color='grey',ls='--',lw=.8)
    for j,name in enumerate(['Continuity','x momentum','y momentum','Energy']):axs[1,0].semilogy([v['time_s'] for v in fs],np.maximum([v['residuals'][j]/limits[j] for v in fs],1e-16),label=name,lw=.9)
    for name in os['targets']:axs[1,1].semilogy([v['time_s'] for v in os['steps']],np.maximum([v['residuals'][name]/os['targets'][name] for v in os['steps']],1e-16),label=name,lw=.9)
    for ax,title in zip(axs[1],['Fluent: end-of-step scaled residuals','OpenFOAM: last-outer initial residuals']):
        ax.set(ylabel='Residual / declared limit',title=title);ax.axhline(1,color='black',ls='--',lw=.8);ax.legend(frameon=False,fontsize=8,ncol=2)
    for ax in axs.flat:ax.set_xlabel('Physical time (s)');ax.grid(alpha=.18);ax.spines[['top','right']].set_visible(False)
    fig.suptitle('Numerical checks for the saved plume calculation',fontsize=18)
    fig.text(.09,.035,'Residual definitions differ between solvers; these ratios assess each solver against its own limits.\nPassing these checks establishes numerical acceptance for this run, not mesh independence or physical validation.',fontsize=10)
    fig.savefig(out/'numerical-checks.png',dpi=160);fig.savefig(out/'numerical-checks.pdf');plt.close(fig)

def main():
    p=argparse.ArgumentParser();p.add_argument('--qualification',action='store_true');p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    assert args.out.is_relative_to(ROOT/'exports/figures') and not args.out.exists()
    assert shutil.disk_usage(ROOT).free>2*2**30
    fields,rows,times,fs,os,limits,audit=load(args)
    args.out.mkdir(parents=True);shutil.copyfile(__file__,args.out/'executed-renderer.py')
    audit['rendering']=render_fields(fields,times,args.out,args.qualification)
    render_checks(rows,fs,os,limits,args.out)
    write(args.out/'native-field-audit.json',audit)
    print(json.dumps(dict(output=str(args.out),accepted_steps=audit['accepted_steps'],rendering=audit['rendering'],final_comparison=audit['comparisons'][-1]),indent=2))

if __name__=='__main__':main()

