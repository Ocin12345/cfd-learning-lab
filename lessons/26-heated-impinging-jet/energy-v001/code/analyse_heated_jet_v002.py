"""Audit and visualize native thermal jet states. Never launches a solver.

Usage: python -B analyse_heated_jet_v002.py ROOT FRESH_EXPORT_DIRECTORY
Requires repo/scripts, lessons/25-impinging-jet and this directory on PYTHONPATH.
"""
import csv, hashlib, io, json, re, sys
from pathlib import Path
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from analyse_nozzle_fluent_v005 import chunk_array, common_mapping
from analyse_nozzle_openfoam_v002 import field, patch_values
from prepare_openfoam_nozzle_v002 import counted
from jet_common_v001 import NX, NY, DX, HALF, DEPTH, RHO, W, ZONES
from heated_jet_common_v001 import CP, K, COLD, HOT, DT, LAST, STRIDE, ENERGY_BALANCE
from plot_heated_jet_fluent_v001 import fluent_temperatures, panel

N=NX*NY
NF=2*N+NX+NY
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
TDIM=[0,0,0,1,0,0,0]
PDIM=[0,3,-1,0,0,0,0]

def write_json(p,v):
    p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def partial(g,n):
    a=np.full(n,np.nan)
    for v in g.values():
        lo=int(v.attrs['minId'][0])-1;hi=int(v.attrs['maxId'][0])
        assert np.isnan(a[lo:hi]).all()
        a[lo:hi]=v[...]
    return a

def csv_write(p,columns,rows):
    with p.open('x',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(columns);w.writerows(rows)

def values(p,name,kind='cells',n=N):
    with h5py.File(p) as h:return chunk_array(h['results/1/phase-1/'+kind+'/'+name],n)

def fluent(run,parent,geo):
    assert json.loads((run/'status.json').read_text())['stage']=='closed'
    frames,mapping=fluent_temperatures(run,geo)
    assert sorted(frames)==sorted({*range(0,LAST+1,STRIDE),LAST-1})
    with h5py.File(run/'heated-jet.cas.h5') as h:
        b=h['meshes/1/faces'];z=b['zoneTopology']
        zones={n:slice(int(lo)-1,int(hi)) for n,lo,hi in zip(z['name'][0].decode().split(';'),z['minId'][...],z['maxId'][...]) if n!='interior'}
        pts=chunk_array(h['meshes/1/nodes/coords'],(NX+1)*(NY+1),(2,))
        assert list(b['nodes'].keys())==['1'] and (b['nodes/1/nnodes'][...]==2).all()
        face_nodes=b['nodes/1/nodes'][...].reshape(-1,2).astype(int)-1
        face_x=pts[face_nodes].mean(axis=1)[:,0]
    assert set(zones)==set(ZONES)
    pm=common_mapping(parent/'jet-1750.cas.h5',geo)
    assert np.array_equal(pm,mapping)
    frozen={}
    for name,kind,n in [('SV_U','cells',N),('SV_V','cells',N),('SV_P','cells',N),('SV_FLUX','faces',NF)]:
        baseline=values(parent/'jet-1750.dat.h5',name,kind,n)
        frozen[name]=float(max(np.max(abs(values(run/f'heat-{s:04d}.dat.h5',name,kind,n)-baseline)) for s in frames))
    data=run/'heat-0600.dat.h5'
    flux=values(data,'SV_FLUX','faces',NF)
    h0=values(run/'heat-0000.dat.h5','SV_H')
    assert np.ptp(h0)<1e-10
    h300=float(h0[0])
    enthalpy=values(data,'SV_H');temp=values(data,'SV_T')
    assert np.max(abs(enthalpy-h300-CP*(temp-COLD)))<1e-8
    with h5py.File(data) as h:heat=partial(h['results/1/phase-1/faces/SV_HEAT_FLUX'],NF)
    mass={n:float(flux[s].sum()) for n,s in zones.items()}
    raw={n:float(heat[s].sum()) for n,s in zones.items()}
    for file,computed in [('final-heat-report.txt',raw),('final-mass-report.txt',{n:-v for n,v in mass.items()})]:
        text=(run/file).read_text()
        for n,v in computed.items():
            found=re.findall(r'^\s*'+n+r'\s+('+NUM+r')\s*$',text,re.M)
            assert len(found)==1 and np.isclose(float(found[0]),v,rtol=1e-7,atol=1e-11),(file,n,v,found)
    heat300={n:raw[n]+h300*mass[n] for n in ZONES}
    storage=float(RHO*DX**2*np.sum(enthalpy-values(run/'heat-0599.dat.h5','SV_H'))/DT)
    x=face_x[zones['plate']];order=np.argsort(x)
    q=heat[zones['plate']][order]/DX
    assert np.allclose(x[order],geo['centres'][:NX,0],rtol=0,atol=1e-12)
    text=(run/'solver-transcript.txt').read_text()
    terminal=[];last=0
    for m in re.finditer(r'Flow time = ('+NUM+r')s, time step = (\d+)',text):
        rows=re.findall(r'^\s*(\d+)\s+('+NUM+r')\s+\d+:\d+:\d+\s+\d+\s*$',text[last:m.start()],re.M)
        assert rows
        terminal.append((int(m[2]),float(m[1]),float(rows[-1][1])))
        last=m.end()
    assert len(terminal)==LAST and np.array_equal([r[0] for r in terminal],np.arange(1,LAST+1))
    assert np.allclose([r[1] for r in terminal],np.arange(1,LAST+1)*DT,rtol=0,atol=1e-10)
    mesh=json.loads((parent/'mesh-quality.json').read_text())
    return frames,q,dict(frozen_field_max_abs_change=frozen,outward_mass_kg_s_per_m=mass,
        native_heat_in_W_per_m=raw,heat_in_relative_to_300K_W_per_m=heat300,
        native_enthalpy_at_300K_J_kg=h300,storage_rate_W_per_m=storage,
        terminal_energy_residual=float(terminal[-1][2]),maximum_step_terminal_energy_residual=max(r[2] for r in terminal),
        residual_limit=1e-8,residual_definition='Fluent scaled energy residual at the last iteration of each time step',
        actual_steps=LAST,mesh_quality=mesh),terminal

def foam(run,geo):
    assert (run/'run-status.txt').read_text().strip()=='complete'
    before=json.loads((run/'frozen-input-sha256.json').read_text())
    assert all(sha(run/p)==h for p,h in before.items())
    mesh=run/'constant/polyMesh'
    n,b=counted(mesh/'boundary')
    patches={}
    for name,body in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',b):
        patches[name]=dict(count=int(re.search(r'nFaces\s+(\d+)',body)[1]),start=int(re.search(r'startFace\s+(\d+)',body)[1]))
    assert len(patches)==n and set(patches)==set(ZONES)|{'frontAndBack'}
    n,b=counted(mesh/'owner');owner=np.fromstring(b,sep=' ',dtype=int);assert len(owner)==n
    ni=min(r['start'] for r in patches.values())
    c=field(run/'1250/C',N,3,[0,1,0,0,0,0,0])
    i=np.rint((c[:,0]+HALF)/DX-.5).astype(int);j=np.rint(c[:,1]/DX-.5).astype(int)
    mapping=j*NX+i
    assert len(np.unique(mapping))==N and mapping.min()==0 and mapping.max()==N-1
    assert np.max(abs(c[:,:2]-geo['centres'][mapping]))<1e-12
    native={};frames={}
    for s in sorted({*range(0,LAST+1,STRIDE),LAST-1}):
        t=f'{s*DT:g}';v=field(run/t/'T',N,1,TDIM);native[s]=v
        out=np.empty(N);out[mapping]=v;frames[s]=out.reshape(NY,NX)
    assert np.allclose(frames[0],COLD,rtol=0,atol=1e-10)
    frozen={}
    for name,width,dim,num in [('U',3,[0,1,-1,0,0,0,0],N),('p',1,[0,2,-2,0,0,0,0],N),('phi',1,PDIM,ni)]:
        frozen[name]=float(np.max(abs(field(run/'3'/name,num,width,dim)-field(run/'1250'/name,num,width,dim))))
    heat={};mass={};T=native[LAST]
    for name in ZONES:
        p=patches[name];o=owner[p['start']:p['start']+p['count']]
        ph=patch_values(run/'3/phi',name,p['count'])
        baseline=patch_values(run/'1250/phi',name,p['count'])
        frozen['phi_'+name]=float(np.max(abs(ph-baseline)))
        mass[name]=float(ph.sum()*RHO/DEPTH)
        upwind=np.where(ph<0,COLD,T[o])
        convection=-RHO*CP*float(np.sum(ph*(upwind-COLD)))/DEPTH
        if name=='plate': delta=HOT-T[o]
        elif name=='inlet':delta=COLD-T[o]
        elif name in ('left','right'):delta=np.where(ph<0,COLD-T[o],0.)
        else:delta=np.zeros_like(o,dtype=float)
        conduction=float(np.sum(K*delta/(DX/2)*DX))
        heat[name]=convection+conduction
        if name=='plate':
            assert np.max(abs(ph))<1e-15
            x=c[o,0];order=np.argsort(x);q=(K*delta/(DX/2))[order]
            assert np.allclose(x[order],geo['centres'][:NX,0],rtol=0,atol=1e-12)
    storage=float(RHO*CP*DX**2*np.sum(T-native[LAST-1])/DT)
    text=(run/'solver.log').read_text()
    times=np.array([float(x) for x in re.findall(r'^Time = ('+NUM+r')s?\s*$',text,re.M)])
    residuals=np.array([float(x) for x in re.findall(r'Solving for T, Initial residual = '+NUM+r', Final residual = ('+NUM+r')',text)])
    assert len(times)==LAST and len(residuals)==LAST and np.allclose(times,np.arange(1,LAST+1)*DT,rtol=0,atol=1e-10)
    mq=(run/'mesh-check.log').read_text();assert 'Mesh OK.' in mq
    mx=float(re.search(r'Mesh non-orthogonality Max:\s*('+NUM+r')',mq)[1])
    return frames,q,dict(frozen_field_max_abs_change=frozen,outward_mass_kg_s_per_m=mass,
        heat_in_relative_to_300K_W_per_m=heat,storage_rate_W_per_m=storage,
        heat_source='Reconstructed from native owner T and phi using the declared Euler/upwind/orthogonal scheme; not a native heat-transfer report',
        terminal_energy_residual=float(residuals[-1]),maximum_step_terminal_energy_residual=float(residuals.max()),
        residual_limit=1e-12,residual_definition='OpenFOAM final linear T residual in each scalar solve; not comparable in scale to Fluent',
        actual_steps=LAST,mesh_quality=dict(checkMesh='Mesh OK.',max_nonorthogonality_deg=mx)),list(zip(range(1,LAST+1),times,residuals))

def main(root,out):
    assert not out.exists();out.mkdir(parents=True)
    f=root/'native/runs/heated-impinging-jet-fluent-v002'
    o=root/'native/runs/heated-impinging-jet-openfoam-v002'
    parent=root/'native/runs/impinging-jet-fluent-v001'
    with np.load(parent/'geometry.npz') as z:geo={k:z[k].copy() for k in z.files}
    results={'ANSYS Fluent':fluent(f,parent,geo),'OpenFOAM':foam(o,geo)}
    checks={};x=geo['centres'][:NX,0]/W
    for name,(frames,q,record,residuals) in results.items():
        mass=record['outward_mass_kg_s_per_m'];heat=record['heat_in_relative_to_300K_W_per_m']
        record['relative_mass_imbalance']=abs(sum(mass.values()))/abs(mass['inlet'])
        record['net_heat_in_W_per_m']=sum(heat.values())
        record['energy_balance_relative_error']=abs(sum(heat.values())-record['storage_rate_W_per_m'])/max(sum(abs(v) for v in heat.values()),abs(record['storage_rate_W_per_m']))
        record['energy_balance_normalization']='abs(net heat in - storage) / max(sum(abs(patch heat in)),abs(storage))'
        record['temperature_bounds_K']=[float(min(v.min() for v in frames.values())),float(max(v.max() for v in frames.values()))]
        record['plate_heat_W_per_m']=float(q.sum()*DX)
        record['plate_heat_flux_centre_W_m2']=float(q[NX//2-1:NX//2+1].mean())
        record['plate_heat_flux_max_W_m2']=float(q.max())
        record['plate_heat_flux_max_x_W']=float(x[np.argmax(q)])
        record['checks']=dict(frozen_flow=max(record['frozen_field_max_abs_change'].values())<1e-12,
            mass_balance=record['relative_mass_imbalance']<1e-4,
            energy_balance=record['energy_balance_relative_error']<ENERGY_BALANCE,
            bounded_temperature=record['temperature_bounds_K'][0]>=COLD-.01 and record['temperature_bounds_K'][1]<=HOT+.01,
            thermal_residual=record['maximum_step_terminal_energy_residual']<=record['residual_limit'])
        checks[name]=record
        csv_write(out/('fluent-residuals.csv' if name=='ANSYS Fluent' else 'openfoam-residuals.csv'),['step','time_s','terminal_energy_residual'],residuals)
    ft,fq,_,_=results['ANSYS Fluent'];ot,oq,_,_=results['OpenFOAM']
    summary=dict(model='heated-impinging-jet',cells_per_solver=N,Reynolds=100,slot_width_m=W,gap_over_width=3,
        temperature_inlet_K=COLD,temperature_plate_K=HOT,dt_s=DT,end_time_s=LAST*DT,actual_temperature_frames_per_solver=31,
        frozen_velocity=True,buoyancy=False,first_order_time_and_advection=True,
        thermal_steady_state_established=False,mesh_time_independence_established=False,physical_validation=False,
        cases=checks,temperature_field_rms_difference_at_3s_K=float(np.sqrt(np.mean((ft[LAST]-ot[LAST])**2))),
        integrated_plate_heat_difference_percent=float(abs(fq.sum()-oq.sum())/abs(fq.sum())*100),
        failed_attempt='OpenFOAM v001 stopped before its first T solve: missing TFinal; corrected in v002 with TFinal inheriting T. MPI FIFO storage uses /dev/shm because the E: mount lacks FIFO support; case data remain on E:.')
    write_json(out/'heated-jet-summary.json',summary)
    csv_write(out/'plate-heat-flux.csv',['x_over_W','fluent_heat_flux_W_m2','openfoam_heat_flux_W_m2'],zip(x,fq,oq))
    csv_write(out/'thermal-history.csv',['step','time_s','fluent_mean_temperature_K','openfoam_mean_temperature_K'],[(s,s*DT,float(ft[s].mean()),float(ot[s].mean())) for s in range(0,LAST+1,STRIDE)])
    fig=plt.figure(figsize=(11.5,7),facecolor='white')
    gs=fig.add_gridspec(2,2,width_ratios=[1,.025],left=.085,right=.92,top=.83,bottom=.16,hspace=.3,wspace=.06)
    axes=[fig.add_subplot(gs[i,0]) for i in range(2)]
    ims=[panel(ax,frames[LAST],name,xlabel=i==1) for i,(ax,(name,(frames,*_))) in enumerate(zip(axes,results.items()))]
    fig.colorbar(ims[0],cax=fig.add_subplot(gs[:,1])).set_label('Air temperature (K)')
    fig.text(.085,.947,'A cold jet meets a hot plate',size=19,weight='bold')
    fig.text(.085,.897,'300 K inlet  |  320 K plate  |  Re = 100  |  fixed velocity fields',color='#426178')
    clock=fig.text(.085,.09,'Actual temperature fields  |  t = 3.00 s',size=11)
    fig.text(.085,.04,'Orange: hot plate. Cyan: inlet. Single mesh; no experimental validation.',size=9,color='#426178')
    fig.savefig(out/'heated-jet-comparison.png',dpi=150)
    images=[]
    for s in range(0,LAST+1,STRIDE):
        for im,frames in zip(ims,[ft,ot]):im.set_data(frames[s])
        clock.set_text(f'Actual temperature fields  |  t = {s*DT:.2f} s')
        fig.canvas.draw();images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()))
    images[0].save(out/'heated-jet-paired.gif',save_all=True,append_images=images[1:],duration=[700]+[110]*29+[1500],loop=0,disposal=2)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
    ax.plot(x,fq,label='ANSYS Fluent: native wall heat flux',color='#1768a6',lw=2)
    ax.plot(x,oq,label='OpenFOAM: reconstructed from temperature',color='#d56629',ls='--',lw=1.8)
    ax.set(xlabel='$x/W$',ylabel='Heat flux from plate to air (W m$^{-2}$)',xlim=(-6,6),ylim=(0,None),title='Plate heat transfer at t = 3 s')
    ax.grid(alpha=.17);ax.legend(frameon=False,fontsize=9)
    ax.text(.5,.97,'Edge peaks are influenced by the cold side backflow boundary.',transform=ax.transAxes,ha='center',va='top',fontsize=9)
    fig.savefig(out/'plate-heat-flux.png',dpi=165);fig.savefig(out/'plate-heat-flux.svg');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10.5,3.8),layout='constrained')
    for ax,(name,(_,_,rec,rows)) in zip(axs,results.items()):
        a=np.array(rows);ax.semilogy(a[:,1],a[:,2],color='#1768a6' if name=='ANSYS Fluent' else '#d56629',lw=1.3)
        ax.axhline(rec['residual_limit'],ls=':',color='#687d8b',label='Configured criterion')
        ax.set(xlabel='Time (s)',ylabel='Terminal energy residual',title=name,xlim=(0,3));ax.grid(alpha=.17);ax.legend(frameon=False,fontsize=8)
    fig.suptitle('Thermal equation convergence at every step (different solver definitions)',fontsize=11)
    fig.savefig(out/'thermal-residuals.png',dpi=150);plt.close(fig)
    with Image.open(out/'heated-jet-paired.gif') as im:
        digests=[]
        for i in range(im.n_frames):im.seek(i);digests.append(hashlib.sha256(im.convert('RGB').tobytes()).hexdigest())
        assert im.n_frames==31 and len(set(digests))==31
    print(json.dumps(summary,indent=2))
    assert all(all(r['checks'].values()) for r in checks.values()),'Numerical checks failed; inspect saved summary'

if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]))
