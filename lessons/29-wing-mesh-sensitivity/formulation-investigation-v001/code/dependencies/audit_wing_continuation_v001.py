"""Independent read-only comparison of the preserved wing solver experiments."""
import csv,hashlib,json,re,sys
from pathlib import Path
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wing_fluent_v001 import force_rows,force_stability
from plot_finite_wing_native_v002 import read_fluent_geometry,indexed

R=Path('E:/CFD-Learning-Lab');PARENT=R/'native/runs/2026-11-01-finite-wing-v003/fluent'
TIGHT=R/'native/runs/finite-wing-continuation-fluent-v001'
COUPLED=R/'native/runs/finite-wing-coupled-fluent-v001'
EQ=('continuity','x-velocity','y-velocity','z-velocity')
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'

def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def native(case,data):
    g=read_fluent_geometry(case);n=len(g['tets']);nf=len(g['faces'])
    with h5py.File(data) as h:
        residual={};scale={};iterations=None
        for name in EQ:
            a=h['results/residuals/phase-1/'+name];it=a['iterations'][...].reshape(-1).astype(int);v=a['data'][...]
            assert np.all(np.diff(it)==1) and np.isfinite(v[:,:2]).all() and (v[:,1]>0).all()
            if iterations is not None:assert np.array_equal(it,iterations)
            iterations=it;residual[name]=v[:,0]/v[:,1];scale[name]=float(v[-1,1])
        fields=h['results/1/phase-1/cells'];U=np.column_stack([indexed(fields[s],n) for s in ('SV_U','SV_V','SV_W')]);p=indexed(fields['SV_P'],n)
        assert np.isfinite(U).all() and np.isfinite(p).all()
        assert np.max(np.linalg.norm(U,axis=1))<10
        assert np.allclose(indexed(fields['SV_DENSITY'],n),1.,rtol=0,atol=1e-12)
        assert np.allclose(indexed(fields['SV_MU_LAM'],n),.0025,rtol=0,atol=1e-12)
        flux=indexed(h['results/1/phase-1/faces/SV_FLUX'],nf)
    c0,c1=g['owner'],g['neighbour'];internal=c1>=0
    balance=np.bincount(c0,weights=flux,minlength=n)-np.bincount(c1[internal],weights=flux[internal],minlength=n)
    boundary=~internal;inflow=-float(flux[boundary & (flux<0)].sum());net=float(flux[boundary].sum())
    return dict(iteration=int(iterations[-1]),final_residuals={k:float(v[-1]) for k,v in residual.items()},
        last100_max_residuals={k:float(v[-100:].max()) for k,v in residual.items()},native_scaling_denominators=scale,
        boundary_mass_in_kg_s=inflow,net_outward_mass_kg_s=net,relative_mass_imbalance=abs(net)/inflow,
        cell_mass_imbalance_abs_sum_kg_s=float(abs(balance).sum()),cell_mass_imbalance_max_kg_s=float(abs(balance).max()),
        maximum_speed_m_s=float(np.linalg.norm(U,axis=1).max()),cells=n),iterations,residual,g

def main(out):
    assert not out.exists();out.mkdir(parents=True)
    assert all(read(r/'status.json')['stage']=='closed' for r in [TIGHT,COUPLED])
    for r in [TIGHT,COUPLED]:
        h=read(r/'source-hashes.json');assert all(sha(PARENT/n)==v==sha(r/'source'/n) for n,v in h.items())
    specs=[('Original SIMPLE',PARENT,'wing-final','force-history.out'),
           ('Extra 100 SIMPLE',TIGHT,'baseline-1300','forces-history.out'),
           ('Tighter pressure AMG',TIGHT,'tight-pressure-1600','forces-history.out'),
           ('Coupled',COUPLED,'wing-final','forces-history.out')]
    cases={};hist={};force_hist={};geometry=None;source_hashes={}
    for name,run,stem,forces in specs:
        cas=run/(stem+'.cas.h5');dat=run/(stem+'.dat.h5');result,it,res,g=native(cas,dat)
        if geometry is None:geometry=g
        else:
            for k in ['points','tets','faces','owner','neighbour']:
                assert np.array_equal(g[k],geometry[k]),'Native topology changed: '+k
        count=result['iteration'];rows=[v for v in force_rows((run/forces).read_text()) if v['iteration']<=count]
        assert rows[-1]['iteration']==count
        result['force_stability']=force_stability(rows,count);result['forces']=rows[-1]
        result['checks']=dict(mesh_and_native_order_unchanged=True,finite_fields_and_materials=True,
            residuals=max(result['final_residuals'].values())<=1e-5,
            mass=result['relative_mass_imbalance']<=1e-4,force_stability=result['force_stability']['passed'])
        cases[name]=result;hist[name]=(it,res);force_hist[name]=rows
        for p in [cas,dat,run/forces]:source_hashes[str(p)]=sha(p)
    mesh=read(COUPLED/'mesh-quality.json');assert mesh['minimum_orthogonal_quality']>=.1 and mesh['minimum_cell_volume_m3']>0
    loaded=read(COUPLED/'loaded-methods.json');final=read(COUPLED/'final-methods.json')
    assert loaded['spatial_discretization']==final['spatial_discretization']
    assert final['p_v_coupling']['flow_scheme']=='Coupled' and final['pseudo_time_method']['formulation']['coupled_solver']=='off'
    setup=read(COUPLED/'loaded-setup.json')
    assert setup['general']['solver']['time']=='steady' and setup['general']['solver']['type']=='pressure-based'
    assert setup['material']['density']['value']==1 and setup['material']['viscosity']['value']==.0025
    assert setup['references']['area']==4 and setup['references']['velocity']==1
    target=read(COUPLED/'declared-test.json');assert target['residual_target']==1e-5
    # Check the reported native total force coefficient against the monitor.
    for name,coef in [('drag','Cd'),('lift','Cl')]:
        text=(COUPLED/(name+'-final.txt')).read_text()
        rows=re.findall(r'^Net\s+((?:'+NUM+r'\s+){5}'+NUM+r')\s*$',text,re.M)
        assert len(rows)==1
        v=float(rows[0].split()[-1]);assert np.isclose(v,cases['Coupled']['forces'][coef],rtol=1e-7,atol=1e-10)
    mtext=(COUPLED/'mass-balance.txt').read_text()
    net=float(re.findall(r'^\s*Net\s+('+NUM+r')\s*$',mtext,re.M)[0])
    assert np.isclose(-net,cases['Coupled']['net_outward_mass_kg_s'],rtol=1e-4,atol=1e-12)
    # Fresh hash check connects the reused OpenFOAM audit to its unchanged native evidence.
    source=read(R/'local/finite-wing-v003/openfoam-checks-v001/source-manifest.json')
    for row in source['sources']:assert sha(Path(row['path']))==row['sha256_after']
    foam=read(R/'repo/cfd-learning-lab/results/finite-wing-v003/openfoam/checks.json')
    assert foam['all_declared_checks_passed'] and not foam['failed_checks']
    compare={n:100*(cases['Coupled']['forces'][n]-foam['final_forces'][n])/abs(foam['final_forces'][n]) for n in ['Cd','Cl']}
    summary=dict(model='finite-wing-solver-continuation',geometry='Unchanged 169800-cell Re400 finite wing at 8 degrees',
        cases=cases,mesh_quality=mesh,openfoam_reused=True,openfoam_new_iterations=0,openfoam_source_hashes_reverified=len(source['sources']),
        openfoam_final_forces=foam['final_forces'],openfoam_checks_passed=True,
        paired_force_difference_percent_relative_to_openfoam=compare,
        coupled_additional_iterations=300,pressure_amg_test_additional_iterations=400,
        physical_validation=False,mesh_independence=False,pressure_amg_hypothesis_sufficient=False,
        conclusion='Coupled iteration meets the original Fluent numerical targets; tighter pressure AMG alone did not. The force discrepancy between codes remains unresolved.',
        residuals_note='Native scaled residuals are computed from stored numerator and denominator, with scaling recorded per algorithm. Residual histories are iteration diagnostics, not physical time.')
    write(out/'wing-continuation-summary.json',summary)
    write(out/'native-source-hashes-private.json',source_hashes)
    with (out/'residual-history.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['case','iteration',*EQ])
        for name,(it,res) in hist.items():
            for k,j in enumerate(it):w.writerow([name,j,*[res[n][k] for n in EQ]])
    with (out/'force-history.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['case','iteration','Cd','Cl'])
        for name,rows in force_hist.items():
            for v in rows:w.writerow([name,v['iteration'],v['Cd'],v['Cl']])
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,3,figsize=(14.5,4.5),layout='constrained')
    colors={'Original SIMPLE':'#87939b','Extra 100 SIMPLE':'#87939b','Tighter pressure AMG':'#d16b25','Coupled':'#1768a6'}
    for name,(it,res) in hist.items():
        first=1000 if name=='Original SIMPLE' else (1200 if name in ['Coupled','Extra 100 SIMPLE'] else 1300)
        valid=it>=first
        axes[0].semilogy(it[valid],res['continuity'][valid],color=colors[name],lw=1.6,label=name)
        rows=[v for v in force_hist[name] if v['iteration']>=first]
        for ax,coef in zip(axes[1:],['Cd','Cl']):ax.plot([v['iteration'] for v in rows],[v[coef] for v in rows],color=colors[name],lw=1.6)
    axes[0].axhline(1e-5,color='#414952',ls=':',lw=1.2,label='Original residual target')
    axes[0].set(ylabel='Scaled continuity residual',title='(a) Convergence');axes[0].legend(frameon=False,fontsize=8)
    for ax,coef,limits,title in zip(axes[1:],['Cd','Cl'],[(.20,.24),(.26,.32)],['(b) Drag coefficient','(c) Lift coefficient']):
        ax.axhline(foam['final_forces'][coef],color='#6c8243',ls='--',lw=1.5,label='Saved OpenFOAM result')
        ax.set(ylabel='$C_D$' if coef=='Cd' else '$C_L$',ylim=limits,title=title);ax.legend(frameon=False,fontsize=8)
    for ax in axes:ax.set(xlabel='Steady solver iteration',xlim=(1000,1600));ax.grid(alpha=.15)
    fig.suptitle('Same wing and mesh: convergence improves, cross-solver force differences remain',fontsize=13)
    fig.savefig(out/'wing-convergence-comparison.png',dpi=145);fig.savefig(out/'wing-convergence-comparison.svg');plt.close(fig)
    assert all(cases['Coupled']['checks'].values())
    print(json.dumps(dict(coupled=cases['Coupled'],force_difference_percent=compare,openfoam_hashes_verified=len(source['sources'])),indent=2))

if __name__=='__main__':main(Path(sys.argv[1]))
