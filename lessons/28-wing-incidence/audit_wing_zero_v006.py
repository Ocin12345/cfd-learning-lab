"""Native numerical checks for the paired zero-incidence wing; no solver launches."""
from pathlib import Path
import csv,hashlib,json,re,sys
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('E:/CFD-Learning-Lab');REPO=R/'repo/cfd-learning-lab'
sys.path.insert(0,str(REPO/'lessons/27-wing-convergence'))
from audit_wing_continuation_v001 import native
from wing_fluent_v001 import force_rows,force_stability
from analyse_finite_wing_openfoam_v003 import residual_history,forces,boundary_mass,field
from prepare_openfoam_finite_wing_v003 import native_poly,native_patches,mesh_quality
from plot_finite_wing_native_v002 import indexed,TetGeometry,read_fluent,read_foam
from wing_zero_common_v001 import U,LIFT
F0=R/'native/runs/finite-wing-zero-fluent-v001';F=R/'native/runs/finite-wing-zero-fluent-v002';O=R/'native/runs/finite-wing-zero-openfoam-v002'
G=R/'native/runs/2026-11-01-finite-wing-v003/mesh-export-v002/geometry.npz'
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csvwrite(p,rows):
    with p.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def mesh_quality(log):
    command=re.search(r'^Exec\s*:\s*(.+)$',log,re.M)
    assert command and set(command[1].split())=={'checkMesh','-allGeometry','-allTopology'}
    assert 'Mesh OK.' in log and not re.search(r'Failed\s+\d+\s+mesh checks|FOAM FATAL',log)
    patterns={'cells':r'\bcells:\s+(\d+)',
        'minimum_volume_m3':r'Min volume\s*=\s*('+NUM+r')',
        'maximum_nonorthogonality_degrees':r'Mesh non-orthogonality Max:\s*('+NUM+r')',
        'maximum_skewness':r'Max skewness\s*=\s*('+NUM+r')'}
    result={k:float(re.search(v,log)[1]) for k,v in patterns.items()};result['cells']=int(result['cells'])
    assert result['minimum_volume_m3']>0 and result['maximum_nonorthogonality_degrees']<=75 and result['maximum_skewness']<=4
    return result

def main(out):
    assert not out.exists();out.mkdir(parents=True)
    assert read(F/'status.json')['stage']=='closed'
    assert (O/'run-status.txt').read_text().strip()=='complete'
    source=read(F/'source-hashes.json');assert all(sha(Path(p))==h for p,h in source.items())
    assert sha(F/'source/wing.cas.h5')==next(iter(source.values()))
    r,it,res,geom=native(F/'wing-final.cas.h5',F/'wing-final.dat.h5')
    with h5py.File(F0/'wing-initial.dat.h5') as h:
        a=h['results/1/phase-1/cells']
        for name,v in zip(['SV_U','SV_V','SV_W'],U):assert np.allclose(indexed(a[name],r['cells']),v,atol=1e-12,rtol=0)
        assert np.max(abs(indexed(a['SV_P'],r['cells'])))<1e-12
    settings=read(F/'actual-settings.json');assert settings['methods']['p_v_coupling']['flow_scheme']=='Coupled'
    assert settings['viscous']['model']=='laminar' and settings['general']['solver']['time']=='steady'
    for name in ['inlet','farfield']:
        got=[v['value'] for v in settings['boundaries']['velocity_inlet'][name]['momentum']['velocity_components']]
        assert np.allclose(got,U,rtol=0,atol=1e-14)
    for kind,name,vec in [('drag','Cd',U),('lift','Cl',LIFT)]:
        assert np.allclose(settings['reports'][kind][name]['force_vector'],vec,rtol=0,atol=1e-14)
    parent_rows=force_rows((F0/'forces-history.out').read_text())
    restart_rows=force_rows((F/'forces-history.out').read_text())
    assert restart_rows[0]==parent_rows[-1] and restart_rows[0]['iteration']==1200
    frows=parent_rows+restart_rows[1:]
    assert [v['iteration'] for v in frows]==list(range(1,r['iteration']+1))
    r['force_stability']=force_stability(frows,r['iteration']);r['forces']=frows[-1]
    mq=read(F/'mesh-quality.json');assert mq['minimum_orthogonal_quality']>=.1 and mq['minimum_cell_volume_m3']>0;r['mesh_quality']=mq
    for name,key in [('drag','Cd'),('lift','Cl')]:
        text=(F/(name+'-final.txt')).read_text()
        lines=re.findall(r'^Net\s+((?:'+NUM+r'\s+){5}'+NUM+r')\s*$',text,re.M)
        assert len(lines)==1 and np.isclose(float(lines[0].split()[-1]),r['forces'][key],rtol=1e-6,atol=1e-10)
    r['checks']=dict(residuals=max(r['final_residuals'].values())<=1e-5,mass=r['relative_mass_imbalance']<=1e-4,
        force_stability=r['force_stability']['passed'],mesh=True,material_and_bc=True,uniform_initialization_verified=True)
    omesh=native_poly(O);patches=native_patches(omesh);quality=mesh_quality((O/'mesh-check.log').read_text());assert quality['cells']==r['cells']
    copied=read(O/'copied-source-hashes.json');assert all(sha(O/p)==h for p,h in copied.items())
    for p,h in copied.items():
        original=p.replace('constant/polyMesh/','0/polyMesh/') if p.startswith('constant/polyMesh/') else p
        assert sha(R/'native/runs/2026-11-01-finite-wing-v003/openfoam'/original)==h
    for path,expr,vec in [('0/U',r'internalField\s+uniform\s*\(([^)]+)\)',U),
        ('system/controlDict',r'liftDir\s*\(([^)]+)\)',LIFT),('system/controlDict',r'dragDir\s*\(([^)]+)\)',U)]:
        got=np.fromstring(re.search(expr,(O/path).read_text())[1],sep=' ');assert np.allclose(got,vec,atol=1e-14,rtol=0)
    # Exact fixedValue boundary vectors are in the initial U dictionary copied to runtime.
    ut=(O/'0/U').read_text()
    for name in ['inlet','farfield']:
        b=re.search(name+r'\s*\{([^}]+)\}',ut)[1]
        assert 'type fixedValue;' in b
        assert np.allclose(np.fromstring(re.search(r'value\s+uniform\s*\(([^)]+)\)',b)[1],sep=' '),U,atol=1e-14,rtol=0)
    n=r['cells'];ni=min(v['start'] for v in patches.values())
    for t in [300,600,900,1200]:
        u=field((O/str(t)/'U').read_text(),n,3,[0,1,-1,0,0,0,0]);p=field((O/str(t)/'p').read_text(),n,1,[0,2,-2,0,0,0,0])
        phi=field((O/str(t)/'phi').read_text(),ni,1,[0,3,-1,0,0,0,0]);assert np.linalg.norm(u,axis=1).max()<10
    mass=boundary_mass((O/'1200/phi').read_text(),patches)
    residuals=residual_history((O/'solver.log').read_text())
    fp=list((O/'postProcessing/forceCoefficients').glob('*/forceCoeffs.dat'));assert len(fp)==1
    orows=forces(fp[0].read_text());assert [v['iteration'] for v in orows]==list(range(1,1201))
    stability=force_stability(orows,1200)
    ofinal={k:v for k,v in residuals[-1].items() if k!='iteration'}
    o=dict(iteration=1200,final_residuals=ofinal,maximum_last100_residuals={k:max(v[k] for v in residuals[-100:]) for k in ofinal},
        forces=orows[-1],force_stability=stability,mass_balance=mass,mesh_quality=quality,
        checks=dict(residuals=max(ofinal.values())<=1e-5,mass=mass['relative_imbalance']<=1e-4,
                    force_stability=stability['passed'],mesh=True,material_and_bc=True))
    # Both native meshes must match the preserved common tetrahedra, independent of ordering.
    geometry=TetGeometry(G)
    for nat in [read_fluent(F/'wing-final.cas.h5',F/'wing-final.dat.h5'),read_foam(O,'1200','wing')]:geometry.match(nat)
    old=read(REPO/'lessons/27-wing-convergence/wing-continuation-summary.json')
    compare=[dict(solver='ANSYS Fluent',angle_deg=0,**{k:r['forces'][k] for k in ['Cd','Cl']}),
        dict(solver='ANSYS Fluent',angle_deg=8,**{k:old['cases']['Coupled']['forces'][k] for k in ['Cd','Cl']}),
        dict(solver='OpenFOAM',angle_deg=0,**{k:o['forces'][k] for k in ['Cd','Cl']}),
        dict(solver='OpenFOAM',angle_deg=8,**{k:old['openfoam_final_forces'][k] for k in ['Cd','Cl']})]
    s=dict(model='finite-wing-zero-incidence',cells_per_solver=n,body_clockwise_pitch_deg=8,
        effective_incidence_deg=0,freestream_U=U,lift_direction=LIFT,Re=400,domain_and_mesh_unchanged=True,
        cases={'ANSYS Fluent':r,'OpenFOAM':o},incidence_comparison=compare,
        zero_incidence_drag_difference_percent_relative_to_openfoam=100*(r['forces']['Cd']-o['forces']['Cd'])/abs(o['forces']['Cd']),
        zero_incidence_absolute_lift_difference=abs(r['forces']['Cl']-o['forces']['Cl']),
        physical_validation=False,mesh_independence=False,domain_orientation_sensitivity=False,
        fluent_continuation=dict(parent_iteration=1200,parent_continuity_residual=read(F0/'summary.json')['residuals']['continuity'],parent_miss_preserved=True,restart_monitor_overlap_verified=True,changed_flow_courant_from=20,to=5),
        residual_definitions='Fluent saved scaled residuals; OpenFOAM maximum printed initial linear residual per equation per steady iteration; magnitudes not directly comparable.')
    write(out/'wing-zero-summary.json',s)
    csvwrite(out/'incidence-forces.csv',compare);csvwrite(out/'fluent-forces.csv',frows);csvwrite(out/'openfoam-forces.csv',orows)
    csvwrite(out/'openfoam-residuals.csv',residuals)
    frr=[dict(iteration=int(t),**{k:float(v[i]) for k,v in res.items()}) for i,t in enumerate(it)];csvwrite(out/'fluent-residuals.csv',frr)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    for name,color,marker in [('ANSYS Fluent','#1768a6','o'),('OpenFOAM','#d56629','s')]:
        a=[v for v in compare if v['solver']==name]
        for ax,key,title in zip(axs,['Cl','Cd'],['(a) Lift','(b) Drag']):
            ax.plot([v['angle_deg'] for v in a],[v[key] for v in a],marker=marker,color=color,label=name,lw=1.6)
            ax.set(xlabel='Incidence (degrees)',ylabel='$C_L$' if key=='Cl' else '$C_D$',title=title,xticks=[0,8]);ax.grid(alpha=.15)
    axs[0].axhline(0,color='#687d8b',lw=.6);axs[1].set_ylim(bottom=0);axs[0].legend(frameon=False)
    fig.suptitle('Re = 400; lines connect two cases\nFluent at 0 deg misses the residual target',fontsize=11)
    fig.savefig(out/'wing-incidence-forces.png',dpi=165);fig.savefig(out/'wing-incidence-forces.svg');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    for ax,name,rows in [(axs[0],'ANSYS Fluent',frr),(axs[1],'OpenFOAM',residuals)]:
        for key in rows[0]:
            if key!='iteration':ax.semilogy([v['iteration'] for v in rows],[v[key] for v in rows],label=key,lw=1)
        ax.axhline(1e-5,ls=':',color='#687d8b');ax.set(title=name,xlabel='Steady iteration',ylabel='Native residual');ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
    fig.suptitle('Zero-incidence convergence: Fluent MISS; OpenFOAM PASS\nSolver residual definitions differ',fontsize=11)
    fig.savefig(out/'wing-zero-residuals.png',dpi=150);plt.close(fig)
    files=[F/n for n in ['wing-final.cas.h5','wing-final.dat.h5','actual-settings.json','forces-history.out','solver-transcript.txt']]
    files += [F0/'wing-initial.dat.h5',F0/'forces-history.out',F0/'summary.json']
    files += [O/n for n in ['1200/U','1200/p','1200/phi','1200/C','solver.log','mesh-check.log']]+[fp[0],G]
    write(out/'native-source-hashes-private.json',{str(p):sha(p) for p in files})
    print(json.dumps(dict(fluent=r,openfoam=o,comparison=compare),indent=2))
    if not (all(r['checks'].values()) and all(o['checks'].values())):
        print('NUMERICAL CHECK MISS: retained in saved summary; not a parsing failure.');sys.exit(2)
if __name__=='__main__':main(Path(sys.argv[1]))
