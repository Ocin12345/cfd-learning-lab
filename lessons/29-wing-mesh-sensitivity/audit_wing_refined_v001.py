"""Audit saved paired wing refinement; no solver launches or criterion changes."""
from pathlib import Path
import csv,hashlib,json,re,sys
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path('E:/CFD-Learning-Lab');REPO=R/'repo/cfd-learning-lab';L=R/'local/finite-wing-mesh-sensitivity-v001'
sys.path.insert(0,str(REPO/'lessons/27-wing-convergence'))
from audit_wing_continuation_v001 import native
from wing_fluent_v001 import force_rows,force_stability
from wing_refined_foam_fields_v001 import residual_history,forces,boundary_mass,field
from prepare_openfoam_refined_v001 import dictionaries,native_poly,native_patches,mesh_quality
from plot_finite_wing_native_v002 import TetGeometry,read_fluent,read_foam,indexed
F=R/'native/runs/finite-wing-refined-v001/fluent';O=F.parent/'openfoam'
M=R/'native/meshes/finite-wing-refined-v002';B=R/'native/runs/2026-11-01-finite-wing-v003/mesh-export-v002'
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def csvwrite(p,rows):
    with p.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main(out):
    assert not out.exists();out.mkdir(parents=True)
    assert read(F/'status.json')['stage'] in ['completed','completed_with_export_errors']
    assert (O/'run-status.txt').read_text().strip()=='complete'
    files=[F/n for n in ['wing-final.cas.h5','wing-final.dat.h5','wing-setup.dat.h5','force-history.out','final-settings.json','geometry-quality.json','code-provenance.json','declared-targets.json']]
    files += [O/n for n in [*dictionaries(),'1500/U','1500/p','1500/phi','1500/C','solver.log','mesh-check.log','provenance/preflight-validation.json','provenance/prepared-input-sha256.json']]
    files += [M/n for n in ['geometry.npz','setup.json','wing.msh','profile.json','node-relocation-witness.json','mesh-source-manifest.json']]
    fp=list((O/'postProcessing/forceCoefficients').glob('*/forceCoeffs.dat'));assert len(fp)==1;files+=fp
    hashes={str(p):sha(p) for p in files}
    for row in read(R/'local/finite-wing-v003/openfoam-checks-v001/source-manifest.json')['sources']:
        assert sha(Path(row['path']))==row['sha256_after'];hashes[row['path']]=row['sha256_after']
    for p,h in read(R/'exports/figures/2026-11-07-wing-convergence-v001/native-source-hashes-private.json').items():
        assert sha(Path(p))==h;hashes[p]=h
    assert sha(M/'profile.json')==sha(B/'profile.json')
    setup=read(M/'setup.json');base=read(B/'setup.json');witness=read(M/'node-relocation-witness.json')
    assert witness['boundary_coordinates_byte_identical'] and witness['candidate_gate_pass']
    assert witness['connectivity_before_sha256']==witness['connectivity_after_sha256']
    for k in ['chord_m','span_m','reference_area_m2','domain_bounds_m','incidence_degrees_clockwise','density_kg_m3','dynamic_viscosity_Pa_s','inlet_velocity_m_s','Re_c']:
        assert setup[k]==base[k],k
    assert setup['cells']==224972>base['cells']==169800
    assert read(F/'code-provenance.json')['mesh_sha256']==sha(M/'wing.msh')==sha(O/'input-wing.msh')
    f,it,res,g=native(F/'wing-final.cas.h5',F/'wing-final.dat.h5');assert f['iteration']==1500 and f['cells']==setup['cells']
    fr=force_rows((F/'force-history.out').read_text());assert [v['iteration'] for v in fr]==list(range(1,1501))
    f['forces']=fr[-1];f['force_stability']=force_stability(fr,1500);f['mesh_quality']=read(F/'geometry-quality.json')
    assert f['mesh_quality']['minimum_orthogonal_quality']>=.1 and f['mesh_quality']['minimum_cell_volume_m3']>0
    st=read(F/'final-settings.json');assert st['methods']['p_v_coupling']['flow_scheme']=='Coupled'
    assert st['methods']['pseudo_time_method']['formulation']['coupled_solver']=='off'
    assert st['controls']['p_v_controls']['flow_courant_number']==20
    for k in ['inlet','farfield']:assert [v['value'] for v in st[k]['momentum']['velocity_components']]==[1,0,0]
    assert st['references']['area']==4 and st['references']['length']==1 and st['references']['density']==1 and st['references']['velocity']==1
    assert st['viscous']['model']=='laminar' and st['general']['solver']['time']=='steady'
    with h5py.File(F/'wing-setup.dat.h5') as h:
        a=h['results/1/phase-1/cells']
        for name,val in [('SV_U',1.),('SV_V',0.),('SV_W',0.),('SV_P',0.)]:assert np.allclose(indexed(a[name],f['cells']),val,rtol=0,atol=1e-12)
    for label,coef in [('drag','Cd'),('lift','Cl')]:
        p=F/'force-reports'/(label+'-final.txt');txt=p.read_text();hashes[str(p)]=sha(p)
        rows=re.findall(r'^Net\s+((?:'+NUM+r'\s+){5}'+NUM+r')\s*$',txt,re.M)
        assert len(rows)==1 and np.isclose(float(rows[0].split()[-1]),f['forces'][coef],rtol=1e-6,atol=1e-10)
    f['checks']=dict(mesh=True,settings=True,uniform_initialization=True,residuals=max(f['final_residuals'].values())<=1e-5,mass=f['relative_mass_imbalance']<=1e-4,force_stability=f['force_stability']['passed'])
    for name,expected in dictionaries().items():assert (O/name).read_text()==expected,name
    for name,h in read(O/'provenance/prepared-input-sha256.json').items():assert sha(O/name)==h
    poly=native_poly(O);patches=native_patches(poly);quality=mesh_quality((O/'mesh-check.log').read_text());assert quality['cells']==setup['cells']
    for name,h in read(O/'provenance/preflight-validation.json')['native_mesh_sha256'].items():assert sha(poly/name)==h
    n=quality['cells'];ni=min(v['start'] for v in patches.values())
    for t in [300,600,900,1200,1500]:
        for name,nc,dim,count in [('U',3,[0,1,-1,0,0,0,0],n),('p',1,[0,2,-2,0,0,0,0],n),('phi',1,[0,3,-1,0,0,0,0],ni)]:
            p=O/str(t)/name;a=field(p.read_text(),count,nc,dim);hashes[str(p)]=sha(p)
            if name=='U':assert np.linalg.norm(a,axis=1).max()<10
    mass=boundary_mass((O/'1500/phi').read_text(),patches);orr=residual_history((O/'solver.log').read_text());of=forces(fp[0].read_text())
    assert [v['iteration'] for v in of]==list(range(1,1501))
    ofs=force_stability(of,1500);last={k:v for k,v in orr[-1].items() if k!='iteration'}
    o=dict(iteration=1500,forces=of[-1],force_stability=ofs,final_residuals=last,maximum_last100_residuals={k:max(v[k] for v in orr[-100:]) for k in last},mass_balance=mass,mesh_quality=quality,
        checks=dict(mesh=True,settings=True,residuals=max(last.values())<=1e-5,mass=mass['relative_imbalance']<=1e-4,force_stability=ofs['passed']))
    geo=TetGeometry(M/'geometry.npz')
    for a in [read_fluent(F/'wing-final.cas.h5',F/'wing-final.dat.h5'),read_foam(O,'1500','wing')]:geo.match(a)
    old=read(REPO/'lessons/27-wing-convergence/wing-continuation-summary.json')
    compare=[];change={};gap={}
    for name,new,baseline in [('ANSYS Fluent',f,old['cases']['Coupled']['forces']),('OpenFOAM',o,old['openfoam_final_forces'])]:
        compare += [dict(solver=name,mesh='baseline',cells=169800,**{k:baseline[k] for k in ['Cd','Cl']}),dict(solver=name,mesh='refined',cells=224972,**{k:new['forces'][k] for k in ['Cd','Cl']})]
        change[name]={k:100*(new['forces'][k]-baseline[k])/abs(baseline[k]) for k in ['Cd','Cl']}
    for k in ['Cd','Cl']:gap[k]=100*(f['forces'][k]-o['forces'][k])/abs(o['forces'][k])
    s=dict(model='finite-wing-mesh-sensitivity',cases={'ANSYS Fluent':f,'OpenFOAM':o},force_comparison=compare,refinement_change_percent=change,refined_cross_solver_difference_percent_relative_to_openfoam=gap,
        baseline_cells=169800,refined_cells=224972,mesh_size_factor=.9,effective_global_h_ratio=(224972/169800)**(1/3),geometry_and_physical_conditions_unchanged=True,
        native_common_geometry_match=True,boundary_preserving_repair_moved_nodes=witness['moved_node_count'],physical_validation=False,mesh_independence=False,
        limitations='Two modestly separated unstructured grids do not establish asymptotic convergence or a Grid Convergence Index; iterative misses remain explicit.')
    write(out/'wing-refinement-summary.json',s);csvwrite(out/'mesh-forces.csv',compare);csvwrite(out/'fluent-forces.csv',fr);csvwrite(out/'openfoam-forces.csv',of)
    frr=[dict(iteration=int(t),**{k:float(v[i]) for k,v in res.items()}) for i,t in enumerate(it)]
    csvwrite(out/'fluent-residuals.csv',frr);csvwrite(out/'openfoam-residuals.csv',orr)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
    for name,color,marker in [('ANSYS Fluent','#1768a6','o'),('OpenFOAM','#d56629','s')]:
        a=[v for v in compare if v['solver']==name]
        for ax,key,title in zip(axs,['Cl','Cd'],['(a) Lift','(b) Drag']):
            ax.plot([v['cells']/1000 for v in a],[v[key] for v in a],marker=marker,color=color,label=name,lw=1.6)
            ax.set(xlabel='Cells (thousands)',ylabel='$C_L$' if key=='Cl' else '$C_D$',title=title,xticks=[169.8,224.972]);ax.grid(alpha=.15)
    axs[0].legend(frameon=False)
    status='; '.join(name+(' PASS' if all(case['checks'].values()) else ' MISS') for name,case in [('Fluent',f),('OpenFOAM',o)])
    fig.suptitle('Same wing at 8 deg, Re = 400; lines connect two meshes\nRefined numerical checks: '+status,fontsize=10)
    fig.savefig(out/'wing-mesh-forces.png',dpi=170);fig.savefig(out/'wing-mesh-forces.svg');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
    for ax,name,rows in [(axs[0],'ANSYS Fluent',frr),(axs[1],'OpenFOAM',orr)]:
        for key in rows[0]:
            if key!='iteration':ax.semilogy([v['iteration'] for v in rows],[v[key] for v in rows],label=key,lw=1)
        ax.axhline(1e-5,ls=':',color='#687d8b');ax.set(title=name,xlabel='Steady iteration',ylabel='Native residual');ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
    fig.suptitle(status+'\nResidual definitions differ between solvers',fontsize=10)
    fig.savefig(out/'wing-refinement-residuals.png',dpi=160);plt.close(fig)
    assert all(sha(Path(p))==h for p,h in hashes.items());write(out/'native-source-hashes-private.json',hashes)
    print(json.dumps(dict(checks={k:v['checks'] for k,v in s['cases'].items()},changes=change,gap=gap),indent=2))
    return 0 if all(all(v['checks'].values()) for v in s['cases'].values()) else 2
if __name__=='__main__':raise SystemExit(main(Path(sys.argv[1])))
