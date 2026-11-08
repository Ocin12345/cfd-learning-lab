from pathlib import Path
import sys,json,re,datetime,shutil
R=Path('E:/CFD-Learning-Lab');W=R/'local/showcase-improvements-2026-11-08-v001/wing-convergence-repair-v001';L=R/'local/finite-wing-mesh-sensitivity-v001';C=R/'local/finite-wing-continuation-v001';repo=R/'repo/cfd-learning-lab';F=R/'native/runs/finite-wing-refined-v001/fluent';O=F.parent/'openfoam';M=R/'native/meshes/finite-wing-refined-v002'
sys.path[:0]=[str(L),str(C),str(repo/'scripts')]
import numpy as np,h5py
from audit_wing_continuation_v001 import native
from wing_fluent_v001 import available_memory,force_stability
from wing_refined_foam_fields_v001 import residual_history,forces,boundary_mass,field
from prepare_openfoam_refined_v001 import dictionaries,native_poly,native_patches,mesh_quality
from plot_finite_wing_native_v002 import TetGeometry,read_fluent,read_foam
import hashlib

def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(n,v):
 p=W/n;assert not p.exists(),p;p.write_text(json.dumps(v,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
assert available_memory()>=4*2**30 and shutil.disk_usage(R).free>=10*2**30
# Explicit persistent scopes, never a private/runtime-cleanup/server glob.
paths=[F/n for n in ['wing-final.cas.h5','wing-final.dat.h5','wing-setup.cas.h5','wing-setup.dat.h5','force-history.out','final-settings.json','geometry-quality.json','code-provenance.json','declared-targets.json','native-force-monitor-settings.json','mass-balance.txt','solver-transcript.txt','status.json']]
paths += list((F/'force-reports').glob('*-final.txt'))
paths += [M/n for n in ['geometry.npz','setup.json','wing.msh','profile.json','node-relocation-witness.json','mesh-source-manifest.json']]
paths += [O/n for n in [*dictionaries(),'input-wing.msh','solver.log','mesh-check.log','provenance/preflight-validation.json','provenance/prepared-input-sha256.json','run-status.txt']]
poly=native_poly(O);paths += [p for p in poly.iterdir() if p.is_file()]
for t in [300,600,900,1200,1500]:paths += [O/str(t)/n for n in ['U','p','phi']]
paths += [O/'1500/C'];fp=list((O/'postProcessing/forceCoefficients').glob('*/forceCoeffs.dat'));assert len(fp)==1;paths+=fp
paths += [L/n for n in ['wing_fluent_refined_v001.py','audit_wing_refined_v001.py','wing_refined_foam_fields_v001.py','prepare_openfoam_refined_v001.py','declared-comparison.json','DESIGN.md']]
paths += [C/'audit_wing_continuation_v001.py',repo/'scripts/wing_fluent_v001.py',repo/'scripts/plot_finite_wing_native_v002.py']
for n in ['27-wing-convergence','28-wing-incidence','29-wing-mesh-sensitivity']:paths.append(repo/'lessons'/n/'README.md')
bindings={str(p):sha(p) for p in sorted(set(paths))};dump('original-persistent-inputs-v001.json',bindings)
baseline={str(p.relative_to(repo)).replace('\\','/'):sha(p) for p in repo.rglob('*') if p.is_file() and not set(p.relative_to(repo).parts)&{'.git','.superpowers','__pycache__'}};assert len(baseline)==2798;dump('repository-before-wing-improvement-v001.json',baseline)
f,it,res,g=native(F/'wing-final.cas.h5',F/'wing-final.dat.h5');assert f['iteration']==1500 and f['cells']==224972
raw=[]
for line in (F/'force-history.out').read_text().splitlines():
 a=line.split()
 if len(a)!=3:continue
 try:i,cd,cl=map(float,a)
 except ValueError:continue
 assert i.is_integer() and np.isfinite([i,cd,cl]).all();raw.append(dict(iteration=int(i),Cd=cd,Cl=cl))
assert [r['iteration'] for r in raw]==list(range(1,1501));f['forces']=raw[-1];f['force_stability']=force_stability(raw,1500);f['continuity_windows']={str(k):dict(min=float(res['continuity'][it>=k].min()),max=float(res['continuity'][it>=k].max()),mean=float(res['continuity'][it>=k].mean())) for k in [900,1200,1401]}
st=json.loads((F/'final-settings.json').read_text());f['pressure_coupled_AMG']=st['controls']['advanced']['multi_grid']['mg_controls']['pressure-coupled'];assert f['pressure_coupled_AMG']['termination_criteria']==.1 and f['pressure_coupled_AMG']['cycle_type']=='f-cycle'
assert st['controls']['p_v_controls']['flow_courant_number']==20 and st['methods']['p_v_coupling']['flow_scheme']=='Coupled'
assert st['methods']['pseudo_time_method']['formulation']['coupled_solver']=='off'
assert st['references']['area']==4 and st['references']['density']==st['references']['velocity']==st['references']['length']==1
for patch in ['inlet','farfield']:assert [v['value'] for v in st[patch]['momentum']['velocity_components']]==[1,0,0]
assert st['material']['density']['value']==1 and st['material']['viscosity']['value']==.0025
assert st['outlet']['momentum']['gauge_pressure']['value']==0 and st['walls']['wing']['momentum']['wall_motion']=='Stationary Wall' and st['walls']['wing']['momentum']['shear_condition']=='No Slip'
for p,t in dictionaries().items():assert (O/p).read_text()==t,p
prepared=json.loads((O/'provenance/prepared-input-sha256.json').read_text());assert all(sha(O/p)==h for p,h in prepared.items())
preflight=json.loads((O/'provenance/preflight-validation.json').read_text());assert preflight['passed'];assert all(sha(poly/p)==h for p,h in preflight['native_mesh_sha256'].items())
patch=native_patches(poly);quality=mesh_quality((O/'mesh-check.log').read_text());assert quality['cells']==224972
mass=boundary_mass((O/'1500/phi').read_text(),patch);orr=residual_history((O/'solver.log').read_text());of=forces(fp[0].read_text());assert [x['iteration'] for x in of]==list(range(1,1501));ofs=force_stability(of,1500);last={k:v for k,v in orr[-1].items() if k!='iteration'};assert max(last.values())<=1e-5 and mass['relative_imbalance']<=1e-4 and ofs['passed']
geo=TetGeometry(M/'geometry.npz');geo.match(read_fluent(F/'wing-final.cas.h5',F/'wing-final.dat.h5'));geo.match(read_foam(O,'1500','wing'))
assert sha(O/'input-wing.msh')==sha(M/'wing.msh')==json.loads((F/'code-provenance.json').read_text())['mesh_sha256']
for p,h in bindings.items():assert sha(p)==h,p
record=dict(observed_at=datetime.datetime.now().astimezone().isoformat(),repository_baseline_files=len(baseline),persistent_native_source_bindings=len(bindings),Fluent=f,OpenFOAM=dict(accepted=True,iteration=1500,final_residuals=last,mass=mass,forces=of[-1],force_stability=ofs,quality=quality),native_common_refined_geometry_match=True,shared_force_denominator_N=2.,physical_conditions_unchanged=True,private_or_transient_metadata_opened_or_copied=False,selected_hypothesis='Same1500checkpoint100iteration baseline versus pressure-coupled AMG termination .1->.001 only; no prediction of success',new_solver_launched=False)
dump('native-evidence-inventory-v001.json',record);print(json.dumps(dict(Fluent_terminal=f['final_residuals'],Fluent_last100=f['last100_max_residuals'],Fluent_cell_mass_abs_sum=f['cell_mass_imbalance_abs_sum_kg_s'],OpenFOAM_pass=True,geometry_match=True,baseline_files=len(baseline),source_bindings=len(bindings))))
