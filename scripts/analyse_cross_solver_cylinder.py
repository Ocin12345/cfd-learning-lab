"""Check a completed matched Fluent/OpenFOAM cylinder calculation.

No CFD is executed. Raw native data stay outside Git; numeric tables are written
to the specified fresh review directory only after all required checks pass.
"""
from pathlib import Path
import argparse,csv,hashlib,json,math,re,shutil

parser=argparse.ArgumentParser()
parser.add_argument('--run-root',type=Path,required=True)
parser.add_argument('--source-mesh',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
a=parser.parse_args()
foam=a.run_root/'openfoam'; fluent=a.run_root/'fluent'
targets=json.loads((foam/'declared-comparison.json').read_text())
source_hash=hashlib.sha256(a.source_mesh.read_bytes()).hexdigest()
assert source_hash==targets['mesh_source_sha256'],'Declared source mesh differs'
assert hashlib.sha256((foam/'input-cylinder.msh').read_bytes()).hexdigest()==source_hash
settings=json.loads((fluent/'loaded-settings.json').read_text())
assert settings['general']['solver']['time']=='steady'
assert settings['general']['solver']['two_dim_space']=='planar'
assert settings['viscous']['model']=='laminar'
rho=settings['material']['density']; mu=settings['material']['viscosity']
assert rho['option']==mu['option']=='value'
assert rho['value']==targets['density_kg_m3']
assert math.isclose(mu['value']/rho['value'],targets['kinematic_viscosity_m2_s'])
refs=settings['references']
for key,target in [('density','density_kg_m3'),('velocity','velocity_m_s'),
                   ('length','diameter_m'),('area','reference_area_m2'),('depth','span_m')]:
    assert refs[key]==targets[target],('Fluent reference mismatch',key)
assert [v['value'] for v in settings['inlet']['momentum']['velocity_components']]==[targets['velocity_m_s'],0]
assert settings['outlet']['momentum']['gauge_pressure']['value']==0
assert not settings['general']['operating_conditions']['gravity']['enable']
def dictionary(path):
    return re.sub(r'/\*.*?\*/|//[^\n]*','',path.read_text(),flags=re.S)
def scalar(text,key):
    return float(re.search(r'\b'+key+r'\s+([-+0-9.eE]+)\s*;',text)[1])
properties=dictionary(foam/'constant/physicalProperties')
assert scalar(properties,'nu')==targets['kinematic_viscosity_m2_s']
assert re.search(r'viscosityModel\s+constant\s*;',properties)
assert re.search(r'simulationType\s+laminar\s*;',dictionary(foam/'constant/momentumTransport'))
control=dictionary(foam/'system/controlDict')
for key,target in [('rhoInf','density_kg_m3'),('magUInf','velocity_m_s'),
                   ('lRef','diameter_m'),('Aref','reference_area_m2')]:
    assert scalar(control,key)==targets[target],('OpenFOAM reference mismatch',key)
assert re.search(r'dragDir\s+\(\s*1\s+0\s+0\s*\)',control)
assert re.search(r'liftDir\s+\(\s*0\s+1\s+0\s*\)',control)
assert re.search(r'patches\s+\(\s*cylinder\s*\)',control)
u=dictionary(foam/'0/U'); p=dictionary(foam/'0/p')
assert re.search(r'inlet\s*\{[^}]*type\s+fixedValue;[^}]*value\s+uniform\s+\(1 0 0\)',u)
assert re.search(r'cylinder\s*\{[^}]*type\s+fixedValue;[^}]*value\s+uniform\s+\(0 0 0\)',u)
assert re.search(r'outlet\s*\{[^}]*type\s+zeroGradient',u)
assert re.search(r'outlet\s*\{[^}]*type\s+fixedValue;[^}]*value\s+uniform\s+0\s*;',p)
assert re.search(r'ddtSchemes\s*\{\s*default\s+steadyState',dictionary(foam/'system/fvSchemes'))
reynolds=targets['velocity_m_s']*targets['diameter_m']/targets['kinematic_viscosity_m2_s']
assert reynolds==targets['reynolds_number']==20
status=json.loads((foam/'run-status.json').read_text())
assert status['stage']=='finished_execution','OpenFOAM execution still active or interrupted'
assert 'Mesh OK.' in (foam/'mesh-check.log').read_text()
mesh_log=(foam/'mesh-check.log').read_text()
cells=int(re.search(r'^\s*cells:\s+(\d+)',mesh_log,re.M)[1])
assert cells==targets['intended_cells']==98304

# Verify that extrusion preserved the original x/y coordinates, not just bounds.
lines=a.source_mesh.read_text().splitlines()
start=next(i for i,s in enumerate(lines) if s.startswith('(10 (1 '))+1
end=next(i for i in range(start,len(lines)) if lines[i]=='))')
original=[tuple(map(float,s.split())) for s in lines[start:end]]
radii=[math.hypot(*xy) for xy in original]
assert math.isclose(min(radii)*2,targets['diameter_m'],abs_tol=1e-9)
assert math.isclose(max(radii),targets['outer_radius_m'],abs_tol=1e-9)
points=[tuple(map(float,m)) for m in re.findall(r'^\(\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s*\)$',(foam/'0/polyMesh/points').read_text(),re.M)]
assert len(points)==2*len(original)
assert set(round(p[2],12) for p in points)=={-.5,.5}
max_xy_error=max(math.hypot(p[0]-xy[0],p[1]-xy[1]) for p,xy in zip(points,original+original))
assert max_xy_error<1e-10,'Coordinate/order preservation failed'

force_by_k={}
for file in sorted((foam/'postProcessing/forceCoefficients').glob('*/forceCoeffs.dat'),key=lambda p:float(p.parent.name)):
    for line in file.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            values=list(map(float,line.split()))
            force_by_k[int(values[0])]=dict(iteration=int(values[0]),Cm=values[1],Cd=values[2],Cl=values[3])
forces=[force_by_k[k] for k in sorted(force_by_k)]
assert len(forces)>=1200
window=forces[-targets['final_force_window']:]
drag=[r['Cd'] for r in window]
cd_range=(max(drag)-min(drag))/abs(sum(drag)/len(drag))
assert cd_range<=targets['cd_relative_range_target'],('OpenFOAM force drift',cd_range)
final_k=forces[-1]['iteration']
assert (foam/str(final_k)/'phi').exists(),'Reconstructed final flux missing'
logs=sorted(foam.glob('solver*.log'))
residual_by_k={}
for file in logs:
    text=file.read_text(errors='replace')
    assert 'FOAM FATAL' not in text
    for k,block in re.findall(r'^Time = ([0-9.eE+\-]+)s?\s*\n(.*?)(?=^Time = |\Z)',text,re.M|re.S):
        values={name:float(value) for name,value in re.findall(r'Solving for (Ux|Uy|p), Initial residual = ([0-9.eE+\-]+)',block)}
        if set(values)=={'Ux','Uy','p'}:residual_by_k[int(float(k))]=values
final_residual=residual_by_k[final_k]
assert max(final_residual.values())<=targets['initial_p_U_residual_learning_target'],final_residual

phi=(foam/str(final_k)/'phi').read_text()
flux={}
for patch in ['inlet','outlet','cylinder']:
    block=re.search(r'^\s*'+patch+r'\s*\{([^}]*)\}',phi,re.M|re.S)[1]
    uniform=re.search(r'value\s+uniform\s+([-+0-9.eE]+)',block)
    if uniform:
        assert float(uniform[1])==0,'Unexpected nonzero uniform flux; requires face count'
        flux[patch]=0.
    else:
        match=re.search(r'value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',block,re.S)
        assert match,('No boundary flux values',patch)
        values=list(map(float,match[2].split()))
        assert len(values)==int(match[1])
        flux[patch]=math.fsum(values)
relative_mass=abs(sum(flux.values()))/abs(flux['inlet'])
assert relative_mass<=targets['relative_boundary_mass_imbalance_target'],relative_mass

fl_summary=json.loads((fluent/'summary.json').read_text())
fl_targets=json.loads((fluent/'declared-targets.json').read_text())
with (fluent/'force-history.csv').open(newline='') as f:
    fl_forces=[{key:float(value) for key,value in row.items()} for row in csv.DictReader(f)]
assert len(fl_forces)==fl_targets['additional_iterations']+1
assert [r['additional_iteration'] for r in fl_forces]==list(range(len(fl_forces)))
assert all(math.isfinite(v) for row in fl_forces for v in row.values())
fl_window=fl_forces[-fl_targets['last_window_samples']:]
fl_drag=[r['drag_coefficient'] for r in fl_window]
fl_lift=[r['lift_coefficient'] for r in fl_window]
fl_cd_range=(max(fl_drag)-min(fl_drag))/abs(math.fsum(fl_drag)/len(fl_drag))
assert fl_cd_range<=fl_targets['drag_relative_range_max']
assert max(fl_lift)-min(fl_lift)<=fl_targets['lift_absolute_range_max']
assert max(map(abs,fl_lift))<=fl_targets['lift_absolute_magnitude_max']
for row in fl_forces:
    n=int(row['additional_iteration'])
    for direction,key in [('drag','drag_coefficient'),('lift','lift_coefficient')]:
        report=(fluent/f'{direction}-{n:03d}.txt').read_text()
        net=re.findall(r'^\s*Net\s+((?:[-+0-9.eE]+\s+){5}[-+0-9.eE]+)\s*$',report,re.M)
        assert net and float(net[-1].split()[-1])==row[key],('Force report mismatch',n,key)
for key,value in fl_forces[-1].items():
    assert value==fl_summary['final'][key],('Stale Fluent summary',key)
assert fl_summary['native_case_saved']
assert (fluent/'cylinder-Re20-replay.cas.h5').is_file()
assert (fluent/'cylinder-Re20-replay.dat.h5').is_file()
fluent_text=(fluent/'solver-transcript.txt').read_text(errors='replace')
fl_rows=re.findall(r'^\s*(\d+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)\s+([0-9.eE+\-]+)',fluent_text,re.M)
fl_residual=list(map(float,fl_rows[-1][1:]))
assert max(fl_residual)<=1e-7
fl_mass=(fluent/'mass-balance.txt').read_text()
fl_net=float(re.search(r'^\s*Net\s+([-+0-9.eE]+)',fl_mass,re.M)[1])
fl_in=float(re.search(r'^\s*inlet\s+([-+0-9.eE]+)',fl_mass,re.M)[1])
assert abs(fl_net/fl_in)<=1e-6
fl_quality=float(re.search(r'Minimum Orthogonal Quality\s*=\s*([-+0-9.eE]+)',fluent_text)[1])
assert fl_quality>.2
fl_cd=fl_forces[-1]['drag_coefficient']
foam_cd=forces[-1]['Cd']
difference=100*abs(foam_cd-fl_cd)/abs(fl_cd)
record=dict(reynolds_number=reynolds,cells=cells,
            outer_radius_over_diameter=max(radii)/targets['diameter_m'],span_m=max(p[2] for p in points)-min(p[2] for p in points),
            mesh_xy_max_difference_m=max_xy_error,source_mesh_sha256=hashlib.sha256(a.source_mesh.read_bytes()).hexdigest(),
            openfoam=dict(final_iteration=final_k,Cd=foam_cd,Cl=forces[-1]['Cl'],
                          final_initial_residuals=final_residual,last100_Cd_relative_range=cd_range,
                          boundary_mass_flux_kg_s=flux,relative_mass_imbalance=relative_mass),
            fluent=dict(final_iteration=int(fl_rows[-1][0]),Cd=fl_cd,Cl=fl_summary['final']['lift_coefficient'],
                        final_scaled_residuals=fl_residual,last20_Cd_relative_range=fl_cd_range,
                        relative_mass_imbalance=abs(fl_net/fl_in),minimum_orthogonal_quality=fl_quality),
            cross_code_Cd_difference_percent=difference,
            cross_code_learning_target_percent=targets['cross_code_cd_difference_learning_target_percent'],
            cross_code_target_pass=difference<=targets['cross_code_cd_difference_learning_target_percent'],
            numerical_checks_passed=True,physical_validation=False,
            note='Distinct residual definitions and initialization; no performance or experimental-validation claim')
assert record['cross_code_target_pass'],('Cross-code difference exceeds declared comparison target',difference)
if a.out.exists():raise FileExistsError('Use a fresh review output directory')
a.out.mkdir(parents=True)
(a.out/'comparison.json').write_text(json.dumps(record,indent=2)+'\n')
with (a.out/'openfoam-force-history.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['iteration','Cm','Cd','Cl']);writer.writeheader();writer.writerows(forces)
with (a.out/'openfoam-initial-residuals.csv').open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=['iteration','Ux','Uy','p']);writer.writeheader()
    writer.writerows(dict(iteration=k,**residual_by_k[k]) for k in sorted(residual_by_k))
shutil.copy2(fluent/'force-history.csv',a.out/'fluent-force-history.csv')
shutil.copy2(foam/'declared-comparison.json',a.out/'declared-comparison.json')
print(json.dumps(record,indent=2))
