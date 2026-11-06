"""Read-only numerical audit of the saved finite-wing Fluent calculation."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re
import shutil

import h5py
import numpy as np

EQUATIONS = ('continuity', 'x-velocity', 'y-velocity', 'z-velocity')
NUMBER = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def close(value, expected, tolerance=1e-10):
    try:
        return math.isfinite(float(value)) and abs(float(value)-expected) <= tolerance
    except (TypeError, ValueError):
        return False


def json_write(path, data):
    with path.open('x', encoding='utf-8', newline='\n') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')


def csv_write(path, rows):
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def native_residuals(path):
    result = {}
    common = None
    with h5py.File(path, 'r') as handle:
        for equation in EQUATIONS:
            group = handle['results/residuals/phase-1/'+equation]
            index = np.asarray(group['iterations'][:]).reshape(-1)
            values = np.asarray(group['data'][:])
            assert len(index) and np.array_equal(index, index.astype(int))
            assert np.all(np.diff(index)>0) and values.shape[0] == len(index)
            assert values.ndim == 2 and values.shape[1]>=2
            assert np.isfinite(values[:,:2]).all() and np.all(values[:,1]>0)
            if common is not None:
                assert np.array_equal(common,index)
            common = index
            scaled = values[:,0]/values[:,1]
            assert np.isfinite(scaled).all() and np.all(scaled>=0)
            for iteration, value in zip(index,scaled):
                result.setdefault(int(iteration),dict(iteration=int(iteration)))[equation] = float(value)
    return [result[i] for i in sorted(result)]


def native_fields(path, count):
    result = {}
    with h5py.File(path,'r') as handle:
        root = handle['results/1/phase-1/cells']
        for name, expected in [('SV_U',None),('SV_V',None),('SV_W',None),
                               ('SV_P',None),('SV_DENSITY',1.),('SV_MU_LAM',.0025)]:
            seen = np.zeros(count,dtype=bool)
            low,high = math.inf,-math.inf
            for dataset in root[name].values():
                first = int(np.asarray(dataset.attrs['minId']).ravel()[0])-1
                last = int(np.asarray(dataset.attrs['maxId']).ravel()[0])
                values = np.asarray(dataset[:]).reshape(-1)
                assert 0<=first<last<=count and len(values)==last-first
                assert not seen[first:last].any() and np.isfinite(values).all()
                if expected is not None:
                    assert np.allclose(values,expected,rtol=1e-10,atol=1e-12)
                seen[first:last]=True
                low=min(low,float(values.min()));high=max(high,float(values.max()))
            assert seen.all()
            result[name]=dict(cells=count,minimum=low,maximum=high)
    return result


def native_mass(text):
    result={}
    for name in ('inlet','outlet','farfield','Net'):
        rows=re.findall(r'^\s*'+name+r'\s+('+NUMBER+r')\s*$',text,re.M)
        assert rows, 'Missing native mass report row '+name
        result[name]=float(rows[-1])
    assert all(map(math.isfinite,result.values())) and result['inlet']>0 and result['outlet']<0
    result['relative_imbalance']=abs(result['Net'])/abs(result['inlet'])
    result['definition']='Native net mass flux / inlet mass flux; farfield row included in native net'
    return result


def force_rows(text):
    output=[]
    for line in text.splitlines():
        tokens=line.split()
        if len(tokens)!=3:
            continue
        try:
            iteration,cd,cl=map(float,tokens)
        except ValueError:
            continue
        assert all(map(math.isfinite,(iteration,cd,cl))) and iteration.is_integer()
        if iteration>0:
            output.append(dict(iteration=int(iteration),Cd=cd,Cl=cl))
    assert output and len({row['iteration'] for row in output})==len(output)
    assert all(b['iteration']==a['iteration']+1 for a,b in zip(output,output[1:]))
    return output


def components(momentum):
    assert momentum['velocity_specification_method']=='Components'
    values=momentum['velocity_components']
    assert len(values)==3 and all(v['option']=='value' for v in values)
    return [v['value'] for v in values]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    run,out=args.run.resolve(),args.out.resolve()
    assert run.is_dir() and out.drive.upper()=='E:' and not out.exists()
    out.mkdir(parents=True)
    checks,evidence,errors,before={},{},{},{}

    def remember(path):
        path=Path(path).resolve()
        assert path.is_file(),str(path)
        before.setdefault(path,sha(path))
        return path

    def load(name):
        return json.loads(remember(run/name).read_text(encoding='utf-8-sig'))

    def attempt(name,operation):
        try:
            value=operation();checks[name]=True;return value
        except Exception as error:
            checks[name]=False;errors[name]=type(error).__name__+': '+str(error)
            return None

    status=attempt('execution_status_read',lambda:load('status.json')) or {}
    summary=attempt('summary_read',lambda:load('summary.json')) or {}
    settings=attempt('actual_settings_read',lambda:load('final-settings.json')) or {}
    quality=attempt('native_mesh_quality_read',lambda:load('geometry-quality.json')) or {}
    targets=attempt('declared_targets_read',lambda:load('declared-targets.json')) or {}
    source=attempt('source_provenance_read',lambda:load('code-provenance.json')) or {}
    mesh_dir=Path(source.get('mesh',str(run.parent/'mesh/wing.msh'))).resolve().parent
    checks['recorded_mesh_folder_within_selected_run']=mesh_dir.parent==run.parent
    setup=attempt('shared_geometry_setup_read',lambda:json.loads(remember(mesh_dir/'setup.json').read_text(encoding='utf-8-sig'))) or {}
    evidence.update(status=status,summary=summary,mesh_quality=quality,declared_targets=targets)
    checks['declared_numerical_targets_unchanged']=(targets.get('maximum_total_iterations')==1200
        and close(targets.get('residual_target'),1e-5,1e-12)
        and close(targets.get('relative_mass_imbalance_target'),1e-4,1e-12)
        and close(targets.get('minimum_orthogonal_quality'),.10)
        and targets.get('force_stability_samples')==100
        and close(targets.get('force_stability_relative_range_max'),1e-3,1e-12)
        and targets.get('solver_processes')==2)
    checks['native_execution_completed']=status.get('stage') in ('completed','completed_with_export_errors') and summary.get('native_final_saved') is True
    checks['bounded_shared_mesh_count']=0<setup.get('cells',0)<=250000
    checks['native_mesh_quality_gate']=quality.get('minimum_orthogonal_quality',0)>=.10 and quality.get('minimum_cell_volume_m3',0)>0 and quality.get('checked_before_flow') is True
    def source_check():
        assert sha(remember(run/'executed-script.py'))==source.get('source_sha256')
        return True
    attempt('executed_source_matches_prelaunch_hash',source_check)
    def mesh_binding_check():
        parent=mesh_dir
        assert parent.parent==run.parent
        assert Path(source['mesh']).resolve()==(parent/'wing.msh').resolve()
        assert sha(remember(parent/'wing.msh'))==source['mesh_sha256']==setup['mesh_sha256']
        assert sha(remember(parent/'setup.json'))==source['mesh_setup_sha256']
        assert sha(remember(parent/'geometry.npz'))==source['geometry_npz_sha256']==setup['geometry_npz_sha256']
        assert sha(remember(run/'mesh-setup.json'))==source['mesh_setup_sha256']
        return True
    attempt('native_mesh_and_geometry_match_prelaunch_hashes',mesh_binding_check)

    def settings_check():
        assert settings['general']['solver']['type']=='pressure-based'
        assert settings['general']['solver']['time']=='steady'
        assert settings['viscous']['model']=='laminar'
        assert close(settings['material']['density']['value'],1.) and close(settings['material']['viscosity']['value'],.0025)
        assert settings['cell_zones']['fluid']['general']['material']=='air'
        for patch in ('inlet','farfield'):
            assert all(close(v,e) for v,e in zip(components(settings[patch]['momentum']),(1.,0.,0.)))
        assert close(settings['outlet']['momentum']['gauge_pressure']['value'],0.)
        wall=settings['walls']['wing']['momentum']
        assert wall['wall_motion']=='Stationary Wall' and wall['shear_condition']=='No Slip'
        ref=settings['references']
        for name,value in dict(density=1.,velocity=1.,length=1.,area=4.,viscosity=.0025).items():
            assert close(ref[name],value)
        assert settings['methods']['p_v_coupling']['flow_scheme']=='SIMPLE'
        discretization=settings['methods']['spatial_discretization']['discretization_scheme']
        assert discretization['mom']=='second-order-upwind' and discretization['pressure']=='second-order'
        for equation in EQUATIONS:
            entry=settings['residual']['equations'][equation]
            assert close(entry['absolute_criteria'],1e-5,1e-12) and entry['check_convergence'] is True
        return True
    attempt('actual_settings_match_declared_physics_and_numerics',settings_check)
    data=attempt('native_final_data_present',lambda:remember(run/'wing-final.dat.h5'))
    attempt('native_final_case_present',lambda:remember(run/'wing-final.cas.h5'))
    residual=attempt('native_residual_history_valid',lambda:native_residuals(data)) if data else None
    if residual:
        csv_write(out/'fluent-residuals.csv',residual)
        evidence['residuals']=dict(final_scaled=residual[-1],native_history_records=len(residual),definition='Saved native HDF residual data column0 / column1; no cross-code norm equivalence')
        checks['final_residuals_le_1e_5']=all(residual[-1][name]<=1e-5 for name in EQUATIONS)
        checks['actual_iteration_cap_1200']=0<residual[-1]['iteration']<=1200
        checks['saved_native_iteration_matches_driver_summary']=residual[-1]['iteration']==summary.get('actual_iterations')
    if data and setup.get('cells'):
        evidence['fields']=attempt('all_native_3d_fields_finite_and_materials_match',lambda:native_fields(data,setup['cells']))
    mass=attempt('native_mass_balance_report_valid',lambda:native_mass(remember(run/'mass-balance.txt').read_text(encoding='utf-8',errors='replace')))
    if mass:
        evidence['mass_balance']=mass;checks['relative_mass_imbalance_le_1e_4']=mass['relative_imbalance']<=1e-4
    forces=attempt('native_per_iteration_force_history_valid',lambda:force_rows(remember(run/'force-history.out').read_text(encoding='utf-8',errors='replace')))
    if forces:
        csv_write(out/'fluent-forces.csv',forces)
        def force_monitor_check():
            native=load('native-force-monitor-settings.json')
            assert native['native_column_order']==['iteration','Cd','Cl']
            assert native['report_files']['forces']['report_defs']==['Cd','Cl']
            definitions=native['report_definitions']
            for kind,name,direction in [('drag','Cd',[1,0,0]),('lift','Cl',[0,1,0])]:
                report=definitions[kind][name]
                assert report['zones']==['wing'] and report['force_vector']==direction
            return True
        attempt('native_force_monitor_order_directions_and_wing_zone',force_monitor_check)
        if residual:
            checks['native_force_final_iteration_matches_residual']=forces[-1]['iteration']==residual[-1]['iteration']
        checks['force_stability_has_100_samples']=len(forces)>=100
        stability={}
        for name in ('Cd','Cl'):
            values=np.array([r[name] for r in forces[-100:]])
            metric=float(np.ptp(values)/max(abs(float(values.mean())),1e-3))
            stability[name]=dict(samples=len(values),mean=float(values.mean()),minimum=float(values.min()),maximum=float(values.max()),relative_range=metric)
            checks['last100_'+name+'_relative_range_le_1e_3']=len(values)==100 and metric<=1e-3
        evidence['forces']=dict(samples=len(forces),last=forces[-1],stability=stability,reference_area_m2=4.,reference_speed_m_s=1.,reference_density_kg_m3=1.)
    for file in ('wing.msh','geometry.npz','executed-mesher.py','profile.json'):
        attempt('shared_mesh_source_'+file,lambda file=file:remember(mesh_dir/file))
    for name in ('solver-transcript.txt','checkpoints.json'):
        attempt('native_'+name,lambda name=name:remember(run/name))
    for path,digest in list(before.items()):
        assert sha(path)==digest,'Native input changed: '+str(path)
    checks['inputs_unchanged_during_analysis']=True
    result=dict(solver='ANSYS Fluent Student 2026 R1',analysis_completed=True,check_results=checks,failed_checks=[k for k,v in checks.items() if not v],evidence=evidence,errors=errors,all_declared_checks_passed=all(checks.values()),physical_validation=False,mesh_independence_established=False,unsteady_stability_tested=False,solver_started_or_modified=False)
    json_write(out/'fluent-checks.json',result)
    shutil.copy2(Path(__file__),out/'executed-fluent-analyser.py')
    shutil.copy2(run/'executed-script.py',out/'executed-fluent.py')
    assert sha(out/'executed-fluent-analyser.py')==sha(Path(__file__))
    assert sha(out/'executed-fluent.py')==sha(run/'executed-script.py')
    manifest=dict(sources=[dict(path=str(p),sha256=digest) for p,digest in before.items()],scripts=[dict(path=str(Path(__file__).resolve()),sha256=sha(Path(__file__)))],outputs=[dict(path=p.name,sha256=sha(p)) for p in sorted(out.iterdir()) if p.is_file()],inputs_unchanged_during_analysis=True)
    json_write(out/'fluent-source-manifest.json',manifest)
    print(json.dumps(dict(output=str(out),all_declared_checks_passed=result['all_declared_checks_passed'],failed_checks=result['failed_checks'],errors=errors)))
    return 0 if result['all_declared_checks_passed'] else 2


if __name__=='__main__':
    raise SystemExit(main())
