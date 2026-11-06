"""Read-only, hash-bound audit of the actual saved planar Euler nozzle Fluent run."""
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

EQUATIONS=('continuity','x-velocity','y-velocity','energy')
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):result.update(block)
    return result.hexdigest()


def write_json(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def write_csv(path,rows):
    with Path(path).open('x',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def near(value,expected,tolerance=1e-9):
    return math.isfinite(float(value)) and abs(float(value)-expected)<=tolerance


def chunk_array(group,count,tail=()):
    result=np.empty((count,*tail));seen=np.zeros(count,bool)
    for dataset in group.values():
        first=int(np.asarray(dataset.attrs['minId']).ravel()[0])-1
        last=int(np.asarray(dataset.attrs['maxId']).ravel()[0]);values=np.asarray(dataset[:])
        assert 0<=first<last<=count and not seen[first:last].any()
        values=values.reshape((last-first,*tail));assert np.isfinite(values).all()
        result[first:last]=values;seen[first:last]=True
    assert seen.all()
    return result


def residual_rows(path):
    rows={};common=None
    with h5py.File(path,'r') as handle:
        for equation in EQUATIONS:
            group=handle['results/residuals/phase-1/'+equation]
            index=np.asarray(group['iterations'][:]).reshape(-1);values=np.asarray(group['data'][:])
            assert len(index) and np.array_equal(index,index.astype(int)) and np.all(np.diff(index)>0)
            assert values.ndim==2 and len(values)==len(index) and values.shape[1]>=2
            assert np.isfinite(values[:,:2]).all() and (values[:,1]>0).all()
            if common is not None:assert np.array_equal(index,common)
            common=index;scaled=values[:,0]/values[:,1]
            assert np.isfinite(scaled).all() and (scaled>=0).all()
            for iteration,value in zip(index,scaled):
                rows.setdefault(int(iteration),dict(iteration=int(iteration)))[equation]=float(value)
    return [rows[key] for key in sorted(rows)]


def mass_rows(text):
    rows=[]
    for line in text.splitlines():
        fields=line.split()
        if len(fields)!=3:continue
        try:iteration,mi,mo=map(float,fields)
        except ValueError:continue
        assert all(map(math.isfinite,(iteration,mi,mo))) and iteration.is_integer()
        if iteration>0:rows.append(dict(iteration=int(iteration),mdot_in_kg_s_per_m_depth=mi,mdot_out_kg_s_per_m_depth=mo))
    assert rows and len({row['iteration'] for row in rows})==len(rows)
    assert all(b['iteration']==a['iteration']+1 for a,b in zip(rows,rows[1:]))
    return rows


def native_mass(text):
    result={}
    for name in ('inlet','outlet','upper_wall','lower_wall','Net'):
        values=re.findall(r'^\s*'+name+r'\s+('+NUMBER+r')\s*$',text,re.M)
        assert values,'Missing actual mass report row '+name
        result[name]=float(values[-1])
    assert all(map(math.isfinite,result.values())) and result['inlet']>0 and result['outlet']<0
    result['reference_depth_m']=1.;result['relative_imbalance']=abs(result['Net'])/abs(result['inlet'])
    result['definition']='Absolute native Net mass flux / absolute inlet mass flux; includes both slip wall rows; depth1m'
    return result


def common_mapping(case,geometry):
    nodes=geometry['nodes'];quads=geometry['quads'];count=len(quads)
    with h5py.File(case,'r') as handle:
        root=handle['meshes/1'];native_nodes=chunk_array(root['nodes/coords'],len(nodes),(2,))
        lookup={tuple(np.round(x,10)):i for i,x in enumerate(nodes)}
        assert len(lookup)==len(nodes)
        node_map=np.array([lookup[tuple(np.round(x,10))] for x in native_nodes])
        assert len(np.unique(node_map))==len(nodes)
        assert np.max(np.abs(nodes[node_map]-native_nodes))<1e-10
        face_count=max(int(np.asarray(group.attrs['maxId']).ravel()[0]) for group in root['faces/nodes'].values())
        c0=chunk_array(root['faces/c0'],face_count).astype(int)
        c1=np.zeros(face_count,int)
        for ds in root['faces/c1'].values():
            first=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;last=int(np.asarray(ds.attrs['maxId']).ravel()[0])
            c1[first:last]=np.asarray(ds[:]).reshape(-1)
        members=[set() for _ in range(count)]
        face_seen=np.zeros(face_count,bool)
        for group in root['faces/nodes'].values():
            first=int(np.asarray(group.attrs['minId']).ravel()[0])-1;last=int(np.asarray(group.attrs['maxId']).ravel()[0])
            nnodes=np.asarray(group['nnodes'][:]).reshape(-1).astype(int)
            flat=np.asarray(group['nodes'][:]).reshape(-1).astype(int)-1
            assert len(nnodes)==last-first and np.all(nnodes==2) and len(flat)==2*len(nnodes)
            assert not face_seen[first:last].any();face_seen[first:last]=True
            for face,native in zip(range(first,last),flat.reshape(-1,2)):
                mapped=node_map[native]
                for cid in (c0[face],c1[face]):
                    if cid:
                        assert 1<=cid<=count;members[cid-1].update(map(int,mapped))
        assert face_seen.all() and all(len(values)==4 for values in members)
    by_key={tuple(sorted(map(int,quad))):i for i,quad in enumerate(quads)}
    assert len(by_key)==count
    mapping=np.array([by_key[tuple(sorted(values))] for values in members])
    assert len(np.unique(mapping))==count
    return mapping


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();run=args.run.resolve();out=args.out.resolve()
    assert run.is_dir() and out.drive.upper()=='E:' and not out.exists()
    out.mkdir(parents=True);checks={};evidence={};errors={};before={}
    def remember(path):
        path=Path(path).resolve();assert path.is_file(),str(path)
        before.setdefault(path,sha(path));return path
    def load(name):return json.loads(remember(run/name).read_text(encoding='utf-8-sig'))
    def attempt(name,operation):
        try:
            value=operation();checks[name]=True;return value
        except Exception as error:
            checks[name]=False;errors[name]=type(error).__name__+': '+str(error);return None
    status=attempt('execution_status_read',lambda:load('status.json')) or {}
    summary=attempt('summary_read',lambda:load('summary.json')) or {}
    settings=attempt('settings_read',lambda:load('final-settings.json')) or {}
    quality=attempt('quality_read',lambda:load('geometry-quality.json')) or {}
    targets=attempt('targets_read',lambda:load('declared-targets.json')) or {}
    source=attempt('provenance_read',lambda:load('code-provenance.json')) or {}
    gas=attempt('gas_readback_read',lambda:load('gas-constant-readback.json')) or {}
    auth=attempt('root_authorization_read',lambda:load('authorization.json')) or {}
    mesh=Path(source.get('mesh',str(run.parent/'mesh/nozzle.msh'))).resolve()
    mesh_dir=mesh.parent
    def provenance():
        assert source['source_sha256']==sha(remember(run/'executed-script.py'))==auth['source_sha256']
        assert source['mesh_sha256']==sha(remember(mesh))==auth['mesh_sha256']
        assert source['mesh_setup_sha256']==sha(remember(mesh_dir/'setup.json'))
        assert source['geometry_npz_sha256']==sha(remember(mesh_dir/'geometry.npz'))
        assert auth['approved'] is True and auth['root']=='root' and auth['allowed_mode']=='Run'
        assert auth['action']=='run_nozzle_fluent_v001' and Path(auth['out']).resolve()==run
        assert source['snapshot_before_solver_launch'] is True
        assert source['available_ram_before_launch_bytes']>=4*1024**3 and source['free_disk_before_launch_bytes']>=10*1024**3
        return dict(source_sha256=source['source_sha256'],mesh_sha256=source['mesh_sha256'],geometry_npz_sha256=source['geometry_npz_sha256'])
    evidence['source_binding']=attempt('source_binding',provenance)
    checks['fixed_run_completed']=status.get('stage') in ('completed','completed_with_export_errors') and summary.get('native_final_saved') is True and summary.get('actual_iterations')==2000
    checks['native_mesh_positive_and_quality']=quality.get('checked_before_flow') is True and quality.get('minimum_cell_volume_m3',0)>0 and quality.get('minimum_orthogonal_quality',-1)>=.15
    evidence['native_quality']=quality
    def actual_settings():
        assert settings['general']['solver']['type']=='density-based-implicit'
        assert settings['general']['solver']['two_dim_space']=='planar' and settings['general']['solver']['time']=='steady'
        assert settings['viscous']['model']=='inviscid' and settings['energy']['enabled'] is True
        assert near(settings['general']['operating_conditions']['operating_pressure'],0.)
        assert near(settings['references']['depth'],1.)
        assert settings['material']['density']['option']=='ideal-gas'
        assert near(settings['material']['specific_heat']['value'],1004.675)
        mw=float(settings['material']['molecular_weight']['value']);assert near(8314.46261815324/mw,287.05)
        assert near(gas['derived_gamma'],1.4) and near(gas['derived_R_J_kg_K'],287.05)
        assert near(settings['inlet']['momentum']['gauge_total_pressure']['value'],200000.)
        assert near(settings['inlet']['thermal']['total_temperature']['value'],300.)
        assert near(settings['outlet']['momentum']['gauge_pressure']['value'],140000.)
        assert near(settings['outlet']['thermal']['backflow_total_temperature']['value'],300.)
        assert settings['methods']['p_v_coupling']['flux_type']=='Roe-FDS'
        assert settings['methods']['spatial_discretization']['discretization_scheme']['flow']=='second-order-upwind'
        assert near(settings['controls']['courant_number'],3.)
        assert all(settings['residual']['equations'][name]['check_convergence'] is False for name in EQUATIONS)
        return gas
    evidence['actual_settings']=attempt('actual_Euler_gas_boundary_settings',actual_settings)
    residuals=attempt('native_residual_history',lambda:residual_rows(remember(run/'nozzle-final.dat.h5')))
    if residuals:
        write_csv(out/'fluent-residuals.csv',residuals)
        checks['actual_iteration_count_and_coverage']=residuals[-1]['iteration']==2000 and [row['iteration'] for row in residuals]==list(range(1,2001))
        evidence['final_scaled_residuals']={name:residuals[-1][name] for name in EQUATIONS}
        checks['final_residual_targets']=all(residuals[-1][name]<=targets['residual_targets'][name] for name in EQUATIONS)
        checks['driver_residual_readback_matches']=all(near(residuals[-1][name],summary['final_scaled_residuals'][name],1e-12) for name in EQUATIONS)
    def checkpoints():
        entries=load('checkpoints.json');assert [row['iteration'] for row in entries]==[300,600,900,1200,1500,1800,2000]
        for row in entries:
            remember(run/row['case']);data=remember(run/row['data']);native=residual_rows(data)
            assert native[-1]['iteration']==row['iteration']
            assert all(near(native[-1][name],row['scaled_residuals'][name],1e-12) for name in EQUATIONS)
        first=load('initialized-settings.json');second=load('second-order-settings.json')
        assert first['methods']['spatial_discretization']['discretization_scheme']['flow']=='first-order-upwind'
        assert near(first['controls']['courant_number'],1.)
        assert second['methods']['spatial_discretization']['discretization_scheme']['flow']=='second-order-upwind'
        assert near(second['controls']['courant_number'],3.)
        return entries
    evidence['checkpoints']=attempt('bounded_schedule_and_native_checkpoints',checkpoints)
    masses=attempt('native_per_iteration_massflows',lambda:mass_rows(remember(run/'massflow-history.out').read_text(encoding='utf-8',errors='replace')))
    if masses:
        write_csv(out/'fluent-massflows.csv',masses)
        checks['massflow_iteration_coverage']=[row['iteration'] for row in masses]==list(range(1,2001))
        final=masses[-100:];stability={}
        for name in ('mdot_in_kg_s_per_m_depth','mdot_out_kg_s_per_m_depth'):
            values=np.array([row[name] for row in final]);relative=float(np.ptp(values)/max(abs(np.mean(values)),1e-9))
            stability[name]=dict(mean=float(np.mean(values)),relative_range=relative)
        evidence['final100_massflow_stability']=stability
        checks['massflow_stability_target']=len(final)==100 and all(item['relative_range']<=.001 for item in stability.values())
    balance=attempt('native_mass_report_read',lambda:native_mass(remember(run/'mass-balance.txt').read_text(encoding='utf-8',errors='replace')))
    if balance:
        evidence['mass_balance']=balance;checks['mass_balance_target']=balance['relative_imbalance']<=1e-4
        if masses:
            checks['native_report_and_monitor_agree']=near(balance['inlet'],masses[-1]['mdot_in_kg_s_per_m_depth'],1e-4) and near(balance['outlet'],masses[-1]['mdot_out_kg_s_per_m_depth'],1e-4)
    def fields():
        geometry=dict(np.load(remember(mesh_dir/'geometry.npz')));count=len(geometry['quads']);assert count==9600
        mapping=common_mapping(remember(run/'nozzle-final.cas.h5'),geometry)
        arrays={}
        with h5py.File(run/'nozzle-final.dat.h5','r') as handle:
            root=handle['results/1/phase-1/cells']
            for name in ('SV_U','SV_V','SV_P','SV_DENSITY','SV_T'):
                values=chunk_array(root[name],count);ordered=np.empty(count);ordered[mapping]=values;arrays[name]=ordered
        p=arrays['SV_P'];rho=arrays['SV_DENSITY'];T=arrays['SV_T']
        assert np.all(p>0) and np.all(rho>0) and np.all(T>0)
        gas_error=np.abs(rho-p/(287.05*T))/np.maximum(rho,1e-12)
        assert np.max(gas_error)<1e-5
        speed=np.hypot(arrays['SV_U'],arrays['SV_V']);mach=speed/np.sqrt(1.4*287.05*T)
        assert np.isfinite(mach).all()
        np.savez_compressed(out/'native-cell-fields.npz',**arrays,Mach=mach,centres=geometry['centres'],
            native_cell_to_common_cell=mapping,common_geometry_sha256=np.array(sha(mesh_dir/'geometry.npz')))
        centre=[];ix=geometry['ix'];iy=geometry['iy']
        for column in range(240):
            selected=(ix==column)&np.isin(iy,[19,20]);assert np.count_nonzero(selected)==2
            centre.append(dict(x_m=float(np.mean(geometry['centres'][selected,0])),
                pressure_Pa=float(np.mean(p[selected])),temperature_K=float(np.mean(T[selected])),
                density_kg_m3=float(np.mean(rho[selected])),Mach=float(np.mean(mach[selected])),
                Ux_m_s=float(np.mean(arrays['SV_U'][selected])),Uy_m_s=float(np.mean(arrays['SV_V'][selected]))))
        write_csv(out/'fluent-centreline.csv',centre)
        x=np.array([row['x_m'] for row in centre]);pc=np.array([row['pressure_Pa'] for row in centre]);mc=np.array([row['Mach'] for row in centre])
        gradient=np.gradient(pc,x);eligible=np.flatnonzero((x>.42)&(x<.98));index=int(eligible[np.argmax(gradient[eligible])])
        left=max(index-4,0);right=min(index+4,239)
        shock=dict(method='Steepest positive centreline dp/dx within .42<x<.98m; arithmetic average of two cells straddling centreline',
            candidate_x_m=float(x[index]),candidate_dp_dx_Pa_m=float(gradient[index]),
            upstream_sample_x_m=float(x[left]),downstream_sample_x_m=float(x[right]),
            upstream_Mach=float(mc[left]),downstream_Mach=float(mc[right]),
            pressure_ratio_samples=float(pc[right]/pc[left]),
            resolved_compressive_supersonic_to_subsonic_transition=bool(mc[left]>1.05 and mc[right]<.95 and pc[right]/pc[left]>1.1),
            mesh_independence=False,benchmark_validation=False)
        return dict(cells=count,field_ranges={name:dict(minimum=float(np.min(value)),maximum=float(np.max(value))) for name,value in arrays.items()},
            Mach_range=[float(np.min(mach)),float(np.max(mach))],maximum_relative_ideal_gas_density_error=float(np.max(gas_error)),
            common_cell_mapping='All native face-derived quad node sets matched bijectively to shared geometry; no ordering assumption',
            native_Mach_derivation='sqrt(Ux^2+Uy^2)/sqrt(gamma*R*T), saved native static T and U; gamma1.4,R287.05',
            shock_diagnostic=shock)
    evidence['native_field_checks']=attempt('native_positive_fields_ideal_gas_and_common_mesh',fields)
    errors_export=summary.get('export_errors',[]);evidence['optional_export_errors']=errors_export
    checks['optional_exports_without_error']=not errors_export
    checks['all_input_hashes_unchanged']=all(sha(path)==digest for path,digest in before.items())
    shutil.copyfile(__file__,out/'executed-analysis.py')
    result=dict(model='Planar full converging-diverging nozzle; steady inviscid ideal-gas Euler',
        checks=checks,failed_checks=[name for name,value in checks.items() if not value],all_declared_checks_passed=all(checks.values()),
        evidence=evidence,errors=errors,actual_final_iteration=summary.get('actual_iterations'),
        physical_validation=False,mesh_independence=False,
        interpretation='A numerical calculation and cross-code illustration; convergence alone does not establish accuracy or physical validation',
        analysis_source_sha256=sha(__file__))
    write_json(out/'fluent-checks.json',result)
    manifest=dict(source_files=[dict(path=str(path),sha256=digest,after_sha256=sha(path),unchanged=sha(path)==digest) for path,digest in before.items()],
        analysis_source=dict(path=str(Path(__file__).resolve()),sha256=sha(__file__)),
        output_files=[dict(path=path.name,sha256=sha(path),bytes=path.stat().st_size) for path in sorted(out.iterdir()) if path.is_file()])
    write_json(out/'source-manifest.json',manifest)
    print(json.dumps(dict(output=str(out),checks=len(checks),all_passed=all(checks.values()),failed_checks=result['failed_checks'],errors=errors,
        final_scaled_residuals=evidence.get('final_scaled_residuals'),mass_balance=evidence.get('mass_balance'),native_field_checks=evidence.get('native_field_checks')),indent=2))


if __name__=='__main__':main()
