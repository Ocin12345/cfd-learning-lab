"""Read-only hash-bound audit of saved native Fluent Re100 cavity fields.

Create a fresh analysis directory; never changes a native solver file.
The mesh mapping helper reconstructs native cell IDs from saved face connectivity.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import h5py
import numpy as np
from analyse_nozzle_fluent_v005 import chunk_array, common_mapping

EQUATIONS=('continuity','x-velocity','y-velocity')
WALLS=('lid','bottom','left','right')
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):digest.update(block)
    return digest.hexdigest()


def write_json(path,value):
    with Path(path).open('x',encoding='utf-8',newline='\n') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def local_mass_sums(c0,c1,flux,cell_count):
    c0=np.asarray(c0);c1=np.asarray(c1);flux=np.asarray(flux,dtype=float)
    if not (c0.ndim==c1.ndim==flux.ndim==1 and len(c0)==len(c1)==len(flux)):
        raise ValueError('Native face index/flux shapes differ')
    if not (np.isfinite(flux).all() and np.array_equal(c0,c0.astype(int)) and np.array_equal(c1,c1.astype(int))):
        raise ValueError('Invalid native face mass flux or cell ID')
    c0=c0.astype(int);c1=c1.astype(int)
    if not (np.all((c0>=1)&(c0<=cell_count)) and np.all((c1>=0)&(c1<=cell_count))):
        raise ValueError('Native face cell index outside approved mesh')
    result=np.zeros(cell_count);np.add.at(result,c0-1,flux)
    inside=c1>0;np.add.at(result,c1[inside]-1,-flux[inside])
    return result


def flux_metrics(c0,c1,flux,volumes,rho=1.,U=1.,L=1.,depth=1.):
    volumes=np.asarray(volumes,dtype=float)
    if not (volumes.ndim==1 and len(volumes) and np.isfinite(volumes).all() and np.all(volumes>0)):
        raise ValueError('Positive finite native-mapped cell volumes required')
    if min(rho,U,L,depth)<=0:raise ValueError('Positive reference scales required')
    sums=local_mass_sums(c0,c1,flux,len(volumes));dimensionless=sums/(rho*volumes)*L/U
    boundary=np.asarray(c1)==0;net=float(np.asarray(flux)[boundary].sum())
    return dict(native_faces=int(len(flux)),native_cells=int(len(volumes)),
        max_cell_absolute_mass_sum_kg_s=float(abs(sums).max()),sum_cell_absolute_mass_sum_kg_s=float(abs(sums).sum()),
        max_dimensionless_divergence=float(abs(dimensionless).max()),
        rms_dimensionless_divergence=float(np.sqrt(np.sum(volumes*dimensionless**2)/volumes.sum())),
        boundary_signed_mass_flux_kg_s=net,boundary_absolute_face_mass_flux_sum_kg_s=float(abs(np.asarray(flux)[boundary]).sum()),
        boundary_net_characteristic_normalized=abs(net)/(rho*U*L*depth),
        incidence_sum_matches_boundary=bool(abs(float(sums.sum())-net)<=1e-12*max(1.,float(abs(flux).sum()))),
        normalization='Cell divergence: sum signed native SV_FLUX/(rho*cell area*depth)*L/U; volume-weighted RMS',
        boundary_limitation='Closed wallNet0 is imposed by impermeability and is not sole conservation proof')


def weighted_velocity_change(baseline,final,areas):
    baseline=np.asarray(baseline,dtype=float);final=np.asarray(final,dtype=float);areas=np.asarray(areas,dtype=float)
    if not (baseline.shape==final.shape==(len(areas),2) and np.isfinite(baseline).all() and np.isfinite(final).all() and np.all(areas>0)):
        raise ValueError('Invalid velocity stability input')
    denominator=float(np.sum(areas[:,None]*baseline**2))
    if denominator<=0:raise ValueError('Zero velocity norm cannot define relative stability')
    return float(np.sqrt(np.sum(areas[:,None]*(final-baseline)**2)/denominator))


def residual_history(path):
    rows={};common=None
    with h5py.File(path,'r') as handle:
        for name in EQUATIONS:
            group=handle['results/residuals/phase-1/'+name]
            iterations=np.asarray(group['iterations'][:]).reshape(-1);data=np.asarray(group['data'][:])
            assert np.array_equal(iterations,np.arange(1,2001)), 'Require complete native1..2000 residual history'
            assert data.ndim==2 and len(data)==len(iterations) and data.shape[1]>=2 and np.all(data[:,1]>0)
            assert np.isfinite(data[:,:2]).all()
            assert common is None or np.array_equal(iterations,common)
            common=iterations;scaled=data[:,0]/data[:,1];assert np.all(scaled>=0)
            for iteration,value in zip(iterations,scaled):rows.setdefault(int(iteration),dict(iteration=int(iteration)))[name]=float(value)
    return [rows[key] for key in sorted(rows)]


def native_faces(case,data):
    with h5py.File(case,'r') as handle,h5py.File(data,'r') as result:
        root=handle['meshes/1/faces']
        count=max(int(np.asarray(ds.attrs['maxId']).ravel()[0]) for ds in root['c0'].values())
        c0=chunk_array(root['c0'],count).astype(int);c1=np.zeros(count,int);seen=np.zeros(count,bool)
        for ds in root['c1'].values():
            lo=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;hi=int(np.asarray(ds.attrs['maxId']).ravel()[0])
            assert 0<=lo<hi<=count and not seen[lo:hi].any()
            c1[lo:hi]=np.asarray(ds[:]).reshape(-1).astype(int);seen[lo:hi]=True
        flux=chunk_array(result['results/1/phase-1/faces/SV_FLUX'],count)
        topology=root['zoneTopology'];names=b''.join(np.asarray(topology['name'][:]).ravel()).decode().split(';')
        minimum=np.asarray(topology['minId'][:]).reshape(-1).astype(int);maximum=np.asarray(topology['maxId'][:]).reshape(-1).astype(int)
        assert len(names)==len(minimum)==len(maximum)
        zones={name:dict(first=int(lo)-1,last=int(hi)) for name,lo,hi in zip(names,minimum,maximum)}
    return c0,c1,flux,zones


def same_mesh(first,last):
    with h5py.File(first,'r') as a,h5py.File(last,'r') as b:
        names=[];a['meshes/1'].visititems(lambda name,item:names.append(name) if isinstance(item,h5py.Dataset) else None)
        for name in names:assert np.array_equal(a['meshes/1/'+name][:],b['meshes/1/'+name][:]), 'Saved native mesh changed: '+name
    return len(names)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();run=args.run.resolve();out=args.out.resolve()
    assert run.is_dir() and out.drive.upper()=='E:' and not out.exists(), 'Fresh E: analysis directory required'
    out.mkdir(parents=True);before={};checks={};evidence={};errors={}
    def remember(path):
        path=Path(path).resolve();assert path.is_file(),str(path)
        before.setdefault(path,sha256(path));return path
    def load(name):return json.loads(remember(run/name).read_text(encoding='utf-8-sig'))
    def attempt(name,operation):
        try:value=operation();checks[name]=True;return value
        except Exception as error:
            checks[name]=False;errors[name]=type(error).__name__+': '+str(error);return None
    status=attempt('execution_status_read',lambda:load('status.json')) or {}
    summary=attempt('summary_read',lambda:load('summary.json')) or {}
    source=attempt('provenance_read',lambda:load('code-provenance.json')) or {}
    authorization=attempt('root_authorization_read',lambda:load('authorization.json')) or {}
    targets=attempt('targets_read',lambda:load('declared-targets.json')) or {}
    quality=attempt('native_quality_read',lambda:load('geometry-quality.json')) or {}
    settings=attempt('native_settings_read',lambda:load('final-settings.json')) or {}
    physical=attempt('native_physical_readback_read',lambda:load('native-physical-readback.json')) or {}
    mesh=Path(source.get('mesh',str(run.parent/'mesh/cavity.msh'))).resolve();geometry_path=mesh.parent/'geometry.npz';setup_path=mesh.parent/'setup.json'
    def provenance():
        assert source['source_sha256']==sha256(remember(run/'executed-script.py'))==authorization['source_sha256']
        assert source['mesh_sha256']==sha256(remember(mesh))==authorization['mesh_sha256']
        assert source['mesh_setup_sha256']==sha256(remember(setup_path))
        assert source['geometry_npz_sha256']==sha256(remember(geometry_path))
        assert authorization['action']=='run_cavity_fluent_v001' and authorization['approved'] is True and authorization['root']=='root' and authorization['allowed_mode']=='Run'
        assert Path(authorization['out']).resolve()==run
        return dict(source_sha256=source['source_sha256'],mesh_sha256=source['mesh_sha256'],geometry_npz_sha256=source['geometry_npz_sha256'])
    evidence['provenance']=attempt('native_inputs_bound_to_reviewed_authorization',provenance)
    checks['native_completed_2000']=status.get('stage') in ('completed','completed_with_export_errors') and summary.get('actual_iterations')==2000 and summary.get('native_final_saved') is True
    checks['native_quality_positive_volume_and_total1']=quality.get('minimum_orthogonal_quality',-1)>=.15 and quality.get('minimum_cell_volume_m3',0)>0 and abs(quality.get('total_cell_volume_m3',0)-1)<1e-8
    def settings_check():
        assert physical['density']==1 and physical['viscosity']==.01 and physical['viscous_model']=='laminar' and physical['energy_enabled'] is False
        assert physical['operating_pressure_Pa']==0 and physical['reference_depth_m']==1
        assert settings['general']['solver']['type']=='pressure-based' and settings['general']['solver']['time']=='steady'
        assert settings['viscous']['model']=='laminar' and settings['energy']['enabled'] is False
        assert settings['general']['operating_conditions']['operating_pressure']==0 and settings['references']['depth']==1
        assert settings['material']['density']['value']==1 and settings['material']['viscosity']['value']==.01
        assert settings['methods']['p_v_coupling']['flow_scheme']=='SIMPLE'
        spatial=settings['methods']['spatial_discretization'];assert spatial['gradient_scheme']=='least-square-cell-based'
        assert spatial['discretization_scheme']['mom']=='second-order-upwind' and spatial['discretization_scheme']['pressure']=='second-order'
        lid=settings['walls']['lid']['momentum']
        assert lid['wall_motion']=='Moving Wall' and lid['relative']=='Absolute' and lid['velocity_spec']=='Components'
        assert lid['shear_condition']=='No Slip' and [item['value'] for item in lid['velocity_components']]==[1,0]
        for name in WALLS[1:]:
            wall=settings['walls'][name]['momentum'];assert wall['wall_motion']=='Stationary Wall' and wall['shear_condition']=='No Slip'
        for name in EQUATIONS:assert settings['residual']['equations'][name]['check_convergence'] is False
        return physical
    evidence['settings']=attempt('native_final_physical_wall_and_numerical_settings',settings_check)
    final_case=attempt('native_final_case_read',lambda:remember(run/'cavity-final.cas.h5'))
    final_data=attempt('native_final_data_read',lambda:remember(run/'cavity-final.dat.h5'))
    baseline_case=attempt('native1800_case_read',lambda:remember(run/'cavity-iter1800.cas.h5'))
    baseline_data=attempt('native1800_data_read',lambda:remember(run/'cavity-iter1800.dat.h5'))
    for endpoint in (500,1000,1500,1800,2000):
        attempt('checkpoint'+str(endpoint)+'_case_read',lambda endpoint=endpoint:remember(run/f'cavity-iter{endpoint:04d}.cas.h5'))
        attempt('checkpoint'+str(endpoint)+'_data_read',lambda endpoint=endpoint:remember(run/f'cavity-iter{endpoint:04d}.dat.h5'))
    rows=attempt('native_residual_history_complete',lambda:residual_history(final_data))
    if rows:
        evidence['final_scaled_residuals']={name:rows[-1][name] for name in EQUATIONS}
        checks['final_residuals_le_1e_6']=all(rows[-1][name]<=1e-6 for name in EQUATIONS)
        evidence['final200_max_scaled_residuals']={name:max(row[name] for row in rows[-200:]) for name in EQUATIONS}
        with (out/'fluent-residuals.csv').open('x',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=['iteration',*EQUATIONS]);writer.writeheader();writer.writerows(rows)
    geometry=attempt('source_geometry_read',lambda:dict(np.load(remember(geometry_path),allow_pickle=False)))
    fields=None
    if geometry is not None and final_case and final_data and baseline_case and baseline_data:
        def mapped_fields():
            assert geometry['quads'].shape==(6400,4) and geometry['nodes'].shape==(6561,2)
            areas=geometry['areas'].reshape(-1);assert np.isfinite(areas).all() and np.all(areas>0) and abs(areas.sum()-1)<1e-10
            mapping=common_mapping(final_case,geometry);assert mapping.shape==(6400,)
            native={};common={}
            with h5py.File(final_data,'r') as handle:
                group=handle['results/1/phase-1/cells']
                for name,key in (('u','SV_U'),('v','SV_V'),('p_raw','SV_P'),('rho','SV_DENSITY'),('mu','SV_MU_LAM')):
                    native[name]=chunk_array(group[key],6400);common[name]=np.empty(6400);common[name][mapping]=native[name]
            assert np.allclose(common['rho'],1,atol=1e-12,rtol=0) and np.allclose(common['mu'],.01,atol=1e-12,rtol=0)
            native_baseline={}
            with h5py.File(baseline_data,'r') as handle:
                for name,key in (('u','SV_U'),('v','SV_V')):native_baseline[name]=chunk_array(handle['results/1/phase-1/cells/'+key],6400)
            baseline=np.empty((6400,2));baseline[mapping]=np.column_stack([native_baseline['u'],native_baseline['v']])
            final=np.column_stack([common['u'],common['v']]);change=weighted_velocity_change(baseline,final,areas)
            common['p']=common['p_raw']-np.sum(areas*common['p_raw'])/areas.sum();common['speed']=np.linalg.norm(final,axis=1)
            return dict(mapping=mapping,areas=areas,common=common,velocity_change=change,native_area=areas[mapping])
        fields=attempt('native_topology_mapping_and_finite_fields',mapped_fields)
        evidence['unchanged_native_mesh_datasets']=attempt('native1800_and2000_mesh_identical',lambda:same_mesh(baseline_case,final_case))
        if fields is not None:
            evidence['velocity_stability']=dict(endpoints=[1800,2000],weighted_velocity_L2_relative_change=fields['velocity_change'],
                definition='sqrt(sum(source cell area*|U2000-U1800|²)/sum(source cell area*|U1800|²))')
            checks['final200_weighted_velocity_L2_relative_change_le_1e_3']=fields['velocity_change']<=.001
            def local_conservation():
                c0,c1,flux,zones=native_faces(final_case,final_data)
                metrics=flux_metrics(c0,c1,flux,fields['native_area'])
                metrics['zones']={}
                for name in WALLS:
                    zone=zones[name];selection=slice(zone['first'],zone['last']);assert np.all(c1[selection]==0)
                    metrics['zones'][name]=dict(faces=zone['last']-zone['first'],native_signed_SV_FLUX_kg_s=float(flux[selection].sum()),
                        absolute_face_flux_sum_kg_s=float(abs(flux[selection]).sum()))
                sums=local_mass_sums(c0,c1,flux,6400);common_sums=np.empty(6400);common_sums[fields['mapping']]=sums
                fields['common']['cell_mass_sum_kg_s']=common_sums
                fields['common']['dimensionless_divergence']=common_sums/fields['areas']
                return metrics
            conservation=attempt('full_native_face_flux_and_incidence_coverage',local_conservation)
            evidence['local_conservation']=conservation
            if conservation:
                checks['native_local_divergence_max_le_1e_5']=conservation['max_dimensionless_divergence']<=1e-5
                checks['native_local_divergence_rms_le_1e_6']=conservation['rms_dimensionless_divergence']<=1e-6
                checks['closed_wall_net_characteristic_normalized_le_1e_10']=conservation['boundary_net_characteristic_normalized']<=1e-10
                checks['native_incidence_sum_equals_boundary_flux']=conservation['incidence_sum_matches_boundary']
            np.savez_compressed(out/'fluent-fields-common.npz',**geometry,**fields['common'])
    def wall_report():
        text=remember(run/'wall-mass-balance.txt').read_text(encoding='utf-8',errors='replace');result={}
        for name in (*WALLS,'Net'):
            matches=re.findall(r'^\s*'+re.escape(name)+r'\s+('+NUMBER+r')\s*$',text,re.M);assert matches,'Missing native wall report '+name
            result[name]=float(matches[-1]);assert math.isfinite(result[name])
        assert abs(result['Net'])<=1e-10
        result['relative_mass_imbalance']=None
        result['relative_mass_imbalance_reason']='No inlet or outlet; dividing wallNet0 by inlet0 would be undefined. Use local native flux divergence.'
        return result
    evidence['closed_wall_native_mass_report']=attempt('closed_wall_report_present_without_zero_over_zero',wall_report)
    def monitors():
        text=remember(run/'field-monitor-history.out').read_text(encoding='utf-8',errors='replace');records=[]
        for line in text.splitlines():
            parts=line.split()
            if len(parts)!=3:continue
            try:iteration,speed,pressure=map(float,parts)
            except ValueError:continue
            assert all(map(math.isfinite,(iteration,speed,pressure))) and iteration.is_integer()
            if iteration>0:records.append(dict(iteration=int(iteration),mean_speed=speed,mean_pressure=pressure))
        assert [row['iteration'] for row in records]==list(range(1,2001))
        with (out/'fluent-field-monitors.csv').open('x',newline='',encoding='utf-8') as stream:
            writer=csv.DictWriter(stream,fieldnames=['iteration','mean_speed','mean_pressure']);writer.writeheader();writer.writerows(records)
        return dict(actual_rows=len(records),last=records[-1])
    evidence['monitors']=attempt('native_per_iteration_monitors_complete',monitors)
    def licence_and_logs():
        startup=list(run.glob('fluent-*.trn'));assert len(startup)==1
        text=remember(startup[0]).read_text(encoding='utf-8',errors='replace')
        transcript=remember(run/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace')
        assert 'student version' in text.lower()
        processes=len(re.findall(r'^n\d+\*?\s+',text,re.M));assert processes==2
        assert re.search(r'^\s*2000\s+',transcript,re.M)
        assert not re.search(r'Divergence detected|floating point exception|fatal error',transcript,re.I)
        return dict(native_compute_processes=processes,student_version_verified=True)
    evidence['live_native_licence_and_log']=attempt('native_student_two_process_licence_and_no_fatal',licence_and_logs)
    def plots():
        os.environ.setdefault('MPLCONFIGDIR','E:/CFD-Learning-Lab/temp/matplotlib-cavity')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.collections import PolyCollection
        for name,key,title,cmap in (('velocity','speed','Fluent cavity speed |U| (m/s)','viridis'),('pressure','p','Fluent cavity pressure minus area mean (Pa)','coolwarm')):
            fig,ax=plt.subplots(figsize=(7,6),layout='constrained')
            image=PolyCollection(geometry['nodes'][geometry['quads']],array=fields['common'][key],cmap=cmap,edgecolors='none')
            ax.add_collection(image);ax.set(xlim=(0,1),ylim=(0,1),xlabel='x (m)',ylabel='y (m)',title=title);ax.set_aspect('equal')
            fig.colorbar(image,ax=ax);fig.savefig(out/('fluent-'+name+'.png'),dpi=160);plt.close(fig)
        if rows:
            fig,ax=plt.subplots(figsize=(8,4),layout='constrained')
            for name in EQUATIONS:ax.semilogy([row['iteration'] for row in rows],[row[name] for row in rows],label=name)
            ax.axhline(1e-6,color='black',linestyle='--',linewidth=1);ax.set(xlabel='Native steady iteration',ylabel='Scaled residual',title='Fluent Re100 cavity residuals');ax.legend();ax.grid(alpha=.2)
            fig.savefig(out/'fluent-residuals.png',dpi=160);plt.close(fig)
        return ['fluent-velocity.png','fluent-pressure.png','fluent-residuals.png']
    if fields is not None:evidence['derived_saved_native_images']=attempt('saved_native_visualizations_created',plots)
    checks['all_hashed_native_inputs_unchanged']=all(sha256(path)==digest for path,digest in before.items())
    shutil.copyfile(__file__,out/'executed-analyser.py')
    helper=Path(__file__).with_name('analyse_nozzle_fluent_v005.py');shutil.copyfile(helper,out/'mapping-helper-source.py')
    report=dict(kind='Fluent Re100 lid-driven cavity saved-native audit',checks=checks,evidence=evidence,errors=errors,
        all_checks_passed=all(checks.values()),failed_checks=[key for key,value in checks.items() if not value],
        physical_validation=False,mesh_independence=False,
        boundary_note='Closed boundaryNet0 alone is not independent conservation; local native face flux divergence is checked.',
        pressure_note='Pressure comparisons and derived visuals subtract source area weighted mean.',
        analyser_sha256=sha256(__file__),mapping_helper_sha256=sha256(helper))
    write_json(out/'fluent-checks.json',report)
    write_json(out/'fluent-source-manifest.json',dict(inputs=[dict(path=str(path),sha256=digest) for path,digest in before.items()],
        all_inputs_unchanged=checks['all_hashed_native_inputs_unchanged']))
    print(json.dumps(dict(all_checks_passed=report['all_checks_passed'],failed_checks=report['failed_checks'],evidence=evidence,errors=errors),allow_nan=False))
    return 0 if report['all_checks_passed'] else 1

if __name__=='__main__':raise SystemExit(main())
