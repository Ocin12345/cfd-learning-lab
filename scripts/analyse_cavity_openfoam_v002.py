"""Read-only native cavity checks; writes only a fresh external E: analysis folder."""
from pathlib import Path
import argparse
import csv
import json
import re
import sys
import numpy as np
from prepare_openfoam_cavity_v001 import (
    PATCHES, EMPTY, FLOAT, require, sha, clean, dictionaries, native_poly,
    native_patches, native_geometry, mesh_quality, json_new, copied, checked_case,
)
class Evidence:
    def __init__(self):
        self.before={}
    def touch(self,path):
        path=Path(path).resolve()
        require(path.is_file(),'Missing actual native evidence '+str(path))
        if path not in self.before:
            self.before[path]=sha(path)
        return path
    def text(self,path):
        return self.touch(path).read_text(encoding='utf-8-sig')
    def read_json(self,path):
        return json.loads(self.text(path))
    def finish(self):
        rows=[dict(path=str(path),bytes=path.stat().st_size,sha256_before=digest,sha256_after=sha(path)) for path,digest in self.before.items()]
        require(all(r['sha256_before']==r['sha256_after'] for r in rows),'Native input changed during read-only analysis')
        return dict(inputs_unchanged_during_analysis=True,sources=rows,solver_started_or_modified=False)
def nonfinite_check(text):
    require(not re.search(r'(?i)(?<![A-Za-z0-9_])[-+]?(?:nan|inf(?:inity)?)(?![A-Za-z0-9_])',text),'Nonfinite native field token')
def parse_value(text,count,components,label='internalField'):
    kind='vector' if components==3 else 'scalar'
    m=re.search(r'\b'+label+r'\s+uniform\s+([^;]+);',text)
    if m:
        row=np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ')
        require(row.size==components and np.isfinite(row).all(),'Invalid uniform native '+label)
        return np.tile(row,(count,1)) if components>1 else np.full(count,row[0])
    m=re.search(r'\b'+label+r'\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    require(m is not None and int(m[1])==count,'Native '+label+' count differs')
    values=np.fromstring(m[2].replace('(',' ').replace(')',' '),sep=' ')
    require(values.size==count*components and np.isfinite(values).all(),'Native '+label+' length/finite gate')
    return values.reshape((count,components)) if components>1 else values
def field(raw,count,components,dimensions):
    text=clean(raw); nonfinite_check(text)
    require(re.search(r'\bformat\s+ascii\s*;',text),'Native ASCII field required')
    d=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',text)
    require(d is not None and list(map(int,d[1].split()))==dimensions,'Native dimensions differ')
    return parse_value(text,count,components)
def boundary_blocks(raw):
    text=clean(raw); nonfinite_check(text)
    require('boundaryField' in text,'Missing native boundaryField')
    rows=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',text.split('boundaryField',1)[1]))
    require(set(rows)==set(PATCHES)|{EMPTY},'Actual saved boundary names differ')
    return rows
def flux_and_divergence(raw,patches,geometry):
    n_internal=len(geometry['neighbour'])
    internal=field(raw,n_internal,1,[0,3,-1,0,0,0,0])
    rows=boundary_blocks(raw)
    phi=np.zeros(len(geometry['faces']))
    phi[:n_internal]=internal
    totals={}; absolutes={}
    for name in PATCHES:
        patch=patches[name]; values=parse_value(rows[name],patch['count'],1,'value')
        phi[patch['start']:patch['start']+patch['count']]=values
        totals[name]=float(values.sum()); absolutes[name]=float(np.abs(values).sum())
    require(re.search(r'\btype\s+empty\s*;',rows[EMPTY]),'Native phi planes must be empty')
    owner=np.asarray(geometry['owner']); neighbour=np.asarray(geometry['neighbour'])
    cell_sums=np.bincount(owner,weights=phi,minlength=6400)-np.bincount(neighbour,weights=internal,minlength=6400)
    volumes=np.asarray(geometry['volumes'])
    local=cell_sums/volumes
    # L/U=1s for this dimensionless scaling; phi is volume flux through0.01m depth.
    metrics=dict(physical_boundary_signed_outward_volume_flux_m3_s=totals,
        physical_boundary_absolute_volume_flux_m3_s=absolutes,net_outward_volume_flux_m3_s=sum(totals.values()),
        closed_wall_net_over_U_L_depth=abs(sum(totals.values()))/.01,
        closed_wall_all_face_absolute_over_U_L_depth=sum(absolutes.values())/.01,
        native_internal_face_count=n_internal,native_physical_boundary_face_count=320,native_empty_face_count=12800,
        actual_native_cell_count=len(volumes),sum_cell_signed_phi_m3_s=float(cell_sums.sum()),
        cell_sum_minus_boundary_net_m3_s=float(cell_sums.sum()-sum(totals.values())),
        local_divergence_scaled_max_abs=float(np.max(np.abs(local))),
        local_divergence_scaled_RMS=float(np.sqrt(np.mean(local**2))),
        local_divergence_scaled_volume_weighted_RMS=float(np.sqrt(np.average(local**2,weights=volumes))),
        local_divergence_definition='sum signed outward native phi on every cell face / oriented native cell volume * L/U; L/U=1s',
        global_boundary_caveat='Closed impermeable walls have zero expected net flux by boundary condition; this alone does not establish interior continuity. No inlet normalization exists.',
        owner_neighbour_orientation='phi positive outward from owner; sign reversed for neighbour')
    require(np.isfinite(local).all(),'Nonfinite actual local continuity')
    return metrics,local
def residual_history(text):
    require(not re.search(r'FOAM FATAL|Segmentation fault|Floating point exception\s*$',text,re.M),'Actual solver fatal error')
    require(re.search(r'^End\s*$',text,re.M),'Actual solver log missing End')
    sections=re.split(r'^Time = ('+FLOAT+r')s?\s*$',text,flags=re.M); rows=[]
    for i in range(1,len(sections),2):
        iteration=float(sections[i]); samples={name:[] for name in ('p','Ux','Uy','Uz')}
        for name,value in re.findall(r'Solving for (p|Ux|Uy|Uz), Initial residual = ('+FLOAT+r')',sections[i+1]):
            samples[name].append(float(value))
        require(all(samples[name] for name in ('p','Ux','Uy')),'Missing actual p/Ux/Uy residual samples')
        row=dict(iteration=int(iteration),**{name:max(samples[name]) for name in ('p','Ux','Uy')})
        if samples['Uz']:
            row['Uz']=max(samples['Uz'])
        require(iteration==row['iteration'] and all(np.isfinite(v) for v in row.values()),'Invalid printed residual')
        rows.append(row)
    require([row['iteration'] for row in rows]==list(range(1,2001)),'Require fixed-cap actual iterations1..2000')
    require(len({tuple(row) for row in rows})==1,'Printed residual variables change across iterations')
    return rows
def cell_mapping(centres,geometry_file,volumes):
    require(centres.shape==(6400,3),'Native C count differs')
    require(np.max(np.abs(centres[:,2]))<=1e-10,'Actual C must be centred in extrusion')
    ix=np.rint(centres[:,0]*80-.5).astype(int); iy=np.rint(centres[:,1]*80-.5).astype(int)
    require(np.all((ix>=0)&(ix<80)&(iy>=0)&(iy<80)),'Actual C outside centre grid')
    expected=np.column_stack(((ix+.5)/80,(iy+.5)/80))
    require(np.max(np.abs(centres[:,:2]-expected))<=1e-10,'Actual C differs from exact80x80 grid')
    keys=iy*80+ix
    require(np.array_equal(np.sort(keys),np.arange(6400)),'Native-to-grid bijection failed')
    with np.load(geometry_file,allow_pickle=False) as source:
        require(all(name in source for name in ('centres','areas','grid_shape','ix','iy')),'Shared source geometry mapping arrays missing')
        source_c=np.asarray(source['centres']).reshape(6400,2)
        source_a=np.asarray(source['areas']).reshape(6400)
        sx=np.asarray(source['ix']).reshape(6400).astype(int); sy=np.asarray(source['iy']).reshape(6400).astype(int)
        require(np.array_equal(np.asarray(source['grid_shape']).reshape(-1),[80,80]),'Shared grid shape differs')
    skeys=sy*80+sx
    require(np.array_equal(np.sort(skeys),np.arange(6400)),'Source-to-grid bijection failed')
    source_by_key=np.argsort(skeys); native_to_source=source_by_key[keys]
    require(np.max(np.abs(centres[:,:2]-source_c[native_to_source]))<=1e-10,'Actual C and shared source centres differ')
    require(np.isfinite(source_a).all() and np.all(source_a>0) and np.max(np.abs(volumes/.01-source_a[native_to_source]))<=1e-12,'Actual native/source area mapping differs')
    native_index_grid=np.argsort(keys).reshape(80,80)
    return dict(native_cell_to_source_cell=native_to_source,native_cell_index_grid=native_index_grid,
        source_cell_index_grid=source_by_key.reshape(80,80),ix=ix,iy=iy),dict(native_cells=6400,source_cells=6400,unique_grid_keys=6400,
        centre_max_mapping_error_m=float(np.max(np.abs(centres[:,:2]-source_c[native_to_source]))),
        native_source_area_max_error_m2=float(np.max(np.abs(volumes/.01-source_a[native_to_source]))),
        native_to_source_bijective=True,native_to_grid_bijective=True,ordering_assumed=False)
def check_wall_fields(u_raw,p_raw,patches):
    ub=boundary_blocks(u_raw);pb=boundary_blocks(p_raw)
    for name in PATCHES:
        expected=np.array([1,0,0]) if name=='lid' else np.zeros(3)
        require(re.search(r'\btype\s+'+('fixedValue' if name=='lid' else 'noSlip')+r'\s*;',ub[name]),'Actual native wall velocity type differs')
        # Foundation14 noSlip writes its type without a redundant value entry.
        # Exact noSlip type on a stationary native mesh implies zero wall U.
        # Validate every explicit value if one is actually present.
        if name=='lid' or re.search(r'\bvalue\s+',ub[name]):
            require(np.max(np.abs(parse_value(ub[name],patches[name]['count'],3,'value')-expected))<=1e-12,'Actual native saved wall velocity differs')
        require(re.search(r'\btype\s+zeroGradient\s*;',pb[name]),'Actual native saved wall pressure type differs')
    require(re.search(r'\btype\s+empty\s*;',ub[EMPTY]) and re.search(r'\btype\s+empty\s*;',pb[EMPTY]),'Actual native U/p empty planes differ')
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();case=checked_case(args.run);out=args.out.resolve()
    require(not out.exists() and out.drive.upper()=='E:','Require fresh E: analysis output')
    require(case not in out.parents and out!=case,'Read-only analysis must write outside native case')
    e=Evidence();checks={}
    for p in (Path(__file__),Path(__file__).with_name('prepare_openfoam_cavity_v001.py')):
        e.touch(p)
    declaration=e.read_json(case/'declared-comparison.json');status=e.read_json(case/'run-status.json')
    require(status.get('stage')=='finished_execution' and status.get('final_iteration')==2000,'Fixed-cap native execution incomplete')
    for name,expected in dictionaries().items():
        require(e.text(case/name)==expected,'Actual native dictionary differs '+name)
    checks['actual_native_dictionaries_match_exact_settings']=True
    prepared=e.read_json(case/'provenance/prepared-input-sha256.json')
    for name,digest in prepared.items():
        e.touch(case/name);require(sha(case/name)==digest,'Prepared source changed '+name)
    require(sha(e.touch(case/'input-cavity.msh'))==declaration['mesh_source_sha256'],'Shared source mesh hash differs')
    require(sha(e.touch(case/'provenance/geometry.npz'))==declaration['geometry_sha256'],'Shared source geometry hash differs')
    checks['prepared_sources_and_shared_mesh_hashes_match']=True
    preflight=e.read_json(case/'provenance/preflight-validation.json');require(preflight.get('passed') is True,'Native preflight not passed')
    poly=native_poly(case)
    for name,digest in preflight['native_mesh_sha256'].items():
        e.touch(poly/name);require(sha(poly/name)==digest,'Native static mesh changed')
    patches=native_patches(poly);geometry=native_geometry(poly);volumes=np.asarray(geometry['volumes'])
    quality=mesh_quality(e.text(case/'mesh-check.log'))
    checks['full_checkMesh_positive_actual_volumes_and_exact_native_counts_passed']=True
    centres=field(e.text(case/'2000/C'),6400,3,[0,1,0,0,0,0,0])
    mapping,mapping_record=cell_mapping(centres,case/'provenance/geometry.npz',volumes)
    checks['actual_C_source_and_grid_mappings_bijective']=True
    data={};stats={};continuity={}
    for iteration in range(100,2001,100):
        u_raw=e.text(case/str(iteration)/'U');p_raw=e.text(case/str(iteration)/'p');phi_raw=e.text(case/str(iteration)/'phi')
        u=field(u_raw,6400,3,[0,1,-1,0,0,0,0]);p=field(p_raw,6400,1,[0,2,-2,0,0,0,0])
        check_wall_fields(u_raw,p_raw,patches)
        metrics,local=flux_and_divergence(phi_raw,patches,geometry)
        p_mean=float(np.average(p,weights=volumes));p0=p-p_mean
        stats[str(iteration)]=dict(U_min_m_s=u.min(axis=0).tolist(),U_max_m_s=u.max(axis=0).tolist(),
            speed_min_m_s=float(np.linalg.norm(u,axis=1).min()),speed_max_m_s=float(np.linalg.norm(u,axis=1).max()),
            raw_kinematic_pressure_min=float(p.min()),raw_kinematic_pressure_max=float(p.max()),
            area_mean_pressure_Pa=p_mean,mean_subtracted_pressure_min_Pa=float(p0.min()),mean_subtracted_pressure_max_Pa=float(p0.max()),
            max_abs_out_of_plane_U_m_s=float(np.abs(u[:,2]).max()))
        continuity[str(iteration)]=metrics
        data[iteration]=dict(U=u,p=p,p_area_mean_subtracted_Pa=p0,local_divergence_scaled=local)
    checks['all20_saved_native_U_p_phi_finite_correct_counts_dimensions_and_wall_values']=True
    final=continuity['2000']
    checks['final_native_local_divergence_max_le_1e_5']=final['local_divergence_scaled_max_abs']<=1e-5
    checks['final_native_local_divergence_RMS_le_1e_6']=final['local_divergence_scaled_RMS']<=1e-6
    checks['closed_wall_net_over_U_L_depth_le_1e_10']=final['closed_wall_net_over_U_L_depth']<=1e-10
    checks['closed_wall_all_face_abs_over_U_L_depth_le_1e_10']=final['closed_wall_all_face_absolute_over_U_L_depth']<=1e-10
    checks['actual_native_2D_Uz_max_le_1e_12']=stats['2000']['max_abs_out_of_plane_U_m_s']<=1e-12
    residuals=residual_history(e.text(case/'solver.log'));terminal={name:residuals[-1][name] for name in ('p','Ux','Uy')}
    checks['terminal_printed_p_Ux_Uy_initial_residuals_le_1e_6']=max(terminal.values())<=1e-6
    delta=data[2000]['U']-data[1800]['U']
    numerator=float(np.sum(volumes[:,None]*delta**2));denominator=float(np.sum(volumes[:,None]*data[2000]['U']**2))
    require(denominator>0,'Missing actual cavity circulation')
    relative=float(np.sqrt(numerator/denominator))
    checks['U_1800_to_2000_volume_weighted_relative_L2_le_1e_3']=relative<=1e-3
    auth=e.read_json(case/'provenance/root-run-authorization.json')
    require(auth.get('approved') is True and auth.get('approved_by')=='root' and auth.get('action')=='run_cavity_openfoam_v001' and auth.get('allowed_mode')=='Run','Missing root run authorization')
    executed=e.touch(case/'provenance/run-openfoam-executed.sh');original=e.touch(case/'provenance/run-openfoam.sh')
    require(sha(executed)==sha(original)==auth.get('wrapper_sha256'),'Executed wrapper hash differs')
    require(sha(case/'provenance/prepare_openfoam_cavity_v001.py')==auth.get('preparation_sha256'),'Executed preparer hash differs')
    require(sha(case/'provenance/Launch-OpenFOAM.ps1')==auth.get('launcher_sha256'),'Reviewed launcher hash differs')
    checks['reviewed_execution_sources_match_root_authorization']=True
    installed=e.text(case/'provenance/installed-source-sha256.txt')
    installed_rows=[line.split(None,1) for line in installed.splitlines() if line.strip()]
    require(len(installed_rows)==6,'Expected six actual installed official sources')
    for digest,source_path in installed_rows:
        require(source_path.startswith('/opt/openfoam14/'),'Installed official source root differs')
        relative=source_path[len('/opt/openfoam14/'):]
        require(sha(e.touch(case/'provenance/installed-sources'/relative))==digest,'Copied installed source hash differs')
    checks['installed_official_sources_match_recorded_hashes']=True
    for path in sorted((case/'provenance/installed-sources').rglob('*')):
        if path.is_file():
            e.touch(path)
    source=e.finish();checks['inputs_unchanged_during_read_only_analysis']=True
    record=dict(solver='OpenFOAM Foundation14 incompressibleFluid',actual_settings=declaration,steady_iteration_indices_are_not_physical_time=True,
        mesh_cells=6400,mesh_quality=quality,actual_native_patches=patches,
        native_geometry=dict(points=13122,faces=25760,internal_faces=12640,cell_volume_min_m3=float(volumes.min()),cell_volume_max_m3=float(volumes.max()),total_volume_m3=float(volumes.sum()),bounds_m=geometry['bounds']),
        actual_cell_mapping=mapping_record,native_fields=stats,native_continuity_all_saved_states=continuity,final_native_continuity=final,
        residuals=dict(final_maximum_printed_initial_residual_per_variable=terminal,
            last100_maximum_printed_initial_residual_per_variable={name:max(row[name] for row in residuals[-100:]) for name in terminal},
            definition='Maximum printed initial linear-system residual per variable in a steady SIMPLE iteration',norm_equivalent_to_fluent=False),
        U_stability=dict(from_iteration=1800,to_iteration=2000,volume_weighted_relative_L2=relative,definition='sqrt(sum(V*|U2000-U1800|^2)/sum(V*|U2000|^2))'),
        final_iteration=2000,saved_iterations=list(range(100,2001,100)),selected_iterations=[500,1000,1500,1800,2000],
        native_wall_value_semantics='Exact native noSlip on three stationary walls fixes zeroU; explicit values checked if present; lid explicit fixedValue(1,0,0)',
        noSlip_official_source='https://cpp.openfoam.org/v14/classFoam_1_1noSlipFvPatchVectorField.html',
        pressure_comparison='rho*p minus its own area-weighted mean; rho=1kg/m3; pressure reference has arbitrary additive gauge',
        check_results=checks,failed_checks=[name for name,value in checks.items() if not value],all_declared_checks_passed=all(checks.values()),
        physical_validation=False,mesh_independence_established=False,benchmark_validation=False,
        interpretation='AI-assisted educational steady one-mesh solve; closed-wall net balance alone is uninformative; numerical target misses retained without extension')
    out.mkdir(parents=True)
    json_new(out/'checks.json',record);json_new(out/'source-manifest.json',source)
    with (out/'residual-history.csv').open('x',encoding='ascii',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(residuals[0]));writer.writeheader();writer.writerows(residuals)
    arrays=dict(centres=centres,volumes=volumes,areas=volumes/.01,grid_shape=np.array([80,80]),**mapping,**data[2000])
    for iteration in (500,1000,1500,1800,2000):
        for key,value in data[iteration].items():
            arrays[key+'_'+str(iteration)]=value
    np.savez_compressed(out/'native-fields.npz',**arrays)
    with (out/'actual-mapped-cell-fields.csv').open('x',encoding='ascii',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['native_cell','source_cell','ix','iy','x_m','y_m','z_m','volume_m3','Ux_m_s','Uy_m_s','Uz_m_s','kinematic_p','p_minus_area_mean_Pa','div_phi_over_V_times_L_over_U'])
        for i in range(6400):
            writer.writerow([i,int(mapping['native_cell_to_source_cell'][i]),int(mapping['ix'][i]),int(mapping['iy'][i]),*centres[i],volumes[i],*data[2000]['U'][i],data[2000]['p'][i],data[2000]['p_area_mean_subtracted_Pa'][i],data[2000]['local_divergence_scaled'][i]])
    for path in (Path(__file__),Path(__file__).with_name('prepare_openfoam_cavity_v001.py')):
        copied(path,out/('executed-'+path.name))
    print(json.dumps(dict(output=str(out),all_declared_checks_passed=record['all_declared_checks_passed'],failed_checks=record['failed_checks']),indent=2))
    return 0 if record['all_declared_checks_passed'] else 2
if __name__=='__main__':
    sys.exit(main())

