"""Read-only native square-obstacle wake checks; writes only a fresh external E: analysis folder."""
from pathlib import Path
import argparse
import csv
import json
import re
import sys
import numpy as np
from prepare_openfoam_square_wake_v003 import (
    PATCHES, EMPTY, FLOAT, CELLS, INTERNAL_FACES, TOTAL_FACES, POINTS, require, sha, clean, dictionaries, native_poly,
    native_patches, native_geometry, mesh_quality, json_new, copied, checked_case, INSTALLED_SOURCES,
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
    phi=np.zeros(len(geometry['faces']));phi[:n_internal]=internal
    totals={};absolutes={};inward={};outward={}
    for name in PATCHES:
        patch=patches[name]
        if name in ('top','bottom') and not re.search(r'\bvalue\s+',rows[name]):
            require(patch['kind']=='symmetry' and re.search(r'\btype\s+symmetry\s*;',rows[name]),'Zero flux requires actual symmetry constraint')
            values=np.zeros(patch['count'])
        else:
            values=parse_value(rows[name],patch['count'],1,'value')
        phi[patch['start']:patch['start']+patch['count']]=values
        totals[name]=float(values.sum());absolutes[name]=float(np.abs(values).sum())
        inward[name]=float(-values[values<0].sum());outward[name]=float(values[values>0].sum())
    require(re.search(r'\btype\s+empty\s*;',rows[EMPTY]),'Native phi planes must be empty')
    owner=np.asarray(geometry['owner']);neighbour=np.asarray(geometry['neighbour'])
    cell_sums=np.bincount(owner,weights=phi,minlength=CELLS)-np.bincount(neighbour,weights=internal,minlength=CELLS)
    volumes=np.asarray(geometry['volumes']);local=cell_sums/volumes*1
    total_inward=sum(inward.values());total_outward=sum(outward.values());net=sum(totals.values())
    total_inlet_inward=inward['inlet']

    metrics=dict(physical_boundary_signed_outward_volume_flux_m3_s=totals,
        physical_boundary_absolute_volume_flux_m3_s=absolutes,physical_boundary_inward_volume_flux_m3_s=inward,
        physical_boundary_outward_volume_flux_m3_s=outward,total_inward_volume_flux_m3_s=total_inward,
        total_outward_volume_flux_m3_s=total_outward,total_inlet_inward_volume_flux_m3_s=total_inlet_inward,net_outward_volume_flux_m3_s=net,
        relative_all_boundary_mass_imbalance=abs(net)/total_inlet_inward if total_inlet_inward>0 else None,
        mass_normalization='abs(sum actual signed outward phi on every physical boundary face)/sum magnitude actual inward phi on inlet; rho1kg/m3',
        patch_signed_flow_per_unit_depth_m2_s={n:v/.01 for n,v in totals.items()},
        walls_absolute_volume_flux_m3_s=sum(absolutes[n] for n in ('body','top','bottom')),
        native_internal_face_count=n_internal,native_physical_boundary_face_count=600,native_empty_face_count=31800,
        actual_native_cell_count=len(volumes),sum_cell_signed_phi_m3_s=float(cell_sums.sum()),
        cell_sum_minus_boundary_net_m3_s=float(cell_sums.sum()-net),
        local_divergence_scaled_max_abs=float(np.max(np.abs(local))),local_divergence_scaled_RMS=float(np.sqrt(np.mean(local**2))),
        local_divergence_scaled_volume_weighted_RMS=float(np.sqrt(np.average(local**2,weights=volumes))),
        local_divergence_definition='sum signed outward native phi on every cell face / actual oriented native cell volume * Lref/Uref; Dref/Uref=1s',
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
def cell_mapping(centres,geometry_file,volumes,native_points=None):
    require(centres.shape==(CELLS,3),'Native C count differs')
    require(np.max(np.abs(centres[:,2]))<=1e-10,'Native C not centred in extrusion')
    ix=np.rint(centres[:,0]/.01-.5).astype(int);iy=np.rint((centres[:,1]+.4)/.01-.5).astype(int)
    require(np.all((ix>=0)&(ix<200)&(iy>=0)&(iy<80)),'Native centres outside sparse background grid')
    expected=np.column_stack(((ix+.5)*.01,-.4+(iy+.5)*.01))
    require(np.max(np.abs(centres[:,:2]-expected))<=1e-10,'Native C off canonical0.01m lattice')
    keys=iy*200+ix
    require(len(np.unique(keys))==CELLS,'Native C has duplicate sparse grid cells')
    with np.load(geometry_file,allow_pickle=False) as source:
        require(all(n in source for n in ('nodes','quads','centres','areas','grid_shape','ix','iy')),'Missing canonical source mapping arrays')
        source_c=np.asarray(source['centres']).reshape(CELLS,2);source_a=np.asarray(source['areas']).reshape(CELLS)
        source_nodes=np.asarray(source['nodes']).reshape(16200,2)
        sx=np.asarray(source['ix']).reshape(CELLS).astype(int);sy=np.asarray(source['iy']).reshape(CELLS).astype(int)
        require(np.array_equal(np.asarray(source['grid_shape']).reshape(-1),[80,200]),'Source sparse background shape differs')
        source_quads=np.asarray(source['quads']).reshape(CELLS,4)
    require(np.all((sx>=0)&(sx<200)&(sy>=0)&(sy<80)),'Source sparse indices out of background')
    skeys=sy*200+sx
    require(len(np.unique(skeys))==CELLS and np.array_equal(np.sort(keys),np.sort(skeys)),'Native C sparse occupancy differs from canonical source')
    canonical={(i,j) for j in range(80) for i in range(200) if not(45<=i<55 and 35<=j<45)}
    require(set(zip(sx.tolist(),sy.tolist()))==canonical,'Source sparse occupancy differs from true square hole mask')
    source_grid=np.full(80*200,-1,dtype=np.int64);source_grid[skeys]=np.arange(CELLS)
    native_to_source=source_grid[keys]
    require(np.array_equal(np.sort(native_to_source),np.arange(CELLS)),'Native-to-source bijection failed')
    require(np.max(np.abs(centres[:,:2]-source_c[native_to_source]))<=1e-10,'Native C differs from actual source centre')
    require(np.isfinite(source_a).all() and np.all(source_a>0) and np.max(np.abs(volumes/.01-source_a[native_to_source]))<=1e-12,'Native/source mapped areas differ')
    native_grid=np.full(80*200,-1,dtype=np.int64);native_grid[keys]=np.arange(CELLS)
    mapping=dict(native_cell_to_source_cell=native_to_source,native_cell_index_grid=native_grid.reshape(80,200),
        source_cell_index_grid=source_grid.reshape(80,200),occupancy_mask=(source_grid>=0).reshape(80,200),ix=ix,iy=iy)
    report=dict(native_cells=CELLS,source_cells=CELLS,unique_sparse_grid_keys=CELLS,background_grid_cells=16000,void_background_cells=100,
        centre_max_mapping_error_m=float(np.max(np.abs(centres[:,:2]-source_c[native_to_source]))),
        native_source_area_max_error_m2=float(np.max(np.abs(volumes/.01-source_a[native_to_source]))),
        native_to_source_bijective=True,native_to_sparse_grid_bijective=True,ordering_assumed=False,void_indices=-1)
    if native_points is not None:
        points=np.asarray(native_points)
        require(points.shape==(POINTS,3),'Native vertex count differs')
        xyints=np.rint(points[:,:2]/.01).astype(int);nodeints=np.rint(source_nodes/.01).astype(int)
        require(np.max(np.abs(points[:,:2]-xyints*.01))<=1e-10 and np.max(np.abs(source_nodes-nodeints*.01))<=1e-10,'Native/source vertices off canonical lattice')
        require(np.all((xyints[:,0]>=0)&(xyints[:,0]<=200)&(xyints[:,1]>=-40)&(xyints[:,1]<=40)),'Native vertex grid bounds differ')
        nkeys=(nodeints[:,1]+40)*201+nodeints[:,0];pkeys=(xyints[:,1]+40)*201+xyints[:,0]
        require(len(np.unique(nkeys))==16200,'Canonical source vertices duplicate')
        lookup={int(key):i for i,key in enumerate(nkeys)}
        require(all(int(key) in lookup for key in pkeys),'Native projected vertex absent from actual source')
        mapped=np.array([lookup[int(key)] for key in pkeys],dtype=np.int64)
        low=np.abs(points[:,2]+.005)<=1e-10;high=np.abs(points[:,2]-.005)<=1e-10
        require(np.all(low^high) and low.sum()==high.sum()==16200,'Native vertices not exactly two extrusion planes')
        require(np.array_equal(np.sort(mapped[low]),np.arange(16200)) and np.array_equal(np.sort(mapped[high]),np.arange(16200)),'Each native plane must biject to all source vertices')
        require(np.max(np.abs(points[:,:2]-source_nodes[mapped]))<=1e-10,'Native-source projected vertex coordinates differ')
        mapping['native_point_to_source_node']=mapped;mapping['native_points']=points
        report.update(projected_native_vertices=POINTS,source_vertices_per_plane=16200,native_to_source_vertex_bijection_each_plane=True,
            vertex_max_mapping_error_m=float(np.max(np.abs(points[:,:2]-source_nodes[mapped]))))
    return mapping,report
def check_wall_fields(u_raw,p_raw,patches):
    ub=boundary_blocks(u_raw);pb=boundary_blocks(p_raw);expected=np.array([.1,0,0])
    require(re.search(r'\btype\s+fixedValue\s*;',ub['inlet']),'Actual native inlet velocity type differs')
    require(np.max(np.abs(parse_value(ub['inlet'],patches['inlet']['count'],3,'value')-expected))<=1e-12,'Actual native inlet velocity differs')
    require(re.search(r'\btype\s+zeroGradient\s*;',pb['inlet']),'Native inlet pressure type differs')
    require(re.search(r'\btype\s+zeroGradient\s*;',ub['outlet']),'Actual native outlet U must be zeroGradient')
    require(re.search(r'\btype\s+fixedValue\s*;',pb['outlet']),'Actual native outlet pressure type differs')
    require(np.max(np.abs(parse_value(pb['outlet'],patches['outlet']['count'],1,'value')))<=1e-12,'Actual native outlet pressure not zero')
    require(re.search(r'\btype\s+noSlip\s*;',ub['body']),'Actual body must be noSlip')
    if re.search(r'\bvalue\s+',ub['body']):
        require(np.max(np.abs(parse_value(ub['body'],patches['body']['count'],3,'value')))<=1e-12,'Explicit native noSlip body value differs from zero')
    require(re.search(r'\btype\s+zeroGradient\s*;',pb['body']),'Native body pressure type differs')
    for name in ('top','bottom'):
        require(re.search(r'\btype\s+symmetry\s*;',ub[name]) and re.search(r'\btype\s+symmetry\s*;',pb[name]),'Native top/bottom symmetry U/p differs')
    require(re.search(r'\btype\s+empty\s*;',ub[EMPTY]) and re.search(r'\btype\s+empty\s*;',pb[EMPTY]),'Actual native empty planes differ')

def numeric_history(raw,columns,label):
    nonfinite_check(raw)
    rows=[]
    for line in raw.splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        require(re.fullmatch(r'[\s()0-9eE+.\-]+',line) is not None,'Unexpected native '+label+' history token')
        row=list(map(float,re.findall(FLOAT,line)))
        require(len(row)==columns and all(np.isfinite(row)),'Native '+label+' history shape/finite error')
        rows.append(row)
    a=np.asarray(rows,dtype=float)
    require(a.ndim==2 and a.shape[1]==columns,'Empty native '+label+' history')
    indices=a[:,0]
    require(np.array_equal(indices,np.rint(indices)) and np.all(indices>=0),'Native '+label+' iteration not nonnegative integer')
    require(len(np.unique(indices))==len(indices),'Duplicate native '+label+' iteration')
    require(np.array_equal(indices[indices>0],np.arange(1,2001)),'Require complete native '+label+' iterations1..2000')
    require(np.count_nonzero(indices==0)<=1,'Duplicate zero probe in actual history')
    return a[indices>0],dict(native_rows=len(a),retained_positive_rows=2000,initial_iteration0_present=bool(np.any(indices==0)),native_initial_row=a[indices==0].tolist())

def force_history(force_raw,coefficient_raw):
    require('pressure' in force_raw and 'viscous' in force_raw,'Native raw force header missing components')
    require(all(re.search(r'\b'+name+r'\b',coefficient_raw) for name in ('Cm','Cd','Cl')),'Native coefficient header differs')
    f,fr=numeric_history(force_raw,13,'forces');c,cr=numeric_history(coefficient_raw,6,'forceCoeffs')
    require(np.array_equal(f[:,0],c[:,0]),'Native force/coefficient iteration mismatch')
    pressure=f[:,1:4];viscous=f[:,4:7];total=pressure+viscous
    derived=total[:,:2]/5e-6;error=np.abs(derived-c[:,2:4])
    window=(c[:,0]>=1801)&(c[:,0]<=2000);require(window.sum()==200,'Incomplete predeclared final200 force window')
    mean_Cd=float(c[window,2].mean());mean_Cl=float(c[window,3].mean());span=float(np.ptp(c[window,2]));ratio=span/abs(mean_Cd) if mean_Cd!=0 else None
    report=dict(force_rows=2000,coefficient_rows=2000,first_iteration=1,last_iteration=2000,raw_force_history=fr,native_coefficient_history=cr,
      density_kg_m3=1,Uref_m_s=.1,Dref_m=.1,Aref_m2=.001,depth_m=.01,coefficient_denominator_N=5e-6,patches=['body'],CofR_m=[.5,0,0],
      raw_force_derived_CdCl_max_absolute_error=float(error.max()),raw_force_derived_Cd_max_absolute_error=float(error[:,0].max()),raw_force_derived_Cl_max_absolute_error=float(error[:,1].max()),
      final_window_first_iteration=1801,final_window_last_iteration=2000,final_window_rows=200,last200_mean_Cd=mean_Cd,last200_mean_Cl=mean_Cl,last200_Cd_peak_to_peak=span,
      last200_Cd_peak_to_peak_over_abs_mean=ratio,last200_Cd_relative_range_defined=mean_Cd!=0,
      last200_Cd_relative_range_undefined_reason=None if mean_Cd!=0 else 'Mean Cd is exactly zero; relative range undefined; declared numerical check false',
      native_force_units='N on body at native extrusion depth .01m',native_pressure='kinematic p times rho gives Pa',
      official_force_formula='pressureF=sum(rho*p*Sf); viscousF=sum(area*rho*devSigma); native wallShearStress=-devSigma',
      paired_mean_Cd_comparison_target=.01,paired_mean_Cd_comparison_performed=False)
    arrays=dict(force_iteration=f[:,0].astype(np.int64),Cm=c[:,1],Cd=c[:,2],Cl=c[:,3],Cl_front=c[:,4],Cl_rear=c[:,5],pressure_force_N=pressure,
      viscous_force_N=viscous,total_force_N=total,pressure_moment_N_m=f[:,7:10],viscous_moment_N_m=f[:,10:13],raw_force_derived_Cd=derived[:,0],raw_force_derived_Cl=derived[:,1])
    return arrays,report

def body_field_forces(p_raw,p,shear_raw,patches,geometry):
    # Read actual native body-face pressure and wallShearStress only; no interior stress inference.
    pb=boundary_blocks(p_raw);wb=boundary_blocks(shear_raw);patch=patches['body'];count=patch['count'];start=patch['start']
    field(shear_raw,CELLS,3,[0,2,-2,0,0,0,0])
    owner=np.asarray(geometry['owner'])[start:start+count]
    require(re.search(r'\btype\s+zeroGradient\s*;',pb['body']),'Body p must retain exact zeroGradient semantics')
    if re.search(r'\bvalue\s+',pb['body']):
        pressure=parse_value(pb['body'],count,1,'value');pressure_source='Actual saved body p value'
        require(np.max(np.abs(pressure-p[owner]))<=1e-10,'Saved zeroGradient body pressure differs from owner')
    else:
        pressure=p[owner];pressure_source='Exact native zeroGradient body p equals actual native owner cell p'
    shear=parse_value(wb['body'],count,3,'value')
    points=np.asarray(geometry['points']);facepoints=points[np.asarray(geometry['faces'])[start:start+count]]
    Sf=.5*np.sum(np.cross(facepoints,np.roll(facepoints,-1,axis=1)),axis=1);area=np.linalg.norm(Sf,axis=1)
    require(np.all(area>0) and np.linalg.norm(Sf.sum(axis=0))<=1e-14,'Actual body surface not closed or positive')
    return (pressure[:,None]*Sf).sum(axis=0),-(shear*area[:,None]).sum(axis=0),dict(actual_body_faces=count,body_area_m2=float(area.sum()),
      actual_outward_face_vector_sum_m2=Sf.sum(axis=0).tolist(),pressure_boundary_source=pressure_source,
      shear_source='Actual body patch value of native wallShearStress; kinematic traction m2/s2; negative rho*area integral gives force on body',interior_wall_stress_claimed=False)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();case=checked_case(args.run);out=args.out.resolve()
    require(not out.exists() and out.drive.upper()=='E:','Require fresh E: analysis output')
    require(case not in out.parents and out!=case,'Read-only analysis writes outside native case')
    e=Evidence();checks={}
    for path in (Path(__file__),Path(__file__).with_name('prepare_openfoam_square_wake_v003.py')):e.touch(path)
    declaration=e.read_json(case/'declared-comparison.json');status=e.read_json(case/'run-status.json')
    require(status.get('stage')=='finished_execution' and status.get('final_iteration')==2000,'Fixed-cap native execution incomplete')
    for name,expected in dictionaries().items():require(e.text(case/name)==expected,'Actual native initial dictionary differs '+name)
    require(declaration['initial_U_m_s']==[0,0,0] and declaration['initial_p']==0,'Matched declared initial state differs')
    require(np.max(np.abs(field(e.text(case/'0/U'),CELLS,3,[0,1,-1,0,0,0,0])))==0 and np.max(np.abs(field(e.text(case/'0/p'),CELLS,1,[0,2,-2,0,0,0,0])))==0,'Actual native initial state differs from zero')
    checks['actual_native_dictionaries_and_zero_initial_U_p_match_exact_settings']=True
    prepared=e.read_json(case/'provenance/prepared-input-sha256.json')
    for name,digest in prepared.items():e.touch(case/name);require(sha(case/name)==digest,'Prepared source changed '+name)
    require(sha(e.touch(case/'input-square-wake.msh'))==declaration['mesh_source_sha256'],'Shared source mesh hash differs')
    require(sha(e.touch(case/'provenance/geometry.npz'))==declaration['geometry_sha256'],'Shared geometry hash differs')
    require(sha(e.touch(case/'provenance/FORCE_TARGETS.md'))==declaration['force_target_document_sha256'],'Predeclared force target document changed')
    checks['prepared_sources_shared_mesh_and_predeclared_force_targets_hashes_match']=True
    preflight=e.read_json(case/'provenance/preflight-validation.json');require(preflight.get('passed') is True,'Native preflight not passed')
    poly=native_poly(case)
    for name,digest in preflight['native_mesh_sha256'].items():e.touch(poly/name);require(sha(poly/name)==digest,'Native static mesh changed')
    for name,digest in preflight['native_initial_C_Vc_sha256'].items():e.touch(case/'0'/name);require(sha(case/'0'/name)==digest,'Initial C/Vc changed')
    patches=native_patches(poly);geometry=native_geometry(poly);oriented=np.asarray(geometry['volumes'])
    quality=mesh_quality(e.text(case/'mesh-check.log'))
    volumes=field(e.text(case/'2000/Vc'),CELLS,1,[0,3,0,0,0,0,0]);initial_volumes=field(e.text(case/'0/Vc'),CELLS,1,[0,3,0,0,0,0,0])
    require(np.all(volumes>0) and np.max(np.abs(volumes-oriented))<=1e-14 and np.max(np.abs(volumes-initial_volumes))<=1e-14,'Actual native Vc differs from oriented/static volumes')
    centres=field(e.text(case/'2000/C'),CELLS,3,[0,1,0,0,0,0,0]);initial_centres=field(e.text(case/'0/C'),CELLS,3,[0,1,0,0,0,0,0])
    require(np.max(np.abs(centres-initial_centres))<=1e-12,'Actual native C changed on static mesh')
    geometry['volumes']=volumes.tolist()
    checks['full_checkMesh_positive_oriented_and_native_Vc_exact_counts_passed']=True
    mapping,mapping_record=cell_mapping(centres,case/'provenance/geometry.npz',volumes,geometry['points'])
    checks['actual_C_source_and_true_square_hole_grid_mappings_bijective']=True
    with np.load(case/'provenance/geometry.npz',allow_pickle=False) as canonical:source_quads=np.asarray(canonical['quads'])
    cell_vertices=[set() for _ in range(CELLS)]
    for face_index,face in enumerate(geometry['faces']):
        cell_vertices[geometry['owner'][face_index]].update(face)
        if face_index<INTERNAL_FACES:cell_vertices[geometry['neighbour'][face_index]].update(face)
    point_map=mapping['native_point_to_source_node'];native_points=np.asarray(geometry['points'])
    for index,vertices in enumerate(cell_vertices):
        require(len(vertices)==8,'Actual native cell must have8vertices')
        expected=set(source_quads[mapping['native_cell_to_source_cell'][index]].tolist())
        for sign in (-1,1):require({int(point_map[v]) for v in vertices if sign*native_points[v,2]>0}==expected,'Actual native cell vertices differ from source quad')
    checks['actual_native_cell_vertices_match_canonical_source_quads_both_planes']=True
    forces,force_record=force_history(e.text(case/'postProcessing/bodyForces/0/forces.dat'),e.text(case/'postProcessing/bodyCoefficients/0/forceCoeffs.dat'))
    checks['complete_actual_native_force_and_coefficient_histories_1_to_2000']=True
    checks['native_raw_force_derived_Cd_Cl_absolute_crosscheck_le_1e_8']=force_record['raw_force_derived_CdCl_max_absolute_error']<=1e-8
    ratio=force_record['last200_Cd_peak_to_peak_over_abs_mean']
    checks['last200_Cd_peak_to_peak_over_abs_mean_le_1e_3']=ratio is not None and ratio<=1e-3
    checks['last200_abs_mean_Cl_le_1e_3']=abs(force_record['last200_mean_Cl'])<=1e-3
    data={};stats={};continuity={};field_force_records=[];field_pressure=[];field_viscous=[]
    for iteration in range(100,2001,100):
        u_raw=e.text(case/str(iteration)/'U');p_raw=e.text(case/str(iteration)/'p');phi_raw=e.text(case/str(iteration)/'phi');shear_raw=e.text(case/str(iteration)/'wallShearStress')
        u=field(u_raw,CELLS,3,[0,1,-1,0,0,0,0]);p=field(p_raw,CELLS,1,[0,2,-2,0,0,0,0]);check_wall_fields(u_raw,p_raw,patches)
        metrics,local=flux_and_divergence(phi_raw,patches,geometry)
        fp,fv,fr=body_field_forces(p_raw,p,shear_raw,patches,geometry);field_pressure.append(fp);field_viscous.append(fv)
        pe=float(np.max(np.abs(fp-forces['pressure_force_N'][iteration-1]))/5e-6);ve=float(np.max(np.abs(fv-forces['viscous_force_N'][iteration-1]))/5e-6)
        te=float(np.max(np.abs(fp+fv-forces['total_force_N'][iteration-1]))/5e-6)
        field_force_records.append(dict(iteration=iteration,pressure_force_N=fp.tolist(),viscous_force_N=fv.tolist(),total_force_N=(fp+fv).tolist(),
          pressure_coefficient_equivalent_max_error=pe,viscous_coefficient_equivalent_max_error=ve,total_coefficient_equivalent_max_error=te,**fr))
        p_mean=float(np.average(p,weights=volumes))
        stats[str(iteration)]=dict(U_min_m_s=u.min(axis=0).tolist(),U_max_m_s=u.max(axis=0).tolist(),speed_min_m_s=float(np.linalg.norm(u,axis=1).min()),speed_max_m_s=float(np.linalg.norm(u,axis=1).max()),
          raw_gauge_pressure_min_Pa=float(p.min()),raw_gauge_pressure_max_Pa=float(p.max()),area_mean_pressure_Pa=p_mean,max_abs_out_of_plane_U_m_s=float(np.abs(u[:,2]).max()))
        continuity[str(iteration)]=metrics;data[iteration]=dict(U=u,p=p,p_gauge_Pa=p,local_divergence_scaled=local)
    checks['all20_saved_native_U_p_phi_wallShearStress_finite_counts_dimensions_and_boundary_types']=True
    max_field_error=max(max(r[n] for n in ('pressure_coefficient_equivalent_max_error','viscous_coefficient_equivalent_max_error','total_coefficient_equivalent_max_error')) for r in field_force_records)
    checks['native_body_pressure_shear_force_integrals_coefficient_equivalent_crosscheck_le_1e_8']=max_field_error<=1e-8
    final=continuity['2000'];flux=final['physical_boundary_signed_outward_volume_flux_m3_s']
    checks['final_native_local_divergence_max_le_1e_5']=final['local_divergence_scaled_max_abs']<=1e-5
    checks['final_native_local_divergence_RMS_le_1e_6']=final['local_divergence_scaled_RMS']<=1e-6
    checks['actual_inlet_inward_and_outlet_outward']=flux['inlet']<0 and flux['outlet']>0
    mass=final['relative_all_boundary_mass_imbalance'];checks['relative_all_boundary_mass_imbalance_le_1e_4']=mass is not None and mass<=1e-4
    checks['native_body_and_symmetry_absolute_normal_flux_le_1e_14_m3_s']=final['walls_absolute_volume_flux_m3_s']<=1e-14
    checks['actual_native_2D_Uz_max_le_1e_12']=stats['2000']['max_abs_out_of_plane_U_m_s']<=1e-12
    residuals=residual_history(e.text(case/'solver.log'));terminal={name:residuals[-1][name] for name in ('p','Ux','Uy')}
    checks['terminal_printed_p_Ux_Uy_initial_residuals_le_1e_6']=max(terminal.values())<=1e-6
    delta=data[2000]['U']-data[1800]['U'];denominator=float(np.sum(volumes[:,None]*data[2000]['U']**2))
    relative=float(np.sqrt(np.sum(volumes[:,None]*delta**2)/denominator)) if denominator>0 else None
    checks['U_1800_to_2000_volume_weighted_relative_L2_le_1e_3']=relative is not None and relative<=1e-3
    auth=e.read_json(case/'provenance/root-run-authorization.json')
    require(auth.get('approved') is True and auth.get('approved_by')=='root' and auth.get('action')=='run_square_wake_openfoam_v002' and auth.get('allowed_mode')=='Run','Missing exact root run authorization')
    executed=e.touch(case/'provenance/run-openfoam-executed.sh');original=e.touch(case/'provenance/run-openfoam-v003.sh')
    require(sha(executed)==sha(original)==auth.get('wrapper_sha256'),'Executed wrapper hash differs')
    require(sha(case/'provenance/prepare_openfoam_square_wake_v003.py')==sha(Path(__file__).with_name('prepare_openfoam_square_wake_v003.py'))==auth.get('preparation_sha256'),'Executed/current preparer hash differs')
    require(sha(case/'provenance/Launch-OpenFOAM-v003.ps1')==auth.get('launcher_sha256'),'Reviewed launcher hash differs')
    checks['reviewed_execution_sources_match_root_authorization']=True
    installed_rows=[line.split(None,1) for line in e.text(case/'provenance/installed-source-sha256.txt').splitlines() if line.strip()]
    require(len(installed_rows)==len(INSTALLED_SOURCES),'Installed source count differs')
    require({source_path for digest,source_path in installed_rows}=={'/opt/openfoam14/'+n for n in INSTALLED_SOURCES},'Installed source path set differs')
    for digest,source_path in installed_rows:
        require(sha(e.touch(case/'provenance/installed-sources'/source_path[len('/opt/openfoam14/'):]))==digest,'Copied installed source hash differs')
    checks['installed_official_force_shear_solver_and_dictionary_sources_match_recorded_hashes']=True
    source=e.finish();checks['inputs_unchanged_during_read_only_analysis']=True
    checks={name:bool(value) for name,value in checks.items()}
    record=dict(solver='OpenFOAM Foundation14 incompressibleFluid',actual_settings=declaration,steady_iteration_indices_are_not_physical_time=True,
      mesh_cells=CELLS,mesh_quality=quality,actual_native_patches=patches,native_geometry=dict(points=POINTS,faces=TOTAL_FACES,internal_faces=INTERNAL_FACES,
        cell_volume_min_m3=float(volumes.min()),cell_volume_max_m3=float(volumes.max()),total_volume_m3=float(volumes.sum()),bounds_m=geometry['bounds']),
      actual_cell_mapping=mapping_record,native_fields=stats,native_continuity_all_saved_states=continuity,final_native_continuity=final,
      force_history=force_record,native_body_field_force_crosscheck=dict(maximum_coefficient_equivalent_error=max_field_error,saved_state_records=field_force_records),
      residuals=dict(final_maximum_printed_initial_residual_per_variable=terminal,last100_maximum_printed_initial_residual_per_variable={name:max(row[name] for row in residuals[-100:]) for name in terminal},
        definition='Maximum printed initial linear-system residual per variable in each steady SIMPLE iteration',norm_equivalent_to_fluent=False),
      U_stability=dict(from_iteration=1800,to_iteration=2000,volume_weighted_relative_L2=relative,definition='sqrt(sum(Vc*|U2000-U1800|^2)/sum(Vc*|U2000|^2))'),
      final_iteration=2000,saved_iterations=list(range(100,2001,100)),selected_iterations=[500,1000,1500,1800,2000],
      native_wall_value_semantics='Native body noSlip fixes stationary U0, explicit value checked if present; top/bottom symmetry imposes zero normal flux',
      pressure_comparison='Raw gauge Pa=rho*p with actual rho1 and common outlet0',
      checks=checks,check_results=checks,failed_checks=[n for n,v in checks.items() if not v],all_checks_passed=all(checks.values()),all_declared_checks_passed=all(checks.values()),
      physical_validation=False,mesh_independence_established=False,benchmark_validation=False,unbounded_freestream_established=False,blockage_ratio=.125,domain_height_over_D=8,
      transient_shedding_simulated=False,interior_wall_stress_claimed=False,interpretation='AI-assisted educational steady one-mesh confined square-obstacle wake; all numerical misses retained at fixed2000 cap')
    out.mkdir(parents=True)
    arrays=dict(centres=centres,volumes=volumes,areas=volumes/.01,grid_shape=np.array([80,200]),**mapping,**data[2000],**forces,
      iterations=np.arange(100,2001,100),saved_iterations=np.arange(100,2001,100),U_history=np.stack([data[i]['U'] for i in range(100,2001,100)]),p_history=np.stack([data[i]['p'] for i in range(100,2001,100)]),
      initial_U=np.zeros((CELLS,3)),initial_p=np.zeros(CELLS),residual_iteration=np.array([r['iteration'] for r in residuals]),residual_names=np.array(['p','Ux','Uy']),
      residual_values=np.array([[r[n] for n in ('p','Ux','Uy')] for r in residuals]),body_field_force_iteration=np.arange(100,2001,100),
      body_field_pressure_force_N=np.asarray(field_pressure),body_field_viscous_force_N=np.asarray(field_viscous),rho_kg_m3=np.array(1.),reference_velocity_m_s=np.array(.1),reference_length_m=np.array(.1),reference_area_m2=np.array(.001),depth_m=np.array(.01))
    for iteration in (500,1000,1500,1800,2000):
        for key,value in data[iteration].items():arrays[key+'_'+str(iteration)]=value
    np.savez_compressed(out/'fields.npz',**arrays)
    header=['iteration','Cm','Cd','Cl','Fx_pressure_N','Fy_pressure_N','Fz_pressure_N','Fx_viscous_N','Fy_viscous_N','Fz_viscous_N','Fx_total_N','Fy_total_N','Fz_total_N']
    with (out/'force-history.csv').open('x',encoding='ascii',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(header)
        for j in range(2000):writer.writerow([int(forces['force_iteration'][j]),forces['Cm'][j],forces['Cd'][j],forces['Cl'][j],*forces['pressure_force_N'][j],*forces['viscous_force_N'][j],*forces['total_force_N'][j]])
    with (out/'residual-history.csv').open('x',encoding='ascii',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(residuals[0]));writer.writeheader();writer.writerows(residuals)
    record['audited_fields_npz_sha256']=sha(out/'fields.npz');record['force_history_csv_sha256']=sha(out/'force-history.csv')
    record['analysis_source_sha256']=sha(Path(__file__));record['preparation_source_sha256']=sha(Path(__file__).with_name('prepare_openfoam_square_wake_v003.py'))
    source['audited_outputs']={name:dict(sha256=sha(out/name),bytes=(out/name).stat().st_size) for name in ('fields.npz','force-history.csv','residual-history.csv')}
    json_new(out/'checks.json',record);json_new(out/'source-manifest.json',source)
    for path in (Path(__file__),Path(__file__).with_name('prepare_openfoam_square_wake_v003.py')):copied(path,out/('executed-'+path.name))
    print(json.dumps(dict(output=str(out),all_checks_passed=record['all_checks_passed'],failed_checks=record['failed_checks'],fields_npz_sha256=record['audited_fields_npz_sha256'],force_history_csv_sha256=record['force_history_csv_sha256']),indent=2))
    return 0 if record['all_checks_passed'] else 2
if __name__=='__main__':
    sys.exit(main())
