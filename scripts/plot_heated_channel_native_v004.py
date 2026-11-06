"""Actual601-state paired heated-channel production fields and histories.

Main61 frames select native steps0,10,...600. Supplementary101 frames select
native0..100 through0.5 s. No time interpolation, particles or forced motion.
Execution requires root-selected complete audits and exact native-source hashes.
OpenFOAM temperatures are reread from both literal native ranks and assembled
with actual cellProcAddressing. Partial slow serial exports stay diagnostics.
"""
from pathlib import Path
import argparse,csv,gzip,hashlib,json,os,re,sys
os.environ['TEMP']=os.environ['TMP']='E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-heated-channel-paired'
import h5py,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw
import heated_channel_fluent_production_common_v001 as fr
import plot_cavity_native_v003 as common
import plot_heated_channel_native_scaffold_v001 as layout
os.environ['MPLCONFIGDIR']='E:/CFD-Learning-Lab/temp/matplotlib-heated-channel-paired'
require,sha=common.require,common.sha
BOUNDS=('inlet','outlet','top','bottom');EXPECTED_STEPS=np.arange(601)
FROZEN={'reader':'b85769c76ca3d992ba2e491865c3e813e8436a2baca8d1d30b4775f775213e09','layout':'d5970e22c16877649fcf23db5986a1139a22194118188a773e9405c3d35bc21f',
        'cavity':'c2425654f684bc8de8d9219dffa872d8b3c80a605bb35fce2b9070362f388445','native':'4596d6a6b975abc6110d80fdb8b5b680878cdff18022ce969866b288aacb2f02'}


def read_audit(name,path,npz,digest):
    record=json.loads(path.read_text(encoding='utf-8-sig'));raw_flags=record['checks'];flags=dict(raw_flags);normalization=None
    key='all601_native_states_and600_step_metric_records'
    if name=='fluent' and type(flags.get(key)) is str and flags[key]=='True':
        normalization=dict(key=key,raw_value=flags[key],raw_type='str',normalized_value=None,independent_native_proof_completed=False)
        flags[key]=True  # Pending independent full native proof below; never a generic truthiness conversion.
    require(flags and all(type(v) is bool for v in flags.values()),'Literal boolean audit flags required except the one admitted Fluent serialization key')
    failed=[k for k,v in flags.items() if not v];require(record['all_checks_passed']==(not failed) and set(record.get('failed_checks',failed))==set(failed),'Inconsistent audit flags/aggregate')
    require(not record.get('errors',{}),'Structural/native parser audit errors block selection')
    manifest=path.parent/('fluent-thermal-source-manifest.json' if name=='fluent' else 'source-manifest.json');m=json.loads(manifest.read_text(encoding='utf-8-sig'))
    rows=m['inputs'] if name=='fluent' else m['sources'];require(m['all_inputs_unchanged' if name=='fluent' else 'inputs_unchanged_during_analysis'] is True,'Audit source hashes were not stable')
    hashes={}
    for row in rows:
        p=Path(row['path']).resolve()
        if 'sha256_before' in row:
            require(row['sha256_before']==row['sha256_after'],'Native before/after digest differs');h=row['sha256_before']
        else:h=row['sha256']
        require(str(p) not in hashes and p.is_file(),'Duplicate/missing audit input');hashes[str(p)]=h
    bound=m['outputs'].get(npz.name) if name=='fluent' else record['audited_fields_npz_sha256']
    if name=='openfoam':require(m['audited_outputs'][npz.name]['sha256']==bound,'OpenFOAM manifest/NPZ binding differs')
    require(bound==digest==sha(npz),'Selected actual NPZ not audit-hash bound')
    require(all(sha(Path(p))==h for p,h in hashes.items()),'Audit-bound native/source input changed')
    status=dict(path=str(path),source_manifest_path=str(manifest),native_source_hashes=hashes,checks=flags,failed_checks=failed,all_declared_checks_passed=not failed,raw_check_flags=raw_flags,nonboolean_check_normalization=normalization)
    return record,status


def geometry(path):
    g=dict(np.load(path,allow_pickle=False));fr.geometry_check(g)
    require(sha(path)=='23209a3b5d858075a85a053618c5ab0005acef3121f2d798b9f0c4c335212920','Canonical geometry changed')
    require(np.allclose(g['x_edges'],np.arange(97)/96,rtol=0,atol=1e-14) and np.allclose(g['y_edges'],np.arange(25)*(.1/24),rtol=0,atol=1e-14),'Channel extent or lattice differs')
    return g


def common_proof():
    return dict(all_native_temperatures_independently_verified=True,actual_all_state_frozen_U_p_and_native_face_flux_verified=True,
        native_boundary_power_semantics_proved=True,initial_heat_unsolved_and_not_padded=True)


def finish_bundle(b,g):
    b.validate();require(np.array_equal(b.times_s[1:],b.solved_times_s),'Power times must be solved states only')
    energy=np.sum(1000*g['areas'][None,:]*(b.temperature_K-300.),axis=1)
    b.energy_recompute_max_abs_difference_J_per_m=float(np.max(abs(energy-b.energy_J_per_m)))
    require(b.energy_recompute_max_abs_difference_J_per_m<1e-8,'Native temperature/volume energy differs from audited history')
    storage=np.diff(b.energy_J_per_m)/.005;b.storage_recompute_max_abs_difference_W_per_m=float(np.max(abs(storage-b.storage)))
    require(b.storage_recompute_max_abs_difference_W_per_m<1e-6,'Actual backward-Euler storage differs')
    normalization=np.maximum.reduce((np.sum(abs(b.native_Qin_W_per_m),axis=1),abs(b.storage),np.full(600,5000.)))
    native=abs(b.storage-b.native_Qin_W_per_m.sum(axis=1))/normalization
    reconstructed=abs(b.storage-b.reconstructed_Qin_W_per_m.sum(axis=1))/normalization
    require(np.max(abs(native-b.native_Euler_relative_defect))<1e-10 and np.max(abs(reconstructed-b.reconstructed_PDE_relative_defect))<1e-10,'Declared native/reconstructed defect histories differ')
    b.normalization=normalization;b.actual_T_delta=np.max(abs(np.diff(b.temperature_K,axis=0)),axis=1)
    return b


def fluent_bundle(run,npz,status,g):
    z=np.load(npz,allow_pickle=False);require(z['physical_boundary_names'].astype(str).tolist()==list(BOUNDS),'Fluent boundary order differs')
    maps={};fieldT=[];fieldH=[];changes=[];raw=[];times=[];steps=[];dts=[]
    flow=fr.snapshot(run/'flow-final.cas.h5',run/'flow-final.dat.h5',g,False)
    require(all(np.array_equal(z['frozen_'+key],flow['common'][key]) for key in ('u','v','p')) and np.array_equal(z['frozen_face_mass_flux_kg_s'],flow['flux']) and np.array_equal(z['native_cell_to_common'],flow['mapping']),'Audited frozen Fluent flow/mapping differs from actual native source')
    needed=[run/'flow-final.cas.h5',run/'flow-final.dat.h5',run/'thermal-initial.cas.h5',run/'inherited-step000001.cas.h5',run/'production-start.cas.h5',run/'thermal-final.cas.h5',run/'thermal-final.dat.h5']
    hashes=status['native_source_hashes'];require(all(str(p.resolve()) in hashes for p in needed),'Selected Fluent native cases not bound')
    for step in range(601):
        if step%40==0:common.resource_guard()
        case=run/('thermal-initial.cas.h5' if step==0 else 'inherited-step000001.cas.h5' if step==1 else 'production-start.cas.h5');data=run/'states'/f'step-{step:06d}.dat.h5'
        require(str(data.resolve()) in hashes,'Actual saved Fluent state absent from audit source binding')
        if case not in maps:
            snap=fr.snapshot(case,data,g,True);fr.frozen_compare(flow,snap);maps[case]=snap['mapping']
        mapping=maps[case]
        with h5py.File(data,'r') as handle:
            cell=handle['results/1/phase-1/cells'];T=np.empty(2304);H=np.empty(2304);T[mapping]=fr.chunk_array(cell['SV_T'],2304);H[mapping]=fr.chunk_array(cell['SV_H'],2304)
            delta=[]
            for key,native_key in (('u','SV_U'),('v','SV_V'),('p','SV_P')):
                value=np.empty(2304);value[mapping]=fr.chunk_array(cell[native_key],2304);delta.append(float(abs(value-flow['common'][key]).max()))
            flux=fr.chunk_array(handle['results/1/phase-1/faces/SV_FLUX'],4728);delta.append(float(abs(flux-flow['flux']).max()))
        fieldT.append(T);fieldH.append(H);changes.append(delta);t=fr.native_time(data);times.append(t['flow-time']);steps.append(t['time-step']);dts.append(t['physical-time-step'])
        if step:
            _,Q,_=fr.boundary_array(data,'SV_HEAT_FLUX',flow['zones']);raw.append([Q[n] for n in BOUNDS])
    T,H=np.array(fieldT),np.array(fieldH);times=np.array(times);steps=np.array(steps);changes=np.array(changes);raw=np.array(raw)
    require(T.shape==H.shape==(601,2304) and times.shape==steps.shape==np.array(dts).shape==(601,) and changes.shape==(601,4) and raw.shape==(600,4),'All actual Fluent native production array shapes required')
    require(all(np.isfinite(a).all() for a in (T,H,times,steps,np.array(dts),changes,raw)),'Actual native Fluent production arrays must be finite')
    require(np.array_equal(T,z['temperature_K']) and np.array_equal(H,z['native_enthalpy_J_kg']),'All actual native Fluent T/h must equal audited arrays')
    require(np.array_equal(times,z['times_s']) and np.array_equal(steps,z['native_step']) and np.array_equal(np.array(dts),z['native_dt_s']),'Actual native Fluent time/step/dt differs')
    require(np.array_equal(changes,z['frozen_field_maximum_absolute_changes_SI']) and np.array_equal(raw,z['native_raw_boundary_Qin_W']),'Actual native frozen-field or raw integrated-W history differs')
    h300=float(H[0].mean());mass=np.array([flow['flux'][slice(*flow['zones'][n])].sum() for n in BOUNDS]);Q=raw+h300*mass[None,:]
    require(np.max(abs(Q-z['native_reference_corrected_boundary_Qin_W_per_m']))<1e-8,'Fluent h300 reference correction differs; no extra convection/area factor allowed')
    b=layout.NativeProduction('fluent',times,steps,T,np.column_stack((z['frozen_u'],z['frozen_v'])),z['energy_J_per_m'],times[1:],z['native_reference_corrected_boundary_Qin_W_per_m'],z['reconstructed_boundary_Qin_W_per_m'],z['native_Euler_relative_defect'],z['reconstructed_PDE_relative_defect'],z['energy_final_scaled_residual'],times,z['outlet_bulk_temperature_K'],common_proof(),status)
    b.storage=z['storage_W_per_m'];b.raw_native_Q_W=raw;b.depth=1.;b.native_diffusion=None;b.native_convection=None;b.frozen_changes=changes;b.inner_iterations=z['actual_inner_iterations'];b.initial_residual=None;b.native_inner_residual_iteration=z['residual_iteration'];b.native_inner_residual=z['residual_energy_scaled']
    b.meta=dict(native_T_values_verified=1384704,native_h_values_verified=1384704,native_h300_J_kg=h300,complete_vertex_set_bijection=True,all601_native_U_V_p_and4728_face_flux_rechecked=True,all600_native_boundary_SV_HEAT_FLUX_integrated_W_rechecked=True,initial_heat_unsolved=True,raw_native_total_heat_reference_correction='Qin_above300 = native_integrated_W + actual_h300 * native_outward_mass_flow; no added convection or area multiplier')
    b.reconstructed_diffusion=np.empty((600,4));b.reconstructed_convection=np.empty((600,4))
    for step in range(1,601):
        path=run/'step-metrics'/f'step-{step:06d}.json';require(str(path.resolve()) in hashes,'Actual Fluent step metric not audit-bound');metric=json.loads(path.read_text())['metrics']
        require(metric['actual_thermal_step']==step and metric['native_time_step1']==dict(zip(('flow-time','time-step','physical-time-step'),(times[step],steps[step],dts[step]))),'Source-bound Fluent metric/native actual time differs')
        component=metric['reconstructed_PDE_components_W_per_m'];b.reconstructed_diffusion[step-1]=[component[n]['diffusion_in_W'] for n in BOUNDS];b.reconstructed_convection[step-1]=[component[n]['convection_in_above300_W'] for n in BOUNDS]
    require(np.isfinite(b.reconstructed_diffusion).all() and np.isfinite(b.reconstructed_convection).all() and np.max(abs(b.reconstructed_diffusion+b.reconstructed_convection-b.reconstructed_Qin_W_per_m))<1e-8,'All600 Fluent source metric reconstructed powers differ from audited arrays')
    if status['nonboolean_check_normalization'] is not None:
        status['nonboolean_check_normalization'].update(normalized_value=True,independent_native_proof_completed=True,proof='All601 finite native T/h/time/dt/frozen-field arrays exactly match audited arrays; all600 source-bound actual metric records match native step/time and reconstructed powers')
    z.close();return finish_bundle(b,g)


def foam_mapping(case,g):
    mesh=common.native.foam_mesh(case);lookup={tuple(np.rint(p/1e-10).astype(np.int64)):i for i,p in enumerate(g['nodes'])}
    ids=np.array([lookup[tuple(np.rint(p[:2]/1e-10).astype(np.int64))] for p in mesh['nodes']]);require(np.max(abs(mesh['nodes'][:,:2]-g['nodes'][ids]))<1e-10,'Native OF vertex differs')
    cells={tuple(v):i for i,v in enumerate(np.sort(g['quads'],axis=1))};mapping=np.array([cells[tuple(np.unique(ids[v]))] for v in mesh['cells']]);require(np.array_equal(np.sort(mapping),np.arange(2304)),'Complete native hex/source-quad mapping not bijective')
    return mapping


def native_labels(path):
    count,body=common.native.foam_list(path);require(re.fullmatch(r'[\s+\-\d]*',body) is not None,'Actual integer native addressing required')
    values=np.fromstring(body,sep=' ',dtype=np.int64);require(values.shape==(count,),'Native addressing count differs');return values


def rank_topology(row,status,g):
    """Independently bind rank cells/vertices to the actual global/source mesh."""
    case=Path(row['native_case']).resolve();global_mesh=common.native.foam_mesh(case);mapping=foam_mapping(case,g)
    hashes=status['native_source_hashes'];require(all(str(p.resolve()) in hashes for p in global_mesh['inputs']),'Actual global mesh is not audit-bound')
    require(len(row['rank_mesh'])==2,'Exactly two literal native rank meshes required')
    node_lookup={tuple(np.rint(p/1e-10).astype(np.int64)):i for i,p in enumerate(g['nodes'])};source_nodes=np.array([node_lookup[tuple(np.rint(p[:2]/1e-10).astype(np.int64))] for p in global_mesh['nodes']])
    source_cells={tuple(v):i for i,v in enumerate(np.sort(g['quads'],axis=1))};ranks=[];proof=[]
    expected={'points','faces','owner','neighbour','boundary','cellProcAddressing','faceProcAddressing','pointProcAddressing'}
    for rank,bindings in enumerate(row['rank_mesh']):
        require(set(bindings)==expected,'Literal rank mesh/addressing schema differs');paths={k:Path(p).resolve() for k,p in bindings.items()}
        poly=case/f'processor{rank}/0/polyMesh';require(all(p==poly/k and str(p) in hashes for k,p in paths.items()),'Actual rank topology/addressing path is not uniquely bound')
        mesh=common.native.foam_mesh(case/f'processor{rank}');require(all(p.resolve()==paths[p.name] for p in mesh['inputs']),'Rank mesh reader selected different native topology')
        cells=native_labels(paths['cellProcAddressing']);points=native_labels(paths['pointProcAddressing']);faces=native_labels(paths['faceProcAddressing'])
        require(len(cells)==len(mesh['cells']) and np.all((cells>=0)&(cells<2304)) and len(np.unique(cells))==len(cells),'Rank/global cell addressing differs')
        require(len(points)==len(mesh['nodes']) and np.all((points>=0)&(points<len(global_mesh['nodes']))) and np.max(abs(mesh['nodes']-global_mesh['nodes'][points]))<1e-12,'Rank/global native point addressing or coordinates differ')
        nfaces,_=common.native.foam_list(paths['owner']);require(faces.shape==(nfaces,) and np.all(faces!=0) and np.all(abs(faces)<=9336),'Actual signed rank face addressing differs')
        projected=np.array([source_cells[tuple(np.unique(source_nodes[points[v]]))] for v in mesh['cells']]);require(np.array_equal(projected,mapping[cells]),'Rank native eight-vertex cells disagree with global/source quad bijection')
        ranks.append(dict(rank=rank,cells=cells,source_cells=projected,processor_case=case/f'processor{rank}'))
        proof.append(dict(rank=rank,native_cell_count=len(cells),native_vertex_sets_match_source_quads=True,point_addressing_and_actual_coordinates_match_global_mesh=True,cellProcAddressing=str(paths['cellProcAddressing']),pointProcAddressing=str(paths['pointProcAddressing']),faceProcAddressing=str(paths['faceProcAddressing'])))
    require(np.array_equal(np.sort(np.concatenate([r['cells'] for r in ranks])),np.arange(2304)),'Two native ranks must cover every global cell exactly once')
    return ranks,mapping,dict(native_case=str(case),rank_meshes=proof,disjoint_full_global_cell_coverage=True)


def openfoam_bundle(npz,record,status,g):
    z=np.load(npz,allow_pickle=False);require(z['boundary_names'].astype(str).tolist()==list(BOUNDS),'OpenFOAM physical boundary order differs')
    rows=record['native_state_paths'];require(len(rows)==601,'Exact601 literal native source rows required')
    mapping=z['native_cell_to_source_cell'].astype(int);require(np.array_equal(np.sort(mapping),np.arange(2304)),'Audited OF mapping not bijective');order=np.argsort(mapping)
    require(np.max(abs(z['centres'][order,:2]-g['centres']))<1e-10 and np.max(abs(z['volumes'][order]/.01-g['areas']))<1e-12,'Native OF C/volumes differ from source')
    cases={};T=[];times=[];topology_proofs=[];time_proofs=0
    for step,row in enumerate(rows):
        if step%40==0:common.resource_guard()
        require(row['index']==step and abs(row['time_s']-step*.005)<1e-10,'Actual OF source row index/time differs')
        require(row['source_kind']=='native_rank_reassembly' and abs(float(row['time_label'])-row['time_s'])<1e-12,'Actual native rank source kind/time differs')
        case=Path(row['native_case']).resolve()
        for key in ('rank_T','rank_U','rank_p','rank_phi','rank_gradT','rank_time_metadata'):
            require(len(row[key])==(0 if step==0 and key in ('rank_gradT','rank_time_metadata') else 2),'Incomplete two-rank native state/metadata paths')
            require(all(str(Path(p).resolve()) in status['native_source_hashes'] for p in row[key]),'Literal rank state/time/gradient input is not audit-bound')
        if case not in cases:
            ranks,actual_mapping,proof=rank_topology(row,status,g);cases[case]=(ranks,row['rank_mesh']);topology_proofs.append(proof)
            require(np.array_equal(actual_mapping,mapping),'Actual rank source case global/source cell order differs from audited mapping')
        ranks,mesh_bindings=cases[case];require(row['rank_mesh']==mesh_bindings,'A native source case changes rank topology/addressing within the history');actual=np.empty(2304)
        for rank in ranks:
            r=rank['rank'];folder=rank['processor_case']/row['time_label']
            require(all(Path(row[key][r]).resolve()==folder/name for key,name in (('rank_T','T'),('rank_U','U'),('rank_p','p'),('rank_phi','phi'))),'Rank field path is not the actual matching rank/time')
            actual[rank['cells']]=common.native.foam_field(Path(row['rank_T'][r]),len(rank['cells']),1,[0,0,0,1,0,0,0])
            if step:
                path=Path(row['rank_time_metadata'][r]).resolve();require(path==folder/'uniform/time' and Path(row['rank_gradT'][r]).resolve()==folder/'gradT','Native solved rank metadata/gradient path differs')
                raw=common.native.foam_text(path);values={key:float(re.search(r'\b'+key+r'\s+('+fr.NUMBER+r')\s*;',raw)[1]) for key in ('value','deltaT','index')}
                require(abs(values['value']-step*.005)<1e-12 and abs(values['deltaT']-.005)<1e-14 and values['index']==step,'Actual rank native time/dt/global step differs');time_proofs+=1
        require(np.array_equal(actual,z['T_history'][step]),'Actual assembled two-rank native T differs from audited global-native-order state')
        T.append(actual[order]);times.append(row['time_s'])
    T=np.array(T);times=np.array(times);require(np.array_equal(T,z['T_history'][:,order]),'All actual native OpenFOAM T values must equal audited arrays')
    require(np.array_equal(times,z['time_s']) and np.array_equal(z['solved_time_s'],times[1:]),'Actual OF saved times differ')
    changeU=np.max(abs(z['U_history']-z['frozen_reference_U'][None,:,:]),axis=(1,2));changeP=np.max(abs(z['p_history']-z['frozen_reference_p'][None,:]),axis=1);changePhi=np.max(abs(z['native_phi_history']-z['frozen_reference_native_phi'][None,:]),axis=1)
    # Audited maxima also cover physical boundaries and both processor-interface
    # sides; internal/global histories cannot replace that broader native proof.
    full_changes=np.column_stack((z['frozen_U_max_abs_change_m_s'],z['frozen_p_max_abs_change'],z['frozen_native_phi_max_abs_change_m3_s']))
    require(np.all(np.column_stack((changeU,changeP,changePhi))<=full_changes),'Internal/global frozen OF changes exceed audited full native maxima')
    b=layout.NativeProduction('openfoam',times,EXPECTED_STEPS.copy(),T,z['frozen_reference_U'][order,:2],z['energy_above300_J_per_m'],z['solved_time_s'],z['native_Qin_W_per_m'],z['reconstructed_Qin_W_per_m'],z['native_Euler_relative_defect'],z['reconstructed_PDE_relative_defect'],z['final_linear_T_residual'],times,z['outlet_bulk_T_K'],common_proof(),status)
    b.storage=z['storage_W_per_m'];b.raw_native_Q_W=None;b.depth=.01;b.native_diffusion=z['native_diffusive_Qin_W_per_m'];b.native_convection=z['native_convective_Qin_W_per_m'];b.reconstructed_diffusion=None;b.reconstructed_convection=None;b.frozen_changes=full_changes;b.inner_iterations=z['linear_T_iterations'];b.initial_residual=z['initial_linear_T_residual'];b.native_inner_residual_iteration=None;b.native_inner_residual=None
    require(np.max(abs(b.native_diffusion+b.native_convection-b.native_Qin_W_per_m))<1e-6,'Actual native OF diffusive/convective Q sum differs')
    require(time_proofs==1200,'All600 solved physical times must have two actual native rank metadata proofs')
    b.meta=dict(native_T_values_verified=1384704,complete_vertex_set_bijection=True,native_temperature_source_kind='native_rank_reassembly',all601_actual_both_rank_T_reread_with_cellProcAddressing=True,native_rank_topology_proofs=topology_proofs,actual_rank_time_metadata_files_verified=time_proofs,accepted_native_rank_reassembly_proof=record['native_rank_reassembly_proof'],actual_execution_completion=record['execution_completion'],actual_all601_U_p_and9336_native_phi_histories_consistent=True,native_physical_gradient_source='Literal retained actual two-rank physical gradT patches mapped by accepted audit; incomplete slow serial field exports are diagnostics only',diagnostic_comparisons=record.get('diagnostic_comparisons',{}),initial_heat_unsolved=True)
    z.close();return finish_bundle(b,g)


def status_notes(statuses,record):
    notes=[line.replace('declared numerical checks passed','declared primary checks passed') for line in common.footer_lines(statuses)];d=record.get('diagnostic_comparisons',{})
    if d and d.get('original_strict_stencil_comparison_passed') is False:
        require('original_strict_stencil_max_relative_error' in d and 'original_strict_stencil_relative_target' in d,'Auxiliary miss requires actual metric/target')
        notes.append(f"OpenFOAM auxiliary strict gradient: FAIL {len(d['original_strict_failed_step_indices'])}/600 states; max {d['original_strict_stencil_max_relative_error']:.5g} > {d['original_strict_stencil_relative_target']:.5g}")
    if d.get('precision_interval_all600_consistent') is False:
        notes.append(f"OpenFOAM auxiliary print-precision interval: FAIL {len(d['precision_interval_inconsistent_step_indices'])}/600 states inconsistent at {d['print_precision_significant_digits']} significant digits")
    return notes


def render_animation(g,bundles,out,notes,indices,name,poster=False):
    fig,artists,clock,axes=layout.temperature_canvas(g,bundles);common.add_footer(fig,notes)
    for text in fig.texts:
        if text.get_text().startswith('Physical time is the native solver time;'):
            text.set_text(f'Actual native states {indices[0]}-{indices[-1]} | playback {len(indices)} frames at 50 ms each. No temporal interpolation.')
        if not poster and text.get_text()=='Heating a channel on a frozen computed flow':text.set_text('Early heating | actual first 0.5 s')
    # A small stationary overlay uses actual frozen cell-centre U once only.
    xc=(g['x_edges'][:-1]+g['x_edges'][1:])/2;yc=(g['y_edges'][:-1]+g['y_edges'][1:])/2
    for key,b in bundles.items():
        U=np.empty((24,96,2));U[g['iy'],g['ix']]=b.frozen_U_m_s
        axes[key].streamplot(xc,yc,U[:,:,0],U[:,:,1],color='#b8cbd2',linewidth=.33,density=.6,arrowsize=.45,zorder=3)
    images=[];actual_times=[];map_hashes={key:[] for key in bundles};native_deltas={key:[] for key in bundles};map_crops={}
    for frame,state in enumerate(indices):
        if frame%10==0:common.resource_guard()
        times=[b.times_s[state] for b in bundles.values()];require(max(times)-min(times)<1e-10,'Matched actual frame times differ');actual_times.append(times[0])
        for key,b in bundles.items():artists[key].set_array(b.temperature_K[state]);native_deltas[key].append(float(abs(b.temperature_K[state]-b.temperature_K[indices[max(0,frame-1)]]).max()))
        clock.set_text(f'Actual physical time {times[0]:.3f} s | native thermal step {state}/600 | playback 50 ms/frame')
        if frame==0:common.verify_text_bounds(fig,out,name)
        fig.canvas.draw();rgb=Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy());images.append(rgb)
        for key,ax in axes.items():
            box=ax.get_window_extent();map_crops[key]=(round(box.x0)+2,round(rgb.height-box.y1)+2,round(box.x1)-2,round(rgb.height-box.y0)-2)
    if poster:images[-1].save(out/'heated-channel-native-temperature.png')
    gif=out/(name+'.gif');images[0].save(gif,save_all=True,append_images=images[1:],duration=50,loop=0,optimize=False,disposal=2);plt.close(fig)
    with Image.open(gif) as decoded:
        require(decoded.n_frames==len(indices),'Actual selected GIF frame count differs');durations=[];full_hashes=[]
        for frame in range(decoded.n_frames):
            decoded.seek(frame);durations.append(decoded.info['duration']);rgb=decoded.convert('RGB');full_hashes.append(hashlib.sha256(rgb.tobytes()).hexdigest())
            for key,crop in map_crops.items():map_hashes[key].append(hashlib.sha256(rgb.crop(crop).tobytes()).hexdigest())
        require(durations==[50]*len(indices),'Decoded GIF delays differ')
        if poster:
            sheet=Image.new('RGB',(images[0].width,(images[0].height+30)*3),'#e9eff2');draw=ImageDraw.Draw(sheet)
            for slot,frame in enumerate((0,len(indices)//2,len(indices)-1)):
                decoded.seek(frame);draw.text((16,slot*(images[0].height+30)+8),f'Decoded actual native time {actual_times[frame]:.3f} s | step {indices[frame]}',fill='#203c4b');sheet.paste(decoded.convert('RGB'),(0,slot*(images[0].height+30)+30))
            sheet.save(out/'heated-channel-decoded-GIF-contact-sheet.png')
    result=dict(frames=len(indices),frame_delay_ms=50,playback_duration_s=len(indices)*.05,actual_native_step_indices=[int(v) for v in indices],actual_physical_times_s=actual_times,actual_solver_temperature_evolution=True,solver_flow_evolution=False,temporal_interpolation=False,particles=False,temperature_limits_K=[300.,350.],distinct_decoded_RGB_frames=len(set(full_hashes)),map_region_unique_frames={key:len(set(v)) for key,v in map_hashes.items()},map_uniqueness_source='Decoded GIF RGB map crops',map_crop_pixel_bounds=map_crops,native_maximum_T_change_since_previous_selected_state_K=native_deltas,map_frame_uniqueness_requirement=False,plateau_policy='Actual settled or quantized map repeats retained without fabricated motion',sha256=sha(gif))
    return result


def write_histories(out,g,bundles):
    with gzip.open(out/'actual-native-temperatures.csv.gz','wt',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','actual_native_step','actual_physical_time_s','source_quad_0based','x_m','y_m','native_temperature_K'])
        for name,b in bundles.items():
            for step in range(601):w.writerows([name,step,b.times_s[step],i,*g['centres'][i],b.temperature_K[step,i]] for i in range(2304))
    with (out/'record-selected-native-temperatures.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','actual_native_step','actual_physical_time_s','source_quad_0based','x_m','y_m','native_temperature_K'])
        for name,b in bundles.items():
            for step in (0,1,10,100,600):w.writerows([name,step,b.times_s[step],i,*g['centres'][i],b.temperature_K[step,i]] for i in range(2304))
    with (out/'actual-time-metrics.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','step','actual_time_s','heat_solved','energy_above300_J_per_m','native_outlet_bulk_T_K','T_min_K','T_max_K','storage_W_per_m','native_Euler_relative_defect','reconstructed_PDE_relative_defect','native_terminal_thermal_residual'])
        for name,b in bundles.items():
            for step in range(601):w.writerow([name,step,b.times_s[step],step>0,b.energy_J_per_m[step],b.outlet_bulk_T_K[step],b.temperature_K[step].min(),b.temperature_K[step].max(),*([b.storage[step-1],b.native_Euler_relative_defect[step-1],b.reconstructed_PDE_relative_defect[step-1],b.terminal_native_thermal_residual[step-1]] if step else ['']*4)])
    with (out/'actual-boundary-powers.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','step','actual_time_s','patch','raw_native_total_report_W_if_available','native_reference_corrected_Qin_W_per_m','reconstructed_Qin_W_per_m','reconstructed_minus_native_W_per_m','native_diffusive_Qin_W_per_m_if_available','native_convective_Qin_W_per_m_if_available','reconstructed_diffusive_Qin_W_per_m_if_available','reconstructed_convective_Qin_W_per_m_if_available'])
        for name,b in bundles.items():
            for row in range(600):
                for j,patch in enumerate(BOUNDS):
                    pick=lambda a:'' if a is None else a[row,j]
                    w.writerow([name,row+1,b.solved_times_s[row],patch,pick(b.raw_native_Q_W),b.native_Qin_W_per_m[row,j],b.reconstructed_Qin_W_per_m[row,j],b.reconstructed_Qin_W_per_m[row,j]-b.native_Qin_W_per_m[row,j],pick(b.native_diffusion),pick(b.native_convection),pick(b.reconstructed_diffusion),pick(b.reconstructed_convection)])
    with (out/'actual-energy-defects.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','step','actual_time_s','storage_W_per_m','native_sum_Qin_W_per_m','reconstructed_sum_Qin_W_per_m','normalization_W_per_m','native_Euler_relative_defect','native_target_1e_5','reconstructed_PDE_relative_defect','reconstructed_target_0p005'])
        for name,b in bundles.items():
            for i in range(600):w.writerow([name,i+1,b.solved_times_s[i],b.storage[i],b.native_Qin_W_per_m[i].sum(),b.reconstructed_Qin_W_per_m[i].sum(),b.normalization[i],b.native_Euler_relative_defect[i],b.native_Euler_relative_defect[i]<=1e-5,b.reconstructed_PDE_relative_defect[i],b.reconstructed_PDE_relative_defect[i]<=.005])
    with (out/'actual-native-thermal-residuals.csv').open('x',encoding='utf-8',newline='') as stream:
        w=csv.writer(stream);w.writerow(['solver','step','actual_time_s','initial_linear_residual_if_available','actual_terminal_native_residual','actual_inner_iterations_if_available','own_norm'])
        for name,b in bundles.items():
            for i in range(600):w.writerow([name,i+1,b.solved_times_s[i],'' if b.initial_residual is None else b.initial_residual[i],b.terminal_native_thermal_residual[i],'' if b.inner_iterations is None else b.inner_iterations[i],'Fluent scaled energy' if name=='fluent' else 'OpenFOAM final linear T'])


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for flag in ('geometry','out','fluent-run','fluent-audit','fluent-fields','openfoam-audit','openfoam-fields'):ap.add_argument('--'+flag,type=Path,required=True)
    for name in ('fluent','openfoam'):ap.add_argument('--'+name+'-fields-sha256',required=True)
    args=ap.parse_args();out=args.out.resolve();require(out.drive.upper()=='E:' and not out.exists(),'Fresh external E: export required')
    for key,module in (('reader',fr),('layout',layout),('cavity',common),('native',common.native)):require(sha(Path(module.__file__))==FROZEN[key],'Frozen helper differs: '+key)
    resources_before=common.resource_guard(render_reserve=True);g=geometry(args.geometry)
    paths={'fluent':args.fluent_fields.resolve(),'openfoam':args.openfoam_fields.resolve()};audits={'fluent':args.fluent_audit.resolve(),'openfoam':args.openfoam_audit.resolve()};digests={'fluent':args.fluent_fields_sha256,'openfoam':args.openfoam_fields_sha256}
    records={};statuses={}
    for name in paths:records[name],statuses[name]=read_audit(name,audits[name],paths[name],digests[name])
    from ansys.fluent.core.offline import _lispy
    inputs=[Path(__file__).resolve(),Path(fr.__file__),Path(layout.__file__),Path(common.__file__),Path(common.native.__file__),Path(sys.modules['analyse_nozzle_fluent_v005'].__file__),Path(_lispy.__file__),args.geometry.resolve(),args.geometry.parent/'setup.json',*paths.values()]
    for s in statuses.values():inputs += [Path(s['path']),Path(s['source_manifest_path']),*(Path(p) for p in s['native_source_hashes'])]
    inputs=sorted(set(inputs),key=str);before={str(p):sha(p) for p in inputs}
    bundles={'fluent':fluent_bundle(args.fluent_run.resolve(),paths['fluent'],statuses['fluent'],g),'openfoam':openfoam_bundle(paths['openfoam'],records['openfoam'],statuses['openfoam'],g)}
    require(all(s['nonboolean_check_normalization'] is None or s['nonboolean_check_normalization']['independent_native_proof_completed'] is True for s in statuses.values()),'Pending raw audit serialization proof blocks rendering')
    require(np.max(abs(bundles['fluent'].times_s-bundles['openfoam'].times_s))<1e-10,'Actual paired times differ')
    notes=status_notes(statuses,records['openfoam']);plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.linewidth':.7});out.mkdir(parents=True)
    animation=render_animation(g,bundles,out,notes,np.arange(0,601,10),'heated-channel-actual-temperature',poster=True)
    early=render_animation(g,bundles,out,notes,np.arange(101),'heated-channel-early-actual-temperature',poster=False)
    fig=layout.outlet_and_energy_canvas(bundles);common.add_footer(fig,notes);common.save(fig,out,'heated-channel-outlet-and-energy')
    fig=layout.power_and_defect_canvas(bundles);common.add_footer(fig,notes);common.save(fig,out,'heated-channel-boundary-power-and-balance')
    write_histories(out,g,bundles);require(before=={str(p):sha(p) for p in inputs},'Native/audit/source input changed during export')
    diff=bundles['fluent'].temperature_K-bundles['openfoam'].temperature_K
    checks=dict(model='Actual transient heated channel on separately frozen native computed flows',actual_native_temperatures=True,actual_states_per_code=601,native_temperature_values_verified_per_code=1384704,paired_native_temperature_values_verified=2769408,native_mapping={k:b.meta for k,b in bundles.items()},numerical_audits=statuses,auxiliary_diagnostics=records['openfoam'].get('diagnostic_comparisons',{}),input_hashes_stable=True,animation=animation,early_animation=early,initial_heat_unsolved_and_not_padded=True,temperature_scale_K=[300.,350.],true_channel_aspect=10.,solver_temperature_evolution=True,solver_flow_evolution=False,particles=False,temporal_interpolation=False,global_energy_conservation_only=True,all_cell_energy_flux_conservation_claim=False,physical_validation=False,mesh_independence=False,time_step_independence=False,
        source_bound_audited_NPZ={k:dict(path=str(paths[k]),sha256=digests[k]) for k in paths},actual_T_ranges_K={k:[float(b.temperature_K.min()),float(b.temperature_K.max())] for k,b in bundles.items()},final_outlet_bulk_T_K={k:float(b.outlet_bulk_T_K[-1]) for k,b in bundles.items()},final_energy_J_per_m={k:float(b.energy_J_per_m[-1]) for k,b in bundles.items()},maximum_native_Euler_relative_defect={k:float(b.native_Euler_relative_defect.max()) for k,b in bundles.items()},maximum_reconstructed_PDE_relative_defect={k:float(b.reconstructed_PDE_relative_defect.max()) for k,b in bundles.items()},maximum_terminal_native_thermal_residual={k:float(b.terminal_native_thermal_residual.max()) for k,b in bundles.items()},maximum_frozen_field_change_SI={k:b.frozen_changes.max(axis=0).tolist() for k,b in bundles.items()},cross_code_temperature_difference_descriptive_K=dict(max_abs_by_actual_step=np.max(abs(diff),axis=1).tolist(),area_weighted_RMS_by_actual_step=np.sqrt(np.sum(g['areas'][None,:]*diff*diff,axis=1)/g['areas'].sum()).tolist(),accuracy_gate=False),record_data_policy='Full time-cell CSV external E only; compact all600 solved-step histories plus root-selected native maps for later reviewed record.',resources=dict(before=resources_before,after=common.resource_guard()))
    (out/'checks.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for path,name in ((Path(__file__),'executed-plotter.py'),(Path(fr.__file__),'reused-Fluent-production-reader.py'),(Path(layout.__file__),'reused-figure-layout.py'),(Path(common.__file__),'reused-cavity-helper.py'),(Path(common.native.__file__),'reused-native-readers.py'),(Path(sys.modules['analyse_nozzle_fluent_v005'].__file__),'reused-native-mapping-helper.py')):(out/name).write_bytes(path.read_bytes())
    (out/'input-hash-manifest.json').write_text(json.dumps(dict(command_arguments=sys.argv[1:],inputs=[dict(path=str(p),sha256=before[str(p)],bytes=p.stat().st_size) for p in inputs]),indent=2)+'\n',encoding='utf-8')
    caption='Actual matched 601 native temperature states per code. Main 61 states 0,10,..600 cover 0..3 s; supplementary 101 native states 0..100 cover 0..0.5 s. Playback 50 ms/frame is distinct from physical time. No interpolated states or particles; true 10:1 maps and fixed 300..350 K scale. Actual plateaus/quantized equal map crops remain unchanged.\n\nNative U,p and all face fluxes stay frozen in each solver. All 2,769,408 paired native temperatures are independently reread and compared to accepted audited arrays. OpenFOAM uses both literal native rank T files and actual cellProcAddressing for every state, with disjoint complete cell coverage and actual vertex/source geometry proof. Both rank time/dt/global step metadata are reread at all 600 solved times. Native thermal execution and slow serial reconstruction are separate statuses; partial serial field exports remain diagnostic crosschecks, without a CFD replay. Accepted raw completion and available closed serial-state crosschecks remain in checks.json. Initial heat is unsolved/unavailable; power and defect plots use only 600 solved steps. Fluent actual h300 reference correction is applied to integrated native total-W reports without double-counting convection or multiplying by area. OpenFOAM direct gradient uses retained native rank physical patches; reconstructed serialized gradT is diagnostic only. Native global Euler and reconstructed PDE flux defects are distinct and use their unchanged targets 1e-5 and 0.005. Native and auxiliary false outcomes remain explicit. No physical validation, mesh/time independence or unobserved all-cell energy claim.\n\nOpenFOAM execution/postprocessing status: '+json.dumps(records['openfoam']['execution_completion'],sort_keys=True)+'\n\n'+'\n'.join(notes)+'\n'
    (out/'CAPTIONS.md').write_text(caption,encoding='utf-8');files=sorted(p for p in out.iterdir() if p.is_file());(out/'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name,sha256=sha(p),bytes=p.stat().st_size) for p in files],indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(out),animation=animation,early_animation=early,native_numerical_failures={k:s['failed_checks'] for k,s in statuses.items()}),indent=2))

if __name__=='__main__':main()
