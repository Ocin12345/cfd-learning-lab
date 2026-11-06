"""Read-only audit of ONE actual heated-channel frozen-flow thermal pilot; never launches a solver."""
from pathlib import Path
import argparse,csv,json,re,sys,math
import numpy as np
import prepare_openfoam_heated_channel_pilot_v003 as prep
from analyse_heated_channel_flow_openfoam_v002 import Evidence,field,parse_value,boundary_blocks,cell_mapping,check_wall_fields
require=prep.require
sha=prep.sha
PATCHES=prep.PATCHES
CELLS=prep.CELLS
EMPTY=prep.EMPTY

def native_T_field(raw,components,dimensions):
    # Count/dimension errors are structural; genuine nonfinite numerical values remain data/check failures.
    text=prep.clean(raw);require(re.search(r'\bformat\s+ascii\s*;',text),'Native ASCII thermal field required')
    dm=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',text);require(dm is not None and list(map(int,dm[1].split()))==dimensions,'Thermal dimensions differ')
    return thermal_value(text,CELLS,components,'internalField')

def thermal_value(text,count,components,label='value'):
    m=re.search(r'\b'+label+r'\s+uniform\s+([^;]+);',text)
    if m:
        row=np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ');require(row.size==components,'Native thermal uniform component count differs')
        return np.tile(row,(count,1)) if components>1 else np.full(count,row[0])
    kind='vector' if components==3 else 'scalar';m=re.search(r'\b'+label+r'\s+nonuniform\s+List<'+kind+r'>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    require(m is not None and int(m[1])==count,'Native thermal count differs')
    v=np.fromstring(m[2].replace('(',' ').replace(')',' '),sep=' ');require(v.size==count*components,'Native thermal value count differs')
    return v.reshape(count,components) if components>1 else v

def thermal_blocks(raw):
    text=prep.clean(raw);require('boundaryField' in text,'Missing native thermal boundaries')
    b=dict(re.findall(r'\b(\w+)\s*\{([^{}]*)\}',text.split('boundaryField',1)[1]));require(set(b)==set(PATCHES)|{EMPTY},'Thermal patch set differs');return b

def face_geometry(geometry,centres,patches):
    points=np.asarray(geometry['points']);faces=np.asarray(geometry['faces']);owner=np.asarray(geometry['owner']);result={}
    for name in PATCHES:
        q=patches[name];ids=np.arange(q['start'],q['start']+q['count']);vertices=points[faces[ids]]
        Sf=.5*np.sum(np.cross(vertices,np.roll(vertices,-1,axis=1)),axis=1);A=np.linalg.norm(Sf,axis=1);n=Sf/A[:,None];Cf=vertices.mean(axis=1);own=owner[ids]
        distance=np.sum((Cf-centres[own])*n,axis=1);require(np.all(A>0) and np.all(distance>0),'Actual native boundary area/distance not positive')
        require(np.max(np.linalg.norm((Cf-centres[own])-distance[:,None]*n,axis=1))<=1e-10,'Boundary mesh not orthogonal to owner centres')
        result[name]=dict(ids=ids,owner=own,Sf=Sf,area=A,normal=n,centre=Cf,distance=distance)
    return result

def phi_values(raw,patches):
    internal=field(raw,prep.INTERNAL_FACES,1,[0,3,-1,0,0,0,0]);b=boundary_blocks(raw);allphi=np.zeros(prep.TOTAL_FACES);allphi[:prep.INTERNAL_FACES]=internal;by_patch={}
    for n in PATCHES:
        q=patches[n];v=parse_value(b[n],q['count'],1,'value');allphi[q['start']:q['start']+q['count']]=v;by_patch[n]=v
    require(re.search(r'\btype\s+empty\s*;',b[EMPTY]),'Native phi empty constraint differs')
    return allphi,by_patch

def flow_boundary_arrays(u_raw,p_raw,U,P,patches,geometry):
    check_wall_fields(u_raw,p_raw,patches);ub=boundary_blocks(u_raw);pb=boundary_blocks(p_raw);owner=np.asarray(geometry['owner']);us=[];ps=[]
    for n in PATCHES:
        q=patches[n];own=owner[q['start']:q['start']+q['count']]
        if re.search(r'\bvalue\s+',ub[n]):uv=parse_value(ub[n],q['count'],3,'value')
        elif n in ('top','bottom'):uv=np.zeros((q['count'],3))
        else:require(n=='outlet','Only native outlet zeroGradient U may omit value');uv=U[own]
        if re.search(r'\bvalue\s+',pb[n]):pv=parse_value(pb[n],q['count'],1,'value')
        else:require(re.search(r'\btype\s+zeroGradient\s*;',pb[n]),'Missing p value without zeroGradient');pv=P[own]
        us.append(uv);ps.append(pv)
    return np.vstack(us),np.concatenate(ps)

def thermal_boundary(Traw,Graw,T,phi,geometry):
    tb=thermal_blocks(Traw);gb=thermal_blocks(Graw);vectors=[];normal_grad=[];boundary_T=[];recon_sn=[];upwind=[];native_Q=[];recon_Q=[];qdiff=[];qconv=[];native_sn={};stencil_sn={};zone_report={}
    for name in PATCHES:
        g=geometry[name];count=len(g['ids']);own=g['owner'];block=tb[name]
        if name=='outlet':
            require(re.search(r'\btype\s+zeroGradient\s*;',block),'Native outlet T not zeroGradient')
            Tb=thermal_value(block,count,1) if re.search(r'\bvalue\s+',block) else T[own]
            recon=np.zeros(count)
        else:
            require(re.search(r'\btype\s+fixedValue\s*;',block),'Native inlet/wall T not fixedValue');Tb=thermal_value(block,count,1)
            expected=300 if name=='inlet' else 350;require(np.all(Tb==expected),'Prescribed native T boundary differs')
            recon=(Tb-T[own])/g['distance']
        G=thermal_value(gb[name],count,3);sn=np.sum(G*g['normal'],axis=1)
        # On this admitted flow inlet phi<0, outlet phi>0, walls0. Boundary field values
        # therefore give the native Gauss-upwind convection; zeroGradient outlet equals owner.
        Tup=Tb if name!='outlet' else T[own]
        reconstructed_upwind=np.where(phi[name]>=0,T[own],Tb)
        Qd=10*sn*g['area']/.01;Qc=-1000*phi[name]*(Tup-300)/.01;Q=Qd+Qc
        Qr=10*recon*g['area']/.01-1000*phi[name]*(reconstructed_upwind-300)/.01
        native_Q.append(float(Q.sum()));recon_Q.append(float(Qr.sum()));qdiff.append(float(Qd.sum()));qconv.append(float(Qc.sum()))
        vectors.append(G);normal_grad.append(sn);boundary_T.append(Tb);recon_sn.append(recon);upwind.append(Tup);native_sn[name]=sn;stencil_sn[name]=recon
        zone_report[name]=dict(actual_faces=count,diffusive_Qin_W_per_m=float(Qd.sum()),convective_reference300_Qin_W_per_m=float(Qc.sum()),total_Qin_W_per_m=float(Q.sum()),reconstructed_Qin_W_per_m=float(Qr.sum()),
          snGradT_native_min_K_m=float(np.min(sn)),snGradT_native_max_K_m=float(np.max(sn)),native_stencil_snGrad_max_abs_K_m=float(np.max(np.abs(sn-recon))))
    arrays=dict(boundary_gradT=np.vstack(vectors),boundary_snGradT=np.concatenate(normal_grad),boundary_T=np.concatenate(boundary_T),reconstructed_boundary_snGradT=np.concatenate(recon_sn),upwind_boundary_T=np.concatenate(upwind),
      native_Qin_W_per_m=np.array(native_Q),reconstructed_Qin_W_per_m=np.array(recon_Q),native_diffusive_Qin_W_per_m=np.array(qdiff),native_convective_Qin_W_per_m=np.array(qconv))
    return arrays,zone_report

def finite_json(value):
    if isinstance(value,dict):return {k:finite_json(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [finite_json(v) for v in value]
    if isinstance(value,(float,np.floating)):return float(value) if math.isfinite(float(value)) else None
    if isinstance(value,np.integer):return int(value)
    if isinstance(value,np.bool_):return bool(value)
    return value

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    case=prep.checked_case(args.run,'pilot');out=args.out.resolve();require(out.drive.upper()=='E:' and not out.exists() and case not in out.parents,'Fresh external E: output required')
    e=Evidence();checks={};e.touch(Path(__file__));e.touch(Path(prep.__file__));e.touch(Path(__file__).with_name('analyse_heated_channel_flow_openfoam_v002.py'))
    status=e.read_json(case/'run-status.json');require(status.get('stage')=='finished_thermal_pilot_execution' and status.get('actual_thermal_steps')==1 and status.get('physical_time_s')==.005 and status.get('thermal_production_admitted') is False,'Only genuinely completed single-step pilot can be audited')
    declaration=e.read_json(case/'declared-comparison.json');require(declaration['thermal_end_time_s']==.005 and declaration['actual_steps_authorized']==1 and declaration['production_admitted'] is False,'Pilot cap/admission changed')
    for n,expected in prep.dictionaries('pilot').items():
        if n!='0/T':require(e.text(case/n)==expected,'Actual thermal dictionary differs '+n)
    prepared=e.read_json(case/'provenance/prepared-input-sha256.json')
    for n,digest in prepared.items():
        e.touch(case/n)
        if n!='0/T':require(sha(case/n)==digest,'Prepared immutable pilot source changed '+n)
    require(sha(e.touch(case/'provenance/TARGETS.md'))==prep.TARGET_SHA==declaration['target_document_sha256'],'Frozen targets changed')
    for folder in ('constant','system'):
        for name in ('fvModels','fvConstraints','fvOptions'):
            candidate=case/folder/name
            if candidate.exists():
                body=re.sub(r'FoamFile\s*\{[^{}]*\}','',prep.clean(e.text(candidate)))
                require(not body.strip(),'Unexpected external scalar source/constraint models configured')
    checks['actual_pilot_dictionaries_match_one_step_Euler_upwind_D_and_frozen_targets']=True
    ref=e.read_json(case/'provenance/frozen-flow-references.json');flow=prep.canonical_native_source_path(ref['actual_flow_case']);require(flow==prep.run_root()/'openfoam-flow-v002','Unexpected actual frozen flow source')
    flow_checks=e.read_json(case/'provenance/flow-checks.json');require(flow_checks.get('all_checks_passed') is True,'Required flow audit not passed')
    require(flow_checks['analysis_source_sha256']==sha(Path(__file__).with_name('analyse_heated_channel_flow_openfoam_v002.py')) and flow_checks['preparation_source_sha256']==prep.FLOW_PREPARATION_SHA,'Current imported native parsers differ from audited frozen flow sources')
    require(sha(e.touch(prep.private_root()/'openfoam-flow-checks-v001/checks.json'))==sha(case/'provenance/flow-checks.json'),'Selected actual flow audit differs')
    require(sha(e.touch(prep.private_root()/'openfoam-flow-checks-v001/checks.json'))==prep.FLOW_CHECKS_SHA,'Pinned original flow audit changed')
    require(sha(e.touch(case/'provenance/prepare_openfoam_heated_channel_v002.py'))==prep.FLOW_PREPARATION_SHA,'Original flow preparation source changed')
    for n,v in ref['fields'].items():
        require(sha(e.touch(prep.canonical_native_source_path(v['flow_path'])))==v['flow_sha256']==sha(e.touch(case/v['pilot_path'])),'Actual frozen source/copy bytes changed '+n)
    for n in ('U','p','phi'):require(sha(e.touch(case/'provenance/baseline-native-inputs'/n))==ref['fields'][n]['flow_sha256'],'Immutable baseline native input bytes changed')
    checks['actual_audited_flow2000_and_both_frozen_reference_copies_byte_bound_unchanged']=True
    poly=prep.native_poly(case);flowpoly=prep.native_poly(flow)
    for n in ('points','faces','owner','neighbour','boundary'):require(sha(e.touch(poly/n))==sha(e.touch(flowpoly/n)),'Pilot mesh differs from audited flow')
    patches=prep.native_patches(poly);geometry=prep.native_geometry(poly);centres=field(e.text(case/'2000/C'),CELLS,3,[0,1,0,0,0,0,0]);volumes=field(e.text(case/'2000/Vc'),CELLS,1,[0,3,0,0,0,0,0])
    require(np.all(volumes>0) and np.max(np.abs(volumes-np.array(geometry['volumes'])))<=1e-14,'Actual native Vc differs from oriented mesh')
    mapping,mapping_record=cell_mapping(centres,e.touch(case/'provenance/geometry.npz'),volumes,geometry['points']);FG=face_geometry(geometry,centres,patches)
    checks['actual_native_mesh_C_Vc_source_cells_vertices_and_boundary_orientation_match_flow']=True
    log=e.text(case/'solver.log');require(re.search(r'^End\s*$',log,re.M) and re.search(r'\bnProcs\s*:\s*2\b',log),'Actual completed two-rank log required')
    actual_times=[float(x) for x in re.findall(r'^Time = ('+prep.FLOAT+r')s?\s*$',log,re.M)];require(actual_times==[.005],'Exactly one real physical thermal step required')
    require(not re.search(r'Solving for (?:Ux|Uy|Uz|p),',log),'Frozen-flow pilot unexpectedly solved momentum/pressure')
    residuals=re.findall(r'Solving for T, Initial residual = ('+prep.FLOAT+r'), Final residual = ('+prep.FLOAT+r'), No Iterations (\d+)',log);require(len(residuals)==1,'Exactly one native scalar matrix solve required')
    initial_residual,final_residual=float(residuals[0][0]),float(residuals[0][1]);iterations=int(residuals[0][2])
    for rank in (0,1):
        for n in ('T','gradT','U','p','phi'):e.touch(case/('processor'+str(rank))/'0.005'/n)
    time_raw=e.text(case/'0.005/uniform/time');time_meta={n:float(re.search(r'\b'+n+r'\s+('+prep.FLOAT+r')\s*;',time_raw)[1]) for n in ('value','deltaT','index')}
    require(time_meta['value']==.005 and time_meta['deltaT']==.005 and time_meta['index']==1,'Native physical time/dt/index differs')
    checks['actual_two_rank_one_step_physical_time_dt_index_and_no_flow_equations_proved']=True
    checks['final_linear_T_residual_le_1e_8']=final_residual<=1e-8
    Uref_raw=e.text(case/'2000/U');Pref_raw=e.text(case/'2000/p');phiref_raw=e.text(case/'2000/phi');Uref=field(Uref_raw,CELLS,3,[0,1,-1,0,0,0,0]);Pref=field(Pref_raw,CELLS,1,[0,2,-2,0,0,0,0]);phiref,phiref_by=phi_values(phiref_raw,patches);Ubref,Pbref=flow_boundary_arrays(Uref_raw,Pref_raw,Uref,Pref,patches,geometry)
    Tstates=[];Ustates=[];Pstates=[];phistates=[];boundary={};zone_reports={};freeze=[];state_stats={};temp_finite=True
    for label in ('0','0.005'):
        ur=e.text(case/label/'U');pr=e.text(case/label/'p');phr=e.text(case/label/'phi');tr=e.text(case/label/'T');gr=e.text(case/label/'gradT')
        U=field(ur,CELLS,3,[0,1,-1,0,0,0,0]);P=field(pr,CELLS,1,[0,2,-2,0,0,0,0]);Phi,Phi_by=phi_values(phr,patches);Ub,Pb=flow_boundary_arrays(ur,pr,U,P,patches,geometry)
        T=native_T_field(tr,1,[0,0,0,1,0,0,0]);G=native_T_field(gr,3,[0,-1,0,1,0,0,0]);b,zr=thermal_boundary(tr,gr,T,Phi_by,FG)
        Tstates.append(T);Ustates.append(U);Pstates.append(P);phistates.append(Phi);boundary[label]=b;zone_reports[label]=zr
        finite=bool(np.isfinite(T).all() and np.isfinite(G).all() and all(np.isfinite(v).all() for v in b.values()));temp_finite=temp_finite and finite
        freeze.append([float(max(np.max(np.abs(U-Uref)),np.max(np.abs(Ub-Ubref)))),float(max(np.max(np.abs(P-Pref)),np.max(np.abs(Pb-Pbref)))),float(np.max(np.abs(Phi-phiref)))])
        state_stats[label]=dict(T_min_K=float(np.min(T)),T_max_K=float(np.max(T)),all_T_and_gradT_finite=finite)
    Tstates=np.asarray(Tstates);freeze=np.asarray(freeze)
    require(np.isfinite(Tstates[0]).all() and np.all(Tstates[0]==300),'Actual native initial thermal field must be300K')
    checks['all_actual_initial_and_solved_native_T_gradT_finite']=temp_finite
    checks['all_actual_T_within_300_350_K_plus_0_01_tolerance']=bool(temp_finite and np.all((Tstates>=299.99)&(Tstates<=350.01)))
    checks['actual_first_step_T_changes_from_initial']=bool(temp_finite and np.max(np.abs(Tstates[1]-Tstates[0]))>0)
    checks['all_actual_runtime_internal_and_boundary_U_p_phi_frozen_le_1e_12']=bool(np.max(freeze)<=1e-12)
    checks['actual_inlet_inward_outlet_outward_no_backflow_and_wall_zero_flux']=bool(np.all(phiref_by['inlet']<0) and np.all(phiref_by['outlet']>=0) and phiref_by['outlet'].sum()>0 and sum(np.abs(phiref_by[n]).sum() for n in ('top','bottom'))<=1e-14)
    energy=1000*np.sum(volumes[None,:]*(Tstates-300),axis=1)/.01;storage=(energy[1]-energy[0])/.005
    native_Q=boundary['0.005']['native_Qin_W_per_m'];recon_Q=boundary['0.005']['reconstructed_Qin_W_per_m'];denom=max(float(np.abs(native_Q).sum()),abs(float(storage)),5000.)
    native_defect=abs(float(storage-native_Q.sum()))/denom;PDE_defect=abs(float(storage-recon_Q.sum()))/denom;power_difference=recon_Q-native_Q
    sn_native=boundary['0.005']['boundary_snGradT'];sn_recon=boundary['0.005']['reconstructed_boundary_snGradT'];sn_error=float(np.max(np.abs(sn_native-sn_recon)/np.maximum(1,np.abs(sn_native))))
    checks['native_boundary_snGradT_matches_selected_orthogonal_stencil_rel_le_1e_10']=bool(temp_finite and sn_error<=1e-10)
    checks['global_native_Euler_relative_energy_defect_le_1e_5']=bool(temp_finite and native_defect<=1e-5)
    checks['separate_reconstructed_PDE_relative_energy_defect_le_0_005']=bool(temp_finite and PDE_defect<=.005)
    gauss_source=e.text(case/'provenance/installed-sources/src/finiteVolume/finiteVolume/gradSchemes/gaussGrad/gaussGrad.C');require('vsf.boundaryField()[patchi].snGrad()' in gauss_source and '(n & gGradbf[patchi])' in gauss_source,'Installed native boundary snGrad source proof differs')
    installed=[line.split(None,1) for line in e.text(case/'provenance/installed-source-sha256.txt').splitlines() if line.strip()];require(len(installed)==len(prep.INSTALLED_SOURCES),'Installed source count differs')
    for digest,path in installed:require(sha(e.touch(case/'provenance/installed-sources'/path[len('/opt/openfoam14/'):]))==digest,'Installed source copy hash differs')
    matrix_source=e.text(case/'provenance/installed-fvMatrixSolve.C');matrix_row=e.text(case/'provenance/installed-fvMatrixSolve-sha256.txt').split()
    require(len(matrix_row)==2 and matrix_row[1]=='/opt/openfoam14/src/finiteVolume/fvMatrices/fvMatrix/fvMatrixSolve.C' and sha(case/'provenance/installed-fvMatrixSolve.C')==matrix_row[0],'Actual installed final-solver source hash differs')
    require('word(name + \"Final\")' in matrix_source and 'solutionControl::finalIteration' in matrix_source,'Actual installed TFinal selection proof differs')
    checks['installed_functions_scalarTransport_grad_and_field_IO_source_hashes_bound']=True
    auth=e.read_json(case/'provenance/root-run-authorization.json');require(auth.get('approved') is True and auth.get('approved_by')=='root' and auth.get('action')=='run_heated_channel_openfoam_pilot_v003' and auth.get('allowed_mode')=='ThermalPilot','Separate root pilot authorization missing')
    for n,key in [('prepare_openfoam_heated_channel_pilot_v003.py','preparation_sha256'),('run-openfoam-thermal-pilot-v003.sh','wrapper_sha256'),('Launch-OpenFOAM-ThermalPilot-v003.ps1','launcher_sha256')]:require(sha(e.touch(case/'provenance'/n))==auth[key],'Reviewed execution source hash differs '+n)
    require(sha(e.touch(case/'provenance/run-openfoam-executed.sh'))==auth['wrapper_sha256'],'Executed wrapper differs')
    require(sha(Path(prep.__file__))==auth['preparation_sha256'],'Current imported preparer differs from reviewed execution source')
    checks['separate_root_admission_and_all_frozen_execution_source_hashes_match']=True
    source=e.finish();checks['all_native_source_inputs_unchanged_during_read_only_analysis']=True;checks={n:bool(v) for n,v in checks.items()}
    physical_ids=np.concatenate([FG[n]['ids'] for n in PATCHES]);owner=np.asarray(geometry['owner']);outowner=FG['outlet']['owner'];weights=phiref_by['outlet'];positive=weights>0
    bulk=np.array([float(np.sum(weights[positive]*T[outowner][positive])/weights[positive].sum()) for T in Tstates]);reference_mass_correction=-1000*300*sum(v.sum() for v in phiref_by.values())/.01
    record=dict(phase='one_step_thermal_pilot',native_case=str(case),actual_settings=declaration,physical_times_s=[0,.005],native_time_metadata=time_meta,actual_steps=1,mpi_ranks=2,
      flow_source=ref,actual_cell_mapping=mapping_record,actual_native_patches=patches,native_T_states=state_stats,
      initial_linear_T_residual=initial_residual,final_linear_T_residual=final_residual,linear_T_iterations=iterations,
      numeric_undefined_reason=None if temp_finite else 'Actual native T or gradT contains nonfinite values; derived thermal metrics undefined and numerical checks false',
      frozen_fields=dict(per_state_times_s=[0,.005],max_abs_change_U_m_s=freeze[:,0].tolist(),max_abs_change_kinematic_p=freeze[:,1].tolist(),max_abs_change_all_native_face_phi_m3_s=freeze[:,2].tolist(),internal_and_physical_boundary_values_checked=True,empty_phi_zero_by_exact_native_empty_constraint=True),
      energy=dict(reference_T_K=300,depth_m=.01,energy_above300_J_per_m=energy.tolist(),Euler_storage_W_per_m=float(storage),boundary_order=list(PATCHES),native_boundary_Qin_W_per_m=native_Q.tolist(),
        native_diffusive_Qin_W_per_m=boundary['0.005']['native_diffusive_Qin_W_per_m'].tolist(),native_convective_reference300_Qin_W_per_m=boundary['0.005']['native_convective_Qin_W_per_m'].tolist(),
        denominator_W_per_m=denom,native_Euler_relative_defect=native_defect,reconstructed_PDE_relative_defect=PDE_defect,reconstructed_boundary_Qin_W_per_m=recon_Q.tolist(),reconstructed_minus_native_Qin_W_per_m=power_difference.tolist(),
        reference300_net_mass_flux_correction_W_per_m=float(reference_mass_correction),absolute_T_discrete_Euler_defect_with_reference_mass_correction=abs(float(storage-native_Q.sum()-reference_mass_correction))/denom,
        normal_gradient_native_stencil_max_relative_error=sn_error,all_cell_energy_balance_claimed=False),
      native_boundary_diagnostics_all_states=zone_reports,outlet_outward_flux_weighted_bulk_T_K=bulk.tolist(),
      checks=checks,check_results=checks,failed_checks=[n for n,v in checks.items() if not v],all_checks_passed=all(checks.values()),all_declared_checks_passed=all(checks.values()),
      analysis_source_sha256=sha(Path(__file__)),preparation_source_sha256=sha(Path(prep.__file__)),physical_validation=False,production_admitted=False,completed_600_step_transient=False,
      interpretation='One actual thermal pilot only; direct native gradT and native phi with upwind boundary values reconstruct the discrete global energy flux. Separate orthogonal PDE reconstruction retained; no thermo wallHeatFlux assumed.')
    out.mkdir(parents=True)
    arrays=dict(time_s=np.array([0,.005]),solved_time_s=np.array([.005]),T_history=Tstates,T=Tstates[-1],initial_T=Tstates[0],centres=centres,volumes=volumes,areas=volumes/.01,grid_shape=np.array([24,96]),**mapping,
      U_history=np.asarray(Ustates),p_history=np.asarray(Pstates),native_phi_history=np.asarray(phistates),frozen_reference_U=Uref,frozen_reference_p=Pref,frozen_reference_native_phi=phiref,
      physical_boundary_face_indices=physical_ids,boundary_names=np.array(PATCHES),boundary_offsets=np.array([0,24,48,144,240]),
      native_boundary_gradT_history=np.stack([boundary[t]['boundary_gradT'] for t in ('0','0.005')]),native_boundary_snGradT_history=np.stack([boundary[t]['boundary_snGradT'] for t in ('0','0.005')]),
      native_boundary_T_history=np.stack([boundary[t]['boundary_T'] for t in ('0','0.005')]),native_upwind_boundary_T_history=np.stack([boundary[t]['upwind_boundary_T'] for t in ('0','0.005')]),
      energy_above300_J_per_m=energy,storage_W_per_m=np.array([storage]),native_Qin_W_per_m=native_Q[None,:],reconstructed_Qin_W_per_m=recon_Q[None,:],native_Euler_relative_defect=np.array([native_defect]),
      reconstructed_PDE_relative_defect=np.array([PDE_defect]),native_diffusive_Qin_W_per_m=boundary['0.005']['native_diffusive_Qin_W_per_m'][None,:],native_convective_Qin_W_per_m=boundary['0.005']['native_convective_Qin_W_per_m'][None,:],
      boundary_power_reconstructed_minus_native_W_per_m=power_difference[None,:],outlet_bulk_T_K=bulk,initial_linear_T_residual=np.array([initial_residual]),final_linear_T_residual=np.array([final_residual]),
      frozen_U_max_abs_change_m_s=freeze[:,0],frozen_p_max_abs_change=freeze[:,1],frozen_native_phi_max_abs_change_m3_s=freeze[:,2],rho_kg_m3=np.array(1.),cp_J_kg_K=np.array(1000.),k_W_m_K=np.array(10.),depth_m=np.array(.01))
    np.savez_compressed(out/'fields.npz',**arrays);record['audited_fields_npz_sha256']=sha(out/'fields.npz')
    with (out/'thermal-step-history.csv').open('x',encoding='ascii',newline='') as f:
        w=csv.writer(f);w.writerow(['time_s','dt_s','E_J_per_m','storage_W_per_m','Qin_inlet_W_per_m','Qin_outlet_W_per_m','Qin_top_W_per_m','Qin_bottom_W_per_m','native_Euler_relative_defect','reconstructed_PDE_relative_defect','initial_linear_T_residual','final_linear_T_residual']);w.writerow([.005,.005,energy[-1],storage,*native_Q,native_defect,PDE_defect,initial_residual,final_residual])
    source['audited_outputs']={n:dict(sha256=sha(out/n),bytes=(out/n).stat().st_size) for n in ('fields.npz','thermal-step-history.csv')};prep.json_new(out/'checks.json',finite_json(record));prep.json_new(out/'source-manifest.json',source)
    for path in (Path(__file__),Path(prep.__file__),Path(__file__).with_name('analyse_heated_channel_flow_openfoam_v002.py')):prep.copied(path,out/('executed-'+path.name))
    print(json.dumps(dict(output=str(out),all_checks_passed=all(checks.values()),failed_checks=record['failed_checks']),indent=2));return 0 if all(checks.values()) else 2
if __name__=='__main__':sys.exit(main())
