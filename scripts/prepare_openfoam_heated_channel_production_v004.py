"""Prepare a separately admitted 599-step continuation of the actual passed one-step thermal pilot; never launches CFD."""
from prepare_openfoam_heated_channel_pilot_v003 import *
import prepare_openfoam_heated_channel_pilot_v003 as pilot
PILOT_PREPARATION_SHA='807c30de1c7c7400ab38d2b4592896db774be4c98ec2b853afd2f9fff2ed90b0'
PILOT_ANALYSIS_SHA='52d494048f6458aa2ccb62d1bb1d1f1bc946e723528f71b4f464292e60820faa'
PILOT_CHECKS_SHA='41ff89b6182b88c28ac4a331437d4303c292b4fbec6173430047242e7d6b9a2a'
PILOT_FIELDS_SHA='568a7a26c27d56f766d1cc2d0a50aa9d4555ab3e66f8b999205940bd9c3e034e'
PILOT_MANIFEST_SHA='7b417f82345b14732889ed161156effd8ae0f65e6b26ce961d807a2db47bd15c'
DEPENDENCIES={'prepare_openfoam_heated_channel_pilot_v003.py':PILOT_PREPARATION_SHA,'analyse_heated_channel_pilot_openfoam_v003.py':PILOT_ANALYSIS_SHA,'prepare_openfoam_heated_channel_v002.py':FLOW_PREPARATION_SHA,'analyse_heated_channel_flow_openfoam_v002.py':FLOW_ANALYSIS_SHA}

def checked_case(path,phase='production'):
    p=Path(path).resolve();require(p==(run_root()/'openfoam-thermal-production-v004').resolve(),'Only fresh dedicated thermal production case allowed');return p

def project_path_suffix(raw):
    # Accept only these two literal absolute project prefixes, never ambiguous E:relative paths.
    text=str(raw).replace('\\','/')
    for prefix in ('E:/CFD-Learning-Lab','/mnt/e/CFD-Learning-Lab'):
        if text==prefix:return ''
        if text.startswith(prefix+'/'):
            suffix=text[len(prefix)+1:];require(not suffix.startswith('/') and ':' not in suffix and chr(0) not in suffix,'Ambiguous project-relative path')
            return suffix
    raise ValueError('Native source path lacks exact Windows/Linux project prefix: '+str(raw))

def canonical_native_source_path(raw):
    root=lab_root().resolve();path=(root/project_path_suffix(raw)).resolve()
    require(path==root or root in path.parents,'Native source path resolves outside exact project root')
    return path

def production_dicts():
    content=pilot.dictionaries('pilot');content.pop('0/T')
    content['system/controlDict']=content['system/controlDict'].replace('startTime 0;','startTime 0.005;').replace('endTime 0.005;','endTime 3;')
    return content

def verify_dependencies():
    for name,digest in DEPENDENCIES.items():require(sha(Path(__file__).with_name(name))==digest,'Pinned production dependency changed '+name)

def copy_tree_files(src,dst,rows):
    for p in sorted(src.rglob('*')):
        if p.is_file():
            q=dst/p.relative_to(src);copied(p,q);rows.append(dict(source=str(p),destination=str(q),sha256=sha(p)))

def verify_inherited(case):
    case=checked_case(case);verify_dependencies();r=json.loads((case/'provenance/inherited-pilot-inputs.json').read_text())
    for row in r['files']:
        require(sha(canonical_native_source_path(row['source']))==row['sha256'],'Original passed pilot input changed '+row['source'])
    refs=json.loads((case/'provenance/frozen-flow-references.json').read_text())
    for n,v in refs['fields'].items():
        require(sha(canonical_native_source_path(v['flow_path']))==v['flow_sha256']==sha(case/'2000'/n),'Original frozen flow/copy changed '+n)
    for n in ('U','p','phi'):require(sha(case/'provenance/baseline-native-inputs'/n)==refs['fields'][n]['flow_sha256'],'Immutable frozen baseline changed '+n)
    for name,text in production_dicts().items():require((case/name).read_text()==text,'Reviewed production dictionary changed '+name)
    return dict(original_pilot_source_files_unchanged=True,original_frozen_flow_and_reference_copies_unchanged=True,production_dictionaries_exact=True)

def prepare_production(case):
    case=checked_case(case);require(not case.exists(),'Fresh thermal production case required');verify_dependencies()
    source=run_root()/'openfoam-thermal-pilot-v003';audit=private_root()/'openfoam-pilot-checks-v004'
    status=json.loads((source/'run-status.json').read_text());require(status.get('stage')=='finished_thermal_pilot_execution' and status.get('actual_thermal_steps')==1,'Exactly one completed pilot required')
    checks=json.loads((audit/'checks.json').read_text());require(checks.get('all_checks_passed') is True and checks.get('actual_steps')==1 and checks.get('completed_600_step_transient') is False,'Passed actual one-step pilot audit required')
    require(sha(audit/'checks.json')==PILOT_CHECKS_SHA and sha(audit/'fields.npz')==PILOT_FIELDS_SHA and sha(audit/'source-manifest.json')==PILOT_MANIFEST_SHA,'Selected corrected pilot audit bytes changed')
    require(checks['analysis_source_sha256']==PILOT_ANALYSIS_SHA and checks['preparation_source_sha256']==PILOT_PREPARATION_SHA,'Pilot source selection differs')
    require(sha(private_root()/'TARGETS.md')==TARGET_SHA,'Root frozen targets changed')
    require(shutil.disk_usage(case.anchor).free>=10*1024**3,'Require10GiB E disk')
    # Every selected original pilot input is re-hashed before continuation admission.
    pm=json.loads((audit/'source-manifest.json').read_text())
    for row in pm['sources']:require(sha(canonical_native_source_path(row['path']))==row['sha256_after'],'Passed pilot source changed '+row['path'])
    scalar=(source/'provenance/installed-sources/src/functionObjects/solvers/scalarTransport/scalarTransport.H').read_text()
    require(re.search(r'bool\s+executeAtStart\(\)\s+const\s*\{\s*return false;',scalar),'Installed scalarTransport startup must not solve a duplicate pilot step')
    case.mkdir(parents=True);inherited=[]
    for name,text in production_dicts().items():write_new(case/name,text)
    for label in ('0','0.005','2000'):copy_tree_files(source/label,case/label,inherited)
    for rank in (0,1):
        for label in ('0','0.005','2000'):copy_tree_files(source/('processor'+str(rank))/label,case/('processor'+str(rank))/label,inherited)
    for label in ('0','0.005'):
        for n in ('T','gradT','U','p','phi'):require((case/label/n).is_file(),'Missing copied actual inherited thermal field')
    require((case/'0.005/uniform/time').read_bytes()==(source/'0.005/uniform/time').read_bytes(),'Actual native restart time/index must be copied unchanged')
    for n in ('setup.json','geometry.npz','TARGETS.md','frozen-flow-references.json','installed-source-sha256.txt','installed-fvMatrixSolve.C','installed-fvMatrixSolve-sha256.txt','prepare_openfoam_heated_channel_v002.py','flow-checks.json','flow-source-manifest.json'):
        copied(source/'provenance'/n,case/'provenance'/n)
    for folder in ('installed-sources','baseline-native-inputs'):copy_tree_files(source/'provenance'/folder,case/'provenance'/folder,inherited)
    for name in ('checks.json','fields.npz','source-manifest.json','thermal-step-history.csv'):copied(audit/name,case/'provenance'/('accepted-pilot-'+name))
    copied(source/'solver.log',case/'provenance/accepted-pilot-solver.log');copied(source/'run-status.json',case/'provenance/accepted-pilot-run-status.json')
    for p in (Path(__file__),private_root()/'run-openfoam-thermal-production-v004.sh',private_root()/'Launch-OpenFOAM-ThermalProduction-v004.ps1',*[Path(__file__).with_name(n) for n in DEPENDENCIES]):copied(p,case/'provenance'/p.name)
    json_new(case/'provenance/inherited-pilot-inputs.json',dict(actual_source_case=str(source),inherited_physical_steps=1,original_pilot_is_authoritative_for_steps_0_and_1=True,files=inherited))
    declaration=dict(phase='thermal_production_continuation',solver='OpenFOAM Foundation14 functions/subSolver incompressibleFluid',frozen_read_label=2000,thermal_start_time_s=.005,thermal_end_time_s=3,deltaT_s=.005,
      inherited_actual_steps=1,additional_actual_steps_authorized=599,total_actual_steps_authorized=600,saved_physical_state_count=601,production_admitted_only_by_separate_root_authorization=True,
      original_pilot_case=str(source),original_pilot_checks_sha256=PILOT_CHECKS_SHA,original_pilot_fields_sha256=PILOT_FIELDS_SHA,original_pilot_manifest_sha256=PILOT_MANIFEST_SHA,
      rho_kg_m3=1,cp_J_kg_K=1000,k_W_m_K=10,alpha_m2_s=.01,depth_m=.01,initial_T_K=300,inlet_T_K=300,wall_T_K=350,
      thermal_advection='Gauss upwind',thermal_ddt='Euler',thermal_laplacian='Gauss linear orthogonal',thermal_nCorrectors=0,thermal_equation_relaxation=1,
      final_linear_T_residual_target=1e-8,T_bounds_tolerance_K=.01,frozen_field_max_abs_change=1e-12,global_Euler_relative_energy_target=1e-5,reconstructed_PDE_boundary_flux_target=.005,
      energy_reference_T_K=300,energy_normalization='per native depth; max(sum(abs(Qin_zone)),abs(storage),5000W/m)',native_gradient_source='Retained rank physical gradT values; native face/cell/pointProcAddressing plus current native T and exact outward geometry',
      reconstructed_serial_gradT_is_diagnostic_only=True,native_frozen_snapshot_source='writeObjects anyWrite actual registered U,p,phi',target_document_sha256=TARGET_SHA,physical_validation=False)
    json_new(case/'declared-comparison.json',declaration)
    # Runtime inherited fields may be written during grad/writeObjects startup. Their original pilot source files remain immutable and authoritative.
    excluded=re.compile(r'^(?:processor[01]/)?(?:0|0\.005)/(?:T|gradT|U|p|phi|uniform/.*)$')
    json_new(case/'provenance/prepared-input-sha256.json',{p.relative_to(case).as_posix():sha(p) for p in sorted(case.rglob('*')) if p.is_file() and not excluded.match(p.relative_to(case).as_posix())})
    verify_inherited(case);json_new(case/'run-status.json',dict(stage='thermal_production_prepared',solver_started=False,inherited_steps=1,additional_steps_authorized=599,total_steps_authorized=600))
    return dict(case=str(case),solver_started=False,inherited_actual_steps=1,additional_actual_steps_authorized=599,total_actual_steps_authorized=600)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True);g=ap.add_mutually_exclusive_group(required=True);g.add_argument('--prepare-production',action='store_true');g.add_argument('--verify-inherited',action='store_true');a=ap.parse_args()
    print(json.dumps(prepare_production(a.case) if a.prepare_production else verify_inherited(a.case),indent=2))
if __name__=='__main__':main()
