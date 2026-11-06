"""SOURCE ONLY: build offline visual adapters from actual root-accepted full audits.
No CFD, native commands, rendering, field creation or acceptance receipts are made.
Generated selection deliberately keeps root_authorizes_offline_render false.
"""
from pathlib import Path
import argparse,datetime,hashlib,json,re,shutil
ROOT=Path('E:/CFD-Learning-Lab');PRIVATE=ROOT/'local/prospective-tracer-v001'
AUDITORS={'fluent':('audit_tracer_tjunction_fluent_production_v002.py','9dff240e8343eaf1a6a97344002983a773641924fcce3bcbf83a5e0cd1dd3aa5'),'openfoam':('audit_openfoam_tracer_production_v005.py','7832c8285d45e7f4dfdc2a40efbff741320b4c265625b568c70c661cf0f01c1e')}
PLOTTER=('plot_tracer_tjunction_native_v004.py','60c1bb85858beedbd6f68c08818bc620a942339ed4640a76d061bf5cabc906d2')
TARGET_SHA='07f5b0745fc71758b9cf5678171ce17fd558f00a01cd773560911920e91318ba'
DT=.025
INSTALL_SRC=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/fluent26.1.0/src')
# Exact read-only installed-header exceptions already pinned by accepted flow qualification.
HEADER_PINS={INSTALL_SRC/'storage/mem.h':'fef5d53f4ea2ae006d2d2609bf80aa5afecc361c0998cdba82a8680a0af5b983',INSTALL_SRC/'storage/storage.h':'936eb947aa98762b31e9d516ffa4aa30ba1f1a782ac1539dc84a76f3e7884269',INSTALL_SRC/'storage/threads.h':'84c5706c64b41004c03bcb2925bc9778abdeb07c6662e0bd8d66706db3a169d9',INSTALL_SRC/'udf/sg_udms.h':'70f76587daec2aeb493ab009f3ed287b4f9d65f75903ab61f7d74bc8e335f6d3',INSTALL_SRC/'udf/udf.h':'4414104220f454bdc0c8cece67a8799a5ee8f8d4d7fc4f86dde55f04fc977a21'}

def require(value,message):
    if not value:raise ValueError(message)
def sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1048576),b''):digest.update(chunk)
    return digest.hexdigest()
def project_path(raw):
    path=Path(raw);require(path.is_absolute(),'Absolute project path required');path=path.resolve()
    require(ROOT.resolve() in path.parents,'Project path required: '+str(path));return path
class Evidence:
    def __init__(self):self.rows={}
    def pin(self,path,expected):
        path=Path(path);require(path.is_absolute(),'Absolute actual evidence path required');path=path.resolve()
        require(ROOT.resolve() in path.parents or HEADER_PINS.get(path)==expected,'Evidence outside E project or exact pinned D headers')
        require(path.is_file(),'Missing actual selected input: '+str(path))
        require(isinstance(expected,str) and re.fullmatch('[0-9a-f]{64}',expected),'Explicit lowercase SHA256 required')
        if path not in self.rows:require(sha(path)==expected,'Selected input changed: '+str(path));self.rows[path]=expected
        require(self.rows[path]==expected,'Conflicting input pins: '+str(path));return path
    def ref(self,path,expected=None):
        path=project_path(path);expected=sha(path) if expected is None else expected;self.pin(path,expected);return dict(path=str(path),sha256=expected)
    def json(self,path,expected=None):
        reference=self.ref(path,expected);return json.loads(Path(reference['path']).read_text(encoding='utf-8-sig'))
    def finish(self):
        result=[]
        for path,expected in self.rows.items():
            actual=sha(path);require(actual==expected,'Input changed during selector build: '+str(path));result.append(dict(path=str(path),sha256_before=expected,sha256_after=actual,bytes=path.stat().st_size))
        return result

def same_file_pin(reference,path,evidence):
    require(isinstance(reference,dict) and set(reference)>={'path','sha256'},'Root receipt must bind exact artifact path/SHA')
    require(project_path(reference['path'])==project_path(path),'Root receipt selected a different artifact')
    return evidence.ref(path,reference['sha256'])
def root_acceptance(solver,args,evidence):
    receipt_path=getattr(args,solver+'_acceptance');pin=getattr(args,solver+'_acceptance_sha256')
    receipt=evidence.json(receipt_path,pin)
    require(receipt.get('reviewer')=='root' and receipt.get('full_audit_accepted') is True,'Actual root full-audit acceptance required')
    require(receipt.get('solver')==solver,'Root receipt solver differs')
    files={key:same_file_pin(receipt[key],getattr(args,solver+'_'+arg),evidence) for key,arg in [('checks','checks'),('source_manifest','manifest'),('fields','fields')]}
    files['acceptance']=evidence.ref(receipt_path,pin)
    return files

def manifest_inputs(solver,manifest,evidence):
    if solver=='fluent':
        require(manifest.get('source_bytes_unchanged_during_read') is True,'Actual Fluent source stability evidence required');rows=manifest['inputs']
        pairs=[(row['path'],row['sha256']) for row in rows]
    else:
        require(manifest.get('inputs_unchanged_during_analysis') is True,'Actual OpenFOAM source stability evidence required');rows=manifest['sources'];pairs=[]
        for row in rows:
            require(row['sha256_before']==row['sha256_after'],'OpenFOAM input changed in full audit');pairs.append((row['path'],row['sha256_after']))
    index={}
    for raw,pin in pairs:
        path=evidence.pin(raw,pin);require(path not in index or index[path]==pin,'Conflicting manifest source pins');index[path]=pin
    return index

def native_ref(raw,manifest_index,evidence,declared=None):
    path=project_path(raw);require(path in manifest_index,'Native path absent from accepted full source manifest: '+str(path));pin=manifest_index[path]
    require(declared is None or declared==pin,'Audit/native declared hash differs');return evidence.ref(path,pin)
def assert_checks(record,flag):
    require(record.get(flag) is True and record['checks'] and all(value is True for value in record['checks'].values()),'Accepted full numerical checks are not all literal true')
def json_reference(files,key_path):return dict(**files['checks'],key_path=key_path)

def fluent_rows(checks,files,index,evidence):
    assert_checks(checks,'all_primary_checks_passed')
    require(checks['solver']=='Fluent261' and checks['actual_solved_step_count']==800 and checks['actual_initial_state_count']==1,'Actual full Fluent count/version differs')
    require(checks['production_targets_sha256']==TARGET_SHA and checks['audited_fields_npz_sha256']==files['fields']['sha256'] and project_path(checks['audited_fields_npz'])==project_path(files['fields']['path']),'Fluent audited NPZ/target binding differs')
    states=checks['states'];paths=checks['native_state_paths'];require(len(states)==len(paths)==801,'Full Fluent actual801-state records required');out=[]
    for step,(state,native) in enumerate(zip(states,paths)):
        require(state['step']==native['index']==step and abs(state['time_s']-step*DT)<=1e-12 and state['time_s']==native['time_s'],'Actual Fluent state/time order differs')
        require(native['source_kind']=='native_fluent_hdf','Actual HDF native source required')
        observed=state['observed_native_scalar_datasets']['SV_UDS_0'];require(observed==native['observed_scalar_datasets']['SV_UDS_0'] and len(observed)==1,'Exactly one actually observed current scalar leaf required; ambiguous chunks need separate root review')
        leaf=observed[0];require(isinstance(leaf,str) and leaf.startswith('results/1/phase-1/cells/SV_UDS_0/'),'Observed scalar leaf is not current cell UDS0')
        for key in ('case','data'):require(project_path(state[key])==project_path(native[key]) and state[key+'_sha256']==native[key+'_sha256'],'Fluent audited path row differs')
        row=dict(step=step,time_s=state['time_s'],native_field_name='SV_UDS_0',native=dict(case=native_ref(native['case'],index,evidence,native['case_sha256']),data=native_ref(native['data'],index,evidence,native['data_sha256']),scalar_dataset=leaf))
        if step:
            require(state['scheduled_top_value']==(1. if step<=80 else 0.) and state['actual_top_boundary_max_abs_error']<=1e-12 and state['actual_left_boundary_max_abs_error']<=1e-12,'Actual accepted Fluent inlet schedule evidence differs')
            row.update(prescribed_top_c=state['scheduled_top_value'],prescribed_left_c=0.,inlet_evidence=dict(path=files['checks']['path'],sha256=files['checks']['sha256'],top_key_path=['states',step,'scheduled_top_value'],top_error_key_path=['states',step,'actual_top_boundary_max_abs_error'],left_zero_error_key_path=['states',step,'actual_left_boundary_max_abs_error']),budget=json_reference(files,['states',step,'RECONSTRUCTED_budget']),residual=json_reference(files,['states',step,'own_new_residual','terminal_scaled_residual']))
            require('RECONSTRUCTED' in state['RECONSTRUCTED_budget']['classification'] and isinstance(state['own_new_residual']['terminal_scaled_residual'],(int,float)),'Actual Fluent budget/residual leaf missing')
        else:row.update(budget=None,residual=None,initial_state='Actual accepted unsolved c0; no solved residual or step budget')
        out.append(row)
    return out
def openfoam_rows(checks,files,index,evidence):
    assert_checks(checks,'all_declared_numerical_checks_passed')
    require(checks['analysis_source_sha256']==AUDITORS['openfoam'][1] and checks['production_targets_sha256']==TARGET_SHA,'OpenFOAM actual auditor/target binding differs')
    require(checks['actual_steps']==800 and checks['actual_state_count']==801 and checks['dt_s']==DT and checks['final_time_s']==20.,'Actual full OpenFOAM schedule differs')
    states=checks['native_state_paths'];budgets=checks['reconstructed_scalar_steps'];residuals=checks['actual_native_residuals']
    require(len(states)==801 and len(budgets)==len(residuals)==800,'Full OpenFOAM actual801/800 records required');out=[]
    source_path=PRIVATE/AUDITORS['openfoam'][0];source_reference=evidence.ref(source_path,AUDITORS['openfoam'][1])
    for step,state in enumerate(states):
        require(state['index']==step and abs(state['time_s']-step*DT)<=1e-12,'Actual OpenFOAM state/time order differs')
        native=dict(field=native_ref(state['global_tracer_path'],index,evidence,state['global_tracer_sha256']),lineage=[])
        require(len(state['rank_tracer_paths'])==len(state['rank_tracer_sha256'])==2,'Two actual native scalar ranks required')
        for raw,pin in zip(state['rank_tracer_paths'],state['rank_tracer_sha256']):native['lineage'].append(native_ref(raw,index,evidence,pin))
        row=dict(step=step,time_s=state['time_s'],native_field_name='tracer',native=native)
        if step:
            clock=state['global_clock'];require(clock['expected_index']==step and clock['index']==step and abs(clock['value']-step*DT)<=1e-12 and abs(clock['deltaT']-DT)<=1e-12,'Actual OpenFOAM clock differs')
            native['time_metadata']=native_ref(clock['path'],index,evidence)
            require(len(state['rank_clock_paths'])==len(state['rank_clock_sha256'])==2,'Both actual native rank clocks required')
            for raw,pin in zip(state['rank_clock_paths'],state['rank_clock_sha256']):native['lineage'].append(native_ref(raw,index,evidence,pin))
            top=state['actual_top_c_min'];require(top==state['actual_top_c_max'] and top==(1. if step<=80 else 0.),'Actual saved OpenFOAM top inlet values differ')
            budget,residual=budgets[step-1],residuals[step-1]
            require(budget['index']==residual['index']==step and budget['time_s']==residual['time_s']==state['time_s'],'Actual OpenFOAM budget/residual alignment differs')
            require(budget['actual_native_scalar_binding'] is True and 'RECONSTRUCTED' in budget['fv_helper_return']['classification'],'Actual OpenFOAM scalar binding/classification missing')
            require(project_path(budget['native_current_scalar_path'])==project_path(state['global_tracer_path']),'OpenFOAM budget current native scalar differs')
            row.update(prescribed_top_c=top,prescribed_left_c=0.,inlet_evidence=dict(top_observed_min_key_path=['native_state_paths',step,'actual_top_c_min'],top_observed_max_key_path=['native_state_paths',step,'actual_top_c_max'],checks=files['checks'],left_zero_evidence='Accepted pinned full auditor reads every actual saved inlet_left fixedValue and checks its values against0 within existing1e-12; no separate persisted left-min/max audit field is claimed',auditor_source=source_reference),budget=json_reference(files,['reconstructed_scalar_steps',step-1,'fv_helper_return']),residual=json_reference(files,['actual_native_residuals',step-1,'final_residual']))
        else:
            require(state['global_clock'] is None and not state['rank_clock_paths'],'Initial actual clock availability differs')
            row.update(budget=None,residual=None,initial_state='Actual accepted saved c0; no invented solved budget/residual or clock metadata')
        out.append(row)
    return out

def load_solver(solver,args,evidence):
    files=root_acceptance(solver,args,evidence);checks=evidence.json(files['checks']['path'],files['checks']['sha256']);manifest=evidence.json(files['source_manifest']['path'],files['source_manifest']['sha256'])
    index=manifest_inputs(solver,manifest,evidence);source_path=PRIVATE/AUDITORS[solver][0];source_sha=AUDITORS[solver][1]
    require(source_path.resolve() in index and index[source_path.resolve()]==source_sha,'Accepted full audit source is not the reviewed source')
    evidence.pin(source_path,source_sha)
    if solver=='fluent':
        output_pins={project_path(row['path']):row['sha256'] for row in manifest['outputs']}
        require(all(output_pins.get(project_path(files[key]['path']))==files[key]['sha256'] for key in ('checks','fields')),'Fluent exporter output manifest disagrees with root acceptance')
        rows=fluent_rows(checks,files,index,evidence)
    else:rows=openfoam_rows(checks,files,index,evidence)
    return files,rows

def write_json(path,value):
    with path.open('x',encoding='utf-8',newline='\n') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for solver in ('fluent','openfoam'):
        for name in ('checks','manifest','fields','acceptance'):parser.add_argument('--'+solver+'-'+name,type=Path,required=True)
        parser.add_argument('--'+solver+'-acceptance-sha256',required=True)
    parser.add_argument('--targets',type=Path,required=True);parser.add_argument('--source-sha256',required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    evidence=Evidence();source_ref=evidence.ref(Path(__file__).resolve(),args.source_sha256)
    out=project_path(args.out);require(out.parent==PRIVATE.resolve() and re.fullmatch(r'visual-input-selection-v\d{3}',out.name) and not out.exists(),'Fresh private versioned visual selector directory required')
    require(shutil.disk_usage(ROOT).free>=10*1024**3,'E free-space floor not met')
    targets=evidence.ref(args.targets,TARGET_SHA);plotter=evidence.ref(PRIVATE/PLOTTER[0],PLOTTER[1])
    accepted={};adapters={}
    for solver in ('fluent','openfoam'):accepted[solver],adapters[solver]=load_solver(solver,args,evidence)
    # All actual acceptance, file, source and row integrity checks finish before any output write.
    before=evidence.finish();out.mkdir()
    solver_selection={}
    for solver in ('fluent','openfoam'):
        path=out/(solver+'-native-state-adapter.json');write_json(path,dict(schema_version=1,solver=solver,source='Actual root-accepted full audit; no new physical states',states=adapters[solver]))
        files=accepted[solver];manifest_ref=dict(path=str(path),sha256=sha(path))
        solver_selection[solver]=dict(full_actual_audit_accepted_by_root=True,accepted_evidence=[files['acceptance'],files['checks'],files['source_manifest'],dict(path=str(PRIVATE/AUDITORS[solver][0]),sha256=AUDITORS[solver][1])],fields=files['fields'],manifest=manifest_ref,residual_label='Native scaled UDS residual' if solver=='fluent' else 'Native final tracer matrix residual')
    selection=dict(schema_version=1,root_authorizes_offline_render=False,render_permission_pending=True,actual_paired_physical_steps=800,dt_s=DT,context_inputs=[source_ref,targets,plotter],solvers=solver_selection,physical_validation=False,mesh_independence=False,time_step_independence=False,complete_washout_established=False)
    selection_path=out/'selection-awaiting-root-render-review.json';write_json(selection_path,selection)
    after=evidence.finish();require(before==after,'Inputs changed while writing adapters')
    write_json(out/'builder-input-manifest.json',dict(all_declared_inputs_hash_stable=True,files=after))
    write_json(out/'builder-checks.json',dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),builder_source=source_ref,actual_root_full_audits_accepted=True,actual_state_rows_each=801,initial_solved_budget_or_residual_invented=False,scalar_fields_created_or_modified=False,native_paths_from_accepted_audits=True,actual_left_zero_evidence_scope='Fluent actual perstate saved-face error against0; OpenFOAM source-bound accepted auditor checks actual saved left fixedValue against0',render_authorized=False,render_started=False,CFD_or_WSL_or_GUI_started=False,public_files_modified=False))
    outputs=[dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path)) for path in sorted(out.iterdir()) if path.is_file()]
    write_json(out/'builder-output-manifest.json',dict(manifest_self_excluded=True,outputs=outputs))
    print(json.dumps(dict(out=str(out),selection=str(selection_path),selection_sha256=sha(selection_path),root_render_review_pending=True,render_started=False)),flush=True)
if __name__=='__main__':main()
