"""Future ROOT-admitted setup-only Fluent UDS probe v012. ZERO physical steps; not a solver driver.

Preparation of this private source grants no native admission. --help performs no
native work. A future authorization must bind this exact source and completed
Magnus evidence. Runtime metadata choices fail closed if absent or ambiguous.
"""
from pathlib import Path
import argparse,ast,ctypes,datetime,hashlib,importlib.util,json,os,re,shutil,subprocess,sys,traceback
ROOT=Path('E:/CFD-Learning-Lab')
PRIVATE=ROOT/'local/prospective-tracer-v001'
QUALIFICATION=PRIVATE/'native-flow-qualification-v001.json'
QUALIFICATION_SHA='dccbb17e813f2552298332dc3fd7ff713e30df71a70508316a30831995ccd312'
PATCH_ALIAS_PROOF=PRIVATE/'fluent-native-patch-UDS0-alias-proof-v001.json'
PATCH_ALIAS_PROOF_SHA='2f7492378f39e7b9c203303e2888a9ce3de5a0d75018ce7df1cd8232fc63bc48'
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
REPO=ROOT/'repo/cfd-learning-lab/scripts'
ACTION='tracer_tjunction_fluent_setup_probe_v012'
sys.dont_write_bytecode=True
FORBIDDEN_CALLS={'iterate','dual_time_iterate','calculate','initialize','hybrid_initialize','standard_initialize','fm_initialize','read_journal','execute_tui','exec','eval','connect_to_fluent'}

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()

def require(condition,message):
    if not condition:raise RuntimeError(message)

def json_default(value):
    if isinstance(value,Path):return str(value)
    if isinstance(value,bytes):return {'bytes_hex':value.hex()}
    if hasattr(value,'item'):return value.item()
    raise TypeError(type(value).__name__)

def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=json_default)+'\n',encoding='utf-8')

def inside(path,parent):return path==parent or parent in path.parents

def no_physical_step_source():
    tree=ast.parse(Path(__file__).read_text(encoding='utf-8'))
    calls=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.Call):
            name=node.func.attr if isinstance(node.func,ast.Attribute) else node.func.id if isinstance(node.func,ast.Name) else ''
            if name in FORBIDDEN_CALLS:calls.append(name)
    require(not calls,'Forbidden physical-step/init/arbitrary-evaluation call in probe source: '+str(calls))

def memory_bytes():
    class MemoryStatus(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong),('total',ctypes.c_ulonglong),('available',ctypes.c_ulonglong),('total_page',ctypes.c_ulonglong),('available_page',ctypes.c_ulonglong),('total_virtual',ctypes.c_ulonglong),('available_virtual',ctypes.c_ulonglong),('available_extended',ctypes.c_ulonglong)]
    state=MemoryStatus();state.length=ctypes.sizeof(state)
    require(bool(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(state))),'Native memory query failed')
    return int(state.available)

def require_idle():
    command="Get-Process -ErrorAction SilentlyContinue | Where-Object { $_.ProcessName -match '^(fluent|cortex|fl[0-9])' } | Select-Object ProcessName,Id | ConvertTo-Json -Compress"
    native=subprocess.run(['powershell.exe','-NoProfile','-NonInteractive','-Command',command],capture_output=True,text=True,check=True).stdout.strip()
    require(not native,'Existing Fluent process; save/close outside this probe first: '+native)
    result=subprocess.run(['wsl.exe','--list','--running','--quiet'],capture_output=True,check=True)
    raw=result.stdout;running=raw.decode('utf-16-le' if b'\x00' in raw else 'utf-8',errors='replace').replace('\x00','').strip()
    require(not running,'Sequential solver gate requires WSL inactive; this source never stops it: '+running)
    return {'existing_fluent_processes':[],'running_wsl_distributions':[]}

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--authorization',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    no_physical_step_source()
    out=args.out.resolve();auth_path=args.authorization.resolve()
    require(out.parent.parent==ROOT/'native/runs' and bool(re.fullmatch(r'\d{4}-\d{2}-\d{2}-tracer-tjunction-v\d{3}',out.parent.name)) and bool(re.fullmatch(r'fluent-setup-probe-v\d{3}',out.name)) and not out.exists(),'Require a fresh versioned tracer native setup-probe output on E:')
    require(inside(auth_path,ROOT/'local') and auth_path.is_file(),'Future ROOT authorization must be an existing private E: file')
    auth=json.loads(auth_path.read_text(encoding='utf-8-sig'))
    for key,value in {'approved':True,'approved_by':'root','action':ACTION,'phase':'setup-only-zero-physical-steps','physical_steps_authorized':0,'solver_processes':2,'magnus_pair_and_local_record_complete':True}.items():
        require(type(auth.get(key)) is type(value) and auth[key]==value,'Authorization mismatch: '+key)
    require(Path(auth.get('output_directory','')).resolve()==out and auth.get('source_sha256')==sha(__file__),'Output/source is not bound to exact future ROOT authorization')
    completion=auth.get('magnus_completion_evidence',[])
    require(isinstance(completion,list) and len(completion)>0,'Require ROOT-bound actual Magnus completion evidence')
    for row in completion:
        p=Path(row['path']).resolve();require(inside(p,ROOT/'local') and p.is_file() and sha(p)==row['sha256'],'Magnus completion evidence path/hash differs')
    require(sha(QUALIFICATION)==QUALIFICATION_SHA,'Accepted read-only qualification changed')
    q=json.loads(QUALIFICATION.read_text());source_rows=q['source_manifest'];watched={}
    for row in source_rows:
        p=Path(row['path']);digest=sha(p);require(digest==row['sha256_before']==row['sha256_after'],'Qualified ancestor/source changed: '+str(p));watched[p]=digest
    selection=json.loads(Path(q['selection_path']).read_text());case=Path(selection['fluent_final_case']);data=Path(selection['fluent_final_data'])
    require(EXECUTABLE.is_file() and case.is_file() and data.is_file(),'Installed Fluent or actual inherited case/data missing')
    require(case==ROOT/'native/runs/2026-11-01-tjunction-v001/fluent/tjunction-final.cas.h5' and data==case.with_name('tjunction-final.dat.h5'),'Actual accepted baseline path differs')
    reader_path=Path(selection['fluent_checks_json']).parent/'executed-analyser.py';helper_path=REPO/'analyse_nozzle_fluent_v005.py'
    original_checks=json.loads(Path(selection['fluent_checks_json']).read_text())
    require(sha(reader_path)==original_checks['analyser_sha256'] and sha(helper_path)==original_checks['mapping_helper_sha256'],'Actual loader/helper differs from accepted audit source')
    interfaces=[]
    for path,name,arguments,defaults in ((reader_path,'native_faces',['case','data'],[]),(helper_path,'chunk_array',['group','count','tail'],['()'])):
        tree=ast.parse(path.read_text(encoding='utf-8'));functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name]
        require(len(functions)==1,'Actual reader interface absent/ambiguous: '+name);function=functions[0]
        actual_arguments=[a.arg for a in function.args.args];actual_defaults=[ast.unparse(x) for x in function.args.defaults]
        require(actual_arguments==arguments and actual_defaults==defaults,'Actual native reader signature differs: '+name)
        interfaces.append({'source_path':str(path),'source_sha256':sha(path),'function':name,'arguments':actual_arguments,'defaults':actual_defaults})
    require(any(isinstance(n,ast.ImportFrom) and n.module=='analyse_nozzle_fluent_v005' and any(a.name=='chunk_array' for a in n.names) for n in ast.parse(reader_path.read_text()).body),'Actual native_faces module import of chunk_array differs')

    ram=memory_bytes();free=shutil.disk_usage(ROOT).free;require(ram>=4*1024**3 and free>=10*1024**3,'Require4GiB available RAM and10GiB E: free')
    idle=require_idle()
    import numpy as np,h5py
    from ansys.fluent.core.offline import _lispy
    sys.path.insert(0,str(REPO))
    reader=load_module('tracer_original_native_reader',reader_path)
    out.mkdir(parents=True)
    for name in ('inputs','examples','reports'):(out/name).mkdir()
    write_json(out/'actual-loader-interface-proof.json',{'interfaces':interfaces,'baseline_case_path':str(case),'baseline_data_path':str(data),'baseline_case_sha256':sha(case),'baseline_data_sha256':sha(data),'reader_and_helper_match_accepted_audit_sources':True})
    temporary=ROOT/'temp/tracer-tjunction'/out.parent.name/out.name
    require(not temporary.exists(),'Require a fresh versioned project temp directory')
    temporary.mkdir(parents=True)
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
    copied_case=out/'inputs'/case.name;copied_data=out/'inputs'/data.name
    for src,dst in ((case,copied_case),(data,copied_data)):
        shutil.copyfile(src,dst);require(sha(src)==sha(dst),'Inherited native byte copy differs')
    shutil.copyfile(auth_path,out/'root-authorization.json');shutil.copyfile(__file__,out/'executed-source.py')
    events=[];discovery={};session=None;loaded=False;report=None
    def stage(name,**values):
        row={'stage':name,'time':datetime.datetime.now().astimezone().isoformat(),**values};events.append(row);write_json(out/'history.json',events);write_json(out/'status.json',row);print(json.dumps(row,default=json_default),flush=True)
    def metadata(label,obj):
        values=obj.allowed_values();row={'allowed_values':values}
        try:row['current']=obj.get_state()
        except Exception as error:row['current_state_unavailable']=str(error)
        discovery[label]=row;write_json(out/'runtime-discovery.json',discovery)
        require(isinstance(values,(list,tuple)) and len(values)>0,'Runtime enumeration unavailable: '+label)
        return list(values)
    def pick(label,obj,canonical):
        values=metadata(label,obj);normal=lambda value:re.sub(r'[^a-z0-9]','',str(value).lower())
        matches=[v for v in values if normal(v) in set(canonical)]
        require(len(matches)==1,'Runtime semantic option absent/ambiguous: '+label+' '+str(values));obj.set_state(matches[0]);return matches[0]
    def native_clock(path):
        with h5py.File(path,'r') as handle:values={x[0]:x[1] for x in _lispy.parse(handle['settings/Data Variables'][0].decode())[1]}
        return {k:float(values[k]) for k in ('flow-time','time-step','physical-time-step')}
    def native_fields(case_path,data_path):
        with h5py.File(data_path,'r') as handle:
            cells=handle['results/1/phase-1/cells'];fields={name:reader.chunk_array(cells[key],4400) for name,key in [('U','SV_U'),('V','SV_V'),('p','SV_P'),('rho','SV_DENSITY'),('mu','SV_MU_LAM')]}
        c0,c1,flux,zones=reader.native_faces(case_path,data_path);fields.update(phi=flux,c0=c0,c1=c1)
        require(all(np.isfinite(a).all() for a in fields.values()),'Nonfinite saved native flow array')
        return fields
    def inventory(path):
        rows=[]
        with h5py.File(path,'r') as handle:
            def visit(name,obj):
                if isinstance(obj,h5py.Dataset) and 'UDS' in name.upper():
                    row={'path':name,'shape':list(obj.shape),'dtype':str(obj.dtype),'attrs':{str(k):np.asarray(v).tolist() for k,v in obj.attrs.items()}}
                    if obj.dtype.kind in 'fiu' and obj.size:
                        a=np.asarray(obj[:]);row.update(finite=bool(np.isfinite(a).all()),minimum=float(a.min()),maximum=float(a.max()))
                    rows.append(row)
            handle.visititems(visit)
        return rows
    def mesh_fingerprint(case_path):
        rows=[]
        with h5py.File(case_path,'r') as handle:
            require('meshes/1' in handle,'Saved native meshes/1 group missing')
            def visit(name,obj):
                if not isinstance(obj,h5py.Dataset):return
                values=np.asarray(obj[()])
                # Variable-length string/object arrays must not hash Python object
                # pointers. Canonical JSON preserves their actual values/bytes.
                payload=(json.dumps(values.tolist(),sort_keys=True,separators=(',',':'),allow_nan=False,default=json_default).encode('utf-8') if values.dtype.kind=='O' else values.tobytes(order='C'))
                row={'path':'meshes/1/'+name,'shape':list(obj.shape),'dtype':str(obj.dtype),'value_sha256':hashlib.sha256(payload).hexdigest(),'attributes':{str(k):np.asarray(v).tolist() for k,v in obj.attrs.items()}}
                rows.append(row)
            handle['meshes/1'].visititems(visit)
        rows.sort(key=lambda row:row['path']);require(len(rows)>0,'No native mesh datasets covered')
        digest=hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':'),allow_nan=False,default=json_default).encode('utf-8')).hexdigest()
        return {'definition':'Every dataset under meshes/1, with path/dtype/shape/value hash and dataset attributes; object values canonicalized instead of pointer bytes. Group-only attributes are outside this dataset fingerprint.','dataset_count':len(rows),'sha256':digest,'datasets':rows}
    baseline=native_fields(copied_case,copied_data);original_clock=native_clock(copied_data);baseline_mesh=mesh_fingerprint(copied_case)
    write_json(out/'copied-ancestor-native-mesh-fingerprint.json',baseline_mesh)
    def compare(label,c,d):
        current=native_fields(c,d);clock=native_clock(d);mesh=mesh_fingerprint(c)
        require(mesh['sha256']==baseline_mesh['sha256'] and mesh['dataset_count']==baseline_mesh['dataset_count'],'Saved native complete mesh dataset fingerprint changed at '+label)
        changes={k:{'bitwise_equal':bool(np.array_equal(current[k],a)),'maximum_absolute_change':float(np.max(abs(current[k]-a)))} for k,a in baseline.items()}
        require(all(current[k].shape==a.shape for k,a in baseline.items()),'Native field topology/order changed')
        require(all(r['maximum_absolute_change']<=1e-12 for r in changes.values()),'Frozen flow changed at '+label)
        require(clock['flow-time']==original_clock['flow-time'] and clock['time-step']==original_clock['time-step'],'Physical time/index advanced in setup-only probe')
        record={'label':label,'clock':clock,'flow_changes':changes,'native_mesh_fingerprint':mesh,'mesh_equals_copied_ancestor':True,'UDS_native_datasets':inventory(d)};write_json(out/(label+'-native-readback.json'),record);return record
    def snapshot(label):
        p=out/(label+'.cas.h5');session.settings.file.write_case_data(file_name=str(p));return compare(label,p,p.with_name(label+'.dat.h5'))
    def settings_state():
        st=session.settings
        return {'general':st.setup.general.get_state(),'models':st.setup.models.get_state(),'material':st.setup.materials.fluid['air'].get_state(),'scalars':st.setup.user_defined.scalars.get_state(),'controls':st.solution.controls.get_state(),'residual':st.solution.monitor.residual.get_state(),'methods':st.solution.methods.get_state(),'run_calculation':st.solution.run_calculation.get_state(),'boundaries':{n:obj.get_state() for n,obj in [('inlet_left',st.setup.boundary_conditions.velocity_inlet['inlet_left']),('inlet_top',st.setup.boundary_conditions.velocity_inlet['inlet_top']),('outlet',st.setup.boundary_conditions.pressure_outlet['outlet']),('walls',st.setup.boundary_conditions.wall['walls'])]}}
    def set_boundary(name,obj,condition,value):
        group=obj.uds;types=group.uds_bc.get_object_names();values=group.uds.get_object_names();require(len(types)==len(values)==1 and set(types)==set(values),'Single runtime scalar BC required: '+name)
        key=types[0];pick(name+'.UDS_condition',group.uds_bc[key],condition);entry=group.uds[key];pick(name+'.UDS_value_method',entry.option,{'value'});entry.value.set_state(value)
        discovery[name+'.UDS_configured']={'native_key':key,'condition':group.uds_bc[key].get_state(),'value_state':entry.get_state(),'semantic_status':'Runtime setting readback only; outlet total-versus-diffusive interpretation remains unqualified.'};write_json(out/'runtime-discovery.json',discovery)
    try:
        stage('admitted_setup_only_no_physical_steps',available_ram_bytes=ram,free_E_bytes=free,idle=idle,root_authorization_sha256=sha(auth_path))
        import ansys.fluent.core as pf
        require(pf.__version__=='0.43.0','PyFluent version changed; source review required')
        pf.config.examples_path=str(out/'examples')
        ram_now=memory_bytes();free_now=shutil.disk_usage(ROOT).free;require(ram_now>=4*1024**3 and free_now>=10*1024**3,'Fresh launch requires4GiB RAM/10GiB E: free')
        stage('fresh_resource_gate_before_launch',available_ram_bytes=ram_now,free_E_bytes=free_now)
        session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
        session.transcript.start(file_name=str(out/'native-transcript.txt'))
        from ansys.fluent.core.utils.fluent_version import FluentVersion
        version_object=session.get_fluent_version();version=str(version_object)
        write_json(out/'actual-connected-version-identity.json',{'type_module':type(version_object).__module__,'type_name':type(version_object).__qualname__,'str':version,'repr':repr(version_object),'value':getattr(version_object,'value',None),'name':getattr(version_object,'name',None),'expected_enum':'FluentVersion.v261','equals_v261':version_object==FluentVersion.v261})
        require(isinstance(version_object,FluentVersion) and version_object==FluentVersion.v261,'Connected Fluent version differs from exact FluentVersion.v261')
        stage('single_hidden_session_connected',version=version,processor_count=2)
        st=session.settings;st.file.read_case_data(file_name=str(copied_case));loaded=True
        before=settings_state();write_json(out/'loaded-settings.json',before);snapshot('loaded-baseline')
        require(before['models']['energy']['enabled'] is False and before['models']['viscous']['model']=='laminar','Inherited laminar/no-energy model changed')
        equations_before=st.solution.controls.equations.get_state();schemes_before=st.solution.methods.spatial_discretization.discretization_scheme.get_state()
        scalars=st.setup.user_defined.scalars;require(len(scalars.scalar)==0,'Inherited case unexpectedly contains scalars')
        scalars.scalar.resize(size=1);scalar=scalars.scalar[0]
        pick('solution_zones',scalar.solution_zones,{'allfluidzones'});pick('advective_flux_function',scalar.flux_fcn,{'massflowrate'})
        scalars.inlet_diffusion.set_state(True);scalars.update();require(len(scalars.scalar)==1 and scalars.inlet_diffusion.get_state() is True,'Native scalar allocation/inlet diffusion readback failed')
        # Capture the fresh native UDS residual key/default criteria/scaling metadata.
        write_json(out/'post-UDS-allocation-settings.json',settings_state())
        diffusivity=st.setup.materials.fluid['air'].uds_diffusivity;pick('UDS_diffusivity_method',diffusivity.option,{'definedperuds'})
        names=diffusivity.uds_diffusivities.get_object_names();require(len(names)==1,'Exactly one actual per-UDS diffusivity property required')
        coefficient=diffusivity.uds_diffusivities[names[0]];pick('UDS0_constant_diffusivity',coefficient.option,{'constant','value'});coefficient.value.set_state(.001);require(abs(float(coefficient.value.get_state())-.001)<=1e-15,'Native Gamma constant readback differs')
        discovery['UDS_diffusivity']={'native_key':names[0],'state':coefficient.get_state(),'units_quantity_metadata':coefficient.value.get_attr('units-quantity'),'intended_Gamma_kg_m_s':.001,'rho_kg_m3':1.,'D_m2_s':.001};write_json(out/'runtime-discovery.json',discovery)
        eq=st.solution.controls.equations.get_state();added=sorted(set(eq)-set(equations_before));require(len(added)==1,'Actual new UDS equation key absent/ambiguous')
        for name in eq:st.solution.controls.equations[name]=(name==added[0])
        active=st.solution.controls.equations.get_state();require(active.get(added[0]) is True and all(v is False for k,v in active.items() if k!=added[0]),'Flow equation freeze failed')
        discovery['equations']={'before':equations_before,'actual_scalar_key':added[0],'active_after':active};write_json(out/'runtime-discovery.json',discovery)
        pick('solver_time',st.setup.general.solver.time,{'transient'});pick('Euler_formulation',st.solution.methods.transient_formulation,{'unsteady1storder'})
        pick('unsteady_function',scalar.unsteady_fcn,{'default'})
        schemes=st.solution.methods.spatial_discretization.discretization_scheme;new_schemes=sorted(set(schemes.get_state())-set(schemes_before));require(len(new_schemes)==1,'Actual new UDS discretization key absent/ambiguous');pick('UDS_advection',schemes[new_schemes[0]],{'firstorderupwind'})
        run=st.solution.run_calculation;run.parameters.time_step_size=.025;run.parameters.max_iter_per_time_step=60
        set_boundary('inlet_left',st.setup.boundary_conditions.velocity_inlet['inlet_left'],{'specifiedvalue','value'},0.)
        set_boundary('inlet_top',st.setup.boundary_conditions.velocity_inlet['inlet_top'],{'specifiedvalue','value'},1.)
        set_boundary('walls',st.setup.boundary_conditions.wall['walls'],{'specifiedflux','flux'},0.)
        # This setup-only candidate is not a verified outlet PDE boundary condition.
        # Native frozen flow has no reversed outlet faces. Its eventual zero-diffusion
        # versus total-flux semantics must be classified before production; a separately
        # ROOT-admitted physical qualification may investigate them.
        set_boundary('outlet',st.setup.boundary_conditions.pressure_outlet['outlet'],{'specifiedflux','flux'},0.)
        write_json(out/'configured-before-patch-settings.json',settings_state());snapshot('activated-before-patch')
        patch=st.solution.initialization.patch.calculate_patch
        def patch_argument(name):
            obj=getattr(patch,name);row={}
            try:row['allowed_values']=obj.allowed_values()
            except Exception as error:row['allowed_values_unavailable']=str(error)
            try:row['current']=obj.get_state()
            except Exception as error:row['current_unavailable']=str(error)
            return row
        patch_before={name:patch_argument(name) for name in patch.argument_names}
        write_json(out/'patch-all-argument-metadata-before-context.json',patch_before)
        domain_row=patch_before.get('domain',{});domain_values=domain_row.get('allowed_values',[])
        domain_values=[v for v in domain_values if isinstance(v,str) and v.strip()] if isinstance(domain_values,(list,tuple)) else []
        domain_current=domain_row.get('current');domain_choice=domain_values[0] if len(domain_values)==1 else domain_current if isinstance(domain_current,str) and domain_current.strip() and (not domain_values or domain_current in domain_values) else None
        domain_context={'actual_candidate':domain_choice,'passed_as_explicit_patch_argument':domain_choice is not None};patch_kwargs={'domain':domain_choice} if domain_choice is not None else {}
        if domain_choice is not None:
            try:
                patch.domain.set_state(domain_choice);domain_context['actual_readback']=patch.domain.get_state()
                if domain_context['actual_readback']==domain_choice:patch_kwargs['domain']=domain_choice;domain_context['passed_as_explicit_patch_argument']=True
            except Exception as error:domain_context['context_set_unavailable']=str(error)
        # Pass an actual discovered domain as a command argument even when
        # storing argument defaults is unsupported; omit only undiscovered context.
        patch_after={name:patch_argument(name) for name in ('domain','variable','cell_zones')}
        write_json(out/'patch-argument-metadata-after-context.json',{'arguments':patch_after,'domain_context':domain_context})
        catalog=session.fields._field_info._get_scalar_fields_info()
        write_json(out/'actual-native-scalar-field-catalog.json',catalog)
        scalar0=lambda name:isinstance(name,str) and bool(re.fullmatch(r'(?:(uds|udsi|userscalar|userdefinedscalar)0|uds0scalar|scalar0)',re.sub(r'[^a-z0-9]','',name.lower())))
        catalog_candidates=[{'identifier':name,'native_metadata':info} for name,info in catalog.items() if isinstance(info,dict) and (scalar0(name) or scalar0(info.get('display_name')))] if isinstance(catalog,dict) else []
        allowed=patch_after['variable'].get('allowed_values',[]);allowed=allowed if isinstance(allowed,(list,tuple)) else []
        matches=[name for name in allowed if scalar0(name)]
        zones=patch_after['cell_zones'].get('allowed_values',[]);zones=zones if isinstance(zones,(list,tuple)) else []
        selection={'variable_allowed_values':allowed,'cell_zones_allowed_values':zones,'unique_native_catalog_candidates':catalog_candidates,'domain_context':domain_context,'empty_variable_enumeration_classification':'Unavailable argument metadata when native catalog proves UDS0; empty enumeration alone does not establish absent scalar.','native_postpatch_c0_and_frozen_flow_proof_required':True,'patch1_identity_control_still_required_before_physical_steps':True}
        write_json(out/'scalar-only-patch-candidate-discovery.json',selection)
        require(len(catalog_candidates)==1,'Actual native field catalog UDS0 identifier/display label absent or ambiguous')
        if len(matches)==1:variable=matches[0];selection['candidate_source']='actual patch-variable enumeration'
        elif not allowed:
            require(sha(PATCH_ALIAS_PROOF)==PATCH_ALIAS_PROOF_SHA,'Pinned native UDS0 patch-alias proof changed')
            alias_proof=json.loads(PATCH_ALIAS_PROOF.read_text(encoding='utf-8'))
            require(alias_proof['accepted_source_case_sha256']==sha(case)==sha(copied_case) and alias_proof['accepted_source_data_sha256']==sha(data)==sha(copied_data),'Observed patch alias is not bound to this accepted native flow case/data')
            require(alias_proof['observed_connected_version_value']==version_object.value=='26.1.0' and alias_proof['observed_connected_version_name']==version_object.name=='v261' and alias_proof['observed_PyFluent_version']==pf.__version__=='0.43.0','Observed patch alias version binding differs')
            for row in alias_proof['bound_native_evidence']:
                observed=Path(row['path']);require(inside(observed,ROOT) and observed.is_file() and sha(observed)==row['sha256'],'Pinned native patch-alias evidence changed')
            require('uds-0-scalar is_not_in (pressure x-velocity y-velocity uds-0)' in alias_proof['raw_native_error'] and alias_proof['solver_returned_allowed_patch_values']==['pressure','x-velocity','y-velocity','uds-0'],'Native allowed patch-token observation differs')
            variable=alias_proof['actual_native_error_UDS0_patch_token'];require(variable=='uds-0','Observed scalar0 patch token differs')
            selection.update(candidate_source='exact uds-0 token observed in pinned native solver error allowed-values list; actual patch execution/postpatch proof still required',native_alias_proof_path=str(PATCH_ALIAS_PROOF),native_alias_proof_sha256=PATCH_ALIAS_PROOF_SHA)
        else:raise RuntimeError('Actual scalar-only patch enumeration excludes or ambiguously identifies UDS0')
        selection['actual_patch_variable_candidate']=variable;write_json(out/'scalar-only-patch-candidate-discovery.json',selection)
        require('fluid' in zones,'Inherited fluid patch zone absent from actual patch cell-zone metadata')
        patch(variable=variable,value=0.,cell_zones=['fluid'],**patch_kwargs)
        write_json(out/'setup-initial-settings.json',settings_state());final=snapshot('setup-initial')
        with h5py.File(out/'setup-initial.dat.h5','r') as handle:
            cells=handle['results/1/phase-1/cells'];names=[n for n in cells if re.fullmatch(r'SV_UDS_0',n)];require(len(names)==1,'Actual native UDS0 cell array absent/ambiguous');values=reader.chunk_array(cells[names[0]],4400);require(np.isfinite(values).all() and np.max(abs(values))<=1e-12,'Native scalar-only patch0 readback failed')
        # Persist fresh actual settings readback after the final native write, then
        # verify the stored JSON rather than intended setter values.
        post_path=out/'actual-post-settings.json';write_json(post_path,settings_state());post=json.loads(post_path.read_text(encoding='utf-8'))
        post_equations=post['controls']['equations'];gamma_key=discovery['UDS_diffusivity']['native_key']
        gamma=float(post['material']['uds_diffusivity']['uds_diffusivities'][gamma_key]['value'])
        post_checks={'only_actual_UDS_equation_active':post_equations.get(added[0]) is True and all(v is False for k,v in post_equations.items() if k!=added[0]),'energy_off':post['models']['energy']['enabled'] is False,'laminar':post['models']['viscous']['model']=='laminar','Gamma0p001':abs(gamma-.001)<=1e-15,'native_physical_dt0p025':abs(final['clock']['physical-time-step']-.025)<=1e-12,'physical_time_and_index_unchanged':final['clock']['flow-time']==original_clock['flow-time'] and final['clock']['time-step']==original_clock['time-step']}
        write_json(out/'actual-post-settings-verification.json',{'checks':post_checks,'all_checks_passed':all(post_checks.values()),'persisted_actual_settings_path':str(post_path),'persisted_actual_settings_sha256':sha(post_path),'native_final_case_sha256':sha(out/'setup-initial.cas.h5'),'native_final_data_sha256':sha(out/'setup-initial.dat.h5'),'native_final_clock':final['clock'],'actual_scalar_equation_key':added[0],'actual_Gamma_kg_m_s':gamma,'evidence_limit':'Settings JSON is actual live readback immediately after the final native save; native HDF separately proves dt/time/index and all flow/mesh arrays. It is not an independent reopen-and-reload of settings.'})
        require(all(post_checks.values()),'Actual persisted postsettings/native dt check failed')
        command=st.results.report.flow.uds_flow;domains=command.domain.allowed_values()
        write_json(out/'UDS-report-command-metadata.json',{'active':bool(command.is_active()),'argument_names':list(command.argument_names),'domain_allowed_values':domains,'report_invoked':False,'total_vs_diffusive_sign_units':'UNKNOWN; requires separately ROOT-admitted native qualification; no initialized report is conservation evidence.'})
        write_json(out/'source-ancestry-final-check.json',{'unchanged':all(sha(p)==v for p,v in watched.items()),'source_files':[{'path':str(p),'sha256_before':v,'sha256_after':sha(p)} for p,v in watched.items()]})
        require(all(sha(p)==v for p,v in watched.items()),'Original ancestor/source changed during probe')
        report={'setup_only_probe_completed':True,'actual_physical_steps':0,'actual_physical_time_advanced_s':0.,'ready_for_five_step_pilot':False,'production_authorized':False,'scalar_material_and_patch_settings_read_back':True,'frozen_flow_compared_all_native_cells_and_faces':True,'complete_mesh_dataset_fingerprint_compared_at_all_saved_setup_states':True,'native_physical_time_step_s':final['clock']['physical-time-step'],'actual_post_settings_checks':post_checks,'native_scalar_dataset':names[0],'native_UDS_inventory':final['UDS_native_datasets'],'unresolved':['Outlet flux0 versus diffusive-zero/backflow semantics','UDS report total-versus-diffusive/sign/units and actual allocated face-flux coverage','Scalar residual behavior during a separately admitted physical solve','Matched future pulse step80/81 schedule; this setup probe has no physical time loop'],'physical_validation':False}
        write_json(out/'summary.json',report);stage('setup_only_probe_saved_for_root_review',actual_physical_steps=0,ready_for_five_step_pilot=False,production_authorized=False)
    except Exception:
        failure={'error':traceback.format_exc(),'actual_physical_steps':0,'production_authorized':False,'runtime_discovery_required':True};write_json(out/'failure.json',failure);stage('setup_probe_stopped_no_physical_steps',**failure)
        if session is not None and loaded:
            try:session.settings.file.write_case_data(file_name=str(out/'failure-state.cas.h5'))
            except Exception:write_json(out/'failure-save-error.json',{'error':traceback.format_exc()})
        raise
    finally:
        if session is not None:
            try:session.transcript.stop()
            except Exception:write_json(out/'transcript-stop-error.json',{'error':traceback.format_exc()})
            try:session.exit()
            except Exception:write_json(out/'session-exit-error.json',{'error':traceback.format_exc()})
    print(json.dumps({'out':str(out),'actual_physical_steps':0,'setup_only_probe_completed':bool(report)},indent=2))

if __name__=='__main__':main()