"""ROOT-authorized bounded Fluent tracer qualification: separate5-step c1/control and c0/pulse.

SOURCE PREPARATION ONLY until an exact root authorization is supplied. One hidden
two-process session, no flow initialization, no production entry. --help does no native work.
"""
from pathlib import Path
import argparse,ast,ctypes,datetime,hashlib,importlib.util,json,os,re,shutil,subprocess,sys,traceback
ROOT=Path('E:/CFD-Learning-Lab')
PRIVATE=ROOT/'local/prospective-tracer-v001'
QUALIFICATION=PRIVATE/'native-flow-qualification-v001.json'
QUALIFICATION_SHA='dccbb17e813f2552298332dc3fd7ff713e30df71a70508316a30831995ccd312'
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
REPO=ROOT/'repo/cfd-learning-lab/scripts'
ACTION='tracer_tjunction_fluent_qualification_v002'
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

def bounded_driver_source():
    tree=ast.parse(Path(__file__).read_text(encoding='utf-8'));solves=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.Call):
            name=node.func.attr if isinstance(node.func,ast.Attribute) else node.func.id if isinstance(node.func,ast.Name) else ''
            require(name not in FORBIDDEN_CALLS-{'dual_time_iterate'},'Forbidden initialization/flow/arbitrary-evaluation call: '+name)
            if name=='dual_time_iterate':
                require(not node.args and {k.arg:ast.literal_eval(k.value) for k in node.keywords}=={'time_step_count':1,'max_iter_per_step':60},'Only one-step/60-inner scalar calls are admitted');solves.append(node)
    require(len(solves)==1,'Exactly one bounded scalar solve site required')

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
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--authorization',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);a=parser.parse_args();bounded_driver_source()
    auth_path=a.authorization.resolve();out=a.out.resolve();auth=json.loads(auth_path.read_text(encoding='utf-8-sig'))
    expected={'approved':True,'approved_by':'root','action':ACTION,'phase':'qualification-5-control-plus-5-pulse','solver_processes':2,'control_steps':5,'pulse_steps':5,'maximum_physical_steps':10,'dt_s':.025,'max_inner_iterations':60,'setup_accepted_by_root':True}
    for key,value in expected.items():require(type(auth.get(key)) is type(value) and auth[key]==value,'Root authorization mismatch: '+key)
    require(inside(auth_path,ROOT/'local') and auth.get('source_sha256')==sha(__file__) and Path(auth['output_directory']).resolve()==out,'Authorization source/output/private path differs')
    require(out.parent==ROOT/'native/runs/2026-11-02-tracer-tjunction-v001' and re.fullmatch(r'fluent-qualification-v\d{3}',out.name) and not out.exists(),'Fresh versioned qualification output required')
    targets=PRIVATE/'QUALIFICATION-TARGETS-v001.md';targets_sha='173fa3a540d36e4dad195f25db5586d6c15527a6df592dca8076b1999ecb41cd';require(sha(targets)==auth.get('targets_sha256')==targets_sha,'Predeclared qualification targets changed')
    setup=ROOT/'native/runs/2026-11-02-tracer-tjunction-v001/fluent-setup-probe-v007';case=setup/'setup-initial.cas.h5';data=setup/'setup-initial.dat.h5';closure=PRIVATE/'fluent-setup-probe-v007-attempt-closure-v001.json'
    require(sha(case)=='511a778e3dfac16a9dcbf533595bcbf2d43022c096ae8d4790887197ec5b3b94' and sha(data)=='c74d1b7aca58b301820e35e18d76445d00d927ec8f6c4d4cd11c75ffa83f084d' and sha(closure)==auth.get('setup_closure_sha256')=='22fc7f439785880f58db05505b9b0c4b5fa90bfe297bbb585b08deb56f250c06','Accepted native c0 setup ancestry differs')
    evidence=auth.get('root_acceptance_evidence',[]);require(isinstance(evidence,list) and evidence,'Root-reviewed setup/WSL closure evidence required')
    for row in evidence:
        f=Path(row['path']).resolve();require(inside(f,ROOT/'local') and f.is_file() and sha(f)==row['sha256'],'Root acceptance/closure evidence changed')
    require(sha(QUALIFICATION)==QUALIFICATION_SHA and EXECUTABLE.is_file(),'Original qualification or installed Fluent missing')
    q=json.loads(QUALIFICATION.read_text(encoding='utf-8'));watched={Path(row['path']):row['sha256_before'] for row in q['source_manifest']}
    pins={PRIVATE/'scalar_fv_audit_v002.py':'75ceec60add2f5c53fb82585d105dcc5cd7c9d8d19e1ea8ac3f2dde14aa455cd',PRIVATE/'scalar_native_geometry_v003.py':'06f6b4089116d8e7cdab8b4050c30d2939d3876381ebbaa878c7908e6f2d1903',PRIVATE/'tracer_tjunction_fluent_setup_probe_v012.py':'6ac86a99d4b6fd858742734a91dbdc78d4982c61ff9b766979d8970d58501af8',ROOT/'environments/fluent/Lib/site-packages/ansys/fluent/core/generated/v261/solver/settings.py':'4a314bc418fdd2a9cfefd21b6c7a7286545e3ba7a032eeb1976e9fb6f5678d8e'}
    watched.update(pins);watched.update({case:sha(case),data:sha(data),closure:sha(closure),targets:targets_sha})
    require(all(sha(f)==digest for f,digest in watched.items()),'Pinned input/source changed before import')
    ram=memory_bytes();free=shutil.disk_usage(ROOT).free;require(ram>=4*1024**3 and free>=10*1024**3,'Require4GiB RAM/10GiB E: free');idle=require_idle()
    import numpy as np,h5py,time
    from ansys.fluent.core.offline import _lispy
    fv=load_module('tracer_qualification_fv_checked',PRIVATE/'scalar_fv_audit_v002.py');geo=load_module('tracer_qualification_geometry_checked',PRIVATE/'scalar_native_geometry_v003.py');geometries,geometry_proof=geo.qualified_geometry();reader=sys.modules['qualified_original_tjunction_fluent'];g=geometries['fluent']['stencil']
    require(g['cell_count']==4400 and g['face_count']==9040,'Actual qualified Fluent geometry differs')
    out.mkdir(parents=True);[(out/name).mkdir() for name in ('inputs','examples','control','pulse','reports')]
    temporary=ROOT/'temp/tracer-tjunction'/out.parent.name/out.name;require(not temporary.exists(),'Fresh project temp required');temporary.mkdir(parents=True)
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
    copied_case=out/'inputs'/case.name;copied_data=out/'inputs'/data.name
    for src,dst in ((case,copied_case),(data,copied_data),(auth_path,out/'root-authorization.json'),(Path(__file__),out/'executed-source.py'),(targets,out/'QUALIFICATION-TARGETS-v001.md')):shutil.copyfile(src,dst);require(sha(src)==sha(dst),'Immutable native/source copy differs')
    write_json(out/'source-provenance.json',{'watched':[{'path':str(f),'sha256':v} for f,v in watched.items()],'root_evidence':evidence,'native_setup_case_sha256':sha(case),'native_setup_data_sha256':sha(data)});write_json(out/'qualified-geometry.json',geometry_proof)
    events=[];discovery={};session=None;loaded=False;counts={'control':0,'pulse':0};phase=None;all_states={};optional_errors={};transcript=out/'native-transcript.txt'
    def stage(name,**values):
        row={'stage':name,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**values};events.append(row);write_json(out/'events.json',events);write_json(out/'status.json',row);print(json.dumps(row,default=json_default),flush=True)
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
    def settings_state():
            st=session.settings
            return {'general':st.setup.general.get_state(),'models':st.setup.models.get_state(),'material':st.setup.materials.fluid['air'].get_state(),'scalars':st.setup.user_defined.scalars.get_state(),'controls':st.solution.controls.get_state(),'residual':st.solution.monitor.residual.get_state(),'methods':st.solution.methods.get_state(),'run_calculation':st.solution.run_calculation.get_state(),'boundaries':{n:obj.get_state() for n,obj in [('inlet_left',st.setup.boundary_conditions.velocity_inlet['inlet_left']),('inlet_top',st.setup.boundary_conditions.velocity_inlet['inlet_top']),('outlet',st.setup.boundary_conditions.pressure_outlet['outlet']),('walls',st.setup.boundary_conditions.wall['walls'])]}}
    def set_boundary(name,obj,condition,value):
            group=obj.uds;types=group.uds_bc.get_object_names();values=group.uds.get_object_names();require(len(types)==len(values)==1 and set(types)==set(values),'Single runtime scalar BC required: '+name)
            key=types[0];pick(name+'.UDS_condition',group.uds_bc[key],condition);entry=group.uds[key];pick(name+'.UDS_value_method',entry.option,{'value'});entry.value.set_state(value)
            discovery[name+'.UDS_configured']={'native_key':key,'condition':group.uds_bc[key].get_state(),'value_state':entry.get_state(),'semantic_status':'Runtime setting readback only; outlet total-versus-diffusive interpretation remains unqualified.'};write_json(out/'runtime-discovery.json',discovery)
    reference=native_fields(copied_case,copied_data);reference_mesh=mesh_fingerprint(copied_case);require(np.array_equal(reference['phi'],geometries['fluent']['native_mass_phi']),'Frozen native phi differs from qualified geometry')
    _,_,_,zones=reader.native_faces(copied_case,copied_data);physical_ids=np.flatnonzero(g['patch']!='interior')+1;require(len(physical_ids)==480,'Actual physical face coverage differs')
    def native_scalar(path):
        with h5py.File(path,'r') as h:
            cells=h['results/1/phase-1/cells'];current=reader.chunk_array(cells['SV_UDS_0'],4400);old=reader.chunk_array(cells['SV_UDS_0_M1'],4400);faces=h['results/1/phase-1/faces/SV_UDS_0'];values=np.full(9040,np.nan);seen=np.zeros(9040,dtype=bool)
            for chunk in faces.values():
                lo=int(np.asarray(chunk.attrs['minId']).ravel()[0])-1;hi=int(np.asarray(chunk.attrs['maxId']).ravel()[0]);require(not seen[lo:hi].any(),'Overlapping native boundary c chunks');values[lo:hi]=np.asarray(chunk[()]).ravel();seen[lo:hi]=True
            require(np.array_equal(np.flatnonzero(seen)+1,physical_ids) and np.isfinite(values[seen]).all(),'Actual boundary c480 coverage differs')
        require(np.isfinite(current).all() and np.isfinite(old).all(),'Nonfinite actual native current/M1 scalar');return current,old,values
    def scalar_controls():
        st=session.settings
        return {'equations':st.solution.controls.equations.get_state(),'UDS_residual':st.solution.monitor.residual.equations['uds-0'].get_state(),'Gamma':st.setup.materials.fluid['air'].uds_diffusivity.get_state(),'inlet_diffusion':st.setup.user_defined.scalars.inlet_diffusion.get_state(),'transient_formulation':st.solution.methods.transient_formulation.get_state(),'dt_s':st.solution.run_calculation.parameters.time_step_size.get_state(),'max_inner':st.solution.run_calculation.parameters.max_iter_per_time_step.get_state(),'scalar_boundaries':{n:o.uds.get_state() for n,o in [('inlet_left',st.setup.boundary_conditions.velocity_inlet['inlet_left']),('inlet_top',st.setup.boundary_conditions.velocity_inlet['inlet_top']),('walls',st.setup.boundary_conditions.wall['walls']),('outlet',st.setup.boundary_conditions.pressure_outlet['outlet'])]}}
    def controls_gate(state,left):
        require(state['equations']=={'flow':False,'uds-0':True} and state['UDS_residual']['check_convergence'] is True and state['UDS_residual']['absolute_criteria']==1e-8,'Actual scalar-only equation/residual controls differ')
        require(state['inlet_diffusion'] is True and state['transient_formulation']=='unsteady-1st-order' and abs(state['dt_s']-.025)<=1e-12 and state['max_inner']==60,'Actual scalar time/diffusion controls differ')
        require(state['Gamma']['option']=='defined-per-uds' and state['Gamma']['uds_diffusivities']['uds-0']=={'option':'value','value':.001},'Actual Gamma differs')
        for name,value in [('inlet_left',left),('inlet_top',1.),('walls',0.),('outlet',0.)]:
            b=state['scalar_boundaries'][name];require(b['uds_bc']['User Scalar 0']==('Specified Value' if name.startswith('inlet') else 'Specified Flux') and b['uds']['User Scalar 0']['value']==value,'Actual scalar boundary differs: '+name)
    def static_state():
        state=settings_state();state.pop('run_calculation');return state
    def residual_history(path):
        with h5py.File(path,'r') as h:r=h['results/residuals/phase-1/uds-0'];indices=np.asarray(r['iterations'][()]).ravel();values=np.asarray(r['data'][()])
        require(indices.size and np.array_equal(indices,indices.astype(int)) and np.all(np.diff(indices)>0) and values.shape[0]==len(indices) and values.shape[1]>=2,'Native UDS residual index/shape invalid');return indices.astype(int),values
    def new_residual(path,previous,start_line):
        indices,values=residual_history(path);mask=indices>previous;new_indices=indices[mask];new_values=values[mask,:2];require(1<=len(new_indices)<=60 and int(new_indices[-1])-previous<=60 and np.isfinite(new_values).all(),'New own residual rows/inner cap invalid')
        iteration=int(new_indices[-1]);numerator=float(new_values[-1,0]);denominator=float(new_values[-1,1]);require(numerator>=0 and denominator>=0,'Negative new UDS residual numerator/denominator')
        text=transcript.read_text(encoding='utf-8',errors='replace');headers=None;printed=[]
        for lineno,line in enumerate(text.splitlines()):
            words=line.split()
            if words and words[0]=='iter' and 'uds-0' in words:headers=words;continue
            if lineno>=start_line and headers and words and words[0].isdigit() and int(words[0])==iteration:
                position=headers.index('uds-0')
                try:printed.append(float(words[position]))
                except (ValueError,IndexError):pass
        if denominator>0:value=numerator/denominator;method='native new numerator/positive native denominator'
        else:
            require(printed and float(format(numerator,'.5g'))==printed[-1],'Zero normalization lacks own printed native residual proof');value=printed[-1];method='actual newly printed solver value; native ratio undefined because denominator0'
        require(np.isfinite(value) and value>=0 and value<=1e-8,'Own new terminal UDS residual exceeds1e-8 or unproved')
        return {'last_iteration':iteration,'new_iterations':new_indices.tolist(),'new_native_residual_rows':new_values.tolist(),'terminal_numerator':numerator,'terminal_denominator':denominator,'own_terminal_scaled_residual':float(value),'method':method,'new_printed_terminal_candidates':printed,'all_inherited_indices_at_or_below':previous,'new_row_count':len(new_indices)}
    def native_report(folder,step):
        command=session.settings.results.report.flow.uds_flow;record={'step':step,'physical_time_s':step*.025,'active':bool(command.is_active()),'argument_names':list(command.argument_names),'domain_allowed_values':command.domain.allowed_values(),'native_flux_semantics':'UNCLASSIFIED; no reconstructed values are labeled native'};before=len(transcript.read_bytes().splitlines()) if transcript.exists() else 0
        try:record['return_value']=command();record['executed']=True
        except Exception:record.update(executed=False,error=traceback.format_exc())
        text=transcript.read_text(encoding='utf-8',errors='replace') if transcript.exists() else '';record['transcript_start_line']=before;record['transcript_end_line']=len(text.splitlines());(folder/f'state-{step:04d}-native-UDS-report.txt').write_text('\n'.join(text.splitlines()[before:])+'\n',encoding='utf-8');write_json(folder/f'state-{step:04d}-native-UDS-report.json',record);return record
    def snapshot(folder,step,left,previous=None):
        st=session.settings;cp=folder/f'state-{step:04d}.cas.h5';dp=cp.with_name(f'state-{step:04d}.dat.h5');st.file.write_case_data(file_name=str(cp));flow=native_fields(cp,dp);clock=native_clock(dp);mesh=mesh_fingerprint(cp);changes=fv.frozen_field_changes(reference,flow)
        require(all(v<=1e-12 for v in changes.values()) and mesh['sha256']==reference_mesh['sha256'],'Frozen native flow or complete mesh changed')
        require(abs(clock['flow-time']-step*.025)<=1e-12 and clock['time-step']==step and abs(clock['physical-time-step']-.025)<=1e-12,'Actual native physical clock differs')
        current,old,face_c=native_scalar(dp);small=scalar_controls();controls_gate(small,left);boundary={}
        for name,value in [('inlet_left',left),('inlet_top',1.),('walls',None),('outlet',None)]:
            pick=g['patch']==name;observed=face_c[pick];boundary[name]={'native_face_ids':g['native_face_ids'][pick].tolist(),'count':int(pick.sum()),'minimum':float(observed.min()),'maximum':float(observed.max()),'prescribed_value':value,'max_abs_difference_from_prescribed':float(np.max(abs(observed-value))) if value is not None else None,'initial_unsolved_faces_are_not_BC_enforcement_proof':step==0}
        record={'step':step,'clock':clock,'case':str(cp),'data':str(dp),'case_sha256':sha(cp),'data_sha256':sha(dp),'native_scalar_dataset':'results/1/phase-1/cells/SV_UDS_0/1','native_old_scalar_dataset':'results/1/phase-1/cells/SV_UDS_0_M1/1','current_min':float(current.min()),'current_max':float(current.max()),'M1_min':float(old.min()),'M1_max':float(old.max()),'M1_equal_to_previous_current_max_abs':float(np.max(abs(old-previous))) if previous is not None else None,'frozen_native_max_abs_changes':changes,'native_mesh_sha256':mesh['sha256'],'scalar_controls':small,'native_boundary_scalar_values':boundary}
        if previous is not None:require(np.max(abs(old-previous))<=1e-12,'Actual native M1 differs from preceding actual current scalar')
        np.savez_compressed(folder/f'state-{step:04d}-native-fields.npz',c=current,M1=old,boundary_c=face_c,**flow);write_json(folder/f'state-{step:04d}-native-readback.json',record);return current,record
    try:
        import ansys.fluent.core as pf
        from ansys.fluent.core.utils.fluent_version import FluentVersion
        require(pf.__version__=='0.43.0','PyFluent changed');pf.config.examples_path=str(out/'examples');ram_now=memory_bytes();free_now=shutil.disk_usage(ROOT).free;require(ram_now>=4*1024**3 and free_now>=10*1024**3,'Fresh launch resource gate failed');require_idle();stage('launching_one_hidden_qualification_session',available_ram_bytes=ram_now,free_E_bytes=free_now,maximum_physical_steps=10)
        session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True);session.transcript.start(file_name=str(transcript));version=session.get_fluent_version();write_json(out/'actual-version.json',{'type':type(version).__qualname__,'str':str(version),'value':version.value,'name':version.name});require(version==FluentVersion.v261,'Exact Fluent261 required')
        for phase,left,initial in [('control',1.,1.),('pulse',0.,0.)]:
            folder=out/phase;st=session.settings;st.file.read_case_data(file_name=str(copied_case));loaded=True;require(counts[phase]==0,'No native replay');write_json(folder/'loaded-settings.json',settings_state())
            st.solution.monitor.residual.equations['uds-0'].absolute_criteria=1e-8;st.solution.monitor.residual.equations['uds-0'].check_convergence=True;st.solution.run_calculation.parameters.max_iter_per_time_step=60
            set_boundary('inlet_left',st.setup.boundary_conditions.velocity_inlet['inlet_left'],{'specifiedvalue'},left);set_boundary('inlet_top',st.setup.boundary_conditions.velocity_inlet['inlet_top'],{'specifiedvalue'},1.)
            if phase=='control':st.solution.initialization.patch.calculate_patch(variable='uds-0',value=1.,cell_zones=['fluid'])
            c0,row=snapshot(folder,0,left);require(np.max(abs(c0-initial))<=1e-12,'Actual initial current scalar not prescribed');initial_M1=row['M1_min']==initial and row['M1_max']==initial;write_json(folder/'initial-old-time-storage-classification.json',{'initial_current_prescribed':True,'initial_M1_matches_initial_current':initial_M1,'initial_M1_min':row['M1_min'],'initial_M1_max':row['M1_max'],'M1_not_assumed':True,'step1_M1_must_match_actual_preceding_current':True})
            static_before=static_state();write_json(folder/'static-settings-before.json',static_before);row['native_report']=native_report(folder,0);rows=[row];previous=c0;previous_iteration=int(residual_history(Path(row['data']))[0][-1]);require(previous_iteration==2000,'Expected inherited unsolved residual index2000 differs')
            if phase=='control':
                try:write_json(out/'autosave-actual-metadata.json',{'canonical_path':'solution.calculation_activity.auto_save','state':st.solution.calculation_activity.auto_save.get_state(),'configuration_changed':False,'filename_alignment_qualified':False,'future_production_still_held':True})
                except Exception:write_json(out/'autosave-actual-metadata.json',{'canonical_path':'solution.calculation_activity.auto_save','metadata_error':traceback.format_exc(),'configuration_changed':False,'manual_five_step_saves_retained':True,'future_production_still_held':True})
            boundary={'inlet_left':{'value':left,'diffusion':'dirichlet','evidence':'Actual prescribed Specified Value and native inlet-diffusion on; reconstruction'},'inlet_top':{'value':1.,'diffusion':'dirichlet','evidence':'Actual prescribed Specified Value and native inlet-diffusion on; reconstruction'},'walls':{'diffusion':'zero_gradient','evidence':'Candidate Specified Flux0 interpreted as zero diffusion; control/pulse qualification hypothesis'},'outlet':{'diffusion':'zero_gradient','evidence':'Candidate advective outlet plus zero diffusion; constant-c1 preservation and reconstructed budget must qualify interpretation'}}
            require(np.all(reference['phi'][g['patch']=='outlet']>=0),'Actual native outlet backflow is unsupported')
            for step in range(1,6):
                require(counts[phase]==step-1 and sum(counts.values())<10,'Bounded source solve budget exceeded');line_start=len(transcript.read_bytes().splitlines());stage('advancing_one_scalar_step',case=phase,step=step,max_inner=60)
                st.solution.run_calculation.dual_time_iterate(time_step_count=1,max_iter_per_step=60);counts[phase]+=1
                current,row=snapshot(folder,step,left,previous);row['own_new_residual']=new_residual(Path(row['data']),previous_iteration,line_start);previous_iteration=row['own_new_residual']['last_iteration'];terms,defect,audit=fv.audit_step(g,previous,current,reference['phi'],rho=reference['rho'],gamma=.001,dt_s=.025,boundary=boundary,normalization_time_s=2.,global_floor_kg_m_s=.02)
                require(audit['global_Euler_defect_relative']<=1e-5,'Reconstructed global Euler budget fails predeclared target')
                if phase=='control':require(audit['constant_c1_max_absolute_deviation']<=1e-6,'Constant-c1 control not preserved')
                else:require(audit['scalar_min']>=-1e-8 and audit['scalar_max']<=1+1e-8,'Pulse scalar bounds fail')
                np.savez_compressed(folder/f'state-{step:04d}-RECONSTRUCTED-fluxes.npz',**terms,reconstructed_cell_defect_kg_s=defect);write_json(folder/f'state-{step:04d}-RECONSTRUCTED-budget.json',audit);row['reconstructed_budget']=audit;row['native_report']=native_report(folder,step);rows.append(row);write_json(folder/'states.json',rows);previous=current;stage('actual_scalar_state_saved',case=phase,step=step,time_s=row['clock']['flow-time'],global_Euler_reconstructed_defect=audit['global_Euler_defect_relative'])
            static_after=static_state();write_json(folder/'static-settings-after.json',static_after);require(static_before==static_after,'Static actual settings changed during scalar phase');all_states[phase]=rows;st.solution.monitor.residual.write(file_name=str(folder/'native-residuals.xy'))
            if phase=='pulse':require(np.max(abs(previous-c0))>1e-12,'Pulse did not genuinely change')
        try:
            graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1600;graphics.picture.y_resolution=1000;graphics.picture.invert_background=False;graphics.contour.create(name='tracer-c');contour=graphics.contour['tracer-c'];contour.field='uds-0-scalar';contour.surfaces_list=[];contour.options.filled=True;contour.options.node_values=False;contour.range_options.auto_range=False;contour.range_options.minimum=0.;contour.range_options.maximum=1.;contour.display();graphics.views.camera.projection='orthographic';graphics.views.camera.position=[.4,.3,10.];graphics.views.camera.target=[.4,.3,0.];graphics.views.camera.field=[1.8,1.125];graphics.picture.save_picture(file_name=str(out/'pulse/native-scalar-c-final.png'));write_json(out/'pulse/native-contour-settings.json',contour.get_state())
        except Exception:optional_errors['native_contour']=traceback.format_exc()
        require(counts=={'control':5,'pulse':5} and all(len(v)==6 for v in all_states.values()),'Incomplete bounded qualification coverage');write_json(out/'summary.json',{'qualification_driver_completed':True,'actual_physical_steps_by_case':counts,'end_time_s_each':.125,'native_states_each':6,'all_frozen_predeclared_numerical_gates_passed':True,'native_report_BC_flux_classification_requires_root_review':True,'reconstructed_budget_not_native_matrix_flux':True,'optional_export_errors':optional_errors,'production_authorized':False,'physical_validation':False});stage('bounded_qualification_saved_for_root_review',physical_steps_by_case=counts,production_authorized=False)
    except Exception:
        write_json(out/'failure.json',{'error':traceback.format_exc(),'completed_step_calls_by_case':counts,'active_case':phase,'native_failure_clock_is_authoritative_if_RPC_interrupted':True,'production_authorized':False})
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
        write_json(out/'immutable-source-final-check.json',{'all_unchanged':all(sha(f)==v for f,v in watched.items()),'watched_file_count':len(watched)})
    print(json.dumps({'out':str(out),'actual_steps':counts,'production_authorized':False}))

if __name__=='__main__':main()
