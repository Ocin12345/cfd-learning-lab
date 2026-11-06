"""Future ROOT-admitted Fluent tracer continuation: existing pulse5 plus exactly795 new steps.

SOURCE PREPARATION ONLY. One hidden2-process session; stages6,80,81,800.
No initialization, scalar patch, replay, control run or production admission by preparation.
"""
from pathlib import Path
import argparse,ast,ctypes,datetime,hashlib,importlib.util,json,os,re,shutil,subprocess,sys,traceback
ROOT=Path('E:/CFD-Learning-Lab')
PRIVATE=ROOT/'local/prospective-tracer-v001'
QUALIFICATION=PRIVATE/'native-flow-qualification-v001.json'
QUALIFICATION_SHA='dccbb17e813f2552298332dc3fd7ff713e30df71a70508316a30831995ccd312'
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
REPO=ROOT/'repo/cfd-learning-lab/scripts'
ACTION='tracer_tjunction_fluent_production_v002'
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
    tree=ast.parse(Path(__file__).read_text(encoding='utf-8'));steps=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.Call):
            name=node.func.attr if isinstance(node.func,ast.Attribute) else node.func.id if isinstance(node.func,ast.Name) else ''
            require(name not in (FORBIDDEN_CALLS-{'dual_time_iterate'})|{'calculate_patch'},'Forbidden initialize/patch/flow/arbitrary-evaluation call: '+name)
            if name=='dual_time_iterate':
                kw={k.arg:ast.literal_eval(k.value) for k in node.keywords};require(not node.args and set(kw)=={'time_step_count','max_iter_per_step'} and kw['max_iter_per_step']==60,'Fixed scalar-step arguments required');steps.append(kw['time_step_count'])
    require(steps==[1,74,1,719] and sum(steps)==795,'Exactly staged795 new steps required')

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



PILOT=ROOT/'native/runs/2026-11-02-tracer-tjunction-v001/fluent-qualification-v001'
PILOT_FILES={'pulse/state-0000.cas.h5': 'cf61ab2c02d9097482eff5de5c04acfec0381b67755a1506a91d8a72798cac43', 'pulse/state-0000.dat.h5': 'ba6cdca5623f341a121b4a451ebe8ad8559bb244020d5bb30453373f2a063517', 'pulse/state-0001.cas.h5': 'a0ce435c533e0b79b74d327213259770f443a27e55dfd7ff6ef41c16c189e722', 'pulse/state-0001.dat.h5': '92afd4f2a3cbfa2a1aad51f950bf8d7f9b246428daa0fe9aa7038d3e15644e30', 'pulse/state-0002.cas.h5': '25f6f0df522954294f9c20ad0e49ad4c84ab8368407625534fc7669bb9a91b7c', 'pulse/state-0002.dat.h5': '7c3c8fc244b2995cfa919b28bd4d902ae90975f5a79b864f53cca996de66b4a5', 'pulse/state-0003.cas.h5': '3b2ffc6627578efdefd26cf7422e120f08afb5c809964099e6190cea3df5f0f2', 'pulse/state-0003.dat.h5': 'bbab5a534d018ce77ae8cfa5e45df97b4ed3391e278f3fb41d62566899831853', 'pulse/state-0004.cas.h5': 'a1e4d068ffd9a34abb541b2af57286d773ababbdb14c4ec782e5a32e45f77e4b', 'pulse/state-0004.dat.h5': '7c0bb0ae47433c47d6bf33a9b42f4fbeb7a4db9c1ef644473cff2cf4e0e869ee', 'pulse/state-0005.cas.h5': 'c0f9f6efbec002147623a2f2faefcc69f127d4a79d8a0ede1db8af175d22b2cf', 'pulse/state-0005.dat.h5': 'da9b75709dbfc1ecb72078d01adc75a80f4133bd272125ffd5ed2f40106a10f8', 'summary.json': '299ec94b84b1bfc957dd66c9aa584ad70233915b0f5d8b1b1b594b97b797710d', 'pulse/states.json': '3e52618e2f8085d2529c99abcaba2600be441c14071338f539f2e5f1deb8c6e4', 'pulse/static-settings-after.json': '5db86f42c9818f9a1c5f39b1e9292c704714260f0b7036e806055aa4d4bc9f03', 'executed-source.py': '87650ec97332de81f285a04661a7788c4cd15e4048190036f0467dab90dffdf6', 'autosave-actual-metadata.json': '95a78f0b783b79096c9e650094b2d0e3e7c787f61525788856bf4d13b27a815d'}

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--authorization',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);a=parser.parse_args();bounded_driver_source();auth_path=a.authorization.resolve();out=a.out.resolve();auth=json.loads(auth_path.read_text(encoding='utf-8-sig'))
    expected={'approved':True,'approved_by':'root','action':ACTION,'phase':'production-continuation-795-steps','solver_processes':2,'pilot_pulse_steps':5,'new_steps':795,'final_step':800,'dt_s':.025,'max_inner_iterations':60,'both_qualifiers_accepted_by_root':True,'BC_flux_classification_accepted_by_root':True}
    for key,value in expected.items():require(type(auth.get(key)) is type(value) and auth[key]==value,'Authorization differs: '+key)
    require(inside(auth_path,ROOT/'local') and auth['source_sha256']==sha(__file__) and Path(auth['output_directory']).resolve()==out,'Source/output/private authorization differs')
    require(out.parent==ROOT/'native/runs/2026-11-02-tracer-tjunction-v001' and re.fullmatch(r'fluent-production-v\d{3}',out.name) and not out.exists(),'Fresh versioned production output required')
    targets=PRIVATE/'PRODUCTION-TARGETS-v001.md';digest=auth.get('production_targets_sha256','');require(re.fullmatch(r'[0-9a-f]{64}',digest) and targets.is_file() and sha(targets)==digest,'Future ROOT-frozen production targets required')
    acceptance=PRIVATE/'root-Fluent-qualification-acceptance-v001.json';require(sha(acceptance)=='93a4a8c894b1ef3241abaa0508ce4f423405249ff57b04ffe15ad2155977fe6d','Accepted Fluent pilot receipt differs');accepted=json.loads(acceptance.read_text(encoding='utf-8'));require(accepted['independent_all10_checks_passed'] is True and accepted['outlet_convection_zero_diffusion_qualified_by_constant_preservation_and_boundary_value_evidence'] is True,'Fluent pilot/BC interpretation not accepted')
    evidence=auth.get('root_acceptance_evidence',[]);require(isinstance(evidence,list) and evidence,'Both-code/BC/WSL closure root evidence required')
    for row in evidence:
        f=Path(row['path']).resolve();require(inside(f,ROOT/'local') and f.is_file() and sha(f)==row['sha256'],'Root acceptance evidence changed')
    require(sha(QUALIFICATION)==QUALIFICATION_SHA and EXECUTABLE.is_file(),'Original qualification/installed Fluent differs');q=json.loads(QUALIFICATION.read_text(encoding='utf-8'));watched={Path(r['path']):r['sha256_before'] for r in q['source_manifest']}
    watched.update({PILOT/name:expected_sha for name,expected_sha in PILOT_FILES.items()});watched.update({acceptance:sha(acceptance),targets:digest,PRIVATE/'scalar_fv_audit_v002.py':'75ceec60add2f5c53fb82585d105dcc5cd7c9d8d19e1ea8ac3f2dde14aa455cd',PRIVATE/'scalar_native_geometry_v003.py':'06f6b4089116d8e7cdab8b4050c30d2939d3876381ebbaa878c7908e6f2d1903',ROOT/'environments/fluent/Lib/site-packages/ansys/fluent/core/generated/v261/solver/settings.py':'4a314bc418fdd2a9cfefd21b6c7a7286545e3ba7a032eeb1976e9fb6f5678d8e'})
    require(all(sha(f)==v for f,v in watched.items()),'Pinned source/pilot/targets changed before imports');ram=memory_bytes();free=shutil.disk_usage(ROOT).free;require(ram>=4*1024**3 and free>=10*1024**3,'Require4GiB RAM/10GiB E: free');require_idle()
    import numpy as np,h5py
    from ansys.fluent.core.offline import _lispy
    fv=load_module('tracer_production_fv_checked',PRIVATE/'scalar_fv_audit_v002.py');geo=load_module('tracer_production_geometry_checked',PRIVATE/'scalar_native_geometry_v003.py');geometries,geometry_proof=geo.qualified_geometry();reader=sys.modules['qualified_original_tjunction_fluent'];g=geometries['fluent']['stencil'];physical_ids=np.flatnonzero(g['patch']!='interior')+1
    out.mkdir(parents=True);[(out/n).mkdir() for n in ('ancestry','autosave','checkpoints','audits','exports','examples','reports')];temporary=ROOT/'temp/tracer-tjunction'/out.parent.name/out.name;require(not temporary.exists(),'Fresh project temp required');temporary.mkdir(parents=True)
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
    ancestry=[]
    for step in range(6):
        row={'step':step,'ancestry':'accepted pulse pilot; control excluded'}
        for kind in ('cas','dat'):
            src=PILOT/f'pulse/state-{step:04d}.{kind}.h5';dst=out/'ancestry'/src.name;shutil.copyfile(src,dst);require(sha(src)==sha(dst),'Immutable pilot copy differs');row[kind+'_path']=str(dst);row[kind+'_sha256']=sha(dst)
        row['case']=row.pop('cas_path');row['data']=row.pop('dat_path');row['case_sha256']=row.pop('cas_sha256');row['data_sha256']=row.pop('dat_sha256');row['clock']=json.loads((PILOT/'pulse/states.json').read_text(encoding='utf-8'))[step]['clock'];row['native_current_dataset']='cells/SV_UDS_0/1';row['native_M1_dataset']='cells/SV_UDS_0_M1/1'
        ancestry.append(row)
    copied_case=out/'ancestry/state-0005.cas.h5';copied_data=out/'ancestry/state-0005.dat.h5'
    require(sha(copied_case)=='c0f9f6efbec002147623a2f2faefcc69f127d4a79d8a0ede1db8af175d22b2cf' and sha(copied_data)=='da9b75709dbfc1ecb72078d01adc75a80f4133bd272125ffd5ed2f40106a10f8','Actual accepted pulse5 copy differs')
    for src,dst in ((auth_path,out/'root-authorization.json'),(Path(__file__),out/'executed-source.py'),(targets,out/'PRODUCTION-TARGETS-v001.md')):shutil.copyfile(src,dst);require(sha(src)==sha(dst),'Source/authorization/targets copy differs')
    write_json(out/'ancestry.json',ancestry);write_json(out/'source-provenance.json',{'watched':[{'path':str(f),'sha256':v} for f,v in watched.items()],'root_acceptance_evidence':evidence});write_json(out/'qualified-geometry.json',geometry_proof)
    events=[];discovery={};session=None;loaded=False;completed_new_steps=0;audited_through=5;rows=list(ancestry);validated={};optional_errors={};transcript=out/'native-transcript.txt';auto_prefix=out/'autosave/pulse'
    def stage(name,**values):
        record={'stage':name,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),**values};events.append(record);write_json(out/'events.json',events);write_json(out/'status.json',record);print(json.dumps(record,default=json_default),flush=True)
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
    def static_state():
            state=settings_state();state.pop('run_calculation');return state
    def native_report(folder,step):
            command=session.settings.results.report.flow.uds_flow;record={'step':step,'physical_time_s':step*.025,'active':bool(command.is_active()),'argument_names':list(command.argument_names),'domain_allowed_values':command.domain.allowed_values(),'native_flux_semantics':'UNCLASSIFIED; no reconstructed values are labeled native'};before=len(transcript.read_bytes().splitlines()) if transcript.exists() else 0
            try:record['return_value']=command();record['executed']=True
            except Exception:record.update(executed=False,error=traceback.format_exc())
            text=transcript.read_text(encoding='utf-8',errors='replace') if transcript.exists() else '';record['transcript_start_line']=before;record['transcript_end_line']=len(text.splitlines());(folder/f'state-{step:04d}-native-UDS-report.txt').write_text('\n'.join(text.splitlines()[before:])+'\n',encoding='utf-8');write_json(folder/f'state-{step:04d}-native-UDS-report.json',record);return record
    reference=native_fields(copied_case,copied_data);reference_mesh=mesh_fingerprint(copied_case);previous,initial_M1,_=native_scalar(copied_data);initial_clock=native_clock(copied_data);require(initial_clock=={'flow-time':.125,'time-step':5.,'physical-time-step':.025} and np.array_equal(reference['phi'],geometries['fluent']['native_mass_phi']),'Accepted pulse5 native clock/flow identity differs')
    require(np.all(reference['phi'][g['patch']=='outlet']>=0),'Unsupported native outlet backflow')
    def controls_gate(state,top):
        require(state['equations']=={'flow':False,'uds-0':True} and state['UDS_residual']['check_convergence'] is True and state['UDS_residual']['absolute_criteria']==1e-8,'Scalar-only equation/residual controls differ')
        require(state['inlet_diffusion'] is True and state['transient_formulation']=='unsteady-1st-order' and abs(state['dt_s']-.025)<=1e-12 and state['max_inner']==60,'Time/scalar diffusion controls differ')
        require(state['Gamma']['option']=='defined-per-uds' and state['Gamma']['uds_diffusivities']['uds-0']=={'option':'value','value':.001},'Gamma differs')
        for name,value in [('inlet_left',0.),('inlet_top',top),('walls',0.),('outlet',0.)]:
            b=state['scalar_boundaries'][name];require(b['uds_bc']['User Scalar 0']==('Specified Value' if name.startswith('inlet') else 'Specified Flux') and b['uds']['User Scalar 0']['value']==value,'Scalar boundary differs: '+name)
        require(session.settings.solution.monitor.residual.options.n_save.get_state()>=50000,'Native residual history retention below50000')
    def configure_autosave(label):
        auto=session.settings.solution.calculation_activity.auto_save;before=auto.get_state();require('data_frequency' in before,'Active native data_frequency metadata required')
        auto.case_frequency='each-time';auto.retain_most_recent_files=False;auto.root_name=str(auto_prefix);auto.save_data_file_every.frequency_type='time-step';auto.append_file_name_with.file_suffix_type='time-step';auto.data_frequency=1
        after=auto.get_state();native_root=Path(after['root_name']);native_root=native_root if native_root.is_absolute() else out/native_root
        require(native_root.resolve()==auto_prefix.resolve() and after['case_frequency']=='each-time' and after['retain_most_recent_files'] is False and after['data_frequency']==1 and after['save_data_file_every']['frequency_type']=='time-step' and after['append_file_name_with']['file_suffix_type']=='time-step','Actual canonical all-state autosave settings differ')
        write_json(out/'reports'/(label+'-autosave-readback.json'),{'before':before,'after':after,'canonical_prefix':str(auto_prefix),'data_frequency_used':True,'inactive_save_frequency_not_set':True});return after
    def history_index(path):
        with h5py.File(path,'r') as h:k=np.asarray(h['results/residuals/phase-1/uds-0/iterations'][()]).ravel()
        require(k.size and np.array_equal(k,k.astype(int)) and np.all(np.diff(k)>0),'Actual native residual indices invalid');return int(k[-1])
    previous_iteration=history_index(copied_data)
    def printed_residuals():
        header=None;printed={}
        for line in transcript.read_text(encoding='utf-8',errors='replace').splitlines():
            words=line.split()
            if words and words[0]=='iter' and 'uds-0' in words:header=words;continue
            if header and words and words[0].isdigit():
                try:printed[int(words[0])]=float(words[header.index('uds-0')])
                except (ValueError,IndexError):pass
        return printed
    def own_residual(path,before,printed):
        with h5py.File(path,'r') as h:
            group=h['results/residuals/phase-1/uds-0'];k=np.asarray(group['iterations'][()]).ravel();require(k.size and np.array_equal(k,k.astype(int)) and np.all(np.diff(k)>0),'Native residual index invalid');start=int(np.searchsorted(k,before+1));indices=k[start:].astype(int);values=np.asarray(group['data'][start:,:2])
        require(1<=len(indices)<=60 and int(indices[-1])-before<=60 and np.isfinite(values).all(),'Own new residual rows/60-inner cap fail');num=float(values[-1,0]);den=float(values[-1,1]);require(num>=0 and den>=0,'Negative native new residual numerator/denominator')
        if den>0:value=num/den;method='native new numerator/positive denominator'
        else:require(int(indices[-1]) in printed and float(format(num,'.5g'))==printed[int(indices[-1])],'Zero normalization lacks own printed native proof');value=printed[int(indices[-1])];method='actual native printed value; ratio undefined because denominator0'
        require(np.isfinite(value) and 0<=value<=1e-8,'Own new terminal UDS residual target fails');return {'last_iteration':int(indices[-1]),'new_indices':indices.tolist(),'new_numerator_denominator_rows':values.tolist(),'terminal_scaled_residual':float(value),'method':method,'terminal_numerator':num,'terminal_denominator':den,'all_older_indices_excluded_at_or_below':before}
    def discover_autosaves(expected):
        found={};all_data=list((out/'autosave').glob('*.dat.h5'))
        for dp in all_data:
            if str(dp) in validated:continue
            cp=dp.with_name(dp.name.removesuffix('.dat.h5')+'.cas.h5');require(cp.is_file(),'Each native autosave must retain matching CFF case/data basename');clock=native_clock(dp);step=clock['time-step'];require(step==int(step) and int(step) in expected and int(step) not in found,'Unexpected/duplicate native autosave clock');found[int(step)]=(cp,dp,clock)
        require(set(found)==set(expected),'Missing actual native per-step autosaved state')
        for name,digest in validated.items():require(sha(name)==digest,'Earlier native autosave overwritten')
        return found
    def audit_batch(expected,top):
        nonlocal previous,previous_iteration,audited_through
        found=discover_autosaves(expected);printed=printed_residuals();require(expected[0]==audited_through+1,'Nonconsecutive continuation audit')
        for step in expected:
            cp,dp,clock=found[step];require(abs(clock['flow-time']-step*.025)<=1e-12 and abs(clock['physical-time-step']-.025)<=1e-12,'Native autosave time/dt mismatch');flow=native_fields(cp,dp);changes=fv.frozen_field_changes(reference,flow);mesh=mesh_fingerprint(cp);require(all(v<=1e-12 for v in changes.values()) and mesh['sha256']==reference_mesh['sha256'],'Frozen native flow/mesh changed')
            current,M1,face_c=native_scalar(dp);require(np.max(abs(M1-previous))<=1e-12,'Native M1 differs from actual preceding current');require(float(current.min())>=-1e-8 and float(current.max())<=1+1e-8,'Native scalar bounds fail')
            top_error=float(np.max(abs(face_c[g['patch']=='inlet_top']-top)));left_error=float(np.max(abs(face_c[g['patch']=='inlet_left'])));require(top_error<=1e-12 and left_error<=1e-12,'Actual solved native inlet scalar boundary fails scheduled value')
            boundary={'inlet_left':{'value':0.,'diffusion':'dirichlet','evidence':'Accepted actual native Specified Value0 and inlet diffusion on; RECONSTRUCTED'},'inlet_top':{'value':top,'diffusion':'dirichlet','evidence':'Actual scheduled native boundary field and Specified Value; RECONSTRUCTED'},'walls':{'diffusion':'zero_gradient','evidence':'Root-qualified constant-c1/early-pulse zero diffusion interpretation; RECONSTRUCTED'},'outlet':{'diffusion':'zero_gradient','evidence':'Root-qualified outward convection/zero diffusion from constant-c1 control; RECONSTRUCTED'}}
            terms,defect,budget=fv.audit_step(g,previous,current,reference['phi'],rho=reference['rho'],gamma=.001,dt_s=.025,boundary=boundary,normalization_time_s=2.,global_floor_kg_m_s=.02);require(budget['global_Euler_defect_relative']<=1e-5,'Reconstructed global Euler budget target fails');residual=own_residual(dp,previous_iteration,printed)
            row={'step':step,'case':str(cp),'data':str(dp),'case_sha256':sha(cp),'data_sha256':sha(dp),'clock':clock,'native_current_dataset':'cells/SV_UDS_0/1','native_M1_dataset':'cells/SV_UDS_0_M1/1','M1_previous_current_max_abs':float(np.max(abs(M1-previous))),'frozen_native_field_max_abs':changes,'all_mesh_datasets_sha256':mesh['sha256'],'scheduled_top_value':top,'actual_top_boundary_max_abs_error':top_error,'actual_left_boundary_max_abs_error':left_error,'own_new_residual':residual,'RECONSTRUCTED_budget':budget}
            np.savez_compressed(out/'audits'/f'state-{step:04d}-native-fields.npz',c=current,M1=M1,boundary_c=face_c);np.savez_compressed(out/'audits'/f'state-{step:04d}-RECONSTRUCTED-fluxes.npz',**terms,reconstructed_cell_defect_kg_s=defect);write_json(out/'audits'/f'state-{step:04d}.json',row)
            with (out/'state-index.jsonl').open('a',encoding='utf-8') as stream:stream.write(json.dumps(row,allow_nan=False)+'\n')
            validated[str(cp)]=row['case_sha256'];validated[str(dp)]=row['data_sha256'];rows.append(row);previous=current;previous_iteration=residual['last_iteration'];audited_through=step
            if step%50==0:stage('native_saved_states_audited',audited_through=step,completed_new_step_calls=completed_new_steps)
        write_json(out/'states.json',rows);write_json(out/'reports'/f'autosave-coverage-through-{expected[-1]:04d}.json',{'actual_steps':list(expected),'native_files':[{'step':s,'case':str(found[s][0]),'data':str(found[s][1]),'clock':found[s][2]} for s in expected],'naming_inferred_from_native_clock_not_padding':True});return found
    def manual_checkpoint(label,step,top):
        cp=out/'checkpoints'/(label+'.cas.h5');dp=cp.with_name(label+'.dat.h5');session.settings.file.write_case_data(file_name=str(cp));clock=native_clock(dp);current,M1,face_c=native_scalar(dp);require(clock['time-step']==step and abs(clock['flow-time']-step*.025)<=1e-12 and np.array_equal(current,previous),'Checkpoint changed physical state');flow=native_fields(cp,dp);require(all(v<=1e-12 for v in fv.frozen_field_changes(reference,flow).values()) and mesh_fingerprint(cp)['sha256']==reference_mesh['sha256'],'Checkpoint changed frozen flow/mesh');small=scalar_controls();controls_gate(small,top);write_json(out/'checkpoints'/(label+'-readback.json'),{'clock':clock,'case_sha256':sha(cp),'data_sha256':sha(dp),'scalar_controls':small,'native_top_face_min':float(face_c[g['patch']=='inlet_top'].min()),'native_top_face_max':float(face_c[g['patch']=='inlet_top'].max()),'unsolved_faces_not_new_BC_enforcement_proof':label=='top-off-before-step81'});configure_autosave('restored-after-'+label)
    try:
        import ansys.fluent.core as pf
        from ansys.fluent.core.utils.fluent_version import FluentVersion
        require(pf.__version__=='0.43.0','PyFluent changed');pf.config.examples_path=str(out/'examples');ram_now=memory_bytes();free_now=shutil.disk_usage(ROOT).free;require(ram_now>=4*1024**3 and free_now>=10*1024**3,'Fresh launch resource gate fails');require_idle();stage('launching_one_hidden_continuation',new_step_budget=795,available_ram_bytes=ram_now,free_E_bytes=free_now)
        session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True);session.transcript.start(file_name=str(transcript));version=session.get_fluent_version();write_json(out/'actual-version.json',{'type':type(version).__qualname__,'str':str(version),'value':version.value,'name':version.name});require(version==FluentVersion.v261,'Exact Fluent261 required');st=session.settings;st.file.read_case_data(file_name=str(copied_case));loaded=True;write_json(out/'loaded-settings.json',settings_state())
        st.solution.monitor.residual.options.n_save=50000;controls_gate(scalar_controls(),1.);configure_autosave('before-step6');static_on=static_state();write_json(out/'static-top-on-before.json',static_on)
        stage('advancing_native_alignment_step6',new_steps=1);st.solution.run_calculation.dual_time_iterate(time_step_count=1,max_iter_per_step=60);completed_new_steps+=1;audit_batch([6],1.);require(static_state()==static_on,'Static settings changed at native step6');write_json(out/'reports/native-step6-autosave-qualified.json',{'case_and_data_native_clock_step6_verified':True,'actual_clock':rows[-1]['clock'],'actual_case':rows[-1]['case'],'actual_data':rows[-1]['data'],'all_numerical_gates_passed':True,'bulk_remaining_not_yet_executed':True})
        require(completed_new_steps==1 and audited_through==6,'Alignment stage coverage differs');configure_autosave('before-bulk74');stage('advancing_bulk74_to_step80',new_steps=74);st.solution.run_calculation.dual_time_iterate(time_step_count=74,max_iter_per_step=60);completed_new_steps+=74;audit_batch(list(range(7,81)),1.);require(static_state()==static_on,'Static top-on settings changed during bulk74');manual_checkpoint('top-on-step80',80,1.)
        set_boundary('inlet_top',st.setup.boundary_conditions.velocity_inlet['inlet_top'],{'specifiedvalue'},0.);controls_gate(scalar_controls(),0.);manual_checkpoint('top-off-before-step81',80,0.);static_off=static_state();write_json(out/'static-top-off-before.json',static_off)
        require(completed_new_steps==75 and audited_through==80,'Top-on stage budget differs');configure_autosave('before-step81');stage('advancing_actual_switched_step81',new_steps=1);st.solution.run_calculation.dual_time_iterate(time_step_count=1,max_iter_per_step=60);completed_new_steps+=1;audit_batch([81],0.);require(static_state()==static_off,'Static settings changed at switched step81');write_json(out/'reports/actual-step80-81-switch-proof.json',{'step80_top_c':1.,'step81_top_c':0.,'step80':rows[80],'step81':rows[81],'actual_boundary_field_proof':True,'future_continuum_edge_not_inferred_from_discrete_schedule':True})
        require(completed_new_steps==76 and audited_through==81,'Switch stage budget differs');configure_autosave('before-bulk719');stage('advancing_bulk719_to_step800',new_steps=719);st.solution.run_calculation.dual_time_iterate(time_step_count=719,max_iter_per_step=60);completed_new_steps+=719;audit_batch(list(range(82,801)),0.);require(static_state()==static_off,'Static top-off settings changed during bulk719');controls_gate(scalar_controls(),0.);manual_checkpoint('final-step800',800,0.)
        try:
            graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1600;graphics.picture.y_resolution=1000;graphics.picture.invert_background=False;graphics.contour.create(name='tracer-c');contour=graphics.contour['tracer-c'];contour.field='uds-0-scalar';contour.surfaces_list=[];contour.options.filled=True;contour.options.node_values=False;contour.range_options.auto_range=False;contour.range_options.minimum=0.;contour.range_options.maximum=1.;contour.display();graphics.views.camera.projection='orthographic';graphics.views.camera.position=[.4,.3,10.];graphics.views.camera.target=[.4,.3,0.];graphics.views.camera.field=[1.8,1.125];graphics.picture.save_picture(file_name=str(out/'exports/native-scalar-c-final.png'));write_json(out/'exports/native-contour-settings.json',contour.get_state())
        except Exception:optional_errors['native_contour']=traceback.format_exc()
        require(completed_new_steps==795 and audited_through==800 and len(rows)==801,'Incomplete exact continuation coverage');write_json(out/'summary.json',{'production_driver_completed':True,'bound_accepted_pilot_steps':5,'new_native_step_calls_completed':795,'actual_final_step':800,'actual_final_time_s':20.,'actual_native_state_count':801,'all_frozen_numerical_gates_passed':True,'actual_step80_top1_step81_top0_native_boundary_checked':True,'native_autosave_all_state_naming_alignment_qualified':True,'residual_history_n_save':50000,'RECONSTRUCTED_fluxes_are_not_native_matrix_flux':True,'optional_errors':optional_errors,'root_final_acceptance_pending':True,'physical_validation':False});stage('full_continuation_saved_for_root_review',new_native_steps=795,final_step=800,root_final_acceptance_pending=True)
    except Exception:
        write_json(out/'failure.json',{'error':traceback.format_exc(),'completed_new_step_calls':completed_new_steps,'audited_through':audited_through,'native_failure_clock_authoritative_if_RPC_interrupted':True,'automatic_replay_or_extension':False});write_json(out/'states.json',rows)
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
    print(json.dumps({'out':str(out),'new_steps':completed_new_steps,'audited_through':audited_through,'root_final_acceptance_pending':True}))

if __name__=='__main__':main()
