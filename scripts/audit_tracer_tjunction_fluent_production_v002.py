"""Read-only actual Fluent tracer 801-state auditor/exporter. No CFD, RPC, WSL or case writes.
Run only after root observes completed production and native session/process closure.
An explicit hash-bound root recovery receipt can accept saved computation after interrupted postprocessing.
The original five pulse steps remain ancestry. All scalar face fluxes below are RECONSTRUCTED.
"""
from pathlib import Path
import argparse,ast,ctypes,datetime,hashlib,importlib.util,json,re,shutil,sys
import h5py,numpy as np
ROOT=Path('E:/CFD-Learning-Lab');PRIVATE=ROOT/'local/prospective-tracer-v001'
DRIVER=PRIVATE/'tracer_tjunction_fluent_production_v002.py'
DRIVER_SHA='9b10b7f3cc016153eec6be7cba08569cbf0ddb56788a25fccca789226210c42c'
TARGET_SHA='07f5b0745fc71758b9cf5678171ce17fd558f00a01cd773560911920e91318ba'
PILOT=ROOT/'native/runs/2026-11-02-tracer-tjunction-v001/fluent-qualification-v001'
ACCEPT_SHA='93a4a8c894b1ef3241abaa0508ce4f423405249ff57b04ffe15ad2155977fe6d'
FV_SHA='75ceec60add2f5c53fb82585d105dcc5cd7c9d8d19e1ea8ac3f2dde14aa455cd'
GEO_SHA='06f6b4089116d8e7cdab8b4050c30d2939d3876381ebbaa878c7908e6f2d1903'
FLOAT=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def require(v,message):
 if not v:raise ValueError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
def inside(path,parent):
 try:Path(path).resolve().relative_to(Path(parent).resolve());return True
 except ValueError:return False
def json_default(value):
 if isinstance(value,np.generic):return value.item()
 if isinstance(value,np.ndarray):return value.tolist()
 if isinstance(value,Path):return str(value)
 if isinstance(value,bytes):return {'bytes_hex':value.hex()}
 raise TypeError(type(value).__name__)
def write_json(path,value):
 with Path(path).open('x',encoding='utf-8',newline='\n') as f:json.dump(value,f,indent=2,allow_nan=False,default=json_default);f.write('\n')
class Evidence:
 def __init__(self):self.before={}
 def touch(self,path,expected=None):
  path=Path(path).resolve();require(path.is_file(),'Missing actual evidence '+str(path));digest=sha(path)
  require(expected is None or digest==expected,'Changed pinned input '+str(path))
  require(path not in self.before or self.before[path]==digest,'Input changed between reads '+str(path));self.before[path]=digest;return path
 def text(self,path,expected=None):return self.touch(path,expected).read_text(encoding='utf-8-sig')
 def json(self,path,expected=None):return json.loads(self.text(path,expected))
 def finish(self):
  rows=[dict(path=str(p),sha256=h,bytes=p.stat().st_size) for p,h in self.before.items()]
  require(all(sha(p)==h for p,h in self.before.items()),'Consumed source/native input changed during read');return rows
def checked_load(name,path,expected,watched):
 path=Path(path).resolve();require(sha(path)==expected,'Pre-import source hash '+str(path));watched[path]=expected
 require(name not in sys.modules,'Checked module namespace already loaded')
 spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);sys.modules[name]=module;spec.loader.exec_module(module);return module
def check_clock(clock,step):
 require(set(clock)=={'flow-time','time-step','physical-time-step'} and all(np.isfinite(list(clock.values()))),'Finite complete native clock')
 require(abs(clock['flow-time']-step*.025)<=1e-12 and clock['time-step']==step and abs(clock['physical-time-step']-.025)<=1e-12,'Actual native clock/time/index/dt differs')
def read_clock(path):
 with h5py.File(path,'r') as h:raw=h['settings/Data Variables'][0].decode()
 result={}
 for key in ('flow-time','time-step','physical-time-step'):
  found=re.findall(r'\('+re.escape(key)+r'\s+('+FLOAT+r')\)',raw);require(len(found)==1,'Unambiguous native clock '+key);result[key]=float(found[0])
 return result
def validate_rows(rows,run):
 require(isinstance(rows,list) and len(rows)==801,'Exactly801 native pulse states required')
 cases=set();data=set()
 for step,row in enumerate(rows):
  require(type(row.get('step')) is int and row['step']==step,'Missing/duplicate/out-of-order actual state')
  check_clock(row['clock'],step)
  for key,seen,suffix in (('case',cases,'.cas.h5'),('data',data,'.dat.h5')):
   path=Path(row[key]).resolve();parent=Path(run)/('ancestry' if step<=5 else 'autosave')
   require(inside(path,parent) and path.name.endswith(suffix) and path not in seen,'Native state path escape/duplicate/kind');seen.add(path)
   require(re.fullmatch(r'[0-9a-f]{64}',row[key+'_sha256']) is not None,'Actual native state hash required')
 return rows
def mesh_equal(reference,current):
 with h5py.File(reference,'r') as a,h5py.File(current,'r') as b:
  names=[];other=[];a['meshes/1'].visititems(lambda n,v:names.append(n) if isinstance(v,h5py.Dataset) else None);b['meshes/1'].visititems(lambda n,v:other.append(n) if isinstance(v,h5py.Dataset) else None)
  require(names==other,'Native mesh dataset coverage differs')
  for name in names:
   x=a['meshes/1/'+name];y=b['meshes/1/'+name]
   require(x.dtype==y.dtype and x.shape==y.shape and np.array_equal(x[()],y[()]),'Native mesh dataset values differ '+name)
   require(set(x.attrs)==set(y.attrs) and all(np.array_equal(x.attrs[k],y.attrs[k]) for k in x.attrs),'Native mesh dataset attributes differ '+name)
 return len(names)
def rounded_tolerance(token):
 token=token.lower();mant,exponent=token.split('e') if 'e' in token else (token,'0');digits=len(mant.split('.')[1]) if '.' in mant else 0
 return .500001*10**(int(exponent)-digits)
def printed_residuals(text):
 header=None;result={}
 for line in text.splitlines():
  words=line.split()
  if words and words[0]=='iter' and 'uds-0' in words:header=words;continue
  if header and words and words[0].isdigit():
   try:token=words[header.index('uds-0')];float(token);result[int(words[0])]=token
   except (ValueError,IndexError):pass
 return result
def new_residual(path,before,printed):
 with h5py.File(path,'r') as h:
  r=h['results/residuals/phase-1/uds-0'];raw=np.asarray(r['iterations'][()]).ravel();require(raw.size and np.array_equal(raw,raw.astype(int)) and np.all(np.diff(raw)>0),'Native residual indices invalid')
  start=int(np.searchsorted(raw,before+1));indices=raw[start:].astype(int);values=np.asarray(r['data'][start:,:2])
 require(1<=len(indices)<=60 and np.array_equal(indices,np.arange(before+1,int(indices[-1])+1)) and values.shape==(len(indices),2) and np.isfinite(values).all(),'Own new contiguous native residual rows/cap')
 num,den=map(float,values[-1]);require(num>=0 and den>=0,'Native residual signs');last=int(indices[-1]);require(last in printed,'Own terminal native printed normalization proof missing');token=str(printed[last]);native=float(token)
 value=num/den if den>0 else native
 require((den>0 and abs(value-native)<=rounded_tolerance(token)) or (den==0 and float(format(num,'.5g'))==native),'Native residual numerator/denominator vs print normalization mismatch')
 require(np.isfinite(value) and 0<=value<=1e-8,'Own new native UDS scaled residual target')
 return dict(last_iteration=last,new_indices=indices.tolist(),terminal_numerator=num,terminal_denominator=den,terminal_scaled_residual=float(value),printed_native_token=token,normalization='native numerator/positive denominator' if den>0 else 'own native printed value; zero denominator ratio undefined',older_indices_excluded_at_or_below=before)
def native_fields(case,data,reader):
 with h5py.File(data,'r') as h:
  cells=h['results/1/phase-1/cells'];fields={key:reader.chunk_array(cells[group],4400) for key,group in [('U','SV_U'),('V','SV_V'),('p','SV_P'),('rho','SV_DENSITY'),('mu','SV_MU_LAM')]}
 c0,c1,phi,_=reader.native_faces(case,data);fields.update(phi=phi,c0=c0,c1=c1)
 require(all(np.isfinite(v).all() for v in fields.values()),'Finite actual native frozen fields');return fields
def native_scalar(data,reader,physical):
 with h5py.File(data,'r') as h:
  cells=h['results/1/phase-1/cells'];c=reader.chunk_array(cells['SV_UDS_0'],4400);old=reader.chunk_array(cells['SV_UDS_0_M1'],4400);face=np.full(9040,np.nan);seen=np.zeros(9040,bool);datasets={}
  for group in ('SV_UDS_0','SV_UDS_0_M1'):datasets[group]=['results/1/phase-1/cells/'+group+'/'+name for name in cells[group]]
  for ds in h['results/1/phase-1/faces/SV_UDS_0'].values():
   lo=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;hi=int(np.asarray(ds.attrs['maxId']).ravel()[0]);v=np.asarray(ds[()]).ravel()
   require(0<=lo<hi<=9040 and len(v)==hi-lo and not seen[lo:hi].any(),'Native physical scalar face chunk cover');face[lo:hi]=v;seen[lo:hi]=True
  require(np.array_equal(seen,physical) and np.isfinite(face[seen]).all(),'Observed480 native physical face scalar coverage')
 require(np.isfinite(c).all() and np.isfinite(old).all(),'Finite actual current/M1 scalar')
 return c,old,face,datasets
def scalar_step_checks(c,old,previous,face,g,step):
 require(np.isfinite(c).all() and np.isfinite(old).all(),'Finite actual current/M1')
 row=dict(c_min=float(c.min()),c_max=float(c.max()),negative_cells=int((c<0).sum()),above_one_cells=int((c>1).sum()),M1_range=[float(old.min()),float(old.max())])
 if step==0:require(np.max(abs(c))<=1e-12,'Actual initial cell scalar must be0');row['initial_M1_is_not_advanced_proof']=True;return row
 require(c.min()>=-1e-8 and c.max()<=1+1e-8,'Actual scalar bounds')
 row['M1_previous_current_max_abs']=float(abs(old-previous).max());require(row['M1_previous_current_max_abs']<=1e-12,'Native M1 vs preceding actual current')
 row['scheduled_top_value']=1. if step<=80 else 0.
 row['actual_top_boundary_max_abs_error']=float(abs(face[g['patch']=='inlet_top']-row['scheduled_top_value']).max());row['actual_left_boundary_max_abs_error']=float(abs(face[g['patch']=='inlet_left']).max())
 require(row['actual_top_boundary_max_abs_error']<=1e-12 and row['actual_left_boundary_max_abs_error']<=1e-12,'Actual solved native inlet scalar schedule')
 row['outlet_wall_zero_gradient_value_error_descriptive']={n:float(abs(face[g['patch']==n]-c[g['owner'][g['patch']==n]]).max()) for n in ('outlet','walls')}
 return row
def close_budget(actual,stored):
 if isinstance(actual,dict):
  require(isinstance(stored,dict) and set(actual)==set(stored),'Producer reconstructed budget keys differ')
  for key in actual:close_budget(actual[key],stored[key])
 elif isinstance(actual,(float,int)) and not isinstance(actual,bool):
  require(isinstance(stored,(float,int)) and not isinstance(stored,bool) and abs(actual-stored)<=1e-12+1e-10*abs(actual),'Producer reconstructed budget stale numeric')
 elif isinstance(actual,list):
  require(isinstance(stored,list) and len(actual)==len(stored),'Producer budget list')
  for a,b in zip(actual,stored):close_budget(a,b)
 elif isinstance(actual,str):require(isinstance(stored,str),'Producer budget text type')
 else:require(actual==stored,'Producer reconstructed budget flag/null differs')
def validate_recovery(path,expected,run,e):
 path=Path(path).resolve();require(inside(path,PRIVATE) and re.fullmatch(r'[0-9a-f]{64}',expected or ''),'Private root recovery receipt and external hash required');record=e.json(path,expected)
 require(record.get('schema_version')==1 and record.get('reviewer')=='root' and record.get('action')=='accept_saved_fluent_tracer_computation_for_offline_audit','Explicit root saved-computation recovery approval')
 require(Path(record['native_run']).resolve()==Path(run).resolve() and record['driver_sha256']==DRIVER_SHA and record['production_targets_sha256']==TARGET_SHA,'Recovery run/source/targets differ')
 for key in ('native_step800_time20_saved_before_close','fluent_session_closed','postprocessor_stopped'):require(type(record.get(key)) is bool and record[key] is True,'Recovery must prove '+key)
 for key in ('closure','postprocessor_stop_proof','interrupted_status','interrupted_states'):
  item=record[key];require(inside(item['path'],ROOT/'local') if key in ('closure','postprocessor_stop_proof') else Path(item['path']).resolve()==(Path(run)/('status.json' if key=='interrupted_status' else 'states.json')).resolve(),'Recovery proof path '+key);e.touch(item['path'],item['sha256'])
 if 'amendment' in record:
  amendment=record['amendment'];require(Path(amendment['path']).resolve()==PRIVATE/'ARCHIVE-RECOVERY-AMENDMENT-v001.md' and amendment['sha256']=='d1fe61441f359a9a0a536bb6713bc047f35eff731b35a88ff4add65663e7edf9','Root archive amendment metadata differs');e.touch(amendment['path'],amendment['sha256'])
 closure=e.json(record['closure']['path'],record['closure']['sha256'])
 for key in ('saved_step800_time20_before_close','native_fields_identical_to_autosave800','session_exit_returned'):require(type(closure.get(key)) is bool and closure[key] is True,'Actual closure proof lacks '+key)
 e.json(record['postprocessor_stop_proof']['path'],record['postprocessor_stop_proof']['sha256'])
 final=record['saved_step800'];check_clock(final['clock'],800)
 for key in ('case','data'):
  require(inside(final[key],run),'Final saved state escaped actual run');e.touch(final[key],final[key+'_sha256'])
 return record
def discover_saved_pairs(run,e):
 found={}
 for directory,first,last in (('ancestry',0,5),('autosave',6,800)):
  for data in sorted((Path(run)/directory).glob('*.dat.h5')):
   data=e.touch(data);case=data.with_name(data.name.removesuffix('.dat.h5')+'.cas.h5');e.touch(case);clock=read_clock(data);step=clock['time-step'];require(step==int(step) and first<=step<=last and int(step) not in found,'Recovery unexpected/duplicate native clock');step=int(step);check_clock(clock,step)
   found[step]=dict(step=step,case=str(case),data=str(data),case_sha256=e.before[case.resolve()],data_sha256=e.before[data.resolve()],clock=clock)
 require(set(found)==set(range(801)),'Recovery requires all801 actual native pairs');return validate_rows([found[i] for i in range(801)],run)
def recovery_producer_records(run,rows,partial,e):
 require(isinstance(partial,list) and 6<=len(partial)<=801,'Preserved partial producer ledger required')
 known={}
 def bind(record):
  step=record['step'];require(type(step) is int and 0<=step<=800,'Partial producer state index');actual=rows[step]
  for key in ('case','data'):require(Path(record[key]).resolve()==Path(actual[key]).resolve() and record[key+'_sha256']==actual[key+'_sha256'],'Partial producer/native path/hash mismatch')
  require(record['clock']==actual['clock'],'Partial producer/native clock mismatch');return step
 for expected,record in enumerate(partial):
  step=bind(record);require(step==expected and step not in known,'Preserve consecutive unmodified partial ledger');known[step]=record
 result={}
 for path in sorted((Path(run)/'audits').glob('state-*.json')):
  record=e.json(path);step=bind(record);require(6<=step<=800 and path.name==f'state-{step:04d}.json' and step not in result,'Unexpected producer audit JSON')
  if step in known:require(record==known[step],'Producer JSON differs from preserved partial ledger')
  require('RECONSTRUCTED_budget' in record,'Available producer audit lacks its actual budget');result[step]=record
 require(set(range(6,len(partial))).issubset(result),'Partial ledger claims missing producer audit JSON');return result
def final_recovery_readback(record,e,reader,g,baseline,basecase,expected_current,expected_M1):
 saved=record['saved_step800'];case=e.touch(saved['case'],saved['case_sha256']);data=e.touch(saved['data'],saved['data_sha256']);clock=read_clock(data);check_clock(clock,800);require(clock==saved['clock'],'Root final saved/native clock differs');count=mesh_equal(basecase,case);require(count==32,'Root final saved actual32 mesh datasets')
 flow=native_fields(case,data,reader);changes={k:float(abs(v-flow[k]).max()) for k,v in baseline.items()};require(max(changes.values())<=1e-12,'Final saved native frozen fields differ')
 c,M1,face,datasets=native_scalar(data,reader,g['patch']!='interior');facts=scalar_step_checks(c,M1,expected_M1,face,g,800);difference=float(abs(c-expected_current).max());require(difference<=1e-12,'New final save differs from actual autosave800')
 return dict(case=str(case),data=str(data),case_sha256=saved['case_sha256'],data_sha256=saved['data_sha256'],clock=clock,mesh_dataset_count=count,current_vs_autosave800_max_abs=difference,M1_vs_actual799_max_abs=facts['M1_previous_current_max_abs'],frozen_native_changes=changes,observed_scalar_datasets=datasets)
def available_ram():
 class MS(ctypes.Structure):_fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ('total','available','total_page','available_page','total_virtual','available_virtual','available_extended')]
 m=MS();m.length=ctypes.sizeof(m);require(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)),'Windows RAM query');return int(m.available)
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--recovery',type=Path);parser.add_argument('--recovery-sha256');args=parser.parse_args();run=args.run.resolve();out=args.out.resolve()
 require(run.parent==ROOT/'native/runs/2026-11-02-tracer-tjunction-v001' and re.fullmatch(r'fluent-production-v\d{3}',run.name),'Versioned actual production path')
 require(inside(out,PRIVATE) and out.parent==PRIVATE and re.fullmatch(r'fluent-production-checks-v\d{3}',out.name) and not out.exists(),'Fresh private audit output required')
 ram=available_ram();require(ram>=2.5*1024**3,'Fresh readonly analysis RAM floor2.5GiB');e=Evidence();e.touch(__file__);e.touch(DRIVER,DRIVER_SHA);e.touch(run/'executed-source.py',DRIVER_SHA)
 e.touch(PRIVATE/'PRODUCTION-TARGETS-v001.md',TARGET_SHA);e.touch(run/'PRODUCTION-TARGETS-v001.md',TARGET_SHA);acceptance=e.json(PRIVATE/'root-Fluent-qualification-acceptance-v001.json',ACCEPT_SHA)
 require(acceptance['independent_all10_checks_passed'] is True and acceptance['outlet_convection_zero_diffusion_qualified_by_constant_preservation_and_boundary_value_evidence'] is True,'Accepted pilot/BC interpretation')
 require((args.recovery is None)==(args.recovery_sha256 is None),'Recovery path/hash must be supplied together');recovery=validate_recovery(args.recovery,args.recovery_sha256,run,e) if args.recovery is not None else None
 if recovery is None:
  status=e.json(run/'status.json');summary=e.json(run/'summary.json');require(status['stage']=='full_continuation_saved_for_root_review','Actual production completion required')
  for key,want in {'production_driver_completed':True,'bound_accepted_pilot_steps':5,'new_native_step_calls_completed':795,'actual_final_step':800,'actual_native_state_count':801,'all_frozen_numerical_gates_passed':True,'actual_step80_top1_step81_top0_native_boundary_checked':True,'residual_history_n_save':50000}.items():require(type(summary.get(key)) is type(want) and summary[key]==want,'Actual production summary differs '+key)
  require(summary['physical_validation'] is False and abs(summary['actual_final_time_s']-20.)<=1e-12,'Actual final schedule/limitations')
  require(not (run/'failure.json').exists() and not (run/'session-exit-error.json').exists(),'Failed/incomplete native session preserved; refuse accepted full export')
  require(e.json(run/'immutable-source-final-check.json')['all_unchanged'] is True,'Producer finally/source preservation incomplete')
 else:
  status=e.json(run/'status.json');summary=e.json(run/'summary.json') if (run/'summary.json').exists() else None
  # Preserve raw interrupted producer status/summary exactly; do not synthesize its completion/final-check.
  if (run/'immutable-source-final-check.json').exists():e.json(run/'immutable-source-final-check.json')
 auth=e.json(run/'root-authorization.json');require(auth['approved'] is True and auth['approved_by']=='root' and auth['source_sha256']==DRIVER_SHA and auth['both_qualifiers_accepted_by_root'] is True and auth['BC_flux_classification_accepted_by_root'] is True and Path(auth['output_directory']).resolve()==run,'Actual root production admission binding')
 provenance=e.json(run/'source-provenance.json')
 for row in provenance['watched']:e.touch(row['path'],row['sha256'])
 for row in provenance['root_acceptance_evidence']:require(inside(row['path'],ROOT/'local'),'Private root evidence path');e.touch(row['path'],row['sha256'])
 tree=ast.parse(e.text(DRIVER));pilot_pins=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PILOT_FILES' for t in n.targets))
 for name,digest in pilot_pins.items():e.touch(PILOT/name,digest)
 partial_states=e.json(run/'states.json');rows=validate_rows(partial_states,run) if recovery is None else discover_saved_pairs(run,e);producer_records={} if recovery is None else recovery_producer_records(run,rows,partial_states,e);missing_producer_budgets=[];pilot_rows=e.json(PILOT/'pulse/states.json',pilot_pins['pulse/states.json'])
 fv=checked_load('production_independent_fv',PRIVATE/'scalar_fv_audit_v002.py',FV_SHA,e.before);geo=checked_load('production_independent_geometry',PRIVATE/'scalar_native_geometry_v003.py',GEO_SHA,e.before);geos,gproof=geo.qualified_geometry();reader=sys.modules['qualified_original_tjunction_fluent'];g=geos['fluent']['stencil']
 for source in gproof['consumed_inputs']:e.touch(source['path'],source['sha256'])
 physical=g['patch']!='interior';mapping=np.asarray(geos['fluent']['native_cell_to_source_cell'],dtype=int);require(np.array_equal(np.sort(mapping),np.arange(4400)),'Native/source cell bijection');source_to_native=np.argsort(mapping)
 basecase=e.touch(PILOT/'pulse/state-0000.cas.h5',pilot_pins['pulse/state-0000.cas.h5']);basedata=e.touch(PILOT/'pulse/state-0000.dat.h5',pilot_pins['pulse/state-0000.dat.h5']);baseline=native_fields(basecase,basedata,reader)
 require(np.array_equal(baseline['phi'],geos['fluent']['native_mass_phi']),'Native all-face phi geometry binding');require(np.all(baseline['phi'][g['patch']=='outlet']>=0),'Unsupported actual outlet backflow')
 pilot_text=e.text(PILOT/'native-transcript.txt');production_printed=printed_residuals(e.text(run/'native-transcript.txt'));pilot_lines=pilot_text.splitlines()
 # Source-order and native-order scalar arrays are preallocated once (~54MiB together).
 c_native=np.empty((801,4400));c_source=np.empty_like(c_native);times=np.empty(801);inventory=np.empty(801);residual=np.empty(800);iteration=np.empty(800,dtype=np.int64)
 scalar_min=np.empty(801);scalar_max=np.empty(801);negatives=np.empty(801,dtype=np.int64);overshoots=np.empty(801,dtype=np.int64);global_defect=np.full(801,np.nan);storage=np.full(801,np.nan);boundary_names=np.array(['inlet_left','inlet_top','outlet','walls']);adv=np.full((801,4),np.nan);diff=np.full_like(adv,np.nan);total=np.full_like(adv,np.nan)
 audited=[];native_paths=[];previous=None;before=2000;mesh_count=None;frozen_max=0.
 for step,row in enumerate(rows):
  cp=e.touch(row['case'],row['case_sha256']);dp=e.touch(row['data'],row['data_sha256'])
  if step<=5:
   require(sha(cp)==pilot_pins[f'pulse/state-{step:04d}.cas.h5'] and sha(dp)==pilot_pins[f'pulse/state-{step:04d}.dat.h5'],'Original pulse ancestry changed')
  elif recovery is None:require(e.json(run/'audits'/f'state-{step:04d}.json')==row,'Per-step producer JSON differs from state ledger')
  native_clock=read_clock(dp);check_clock(native_clock,step);require(native_clock==row['clock'],'Ledger/native actual clock differs')
  mesh_count=mesh_equal(basecase,cp);require(mesh_count==32,'Actual accepted32 native mesh datasets expected')
  flow=native_fields(cp,dp,reader);changes=fv.frozen_field_changes(baseline,flow);frozen_max=max(frozen_max,max(changes.values()));require(max(changes.values())<=1e-12,'Frozen U/V/p/rho/mu/9040phi/c0/c1')
  c,old,face,datasets=native_scalar(dp,reader,physical);facts=scalar_step_checks(c,old,previous,face,g,step)
  item=dict(step=step,time_s=native_clock['flow-time'],clock=native_clock,case=str(cp),data=str(dp),case_sha256=row['case_sha256'],data_sha256=row['data_sha256'],observed_native_scalar_datasets=datasets,mesh_dataset_count=mesh_count,frozen_changes=changes,**facts)
  c_native[step]=c;c_source[step]=c[source_to_native];times[step]=native_clock['flow-time'];inventory[step]=np.sum(flow['rho']*g['volumes']*c);scalar_min[step]=c.min();scalar_max[step]=c.max();negatives[step]=(c<0).sum();overshoots[step]=(c>1).sum()
  if step:
   if step<=5:
    starts=[i for i,line in enumerate(pilot_lines) if f'pulse/state-{step-1:04d}.dat.h5' in line.replace('\\','/')];ends=[i for i,line in enumerate(pilot_lines) if f'pulse/state-{step:04d}.cas.h5' in line.replace('\\','/')];require(len(starts)==len(ends)==1 and starts[0]<ends[0],'Own pilot transcript phase/write anchors');printed=printed_residuals('\n'.join(pilot_lines[starts[0]:ends[0]]))
   else:printed=production_printed
   own=new_residual(dp,before,printed);before=own['last_iteration'];residual[step-1]=own['terminal_scaled_residual'];iteration[step-1]=before;item['own_new_residual']=own
   boundary={'inlet_left':dict(value=0.,diffusion='dirichlet',evidence='Actual solved native inlet field/accepted Specified Value0'),'inlet_top':dict(value=facts['scheduled_top_value'],diffusion='dirichlet',evidence='Actual solved native inlet field/source-bound scheduled Specified Value'),'outlet':dict(diffusion='zero_gradient',evidence='Root-qualified native boundary interpretation; outward convection/zero diffusion'),'walls':dict(diffusion='zero_gradient',evidence='Root-qualified stationary impermeable wall zero diffusion')}
   _,_,budget=fv.audit_step(g,previous,c,flow['phi'],rho=flow['rho'],gamma=.001,dt_s=.025,boundary=boundary,normalization_time_s=2.,global_floor_kg_m_s=.02)
   require(budget['global_Euler_defect_relative']<=1e-5,'Reconstructed global Euler target');stored=pilot_rows[step]['reconstructed_budget'] if step<=5 else (row['RECONSTRUCTED_budget'] if recovery is None else producer_records.get(step,{}).get('RECONSTRUCTED_budget'))
   if stored is not None:close_budget(budget,stored)
   else:missing_producer_budgets.append(step)
   item['RECONSTRUCTED_budget']=budget;global_defect[step]=budget['global_Euler_defect_relative'];storage[step]=budget['storage_rate_kg_m_s']
   for j,name in enumerate(boundary_names):
    b=budget['reconstructed_boundary'][str(name)];adv[step,j]=b['advective_inward_kg_m_s'];diff[step,j]=b['diffusive_inward_kg_m_s'];total[step,j]=b['total_inward_kg_m_s']
  else:
   with h5py.File(dp,'r') as h:k=np.asarray(h['results/residuals/phase-1/uds-0/iterations'][()]).ravel();require(k.size and k[-1]==2000,'Inherited unsolved2000 row exclusion')
   item['initial_residual_and_Euler_budget_unavailable_not_fabricated']=True
  audited.append(item);native_paths.append(dict(index=step,time_s=float(times[step]),source_kind='native_fluent_hdf',case=str(cp),data=str(dp),case_sha256=row['case_sha256'],data_sha256=row['data_sha256'],observed_scalar_datasets=datasets,budget_json=str(out/'checks.json'),budget_key_path=['states',step,'RECONSTRUCTED_budget'] if step else None))
  previous=c
  if step%100==0:print(json.dumps(dict(readonly_audited_through=step,actual_time_s=float(times[step]))),flush=True)
 missing_producer_postprocessing_files=[];final_saved_readback=None
 if recovery is None:
  for name in ('reports/native-step6-autosave-qualified.json','reports/actual-step80-81-switch-proof.json','checkpoints/top-on-step80-readback.json','checkpoints/top-off-before-step81-readback.json','checkpoints/final-step800-readback.json','loaded-settings.json','static-top-on-before.json','static-top-off-before.json','runtime-discovery.json'):e.json(run/name)
  switch=e.json(run/'reports/actual-step80-81-switch-proof.json');require(switch['actual_boundary_field_proof'] is True and switch['step80']==rows[80] and switch['step81']==rows[81],'Actual80/81 producer/native switch binding')
 else:
  for name in ('reports/native-step6-autosave-qualified.json','reports/actual-step80-81-switch-proof.json','checkpoints/top-on-step80-readback.json','checkpoints/top-off-before-step81-readback.json','checkpoints/final-step800-readback.json','loaded-settings.json','static-top-on-before.json','static-top-off-before.json','runtime-discovery.json'):
   if (run/name).exists():e.json(run/name)
   else:missing_producer_postprocessing_files.append(name)
  final_saved_readback=final_recovery_readback(recovery,e,reader,g,baseline,basecase,previous,old)
 inputs=e.finish();out.mkdir();shutil.copyfile(Path(__file__),out/'executed-auditor.py')
 np.savez_compressed(out/'fields.npz',step=np.arange(801,dtype=np.int64),time_s=times,c=c_source,c_native=c_native,native_cell_to_source_cell=mapping,source_native_cell_index=source_to_native,source_centres_m=g['cell_centres'][source_to_native],source_areas_m2=g['volumes'][source_to_native],solved_time_s=times[1:],terminal_scaled_residual=residual,terminal_native_iteration=iteration,inventory_kg_per_m=inventory,scalar_min=scalar_min,scalar_max=scalar_max,negative_cell_count=negatives,above_one_cell_count=overshoots,global_Euler_defect_relative=global_defect,storage_rate_kg_m_s=storage,boundary_names=boundary_names,reconstructed_boundary_advective_inward_kg_m_s=adv,reconstructed_boundary_diffusive_inward_kg_m_s=diff,reconstructed_boundary_total_inward_kg_m_s=total)
 checks={k:True for k in ('exact801_native_states_and_800_solved_steps','original5_pulse_steps_preserved_and795_new','all_native_clocks_and_state_hashes','all32_native_mesh_datasets_unchanged_each_state','all_frozen_native_U_V_p_rho_mu_9040phi_c0_c1','actual_initial0_and_solved_bounds','all800_native_M1_matches_previous_actual_current','all800_new_native_UDS_residual_targets_and_norms','actual_native_top1_through80_top0_from81','all800_RECONSTRUCTED_global_Euler_targets','all_source_and_native_bytes_unchanged_during_read')}
 record=dict(schema_version=1,kind='Independent READ-ONLY actual Fluent full tracer audit',solver='Fluent261',checks=checks,all_primary_checks_passed=all(checks.values()),source_bytes_unchanged_during_read=True,audited_fields_npz=str(out/'fields.npz'),audited_fields_npz_sha256=sha(out/'fields.npz'),native_run=str(run),native_state_paths=native_paths,states=audited,raw_producer_status=status,raw_producer_summary=summary,accepted_pilot_receipt_sha256=ACCEPT_SHA,production_driver_sha256=DRIVER_SHA,production_targets_sha256=TARGET_SHA,actual_initial_state_count=1,actual_solved_step_count=800,actual_final_native_clock=audited[-1]['clock'],native_mesh_dataset_count=mesh_count,maximum_frozen_field_change=frozen_max,available_ram_before_read_bytes=ram,memory_strategy='Two preallocated801x4400 scalar arrays plus one native state; no801 native object list',native_scalar_flux_classification='Prior actual root-accepted5+5 integrated UDS report total-inward interpretation retained; this exporter does not invent missing production integrated reports',reconstructed_flux_classification='Orthogonal Euler/upwind/Gamma reconstruction, never directly saved native matrix flux',initial_diagnostic_availability='Initial NPZ Euler/storage/boundary-rate NaNs mean unavailable; no initial residual row',root_final_acceptance_pending=True,physical_validation=False,cross_code_agreement_claimed=False,public_record_complete=False,CFD_or_WSL_or_native_mutation_performed=False)
 record['postprocessing_interruption']={'recovery_used':recovery is not None,'root_recovery_receipt':str(args.recovery.resolve()) if recovery is not None else None,'root_recovery_receipt_sha256':args.recovery_sha256,'raw_producer_status_retained_unchanged':True,'raw_producer_summary_absent':summary is None,'producer_budget_steps_available':sorted(producer_records) if recovery is not None else list(range(6,801)),'producer_budget_steps_missing':missing_producer_budgets,'producer_postprocessing_files_missing':missing_producer_postprocessing_files,'final_saved_step800_direct_readback':final_saved_readback,'root_approval':recovery}
 record['producer_missing_audit_budget_steps']=missing_producer_budgets
 record['raw_partial_producer_states']=partial_states if recovery is not None else None
 write_json(out/'checks.json',record);inputs=e.finish();outputs=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in (out/'fields.npz',out/'checks.json',out/'executed-auditor.py')];write_json(out/'source-manifest.json',dict(schema_version=1,kind='Read-only actual native full801-state inputs and exporter outputs',inputs=inputs,outputs=outputs,source_bytes_unchanged_during_read=True,physical_validation=False,root_final_acceptance_pending=True))
 print(json.dumps(dict(out=str(out),all_primary_checks_passed=True,checks_sha256=sha(out/'checks.json'),fields_sha256=sha(out/'fields.npz'),manifest_sha256=sha(out/'source-manifest.json'),source_inputs=len(inputs))),flush=True)
if __name__=='__main__':main()

