"""Read-only five-step OpenFOAM tracer qualification audit. RECONSTRUCTED scalar budgets; no native commands/admission."""
from pathlib import Path
import argparse,csv,ctypes,datetime,hashlib,importlib.util,json,re,shutil,sys
import numpy as np
ROOT=Path('E:/CFD-Learning-Lab');PRIVATE=ROOT/'local/prospective-tracer-v001';RUN=ROOT/'native/runs/2026-11-02-tracer-tjunction-v001';ANCESTOR=ROOT/'native/runs/2026-11-01-tjunction-v001/openfoam'
CASES={'Control':'openfoam-control-v003','Pulse':'openfoam-pulse-pilot-v003'}
PINS={'scalar_fv_audit_v002.py':'75ceec60add2f5c53fb82585d105dcc5cd7c9d8d19e1ea8ac3f2dde14aa455cd','scalar_native_geometry_v003.py':'06f6b4089116d8e7cdab8b4050c30d2939d3876381ebbaa878c7908e6f2d1903','prepare_openfoam_tracer_qualification_v004.py':'278e6948a244e82887d8eb8dc2da90b502a14a7165de2238c3ba1e01b2890888','run-openfoam-tracer-qualification-v004.sh':'a4e4437d8ccd241c724d99148a9a9b3b4c4764785ad01e6b516ea08187082c41','Launch-OpenFOAM-TracerQualification-v004.ps1':'e4c7f96819c3a4bdf25681c68c517e5ec7a107866529f2b171c5540096ec268c'}
TARGET_SHA='173fa3a540d36e4dad195f25db5586d6c15527a6df592dca8076b1999ecb41cd';QUAL_SHA='dccbb17e813f2552298332dc3fd7ff713e30df71a70508316a30831995ccd312';CAPTURE_SHA='79227f2f31fa6dee3fbe03a96ac2429601ee8bf25beb377ea37eddc266b3dcd6';PARSER_SHA='8d87c1aab5a92cfd4c86d499ce83d2317965df4f490ad830d01f37934f674f00'
LINEAR_CAPTURE_SHA='16a0ade334dd6e6eb271180509bc6962c62f83d7106415ce1f32166c43b483c9'
DT=.025;D=.001;RHO=1.;GAMMA=.001;DEPTH=.01

def require(v,m):
 if not v:raise ValueError(m)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
class Evidence:
 def __init__(self):self.rows={}
 def touch(self,p,pin=None):
  p=Path(p).resolve();require(p.is_file(),'Missing actual evidence '+str(p));digest=sha(p);require(pin is None or digest==pin,'Pinned evidence changed '+str(p))
  if p in self.rows:require(self.rows[p]==digest,'Consumed input changed '+str(p))
  else:self.rows[p]=digest
  return p
 def text(self,p,pin=None):return self.touch(p,pin).read_text(encoding='utf-8-sig')
 def readj(self,p,pin=None):return json.loads(self.text(p,pin))
 def module(self,name,p,pin):
  p=self.touch(p,pin);require(name not in sys.modules,'Fresh audit import namespace');spec=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(spec);sys.modules[name]=m;spec.loader.exec_module(m);return m
 def finish(self):
  rows=[dict(path=str(p),bytes=p.stat().st_size,sha256_before=digest,sha256_after=sha(p)) for p,digest in self.rows.items()];require(all(r['sha256_before']==r['sha256_after'] for r in rows),'Native/source input changed during read-only audit')
  return dict(inputs_unchanged_during_analysis=True,source_count=len(rows),sources=rows,native_or_WSL_started=False,case_inputs_modified=False)
def json_new(p,v):
 with Path(p).open('x',encoding='utf-8',newline='\n') as f:json.dump(v,f,indent=2,allow_nan=False);f.write('\n')
def memory_guard():
 class M(ctypes.Structure):_fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ('total','available','pageTotal','pageAvailable','virtualTotal','virtualAvailable','extended')]
 require(sys.platform=='win32','Audit uses the existing Windows E: Python environment');m=M();m.length=ctypes.sizeof(m);require(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))!=0,'Read actual memory');require(m.available>=2.5*1024**3,'Require2.5GiB available RAM for native audit');require(shutil.disk_usage(ROOT).free>=10*1024**3,'Require10GiB E: disk');return int(m.available)
def clock(raw,index,label):
 values={k:re.search(r'\b'+k+r'\s+([-+0-9.eE]+)\s*;',raw) for k in ('value','index','deltaT')};require(all(values.values()),'Native clock value/index/dt missing '+label)
 result={k:float(v[1]) for k,v in values.items()};require(all(np.isfinite(list(result.values()))) and abs(result['value']-index*DT)<=1e-12 and result['index']==index and abs(result['deltaT']-DT)<=1e-12,'Native clock/index/dt differs '+label);return result

def audit(phase,out):
 require(phase in CASES,'Only Control or Pulse');case=RUN/CASES[phase];require(case.resolve()==case,'Native case alias/junction is unsupported')
 out=Path(out).resolve();prefix='openfoam-control-checks-' if phase=='Control' else 'openfoam-pulse-pilot-checks-';require(out.parent==PRIVATE.resolve() and re.fullmatch(re.escape(prefix)+r'v\d{3}',out.name) and not out.exists(),'Fresh versioned private output required')
 available=memory_guard();e=Evidence();e.touch(Path(__file__));e.touch(PRIVATE/'QUALIFICATION-TARGETS-v001.md',TARGET_SHA);e.touch(PRIVATE/'native-flow-qualification-v001.json',QUAL_SHA)
 for name,pin in PINS.items():e.touch(PRIVATE/name,pin)
 # Geometry must instantiate its authenticated original parser namespaces BEFORE the preparer import.
 geo=e.module('tracer_qualification_geometry',PRIVATE/'scalar_native_geometry_v003.py',PINS['scalar_native_geometry_v003.py']);fv=e.module('tracer_qualification_fv',PRIVATE/'scalar_fv_audit_v002.py',PINS['scalar_fv_audit_v002.py'])
 data,geometry_report=geo.qualified_geometry()
 for row in geometry_report['consumed_inputs']:e.touch(row['path'],row['sha256'])
 native=sys.modules['qualified_original_tjunction_openfoam'];e.touch(Path(native.__file__),PARSER_SHA)
 prep=e.module('tracer_qualification_preparer',PRIVATE/'prepare_openfoam_tracer_qualification_v004.py',PINS['prepare_openfoam_tracer_qualification_v004.py'])
 old=prep.old;O=data['openfoam'];g=O['stencil'];qref=O['native_mass_phi'];mapping=O['native_cell_to_source_cell'];checks={};history=[];native_states=[];clock_rows=[]
 require(g['cell_count']==4400 and g['face_count']==9040 and g['depth_m']==DEPTH and np.array_equal(np.sort(mapping),np.arange(4400)),'Qualified native geometry/source mapping differs')
 capture=e.readj(PRIVATE/'openfoam-installed-source-capture-v002/installed-source-capture.json',CAPTURE_SHA)
 auth=e.readj(case/'provenance/root-run-authorization.json');require(auth.get('approved') is True and auth.get('approved_by')=='root' and auth.get('action')=='run_tracer_openfoam_qualification_v004' and auth.get('allowed_phase')==phase and auth.get('actual_steps_authorized')==5 and auth.get('final_time_s')==.125 and auth.get('production_admitted') is False,'Actual exact five-step authorization differs')
 require(auth.get('qualification_targets_sha256')==TARGET_SHA and auth.get('native_flow_qualification_sha256')==QUAL_SHA,'Frozen targets/ancestry binding differs')
 for key,name in zip(prep.KEYS,prep.SOURCES):require(auth.get(key)==PINS[name] and sha(e.touch(case/'provenance'/name))==PINS[name],'Actual reviewed execution source differs '+name)
 linear_capture=e.readj(PRIVATE/'openfoam-linear-solver-sources-v001/capture.json',LINEAR_CAPTURE_SHA);installed=dict(capture['installed_source_sha256'])
 for row in linear_capture['rows']:
  name=row['source'].removeprefix('/opt/openfoam14/');require(name in prep.LINEAR_SOURCE_SHA256 and row['sha256']==prep.LINEAR_SOURCE_SHA256[name],'Actual linear-source capture pin differs');installed[name]=row['sha256'];e.touch(Path(row['copy'].replace('/mnt/e/','E:/')),row['sha256'])
 for name,pin in capture['installed_source_sha256'].items():e.touch(PRIVATE/'openfoam-installed-source-capture-v002/installed-sources'/name,pin)
 require(len(installed)==19 and set(installed)==set(prep.INSTALLED) and auth.get('installed_source_sha256')==installed,'Actual installed-source selection differs')
 for name,pin in prep.HELPERS.items():e.touch(prep.SCRIPTS/name,pin)
 for name,pin in auth['installed_source_sha256'].items():e.touch(case/'provenance/installed-sources'/name,pin)
 for name,pin in e.readj(case/'provenance/immutable-input-sha256.json').items():e.touch(case/name,pin)
 for name,text in prep.dictionaries(phase).items():
  if name!='0/tracer':require(e.text(case/name)==text,'Actual scalar dictionary differs '+name)
 require(e.text(case/'provenance/initial-tracer')==prep.dictionaries(phase)['0/tracer'],'Preserved initial tracer input differs')
 for name in ('physicalProperties','momentumTransport'):require(e.text(case/'constant'/name)==e.text(ANCESTOR/'constant'/name),'Native material/model changed')
 status=e.readj(case/'run-status.json');execution=e.readj(case/'qualification-execution-checks.json');require(status.get('stage')=='finished_five_step_qualification_execution' and status.get('actual_steps')==5 and status.get('native_states')==6 and status.get('phase')==phase,'Actual bounded execution incomplete')
 require(execution.get('completed_five_steps') is True and execution.get('states_verified')==6 and execution.get('phase')==phase,'Actual execution readback incomplete')
 poly=old.native_poly(case);original_poly=old.native_poly(ANCESTOR)
 for name in ('points','faces','owner','neighbour','boundary'):e.touch(poly/name);require(sha(poly/name)==sha(e.touch(original_poly/name)),'Native mesh source copy changed')
 mesh_quality=old.mesh_quality(e.text(case/'mesh-check.log'));patches=old.native_patches(poly)
 # Exact pinned v004 native parser retains []==dimensionless zero; original parsers stay unchanged.
 field=lambda p,n,k,d:np.asarray(prep.field(e.touch(p),n,k,d),dtype=float)
 require(np.max(abs(field(case/'0/Vc',4400,1,[0,3,0,0,0,0,0])-g['volumes']))<=1e-14,'Native Vc differs from oriented geometry')
 require(np.max(abs(field(case/'0/C',4400,3,[0,1,0,0,0,0,0])[:,:2]-g['cell_centres']))<=1e-12,'Native C differs from authenticated native geometry')
 baseU=field(ANCESTOR/'2000/U',4400,3,[0,1,-1,0,0,0,0]);basep=field(ANCESTOR/'2000/p',4400,1,[0,2,-2,0,0,0,0]);cstates=[];qstates=[];qfullstates=[];frozen=[];bcstates=[]
 for i in range(6):
  label=f'{i*DT:.10g}';folder=case/label;raw=e.text(folder/'tracer');require(re.search(r'\bclass\s+volScalarField\s*;',raw) and re.search(r'\bobject\s+tracer\s*;',raw),'Native tracer class/object differs')
  c=field(folder/'tracer',4400,1,[0,0,0,0,0,0,0]);bs=prep.blocks(raw);require(set(bs)==set(patches),'Native scalar boundary coverage differs');boundary={}
  for name in ('inlet_left','inlet_top'):
   expected_type='uniformFixedValue' if phase=='Pulse' and name=='inlet_top' else 'fixedValue';require(re.search(r'\btype\s+'+expected_type+r'\s*;',bs[name]),'Native scalar inlet type differs')
   values=np.asarray(prep.patch_value(bs[name],20));expected=1 if phase=='Control' or name=='inlet_top' else 0;require(np.max(abs(values-expected))<=1e-12,'Actual scalar inlet values differ')
   boundary[name]=dict(diffusion='dirichlet',value=values,evidence=str(folder/'tracer')+' actual saved '+expected_type+' prescribed value; initial0 is unsolved input')
  for name in ('outlet','walls'):
   require(re.search(r'\btype\s+zeroGradient\s*;',bs[name]),'Native scalar zeroGradient differs');boundary[name]=dict(diffusion='zero_gradient',evidence=str(folder/'tracer')+' native zeroGradient type plus authenticated installed zeroGradient source; owner face c and diffusive zero are reconstructed BC consequences')
  require(re.search(r'\btype\s+empty\s*;',bs[old.EMPTY]),'Native scalar empty planes differ')
  e.touch(folder/'phi');qfull=np.asarray(prep.all_phi(folder/'phi',poly));q=qfull[g['native_face_ids']]*RHO;u=field(folder/'U',4400,3,[0,1,-1,0,0,0,0]);p=field(folder/'p',4400,1,[0,2,-2,0,0,0,0])
  frozen.append(dict(index=i,U_max_abs_change=float(abs(u-baseU).max()),p_max_abs_change=float(abs(p-basep).max()),phi_max_abs_change=float(abs(q-qref).max())))
  require(np.all(q[(g['patch']=='outlet')&~g['inside']]>=0) and np.all(q[np.isin(g['patch'],['inlet_left','inlet_top'])]<=0),'Unsupported native boundary flow reversal')
  cstates.append(c);qstates.append(q);qfullstates.append(qfull);bcstates.append(boundary);native_states.append(dict(index=i,time_s=i*DT,folder=label,global_tracer_path=str(folder/'tracer'),global_tracer_sha256=e.rows[(folder/'tracer').resolve()],rank_tracer_paths=[],rank_clock_paths=[]))
  if (folder/'uniform/time').is_file():clock_rows.append(dict(scope='reconstructed',expected_index=i,path=str(folder/'uniform/time'),**clock(e.text(folder/'uniform/time'),i,label)))
 rank_cells=[]
 for rank in (0,1):
  rp=case/('processor'+str(rank));rpoly=old.native_poly(rp)
  for name in ('points','faces','owner','neighbour','boundary','cellProcAddressing','faceProcAddressing','pointProcAddressing'):e.touch(rpoly/name)
  n,body=old.counted(rpoly/'cellProcAddressing');cellmap=np.asarray(list(map(int,body.split())),dtype=int);require(len(cellmap)==n and len(np.unique(cellmap))==n and np.all((cellmap>=0)&(cellmap<4400)),'Native rank cell addressing differs')
  rank_cells.extend(cellmap.tolist());nf,body=old.counted(rpoly/'faceProcAddressing');facemap=np.asarray(list(map(int,body.split())),dtype=int);faceids=abs(facemap)-1
  require(len(facemap)==nf and len(np.unique(faceids))==nf and np.all(facemap!=0) and np.all((faceids>=0)&(faceids<len(qfullstates[0]))),'Native signed rank face addressing differs')
  expected_phi=qfullstates[0][faceids]*np.sign(facemap)
  for i in range(6):
   label=f'{i*DT:.10g}';folder=rp/label
   require(np.max(abs(field(folder/'U',n,3,[0,1,-1,0,0,0,0])-baseU[cellmap]))<=1e-12 and np.max(abs(field(folder/'p',n,1,[0,2,-2,0,0,0,0])-basep[cellmap]))<=1e-12,'Native rank frozen U/p differs')
   e.touch(folder/'phi');require(np.max(abs(np.asarray(prep.all_phi(folder/'phi',rpoly))-expected_phi))<=1e-12,'Native signed rank all-face phi differs')
   cr=field(folder/'tracer',n,1,[0,0,0,0,0,0,0]);require(np.max(abs(cr-cstates[i][cellmap]))<=1e-12,'Native rank/global scalar differs');native_states[i]['rank_tracer_paths'].append(str(folder/'tracer'))
   if i:
    cp=folder/'uniform/time';clock_rows.append(dict(scope='processor'+str(rank),expected_index=i,path=str(cp),**clock(e.text(cp),i,label)));native_states[i]['rank_clock_paths'].append(str(cp))
 require(sorted(rank_cells)==list(range(4400)),'Native two-rank cells do not biject to reconstructed cells')
 require(len([r for r in clock_rows if r['scope'].startswith('processor')])==10,'Require two ranks times/index/dt for five actual steps')
 log=e.text(case/'solver.log');require(re.search(r'^End\s*$',log,re.M) and not re.search(r'FOAM FATAL|Floating point exception|Segmentation fault',log),'Actual scalar log incomplete/fatal')
 times=[float(t) for t in re.findall(r'^Time = ('+old.FLOAT+r')s?\s*$',log,re.M)];require(len(times)==5 and np.max(abs(np.asarray(times)-np.arange(1,6)*DT))<=1e-12,'Actual five scalar time endpoints differ')
 sections=re.split(r'^Time = ('+old.FLOAT+r')s?\s*$',log,flags=re.M);residuals=[]
 for index in range(1,6):
  body=sections[2*index];matches=re.findall(r'Solving for tracer, Initial residual = ('+old.FLOAT+r'), Final residual = ('+old.FLOAT+r'), No Iterations (\d+)',body);require(len(matches)==1,'One scalar matrix solve per physical step required');a,b,n=matches[0];require(np.isfinite([float(a),float(b)]).all(),'Nonfinite native scalar residual');residuals.append(dict(index=index,time_s=index*DT,initial_residual=float(a),final_residual=float(b),iterations=int(n)))
 require(not re.search(r'Solving for (?:Ux|Uy|Uz|p),',log),'Frozen momentum/pressure was recomputed')
 checks['two_rank_native_scalar_and_frozen_fields_addressing_and_clocks_verified']=True;checks['native_mesh_quality_and_scalar_schema_verified']=True;checks['all_saved_U_p_and_every_active_native_phi_frozen']=all(max(r['U_max_abs_change'],r['p_max_abs_change'],r['phi_max_abs_change'])<=1e-12 for r in frozen)
 checks['initial_scalar_matches_declared_value']=bool(np.max(abs(cstates[0]-(1 if phase=='Control' else 0)))<=1e-12);checks['own_native_final_matrix_residuals_le_1e_8']=all(r['final_residual']<=1e-8 for r in residuals)
 checks['all_solved_control_constancy_or_pulse_bounds']=bool(all(np.max(abs(c-1))<=1e-6 for c in cstates[1:])) if phase=='Control' else bool(all(c.min()>=-1e-8 and c.max()<=1+1e-8 for c in cstates[1:]));checks['pulse_genuinely_changes_or_control_exempt']=phase=='Control' or any(not np.array_equal(c,cstates[0]) for c in cstates[1:])
 face_terms=[];local_defects=[];I0=float(np.dot(g['volumes'],cstates[0])*RHO/DEPTH);co=[]
 for i in range(1,6):
  terms,defect,summary=fv.audit_step(g,cstates[i-1],cstates[i],qstates[i],rho=RHO,gamma=GAMMA,dt_s=DT,boundary=bcstates[i],normalization_time_s=2.,global_floor_kg_m_s=.02)
  history.append(dict(index=i,time_s=i*DT,actual_native_scalar_binding=True,native_current_scalar_path=native_states[i]['global_tracer_path'],native_previous_scalar_path=native_states[i-1]['global_tracer_path'],fv_helper_return=summary,global_reconstructed_budget_pass=summary['global_Euler_defect_relative']<=1e-5));face_terms.append(terms);local_defects.append(defect)
  q=qstates[i];inside=g['inside'];sums=np.bincount(g['owner'],weights=abs(q),minlength=4400)+np.bincount(g['neighbour'][inside],weights=abs(q[inside]),minlength=4400);co.append(float((DT*.5*sums/(RHO*g['volumes'])).max()))
 checks['all_solved_reconstructed_global_Euler_defects_le_1e_5']=all(r['global_reconstructed_budget_pass'] for r in history)
 if phase=='Control':checks['control_relative_inventory_change_le_1e_6']=all(abs(r['fv_helper_return']['inventory_kg_per_m']/I0-1)<=1e-6 for r in history)
 manifest=e.finish();checks['all_consumed_native_and_source_bytes_unchanged']=True
 c_native=np.asarray(cstates);c_source=np.empty_like(c_native);c_source[:,mapping]=c_native
 cumulative={name:dict(advective_inward_kg_per_m=0.,diffusive_inward_kg_per_m=0.) for name in ('inlet_left','inlet_top','outlet','walls')}
 for row in history:
  for name,b in row['fv_helper_return']['reconstructed_boundary'].items():
   cumulative[name]['advective_inward_kg_per_m']+=DT*b['advective_inward_kg_m_s'];cumulative[name]['diffusive_inward_kg_per_m']+=DT*b['diffusive_inward_kg_m_s']
 result=dict(phase=phase,native_case=str(case),native_scalar_field_name='tracer',mathematical_scalar_symbol='c',actual_steps=5,native_state_count=6,dt_s=DT,final_time_s=.125,checks=checks,failed_checks=[k for k,v in checks.items() if v is not True],all_declared_numerical_checks_passed=all(v is True for v in checks.values()),root_acceptance_performed=False,production_admitted=False,physical_validation=False,analysis_source_sha256=sha(Path(__file__)),qualification_targets_sha256=TARGET_SHA,native_state_paths=native_states,native_clocks=clock_rows,initial_clock_classification='Native initial tracer fields plus explicit startFrom/startTime0/dt input; initial0 is unsolved and its missing clock file is not invented',raw_native_run_status=status,raw_execution_checker=execution,mesh_quality=mesh_quality,frozen_native_field_changes=frozen,actual_native_residuals=residuals,residual_norm='Actual PBiCGStab recursive residual with native uniform-reference lduMatrix normFactor; original normalization unchanged and no equivalence to Fluent asserted',reconstructed_scalar_steps=history,reconstructed_cumulative_boundary_inventory=cumulative,native_cell_Co_max_each_step=co,Co_is_descriptive_not_accuracy_gate=True,local_reconstructed_defect_has_no_pass_threshold=True,scalar_evidence=dict(current_fields='Actual native volScalarField tracer at all six saved times with rank/global agreement',density='Declared constant rho1 mass conversion; incompressible OpenFOAM has no persisted native rho field',diffusivity='Exact source-bound runtime constant D0.001; Gamma=rhoD0.001; no saved diffusivity-array claim',old_time='Consecutive actual saved scalar states; no persisted internal old-time array claim',direct_native_total_scalar_flux='Unavailable from ordinary scalarTransport; all scalar transfers/defects here are RECONSTRUCTED',FV_helper_binding='Raw pure-helper metadata is retained; this caller authenticates actual scalar paths/hashes separately',initial_flux='Unsolved initial scalar state is not a solved conservation check'),available_RAM_at_audit_bytes=available,geometry_identity_report=geometry_report)
 out.mkdir();json_new(out/'checks.json',result);json_new(out/'source-manifest.json',manifest)
 np.savez_compressed(out/'fields.npz',step=np.arange(6),time_s=np.arange(6)*DT,c=c_source,c_native=c_native,native_cell_to_source_cell=mapping,native_mass_phi=np.asarray(qstates),reconstructed_advective_outward_kg_s=np.asarray([t['reconstructed_advective_outward_kg_s'] for t in face_terms]),reconstructed_diffusive_outward_kg_s=np.asarray([t['reconstructed_diffusive_outward_kg_s'] for t in face_terms]),reconstructed_cell_defect_kg_s=np.asarray(local_defects))
 rows=[]
 for row,res in zip(history,residuals):
  s=row['fv_helper_return'];r=dict(index=row['index'],time_s=row['time_s'],scalar_min=s['scalar_min'],scalar_max=s['scalar_max'],inventory_kg_per_m=s['inventory_kg_per_m'],storage_rate_kg_m_s=s['storage_rate_kg_m_s'],global_Euler_defect_relative=s['global_Euler_defect_relative'],local_reconstructed_defect_max=s['reconstructed_local_dimensionless_defect_max'],local_reconstructed_defect_RMS=s['reconstructed_local_dimensionless_defect_volume_RMS'],final_native_matrix_residual=res['final_residual'])
  for name,b in s['reconstructed_boundary'].items():
   for key,v in b.items():r[name+'_'+key]=v
  rows.append(r)
 with (out/'step-history.csv').open('x',encoding='ascii',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 for name in ('scalar_fv_audit_v002.py','scalar_native_geometry_v003.py','prepare_openfoam_tracer_qualification_v004.py'):shutil.copy2(PRIVATE/name,out/('executed-'+name))
 shutil.copy2(Path(__file__),out/'executed-auditor.py');json_new(out/'output-manifest.json',dict(outputs=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file()]))
 return dict(out=str(out),checks_passed=sum(v is True for v in checks.values()),check_count=len(checks),all_declared_numerical_checks_passed=result['all_declared_numerical_checks_passed'],maximum_reconstructed_global_defect=max(r['fv_helper_return']['global_Euler_defect_relative'] for r in history),actual_steps=5,root_acceptance_performed=False,production_admitted=False)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--phase',choices=tuple(CASES),required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(audit(a.phase,a.out),indent=2,allow_nan=False))