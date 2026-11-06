"""SOURCE ONLY until root review and closed actual production: read-only 801-state OpenFOAM tracer audit/export.
No CFD/WSL/RPC/case writes. All scalar transfers and defects are RECONSTRUCTED, never native matrix flux.
"""
from pathlib import Path,PurePosixPath
import argparse,csv,hashlib,importlib.util,json,re,shutil,sys
import numpy as np
ROOT=Path('E:/CFD-Learning-Lab');PRIVATE=ROOT/'local/prospective-tracer-v001';RUN=ROOT/'native/runs/2026-11-02-tracer-tjunction-v001'
ORIGINAL_CASE=RUN/'openfoam-production-v001';CASE=RUN/'openfoam-production-reconstructed-v002';PILOT=RUN/'openfoam-pulse-pilot-v003';ANCESTOR=ROOT/'native/runs/2026-11-01-tjunction-v001/openfoam'
QA_NAME='audit_openfoam_tracer_qualification_v003.py';QA_SHA='98bbb0764bd1ce9b8da5c325a296044689fddbe2a3a82e1250abab177fca173f'
PROD_PINS={'prepare_openfoam_tracer_production_v001.py':'eab29753a8a63047c1bdf21953630d418cde5f68a46b81846937cd8a52b980ac','run-openfoam-tracer-production-v001.sh':'bb9ba176ef243fb8fa3baaf707249456dcb6ae558927c9e6a8df46d90125e8aa','Launch-OpenFOAM-TracerProduction-v001.ps1':'64ea9e1fff0e4a63331fc019b9bedad8cce5991775eb8c2e6e7f0840f7983ec7'}
TARGET_SHA='07f5b0745fc71758b9cf5678171ce17fd558f00a01cd773560911920e91318ba';DT=.025;RHO=1.;GAMMA=.001;DEPTH=.01

def require(v,message):
 if not v:raise ValueError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1048576),b''):h.update(block)
 return h.hexdigest()
require(sha(PRIVATE/QA_NAME)==QA_SHA,'Reviewed qualifier auditor changed before import')
spec=importlib.util.spec_from_file_location('production_audit_helpers',PRIVATE/QA_NAME);qa=importlib.util.module_from_spec(spec);spec.loader.exec_module(qa)
Evidence=qa.Evidence;clock=qa.clock;json_new=qa.json_new

def log_rows(raw,steps,float_pattern):
 """Actual native text; one selected scalar matrix per exact prescribed endpoint."""
 steps=tuple(steps);require(re.search(r'^End\s*$',raw,re.M) and not re.search(r'FOAM FATAL|Floating point exception|Segmentation fault',raw),'Incomplete/fatal scalar log')
 require(not re.search(r'Solving for (?:Ux|Uy|Uz|p),',raw),'Frozen momentum/pressure was recomputed')
 split=re.split(r'^Time = ('+float_pattern+r')s?\s*$',raw,flags=re.M)
 times=[float(v) for v in split[1::2]];require(len(times)==len(steps) and all(abs(t-i*DT)<=1e-12 for t,i in zip(times,steps)),'Missing/replayed/extra native scalar endpoints')
 rows=[]
 for k,i in enumerate(steps):
  body=split[2*k+2];all_solves=re.findall(r'Solving for tracer,',body)
  found=re.findall(r'DILUPBiCGStab:\s+Solving for tracer, Initial residual = ('+float_pattern+r'), Final residual = ('+float_pattern+r'), No Iterations (\d+)',body)
  require(len(all_solves)==len(found)==1,'Exactly one qualified DILUPBiCGStab scalar matrix per endpoint required')
  a,b,n=found[0];a,b,n=float(a),float(b),int(n);require(np.isfinite([a,b]).all() and a>=0 and 0<=b<=1e-8 and 0<=n<=1000,'Native scalar matrix residual/iteration gate missed')
  rows.append(dict(index=i,time_s=i*DT,initial_residual=a,final_residual=b,iterations=n))
 return rows

def validate_recovery_maps(original_before,original_after,replica_before,replica_after,recovered_final):
 """Whole-file byte lineage; only regenerated global states6..800/new recovery metadata may differ."""
 for rows in (original_before,original_after,replica_before,replica_after,recovered_final):
  require(isinstance(rows,dict) and rows,'Complete nonempty recovery tree SHA map required');folded=set()
  for name,pin in rows.items():
   require(isinstance(name,str),'Native relative path string required');path=PurePosixPath(name)
   require(name==path.as_posix() and not path.is_absolute() and name and all(v not in ('..','.') for v in path.parts),'Unsafe recovery relative path')
   require(name.casefold() not in folded,'Case-folding recovery path collision');folded.add(name.casefold())
   require(isinstance(pin,str) and re.fullmatch('[0-9a-f]{64}',pin),'Native tree SHA required')
 require(original_before==original_after==replica_before,'Original tree changed or initial full replica copy differs')
 require(replica_after==recovered_final,'Exported recovered case differs from complete reconstructed replica')
 require(set(original_before)<=set(recovered_final),'Recovery deleted an original file')
 require(not any(name.startswith('recovery-provenance/') for name in original_before),'Original already contains reserved recovery metadata')
 labels={f'{i*DT:.10g}' for i in range(6,801)};changed=[]
 for name,pin in recovered_final.items():
  if original_before.get(name)!=pin:
   require(name.split('/',1)[0] in labels and '/' in name or name.startswith('recovery-provenance/'),'Changed nonderived original/rank/pilot/producer file '+name);changed.append(name)
 return sorted(changed)

def recovery_evidence(e,prod,authpath,permitpath):
 permitpath=Path(permitpath).resolve();require(PRIVATE.resolve() in permitpath.parents,'Exact private root recovery-audit permit required');permit=e.readj(permitpath)
 require(permit.get('approved') is True and permit.get('approved_by')=='root' and permit.get('action')=='audit_tracer_openfoam_reconstruction_recovery_v001','Exact root recovery-audit authorization required')
 require(prod.project_path(permit['original_case'])==ORIGINAL_CASE and prod.project_path(permit['recovered_case'])==CASE,'Exact original/recovered cases required')
 require(permit.get('original_production_authorization_sha256')==e.rows[authpath]=='377fd5d23cc7c557eaa258db4d97073be86a4c61903ef0002f813a7595194f64','Original795-step production authorization changed')
 require(permit.get('native_commands_closed') is True and permit.get('project_distro_closed') is True and permit.get('full_results_accepted') is False,'Root must observe closed native postprocessing before independent audit')
 def ref(binding):
  require(isinstance(binding,dict) and isinstance(binding.get('sha256'),str) and re.fullmatch('[0-9a-f]{64}',binding['sha256']),'Exact recovery path/SHA binding required')
  return e.touch(prod.project_path(binding['path']),binding['sha256'])
 receiptpath=ref(permit['recovery_execution_receipt']);receipt=e.readj(receiptpath);closurepath=ref(permit['closure_proof']);closure=e.readj(closurepath)
 require(receipt.get('schema')=='openfoam_reconstruction_recovery_v001' and prod.project_path(receipt['original_case'])==ORIGINAL_CASE and prod.project_path(receipt['recovered_case'])==CASE,'Recovery execution schema/cases differ')
 original_auth=ref(receipt['original_production_authorization']);require(original_auth==ORIGINAL_CASE/'production-provenance/root-production-authorization.json' and e.rows[original_auth]==e.rows[authpath] and receipt.get('original_wrapper_windows_pid')==117552,'Exact actual original producer lineage required')
 worker=receipt['recovery_source_sha256'];require(worker and worker==permit['recovery_source_sha256'],'Reviewed recovery worker source bindings differ')
 for raw,pin in worker.items():ref(dict(path=raw,sha256=pin))
 keys=('original_before','original_after','replica_before','replica_after','recovered_final');require(set(receipt['manifests'])==set(keys),'All5 source/copy/final tree manifests required');maps={key:e.readj(ref(receipt['manifests'][key])) for key in keys};changed=validate_recovery_maps(*(maps[key] for key in keys))
 require(prod.tree_sha(ORIGINAL_CASE)==maps['original_before'] and prod.tree_sha(CASE)==maps['recovered_final'],'Current original/final tree bytes differ from native recovery manifests')
 for name,pin in maps['original_before'].items():e.touch(ORIGINAL_CASE/name,pin)
 for name,pin in maps['recovered_final'].items():e.touch(CASE/name,pin)
 for key,name in [('original_status','run-status.json'),('original_sync','production-sync-completed.json'),('original_solver_log','production-solver.log')]:require(ref(receipt[key])==ORIGINAL_CASE/name,'Original producer evidence path differs')
 status=e.readj(ORIGINAL_CASE/'run-status.json');sync=e.readj(ORIGINAL_CASE/'production-sync-completed.json')
 require(status.get('stage')=='failed' and status.get('failed_stage')=='reconstructing_exactly795_new_native_states' and status.get('full_signed_scalar_budget_analysis_pending') is True,'Preserve actual original postprocessing failure/pending flags')
 require(sync.get('native_commands_returned') is True and sync.get('explicit_sync_exit_code')==0 and isinstance(sync.get('wrapper_exit_code'),int) and sync['wrapper_exit_code']!=0,'Preserve actual failed wrapper plus successful sync')
 post=receipt['postprocess'];require(post.get('command')==['reconstructPar','-time','0.15:20'] and post.get('exit_code')==post.get('sync_exit_code')==0 and post.get('native_commands_returned') is True,'Only complete native reconstructPar postprocessing admitted')
 logpath=ref(post['log']);require(logpath==CASE/'recovery-provenance/reconstruction.log','Fresh recovered reconstruction log path required');log=e.text(logpath)
 require(re.search(r'^End\s*$',log,re.M) and not re.search(r'FOAM FATAL|Floating point exception|Segmentation fault',log),'Recovery native reconstruction incomplete/fatal')
 require(all(receipt.get(key) is True for key in ('copied_rank_files_unchanged','pilot183_preserved','original_all_bytes_unchanged','no_CFD_replay')) and receipt.get('full_results_accepted') is False and receipt.get('full_signed_scalar_budget_analysis_pending') is True,'Byte preservation/no-CFD-replay evidence required; independent numerical audit still pending')
 execution_path=ORIGINAL_CASE/'production-execution-checks.json';execution=e.readj(execution_path) if execution_path.is_file() else None
 summary=dict(original_postprocessing_failed=True,original_failed_stage=status['failed_stage'],original_budget_pending=status['full_signed_scalar_budget_analysis_pending'],native_reconstruction_completed=True,reconstruction_exit_code=post['exit_code'],no_CFD_replay=True,original_all_bytes_preserved=True,raw_rank_fields_byte_identical=True,pilot183_preserved=True,full_results_accepted=receipt['full_results_accepted'],raw_recovery_budget_pending=receipt['full_signed_scalar_budget_analysis_pending'])
 evidence=dict(root_recovery_audit_permit=permit,root_recovery_audit_permit_sha256=e.rows[permitpath],worker_execution_receipt=receipt,worker_execution_receipt_sha256=e.rows[receiptpath],root_native_closure_proof=closure,root_native_closure_proof_sha256=e.rows[closurepath],actual_changed_or_added_relative_paths=changed,all_original_files=len(maps['original_before']),all_final_files=len(maps['recovered_final']))
 return status,execution,sync,summary,evidence

def audit(authpath,out,recovery_authorization):
 require(CASE.resolve()==CASE and ORIGINAL_CASE.resolve()==ORIGINAL_CASE and PILOT.resolve()==PILOT,'Native case aliases are unsupported')
 out=Path(out).resolve();require(out.parent==PRIVATE.resolve() and re.fullmatch(r'openfoam-production-checks-v\d{3}',out.name) and not out.exists(),'Fresh private versioned output required')
 available=qa.memory_guard();e=Evidence();e.touch(Path(__file__));e.touch(PRIVATE/QA_NAME,QA_SHA);e.touch(PRIVATE/'PRODUCTION-TARGETS-v001.md',TARGET_SHA);e.touch(PRIVATE/'QUALIFICATION-TARGETS-v001.md',qa.TARGET_SHA)
 e.touch(PRIVATE/'native-flow-qualification-v001.json',qa.QUAL_SHA)
 for name,pin in {**qa.PINS,**PROD_PINS}.items():e.touch(PRIVATE/name,pin)
 # Instantiate the original parser namespaces through the pinned geometry adapter first.
 geo=e.module('openfoam_production_geometry',PRIVATE/'scalar_native_geometry_v003.py',qa.PINS['scalar_native_geometry_v003.py']);fv=e.module('openfoam_production_fv',PRIVATE/'scalar_fv_audit_v002.py',qa.PINS['scalar_fv_audit_v002.py'])
 data,geometry_report=geo.qualified_geometry()
 for row in geometry_report['consumed_inputs']:e.touch(row['path'],row['sha256'])
 prod=e.module('openfoam_production_preparer',PRIVATE/'prepare_openfoam_tracer_production_v001.py',PROD_PINS['prepare_openfoam_tracer_production_v001.py']);q=prod.q;old=q.old
 require(prod.CASE==ORIGINAL_CASE and prod.PILOT==PILOT and prod.TARGET_SHA==TARGET_SHA,'Pinned production path/targets differ')
 authpath=e.touch(Path(authpath).resolve());auth=e.readj(authpath);prod.authorization(authpath)
 e.touch(CASE/'production-provenance/root-production-authorization.json',e.rows[authpath]);e.touch(CASE/'production-provenance/PRODUCTION-TARGETS-v001.md',TARGET_SHA)
 for name,pin in PROD_PINS.items():e.touch(CASE/'production-provenance'/name,pin)
 for name,pin in q.HELPERS.items():e.touch(q.SCRIPTS/name,pin);e.touch(CASE/'provenance'/name,pin)
 for name,pin in prod.QUAL_PINS.items():e.touch(CASE/'provenance'/name,pin)
 for raw,pin in auth['selected_ancestry_sha256'].items():e.touch(prod.project_path(raw),pin)
 for row in auth['acceptance_evidence'].values():e.touch(prod.project_path(row['path']),row['sha256'])
 e.touch(prod.project_path(auth['pilot_run_authorization']['path']),auth['pilot_run_authorization']['sha256'])
 acceptances={name:e.readj(PRIVATE/name,pin) for name,pin in prod.ACCEPTANCE_PINS.items()};accepted=acceptances['root-OpenFOAM-qualification-acceptance-v001.json']
 pins={r['relative_path']:r['sha256'] for r in accepted['pulse_ancestry_file_manifest']};require(len(pins)==183 and pins==auth['pilot_file_sha256'],'Exact accepted183-file pulse lineage differs');prod.validate_relative_paths(pins)
 require(e.readj(CASE/'production-provenance/accepted-pilot-all-file-sha256.json')==pins,'Copied accepted pulse manifest differs')
 require(prod.tree_sha(PILOT)==pins,'Original accepted pilot tree changed')
 for name,pin in pins.items():
  e.touch(PILOT/name,pin);copy=CASE/'production-provenance/pilot-overridden'/name if name in prod.OVERRIDDEN else CASE/name;e.touch(copy,pin)
 for name,pin in e.readj(CASE/'provenance/immutable-input-sha256.json').items():e.touch(CASE/'production-provenance/pilot-overridden'/name if name in prod.OVERRIDDEN else CASE/name,pin)
 require(e.text(CASE/'system/controlDict')==prod.production_control(),'Exact795-step restart controls differ')
 for name,text in q.dictionaries('Pulse').items():
  if name not in ('system/controlDict','0/tracer'):require(e.text(CASE/name)==text,'Qualified scalar dictionary changed '+name)
 capture=e.readj(PRIVATE/'openfoam-installed-source-capture-v002/installed-source-capture.json',qa.CAPTURE_SHA);linear=e.readj(PRIVATE/'openfoam-linear-solver-sources-v001/capture.json',qa.LINEAR_CAPTURE_SHA);installed=dict(capture['installed_source_sha256'])
 for row in linear['rows']:
  name=row['source'].removeprefix('/opt/openfoam14/');require(q.LINEAR_SOURCE_SHA256.get(name)==row['sha256'],'Linear source capture differs');installed[name]=row['sha256'];e.touch(Path(row['copy'].replace('/mnt/e/','E:/')),row['sha256'])
 require(len(installed)==19 and set(installed)==set(q.INSTALLED) and installed==auth['installed_source_sha256'],'All19 installed-source bindings differ')
 for name,pin in installed.items():
  e.touch(CASE/'provenance/installed-sources'/name,pin)
  if name in capture['installed_source_sha256']:e.touch(PRIVATE/'openfoam-installed-source-capture-v002/installed-sources'/name,pin)
 qualification=e.readj(PRIVATE/'native-flow-qualification-v001.json',qa.QUAL_SHA)
 for row in qualification['source_manifest']:e.touch(Path(row['path']),row['sha256_after'])
 q.verify_sources()
 status,execution,sync,recovery_summary,recovery_details=recovery_evidence(e,prod,authpath,recovery_authorization)
 O=data['openfoam'];g=O['stencil'];mapping=np.asarray(O['native_cell_to_source_cell'],dtype=int);require(g['cell_count']==4400 and g['face_count']==9040 and g['depth_m']==DEPTH and np.array_equal(np.sort(mapping),np.arange(4400)),'Actual source-order geometry differs')
 poly=old.native_poly(CASE);original_poly=old.native_poly(ANCESTOR);patches=old.native_patches(poly)
 for name in ('points','faces','owner','neighbour','boundary'):require(sha(e.touch(poly/name))==sha(e.touch(original_poly/name)),'Actual native mesh changed')
 mesh_quality=old.mesh_quality(e.text(CASE/'production-mesh-check.log'))
 field=lambda path,n,k,d:np.asarray(q.field(e.touch(path),n,k,d),dtype=float)
 require(np.max(abs(field(CASE/'0/Vc',4400,1,[0,3,0,0,0,0,0])-g['volumes']))<=1e-14,'Actual native Vc differs from geometry')
 require(np.max(abs(field(CASE/'0/C',4400,3,[0,1,0,0,0,0,0])[:,:2]-g['cell_centres']))<=1e-12,'Actual native C/source ordering differs')
 baseU=field(ANCESTOR/'2000/U',4400,3,[0,1,-1,0,0,0,0]);basep=field(ANCESTOR/'2000/p',4400,1,[0,2,-2,0,0,0,0]);qfull=np.asarray(q.all_phi(e.touch(ANCESTOR/'2000/phi'),original_poly));phi=qfull[g['native_face_ids']]*RHO
 require(np.max(abs(phi-O['native_mass_phi']))<=1e-12 and np.all(phi[(g['patch']=='outlet')&~g['inside']]>=0) and np.all(phi[np.isin(g['patch'],['inlet_left','inlet_top'])]<=0),'Changed/unsupported native frozen flux')
 inside=g['inside'];sums=np.bincount(g['owner'],weights=abs(phi),minlength=4400)+np.bincount(g['neighbour'][inside],weights=abs(phi[inside]),minlength=4400);Co_max=float((DT*.5*sums/(RHO*g['volumes'])).max())
 old_log=e.text(CASE/'solver.log');new_log=e.text(CASE/'production-solver.log');residuals=log_rows(old_log,tuple(range(1,6)),old.FLOAT)+log_rows(new_log,tuple(range(6,801)),old.FLOAT);prod.native_times(new_log)
 require([r['index'] for r in residuals]==list(range(1,801)),'Actual log lineage does not cover800 solved states')
 c_native=np.empty((801,4400),dtype=float);states=[];budgets=[];clocks=[];frozen=[];edge=[];checks={};cumulative={name:dict(advective_inward_kg_per_m=0.,diffusive_inward_kg_per_m=0.) for name in ('inlet_left','inlet_top','outlet','walls')}
 expected_dirs={prod.label(i) for i in range(801)}|{'2000'};require({p.name for p in CASE.iterdir() if p.is_dir() and re.fullmatch(old.FLOAT,p.name)}==expected_dirs,'Missing/extra reconstructed state directories')
 for i in range(801):
  label=prod.label(i);folder=CASE/label;scalar_path=folder/'tracer';raw=e.text(scalar_path);require(re.search(r'\bclass\s+volScalarField\s*;',raw) and re.search(r'\bobject\s+tracer\s*;',raw),'Native tracer class/object differs')
  c=field(scalar_path,4400,1,[0]*7);c_native[i]=c;bs=q.blocks(raw);require(set(bs)==set(patches),'Native scalar boundary coverage differs');boundary={};top=prod.expected_top(i)
  for name,expected,kind in [('inlet_left',0,'fixedValue'),('inlet_top',top,'uniformFixedValue')]:
   require(re.search(r'\btype\s+'+kind+r'\s*;',bs[name]),'Native scalar inlet type differs');values=np.asarray(q.patch_value(bs[name],20));require(np.max(abs(values-expected))<=1e-12,'Actual saved inlet value differs at'+label)
   boundary[name]=dict(diffusion='dirichlet',value=values,evidence=str(scalar_path)+' actual native saved '+kind+' value; source-bound inlet diffusion retained')
  for name in ('outlet','walls'):
   require(re.search(r'\btype\s+zeroGradient\s*;',bs[name]),'Native scalar zeroGradient differs');boundary[name]=dict(diffusion='zero_gradient',evidence=str(scalar_path)+' native zeroGradient and authenticated installed source; reconstructed owner c/zero diffusive gradient')
  require(re.search(r'\btype\s+empty\s*;',bs[old.EMPTY]),'Native empty planes differ')
  u=field(folder/'U',4400,3,[0,1,-1,0,0,0,0]);p=field(folder/'p',4400,1,[0,2,-2,0,0,0,0]);nowphi=np.asarray(q.all_phi(e.touch(folder/'phi'),poly));frozen.append(dict(index=i,U_max_abs_change=float(abs(u-baseU).max()),p_max_abs_change=float(abs(p-basep).max()),all_native_phi_max_abs_change=float(abs(nowphi-qfull).max())))
  current_clock=None
  if i:
   cp=folder/'uniform/time';current_clock=dict(scope='reconstructed',expected_index=i,path=str(cp),**clock(e.text(cp),i,label));clocks.append(current_clock)
  state=dict(index=i,time_s=i*DT,native_folder=label,global_tracer_path=str(scalar_path),global_tracer_sha256=e.rows[scalar_path.resolve()],rank_tracer_paths=[],rank_tracer_sha256=[],rank_clock_paths=[],rank_clock_sha256=[],global_clock=current_clock,scalar_min=float(c.min()),scalar_max=float(c.max()),negative_cells=int((c<0).sum()),above_one_cells=int((c>1).sum()),inventory_kg_per_m=float(np.dot(g['volumes'],c)*RHO/DEPTH),actual_top_c_min=float(min(boundary['inlet_top']['value'])),actual_top_c_max=float(max(boundary['inlet_top']['value'])))
  states.append(state)
  if i in (80,81):edge.append(dict(index=i,time_s=i*DT,tracer_path=str(scalar_path),tracer_sha256=state['global_tracer_sha256'],top_c_min=state['actual_top_c_min'],top_c_max=state['actual_top_c_max'],native_clock=current_clock))
  if i:
   _,_,summary=fv.audit_step(g,c_native[i-1],c,phi,rho=RHO,gamma=GAMMA,dt_s=DT,boundary=boundary,normalization_time_s=2.,global_floor_kg_m_s=.02)
   budgets.append(dict(index=i,time_s=i*DT,actual_native_scalar_binding=True,native_previous_scalar_path=states[i-1]['global_tracer_path'],native_current_scalar_path=str(scalar_path),fv_helper_return=summary,global_reconstructed_budget_pass=summary['global_Euler_defect_relative']<=1e-5))
   for name,b in summary['reconstructed_boundary'].items():cumulative[name]['advective_inward_kg_per_m']+=DT*b['advective_inward_kg_m_s'];cumulative[name]['diffusive_inward_kg_per_m']+=DT*b['diffusive_inward_kg_m_s']
 rank_cells=[];rank_edges=[]
 for rank in (0,1):
  rp=CASE/('processor'+str(rank));rpoly=old.native_poly(rp);pilot_poly=old.native_poly(PILOT/('processor'+str(rank)))
  for name in ('points','faces','owner','neighbour','boundary','cellProcAddressing','faceProcAddressing','pointProcAddressing'):require(sha(e.touch(rpoly/name))==sha(e.touch(pilot_poly/name)),'Native rank mesh/addressing changed')
  n,body=old.counted(rpoly/'cellProcAddressing');cellmap=np.asarray(list(map(int,body.split())),dtype=int);require(len(cellmap)==n and len(np.unique(cellmap))==n and np.all((cellmap>=0)&(cellmap<4400)),'Invalid native rank cell addressing');rank_cells.extend(cellmap.tolist())
  nf,body=old.counted(rpoly/'faceProcAddressing');facemap=np.asarray(list(map(int,body.split())),dtype=int);faceids=abs(facemap)-1;require(len(facemap)==nf and len(np.unique(faceids))==nf and np.all(facemap!=0) and np.all((faceids>=0)&(faceids<len(qfull))),'Invalid signed native rank face addressing');expected_phi=qfull[faceids]*np.sign(facemap)
  require({p.name for p in rp.iterdir() if p.is_dir() and re.fullmatch(old.FLOAT,p.name)}==expected_dirs,'Missing/extra rank native states')
  n_top=prod.rank_patch_counts(rpoly)['inlet_top']
  for i in range(801):
   label=prod.label(i);folder=rp/label;path=folder/'tracer';raw=e.text(path);require(re.search(r'\bclass\s+volScalarField\s*;',raw) and re.search(r'\bobject\s+tracer\s*;',raw),'Native rank tracer object differs');cr=field(path,n,1,[0]*7);require(np.max(abs(cr-c_native[i,cellmap]))<=1e-12,'Native rank/reconstructed scalar differs')
   require(np.max(abs(field(folder/'U',n,3,[0,1,-1,0,0,0,0])-baseU[cellmap]))<=1e-12 and np.max(abs(field(folder/'p',n,1,[0,2,-2,0,0,0,0])-basep[cellmap]))<=1e-12,'Native rank frozen U/p differs')
   require(np.max(abs(np.asarray(q.all_phi(e.touch(folder/'phi'),rpoly))-expected_phi))<=1e-12,'Native rank signed all-face phi differs')
   states[i]['rank_tracer_paths'].append(str(path));states[i]['rank_tracer_sha256'].append(e.rows[path.resolve()])
   if i:
    cp=folder/'uniform/time';row=dict(scope='processor'+str(rank),expected_index=i,path=str(cp),**clock(e.text(cp),i,label));clocks.append(row);states[i]['rank_clock_paths'].append(str(cp));states[i]['rank_clock_sha256'].append(e.rows[cp.resolve()])
   if i in (80,81) and n_top:
    values=q.patch_value(q.blocks(raw)['inlet_top'],n_top);require(max(abs(v-prod.expected_top(i)) for v in values)<=1e-12,'Actual rank80/81 inlet differs');rank_edges.append(dict(rank=rank,index=i,time_s=i*DT,top_faces=n_top,tracer_path=str(path),tracer_sha256=e.rows[path.resolve()],top_c_min=min(values),top_c_max=max(values)))
 require(sorted(rank_cells)==list(range(4400)) and len(clocks)==2400,'Actual rank bijection/global and rank clocks incomplete')
 for i in (80,81):require(sum(row['top_faces'] for row in rank_edges if row['index']==i)==20,'Actual rank80/81 top boundary coverage incomplete')
 checks['exact801_native_states_and795_new_solve_endpoints']=True;checks['accepted183_pilot_files_and_original_flow_preserved']=True;checks['native_mesh_and_both_rank_fields_addressing_clocks_agree']=True
 checks['all_frozen_U_p_all9040_active_phi_within1e_12']=all(max(row['U_max_abs_change'],row['p_max_abs_change'],row['all_native_phi_max_abs_change'])<=1e-12 for row in frozen)
 checks['initial_scalar_zero_within1e_12']=bool(np.max(abs(c_native[0]))<=1e-12);checks['all_native_scalars_finite_and_bounds1e_8']=bool(np.isfinite(c_native).all() and c_native.min()>=-1e-8 and c_native.max()<=1+1e-8);checks['genuine_transport']=bool(np.any(c_native[1:]!=c_native[0]))
 checks['all800_own_native_matrix_residuals_le1e_8']=all(row['final_residual']<=1e-8 for row in residuals);checks['actual_global_and_rank80top1_81top0']=True;checks['all800_reconstructed_global_Euler_defects_le1e_5']=all(row['global_reconstructed_budget_pass'] for row in budgets)
 geometry_path=e.touch(ROOT/'native/runs/2026-11-01-tjunction-v001/mesh/geometry.npz')
 with np.load(geometry_path,allow_pickle=False) as z:canonical={name:z[name].copy() for name in ('nodes','quads','centres','areas')}
 c_source=np.empty_like(c_native);c_source[:,mapping]=c_native
 source_to_native=np.argsort(mapping);require(np.array_equal(mapping[source_to_native],np.arange(4400)),'Actual inverse native/source mapping differs')
 source_centres_m=g['cell_centres'][source_to_native];source_areas_m2=(g['volumes']/g['depth_m'])[source_to_native]
 require(np.max(abs(source_centres_m-canonical['centres']))<=1e-12 and np.max(abs(source_areas_m2-canonical['areas']))<=1e-14,'Authenticated native/source polygon geometry differs')
 manifest=e.finish();checks['all_consumed_source_native_bytes_unchanged']=True
 result=dict(native_case=str(CASE),original_native_case=str(ORIGINAL_CASE),postprocessing_recovery=recovery_summary,postprocessing_recovery_evidence=recovery_details,raw_producer_flags=dict(run_status=status,execution_checks=execution,sync_receipt=sync),raw_native_execution_checks_available=execution is not None,raw_native_execution_checks_availability='Original producer file if present; otherwise unavailable because wrapper stopped during reconstruction, never fabricated',native_scalar_field_name='tracer',mathematical_scalar_symbol='c',actual_steps=800,actual_state_count=801,copied_accepted_pulse_steps=5,new_production_steps=795,dt_s=DT,final_time_s=20.,production_targets_sha256=TARGET_SHA,analysis_source_sha256=sha(Path(__file__)),checks=checks,failed_checks=[k for k,v in checks.items() if v is not True],all_declared_numerical_checks_passed=all(v is True for v in checks.values()),root_acceptance_performed=False,additional_CFD_authorized_by_this_audit=False,physical_validation=False,raw_native_run_status=status,raw_native_execution_checks=execution,raw_native_sync_receipt=sync,raw_root_pilot_acceptances=acceptances,raw_pending_flags_preserved=True,mesh_quality=mesh_quality,actual_native_residuals=residuals,residual_norm='Qualified DILUPBiCGStab recursive residual and native uniform-reference lduMatrix normFactor; no cross-code residual ratio',native_clocks=clocks,native_state_paths=states,actual_pulse_edge_evidence=edge,actual_rank_pulse_edge_evidence=rank_edges,reconstructed_scalar_steps=budgets,reconstructed_cumulative_boundary_inventory=cumulative,frozen_native_field_changes=frozen,native_cell_Co_max=Co_max,Co_is_descriptive_not_accuracy_gate=True,local_reconstructed_defect_has_no_pass_threshold=True,scalar_evidence=dict(observed_native_active_phi_faces=9040,excluded_empty_phi_faces=8800,helper_empty_phi_slots='q.all_phi fills8800 empty-face slots with zero placeholders by convention; those are not observed native flux. Scalar budgets use only9040 active faces.',current='Actual native tracer at801 global and both rank states; no clipping',density='Declared constant rho1 mass conversion; no persisted native incompressible rho field claimed',diffusivity='Source-bound D0.001 and Gamma=rhoD0.001; no saved diffusivity-array claim',native_M1_available=False,old_time='Consecutive actually saved native current states; internal solver old-time array not persisted or invented',native_total_scalar_matrix_flux_available=False,scalar_fluxes='RECONSTRUCTED orthogonal Euler/upwind advection and diffusion from native fields, geometry and phi',inventory='kg/m of dimensionless mass-weighted scalar; not calibrated dye mass',initial_state='Accepted unsolved initial c0; no solved conservation check at step0',source_endpoint_fixture='Actual80/81 boundary values and clocks checked here; Table source fixture remains distinct from executed fields',washout='20 seconds is arrival/early washout; no complete-washout, periodicity, mesh/time-independence or validation claim'),available_RAM_at_audit_bytes=available,geometry_identity_report=geometry_report)
 out.mkdir();json_new(out/'checks.json',result);json_new(out/'source-manifest.json',manifest);json_new(out/'native-state-manifest.json',dict(source_order='Canonical original T-junction quad order',native_scalar='tracer',states=states));json_new(out/'reconstructed-step-budgets.json',dict(classification='RECONSTRUCTED, not native matrix flux',steps=budgets,cumulative_boundary_inventory=cumulative))
 np.savez_compressed(out/'fields.npz',step=np.arange(801,dtype=np.int64),time_s=np.arange(801)*DT,c=c_source,c_native=c_native,native_cell_to_source_cell=mapping,source_to_native_cell=source_to_native,source_native_cell_index=source_to_native,source_centres_m=source_centres_m,source_areas_m2=source_areas_m2,native_cell_centres=g['cell_centres'],native_mass_phi=phi,**canonical)
 rows=[]
 for state,budget,residual in zip(states[1:],budgets,residuals):
  s=budget['fv_helper_return'];row=dict(index=state['index'],time_s=state['time_s'],scalar_min=state['scalar_min'],scalar_max=state['scalar_max'],negative_cells=state['negative_cells'],above_one_cells=state['above_one_cells'],actual_top_c_min=state['actual_top_c_min'],actual_top_c_max=state['actual_top_c_max'],inventory_kg_per_m=s['inventory_kg_per_m'],storage_rate_kg_m_s=s['storage_rate_kg_m_s'],global_Euler_defect_raw_kg_m_s=s['global_Euler_defect_raw_kg_m_s'],global_Euler_defect_denominator_kg_m_s=s['global_Euler_defect_denominator_kg_m_s'],global_Euler_defect_relative=s['global_Euler_defect_relative'],local_reconstructed_defect_max=s['reconstructed_local_dimensionless_defect_max'],local_reconstructed_defect_RMS=s['reconstructed_local_dimensionless_defect_volume_RMS'],outlet_outflow_weighted_mean_c=s['reconstructed_outlet_outflow_weighted_mean_c'],final_native_matrix_residual=residual['final_residual'])
  for name,b in s['reconstructed_boundary'].items():
   for key,value in b.items():row[name+'_'+key]=value
  rows.append(row)
 with (out/'step-history.csv').open('x',encoding='ascii',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 for name in (QA_NAME,'scalar_fv_audit_v002.py','scalar_native_geometry_v003.py','prepare_openfoam_tracer_qualification_v004.py','prepare_openfoam_tracer_production_v001.py'):shutil.copy2(PRIVATE/name,out/('executed-'+name))
 shutil.copy2(Path(__file__),out/'executed-auditor.py');json_new(out/'output-manifest.json',dict(outputs=[dict(name=p.name,bytes=p.stat().st_size,sha256=sha(p)) for p in out.iterdir() if p.is_file()]))
 return dict(out=str(out),checks_passed=sum(v is True for v in checks.values()),check_count=len(checks),all_declared_numerical_checks_passed=result['all_declared_numerical_checks_passed'],maximum_reconstructed_global_defect=max(row['fv_helper_return']['global_Euler_defect_relative'] for row in budgets),actual_steps=800,actual_states=801,new_production_steps=795,root_acceptance_performed=False,physical_validation=False)
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--authorization',type=Path,required=True);p.add_argument('--recovery-authorization',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();print(json.dumps(audit(a.authorization,a.out,a.recovery_authorization),indent=2,allow_nan=False))
