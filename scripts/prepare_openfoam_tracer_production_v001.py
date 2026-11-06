"""SOURCE ONLY until exact root admission: continue the accepted pulse5 by795 native scalar steps."""
from pathlib import Path, PurePosixPath
import argparse, hashlib, importlib.util, json, math, os, re, shutil, sys
ROOT=Path('E:/CFD-Learning-Lab' if os.name=='nt' else '/mnt/e/CFD-Learning-Lab')
PRIVATE=ROOT/'local/prospective-tracer-v001'
RUN=ROOT/'native/runs/2026-11-02-tracer-tjunction-v001'
PILOT=RUN/'openfoam-pulse-pilot-v003'
CASE=RUN/'openfoam-production-v001'
TARGET_SHA='07f5b0745fc71758b9cf5678171ce17fd558f00a01cd773560911920e91318ba'
QUAL_PINS={'prepare_openfoam_tracer_qualification_v004.py':'278e6948a244e82887d8eb8dc2da90b502a14a7165de2238c3ba1e01b2890888','run-openfoam-tracer-qualification-v004.sh':'a4e4437d8ccd241c724d99148a9a9b3b4c4764785ad01e6b516ea08187082c41','Launch-OpenFOAM-TracerQualification-v004.ps1':'e4c7f96819c3a4bdf25681c68c517e5ec7a107866529f2b171c5540096ec268c'}
FV_PINS={'scalar_fv_audit_v002.py':'75ceec60add2f5c53fb82585d105dcc5cd7c9d8d19e1ea8ac3f2dde14aa455cd','scalar_native_geometry_v003.py':'06f6b4089116d8e7cdab8b4050c30d2939d3876381ebbaa878c7908e6f2d1903'}
ACCEPTANCE_PINS={'root-Fluent-qualification-acceptance-v001.json':'93a4a8c894b1ef3241abaa0508ce4f423405249ff57b04ffe15ad2155977fe6d','root-OpenFOAM-qualification-acceptance-v001.json':'7ea48c3ab6e451314a96bef8b40df9857b1a32162d39e0f504a5e80cbca89992'}
PILOT_AUTH_NAME='root-OpenFOAM-Pulse-qualification-authorization-v001.json'
PILOT_AUTH_PIN='1d960967f83b647a00ece6ebf2debf0fccfa03cdadea3f30110c371368ca3a8e'
SOURCES=('prepare_openfoam_tracer_production_v001.py','run-openfoam-tracer-production-v001.sh','Launch-OpenFOAM-TracerProduction-v001.ps1')
KEYS=('preparation_sha256','wrapper_sha256','launcher_sha256')
NEW_STEPS=tuple(range(6,801)); STATE_STEPS=tuple(range(801)); DT=.025
OVERRIDDEN={'system/controlDict','run-status.json'}

def require(v,m):
 if not v:raise ValueError(m)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def readj(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write_new(p,text):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('x',encoding='ascii',newline='\n') as f:f.write(text.rstrip()+'\n')
def json_new(p,v):write_new(p,json.dumps(v,indent=2,allow_nan=False))
def project_path(raw):
 text=str(raw).replace('\\','/')
 for prefix in ('E:/CFD-Learning-Lab','/mnt/e/CFD-Learning-Lab'):
  if text.startswith(prefix+'/'):
   p=(ROOT/text[len(prefix)+1:]).resolve();require(ROOT.resolve() in p.parents,'Path escapes project');return p
 raise ValueError('Expected exact project prefix '+text)
for name,pin in QUAL_PINS.items():require(sha(PRIVATE/name)==pin,'Reviewed qualifier source changed '+name)
spec=importlib.util.spec_from_file_location('frozen_tracer_qualifier_v004',PRIVATE/'prepare_openfoam_tracer_qualification_v004.py')
q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)

def label(i):return f'{i*DT:.10g}'
def expected_top(i):return 1 if i<=80 else 0

def validate_relative_paths(paths):
 seen=set()
 for text in paths:
  p=PurePosixPath(text);require(text==p.as_posix() and not p.is_absolute() and text and all(v not in ('','..','.') for v in p.parts),'Unsafe relative path '+text)
  folded=text.casefold();require(folded not in seen,'Case-folding path collision '+text);seen.add(folded)

def tree_sha(folder):
 folder=Path(folder);require(folder.is_dir() and not folder.is_symlink(),'Real source directory required')
 files=[]
 for p in folder.rglob('*'):
  require(not p.is_symlink(),'Native tree symlinks unsupported')
  if p.is_file():files.append(p)
 result={p.relative_to(folder).as_posix():sha(p) for p in sorted(files)};validate_relative_paths(result);return result

def verify_tree(folder,pins):require(tree_sha(folder)==pins,'Complete accepted pilot tree changed '+str(folder))
def copy_pilot_tree(pilot,case,pins):
 require(not Path(case).exists(),'Fresh production directory required');verify_tree(pilot,pins)
 Path(case).mkdir(parents=True)
 for name in sorted(pins):q.copy_new(Path(pilot)/name,Path(case)/name)
 verify_tree(pilot,pins);verify_tree(case,pins)

def archive_then_replace(case,relative,text):
 path=Path(case)/relative;archive=Path(case)/'production-provenance/pilot-overridden'/relative
 q.copy_new(path,archive)
 with path.open('w',encoding='ascii',newline='\n') as f:f.write(text.rstrip()+'\n')

def production_control():
 raw=q.dictionaries('Pulse')['system/controlDict'];require(raw.count('startTime 0;')==1 and raw.count('endTime 0.125;')==1,'Pinned restart literals differ')
 return raw.replace('startTime 0;','startTime 0.125;').replace('endTime 0.125;','endTime 20;')

def native_times(raw):
 values=[float(v) for v in re.findall(r'^Time = ('+q.old.FLOAT+r')s?\s*$',raw,re.M)]
 require(len(values)==len(NEW_STEPS) and all(abs(t-i*DT)<=1e-12 for i,t in zip(NEW_STEPS,values)),'Exactly795 new endpoints6..800 required; no pilot replay/extension')
 return values

def accepted_pilot_pins():
 for name,pin in ACCEPTANCE_PINS.items():require(sha(PRIVATE/name)==pin,'Actual root acceptance changed '+name)
 f=readj(PRIVATE/'root-Fluent-qualification-acceptance-v001.json');o=readj(PRIVATE/'root-OpenFOAM-qualification-acceptance-v001.json')
 require(f.get('reviewer')=='root' and f.get('independent_all10_checks_passed') is True and f.get('qualified_control_steps')==5 and f.get('qualified_pulse_steps')==5,'Actual both-phase Fluent qualification acceptance differs')
 require(o.get('reviewer')=='root' and o.get('both_phases_accepted') is True and o.get('selected_control')=='openfoam-control-v003' and o.get('selected_pulse')=='openfoam-pulse-pilot-v003' and o.get('qualified_control_steps')==o.get('qualified_pulse_steps')==5,'Actual new separate OF control/pulse acceptance differs')
 rows=o['pulse_ancestry_file_manifest'];pins={r['relative_path']:r['sha256'] for r in rows};require(len(rows)==len(pins)==183,'Exact actual accepted183-file pulse ancestry required')
 for row in rows:require((PILOT/row['relative_path']).stat().st_size==row['bytes'],'Accepted pilot file byte count differs')
 require(sha(PRIVATE/PILOT_AUTH_NAME)==PILOT_AUTH_PIN,'Actual root pilot run authorization changed')
 return pins

def authorization(path):
 path=Path(path).resolve();require(PRIVATE.resolve() in path.parents,'Authorization outside exact candidate folder');a=readj(path)
 require(a.get('approved') is True and a.get('approved_by')=='root' and a.get('action')=='run_tracer_openfoam_production_v001','Exact root production authorization required')
 require(a.get('additional_steps')==795 and a.get('pilot_steps')==5 and a.get('total_steps')==800 and a.get('total_states')==801 and a.get('start_time_s')==.125 and a.get('final_time_s')==20 and a.get('deltaT_s')==DT and a.get('processes')==2,'Frozen production budget differs')
 require(a.get('fluent_sessions_closed') is True and a.get('fluent_qualification_accepted_by_root') is True and a.get('openfoam_control_accepted_by_root') is True and a.get('openfoam_pulse_accepted_by_root') is True,'Both codes and new separate OF controls/pilots must be root-accepted before production')
 require(a.get('production_targets_sha256')==TARGET_SHA and sha(PRIVATE/'PRODUCTION-TARGETS-v001.md')==TARGET_SHA,'Frozen production targets changed')
 require(a.get('native_flow_qualification_sha256')==q.QUAL_SHA,'Original flow ancestry binding differs')
 require(a.get('frozen_field_absolute_tolerance')==1e-12 and a.get('scalar_bounds_tolerance')==1e-8 and a.get('final_matrix_residual_target')==1e-8 and a.get('reconstructed_global_budget_target')==1e-5,'Scientific gates must remain unchanged')
 require(a.get('pulse_continuum_edge_s')==2.0125 and a.get('pulse_step80_value')==1 and a.get('pulse_step81_value')==0,'Frozen sampled pulse schedule required')
 for key,name in zip(KEYS,SOURCES):require(a.get(key)==sha(PRIVATE/name),'Reviewed production source changed '+name)
 require(a.get('qualification_source_sha256')==QUAL_PINS and a.get('fv_helper_sha256')==FV_PINS,'Full reviewed qualifier/FV source binding differs')
 for name,pin in FV_PINS.items():require(sha(PRIVATE/name)==pin,'Reviewed FV helper changed '+name)
 require(set(a.get('installed_source_sha256',{}))==set(q.INSTALLED),'Exact installed19 source coverage required')
 require(all(re.fullmatch('[0-9a-f]{64}',v) for v in a['installed_source_sha256'].values()),'Invalid installed-source hashes')
 require(all(a['installed_source_sha256'][n]==d for n,d in q.LINEAR_SOURCE_SHA256.items()),'Actual linear solver source differs')
 require(a.get('pilot_case')=='native/runs/2026-11-02-tracer-tjunction-v001/openfoam-pulse-pilot-v003' and a.get('production_case')=='native/runs/2026-11-02-tracer-tjunction-v001/openfoam-production-v001','Exact fresh model/pulse lineage required')
 pilot_auth=a.get('pilot_run_authorization',{});pa=project_path(pilot_auth['path']);require(pa==PRIVATE/PILOT_AUTH_NAME and pilot_auth['sha256']==PILOT_AUTH_PIN and sha(pa)==PILOT_AUTH_PIN,'Exact actual private pilot authorization binding required')
 pins=a.get('pilot_file_sha256',{});require(pins and all(re.fullmatch('[0-9a-f]{64}',v) for v in pins.values()),'Actual complete accepted pilot file manifest required');validate_relative_paths(pins);require(pins==accepted_pilot_pins(),'Production must bind exact actual root-accepted183 pilot files')
 evidence=a.get('acceptance_evidence',{});require(set(evidence)=={'Fluent','OpenFOAMControl','OpenFOAMPulse'},'Three actual root acceptance evidence bindings required')
 for role,row in evidence.items():
  p=project_path(row['path']);require(p.is_file() and re.fullmatch('[0-9a-f]{64}',row['sha256']) and sha(p)==row['sha256'],'Actual acceptance evidence changed '+role)
  require(PRIVATE.resolve() in p.parents,'Acceptance evidence must remain in candidate private folder')
  expected='root-Fluent-qualification-acceptance-v001.json' if role=='Fluent' else 'root-OpenFOAM-qualification-acceptance-v001.json'
  require(p==PRIVATE/expected and row['sha256']==ACCEPTANCE_PINS[expected],'Actual selected acceptance role differs '+role)
 ancestry=a.get('selected_ancestry_sha256',{});require(ancestry and all(re.fullmatch('[0-9a-f]{64}',v) for v in ancestry.values()),'Actual full selected ancestry/evidence hashes required')
 for raw,pin in ancestry.items():require(sha(project_path(raw))==pin,'Selected ancestry changed '+raw)
 return path,a

def verify_pilot(a):
 verify_tree(PILOT,a['pilot_file_sha256']);q.verify_sources()
 require(not any(n=='production-provenance' or n.startswith('production-provenance/') or '/' not in n and n.startswith('production-') for n in a['pilot_file_sha256']),'Pilot contains reserved production output names')
 p=project_path(a['pilot_run_authorization']['path']);require(sha(p)==sha(PILOT/'provenance/root-run-authorization.json')==a['pilot_run_authorization']['sha256'],'Actual pilot private/copied run authorization differs');pa=q.readj(p)
 require(pa['installed_source_sha256']==a['installed_source_sha256'],'Selected pilot installed19 source map differs')
 require((PILOT/'system/fvSolution').read_text()==q.dictionaries('Pulse')['system/fvSolution'],'Selected pilot did not qualify exact PBiCGStab/DILU configuration')
 result=q.verify('Pulse',p,True);require(result['completed_five_steps'] is True and result['states_verified']==6,'Actual completed fresh pulse5 required')
 verify_tree(PILOT,a['pilot_file_sha256']);return result

def prepare(authpath):
 authpath,a=authorization(authpath);require(shutil.disk_usage(ROOT).free>=10*1024**3,'Need10GiB E disk');require(not CASE.exists(),'Fresh versioned production case required')
 pilot_result=verify_pilot(a);copy_pilot_tree(PILOT,CASE,a['pilot_file_sha256'])
 archive_then_replace(CASE,'system/controlDict',production_control())
 archive_then_replace(CASE,'run-status.json',json.dumps(dict(stage='prepared_only',pilot_steps=5,additional_steps_authorized=795,total_steps=800,solver_started=False,full_signed_scalar_budget_analysis_pending=True,physical_validation=False),indent=2))
 for name in SOURCES:q.copy_new(PRIVATE/name,CASE/'production-provenance'/name)
 q.copy_new(authpath,CASE/'production-provenance/root-production-authorization.json');q.copy_new(PRIVATE/'PRODUCTION-TARGETS-v001.md',CASE/'production-provenance/PRODUCTION-TARGETS-v001.md')
 json_new(CASE/'production-provenance/accepted-pilot-all-file-sha256.json',a['pilot_file_sha256'])
 json_new(CASE/'production-provenance/actual-five-step-pilot-preparation-readback.json',pilot_result)
 preservation(a);return dict(case=str(CASE),prepared_only=True,pilot_steps=5,additional_steps_authorized=795,total_steps=800,solver_started=False)

def preservation(a):
 verify_tree(PILOT,a['pilot_file_sha256'])
 for name,pin in a['pilot_file_sha256'].items():
  p=CASE/'production-provenance/pilot-overridden'/name if name in OVERRIDDEN else CASE/name
  require(sha(p)==pin,'Preserved pilot byte changed '+name)
 require((CASE/'system/controlDict').read_text()==production_control(),'Reviewed restart control changed')
 for name,text in q.dictionaries('Pulse').items():
  if name not in ('system/controlDict','0/tracer'):require((CASE/name).read_text()==text,'Reviewed scalar dictionary changed '+name)
 for name,pin in a['installed_source_sha256'].items():require(sha(CASE/'provenance/installed-sources'/name)==pin,'Preserved installed source changed '+name)
 for name,pin in QUAL_PINS.items():require(sha(CASE/'provenance'/name)==pin,'Preserved qualifier source changed '+name)
 for name in SOURCES:require(sha(CASE/'production-provenance'/name)==sha(PRIVATE/name),'Executed production source changed '+name)
 q.verify_sources()

def rank_patch_counts(poly):
 n,body=q.old.counted(poly/'boundary');rows={name:int(re.search(r'\bnFaces\s+(\d+)\s*;',b)[1]) for name,b in re.findall(r'([A-Za-z_][A-Za-z_0-9]*)\s*\{([^{}]*)\}',body)}
 require(len(rows)==n and 'inlet_top' in rows,'Rank boundary coverage differs');return rows

def native_clock(path,index):
 raw=path.read_text();v={n:re.search(r'\b'+n+r'\s+('+q.old.FLOAT+r')\s*;',raw) for n in ('value','index','deltaT')};require(all(v.values()),'Missing native clock '+str(path));values={n:float(m[1]) for n,m in v.items()}
 require(all(map(math.isfinite,values.values())) and abs(values['value']-index*DT)<=1e-12 and values['index']==index and abs(values['deltaT']-DT)<=1e-12,'Native clock/index/dt differs '+str(path));return values

def verify(authpath,completed=False):
 _,a=authorization(authpath);preservation(a);q.old.mesh_quality((CASE/'production-mesh-check.log').read_text());poly=q.old.native_poly(CASE);patches=q.old.native_patches(poly);geom=q.old.native_geometry(poly)
 for name in ('points','faces','owner','neighbour','boundary'):require(sha(poly/name)==sha(q.old.native_poly(PILOT)/name),'Native production mesh changed '+name)
 require(q.difference(q.field(CASE/'0/Vc',4400,1,[0,3,0,0,0,0,0]),geom['volumes'])<=1e-14,'Native Vc differs from actual mesh')
 require(q.difference(q.field(CASE/'0/C',4400,3,[0,1,0,0,0,0,0]),q.field(PILOT/'0/C',4400,3,[0,1,0,0,0,0,0]))<=1e-12,'Native source cell centre order changed')
 dims={'U':(3,[0,1,-1,0,0,0,0]),'p':(1,[0,2,-2,0,0,0,0])};base={n:q.field(CASE/'provenance/baseline-native-inputs'/n,4400,k,d) for n,(k,d) in dims.items()};q0=q.all_phi(CASE/'provenance/baseline-native-inputs/phi',poly)
 outlet=patches['outlet'];require(all(v>=0 for v in q0[outlet['start']:outlet['start']+outlet['count']]),'Unsupported frozen outlet backflow')
 steps=STATE_STEPS if completed else tuple(range(6));states={};history=[]
 for i in steps:
  folder=CASE/label(i)
  for name,(k,d) in dims.items():require(q.difference(q.field(folder/name,4400,k,d),base[name])<=1e-12,'Frozen '+name+' changed at'+label(i))
  require(q.difference(q.all_phi(folder/'phi',poly),q0)<=1e-12,'Frozen global all-face phi changed at'+label(i))
  raw=(folder/'tracer').read_text();require(re.search(r'\bclass\s+volScalarField\s*;',raw) and re.search(r'\bobject\s+tracer\s*;',raw),'Native tracer class/object differs')
  c=q.field(folder/'tracer',4400,1,[0,0,0,0,0,0,0]);require(min(c)>=-1e-8 and max(c)<=1+1e-8,'Actual scalar bounds missed at'+label(i));require(i==0 or max(map(abs,c))>0,'Actual transported scalar unchanged');states[i]=c
  bs=q.blocks(raw);require(set(bs)==set(patches),'Actual scalar boundary coverage differs')
  for name in ('outlet','walls'):require(re.search(r'\btype\s+zeroGradient\s*;',bs[name]),'Scalar zero diffusion changed '+name)
  require(re.search(r'\btype\s+empty\s*;',bs[q.old.EMPTY]),'Empty planes changed')
  require(re.search(r'\btype\s+fixedValue\s*;',bs['inlet_left']) and re.search(r'\btype\s+uniformFixedValue\s*;',bs['inlet_top']),'Native inlet scalar type changed')
  for name,expected in [('inlet_left',0),('inlet_top',expected_top(i))]:require(max(abs(v-expected) for v in q.patch_value(bs[name],20))<=1e-12,'Actual executed scalar inlet differs '+name+' at'+label(i))
  if i==0:require(max(map(abs,c))<=1e-12,'Actual pulse initial state changed')
  if (folder/'uniform/time').is_file() and i:native_clock(folder/'uniform/time',i)
  top_values=q.patch_value(bs['inlet_top'],20)
  history.append(dict(index=i,time_s=i*DT,actual_top_c_min=min(top_values),actual_top_c_max=max(top_values),min_c=min(c),max_c=max(c),negative_count=sum(v<0 for v in c),above_one_count=sum(v>1 for v in c),inventory_kg_per_m=sum(v*x for v,x in zip(geom['volumes'],c))/.01,tracer_sha256=sha(folder/'tracer')))
 if completed:
  raw=(CASE/'production-solver.log').read_text();require(re.search(r'^End\s*$',raw,re.M) and not re.search(r'FOAM FATAL|Floating point exception|Segmentation fault',raw),'Native production solver did not close cleanly');native_times(raw)
  require(not re.search(r'Solving for (?:Ux|Uy|Uz|p),',raw),'Unexpected flow update')
  residuals=re.findall(r'DILUPBiCGStab:\s+Solving for tracer, Initial residual = ('+q.old.FLOAT+r'), Final residual = ('+q.old.FLOAT+r'), No Iterations (\d+)',raw)
  require(len(residuals)==795 and all(math.isfinite(float(v[0])) and math.isfinite(float(v[1])) and 0<=float(v[1])<=1e-8 and 0<=int(v[2])<=1000 for v in residuals),'Exactly795 native scalar matrix residual gates required')
 rank_cells=[];clock_count=0;rank_edge_evidence=[]
 for rank in (0,1):
  rp=CASE/('processor'+str(rank));rpoly=q.old.native_poly(rp);nc,b=q.old.counted(rpoly/'cellProcAddressing');cellmap=list(map(int,b.split()));require(len(cellmap)==nc and len(set(cellmap))==nc and all(0<=v<4400 for v in cellmap),'Invalid rank cell mapping');rank_cells+=cellmap
  nf,b=q.old.counted(rpoly/'faceProcAddressing');facemap=list(map(int,b.split()));require(len(facemap)==nf and all(v and abs(v)<=len(q0) for v in facemap),'Invalid signed rank face mapping');expected=[q0[abs(v)-1]*(1 if v>0 else -1) for v in facemap]
  for name in ('points','faces','owner','neighbour','boundary','cellProcAddressing','faceProcAddressing','pointProcAddressing'):require(sha(rpoly/name)==sha(q.old.native_poly(PILOT/('processor'+str(rank)))/name),'Rank native mesh/addressing changed '+name)
  for i in steps:
   folder=rp/label(i)
   for name,(k,d) in dims.items():require(q.difference(q.field(folder/name,nc,k,d),[base[name][g] for g in cellmap])<=1e-12,'Actual rank frozen '+name+' changed')
   require(q.difference(q.all_phi(folder/'phi',rpoly),expected)<=1e-12,'Actual rank frozen all-face phi changed')
   require(q.difference(q.field(folder/'tracer',nc,1,[0,0,0,0,0,0,0]),[states[i][g] for g in cellmap])<=1e-12,'Actual rank/reconstructed tracer differs')
   if i in (80,81):
    rpatch=q.blocks((folder/'tracer').read_text());n_top=rank_patch_counts(rpoly)['inlet_top']
    if n_top:
     top=q.patch_value(rpatch['inlet_top'],n_top);require(max(abs(v-expected_top(i)) for v in top)<=1e-12,'Actual rank executed80/81 inlet differs')
     rank_edge_evidence.append(dict(rank=rank,index=i,time_s=i*DT,top_faces=n_top,tracer_sha256=sha(folder/'tracer'),top_c_min=min(top),top_c_max=max(top)))
   if i:native_clock(folder/'uniform/time',i);clock_count+=1
  if completed:
   numeric={p.name for p in rp.iterdir() if p.is_dir() and re.fullmatch(q.old.FLOAT,p.name)};require(numeric=={label(i) for i in STATE_STEPS}|{'2000'},'Unexpected/missing native rank state directories')
 require(sorted(rank_cells)==list(range(4400)),'Two ranks must biject to all4400 cells')
 if completed:require(clock_count==1600,'Require both rank native clocks1..800');preservation(a)
 return dict(completed_production_execution=completed,states_verified=len(steps),pilot_steps=5,additional_steps_verified=795 if completed else 0,total_steps=800,deltaT_s=DT,physical_time_s=20 if completed else .125,processes=2,native_rank_clock_count=clock_count,frozen_U_p_all_phi_unchanged=True,native_mesh_unchanged=True,accepted_pilot_all_file_bytes_preserved=True,actual_step80_top_c=history[80]['actual_top_c_min'] if completed else None,actual_step81_top_c=history[81]['actual_top_c_max'] if completed else None,actual_rank_pulse_edge_evidence=rank_edge_evidence,actual_pulse_edge_evidence=[dict(index=i,time_s=i*DT,tracer_path=str(CASE/label(i)/'tracer'),tracer_sha256=history[i]['tracer_sha256'],top_c_min=history[i]['actual_top_c_min'],top_c_max=history[i]['actual_top_c_max']) for i in (80,81)] if completed else [],history=history,full_signed_scalar_budget_analysis_pending=True,root_acceptance_performed=False,physical_validation=False)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--authorization',type=Path,required=True);mode=p.add_mutually_exclusive_group(required=True);mode.add_argument('--prepare',action='store_true');mode.add_argument('--verify-sources',action='store_true');mode.add_argument('--verify-baseline',action='store_true');mode.add_argument('--verify-completed',action='store_true');args=p.parse_args();_,a=authorization(args.authorization)
 if args.prepare:r=prepare(args.authorization)
 elif args.verify_sources:verify_pilot(a);r=dict(pilot_and_selected_sources_verified=True,solver_started=False)
 else:r=verify(args.authorization,args.verify_completed)
 print(json.dumps(r,indent=2,allow_nan=False))
if __name__=='__main__':main()
