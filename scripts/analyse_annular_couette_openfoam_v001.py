"""Read-only rank-native five-step annular Couette audit. Never launches CFD or fabricates native global fields."""
from pathlib import Path
import argparse,csv,json,re,math
import numpy as np
import prepare_openfoam_annular_couette_v002 as prep
require=prep.require;sha=prep.sha
PREP_SHA='9e07e2982beca3a34242478bb23bf98e7d8ccb14b00836e69831ffc713bdfe34'
VALIDATOR_SHA='84bdaab73ca3205b07ec0aaa0162052659fbb386c265dc31f88b4861de0e2276'
RECOVERY_SHA='f3f551767965f01ada7b454029392a2cdfea20ea66061f03670c43d2a1a82796'
PREFLIGHT_SHA='ea061f0484e7956db7699ed34b8a9222d18845cafa02dc979c5a6287e26f0ded'
MAPPING_SHA='81034e2d28081c15cb16348dfe9ef426ab04fa04b281d5a7b8f400e5c8e2ebe4'
class Evidence:
 def __init__(self):self.before={}
 def touch(self,p):
  p=Path(p).resolve();require(p.is_file(),'Missing native evidence '+str(p))
  if p not in self.before:self.before[p]=sha(p)
  return p
 def text(self,p):return self.touch(p).read_text(encoding='utf-8-sig')
 def json(self,p):return json.loads(self.text(p))
 def finish(self):
  rows=[dict(path=str(p),bytes=p.stat().st_size,sha256_before=h,sha256_after=sha(p)) for p,h in self.before.items()];require(all(r['sha256_before']==r['sha256_after'] for r in rows),'Read-only input changed');return dict(sources=rows,inputs_unchanged_during_analysis=True,solver_started_or_modified=False)
def field(raw,count,components,dims):
 t=prep.clean(raw);dm=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',t);require(dm is not None and list(map(int,dm[1].split()))==dims and re.search(r'\bformat\s+ascii\s*;',t),'Actual ASCII field dimensions');return np.asarray(prep.values(t,count,components),dtype=float)
def blocks(raw):
 t=prep.clean(raw);k=t.index('boundaryField');start=t.index('{',k)+1;out={};i=start
 while i<len(t):
  m=re.match(r'\s*(\"[^\" ]+\"|\w+)\s*\{',t[i:])
  if not m:break
  name=m[1].strip('"');a=i+m.end();j=a;depth=1
  while depth:
   require(j<len(t),'Unclosed native patch dictionary');depth+=(t[j]=='{')-(t[j]=='}');j+=1
  out[name]=t[a:j-1];i=j
 return out
def top(block):
 chars=[];depth=0
 for c in block:
  if c=='{':depth+=1
  elif c=='}':depth-=1
  elif not depth:chars.append(c)
 return ''.join(chars)
def time_meta(e,folder,index):
 t=e.text(folder/'uniform/time');d={n:float(re.search(r'\b'+n+r'\s+('+prep.FLOAT+r')\s*;',t)[1]) for n in ('value','deltaT','index')};require(abs(d['value']-index*.05)<=1e-12 and abs(d['deltaT']-.05)<=1e-14 and d['index']==index,'Native actual time/dt/index mismatch');return d
class Ranks:
 def __init__(self,e,case,g,patches):
  self.e=e;self.case=case;self.g=g;self.rows=[];self.mesh_paths=[];cells_seen=[];faces_seen=np.zeros(prep.FACES,dtype=int);physical_seen=[];points=np.asarray(g['points']);faces=np.asarray(g['faces']);owner=np.asarray(g['owner']);neighbour=np.asarray(g['neighbour']);Sf=np.asarray(g['Sf']);physical=np.concatenate([np.arange(patches[n]['start'],patches[n]['start']+patches[n]['count']) for n in prep.PATCHES]);self.physical=physical;physical_lookup={int(fi):j for j,fi in enumerate(physical)}
  for rank in (0,1):
   poly=case/f'processor{rank}/0/polyMesh';paths={n:str(poly/n) for n in ('points','faces','owner','neighbour','boundary','cellProcAddressing','faceProcAddressing','pointProcAddressing')};self.mesh_paths.append(paths)
   def labels(n):
    count,b=prep.counted(e.touch(poly/n));a=np.array(list(map(int,b.split())),dtype=np.int64);require(len(a)==count,'Rank label length');return a
   cells=labels('cellProcAddressing');addr=labels('faceProcAddressing');own=labels('owner');nei=labels('neighbour');pm=labels('pointProcAddressing');cells_seen.extend(cells.tolist());require(np.all((cells>=0)&(cells<prep.CELLS)) and len(set(cells))==len(cells),'Native rank cell addressing')
   n,b=prep.counted(e.touch(poly/'points'));rp=np.array([list(map(float,s.split())) for s in re.findall(r'\(([^()]*)\)',b)]);require(rp.shape==(n,3) and len(pm)==n and np.all((pm>=0)&(pm<prep.POINTS)) and np.max(abs(rp-points[pm]))<=1e-12,'Native rank points/coordinates')
   n,b=prep.counted(e.touch(poly/'faces'));rf=np.array([list(map(int,s.split())) for size,s in re.findall(r'(\d+)\s*\(([^()]*)\)',b)]);require(rf.shape==(n,4) and len(own)==len(addr)==n and np.all(addr!=0) and np.all(abs(addr)<=prep.FACES) and len(set(abs(addr)))==len(addr),'Rank native face addressing')
   for lf,a in enumerate(addr):
    fi=abs(int(a))-1;sign=1 if a>0 else -1;require(set(pm[rf[lf]].tolist())==set(faces[fi].tolist()),'Every rank face native vertex set');v=rp[rf[lf]];sf=.5*np.sum(np.cross(v,np.roll(v,-1,axis=0)),axis=0);require(np.max(abs(sf-sign*Sf[fi]))<=1e-14,'Native signed Sf orientation');require(cells[own[lf]]==(owner[fi] if sign>0 else neighbour[fi] if fi<prep.INTERNAL else -1),'Rank/global owner mapping')
    if lf<len(nei):require(cells[nei[lf]]==(neighbour[fi] if sign>0 else owner[fi]),'Rank/global neighbour mapping')
   faces_seen[abs(addr)-1]+=1;n,b=prep.counted(e.touch(poly/'boundary'));bp={name:dict(kind=re.search(r'\btype\s+(\w+)\s*;',q)[1],count=int(re.search(r'\bnFaces\s+(\d+)\s*;',q)[1]),start=int(re.search(r'\bstartFace\s+(\d+)\s*;',q)[1])) for name,q in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',b)};maps=[]
   for name in prep.PATCHES:
    q=bp[name];local=np.arange(q['start'],q['start']+q['count']);a=addr[local];ids=abs(a)-1;require(np.all(a>0) and set(ids).issubset(set(np.arange(patches[name]['start'],patches[name]['start']+64))),'Physical patch signed mapping');slots=np.array([physical_lookup[int(fi)] for fi in ids],dtype=int);physical_seen.extend(ids.tolist());maps.append((name,q,slots,own[local]))
   self.rows.append(dict(rank=rank,cells=cells,addr=addr,ninternal=len(nei),bp=bp,maps=maps))
  require(sorted(cells_seen)==list(range(prep.CELLS)) and np.all((faces_seen>=1)&(faces_seen<=2)),'Complete native cells/faces rank cover');require(sorted(physical_seen)==sorted(physical.tolist()),'Every128physical faces exactlyonce');self.report=dict(rank_cell_counts=[len(r['cells']) for r in self.rows],all_cells_disjoint_complete=True,all_faces_signed_addressed=True,all_native_face_vertex_Sf_owner_neighbour_maps=True,all128_physical_faces_once=True)
 def read(self,label):
  U=np.full((prep.CELLS,3),np.nan);P=np.full(prep.CELLS,np.nan);W=np.full((prep.CELLS,3),np.nan);Phi=np.full(prep.FACES,np.nan);seen=np.zeros(prep.FACES,dtype=bool);Ub=np.full((128,3),np.nan);Pb=np.full(128,np.nan);duplicate_error=0.;rank_finite=True
  for r in self.rows:
   folder=self.case/f'processor{r["rank"]}'/label;nc=len(r['cells']);ur=self.e.text(folder/'U');pr=self.e.text(folder/'p');fr=self.e.text(folder/'phi');wr=self.e.text(folder/'vorticity');u=field(ur,nc,3,[0,1,-1,0,0,0,0]);p=field(pr,nc,1,[0,2,-2,0,0,0,0]);w=field(wr,nc,3,[0,0,-1,0,0,0,0]);U[r['cells']]=u;P[r['cells']]=p;W[r['cells']]=w;phi=np.zeros(len(r['addr']));phi[:r['ninternal']]=field(fr,r['ninternal'],1,[0,3,-1,0,0,0,0]);fb=blocks(fr);ub=blocks(ur);pb=blocks(pr)
   for name,q in r['bp'].items():
    body=top(fb[name])
    if re.search(r'\btype\s+empty\s*;',body):continue
    phi[q['start']:q['start']+q['count']]=np.asarray(prep.values(body,q['count'],1,'value'))
   rank_finite=rank_finite and bool(all(np.isfinite(a).all() for a in (u,p,w,phi)));ids=abs(r['addr'])-1;oriented=np.sign(r['addr'])*phi;overlap=seen[ids]
   if overlap.any():duplicate_error=max(duplicate_error,float(np.max(abs(Phi[ids[overlap]]-oriented[overlap]))))
   Phi[ids[~overlap]]=oriented[~overlap];seen[ids]=True
   for name,q,slots,own in r['maps']:
    b=top(ub[name]);b_p=top(pb[name]);kind='rotatingWallVelocity' if name=='inner' else 'noSlip';require(re.search(r'\btype\s+'+kind+r'\s*;',b) and re.search(r'\btype\s+fixedFluxPressure\s*;',b_p),'Native wall U/p kinds');count=q['count']
    if name=='inner':require(re.search(r'\bvalue\s+',b),'Actual rotating wall numerical value required')
    v=np.asarray(prep.values(b,count,3,'value')) if re.search(r'\bvalue\s+',b) else np.zeros((count,3));Ub[slots]=v;Pb[slots]=np.asarray(prep.values(b_p,count,1,'value'))
  require(seen.all(),'Native allface coverage');return U,P,W,Phi,Ub,Pb,rank_finite,duplicate_error

def residual_history(e,case):
 text=e.text(case/'solver.log');require(re.search(r'^End\s*$',text,re.M) and 'Finalising parallel run' in text and re.search(r'\bnProcs\s*:\s*2\b',text) and 'FOAM FATAL' not in text,'Genuine completed2rank log');parts=re.split(r'^Time = ('+prep.FLOAT+r')s?\s*$',text,flags=re.M);require(len(parts)==11,'Exactlyfive actualpilotsteps');rows=[];final=[];counts=[];courant=[];components=set()
 number='(?:'+prep.FLOAT+'|[-+]?(?:nan|inf))'
 for i in range(5):
  time=float(parts[i*2+1]);require(abs(time-(i+1)*.05)<=1e-12,'Actualphysicalsteptime');body=parts[i*2+2];outers=re.split(r'PIMPLE:\s+Iteration\s+(\d+)\s*',body);require(len(outers)>=5,'Atleast2actualoutercorrectors');outer_numbers=[int(outers[k]) for k in range(1,len(outers),2)];require(outer_numbers==list(range(1,len(outer_numbers)+1)) and 2<=len(outer_numbers)<=10,'Actualoutercorrectorcount');counts.append(len(outer_numbers));last=[]
  for j in range(1,len(outers),2):
   outer=int(outers[j]);samples=re.findall(r'Solving for (p|Ux|Uy|Uz), Initial residual = ('+number+r'), Final residual = ('+number+r'), No Iterations (\d+)',outers[j+1],re.I);require(all(any(q[0]==name for q in samples) for name in ('p','Ux','Uy')),'Actualcomponentpressure solves required')
   for k,(name,initial,terminal,nit) in enumerate(samples):
    components.add(name);row=dict(step_index=i+1,time_s=time,outer_iteration=outer,solve_index_in_outer=k+1,component=name,initial_residual=float(initial),final_residual=float(terminal),linear_iterations=int(nit));rows.append(row)
    if outer==outer_numbers[-1]:last.append(row)
  final.append(last);co=re.findall(r'Courant Number mean:\s*('+number+r')\s+max:\s*('+number+r')',body,re.I);courant.append(float(co[-1][1]) if co else float('nan'))
 names=[n for n in ('p','Ux','Uy','Uz') if n in components];initial=np.array([[max(q['initial_residual'] for q in data if q['component']==n) for n in names] for data in final]);terminal=np.array([[max(q['final_residual'] for q in data if q['component']==n) for n in names] for data in final]);return rows,names,initial,terminal,np.array(counts),np.array(courant)
def finite_json(v):
 if isinstance(v,dict):return {k:finite_json(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [finite_json(x) for x in v]
 if isinstance(v,(float,np.floating)):return float(v) if math.isfinite(float(v)) else None
 if isinstance(v,(bool,np.bool_)):return bool(v)
 if isinstance(v,(int,np.integer)):return int(v)
 return v

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();case=prep.checked_case(a.run);out=a.out.resolve();require(out.drive.upper()=='E:' and not out.exists() and case not in out.parents,'Fresh external E auditfolder');require(sha(Path(prep.__file__))==PREP_SHA,'Frozenpreparerhash');e=Evidence();checks={};e.touch(Path(__file__));e.touch(Path(prep.__file__));status=e.json(case/'run-status.json');require(status.get('stage')=='finished_pilot_execution' and status.get('actual_steps')==5 and status.get('physical_time_s')==.25 and status.get('production_admitted') is False,'Onlygenuinecompletedfive-step pilot')
 for n,h in e.json(case/'provenance/prepared-input-sha256.json').items():require(sha(e.touch(case/n))==h,'Preparedsourcechanged '+n)
 require(sha(e.touch(case/'provenance/TARGETS.md'))==prep.TARGET_SHA,'Frozen targets changed')
 for n,t in prep.dictionaries().items():
  if n not in ('0/U','0/p'):require(e.text(case/n)==t,'Actualdictionarychanged '+n)
 for folder in ('constant','system'):
  for n in ('MRFProperties','dynamicMeshDict','fvModels','fvConstraints','fvOptions'):
   path=case/folder/n
   if path.exists():require(not re.sub(r'FoamFile\s*\{[^{}]*\}','',prep.clean(e.text(path))).strip(),'Unexpectedmesh/frame/source physics')
 checks['frozen_Euler_upwind_material_PIMPLE_closed_pressure_reference_configuration']=True
 recovery=e.json(case/'provenance/preflight-readonly-recovery.json');require(sha(case/'provenance/preflight-readonly-recovery.json')==RECOVERY_SHA and sha(e.touch(case/'provenance/validate_annular_openfoam_preflight_v001.py'))==VALIDATOR_SHA,'Exact preserved read-only parser recovery source/receipt');require(recovery['readonly_validation_passed'] is True and recovery['native_preflight_commands_completed'] is True and recovery['all_original_input_files_unchanged'] is True and recovery['original_input_file_count']==101 and recovery['CFD_replayed'] is False and recovery['native_commands_replayed'] is False and recovery['solver_started'] is False and recovery['original_preparation_source_sha256']==PREP_SHA and recovery['readonly_validator_source_sha256']==VALIDATOR_SHA and recovery['preflight_receipt_sha256']==PREFLIGHT_SHA and recovery['native_mapping_sha256']==MAPPING_SHA,'Native preflight completed before narrow parser-only continuation');require(recovery['raw_failed_run_status_preserved']['solver_started'] is False,'Original preflight parser failure preserved separately');require(sha(e.touch(case/'provenance/preflight-validation.json'))==PREFLIGHT_SHA and sha(e.touch(case/'provenance/native-mapping.json'))==MAPPING_SHA,'Exact recovered native preflight/mapping');checks['native_preflight_completed_parser_only_recovery_preserved_no_CFD_replay']=True
 pre=e.json(case/'provenance/preflight-validation.json');require(pre['passed'] is True and pre['full_source_cell_point_face_vertex_mapping'] is True,'Actualnativepreflightrequired');poly=prep.native_poly(case)
 for n,h in pre['native_mesh_sha256'].items():require(sha(e.touch(poly/n))==h,'Native mesh changed')
 for n,h in pre['actual_native_C_Vc_sha256'].items():require(sha(e.touch(case/'0'/n))==h,'ActualC/Vcchanged')
 g=prep.native_geometry(poly);pa=prep.patches(poly);C=field(e.text(case/'0/C'),prep.CELLS,3,[0,1,0,0,0,0,0]);V=field(e.text(case/'0/Vc'),prep.CELLS,1,[0,3,0,0,0,0,0]);require(np.all(V>0) and np.max(abs(V-np.asarray(g['volumes'])))<=1e-14,'Positive native orientedVc');mapping=e.json(case/'provenance/native-mapping.json');n2s=np.asarray(mapping['native_cell_to_source_cell']);s2n=np.asarray(mapping['source_native_cell_index']);require(np.array_equal(np.sort(n2s),np.arange(prep.CELLS)) and np.array_equal(n2s[s2n],np.arange(prep.CELLS)),'Source/nativecellbijection');src=e.json(case/'provenance/source-geometry.json');pm=np.asarray(mapping['native_point_to_source_node'])
 for n,vertices in enumerate(g['cell_vertices']):require(set(pm[vertices].tolist())==set(src['quads'][n2s[n]]),'Independent native eightvertex/sourcequad check')
 require(np.max(abs(C[s2n,:2]-np.asarray(src['centres'])))<=1e-10 and np.max(abs(V[s2n]/.01-np.asarray(src['areas'])))<=1e-12,'Native/sourcecentresareas');ranks=Ranks(e,case,g,pa);checks['actual_native_quality_positive_volumes_full_source_and_rank_vertex_mapping']=True
 ref=pre['source_cell0_pressure_reference_unique_rank'];require(n2s[ref['global_native_cell']]==0 and ref['global_native_cell'] in ranks.rows[ref['rank']]['cells'] and np.max(abs(C[ref['global_native_cell']]-np.array([.05092584598549641,.00250182638511251,0])))<=1e-10,'Actualsamefluidreferencepointuniqueownerrank');checks['same_sourcecell0_pressure_reference_unique_actual_decomposed_owner']=True
 registration=e.json(case/'provenance/native-registration-process-status.json')
 for label,r in registration.items():
  text=e.text(case/r['log']);require(r['return_code']==0 and sha(case/r['log'])==r['log_sha256'] and re.search(r'^End\s*$',text,re.M),'Native initialregistrationcompletedlog')
  if label=='parallel':require(re.search(r'\bnProcs\s*:\s*2\b',text),'Actualparallelregistration')
 checks['actual_initial_rotating_BC_fixedFluxPressure_registration_returncode_logs_bound']=True
 rows,names,outer_initial,outer_final,counts,courant=residual_history(e,case);checks['five_actual_Euler_step_logs_at_least2_at_most10_outer_correctors']=True
 checks['every_step_last_outer_initial_p_U_residuals_le1e_6']=bool(np.isfinite(outer_initial).all() and np.all(outer_initial<=1e-6));checks['every_step_last_outer_final_linear_p_U_residuals_le1e_8']=bool(np.isfinite(outer_final).all() and np.all(outer_final<=1e-8))
 Uhist=[];Phist=[];Whist=[];Fhist=[];divhist=[];divmax=[];divrms=[];Lz=[];wallmax=[];wallnet=[];wallabs=[];innererror=[];outererror=[];interface=[];paths=[];native_times=[];finite=[];sf=np.asarray(g['Sf'])[ranks.physical];normal=sf/np.linalg.norm(sf,axis=1)[:,None];cf=np.asarray(g['Cf'])[ranks.physical];target=np.column_stack((-cf[:,1],cf[:,0],np.zeros(128)));target-=normal*np.sum(normal*target,axis=1)[:,None];target[64:]=0;own=np.asarray(g['owner']);nei=np.asarray(g['neighbour'])
 for i in range(6):
  label=f'{i*.05:.10g}';U,P,W,Phi,Ub,Pb,rankfinite,dup=ranks.read(label);statefinite=bool(rankfinite and all(np.isfinite(v).all() for v in (U,P,W,Phi,Ub,Pb)));finite.append(statefinite);interface.append(dup)
  if i:
   for rank in (0,1):native_times.append(dict(rank=rank,step_index=i,**time_meta(e,case/f'processor{rank}'/label,i)))
  net=np.bincount(own,weights=Phi,minlength=prep.CELLS)-np.bincount(nei,weights=Phi[:prep.INTERNAL],minlength=prep.CELLS);div=net/V;divhist.append(div);divmax.append(float(np.max(abs(div))));divrms.append(float(np.sqrt(np.sum(V*div*div)/V.sum())));wall=Phi[ranks.physical];wallmax.append(float(np.max(abs(wall))));wallnet.append(float(wall.sum()));wallabs.append(float(abs(wall).sum()));innererror.append(float(np.max(abs(Ub[:64]-target[:64]))));outererror.append(float(np.max(abs(Ub[64:]))));Lz.append(float(np.sum(V*(C[:,0]*U[:,1]-C[:,1]*U[:,0]))/.01))
  Uhist.append(U);Phist.append(P);Whist.append(W);Fhist.append(Phi);paths.append(dict(index=i,step_index=i,time_s=i*.05,native_case=str(case),time_label=label,source_kind='native_rank_reassembly',initial_state_kind='prescribed_cell_rest_with_actual_native_runtime_BC_flux_readback' if i==0 else None,rank_U=[str(case/f'processor{r}'/label/'U') for r in (0,1)],rank_p=[str(case/f'processor{r}'/label/'p') for r in (0,1)],rank_phi=[str(case/f'processor{r}'/label/'phi') for r in (0,1)],rank_vorticity=[str(case/f'processor{r}'/label/'vorticity') for r in (0,1)],rank_time_metadata=[str(case/f'processor{r}'/label/'uniform/time') for r in (0,1)] if i else [],rank_mesh=ranks.mesh_paths))
 Uhist=np.asarray(Uhist);Phist=np.asarray(Phist);Whist=np.asarray(Whist);Fhist=np.asarray(Fhist);Lz=np.array(Lz);checks['all_saved_native_U_p_phi_vorticity_boundary_values_finite']=bool(all(finite));checks['actual_initial_cell_velocity_rest_le1e_12']=bool(np.max(abs(Uhist[0]))<=1e-12);checks['all_actual_inner_wall_U_matches_projected_Omega_cross_r_le1e_10']=bool(np.isfinite(innererror).all() and max(innererror)<=1e-10);checks['all_actual_outer_wall_noSlip_U_le1e_12']=bool(np.isfinite(outererror).all() and max(outererror)<=1e-12)
 checks['every_saved_wall_phi_max_and_signed_total_abs_le1e_12_m3_s']=bool(np.isfinite(wallmax).all() and max(wallmax)<=1e-12 and max(abs(np.array(wallnet)))<=1e-12);checks['all_native_rank_interface_phi_agrees_le1e_12']=bool(np.isfinite(interface).all() and max(interface)<=1e-12);checks['every_saved_allcell_native_phi_divergence_max_le1e_5_RMS_le1e_6']=bool(np.isfinite(divmax).all() and max(divmax)<=1e-5 and max(divrms)<=1e-6);checks['every_solved_native_Lz_per_depth_positive']=bool(np.isfinite(Lz).all() and np.all(Lz[1:]>0))
 source_rows=[line.split(None,1) for line in e.text(case/'provenance/installed-source-sha256.txt').splitlines() if line.strip()];require(len(source_rows)==len(prep.INSTALLED_SOURCES),'Actualinstalledsourcecount')
 for h,path in source_rows:require(sha(e.touch(case/'provenance/installed-sources'/path[len('/opt/openfoam14/'):]))==h,'Installedsourcehashchanged')
 require('units::radians/dimTime' in e.text(case/'provenance/installed-sources/src/OpenFOAM/primitives/functions/Function1/omega1/omega1.C'),'ActualinstalledomegaSIunitsproof');checks['actual_installed_BC_pressure_PIMPLE_omega_rad_per_s_vorticity_sources_hash_bound']=True
 auth=e.json(case/'provenance/root-pilot-authorization.json');require(auth.get('approved') is True and auth.get('approved_by')=='root' and auth.get('action')=='run_annular_couette_openfoam_v002' and auth.get('allowed_mode')=='PilotRun' and auth['preflight_receipt_sha256']==sha(case/'provenance/preflight-validation.json'),'Separateexactrootpilotadmission')
 for n,key in [('prepare_openfoam_annular_couette_v002.py','preparation_sha256'),('run-openfoam-v002.sh','wrapper_sha256'),('Launch-OpenFOAM-v002.ps1','launcher_sha256')]:require(sha(e.touch(case/'provenance'/n))==auth[key],'Reviewedexecutionsourcechanged '+n)
 require(auth['preparation_sha256']==PREP_SHA and sha(e.touch(case/'provenance/run-openfoam-pilot-executed.sh'))==auth['wrapper_sha256'],'Executedwrapperbinding');checks['separate_five_step_root_admission_execution_source_hashes_bound']=True;source=e.finish();checks['all_consumed_native_inputs_unchanged_during_read_only_audit']=True;checks={n:bool(v) for n,v in checks.items()}
 area=V[s2n]/.01;Psource=Phist[:,s2n];Pmean=np.sum(Psource*area,axis=1)/area.sum();record=dict(phase='five_step_annular_couette_pilot',actual_steps=5,native_states=6,native_case=str(case),native_state_paths=paths,physical_times_s=(np.arange(6)*.05).tolist(),native_time_metadata=native_times,mpi_ranks=2,initial_state_kind='prescribed cell rest; native preflight BC/flux readback, no transient advance',actual_native_patches=pa,actual_native_rank_mapping=ranks.report,pressure_reference=ref,pressure_units='kinematic m2/s2 timesrho1 toPa; comparevolume-mean-subtracted gauges',native_quality=pre,preflight_parser_recovery=recovery,outer_corrector_count=counts.tolist(),residual_component_names=names,last_outer_initial_residual_max=float(np.max(outer_initial)),last_outer_final_linear_residual_max=float(np.max(outer_final)),all_pressure_component_correction_solves_recorded=True,dimensionless_divergence_max=float(max(divmax)),dimensionless_divergence_RMS_max=float(max(divrms)),absolute_wall_phi_max_m3_s=float(max(wallmax)),absolute_signed_boundary_phi_total_max_m3_s=float(max(abs(np.array(wallnet)))),sum_absolute_wall_phi_max_m3_s=float(max(wallabs)),inner_wall_velocity_error_max_m_s=float(max(innererror)),outer_wall_velocity_max_m_s=float(max(outererror)),angular_momentum_per_depth=Lz.tolist(),native_phi_interface_max_difference=float(max(interface)),closed_net_flux_alone_uninformative=True,inlet_relative_mass_metric=None,no_inlet_normalization_reason='Closeddomain; noinlet/no0over0 metric',torque_reported=False,analysis_source_sha256=sha(Path(__file__)),preparation_source_sha256=PREP_SHA,checks=checks,failed_checks=[n for n,v in checks.items() if not v],all_checks_passed=all(checks.values()),physical_validation=False,production_admitted=False,completed_400_step_transient=False)
 out.mkdir(parents=True);np.savez_compressed(out/'fields.npz',time_s=np.arange(6)*.05,step_index=np.arange(6),solved_time_s=np.arange(1,6)*.05,U_source_history_m_s=Uhist[:,s2n,:2],p_source_history_Pa=Psource,p_source_mean_subtracted_history_Pa=Psource-Pmean[:,None],source_centres_m=C[s2n,:2],source_areas_m2=area,source_native_cell_index=s2n,dimensionless_divergence_max=np.array(divmax),dimensionless_divergence_RMS=np.array(divrms),angular_momentum_per_depth=Lz,U_history=Uhist,p_history=Phist,native_phi_history=Fhist,native_vorticity_history=Whist,centres=C,volumes=V,native_cell_to_source_cell=n2s,native_face_to_source_edge=np.asarray(mapping['native_face_to_source_edge']),native_face_source_orientation=np.asarray(mapping['native_face_source_orientation']),native_point_to_source_node=pm,local_divergence_scaled_history=np.asarray(divhist),residual_component_names=np.array(names),last_outer_initial_residuals=outer_initial,last_outer_final_linear_residuals=outer_final,outer_corrector_count=counts,courant_max=courant,courant_available=np.isfinite(courant),inner_wall_velocity_error_m_s=np.array(innererror),outer_wall_velocity_error_m_s=np.array(outererror),wall_native_phi_max_abs_m3_s=np.array(wallmax),signed_boundary_phi_m3_s=np.array(wallnet),rho_kg_m3=np.array(1.),depth_m=np.array(.01))
 with (out/'step-history.csv').open('x',encoding='ascii',newline='') as f:
  w=csv.writer(f);w.writerow(['step_index','time_s','divergence_max_scaled','divergence_RMS_scaled','angular_momentum_per_depth','wall_phi_max_abs_m3_s','signed_boundary_phi_m3_s','inner_U_max_error_m_s','outer_U_max_error_m_s'])
  for i in range(6):w.writerow([i,i*.05,divmax[i],divrms[i],Lz[i],wallmax[i],wallnet[i],innererror[i],outererror[i]])
 with (out/'residual-history.csv').open('x',encoding='ascii',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 record['audited_fields_npz_sha256']=sha(out/'fields.npz');record['step_history_csv_sha256']=sha(out/'step-history.csv');record['residual_history_csv_sha256']=sha(out/'residual-history.csv');source['audited_outputs']={n:dict(sha256=sha(out/n),bytes=(out/n).stat().st_size) for n in ('fields.npz','step-history.csv','residual-history.csv')};prep.json_new(out/'checks.json',finite_json(record));prep.json_new(out/'source-manifest.json',source)
 for p in (Path(__file__),Path(prep.__file__)):prep.copied(p,out/('executed-'+p.name))
 print(json.dumps(dict(output=str(out),all_checks_passed=all(checks.values()),failed_checks=record['failed_checks']),indent=2));return 0 if all(checks.values()) else 2
if __name__=='__main__':raise SystemExit(main())
