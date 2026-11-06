"""Read actual nozzle OpenFOAM fields and native mass flux; never start solvers.

LTS indices are pseudo-iterations, not physical time. Diagonal equation residuals
are retained as raw evidence and deliberately excluded from convergence gates.
"""
from pathlib import Path
import argparse, csv, itertools, json, math, re, shutil
import numpy as np
from prepare_openfoam_nozzle_v001 import ROOT,RUN,F,require,sha,write,jwrite,copy,poly,patches,mesh_quality,clean
FIELDS={'U':(3,[0,1,-1,0,0,0,0]),'p':(1,[1,-1,-2,0,0,0,0]),'T':(1,[0,0,0,1,0,0,0]),'rho':(1,[1,-3,0,0,0,0,0]),'C':(3,[0,1,0,0,0,0,0]),'phi':(1,[1,0,-1,0,0,0,0])}

def field(path,count,width,dims):
 s=clean(path.read_text(encoding='ascii'));h=re.search(r'FoamFile\s*\{([^{}]*)\}',s)
 require(h is not None and re.search(r'\bformat\s+ascii\s*;',h[1]),'Require native ASCII field '+str(path))
 d=re.search(r'\bdimensions\s*\[([^]]+)\]\s*;',s);require(d is not None and list(map(float,d[1].split()))==dims,'Wrong field dimensions '+str(path))
 un=re.search(r'\binternalField\s+uniform\s+([^;]+);',s)
 if un:
  v=np.fromstring(un[1].strip(' ()\n\t'),sep=' ');require(len(v)==width,'Uniform width differs '+str(path));a=np.tile(v,(count,1))
 else:
  m=re.search(r'\binternalField\s+nonuniform\s+List<(?:scalar|vector)>\s+(\d+)\s*\((.*?)\)\s*;',s,re.S)
  require(m is not None and int(m[1])==count,'Nonuniform count differs '+str(path))
  a=np.fromstring(m[2].replace('(',' ').replace(')',' '),sep=' ')
  require(len(a)==count*width,'Nativefield length differs '+str(path));a=a.reshape(count,width)
 require(np.isfinite(a).all(),'Nonfinite field '+str(path))
 return a[:,0] if width==1 else a

def centre_map(native,reference):
 require(native.shape==(9600,3) and reference.shape==(9600,2),'Centre shape differs')
 require(np.max(np.abs(native[:,2]))<1e-10,'Nonplanar native cell centres')
 step=1e-8; keys=np.floor(reference/step).astype(np.int64);lookup={}
 for i,k in enumerate(keys):
  key=tuple(k);require(key not in lookup,'Duplicate reference coordinatebins');lookup[key]=i
 order=np.full(9600,-1,dtype=np.int64);worst=0.
 for native_id,point in enumerate(native[:,:2]):
  key=np.floor(point/step).astype(np.int64);hits=[]
  for dx,dy in itertools.product((-1,0,1),repeat=2):
   cid=lookup.get((key[0]+dx,key[1]+dy))
   if cid is not None and np.linalg.norm(reference[cid]-point)<1e-9:hits.append(cid)
  require(len(hits)==1 and order[hits[0]]<0,'Nativecentres do notbiject exactsource polygoncentres')
  order[hits[0]]=native_id;worst=max(worst,float(np.linalg.norm(reference[hits[0]]-point)))
 require(np.all(order>=0),'Incomplete nativecentre mapping')
 return order,dict(bijective=True,maximum_coordinate_error_m=worst,tolerance_m=1e-9,method='Unique nativeC to true source polygoncentroid bijection, independent of cell order')

def patch_values(path,patch,count):
 s=clean(path.read_text(encoding='ascii'));require('boundaryField' in s,'Missing native boundary field')
 blocks={n.strip('"'):b for n,b in re.findall(r'([^\s{}]+)\s*\{([^{}]*)\}',s.split('boundaryField',1)[1])}
 require(patch in blocks,'Missing phi patch '+patch);b=blocks[patch]
 un=re.search(r'\bvalue\s+uniform\s+([^;]+);',b)
 if un:a=np.full(count,float(un[1]))
 else:
  m=re.search(r'\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',b,re.S)
  require(m is not None and int(m[1])==count,'Missing/incorrect boundary phi value '+patch);a=np.fromstring(m[2],sep=' ')
 require(len(a)==count and np.isfinite(a).all(),'Nonfinite/incomplete patch phi');return a

def native_flux(case,mesh,iteration):
 rows=patches(mesh);p=case/str(iteration)/'phi'; ninternal=min(r['start'] for r in rows.values())
 field(p,ninternal,*FIELDS['phi'])
 sums={name:float(patch_values(p,name,row['count']).sum()) for name,row in rows.items() if row['kind']!='empty'}
 require(sums['inlet']<0<sums['outlet'],'Expected actual inflow/outflow')
 return dict(signed_mass_flux_kg_s=sums,signed_mass_flux_kg_s_per_m_depth={k:v/.01 for k,v in sums.items()},relative_imbalance=abs(sum(sums.values()))/abs(sums['inlet']),definition='abs(sum signed outward native phi over4 physicalpatches)/abs(inlet sum); emptyplanes excluded')

def flux_history(case):
 sources=[];hist={}
 for patch in ('inlet','outlet','upper_wall','lower_wall'):
  files=sorted((case/'postProcessing'/('flux_'+patch)).glob('*/surfaceFieldValue.dat'))
  require(files,'Missing native phi history '+patch);hist[patch]={}
  for path in files:
   sources.append(path)
   for line in path.read_text().splitlines():
    if not line.strip() or line.lstrip().startswith('#'):continue
    row=[float(v) for v in line.split()];require(len(row)==2 and all(map(math.isfinite,row)),'Unexpected native phi history row')
    t=int(round(row[0]));require(abs(row[0]-t)<1e-8 and t not in hist[patch],'Duplicate/noninteger pseudoiteration')
    hist[patch][t]=row[1]
  require(all(t in hist[patch] for t in range(1,6001)),'Incomplete1..6000masshistory')
 rows=[]
 for t in range(1,6001):
  values=[hist[p][t] for p in ('inlet','outlet','upper_wall','lower_wall')];den=abs(values[0]);relative=abs(sum(values))/den if den>1e-30 else float('inf')
  rows.append([t,*values,relative])
 require(all(math.isfinite(v) for row in rows for v in row),'Nonfinite massbalance history')
 last=rows[-100:]
 return rows,sources,dict(samples=6000,last100_maximum_relative_imbalance=max(r[-1] for r in last),final_relative_imbalance=rows[-1][-1],mass_flux_units='kg/s through0.01m depth',history_is_pseudo_iteration_not_time=True)

def shock(fields,centres):
 mach=fields['Mach'].reshape(40,240);p=fields['p'].reshape(40,240);x=centres[:,0].reshape(40,240)
 m=mach[19:21].mean(axis=0);pr=p[19:21].mean(axis=0);xx=x[19:21].mean(axis=0)
 candidates=np.flatnonzero((xx[:-1]>.4)&(m[:-1]>=1)&(m[1:]<1)&(pr[1:]>pr[:-1]))
 if not len(candidates):return dict(detected=False,description='No postthroat centerline supersonic-to-subsonic crossing with local pressure rise detected')
 i=int(candidates[np.argmax((np.diff(pr)/np.diff(xx))[candidates])]);fraction=float((m[i]-1)/(m[i]-m[i+1]));location=float(xx[i]+fraction*(xx[i+1]-xx[i]))
 return dict(detected=True,x_m=location,bracket_x_m=[float(xx[i]),float(xx[i+1])],bracket_Mach=[float(m[i]),float(m[i+1])],bracket_static_p_Pa=[float(pr[i]),float(pr[i+1])],description='Centerline compression crossing consistent with a numerically smeared shock; two central cellrows averaged; no exact normalshock jump validation')

def rel_l2(a,b,w):
 if a.ndim==2:w=w[:,None]
 return float(np.sqrt(np.sum(w*(a-b)**2)/max(np.sum(w*b*b),1e-30)))

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
 case=args.case.resolve();out=args.out.resolve();run=case.parent;meshdir=run/'mesh'
 require(case==(ROOT/'native/runs'/RUN/'openfoam').resolve(),'Wrong native case');require(not out.exists() and out.is_relative_to(ROOT/'local/nozzle-v001'),'Fresh private E:analysis folder required')
 initial=poly(case);target_path=case/'declared-comparison.json';status_path=case/'run-status.json';setup_path=meshdir/'setup.json';npz_path=meshdir/'geometry.npz'
 target=json.loads(target_path.read_text());setup=json.loads(setup_path.read_text());status=json.loads(status_path.read_text())
 require(status.get('stage')=='finished_execution' and status.get('final_pseudo_iteration')==6000,'Require completed fixedcap6000nativeexecution')
 require(target['max_pseudo_iterations']==6000 and target['iteration_indices_are_not_physical_time'] is True,'Wrong declared cap/time interpretation')
 require(sha(meshdir/'nozzle.msh')==sha(case/'input-nozzle.msh')==setup['mesh_sha256']==target['mesh_source_sha256'],'Sharedmesh hashes differ')
 require(sha(npz_path)==setup['geometry_npz_sha256'],'SourcegeometryNPZ hash differs')
 require(sha(Path(__file__).with_name('prepare_openfoam_nozzle_v001.py'))==sha(case/'provenance/prepare_openfoam_nozzle_v001.py'),'Preparation source differs fromfreeze')
 inputs=[Path(__file__),Path(__file__).with_name('prepare_openfoam_nozzle_v001.py'),target_path,status_path,setup_path,npz_path,meshdir/'nozzle.msh',case/'input-nozzle.msh',case/'solver.log',case/'mesh-check.log',case/'reconstruction.log',case/'cell-centres.log',case/'6000/C',*(initial/n for n in ('points','faces','owner','neighbour','boundary')),*sorted((case/'system').glob('*')),*sorted((case/'constant').glob('*')),*sorted((case/'provenance').rglob('*'))]
 inputs=[p for p in inputs if p.is_file()]
 for t in (5000,5500,6000):inputs.extend(case/str(t)/n for n in ('U','p','T','rho','phi'))
 hrows,hsources,hsummary=flux_history(case);inputs.extend(hsources);inputs=list(dict.fromkeys(inputs));before={str(p):sha(p) for p in inputs}
 prepared=json.loads((case/'provenance/prepared-input-sha256.json').read_text());require(all(sha(case/n)==v for n,v in prepared.items()),'Frozen preparedinput differs')
 preflight=json.loads((case/'provenance/preflight-validation.json').read_text());require(preflight['passed'] and all(sha(initial/n)==v for n,v in preflight['native_mesh_sha256'].items()),'Native mesh changed aftercheckMesh')
 q=mesh_quality((case/'mesh-check.log').read_text());nativeC=field(case/'6000/C',9600,*FIELDS['C'])
 with np.load(npz_path) as data:centres=data['centres'].copy();areas=data['areas'].copy()
 order,mapping=centre_map(nativeC,centres); require(np.all(areas>0),'Invalid sourceareas')
 states={};fluxes={};shocks={};rhoerrors={}
 for t in (5000,5500,6000):
  fields={name:field(case/str(t)/name,9600,*FIELDS[name])[order] for name in ('U','p','T','rho')}
  require(min(float(fields[n].min()) for n in ('p','T','rho'))>0,'Nonpositive actualgasproperty')
  require(np.max(np.abs(fields['U'][:,2]))<1e-9,'Nonplanar savedvelocity')
  fields['Mach']=np.linalg.norm(fields['U'],axis=1)/np.sqrt(1.4*287.05*fields['T'])
  rhoerrors[t]=float(np.max(np.abs(fields['rho']-fields['p']/(287.05*fields['T']))/fields['rho']))
  states[t]=fields;fluxes[t]=native_flux(case,initial,t);shocks[t]=shock(fields,centres)
 stability={}
 for first,last in ((5000,5500),(5500,6000),(5000,6000)):
  stability[str(first)+'_to_'+str(last)]={n:rel_l2(states[first][n],states[last][n],areas) for n in ('U','p','T')}
 allshock=all(v['detected'] for v in shocks.values());shockspan=max(v['x_m'] for v in shocks.values())-min(v['x_m'] for v in shocks.values()) if allshock else None
 log=(case/'solver.log').read_text();times=[int(float(t)) for t in re.findall(r'^Time = ('+F+r')s?\s*$',log,re.M)]
 require(times==list(range(1,6001)) and re.search(r'^End\s*$',log,re.M),'Require exact6000actualnativepseudoiterations andEnd')
 require('Using LTS' in log and 'localEuler' in (case/'system/fvSchemes').read_text(),'NativeLTS evidence missing')
 require(not re.search(r'FOAM FATAL|Floating point exception(?! trapping)|\bnan\b|\binf\b',log,re.I),'Native solver failure/nonfinite log')
 physical=(case/'constant/physicalProperties').read_text();control=(case/'system/controlDict').read_text();schemes=(case/'system/fvSchemes').read_text();solution=(case/'system/fvSolution').read_text()
 expected_settings=all(x in physical for x in ('perfectGas','sensibleInternalEnergy','Cp 1004.675','mu 0')) and all(x in control for x in ('solver shockFluid;','endTime 6000;','writeInterval 500;','adjustTimeStep false;')) and 'maxCo 0.4;' in solution and 'default localEuler;' in schemes and 'fluxScheme Kurganov;' in schemes
 checks=dict(actual_fixedcap_saved_execution_present=True,native_mesh_check_passed=True,shared_geometry_and_prepared_inputs_unchanged=True,actual_cellcentre_bijection_passed=True,expected_inviscid_gas_LTS_settings_match=expected_settings,finite_positive_gas_fields=True,gas_EOS_max_relative_error_le_1e_5=max(rhoerrors.values())<=1e-5,final_native_relative_mass_imbalance_le_1e_4=fluxes[6000]['relative_imbalance']<=1e-4,last100_phi_history_max_relative_imbalance_le_1e_4=hsummary['last100_maximum_relative_imbalance']<=1e-4,last_three_saved_states_relative_L2_U_p_T_le_1e_3=max(v for row in stability.values() for v in row.values())<=1e-3,compression_crossing_location_stable_if_detected=bool(not any(v['detected'] for v in shocks.values()) or (allshock and shockspan<=.01)))
 require(all(sha(p)==before[str(p)] for p in inputs),'Sourcechanged duringreadonlyanalysis')
 stats={n:dict(minimum=float(a.min()),maximum=float(a.max()),area_weighted_mean=float(np.average(a,weights=areas))) for n,a in states[6000].items() if a.ndim==1}
 record=dict(solver='OpenFOAM Foundation14 shockFluid',analysis_completed=True,mesh_cells=9600,final_pseudo_iteration=6000,physical_time=False,diagonal_solver_residuals_are_not_convergence_proof=True,mesh_quality=q,centre_mapping=mapping,boundary_mass_balance_by_pseudo_iteration=fluxes,mass_flux_history=hsummary,final_field_statistics=stats,gas_EOS_relative_error_by_iteration=rhoerrors,last_three_saved_state_relative_L2_change=stability,field_change_definition='sqrt(sum(area*(old-new)^2)/sum(area*new^2)); vectorcomponents summed; savedendpoint comparisons not everyiteration stability',centerline_compression_crossing_by_iteration=shocks,centerline_compression_location_range_m=shockspan,compression_crossing_detected_all_three_saved_states=allshock,shock_location_check_applicable=any(v['detected'] for v in shocks.values()),shock_detection_is_observation_not_convergence_requirement=True,check_results=checks,failed_checks=[k for k,v in checks.items() if not v],all_declared_checks_passed=all(checks.values()),physical_validation=False,benchmark_validation=False,mesh_independence_established=False,limitation='Bounded steadyseeking Euler calculation onone modestmesh. LTSindices andtracers are not physical transient CFD. A centerline compression crossing does not validate exact shock position/jump.')
 out.mkdir(parents=True)
 jwrite(out/'openfoam-checks.json',record)
 with (out/'openfoam-mass-history.csv').open('x',encoding='ascii',newline='') as f:
  writer=csv.writer(f);writer.writerow(['pseudo_iteration','inlet_kg_s','outlet_kg_s','upper_wall_kg_s','lower_wall_kg_s','relative_imbalance']);writer.writerows(hrows)
 with (out/'openfoam-centerline-final.csv').open('x',encoding='ascii',newline='') as f:
  writer=csv.writer(f);writer.writerow(['x_m','static_p_Pa','temperature_K','Mach'])
  for i in range(240):writer.writerow([float(centres[:,0].reshape(40,240)[19:21,i].mean()),*(float(states[6000][n].reshape(40,240)[19:21,i].mean()) for n in ('p','T','Mach'))])
 copy(Path(__file__),out/'executed-analyser.py');copy(Path(__file__).with_name('prepare_openfoam_nozzle_v001.py'),out/'executed-preparation-helper.py')
 outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()}
 jwrite(out/'source-manifest.json',dict(inputs_unchanged_during_analysis=True,sources=[dict(path=p,sha256=v) for p,v in before.items()],outputs=outputs))
 print(json.dumps(record,indent=2,allow_nan=False));return 0 if all(checks.values()) else 2
if __name__=='__main__':raise SystemExit(main())
