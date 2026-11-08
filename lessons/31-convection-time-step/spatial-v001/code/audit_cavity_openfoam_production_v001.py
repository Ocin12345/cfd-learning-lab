from pathlib import Path
import hashlib,json,re,math
import numpy as np
from cavity_spatial_common_v001 import *
from cavity_native_fields_v001 import foam,metrics
from cavity_openfoam_heat_v001 import native_foam_time,foam_walls
from analyse_nozzle_openfoam_v002 import field,FIELDS,patch_values
S=Path(__file__).parent;N=ROOT/'native/runs/showcase-convection-spatial-v001';P=N/'openfoam-qualification-v001'
def audit(o):
 assert (o/'run-status.txt').read_text().strip()=='complete';steps=json.loads((o/'all-step-checks.json').read_text());assert steps['all_passed'] and steps['steps_found']==1900
 assert max(abs(row['time_s']-(i+101)*.01) for i,row in enumerate(steps['steps']))<1e-8
 prefix=json.loads((P/'all-step-checks.json').read_text());assert prefix['all_passed'] and prefix['steps_found']==100
 merged=prefix['steps']+steps['steps'];assert len(merged)==2000 and max(abs(r['time_s']-(i+1)*.01) for i,r in enumerate(merged))<1e-8
 log=(o/'solver.log').read_text();co=[];current=1.
 for line in log.splitlines():
  m=re.match(r'Time = ([0-9.eE+-]+)s?$',line)
  if m:current=float(m[1])
  m=re.match(r'Courant Number mean: [^ ]+ max: ([-+0-9.eE]+)',line)
  if m:co.append(dict(actual_pre_advance_time_s=current,maximum_native_courant=float(m[1])))
 assert len(co) in (1900,1901) and max(v['maximum_native_courant'] for v in co)<=.5
 reference=foam(P,1.,96);restart=foam(o,1.,96);difference={k:float(abs(reference[k]-restart[k]).max()) for k in reference};assert all(v==0 for v in difference.values())
 m0=metrics(foam(P,0.,96),96)['mass_per_depth_kg_m'];rows=[];scale=RHO0*L*math.sqrt(G*(HOT-COLD)*L/T0)
 for step in [100]+list(range(150,2001,50)):
  t=step*.01;folder=o/f'{t:g}';clock=native_foam_time(folder/'uniform/time',t);assert clock['step']==step
  q=foam(o,t,96);m=metrics(q,96);assert all(np.isfinite(v).all() for v in q.values()) and all((q[k]>0).all() for k in ['p','rho','T']);m['relative_mass_drift']=abs(m['mass_per_depth_kg_m']/m0-1)
  assert COLD-1e-4<=m['Tmin_K']<=m['Tmax_K']<=HOT+1e-4 and m['eos_relative_error']<=1e-6 and m['relative_mass_drift']<=1e-5
  U=field(folder/'U',9216,*FIELDS['U']);assert abs(U[:,2]).max()<1e-12
  leaks={name:float(abs(patch_values(folder/'phi',name,96)).sum())/DEPTH for name in WALLS};leak=sum(leaks.values());assert leak<=1e-12 and leak/scale<=1e-10
  rows.append(dict(step=step,time_s=t,clock=clock,metrics=m,total_absolute_wall_mass_flux_kg_ms=leak))
 for part in ['processor0','processor1']:native_foam_time(o/part/'20/uniform/time',20.)
 walls,_=foam_walls(o,WALLS,20.);assert all(v['face_count']==96 and abs(v['wall_length_m']-L)<1e-10 for v in walls.values())
 tolerance=1e-10+1e-7*sum(abs(v['report_heat_per_depth_W_m']) for v in walls.values());assert all(abs(v['heat_per_depth_W_m']-v['report_heat_per_depth_W_m'])<=tolerance for v in walls.values())
 expected=(P/'system/controlDict').read_text();expected=re.sub(r'\bstartTime\s+0;', 'startTime 1;',expected);expected=re.sub(r'\bendTime\s+1;', 'endTime 20;',expected);expected=re.sub(r'\bwriteInterval\s+1;', 'writeInterval 50;',expected);assert (o/'system/controlDict').read_text()==expected
 for d in ['constant','system']:
  for p in (P/d).iterdir():
   if p.is_file() and p.name!='controlDict':assert p.read_bytes()==(o/d/p.name).read_bytes(),p
 result=dict(numerical_pass=True,new_steps=1900,merged_steps=2000,continued_native_saved_states=len(rows),selected_common_states=41,native_checkpoint_source=str(P),restart_field_max_absolute_differences=difference,fixed_dt_s=.01,maximum_native_courant=max(r['maximum_native_courant'] for r in co),native_courant_records=co,saved_native_checks=rows,native_wall_heat_at20s=walls,Nu_hot_at20s=walls['left']['heat_per_depth_W_m']/(K*(HOT-COLD)),no_independent_full_rerun=True,physical_validation=False,energy_scope='Native residual/EOS/gravity formulation and instantaneous face heat; no full transient storage closure',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in o.rglob('*') if p.is_file()})
 p=o/'production-audit-v001.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');(S/'openfoam-combined-2000-step-checks-v001.json').write_text(json.dumps(merged,indent=2)+'\n');return result