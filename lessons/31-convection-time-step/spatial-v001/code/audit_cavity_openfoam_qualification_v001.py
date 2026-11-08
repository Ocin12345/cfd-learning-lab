from pathlib import Path
import json,hashlib,re,math
import numpy as np
from cavity_spatial_common_v001 import *
from cavity_native_fields_v001 import foam,metrics
from cavity_restriction_v001 import compare_fields
from cavity_openfoam_heat_v001 import foam_walls,native_foam_time
from analyse_nozzle_openfoam_v002 import field,FIELDS,patch_values
S=Path(__file__).parent;BASE=ROOT/'native/runs/natural-convection-dt001-v001/openfoam'
def audit(o):
 assert (o/'run-status.txt').read_text().strip()=='complete'
 steps=json.loads((o/'all-step-checks.json').read_text());assert steps['all_passed'] and steps['steps_found']==100
 assert max(abs(row['time_s']-(i+1)*.01) for i,row in enumerate(steps['steps']))<1e-9
 log=(o/'solver.log').read_text();co=[];current=0.
 for line in log.splitlines():
  m=re.match(r'Time = ([0-9.eE+-]+)s?$',line)
  if m:current=float(m[1])
  m=re.match(r'Courant Number mean: [^ ]+ max: ([-+0-9.eE]+)',line)
  if m:co.append(dict(actual_pre_advance_time_s=current,maximum_native_courant=float(m[1])))
 assert len(co) in (100,101) and max(v['maximum_native_courant'] for v in co)<=.5
 initial=foam(o,0,96);coarse=foam(BASE,0,48);initial_comparison=compare_fields(coarse,initial,48,relative=False)
 m0=metrics(initial,96)['mass_per_depth_kg_m'];initial_comparison['relative_total_mass_difference']=abs(m0/metrics(coarse,48)['mass_per_depth_kg_m']-1)
 assert initial_comparison['max_temperature_difference_K']<=2e-4 and initial_comparison['max_velocity_difference_m_s']<=2e-5 and initial_comparison['relative_total_mass_difference']<=2e-5
 data=[];scale=RHO0*L*math.sqrt(G*(HOT-COLD)*L/T0)
 for step in range(1,101):
  t=step*.01;folder=o/f'{t:g}';clock=native_foam_time(folder/'uniform/time',t);assert clock['step']==step
  U=field(folder/'U',9216,*FIELDS['U']);assert abs(U[:,2]).max()<1e-12
  q=foam(o,t,96);m=metrics(q,96);m['relative_mass_drift']=abs(m['mass_per_depth_kg_m']/m0-1)
  assert COLD-1e-4<=m['Tmin_K']<=m['Tmax_K']<=HOT+1e-4 and m['eos_relative_error']<=1e-6 and m['relative_mass_drift']<=1e-5
  flux={name:float(abs(patch_values(folder/'phi',name,96)).sum())/DEPTH for name in WALLS};leak=sum(flux.values());assert leak<=1e-12 and leak/scale<=1e-10
  data.append(dict(step=step,time_s=t,clock=clock,metrics=m,total_absolute_wall_mass_flux_kg_ms=leak))
 for part in ['processor0','processor1']:native_foam_time(o/part/'1/uniform/time',1.)
 walls,_=foam_walls(o,WALLS,1.);assert all(w['face_count']==96 and abs(w['wall_length_m']-L)<1e-10 for w in walls.values())
 tolerance=1e-10+1e-7*sum(abs(w['report_heat_per_depth_W_m']) for w in walls.values());assert all(abs(w['heat_per_depth_W_m']-w['report_heat_per_depth_W_m'])<=tolerance for w in walls.values())
 for folder in ['constant','system']:
  for p in (BASE/folder).iterdir():
   if p.is_file() and p.name!='controlDict':assert (o/folder/p.name).read_bytes()==p.read_bytes(),p
 expected=(BASE/'system/controlDict').read_text();expected=re.sub(r'\bendTime\s+20;', 'endTime 1;',expected);expected=re.sub(r'\bwriteInterval\s+50;', 'writeInterval 1;',expected);assert (o/'system/controlDict').read_text()==expected
 result=dict(numerical_pass=True,heated_steps=100,actual_native_saved_states=101,selected_matched_times_s=[0,.01,.1,.25,.5,.75,1],cells=9216,fixed_dt_s=.01,initial_comparison=initial_comparison,maximum_native_courant=max(v['maximum_native_courant'] for v in co),native_courant_records=co,native_state_checks=data,native_wall_heat=walls,Nu_hot_at_1s=walls['left']['heat_per_depth_W_m']/(K*(HOT-COLD)),identical_physical_numerical_dictionaries_to_baseline_except_end_write=True,energy_scope='Native energy residuals/EOS/gravity formulation and independently integrated instantaneous wall heat, no full transient energy closure',physical_validation=False,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in o.rglob('*') if p.is_file()})
 p=o/'qualification-audit-v001.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');return result