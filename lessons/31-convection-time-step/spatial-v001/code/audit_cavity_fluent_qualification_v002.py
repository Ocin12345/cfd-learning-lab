from pathlib import Path
import json,hashlib,math
import numpy as np
from cavity_spatial_common_v001 import *
from cavity_native_fields_v001 import fluent,metrics
from cavity_restriction_v001 import compare_fields
from cavity_native_heat_v001 import fluent_wall_flux,native_fluent_time
from cavity_step_audit_v001 import audit_steps
S=Path(__file__).parent;BASE=ROOT/'native/runs/natural-convection-dt001-v001/fluent-v002'

def check_settings(actual):
 base=json.loads((BASE/'heated-settings.json').read_text());roundoff=[]
 def same(a,b,path):
  if isinstance(a,dict) and isinstance(b,dict):
   assert a.keys()==b.keys(),(path,a.keys(),b.keys())
   for k in a:same(a[k],b[k],path+'/'+k)
  elif isinstance(a,list) and isinstance(b,list):
   assert len(a)==len(b),path
   for i,(x,y) in enumerate(zip(a,b)):same(x,y,path+'/'+str(i))
  elif isinstance(a,(int,float)) and not isinstance(a,bool) and isinstance(b,(int,float)) and not isinstance(b,bool):
   assert math.isclose(a,b,rel_tol=1e-14,abs_tol=1e-30),(path,a,b)
   if a!=b:roundoff.append(dict(path=path,actual=a,baseline=b))
  else:assert a==b,(path,a,b)
 for key in ['general','energy','viscous','material','boundaries','references','methods','controls']:same(actual[key],base[key],key)
 assert actual['methods']['p_v_coupling']['flow_scheme']=='PISO'
 return dict(passed=True,numeric_relative_tolerance=1e-14,roundoff_only=roundoff,baseline=str(BASE/'heated-settings.json'),compared_keys=['general','energy','viscous','material','boundaries','references','methods','controls'])

def check_initial(o):
 coarse=fluent(BASE/'heated-start.cas.h5',BASE/'heated-0000.dat.h5',48);fine=fluent(o/'heated-start.cas.h5',o/'heated-start.dat.h5',96)
 a=compare_fields(coarse,fine,48,relative=False);a['relative_total_mass_difference']=float(abs(metrics(fine,96)['mass_per_depth_kg_m']/metrics(coarse,48)['mass_per_depth_kg_m']-1))
 assert a['max_temperature_difference_K']<=2e-4 and a['max_velocity_difference_m_s']<=2e-5 and a['relative_total_mass_difference']<=2e-5,a
 return a

def audit(o):
 assert json.loads((o/'status.json').read_text())['stage']=='closed'
 rows=audit_steps((o/'solver-transcript.txt').read_text(errors='replace'),100,.01);assert all(r['passed'] for r in rows)
 actual=json.loads((o/'heated-settings.json').read_text());settings=check_settings(actual);initial=check_initial(o)
 assert json.loads((o/'heated-gravity-work-source.json').read_text())['readback']==json.loads((BASE/'heated-gravity-work-source.json').read_text())['readback']
 assert json.loads((o/'mesh-quality.json').read_text())['minimum_orthogonal_quality']>=.2
 frames=json.loads((o/'frames.json').read_text());assert [r['step'] for r in frames]==[0,1,10,25,50,75,100]
 data=[];m0=None;scale=RHO0*L*math.sqrt(G*(HOT-COLD)*L/T0)
 for row in frames:
  p=o/row['file'];clock=native_fluent_time(p);assert abs(clock['time_s']-row['step']*.01)<1e-9 and clock['step']==row['step'] and abs(clock['dt_s']-.01)<1e-12
  q=fluent(o/'qualification-final.cas.h5',p,96);m=metrics(q,96)
  if m0 is None:m0=m['mass_per_depth_kg_m']
  m['relative_mass_drift']=abs(m['mass_per_depth_kg_m']/m0-1);assert COLD-1e-4<=m['Tmin_K']<=m['Tmax_K']<=HOT+1e-4 and m['relative_mass_drift']<=1e-5 and m['eos_relative_error']<=1e-6
  report=(o/f"wall-mass-{row['step']:04d}.txt").read_text();import re
  mass=[float(re.findall(r'^\s*'+name+r'\s+([-+0-9.eE]+)\s*$',report,re.M)[0]) for name in WALLS];assert sum(abs(v) for v in mass)<=1e-12 and sum(abs(v) for v in mass)/scale<=1e-10
  data.append(dict(step=row['step'],clock=clock,metrics=m,native_wall_mass_kg_ms=dict(zip(WALLS,mass))))
 walls,paths=fluent_wall_flux(o/'qualification-final.cas.h5',o/'qualification-final.dat.h5',o/'wall-heat-final.txt',WALLS)
 assert set(walls)==set(WALLS) and all(v['face_count']==96 and abs(v['wall_length_m']-L)<1e-10 for v in walls.values())
 assert sum(w['absolute_mass_flux_per_depth_kg_ms'] for w in walls.values())<=1e-12
 result=dict(numerical_pass=True,time_steps=100,saved_states=7,cells=9216,dt_s=.01,settings=settings,initial_fields=initial,saved_native_states=data,native_wall_heat=walls,Nu_hot_at_1s=walls['left']['heat_per_depth_W_m']/(K*(HOT-COLD)),energy_scope='Native energy residuals, thermodynamic/source checks and independently integrated instantaneous wall heat; no full transient global energy closure',physical_validation=False,source_binding_scope='Explicit persistent native fields, settings, reports and executed source; asynchronous cleanup script is not a preserved simulation input',source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in o.iterdir() if p.is_file() and (p.suffix=='.h5' or p.name in ['executed-driver.py','launch-admission.json','declared-settings.json','gravity-work-source.json','control-settings.json','control-checks.json','heated-gravity-work-source.json','heated-settings.json','matched-baseline-settings.json','initial-comparison.json','frames.json','time-step-checks.json','final-settings.json','status.json','mesh-quality.json','solver-transcript.txt'] or p.name.startswith(('wall-mass-','wall-heat-'))) })
 p=o/'qualification-audit-v002.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');return result