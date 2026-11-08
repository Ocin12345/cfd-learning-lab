from pathlib import Path
import json,re,hashlib
import numpy as np,h5py
from ansys.fluent.core.offline import _lispy
from analyse_nozzle_fluent_v005 import common_mapping
from sloshing_dt_common_v001 import *
from sloshing_courant_v001 import courant_records
from sloshing_fluent_steps_v001 import parse_steps
from audit_sloshing_dt_qualification_v003 import field,saved_check
from sloshing_phi_v001 import phi_wall_flux
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fclock(p):
 with h5py.File(p) as f:d={x[0]:x[1] for x in _lispy.parse(f['settings/Data Variables'][0].decode())[1]}
 return dict(step=int(d['time-step']),time_s=float(d['flow-time']),dt_s=float(d['physical-time-step']))
def ffields(case,data,geometry):
 with np.load(geometry) as g:m=common_mapping(case,g)
 with h5py.File(data) as f:
  q={}
  for k,name in [('u','SV_U'),('v','SV_V'),('p','SV_P')]:
   x=np.empty(NX*NY);x[m]=cell_values(f['results/1/phase-1/cells/'+name]);q[k]=x
  x=np.empty(NX*NY);x[m]=cell_values(alpha_group(f));q['a']=x
 return q
def guardF(case,data,geometry,step):
 c=fclock(data);assert c['step']==step and abs(c['time_s']-step*DT)<1e-10 and abs(c['dt_s']-DT)<1e-14
 q=ffields(case,data,geometry);r=saved_check(step,c['time_s'],q['a'],np.column_stack([q['u'],q['v']]),q['a'],'continuation',q['p']);assert all(r[k] for k in ['alpha_bounds_pass','volume_pass','quiet_pass']);return r

def auditF(O):
 assert not (O/'production-audit-v001.json').exists();P=N/'fluent-qualification-v001';prefix=json.loads((P/'qualification-audit-v001.json').read_text());assert prefix['numerical_pass']
 rows=parse_steps((O/'solver-transcript.txt').read_text(errors='replace'),80,4800,DT);merged=prefix['step_records']+rows;assert [r['step'] for r in merged]==list(range(1,4801)) and all(r['passed'] for r in merged)
 frames=json.loads((O/'frames.json').read_text());assert [r['step'] for r in frames]==[80,81]+list(range(120,4801,40))
 records=[guardF(O/'active.cas.h5',O/f'state-{r["step"]:04d}.dat.h5',O/'geometry.npz',r['step']) for r in frames]
 q=ffields(P/'tank-final.cas.h5',P/'tank-final.dat.h5',P/'geometry.npz');z=ffields(O/'active.cas.h5',O/'state-0080.dat.h5',O/'geometry.npz');assert all(np.array_equal(q[k],z[k]) for k in q)
 co=courant_records((O/'solver-transcript.txt').read_text(errors='replace'),'fluent',DT,80,4800);assert max(r['maximum_native_interface_courant'] for r in co)<=.25
 sources=[p for p in O.iterdir() if p.is_file() and (p.suffix in ['.h5','.npz','.csv'] or p.name in ['solver-transcript.txt','executed-driver.py','frames.json','loaded-fixed-settings.json','fixed-settings-comparison.json','declared-settings.json','all-step-checks.json'])]
 assert json.loads((O/'fixed-settings-comparison.json').read_text())['passed'];a=dict(numerical_pass=True,new_steps=4720,merged_steps=4800,continued_states=len(records),selected_matched_states=121,records=records,step_records=rows,courant_records=co,exact_checkpoint_restart_fields=True,physical_validation=False,independent_full_rerun=False,source_sha256={str(p):sha(p) for p in sources})
 (O/'production-audit-v001.json').write_text(json.dumps(a,indent=2)+'\n');(S/'fluent-combined-4800-step-checks-v001.json').write_text(json.dumps(merged,indent=2)+'\n');return a

def auditO(O):
 assert not (O/'production-audit-v001.json').exists();P=N/'openfoam-qualification-v001';prefix=json.loads((P/'qualification-audit-v001.json').read_text());assert prefix['numerical_pass'];text=(O/'solver.log').read_text();assert re.search(r'^End\s*$',text,re.M) and 'FOAM FATAL' not in text
 parts=re.split(r'^Time = ([0-9.eE+-]+)s?\s*$',text,flags=re.M);rows=[]
 for i in range(1,len(parts),2):
  t=float(parts[i]);chunk=parts[i+1];assert 'PIMPLE: Iteration 8' in chunk;last=chunk.rsplit('PIMPLE: Iteration 8',1)[1];res={}
  for m in re.finditer(r'Solving for (\w+), Initial residual = ('+NUM+r'), Final residual = ('+NUM+')',last):
   if m[1] in ['Ux','Uy','p_rgh']:res[m[1]]=max(res.get(m[1],0),float(m[2]))
  assert set(res)=={'Ux','Uy','p_rgh'};a=re.findall(r'Min\(alpha.water\)\s*=\s*('+NUM+r').*?Max\(alpha.water\)\s*=\s*('+NUM+')',last);assert a;lo,hi=map(float,a[-1]);assert np.isfinite([lo,hi,*res.values()]).all()
  rows.append(dict(step=81+len(rows),time_s=t,residuals=res,alpha_min=lo,alpha_max=hi,passed=all(v<=1e-4 for v in res.values()) and lo>=-1e-6 and hi<=1+1e-6))
 assert len(rows)==4720 and all(abs(r['time_s']-r['step']*DT)<1e-10 for r in rows);merged=prefix['step_records']+rows;assert all(r['passed'] for r in merged)
 co=courant_records(text,'openfoam',DT,80,4800);assert max(r['maximum_native_interface_courant'] for r in co)<=.25
 states=[80]+list(range(120,4801,40));records=[]
 for step in states:
  folder=O/f'{step*DT:.8g}';q=[field(folder/k,c) for k,c in [('alpha.water',1),('U',3),('p_rgh',1)]];r=saved_check(step,step*DT,q[0].ravel(),q[1],q[0].ravel(),'continuation',q[2].ravel());assert r['alpha_bounds_pass'] and r['volume_pass']
  raw=(folder/'uniform/time').read_text()
  for key,val in [('value',step*DT),('index',step),('deltaT',DT)]:assert abs(float(re.search(r'^\s*'+key+r'\s+('+NUM+r');',raw,re.M)[1])-val)<1e-10
  f=phi_wall_flux((folder/'phi').read_text(),480,25360)/.01;assert f<=1e-12;r['absolute_wall_volume_flux_per_depth_m2_s']=f;records.append(r)
 for part in ['processor0','processor1']:
  raw=(O/part/'1.2/uniform/time').read_text()
  for key,val in [('value',1.2),('index',4800),('deltaT',DT)]:assert abs(float(re.search(r'^\s*'+key+r'\s+('+NUM+r');',raw,re.M)[1])-val)<1e-10
 for key,c in [('alpha.water',1),('U',3),('p_rgh',1)]:assert np.array_equal(field(O/'0.02'/key,c),field(P/'0.02'/key,c))
 sources=[p for p in O.rglob('*') if p.is_file() and p.name not in ['production-audit-v001.json']];a=dict(numerical_pass=True,new_steps=4720,merged_steps=4800,continued_states=len(records),selected_matched_states=121,records=records,step_records=rows,courant_records=co,exact_checkpoint_restart_fields=True,physical_validation=False,independent_full_rerun=False,source_sha256={str(p):sha(p) for p in sources})
 (O/'production-audit-v001.json').write_text(json.dumps(a,indent=2)+'\n');(S/'openfoam-combined-4800-step-checks-v001.json').write_text(json.dumps(merged,indent=2)+'\n');return a
