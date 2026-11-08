import re,math
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def courant_records(text,solver,dt,start_step,end_step):
 n=end_step-start_step;assert n>0
 if solver=='fluent':
  reports=list(re.finditer(r'Global Courant Number \[Explicit VOF Criteria\]\s*:\s*('+NUM+')',text))
  clocks=list(re.finditer(r'Flow time = ('+NUM+r')s, time step = (\d+)',text));assert len(reports)==len(clocks)==n
 else:
  assert solver=='openfoam'
  reports=list(re.finditer(r'^Interface Courant Number mean: '+NUM+r' max: ('+NUM+')',text,re.M))
  clocks=list(re.finditer(r'^Time = ('+NUM+r')s?\s*$',text,re.M));assert len(reports)==n+1 and len(clocks)==n
 rows=[]
 for i,m in enumerate(reports):
  preceding=[c for c in clocks if c.start()<m.start()]
  time=float(preceding[-1][1]) if preceding else start_step*dt
  expected=(start_step+(i if solver=='fluent' else max(i-1,0)))*dt
  assert abs(time-expected)<1e-10,(i,time,expected)
  value=float(m[1]);assert math.isfinite(value) and value>=0
  rows.append(dict(actual_pre_advance_time_s=time,maximum_native_interface_courant=value,initial_duplicate=solver=='openfoam' and i==0))
 return rows
