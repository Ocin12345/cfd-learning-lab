import re,numpy as np
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def parse_steps(text,start_step,end_step,dt):
 parts=re.split(r'^Time = ([0-9.eE+-]+)s?\s*$',text,flags=re.M);rows=[]
 for i in range(1,len(parts),2):
  t=float(parts[i]);chunk=parts[i+1];assert 'PIMPLE: Iteration 8' in chunk;last=chunk.rsplit('PIMPLE: Iteration 8',1)[1];res={}
  for m in re.finditer(r'Solving for (\w+), Initial residual = ('+NUM+r'), Final residual = ('+NUM+')',last):
   if m[1] in ['Ux','Uy','p_rgh']:res[m[1]]=max(res.get(m[1],0),float(m[2]))
  assert set(res)=={'Ux','Uy','p_rgh'};a=re.findall(r'Min\(alpha.water\)\s*=\s*('+NUM+r').*?Max\(alpha.water\)\s*=\s*('+NUM+')',last);assert a;lo,hi=map(float,a[-1]);assert np.isfinite([lo,hi,*res.values()]).all()
  step=start_step+1+len(rows);assert abs(t-step*dt)<1e-10
  rows.append(dict(step=step,time_s=t,residuals=res,alpha_min=lo,alpha_max=hi,passed=all(v<=1e-4 for v in res.values()) and lo>=-1e-6 and hi<=1+1e-6))
 assert len(rows)==end_step-start_step;return rows
