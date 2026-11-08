import re,numpy as np
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def parse_steps(text,start_step,end_step,dt):
 boundaries=list(re.finditer(r'Flow time = ('+NUM+r')s, time step = (\d+)',text));assert len(boundaries)==end_step-start_step
 pattern=re.compile(r'^\s*(\d+)\s+('+NUM+r')\s+('+NUM+r')\s+('+NUM+r')\s+\d+:\d+:\d+\s+(\d+)\s*$',re.M)
 result=[];start=0;previous=None
 for step,b in zip(range(start_step+1,end_step+1),boundaries):
  assert int(b[2])==step and abs(float(b[1])-step*dt)<1e-10
  rows=list(pattern.finditer(text,start,b.start()));assert rows
  total=int(rows[0][1])+int(rows[0][5]);assert all(int(r[1])+int(r[5])==total for r in rows)
  r=rows[-1];inner=400-int(r[5]);assert 1<=inner<=400
  if previous is not None:assert int(r[1])-previous==inner
  previous=int(r[1]);values=[float(r[k]) for k in [2,3,4]];assert np.isfinite(values).all()
  result.append(dict(step=step,time_s=float(b[1]),residuals=dict(zip(['continuity','x-velocity','y-velocity'],values)),inner_iterations=inner,passed=all(v<=1e-4 for v in values)))
  start=b.end()
 return result
