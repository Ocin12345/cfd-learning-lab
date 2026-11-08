import re,math

def audit_steps(text,expected,dt):
 passed=[]
 for block in re.split(r'Updating solution at time level N',text)[1:]:
  event=re.search(r'Flow time = ([0-9.eE+-]+)s, time step = (\d+)',block)
  if not event:continue
  rows=re.findall(r'^\s*(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+\d+:',block,re.M)
  if not rows:raise ValueError('No terminal scaled residual row for actual physical step')
  values=[float(v) for v in rows[-1][1:]];step=int(event[2]);time=float(event[1])
  if not all(math.isfinite(v) and v>=0 for v in values):raise ValueError('Invalid native residual values')
  passed.append(dict(step=step,time_s=time,residuals=values,passed=all(v<=lim for v,lim in zip(values,[1e-6,1e-6,1e-6,1e-9]))))
 if [r['step'] for r in passed]!=list(range(1,expected+1)) or any(abs(r['time_s']-r['step']*dt)>1e-8 for r in passed):raise ValueError('Missing, duplicate or wrong-time native physical step')
 return passed