"""Strict actual-iteration and restart evidence helpers, no solver launch."""
import numpy as np
EQ=('continuity','x-velocity','y-velocity','z-velocity')
def force_window(text,start,end,parent):
 rows=[]
 for line in text.splitlines():
  a=line.split()
  if len(a)!=3:continue
  try:i=float(a[0])
  except ValueError:continue
  try:cd,cl=map(float,a[1:])
  except ValueError:raise ValueError('Invalid numeric monitor row')
  if not np.isfinite([i,cd,cl]).all() or not i.is_integer():raise ValueError('Nonfinite/fractional monitor iteration')
  rows.append(dict(iteration=int(i),Cd=cd,Cl=cl))
 if rows and rows[0]['iteration']==start-1:
  if parent is None or any(rows[0][k]!=parent[k] for k in ['iteration','Cd','Cl']):raise ValueError('Changed restart monitor echo')
  rows=rows[1:]
 if [r['iteration'] for r in rows]!=list(range(start,end+1)):raise ValueError('Missing/duplicate/unexpected actual force iteration')
 return rows

def residual_window(iterations,values,start,end):
 i=np.asarray(iterations,dtype=float)
 if i.ndim!=1 or not len(i) or not np.isfinite(i).all() or not np.array_equal(i,i.astype(int)) or not np.all(np.diff(i)==1):raise ValueError('Invalid residual iteration coverage')
 if set(values)!=set(EQ) or int(i[-1])!=end:raise ValueError('Residual equations/terminal mismatch')
 q={n:np.asarray(values[n],dtype=float) for n in EQ}
 if any(v.shape!=i.shape or not np.isfinite(v).all() or (v<0).any() for v in q.values()):raise ValueError('Invalid native residual array')
 mask=(i>=start)&(i<=end)
 if not np.array_equal(i[mask],np.arange(start,end+1)):raise ValueError('Missing new global residual iteration')
 return [dict(iteration=int(it),**{n:float(q[n][k]) for n in EQ}) for k,it in enumerate(i) if mask[k]]

def same_residual_prefix(parent_it,parent_raw,new_it,new_raw,last):
 pi=np.asarray(parent_it);ni=np.asarray(new_it);pm=pi<=last;nm=ni<=last
 if set(parent_raw)!=set(new_raw) or not np.array_equal(pi[pm],ni[nm]):raise ValueError('Native residual prefix coverage changed')
 for n in parent_raw:
  a=np.asarray(parent_raw[n])[pm];b=np.asarray(new_raw[n])[nm]
  if not np.isfinite(a).all() or not np.array_equal(a,b):raise ValueError('Native residual numerator/scaling history changed: '+n)
 return True

def original_target(equation):
 old=equation.check_convergence()
 try:
  equation.check_convergence.set_state(True)
  return float(equation.absolute_criteria())
 finally:equation.check_convergence.set_state(old)
