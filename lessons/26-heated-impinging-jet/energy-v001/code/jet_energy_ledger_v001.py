"""Independent heat/storage arithmetic. No solver launch or field interpolation."""
import re
import numpy as np
NUM=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'

def finite(v):
 a=np.asarray(v,dtype=float)
 if not np.isfinite(a).all():raise ValueError('Nonfinite energy evidence')
 return a

def fluent_heat300(raw_heat,mass_out,h300):
 if set(raw_heat)!=set(mass_out):raise ValueError('Patch names differ')
 finite([*raw_heat.values(),*mass_out.values(),h300])
 return {n:float(raw_heat[n]+h300*mass_out[n]) for n in raw_heat}

def foam_patch_heat300(name,T,phi,*,dx,depth,rho,cp,k):
 T=finite(T);phi=finite(phi)
 if T.shape!=phi.shape or T.ndim!=1:raise ValueError('Patch arrays differ')
 if min(dx,depth,rho,cp,k)<=0:raise ValueError('Invalid geometry/material')
 delta={'plate':lambda:320.-T,'inlet':lambda:300.-T,'left':lambda:np.where(phi<0,300.-T,0.),'right':lambda:np.where(phi<0,300.-T,0.),'ceiling':lambda:np.zeros_like(T)}
 if name not in delta:raise ValueError('Unknown patch')
 upwind=np.where(phi<0,300.,T)
 conv=-rho*cp*np.sum(phi*(upwind-300.))/depth
 cond=np.sum(k*delta[name]()/ (dx/2)*dx)
 return float(conv+cond)

def ledger(steps,times,energy,patch_heat,dt):
 s=finite(steps);t=finite(times);E=finite(energy);n=len(s)
 if not np.isfinite(dt) or dt<=0 or n==0:raise ValueError('Invalid time interval')
 if s.shape!=(n,) or not np.array_equal(s,np.arange(1,n+1)):raise ValueError('Missing or duplicate actual steps')
 if t.shape!=(n,) or not np.allclose(t,s*dt,rtol=0,atol=1e-10):raise ValueError('Native clock mismatch')
 if E.shape!=(n+1,):raise ValueError('Energy requires initial state plus every step')
 if not patch_heat:raise ValueError('No independent boundary evidence')
 Q=np.array([finite(q) for q in patch_heat.values()])
 if Q.shape!=(len(patch_heat),n):raise ValueError('Incomplete patch heat series')
 net=Q.sum(axis=0);absolute=np.abs(Q).sum(axis=0);storage=np.diff(E)/dt
 def ratio(a,b):
  return np.divide(a,b,out=np.where(a==0,0.,np.inf),where=b!=0)
 local=ratio(abs(net-storage),np.maximum(absolute,abs(storage)))
 cumulative=np.cumsum(dt*net);throughput=np.cumsum(dt*absolute);change=E[1:]-E[0]
 err=ratio(abs(cumulative-change),np.maximum(throughput,abs(change)))
 return dict(steps=s.astype(int).tolist(),time_s=t.tolist(),stored_energy_change_J_per_m=change.tolist(),storage_rate_W_per_m=storage.tolist(),net_input_W_per_m=net.tolist(),step_relative_error=local.tolist(),cumulative_net_input_J_per_m=cumulative.tolist(),cumulative_patch_absolute_J_per_m=throughput.tolist(),cumulative_relative_error=err.tolist(),limit=.005,passed=bool((local<=.005).all() and (err<=.005).all()),definition='Backward-Euler right-endpoint discrete ledger; not exact continuous integration')

def fluent_steps(text,start,end,dt,cap=30):
 out=[];last=0
 for m in re.finditer(r'Flow time = ('+NUM+r')s, time step = (\d+)',text):
  rows=re.findall(r'^\s*(\d+)\s+('+NUM+r')\s+\d+:\d+:\d+\s+(\d+)\s*$',text[last:m.start()],re.M);last=m.end()
  advancing=[r for r in rows if int(r[2])<cap]
  if not advancing:raise ValueError('No actual new inner iteration')
  remains=[int(r[2]) for r in advancing]
  if min(remains)<0 or remains!=list(range(cap-1,remains[-1]-1,-1)):raise ValueError('Unexpected inner iteration budget')
  r=advancing[-1];res=float(r[1]);clock=float(m[1]);step=int(m[2])
  if not np.isfinite(res) or abs(clock-step*dt)>1e-10:raise ValueError('Residual or native clock invalid')
  out.append(dict(step=step,time_s=clock,residual=res,inner_iterations=cap-int(r[2])))
 if [r['step'] for r in out]!=list(range(start,end+1)):raise ValueError('Incomplete native step history')
 return out
