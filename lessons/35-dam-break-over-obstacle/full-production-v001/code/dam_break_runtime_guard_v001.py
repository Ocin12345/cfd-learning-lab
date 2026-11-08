"""Fail early on actual terminal residual, VOF Courant or open-mass misses."""
import math,re
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
class NumericalMiss(RuntimeError):pass
def check_chunk(text,flux_text,*,completed_steps,start_step,dt,current_mass,initial_mass,targets):
 ends=list(re.finditer(r'Flow time = ('+NUMBER+r')s, time step = (\d+)',text))
 assert len(ends)==completed_steps and 1<=start_step<=completed_steps,'Incomplete native time-step coverage'
 pattern=re.compile(r'^\s*(\d+)\s+('+NUMBER+r')\s+('+NUMBER+r')\s+('+NUMBER+r')\s+\d+:\d+:\d+\s+\d+\s*$',re.M)
 misses=[];ratios=[]
 for k in range(start_step,completed_steps+1):
  end=ends[k-1];assert int(end[2])==k and abs(float(end[1])-k*dt)<1e-10
  begin=0 if k==1 else ends[k-2].end()
  rows=list(pattern.finditer(text,begin,end.start()));assert rows,'Missing terminal residuals'
  values=[float(rows[-1][i]) for i in [2,3,4]];assert all(math.isfinite(v) for v in values)
  ratio=max(v/limit for v,limit in zip(values,targets.values()));ratios.append(ratio)
  if ratio>1:misses.append(dict(step=k,residuals=dict(zip(targets,values))))
 co=[float(v) for v in re.findall(r'Global Courant Number \[Explicit VOF Criteria\]\s*:\s*('+NUMBER+')',text)]
 assert len(co)==completed_steps and all(math.isfinite(v) for v in co)
 peak=max(co[start_step-1:])
 flux=[re.fullmatch(r'\s*(\d+)\s+('+NUMBER+r')\s*',line) for line in flux_text.splitlines()];flux=[m for m in flux if m]
 assert len(flux)==completed_steps and [int(m[1]) for m in flux]==list(range(1,completed_steps+1))
 outgoing=[-float(m[2]) for m in flux];assert all(math.isfinite(v) for v in outgoing)
 closure=(current_mass-initial_mass+sum(outgoing)*dt)/initial_mass
 result=dict(start_step=start_step,completed_steps=completed_steps,steps_checked=completed_steps-start_step+1,maximum_residual_limit_ratio=max(ratios),peak_reported_interface_courant=peak,relative_mixture_mass_closure=closure,passed=not misses and peak<=.25 and abs(closure)<=.005)
 if not result['passed']:raise NumericalMiss(str(dict(**result,residual_misses=misses[:5])))
 return result