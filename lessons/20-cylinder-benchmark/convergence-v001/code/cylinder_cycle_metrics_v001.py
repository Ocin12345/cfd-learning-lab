import numpy as np

def cycle_metrics(t,cd,cl):
 t,cd,cl=[np.asarray(a,dtype=float) for a in [t,cd,cl]]
 if t.ndim!=1 or len(t)<2 or cd.shape!=t.shape or cl.shape!=t.shape or not np.isfinite([t,cd,cl]).all() or not (np.diff(t)>0).all():raise ValueError('Finite force samples and strictly increasing native times required')
 j=np.where((cl[:-1]<0)&(cl[1:]>=0))[0];roots=t[j]-cl[j]*(t[j+1]-t[j])/(cl[j+1]-cl[j]);cycles=[]
 for a,b in zip(roots[:-1],roots[1:]):
  inside=(t>a)&(t<b);tc=np.r_[a,t[inside],b];dc=np.interp(tc,t,cd);lc=np.interp(tc,t,cl)
  if inside.sum()<2:raise ValueError('Insufficient native samples in a cycle')
  native_lift=cl[inside];period=b-a
  cycles.append(dict(start_time_s=float(a),end_time_s=float(b),period_s=float(period),St=float(.1/period),sampled_Cd_min=float(cd[inside].min()),sampled_Cd_max=float(cd[inside].max()),sampled_Cl_min=float(native_lift.min()),sampled_Cl_max=float(native_lift.max()),sampled_Cl_peak_to_peak=float(np.ptp(native_lift)),time_mean_Cd=float(np.trapezoid(dc,tc)/period),native_samples=int(inside.sum())))
 result=dict(method='Successive rising lift zero crossings, linearly estimated between native force samples. Extrema use native samples only. Cd time means use trapezoidal integration with linearly interpolated force values only at cycle boundaries; no native fields or video frames are interpolated.',criteria=dict(required_complete_cycles=6,maximum_relative_range=.005,measures=['period','sampled lift peak-to-peak','time-mean Cd']),complete_cycles=cycles,last_six_cycles_stable_under_declared_criterion=False,physical_validation=False,benchmark_pass=False)
 if len(cycles)<6:result['reason']='Fewer than six complete cycles';return result
 last=cycles[-6:];ranges={}
 for key in ['period_s','sampled_Cl_peak_to_peak','time_mean_Cd']:
  values=np.array([r[key] for r in last]);denom=abs(values.mean());ranges[key]=float(np.ptp(values)/denom) if denom>0 else None
 result.update(last_six_cycles=last,last_six_relative_ranges=ranges,last_six_mean_period_s=float(np.mean([r['period_s'] for r in last])),last_six_mean_Cd=float(np.mean([r['time_mean_Cd'] for r in last])),last_six_cycles_stable_under_declared_criterion=all(v is not None and v<=.005 for v in ranges.values()))
 return result