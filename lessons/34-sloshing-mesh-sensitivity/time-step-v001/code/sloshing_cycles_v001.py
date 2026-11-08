import numpy as np

def cycle_evidence(t,depth,mean):
 t=np.asarray(t,float);y=np.asarray(depth,float)-mean
 if t.ndim!=1 or y.shape!=t.shape or len(t)<2 or not np.isfinite(t).all() or not np.isfinite(y).all() or not np.all(np.diff(t)>0):raise ValueError('Finite strictly increasing paired samples required')
 nz=np.flatnonzero(y!=0);cross=[]
 for a,b in zip(nz[:-1],nz[1:]):
  if y[a]<0<y[b]:
   if b-a>2:raise ValueError('A flat zero interval does not define a unique crossing time')
   c=float(t[a+1]) if b-a==2 else float(t[a]-y[a]*(t[b]-t[a])/(y[b]-y[a]));cross.append(c)
 cycles=[]
 for a,b in zip(cross[:-1],cross[1:]):
  samples=y[(t>a)&(t<b)];values=np.r_[0,samples,0];cycles.append(dict(start_s=a,end_s=b,period_s=b-a,peak_to_peak=float(values.max()-values.min()),sampled_peak=float(mean+values.max()),sampled_trough=float(mean+values.min())))
 return dict(positive_crossing_times_s=cross,complete_cycles=cycles,method='Positive-going mean crossings linearly bracketed by saved samples; only complete consecutive-crossing cycles; extrema sampled without extrapolation',maximum_sampling_interval_s=float(np.diff(t).max()),phase_alignment_applied=False)
