import numpy as np
from sloshing_cycles_v001 import cycle_evidence

def depth_comparison(t,base,new,mean):
 t=np.asarray(t,float);a=np.asarray(base,float);b=np.asarray(new,float)
 if a.ndim!=2 or b.shape!=a.shape or a.shape[0]!=len(t) or a.shape[1]%2 or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('Finite matching time/column arrays with even native columns required')
 paired_a=a.reshape(len(t),-1,2).mean(axis=2);paired_b=b.reshape(len(t),-1,2).mean(axis=2);d=b-a;paired=paired_b-paired_a
 c1=cycle_evidence(t,paired_a[:,-1],mean);c2=cycle_evidence(t,paired_b[:,-1],mean)
 return dict(maximum_absolute_native_column_difference_m=float(abs(d).max()),native_column_sample_RMS_m=float(np.sqrt(np.mean(d*d))),endpoint_native_column_RMS_m=float(np.sqrt(np.mean(d[-1]**2))),maximum_absolute_paired_strip_difference_m=float(abs(paired).max()),paired_strip_sample_RMS_m=float(np.sqrt(np.mean(paired*paired))),endpoint_paired_strip_RMS_m=float(np.sqrt(np.mean(paired[-1]**2))),right_strip_base_depth_m=paired_a[:,-1].tolist(),right_strip_new_depth_m=paired_b[:,-1].tolist(),right_strip_difference_m=paired[:,-1].tolist(),base_cycle_evidence=c1,new_cycle_evidence=c2,ordinal_matched_crossing_time_difference_s=[float(y-x) for x,y in zip(c1['positive_crossing_times_s'],c2['positive_crossing_times_s'])],crossing_counts=[len(c1['positive_crossing_times_s']),len(c2['positive_crossing_times_s'])],right_strip_definition='Mean of last two fine native columns; same physical 2.5mm strip',RMS_definition='Unweighted uniformly spaced saved-time/native-column samples, not continuous time integration',phase_alignment_applied=False)
