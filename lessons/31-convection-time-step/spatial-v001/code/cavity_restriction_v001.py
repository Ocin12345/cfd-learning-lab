import numpy as np

def restrict_square(fine,n):
 q=np.asarray(fine,dtype=float)
 if not isinstance(n,int) or n<=0 or q.ndim<1 or q.shape[0]!=4*n*n or not np.isfinite(q).all():raise ValueError('Finite square2x2 fine-cell field and positive integer coarse n required')
 return q.reshape(n,2,n,2,*q.shape[1:]).mean(axis=(1,3)).reshape(n*n,*q.shape[1:])

def compare_fields(coarse,fine,n,relative=True):
 reduced={k:restrict_square(v,n) for k,v in fine.items()}
 if any(reduced[k].shape!=np.asarray(coarse[k]).shape for k in coarse):raise ValueError('Coarse/fine components mismatch')
 du=reduced['U']-coarse['U'];dT=reduced['T']-coarse['T'];u=np.linalg.norm(coarse['U']);thermal=np.linalg.norm(coarse['T']-300)
 if relative and (u<=0 or thermal<=0):raise ValueError('Relative heated-field denominator is zero; use explicit absolute-only initial comparison')
 return dict(velocity_relative_L2=float(np.linalg.norm(du)/u) if relative else None,temperature_departure_relative_L2=float(np.linalg.norm(dT)/thermal) if relative else None,max_velocity_difference_m_s=float(np.linalg.norm(du,axis=-1).max()),max_temperature_difference_K=float(abs(dT).max()),method='Conservative equal-area2x2 average onto matching coarse cells; relative denominators use coarse U and coarse(T-300)',relative_metrics_enabled=relative)