"""Canonical native-cell fields for cavity postprocessing; no solver launch."""
from pathlib import Path
import numpy as np
from plot_nozzle_native_v004 import fluent_native,foam_mesh,foam_field
L=.05;R=287.0028467801174;P0=101325.
def order(q,n):
 xy=q['nodes'][q['cells'],:2].mean(axis=1);ij=np.rint(xy/(L/n)-.5).astype(int)
 assert np.all((ij>=0)&(ij<n)) and np.max(abs(xy-(ij+.5)*L/n))<1e-9
 slots=ij[:,1]*n+ij[:,0];assert np.array_equal(np.sort(slots),np.arange(n*n))
 assert np.max(abs(np.ptp(q['nodes'][q['cells'],:2],axis=1)-L/n))<1e-9
 idx=np.argsort(slots)
 return {k:q[k][idx] for k in ['U','T','p','rho']}
def fluent(case,data,n):return order(fluent_native(case,data,P0),n)
def foam(case,time,n):
 q=foam_mesh(case);p=case/f'{time:g}'
 for key,w,dim in [('U',3,[0,1,-1,0,0,0,0]),('T',1,[0,0,0,1,0,0,0]),('p',1,[1,-1,-2,0,0,0,0])]:q[key]=foam_field(p/key,n*n,w,dim)
 q['U']=q['U'][:,:2];q['rho']=foam_field(p/'rho',n*n,1,[1,-3,0,0,0,0,0]) if time else q['p']/(R*q['T'])
 return order(q,n)
def metrics(q,n):
 assert all(np.isfinite(v).all() for v in q.values())
 return dict(Tmin_K=float(q['T'].min()),Tmax_K=float(q['T'].max()),Umax_m_s=float(np.linalg.norm(q['U'],axis=1).max()),mass_per_depth_kg_m=float(q['rho'].sum()*(L/n)**2),eos_relative_error=float(np.max(abs(q['rho']*R*q['T']/q['p']-1))))
def mesh_comparison(coarse,fine):
 reduced={k:v.reshape(24,2,24,2,*v.shape[1:]).mean(axis=(1,3)).reshape(coarse[k].shape) for k,v in fine.items()}
 return dict(method='Area-average each 2x2 group of uniform fine cells into its coarse cell; compare at same physical time',velocity_relative_L2=float(np.linalg.norm(reduced['U']-coarse['U'])/np.linalg.norm(coarse['U'])),maximum_T_difference_K=float(abs(reduced['T']-coarse['T']).max()),temperature_departure_relative_L2=float(np.linalg.norm(reduced['T']-coarse['T'])/np.linalg.norm(coarse['T']-300)),coarse_max_speed_m_s=float(np.linalg.norm(coarse['U'],axis=1).max()),fine_max_speed_m_s=float(np.linalg.norm(fine['U'],axis=1).max()),mesh_independence_established=False)
