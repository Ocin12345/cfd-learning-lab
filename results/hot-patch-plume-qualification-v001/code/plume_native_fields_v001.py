from pathlib import Path
import numpy as np
import plot_nozzle_native_v004 as native
from plume_common_v001 import *
def order(q):
 xy=q['nodes'][q['cells'],:2].mean(axis=1);ij=np.rint(xy/[DX,DY]-.5).astype(int)
 assert np.all((ij[:,0]>=0)&(ij[:,0]<NX)) and np.all((ij[:,1]>=0)&(ij[:,1]<NY))
 assert np.max(abs(xy-(ij+.5)*[DX,DY]))<1e-9
 slots=ij[:,1]*NX+ij[:,0];assert np.array_equal(np.sort(slots),np.arange(NC))
 assert np.max(abs(np.ptp(q['nodes'][q['cells'],:2],axis=1)-[DX,DY]))<1e-9
 idx=np.argsort(slots);return {k:q[k][idx] for k in ['U','T','p','rho']}
def fluent(case,data):return order(native.fluent_native(case,data,P0))
def foam(case,time):
 q=native.foam_mesh(case);p=case/f'{time:g}'
 for key,w,dim in [('U',3,[0,1,-1,0,0,0,0]),('T',1,[0,0,0,1,0,0,0]),('p',1,[1,-1,-2,0,0,0,0])]:q[key]=native.foam_field(p/key,NC,w,dim)
 q['U']=q['U'][:,:2];q['rho']=native.foam_field(p/'rho',NC,1,[1,-3,0,0,0,0,0]) if time else q['p']/(R*q['T'])
 return order(q)
def metrics(q):
 assert all(np.isfinite(v).all() for v in q.values())
 U=q['U'].reshape(NY,NX,2);T=q['T'].reshape(NY,NX)
 return dict(Tmin_K=float(T.min()),Tmax_K=float(T.max()),Umax_m_s=float(np.linalg.norm(q['U'],axis=1).max()),mass_per_depth_kg_m=float(q['rho'].sum()*DX*DY),eos_relative_error=float(np.max(abs(q['rho']*R*q['T']/q['p']-1))),T_reflection_max_difference_K=float(np.max(abs(T-T[:,::-1]))),Uy_reflection_max_difference_m_s=float(np.max(abs(U[:,:,1]-U[:,::-1,1]))),Ux_antisymmetry_max_difference_m_s=float(np.max(abs(U[:,:,0]+U[:,::-1,0]))))
