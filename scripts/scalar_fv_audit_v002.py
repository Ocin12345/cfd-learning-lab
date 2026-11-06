"""Source-only scalar FV diagnostics. No solver or dataset discovery; no acceptance gates.

All c arrays must be actually observed native cell fields in stencil cell order.
Mass phi is native outward-from-owner kg/s (rho*native volume phi for OF).
Every calculated scalar face flux/PDE defect is RECONSTRUCTED, never native matrix flux.
Dirichlet outflow convection is outside this qualified scope and is refused.
"""
import numpy as np

def require(v,message):
 if not v:raise ValueError(message)
def array(v,shape,label):
 x=np.asarray(v,dtype=float);require(x.shape==shape and np.isfinite(x).all(),'Finite '+label+' shape '+str(shape));return x

def make_stencil(*,owner,neighbour,cell_centres,cell_volumes,face_centres,area_vectors,patch,depth_m,native_face_ids):
 C=np.asarray(cell_centres,dtype=float);require(C.ndim==2 and C.shape[1]==2 and np.isfinite(C).all(),'Finite planar cell centres');nc=len(C);nf=len(owner)
 o=array(owner,(nf,),'owner');n=array(neighbour,(nf,),'neighbour');require(np.array_equal(o,o.astype(int)) and np.array_equal(n,n.astype(int)),'Integral native cell indices');o=o.astype(int);n=n.astype(int)
 require(np.all((o>=0)&(o<nc)) and np.all((n>=-1)&(n<nc)) and np.all(o!=n),'Owner/neighbour ranges');inside=n>=0
 V=array(cell_volumes,(nc,),'native volumes');F=array(face_centres,(nf,2),'face centres');S=array(area_vectors,(nf,2),'outward face area vectors');A=np.linalg.norm(S,axis=1);require(np.all(V>0) and np.all(A>0) and np.isfinite(depth_m) and depth_m>0,'Positive volume/area/depth')
 names=np.asarray(patch,dtype=str);require(names.shape==(nf,) and np.all(names[inside]=='interior') and np.all(names[~inside]!='interior'),'Explicit physical patch cover')
 ids=array(native_face_ids,(nf,),'native face IDs');require(np.array_equal(ids,ids.astype(int)) and len(np.unique(ids))==nf,'Unique integral face IDs')
 normal=S/A[:,None];delta=F-C[o];delta[inside]=C[n[inside]]-C[o[inside]];distance=np.sum(delta*normal,axis=1);tangent=delta-distance[:,None]*normal
 require(np.all(distance>0) and np.max(np.linalg.norm(tangent,axis=1))<=1e-10,'Strict orthogonal owner-to-neighbour/boundary stencil')
 closure=np.zeros((nc,2));np.add.at(closure,o,S);np.add.at(closure,n[inside],-S[inside]);require(np.max(abs(closure))<=1e-10*max(1.,float(A.max())),'Every cell oriented face-vector closure')
 return dict(owner=o,neighbour=n,inside=inside,cell_centres=C,volumes=V,face_centres=F,area_vectors=S,area=A,distance=distance,patch=names,depth_m=float(depth_m),native_face_ids=ids.astype(int),cell_count=nc,face_count=nf)

def face_terms(stencil,c,mass_phi,*,gamma,boundary):
 g=stencil;nc=g['cell_count'];nf=g['face_count'];c=array(c,(nc,),'native c');q=array(mass_phi,(nf,),'native mass phi');require(np.isfinite(gamma) and gamma>=0,'Nonnegative constant Gamma kg/(m s)')
 o=g['owner'];n=g['neighbour'];inside=g['inside'];up=c[o].copy();grad=np.zeros(nf);reverse=inside&(q<0);up[reverse]=c[n[reverse]];grad[inside]=(c[n[inside]]-c[o[inside]])/g['distance'][inside]
 physical=~inside;required=set(g['patch'][physical]);require(set(boundary)==required,'Exact boundary declarations, no missing/extra patches')
 for name in sorted(required):
  pick=physical&(g['patch']==name);b=boundary[name];require(isinstance(b.get('evidence'),str) and bool(b['evidence']),'Explicit boundary evidence/declaration label required')
  v=np.asarray(b.get('value',float('nan')),dtype=float);v=np.broadcast_to(v,(int(pick.sum()),));inflow=q[pick]<0
  require(np.isfinite(v[inflow]).all(),'Observed inflow boundary c required');local=up[pick];local[inflow]=v[inflow];up[pick]=local
  mode=b.get('diffusion')
  if mode=='dirichlet':require(np.isfinite(v).all(),'Observed Dirichlet boundary c required');require(not np.any(q[pick]>0),'Dirichlet outflow convection is not source-qualified');grad[pick]=(v-c[o[pick]])/g['distance'][pick]
  elif mode=='zero_gradient':grad[pick]=0.
  else:raise ValueError('Unqualified/unsupported diffusion mode; a setting named flux0 is insufficient')
 adv=q*up;diff=-float(gamma)*g['area']*grad
 return dict(native_mass_phi_outward_kg_s=q,reconstructed_upwind_face_c=up,reconstructed_outward_snGrad_c_m_inverse=grad,reconstructed_advective_outward_kg_s=adv,reconstructed_diffusive_outward_kg_s=diff,reconstructed_total_outward_kg_s=adv+diff)

def audit_step(stencil,c_prev,c_now,mass_phi,*,rho,gamma,dt_s,boundary,normalization_time_s,global_floor_kg_m_s):
 g=stencil;nc=g['cell_count'];old=array(c_prev,(nc,),'previous native c');c=array(c_now,(nc,),'current native c');density=np.broadcast_to(np.asarray(rho,dtype=float),(nc,))
 require(np.isfinite(density).all() and np.all(density>0) and np.isfinite(dt_s) and dt_s>0 and np.isfinite(normalization_time_s) and normalization_time_s>0 and np.isfinite(global_floor_kg_m_s) and global_floor_kg_m_s>0,'Positive actual density/dt and explicit diagnostic normalizations')
 # This rho-weighted storage is for unchanged density/volumes; verify frozen flow separately.
 f=face_terms(g,c,mass_phi,gamma=gamma,boundary=boundary);storage=density*g['volumes']*(c-old)/dt_s;o=g['owner'];inside=g['inside'];n=g['neighbour'];net=np.zeros(nc);np.add.at(net,o,f['reconstructed_total_outward_kg_s']);np.add.at(net,n[inside],-f['reconstructed_total_outward_kg_s'][inside]);defect=storage+net;depth=g['depth_m'];budgets={}
 for name in sorted(boundary):
  pick=(~inside)&(g['patch']==name);a=float(-f['reconstructed_advective_outward_kg_s'][pick].sum()/depth);d=float(-f['reconstructed_diffusive_outward_kg_s'][pick].sum()/depth);budgets[name]=dict(advective_inward_kg_m_s=a,diffusive_inward_kg_m_s=d,total_inward_kg_m_s=a+d)
 oldI=float(np.sum(density*g['volumes']*old)/depth);I=float(np.sum(density*g['volumes']*c)/depth);S=float(storage.sum()/depth);B=sum(v['total_inward_kg_m_s'] for v in budgets.values());raw=S-B;den=max(abs(S),sum(abs(v['total_inward_kg_m_s']) for v in budgets.values()),float(global_floor_kg_m_s));rate=defect/(density*g['volumes'])*normalization_time_s
 outlet=(~inside)&(g['patch']=='outlet');weights=np.maximum(f['native_mass_phi_outward_kg_s'][outlet],0);w=float(weights.sum());mean=float(np.dot(weights,f['reconstructed_upwind_face_c'][outlet])/w) if w>0 else None
 summary=dict(classification='RECONSTRUCTED orthogonal Euler/upwind/Gamma operator; not native matrix flux',scalar_field_provenance='Caller must authenticate actual native scalar arrays; this pure function does not bind datasets',native_scalar_dataset_binding_not_performed=True,acceptance_gates_applied=False,inventory_kg_per_m=I,previous_inventory_kg_per_m=oldI,storage_rate_kg_m_s=S,reconstructed_boundary=budgets,reconstructed_total_inward_kg_m_s=B,global_Euler_defect_raw_kg_m_s=raw,global_Euler_defect_denominator_kg_m_s=den,global_Euler_defect_relative=abs(raw)/den,reconstructed_local_dimensionless_defect_max=float(abs(rate).max()),reconstructed_local_dimensionless_defect_volume_RMS=float(np.sqrt(np.average(rate**2,weights=g['volumes']))),local_normalization_time_s=float(normalization_time_s),reconstructed_outlet_outflow_weighted_mean_c=mean,outlet_outward_mass_flow_kg_s=w,scalar_min=float(c.min()),scalar_max=float(c.max()),negative_cell_count=int(np.count_nonzero(c<0)),above_one_cell_count=int(np.count_nonzero(c>1)),minimum_undershoot_below_zero=float(max(0.,-c.min())),maximum_overshoot_above_one=float(max(0.,c.max()-1)),constant_c1_max_absolute_deviation=float(abs(c-1).max()),relative_inventory_change=(I-oldI)/oldI if oldI!=0 else None,cell_defect_sum_minus_global_raw_kg_m_s=float(defect.sum()/depth-raw),boundary_diffusion_declarations={k:dict(mode=v['diffusion'],evidence=v['evidence']) for k,v in boundary.items()})
 return f,defect,summary

def frozen_field_changes(reference,current):
 require(set(reference)==set(current),'All declared frozen U/p/phi/rho/mu components required');out={}
 for name,ref in reference.items():
  a=np.asarray(ref,dtype=float);b=array(current[name],a.shape,'current frozen '+name);require(np.isfinite(a).all(),'Finite reference '+name);out[name]=float(abs(a-b).max())
 return out
