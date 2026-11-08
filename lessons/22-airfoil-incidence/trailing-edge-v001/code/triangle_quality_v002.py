import numpy as np

def cell_quality(xy,tri):
 xy=np.asarray(xy,dtype=float);tri=np.asarray(tri,dtype=int);p=xy[tri];centres=p.mean(axis=1);n=len(tri);edges=np.stack([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]],axis=1).reshape(-1,2);owners=np.repeat(np.arange(n),3);unique,inv,counts=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
 if not np.isfinite(p).all() or any(counts>2):raise ValueError('Invalid/nonmanifold mesh')
 v=xy[edges[:,1]]-xy[edges[:,0]];normal=np.c_[v[:,1],-v[:,0]];face=xy[edges].mean(axis=1)-centres[owners];norm=lambda x:np.linalg.norm(x,axis=1)
 qf=np.abs(np.sum(normal*face,axis=1)/(norm(normal)*norm(face)));summed=np.zeros((len(unique),2));np.add.at(summed,inv,centres[owners]);cc=summed[inv]-2*centres[owners];internal=counts[inv]==2;qc=np.ones(len(edges));qc[internal]=np.abs(np.sum(normal[internal]*cc[internal],axis=1)/(norm(normal[internal])*norm(cc[internal])));ortho=np.minimum(qf,qc).reshape(-1,3).min(axis=1)
 a=p[:,1]-p[:,0];b=p[:,2]-p[:,0];area=np.abs(a[:,0]*b[:,1]-a[:,1]*b[:,0])/2;lengths=np.linalg.norm(p-np.roll(p,1,axis=1),axis=2)
 # Area divided by the equilateral triangle area at the same circumradius.
 one_minus_equivolume_skew=64*area**3/(3*np.sqrt(3)*np.prod(lengths**2,axis=1));q=np.minimum(ortho,one_minus_equivolume_skew)
 if not np.isfinite(q).all():raise ValueError('Invalid quality metric')
 return q