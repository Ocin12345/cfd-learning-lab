"""Gradient on native curvilinear cell centres using the full map Jacobian."""
import numpy as np

def gradient_magnitude(f,x,y):
 f,x,y=map(lambda v:np.asarray(v,dtype=float),(f,x,y))
 if f.ndim!=2 or f.shape!=x.shape or x.shape!=y.shape or min(f.shape)<3 or not all(np.isfinite(v).all() for v in (f,x,y)):raise ValueError('Finite matching 2D native fields required')
 fj,fi=np.gradient(f,edge_order=2);xj,xi=np.gradient(x,edge_order=2);yj,yi=np.gradient(y,edge_order=2)
 det=xi*yj-xj*yi
 if np.any(abs(det)<1e-14) or not (np.all(det>0) or np.all(det<0)):raise ValueError('Degenerate or folded centroid map')
 gx=(fi*yj-fj*yi)/det;gy=(xi*fj-xj*fi)/det
 return np.hypot(gx,gy)