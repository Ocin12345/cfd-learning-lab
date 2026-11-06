"""Exact 1D constant-property Dirichlet slab solution; no solver data generated.

Initial interior300K, left400K for t>0, right300K; L.1m, alpha1e-5m2/s.
The cell-average Fourier solution is a numerical reference, not simulated output.
"""
import numpy as np
L=.1
ALPHA=1e-5
COLD=300.
DELTA=100.

def cell_average_temperature(x_left,x_right,time_s,terms=80):
 left,right=np.broadcast_arrays(np.asarray(x_left,float),np.asarray(x_right,float))
 if not(np.isfinite(left).all() and np.isfinite(right).all() and np.isfinite(time_s) and time_s>=0 and terms>=1 and np.all(left>=0) and np.all(right<=L+1e-14) and np.all(right>left)):
  raise ValueError('Finite ordered cells within slab and nonnegative physical time required')
 if time_s==0:return np.full(left.shape,COLD)
 n=np.arange(1,int(terms)+1,dtype=float)
 decay=np.exp(-ALPHA*(n*np.pi/L)**2*time_s)
 weights=decay/n**2
 c=np.cos(left[...,None]*n*np.pi/L)-np.cos(right[...,None]*n*np.pi/L)
 theta=1-(left+right)/(2*L)-2*L/(np.pi**2*(right-left))*np.sum(c*weights,axis=-1)
 return COLD+DELTA*theta

def mean_temperature(time_s,terms=80):
 if not np.isfinite(time_s) or time_s<0:raise ValueError('Nonnegative finite physical time required')
 if time_s==0:return COLD
 odd=np.arange(1,int(terms)+1,dtype=float);odd=odd[(odd.astype(int)%2)==1]
 return float(COLD+DELTA*(.5-4/np.pi**2*np.sum(np.exp(-ALPHA*(odd*np.pi/L)**2*time_s)/odd**2)))

def point_temperature(x,time_s,terms=80):
 x=np.asarray(x,float)
 if time_s<0 or not np.isfinite(time_s) or np.any(x<0) or np.any(x>L):raise ValueError('Point outside slab/time invalid')
 if time_s==0:return np.full(x.shape,COLD)
 n=np.arange(1,int(terms)+1,dtype=float)
 theta=1-x/L-2/np.pi*np.sum(np.sin(x[...,None]*n*np.pi/L)*np.exp(-ALPHA*(n*np.pi/L)**2*time_s)/n,axis=-1)
 return COLD+DELTA*theta
