import numpy as np

def wall_gradient_heat(T,k,hot=300.4,cold=299.6,length=.05):
 T=np.asarray(T,dtype=float)
 if T.ndim!=2 or T.shape[0]!=T.shape[1] or min(T.shape)<2 or not np.isfinite(T).all() or not np.isfinite(k) or k<=0 or length<=0:raise ValueError('Finite uniform square native-cell temperature and positive conductivity required')
 n=T.shape[0];face_length=length/n;normal_distance=length/(2*n)
 # Positive INTO fluid; top/bottom use the imposed adiabatic boundary condition.
 return dict(left=float(k*face_length/normal_distance*np.sum(hot-T[:,0])),right=float(k*face_length/normal_distance*np.sum(cold-T[:,-1])),top=0.,bottom=0.)