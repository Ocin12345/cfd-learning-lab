import numpy as np
from plot_openfoam_airfoil_fields_v002 import coordinate_bijection

def verify_extrusion(C,points,reference):
 C,points,reference=map(lambda a:np.asarray(a,float),(C,points,reference))
 if C.shape!=(len(reference),3) or points.ndim!=2 or points.shape[1]!=3 or not all(np.isfinite(a).all() for a in (C,points,reference)):raise ValueError('Finite native centre/point/reference arrays required')
 lower,upper=float(points[:,2].min()),float(points[:,2].max());middle=(lower+upper)/2
 if max(abs(lower+.0025),abs(upper-.0025))>1e-10 or abs(C[:,2]-middle).max()>1e-10 or np.any(np.minimum(abs(points[:,2]-lower),abs(points[:,2]-upper))>1e-10):raise ValueError('Native symmetric extrusion bounds/centroid plane mismatch')
 if max(abs(points[:,:2].min(axis=0)).max(),abs(points[:,:2].max(axis=0)-.05).max())>1e-10:raise ValueError('Square domain changed')
 order,mapping=coordinate_bijection(reference,C[:,:2],1e-10)
 return dict(native_depth_m=upper-lower,expected_z_bounds_m=[-.0025,.0025],native_centroid_plane_z_m=middle,maximum_centroid_plane_error_m=float(abs(C[:,2]-middle).max()),cell_bijection=mapping)