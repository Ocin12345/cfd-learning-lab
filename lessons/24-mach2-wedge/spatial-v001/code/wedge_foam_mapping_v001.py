import numpy as np
from plot_openfoam_airfoil_fields_v002 import coordinate_bijection
def native_mapping(native,reference,depth):
 n=len(reference)
 if native.shape!=(n,3) or not np.isfinite(native).all() or np.max(abs(native[:,2]-depth/2))>1e-9:raise ValueError("Native cell count/finite coordinates/extrusion plane mismatch")
 mapping,metrics=coordinate_bijection(reference,native[:,:2],1e-9)
 return np.argsort(mapping),dict(**metrics,extrusion_depth_m=depth,native_z_plane_m=depth/2)
