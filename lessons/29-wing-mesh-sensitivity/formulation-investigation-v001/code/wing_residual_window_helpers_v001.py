import numpy as np
from wing_cell_identity_v001 import centroid_bijection
from wing_retained_residual_helpers_v001 import classify_stage
EQ=('continuity','x-velocity','y-velocity','z-velocity')
def declared_iterations(parent,end):
 if parent!=1502 or end!=1510:raise ValueError('Only declared1502->1510 window admitted')
 return list(range(parent+1,end+1))
def verify_runtime(centres,native_u,native_p,arrays,raw,corrected):
 n=len(centres);api=np.asarray(arrays['SV_CENTROID']).reshape(n,3);index,distance=centroid_bijection(centres,api)
 u=np.column_stack([arrays[k] for k in ['SV_U','SV_V','SV_W']]);p=np.asarray(arrays['SV_P'])
 if not np.isfinite(u).all() or not np.isfinite(p).all() or not np.allclose(u,np.asarray(native_u)[index],rtol=1e-10,atol=1e-12) or not np.allclose(p,np.asarray(native_p)[index],rtol=1e-10,atol=1e-12):raise ValueError('Runtime/nativeU,p mismatch')
 values={};sums={}
 for key,eq in zip(['SV_P_R','SV_U_R','SV_V_R','SV_W_R'],EQ):
  v=np.asarray(arrays[key]);
  if v.shape!=(n,) or not np.isfinite(v).all():raise ValueError('Invalid retained residual array')
  sums[key]=float(np.abs(v).sum())
  if not np.isclose(sums[key],raw[eq],rtol=1e-6,atol=1e-12):raise ValueError('Retained sum does not match rawHDF')
  values[key]=np.empty(n);values[key][index]=v
 stage=classify_stage(values['SV_P_R'],raw['continuity'],corrected)
 if stage!='assembled_initial':raise ValueError('Wrong retained residual stage')
 return dict(stage=stage,residuals=values,absolute_sums=sums,maximum_centroid_distance_m=float(distance.max()))
