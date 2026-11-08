import numpy as np
from wing_cell_identity_v001 import centroid_bijection
EQ=('continuity','x-velocity','y-velocity','z-velocity')
def verify_simplec_runtime(centres,u,p,arrays,raw,corrected):
 n=len(centres);index,distance=centroid_bijection(centres,np.asarray(arrays['SV_CENTROID']).reshape(n,3))
 liveu=np.column_stack([arrays[k] for k in ['SV_U','SV_V','SV_W']]);livep=np.asarray(arrays['SV_P'])
 if not np.isfinite(liveu).all() or not np.isfinite(livep).all() or not np.allclose(liveu,u[index],rtol=1e-10,atol=1e-12) or not np.allclose(livep,p[index],rtol=1e-10,atol=1e-12):raise ValueError('Runtime/native field mismatch')
 sums={};matches={}
 for k,e in zip(['SV_P_R','SV_U_R','SV_V_R','SV_W_R'],EQ):
  a=np.asarray(arrays[k])
  if a.shape!=(n,) or not np.isfinite(a).all():raise ValueError('Invalid retained array')
  sums[k]=float(np.abs(a).sum());matches[k]=bool(np.isclose(sums[k],raw[e],rtol=1e-6,atol=1e-12))
 c=bool(np.isclose(sums['SV_P_R'],corrected,rtol=1e-6,atol=1e-12));h=matches['SV_P_R'];stage='matches_both_nonunique' if h and c else 'matches_HDF_numerator_only' if h else 'matches_corrected_native_phi_only' if c else 'unidentified'
 return dict(capture_identity_pass=True,maximum_centroid_distance_m=float(distance.max()),absolute_sums=sums,raw_numerator_matches=matches,pressure_matches_corrected_phi=c,pressure_stage=stage,no_Coupled_linear_block_stage_inferred=True)
