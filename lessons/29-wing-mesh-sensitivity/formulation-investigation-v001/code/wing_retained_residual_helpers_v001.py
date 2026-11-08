import numpy as np
from wing_repair_native_v003 import expected_settings
def expected_probe_settings():
 a=expected_settings(.1,20.,.5);a['controls']['advanced']['expert']['save_cell_residuals']=True;return a
def classify_stage(values,initial,corrected):
 a=np.asarray(values,dtype=float)
 if a.ndim!=1 or not len(a) or not np.isfinite(a).all():raise ValueError('Invalid retained residual array')
 s=float(np.abs(a).sum());matches=[name for name,v in [('assembled_initial',initial),('corrected_final',corrected)] if np.isclose(s,v,rtol=1e-6,atol=1e-12)]
 if len(matches)!=1:raise ValueError('Unknown or ambiguous retained residual stage')
 return matches[0]
