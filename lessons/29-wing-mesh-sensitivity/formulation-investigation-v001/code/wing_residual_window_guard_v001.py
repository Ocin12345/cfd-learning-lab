import h5py,numpy as np,json
from wing_repair_common_v001 import *
from wing_repair_native_v003 import raw_residual
from audit_wing_continuation_v001 import native
from plot_finite_wing_native_v002 import indexed
from wing_repair_helpers_v001 import same_residual_prefix
from wing_residual_window_helpers_v001 import verify_runtime

def checkpoint_guard(out,end):
 capture=out/'runtime'/str(end);case=out/f'iter-{end}.cas.h5';data=out/f'iter-{end}.dat.h5';result,_,_,g=native(case,data);assert result['iteration']==end and result['cells']==224972
 assert max(result['final_residuals'].values())<.1 and result['relative_mass_imbalance']<=1e-4
 n=result['cells'];pi,pr=raw_residual(out/'source.dat.h5');it,raw=raw_residual(data);same_residual_prefix(pi,pr,it,raw,1502)
 with h5py.File(data) as h:
  cells=h['results/1/phase-1/cells'];u=np.column_stack([indexed(cells[k],n) for k in ['SV_U','SV_V','SV_W']]);p=indexed(cells['SV_P'],n)
 with np.load(capture/'runtime-residual-arrays.npz') as a:proof=verify_runtime(g['points'][g['tets']].mean(axis=1),u,p,a,{k:float(raw[k][-1,0]) for k in EQ},result['cell_mass_imbalance_abs_sum_kg_s'])
 record=dict(iteration=end,stage=proof['stage'],maximum_centroid_distance_m=proof['maximum_centroid_distance_m'],absolute_sums=proof['absolute_sums'],native=result,capture_pass=True,original_residual_pass=max(result['final_residuals'].values())<=1e-5,files_sha256={str(p):sha(p) for p in [case,data,capture/'runtime-residual-arrays.npz',capture/'runtime-residual-catalogue.json']})
 dump(capture/'native-runtime-guard.json',record);return record
