from pathlib import Path
import numpy as np,json
from plot_nozzle_native_v004 import foam_field
from plume_common_v001 import *
p=ROOT/'native/runs/hot-patch-plume-qualification-v001/openfoam-control/600'
T=foam_field(p/'T',NC,1,[0,0,0,1,0,0,0]);U=foam_field(p/'U',NC,3,[0,1,-1,0,0,0,0]);pres=foam_field(p/'p',NC,1,[1,-1,-2,0,0,0,0]);rho=foam_field(p/'rho',NC,1,[1,-3,0,0,0,0,0])
s=dict(max_speed_m_s=float(np.linalg.norm(U,axis=1).max()),min_T_K=float(T.min()),max_T_K=float(T.max()),mass_per_depth_kg_m=float(rho.sum()*DX*DY),mass_relative_drift=float(abs(rho.mean()/RHO0-1)),eos_relative_error=float(np.max(abs(rho*R*T/pres-1))))
s['checks']=dict(quiet=s['max_speed_m_s']<=1e-5,temperature=bool(np.max(abs(T-T0))<=1e-4),mass=s['mass_relative_drift']<=1e-5,eos=s['eos_relative_error']<=1e-6)
write(p.parent/'control-checks.json',s);print(json.dumps(s,default=str,indent=2));assert all(s['checks'].values())
