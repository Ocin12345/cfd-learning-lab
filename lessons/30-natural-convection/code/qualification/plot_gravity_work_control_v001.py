from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from vertical_energy_reader_v002 import fluent_native
from cavity_buoyancy_common_v003 import *
out=ROOT/'exports/figures/2026-11-07-gravity-work-control-v001';assert not out.exists();out.mkdir(parents=True)
fig,axs=plt.subplots(1,2,figsize=(11,4.5),layout='constrained');record=[]
for folder,stem,label,color in [('vertical-energy-control-fluent-v002','zero-gravity','Gravity off','#808d9e'),('vertical-energy-control-fluent-v003','with-gravity','Gravity on, original energy model','#cf652f'),('vertical-energy-control-fluent-v005','corrected','Gravity on, explicit work source','#007e83')]:
 p=ROOT/'native/runs'/folder;case=p/(stem+'.cas.h5');data=p/(stem+'.dat.h5');q=fluent_native(case,data,P0)
 xy=q['nodes'][q['cells'],:2].mean(axis=1);ys=np.unique(np.round(xy[:,1],10));T=[];H=[]
 for y in ys:
  idx=abs(xy[:,1]-y)<1e-8;assert idx.sum()==24;T.append(float(q['T'][idx].mean()));H.append(float((CP*q['T'][idx]+.5*(q['U'][idx]**2).sum(axis=1)).mean()))
 slope=float(np.polyfit(ys,H,1)[0]);axs[0].plot(ys,(np.array(T)-T0)*1000,'o-',ms=3,label=label,color=color);axs[1].plot(ys,np.array(H)-H[0],'o-',ms=3,label=label,color=color)
 record.append(dict(label=label,height_m=ys.tolist(),temperature_K=T,enthalpy_plus_kinetic_J_kg=H,slope_J_kg_m=slope,native_residuals=q['solve_evidence'],source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [case,data]}))
axs[1].plot(ys,-G*(ys-ys[0]),'k--',lw=1.5,label='Expected gravitational work')
axs[0].set(xlabel='Height (m)',ylabel='Temperature relative to 300 K (mK)',title='Actual Fluent cell-row means')
axs[1].set(xlabel='Height (m)',ylabel='Change in h + |U|²/2 (J/kg)',title='Energy-source verification')
for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('Check gravitational work before comparing buoyant flow',fontsize=16)
fig.savefig(out/'gravity-work-check.png',dpi=180);fig.savefig(out/'gravity-work-check.pdf');plt.close(fig)
write(out/'evidence.json',dict(control='adiabatic vertical slip duct, inlet speed 1 m/s, inlet temperature 300 K, 1 m x 10 m; diagnostic only',records=record,corrected_slope_relative_error=abs(record[-1]['slope_J_kg_m']/G+1),source_expression='-9.81 [m s^-2] * Density * Velocity.y',physical_validation=False,source_scope='tested pressure-based ideal-gas laminar settings only'))
print('Energy diagnostic plotted; fine OpenFOAM control prepared')
