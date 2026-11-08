"""Keep native pre-advance Courant context; separately scale by new/old step size."""
import numpy as np

def record_courant(raw,dt,times,initial_dt):
 raw,dt,times=map(lambda a:np.asarray(a,float),(raw,dt,times))
 if len(raw)!=len(times)+1 or len(dt)!=len(times) or len(dt)<1 or not all(np.isfinite(a).all() for a in (raw,dt,times)) or np.any(dt<=0) or np.any(np.diff(times)<=0):raise ValueError('Complete monotonic native Courant/deltaT/Time records required')
 if abs(raw[0]-raw[1])>1e-8 or raw.max()>.306:raise ValueError('Original raw Courant limit or duplicate startup mismatch')
 if np.any(raw[2:]>.153):raise ValueError('Reduced-CFL established native records exceed unchanged 1.02 margin')
 old=np.r_[initial_dt,dt[:-1]];scaled=raw[1:]*dt/old
 if np.any(scaled>.153):raise ValueError('An actual advance exceeds the reduced-CFL admission')
 records=[dict(role='constructor, inherited checkpoint deltaT',pre_advance_time_s=.007,native_raw_courant=float(raw[0]))]
 for i in range(len(dt)):records.append(dict(role='pre-advance',pre_advance_time_s=.007 if i==0 else float(times[i-1]),end_time_s=float(times[i]),old_deltaT_s=float(old[i]),new_deltaT_s=float(dt[i]),native_raw_courant=float(raw[i+1]),advance_scaled_courant=float(scaled[i])))
 return dict(raw_startup_courant=raw[:2].tolist(),maximum_raw_courant=float(raw.max()),maximum_established_raw_courant=float(raw[2:].max()) if len(raw)>2 else None,maximum_advance_scaled_courant=float(scaled.max()),interpretation='Native Co uses current pre-advance state and inherited/preceding deltaT; scaled Co multiplies raw by newly selected deltaT / old deltaT at that state. No timestamp shift or removal of startup evidence.',records=records)