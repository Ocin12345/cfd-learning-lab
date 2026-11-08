from pathlib import Path
import math,re
import h5py,numpy as np
from wedge_common_fine_v001 import NX,NY,RESIDUALS,LIMITS
from plot_wedge_fluent_fine_v001 import fields,measure
from analyse_nozzle_fluent_v005 import common_mapping,chunk_array
def accepted(residual,metrics,mass):
 if set(residual)!=set(RESIDUALS) or not all(math.isfinite(float(v)) and v>=0 for v in residual.values()):raise ValueError("Four finite native scaled residuals required")
 return dict(residuals=all(residual[k]<=v for k,v in RESIDUALS.items()),mass_balance=math.isfinite(mass) and mass<=LIMITS['mass_imbalance'],field_stability=max(metrics['relative_field_change_l2'].values())<=LIMITS['field_stability_relative'],supersonic_outlet=metrics['outlet_minimum_cell_mach']>1,ideal_gas=metrics['maximum_ideal_gas_relative_error']<1e-5,analytic_angle=bool(metrics['analytical_angle_check']),analytic_pressure=bool(metrics['analytical_pressure_check']),column_coverage=metrics['shock_fit_columns']>=160)
def load(case,data,geometry):
 mapping=common_mapping(case,geometry)
 with h5py.File(data) as h:
  cells=h['results/1/phase-1/cells'];out={}
  for key,native in [('p','SV_P'),('rho','SV_DENSITY'),('T','SV_T'),('u','SV_U'),('v','SV_V')]:
   a=np.empty(NX*NY);a[mapping]=chunk_array(cells[native],NX*NY);out[key]=a.reshape(NY,NX)
 return fields(out)
def evaluate(run,checkpoint,previous,geometry):
 now=load(run/checkpoint['case'],run/checkpoint['data'],geometry);old=load(run/previous['case'],run/previous['data'],geometry);metrics=measure(now,old,geometry)
 text=(run/checkpoint['mass_report']).read_text();rows={n:float(re.search(r'^\s*'+n+r'\s+([-+0-9.eE]+)\s*$',text,re.M)[1]) for n in ('inlet','outlet','ramp','top','Net')}
 assert all(math.isfinite(v) for v in rows.values()) and rows['inlet']>0>rows['outlet'];assert abs(rows['Net']-sum(v for k,v in rows.items() if k!='Net'))<1e-4
 imbalance=abs(rows['Net']/rows['inlet']);checks=accepted(checkpoint['residuals'],metrics,imbalance)
 return dict(check_results=checks,numerical_endpoint_pass=all(checks.values()),measurements=metrics,native_mass=rows,mass_imbalance=imbalance,field_values_positive_finite=True,actual_iteration=checkpoint['iteration'],previous_iteration=previous['iteration'],physical_validation=False)
