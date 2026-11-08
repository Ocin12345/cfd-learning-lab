"""Source-native column profiles, connected projected-normal numerical widths."""
import numpy as np
from shock_width_v001 import connected_width

def measure_columns(p,geo,metrics):
 p=np.asarray(p);ny,nx=p.shape
 xy=np.asarray(geo['centres']).reshape(ny,nx,2);x,y=xy[...,0],xy[...,1]
 assert np.asarray(geo['y_nodes']).shape==(ny+1,nx+1)
 slope=metrics['shock_fit_slope'];b=metrics['shock_fit_intercept_m']
 p1=metrics['pressure_upstream_pa'];p2=metrics['pressure_downstream_pa']
 if not p2>p1>0:raise ValueError('Positive compression plateau values required')
 d=(slope*x+b-y)/np.sqrt(1+slope*slope);q=(p-p1)/(p2-p1)
 selected=np.flatnonzero((x.mean(axis=0)>=.36)&(x.mean(axis=0)<=.84))
 if not len(selected):raise ValueError('No prescribed columns')
 rows=[]
 for j in selected:
  # Retain the complete column. No smoothing or skipped invalid profiles.
  try:r=connected_width(d[:,j],q[:,j])
  except ValueError as e:raise ValueError(f'column {j}: {e}') from e
  samples=[]
  for name in ('lower_bracket_indices','mid_bracket_indices','upper_bracket_indices'):
   for i in r[name]:samples.append(dict(bracket=name,row=i,common_cell=int(i*nx+j),x_m=float(x[i,j]),y_m=float(y[i,j]),signed_normal_m=float(d[i,j]),p_pa=float(p[i,j]),normalized_jump=float(q[i,j])))
  rows.append(dict(column=int(j),x_mean_m=float(x[:,j].mean()),**r,bracket_samples=samples))
 widths=np.array([r['width_m'] for r in rows]);intervals=np.array([r['width_in_projected_column_intervals'] for r in rows])
 return dict(prescribed_x_band_m=[.36,.84],prescribed_columns=len(selected),valid_columns=len(rows),full_coverage=True,pressure_plateaus_pa=[p1,p2],fit_slope=float(slope),fit_intercept_m=float(b),median_width_m=float(np.median(widths)),p05_width_m=float(np.percentile(widths,5)),p95_width_m=float(np.percentile(widths,95)),median_projected_column_intervals=float(np.median(intervals)),cell_unit_definition='Width / median normal-projected vertical-column centroid increment near threshold brackets; not isotropic cell size',interpolation='Linear between actual native cell centres, connected 10 and 90 percent brackets around strongest 50 percent compression crossing, no smoothing/extrapolation',physical_shock_thickness=False,columns=rows)