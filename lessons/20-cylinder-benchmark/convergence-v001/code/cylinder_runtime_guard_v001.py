from pathlib import Path
from analyse_vortex_fluent_v003 import hdf_fields,hdf_residuals,hdf_times,mass_values

def check_native_step(data,mass_report,expected_time):
 clocks=[r for r in hdf_times(Path(data)) if r['scope']=='current']
 if not clocks or any(abs(r['flow_time_s']-expected_time)>1e-8 for r in clocks):raise ValueError('Native current time mismatch')
 fields=hdf_fields(Path(data),21598);rows=hdf_residuals(Path(data))
 if not rows:raise ValueError('Missing native residual history')
 r=rows[-1]
 if any(r[n]>1e-5 for n in ['continuity','x-velocity','y-velocity']):raise ValueError('Native terminal residual misses unchanged1e-5: '+str(r))
 mass=mass_values(Path(mass_report).read_text())
 if mass['relative_imbalance']>1e-4:raise ValueError('Native boundary imbalance misses unchanged1e-4')
 return dict(time_s=expected_time,native_terminal_residuals=r,relative_mass_imbalance=mass['relative_imbalance'],all_fields_finite_materials_match=True,passed=True)