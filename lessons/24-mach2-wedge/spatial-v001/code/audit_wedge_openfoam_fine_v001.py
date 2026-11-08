from pathlib import Path
import hashlib,json
import numpy as np
from wedge_openfoam_fine_io_v001 import read_openfoam
from plot_wedge_fluent_fine_v001 import measure
from wedge_common_fine_v001 import LIMITS
S=Path(__file__).parent;R=Path('E:/CFD-Learning-Lab')
def audit(case):
 with np.load(R/'native/meshes/showcase-mach2-wedge-spatial-v001/fine-v001/geometry.npz') as geo:
  times,frames,masses,native=read_openfoam(case,geo);m=measure(frames[-1],frames[-5],geo)
 checks=dict(mesh_quality=True,finite_positive_frames=True,complete_saved_times=len(times)==33,field_stability=bool(m['field_stability_check']),mass_balance=native['mass_imbalance']<=LIMITS['mass_imbalance'],late_five_mass=native['last_five_frames_max_mass_imbalance']<=LIMITS['mass_imbalance'],supersonic_outlet=m['outlet_minimum_cell_mach']>1,ideal_gas=m['maximum_ideal_gas_relative_error']<1e-5,analytic_angle=bool(m['analytical_angle_check']),analytic_pressure=bool(m['analytical_pressure_check']),column_coverage=m['shock_fit_columns']>=160)
 paths=[p for p in case.rglob('*') if p.is_file() and 'processor' not in p.relative_to(case).parts[0]]
 record=dict(source_case=str(case),cells=89600,native=native,measurements=m,check_results=checks,failed_checks=[k for k,v in checks.items() if not v],numerical_endpoint_pass=all(checks.values()),spatial_and_adaptive_time_changes_confounded=True,physical_validation=False,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths})
 p=case/'endpoint-audit.json';assert not p.exists();p.write_text(json.dumps(record,indent=2)+'\n');return record
