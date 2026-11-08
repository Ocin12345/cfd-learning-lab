"""Check all fresh half-step native copies and exact initial fields."""
from pathlib import Path
import hashlib,json,re
import numpy as np
from dam_break_common_v002 import *
R=Path('E:/CFD-Learning-Lab');P=R/'local/dam-break-over-obstacle-v001';N=R/'native/runs/dam-break-over-obstacle-v001'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 g=geometry();order=g['openfoam_to_common'];records=[];sources={}
 for mode,rel in [('flat','qualification/openfoam-flat-dt005ms-v001'),('collapse','qualification/openfoam-collapse-dt005ms-v001'),('production','production/openfoam-dt005ms-v001')]:
  case=N/rel;mesh=(case/'mesh-check.log').read_text();assert 'Mesh OK.' in mesh and re.search(r'\bcells:\s+9072\b',mesh)
  c=native_field(case/'0/C',3);v=native_field(case/'0/Vc',1).ravel();a=native_field(case/'0/alpha.water',1).ravel();u=native_field(case/'0/U',3);pr=native_field(case/'0/p_rgh',1).ravel();expected=initial_fields(g,mode)
  assert np.max(abs(c[:,:2]-g['centres'][order]))<1e-12 and np.max(abs(c[:,2]-DEPTH/2))<1e-12 and np.max(abs(v/DEPTH-g['areas'][order]))<1e-15
  assert np.max(abs(a-expected['alpha'][order]))<1e-14 and np.max(abs(pr-expected['p_rgh'][order]))<1e-10 and np.max(abs(u))==0
  control=(case/'system/controlDict').read_text();declared=json.loads((case/'declared-settings.json').read_text())
  steps=10000 if mode=='production' else 400;interval=100 if mode=='production' else 20
  assert declared['dt_s']==DT==.00005 and declared['steps']==steps and declared['save_interval']==interval
  assert 'deltaT 5e-05;' in control and f'writeInterval {interval};' in control and 'adjustTimeStep no;' in control
  for name in ['system/fvSchemes','system/fvSolution','system/decomposeParDict']:
   assert sha(case/name)==sha(N/'production/openfoam-v002'/name),name
  for path in case.rglob('*'):
   if path.is_file():sources[str(path)]=sha(path)
  records.append(dict(mode=mode,cells=NCELLS,native_mesh_pass=True,native_mapping_and_cell_areas_pass=True,exact_initial_fields_pass=True,dt_s=DT,steps=steps,water_area_m2=float(a@v/DEPTH)))
 out=P/'openfoam-native-input-check-v003.json';assert not out.exists();out.write_text(json.dumps(dict(checks=records,native_source_sha256=sources,physical_validation=False,solvers_launched=False),indent=2)+'\n');print(json.dumps(records,indent=2))
if __name__=='__main__':main()