from pathlib import Path
import hashlib,json
from analyse_vortex_fluent_v003 import hdf_times
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');p=Path('E:/CFD-Learning-Lab/native/runs/2026-11-01-vortex-shedding-v001/fluent-corrected-v003/vortex-final.dat.h5');rows=hdf_times(p)
try:old=[r for r in rows if r['current']]
except KeyError as e:red=str(e)
else:raise AssertionError('Expected schema mismatch RED')
current=[r for r in rows if r['scope']=='current'];assert current and all(abs(r['flow_time_s']-4)<1e-8 for r in current)
s=(S/'audit_cylinder_window_v001.py').read_text().replace("r['current']","r['scope']=='current'").replace("r['value']-t","r['flow_time_s']-t");compile(s,'auditor002','exec');(S/'audit_cylinder_window_v002.py').write_text(s)
s=(S/'run_cylinder_pressure_pair_v001.py').read_text().replace('audit_cylinder_window_v001','audit_cylinder_window_v002');(S/'run_cylinder_pressure_pair_v002.py').write_text(s)
(S/'native-clock-schema-regression-v001.json').write_text(json.dumps(dict(actual_native_records=len(rows),old_accessor_red=red,new_current_time_records=current,passed=True),indent=2)+'\n')
files=json.loads((S/'cylinder-inputs-v002.json').read_text())
for name in ['audit_cylinder_window_v002.py','run_cylinder_pressure_pair_v002.py','test_native_coverage_v001.py','native-clock-schema-regression-v001.json','native-force-schema-regression-v001.json','guard-green-v001.log','coverage-green-v001.log']:
 p=S/name;files[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
(S/'diagnostic-controller-inputs-v001.json').write_text(json.dumps(files,indent=2)+'\n')
with (S/'progress.md').open('a') as f:f.write('\nPrelaunchnativeclockRED/GREEN: helperuses scope/flow_time_s, not current/value. RealparentHDF reproducesKeyError; auditor002usesactualschema andcurrent4s verified. Preservedunlaunchedauditor001/controller001. Controller002/frozencontrollerinputs001selected.\n')
print('Native clock schema verified; corrected frozen controller ready')