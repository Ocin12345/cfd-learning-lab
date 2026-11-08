"""Immutable-source gate for the predeclared warm continuation."""
from pathlib import Path
import json,sys
from dam_break_admission_v004 import P,project_path,sha,verify
verify('openfoam-production')
receipt=P/'openfoam-scheduled-inputs-v001.json';a=json.loads(receipt.read_text())
for name,h in a['source_sha256'].items():assert sha(project_path(name))==h,name
parent=json.loads((P/'openfoam-production-dt005ms-audit.json').read_text())
assert all(r['passed'] for r in parent['step_records'][:5000])
assert all(r['fields_pass'] for r in parent['records'] if r['time_s']<=.25+1e-12)
fa=json.loads((P/'fluent-production-dt005ms-audit.json').read_text());assert fa['numerical_checks_pass'] and fa['production_complete']
window=project_path(a['damping_window_audit']);w=json.loads(window.read_text());assert w['window_checks_pass'] and w['steps']==100
for record in [parent,fa,w]:
 for name,h in record['native_source_sha256'].items():assert sha(project_path(name))==h,name
print('Warm-continuation source, passed prefix and damping-window admission PASS')