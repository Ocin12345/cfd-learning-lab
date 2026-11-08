"""Verify immutable pressure-cycle diagnostic inputs on either host."""
import json,sys
from pathlib import Path
from dam_break_admission_v004 import verify,project_path,sha,P
verify('openfoam-production')
a=json.loads((P/'openfoam-pressure-cycle-inputs-v002.json').read_text())
for name,h in a['source_sha256'].items():assert sha(project_path(name))==h,'Pressure-cycle source/input changed: '+name
print('Pressure-cycle native admission PASS')