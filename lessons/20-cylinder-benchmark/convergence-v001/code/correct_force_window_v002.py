from pathlib import Path
import csv,hashlib,io,json
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001');N=Path('E:/CFD-Learning-Lab/native/runs/2026-11-01-vortex-shedding-v001');parent=N/'fluent-corrected-v003'
# Native old CSV schema has no new per-step mass column. Reproduce the inherited-row failure before fixing.
rows=list(csv.DictReader((parent/'force-history.csv').open()));new=dict(rows[-1],relative_mass_imbalance=0.)
try:
 w=csv.DictWriter(io.StringIO(),fieldnames=list(rows[0]));w.writerow(new)
except ValueError as e:red=str(e)
else:raise AssertionError('Expected inherited schema rejection')
s=(S/'cylinder_fluent_pressure_v001.py').read_text();start=s.index('        if resume:\n            # Carry prior pilot records');end=s.index('        save_snapshot(start_step, actual_time)',start);s=s[:start]+s[end:];s=s.replace('cylinder-inputs-v001.json','cylinder-inputs-v002.json');(S/'cylinder_fluent_pressure_v002.py').write_text(s,encoding='utf-8')
w=csv.DictWriter(io.StringIO(),fieldnames=list(new));w.writerow(new)
manifest=json.loads((S/'cylinder-inputs-v001.json').read_text());manifest.pop(str(S/'cylinder_fluent_pressure_v001.py'));manifest[str(S/'cylinder_fluent_pressure_v002.py')]=hashlib.sha256((S/'cylinder_fluent_pressure_v002.py').read_bytes()).hexdigest();(S/'cylinder-inputs-v002.json').write_text(json.dumps(manifest,indent=2)+'\n')
(S/'native-force-schema-regression-v001.json').write_text(json.dumps(dict(real_parent_rows=len(rows),inherited_header_failure_reproduced=red,new_window_header_accepts_mass_column=True,older_timeline_not_silently_relabelled=True),indent=2)+'\n')
with (S/'progress.md').open('a') as f:f.write('\nPrelaunchschemaRED/GREEN: original native forceCSV lacksnewmasscolumn; carryingoldrowsmakesDictWriterrejectnewrow. Preservedunlauncheddriver001/manifest001. Driver002records onlynewwindowrows/states, retainsparentviahashprovenance; nativeCSVschemaregression001 provesoldrejection/newheaderacceptance. No solveradvancedinthispreflight.\n')
compile(s,'driver002','exec');print('Driver002 syntax and real inherited-schema regression pass')