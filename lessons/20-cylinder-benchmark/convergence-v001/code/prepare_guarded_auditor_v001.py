from pathlib import Path
S=Path('E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/cylinder-convergence-v001')
s=(S/'audit_cylinder_window_v002.py').read_text()
needle=' result=dict(start_time_s=start'
assert needle in s
s=s.replace(needle," guards=json.loads((out/'runtime-guard-records.json').read_text());assert len(guards)==steps\n for k,g in enumerate(guards,1):\n  assert g['passed'] and g['all_fields_finite_materials_match'] and abs(g['time_s']-(start+k*.005))<1e-8\n  assert g['relative_mass_imbalance']<=1e-4 and max(g['native_terminal_residuals'][n] for n in ['continuity','x-velocity','y-velocity'])<=1e-5\n  for n in ['continuity','x-velocity','y-velocity']:assert math.isclose(g['native_terminal_residuals'][n],rows[k-1][n],rel_tol=5e-5,abs_tol=1e-12)\n hashes[str(out/'runtime-guard-records.json')]=hashlib.sha256((out/'runtime-guard-records.json').read_bytes()).hexdigest()\n result=dict(runtime_guard_records=steps,all_runtime_guards_pass=True,start_time_s=start")
compile(s,'guardedauditor','exec');(S/'audit_cylinder_guarded_window_v001.py').write_text(s)