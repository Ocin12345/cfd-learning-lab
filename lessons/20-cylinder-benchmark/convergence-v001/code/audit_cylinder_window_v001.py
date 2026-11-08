from pathlib import Path
import csv,hashlib,json,math
from analyse_vortex_fluent_v003 import terminal_rows,mass_values,hdf_fields,hdf_times

def verify_step_rows(rows,start,steps,dt=.005):
 if len(rows)!=steps:raise ValueError('Missing/extra native step records')
 first=round(start/dt)
 for k,r in enumerate(rows,1):
  if r['step']!=first+k or abs(r['time_s']-(start+k*dt))>1e-8:raise ValueError('Native time/step mismatch')
 return rows

def audit_window(out,start=4.,steps=10):
 out=Path(out);summary=json.loads((out/'summary.json').read_text());status=json.loads((out/'status.json').read_text());assert status['stage']=='completed' and not summary['export_errors']
 assert abs(summary['start_time_s']-start)<1e-8 and abs(summary['end_time_s']-(start+steps*.005))<1e-8
 targets=json.loads((out/'declared-targets.json').read_text());assert targets['residual_target']==1e-5 and targets['relative_mass_imbalance_target']==1e-4 and targets['max_inner_iterations_per_time_step']==120
 rows=verify_step_rows(terminal_rows((out/'solver-transcript.txt').read_text()),start,steps)
 misses=[r for r in rows if max(r[n] for n in ['continuity','x-velocity','y-velocity'])>1e-5]
 states=json.loads((out/'snapshots.json').read_text());assert len(states)==steps+1
 checked=[]
 for k,s in enumerate(states):
  t=start+k*.005;assert abs(s['time_s']-t)<1e-8
  path=out/s['data'];fields=hdf_fields(path,21598);clocks=[r for r in hdf_times(path) if r['current']];assert clocks and all(abs(r['value']-t)<1e-8 for r in clocks)
  mass=mass_values((out/'force-reports'/('mass-%06d.txt'%s['step'])).read_text());assert mass['relative_imbalance']<=1e-4
  checked.append(dict(time_s=t,fields=fields,relative_mass_imbalance=mass['relative_imbalance']))
 forces=list(csv.DictReader((out/'force-history.csv').open()));assert len(forces)==steps+1
 for k,r in enumerate(forces):
  assert abs(float(r['time_s'])-(start+k*.005))<1e-8 and all(math.isfinite(float(v)) for v in r.values())
  for coef,key in [('Cd','drag_force_N'),('Cl','lift_force_N')]:assert math.isclose(float(r[coef]),float(r[key])/.05,rel_tol=2e-4,abs_tol=1e-8)
 controls=json.loads((out/'pressure-controls.json').read_text());quality=json.loads((out/'geometry-quality.json').read_text());assert quality['native_numeric_gate_pass']
 files=[p for p in out.iterdir() if p.is_file() and (p.suffix=='.h5' or p.name in ['solver-transcript.txt','executed-script.py','force-history.csv','pressure-controls.json','declared-targets.json','immediately-prelaunch-resources.json','geometry-quality.json'])]+list((out/'snapshots').glob('*.h5'))+list((out/'force-reports').glob('*.txt'))
 hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
 result=dict(start_time_s=start,end_time_s=start+steps*.005,steps=steps,saved_states=len(states),residual_miss_steps=len(misses),residual_misses=misses,step_records=rows,fields_and_mass=checked,maximum_relative_mass_imbalance=max(r['relative_mass_imbalance'] for r in checked),pressure_termination=controls['actual'],original_pressure_termination=controls['original'],numerical_window_pass=not misses,diagnostic_only=True,inherited_startup_not_repaired=True,physical_validation=False,benchmark_pass=False,native_source_sha256=hashes)
 (out/'window-audit.json').write_text(json.dumps(result,indent=2)+'\n');return result