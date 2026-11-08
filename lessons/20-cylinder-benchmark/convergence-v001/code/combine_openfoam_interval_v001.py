from pathlib import Path
import csv,hashlib,json,math
import numpy as np
from plot_vortex_shedding_v002 import openfoam_mesh,native_field
from cylinder_cycle_metrics_v001 import cycle_metrics
S=Path(__file__).parent;N=Path('E:/CFD-Learning-Lab/native/runs/showcase-cylinder-convergence-v001');O=S/'analysis-openfoam-v001';assert not O.exists();O.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((S/'cycle-development-protocol-v001.json').read_text())
for p,h in protocol['source_sha256'].items():assert sha(Path(p))==h
segments=[('openfoam-qualification-v001',4.,4.5,100),('openfoam-developed-v001',4.5,8.,700)];forces={};states={};steps=[];sources={};boundaries=[];courants=[];massmax=0.;mesh0=None
for name,start,end,count in segments:
 p=N/name;a=json.loads((p/'window-audit.json').read_text());assert a['numerical_window_pass'] and a['steps']==count and a['saved_states']==count+1 and a['residual_miss_steps']==0
 for n,h in a['native_source_sha256'].items():assert sha(Path(n))==h,n;sources[n]=h
 sources[str(p/'window-audit.json')]=sha(p/'window-audit.json');mesh=openfoam_mesh(p)
 if mesh0 is None:mesh0=mesh
 else:assert np.array_equal(mesh['xy'],mesh0['xy']) and np.array_equal(mesh['triangles'],mesh0['triangles'])
 steps += [dict(segment=name,**r) for r in a['step_records']];courants +=[dict(segment=name,**r) for r in a['courant_records']];massmax=max(massmax,a['maximum_relative_mass_imbalance'])
 for r in a['fields_and_mass']:
  t=round(r['time_s'],8);folder=p/f'{t:g}';assert folder.is_dir();entry=dict(time_s=t,step=round(t/.005),data=str(folder),case=str(p),segment=name)
  if t in states:
   previous=Path(states[t]['data']);identical={n:sha(previous/n)==sha(folder/n) for n in ['U','p','phi']};assert all(identical.values());boundaries.append(dict(time_s=t,scientific_fields_byte_identical=identical));continue
  states[t]=entry
 for raw in csv.DictReader((p/'force-history.csv').open()):
  r={k:float(v) for k,v in raw.items()};t=round(r['time_s'],8)
  if t in forces:assert all(math.isclose(r[k],forces[t][k],rel_tol=1e-10,abs_tol=1e-11) for k in ['Cd','Cl']);continue
  forces[t]=r
assert len(steps)==800 and len(states)==len(forces)==801 and len(courants)==802
np.testing.assert_allclose([r['time_s'] for r in steps],4.005+np.arange(800)*.005,rtol=0,atol=1e-8);np.testing.assert_allclose(sorted(states),4+np.arange(801)*.005,rtol=0,atol=1e-8)
fr=[forces[t] for t in sorted(forces)];st=[states[t] for t in sorted(states)];m=cycle_metrics([r['time_s'] for r in fr],[r['Cd'] for r in fr],[r['Cl'] for r in fr])
with (O/'openfoam-force-history.csv').open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=list(fr[0]));w.writeheader();w.writerows(fr)
(O/'openfoam-snapshots.json').write_text(json.dumps(st,indent=2)+'\n');(O/'openfoam-step-records.json').write_text(json.dumps(steps,indent=2)+'\n');(O/'openfoam-courant-records.json').write_text(json.dumps(courants,indent=2)+'\n');(O/'openfoam-cycle-metrics.json').write_text(json.dumps(m,indent=2)+'\n')
record=dict(steps=800,unique_native_states=801,start_time_s=4.,end_time_s=8.,all_new_steps_pass=True,native_topology_and_global_cell_order_identical=True,duplicate_boundary_states_checked=boundaries,initial_four_second_state_inherited_not_residually_repaired=True,maximum_relative_mass_imbalance=massmax,courant_reports=802,courant_unique_pre_advance_state_times=len(set(round(r['native_pre_advance_state_time_s'],8) for r in courants)),courant_time_convention='Actual native pre-advance state times, including both initial reports at each restart; not interpolated state times or an imposed CFL acceptance target.',cycle_protocol_sha256=sha(S/'cycle-development-protocol-v001.json'),source_sha256=sources,analysis_source_sha256=sha(Path(__file__)),physical_validation=False,benchmark_pass=False)
(O/'combined-openfoam-audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(dict(steps=800,states=801,complete_cycles=len(m['complete_cycles']),last_six_stable=m['last_six_cycles_stable_under_declared_criterion'],last_six_relative_ranges=m.get('last_six_relative_ranges'),period=m.get('last_six_mean_period_s'),mean_Cd=m.get('last_six_mean_Cd'),mass_max=massmax),indent=2))