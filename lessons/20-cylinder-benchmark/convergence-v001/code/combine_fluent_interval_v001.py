from pathlib import Path
import csv,hashlib,json,math
import numpy as np
from plot_vortex_shedding_v002 import fluent_mesh,native_field
from cylinder_cycle_metrics_v001 import cycle_metrics
S=Path(__file__).parent;N=Path('E:/CFD-Learning-Lab/native/runs/showcase-cylinder-convergence-v001');O=S/'analysis-fluent-v001';assert not O.exists();O.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((S/'cycle-development-protocol-v001.json').read_text())
for p,h in protocol['source_sha256'].items():assert sha(Path(p))==h
segments=[('fluent-coupled-pilot-v001',4.,4.05,10),('fluent-coupled-window-v001',4.05,4.55,100),('fluent-developed-v001',4.55,8.,690)]
allforces={};allstates={};allsteps=[];sources={};boundaries=[];mesh0=None
for name,start,end,steps in segments:
 p=N/name;a=json.loads((p/'window-audit.json').read_text());assert a['numerical_window_pass'] and a['steps']==steps and a['saved_states']==steps+1 and abs(a['start_time_s']-start)<1e-8 and abs(a['end_time_s']-end)<1e-8
 if steps>10:assert a['all_runtime_guards_pass'] and a['runtime_guard_records']==steps
 for source,h in a['native_source_sha256'].items():assert sha(Path(source))==h,source;sources[source]=h
 sources[str(p/'window-audit.json')]=sha(p/'window-audit.json');mesh=fluent_mesh(p/'vortex-final.cas.h5')
 if mesh0 is None:mesh0=mesh
 else:assert np.array_equal(mesh['xy'],mesh0['xy']) and np.array_equal(mesh['triangles'],mesh0['triangles']),'Native global cell/mesh order changed across segments'
 for r in a['step_records']:allsteps.append(dict(segment=name,**r))
 for r in json.loads((p/'snapshots.json').read_text()):
  t=round(r['time_s'],8);data=p/r['data'];entry=dict(time_s=t,step=r['step'],data=str(data),case=str(p/'vortex-final.cas.h5'),segment=name)
  if t in allstates:
   old=allstates[t];f0=native_field('fluent',dict(path=Path(old['data'])),mesh);f1=native_field('fluent',dict(path=data),mesh)
   diffs={k:float(np.max(np.abs(f0[k]-f1[k]))) for k in ['u','v','p']};assert max(diffs.values())<1e-12,(t,diffs);boundaries.append(dict(time_s=t,maximum_field_differences=diffs,initial_copy_matches_previous_endpoint=True));continue
  allstates[t]=entry
 for raw in csv.DictReader((p/'force-history.csv').open()):
  r={k:float(v) for k,v in raw.items()};t=round(r['time_s'],8)
  if t in allforces:
   assert all(math.isclose(r[k],allforces[t][k],rel_tol=5e-6,abs_tol=1e-7) for k in ['Cd','Cl']);continue
  allforces[t]=r
assert len(allsteps)==800 and len(allstates)==len(allforces)==801
expected=np.arange(800)*.005+4.005;assert np.allclose([r['time_s'] for r in allsteps],expected,atol=1e-8,rtol=0)
assert np.allclose(sorted(allstates),np.arange(801)*.005+4,atol=1e-8,rtol=0)
forces=[allforces[t] for t in sorted(allforces)];states=[allstates[t] for t in sorted(allstates)]
with (O/'fluent-force-history.csv').open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=list(forces[0]));w.writeheader();w.writerows(forces)
with (O/'fluent-step-residuals.csv').open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=list(allsteps[0]));w.writeheader();w.writerows(allsteps)
(O/'fluent-snapshots.json').write_text(json.dumps(states,indent=2)+'\n')
m=cycle_metrics([r['time_s'] for r in forces],[r['Cd'] for r in forces],[r['Cl'] for r in forces]);(O/'fluent-cycle-metrics.json').write_text(json.dumps(m,indent=2)+'\n')
record=dict(steps=800,unique_native_states=801,start_time_s=4.,end_time_s=8.,all_new_steps_pass=True,native_topology_and_global_cell_order_identical=True,duplicate_boundary_states_checked=boundaries,initial_four_second_state_inherited_not_residually_repaired=True,maximum_relative_mass_imbalance=max(r['relative_mass_imbalance'] for r in forces),cycle_protocol=str(S/'cycle-development-protocol-v001.json'),cycle_protocol_sha256=sha(S/'cycle-development-protocol-v001.json'),source_sha256=sources,analysis_source_sha256=sha(Path(__file__)),physical_validation=False,benchmark_pass=False)
(O/'combined-fluent-audit.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(dict(steps=800,states=801,complete_cycles=len(m['complete_cycles']),last_six_stable=m['last_six_cycles_stable_under_declared_criterion'],last_six_relative_ranges=m.get('last_six_relative_ranges'),period=m.get('last_six_mean_period_s'),mean_Cd=m.get('last_six_mean_Cd'),mass_max=record['maximum_relative_mass_imbalance']),indent=2))