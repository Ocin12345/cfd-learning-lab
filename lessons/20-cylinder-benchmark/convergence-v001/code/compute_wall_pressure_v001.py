from pathlib import Path
import csv,datetime,hashlib,json
import numpy as np
from cylinder_wall_pressure_v001 import fluent_wall_geometry,openfoam_wall_geometry,wall_pressures
from wing_fluent_v001 import available_memory
S=Path(__file__).parent;O=S/'analysis-pressure-v001';assert not O.exists();assert available_memory()>=4*2**30;O.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
states={};audits={};geoms={};metrics={};sources={};geometry_record={}
for solver in ['fluent','openfoam']:
 folder=S/('analysis-'+solver+'-v001');states[solver]=json.loads((folder/(solver+'-snapshots.json')).read_text());audits[solver]=json.loads((folder/('combined-'+solver+'-audit.json')).read_text());assert audits[solver]['all_new_steps_pass'];metrics[solver]=json.loads((folder/(solver+'-cycle-metrics.json')).read_text())
 geoms[solver]=(fluent_wall_geometry if solver=='fluent' else openfoam_wall_geometry)(Path(states[solver][0]['case']));geometry_record[solver]=dict(stencils=geoms[solver]['stencils'],meaning=geoms[solver]['meaning'],native_face_indices_zero_based=geoms[solver]['indices'].tolist())
for name in ['front','rear']:np.testing.assert_allclose(geoms['fluent']['stencils'][name]['face_centres_m'],geoms['openfoam']['stencils'][name]['face_centres_m'],rtol=0,atol=1e-9)
summary={}
for solver in ['fluent','openfoam']:
 rows=[]
 for r in states[solver]:
  path=Path(r['data']);consumed=path if solver=='fluent' else path/'p';expected=audits[solver]['source_sha256'][str(consumed)];assert sha(consumed)==expected;v=wall_pressures(solver,path,geoms[solver]);assert sha(consumed)==expected;sources[str(consumed)]=expected;rows.append(dict(time_s=r['time_s'],front_wall_estimate_Pa=v['front'],rear_wall_estimate_Pa=v['rear'],difference_Pa=v['difference_Pa']))
 assert len(rows)==801
 with (O/(solver+'-wall-pressure-estimates.csv')).open('x',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 t=np.array([r['time_s'] for r in rows]);y=np.array([r['difference_Pa'] for r in rows]);cycles=[]
 for c in metrics[solver]['last_six_cycles']:
  a,b=c['start_time_s'],c['end_time_s'];inside=(t>a)&(t<b);tx=np.r_[a,t[inside],b];v=np.interp(tx,t,y);cycles.append(dict(start_time_s=a,end_time_s=b,time_mean_difference_Pa=float(np.trapezoid(v,tx)/(b-a)),sampled_min_difference_Pa=float(y[inside].min()),sampled_max_difference_Pa=float(y[inside].max())))
 summary[solver]=dict(last_six_cycles=cycles,last_six_mean_difference_Pa=float(np.mean([c['time_mean_difference_Pa'] for c in cycles])),complete_trace_min_Pa=float(y.min()),complete_trace_max_Pa=float(y.max()))
record=dict(created_at=datetime.datetime.now().astimezone().isoformat(),method='Linear interpolation in polar angle of two native evaluated wall-face pressures bracketing each specified cylinder vertex (.15,.2) and(.25,.2). Fluent uses face SV_P; OpenFOAM zeroGradient face value is native owner-cell p. rho=1 converts kinematic p numerically to Pa. These are discretized wall-vertex estimates, not exact continuum point pressures.',cycle_mean_method='Trapezoidal integration with linearly interpolated temporal values only at cycle boundaries; all interior trace samples come from actual saved .005s states.',geometry=geometry_record,matching_native_wall_stencils=True,summary=summary,source_sha256=sources,code_sha256={str(S/n):sha(S/n) for n in ['cylinder_wall_pressure_v001.py','compute_wall_pressure_v001.py','test_wall_pressure_v001.py']},primary_sources=['https://github.com/OpenFOAM/OpenFOAM-14/blob/master/src/finiteVolume/fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.C','https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark2_re100.html'],physical_validation=False,benchmark_pass=False,mesh_independence=False)
(O/'wall-pressure-audit.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(summary,indent=2))