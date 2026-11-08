from pathlib import Path
import hashlib,json,re
import numpy as np
from analyse_nozzle_openfoam_v002 import field,FIELDS,patch_values
from wedge_foam_mapping_v002 import native_mapping
from plot_wedge_fluent_fine_v001 import fields,measure
from wedge_width_columns_v001 import measure_columns
from wedge_late_courant_v001 import record_courant
from wedge_common_fine_v001 import DEPTH,NX,NY,LIMITS
R=Path('E:/CFD-Learning-Lab');N=R/'native/runs/showcase-mach2-wedge-spatial-v001';G=R/'native/meshes/showcase-mach2-wedge-spatial-v001/fine-v001/geometry.npz';NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def audit(o):
 assert (o/'run-status.txt').read_text().strip()=='complete'
 log=(o/'solver.log').read_text();assert re.search(r'^End\s*$',log,re.M) and 'Finalising parallel run' in log and 'FOAM FATAL' not in log
 ts=np.array([float(t) for t in re.findall(r'^Time = ('+NUMBER+r')s?\s*$',log,re.M)]);assert abs(ts[-1]-.008)<1e-10 and ts[0]>.007
 dt=np.array([float(t) for t in re.findall(r'^deltaT = ('+NUMBER+r')\s*$',log,re.M)]);assert len(ts)==len(dt);assert max(abs(np.diff(np.r_[.007,ts])-dt))<2e-11
 co=[float(v) for v in re.findall(r'^Courant Number mean: '+NUMBER+r' max: ('+NUMBER+r')',log,re.M)]
 initial=(o/'0.007/uniform/time').read_text();old_dt=float(re.search(r'\bdeltaT\s+('+NUMBER+');',initial)[1]);cfl=record_courant(co,dt,ts,old_dt)
 for part in ('processor0','processor1'):
  t=(o/part/'0.008/uniform/time').read_text();assert re.search(r'name\s+"0.008";',t),part
 with np.load(G) as geo:
  C=field(o/'0/C',NX*NY,*FIELDS['C']);order,mapping=native_mapping(C,geo['centres'],DEPTH)
  def load(folder):
   arrays={k:field(folder/k,NX*NY,*FIELDS[k])[order].reshape(NY,NX) for k in ('p','rho','T')};v=field(folder/'U',NX*NY,*FIELDS['U'])[order].reshape(NY,NX,3);assert abs(v[...,2]).max()<1e-9
   return fields(dict(**arrays,u=v[...,0],v=v[...,1]))
  times=[.007,.00725,.0075,.00775,.008];frames=[];mass=[]
  for t in times:
   folder=o/f'{t:.9g}';frames.append(load(folder));flux={name:float(patch_values(folder/'phi',name,n).sum())/DEPTH for name,n in [('inlet',NY),('outlet',NY),('top',NX),('ramp',NX)]};assert flux['inlet']<0<flux['outlet'];mass.append(dict(time_s=t,**flux,relative_imbalance=abs(sum(flux.values())/flux['inlet'])))
  metrics=measure(frames[-1],frames[0],geo);width=measure_columns(frames[-1]['p'],geo,metrics)
  parent=load(N/'openfoam-fine-v001/0.008');areas=geo['areas'].reshape(NY,NX)
  difference={k:float(np.sqrt(np.sum(areas*(frames[-1][k]-parent[k])**2)/np.sum(areas*parent[k]**2))) for k in ('p','mach')}
 checks=dict(complete_native_times_and_partition_clocks=True,exact_checkpoint=True,finite_positive_native_fields=True,mesh_mapping=True,courant=True,mass_balance=max(r['relative_imbalance'] for r in mass)<=LIMITS['mass_imbalance'],field_stability=metrics['field_stability_check'],supersonic_outlet=metrics['outlet_minimum_cell_mach']>1,ideal_gas=metrics['maximum_ideal_gas_relative_error']<1e-5,analytic_angle=metrics['analytical_angle_check'],analytic_pressure=metrics['analytical_pressure_check'],column_coverage=width['valid_columns']==192)
 failed=[k for k,v in checks.items() if not v]
 result=dict(native_case=str(o),saved_times_s=times,time_steps=len(ts),initial_deltaT_s=old_dt,courant=cfl,mass=mass,measurements=metrics,width=width,fine_main_to_late_endpoint_relative_L2=difference,check_results=checks,failed_checks=failed,numerical_endpoint_pass=not failed,independent_full_rerun=False,physical_validation=False,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in o.rglob('*') if p.is_file()})
 p=o/'endpoint-audit.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');return result