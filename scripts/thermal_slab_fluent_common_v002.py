"""Read-only native thermal slab primitives; no solver launch or synthetic CFD fields."""
from pathlib import Path
import hashlib,json,math,re
import h5py
import numpy as np
from analyse_nozzle_fluent_v005 import chunk_array,common_mapping
from thermal_slab_reference_v001 import cell_average_temperature
WALLS=('hot','cold','insulated_bottom','insulated_top')
COUNTS=dict(interior=3880,hot=20,cold=20,insulated_bottom=100,insulated_top=100)
REFERENCE_SHA='550c050d4294c08d0fb3905d06a2eeda7c82bb44397072a3dd261db1a6d9846a'
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
TARGETS=dict(time_step_s=2.,time_steps=100,final_time_s=200.,maximum_inner_iterations_per_step=20,energy_residual_max=1e-10,
 temperature_min_K=300.,temperature_max_K=400.,temperature_bound_tolerance_K=1e-6,y_spread_max_K=1e-6,
 global_Euler_energy_relative_max=1e-5,final_analytic_max_error_K=.5,final_analytic_weighted_rms_error_K=.25,
 series80_vs160_max_K=1e-9,minimum_orthogonal_quality=.15,reference_depth_m=1.,reference_heat_power_W=200.,
 rho_kg_m3=1000.,specific_heat_J_kg_K=1000.,thermal_conductivity_W_m_K=10.,source_area_m2=.002,
 energy_definition='sum(rho*cp*source cell area*depth*(T-300K))',
 energy_balance_definition='((E_n-E_prev)/2-sum(Q_in_n))/max(sum(abs(Q_each_wall_n)),abs((E_n-E_prev)/2),200W)',
 heat_flux_interpretation='Native SV_HEAT_FLUX is integrated face power W, inward-positive only after native wall-report and Euler-storage crosschecks; do not multiply by area again.',
 physical_validation=False,mesh_independence=False,time_step_independence=False,all_cell_thermal_flux_conservation_claim=False)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1024**2),b''):h.update(block)
 return h.hexdigest()
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
def source_files():
 p=Path(__file__).resolve().parent
 return {name:p/name for name in ('thermal_slab_fluent_common_v002.py','analyse_nozzle_fluent_v005.py','thermal_slab_reference_v001.py')}
def verify_sources(expected):
 actual={name:sha(path) for name,path in source_files().items()}
 assert actual==expected,'Reviewed helper/reference source hashes differ'
 assert actual['thermal_slab_reference_v001.py']==REFERENCE_SHA
 return actual
def geometry_check(g):
 assert g['quads'].shape==(2000,4) and g['nodes'].shape==(2121,2) and np.array_equal(g['grid_shape'],[20,100])
 assert np.all(g['areas']>0) and np.allclose(g['areas'],1e-6,rtol=0,atol=1e-18) and abs(float(g['areas'].sum())-.002)<1e-14
 assert np.array_equal(g['ix'],np.tile(np.arange(100),20)) and np.array_equal(g['iy'],np.repeat(np.arange(20),100))
 assert np.allclose(g['centres'],np.column_stack([(g['ix']+.5)*.001,(g['iy']+.5)*.001]),rtol=0,atol=1e-12)
 assert g['quads'].min()==0 and g['quads'].max()==2120
 return True
def native_time(path):
 # Same public reader implementation as DataFile.variables(); explicitly readonly and closed.
 from ansys.fluent.core.offline import _lispy
 with h5py.File(path,'r') as f:
  raw=f['settings/Data Variables'][0].decode();values={v[0]:v[1] for v in _lispy.parse(raw)[1]}
 out={key:float(values[key]) for key in ('flow-time','time-step','physical-time-step')}
 assert all(map(math.isfinite,out.values())) and out['time-step']==int(out['time-step'])
 return out
def topology(case,g):
 mapping=common_mapping(case,g);assert mapping.shape==(2000,) and len(np.unique(mapping))==2000
 with h5py.File(case,'r') as f:
  root=f['meshes/1/faces'];zt=root['zoneTopology']
  names=b''.join(np.asarray(zt['name'][:]).ravel()).decode().split(';')
  lo=np.asarray(zt['minId'][:]).reshape(-1).astype(int)-1;hi=np.asarray(zt['maxId'][:]).reshape(-1).astype(int)
  assert len(names)==len(lo)==len(hi) and set(names)==set(COUNTS)
  zones={name:(int(a),int(b)) for name,a,b in zip(names,lo,hi)}
  assert all(b-a==COUNTS[name] for name,(a,b) in zones.items())
  assert sorted(i for a,b in zones.values() for i in range(a,b))==list(range(4120))
  c0=chunk_array(root['c0'],4120).astype(int);c1=np.zeros(4120,int);seen=np.zeros(4120,bool)
  for ds in root['c1'].values():
   a=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;b=int(np.asarray(ds.attrs['maxId']).ravel()[0]);v=np.asarray(ds[:]).reshape(-1)
   assert 0<=a<b<=4120 and not seen[a:b].any() and np.array_equal(v,v.astype(int))
   c1[a:b]=v;seen[a:b]=True
  assert np.all((c0>=1)&(c0<=2000)) and np.all((c1>=0)&(c1<=2000))
  a,b=zones['interior'];assert np.all(c1[a:b]>0)
  for name in WALLS:
   a,b=zones[name];assert np.all(c1[a:b]==0)
 return mapping,zones
def temperatures(path,mapping):
 with h5py.File(path,'r') as f:native=chunk_array(f['results/1/phase-1/cells/SV_T'],2000)
 result=np.empty(2000);result[mapping]=native
 return result
def heat_power(path,zones):
 result=np.zeros(4120);seen=np.zeros(4120,bool)
 with h5py.File(path,'r') as f:
  for ds in f['results/1/phase-1/faces/SV_HEAT_FLUX'].values():
   a=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;b=int(np.asarray(ds.attrs['maxId']).ravel()[0]);v=np.asarray(ds[:]).reshape(-1)
   assert 0<=a<b<=4120 and len(v)==b-a and np.isfinite(v).all() and not seen[a:b].any()
   result[a:b]=v;seen[a:b]=True
 expected=np.zeros(4120,bool)
 for name in WALLS:
  a,b=zones[name];expected[a:b]=True
 assert np.array_equal(seen,expected),'Require exactly all240 boundary face heat powers; no invented internal flux'
 return {name:float(result[a:b].sum()) for name,(a,b) in zones.items() if name in WALLS}
def heat_report(path):
 text=Path(path).read_text(encoding='utf-8',errors='replace');result={}
 assert re.search(r'heat\s+transfer',text,re.I) and re.search(r'\(\s*[wW]\s*\)',text),'Native report must label heat transfer in W'
 for name in (*WALLS,'Net'):
  matches=re.findall(r'^\s*'+re.escape(name)+r'\s+('+NUMBER+r')\s*$',text,re.M)
  assert matches,'Missing native heat report row '+name
  result[name]=float(matches[-1]);assert math.isfinite(result[name])
 return result
def energy_residuals(path):
 with h5py.File(path,'r') as f:
  root=f['results/residuals/phase-1'];assert 'energy' in root,'Native energy residual group absent'
  group=root['energy'];iterations=np.asarray(group['iterations'][:]).reshape(-1);data=np.asarray(group['data'][:])
 assert len(iterations) and np.array_equal(iterations,iterations.astype(int)) and np.all(np.diff(iterations)>0)
 assert data.ndim==2 and data.shape[0]==len(iterations) and data.shape[1]>=2 and np.isfinite(data[:,:2]).all()
 assert np.all(data[:,1]>0),'No undefined thermal normalization or invented0/0 residual'
 scaled=data[:,0]/data[:,1];assert np.isfinite(scaled).all() and np.all(scaled>=0)
 return iterations.astype(int),scaled

def cell_zone_snapshot(parent):
 # Whole parent state omits inactive collections; observed v001 contains only solid.slab.
 assert isinstance(parent,dict) and not parent.get('fluid',{}) and set(parent.get('solid',{}))=={'slab'}
 return dict(cell_zone_conditions=parent,fluid_zones=parent.get('fluid',{}),solid_zones=parent['solid'],fluid_collection_present='fluid' in parent)
def physical_check(s):
 assert s['general']['solver']['time']=='transient' and s['general']['solver']['two_dim_space']=='planar'
 assert s['energy']['enabled'] is True and not s['fluid_zones'] and set(s['solid_zones'])=={'slab'}
 zones=cell_zone_snapshot(s['cell_zone_conditions']);assert all(s[key]==value for key,value in zones.items())
 assert s['solid_zones']['slab']['general']['material']=='slab-solid'
 for key,value in (('density',1000.),('specific_heat',1000.),('thermal_conductivity',10.)):
  assert s['material'][key]['option']=='constant' and s['material'][key]['value']==value
 assert s['references']['depth']==1. and s['methods']['transient_formulation']=='unsteady-1st-order'
 eq=s['controls']['equations'];assert isinstance(eq,dict) and eq
 active=[str(k).lower() for k,v in eq.items() if v is True]
 assert active and all(any(t in k for t in ('energy','temperature')) for k in active),'Only thermal equations may be active'
 run=s['run_calculation'];assert run['transient_controls']['type']=='Fixed' and run['transient_controls']['method']=='User-Specified'
 assert run['parameters']['time_step_size']==2 and run['parameters']['max_iter_per_time_step']==20
 solid_options=run['options'].get('solid_time_step_size',{})
 assert solid_options.get('enable_solid_time_step') is False
 for name in WALLS:
  wall=s['walls'][name]['thermal']
  if name in ('hot','cold'):
   assert wall['thermal_condition']=='Temperature' and wall['temperature']['value']==(400. if name=='hot' else 300.)
  else:assert wall['thermal_condition']=='Heat Flux' and wall['heat_flux']['value']==0.
 residual=s['residual'];assert residual['options']['criterion_type']=='absolute'
 eq=residual['equations']['energy'];assert eq['check_convergence'] is True and eq['absolute_criteria']==1e-10
 return True

def state_metrics(T,g,time_s):
 assert T.shape==(2000,) and np.isfinite(T).all()
 assert T.min()>=300-1e-6 and T.max()<=400+1e-6,'Native temperature bounds failed'
 spread=float(np.ptp(T.reshape(20,100),axis=0).max());assert spread<=1e-6,'Native y-uniformity failed'
 left=g['x_edges'][g['ix']];right=g['x_edges'][g['ix']+1]
 analytic=cell_average_temperature(left,right,time_s,80);reference160=cell_average_temperature(left,right,time_s,160)
 series=float(abs(analytic-reference160).max());assert time_s==0 or series<1e-9
 error=T-analytic;areas=g['areas']
 return dict(temperature_min_K=float(T.min()),temperature_max_K=float(T.max()),max_y_spread_K=spread,
  energy_J=float(np.sum(1000*1000*areas*(T-300))),analytic_max_error_K=float(abs(error).max()),
  analytic_rms_error_K=float(np.sqrt(np.sum(areas*error**2)/areas.sum())),series80_vs160_max_K=series)

def step_metrics(data,report,mapping,zones,g,step,previous_energy,previous_iteration):
 time=native_time(data);assert time['time-step']==step and abs(time['flow-time']-2*step)<1e-9 and abs(time['physical-time-step']-2)<1e-12
 T=temperatures(data,mapping);metrics=state_metrics(T,g,2*step);q=heat_power(data,zones);native_report=heat_report(report)
 # Matching two independent native representations establishes integrated W and inward sign.
 for name in WALLS:
  assert abs(q[name]-native_report[name])<=1e-6*max(1.,abs(q[name])), 'Native report/HDF integrated W mismatch: '+name
 assert abs(sum(q.values())-native_report['Net'])<=1e-6*max(1.,abs(sum(q.values())))
 assert q['hot']>0 and q['cold']<=1e-7 and all(abs(q[k])<=1e-7 for k in WALLS[2:]),'Native heat powers do not establish expected inward-positive convention'
 storage=(metrics['energy_J']-previous_energy)/2;defect=storage-sum(q.values());relative=abs(defect)/max(sum(abs(v) for v in q.values()),abs(storage),200.)
 assert storage>0 and sum(q.values())>0,'Heating storage and native net inward power must both be positive'
 assert relative<=1e-5,'Global implicit Euler energy balance target missed'
 iterations,residuals=energy_residuals(data);last=int(iterations[-1]);inner=last-previous_iteration
 assert 1<=inner<=20,'Actual inner iteration count outside declaredcap'
 assert residuals[-1]<=1e-10,'Native terminal energy residual target missed'
 return T,dict(step=step,time_s=time['flow-time'],native_dt_s=time['physical-time-step'],native_time_index=int(time['time-step']),
  native_iteration=last,inner_iterations=inner,energy_residual=float(residuals[-1]),wall_heat_in_W=q,native_wall_report_W=native_report,
  storage_power_W=storage,net_inward_heat_W=sum(q.values()),energy_defect_W=defect,energy_balance_relative=relative,**metrics)
