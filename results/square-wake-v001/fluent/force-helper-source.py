"""Square-obstacle geometry/force audit primitives; no solver launch."""
from pathlib import Path
import hashlib,re,math
import numpy as np,h5py
from analyse_nozzle_fluent_v005 import chunk_array,common_mapping
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
FORCE_CRITERIA=dict(window_first=1801,window_last=2000,samples=200,Cd_peak_to_peak_over_abs_mean_max=.001,absolute_mean_Cl_max=.001,comparison_coefficient_absolute_tolerance_scale=1e-6,force_normalization_N=.0005)
BOUNDARIES=('inlet','outlet','top','bottom','body')
COUNTS=dict(interior=31500,inlet=80,outlet=80,top=200,bottom=200,body=40)
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def source_files():
 root=Path(__file__).resolve().parent
 return {name:root/name for name in ('square_wake_fluent_common_v001.py','analyse_nozzle_fluent_v005.py')}
def geometry_check(g):
 assert g['quads'].shape==(15900,4) and g['nodes'].shape==(16200,2) and np.array_equal(g['grid_shape'],[80,200])
 assert np.isfinite(g['areas']).all() and np.all(g['areas']>0) and np.allclose(g['areas'],.0001,atol=1e-16,rtol=0) and abs(g['areas'].sum()-1.59)<1e-12
 occupied=np.zeros((80,200),bool);occupied[g['iy'],g['ix']]=True;expected=np.ones((80,200),bool);expected[35:45,45:55]=False
 assert np.array_equal(occupied,expected) and occupied.sum()==15900 and len(np.unique(g['quads']))==16200
 assert np.allclose(g['centres'],np.column_stack([(g['ix']+.5)*.01,(g['iy']-40+.5)*.01]),atol=1e-12,rtol=0)
 return True
def force_rows(text):
 rows=[]
 for line in text.splitlines():
  p=line.split()
  if len(p)!=3:continue
  try:i,cd,cl=map(float,p)
  except ValueError:continue
  assert all(map(math.isfinite,(i,cd,cl))) and i==int(i)
  if i>0:rows.append(dict(iteration=int(i),Cd=cd,Cl=cl))
 assert [r['iteration'] for r in rows]==list(range(1,2001)),'Require actual every-iteration native force history1..2000'
 return rows
def force_statistics(rows):
 # Completeness is a parser/source check; target misses remain independent false flags.
 assert [r['iteration'] for r in rows]==list(range(1,2001))
 window=rows[1800:2000];assert len(window)==200 and window[0]['iteration']==1801 and window[-1]['iteration']==2000
 cd=np.array([r['Cd'] for r in window]);cl=np.array([r['Cl'] for r in window]);mean=float(cd.mean());ratio=float(np.ptp(cd)/abs(mean)) if mean!=0 else None
 return dict(window=[1801,2000],samples=200,mean_Cd=mean,Cd_peak_to_peak=float(np.ptp(cd)),Cd_peak_to_peak_over_abs_mean=ratio,absolute_mean_Cl=abs(float(cl.mean())),mean_Cl=float(cl.mean()),
  Cd_stability_pass=ratio is not None and ratio<=.001,Cl_symmetry_pass=abs(float(cl.mean()))<=.001,zero_mean_Cd_diagnostic='Undefined normalized Cd variation; stability flag false' if mean==0 else None)
def force_report(path,direction,wall_name='body',normalization=.0005):
 text=Path(path).read_text(encoding='utf-8',errors='replace');assert 'Forces [N]' in text
 assert re.search(r'^\s*Forces - Direction Vector\s*\('+r'\s+'.join(map(str,direction))+r'\)\s*$',text,re.M)
 rows=re.findall(r'^\s*Net\s+((?:'+NUMBER+r'\s+){5}'+NUMBER+r')\s*$',text,re.M);assert len(rows)==1
 values=np.array(list(map(float,rows[0].split())));assert len(values)==6 and np.isfinite(values).all()
 body=re.findall(r'^\s*'+re.escape(wall_name)+r'\s+((?:'+NUMBER+r'\s+){5}'+NUMBER+r')\s*$',text,re.M);assert len(body)==1 and np.array_equal(values,np.array(list(map(float,body[0].split()))))
 # No other wall is admitted; vector table supplies an independent cross-direction check.
 vector_rows=re.findall(r'^\s*Net\s+((?:\([^\n()]+\)\s*){6})$',text,re.M);assert len(vector_rows)==1
 vectors=np.array([list(map(float,v.split())) for v in re.findall(r'\(([^()]*)\)',vector_rows[0])]);assert vectors.shape==(6,3) and np.isfinite(vectors).all()
 for force,coef in zip(values[:3],values[3:]):assert abs(force/normalization-coef)<=1e-6*max(1.,abs(coef))
 assert abs(values[0]+values[1]-values[2])/normalization<=1e-6*max(1.,abs(values[5]))
 return dict(pressure_force_N=float(values[0]),viscous_force_N=float(values[1]),total_force_N=float(values[2]),pressure_coefficient=float(values[3]),viscous_coefficient=float(values[4]),total_coefficient=float(values[5]),vectors=vectors.tolist())
def coefficient_close(a,b):return abs(a-b)<=1e-6*max(1.,abs(b))
def native_body_force(case,data):
 count=32100
 with h5py.File(case,'r') as c,h5py.File(data,'r') as d:
  root=c['meshes/1'];zt=root['faces/zoneTopology'];names=b''.join(np.asarray(zt['name'][:]).ravel()).decode().split(';')
  lo=np.asarray(zt['minId'][:]).ravel().astype(int)-1;hi=np.asarray(zt['maxId'][:]).ravel().astype(int);types=np.asarray(zt['zoneType'][:]).ravel().astype(int)
  assert set(names)==set(COUNTS);zones={name:(int(a),int(b)) for name,a,b in zip(names,lo,hi)}
  assert all(b-a==COUNTS[name] for name,(a,b) in zones.items())
  actual_types=dict(zip(names,map(int,types)));assert actual_types==dict(interior=2,inlet=10,outlet=5,top=7,bottom=7,body=3)
  nodes=chunk_array(root['nodes/coords'],16200,(2,));c0=chunk_array(root['faces/c0'],count).astype(int);face_area=np.empty((count,2));seen=np.zeros(count,bool)
  for ds in root['faces/nodes'].values():
   a=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;b=int(np.asarray(ds.attrs['maxId']).ravel()[0]);ids=np.asarray(ds['nodes'][:]).reshape(-1,2).astype(int)-1
   assert len(ids)==b-a and not seen[a:b].any() and np.all(np.asarray(ds['nnodes'][:])==2)
   edge=nodes[ids[:,1]]-nodes[ids[:,0]];face_area[a:b]=np.column_stack([edge[:,1],-edge[:,0]]);seen[a:b]=True
  assert seen.all();a,b=zones['body'];assert np.linalg.norm(face_area[a:b].sum(axis=0))<1e-12 and abs(np.linalg.norm(face_area[a:b],axis=1).sum()-.4)<1e-12
  p=chunk_array(d['results/1/phase-1/faces/SV_P'],count);pressure=np.sum(p[a:b,None]*face_area[a:b],axis=0)
  shear=np.zeros((count,3));covered=np.zeros(count,bool)
  for ds in d['results/1/phase-1/faces/SV_WALL_SHEAR'].values():
   first=int(np.asarray(ds.attrs['minId']).ravel()[0])-1;last=int(np.asarray(ds.attrs['maxId']).ravel()[0]);values=np.asarray(ds[:]).reshape(last-first,3)
   assert 0<=first<last<=count and not covered[first:last].any() and np.isfinite(values).all();shear[first:last]=values;covered[first:last]=True
  assert covered[a:b].all(),'All40 actual body shear vectors required'
  viscous=-shear[a:b,:2].sum(axis=0);total=pressure+viscous
 return dict(pressure_force_N=pressure.tolist(),viscous_force_N=viscous.tolist(),total_force_N=total.tolist(),Cd=float(total[0]/.0005),Cl=float(total[1]/.0005),body_faces=b-a,body_closed_oriented_edgearea_sum_m2=face_area[a:b].sum(axis=0).tolist(),face_zone_types=actual_types,
  interpretation='Pressure=sum(nativefaceSV_P*outward-from-fluid edgearea*depth1); viscous=-sum(nativeSV_WALL_SHEAR), integratedN. Both require matching this run native force report; no extra area multiplier on shear.')
