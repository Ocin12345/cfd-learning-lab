"""Tandem native polygon, flux and separate-body force primitives; no CFD on import."""
from pathlib import Path
import hashlib,importlib.util,json,math,re,sys
import h5py,numpy as np
ROOT=Path('E:/CFD-Learning-Lab');LOCAL=ROOT/'local/tandem-offset-wake-v001';RUN=ROOT/'native/runs/2026-11-02-tandem-offset-wake-v001'
MESH=RUN/'mesh/tandem-offset-wake.msh'
TARGETS_SHA='426448b7e6909e75dbf08f065d93f90699b98f12be007891631d9ed9b10f71d7'
PLAN_SHA='79f894cdddbf65209b91b6cfd0121909eb69428e9dcf101ae267c0112ec5e055'
MESH_SHA='f5db9ca48d6e5bcc27d387febbc36fe1e2239f1deac62061362ce902cfe0583d'
GEOMETRY_SHA='099165a86f06dd6566c6d76860f5a0935fba4b956e16d6908f6f1b735270b178'
SETUP_SHA='263f87d0068959e316215e85c6e1391627bfd5d6d87e68b5d4c6b13fa0324bb4'
EQUATIONS=('continuity','x-velocity','y-velocity');BODIES=('body1','body2');BOUNDARIES=('inlet','outlet','top','bottom',*BODIES)
COUNTS=dict(interior=18692,inlet=60,outlet=60,top=160,bottom=160,body1=32,body2=32)
TYPES=dict(interior=2,inlet=10,outlet=5,top=7,bottom=7,body1=3,body2=3)
N=9472;NN=9723;NF=19196;NORMALIZATION=.0005
NUMBER=r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
def require(ok,message):
 if not ok:raise ValueError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
def checked_module(path,digest):
 path=Path(path);require(sha(path)==digest,'Changed source before import: '+str(path));spec=importlib.util.spec_from_file_location(path.stem+'_'+digest[:12],path);m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=m;spec.loader.exec_module(m);return m
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def json_default(v):
 if isinstance(v,np.generic):return v.item()
 if isinstance(v,np.ndarray):return v.tolist()
 if isinstance(v,Path):return str(v)
 raise TypeError('Unsupported JSON value '+type(v).__name__)
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=json_default)+'\n',encoding='utf-8')
def geometry_check(g):
 require(g['nodes'].shape==(NN,2) and g['quads'].shape==(N,4) and g['face_nodes'].shape==(NF,2),'Wrong actual shared mesh counts')
 require(np.isfinite(g['nodes']).all() and len(np.unique(g['nodes'],axis=0))==NN and np.array_equal(g['grid_shape'],[60,160]),'Invalid vertices/grid')
 mask=np.ones((60,160),bool);mask[24:32,36:44]=False;mask[28:36,68:76]=False;actual=np.zeros_like(mask);actual[g['iy'],g['ix']]=True
 require(np.array_equal(actual,mask) and np.allclose(g['areas'],.00015625,atol=1e-15,rtol=0) and abs(float(g['areas'].sum())-1.48)<1e-12,'Two-hole occupancy/area differs')
 require(set(map(int,np.unique(g['face_zone'])))==set(range(3,10)),'Wrong actual source zones');return True
def chunk(group,count,tail=(),partial=False):
 out=np.full((count,*tail),np.nan);seen=np.zeros(count,bool)
 for d in group.values():
  a=int(np.asarray(d.attrs['minId']).ravel()[0])-1;b=int(np.asarray(d.attrs['maxId']).ravel()[0]);values=np.asarray(d[:]).reshape((b-a,*tail))
  require(0<=a<b<=count and not seen[a:b].any() and np.isfinite(values).all(),'Invalid/missing/overlapping HDF chunk')
  out[a:b]=values;seen[a:b]=True
 require(partial or seen.all(),'Incomplete native HDF field')
 return (out,seen) if partial else out
def polygon(xy):
 cross=xy[:,0]*np.roll(xy[:,1],-1)-np.roll(xy[:,0],-1)*xy[:,1];area=float(cross.sum()/2)
 require(area>0 and math.isfinite(area),'Nonpositive native polygon area')
 centre=np.array([np.sum((xy[:,k]+np.roll(xy[:,k],-1))*cross)/(6*area) for k in [0,1]])
 return area,centre
def mesh_fingerprint(case):
 h=hashlib.sha256()
 with h5py.File(case,'r') as f:
  datasets=[];f['meshes/1'].visititems(lambda name,item:datasets.append(name) if isinstance(item,h5py.Dataset) else None)
  for name in sorted(datasets):
   d=f['meshes/1/'+name];h.update(name.encode());h.update(str((d.shape,str(d.dtype))).encode());h.update(np.asarray(d[:]).tobytes())
 return h.hexdigest()
def native_time(data):
 from ansys.fluent.core.offline import _lispy
 with h5py.File(data,'r') as f:values={v[0]:v[1] for v in _lispy.parse(f['settings/Data Variables'][0].decode())[1]}
 result={key:float(values[key]) for key in ['flow-time','time-step','physical-time-step']};require(all(map(math.isfinite,result.values())) and result['time-step'].is_integer(),'Invalid native clock');return result
def map_native_mesh(case,g,counts=COUNTS,types=TYPES):
 n=len(g['quads']);nn=len(g['nodes']);nf=len(g['face_nodes'])
 with h5py.File(case,'r') as f:
  mesh=f['meshes/1'];faces=mesh['faces'];nodes=chunk(mesh['nodes/coords'],nn,(2,))
  lookup={tuple(np.round(x,12)):i for i,x in enumerate(g['nodes'])};require(len(lookup)==nn,'Duplicate source node')
  node_map=np.array([lookup[tuple(np.round(x,12))] for x in nodes]);require(len(np.unique(node_map))==nn and np.max(abs(nodes-g['nodes'][node_map]))<1e-12,'Native/source vertices differ')
  owner=chunk(faces['c0'],nf);require(np.array_equal(owner,owner.astype(int)),'Noninteger native owner');owner=owner.astype(int)-1
  v,seen=chunk(faces['c1'],nf,partial=True);neighbour=np.full(nf,-1,int);require(np.array_equal(v[seen],v[seen].astype(int)),'Noninteger neighbour');neighbour[seen]=v[seen].astype(int)-1
  require(np.all((owner>=0)&(owner<n)) and np.all((neighbour>=-1)&(neighbour<n)),'Native adjacency outside cell range')
  fnode=np.empty((nf,2),int);seen=np.zeros(nf,bool)
  for z in faces['nodes'].values():
   a=int(np.asarray(z.attrs['minId']).ravel()[0])-1;b=int(np.asarray(z.attrs['maxId']).ravel()[0]);ids=np.asarray(z['nodes'][:]).astype(int).reshape(-1,2)-1
   require(0<=a<b<=nf and not seen[a:b].any() and np.all(z['nnodes'][:]==2) and len(ids)==b-a,'Native edge coverage differs');fnode[a:b]=ids;seen[a:b]=True
  require(seen.all() and np.all((fnode>=0)&(fnode<nn)),'Missing native edge vertices')
  members=[set() for _ in range(n)]
  for face,(a,b) in enumerate(zip(owner,neighbour)):
   members[a].update(map(int,node_map[fnode[face]]))
   if b>=0:members[b].update(map(int,node_map[fnode[face]]))
  by_cell={tuple(sorted(map(int,q))):i for i,q in enumerate(g['quads'])};require(len(by_cell)==n and all(len(v)==4 for v in members),'Native quad membership differs')
  mapping=np.array([by_cell[tuple(sorted(v))] for v in members]);require(len(np.unique(mapping))==n,'Cell map not bijective');inverse_nodes=np.argsort(node_map)
  areas=np.empty(n);centres=np.empty((n,2))
  for i,m in enumerate(members):
   xy=nodes[inverse_nodes[np.array(sorted(m))]];mean=xy.mean(axis=0);xy=xy[np.argsort(np.arctan2(xy[:,1]-mean[1],xy[:,0]-mean[0]))];areas[i],centres[i]=polygon(xy)
  require(np.max(abs(areas-g['areas'][mapping]))<1e-14 and np.max(abs(centres-g['centres'][mapping]))<1e-12,'Actual native polygon area/centroid differs')
  midpoint=nodes[fnode].mean(axis=1);edge=nodes[fnode[:,1]]-nodes[fnode[:,0]];raw=np.column_stack([edge[:,1],-edge[:,0]]);align=np.sum(raw*(midpoint-centres[owner]),axis=1);require(np.all(abs(align)>1e-20),'Indeterminate geometric orientation');av=raw*np.where(align<0,-1,1)[:,None]
  lookup={tuple(sorted(map(int,e))):i for i,e in enumerate(g['face_nodes'])};face_map=np.array([lookup[tuple(sorted(map(int,node_map[e])))] for e in fnode]);require(len(np.unique(face_map))==nf,'Face map not bijective')
  mo=mapping[owner];mn=np.where(neighbour>=0,mapping[np.maximum(neighbour,0)],-1);so=g['face_owner'][face_map];sn=g['face_neighbour'][face_map];forward=(mo==so)&(mn==sn);reverse=(mo==sn)&(mn==so);require(np.all(forward|reverse),'Native/source face adjacency differs')
  z=faces['zoneTopology'];names=b''.join(np.asarray(z['name'][:]).ravel()).decode().split(';');lo=np.asarray(z['minId'][:]).ravel().astype(int)-1;hi=np.asarray(z['maxId'][:]).ravel().astype(int);zt=np.asarray(z['zoneType'][:]).ravel().astype(int)
  require(set(names)==set(counts) and dict(zip(names,map(int,zt)))==types,'Native zones/types differ');zones={name:(int(a),int(b)) for name,a,b in zip(names,lo,hi)};require(all(b-a==counts[name] for name,(a,b) in zones.items()),'Native zone face counts differ')
  ids=np.asarray(z['id'][:]).ravel().astype(int);require(all(np.all(g['face_zone'][face_map[a:b]]==zone_id) for zone_id,a,b in zip(ids,lo,hi)),'Native/source face-zone identity differs')
 return dict(mapping=mapping,node_mapping=node_map,face_mapping=face_map,owner=owner,neighbour=neighbour,flux_to_source_sign=np.where(forward,1,-1),face_nodes=fnode,face_mid=midpoint,face_area_vector=av,native_cell_areas=areas,native_cell_centres=centres,zones=zones,mesh_fingerprint=mesh_fingerprint(case),volume_definition='Actual saved native vertices/connectivity; positive polygon area times actual planar reference depth1m')
def read_fields(data,mesh):
 n=len(mesh['native_cell_areas']);nf=len(mesh['owner'])
 with h5py.File(data,'r') as f:
  cells=f['results/1/phase-1/cells'];faces=f['results/1/phase-1/faces'];native={k:chunk(cells[v],n) for k,v in [('u','SV_U'),('v','SV_V'),('p','SV_P'),('rho','SV_DENSITY'),('mu','SV_MU_LAM')]};flux=chunk(faces['SV_FLUX'],nf)
  u,uc=chunk(faces['SV_U'],nf,partial=True);v,vc=chunk(faces['SV_V'],nf,partial=True);p,pc=chunk(faces['SV_P'],nf,partial=True)
  if 'SV_WALL_SHEAR' in faces:shear,sc=chunk(faces['SV_WALL_SHEAR'],nf,(3,),partial=True)
  else:shear=np.full((nf,3),np.nan);sc=np.zeros(nf,bool)
 common={k:np.empty(n) for k in native}
 for k in native:common[k][mesh['mapping']]=native[k]
 return dict(**mesh,native=native,common=common,flux=flux,wall_u=u,wall_v=v,wall_coverage=uc&vc,face_pressure=p,pressure_coverage=pc,wall_shear=shear,shear_coverage=sc,native_time=native_time(data))
def read_state(case,data,g,mesh=None):
 geometry_check(g)
 if mesh is None:mesh=map_native_mesh(case,g)
 else:require(mesh_fingerprint(case)==mesh['mesh_fingerprint'],'Saved native mesh changed')
 return read_fields(data,mesh)
def flux_metrics(s):
 q=s['flux'];own=s['owner'];nbr=s['neighbour'];rho=s['native']['rho'];vol=s['native_cell_areas'];n=len(vol);require(np.isfinite(q).all() and np.all(vol>0) and np.all(rho>0),'Invalid flux/material/volume')
 sums=np.zeros(n);absolute=np.zeros(n);np.add.at(sums,own,q);np.add.at(absolute,own,abs(q));mask=nbr>=0;np.add.at(sums,nbr[mask],-q[mask]);np.add.at(absolute,nbr[mask],abs(q[mask]));div=sums/(rho*vol);co=.05*absolute/(2*rho*vol)
 return dict(maximum_dimensionless_divergence=float(abs(div).max()),volume_weighted_RMS_dimensionless_divergence=float(np.sqrt(np.sum(vol*div*div)/vol.sum())),maximum_matched_cell_Courant=float(co.max()),volume_weighted_mean_matched_cell_Courant=float(np.sum(vol*co)/vol.sum()),native_total_boundary_mass_flux_kg_s=float(math.fsum(q[nbr<0])),cell_boundary_incidence_error_kg_s=float(abs(math.fsum(sums)-math.fsum(q[nbr<0]))))
def state_gates(s,step):
 t=s['native_time'];require(t['time-step']==step and abs(t['flow-time']-step*.05)<1e-12 and abs(t['physical-time-step']-.05)<1e-12,'Actual native step/time/dt differs');native=s['native']
 require(all(np.isfinite(v).all() for v in native.values()) and np.isfinite(s['flux']).all(),'Nonfinite native primary field');require(np.max(abs(native['rho']-1))<=1e-12 and np.max(abs(native['mu']-.0001))<=1e-14,'Actual native material differs')
 metrics=flux_metrics(s);q=s['flux'];zones=s['zones'];a,b=zones['inlet'];influx=float(math.fsum(q[a:b]));den=abs(influx);a,b=zones['outlet'];reverse=float(np.maximum(-q[a:b],0).sum());metrics.update(inlet_outward_mass_flux_kg_s=influx,outlet_minimum_face_mass_flux_kg_s=float(q[a:b].min()),outlet_reverse_mass_flux_kg_s=reverse,outlet_reverse_flow_ratio=reverse/den if den>0 else None,relative_mass_imbalance=abs(metrics['native_total_boundary_mass_flux_kg_s'])/den if den>0 else None)
 checks=dict(native_fields_material_finite=True,actual_native_clock=True)
 if step==0:checks['initial_native_cells_at_rest']=float(np.max(np.hypot(native['u'],native['v'])))<=1e-12 and float(np.max(abs(native['p'])))<=1e-12
 else:
  a,b=zones['inlet'];coverage=bool(s['wall_coverage'][a:b].all());inlet_error=float(np.max(abs(s['wall_u'][a:b]-.1))) if coverage else None;iv_error=float(np.max(abs(s['wall_v'][a:b]))) if coverage else None
  wall_ok=True;closed=0.
  for name in ['top','bottom',*BODIES]:
   a,b=zones[name];closed=max(closed,float(np.max(abs(q[a:b]))))
   if name in BODIES:wall_ok=wall_ok and bool(s['wall_coverage'][a:b].all()) and bool(np.max(np.hypot(s['wall_u'][a:b],s['wall_v'][a:b]))<=1e-12)
  metrics.update(maximum_closed_face_mass_flux_kg_s=closed,inlet_u_error_m_s=inlet_error,inlet_v_error_m_s=iv_error)
  checks.update(native_inlet_flux=influx<0 and abs(influx+.075)/.075<=1e-10,native_inlet_velocity=coverage and inlet_error<=1e-12 and iv_error<=1e-12,stationary_no_slip_bodies=wall_ok,closed_body_symmetry_flux=closed<=1e-12,relative_mass_imbalance=den>0 and metrics['relative_mass_imbalance']<=1e-5,outlet_reverse_flow=den>0 and metrics['outlet_reverse_flow_ratio']<=1e-8,local_divergence_max=metrics['maximum_dimensionless_divergence']<=1e-5,local_divergence_RMS=metrics['volume_weighted_RMS_dimensionless_divergence']<=1e-6)
 checks={k:bool(v) for k,v in checks.items()};return dict(step=step,time_s=t['flow-time'],checks=checks,all_checks_passed=all(checks.values()),metrics=metrics)
def integrate_body(pressure,shear,area_vector):
 require(np.isfinite(pressure).all() and np.isfinite(shear).all() and np.isfinite(area_vector).all(),'Nonfinite native force input');p=np.sum(pressure[:,None]*area_vector,axis=0);v=-np.sum(shear[:,:2],axis=0);total=p+v
 return dict(pressure_force_N=p.tolist(),viscous_force_N=v.tolist(),total_force_N=total.tolist(),pressure_coefficient=(p/NORMALIZATION).tolist(),viscous_coefficient=(v/NORMALIZATION).tolist(),total_coefficient=(total/NORMALIZATION).tolist())
def body_forces(s,bodies=BODIES):
 result={}
 for name in bodies:
  a,b=s['zones'][name];require(s['pressure_coverage'][a:b].all() and s['shear_coverage'][a:b].all(),'Missing native pressure/shear '+name);av=s['face_area_vector'][a:b];require(np.linalg.norm(av.sum(axis=0))<1e-12 and abs(np.linalg.norm(av,axis=1).sum()-.4)<1e-12,'Native body not closed0.1m square')
  result[name]=dict(**integrate_body(s['face_pressure'][a:b],s['wall_shear'][a:b],av),faces=b-a,normalization_N=NORMALIZATION,native_shear_already_integrated_N=True,pressure_area_vector_points_outward_from_fluid=True)
 return result
def terminal_residuals(data,transcript):
 output={};indices={};text=Path(transcript).read_text(encoding='utf-8',errors='replace')
 printed={int(m[1]):[float(m[i]) for i in [2,3,4]] for m in re.finditer(r'^\s*(\d+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+\d+:\d+:\d+\s+\d+\s*$',text,re.M)}
 with h5py.File(data,'r') as f:
  for i,name in enumerate(EQUATIONS):
   group=f['results/residuals/phase-1/'+name];k=np.asarray(group['iterations'][:]).ravel();a=np.asarray(group['data'][:]);require(len(k) and np.array_equal(k,k.astype(int)) and np.all(np.diff(k)>0) and a.shape[0]==len(k) and a.shape[1]>=2 and np.isfinite(a[:,:2]).all(),'Invalid native residual history')
   require(a[-1,1]>=0,'Negative residual denominator');iteration=int(k[-1])
   if a[-1,1]>0:value=float(a[-1,0]/a[-1,1])
   else:
    require(iteration in printed and float(format(float(a[-1,0]),'.5g'))==printed[iteration][i],'Zero residual normalization lacks printed native numerator proof');value=printed[iteration][i]
   require(math.isfinite(value) and value>=0,'Invalid final native residual');output[name]=value;indices[name]=iteration
 require(len(set(indices.values()))==1,'Residual equation iteration mismatch')
 return output,next(iter(indices.values()))
def coefficient_close(a,b):return abs(a-b)<=1e-6*max(1.,abs(b))
def force_report(path,direction,body):
 text=Path(path).read_text(encoding='utf-8',errors='replace');require('Forces [N]' in text and re.search(r'^\s*Forces - Direction Vector\s*\('+r'\s+'.join(map(str,direction))+r'\)\s*$',text,re.M),'Wrong native force report units/direction')
 pattern=r'^\s*Net\s+((?:'+NUMBER+r'\s+){5}'+NUMBER+r')\s*$';rows=re.findall(pattern,text,re.M);require(len(rows)==1,'Missing/duplicate Net force');values=np.array(list(map(float,rows[0].split())));rows=re.findall(r'^\s*'+re.escape(body)+r'\s+((?:'+NUMBER+r'\s+){5}'+NUMBER+r')\s*$',text,re.M);require(len(rows)==1 and np.array_equal(values,np.array(list(map(float,rows[0].split())))),'Native report contains wrong body/Net')
 rows=re.findall(r'^\s*Net\s+((?:\([^\n()]+\)\s*){6})$',text,re.M);require(len(rows)==1,'Missing native force vectors');vectors=np.array([list(map(float,v.split())) for v in re.findall(r'\(([^()]*)\)',rows[0])]);require(vectors.shape==(6,3) and np.isfinite(vectors).all() and np.isfinite(values).all(),'Nonfinite native force report')
 for i in range(6):require(coefficient_close(float(vectors[i]@direction)/(NORMALIZATION if i<3 else 1),float(values[i])/(NORMALIZATION if i<3 else 1)),'Native force vector/directional table differs')
 for f,k in zip(values[:3],values[3:]):require(coefficient_close(f/NORMALIZATION,k),'Native force normalization differs')
 require(coefficient_close((values[0]+values[1])/NORMALIZATION,values[5]),'Native pressure+viscous != total')
 return dict(pressure_force_N=float(values[0]),viscous_force_N=float(values[1]),total_force_N=float(values[2]),pressure_coefficient=float(values[3]),viscous_coefficient=float(values[4]),total_coefficient=float(values[5]),vectors=vectors.tolist())
def force_report_agreement(forces,body,direction,path):
 report=force_report(path,direction,body);axis=0 if direction==[1,0,0] else 1
 for term in ['pressure','viscous','total']:require(coefficient_close(forces[body][term+'_coefficient'][axis],report[term+'_coefficient']),'Native HDF/native report force disagrees '+body+' '+term)
 return report
def settings_gate(st):
 general=st['general'];solver=general['solver'];require(solver['type']=='pressure-based' and solver['two_dim_space']=='planar' and solver['time']=='transient','Wrong native solver formulation');require(general['operating_conditions']['operating_pressure']==0 and general['operating_conditions']['gravity']['enable'] is False,'Wrong pressure/gravity')
 require(st['energy']['enabled'] is False and st['viscous']['model']=='laminar' and st['dynamic_mesh']['enabled'] is False,'Unexpected native physical model');require(st['material']['density']['value']==1 and st['material']['viscosity']['value']==.0001,'Wrong actual native material')
 zone=st['cell_zones']['fluid'];require(zone['reference_frame']['frame_motion'] is False and zone['mesh_motion']['enable'] is False,'Moving/rotating cell zone')
 for name in BODIES:
  wall=st['walls'][name]['momentum'];require(wall['wall_motion']=='Stationary Wall' and wall['shear_condition']=='No Slip','Wrong stationary no-slip body '+name)
 inlet=st['inlets']['inlet']['momentum'];require(inlet['velocity_specification_method']=='Components' and [v['value'] for v in inlet['velocity_components']]==[.1,0.],'Wrong actual inlet');require(st['outlet']['momentum']['gauge_pressure']['value']==0 and set(st['symmetries'])=={'top','bottom'},'Wrong outlet/symmetry')
 refs=st['references']
 for k,v in dict(density=1.,viscosity=.0001,velocity=.1,length=.1,area=.1,depth=1.,pressure=0.).items():require(refs[k]==v,'Wrong native force reference '+k)
 methods=st['methods'];schemes=methods['spatial_discretization'];require(methods['p_v_coupling']['flow_scheme']=='SIMPLE' and methods['transient_formulation']=='unsteady-1st-order' and schemes['gradient_scheme']=='least-square-cell-based' and schemes['discretization_scheme']['mom']=='first-order-upwind' and schemes['discretization_scheme']['pressure']=='standard','Wrong temporal/spatial method')
 params=st['run_calculation']['parameters'];require(params['time_step_size']==.05 and params['max_iter_per_time_step']==100 and st['run_calculation']['transient_controls']['type']=='Fixed','Wrong native time controls')
 require(all(st['residual']['equations'][name]['absolute_criteria']==1e-6 and st['residual']['equations'][name]['check_convergence'] is True for name in EQUATIONS),'Wrong own residual controls');return True
