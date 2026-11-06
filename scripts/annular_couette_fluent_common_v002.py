"""Annular Couette native HDF polygon mapping and pilot gates; no solver launch."""
from pathlib import Path
import hashlib,json,math,re
import h5py,numpy as np
TARGETS_SHA='4c3f6bdb49ffb208230a7c89f574d9c030fd0d97de800e256f791ecf64abfff9'
MESH_SHA='13d68cc96c8ea8f55f82f212315245c6e7a377223f3fcc87813328912d13f752'
GEOMETRY_SHA='d9fed1a699d1558e32de2fce7cb81776bf7294da8e4f405274d930c2625b1a97'
SETUP_SHA='bcb50febd503aaaf9f951fce98405a3cfa97ea47c2454d76e05fdd503a4cbfee'
EQUATIONS=('continuity','x-velocity','y-velocity');N=1536;NN=1600;NF=3136
COUNTS=dict(interior=3008,inner=64,outer=64)
TARGETS=dict(dt_s=.05,pilot_steps=5,pilot_end_time_s=.25,max_inner_iterations_per_step=100,residual_max=1e-6,
 minimum_orthogonal_quality=.15,wall_normal_flux_max_SI=1e-12,dimensionless_divergence_max=1e-5,
 dimensionless_divergence_RMS=1e-6,divergence_reference_time_s=1.,initial_cell_rest_max_m_s=1e-12,
 inner_wall_velocity_error_max_m_s=1e-10,outer_wall_velocity_max_m_s=1e-12,rho_kg_m3=1.,mu_Pa_s=.000125,
 omega_rad_s=1.,depth_m=1.,physical_validation=False,production_authorized=False)
def require(ok,message):
 if not ok:raise ValueError(message)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for part in iter(lambda:stream.read(1024**2),b''):h.update(part)
 return h.hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def write_json(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')
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
def geometry_check(g):
 require(g['nodes'].shape==(NN,2) and g['quads'].shape==(N,4) and g['face_nodes'].shape==(NF,2),'Wrong shared annular mesh size')
 require(np.isfinite(g['nodes']).all() and np.all(g['areas']>0) and len(np.unique(g['nodes'],axis=0))==NN,'Invalid source geometry/seam')
 require(np.array_equal(g['grid_shape'],[24,64]) and set(np.unique(g['face_zone']))=={3,4,5},'Wrong annular zones/order')
 return True
def native_time(data):
 from ansys.fluent.core.offline import _lispy
 with h5py.File(data,'r') as f:values={v[0]:v[1] for v in _lispy.parse(f['settings/Data Variables'][0].decode())[1]}
 result={key:float(values[key]) for key in ['flow-time','time-step','physical-time-step']}
 require(all(map(math.isfinite,result.values())) and result['time-step'].is_integer(),'Invalid actual native time')
 return result
def read_state(case,data,g):
 geometry_check(g)
 with h5py.File(case,'r') as c,h5py.File(data,'r') as d:
  mesh=c['meshes/1'];faces=mesh['faces'];nodes=chunk(mesh['nodes/coords'],NN,(2,))
  lookup={tuple(np.round(x,12)):i for i,x in enumerate(g['nodes'])};require(len(lookup)==NN,'Duplicate source node')
  node_map=np.array([lookup[tuple(np.round(x,12))] for x in nodes]);require(len(np.unique(node_map))==NN and np.max(abs(nodes-g['nodes'][node_map]))<1e-12,'Native/source vertex mapping differs')
  owner=chunk(faces['c0'],NF);require(np.array_equal(owner,owner.astype(int)),'Noninteger native owner');owner=owner.astype(int)-1
  nbr,nbr_seen=chunk(faces['c1'],NF,partial=True);neighbour=np.full(NF,-1,int);neighbour[nbr_seen]=nbr[nbr_seen].astype(int)-1
  require(np.all((owner>=0)&(owner<N)) and np.all((neighbour>=-1)&(neighbour<N)),'Native face owner/neighbour range')
  fnode=np.empty((NF,2),int);seen=np.zeros(NF,bool)
  for group in faces['nodes'].values():
   a=int(np.asarray(group.attrs['minId']).ravel()[0])-1;b=int(np.asarray(group.attrs['maxId']).ravel()[0]);ids=np.asarray(group['nodes'][:]).astype(int).reshape(-1,2)-1
   require(0<=a<b<=NF and not seen[a:b].any() and np.all(group['nnodes'][:]==2) and len(ids)==b-a,'Native edge coverage differs');fnode[a:b]=ids;seen[a:b]=True
  require(seen.all() and np.all((fnode>=0)&(fnode<NN)),'Missing native face vertices')
  members=[set() for _ in range(N)]
  for face,(a,b) in enumerate(zip(owner,neighbour)):
   members[a].update(map(int,node_map[fnode[face]]))
   if b>=0:members[b].update(map(int,node_map[fnode[face]]))
  by_cell={tuple(sorted(map(int,q))):i for i,q in enumerate(g['quads'])};require(len(by_cell)==N and all(len(v)==4 for v in members),'Native cell polygon membership differs')
  mapping=np.array([by_cell[tuple(sorted(v))] for v in members]);require(len(np.unique(mapping))==N,'Native cell map not bijective')
  areas=np.empty(N);centres=np.empty((N,2))
  for i,m in enumerate(members):
   ids=np.array(sorted(m));xy=g['nodes'][ids];mean=xy.mean(axis=0);xy=xy[np.argsort(np.arctan2(xy[:,1]-mean[1],xy[:,0]-mean[0]))];areas[i],centres[i]=polygon(xy)
  require(np.max(abs(areas-g['areas'][mapping]))<1e-14 and np.max(abs(centres-g['centres'][mapping]))<1e-12,'Actual native polygon volume/centroid differs')
  mid=nodes[fnode].mean(axis=1);edge=nodes[fnode[:,1]]-nodes[fnode[:,0]];raw=np.column_stack([edge[:,1],-edge[:,0]])
  align=np.sum(raw*(mid-centres[owner]),axis=1);require(np.all(abs(align)>1e-20),'Indeterminate native geometric orientation')
  area_vector=raw*np.where(align<0,-1,1)[:,None]
  edge_lookup={tuple(sorted(map(int,e))):i for i,e in enumerate(g['face_nodes'])};face_map=np.array([edge_lookup[tuple(sorted(map(int,node_map[e])))] for e in fnode])
  require(len(np.unique(face_map))==NF,'Native/source face map not bijective')
  mapped_owner=mapping[owner];mapped_neighbour=np.where(neighbour>=0,mapping[np.maximum(neighbour,0)],-1)
  source_owner=g['face_owner'][face_map];source_neighbour=g['face_neighbour'][face_map]
  forward=(mapped_owner==source_owner)&(mapped_neighbour==source_neighbour);reverse=(mapped_owner==source_neighbour)&(mapped_neighbour==source_owner)
  require(np.all(forward|reverse),'Native/source cell adjacency differs');flux_to_source_sign=np.where(forward,1,-1)
  topology=faces['zoneTopology'];names=b''.join(np.asarray(topology['name'][:]).ravel()).decode().split(';')
  lo=np.asarray(topology['minId'][:]).ravel().astype(int)-1;hi=np.asarray(topology['maxId'][:]).ravel().astype(int);types=np.asarray(topology['zoneType'][:]).ravel().astype(int)
  require(set(names)==set(COUNTS),'Unexpected native zone names');zones={n:(int(a),int(b)) for n,a,b in zip(names,lo,hi)}
  require(all(b-a==COUNTS[n] for n,(a,b) in zones.items()) and dict(zip(names,map(int,types)))==dict(interior=2,inner=3,outer=3),'Native face zone types/counts differ')
  cell=d['results/1/phase-1/cells'];face=d['results/1/phase-1/faces'];native={key:chunk(cell[name],N) for key,name in [('u','SV_U'),('v','SV_V'),('p','SV_P'),('rho','SV_DENSITY'),('mu','SV_MU_LAM')]}
  flux=chunk(face['SV_FLUX'],NF);wu,wseen=chunk(face['SV_U'],NF,partial=True);wv,vseen=chunk(face['SV_V'],NF,partial=True)
  boundary_mask=neighbour<0;require(np.all(wseen[boundary_mask]&vseen[boundary_mask]),'Missing actual saved physical wall velocities')
  common={key:np.empty(N) for key in native}
  for key in native:common[key][mapping]=native[key]
 return dict(native=native,common=common,mapping=mapping,node_mapping=node_map,face_mapping=face_map,owner=owner,neighbour=neighbour,
  flux=flux,flux_to_source_sign=flux_to_source_sign,wall_u=wu,wall_v=wv,wall_coverage=wseen&vseen,face_nodes=fnode,face_mid=mid,face_area_vector=area_vector,
  native_cell_areas=areas,native_cell_centres=centres,zones=zones,native_time=native_time(data),geometry_normal_flips=int(np.sum(align<0)),
  volume_definition='Actual saved native quad coordinates/connectivity, polygon shoelace area * proved planar reference depth1m; no ideal circle/Cartesian volumes')
def state_gates(s,step):
 t=s['native_time'];require(t['time-step']==step and abs(t['flow-time']-step*.05)<1e-12 and abs(t['physical-time-step']-.05)<1e-12,'Actual native step/time/dt differs')
 native=s['native'];require(np.max(abs(native['rho']-1))<1e-12 and np.max(abs(native['mu']-.000125))<1e-14,'Actual material differs')
 q=s['flux'];own=s['owner'];nbr=s['neighbour'];sums=np.zeros(N);np.add.at(sums,own,q);mask=nbr>=0;np.add.at(sums,nbr[mask],-q[mask])
 div=sums/(native['rho']*s['native_cell_areas']);vol=s['native_cell_areas'];maximum=float(abs(div).max());rms=float(np.sqrt(np.sum(vol*div*div)/vol.sum()))
 walls={};boundary_flux=[]
 for name in ['inner','outer']:
  a,b=s['zones'][name];uv=np.column_stack([s['wall_u'][a:b],s['wall_v'][a:b]]);normal=s['face_area_vector'][a:b];normal=normal/np.linalg.norm(normal,axis=1)[:,None]
  expected=np.column_stack([-s['face_mid'][a:b,1],s['face_mid'][a:b,0]]) if name=='inner' else np.zeros((b-a,2))
  expected-=np.sum(expected*normal,axis=1)[:,None]*normal
  err=float(np.max(abs(uv-expected)));wallflux=float(np.max(abs(q[a:b])));boundary_flux.extend(q[a:b].tolist())
  walls[name]=dict(maximum_velocity_component_error_m_s=err,maximum_normal_velocity_m_s=float(abs(np.sum(uv*normal,axis=1)).max()),maximum_absolute_native_face_mass_flux_kg_s=wallflux)
 Lz=float(np.sum(native['rho']*vol*(s['native_cell_centres'][:,0]*native['v']-s['native_cell_centres'][:,1]*native['u'])))
 checks=dict(native_materials_and_fields_finite=True,actual_step_time_dt=True,closed_wall_fluxes=all(w['maximum_absolute_native_face_mass_flux_kg_s']<=1e-12 for w in walls.values()) and abs(math.fsum(boundary_flux))<=1e-12,
  local_divergence_max=maximum<=1e-5,local_divergence_RMS=rms<=1e-6,initial_cell_rest_or_positive_angular_momentum=(float(np.max(np.hypot(native['u'],native['v'])))<=1e-12 and np.max(abs(native['p']))<=1e-12) if step==0 else Lz>0)
 if step>0:checks.update(inner_wall_velocity=walls['inner']['maximum_velocity_component_error_m_s']<=1e-10,outer_wall_velocity=walls['outer']['maximum_velocity_component_error_m_s']<=1e-12)
 return dict(step=step,time_s=t['flow-time'],checks=checks,all_checks_passed=all(checks.values()),metrics=dict(maximum_dimensionless_divergence=maximum,volume_weighted_RMS_dimensionless_divergence=rms,total_boundary_mass_flux_kg_s=math.fsum(boundary_flux),walls=walls,angular_momentum_z_per_unit_depth_kg_m_s=Lz))
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
def save_arrays(path,s):
 arrays={**{'native_'+k:v for k,v in s['native'].items()},**{'source_'+k:v for k,v in s['common'].items()}}
 for key in ['mapping','node_mapping','face_mapping','owner','neighbour','flux','flux_to_source_sign','wall_u','wall_v','wall_coverage','face_nodes','face_mid','face_area_vector','native_cell_areas','native_cell_centres']:arrays[key]=s[key]
 np.savez_compressed(path,**arrays)
