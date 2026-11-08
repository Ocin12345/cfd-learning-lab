from pathlib import Path
import re
import h5py
import numpy as np
from plot_vortex_shedding_v002 import indexed,fluent_mesh,foam_list,foam_labels,foam_field,foam_text

def angular_stencil(centres,target_angle):
 xy=np.asarray(centres,dtype=float);theta=np.arctan2(xy[:,1]-.2,xy[:,0]-.2);delta=(theta-target_angle+np.pi)%(2*np.pi)-np.pi
 left=np.where(delta<0)[0];right=np.where(delta>0)[0]
 if not len(left) or not len(right):raise ValueError('Two-sided angular wall-face stencil required')
 i=int(left[np.argmax(delta[left])]);j=int(right[np.argmin(delta[right])]);a,b=float(delta[i]),float(delta[j])
 if max(abs(a),abs(b))>.1:raise ValueError('Wall faces too far from requested reference point')
 weights=[b/(b-a),-a/(b-a)]
 return dict(indices=[i,j],weights=weights,angle_offsets_rad=[a,b],face_centres_m=xy[[i,j]].tolist(),target_m=[float(.2+.05*np.cos(target_angle)),float(.2+.05*np.sin(target_angle))])

def fluent_wall_geometry(case):
 mesh=fluent_mesh(case);xy=mesh['xy']
 with h5py.File(case) as h:
  m=h['meshes/1'];z=m['faces/zoneTopology'];names=z['name'][0].decode().split(';');k=names.index('cylinder');first,last=int(z['minId'][k])-1,int(z['maxId'][k]);count=int(np.max(z['maxId'][:]));edges=np.full((count,2),-1,dtype=int)
  for b in m['faces/nodes'].values():
   d=b['nnodes'];start=int(np.asarray(d.attrs.get('minId',b.attrs.get('minId',[1]))).ravel()[0])-1;assert (d[:]==2).all();values=np.asarray(b['nodes'][:],dtype=int).reshape(-1,2)-1;edges[start:start+len(values)]=values
  assert (edges>=0).all();ids=np.arange(first,last);centres=xy[edges[ids]].mean(axis=1)
 for target in [[.15,.2],[.25,.2]]:assert np.min(np.linalg.norm(xy-target,axis=1))<1e-10
 return dict(centres=centres,indices=ids,total_faces=count,stencils={'front':angular_stencil(centres,np.pi),'rear':angular_stencil(centres,0.)},meaning='Native Fluent face SV_P in gauge Pa; interpolated in polar angle to the specified cylinder vertices.')

def openfoam_wall_geometry(case):
 case=Path(case);d=next(p for p in [case/'constant/polyMesh',case/'0/polyMesh'] if p.is_dir());n,raw=foam_list(d/'points');points=np.fromstring(raw.replace('(',' ').replace(')',' '),sep=' ').reshape(n,3);n,raw=foam_list(d/'faces');matches=list(re.finditer(r'(\d+)\s*\(([^()]*)\)',raw));assert len(matches)==n
 raw=foam_text(d/'boundary');b=re.search(r'\bcylinder\s*\{([^{}]*)\}',raw).group(1);first=int(re.search(r'startFace\s+(\d+)',b)[1]);count=int(re.search(r'nFaces\s+(\d+)',b)[1]);ids=np.arange(first,first+count);centres=np.array([points[np.fromstring(matches[k][2],sep=' ',dtype=int),:2].mean(axis=0) for k in ids]);owners=foam_labels(d/'owner')[ids]
 for target in [[.15,.2],[.25,.2]]:assert np.min(np.linalg.norm(points[:,:2]-target,axis=1))<1e-10
 return dict(centres=centres,indices=ids,owners=owners,stencils={'front':angular_stencil(centres,np.pi),'rear':angular_stencil(centres,0.)},meaning='Native zeroGradient wall face value equals owner-cell kinematic p; rho=1 converts numerically to Pa. This is evaluated FV wall pressure followed by the same angular interpolation, not an exact continuum point measurement.')

def wall_pressures(solver,path,geometry):
 if solver=='fluent':
  with h5py.File(path) as h:p=indexed(h['results/1/phase-1/faces/SV_P'],geometry['total_faces'])[geometry['indices']]
 else:
  raw=foam_text(Path(path)/'p');block=re.search(r'\bcylinder\s*\{([^{}]*)\}',raw).group(1);assert re.search(r'type\s+zeroGradient\s*;',block)
  p=foam_field(Path(path)/'p',21598,1,[0,2,-2,0,0,0,0])[geometry['owners']]
 assert np.isfinite(p).all()
 values={name:float(np.dot(p[st['indices']],st['weights'])) for name,st in geometry['stencils'].items()};values['difference_Pa']=values['front']-values['rear'];return values