"""Read actual native cell fields and explicit wall/adjacent pressure samples."""
from pathlib import Path
import json,re
import h5py,numpy as np
import plot_openfoam_airfoil_fields_v002 as top
import analyse_openfoam_airfoil_v002 as av
from plot_vortex_shedding_v002 import indexed

def fluent(run):
 case=run/'airfoil-Re100-final.cas.h5';data=run/'airfoil-Re100-final.dat.h5'
 with h5py.File(case) as h:
  m=h['meshes/1'];xy=m['nodes/coords/1'][:];n=int(m['cells/zoneTopology/maxId'][0]);z=m['faces/zoneTopology'];nf=int(np.max(z['maxId'][:]));names=z['name'][0].decode().split(';');w=names.index('airfoil');wallids=np.arange(int(z['minId'][w])-1,int(z['maxId'][w]))
  edges=m['faces/nodes/1/nodes'][:].reshape(-1,2).astype('int64')-1;assert len(edges)==nf
  c0=indexed(m['faces/c0'],nf).astype('int64')-1
  d=m['faces/c1/1'];first=int(d.attrs['minId'][0])-1;c1=d[:].astype('int64')-1
  keys=np.unique(np.r_[c0[:,None]*len(xy)+edges,c1[:,None]*len(xy)+edges[first:first+len(c1)]])
  ids,count=np.unique(keys//len(xy),return_counts=True);assert np.array_equal(ids,np.arange(n)) and np.all(count==3)
  tri=(keys%len(xy)).reshape(-1,3);wall=edges[wallids];owners=c0[wallids]
 with h5py.File(data) as h:
  f=h['results/1/phase-1'];u=indexed(f['cells/SV_U'],n);v=indexed(f['cells/SV_V'],n);cp=2*indexed(f['cells/SV_P'],n);wallcp=2*indexed(f['faces/SV_P'],nf)[wallids]
 return dict(xy=xy,tri=tri,U=np.column_stack((u,v)),Cp=cp,wall=xy[wall],wall_midpoint=xy[wall].mean(axis=1),wall_Cp=wallcp,adjacent_Cp=cp[owners],inputs=[case,data],wall_pressure_definition='Native Fluent face SV_P, gauge Pa divided by rhoU^2/2=.5Pa; adjacent-cell p shown separately.')

def openfoam(case,endpoint):
 g=top.load_case_geometry(case);n=len(g['triangles']);end=case/str(endpoint);u=top.internal_field(end/'U','vector',n);cp=2*top.internal_field(end/'p','scalar',n);mesh=case/'0/polyMesh';points=top.read_points(mesh/'points');sizes,vertices=top.read_faces(mesh/'faces');owner=top.read_labels(mesh/'owner');raw=av.block(av.dictionary(mesh/'boundary'),'airfoil');first=int(av.scalar(raw,'startFace'));count=int(av.scalar(raw,'nFaces'));ids=np.arange(first,first+count);offset=np.r_[0,np.cumsum(sizes)];facecentres=np.array([points[vertices[offset[k]:offset[k+1]],:2].mean(axis=0) for k in ids]);assert np.all(sizes[ids]==4)
 source=top.read_source_mesh(case/'input-airfoil.msh');wall=source['xy'][source['edges'][source['zones']==4]];mid=wall.mean(axis=1)
 native_to_source,metrics=top.coordinate_bijection(mid,facecentres,1e-10)
 values=np.empty(count);values[native_to_source]=cp[owner[ids]]
 return dict(xy=g['xy'],tri=g['triangles'],U=u[:,:2],Cp=cp,wall=wall,wall_midpoint=mid,wall_Cp=values,adjacent_Cp=values.copy(),inputs=g['inputs']+[end/f for f in ('U','p','phi')],wall_pressure_definition='Native zeroGradient evaluated FV wall pressure equals owner-cell kinematic p; rho=1. Therefore its wall and adjacent values coincide, with no independent wall-pressure reconstruction.')

def airfoil_frame(xy):
 a=np.deg2rad(5);return (xy-[.5,0])@np.array([[np.cos(a),np.sin(a)],[-np.sin(a),np.cos(a)]])+[.5,0]

def check_fields(f):
 n=len(f['tri']);assert f['U'].shape==(n,2) and f['Cp'].shape==(n,)
 for k in ('xy','U','Cp','wall','wall_midpoint','wall_Cp','adjacent_Cp'):assert np.isfinite(f[k]).all(),k
 pts=f['xy'][f['tri']];a=pts[:,1]-pts[:,0];b=pts[:,2]-pts[:,0];cross=a[:,0]*b[:,1]-a[:,1]*b[:,0];assert np.all(cross!=0)
 f['tri']=f['tri'].copy();f['tri'][cross<0]=f['tri'][cross<0][:,[0,2,1]]
 f['speed']=np.linalg.norm(f['U'],axis=1);f['wall_frame']=airfoil_frame(f['wall_midpoint']);return f
