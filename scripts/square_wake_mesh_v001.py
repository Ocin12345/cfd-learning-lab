"""Shared15900quad square-obstacle channel source mesh, one genuine hole; no solver launch."""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 out=a.out.resolve();repo=Path(__file__).resolve().parents[1]
 if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo):raise ValueError('Fresh external E: mesh required')
 if shutil.disk_usage(out.anchor).free<10*1024**3:raise RuntimeError('Require10GiB free E:')
 cells=[(i,j) for j in range(-40,40) for i in range(200) if not(45<=i<55 and -5<=j<5)]
 vertices=sorted({v for i,j in cells for v in ((i,j),(i+1,j),(i+1,j+1),(i,j+1))},key=lambda v:(v[1],v[0]))
 number={v:k for k,v in enumerate(vertices)};nodes=np.array(vertices,dtype=float)*.01
 quads=np.array([[number[v] for v in ((i,j),(i+1,j),(i+1,j+1),(i,j+1))] for i,j in cells],dtype=np.int64)
 polygons=nodes[quads];nxt=np.roll(polygons,-1,axis=1)
 cross=polygons[:,:,0]*nxt[:,:,1]-nxt[:,:,0]*polygons[:,:,1]
 areas=.5*cross.sum(axis=1);centres=((polygons+nxt)*cross[:,:,None]).sum(axis=1)/(6*areas[:,None])
 assert len(cells)==15900 and len(nodes)==16200 and np.all(areas>0)
 assert np.allclose(areas,.0001,atol=1e-16,rtol=0) and abs(areas.sum()-1.59)<1e-14
 assert np.allclose(centres,(np.array(cells)+.5)*.01,atol=1e-12,rtol=0)
 faces={}
 for cid,q in enumerate(quads,1):
  for aa,bb in zip(q,np.roll(q,-1)):
   key=tuple(sorted((int(aa),int(bb))))
   if key in faces:assert faces[key][3]==0;faces[key][3]=cid
   else:faces[key]=[int(aa)+1,int(bb)+1,cid,0]
 groups={z:[] for z in (3,4,5,6,7,8)};orth=[]
 for aa,bb,c0,c1 in faces.values():
  v0,v1=nodes[aa-1],nodes[bb-1];mid=.5*(v0+v1);edge=v1-v0;norm=np.array([edge[1],-edge[0]])
  if c1:zone=3;delta=centres[c1-1]-centres[c0-1]
  else:
   delta=mid-centres[c0-1];iv0,iv1=vertices[aa-1],vertices[bb-1]
   if iv0[0]==iv1[0]==0:zone=4
   elif iv0[0]==iv1[0]==200:zone=5
   elif iv0[1]==iv1[1]==40:zone=6
   elif iv0[1]==iv1[1]==-40:zone=7
   else:
    assert (iv0[0]==iv1[0] and iv0[0] in (45,55) and -5<=iv0[1]<=5 and -5<=iv1[1]<=5) or (iv0[1]==iv1[1] and iv0[1] in (-5,5) and 45<=iv0[0]<=55 and 45<=iv1[0]<=55)
    zone=8
  cosine=float(np.dot(norm,delta)/(np.linalg.norm(norm)*np.linalg.norm(delta)));assert cosine>0;orth.append(cosine)
  groups[zone].append((aa,bb,c0,c1))
 assert len(faces)==32100 and len(nodes)-len(faces)+len(cells)==0
 assert [len(groups[z]) for z in (3,4,5,6,7,8)]==[31500,80,80,200,200,40]
 out.mkdir(parents=True);mesh=out/'square-wake.msh'
 with mesh.open('w',encoding='ascii',newline='\n') as f:
  f.write('(0 "Square-obstacle Re20 channel;15900 uniform fluid quads; genuine square hole; source geometry only")\n(2 2)\n')
  f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
  for xx,yy in nodes:f.write(f'{xx:.17g} {yy:.17g}\n')
  f.write('))\n');f.write(f'(12 (0 1 {len(cells):x} 0))\n(12 (2 1 {len(cells):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
  first=1
  for zone in (3,4,5,6,7,8):
   rows=groups[zone];bc={3:2,4:10,5:5,6:7,7:7,8:3}[zone]
   f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
   for row in rows:f.write(' '.join(f'{v:x}' for v in row)+'\n')
   f.write('))\n');first+=len(rows)
  for z,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'velocity-inlet','inlet'),(5,'pressure-outlet','outlet'),(6,'symmetry','top'),(7,'symmetry','bottom'),(8,'wall','body')]:f.write(f'(45 ({z} {kind} {name})())\n')
 np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=quads,centres=centres,areas=areas,grid_shape=np.array([80,200]),ix=np.array([i for i,j in cells]),iy=np.array([j+40 for i,j in cells]),x_edges=np.arange(201)*.01,y_edges=np.arange(-40,41)*.01)
 setup=dict(model='Steady laminar square-obstacle confined wake',cells=15900,nx=200,ny=80,dx_m=.01,dy_m=.01,
  density_kg_m3=1.,dynamic_viscosity_Pa_s=.0005,kinematic_viscosity_m2_s=.0005,inlet_velocity_m_s=[.1,0.],outlet_gauge_pressure_Pa=0.,
  Re_body=20.,body_side_m=.1,body_x_bounds_m=[.45,.55],body_y_bounds_m=[-.05,.05],domain_x_bounds_m=[0.,2.],domain_y_bounds_m=[-.4,.4],
  reference_length_m=.1,reference_velocity_m_s=.1,reference_area_m2=.1,reference_depth_m=1.,openfoam_extrusion_depth_m=.01,
  operating_pressure_Pa=0.,initial_velocity_m_s=[0.,0.],initial_pressure_Pa=0.,expected_flow_per_unit_depth_m2_s=dict(inlet=.08,outlet=.08),
  boundary_names=['inlet','outlet','top','bottom','body'],boundary_types=dict(inlet='velocity-inlet',outlet='pressure-outlet',top='symmetry',bottom='symmetry',body='wall'),
  boundary_face_counts=dict(inlet=80,outlet=80,top=200,bottom=200,body=40),
  geometry_checks=dict(nodes=len(nodes),faces=len(faces),interior_faces=31500,positive_areas=True,total_area_m2=float(areas.sum()),
   minimum_area_m2=float(areas.min()),euler_characteristic=0,holes=1,orthogonality_proxy_min=float(min(orth))),
  mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),generated_by='Codex',physical_validation=False,mesh_independence=False,
  confined_domain=True,domain_height_over_bodyD=8.,blockage_ratio=.125,transient_shedding_simulated=False)
 (out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n',encoding='utf-8')
 (out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
 (out/'source-manifest.json').write_text(json.dumps({p.name:sha(p) for p in out.iterdir() if p.is_file()},indent=2)+'\n',encoding='utf-8')
 print(json.dumps(dict(out=str(out),cells=15900,nodes=len(nodes),faces=len(faces),area=float(areas.sum()),mesh_sha256=sha(mesh),geometry_sha256=sha(out/'geometry.npz')),indent=2))
if __name__=='__main__':main()
