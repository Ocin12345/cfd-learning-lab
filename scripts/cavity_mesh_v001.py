"""Shared uniform full-square cavity mesh. Does not launch a solver."""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
 out=args.out.resolve();repo=Path(__file__).resolve().parents[1]
 if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo):raise ValueError('Fresh external E: mesh output required')
 if shutil.disk_usage(out.anchor).free<10*1024**3:raise RuntimeError('Require10GiB free E:')
 nx=ny=80;x=np.linspace(0.,1.,nx+1);y=np.linspace(0.,1.,ny+1)
 nodes=np.array([(xx,yy) for yy in y for xx in x])
 def nid(i,j):return j*(nx+1)+i
 quads=np.array([[nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)] for j in range(ny) for i in range(nx)],dtype=np.int64)
 polygons=nodes[quads];after=np.roll(polygons,-1,axis=1)
 cross=polygons[:,:,0]*after[:,:,1]-after[:,:,0]*polygons[:,:,1]
 areas=.5*cross.sum(axis=1);centres=((polygons+after)*cross[:,:,None]).sum(axis=1)/(6*areas[:,None])
 assert len(nodes)==6561 and len(quads)==6400 and np.all(areas>0)
 assert np.allclose(areas,1/6400,rtol=0,atol=2e-16) and abs(areas.sum()-1)<1e-14
 faces={}
 for cid,q in enumerate(quads,1):
  for a,b in zip(q,np.roll(q,-1)):
   key=tuple(sorted((int(a),int(b))))
   if key in faces:assert faces[key][3]==0;faces[key][3]=cid
   else:faces[key]=[int(a)+1,int(b)+1,cid,0]
 groups={i:[] for i in (3,4,5,6,7)};orth=[]
 for a,b,c0,c1 in faces.values():
  v0,v1=nodes[a-1],nodes[b-1];mid=.5*(v0+v1);edge=v1-v0;norm=np.array([edge[1],-edge[0]])
  if c1:zone=3;delta=centres[c1-1]-centres[c0-1]
  else:
   delta=mid-centres[c0-1]
   if v0[1]==v1[1]==1:zone=4
   elif v0[1]==v1[1]==0:zone=5
   elif v0[0]==v1[0]==0:zone=6
   else:assert v0[0]==v1[0]==1;zone=7
  cos=float(np.dot(norm,delta)/(np.linalg.norm(norm)*np.linalg.norm(delta)));assert cos>0;orth.append(cos)
  groups[zone].append((a,b,c0,c1))
 assert len(faces)==12960 and len(nodes)-len(faces)+len(quads)==1
 assert [len(groups[z]) for z in (3,4,5,6,7)]==[12640,80,80,80,80]
 out.mkdir(parents=True);mesh=out/'cavity.msh'
 with mesh.open('w',encoding='ascii',newline='\n') as f:
  f.write('(0 "Full square lid-driven cavity; uniform80x80quads; no computed flow")\n(2 2)\n')
  f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
  for xx,yy in nodes:f.write(f'{xx:.17g} {yy:.17g}\n')
  f.write('))\n');f.write(f'(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
  first=1
  for zone in (3,4,5,6,7):
   rows=groups[zone];bc=2 if zone==3 else 3;f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
   for row in rows:f.write(' '.join(f'{v:x}' for v in row)+'\n')
   f.write('))\n');first+=len(rows)
  for zone,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','lid'),(5,'wall','bottom'),(6,'wall','left'),(7,'wall','right')]:f.write(f'(45 ({zone} {kind} {name})())\n')
 np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=quads,centres=centres,areas=areas,grid_shape=np.array([ny,nx]),ix=np.tile(np.arange(nx),ny),iy=np.repeat(np.arange(ny),nx),x_edges=x,y_edges=y)
 setup=dict(model='Full-square lid-driven cavity, steady incompressible laminar Re100',nx=nx,ny=ny,cells=6400,length_m=1.,height_m=1.,density_kg_m3=1.,dynamic_viscosity_Pa_s=.01,kinematic_viscosity_m2_s=.01,lid_velocity_m_s=[1.,0.],Re=100.,reference_depth_m=1.,openfoam_extrusion_depth_m=.01,operating_pressure_Pa=0.,initial_velocity_m_s=[0.,0.],initial_pressure_Pa=0.,boundary_names=['lid','bottom','left','right'],corner_lid_velocity_discontinuity=True,geometry_checks=dict(nodes=6561,faces=12960,interior_faces=12640,wall_faces={k:80 for k in ('lid','bottom','left','right')},positive_areas=True,total_area_m2=float(areas.sum()),minimum_area_m2=float(areas.min()),euler_characteristic=1,orthogonality_proxy_min=float(min(orth))),mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),physical_validation=False,generated_by='Codex')
 (out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n',encoding='utf-8')
 (out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
 manifest={p.name:sha(p) for p in out.iterdir() if p.is_file()};(out/'source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
 print(json.dumps(dict(out=str(out),cells=6400,nodes=len(nodes),faces=len(faces),area=float(areas.sum()),mesh_sha256=sha(mesh),geometry_sha256=sha(out/'geometry.npz')),indent=2))
if __name__=='__main__':main()