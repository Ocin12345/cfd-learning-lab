"""Shared 2000-quad rectangular solid-conduction slab source mesh; no solver launch."""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
 out=a.out.resolve();repo=Path(__file__).resolve().parents[1]
 if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo):raise ValueError('Fresh external E: mesh required')
 if shutil.disk_usage(out.anchor).free<10*1024**3:raise RuntimeError('Require10GiB free E:')
 cells=[(i,j) for j in range(20) for i in range(100)]
 vertices=sorted({v for i,j in cells for v in ((i,j),(i+1,j),(i+1,j+1),(i,j+1))},key=lambda v:(v[1],v[0]))
 number={v:k for k,v in enumerate(vertices)};nodes=np.array(vertices,dtype=float)*.001
 quads=np.array([[number[v] for v in ((i,j),(i+1,j),(i+1,j+1),(i,j+1))] for i,j in cells],dtype=np.int64)
 polygons=nodes[quads];nxt=np.roll(polygons,-1,axis=1)
 cross=polygons[:,:,0]*nxt[:,:,1]-nxt[:,:,0]*polygons[:,:,1]
 areas=.5*cross.sum(axis=1);centres=((polygons+nxt)*cross[:,:,None]).sum(axis=1)/(6*areas[:,None])
 assert len(cells)==2000 and len(nodes)==2121 and np.all(areas>0)
 assert np.allclose(areas,.000001,atol=1e-16,rtol=0) and abs(areas.sum()-.002)<1e-14
 assert np.allclose(centres,(np.array(cells)+.5)*.001,atol=1e-12,rtol=0)
 faces={}
 for cid,q in enumerate(quads,1):
  for aa,bb in zip(q,np.roll(q,-1)):
   key=tuple(sorted((int(aa),int(bb))))
   if key in faces:assert faces[key][3]==0;faces[key][3]=cid
   else:faces[key]=[int(aa)+1,int(bb)+1,cid,0]
 groups={z:[] for z in (3,4,5,6,7)};orth=[]
 for aa,bb,c0,c1 in faces.values():
  v0,v1=nodes[aa-1],nodes[bb-1];mid=.5*(v0+v1);edge=v1-v0;norm=np.array([edge[1],-edge[0]])
  if c1:zone=3;delta=centres[c1-1]-centres[c0-1]
  else:
   delta=mid-centres[c0-1];iv0,iv1=vertices[aa-1],vertices[bb-1]
   if iv0[0]==iv1[0]==0:zone=4
   elif iv0[0]==iv1[0]==100:zone=5
   elif iv0[1]==iv1[1]==0:zone=6
   elif iv0[1]==iv1[1]==20:zone=7
   else:raise ValueError('Unclassified slab boundary')
  cosine=float(np.dot(norm,delta)/(np.linalg.norm(norm)*np.linalg.norm(delta)));assert cosine>0;orth.append(cosine)
  groups[zone].append((aa,bb,c0,c1))
 assert len(faces)==4120 and len(nodes)-len(faces)+len(cells)==1
 assert [len(groups[z]) for z in (3,4,5,6,7)]==[3880,20,20,100,100]
 out.mkdir(parents=True);mesh=out/'slab.msh'
 with mesh.open('w',encoding='ascii',newline='\n') as f:
  f.write('(0 "Solid-conduction slab; 2000 uniform quads; import zone converted to solid before solving")\n(2 2)\n')
  f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
  for xx,yy in nodes:f.write(f'{xx:.17g} {yy:.17g}\n')
  f.write('))\n');f.write(f'(12 (0 1 {len(cells):x} 0))\n(12 (2 1 {len(cells):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
  first=1
  for zone in (3,4,5,6,7):
   rows=groups[zone];bc=2 if zone==3 else 3
   f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
   for row in rows:f.write(' '.join(f'{v:x}' for v in row)+'\n')
   f.write('))\n');first+=len(rows)
  for z,kind,name in [(2,'fluid','slab'),(3,'interior','interior'),(4,'wall','hot'),(5,'wall','cold'),(6,'wall','insulated_bottom'),(7,'wall','insulated_top')]:f.write(f'(45 ({z} {kind} {name})())\n')
 np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=quads,centres=centres,areas=areas,grid_shape=np.array([20,100]),ix=np.array([i for i,j in cells]),iy=np.array([j for i,j in cells]),x_edges=np.arange(101)*.001,y_edges=np.arange(21)*.001)
 setup=dict(model='Transient solid conduction slab',cells=2000,nx=100,ny=20,length_m=.1,height_m=.02,dx_m=.001,dy_m=.001,
  density_kg_m3=1000.,specific_heat_J_kg_K=1000.,thermal_conductivity_W_m_K=10.,thermal_diffusivity_m2_s=1e-5,
  initial_temperature_K=300.,hot_temperature_K=400.,cold_temperature_K=300.,time_step_s=2.,time_steps=100,final_time_s=200.,
  reference_depth_m=1.,openfoam_extrusion_depth_m=.01,cell_zone_name='slab',import_zone_type='fluid',solved_zone_type='solid',
  boundary_names=['hot','cold','insulated_bottom','insulated_top'],boundary_face_counts=dict(hot=20,cold=20,insulated_bottom=100,insulated_top=100),
  geometry_checks=dict(nodes=len(nodes),faces=len(faces),interior_faces=3880,positive_areas=True,total_area_m2=float(areas.sum()),
   minimum_area_m2=float(areas.min()),euler_characteristic=1,orthogonality_proxy_min=float(min(orth))),
  mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),generated_by='Codex',physical_validation=False,mesh_independence=False,time_step_independence=False)
 (out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n',encoding='utf-8')
 (out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
 (out/'source-manifest.json').write_text(json.dumps({p.name:sha(p) for p in out.iterdir() if p.is_file()},indent=2)+'\n',encoding='utf-8')
 print(json.dumps(dict(out=str(out),cells=2000,nodes=len(nodes),faces=len(faces),area=float(areas.sum()),mesh_sha256=sha(mesh),geometry_sha256=sha(out/'geometry.npz')),indent=2))
if __name__=='__main__':main()
