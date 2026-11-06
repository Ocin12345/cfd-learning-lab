"""Shared 2304-quad heated-channel geometry; no solver launch."""
from pathlib import Path
import argparse,hashlib,json,shutil
import numpy as np
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);args=ap.parse_args();out=args.out.resolve()
 if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(Path(__file__).resolve().parents[1]):raise ValueError('Fresh external E: mesh required')
 if shutil.disk_usage(out.anchor).free<10*1024**3:raise RuntimeError('Require10GiB free E:')
 nx,ny=96,24;L,H=1.,.1;dx,dy=L/nx,H/ny
 cells=[(i,j) for j in range(ny) for i in range(nx)]
 nodes=np.array([(i*dx,j*dy) for j in range(ny+1) for i in range(nx+1)])
 def node(i,j):return j*(nx+1)+i
 quads=np.array([[node(i,j),node(i+1,j),node(i+1,j+1),node(i,j+1)] for i,j in cells],dtype=np.int64)
 poly=nodes[quads];origin=poly[:,:1];loc=poly-origin;nxt=np.roll(loc,-1,axis=1)
 cross=loc[:,:,0]*nxt[:,:,1]-nxt[:,:,0]*loc[:,:,1];areas=.5*cross.sum(axis=1)
 centres=origin[:,0]+((loc+nxt)*cross[:,:,None]).sum(axis=1)/(6*areas[:,None])
 assert nodes.shape==(2425,2) and quads.shape==(2304,4) and np.all(areas>0)
 assert np.max(np.abs(areas-dx*dy))<1e-17 and abs(areas.sum()-.1)<1e-14
 assert np.max(np.abs(centres-(np.array(cells)+.5)*[dx,dy]))<1e-12
 faces={}
 for cid,q in enumerate(quads,1):
  for aa,bb in zip(q,np.roll(q,-1)):
   key=tuple(sorted((int(aa),int(bb))))
   if key in faces:assert faces[key][3]==0;faces[key][3]=cid
   else:faces[key]=[int(aa)+1,int(bb)+1,cid,0]
 groups={z:[] for z in (3,4,5,6,7)};orth=[]
 for aa,bb,c0,c1 in faces.values():
  v0,v1=nodes[aa-1],nodes[bb-1];mid=.5*(v0+v1);edge=v1-v0;n=np.array([edge[1],-edge[0]])
  if c1:zone=3;delta=centres[c1-1]-centres[c0-1]
  else:
   delta=mid-centres[c0-1]
   if aa-1 in range(0,2425,97) and bb-1 in range(0,2425,97):zone=4
   elif (aa-1)%97==(bb-1)%97==96:zone=5
   elif (aa-1)//97==(bb-1)//97==24:zone=6
   elif (aa-1)//97==(bb-1)//97==0:zone=7
   else:raise ValueError('Unclassified channel boundary')
  cosine=float(np.dot(n,delta)/(np.linalg.norm(n)*np.linalg.norm(delta)));assert cosine>0;orth.append(cosine)
  groups[zone].append((aa,bb,c0,c1))
 assert len(faces)==4728 and len(nodes)-len(faces)+len(cells)==1
 assert [len(groups[z]) for z in groups]==[4488,24,24,96,96]
 out.mkdir(parents=True);mesh=out/'heated-channel.msh'
 with mesh.open('x',encoding='ascii',newline='\n') as f:
  f.write('(0 "Heated channel L1 H0.1;2304 shared quads; source geometry only")\n(2 2)\n')
  f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
  for x,y in nodes:f.write(f'{x:.17g} {y:.17g}\n')
  f.write('))\n');f.write(f'(12 (0 1 {len(cells):x} 0))\n(12 (2 1 {len(cells):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
  first=1
  for zone,rows in groups.items():
   bc={3:2,4:10,5:5,6:3,7:3}[zone];f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
   for row in rows:f.write(' '.join(f'{v:x}' for v in row)+'\n')
   f.write('))\n');first+=len(rows)
  for z,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'velocity-inlet','inlet'),(5,'pressure-outlet','outlet'),(6,'wall','top'),(7,'wall','bottom')]:f.write(f'(45 ({z} {kind} {name})())\n')
 np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=quads,centres=centres,areas=areas,grid_shape=np.array([24,96]),ix=np.array([i for i,j in cells]),iy=np.array([j for i,j in cells]),x_edges=np.arange(97)*dx,y_edges=np.arange(25)*dy)
 setup=dict(model='Transient heated-wall laminar channel on frozen computed flow',cells=2304,nx=nx,ny=ny,length_m=L,height_m=H,dx_m=dx,dy_m=dy,domain_x_bounds_m=[0.,1.],domain_y_bounds_m=[0.,.1],density_kg_m3=1.,dynamic_viscosity_Pa_s=.01,kinematic_viscosity_m2_s=.01,specific_heat_J_kg_K=1000.,thermal_conductivity_W_m_K=10.,thermal_diffusivity_m2_s=.01,inlet_velocity_m_s=[1.,0.],Re_H=10.,Re_Dh=20.,Pr=1.,hydraulic_diameter_m=.2,reference_length_m=.2,reference_velocity_m_s=1.,reference_depth_m=1.,openfoam_extrusion_depth_m=.01,operating_pressure_Pa=0.,outlet_gauge_pressure_Pa=0.,initial_velocity_m_s=[0.,0.],initial_pressure_Pa=0.,initial_temperature_K=300.,inlet_temperature_K=300.,wall_temperature_K=350.,time_step_s=.005,time_steps=600,final_time_s=3.,boundary_names=['inlet','outlet','top','bottom'],boundary_face_counts=dict(inlet=24,outlet=24,top=96,bottom=96),boundary_types=dict(inlet='velocity-inlet',outlet='pressure-outlet',top='wall',bottom='wall'),geometry_checks=dict(nodes=2425,faces=4728,interior_faces=4488,positive_areas=True,total_area_m2=float(areas.sum()),minimum_area_m2=float(areas.min()),euler_characteristic=1,orthogonality_proxy_min=float(min(orth))),mesh_sha256=sha(mesh),geometry_npz_sha256=sha(out/'geometry.npz'),generated_by='Codex',physical_validation=False,mesh_independence=False,time_step_independence=False)
 (out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n',encoding='utf-8');(out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
 (out/'source-manifest.json').write_text(json.dumps({p.name:sha(p) for p in out.iterdir() if p.is_file()},indent=2)+'\n',encoding='utf-8')
 print(json.dumps(dict(output=str(out),cells=2304,nodes=2425,edges=4728,area_m2=float(areas.sum()),mesh_sha256=sha(mesh),geometry_sha256=sha(out/'geometry.npz')),indent=2))
if __name__=='__main__':main()
