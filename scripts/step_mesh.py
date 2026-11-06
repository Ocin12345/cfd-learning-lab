"""Small conforming stair-step mesh, shared by Fluent and OpenFOAM."""
from pathlib import Path
import argparse,hashlib,json,shutil
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
out=a.out.resolve();repo=Path(__file__).resolve().parents[1]
if out.exists() or out.is_relative_to(repo):raise ValueError('Choose a fresh external mesh directory')
if shutil.disk_usage(out.anchor).free<10*1024**3:raise RuntimeError('Require10GiB free')
nx,ny,dx,dy,x0=240,40,.01,.005,-.4
slots=[(i,j) for j in range(ny) for i in range(nx) if i>=40 or j>=20]
used=sorted({(i+di,j+dj) for i,j in slots for di,dj in [(0,0),(1,0),(1,1),(0,1)]},key=lambda q:(q[1],q[0]))
ids={q:k+1 for k,q in enumerate(used)}
nodes=[(x0+i*dx,j*dy) for i,j in used]
faces={}
for cid,(i,j) in enumerate(slots,1):
 v=[ids[(i,j)],ids[(i+1,j)],ids[(i+1,j+1)],ids[(i,j+1)]]
 xy=[nodes[k-1] for k in v]
 area=.5*sum(xy[k][0]*xy[(k+1)%4][1]-xy[(k+1)%4][0]*xy[k][1] for k in range(4))
 assert area>0 and abs(area-dx*dy)<1e-15
 for n0,n1 in zip(v,v[1:]+v[:1]):
  key=tuple(sorted((n0,n1)))
  if key in faces:assert faces[key][3]==0;faces[key][3]=cid
  else:faces[key]=[n0,n1,cid,0]
groups={3:[],4:[],5:[],6:[]}
for n0,n1,c0,c1 in faces.values():
 xa,ya=nodes[n0-1];xb,yb=nodes[n1-1]
 if c1:zone=3
 elif abs(xa-x0)<1e-12 and abs(xb-x0)<1e-12:zone=5
 elif abs(xa-2)<1e-12 and abs(xb-2)<1e-12:zone=6
 else:zone=4
 groups[zone].append((n0,n1,c0,c1))
assert len(slots)==8800 and len(nodes)-len(faces)+len(slots)==1
assert len(groups[5])==20 and len(groups[6])==40 and len(groups[4])==500
assert len(nodes)==9081 and len(faces)==17880
out.mkdir(parents=True);mesh=out/'step.msh'
with mesh.open('w',encoding='ascii',newline='\n') as f:
 f.write('(0 "Laminar backward-facing step; h0.1m, expansion2:1")\n(2 2)\n')
 f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
 for x,y in nodes:f.write(f'{x:.16g} {y:.16g}\n')
 f.write('))\n');f.write(f'(12 (0 1 {len(slots):x} 0))\n(12 (2 1 {len(slots):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
 first=1
 for zone,bc in [(3,2),(4,3),(5,10),(6,5)]:
  rows=groups[zone];f.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
  for row in rows:f.write(' '.join(f'{k:x}' for k in row)+'\n')
  f.write('))\n');first+=len(rows)
 for zone,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','walls'),(5,'velocity-inlet','inlet'),(6,'pressure-outlet','outlet')]:f.write(f'(45 ({zone} {kind} {name})())\n')
setup=dict(model='2D laminar backward-facing step; uniform developing inlet',nx=nx,ny=ny,cells=len(slots),x_min_m=x0,x_max_m=2.,y_min_m=0.,length_m=2.4,height_m=.2,step_height_m=.1,upstream_height_m=.1,upstream_length_m=.4,downstream_length_m=2.,step_column=40,step_row=20,dx_m=dx,dy_m=dy,reference_depth_m=1.,density_kg_m3=1.,dynamic_viscosity_Pa_s=.001,kinematic_viscosity_m2_s=.001,inlet_velocity_m_s=1.,Re_h=100.,Re_upstream_Dh=200.,outlet_gauge_pressure_Pa=0.,positive_cell_areas=True,euler_characteristic=1,unused_nodes=0,boundary_faces=dict(walls=500,inlet=20,outlet=40),mesh_sha256=hashlib.sha256(mesh.read_bytes()).hexdigest(),physical_validation=False,generated_by='Codex')
(out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n',encoding='utf-8')
(out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
print(json.dumps(setup,indent=2))
