"""Small shared planar-channel mesh; native outputs remain outside Git."""
from pathlib import Path
import argparse, hashlib, json, shutil

p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
out=a.out.resolve();repo=Path(__file__).resolve().parents[1]
if out.exists() or out.is_relative_to(repo):raise ValueError('Use a fresh external mesh directory')
if shutil.disk_usage(out.anchor).free<10*1024**3:raise RuntimeError('Require10GiB free output space')
nx,ny,L,H=240,32,1.,.1
nodes=[(L*i/nx,H*j/ny) for j in range(ny+1) for i in range(nx+1)]
def nid(i,j):return j*(nx+1)+i+1
faces={}
for j in range(ny):
    for i in range(nx):
        cid=j*nx+i+1;v=[nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)]
        xy=[nodes[k-1] for k in v]
        area=.5*sum(xy[k][0]*xy[(k+1)%4][1]-xy[(k+1)%4][0]*xy[k][1] for k in range(4))
        assert area>0 and abs(area-L*H/(nx*ny))<1e-15
        for v0,v1 in zip(v,v[1:]+v[:1]):
            key=tuple(sorted((v0,v1)))
            if key in faces:
                assert faces[key][3]==0;faces[key][3]=cid
            else:faces[key]=[v0,v1,cid,0]
groups={3:[],4:[],5:[],6:[]}
for v0,v1,c0,c1 in faces.values():
    x0,y0=nodes[v0-1];x1,y1=nodes[v1-1]
    if c1:zone=3
    elif x0==x1==0:zone=5
    elif x0==x1==L:zone=6
    else:
        assert y0==y1 and y0 in (0,H);zone=4
    groups[zone].append((v0,v1,c0,c1))
assert len(nodes)-len(faces)+nx*ny==1
assert [len(groups[z]) for z in (4,5,6)]==[2*nx,ny,ny]
out.mkdir(parents=True)
mesh=out/'channel.msh'
with mesh.open('w',encoding='ascii',newline='\n') as f:
    f.write('(0 "Planar channel L1m H0.1m; uniform quad mesh")\n(2 2)\n')
    f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
    for x,y in nodes:f.write(f'{x:.16g} {y:.16g}\n')
    f.write('))\n');f.write(f'(12 (0 1 {nx*ny:x} 0))\n(12 (2 1 {nx*ny:x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
    start=1
    for zone,bc in [(3,2),(4,3),(5,10),(6,5)]:
        rows=groups[zone];f.write(f'(13 ({zone:x} {start:x} {start+len(rows)-1:x} {bc:x} 2)(\n')
        for row in rows:f.write(' '.join(f'{k:x}' for k in row)+'\n')
        f.write('))\n');start+=len(rows)
    for zone,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','walls'),(5,'velocity-inlet','inlet'),(6,'pressure-outlet','outlet')]:
        f.write(f'(45 ({zone} {kind} {name})())\n')
setup=dict(model='Developing laminar flow between parallel plates; 2D planar, not circular pipe',nx=nx,ny=ny,cells=nx*ny,length_m=L,height_m=H,reference_depth_m=1.,density_kg_m3=1.,dynamic_viscosity_Pa_s=.01,kinematic_viscosity_m2_s=.01,inlet_velocity_m_s=1.,Re_H=10.,Re_Dh=20.,outlet_gauge_pressure_Pa=0.,geometry_checks=dict(positive_areas=True,euler_characteristic=1,boundary_face_counts={str(z):len(groups[z]) for z in (4,5,6)}),mesh_sha256=hashlib.sha256(mesh.read_bytes()).hexdigest(),generated_by='Codex',physical_validation=False)
(out/'setup.json').write_text(json.dumps(setup,indent=2)+'\n',encoding='utf-8')
(out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
print(json.dumps(setup,indent=2))
