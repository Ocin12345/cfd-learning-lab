"""Shared square quadrilateral cavity mesh for both native solvers."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from cavity_buoyancy_common_v001 import L,N,ROOT,specification,write

def main(out):
    assert not out.exists();out.mkdir(parents=True)
    x=np.linspace(0,L,N+1); nodes=np.array([(xx,yy) for yy in x for xx in x])
    def nid(i,j):return j*(N+1)+i
    cells=np.array([[nid(i,j),nid(i+1,j),nid(i+1,j+1),nid(i,j+1)] for j in range(N) for i in range(N)])
    centres=nodes[cells].mean(axis=1);areas=np.full(N*N,(L/N)**2);faces={}
    for c,q in enumerate(cells,1):
        for a,b in zip(q,np.roll(q,-1)):
            key=tuple(sorted((int(a),int(b))))
            if key in faces:assert faces[key][3]==0;faces[key][3]=c
            else:faces[key]=[int(a)+1,int(b)+1,c,0]
    groups={z:[] for z in [3,4,5,6,7]}
    for a,b,c0,c1 in faces.values():
        v0,v1=nodes[a-1],nodes[b-1];mid=(v0+v1)/2;edge=v1-v0;normal=np.array([edge[1],-edge[0]])
        if c1:z=3;delta=centres[c1-1]-centres[c0-1]
        else:
            delta=mid-centres[c0-1]
            if v0[1]==v1[1]==L:z=4
            elif v0[1]==v1[1]==0:z=5
            elif v0[0]==v1[0]==0:z=6
            else:assert v0[0]==v1[0]==L;z=7
        assert np.dot(normal,delta)>0;groups[z].append((a,b,c0,c1))
    assert len(nodes)-len(faces)+len(cells)==1
    with (out/'cavity.msh').open('x',encoding='ascii',newline='\n') as f:
        f.write('(0 "Shared stationary ideal-gas cavity; no computed flow")\n(2 2)\n')
        f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
        for xx,yy in nodes:f.write(f'{xx:.17g} {yy:.17g}\n')
        f.write(f'))\n(12 (0 1 {len(cells):x} 0))\n(12 (2 1 {len(cells):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
        first=1
        for z,rows in groups.items():
            f.write(f'(13 ({z:x} {first:x} {first+len(rows)-1:x} {2 if z==3 else 3:x} 2)(\n')
            for row in rows:f.write(' '.join(f'{v:x}' for v in row)+'\n')
            f.write('))\n');first+=len(rows)
        for z,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','top'),(5,'wall','bottom'),(6,'wall','left'),(7,'wall','right')]:f.write(f'(45 ({z} {kind} {name})())\n')
    np.savez_compressed(out/'geometry.npz',nodes=nodes,quads=cells,centres=centres,areas=areas,grid_shape=[N,N],x_edges=x,y_edges=x,ix=np.tile(np.arange(N),N),iy=np.repeat(np.arange(N),N))
    write(out/'setup.json',dict(**specification(),mesh_sha256=hashlib.sha256((out/'cavity.msh').read_bytes()).hexdigest(),positive_cell_areas=True,minimum_orthogonal_quality_expected=1.,total_area_m2=float(areas.sum()),boundary_faces={str(z):len(v) for z,v in groups.items()}))
    (out/'executed-mesher.py').write_bytes(Path(__file__).read_bytes());print(out)
if __name__=='__main__':main(Path(sys.argv[1]))
