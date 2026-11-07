"""Tilt the 109,600-triangle airfoil mesh while preserving connectivity and the far field."""
import argparse, hashlib, json, re
from pathlib import Path
import numpy as np

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--angle",type=float,choices=[5.],default=5.)
    a=parser.parse_args()
    assert not a.out.exists() and a.out.drive.upper()=="E:"
    raw=a.source.read_bytes()
    src_hash=hashlib.sha256(raw).hexdigest()
    lines=raw.decode("ascii").splitlines()
    start=next(i for i,x in enumerate(lines) if x.startswith("(10 (1 "))+1
    stop=lines.index("))",start)
    xy=np.array([list(map(float,x.split())) for x in lines[start:stop]])
    setup=json.loads(a.source.with_name("setup.json").read_text())
    assert setup["angle_degrees"]==0 and setup["cells"]==109600 and len(xy)==setup["nodes"]
    sections={}
    for j,line in enumerate(lines):
        m=re.match(r"^\(13 \(([3456]) ",line)
        if m:
            k=lines.index("))",j+1)
            sections[int(m[1])]=np.array([[int(v,16) for v in s.split()] for s in lines[j+1:k]])
    assert set(sections)=={3,4,5,6}
    walls=np.unique(sections[4][:,:2])-1
    far=np.unique(np.vstack((sections[5],sections[6]))[:,:2])-1
    p=xy-[.5,0.]
    r=np.linalg.norm(p,axis=1)
    z=np.clip((r-1.25)/(6.-1.25),0,1)
    weight=1-(3*z**2-2*z**3)
    theta=-np.deg2rad(a.angle)*weight
    assert np.all(weight[walls]==1) and np.all(weight[far]==0)
    q=xy.copy()
    moving=weight>0
    c,s=np.cos(theta[moving]),np.sin(theta[moving])
    q[moving]=np.column_stack((p[moving,0]*c-p[moving,1]*s+.5,p[moving,0]*s+p[moving,1]*c))
    q[moving]=np.array([[float(f"{x:.14g}"),float(f"{y:.14g}")] for x,y in q[moving]])
    assert np.array_equal(xy[far],q[far])
    edges=np.vstack(list(sections.values()))
    inside=edges[:,3]>0
    n=len(q)
    node_pair=edges[:,:2]-1
    keys=np.unique(np.r_[(edges[:,2,None]-1)*n+node_pair, (edges[inside,3,None]-1)*n+node_pair[inside]])
    ids,counts=np.unique(keys//n,return_counts=True)
    assert np.array_equal(ids,np.arange(109600)) and (counts==3).all()
    tri=(keys%n).reshape(-1,3)
    def twice_area(nodes):
        b,c=nodes[tri[:,1]]-nodes[tri[:,0]],nodes[tri[:,2]]-nodes[tri[:,0]]
        return b[:,0]*c[:,1]-b[:,1]*c[:,0]
    old,new=twice_area(xy),twice_area(q)
    assert np.isfinite(new).all() and np.min(new/old)>.5
    edge=sections[4][:,:2]-1
    length=lambda q:np.linalg.norm(q[edge[:,1]]-q[edge[:,0]],axis=1)
    assert np.max(abs(length(xy)-length(q)))<1e-12
    for i in np.where(moving)[0]:
        lines[start+i]=" ".join(f"{v:.14g}" for v in q[i])
    a.out.mkdir(parents=True)
    mesh=a.out/"airfoil.msh"
    mesh.write_text("\n".join(lines)+"\n",encoding="ascii",newline="\n")
    setup.update(angle_degrees=a.angle, geometry="scaled closed NACA0012 polygon rotated about (0.5,0) m",
        base_mesh_sha256=src_hash, mesh_sha256=hashlib.sha256(mesh.read_bytes()).hexdigest(),
        minimum_triangle_area_m2=float(abs(new).min()/2), min_area_ratio=float(np.min(new/old)),
        max_area_ratio=float(np.max(new/old)), far_boundary_unchanged=True, connectivity_unchanged=True,
        morph_radii_m=[1.25,6.], rotated_wall_nodes=len(walls))
    assert raw==a.source.read_bytes()
    (a.out/"setup.json").write_text(json.dumps(setup,indent=2)+"\n",encoding="ascii")
    print(json.dumps({k:setup[k] for k in ("angle_degrees","cells","min_area_ratio","max_area_ratio","mesh_sha256")}))

if __name__=="__main__":
    main()

