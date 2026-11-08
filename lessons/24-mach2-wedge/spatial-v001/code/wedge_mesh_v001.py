"""Original mesh emitter with explicit counts; no Fluent import or launch."""
from pathlib import Path
from math import tan
import numpy as np
from wedge_common_fine_v001 import NX as DEFAULT_NX,NY as DEFAULT_NY,CORNER,H,L,THETA
def make_mesh(out,nx=DEFAULT_NX,ny=DEFAULT_NY):
    NX,NY=int(nx),int(ny)
    if NX!=nx or NY!=ny or NX<5 or NY<2 or abs(NX*CORNER/L-round(NX*CORNER/L))>1e-12:raise ValueError("Positive grid counts and an exact corner grid node required")
    if any(out.iterdir()):raise FileExistsError("Use a fresh empty mesh directory")
    x = np.linspace(0, L, NX + 1)
    lower = np.maximum(x - CORNER, 0.) * tan(THETA)
    eta = np.linspace(0, 1, NY + 1)[:, None]
    yy = lower + eta * (H - lower)
    nodes = np.column_stack([np.broadcast_to(x, yy.shape).ravel(), yy.ravel()])
    ix = np.arange(NX)[None, :]
    iy = np.arange(NY)[:, None]
    k = (iy * (NX + 1) + ix).ravel()
    quads = np.column_stack([k, k + 1, k + NX + 2, k + NX + 1])
    polys = nodes[quads]
    nxt = np.roll(polys, -1, axis=1)
    cross = polys[..., 0] * nxt[..., 1] - nxt[..., 0] * polys[..., 1]
    areas = cross.sum(axis=1) / 2
    centres = ((polys + nxt) * cross[..., None]).sum(axis=1) / (6 * areas[:, None])
    assert np.all(areas > 0.)
    faces = {}
    for cid, q in enumerate(quads, 1):
        for a, b in zip(q, np.roll(q, -1)):
            key = tuple(sorted((int(a), int(b))))
            if key in faces:
                assert faces[key][3] == 0
                faces[key][3] = cid
            else:
                faces[key] = [int(a) + 1, int(b) + 1, cid, 0]
    zones = {zone: [] for zone in (3, 4, 5, 6, 7)}
    for a, b, c0, c1 in faces.values():
        n0, n1 = nodes[[a - 1, b - 1]]
        if c1:
            zone = 3
        elif n0[0] == n1[0] == 0:
            zone = 6
        elif n0[0] == n1[0] == L:
            zone = 7
        elif np.allclose([n0[1], n1[1]], H, atol=1e-13, rtol=0):
            zone = 4
        else:
            zone = 5
        zones[zone].append((a, b, c0, c1))
    assert [len(zones[z]) for z in (4, 5, 6, 7)] == [NX, NX, NY, NY]
    mesh = out / "wedge.msh"
    with mesh.open("x", encoding="ascii", newline="\n") as f:
        f.write('(0 "Mach 2 compression corner")\n(2 2)\n')
        f.write(f"(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n")
        for xx, y in nodes:
            f.write(f"{xx:.17g} {y:.17g}\n")
        f.write(f"))\n(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n")
        f.write(f"(13 (0 1 {len(faces):x} 0))\n")
        start = 1
        for zone, bc in [(3, 2), (4, 3), (5, 3), (6, 4), (7, 5)]:
            rows = zones[zone]
            f.write(f"(13 ({zone:x} {start:x} {start+len(rows)-1:x} {bc:x} 2)(\n")
            for row in rows:
                f.write(" ".join(f"{v:x}" for v in row) + "\n")
            f.write("))\n")
            start += len(rows)
        for zone, kind, name in ((2, "fluid", "fluid"), (3, "interior", "interior"),
                                  (4, "wall", "top"), (5, "wall", "ramp"),
                                  (6, "pressure-inlet", "inlet"), (7, "pressure-outlet", "outlet")):
            f.write(f"(45 ({zone} {kind} {name})())\n")
    np.savez_compressed(out / "geometry.npz", nodes=nodes, quads=quads, areas=areas,
                        centres=centres, x_nodes=x, y_nodes=yy)
    return mesh

def compare_refinement(coarse_path,fine_path):
    with np.load(coarse_path) as c,np.load(fine_path) as f:
        nx,ny=len(c['x_nodes'])-1,c['y_nodes'].shape[0]-1
        assert len(f['x_nodes'])-1==2*nx and f['y_nodes'].shape[0]-1==2*ny
        old=c['nodes'].reshape(ny+1,nx+1,2);new=f['nodes'].reshape(2*ny+1,2*nx+1,2)
        maximum=float(np.linalg.norm(old-new[::2,::2],axis=-1).max())
        actual=[]
        errors=[]
        for g in (c,f):
            p=g['nodes'][g['quads']];q=np.roll(p,-1,axis=1);a=.5*np.sum(p[:,:,0]*q[:,:,1]-q[:,:,0]*p[:,:,1],axis=1)
            assert np.all(a>0) and np.allclose(a,g['areas'],rtol=0,atol=2e-14)
            x=g['x_nodes'];y=g['y_nodes'];errors.extend([float(np.max(abs(y[0]-np.maximum(x-CORNER,0)*tan(THETA)))),float(np.max(abs(y[-1]-H)))])
            assert abs(x[0])<2e-14 and abs(x[-1]-L)<2e-14 and np.min(abs(x-CORNER))<2e-14
            actual.append(a)
        partition=actual[1].reshape(2*ny,2*nx).reshape(ny,2,nx,2).sum(axis=(1,3))-actual[0].reshape(ny,nx)
        err=float(np.max(abs(partition)));area=L*H-.5*(L-CORNER)**2*tan(THETA)
        assert maximum<=2e-14 and max(errors)<=2e-14 and err<=2e-14
        assert all(abs(a.sum()-area)<2e-14 for a in actual)
        return dict(coarse_cells=nx*ny,fine_cells=4*nx*ny,old_nodes_preserved=True,maximum_old_node_error_m=maximum,same_boundary_polygon=True,maximum_boundary_error_m=max(errors),maximum_parent_area_partition_error_m2=err,domain_area_m2=float(area),minimum_fine_area_m2=float(actual[1].min()),geometry_tolerance_m=2e-14,physical_geometry_changed=False)
