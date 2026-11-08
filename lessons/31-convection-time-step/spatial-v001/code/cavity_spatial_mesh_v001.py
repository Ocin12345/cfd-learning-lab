from pathlib import Path
import numpy as np,hashlib,json
from cavity_spatial_common_v001 import L,ROOT,specification,write

def generate(out, n):
    N = int(n)
    if N != n or N < 2: raise ValueError("Integer square count >=2 required")
    assert not out.exists()
    out.mkdir(parents=True)
    x = np.linspace(0, L, N + 1)
    nodes = np.array([(xx, yy) for yy in x for xx in x])

    def nid(i, j):
        return j * (N + 1) + i
    cells = np.array([[nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)] for j in range(N) for i in range(N)])
    centres = nodes[cells].mean(axis=1)
    areas = np.full(N * N, (L / N) ** 2)
    faces = {}
    for c, q in enumerate(cells, 1):
        for a, b in zip(q, np.roll(q, -1)):
            key = tuple(sorted((int(a), int(b))))
            if key in faces:
                assert faces[key][3] == 0
                faces[key][3] = c
            else:
                faces[key] = [int(a) + 1, int(b) + 1, c, 0]
    groups = {z: [] for z in [3, 4, 5, 6, 7]}
    for a, b, c0, c1 in faces.values():
        v0, v1 = (nodes[a - 1], nodes[b - 1])
        mid = (v0 + v1) / 2
        edge = v1 - v0
        normal = np.array([edge[1], -edge[0]])
        if c1:
            z = 3
            delta = centres[c1 - 1] - centres[c0 - 1]
        else:
            delta = mid - centres[c0 - 1]
            if v0[1] == v1[1] == L:
                z = 4
            elif v0[1] == v1[1] == 0:
                z = 5
            elif v0[0] == v1[0] == 0:
                z = 6
            else:
                assert v0[0] == v1[0] == L
                z = 7
        assert np.dot(normal, delta) > 0
        groups[z].append((a, b, c0, c1))
    assert len(nodes) - len(faces) + len(cells) == 1
    with (out / 'cavity.msh').open('x', encoding='ascii', newline='\n') as f:
        f.write('(0 "Shared stationary ideal-gas cavity; no computed flow")\n(2 2)\n')
        f.write(f'(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n')
        for xx, yy in nodes:
            f.write(f'{xx:.17g} {yy:.17g}\n')
        f.write(f'))\n(12 (0 1 {len(cells):x} 0))\n(12 (2 1 {len(cells):x} 1 3))\n(13 (0 1 {len(faces):x} 0))\n')
        first = 1
        for z, rows in groups.items():
            f.write(f'(13 ({z:x} {first:x} {first + len(rows) - 1:x} {(2 if z == 3 else 3):x} 2)(\n')
            for row in rows:
                f.write(' '.join((f'{v:x}' for v in row)) + '\n')
            f.write('))\n')
            first += len(rows)
        for z, kind, name in [(2, 'fluid', 'fluid'), (3, 'interior', 'interior'), (4, 'wall', 'top'), (5, 'wall', 'bottom'), (6, 'wall', 'left'), (7, 'wall', 'right')]:
            f.write(f'(45 ({z} {kind} {name})())\n')
    np.savez_compressed(out / 'geometry.npz', nodes=nodes, quads=cells, centres=centres, areas=areas, grid_shape=[N, N], x_edges=x, y_edges=x, ix=np.tile(np.arange(N), N), iy=np.repeat(np.arange(N), N))
    write(out / 'setup.json', dict(**dict(specification(),n=N,cells=N*N), mesh_sha256=hashlib.sha256((out / 'cavity.msh').read_bytes()).hexdigest(), positive_cell_areas=True, minimum_orthogonal_quality_expected=1.0, total_area_m2=float(areas.sum()), boundary_faces={str(z): len(v) for z, v in groups.items()}))
    (out / 'executed-mesher.py').write_bytes(Path(__file__).read_bytes())
    print(out)

def refinement_proof(coarse,fine):
 from cavity_restriction_v001 import restrict_square
 c=np.asarray(coarse['grid_shape']);f=np.asarray(fine['grid_shape']);assert np.array_equal(f,2*c) and c[0]==c[1]
 n=int(c[0]);old=coarse['nodes'].reshape(n+1,n+1,2);new=fine['nodes'].reshape(2*n+1,2*n+1,2)
 error=float(abs(old-new[::2,::2]).max());area_error=float(abs(restrict_square(fine['areas'],n)*4-coarse['areas']).max())
 assert error<2e-14 and area_error<1e-17
 assert np.all(fine['areas']>0) and abs(fine['areas'].sum()-.05**2)<1e-16
 for axis in ('x_edges','y_edges'):assert abs(fine[axis][0])<1e-16 and abs(fine[axis][-1]-.05)<1e-16
 return dict(coarse_cells=n*n,fine_cells=4*n*n,old_node_max_error_m=error,maximum_parent_area_error_m2=area_error,total_area_m2=float(fine['areas'].sum()),same_square_boundary=True)
