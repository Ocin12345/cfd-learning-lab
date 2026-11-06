"""Extend a preserved triangular Fluent airfoil mesh by an annulus only.

The added annulus uses Gmsh by default or opt-in deterministic radial layers;
original node IDs, coordinates and triangle cell IDs are retained. This
prepares a mesh, never starts a CFD solver.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import shutil
import sys

import numpy as np


SOURCE_READER = Path(__file__).with_name('plot_openfoam_airfoil_fields.py')
if not SOURCE_READER.exists():
    SOURCE_READER = Path(__file__).with_name('executed-source-reader.py')
_spec = importlib.util.spec_from_file_location('preserved_airfoil_source_reader', SOURCE_READER)
_reader = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_reader)
read_source_mesh = _reader.read_source_mesh


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def polygon_area(xy):
    return float(np.sum(xy[:, 0] * np.roll(xy[:, 1], -1)
                        - xy[:, 1] * np.roll(xy[:, 0], -1)) / 2)


def ordered_loop(xy, edges):
    """Order a single closed degree-two boundary, returning CCW node IDs."""
    adjacency = {}
    unique = set()
    for a, b in np.asarray(edges, dtype=np.int64):
        a, b = int(a), int(b)
        key = tuple(sorted((a, b)))
        require(a != b and key not in unique, 'Repeated/degenerate boundary edge')
        unique.add(key)
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    require(len(adjacency) >= 3 and all(len(v) == 2 for v in adjacency.values()),
            'Boundary must be one closed degree-two loop')
    start = min(adjacency)
    nodes, previous, current = [start], None, start
    while True:
        candidate = next(v for v in adjacency[current] if v != previous)
        if candidate == start:
            break
        require(candidate not in nodes, 'Boundary self-cycle')
        nodes.append(candidate)
        previous, current = current, candidate
    require(len(nodes) == len(adjacency) == len(unique), 'Disconnected boundary loops')
    if polygon_area(xy[nodes]) < 0:
        nodes = [nodes[0], *nodes[:0:-1]]
    require(polygon_area(xy[nodes]) > 0, 'Zero boundary polygon area')
    return np.asarray(nodes, dtype=np.int64)


def signed_triangle_areas(xy, triangles):
    a, b, c = [xy[triangles[:, i]] for i in range(3)]
    return ((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1])
            - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])) / 2


def validate_radial_polygon(xy, ids, radius):
    """A radius-declared polygon may contain subdivisions along its chords."""
    poly = xy[ids]
    radii = np.linalg.norm(poly - [.5, 0], axis=1)
    vertices = np.flatnonzero(abs(radii - radius) <= 1e-10)
    require(len(vertices) >= 3 and np.all(radii <= radius + 1e-10),
            'Polygon has wrong declared radius')
    for j, first in enumerate(vertices):
        last = int(vertices[(j + 1) % len(vertices)])
        stop = last if last > first else last + len(poly)
        interior = poly[np.arange(first, stop + 1) % len(poly)]
        a, b = poly[first], poly[last]
        v = b - a
        length = float(np.linalg.norm(v))
        cross = v[0] * (interior[:, 1] - a[1]) - v[1] * (interior[:, 0] - a[0])
        t = (interior - a) @ v / (length * length)
        require(np.max(abs(cross)) <= 1e-10 * length
                and np.all(t >= -1e-10) and np.all(t <= 1 + 1e-10),
                'Boundary subdivision is not on original polygon chord')
    return dict(radius_vertices=len(vertices), chord_subdivision_nodes=len(poly)-len(vertices),
                node_radius_min_m=float(radii.min()), node_radius_max_m=float(radii.max()))


def orient_triangles(xy, triangles):
    triangles = np.asarray(triangles, dtype=np.int64).copy()
    require(triangles.ndim == 2 and triangles.shape[1] == 3
            and np.all(triangles >= 0) and np.all(triangles < len(xy)), 'Invalid triangles')
    areas = signed_triangle_areas(xy, triangles)
    require(np.isfinite(areas).all() and np.all(areas != 0), 'Zero/nonfinite triangle area')
    reverse = areas < 0
    triangles[reverse] = triangles[reverse][:, [0, 2, 1]]
    return triangles


def edge_codes(edges, nodes):
    edges = np.sort(np.asarray(edges, dtype=np.int64), axis=1)
    return edges[:, 0] * nodes + edges[:, 1]


def build_faces(xy, triangles, wall_edges, outer_edges):
    """Fluent c0 is left of the oriented edge, c1 right (CCW cell edges)."""
    require(np.all(signed_triangle_areas(xy, triangles) > 0), 'Cells must be CCW')
    edges = triangles[:, [[0, 1], [1, 2], [2, 0]]].reshape(-1, 2)
    cells = np.repeat(np.arange(len(triangles), dtype=np.int64), 3)
    codes = edge_codes(edges, len(xy))
    order = np.argsort(codes, kind='stable')
    sorted_codes = codes[order]
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_codes)) + 1]
    counts = np.diff(np.r_[starts, len(order)])
    require(np.all((counts == 1) | (counts == 2)), 'Nonmanifold edge')
    first = order[starts]
    paired = counts == 2
    second = order[starts[paired] + 1]
    require(np.array_equal(edges[first[paired]], edges[second, ::-1]),
            'Adjacent cell edges have inconsistent orientation/overlapping cells')
    faces = np.column_stack([edges[first], cells[first], np.full(len(first), -1)])
    faces[paired, 3] = cells[second]
    single = ~paired
    wall_codes, outer_codes = [edge_codes(e, len(xy)) for e in [wall_edges, outer_edges]]
    require(len(np.unique(wall_codes)) == len(wall_codes)
            and len(np.unique(outer_codes)) == len(outer_codes)
            and not np.isin(wall_codes, outer_codes).any(), 'Invalid boundary edge declarations')
    boundary_codes = sorted_codes[starts[single]]
    is_wall = np.isin(boundary_codes, wall_codes)
    is_outer = np.isin(boundary_codes, outer_codes)
    require(np.all(is_wall ^ is_outer)
            and int(is_wall.sum()) == len(wall_codes)
            and int(is_outer.sum()) == len(outer_codes), 'Missing/unexpected boundary edge')
    boundary = faces[single]
    zone = np.where(is_wall, 4, np.where(xy[boundary[:, :2], 0].mean(axis=1) < .5, 5, 6))
    return {3: faces[paired], **{z: boundary[zone == z] for z in [4, 5, 6]}}


def validate_merged(source, xy, triangles, groups, interface_edges, outer_edges, radius):
    """Assert ID-exact preservation, conformal connectivity and polygon areas."""
    n, c = len(source['xy']), source['cells']
    require(np.array_equal(xy[:n], source['xy']), 'Original node coordinates/IDs changed')
    require(np.array_equal(np.sort(triangles[:c], axis=1), source['triangles']),
            'Original cell IDs/triangles changed')
    require(len(triangles) > c and len(xy) > n, 'No annulus added')
    source_wall = source['edges'][source['zones'] == 4]
    require(np.array_equal(np.sort(edge_codes(groups[4][:, :2], len(xy))),
                           np.sort(edge_codes(source_wall, len(xy)))), 'Original wall edges changed')
    interior_codes = edge_codes(groups[3][:, :2], len(xy))
    wanted = edge_codes(interface_edges, len(xy))
    interface = groups[3][np.isin(interior_codes, wanted)]
    require(len(interface) == len(wanted)
            and np.all((interface[:, 2] < c) ^ (interface[:, 3] < c)),
            'Former farfield is not a conformal inner/annulus interface')
    outer_ids = ordered_loop(xy, outer_edges)
    inner_ids = ordered_loop(xy, interface_edges)
    wall_ids = ordered_loop(xy, source_wall)
    outer_geometry = validate_radial_polygon(xy, outer_ids, radius)
    source_radius = float(np.linalg.norm(xy[inner_ids] - [.5, 0], axis=1).max())
    inner_geometry = validate_radial_polygon(xy, inner_ids, source_radius)
    scaled = (xy[inner_ids] - [.5, 0]) * (radius / source_radius) + [.5, 0]
    actual = xy[outer_ids]
    distances = np.linalg.norm(actual[:, None, :] - scaled[None, :, :], axis=2)
    nearest = distances.argmin(axis=1)
    require(len(actual) == len(scaled) and len(np.unique(nearest)) == len(scaled)
            and np.max(distances[np.arange(len(actual)), nearest]) <= 1e-10,
            'Outer polygon is not uniform scaling of original polygon/subdivisions')
    expected = polygon_area(xy[outer_ids]) - polygon_area(xy[wall_ids])
    inner_area = polygon_area(xy[inner_ids]) - polygon_area(xy[wall_ids])
    annulus_area = polygon_area(xy[outer_ids]) - polygon_area(xy[inner_ids])
    areas = signed_triangle_areas(xy, triangles)
    require(np.all(areas > 0) and np.isfinite(areas).all(), 'Invalid merged cell areas')
    require(math.isclose(math.fsum(areas[:c]), inner_area, rel_tol=1e-11, abs_tol=1e-9)
            and math.isclose(math.fsum(areas[c:]), annulus_area, rel_tol=1e-11, abs_tol=1e-9)
            and math.isclose(math.fsum(areas), expected, rel_tol=1e-11, abs_tol=1e-9),
            'Triangle and boundary polygon areas differ')
    # The original farfield is a convex polygon: every added centroid must be outside it.
    centres = xy[triangles[c:]].mean(axis=1)
    outside = np.zeros(len(centres), dtype=bool)
    poly = xy[inner_ids]
    for a, b in zip(poly, np.roll(poly, -1, axis=0)):
        v = b - a
        outside |= v[0] * (centres[:, 1] - a[1]) - v[1] * (centres[:, 0] - a[0]) < -1e-12
    require(outside.all(), 'Added annulus cell centroid inside original farfield')
    # A single connected cell graph is required; no unattached annular islands.
    parent = np.arange(len(triangles))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for a, b in groups[3][:, 2:]:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    require(all(find(i) == find(0) for i in range(len(triangles))), 'Disconnected cell graph')
    def perimeter(ids):
        return float(np.linalg.norm(xy[ids] - np.roll(xy[ids], -1, axis=0), axis=1).sum())
    return dict(original_node_coordinates_exact=True, original_node_ids_preserved=True,
                original_cell_ids_and_triangles_exact=True, original_wall_edge_ids_exact=True,
                former_boundary_edges_now_interior=len(interface), manifold_oriented_faces=True,
                connected_cell_graph=True, all_added_centroids_outside_source_polygon=True,
                original_nodes=n, original_cells=c, added_nodes=len(xy)-n, added_cells=len(triangles)-c,
                total_area_m2=math.fsum(areas), original_fluid_area_m2=math.fsum(areas[:c]),
                annulus_area_m2=math.fsum(areas[c:]), polygon_expected_area_m2=expected,
                minimum_triangle_area_m2=float(areas.min()),
                wall_perimeter_m=perimeter(wall_ids), original_farfield_perimeter_m=perimeter(inner_ids),
                outer_perimeter_m=perimeter(outer_ids), original_farfield_geometry=inner_geometry,
                outer_farfield_geometry=outer_geometry, outer_polygon_uniform_radial_scale=True,
                boundary_faces={str(z):len(groups[z]) for z in [4,5,6]})


def write_fluent(path, xy, triangles, groups, original_node_lines=None):
    with path.open('x', encoding='ascii', newline='\n') as stream:
        stream.write('(0 "Preserved inner airfoil triangles; outer annulus extension only")\n(2 2)\n')
        stream.write(f'(10 (0 1 {len(xy):x} 1 2))\n(10 (1 1 {len(xy):x} 1 2)(\n')
        for i,(x,y) in enumerate(xy):
            stream.write((original_node_lines[i] if original_node_lines is not None and i<len(original_node_lines)
                          else f'{x:.17g} {y:.17g}')+'\n')
        stream.write('))\n')
        stream.write(f'(12 (0 1 {len(triangles):x} 0))\n(12 (2 1 {len(triangles):x} 1 1))\n')
        stream.write(f'(13 (0 1 {sum(map(len,groups.values())):x} 0))\n')
        first=1
        for zone,bc in [(3,2),(4,3),(5,10),(6,5)]:
            rows=groups[zone]
            require(len(rows)>0,'Empty Fluent face zone')
            stream.write(f'(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n')
            for row in rows:
                stream.write(' '.join(f'{int(v)+1:x}' for v in row)+'\n')
            stream.write('))\n');first+=len(rows)
        for zone,kind,name in [(2,'fluid','fluid'),(3,'interior','interior'),(4,'wall','airfoil'),
                               (5,'velocity-inlet','inlet'),(6,'pressure-outlet','outlet')]:
            stream.write(f'(45 ({zone} {kind} {name})())\n')


def generate_radial(source_xy, inner_ids, radius, layers, source_radius):
    """Append scaled copies of the actual inner polygon and split each strip."""
    require(isinstance(layers, int) and layers >= 1, 'At least one radial layer required')
    require(math.isfinite(radius) and radius > source_radius > 0, 'Invalid radial radii')
    inner = source_xy[inner_ids]
    factors = (radius / source_radius) ** (np.arange(1, layers + 1) / layers)
    rings = (inner[None, :, :] - [.5, 0]) * factors[:, None, None] + [.5, 0]
    xy = np.vstack([source_xy, rings.reshape(-1, 2)])
    count = len(inner_ids)
    previous = np.asarray(inner_ids)
    blocks = []
    for j in range(layers):
        current = len(source_xy) + j * count + np.arange(count)
        next_previous, next_current = np.roll(previous, -1), np.roll(current, -1)
        # Quad order is inner(i), outer(i), outer(i+1), inner(i+1).
        # Alternate the diagonal by both angular segment and radial strip.
        even = (np.arange(count) + j) % 2 == 0
        a = np.column_stack([previous, current, next_current])
        b = np.column_stack([previous, next_current, next_previous])
        a[~even] = np.column_stack([previous, current, next_previous])[~even]
        b[~even] = np.column_stack([current, next_current, next_previous])[~even]
        blocks.append(np.stack([a, b], axis=1).reshape(-1, 3))
        previous = current
    triangles = orient_triangles(xy, np.vstack(blocks))
    info = dict(method='deterministic_radial_triangles', gmsh_used=False,
                radial_layers=layers, new_ring_vertices_per_layer=count,
                radius_schedule='R_j = R0 * (R/R0)**(j/layers)',
                radial_growth_factor=(radius/source_radius)**(1/layers),
                diagonal='Alternate by parity of angular segment plus radial strip',
                inner_segments=count, outer_segments=count)
    return xy, triangles, previous, info


def centroid_cosine_diagnostic(xy, triangles, groups, source_cells):
    """Geometric cosine diagnostic, explicitly not Fluent orthogonal quality."""
    centres = xy[triangles].mean(axis=1)
    all_faces = np.vstack(list(groups.values()))
    edges, owner, neighbour = all_faces[:, :2], all_faces[:, 2], all_faces[:, 3]
    v = xy[edges[:, 1]] - xy[edges[:, 0]]
    normal = np.column_stack([v[:, 1], -v[:, 0]])
    face_centres = xy[edges].mean(axis=1)
    def cosine(delta, n):
        return abs(np.sum(delta*n, axis=1)) / (np.linalg.norm(delta, axis=1)*np.linalg.norm(n, axis=1))
    quality = cosine(face_centres-centres[owner], normal)
    interior = neighbour >= 0
    quality[interior] = np.minimum(quality[interior], cosine(
        centres[neighbour[interior]]-centres[owner[interior]], normal[interior]))
    # Include the neighbour-cell-to-face vector too, independent of c0/c1 ordering.
    quality[interior] = np.minimum(quality[interior], cosine(
        face_centres[interior]-centres[neighbour[interior]], normal[interior]))
    require(np.isfinite(quality).all() and np.all(quality > 0), 'Invalid centroid cosine diagnostic')
    interface = interior & ((owner < source_cells) ^ (neighbour < source_cells))
    annulus = (owner >= source_cells) | (interior & (neighbour >= source_cells))
    return dict(definition='Minimum absolute cosine of oriented face normal with each incident cell-centroid-to-face-centre vector and, for interior faces, the cell-centroid connection',
                scope='Mesh-only geometric diagnostic; not the Fluent orthogonal-quality metric or solver admission',
                minimum_all_faces_cosine=float(quality.min()),
                minimum_former_interface_cosine=float(quality[interface].min()),
                minimum_faces_incident_to_annulus_cosine=float(quality[annulus].min()))


def generate_annulus(source_xy, inner_ids, radius, output, source_radius):
    import gmsh
    gmsh.initialize()
    try:
        options={'General.NumThreads':2,'Mesh.MaxNumThreads1D':2,'Mesh.MaxNumThreads2D':2,
                 'Mesh.Algorithm':6,'Mesh.ElementOrder':1,'Mesh.MeshSizeFromPoints':0,
                 'Mesh.MeshSizeFromCurvature':0,'Mesh.MeshSizeExtendFromBoundary':0,
                 'Mesh.MeshSizeMin':.02,'Mesh.MeshSizeMax':3.,'Mesh.RandomSeed':1,
                 'Mesh.MshFileVersion':4.1}
        for name,value in options.items():gmsh.option.setNumber(name,value)
        gmsh.model.add('outer-annulus-preserved-inner-interface')
        inner=source_xy[inner_ids]
        r=np.linalg.norm(inner-[.5,0],axis=1)
        outer=(inner-[.5,0])*(radius/source_radius)+[.5,0]
        inner_points=[gmsh.model.geo.addPoint(float(x),float(y),0) for x,y in inner]
        outer_points=[gmsh.model.geo.addPoint(float(x),float(y),0) for x,y in outer]
        def curves(points):
            return [gmsh.model.geo.addLine(a,b) for a,b in zip(points,points[1:]+points[:1])]
        inner_curves,outer_curves=curves(inner_points),curves(outer_points)
        il=gmsh.model.geo.addCurveLoop([-tag for tag in inner_curves[::-1]])
        ol=gmsh.model.geo.addCurveLoop(outer_curves)
        surface=gmsh.model.geo.addPlaneSurface([ol,il]);gmsh.model.geo.synchronize()
        for tag in inner_curves+outer_curves:gmsh.model.mesh.setTransfiniteCurve(tag,2)
        gmsh.model.addPhysicalGroup(1,inner_curves,7,'former-farfield-interface')
        for zone,name,test in [(5,'inlet',True),(6,'outlet',False)]:
            lines=[tag for i,tag in enumerate(outer_curves)
                   if bool((outer[i,0]+outer[(i+1)%len(outer),0])/2<.5)==test]
            gmsh.model.addPhysicalGroup(1,lines,zone,name)
        gmsh.model.addPhysicalGroup(2,[surface],2,'annulus')
        # Only added-region background sizes vary; source triangles are never loaded into Gmsh.
        field=gmsh.model.mesh.field.add('MathEval')
        r0=source_radius
        formula=f'1.5 + 1.5*(Sqrt((x-0.5)^2+y^2)-{r0:.17g})/{radius-r0:.17g}'
        gmsh.model.mesh.field.setString(field,'F',formula)
        gmsh.model.mesh.field.setAsBackgroundMesh(field)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(output/'source-gmsh-annulus.msh'))
        tags,coords,_=gmsh.model.mesh.getNodes();coords=np.asarray(coords).reshape(-1,3)
        lookup={int(t):i for i,t in enumerate(tags)}
        mapping={}
        for node,entity in zip(inner_ids,inner_points):
            values,_,_=gmsh.model.mesh.getNodes(0,entity,True)
            require(len(values)==1,'Shared Gmsh geometry vertex has wrong node count')
            tag=int(values[0]);require(np.max(abs(coords[lookup[tag],:2]-source_xy[node]))<=1e-12,
                                       'Gmsh moved shared interface vertex')
            mapping[tag]=int(node)
        merged=[tuple(v) for v in source_xy]
        for tag,xyz in zip(tags,coords):
            if int(tag) not in mapping:
                require(abs(xyz[2])<1e-12,'Nonplanar annulus node')
                mapping[int(tag)]=len(merged);merged.append(tuple(xyz[:2]))
        outer_ids=[]
        for entity in outer_points:
            values,_,_=gmsh.model.mesh.getNodes(0,entity,True)
            require(len(values)==1,'Outer geometry vertex has wrong node count')
            outer_ids.append(mapping[int(values[0])])
        triangles=[];element_tags=[]
        types,tag_blocks,blocks=gmsh.model.mesh.getElements(2,surface)
        for kind,etags,block in zip(types,tag_blocks,blocks):
            require(kind==2,'Nonlinear/nontriangular annulus element')
            triangles.extend([[mapping[int(t)] for t in row] for row in np.asarray(block).reshape(-1,3)])
            element_tags.extend(map(int,etags))
        for curve in inner_curves+outer_curves:
            types,_,blocks=gmsh.model.mesh.getElements(1,curve)
            require(list(types)==[1] and len(blocks[0])==2,'Shared/outer segment was subdivided')
        quality=gmsh.model.mesh.getElementQualities(element_tags,'minSICN')
        require(min(quality)>0,'Nonpositive annulus Gmsh minSICN quality')
        info=dict(method='gmsh_unstructured_annulus',gmsh_used=True,
                  gmsh_version=gmsh.__version__,options=options,background_size_formula=formula,
                  annulus_minimum_gmsh_minSICN=float(min(quality)),
                  inner_segments=len(inner_curves),outer_segments=len(outer_curves),
                  transfinite_nodes_per_boundary_segment=2)
        return np.asarray(merged),np.asarray(triangles,dtype=np.int64),np.asarray(outer_ids),info
    finally:
        gmsh.finalize()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--radius',type=float,required=True,help='New outer radius in metres about (0.5,0)')
    parser.add_argument('--annulus-mode',choices=['gmsh','radial'],default='gmsh',
                        help='Gmsh annulus or deterministic scaled polygon layers')
    parser.add_argument('--radial-layers',type=int,default=96,help='Radial strips for radial mode only')
    args=parser.parse_args()
    source,out=args.source.resolve(),args.out.resolve()
    project=Path('E:/CFD-Learning-Lab').resolve();repo=Path(__file__).resolve().parents[1]
    require(source.is_relative_to(project) and out.is_relative_to(project), 'Inputs/outputs must stay under project root')
    require(not out.is_relative_to(repo) and not out.exists(),'Use a fresh external native output directory')
    setup_path=source.parent/'setup.json';setup=json.loads(setup_path.read_text())
    old_radius=float(setup['outer_radius_m'])
    require(math.isfinite(args.radius) and args.radius>old_radius,'Radius must exceed source radius')
    require(args.radial_layers>=1,'At least one radial layer required')
    require(setup['mesh_type']=='triangles','Source must be triangular')
    src=read_source_mesh(source)
    require(src['cells']==setup['cells'] and len(src['xy'])==setup['nodes'],'Source setup counts differ')
    interface_edges=src['edges'][np.isin(src['zones'],[5,6])]
    inner_ids=ordered_loop(src['xy'],interface_edges)
    if args.annulus_mode=='radial':
        require(src['cells']+2*args.radial_layers*len(inner_ids)<950000
                and len(src['xy'])+args.radial_layers*len(inner_ids)<950000,
                'Radial mesh exceeds 950k cell/node admission ceiling')
    validate_radial_polygon(src['xy'],inner_ids,old_radius)
    source_hash,setup_hash,reader_hash=sha256(source),sha256(setup_path),sha256(SOURCE_READER)
    raw=source.read_text().splitlines();start=next(i for i,s in enumerate(raw) if s.startswith('(10 (1 '))+1
    original_lines=raw[start:start+len(src['xy'])]
    out.mkdir(parents=True)
    shutil.copy2(Path(__file__),out/'executed-mesher.py')
    shutil.copy2(SOURCE_READER,out/'executed-source-reader.py')
    if args.annulus_mode=='radial':
        xy,added,outer_ids,mesher_info=generate_radial(src['xy'],inner_ids,args.radius,args.radial_layers,old_radius)
        np.savez_compressed(out/'radial-annulus-topology.npz',added_nodes_xy=xy[len(src['xy']):],
                            added_triangles=added,
                            radius_factors=(args.radius/old_radius)**(np.arange(1,args.radial_layers+1)/args.radial_layers))
        artifact=dict(radial_topology_file='radial-annulus-topology.npz',
                      radial_topology_sha256=sha256(out/'radial-annulus-topology.npz'))
    else:
        xy,added,outer_ids,mesher_info=generate_annulus(src['xy'],inner_ids,args.radius,out,old_radius)
        artifact=dict(annulus_mesh_sha256=sha256(out/'source-gmsh-annulus.msh'))
    triangles=orient_triangles(xy,np.vstack([src['triangles'],added]))
    cell_cap=950000 if args.annulus_mode=='radial' else 450000
    require(len(triangles)<cell_cap and len(xy)<950000,'Mesh exceeds mode-specific cell/node admission limit')
    outer_edges=np.column_stack([outer_ids,np.roll(outer_ids,-1)])
    groups=build_faces(xy,triangles,src['edges'][src['zones']==4],outer_edges)
    proof=validate_merged(src,xy,triangles,groups,interface_edges,outer_edges,args.radius)
    mesh_diagnostic=centroid_cosine_diagnostic(xy,triangles,groups,src['cells'])
    write_fluent(out/'airfoil.msh',xy,triangles,groups,original_lines)
    reread=read_source_mesh(out/'airfoil.msh')
    require(np.array_equal(reread['xy'],xy) and np.array_equal(reread['triangles'],np.sort(triangles,axis=1)),
            'Written Fluent node/cell readback differs')
    proof['written_fluent_readback_exact']=True
    require(sha256(source)==source_hash and sha256(setup_path)==setup_hash and sha256(SOURCE_READER)==reader_hash,
            'Source inputs changed during preparation')
    provenance=dict(operation='annulus-only domain extension; no CFD',source=str(source),source_sha256=source_hash,
                    source_setup_sha256=setup_hash,source_reader_sha256=reader_hash,script_sha256=sha256(Path(__file__)),
                    old_radius_m=old_radius,new_radius_m=args.radius,centre_m=[.5,0],
                    nodes=len(xy),cells=len(triangles),preservation_checks=proof,mesher=mesher_info,
                    mesh_sha256=sha256(out/'airfoil.msh'),**artifact,
                    centroid_cosine_diagnostic=mesh_diagnostic,
                    boundary_assignment='Original wall unchanged; former inlet/outlet edges interior; new outer xmid<0.5 inlet, other outlet',
                    scientific_scope='Exact inner mesh held fixed; outer domain and added annulus discretization change. This is not mesh independence or physical validation.')
    setup.update(outer_radius_m=args.radius,cells=len(triangles),nodes=len(xy),
                 source_mesh_sha256=source_hash,minimum_triangle_area_m2=proof['minimum_triangle_area_m2'],
                 domain_extension=provenance,preserved_source_node_coordinates=True,
                 preserved_source_triangle_cell_ids=True)
    for name,value in [('setup.json',setup),('domain-extension-provenance.json',provenance)]:
        with (out/name).open('x',encoding='utf-8') as stream:json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps(provenance,indent=2,allow_nan=False),flush=True)


if __name__=='__main__':
    main()
