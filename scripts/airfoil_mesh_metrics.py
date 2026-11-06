"""Audit a saved Gmsh triangular airfoil mesh without meshing, GUI or CFD.

Requires numpy and gmsh. Source geometry is the unrotated wall/outer polygon in
--source, with its original circumferential count in source.parent/setup.json.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

import gmsh
import numpy as np

TOLERANCE_M = 1e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quantiles(values):
    values = np.asarray(values)
    require(values.size > 0, 'Empty region: cannot report resolution statistics')
    return dict(count=int(values.size), minimum=float(values.min()),
                p10=float(np.quantile(values, .1)), median=float(np.median(values)),
                p90=float(np.quantile(values, .9)), maximum=float(values.max()))


def polygon_area(points):
    return abs(float(.5 * np.sum(points[:, 0] * np.roll(points[:, 1], -1)
                                - points[:, 1] * np.roll(points[:, 0], -1))))


def source_polygons(source):
    setup_path = source.parent / 'setup.json'
    setup = json.loads(setup_path.read_text())
    require(setup['chord_m'] == 1 and setup['angle_degrees'] == 0,
            'Fixed metric regions require an unrotated unit-chord source')
    nt = int(setup['circumferential_cells'])
    lines = source.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('(10 (1 ')) + 1
    end = next(i for i in range(start, len(lines)) if lines[i] == '))')
    points = np.array([list(map(float, line.split())) for line in lines[start:end]])
    require(points.ndim == 2 and points.shape[1] == 2 and len(points) >= 2 * nt,
            'Expected original 2D structured ring mesh')
    require(np.isfinite(points).all(), 'Nonfinite source coordinates')
    return points[:nt], points[-nt:], setup_path


def polygon_preservation(edges, reference, xy):
    """Check every new edge lies in exactly one original segment, without gaps/overlap."""
    vectors = np.roll(reference, -1, axis=0) - reference
    squared = np.sum(vectors * vectors, axis=1)
    require(np.all(squared > 0), 'Zero-length original polygon segment')
    lengths = np.sqrt(squared)
    coverage = [[] for _ in reference]
    maximum_distance = 0.
    for edge in edges:
        p, q = xy[edge]
        wp, wq = p - reference, q - reference
        dp = np.abs(vectors[:, 0] * wp[:, 1] - vectors[:, 1] * wp[:, 0]) / lengths
        dq = np.abs(vectors[:, 0] * wq[:, 1] - vectors[:, 1] * wq[:, 0]) / lengths
        tp = np.sum(wp * vectors, axis=1) / squared
        tq = np.sum(wq * vectors, axis=1) / squared
        param_tolerance = TOLERANCE_M / lengths
        valid = ((dp <= TOLERANCE_M) & (dq <= TOLERANCE_M)
                 & (tp >= -param_tolerance) & (tp <= 1 + param_tolerance)
                 & (tq >= -param_tolerance) & (tq <= 1 + param_tolerance))
        require(np.count_nonzero(valid) == 1,
                'Boundary edge is outside or crosses original polygon segments')
        j = int(np.flatnonzero(valid)[0])
        coverage[j].append(tuple(sorted((float(tp[j]), float(tq[j])))))
        maximum_distance = max(maximum_distance, float(dp[j]), float(dq[j]))
    maximum_coverage_error = 0.
    for j, intervals in enumerate(coverage):
        require(intervals, 'An original polygon segment is missing')
        previous = 0.
        for low, high in sorted(intervals):
            gap_or_overlap = abs(low - previous) * lengths[j]
            require(gap_or_overlap <= TOLERANCE_M, 'Boundary segment gap or overlap')
            maximum_coverage_error = max(maximum_coverage_error, gap_or_overlap)
            require(high > low, 'Zero-length boundary subdivision')
            previous = high
        error = abs(previous - 1) * lengths[j]
        require(error <= TOLERANCE_M, 'Original polygon segment is incompletely covered')
        maximum_coverage_error = max(maximum_coverage_error, error)
    boundary_points = xy[np.unique(edges)]
    vertex_distance = max(float(np.linalg.norm(boundary_points - p, axis=1).min())
                          for p in reference)
    require(vertex_distance <= TOLERANCE_M, 'Original polygon vertex missing')
    return dict(passed=True, original_segments=len(reference),
                maximum_original_vertex_distance_m=vertex_distance,
                maximum_edge_to_original_segment_distance_m=maximum_distance,
                maximum_segment_gap_overlap_or_endpoint_error_m=maximum_coverage_error)


def audit(mesh_path, source):
    wall_reference, outer_reference, source_setup = source_polygons(source)
    gmsh.initialize([], readConfigFiles=False)
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.open(str(mesh_path))
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        order = np.argsort(tags)
        ids = tags[order]
        require(len(ids) > 0 and len(np.unique(ids)) == len(ids), 'Invalid node IDs')
        xyz = np.asarray(xyz).reshape(-1, 3)[order]
        require(np.isfinite(xyz).all() and np.max(np.abs(xyz[:, 2])) <= TOLERANCE_M,
                'Expected finite planar z=0 mesh')
        xy = xyz[:, :2]

        def indices(node_tags):
            found = np.searchsorted(ids, node_tags)
            require(np.all(found < len(ids)), 'Unknown node reference')
            require(np.array_equal(ids[found], node_tags), 'Unknown node reference')
            return found

        surfaces = gmsh.model.getEntities(2)
        require(len(surfaces) == 1, 'Only a single fluid surface is supported')
        fluid_groups = gmsh.model.getPhysicalGroups(2)
        require(len(fluid_groups) == 1, 'Expected one fluid physical group')
        dimension, fluid_tag = fluid_groups[0]
        require(gmsh.model.getPhysicalName(dimension, fluid_tag) == 'fluid'
                and list(gmsh.model.getEntitiesForPhysicalGroup(dimension, fluid_tag))
                == [surfaces[0][1]], 'Fluid physical group does not cover the mesh surface')
        types, _, blocks = gmsh.model.mesh.getElements(2)
        require(list(types) == [2], 'Only first-order triangular meshes are supported')
        triangles = indices(blocks[0].reshape(-1, 3))
        count = len(triangles)
        vertices = xy[triangles]
        ab, ac = vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]
        areas = np.abs(ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]) / 2
        require(count > 0 and np.all(areas > 0), 'Degenerate or empty triangle mesh')
        groups = {}
        for dim, tag in gmsh.model.getPhysicalGroups(1):
            name = gmsh.model.getPhysicalName(dim, tag)
            require(name not in groups, 'Duplicate boundary group name')
            edges = []
            for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                kinds, _, values = gmsh.model.mesh.getElements(1, int(entity))
                require(list(kinds) == [1], 'Only straight first-order boundary edges supported')
                edges.append(indices(values[0].reshape(-1, 2)))
            require(edges, 'Empty physical boundary group')
            groups[name] = np.concatenate(edges)
        require(set(groups) == {'airfoil', 'inlet', 'outlet'}, 'Unexpected boundary groups')
        wall = groups['airfoil']
        outer = np.concatenate([groups['inlet'], groups['outlet']])
        nn = len(ids)

        def edge_keys(edges):
            return np.min(edges, axis=1) * nn + np.max(edges, axis=1)

        all_keys = np.concatenate([edge_keys(triangles[:, [i, j]])
                                   for i, j in [(0, 1), (1, 2), (2, 0)]])
        sorter = np.argsort(all_keys)
        sorted_keys = all_keys[sorter]
        unique, occurrences = np.unique(sorted_keys, return_counts=True)
        require(np.all(occurrences <= 2), 'Nonmanifold triangle edges')
        boundary_keys = edge_keys(np.concatenate(list(groups.values())))
        require(len(np.unique(boundary_keys)) == len(boundary_keys), 'Duplicate boundary edge')
        require(np.array_equal(np.sort(boundary_keys), unique[occurrences == 1]),
                'Physical boundaries do not exactly cover triangle boundary edges')
        wall_keys = edge_keys(wall)
        positions = np.searchsorted(sorted_keys, wall_keys)
        adjacent = sorter[positions] % count
        lengths = np.linalg.norm(xy[wall[:, 1]] - xy[wall[:, 0]], axis=1)
        heights = 2 * areas[adjacent] / lengths
        midpoint = xy[wall].mean(axis=1)
        centroids = vertices.mean(axis=1)
        wake = ((centroids[:, 0] >= 1) & (centroids[:, 0] <= 3)
                & (np.abs(centroids[:, 1]) < .25))
        split_x = .5
        for name, sign in [('inlet', -1), ('outlet', 1)]:
            edge_x = xy[groups[name]].mean(axis=1)[:, 0]
            require(np.all(sign * (edge_x - split_x) > 0), 'Incorrect inlet/outlet assignment')
        wall_check = polygon_preservation(wall, wall_reference, xy)
        outer_check = polygon_preservation(outer, outer_reference, xy)
        area = float(areas.sum())
        reference_area = polygon_area(outer_reference) - polygon_area(wall_reference)
        require(math.isclose(area, reference_area, rel_tol=1e-12, abs_tol=1e-10),
                'Triangle area differs from source polygon fluid area')
        regions = {}
        for name, mask in [('all', np.ones(len(wall), dtype=bool)),
                           ('leading_edge_x_le_0.05', midpoint[:, 0] <= .05),
                           ('trailing_edge_x_ge_0.95', midpoint[:, 0] >= .95)]:
            regions[name] = dict(edge_length_m=quantiles(lengths[mask]),
                                 adjacent_triangle_normal_height_m=quantiles(heights[mask]))
        perimeters = {name: float(np.linalg.norm(xy[edges[:, 1]] - xy[edges[:, 0]], axis=1).sum())
                      for name, edges in groups.items()}
        return dict(schema_version=1, checks_passed=True, cells=count, nodes=nn,
                    wall_face_count=len(wall), boundary_face_counts={k: len(v) for k, v in groups.items()},
                    wall_regions=regions, wake=dict(triangle_sqrt_area_m=quantiles(np.sqrt(areas[wake])),
                                                   triangle_area_m2=quantiles(areas[wake])),
                    triangle_area_m2=quantiles(areas), fluid_area_m2=area,
                    source_polygon_fluid_area_m2=reference_area, boundary_perimeters_m=perimeters,
                    outer_perimeter_m=perimeters['inlet'] + perimeters['outlet'],
                    polygon_preservation=dict(wall=wall_check, outer=outer_check),
                    patch_assignment_passed=True, topology_boundary_coverage_passed=True,
                    definitions=dict(coordinates='Source-coordinate metres; requires unrotated unit-chord source polygon',
                        normal_height='Triangle altitude 2*area/wall_edge_length; centroid normal distance is one third of this altitude',
                        leading_edge='Wall-edge midpoint x <= 0.05 m', trailing_edge='Wall-edge midpoint x >= 0.95 m',
                        wake='Triangle centroid 1 <= x <= 3 m and abs(y) < 0.25 m',
                        quantiles='Unweighted by edge length or triangle area; count, minimum, p10, median, p90, maximum',
                        wake_length='sqrt(triangle area), not an isotropic spacing or wall-normal spacing',
                        patch_split='Outer edge midpoint x < 0.5 m is inlet; x > 0.5 m is outlet',
                        preservation_tolerance_m=TOLERANCE_M,
                        scope='Measured nonuniform mesh sensitivity; these metrics do not establish GCI, mesh independence or physical validation'),
                    software=dict(gmsh=gmsh.__version__, numpy=np.__version__),
                    inputs={mesh_path.name: dict(path=str(mesh_path), sha256=sha256(mesh_path)),
                            'original_source_mesh': dict(path=str(source), sha256=sha256(source)),
                            'original_source_setup': dict(path=str(source_setup), sha256=sha256(source_setup))},
                    script_sha256=sha256(Path(__file__)))
    finally:
        gmsh.finalize()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    mesh, source, out = (p.resolve() for p in (args.mesh, args.source, args.out))
    if out.exists():
        raise FileExistsError('Use a fresh metrics JSON path')
    require(out.parent.is_dir(), 'Metrics output parent directory must already exist')
    result = audit(mesh, source)
    with out.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(out=str(out), cells=result['cells'], nodes=result['nodes'],
                          wall_face_count=result['wall_face_count'], checks_passed=True)))


if __name__ == '__main__':
    main()
