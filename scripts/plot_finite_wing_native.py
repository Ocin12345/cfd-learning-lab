"""Native finite-wing figures and steady-field tracer/camera animation.

No solver is launched. NumPy/h5py/Matplotlib/Pillow only; no SciPy, Delaunay,
analytic vortex, invented wall value or transient-flow interpretation is used.
Fluent and OpenFOAM cell order is independently reconstructed and matched to the
common source tetrahedra. Output must be a fresh directory on E:.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys

os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-finite-wing'
os.environ['TEMP'] = 'E:/CFD-Learning-Lab/temp'
os.environ['TMP'] = 'E:/CFD-Learning-Lab/temp'
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from PIL import Image

NAMES = {'fluent': 'ANSYS Fluent', 'openfoam': 'OpenFOAM'}
TOL = 1e-8


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def indexed(group, count, width=1, partial=False):
    values = np.full((count, width), np.nan)
    seen = np.zeros(count, bool)
    for dataset in group.values():
        require(isinstance(dataset, h5py.Dataset), 'Unsupported nested native HDF block')
        first = int(np.asarray(dataset.attrs['minId']).ravel()[0]) - 1
        last = int(np.asarray(dataset.attrs['maxId']).ravel()[0])
        block = np.asarray(dataset[:]).reshape(-1, width)
        require(0 <= first < last <= count and len(block) == last - first
                and not seen[first:last].any(), 'Invalid/overlapping native HDF IDs')
        values[first:last], seen[first:last] = block, True
    require(partial or seen.all(), 'Incomplete native HDF ID coverage')
    return values[:, 0] if width == 1 else values


def foam_text(path):
    return re.sub(r'/\*.*?\*/|//[^\n]*', '', path.read_text(encoding='utf-8'), flags=re.S)


def foam_list(path):
    raw = foam_text(path)
    header = re.search(r'\bFoamFile\s*\{([^}]*)\}', raw, re.S)
    require(header is not None and re.search(r'\bformat\s+ascii\s*;', header[1]),
            'Only saved ASCII OpenFOAM lists are supported: ' + str(path))
    match = re.fullmatch(r'\s*(\d+)\s*\((.*)\)\s*;?\s*', raw[header.end():], re.S)
    require(match is not None, 'Unsupported native OpenFOAM list: ' + str(path))
    return int(match[1]), match[2]


def foam_values(body, width, count, keyword='internalField'):
    kind = 'vector' if width == 3 else 'scalar'
    match = re.search(r'\b' + keyword + r'\s+nonuniform\s+List<' + kind + r'>\s+(\d+)\s*\((.*?)\)\s*;', body, re.S)
    if match:
        require(int(match[1]) == count, 'Native OpenFOAM field count differs')
        values = np.fromstring(match[2].replace('(', ' ').replace(')', ' '), sep=' ')
    else:
        match = re.search(r'\b' + keyword + r'\s+uniform\s+([^;]+);', body)
        require(match is not None, 'Native field has no supported saved ' + keyword)
        value = np.fromstring(match[1].replace('(', ' ').replace(')', ' '), sep=' ')
        require(value.size == width, 'Invalid native uniform field')
        values = np.tile(value, count)
    require(values.size == count * width and np.isfinite(values).all(), 'Invalid native field values')
    return values.reshape(count, width) if width > 1 else values


def foam_field(path, width, count, dimensions):
    raw = foam_text(path)
    require(re.search(r'\bformat\s+ascii\s*;', raw), 'OpenFOAM fields must be ASCII')
    units = re.search(r'\bdimensions\s*\[([^]]+)\]\s*;', raw)
    require(units is not None and np.array_equal(np.fromstring(units[1], sep=' '), dimensions),
            'Unexpected native field dimensions: ' + str(path))
    return foam_values(raw, width, count), raw


def dictionary_block(raw, name):
    match = re.search(r'(?<![\w])"?' + re.escape(name) + r'"?\s*\{', raw)
    require(match is not None, 'Missing native dictionary block: ' + name)
    start, depth = match.end(), 1
    for index in range(start, len(raw)):
        depth += (raw[index] == '{') - (raw[index] == '}')
        if depth == 0:
            return raw[start:index]
    raise ValueError('Unterminated native dictionary block')


def reconstruct_cells(points, faces, owner, neighbour, count):
    require(faces.ndim == 2 and faces.shape[1] == 3 and len(owner) == len(faces)
            and len(neighbour) == len(faces), 'Expected ordinary triangular native tet faces')
    require(np.all((faces >= 0) & (faces < len(points)))
            and np.all((owner >= 0) & (owner < count))
            and np.all((neighbour >= -1) & (neighbour < count)), 'Invalid native connectivity IDs')
    internal = neighbour >= 0
    incidence = np.bincount(owner, minlength=count) + np.bincount(neighbour[internal], minlength=count)
    require(np.all(incidence == 4), 'Each native tetrahedron must have four incident faces')
    keys = np.unique(np.concatenate([np.repeat(owner, 3) * len(points) + faces.ravel(),
          np.repeat(neighbour[internal], 3) * len(points) + faces[internal].ravel()]))
    ids, sizes = np.unique(keys // len(points), return_counts=True)
    require(np.array_equal(ids, np.arange(count)) and np.all(sizes == 4),
            'Each native tetrahedron must have four distinct vertices')
    return dict(points=points, tets=(keys % len(points)).reshape(count, 4), faces=faces,
                owner=owner, neighbour=neighbour)


def read_fluent_geometry(case):
    with h5py.File(case, 'r') as handle:
        mesh = handle['meshes/1']
        topology = mesh['cells/zoneTopology']
        require(np.asarray(topology['minId']).tolist() == [1]
                and len(np.asarray(topology['maxId'])) == 1, 'Expected one saved native fluid-cell zone')
        count = int(topology['maxId'][0])
        nodes = mesh['nodes/coords']
        total = max(int(np.asarray(d.attrs['maxId']).ravel()[0]) for d in nodes.values())
        points = indexed(nodes, total, 3)
        blocks = mesh['faces/nodes']
        nf = max(int(np.asarray(b['nnodes'].attrs.get('maxId', b.attrs.get('maxId'))).ravel()[0]) for b in blocks.values())
        faces = np.full((nf, 3), -1, np.int64)
        for block in blocks.values():
            sizes = block['nnodes']
            first = int(np.asarray(sizes.attrs.get('minId', block.attrs.get('minId', [1]))).ravel()[0]) - 1
            require(np.all(sizes[:] == 3), 'Native Fluent faces must be triangular')
            values = np.asarray(block['nodes'][:], dtype=np.int64).reshape(-1, 3) - 1
            require(np.all(faces[first:first + len(values)] == -1), 'Overlapping native face blocks')
            faces[first:first + len(values)] = values
        owner = indexed(mesh['faces/c0'], nf).astype(np.int64) - 1
        other = indexed(mesh['faces/c1'], nf, partial=True)
        neighbour = np.full(nf, -1, np.int64)
        internal = np.isfinite(other)
        neighbour[internal] = other[internal].astype(np.int64) - 1
        native = reconstruct_cells(points, faces, owner, neighbour, count)
    native['inputs'] = [case]
    return native


def read_fluent(case, data):
    native = read_fluent_geometry(case)
    count, nf = len(native['tets']), len(native['faces'])
    with h5py.File(data, 'r') as handle:
        require('results/residuals/phase-1' in handle, 'No actual Fluent residual iteration history; initialized data is not a computed-field source')
        terminals, residuals = [], {}
        for name in ('continuity', 'x-velocity', 'y-velocity', 'z-velocity'):
            base = 'results/residuals/phase-1/' + name
            require(base in handle and len(handle[base + '/iterations']) > 0,
                    'Missing actual Fluent residual equation history: ' + name)
            iteration = float(handle[base + '/iterations'][-1])
            row = handle[base + '/data'][-1]
            require(np.isfinite(iteration) and iteration == int(iteration) and iteration > 0
                    and len(row) >= 2 and np.isfinite(row[:2]).all() and row[1] > 0,
                    'Invalid or unsolved Fluent residual iteration history')
            value = float(row[0] / row[1])
            require(np.isfinite(value) and value >= 0, 'Nonfinite native Fluent scaled residual')
            terminals.append(int(iteration))
            residuals[name] = value
        require(len(set(terminals)) == 1, 'Native Fluent equation iteration endpoints differ')
        native['solve_evidence'] = dict(actual_terminal_iteration=terminals[0], native_terminal_scaled_residuals=residuals,
                                       initialized_only_data_refused=True, convergence_pass_required_for_render=False)
        fields = handle['results/1/phase-1/cells']
        native['U'] = np.column_stack([indexed(fields[name], count) for name in ('SV_U', 'SV_V', 'SV_W')])
        native['p'] = indexed(fields['SV_P'], count)
        for name, expected in [('SV_DENSITY', 1.), ('SV_MU_LAM', .0025)]:
            if name in fields:
                require(np.allclose(indexed(fields[name], count), expected, rtol=1e-8, atol=1e-12),
                        'Native Fluent material differs: ' + name)
        phase = handle['results/1/phase-1']
        native['face_p'] = (indexed(phase['faces/SV_P'], nf, partial=True)
                            if 'faces' in phase and 'SV_P' in phase['faces'] else None)
    native['inputs'] = [case, data]
    return native


def read_foam(case, time_name, wing_patch):
    mesh = next((p for p in (case / 'constant/polyMesh', case / '0/polyMesh') if (p / 'points').is_file()), None)
    require(mesh is not None, 'No saved native OpenFOAM polyMesh')
    files = [mesh / name for name in ('points', 'faces', 'owner', 'neighbour', 'boundary')]
    n, body = foam_list(files[0])
    points = np.fromstring(body.replace('(', ' ').replace(')', ' '), sep=' ').reshape(-1, 3)
    require(len(points) == n, 'Invalid native point count')
    nf, body = foam_list(files[1])
    rows = list(re.finditer(r'(\d+)\s*\(([^()]*)\)', body))
    require(len(rows) == nf and all(int(m[1]) == 3 for m in rows), 'Expected triangular native tet faces')
    faces = np.array([np.fromstring(m[2], sep=' ', dtype=np.int64) for m in rows])
    labels = []
    for path in files[2:4]:
        n, body = foam_list(path)
        values = np.fromstring(body, sep=' ', dtype=np.int64)
        require(len(values) == n, 'Native label count differs')
        labels.append(values)
    owner, short = labels
    neighbour = np.full(nf, -1, np.int64)
    neighbour[:len(short)] = short
    count = max(int(owner.max()), int(short.max(initial=-1))) + 1
    native = reconstruct_cells(points, faces, owner, neighbour, count)
    if time_name is None:
        candidates = []
        for p in case.iterdir():
            try:
                t = float(p.name)
            except ValueError:
                continue
            if p.is_dir() and np.isfinite(t) and t > 0 and (p / 'U').is_file() and (p / 'p').is_file():
                candidates.append((t, p.name))
        require(bool(candidates), 'No positive saved OpenFOAM iteration; prepared time-0 fields are not computed CFD')
        time_name = max(candidates)[1]
    iteration = float(time_name)
    require(np.isfinite(iteration) and iteration > 0, 'Prepared OpenFOAM time-0 fields are not computed CFD')
    evidence = None
    for log in (case / 'solver.log', case / 'log.solver', case / 'log.foamRun', case / 'log.simpleFoam'):
        if not log.is_file():
            continue
        raw_log = log.read_text(encoding='utf-8', errors='replace')
        events = list(re.finditer(r'^\s*Time\s*=\s*([-+0-9.eE]+)\s*s?\s*$', raw_log, re.M))
        for index, event in enumerate(events):
            if abs(float(event[1]) - iteration) > TOL:
                continue
            block = raw_log[event.start():events[index + 1].start() if index + 1 < len(events) else len(raw_log)]
            if re.search(r'Solving for U(?:x|y|z)?\s*,', block) and re.search(r'Solving for p\s*,', block):
                evidence = dict(actual_saved_iteration=iteration, observed_solver_log=str(log),
                    observed_solve_block=block, initialized_only_data_refused=True,
                    convergence_pass_required_for_render=False)
                break
        if evidence is not None:
            break
    require(evidence is not None, 'Missing native OpenFOAM velocity/pressure solve evidence at the saved iteration')
    native['solve_evidence'] = evidence
    folder = case / time_name
    native['U'], _ = foam_field(folder / 'U', 3, count, [0, 1, -1, 0, 0, 0, 0])
    native['p'], pressure_raw = foam_field(folder / 'p', 1, count, [0, 2, -2, 0, 0, 0, 0])
    native['inputs'] = [*files, folder / 'U', folder / 'p']
    if (folder / 'C').is_file():
        centres, _ = foam_field(folder / 'C', 3, count, [0, 1, 0, 0, 0, 0, 0])
        error = float(np.max(np.linalg.norm(centres - points[native['tets']].mean(axis=1), axis=1)))
        require(error <= TOL, 'Native C differs from reconstructed tet centroids')
        native['native_C_check_max_error_m'] = error
        native['inputs'].append(folder / 'C')
    boundary_raw = foam_list(files[4])[1]
    patch = dictionary_block(boundary_raw, wing_patch)
    first = int(re.search(r'\bstartFace\s+(\d+)', patch)[1])
    length = int(re.search(r'\bnFaces\s+(\d+)', patch)[1])
    native['wing_native_face_ids'] = np.arange(first, first + length)
    patch_field = dictionary_block(dictionary_block(pressure_raw, 'boundaryField'), wing_patch)
    native['wing_pressure_boundary_type'] = re.search(r'\btype\s+([^;]+);', patch_field)[1].strip()
    native['saved_wing_p'] = None
    if re.search(r'\bvalue\s+(?:uniform|nonuniform)', patch_field):
        native['saved_wing_p'] = foam_values(patch_field, 1, length, keyword='value')
    return native


class TetGeometry:
    """Actual common tetrahedra, native-face adjacency and barycentric walking."""
    def __init__(self, path):
        with np.load(path, allow_pickle=False) as stored:
            def get(*names):
                key = next((name for name in names if name in stored), None)
                require(key is not None, 'Missing common geometry array: ' + '/'.join(names))
                return stored[key].copy()
            self.points = get('points_xyz_m', 'nodes', 'points').astype(float)
            self.tets = get('tets', 'tetrahedra').astype(np.int64)
            self.wing_faces = get('wing_faces', 'wing_surface_triangles').astype(np.int64)
            self.wing_owner = get('wing_face_owner_cells', 'wing_owners').astype(np.int64)
        require(self.points.ndim == 2 and self.points.shape[1] == 3 and np.isfinite(self.points).all(), 'Invalid source nodes')
        require(self.tets.ndim == 2 and self.tets.shape[1] == 4 and 100 <= len(self.tets) <= 250000,
                'Unsupported source tet count/type')
        require(np.all((self.tets >= 0) & (self.tets < len(self.points))), 'Invalid source tet IDs')
        xyz = self.points[self.tets]
        matrix = np.stack([xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0], xyz[:, 3] - xyz[:, 0]], axis=-1)
        det = np.linalg.det(matrix)
        require(np.all(det > 1e-16), 'Source tets must have positive nonzero signed volume')
        self.volumes, self.inverse = det / 6, np.linalg.inv(matrix)
        self.centres, self.origins = xyz.mean(axis=1), xyz[:, 0]
        raw = self.tets[:, [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]]].reshape(-1, 3)
        unique, inverse, counts = np.unique(np.sort(raw, axis=1), axis=0, return_inverse=True, return_counts=True)
        require(np.all(np.isin(counts, [1, 2])), 'Nonmanifold common tet mesh')
        order = np.argsort(inverse, kind='stable')
        starts = np.r_[0, np.cumsum(counts[:-1])]
        internal = counts == 2
        a, b = order[starts[internal]], order[starts[internal] + 1]
        adjacent = np.full(len(raw), -1, np.int64)
        adjacent[a], adjacent[b] = b // 4, a // 4
        self.adjacent = adjacent.reshape(-1, 4)
        boundary = {tuple(row): int(order[start] // 4) for row, start, count in zip(unique, starts, counts) if count == 1}
        require(self.wing_faces.ndim == 2 and self.wing_faces.shape[1] == 3
                and self.wing_owner.shape == (len(self.wing_faces),), 'Invalid source wing surface contract')
        for face, owner in zip(np.sort(self.wing_faces, axis=1), self.wing_owner):
            require(boundary.get(tuple(face)) == owner, 'Source wing surface must match actual fluid boundary/owner')
        self.wing_vertices = self.points[self.wing_faces]
        self.bounds = np.array([self.points.min(axis=0), self.points.max(axis=0)])

    def map_native(self, native):
        require(len(native['tets']) == len(self.tets), 'Native/common tet count differs')
        lookup = {tuple(row): index for index, row in enumerate(np.rint(self.points / TOL).astype(np.int64))}
        require(len(lookup) == len(self.points), 'Source nodes collapse at mapping tolerance')
        try:
            point_map = np.array([lookup[tuple(row)] for row in np.rint(native['points'] / TOL).astype(np.int64)])
        except KeyError as error:
            raise ValueError('Native node not found in common source geometry') from error
        require(np.max(np.linalg.norm(native['points'] - self.points[point_map], axis=1)) <= TOL,
                'Native/source coordinate correspondence exceeds tolerance')
        native_tets = np.sort(point_map[native['tets']], axis=1)
        source_tets = np.sort(self.tets, axis=1)
        a, b = np.lexsort(native_tets[:, ::-1].T), np.lexsort(source_tets[:, ::-1].T)
        require(np.array_equal(native_tets[a], source_tets[b])
                and np.all(np.any(np.diff(source_tets[b], axis=0) != 0, axis=1)), 'Native/common tetrahedra are not bijective')
        mapping = np.empty(len(a), np.int64)
        mapping[a] = b
        wing_lookup = {tuple(row): index for index, row in enumerate(np.sort(self.wing_faces, axis=1))}
        wing_ids = np.full(len(self.wing_faces), -1, np.int64)
        for face_id in np.flatnonzero(native['neighbour'] < 0):
            key = tuple(sorted(point_map[native['faces'][face_id]]))
            index = wing_lookup.get(key)
            if index is not None:
                require(wing_ids[index] == -1 and mapping[native['owner'][face_id]] == self.wing_owner[index],
                        'Native wing face/owner correspondence differs')
                wing_ids[index] = face_id
        require(np.all(wing_ids >= 0), 'Incomplete saved native wing surface')
        if 'wing_native_face_ids' in native:
            require(set(native['wing_native_face_ids']) == set(wing_ids),
                    'OpenFOAM named wing patch differs from common source wing geometry')
        return mapping, wing_ids

    def match(self, native):
        mapping, wing_ids = self.map_native(native)
        u, p = np.empty_like(native['U']), np.empty_like(native['p'])
        u[mapping], p[mapping] = native['U'], native['p']
        require(np.isfinite(u).all() and np.isfinite(p).all(), 'Nonfinite native velocity/pressure')
        method = 'owner-cell pressure proxy; exact wall pressure is not claimed'
        wing_p = p[self.wing_owner].copy()
        if native.get('face_p') is not None and np.isfinite(native['face_p'][wing_ids]).all():
            wing_p, method = native['face_p'][wing_ids], 'saved native wing-face pressure'
        elif native.get('saved_wing_p') is not None:
            patch_ids = native['wing_native_face_ids']
            patch_lookup = {int(face): index for index, face in enumerate(patch_ids)}
            wing_p = native['saved_wing_p'][[patch_lookup[int(face)] for face in wing_ids]]
            method = 'saved native wing boundary pressure values'
        return dict(U=u, p=p, wing_p=wing_p, pressure_method=method,
                    native_to_common=mapping, native_wing_face_ids=wing_ids, inputs=native['inputs'])

    def bary(self, point, cell):
        value = self.inverse[cell] @ (point - self.origins[cell])
        return np.r_[1 - value.sum(), value]

    def locate(self, point):
        cell = int(np.argmin(np.sum((self.centres - point) ** 2, axis=1)))
        visited = set()
        for _ in range(10000):
            weights = self.bary(point, cell)
            if weights.min() >= -1e-9:
                return cell
            if cell in visited:
                return None
            visited.add(cell)
            cell = int(self.adjacent[cell, np.argmin(weights)])
            if cell < 0:
                return None
        return None

    def segment(self, start, end, cell):
        """Walk the straight segment through adjacent actual fluid tetrahedra.

        Crossing any external face stops the sample, including the closed wing.
        End-point location alone is insufficient: a step cannot jump a solid.
        """
        point = np.asarray(start, float).copy()
        for _ in range(2000):
            a, b = self.bary(point, cell), self.bary(end, cell)
            if b.min() >= -1e-9:
                return cell
            negative = b < -1e-9
            change = b - a
            fractions = np.full(4, np.inf)
            valid = negative & (change < 0)
            fractions[valid] = np.maximum(0., -a[valid] / change[valid])
            face = int(np.argmin(fractions))
            require(np.isfinite(fractions[face]), 'Native segment walk lost a valid face crossing')
            cell_next = int(self.adjacent[cell, face])
            if cell_next < 0:
                return None
            fraction = min(1., fractions[face] + 1e-10)
            point += fraction * (end - point)
            cell = cell_next
        return None

    def reconstruct(self, field):
        ids = self.tets.ravel()
        weights = np.bincount(ids, weights=np.repeat(self.volumes, 4), minlength=len(self.points))
        require(np.all(weights > 0), 'Unused common mesh nodes')
        nodal = np.column_stack([np.bincount(ids, weights=np.repeat(self.volumes * field['U'][:, k], 4),
                    minlength=len(self.points)) / weights for k in range(3)])
        delta = nodal[self.tets[:, 1:]] - nodal[self.tets[:, :1]]
        gradient = np.einsum('nvc,nvx->ncx', delta, self.inverse)
        field['nodal_U'], field['omega_x'] = nodal, gradient[:, 2, 1] - gradient[:, 1, 2]

    def value(self, point, cell, nodal):
        return self.bary(point, cell) @ nodal[self.tets[cell]]

    def trace(self, seed, nodal, duration, max_step=.02):
        cell = self.locate(seed)
        if cell is None:
            return dict(seed=seed.tolist(), points=np.empty((0, 3)), times=np.array([]), speeds=np.array([]), stop='seed outside located fluid')
        point, time = np.array(seed, float), 0.
        points, times, speeds = [point.copy()], [time], []
        stop = 'integration duration reached'
        for _ in range(10000):
            k1 = self.value(point, cell, nodal)
            speed = float(np.linalg.norm(k1))
            speeds.append(speed)
            if speed < 1e-7:
                stop = 'saved steady-field stagnation'
                break
            dt = min(max_step, .3 * np.cbrt(6 * self.volumes[cell]) / speed, duration - time)
            if dt <= 1e-10:
                break
            samples = []
            for target in (point + .5 * dt * k1,):
                c2 = self.segment(point, target, cell)
                if c2 is None:
                    break
                k2 = self.value(target, c2, nodal)
                target3 = point + .5 * dt * k2
                c3 = self.segment(point, target3, cell)
                if c3 is None:
                    break
                k3 = self.value(target3, c3, nodal)
                target4 = point + dt * k3
                c4 = self.segment(point, target4, cell)
                if c4 is None:
                    break
                k4 = self.value(target4, c4, nodal)
                end = point + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
                next_cell = self.segment(point, end, cell)
                if next_cell is not None:
                    samples = [end, next_cell]
            if not samples:
                stop = 'native fluid boundary reached; integration stopped without crossing'
                break
            point, cell = samples
            time += dt
            points.append(point.copy())
            times.append(time)
            if time >= duration - 1e-10:
                speeds.append(float(np.linalg.norm(self.value(point, cell, nodal))))
                break
        if len(speeds) < len(points):
            speeds.append(float(np.linalg.norm(self.value(point, cell, nodal))))
        return dict(seed=seed.tolist(), points=np.asarray(points), times=np.asarray(times),
                    speeds=np.asarray(speeds[:len(points)]), stop=stop)

    def slice(self, x_plane, field):
        xyz = self.points[self.tets]
        chosen = np.flatnonzero((xyz[:, :, 0].min(axis=1) <= x_plane)
                               & (xyz[:, :, 0].max(axis=1) >= x_plane))
        triangles, omega, vectors, centres, owners = [], [], [], [], []
        for cell in chosen:
            vertices, velocities = [], []
            uv = field['nodal_U'][self.tets[cell]]
            for vertex in np.flatnonzero(np.abs(xyz[cell, :, 0] - x_plane) <= 1e-12):
                vertices.append(xyz[cell, vertex])
                velocities.append(uv[vertex])
            for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
                xa, xb = xyz[cell, a, 0], xyz[cell, b, 0]
                if (xa - x_plane) * (xb - x_plane) < 0:
                    fraction = (x_plane - xa) / (xb - xa)
                    vertices.append(xyz[cell, a] + fraction * (xyz[cell, b] - xyz[cell, a]))
                    velocities.append(uv[a] + fraction * (uv[b] - uv[a]))
            if vertices:
                _, unique_ids = np.unique(np.rint(np.asarray(vertices) / 1e-10).astype(np.int64), axis=0, return_index=True)
                vertices = [vertices[index] for index in unique_ids]
                velocities = [velocities[index] for index in unique_ids]
            if len(vertices) not in (3, 4):
                continue
            vertices, velocities = np.asarray(vertices), np.asarray(velocities)
            centre = vertices.mean(axis=0)
            order = np.argsort(np.arctan2(vertices[:, 1] - centre[1], vertices[:, 2] - centre[2]))
            vertices, velocities = vertices[order], velocities[order]
            for j in range(1, len(vertices) - 1):
                triangles.append(vertices[[0, j, j + 1]][:, [2, 1]])
                omega.append(field['omega_x'][cell])
                owners.append(cell)
            centres.append(centre)
            vectors.append(velocities.mean(axis=0))
        require(len(triangles) > 10, 'Downstream plane has insufficient native-tet intersections')
        return dict(triangles=np.asarray(triangles), omega=np.asarray(omega), owners=np.asarray(owners),
                    centres=np.asarray(centres), U=np.asarray(vectors), x_plane=x_plane)


def display(points):
    """Display physical y vertically; labels preserve x,z,y coordinate names."""
    return np.asarray(points)[..., [0, 2, 1]]


def decorate_3d(ax, title):
    ax.set(xlim=(-.4, 3.7), ylim=(-2.7, 2.7), zlim=(-.8, .8), xlabel='$x$ (m)', ylabel='$z$ (m)', zlabel='$y$ (m)')
    ax.set_box_aspect((4.1, 5.4, 1.6))
    ax.set_title(title, fontsize=11, pad=4)
    ax.tick_params(labelsize=7, pad=0)
    ax.view_init(elev=24, azim=-55)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis._axinfo['grid']['color'] = (.5, .55, .6, .12)


def render_cp(geometry, fields, out):
    maximum = max(float(np.max(np.abs(2 * f['wing_p']))) for f in fields.values())
    maximum = max(maximum, 1e-12)
    norm, cmap = TwoSlopeNorm(vmin=-maximum, vcenter=0, vmax=maximum), plt.get_cmap('RdBu_r')
    count = len(fields)
    fig = plt.figure(figsize=(12.8 if count == 2 else 7.6, 6.5))
    for i, name in enumerate(fields):
        ax = fig.add_axes([.02 + i * .98 / count, .16, .96 / count, .73], projection='3d')
        cp = 2 * fields[name]['wing_p']
        surface = Poly3DCollection(display(geometry.wing_vertices), facecolors=cmap(norm(cp)), edgecolor='none', linewidth=0)
        ax.add_collection3d(surface)
        decorate_3d(ax, NAMES[name])
        ax.set(xlim=(-.1, 1.25), ylim=(-2.3, 2.3), zlim=(-.35, .35))
        ax.set_box_aspect((1.35, 4.6, .7))
        proxy = 'proxy' in fields[name]['pressure_method']
        ax.text2D(.5, .01, 'Owner-cell pressure proxy' if proxy else 'Saved native wing-face pressure',
                  transform=ax.transAxes, ha='center', fontsize=8)
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([.2, .09, .6, .028]), orientation='horizontal')
    cb.set_label('$C_p=(p-0)/(\\rho U_\\infty^2/2)$')
    fig.suptitle('Finite-span wing · pressure coefficient on the saved surface', y=.98, fontsize=15)
    fig.savefig(out / 'wing-surface-Cp.png', dpi=180)
    plt.close(fig)
    return [-maximum, maximum]


def render_slices(slices, out, omega_limit=None):
    maximum = max(float(np.abs(s['omega']).max()) for s in slices.values())
    limit = omega_limit if omega_limit is not None else max(maximum, 1e-12)
    require(np.isfinite(limit) and limit > 0, 'Positive finite streamwise-vorticity display limit required')
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    crossflow = max(float(np.hypot(s['U'][:, 1], s['U'][:, 2]).max()) for s in slices.values())
    scale = max(crossflow / .18, 1e-8)
    fig, axes = plt.subplots(len(slices), 1, figsize=(10.5, 3.6 + 1.95 * len(slices)), layout='constrained')
    axes = np.atleast_1d(axes)
    checks = {}
    for ax, (name, section) in zip(axes, slices.items()):
        polygons = matplotlib.collections.PolyCollection(section['triangles'], array=section['omega'],
                     cmap='RdBu_r', norm=norm, edgecolors='none', rasterized=True)
        ax.add_collection(polygons)
        centres = section['centres']
        bins = np.floor((centres[:, [2, 1]] - [-2.6, -1.]) / [.35, .25]).astype(int)
        _, indices = np.unique(bins, axis=0, return_index=True)
        select = indices[(np.abs(centres[indices, 2]) <= 2.6) & (np.abs(centres[indices, 1]) <= 1.)]
        q = ax.quiver(centres[select, 2], centres[select, 1], section['U'][select, 2], section['U'][select, 1],
                      color='#253540', angles='xy', scale_units='xy', scale=scale, width=.0023)
        if crossflow > 0:
            ax.quiverkey(q, .86, 1.04, crossflow / 2, f'{crossflow / 2:.3g} m/s', coordinates='axes', labelpos='E')
        ax.set(xlim=(-2.6, 2.6), ylim=(-1, 1), aspect='equal', xlabel='$z$ (m)', ylabel='$y$ (m)')
        ax.set_title(NAMES[name], loc='left', fontsize=11)
        ax.spines[['top', 'right']].set_visible(False)
        values = section['omega']
        checks[name] = dict(native_intersection_triangles=len(values), full_omega_x_range_per_s=[float(values.min()), float(values.max())],
                            display_saturated_triangle_fraction=float(np.mean(np.abs(values) > limit)),
                            quiver_native_samples=len(select), crossflow_speed_max_m_s=float(np.hypot(section['U'][:, 1], section['U'][:, 2]).max()))
    extend = 'both' if maximum > limit else 'neither'
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap='RdBu_r'), ax=axes, location='bottom', fraction=.065, pad=.06, aspect=55, extend=extend)
    label = 'Reconstructed streamwise vorticity $\\omega_x$ (s$^{-1}$)'
    if extend != 'neither':
        label += f' · colours saturate outside ±{limit:g}'
    cb.set_label(label)
    x = next(iter(slices.values()))['x_plane']
    fig.suptitle(f'Actual downstream cross-section at $x={x:g}$ m · steady saved field', fontsize=14)
    fig.savefig(out / 'wing-downstream-vorticity-and-crossflow.png', dpi=180)
    plt.close(fig)
    return dict(display_limits_per_s=[-limit, limit], per_solver=checks, vector_scale_shared=True)


def seed_points(geometry):
    wing = geometry.wing_vertices.reshape(-1, 3)
    lo, hi = wing.min(axis=0), wing.max(axis=0)
    z = np.r_[lo[2] + np.array([-.2, -.06, .06, .2]), [0.], hi[2] + np.array([-.2, -.06, .06, .2])]
    return np.array([[lo[0] - .25, y, span] for y in (lo[1] - .06, hi[1] + .06) for span in z])


def render_animation(geometry, fields, traces, out, frames):
    maximum = max(float(np.linalg.norm(field['U'], axis=1).max()) for field in fields.values())
    norm, cmap = Normalize(0, max(maximum, 1e-12)), plt.get_cmap('viridis')
    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    axes, dots = {}, {}
    for i, name in enumerate(fields):
        ax = fig.add_axes([.015 + i * .99 / len(fields), .14, .96 / len(fields), .74], projection='3d')
        ax.add_collection3d(Poly3DCollection(display(geometry.wing_vertices), facecolors='#bcc4ca', edgecolor='none', linewidth=0))
        for trace in traces[name]:
            if len(trace['points']) < 2:
                continue
            path = display(trace['points'])
            line = Line3DCollection(np.stack([path[:-1], path[1:]], axis=1), cmap=cmap, norm=norm, linewidths=.65, alpha=.6)
            line.set_array((trace['speeds'][:-1] + trace['speeds'][1:]) / 2)
            ax.add_collection3d(line)
        dots[name] = ax.scatter([], [], [], s=12, c=[], cmap=cmap, norm=norm, depthshade=False, edgecolor='none')
        decorate_3d(ax, NAMES[name])
        axes[name] = ax
    cb = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), cax=fig.add_axes([.23, .075, .54, .025]), orientation='horizontal')
    cb.set_label('Actual saved speed $|\\mathbf{u}|$ (m/s) · shared scale')
    fig.suptitle('Finite-span wing · steady velocity-field visualization', fontsize=16, y=.975)
    fig.text(.5, .018, 'RK4 tracers and a camera orbit through fixed computed fields · no transient-flow claim', ha='center', fontsize=9, color='#56636d')
    images = []
    for frame in range(frames):
        phase = frame / frames
        for name in fields:
            positions, speeds = [], []
            good = [trace for trace in traces[name] if len(trace['points']) >= 2 and trace['times'][-1] > 0]
            for index, trace in enumerate(good):
                travel = ((phase + index / max(1, len(good))) % 1) * trace['times'][-1]
                positions.append([np.interp(travel, trace['times'], trace['points'][:, k]) for k in range(3)])
                speeds.append(float(np.interp(travel, trace['times'], trace['speeds'])))
            if positions:
                p = display(np.array(positions))
                dots[name]._offsets3d = (p[:, 0], p[:, 1], p[:, 2])
                dots[name].set_array(np.array(speeds))
            axes[name].view_init(elev=24 + 3 * np.sin(2 * np.pi * phase), azim=-55 + 16 * np.sin(2 * np.pi * phase))
        fig.canvas.draw()
        rgb = np.ascontiguousarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3])
        images.append(Image.fromarray(rgb).quantize(colors=256))
        if frame == 0:
            fig.savefig(out / 'wing-steady-tracers.png', dpi=180)
    plt.close(fig)
    path = out / 'wing-steady-field-tracers.gif'
    images[0].save(path, save_all=True, append_images=images[1:], duration=50, loop=0, optimize=False, disposal=2)
    with Image.open(path) as image:
        delays = []
        for index in range(image.n_frames):
            image.seek(index)
            delays.append(image.info.get('duration'))
        require(image.n_frames == frames and delays == [50] * frames, 'GIF frame/timing readback failed')
    return dict(frames=frames, frame_delay_ms=50, playback_fps=20, playback_duration_s=frames / 20,
                steady_field=True, solver_transient_time=False, camera_orbit=True, speed_limits_m_s=[0, maximum],
                tracer_method='RK4 through volume-weighted native nodal velocity; every substage/final segment walks actual adjacent fluid tets and stops at a boundary',
                tracer_loop='Each tracer has a fixed phase offset along its own actual integrated trajectory; visual recycling is not flow unsteadiness')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--geometry', type=Path)
    parser.add_argument('--fluent-case', type=Path)
    parser.add_argument('--fluent-data', type=Path)
    parser.add_argument('--openfoam-case', type=Path)
    parser.add_argument('--openfoam-time')
    parser.add_argument('--wing-patch', default='wing')
    parser.add_argument('--solver', choices=('both', 'fluent', 'openfoam'), default='both')
    parser.add_argument('--slice-x', type=float, default=1.5)
    parser.add_argument('--omega-limit', type=float)
    parser.add_argument('--trace-duration', type=float, default=5.)
    parser.add_argument('--frames', type=int, default=48)
    parser.add_argument('--no-animation', action='store_true')
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    script = Path(__file__).resolve()
    executed_source, source_hash = script.read_bytes(), sha(script)
    require(out.drive.upper() == 'E:' and not out.exists(), 'Use a fresh output directory on E:')
    require(40 <= args.frames <= 60 and np.isfinite(args.trace_duration) and args.trace_duration > 0,
            'Use 40..60 frames and a positive finite tracer integration duration')
    geometry_path = (args.geometry or run / 'mesh/geometry.npz').resolve()
    setup = geometry_path.parent / 'setup.json'
    metadata = json.loads(setup.read_text(encoding='utf-8-sig'))
    for key, expected in dict(chord_m=1., span_m=4., incidence_degrees_clockwise=8.,
                             density_kg_m3=1., dynamic_viscosity_Pa_s=.0025, Re_c=400).items():
        require(metadata.get(key) == expected, 'Unsupported declared finite-wing setup: ' + key)
    require(np.array_equal(metadata.get('inlet_velocity_m_s'), [1., 0., 0.]), 'Unsupported declared inlet vector')
    geometry = TetGeometry(geometry_path)
    case = (args.fluent_case or run / 'fluent/wing-final.cas.h5').resolve()
    data = (args.fluent_data or run / 'fluent/wing-final.dat.h5').resolve()
    foam = (args.openfoam_case or run / 'openfoam').resolve()
    selected = tuple(NAMES) if args.solver == 'both' else (args.solver,)
    natives = {name: (read_fluent(case, data) if name == 'fluent'
                     else read_foam(foam, args.openfoam_time, args.wing_patch)) for name in selected}
    fields = {name: geometry.match(native) for name, native in natives.items()}
    inputs = [geometry_path, *(p for field in fields.values() for p in field['inputs'])]
    inputs.append(setup)
    before = {str(p): sha(p) for p in sorted(set(inputs), key=str)}
    for field in fields.values():
        geometry.reconstruct(field)
    slices = {name: geometry.slice(args.slice_x, field) for name, field in fields.items()}
    seeds = seed_points(geometry)
    traces = {} if args.no_animation else {name: [geometry.trace(seed, field['nodal_U'], args.trace_duration) for seed in seeds]
                                         for name, field in fields.items()}
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.linewidth': .65})
    out.mkdir(parents=True)
    cp_limits = render_cp(geometry, fields, out)
    slice_checks = render_slices(slices, out, args.omega_limit)
    animation = None if args.no_animation else render_animation(geometry, fields, traces, out, args.frames)
    with (out / 'actual-wing-pressure.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['solver', 'source_wing_face_0based', 'owner_common_cell_0based', 'x_m', 'y_m', 'z_m', 'pressure_Pa', 'Cp', 'pressure_method'])
        centres = geometry.wing_vertices.mean(axis=1)
        for name, field in fields.items():
            for index, p in enumerate(field['wing_p']):
                writer.writerow([name, index, geometry.wing_owner[index], *centres[index], p, 2 * p, field['pressure_method']])
    for name, section in slices.items():
        np.savez_compressed(out / (name + '-actual-native-tet-slice.npz'), **section)
    if traces:
        with gzip.open(out / 'actual-steady-RK4-traces.csv.gz', 'wt', encoding='utf-8', compresslevel=1, newline='') as stream:
            writer = csv.writer(stream)
            writer.writerow(['solver', 'seed_0based', 'integration_travel_time_s', 'x_m', 'y_m', 'z_m', 'speed_m_s'])
            for name, paths in traces.items():
                for index, path in enumerate(paths):
                    for time, point, speed in zip(path['times'], path['points'], path['speeds']):
                        writer.writerow([name, index, time, *point, speed])
    after = {str(p): sha(p) for p in sorted(set(inputs), key=str)}
    require(before == after, 'Native source changed during plotting; refuse a stability claim')
    checks = dict(model='Full-span steady laminar finite wing, nominal Re_c=400; independent numerical checks required',
        solvers=list(selected), paired=len(selected) == 2,
        native_solve_evidence={name: {key: value for key, value in native['solve_evidence'].items() if key != 'observed_solve_block'}
                               for name, native in natives.items()},
        source_tets=len(geometry.tets), source_nodes=len(geometry.points), source_wing_faces=len(geometry.wing_faces),
        minimum_source_tet_volume_m3=float(geometry.volumes.min()), native_to_source_tets_bijective=True,
        node_matching_tolerance_m=TOL, actual_native_fields=True, input_hashes_stable_after_loading=True,
        pressure_methods={name: field['pressure_method'] for name, field in fields.items()}, Cp_shared_limits=cp_limits,
        native_field_ranges={name: dict(speed_m_s=[float(np.linalg.norm(field['U'], axis=1).min()), float(np.linalg.norm(field['U'], axis=1).max())],
            pressure_Pa=[float(field['p'].min()), float(field['p'].max())]) for name, field in fields.items()},
        downstream_slice=slice_checks, animation=animation,
        traces={name: [dict(seed=path['seed'], samples=len(path['points']), integration_time_s=float(path['times'][-1]) if len(path['times']) else 0,
                      stop=path['stop']) for path in paths] for name, paths in traces.items()},
        velocity_reconstruction='Volume-weighted averages of actual incident fluid-cell vectors to original native mesh vertices, then piecewise-linear barycentric interpolation; no boundary vectors appended',
        gradient_reconstruction='Piecewise linear curl on original saved fluid tetrahedra; no cross-solid Delaunay connectivity',
        synthetic_vortex=False, physical_validation=False, mesh_independence=False, steady_convergence_acceptance='Separate native analyser evidence; not inferred from rendered traces')
    (out / 'checks.json').write_text(json.dumps(checks, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    if 'openfoam' in natives:
        (out / 'observed-openfoam-solve-block.txt').write_text(natives['openfoam']['solve_evidence']['observed_solve_block'], encoding='utf-8')
    require(sha(script) == source_hash, 'Plotter source changed during execution')
    (out / 'executed-plotter.py').write_bytes(executed_source)
    (out / 'input-hash-manifest.json').write_text(json.dumps(dict(command_arguments=sys.argv[1:], plotter_sha256=source_hash,
        inputs=[dict(path=str(p), bytes=p.stat().st_size, sha256=before[str(p)]) for p in sorted(set(inputs), key=str)]), indent=2) + '\n', encoding='utf-8')
    caption = ('# Actual native finite-wing visualization\n\n'
        'Full-span rectangular wing, chord 1 m, span 4 m, nominal AR=4, clockwise 8 degree pitch about the span axis; '
        'NACA0012-like closed rounded-tail profile, not an exact NACA0012 claim. Declared steady incompressible laminar '
        'conditions are U=(1,0,0) m/s, rho=1 kg/m3 and mu=0.0025 Pa.s, nominal Re_c=400. Native setup, residual, mass and '
        'force-stability evidence is assessed separately; rendered appearance does not establish convergence.\n\n'
        'Every native tetrahedron is reconstructed from saved face/owner/neighbour IDs and matched bijectively to '
        'the common source geometry; imported ordering is never assumed. Saved velocity and pressure are finite. '
        'For traces, actual cell vectors are volume-weighted to incident original vertices and interpolated linearly '
        'inside original fluid tetrahedra. RK4 substage and final-point segments walk adjacent tetrahedra; crossing '
        'a fluid boundary stops integration. No streamline is permitted to jump through the wing. No analytic vortex '
        'or imposed boundary sample enters reconstruction. These are sparse-mesh interpolation/gradient approximations.\n\n'
        'When generated, the GIF depicts fixed computed steady fields: particle phase progression and a quiet common camera orbit '
        'provide motion. This is steady-field visualization, not simulated transient vortex development. Recycled '
        'tracers use fixed offsets along their own actual integrated paths; integration travel time and GIF playback '
        'time are separate from any solver iteration counter. Actual weak or absent tip roll-up is retained. Both '
        'codes in any paired view share speed limits, camera, axes in metres, and visible-domain crop. Physical y is displayed vertically '
        'using the labelled plotting-coordinate permutation (x,z,y). Generated RK4 paths/time samples are exported.\n\n'
        'Surface Cp=(p-0)/(rho*U^2/2) uses saved wing-face/boundary values when available. Each panel labels an '
        'owner-cell pressure proxy if native face values were unavailable; proxy pressure is not exact wall pressure. '
        'The downstream y-z section is made by intersecting actual fluid tetrahedra, splitting only actual intersection '
        'quadrilaterals into triangles. Signed streamwise vorticity is dw/dy-dv/dz of the linear reconstructed field. '
        'Sparse arrows show actual interpolated crossflow at native intersection centroids, with one shared arrow '
        'scale/key. Any user-requested vorticity colour saturation is explicit in the colorbar and checks.\n\n'
        'Native source hashes, executed source, interpolation/trace diagnostics, actual surface pressures and slice '
        'arrays accompany these figures. No aircraft-performance, benchmark-accuracy, mesh/time independence or '
        'physical-validation claim is made. Fluent and OpenFOAM computed the fields; Codex prepared this workflow.\n')
    (out / 'CAPTIONS.md').write_text(caption, encoding='utf-8')
    files = sorted([p for p in out.rglob('*') if p.is_file()], key=str)
    (out / 'output-hash-manifest.json').write_text(json.dumps([dict(path=str(p.relative_to(out)), bytes=p.stat().st_size,
        sha256=sha(p)) for p in files], indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(output=str(out), cells=len(geometry.tets), pressure_methods=checks['pressure_methods'],
                          animation=animation, downstream_slice=slice_checks), indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
