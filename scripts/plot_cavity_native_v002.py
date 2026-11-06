"""Paired lid-driven cavity pictures from explicitly selected native solved fields.

No solver launches. Native connectivity is mapped by complete vertex sets, not
import ordering. Tracers integrate the fixed solved velocity; no vortex is added.
"""
from pathlib import Path
import argparse
import csv
import ctypes
import gzip
import hashlib
import json
import os
import re
import shutil
import sys
import textwrap

os.environ['TEMP'] = os.environ['TMP'] = 'E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-cavity'
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import Normalize
from PIL import Image, ImageDraw

# Existing tested ASCII mesh/field readers and ID-labelled HDF block assembler.
import plot_nozzle_native_v004 as native
os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-cavity'
require, sha = native.require, native.sha
NAMES = {'fluent': 'ANSYS Fluent', 'openfoam': 'OpenFOAM'}
COLOURS = {'fluent': '#007f88', 'openfoam': '#c95d25'}
TOL = 1e-8


def resource_guard(render_reserve=False):
    class MemoryStatus(ctypes.Structure):
        _fields_ = [('length', ctypes.c_ulong), ('load', ctypes.c_ulong),
                    ('total_physical', ctypes.c_ulonglong), ('available_physical', ctypes.c_ulonglong),
                    ('total_page', ctypes.c_ulonglong), ('available_page', ctypes.c_ulonglong),
                    ('total_virtual', ctypes.c_ulonglong), ('available_virtual', ctypes.c_ulonglong),
                    ('available_extended_virtual', ctypes.c_ulonglong)]
    memory = MemoryStatus(); memory.length = ctypes.sizeof(memory)
    require(bool(ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(memory))), 'Cannot read Windows available RAM')
    ram = memory.available_physical / 1024 ** 3
    disk = shutil.disk_usage('E:/').free / 1024 ** 3
    require(ram >= (3.0 if render_reserve else 2.5) and disk >= 10.,
            f'Postprocess resource floor: available RAM {ram:.3f} GiB, E free {disk:.3f} GiB')
    return dict(available_RAM_GiB=ram, E_free_GiB=disk, required_floor_RAM_GiB=2.5, required_floor_E_GiB=10.)


class CavityGeometry:
    def __init__(self, path):
        with np.load(path, allow_pickle=False) as source:
            self.nodes = source['nodes'].copy()
            self.quads = source['quads'].astype(np.int64)
            self.centres = source['centres'].copy()
            self.areas = source['areas'].copy()
            self.ny, self.nx = map(int, source['grid_shape'])
            self.ix, self.iy = source['ix'].astype(int), source['iy'].astype(int)
            self.x_edges, self.y_edges = source['x_edges'].copy(), source['y_edges'].copy()
        require((self.nx, self.ny) == (80, 80), 'Expected shared 80x80 cavity quads')
        count = self.nx * self.ny
        require(self.nodes.shape == (81 * 81, 2) and self.quads.shape == (count, 4)
                and self.centres.shape == (count, 2) and self.areas.shape == (count,), 'Invalid geometry shapes')
        require(all(np.isfinite(a).all() for a in (self.nodes, self.centres, self.areas))
                and np.all((self.quads >= 0) & (self.quads < len(self.nodes))) and np.all(self.areas > 0), 'Invalid geometry values')
        require(np.allclose(self.x_edges, np.linspace(0, 1, 81), rtol=0, atol=TOL)
                and np.allclose(self.y_edges, np.linspace(0, 1, 81), rtol=0, atol=TOL), 'Expected full unit-square uniform mesh')
        self.polygons = self.nodes[self.quads]
        x, y = self.polygons[:, :, 0], self.polygons[:, :, 1]
        cross = x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y
        signed = cross.sum(axis=1) / 2
        centres = np.column_stack((((x + np.roll(x, -1, axis=1)) * cross).sum(axis=1),
                                   ((y + np.roll(y, -1, axis=1)) * cross).sum(axis=1))) / (6 * signed[:, None])
        require(np.all(signed != 0) and np.max(np.abs(np.abs(signed) - self.areas)) < 1e-12
                and np.max(np.abs(centres - self.centres)) < TOL, 'Canonical polygon area or centroid differs')
        slots = self.iy * self.nx + self.ix
        require(np.array_equal(np.sort(slots), np.arange(count)), 'Canonical cell indices not bijective')
        self.order = np.argsort(slots)
        gridc = self.centres[self.order].reshape(self.ny, self.nx, 2)
        self.xc, self.yc = gridc[0, :, 0], gridc[:, 0, 1]
        require(np.allclose(gridc[:, :, 0], self.xc[None, :], rtol=0, atol=TOL)
                and np.allclose(gridc[:, :, 1], self.yc[:, None], rtol=0, atol=TOL), 'Centres not rectilinear')
        self.node_lookup = {tuple(np.rint(p / TOL).astype(np.int64)): i for i, p in enumerate(self.nodes)}
        self.cell_lookup = {tuple(row): i for i, row in enumerate(np.sort(self.quads, axis=1))}
        require(len(self.node_lookup) == len(self.nodes) and len(self.cell_lookup) == count, 'Repeated source entities')

    def match(self, field):
        node_ids, error = [], 0.
        for point in field['nodes'][:, :2]:
            node = self.node_lookup.get(tuple(np.rint(point / TOL).astype(np.int64)))
            if node is None or np.linalg.norm(self.nodes[node] - point) > TOL:
                distance = np.linalg.norm(self.nodes - point, axis=1)
                node = int(np.argmin(distance))
                require(distance[node] <= TOL, 'Native node absent from canonical mesh')
            node_ids.append(node)
            error = max(error, float(np.linalg.norm(self.nodes[node] - point)))
        node_ids = np.array(node_ids, dtype=np.int64)
        cell_ids = []
        for vertices in field['cells']:
            key = tuple(np.unique(node_ids[vertices]))
            require(len(key) == 4 and key in self.cell_lookup, 'Native cell does not project to one canonical quad')
            cell_ids.append(self.cell_lookup[key])
        require(np.array_equal(np.sort(cell_ids), np.arange(len(self.quads))), 'Native vertex-set cell mapping not bijective')
        order = np.argsort(cell_ids)
        U, p = field['U'][order], field['p'][order]
        require(U.shape == (len(self.quads), 2) and p.shape == (len(self.quads),)
                and np.isfinite(U).all() and np.isfinite(p).all(), 'Invalid native saved fields')
        mean = float(np.sum(p * self.areas) / np.sum(self.areas))
        gridU = U[self.order].reshape(self.ny, self.nx, 2)
        du_dy = np.gradient(gridU[:, :, 0], self.yc, axis=0, edge_order=2)
        dv_dx = np.gradient(gridU[:, :, 1], self.xc, axis=1, edge_order=2)
        omega = dv_dx - du_dy
        mapping = dict(vertex_set_bijection=True, maximum_native_node_error_m=error,
                       source_quad_cells=len(self.quads), native_layer=field['layer'])
        if 'C' in field:
            centre_error = float(np.max(np.linalg.norm(field['C'][order, :2] - self.centres, axis=1)))
            require(centre_error <= TOL, 'Saved OpenFOAM C differs from true source polygon centres')
            mapping['maximum_saved_C_error_m'] = centre_error
        return dict(U=U, p=p, relative_p=p - mean, removed_pressure_area_mean_Pa=mean,
                    speed=np.linalg.norm(U, axis=1), gridU=gridU, omega=omega,
                    mapping=mapping, solve_evidence=field['solve_evidence'])

    def interpolate(self, point, grid):
        x, y = point
        if not self.xc[0] <= x <= self.xc[-1] or not self.yc[0] <= y <= self.yc[-1]:
            return None
        i = min(max(int(np.searchsorted(self.xc, x, side='right')) - 1, 0), self.nx - 2)
        j = min(max(int(np.searchsorted(self.yc, y, side='right')) - 1, 0), self.ny - 2)
        a, b = (x - self.xc[i]) / (self.xc[i + 1] - self.xc[i]), (y - self.yc[j]) / (self.yc[j + 1] - self.yc[j])
        return ((1 - a) * (1 - b) * grid[j, i] + a * (1 - b) * grid[j, i + 1]
                + (1 - a) * b * grid[j + 1, i] + a * b * grid[j + 1, i + 1])

    def trace(self, seed, grid, total=16.45):
        p, time = np.array(seed, float), 0.
        points, times, speeds = [p.copy()], [0.], []
        stop = 'requested travel time reached'
        for _ in range(12000):
            k1 = self.interpolate(p, grid)
            if k1 is None:
                stop = 'outside native cell-centre interpolation envelope'; break
            speed = float(np.linalg.norm(k1))
            if not speeds:
                speeds.append(speed)
            if speed < 1e-8:
                stop = 'stagnation below 1e-8 m/s'; break
            if time >= total - 1e-12:
                break
            dt = min(.035, .3 * min(np.diff(self.x_edges).min(), np.diff(self.y_edges).min()) / speed, total - time)
            k2 = self.interpolate(p + dt * k1 / 2, grid)
            k3 = None if k2 is None else self.interpolate(p + dt * k2 / 2, grid)
            k4 = None if k3 is None else self.interpolate(p + dt * k3, grid)
            if k4 is None:
                stop = 'RK4 stage reached native interpolation envelope'; break
            next_p = p + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            value = self.interpolate(next_p, grid)
            if value is None or not np.all((next_p > 0) & (next_p < 1)):
                stop = 'trajectory reached interpolation envelope or actual wall'; break
            time += dt
            p = next_p
            points.append(p.copy()); times.append(time); speeds.append(float(np.linalg.norm(value)))
        else:
            stop = 'maximum integration steps reached'
        require(len(points) >= 2 and len(points) == len(times) == len(speeds), 'Unusable native-field tracer')
        return dict(seed=list(seed), points=np.array(points), times=np.array(times), speeds=np.array(speeds), stop=stop)


def fluent_native(case, data, rho=1., mu=.01):
    with h5py.File(case, 'r') as handle:
        mesh = handle['meshes/1']; topology = mesh['cells/zoneTopology']
        require(np.asarray(topology['minId']).tolist() == [1] and len(np.asarray(topology['maxId'])) == 1, 'One Fluent fluid zone required')
        count = int(topology['maxId'][0]); groups = mesh['nodes/coords']
        total = max(int(np.asarray(d.attrs['maxId']).ravel()[0]) for d in groups.values())
        width = next(iter(groups.values())).shape[1]
        nodes = native.indexed(groups, total, width)
        require(width in (2, 3) and (width == 2 or np.ptp(nodes[:, 2]) <= TOL), 'Fluent nodes not planar')
        blocks = mesh['faces/nodes']
        nf = max(int(np.asarray(b['nnodes'].attrs.get('maxId', b.attrs.get('maxId'))).ravel()[0]) for b in blocks.values())
        faces = np.full((nf, 2), -1, np.int64)
        for block in blocks.values():
            sizes = block['nnodes']; first = int(np.asarray(sizes.attrs.get('minId', block.attrs.get('minId', [1]))).ravel()[0]) - 1
            values = np.asarray(block['nodes'][:], np.int64).reshape(-1, 2) - 1
            require(np.all(sizes[:] == 2) and np.all(faces[first:first + len(values)] == -1), 'Unsupported Fluent face block')
            faces[first:first + len(values)] = values
        require(np.all((faces >= 0) & (faces < total)), 'Invalid Fluent face nodes')
        owner = native.indexed(mesh['faces/c0'], nf).astype(np.int64) - 1
        other = native.indexed(mesh['faces/c1'], nf, partial=True); valid = np.isfinite(other)
        neighbour = other[valid].astype(np.int64) - 1
        require(np.all((owner >= 0) & (owner < count)) and np.all((neighbour >= 0) & (neighbour < count)), 'Invalid Fluent cell IDs')
        keys = np.unique(np.concatenate((np.repeat(owner, 2) * total + faces.ravel(),
                                         np.repeat(neighbour, 2) * total + faces[valid].ravel())))
        ids, sizes = np.unique(keys // total, return_counts=True)
        require(np.array_equal(ids, np.arange(count)) and np.all(sizes == 4), 'Fluent cells require four distinct vertices')
        cells = (keys % total).reshape(count, 4)
    with h5py.File(data, 'r') as handle:
        phase = handle['results/1/phase-1/cells']
        U = np.column_stack([native.indexed(phase[k], count) for k in ('SV_U', 'SV_V')])
        p = native.indexed(phase['SV_P'], count)
        for key, value in (('SV_DENSITY', rho), ('SV_MU_LAM', mu)):
            if key in phase:
                require(np.allclose(native.indexed(phase[key], count), value, rtol=1e-9, atol=1e-12), 'Saved Fluent material differs: ' + key)
        require('results/residuals/phase-1' in handle, 'Fluent initialized-only data refused')
        residuals, terminals = {}, []
        for name, group in handle['results/residuals/phase-1'].items():
            if 'iterations' not in group or 'data' not in group or not len(group['iterations']):
                continue
            iteration, row = float(group['iterations'][-1]), np.asarray(group['data'][-1]).ravel()
            require(iteration > 0 and iteration == int(iteration) and len(row) >= 2 and np.isfinite(row[:2]).all() and row[1] > 0, 'Invalid Fluent solve evidence')
            residuals[name] = float(row[0] / row[1]); terminals.append(int(iteration))
        require(len(residuals) >= 3 and len(set(terminals)) == 1, 'Incomplete Fluent actual residual history')
    return dict(nodes=nodes, cells=cells, U=U, p=p, layer='native planar two-dimensional quads',
                solve_evidence=dict(actual_terminal_iteration=terminals[0], terminal_scaled_residuals=residuals))


def foam_native(case, folder, log, rho=1.):
    field = native.foam_mesh(case); count = len(field['cells'])
    require(np.allclose(np.unique(field['nodes'][:, 2]), [-.005, .005], rtol=0, atol=TOL), 'Expected centred .01 m one-layer extrusion')
    U = native.foam_field(folder / 'U', count, 3, [0, 1, -1, 0, 0, 0, 0])
    require(np.max(np.abs(U[:, 2])) < 1e-10, 'Unexpected spanwise velocity')
    field.update(U=U[:, :2], p=rho * native.foam_field(folder / 'p', count, 1, [0, 2, -2, 0, 0, 0, 0]),
                 C=native.foam_field(folder / 'C', count, 3, [0, 1, 0, 0, 0, 0, 0]),
                 layer='one native .01 m hexahedron layer projected to quads', solve_evidence=native.foam_solve_evidence(log, folder))
    return field


def audit_status(name, path, required):
    path = path.resolve(); record = json.loads(path.read_text(encoding='utf-8-sig'))
    checks = record.get('checks', record.get('check_results'))
    require(isinstance(checks, dict) and checks and all(isinstance(v, bool) for v in checks.values()), 'Audit requires actual boolean checks')
    failed = [key for key, value in checks.items() if not value]
    aggregates = [record[key] for key in ('all_declared_checks_passed', 'all_checks_passed') if key in record]
    require(set(record.get('failed_checks', failed)) == set(failed) and aggregates
            and all(isinstance(value, bool) and value == (not failed) for value in aggregates), 'Audit aggregate or failure list inconsistent')
    manifest_path = path.parent / ('fluent-source-manifest.json' if name == 'fluent' else 'source-manifest.json')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8-sig'))
    rows = manifest.get('inputs') if name == 'fluent' else manifest.get('sources')
    require(isinstance(rows, list) and rows, 'Audit source manifest required')
    unchanged_key = 'all_inputs_unchanged' if name == 'fluent' else 'inputs_unchanged_during_analysis'
    require(manifest.get(unchanged_key) is True, 'Audit source manifest does not prove stable input hashes')
    hashes = {}
    for row in rows:
        source = str(Path(row['path']).resolve())
        if 'sha256_before' in row or 'sha256_after' in row:
            require('sha256_before' in row and 'sha256_after' in row
                    and row['sha256_before'] == row['sha256_after'], 'Audit before/after source hash differs')
            digest = row['sha256_before']
            require('sha256' not in row or row['sha256'] == digest, 'Conflicting audit source digest aliases')
        else:
            digest = row.get('sha256')
        require(isinstance(digest, str) and re.fullmatch(r'[0-9a-f]{64}', digest) and source not in hashes,
                'Invalid digest or duplicate audit source path')
        hashes[source] = digest
    require(all(str(p.resolve()) in hashes for p in required), 'Selected native fields/settings/log absent from audit source manifest')
    require(all(Path(p).is_file() and sha(Path(p)) == digest for p, digest in hashes.items()), 'Audit-bound source changed')
    return dict(path=str(path), source_manifest_path=str(manifest_path), native_source_hashes=hashes,
                all_declared_checks_passed=not failed, failed_checks=failed, checks=checks)


def footer_lines(statuses):
    lines = []
    for name, status in statuses.items():
        detail = 'declared numerical checks passed' if not status['failed_checks'] else 'MISSED: ' + '; '.join(status['failed_checks'])
        lines.extend(textwrap.wrap(NAMES[name] + ': ' + detail, width=145))
    return lines


def collection(ax, g, values, limits, cmap):
    artist = PolyCollection(g.polygons, array=values, cmap=cmap, norm=Normalize(*limits),
                            edgecolors='none', linewidths=0, antialiased=False, rasterized=True)
    ax.add_collection(artist)
    return artist


def decorate(ax, title):
    ax.plot([0, 0, 1, 1], [1, 0, 0, 1], color='#263f4c', lw=1.2)
    ax.plot([0, 1], [1, 1], color='#c15e25', lw=2.1)
    ax.annotate('', xy=(.68, 1.045), xytext=(.32, 1.045), arrowprops=dict(arrowstyle='->', lw=1.25, color='#c15e25'), annotation_clip=False)
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel='x (m)', ylabel='y (m)')
    ax.set_aspect('equal'); ax.set_title(title, loc='left', fontsize=11, pad=21, fontweight='semibold')
    ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=9)


def add_footer(fig, lines, bottom=.025):
    for i, line in enumerate(reversed(lines)):
        fig.text(.07, bottom + i * .019, line, fontsize=7.4, color='#753f25')


def verify_text_bounds(fig, out, name):
    fig.canvas.draw(); renderer = fig.canvas.get_renderer()
    figure_text = [item for item in fig.texts if item.get_visible() and item.get_text()]
    axis_labels = []
    for ax in fig.axes:
        axis_labels.extend(item for item in (ax.title, ax._left_title, ax._right_title, ax.xaxis.label, ax.yaxis.label)
                           if item.get_visible() and item.get_text())
    collisions = []
    for i, first in enumerate(figure_text):
        a = first.get_window_extent(renderer=renderer)
        for second in figure_text[i + 1:] + axis_labels:
            if a.overlaps(second.get_window_extent(renderer=renderer)):
                collisions.append([first.get_text(), second.get_text()])
    records = [dict(text=item.get_text(), bbox_figure_fraction=item.get_window_extent(renderer=renderer).transformed(fig.transFigure.inverted()).bounds)
               for item in figure_text + axis_labels]
    outside = [row['text'] for row in records if row['bbox_figure_fraction'][0] < 0 or row['bbox_figure_fraction'][1] < 0
               or row['bbox_figure_fraction'][0] + row['bbox_figure_fraction'][2] > 1
               or row['bbox_figure_fraction'][1] + row['bbox_figure_fraction'][3] > 1]
    (out / (name + '-text-bounds.json')).write_text(json.dumps(dict(figure_text_vs_titles_and_axis_labels_collisions=collisions,
        text_outside_figure=outside, bounds=records), indent=2) + '\n', encoding='utf-8')
    require(not collisions and not outside, 'Rendered header/footer/axis text collision or crop: ' + name)


def save(fig, out, name):
    verify_text_bounds(fig, out, name)
    fig.savefig(out / (name + '.png'), dpi=170, facecolor='white')
    fig.savefig(out / (name + '.svg'), dpi=170, facecolor='white')
    plt.close(fig)


def render_fields(g, fields, out, status_lines, limits):
    # Rows carry quantities; columns carry solvers with identical spatial scales.
    foot = .12 + .019 * len(status_lines)
    fig, axes = plt.subplots(3, 2, figsize=(12, 15 + .25 * len(status_lines)))
    fig.subplots_adjust(left=.07, right=.90, top=.875, bottom=foot, wspace=.28, hspace=.54)
    fig.suptitle('Lid-driven cavity | actual native CFD fields', x=.07, y=.974, ha='left', fontsize=19, fontweight='semibold')
    fig.text(.07, .945, 'Re = 100 | 1 m square | lid 1 m/s | laminar, steady SIMPLE | common 80 × 80 quads', fontsize=10, color='#425e6d')
    fig.text(.07, .925, 'Selected steady snapshots: Fluent ' + str(fields['fluent']['selected_iteration'])
             + ' / OpenFOAM ' + str(fields['openfoam']['selected_iteration']) + ' iterations; solver iteration counts are independent.', fontsize=8.7, color='#425e6d')
    for col, (name, field) in enumerate(fields.items()):
        speed = collection(axes[0, col], g, field['speed'], limits['speed'], 'viridis')
        axes[0, col].streamplot(g.xc, g.yc, field['gridU'][:, :, 0], field['gridU'][:, :, 1], color='white', density=1.2, linewidth=.65, arrowsize=.7)
        omega = collection(axes[1, col], g, field['omega'].ravel()[np.argsort(g.order)], limits['omega'], 'RdBu_r')
        pressure = collection(axes[2, col], g, field['relative_p'], limits['pressure'], 'RdBu_r')
        for row, label in enumerate(('Speed + streamlines', 'Vorticity', 'Pressure minus area mean')):
            decorate(axes[row, col], NAMES[name] + ' | ' + label)
    for row, (artist, label) in enumerate(((speed, 'Speed |U| (m/s)'), (omega, 'ωz = ∂v/∂x − ∂u/∂y (s⁻¹)'), (pressure, 'p − ⟨p⟩area (Pa)'))):
        position = axes[row, 1].get_position(); cb = fig.add_axes([.925, position.y0, .018, position.height])
        fig.colorbar(artist, cax=cb, label=label)
    fig.text(.07, .075 + .019 * len(status_lines), 'Pressure gauges differ; each field has its own area mean removed. Unclipped common colour scales.', fontsize=8.5, color='#425e6d')
    fig.text(.07, .052 + .019 * len(status_lines), 'Native values; streamlines interpolate saved U. Vorticity uses second-order centre differences, one-sided at the envelope.', fontsize=8.2, color='#425e6d')
    add_footer(fig, status_lines)
    save(fig, out, 'cavity-native-field-comparison')


def render_centrelines(g, fields, out, status_lines):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8 + .15 * len(status_lines)))
    fig.subplots_adjust(left=.08, right=.96, top=.79, bottom=.24 + .018 * len(status_lines), wspace=.28)
    fig.suptitle('Cavity centreline velocities', x=.08, y=.95, ha='left', fontsize=18, fontweight='semibold')
    fig.text(.08, .88, 'Linear interpolation between adjacent native cell-centre columns / rows; boundary values are not appended.', fontsize=9, color='#425e6d')
    curves = {}
    for name, field in fields.items():
        grid = field['gridU']
        u = np.array([np.interp(.5, g.xc, row) for row in grid[:, :, 0]])
        v = np.array([np.interp(.5, g.yc, grid[:, j, 1]) for j in range(g.nx)])
        curves[name] = dict(u_at_x_half=u, v_at_y_half=v)
        style = '--' if name == 'openfoam' else '-'
        axes[0].plot(u, g.yc, color=COLOURS[name], label=NAMES[name], lw=1.8, ls=style)
        axes[1].plot(g.xc, v, color=COLOURS[name], label=NAMES[name], lw=1.8, ls=style)
    axes[0].set(xlabel='u at x = 0.5 m (m/s)', ylabel='y (m)', ylim=(0, 1))
    axes[1].set(xlabel='x (m)', ylabel='v at y = 0.5 m (m/s)', xlim=(0, 1))
    for ax in axes:
        ax.grid(alpha=.17); ax.spines[['top', 'right']].set_visible(False); ax.legend(frameon=False, fontsize=9)
    fig.text(.08, .078 + .019 * len(status_lines), 'Solver agreement is a comparison, not benchmark, mesh-independence or physical validation.', fontsize=8.5, color='#425e6d')
    add_footer(fig, status_lines)
    save(fig, out, 'cavity-centreline-velocities')
    return curves


def render_animation(g, fields, paths, out, status_lines, speed_limits, frames=48, total=16.45):
    fig, axes = plt.subplots(1, 2, figsize=(12, 7.6 + .16 * len(status_lines)), dpi=100)
    fig.subplots_adjust(left=.07, right=.96, top=.80, bottom=.27 + .018 * len(status_lines), wspace=.25)
    fig.suptitle('Inside a lid-driven cavity', x=.07, y=.966, ha='left', fontsize=20, fontweight='semibold')
    fig.text(.07, .913, 'Re = 100 | ANSYS Fluent + OpenFOAM | particles follow each fixed computed velocity field', fontsize=10, color='#425e6d')
    clock = fig.text(.07, .862, '', fontsize=10, color='#263f4c')
    tails, heads = {}, {}
    for ax, (name, field) in zip(axes, fields.items()):
        artist = collection(ax, g, field['speed'], speed_limits, 'viridis')
        decorate(ax, NAMES[name] + ' | steady iteration ' + str(field['selected_iteration']))
        for path in paths[name]:
            ax.plot(path['points'][:, 0], path['points'][:, 1], color='white', lw=.45, alpha=.10)
        tails[name] = LineCollection([], colors='#fffbe6', linewidths=1.5, alpha=.92, zorder=6)
        ax.add_collection(tails[name]); heads[name] = ax.scatter([], [], s=15, c='#fffce7', edgecolors='#254452', linewidths=.28, zorder=7)
    cb = fig.add_axes([.28, .215 + .018 * len(status_lines), .5, .019])
    fig.colorbar(artist, cax=cb, orientation='horizontal', label='Speed |U| (m/s) | common scale')
    fig.text(.07, .12 + .019 * len(status_lines), 'Tracer travel time is integration through a fixed solved field; this is not time-resolved CFD evolution.', fontsize=9, color='#304b5a')
    fig.text(.07, .092 + .019 * len(status_lines), 'Common integration clock for both solvers. GIF replay restarts the seeds; no trajectory is forced to close.', fontsize=8.5, color='#425e6d')
    add_footer(fig, status_lines)
    images = []
    for frame in range(frames):
        if frame % 8 == 0:
            resource_guard()
        time = total * frame / (frames - 1)
        clock.set_text('Tracer travel time: ' + f'{time:05.2f}' + ' s  |  2.40 s playback; fixed native solution')
        if frame == 0:
            verify_text_bounds(fig, out, 'cavity-steady-tracers')
        for name in fields:
            segments, dots = [], []
            for path in paths[name]:
                if time > path['times'][-1] + 1e-10:
                    continue
                dots.append([np.interp(time, path['times'], path['points'][:, k]) for k in (0, 1)])
                sample = np.linspace(max(0, time - .7), time, 14)
                segments.append(np.column_stack([np.interp(sample, path['times'], path['points'][:, k]) for k in (0, 1)]))
            tails[name].set_segments(segments); heads[name].set_offsets(np.array(dots).reshape(-1, 2))
        fig.canvas.draw()
        images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()))
    gif = out / 'cavity-steady-field-tracers.gif'
    images[0].save(gif, save_all=True, append_images=images[1:], duration=50, loop=0, optimize=False, disposal=2)
    images[frames // 2].save(out / 'cavity-steady-tracers.png')
    plt.close(fig)
    with Image.open(gif) as decoded:
        require(decoded.n_frames == frames, 'Decoded GIF frame count differs')
        durations, hashes = [], []
        contact = Image.new('RGB', (images[0].width, (images[0].height + 30) * 3), '#e9eff2')
        draw = ImageDraw.Draw(contact)
        for frame in range(frames):
            decoded.seek(frame); rgb = decoded.convert('RGB')
            durations.append(decoded.info.get('duration')); hashes.append(hashlib.sha256(rgb.tobytes()).hexdigest())
            if frame in (0, frames // 2, frames - 1):
                row = (0, frames // 2, frames - 1).index(frame); top = row * (images[0].height + 30)
                draw.text((16, top + 8), f'Decoded GIF frame {frame + 1} / {frames}', fill='#203c4b')
                contact.paste(rgb, (0, top + 30))
        require(durations == [50] * frames and len(set(hashes)) == frames, 'Decoded timing or distinct-frame verification failed')
        contact.save(out / 'cavity-decoded-GIF-contact-sheet.png')
    return dict(frames=frames, frame_delay_ms=50, distinct_decoded_RGB_frames=frames, playback_duration_s=frames * .05,
                common_tracer_integration_duration_s=total, method='RK4, CFL limited, bilinear native-centre velocity interpolation',
                fixed_native_field=True, solver_time_evolution=False, forced_path_closure=False,
                boundary_values_appended=False, loop='Replay restarts identical actual-velocity seeds')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('run', 'out', 'geometry', 'fluent-case', 'fluent-data', 'fluent-audit', 'openfoam-case', 'openfoam-log', 'openfoam-audit'):
        parser.add_argument('--' + flag, type=Path, required=True)
    parser.add_argument('--openfoam-time', required=True)
    parser.add_argument('--fluent-iteration', type=int, required=True,
                        help='Explicit actual terminal native Fluent iteration; independent of OpenFOAM iteration')
    args = parser.parse_args()
    run, out, script = args.run.resolve(), args.out.resolve(), Path(__file__).resolve()
    require(out.drive.upper() == 'E:' and not out.exists(), 'Fresh output directory on E: required')
    g = CavityGeometry(args.geometry)
    setup_path = args.geometry.parent / 'setup.json'
    setup = json.loads(setup_path.read_text(encoding='utf-8-sig'))
    # Root-provided exact physics is also bound by the audit source manifests.
    rho = setup.get('density_kg_m3', setup.get('rho_kg_m3'))
    mu = setup.get('dynamic_viscosity_Pa_s', setup.get('mu_Pa_s'))
    lid = setup.get('lid_velocity_m_s', setup.get('lid_speed_m_s'))
    reynolds = setup.get('Re', setup.get('reynolds_number'))
    require(rho == 1. and mu == .01 and np.array_equal(np.atleast_1d(lid), [1., 0.]) and reynolds == 100,
            'Expected rho1, mu.01, lid(1,0), Re100 setup')
    require(setup.get('length_m') == 1. and setup.get('height_m') == 1.
            and setup.get('openfoam_extrusion_depth_m') == .01 and setup.get('operating_pressure_Pa') == 0., 'Unexpected cavity dimensions or pressure setting')
    require(setup.get('geometry_npz_sha256') == sha(args.geometry), 'Canonical geometry differs from setup hash')
    case, data = args.fluent_case.resolve(), args.fluent_data.resolve()
    foam, log = args.openfoam_case.resolve(), args.openfoam_log.resolve()
    require(float(args.openfoam_time) > 0, 'OpenFOAM initialized-only field refused')
    folder = foam / args.openfoam_time
    mesh = next((p for p in (foam / 'constant/polyMesh', foam / '0/polyMesh') if (p / 'points').is_file()), None)
    require(mesh is not None, 'Native OpenFOAM mesh required')
    required = {'fluent': [case, data, case.parent / 'final-settings.json', case.parent / 'native-physical-readback.json', case.parent / 'solver-transcript.txt'],
                'openfoam': [*(folder / k for k in ('U', 'p', 'C', 'phi')), *(mesh / k for k in ('points', 'faces', 'owner', 'neighbour', 'boundary')),
                             foam / 'system/fvSchemes', foam / 'system/fvSolution', log]}
    statuses = {name: audit_status(name, args.fluent_audit if name == 'fluent' else args.openfoam_audit, paths) for name, paths in required.items()}
    inputs = [script, Path(native.__file__).resolve(), args.geometry.resolve(), setup_path.resolve()]
    for status in statuses.values():
        inputs.extend([Path(status['path']), Path(status['source_manifest_path']), *(Path(p) for p in status['native_source_hashes'])])
    inputs = sorted(set(inputs), key=str)
    before = {str(path): sha(path) for path in inputs}
    sources = dict(fluent=fluent_native(case, data), openfoam=foam_native(foam, folder, log))
    require(args.fluent_iteration > 0 and sources['fluent']['solve_evidence']['actual_terminal_iteration'] == args.fluent_iteration,
            'Selected Fluent iteration differs from actual terminal native residual history')
    fields = {name: g.match(value) for name, value in sources.items()}
    fields['fluent']['selected_iteration'] = args.fluent_iteration
    fields['openfoam']['selected_iteration'] = int(float(args.openfoam_time))
    omega_max = max(float(np.max(np.abs(field['omega']))) for field in fields.values())
    pressure_max = max(float(np.max(np.abs(field['relative_p']))) for field in fields.values())
    limits = dict(speed=[0., max(float(f['speed'].max()) for f in fields.values())],
                  omega=[-omega_max, omega_max], pressure=[-pressure_max, pressure_max])
    require(all(np.isfinite(v).all() and v[1] > v[0] for v in limits.values()), 'Nonfinite or constant solved field')
    resources_before = resource_guard(render_reserve=True)
    seeds = [[.5, float(y)] for y in np.linspace(.16, .90, 12)] + [[.2, .25], [.8, .25], [.2, .75], [.8, .75]]
    paths = {name: [g.trace(seed, field['gridU']) for seed in seeds] for name, field in fields.items()}
    status_lines = footer_lines(statuses)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.linewidth': .7})
    out.mkdir(parents=True)
    render_fields(g, fields, out, status_lines, limits)
    curves = render_centrelines(g, fields, out, status_lines)
    animation = render_animation(g, fields, paths, out, status_lines, limits['speed'])
    with (out / 'actual-native-cell-fields.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'common_quad_0based', 'x_m', 'y_m', 'u_m_s', 'v_m_s', 'p_native_Pa', 'p_minus_area_mean_Pa', 'vorticity_z_centre_gradient_s-1'])
        for name, field in fields.items():
            omega = field['omega'].ravel()[np.argsort(g.order)]
            for i in range(len(g.quads)):
                writer.writerow([name, i, *g.centres[i], *field['U'][i], field['p'][i], field['relative_p'][i], omega[i]])
    with (out / 'actual-centrelines.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'x_or_y_m', 'u_at_x_half_m_s', 'v_at_y_half_m_s'])
        for name, curve in curves.items():
            writer.writerows(zip([name] * g.ny, g.yc, curve['u_at_x_half'], curve['v_at_y_half']))
    with gzip.open(out / 'actual-RK4-tracer-paths.csv.gz', 'wt', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'seed_0based', 'travel_time_s', 'x_m', 'y_m', 'speed_m_s'])
        for name, trajectories in paths.items():
            for index, path in enumerate(trajectories):
                writer.writerows([name, index, time, *point, speed] for time, point, speed in zip(path['times'], path['points'], path['speeds']))
    after = {str(path): sha(path) for path in inputs}
    require(before == after, 'Native source, audit or executed source changed during visualization')
    checks = dict(model='Closed unit-square lid-driven cavity', Re=100, density_kg_m3=1., dynamic_viscosity_Pa_s=.01,
                  lid_velocity_m_s=1., steady=True, actual_native_fields=True, synthetic_flow=False,
                  common_cell_count=len(g.quads), native_mapping={name: f['mapping'] for name, f in fields.items()},
                  explicitly_selected_steady_iterations={name: f['selected_iteration'] for name, f in fields.items()},
                  equal_solver_iteration_counts_required=False,
                  removed_pressure_area_means_Pa={name: f['removed_pressure_area_mean_Pa'] for name, f in fields.items()},
                  common_unclipped_colour_limits=limits, animation=animation, input_hashes_stable=True,
                  resource_floor_checks=dict(before=resources_before, after=resource_guard()),
                  native_solve_evidence={name: {k: v for k, v in f['solve_evidence'].items() if k != 'observed_solve_block'} for name, f in fields.items()},
                  numerical_audits=statuses, actual_field_ranges={name: {k: [float(f[k].min()), float(f[k].max())] for k in ('speed', 'p', 'relative_p', 'omega')} for name, f in fields.items()},
                  tracer_diagnostics={name: [dict(seed=p['seed'], samples=len(p['times']), duration_s=float(p['times'][-1]),
                       endpoint_m=p['points'][-1].tolist(), maximum_excursion_from_seed_m=float(np.max(np.linalg.norm(p['points'] - p['points'][0], axis=1))),
                       integrated_path_length_m=float(np.linalg.norm(np.diff(p['points'], axis=0), axis=1).sum()), stop=p['stop']) for p in traces] for name, traces in paths.items()},
                  quantitative_acceptance='Separate bound audits, including every missed criterion in figure footers', physical_validation=False, mesh_independence=False)
    (out / 'checks.json').write_text(json.dumps(checks, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (out / 'executed-plotter.py').write_bytes(script.read_bytes())
    (out / 'reused-native-readers.py').write_bytes(Path(native.__file__).read_bytes())
    (out / 'observed-openfoam-solve-block.txt').write_text(sources['openfoam']['solve_evidence']['observed_solve_block'], encoding='utf-8')
    (out / 'input-hash-manifest.json').write_text(json.dumps(dict(command_arguments=sys.argv[1:], inputs=[dict(path=str(p), bytes=p.stat().st_size, sha256=before[str(p)]) for p in inputs]), indent=2) + '\n', encoding='utf-8')
    caption = ('# Lid-driven cavity: actual native field visualization\n\n'
               'Closed 1 m square, lid velocity (1,0) m/s, stationary side/bottom walls, density 1 kg/m3, dynamic viscosity .01 Pa s, Re100. Both solvers compute a steady laminar field on the same 80x80 source quads; OpenFOAM uses one .01 m extruded layer.\n\n'
               'Native Fluent face connectivity and native OpenFOAM faces/owner/neighbour reconstruct complete cell vertex sets. A verified bijection maps them to canonical quads. Saved OpenFOAM C independently verifies true polygon centres. Native pressure gauges differ; each pressure field has its own spatial area mean removed. OpenFOAM kinematic pressure is multiplied by rho. Shared scales include all values.\n\n'
               'Streamlines and RK4 tracers use only bilinearly interpolated actual saved cell-centre velocity. No boundary values, synthetic vortex, forced closed orbit or solver-time evolution is inserted. Tracers stop at the interpolation envelope or wall, and use a common travel-time clock. GIF replay restarts seeds. Integration travel time differs from playback time. Vorticity is a derived second-order cell-centre gradient, one-sided at the outermost centres.\n\n'
               'Every failed numerical criterion from the hash-bound audits remains visible in the figures. A closed-domain zero boundary flux alone does not demonstrate local mass conservation. Solver agreement and attractive flow pictures do not establish benchmark accuracy, mesh independence or physical validation.\n\n'
               'Selected final steady snapshots: Fluent iteration ' + str(args.fluent_iteration) + ', OpenFOAM iteration ' + args.openfoam_time
               + '. Solver iteration counts are independent convergence indices, not a common physical clock.\n\n'
               'ANSYS Fluent and OpenFOAM computed the fields; Codex prepared the workflow and visualization. No independent user execution or reflection is claimed.\n\n' + '\n'.join(status_lines) + '\n')
    (out / 'CAPTIONS.md').write_text(caption, encoding='utf-8')
    outputs = sorted((p for p in out.iterdir() if p.is_file()), key=str)
    (out / 'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name, bytes=p.stat().st_size, sha256=sha(p)) for p in outputs], indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(output=str(out), animation=animation, colour_limits=limits, numerical_failures={k: v['failed_checks'] for k, v in statuses.items()}), indent=2))


if __name__ == '__main__':
    main()
