"""Actual paired T-junction native fields and fixed-field RK4 entry tracers.

No solver launches or fabricated values in the sparse T-domain void. Colours
identify seed origin only: no scalar concentration or mixing equation is solved.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, os, sys
os.environ['TEMP'] = os.environ['TMP'] = 'E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-tjunction'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from PIL import Image, ImageDraw
import plot_cavity_native_v003 as common
os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-tjunction'
require, sha, TOL, NAMES = common.require, common.sha, common.TOL, common.NAMES
ENTRY = {'left': '#72e9ff', 'top': '#ff86c6'}
COHORT_STARTS = (0., 4., 8., 12.)
TRAVEL_TIME = 16.45


class TGeometry:
    def __init__(self, path):
        with np.load(path, allow_pickle=False) as source:
            for key in ('nodes', 'quads', 'centres', 'areas', 'ix', 'iy', 'x_edges', 'y_edges'):
                setattr(self, key, source[key].copy())
            self.ny, self.nx = map(int, source['grid_shape'])
        self.quads = self.quads.astype(np.int64); self.ix = self.ix.astype(int); self.iy = self.iy.astype(int)
        require((self.ny, self.nx) == (80, 160) and self.nodes.shape == (4641, 2)
                and self.quads.shape == (4400, 4) and self.centres.shape == (4400, 2)
                and self.areas.shape == (4400,), 'Unexpected sparse T-junction geometry shapes')
        require(all(np.isfinite(v).all() for v in (self.nodes, self.centres, self.areas))
                and np.all((self.quads >= 0) & (self.quads < len(self.nodes))) and np.all(self.areas > 0), 'Invalid source geometry')
        require(np.allclose(self.x_edges, np.arange(-40, 121) * .01, rtol=0, atol=TOL)
                and np.allclose(self.y_edges, np.arange(-10, 71) * .01, rtol=0, atol=TOL), 'Unexpected physical T-junction extent')
        self.polygons = self.nodes[self.quads]
        x, y = self.polygons[:, :, 0], self.polygons[:, :, 1]
        cross = x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y
        signed = cross.sum(axis=1) / 2
        centre = np.column_stack((((x + np.roll(x, -1, axis=1)) * cross).sum(axis=1),
                                  ((y + np.roll(y, -1, axis=1)) * cross).sum(axis=1))) / (6 * signed[:, None])
        require(np.all(signed > 0) and np.max(np.abs(signed - self.areas)) < 1e-12
                and abs(self.areas.sum() - .44) < 1e-12 and np.max(np.abs(centre - self.centres)) < TOL, 'True polygon geometry differs')
        require(np.all((self.ix >= 0) & (self.ix < self.nx) & (self.iy >= 0) & (self.iy < self.ny))
                and len(np.unique(self.iy * self.nx + self.ix)) == 4400, 'Repeated or invalid sparse grid indices')
        self.xc, self.yc = (self.x_edges[:-1] + self.x_edges[1:]) / 2, (self.y_edges[:-1] + self.y_edges[1:]) / 2
        require(np.max(np.abs(self.centres - np.column_stack((self.xc[self.ix], self.yc[self.iy])))) < TOL, 'Native source centres differ from lattice')
        self.valid = np.zeros((self.ny, self.nx), bool); self.valid[self.iy, self.ix] = True
        expected = np.zeros_like(self.valid); expected[:20, :] = True; expected[20:, 30:50] = True
        require(np.array_equal(self.valid, expected), 'Sparse occupancy differs from exact T union')
        self.node_lookup = {tuple(np.rint(p / TOL).astype(np.int64)): i for i, p in enumerate(self.nodes)}
        self.cell_lookup = {tuple(row): i for i, row in enumerate(np.sort(self.quads, axis=1))}
        require(len(self.node_lookup) == 4641 and len(self.cell_lookup) == 4400, 'Repeated source vertices or quads')

    @staticmethod
    def inside(points):
        p = np.asarray(points); x, y = p[..., 0], p[..., 1]
        horizontal = (x >= -.4) & (x <= 1.2) & (y >= -.1) & (y <= .1)
        stem = (x >= -.1) & (x <= .1) & (y >= .1) & (y <= .7)
        return horizontal | stem

    def stencil(self, point):
        x, y = point
        if not self.inside(point) or not self.xc[0] <= x <= self.xc[-1] or not self.yc[0] <= y <= self.yc[-1]:
            return None
        i = min(max(int(np.searchsorted(self.xc, x, side='right')) - 1, 0), self.nx - 2)
        j = min(max(int(np.searchsorted(self.yc, y, side='right')) - 1, 0), self.ny - 2)
        if not self.valid[j:j + 2, i:i + 2].all():
            return None
        return j, i, (x - self.xc[i]) / (self.xc[i + 1] - self.xc[i]), (y - self.yc[j]) / (self.yc[j + 1] - self.yc[j])

    def interpolate(self, point, grid):
        stencil = self.stencil(point)
        if stencil is None:
            return None
        j, i, a, b = stencil
        result = ((1 - a) * (1 - b) * grid[j, i] + a * (1 - b) * grid[j, i + 1]
                  + (1 - a) * b * grid[j + 1, i] + a * b * grid[j + 1, i + 1])
        require(np.isfinite(result).all(), 'Interpolation touched a void or nonfinite native cell')
        return result

    def match(self, field):
        node_ids, error = [], 0.
        for point in field['nodes'][:, :2]:
            node = self.node_lookup.get(tuple(np.rint(point / TOL).astype(np.int64)))
            if node is None or np.linalg.norm(self.nodes[node] - point) > TOL:
                distance = np.linalg.norm(self.nodes - point, axis=1); node = int(np.argmin(distance))
                require(distance[node] <= TOL, 'Native node absent from source T geometry')
            node_ids.append(node); error = max(error, float(np.linalg.norm(self.nodes[node] - point)))
        node_ids = np.array(node_ids, np.int64); cell_ids = []
        for vertices in field['cells']:
            key = tuple(np.unique(node_ids[vertices]))
            require(len(key) == 4 and key in self.cell_lookup, 'Native cell does not project to one source T quad')
            cell_ids.append(self.cell_lookup[key])
        require(np.array_equal(np.sort(cell_ids), np.arange(4400)), 'Native T vertex-set mapping not bijective')
        order = np.argsort(cell_ids); U, p = field['U'][order], field['p'][order]
        require(U.shape == (4400, 2) and p.shape == (4400,) and np.isfinite(U).all() and np.isfinite(p).all(), 'Invalid actual saved fields')
        grid = np.full((self.ny, self.nx, 2), np.nan); grid[self.iy, self.ix] = U
        require(np.isfinite(grid).all(axis=2).sum() == 4400 and np.isnan(grid[~self.valid]).all(), 'Void must remain masked and unfilled')
        mapping = dict(vertex_set_bijection=True, maximum_native_node_error_m=error, source_quad_cells=4400,
                       sparse_lattice_shape=[80, 160], occupied_cells=4400, masked_void_cells=8400, native_layer=field['layer'])
        if 'C' in field:
            centre_error = float(np.max(np.linalg.norm(field['C'][order, :2] - self.centres, axis=1)))
            require(centre_error <= TOL, 'Saved C differs from true source T polygon centres')
            mapping['maximum_saved_C_error_m'] = centre_error
        return dict(U=U, p=p, speed=np.linalg.norm(U, axis=1), gridU=grid, mapping=mapping, solve_evidence=field['solve_evidence'])

    def trace(self, seed, grid, entry):
        p, time = np.array(seed, float), 0.
        points, times, speeds = [p.copy()], [0.], []
        stop = 'requested travel time reached'
        for _ in range(12000):
            k1 = self.interpolate(p, grid)
            if k1 is None:
                stop = 'outside all-native-centre interpolation envelope'; break
            speed = float(np.linalg.norm(k1))
            if not speeds:
                speeds.append(speed)
            if speed < 1e-8:
                stop = 'stagnation below 1e-8 m/s'; break
            if time >= TRAVEL_TIME - 1e-12:
                break
            dt = min(.04, .25 * .01 / speed, TRAVEL_TIME - time)
            k2 = self.interpolate(p + dt * k1 / 2, grid)
            k3 = None if k2 is None else self.interpolate(p + dt * k2 / 2, grid)
            k4 = None if k3 is None else self.interpolate(p + dt * k3, grid)
            if k4 is None:
                stop = 'RK4 stage reached void, wall or centre envelope'; break
            next_p = p + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            value = self.interpolate(next_p, grid)
            if value is None or not self.inside(np.linspace(p, next_p, 7)).all():
                stop = 'trajectory reached void, actual wall or centre envelope'; break
            time += dt; p = next_p
            points.append(p.copy()); times.append(time); speeds.append(float(np.linalg.norm(value)))
        else:
            stop = 'maximum integration steps reached'
        require(len(points) >= 8 and len(points) == len(times) == len(speeds), 'Unusable actual-field entry tracer')
        return dict(seed=list(seed), entry=entry, points=np.array(points), times=np.array(times), speeds=np.array(speeds), stop=stop)


def seeds():
    return [('left', [-.385, float(v)]) for v in np.linspace(-.075, .075, 6)] + [('top', [float(v), .685]) for v in np.linspace(-.075, .075, 6)]


def decorate(ax, title):
    walls = [[[-.4, -.1], [1.2, -.1]], [[-.4, .1], [-.1, .1]], [[-.1, .1], [-.1, .7]], [[.1, .1], [.1, .7]], [[.1, .1], [1.2, .1]]]
    ax.add_collection(LineCollection(walls, colors='#294653', linewidths=1.2, zorder=6))
    ax.plot([-.4, -.4], [-.1, .1], color=ENTRY['left'], lw=2, ls='--', zorder=7)
    ax.plot([-.1, .1], [.7, .7], color=ENTRY['top'], lw=2, ls='--', zorder=7)
    ax.plot([1.2, 1.2], [-.1, .1], color='#6d8793', lw=1, ls=':', zorder=7)
    ax.set(xlim=(-.4, 1.2), ylim=(-.1, .7), xlabel='x (m)', ylabel='y (m)')
    ax.set_aspect('equal'); ax.set_title(title, loc='left', fontsize=11, pad=10, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=9); ax.set_facecolor('#f6f8f9')


def render_fields(g, fields, out, status_lines, limits):
    fig, axes = plt.subplots(2, 2, figsize=(14, 9.3 + .15 * len(status_lines)))
    fig.subplots_adjust(left=.06, right=.90, top=.81, bottom=.17 + .018 * len(status_lines), hspace=.50, wspace=.22)
    fig.suptitle('Two inlet streams merge | actual native T-junction fields', x=.06, y=.965, ha='left', fontsize=19, fontweight='bold')
    fig.text(.06, .919, 'Steady laminar flow | Re 20 at each inlet | 4,400 common quads | 0.1 m/s inlet speeds', fontsize=10, color='#425e6d')
    fig.text(.06, .882, 'Selected solver iterations: Fluent ' + str(fields['fluent']['selected_iteration']) + ' / OpenFOAM ' + str(fields['openfoam']['selected_iteration']) + '; independent steady convergence indices.', fontsize=9, color='#425e6d')
    for col, (name, field) in enumerate(fields.items()):
        speed = common.collection(axes[0, col], g, field['speed'], limits['speed'], 'viridis')
        grid = np.ma.masked_invalid(field['gridU'])
        axes[0, col].streamplot(g.xc, g.yc, grid[:, :, 0], grid[:, :, 1], color='white', density=1.15, linewidth=.7, arrowsize=.8)
        pressure = common.collection(axes[1, col], g, field['p'], limits['pressure'], 'magma')
        decorate(axes[0, col], NAMES[name] + ' | speed + streamlines')
        decorate(axes[1, col], NAMES[name] + ' | gauge pressure')
    for row, (artist, label) in enumerate(((speed, 'Speed |U| (m/s)'), (pressure, 'Gauge pressure (Pa), outlet reference 0'))):
        box = axes[row, 1].get_position(); cax = fig.add_axes([.925, box.y0, .014, box.height])
        fig.colorbar(artist, cax=cax, label=label)
    fig.text(.06, .102 + .019 * len(status_lines), 'Native cell values with shared unclipped scales. Void cells are masked; streamlines interpolate actual saved U.', fontsize=8.5, color='#425e6d')
    fig.text(.06, .076 + .019 * len(status_lines), 'Gauge pressure retains the specified common outlet reference; no spatial-mean subtraction. No scalar mixing is simulated.', fontsize=8.5, color='#425e6d')
    common.add_footer(fig, status_lines)
    common.save(fig, out, 'tjunction-native-field-comparison')


def render_animation(g, fields, paths, out, status_lines, limits):
    fig, axes = plt.subplots(1, 2, figsize=(12, 7.2 + .16 * len(status_lines)), dpi=100)
    fig.subplots_adjust(left=.07, right=.96, top=.79, bottom=.33 + .018 * len(status_lines), wspace=.25)
    fig.suptitle('Where the two streams turn', x=.07, y=.956, ha='left', fontsize=20, fontweight='bold')
    fig.text(.07, .904, 'T-junction | particles follow actual fixed Fluent / OpenFOAM velocity fields', fontsize=10, color='#425e6d')
    clock = fig.text(.07, .858, '', fontsize=10, color='#304b5a')
    tails, heads = {}, {}
    for ax, (name, field) in zip(axes, fields.items()):
        artist = common.collection(ax, g, field['speed'], limits, 'viridis')
        decorate(ax, NAMES[name] + ' | steady iteration ' + str(field['selected_iteration']))
        for path in paths[name]:
            ax.plot(path['points'][:, 0], path['points'][:, 1], color=ENTRY[path['entry']], lw=.5, alpha=.16)
        tails[name] = LineCollection([], linewidths=1.4, alpha=.94, zorder=8); ax.add_collection(tails[name])
        heads[name] = ax.scatter([], [], s=13, edgecolors='#203e4c', linewidths=.28, zorder=9)
    handles = [Line2D([0], [0], marker='o', color=ENTRY[key], lw=1.5, label=('Left' if key == 'left' else 'Top') + ' entry seed origin') for key in ENTRY]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.065, .319), ncol=2, frameon=False, fontsize=9)
    cax = fig.add_axes([.28, .236, .5, .019]); fig.colorbar(artist, cax=cax, orientation='horizontal', label='Speed |U| (m/s) | common scale')
    fig.text(.07, .153, 'Seed colours identify entry location only; they are not concentration or a scalar-mixing result.', fontsize=9, color='#304b5a')
    fig.text(.07, .123, 'Common tracer travel clock; repeated seed injections; fixed solved fields, not time-resolved CFD evolution.', fontsize=8.5, color='#425e6d')
    fig.text(.07, .095, 'All four interpolation centres must be native fluid cells. Paths stop at walls, voids or the native-centre envelope.', fontsize=8.2, color='#425e6d')
    common.add_footer(fig, status_lines)
    images = []
    for frame in range(48):
        if frame % 8 == 0:
            common.resource_guard()
        time = TRAVEL_TIME * frame / 47
        clock.set_text(f'Tracer travel clock: {time:05.2f} s | 2.40 s playback; not solver evolution')
        for name in fields:
            segments, dots, colours = [], [], []
            for path in paths[name]:
                for start in COHORT_STARTS:
                    age = time - start
                    if not 0 <= age <= path['times'][-1]:
                        continue
                    dots.append([np.interp(age, path['times'], path['points'][:, k]) for k in (0, 1)]); colours.append(ENTRY[path['entry']])
                    sample = np.linspace(max(0, age - .55), age, 12)
                    segments.append(np.column_stack([np.interp(sample, path['times'], path['points'][:, k]) for k in (0, 1)]))
            tails[name].set_segments(segments); tails[name].set_colors(colours)
            heads[name].set_offsets(np.array(dots).reshape(-1, 2)); heads[name].set_facecolors(colours)
        if frame == 0:
            common.verify_text_bounds(fig, out, 'tjunction-steady-tracers')
        fig.canvas.draw(); images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()))
    gif = out / 'tjunction-steady-field-tracers.gif'
    images[0].save(gif, save_all=True, append_images=images[1:], duration=50, loop=0, optimize=False, disposal=2)
    images[24].save(out / 'tjunction-steady-tracers.png'); plt.close(fig)
    with Image.open(gif) as readback:
        require(readback.n_frames == 48, 'Decoded GIF frame count differs')
        durations, hashes = [], []
        contact = Image.new('RGB', (images[0].width, (images[0].height + 30) * 3), '#e9eff2'); draw = ImageDraw.Draw(contact)
        for index in range(48):
            readback.seek(index); rgb = readback.convert('RGB'); durations.append(readback.info.get('duration'))
            hashes.append(hashlib.sha256(rgb.tobytes()).hexdigest())
            if index in (0, 24, 47):
                top = (0, 24, 47).index(index) * (images[0].height + 30)
                draw.text((16, top + 8), f'Decoded GIF frame {index + 1} / 48', fill='#203c4b'); contact.paste(rgb, (0, top + 30))
        require(durations == [50] * 48 and len(set(hashes)) == 48, 'GIF timing or distinct decoded frames differ')
        contact.save(out / 'tjunction-decoded-GIF-contact-sheet.png')
    return dict(frames=48, frame_delay_ms=50, distinct_decoded_RGB_frames=48, playback_duration_s=2.4,
                common_tracer_integration_duration_s=TRAVEL_TIME, method='RK4 CFL limited, bilinear interpolation only across four actual fluid cell centres',
                fixed_native_field=True, solver_time_evolution=False, forced_path_closure=False, boundary_values_appended=False,
                entry_seed_colours_are_not_scalar_mixing=True, seed_injection_travel_times_s=list(COHORT_STARTS),
                loop='Replay restarts the actual-velocity inlet seed injections')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('run', 'out', 'geometry', 'fluent-case', 'fluent-data', 'fluent-audit', 'openfoam-case', 'openfoam-log', 'openfoam-audit'):
        parser.add_argument('--' + flag, type=Path, required=True)
    parser.add_argument('--fluent-iteration', type=int, required=True); parser.add_argument('--openfoam-time', required=True)
    args = parser.parse_args(); out, script = args.out.resolve(), Path(__file__).resolve()
    require(out.drive.upper() == 'E:' and not out.exists(), 'Fresh E: output required')
    g = TGeometry(args.geometry); setup_path = args.geometry.parent / 'setup.json'; setup = json.loads(setup_path.read_text(encoding='utf-8-sig'))
    require(setup.get('density_kg_m3') == 1. and setup.get('dynamic_viscosity_Pa_s') == .001
            and setup.get('inlet_left_velocity_m_s') == [.1, 0.] and setup.get('inlet_top_velocity_m_s') == [0., -.1]
            and setup.get('outlet_gauge_pressure_Pa') == 0. and setup.get('operating_pressure_Pa') == 0.
            and setup.get('openfoam_extrusion_depth_m') == .01 and setup.get('geometry_npz_sha256') == sha(args.geometry), 'Exact T-junction setup or geometry hash differs')
    require(all(g.stencil(seed) is not None for _, seed in seeds()), 'Entry seed outside all-native interpolation envelope')
    case, data, foam, log = args.fluent_case.resolve(), args.fluent_data.resolve(), args.openfoam_case.resolve(), args.openfoam_log.resolve()
    require(float(args.openfoam_time) > 0 and args.fluent_iteration > 0, 'Initialized-only data refused')
    folder = foam / args.openfoam_time
    mesh = next((p for p in (foam / 'constant/polyMesh', foam / '0/polyMesh') if (p / 'points').is_file()), None)
    require(mesh is not None, 'Native OpenFOAM mesh missing')
    required = {'fluent': [case, data, case.parent / 'final-settings.json', case.parent / 'native-physical-readback.json', case.parent / 'solver-transcript.txt'],
                'openfoam': [*(folder / k for k in ('U', 'p', 'C', 'phi')), *(mesh / k for k in ('points', 'faces', 'owner', 'neighbour', 'boundary')),
                             foam / 'system/fvSchemes', foam / 'system/fvSolution', log]}
    statuses = {name: common.audit_status(name, args.fluent_audit if name == 'fluent' else args.openfoam_audit, value) for name, value in required.items()}
    inputs = [script, Path(common.__file__).resolve(), Path(common.native.__file__).resolve(), args.geometry.resolve(), setup_path.resolve()]
    for status in statuses.values():
        inputs += [Path(status['path']), Path(status['source_manifest_path']), *(Path(p) for p in status['native_source_hashes'])]
    inputs = sorted(set(inputs), key=str); before = {str(p): sha(p) for p in inputs}
    sources = dict(fluent=common.fluent_native(case, data, rho=1., mu=.001), openfoam=common.foam_native(foam, folder, log, rho=1.))
    require(sources['fluent']['solve_evidence']['actual_terminal_iteration'] == args.fluent_iteration, 'Selected Fluent terminal iteration differs')
    fields = {name: g.match(value) for name, value in sources.items()}
    fields['fluent']['selected_iteration'] = args.fluent_iteration; fields['openfoam']['selected_iteration'] = int(float(args.openfoam_time))
    limits = dict(speed=[0., max(float(f['speed'].max()) for f in fields.values())],
                  pressure=[min(float(f['p'].min()) for f in fields.values()), max(float(f['p'].max()) for f in fields.values())])
    require(all(np.isfinite(v).all() and v[1] > v[0] for v in limits.values()), 'Invalid actual field colour limits')
    resources_before = common.resource_guard(render_reserve=True)
    paths = {name: [g.trace(seed, f['gridU'], entry) for entry, seed in seeds()] for name, f in fields.items()}
    status_lines = common.footer_lines(statuses)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.linewidth': .7}); out.mkdir(parents=True)
    render_fields(g, fields, out, status_lines, limits); animation = render_animation(g, fields, paths, out, status_lines, limits['speed'])
    with (out / 'actual-native-cell-fields.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'common_quad_0based', 'x_m', 'y_m', 'u_m_s', 'v_m_s', 'gauge_pressure_Pa', 'speed_m_s'])
        for name, field in fields.items():
            for i in range(4400):
                writer.writerow([name, i, *g.centres[i], *field['U'][i], field['p'][i], field['speed'][i]])
    with gzip.open(out / 'actual-RK4-tracer-paths.csv.gz', 'wt', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'entry_seed_origin', 'seed_0based', 'travel_time_s', 'x_m', 'y_m', 'speed_m_s'])
        for name, trajectories in paths.items():
            for i, path in enumerate(trajectories):
                writer.writerows([name, path['entry'], i, time, *point, speed] for time, point, speed in zip(path['times'], path['points'], path['speeds']))
    require(before == {str(p): sha(p) for p in inputs}, 'Native input, audit or executed source changed during visualization')
    checks = dict(model='Laminar T-junction merge', actual_native_fields=True, synthetic_flow=False, common_cell_count=4400,
                  native_mapping={name: f['mapping'] for name, f in fields.items()}, input_hashes_stable=True,
                  explicitly_selected_steady_iterations={name: f['selected_iteration'] for name, f in fields.items()},
                  common_unclipped_colour_limits=limits, animation=animation, numerical_audits=statuses,
                  gauge_pressure_reference='Common specified outlet gauge pressure 0Pa; raw native pressure, rho1; no mean subtraction', pressure_mean_subtracted=False,
                  void_cells_masked=True, void_field_values_fabricated=False, four_native_fluid_centres_required_for_interpolation=True,
                  scalar_mixing_simulated=False, physical_validation=False, mesh_independence=False,
                  resource_floor_checks=dict(before=resources_before, after=common.resource_guard()),
                  actual_field_ranges={name: {k: [float(f[k].min()), float(f[k].max())] for k in ('speed', 'p')} for name, f in fields.items()},
                  native_solve_evidence={name: {k: v for k, v in f['solve_evidence'].items() if k != 'observed_solve_block'} for name, f in fields.items()},
                  tracer_diagnostics={name: [dict(seed=p['seed'], entry=p['entry'], samples=len(p['times']), duration_s=float(p['times'][-1]), endpoint_m=p['points'][-1].tolist(),
                       integrated_path_length_m=float(np.linalg.norm(np.diff(p['points'], axis=0), axis=1).sum()), stop=p['stop']) for p in traces] for name, traces in paths.items()})
    (out / 'checks.json').write_text(json.dumps(checks, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (out / 'executed-plotter.py').write_bytes(script.read_bytes()); (out / 'reused-cavity-helper.py').write_bytes(Path(common.__file__).read_bytes())
    (out / 'reused-native-readers.py').write_bytes(Path(common.native.__file__).read_bytes())
    (out / 'observed-openfoam-solve-block.txt').write_text(sources['openfoam']['solve_evidence']['observed_solve_block'], encoding='utf-8')
    (out / 'input-hash-manifest.json').write_text(json.dumps(dict(command_arguments=sys.argv[1:], inputs=[dict(path=str(p), bytes=p.stat().st_size, sha256=before[str(p)]) for p in inputs]), indent=2) + '\n', encoding='utf-8')
    caption = ('# Actual native T-junction merge visualization\n\n'
               'Two equal uniform .1m/s inlet streams enter the horizontal left arm and the upper stem. Density1kg/m3, dynamic viscosity.001Pa.s, inlet-width Re20; no scalar-mixing equation. Horizontal channel x[-.4,1.2],y[-.1,.1]m, upper stem x[-.1,.1],y[.1,.7]m. Shared4400quads with area.44m2; OpenFOAM one.01m centred extrusion.\n\n'
               'Complete native vertex sets map bijectively to common source quads; saved OpenFOAM C verifies true centres. Raw pressure retains the shared outlet gauge reference0Pa. No area mean is subtracted. Common scales include every native cell value. Sparse-grid voids remain masked; there are no invented field samples in the void.\n\n'
               'RK4 uses actual bilinearly interpolated native-centre velocity only when all four surrounding fluid centres exist. Every stage and segment remains within the true union of rectangles. Paths stop at the native-centre envelope, walls or voids. Seeds lie atx=-.385m ory=.685m inside the actual interpolation envelope. Repeated injections at travel-clock times0/4/8/12s reuse these same native-velocity trajectories. Entry colours are origin labels only, not concentration, composition or mixing. The common travel clock spans16.45s in2.4s playback; the fields stay fixed. No trajectory is forced closed.\n\n'
               'Native numerical audits and every missed criterion remain explicit. Pictures and solver agreement do not establish accuracy, benchmark agreement, mesh independence or physical validation. ANSYS Fluent/OpenFOAM computed the fields; Codex prepared code/operations/visualization; no independent user execution or reflection is claimed.\n\n' + '\n'.join(status_lines) + '\n')
    (out / 'CAPTIONS.md').write_text(caption, encoding='utf-8')
    outputs = sorted((p for p in out.iterdir() if p.is_file()), key=str)
    (out / 'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name, bytes=p.stat().st_size, sha256=sha(p)) for p in outputs], indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(output=str(out), animation=animation, numerical_failures={k: v['failed_checks'] for k, v in statuses.items()}), indent=2))


if __name__ == '__main__':
    main()
