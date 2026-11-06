"""Paired actual native square-obstacle steady wake; fixed-field RK4 tracers.

Finite confined Re 20 model with symmetry/slip top and bottom. No transient
vortex shedding, scalar mixing, published drag benchmark or validation claim.
Tracer travel/injection parameters require a separate actual-field timing probe.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, json, os, sys
os.environ['TEMP'] = os.environ['TMP'] = 'E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-square-wake'
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from PIL import Image, ImageDraw
import plot_elbow_native_v001 as elbow
common = elbow.common
require, sha, TOL, NAMES = common.require, common.sha, common.TOL, common.NAMES
os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-square-wake'
COLOURS = {'fluent': '#00879c', 'openfoam': '#e06421'}
FROZEN_HELPERS = {'elbow': 'fc8fc0272c6f539cbe19533ef7efbef265009435643efbb08a3685fb0773b7a1',
                  'cavity': 'c2425654f684bc8de8d9219dffa872d8b3c80a605bb35fce2b9070362f388445',
                  'readers': '4596d6a6b975abc6110d80fdb8b5b680878cdff18022ce969866b288aacb2f02'}
TRAVEL_TIME = None
COHORT_STARTS = ()


class WakeGeometry(elbow.ElbowGeometry):
    """Use the frozen elbow's four-native-centre stencil/interpolation methods."""
    def __init__(self, path):
        with np.load(path, allow_pickle=False) as source:
            for key in ('nodes', 'quads', 'centres', 'areas', 'ix', 'iy', 'x_edges', 'y_edges'):
                setattr(self, key, source[key].copy())
            self.ny, self.nx = map(int, source['grid_shape'])
        self.quads = self.quads.astype(np.int64); self.ix = self.ix.astype(int); self.iy = self.iy.astype(int)
        require((self.ny, self.nx) == (80, 200) and self.nodes.shape == (16200, 2)
                and self.quads.shape == (15900, 4) and self.centres.shape == (15900, 2)
                and self.areas.shape == (15900,), 'Unexpected square-hole geometry shapes')
        require(all(np.isfinite(v).all() for v in (self.nodes, self.centres, self.areas))
                and np.all((self.quads >= 0) & (self.quads < len(self.nodes))) and np.all(self.areas > 0), 'Invalid source geometry')
        require(np.allclose(self.x_edges, np.arange(201) * .01, rtol=0, atol=TOL)
                and np.allclose(self.y_edges, np.arange(-40, 41) * .01, rtol=0, atol=TOL), 'Unexpected physical extent')
        self.polygons = self.nodes[self.quads]
        x, y = self.polygons[:, :, 0], self.polygons[:, :, 1]
        cross = x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y; signed = cross.sum(axis=1) / 2
        centre = np.column_stack((((x + np.roll(x, -1, axis=1)) * cross).sum(axis=1),
                                  ((y + np.roll(y, -1, axis=1)) * cross).sum(axis=1))) / (6 * signed[:, None])
        require(np.all(signed > 0) and np.max(np.abs(signed - self.areas)) < 1e-12
                and abs(self.areas.sum() - 1.59) < 1e-12 and np.max(np.abs(centre - self.centres)) < TOL, 'True polygon geometry differs')
        require(np.all((self.ix >= 0) & (self.ix < self.nx) & (self.iy >= 0) & (self.iy < self.ny))
                and len(np.unique(self.iy * self.nx + self.ix)) == 15900, 'Repeated or invalid sparse indices')
        self.xc, self.yc = (self.x_edges[:-1] + self.x_edges[1:]) / 2, (self.y_edges[:-1] + self.y_edges[1:]) / 2
        require(np.max(np.abs(self.centres - np.column_stack((self.xc[self.ix], self.yc[self.iy])))) < TOL, 'Source centres differ from lattice')
        self.valid = np.zeros((80, 200), bool); self.valid[self.iy, self.ix] = True
        expected = np.ones((80, 200), bool); expected[35:45, 45:55] = False
        require(np.array_equal(self.valid, expected), 'Occupancy must contain the exact 10 by 10 square hole')
        self.node_lookup = {tuple(np.rint(p / TOL).astype(np.int64)): i for i, p in enumerate(self.nodes)}
        self.cell_lookup = {tuple(row): i for i, row in enumerate(np.sort(self.quads, axis=1))}
        require(len(self.node_lookup) == 16200 and len(self.cell_lookup) == 15900, 'Repeated source vertices or cells')

    @staticmethod
    def inside(points):
        p = np.asarray(points); x, y = p[..., 0], p[..., 1]
        body = (x >= .45) & (x <= .55) & (y >= -.05) & (y <= .05)
        return (x >= 0.) & (x <= 2.) & (y >= -.4) & (y <= .4) & ~body

    def segment_fluid(self, a, b):
        """Exact segment/closed-body-box rejection, including boundary contact."""
        if not self.inside(np.array([a, b])).all():
            return False
        lower, upper = np.array([.45, -.05]), np.array([.55, .05])
        t0, t1 = 0., 1.
        for q, delta, lo, hi in zip(a, np.asarray(b) - a, lower, upper):
            if abs(delta) < 1e-15:
                if q < lo or q > hi:
                    return True
            else:
                r0, r1 = sorted(((lo - q) / delta, (hi - q) / delta))
                t0, t1 = max(t0, r0), min(t1, r1)
                if t0 > t1:
                    return True
        return False

    def match(self, field):
        ids, error = [], 0.
        for point in field['nodes'][:, :2]:
            node = self.node_lookup.get(tuple(np.rint(point / TOL).astype(np.int64)))
            if node is None or np.linalg.norm(self.nodes[node] - point) > TOL:
                distance = np.linalg.norm(self.nodes - point, axis=1); node = int(np.argmin(distance))
                require(distance[node] <= TOL, 'Native node absent from source square-hole geometry')
            ids.append(node); error = max(error, float(np.linalg.norm(self.nodes[node] - point)))
        ids = np.array(ids, np.int64); cells = []
        for vertices in field['cells']:
            key = tuple(np.unique(ids[vertices])); require(len(key) == 4 and key in self.cell_lookup, 'Native cell does not project to one source quad')
            cells.append(self.cell_lookup[key])
        require(np.array_equal(np.sort(cells), np.arange(15900)), 'Native complete vertex-set mapping not bijective')
        order = np.argsort(cells); U, p = field['U'][order], field['p'][order]
        require(U.shape == (15900, 2) and p.shape == (15900,) and np.isfinite(U).all() and np.isfinite(p).all(), 'Invalid actual native fields')
        grid = np.full((80, 200, 2), np.nan); grid[self.iy, self.ix] = U
        require(np.isfinite(grid).all(axis=2).sum() == 15900 and np.isnan(grid[~self.valid]).all(), 'Hole must remain masked and unfilled')
        mapping = dict(vertex_set_bijection=True, maximum_native_node_error_m=error, source_quad_cells=15900,
                       sparse_lattice_shape=[80, 200], occupied_cells=15900, masked_void_cells=100, native_layer=field['layer'])
        if 'C' in field:
            centre_error = float(np.max(np.linalg.norm(field['C'][order, :2] - self.centres, axis=1)))
            require(centre_error <= TOL, 'Saved C differs from actual quad centroids'); mapping['maximum_saved_C_error_m'] = centre_error
        return dict(U=U, p=p, speed=np.linalg.norm(U, axis=1), gridU=grid, mapping=mapping, solve_evidence=field['solve_evidence'])

    def trace(self, seed, grid, origin):
        p, time = np.array(seed, float), 0.; points, times, speeds = [p.copy()], [0.], []
        stop = 'requested travel time reached'
        require(TRAVEL_TIME is not None, 'Travel clock requires actual-field timing selection')
        for _ in range(12000):
            k1 = self.interpolate(p, grid)
            if k1 is None:
                stop = 'outside four-native-fluid-centre interpolation envelope'; break
            speed = float(np.linalg.norm(k1))
            if not speeds:
                speeds.append(speed)
            if speed < 1e-8:
                stop = 'actual speed below 1e-8 m/s'; break
            if time >= TRAVEL_TIME - 1e-12:
                break
            dt = min(.04, .25 * .01 / speed, TRAVEL_TIME - time)
            p2 = p + dt * k1 / 2; k2 = self.interpolate(p2, grid) if self.segment_fluid(p, p2) else None
            p3 = None if k2 is None else p + dt * k2 / 2
            k3 = None if p3 is None or not self.segment_fluid(p, p3) else self.interpolate(p3, grid)
            p4 = None if k3 is None else p + dt * k3
            k4 = None if p4 is None or not self.segment_fluid(p, p4) else self.interpolate(p4, grid)
            if k4 is None:
                stop = 'RK4 stage reached body, wall or native-centre envelope'; break
            next_p = p + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
            value = self.interpolate(next_p, grid)
            if value is None or not self.segment_fluid(p, next_p):
                stop = 'trajectory reached body, wall or native-centre envelope'; break
            time += dt; p = next_p; points.append(p.copy()); times.append(time); speeds.append(float(np.linalg.norm(value)))
        else:
            stop = 'maximum integration steps reached'
        require(len(points) == len(times) == len(speeds) and len(points) >= 1, 'Native-field trajectory arrays differ')
        return dict(seed=list(seed), origin=origin, points=np.array(points), times=np.array(times), speeds=np.array(speeds), stop=stop)


def inlet_seeds():
    return [('inlet', [.015, y]) for y in (-.32, -.22, -.12, -.065, -.025, .025, .065, .12, .22, .32)]


def wake_seeds():
    return [('wake', [x, y]) for x in (.585, .65, .75) for y in (-.025, .025)]


def decorate(ax, title, closeup=False):
    ax.add_patch(Rectangle((.45, -.05), .1, .1, facecolor='#e3e7e9', edgecolor='#263e4d', linewidth=1.1, zorder=8))
    if closeup:
        ax.set(xlim=(.35, 1.15), ylim=(-.2, .2))
    else:
        ax.plot([0, 2], [.4, .4], color='#5c7684', ls='--', lw=.9); ax.plot([0, 2], [-.4, -.4], color='#5c7684', ls='--', lw=.9)
        ax.set(xlim=(0, 2), ylim=(-.4, .4))
    ax.set(xlabel='x (m)', ylabel='y (m)'); ax.set_aspect('equal'); ax.set_title(title, loc='left', fontsize=11, pad=10, fontweight='bold')
    ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=9); ax.set_facecolor('#f6f8f9')


def render_fields(g, fields, paths, out, notes, limits):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10.2 + .16 * len(notes)))
    fig.subplots_adjust(left=.075, right=.89, top=.81, bottom=.19 + .018 * len(notes), hspace=.44, wspace=.28)
    fig.suptitle('A square obstacle in a confined laminar stream', x=.075, y=.965, ha='left', fontsize=19, fontweight='bold')
    fig.text(.075, .919, 'Actual steady fields | Re 20 | 15,900 fluid quads | D = 0.1 m | inlet 0.1 m/s', fontsize=10, color='#425e6d')
    fig.text(.075, .878, 'Symmetry/slip top and bottom | 12.5% blockage | selected steady indices: ' + str(fields['fluent']['selected_iteration']) + ' / ' + str(fields['openfoam']['selected_iteration']), fontsize=9, color='#425e6d')
    for col, (name, f) in enumerate(fields.items()):
        speed = common.collection(axes[0, col], g, f['speed'], limits['speed'], 'viridis')
        axes[0, col].add_collection(LineCollection([p['points'] for p in paths[name] if len(p['points']) > 1], colors='white', linewidths=.6, alpha=.8, zorder=6))
        pressure = common.collection(axes[1, col], g, f['p'], limits['pressure'], 'magma')
        decorate(axes[0, col], NAMES[name] + ' | speed + native-field paths'); decorate(axes[1, col], NAMES[name] + ' | gauge pressure')
    for row, (artist, label) in enumerate(((speed, 'Speed |U| (m/s)'), (pressure, 'Gauge pressure (Pa)\noutlet reference 0'))):
        box = axes[row, 1].get_position(); cax = fig.add_axes([.917, box.y0, .014, box.height]); fig.colorbar(artist, cax=cax, label=label)
    fig.text(.075, .125 + .019 * len(notes), 'Shared unclipped native scales; the actual body hole contains no flow samples. Raw outlet-gauge pressure is retained.', fontsize=8.4, color='#425e6d')
    fig.text(.075, .097 + .019 * len(notes), 'Paths integrate the fixed native velocity field. These are steady solutions; no transient vortex shedding is simulated.', fontsize=8.4, color='#425e6d')
    common.add_footer(fig, notes); common.save(fig, out, 'square-wake-native-field-comparison')


def centreline(g, f):
    xs = g.xc[g.xc > .55]; samples = np.array([g.interpolate([x, 0.], f['gridU']) for x in xs])
    require(samples.shape == (len(xs), 2) and np.isfinite(samples).all(), 'Centreline requires four native fluid centres')
    pgrid = np.full((80, 200), np.nan); pgrid[g.iy, g.ix] = f['p']
    pressure = np.array([g.interpolate([x, 0.], pgrid) for x in xs])
    return xs, samples, pressure


def render_wake(g, fields, out, notes, limits):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10.7 + .16 * len(notes)))
    fig.subplots_adjust(left=.085, right=.91, top=.81, bottom=.19 + .018 * len(notes), hspace=.53, wspace=.33)
    fig.suptitle('The saved steady wake | native vectors and centreline', x=.085, y=.965, ha='left', fontsize=18, fontweight='bold')
    fig.text(.085, .914, 'Arrows retain native velocity magnitudes; cyan dashes mark u = 0 where present. Grey square is solid.', fontsize=9.5, color='#425e6d')
    diagnostics = {}
    for col, (name, f) in enumerate(fields.items()):
        ax = axes[0, col]; artist = common.collection(ax, g, f['speed'], limits['speed'], 'viridis')
        sample = (g.centres[:, 0] >= .35) & (g.centres[:, 0] <= 1.15) & (np.abs(g.centres[:, 1]) <= .2) & (g.ix % 3 == 0) & (g.iy % 3 == 0)
        q = ax.quiver(*g.centres[sample].T, *f['U'][sample].T, color='white', angles='xy', scale_units='xy', scale=1.5, width=.0034, zorder=7)
        ax.quiverkey(q, .71, 1.10, .1, '0.1 m/s', coordinates='axes', labelpos='E')
        grid = np.ma.masked_invalid(f['gridU'][:, :, 0])
        if np.min(f['U'][:, 0]) < 0 < np.max(f['U'][:, 0]):
            ax.contour(g.xc, g.yc, grid, levels=[0.], colors=['#90ffff'], linestyles='--', linewidths=.8)
        decorate(ax, NAMES[name] + ' | actual wake vectors', closeup=True)
        x, U, p = centreline(g, f); style = '--' if name == 'openfoam' else '-'
        axes[1, 0].plot((x - .55) / .1, U[:, 0] / .1, color=COLOURS[name], ls=style, lw=1.8, label=NAMES[name])
        axes[1, 1].plot((x - .55) / .1, p, color=COLOURS[name], ls=style, lw=1.8, label=NAMES[name])
        reversed_cells = (g.centres[:, 0] > .55) & (np.abs(g.centres[:, 1]) < .15) & (f['U'][:, 0] < 0)
        crossing = np.flatnonzero((U[:-1, 0] < 0) & (U[1:, 0] >= 0))
        zeros = [float(x[i] - U[i, 0] * (x[i + 1] - x[i]) / (U[i + 1, 0] - U[i, 0])) for i in crossing]
        diagnostics[name] = dict(downstream_negative_u_cells=int(reversed_cells.sum()), downstream_negative_u_area_m2=float(g.areas[reversed_cells].sum()),
            negative_to_nonnegative_centreline_zero_crossings_x_m=zeros, definition='Native u<0 for x>0.55 m, |y|<0.15 m; no benchmark or time-dependent wake claim')
    box = axes[0, 1].get_position(); cax = fig.add_axes([.932, box.y0, .012, box.height]); fig.colorbar(artist, cax=cax, label='Speed (m/s)')
    axes[1, 0].set(xlabel='Distance from rear face / D', ylabel='Centreline u / inlet U'); axes[1, 0].axhline(0, color='#566d7a', lw=.8)
    axes[1, 1].set(xlabel='Distance from rear face / D', ylabel='Centreline gauge\npressure (Pa)')
    for ax in axes[1]:
        ax.grid(alpha=.17); ax.legend(frameon=False, fontsize=9); ax.spines[['top', 'right']].set_visible(False)
    fig.text(.085, .120 + .019 * len(notes), 'Centreline profiles use bilinear interpolation across four actual fluid centres. Zero crossings are exported as diagnostics.', fontsize=8.3, color='#425e6d')
    fig.text(.085, .092 + .019 * len(notes), 'Finite domain: body centre 5 D from inlet and 15 D from outlet, domain height 8 D. No unbounded-flow accuracy claim.', fontsize=8.3, color='#425e6d')
    common.add_footer(fig, notes); common.save(fig, out, 'square-wake-closeup-and-profiles')
    return diagnostics


def render_diagnostics(histories, out, notes):
    fig, axes = plt.subplots(2, 2, figsize=(12, 10.6 + .16 * len(notes)))
    fig.subplots_adjust(left=.095, right=.96, top=.82, bottom=.19 + .018 * len(notes), hspace=.45, wspace=.3)
    fig.suptitle('Actual native force and convergence histories', x=.095, y=.965, ha='left', fontsize=18, fontweight='bold')
    fig.text(.095, .915, 'Steady iteration indices, not physical time | Cd and Cl use the same rho, inlet speed and D; native depth is retained.', fontsize=9, color='#425e6d')
    for col, (name, h) in enumerate(histories.items()):
        style = '--' if name == 'openfoam' else '-'
        for ax, key in zip(axes[0], ('Cd', 'Cl')):
            ax.plot(h['force_iteration'], h[key], color=COLOURS[name], ls=style, lw=1.6, label=NAMES[name])
        for index, equation in enumerate(h['residual_names']):
            vals = h['residual_values'][:, index]; axes[1, col].semilogy(h['residual_iteration'], np.ma.masked_less_equal(vals, 0), lw=1.1, label=str(equation))
        axes[1, col].axhline(1e-6, color='#5f7180', ls=':', lw=.9, label='Declared terminal target')
        axes[1, col].set(title=NAMES[name] + ' | own native residual norm', xlabel='Steady iteration', ylabel='Native residual\n(own solver norm)')
    axes[0, 0].set(xlabel='Steady iteration', ylabel='Drag coefficient Cd'); axes[0, 1].set(xlabel='Steady iteration', ylabel='Lift coefficient Cl')
    for ax in axes.ravel():
        ax.grid(alpha=.17); ax.spines[['top', 'right']].set_visible(False); ax.legend(frameon=False, fontsize=8.5)
    fig.text(.095, .120 + .019 * len(notes), 'Force history stability and lift symmetry are numerical diagnostics. No published Cd is assumed; solver agreement is not validation.', fontsize=8.3, color='#425e6d')
    fig.text(.095, .092 + .019 * len(notes), 'Residual norms differ between solvers. Any native zero residual is retained in the CSV and omitted from the logarithmic curve.', fontsize=8.3, color='#425e6d')
    common.add_footer(fig, notes); common.save(fig, out, 'square-wake-force-and-residuals')


def render_animation(g, fields, paths, out, notes, limits):
    fig, axes = plt.subplots(2, 1, figsize=(12, 10.7 + .16 * len(notes)), dpi=100)
    fig.subplots_adjust(left=.10, right=.95, top=.80, bottom=.32 + .018 * len(notes), hspace=.49)
    fig.suptitle('A steady stream splits around a square', x=.07, y=.965, ha='left', fontsize=21, fontweight='bold')
    fig.text(.07, .918, 'Actual fixed Fluent / OpenFOAM velocity fields | Re 20 | 12.5% blockage | no transient shedding', fontsize=10, color='#425e6d')
    clock = fig.text(.07, .872, '', fontsize=10, color='#304b5a'); tails, heads = {}, {}
    for ax, (name, f) in zip(axes, fields.items()):
        artist = common.collection(ax, g, f['speed'], limits, 'viridis'); decorate(ax, NAMES[name] + ' | selected steady index ' + str(f['selected_iteration']))
        tails[name] = LineCollection([], linewidths=1.3, colors='#ff9ecf', zorder=9); ax.add_collection(tails[name]); heads[name] = ax.scatter([], [], s=12, color='#ff9ecf', edgecolors='#203e4c', linewidths=.25, zorder=10)
    cax = fig.add_axes([.28, .229, .5, .017]); fig.colorbar(artist, cax=cax, orientation='horizontal', label='Speed |U| (m/s) | shared native scale')
    fig.text(.07, .146, 'RK4 particles follow actual native U. Pink identifies inlet seed origin; no concentration or scalar mixing is simulated.', fontsize=8.6, color='#304b5a')
    fig.text(.07, .116, 'One common tracer travel clock; repeated injections. Fixed solved fields, not physical time-resolved CFD evolution.', fontsize=8.5, color='#425e6d')
    fig.text(.07, .088, 'Four native fluid centres are required at every stage. Paths stop at the solid body, outer envelope or stagnation.', fontsize=8.2, color='#425e6d')
    common.add_footer(fig, notes); images = []
    for frame in range(48):
        if frame % 8 == 0:
            common.resource_guard()
        time = TRAVEL_TIME * frame / 47; clock.set_text(f'Tracer travel clock: {time:05.2f} s | 2.40 s playback | fixed steady fields')
        for name in fields:
            segments, dots = [], []
            for path in paths[name]:
                for start in COHORT_STARTS:
                    age = time - start
                    if not 0 <= age <= path['times'][-1]:
                        continue
                    dots.append([np.interp(age, path['times'], path['points'][:, k]) for k in (0, 1)])
                    samples = np.linspace(max(0, age - .5), age, 12)
                    segments.append(np.column_stack([np.interp(samples, path['times'], path['points'][:, k]) for k in (0, 1)]))
            tails[name].set_segments(segments); heads[name].set_offsets(np.array(dots).reshape(-1, 2))
        if frame == 0:
            common.verify_text_bounds(fig, out, 'square-wake-steady-tracers')
        fig.canvas.draw(); images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()))
    gif = out / 'square-wake-steady-tracers.gif'; images[0].save(gif, save_all=True, append_images=images[1:], duration=50, loop=0, optimize=False, disposal=2)
    images[-1].save(out / 'square-wake-steady-tracers.png'); plt.close(fig)
    with Image.open(gif) as decoded:
        require(decoded.n_frames == 48, 'GIF must decode to 48 actual-field tracer frames'); durations, hashes = [], []
        sheet = Image.new('RGB', (images[0].width, (images[0].height + 30) * 3), '#e9eff2'); draw = ImageDraw.Draw(sheet)
        for index in range(48):
            decoded.seek(index); rgb = decoded.convert('RGB'); durations.append(decoded.info.get('duration')); hashes.append(hashlib.sha256(rgb.tobytes()).hexdigest())
            if index in (0, 24, 47):
                top = (0, 24, 47).index(index) * (images[0].height + 30); draw.text((16, top + 8), f'Decoded tracer frame {index + 1}/48; fixed native fields', fill='#203c4b'); sheet.paste(rgb, (0, top + 30))
        require(durations == [50] * 48 and len(set(hashes)) == 48, 'Decoded GIF timing or distinct frames differ'); sheet.save(out / 'square-wake-decoded-GIF-contact-sheet.png')
    return dict(frames=48, frame_delay_ms=50, distinct_decoded_RGB_frames=48, playback_duration_s=2.4, fixed_native_field=True,
                solver_time_evolution=False, transient_vortex_shedding=False, common_tracer_integration_duration_s=TRAVEL_TIME,
                seed_injection_travel_times_s=list(COHORT_STARTS), method='CFL-limited RK4; four native fluid-centre bilinear interpolation; exact segment/body rejection',
                forced_path_closure=False, boundary_values_appended=False, entry_seed_colours_are_not_scalar_mixing=True)


def audited_fields_and_histories(name, data_path, force_path, expected_sha, f, audit):
    data_path, force_path = data_path.resolve(), force_path.resolve()
    if name == 'fluent':
        manifest = json.loads(Path(audit['source_manifest_path']).read_text(encoding='utf-8-sig'))
        field_digest = manifest.get('outputs', {}).get(data_path.name)
        force_digest = manifest.get('outputs', {}).get(force_path.name)
    else:
        record = json.loads(Path(audit['path']).read_text(encoding='utf-8-sig'))
        field_digest, force_digest = record.get('audited_fields_npz_sha256'), record.get('force_history_csv_sha256')
    require(field_digest == expected_sha == sha(data_path) and force_digest == sha(force_path), 'Selected fields or force CSV not hash-bound to native audit: ' + name)
    with np.load(data_path, allow_pickle=False) as z:
        if name == 'fluent':
            U, p = np.column_stack((z['u'], z['v'])), z['p_raw']
        else:
            source = z['native_cell_to_source_cell'].astype(int); require(np.array_equal(np.sort(source), np.arange(15900)), 'Audited OF source mapping not bijective')
            order = np.argsort(source); U, p = z['U'][order, :2], z['p_gauge_Pa'][order]
        require(np.array_equal(f['U'], U) and np.array_equal(f['p'], p), 'Independently reread native fields differ from audited NPZ: ' + name)
        h = {key: z[key].copy() for key in ('force_iteration', 'Cd', 'Cl', 'residual_iteration', 'residual_names', 'residual_values')}
    with force_path.open(encoding='utf-8-sig', newline='') as stream:
        reader = csv.DictReader(stream); require(set(('iteration', 'Cd', 'Cl')) <= set(reader.fieldnames or []), 'Actual force CSV missing iteration/Cd/Cl')
        rows = list(reader)
    iteration = np.array([float(row['iteration']) for row in rows]); cd = np.array([float(row['Cd']) for row in rows]); cl = np.array([float(row['Cl']) for row in rows])
    require(np.array_equal(iteration, np.arange(1, 2001)) and np.isfinite(cd).all() and np.isfinite(cl).all(), 'Require actual complete force iterations 1..2000')
    require(np.array_equal(h['force_iteration'], iteration) and np.array_equal(h['Cd'], cd) and np.array_equal(h['Cl'], cl), 'Force CSV differs from hash-bound NPZ history')
    expected_names = ['continuity', 'x-velocity', 'y-velocity'] if name == 'fluent' else ['p', 'Ux', 'Uy']
    require(h['residual_names'].astype(str).tolist() == expected_names and np.array_equal(h['residual_iteration'], np.arange(1, 2001))
            and h['residual_values'].shape == (2000, 3) and np.isfinite(h['residual_values']).all() and np.all(h['residual_values'] >= 0), 'Incomplete actual native residual history')
    h['force_CSV_path'] = str(force_path); h['force_CSV_sha256'] = force_digest; h['fields_sha256'] = field_digest
    return h


def force_comparison(histories):
    window = {}
    for name, h in histories.items():
        mask = (h['force_iteration'] >= 1801) & (h['force_iteration'] <= 2000)
        require(np.array_equal(h['force_iteration'][mask], np.arange(1801, 2001)), 'Last force window must be exactly 1801..2000')
        cd, cl = h['Cd'][mask], h['Cl'][mask]; mean = float(cd.mean()); mean_cl = float(cl.mean())
        rel = float(np.ptp(cd) / abs(mean)) if mean != 0 else None
        window[name] = dict(first_iteration=1801, last_iteration=2000, samples=200, mean_Cd=mean, mean_Cl=mean_cl,
            Cd_peak_to_peak=float(np.ptp(cd)), Cd_peak_to_peak_over_abs_mean=rel, Cd_mean_nonzero=mean != 0,
            Cd_stability_target=.001, lift_symmetry_target=.001, Cd_stability_pass=rel is not None and rel <= .001,
            lift_symmetry_pass=abs(mean_cl) <= .001, source_force_CSV_path=h['force_CSV_path'], source_force_CSV_sha256=h['force_CSV_sha256'])
    denom = abs(window['fluent']['mean_Cd']); metric = abs(window['fluent']['mean_Cd'] - window['openfoam']['mean_Cd']) / denom if denom != 0 else None
    return dict(window_metrics=window, definition='abs(mean_Cd_F_1801_2000 - mean_Cd_OF_1801_2000) / abs(mean_Cd_F_1801_2000)',
                relative_mean_Cd_difference=metric, target=.01, defined=metric is not None,
                declared_paired_Cd_target_passed=metric is not None and metric <= .01, benchmark_comparison=False, physical_validation=False)


def main():
    global TRAVEL_TIME, COHORT_STARTS
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('out', 'geometry', 'fluent-case', 'fluent-data', 'fluent-audit', 'fluent-fields', 'fluent-force-history',
                 'openfoam-case', 'openfoam-log', 'openfoam-audit', 'openfoam-fields', 'openfoam-force-history'):
        parser.add_argument('--' + flag, type=Path, required=True)
    parser.add_argument('--fluent-fields-sha256', required=True); parser.add_argument('--openfoam-fields-sha256', required=True)
    parser.add_argument('--fluent-iteration', type=int, required=True); parser.add_argument('--openfoam-time', required=True)
    parser.add_argument('--tracer-duration', type=float, required=True, help='Common travel clock in seconds; select only after an actual native-field timing probe')
    parser.add_argument('--seed-injection-period', type=float, required=True, help='Native-field travel-clock spacing between inlet cohorts; no default')
    args = parser.parse_args(); out, script = args.out.resolve(), Path(__file__).resolve()
    require(0 < args.seed_injection_period <= args.tracer_duration <= 60 and args.tracer_duration / args.seed_injection_period <= 12, 'Invalid bounded tracer clock')
    TRAVEL_TIME = args.tracer_duration; COHORT_STARTS = tuple(float(t) for t in np.arange(0., TRAVEL_TIME, args.seed_injection_period))
    require(out.drive.upper() == 'E:' and not out.exists(), 'Fresh E: square-wake export required')
    for key, module in (('elbow', elbow), ('cavity', common), ('readers', common.native)):
        require(sha(Path(module.__file__)) == FROZEN_HELPERS[key], 'Frozen reused helper changed: ' + key)
    resources_before = common.resource_guard(render_reserve=True)
    g = WakeGeometry(args.geometry); setup_path = args.geometry.parent / 'setup.json'; setup = json.loads(setup_path.read_text(encoding='utf-8-sig'))
    expected = dict(density_kg_m3=1., dynamic_viscosity_Pa_s=.0005, inlet_velocity_m_s=[.1, 0.], Re_body=20., body_side_m=.1,
                    body_x_bounds_m=[.45, .55], body_y_bounds_m=[-.05, .05], domain_x_bounds_m=[0., 2.], domain_y_bounds_m=[-.4, .4],
                    outlet_gauge_pressure_Pa=0., operating_pressure_Pa=0., openfoam_extrusion_depth_m=.01, blockage_ratio=.125,
                    confined_domain=True, transient_shedding_simulated=False)
    require(all(setup.get(key) == value for key, value in expected.items()) and setup.get('geometry_npz_sha256') == sha(args.geometry), 'Exact confined-wake setup or geometry hash differs')
    require(all(g.stencil(seed) is not None for _, seed in inlet_seeds() + wake_seeds()), 'Tracer seed lacks four native fluid centres')
    case, data, foam, log = args.fluent_case.resolve(), args.fluent_data.resolve(), args.openfoam_case.resolve(), args.openfoam_log.resolve()
    require(args.fluent_iteration == 2000 and float(args.openfoam_time) == 2000, 'Current declared steady selection is exactly 2000 each')
    folder = foam / args.openfoam_time; mesh = next((p for p in (foam / 'constant/polyMesh', foam / '0/polyMesh') if (p / 'points').is_file()), None)
    require(mesh is not None, 'Actual native OpenFOAM mesh missing')
    required = {'fluent': [case, data, case.parent / 'final-settings.json', case.parent / 'native-physical-readback.json', case.parent / 'solver-transcript.txt'],
                'openfoam': [*(folder / key for key in ('U', 'p', 'C', 'phi')), *(mesh / key for key in ('points', 'faces', 'owner', 'neighbour', 'boundary')),
                             foam / 'system/fvSchemes', foam / 'system/fvSolution', foam / 'system/controlDict', log]}
    statuses = {name: common.audit_status(name, args.fluent_audit if name == 'fluent' else args.openfoam_audit, paths) for name, paths in required.items()}
    data_paths = {'fluent': args.fluent_fields.resolve(), 'openfoam': args.openfoam_fields.resolve()}
    force_paths = {'fluent': args.fluent_force_history.resolve(), 'openfoam': args.openfoam_force_history.resolve()}
    data_sha = {'fluent': args.fluent_fields_sha256, 'openfoam': args.openfoam_fields_sha256}
    inputs = [script, Path(elbow.__file__).resolve(), Path(common.__file__).resolve(), Path(common.native.__file__).resolve(),
              args.geometry.resolve(), setup_path.resolve(), *data_paths.values(), *force_paths.values()]
    for status in statuses.values():
        inputs += [Path(status['path']), Path(status['source_manifest_path']), *(Path(p) for p in status['native_source_hashes'])]
    inputs = sorted(set(inputs), key=str); before = {str(p): sha(p) for p in inputs}
    sources = dict(fluent=common.fluent_native(case, data, rho=1., mu=.0005), openfoam=common.foam_native(foam, folder, log, rho=1.))
    require(sources['fluent']['solve_evidence']['actual_terminal_iteration'] == args.fluent_iteration, 'Actual Fluent terminal iteration differs')
    fields = {name: g.match(value) for name, value in sources.items()}
    for name in fields:
        fields[name]['selected_iteration'] = 2000
    histories = {name: audited_fields_and_histories(name, data_paths[name], force_paths[name], data_sha[name], fields[name], statuses[name]) for name in fields}
    paired = force_comparison(histories); notes = common.footer_lines(statuses)
    metric = paired['relative_mean_Cd_difference']; shown = 'undefined (zero Fluent mean)' if metric is None else f'{100 * metric:.5g}%'
    notes.append('Paired Cd: ' + ('PASS' if paired['declared_paired_Cd_target_passed'] else 'FAIL') + ' | actual final-window mean difference ' + shown + ' | target <= 1%')
    limits = dict(speed=[0., max(float(f['speed'].max()) for f in fields.values())], pressure=[min(float(f['p'].min()) for f in fields.values()), max(float(f['p'].max()) for f in fields.values())])
    require(all(np.isfinite(v).all() and v[1] > v[0] for v in limits.values()), 'Invalid actual native colour limits')
    paths = {name: [g.trace(seed, f['gridU'], origin) for origin, seed in inlet_seeds()] for name, f in fields.items()}
    steady_paths = {name: values + [g.trace(seed, fields[name]['gridU'], origin) for origin, seed in wake_seeds()] for name, values in paths.items()}
    require(all(sum(len(p['times']) >= 8 for p in values) >= 2 for values in paths.values()), 'Too few usable actual inlet trajectories')
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.linewidth': .7}); out.mkdir(parents=True)
    render_fields(g, fields, steady_paths, out, notes, limits); wake = render_wake(g, fields, out, notes, limits); render_diagnostics(histories, out, notes)
    animation = render_animation(g, fields, paths, out, notes, limits['speed'])
    with (out / 'actual-native-cell-fields.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'common_quad_0based', 'x_m', 'y_m', 'u_m_s', 'v_m_s', 'gauge_pressure_Pa', 'speed_m_s'])
        for name, f in fields.items():
            writer.writerows([name, i, *g.centres[i], *f['U'][i], f['p'][i], f['speed'][i]] for i in range(15900))
    with gzip.open(out / 'actual-RK4-tracer-paths.csv.gz', 'wt', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'seed_origin', 'seed_0based', 'travel_time_s', 'x_m', 'y_m', 'speed_m_s'])
        for name, values in steady_paths.items():
            for index, path in enumerate(values):
                writer.writerows([name, path['origin'], index, t, *point, speed] for t, point, speed in zip(path['times'], path['points'], path['speeds']))
    with (out / 'actual-native-residual-history.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'steady_iteration', 'equation', 'own_native_residual'])
        for name, h in histories.items():
            for row, it in enumerate(h['residual_iteration']):
                writer.writerows([name, it, str(key), h['residual_values'][row, col]] for col, key in enumerate(h['residual_names']))
    for name, path in force_paths.items():
        (out / (name + '-force-history.csv')).write_bytes(path.read_bytes())
    (out / 'paired-force-comparison.json').write_text(json.dumps(paired, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    require(before == {str(p): sha(p) for p in inputs}, 'Native source/audit/force history or executed code changed during export')
    checks = dict(model='Steady finite confined square-obstacle wake Re20', actual_native_fields=True, synthetic_flow=False, common_cell_count=15900,
        native_mapping={name: f['mapping'] for name, f in fields.items()}, numerical_audits=statuses, paired_force_comparison=paired,
        all_native_and_paired_declared_checks_passed=all(s['all_declared_checks_passed'] for s in statuses.values()) and paired['declared_paired_Cd_target_passed'],
        input_hashes_stable=True, explicitly_selected_steady_iterations={name: 2000 for name in fields}, common_unclipped_colour_limits=limits,
        source_bound_audited_NPZ={name: dict(path=str(data_paths[name]), sha256=data_sha[name], all15900_native_U_and_p_identical=True) for name in fields},
        animation=animation, wake_diagnostics=wake, gauge_pressure_reference='Native outlet gauge 0 Pa; Fluent static p and rho*OpenFOAM kinematic p, rho=1 kg/m3',
        pressure_mean_subtracted=False, void_cells_masked=True, void_field_values_fabricated=False, four_native_fluid_centres_required_for_interpolation=True,
        exact_segment_body_intersection_rejected=True, finite_confined_domain=True, blockage_ratio=.125, domain_height_over_bodyD=8.,
        body_centre_distance_from_inlet_D=5., body_centre_distance_to_outlet_D=15., top_bottom_boundary='symmetry/slip',
        transient_vortex_shedding_simulated=False, scalar_mixing_simulated=False, three_dimensional_flow_simulated=False,
        published_Cd_assumed=False, physical_validation=False, mesh_independence=False,
        actual_field_ranges={name: {key: [float(f[key].min()), float(f[key].max())] for key in ('speed', 'p')} for name, f in fields.items()},
        native_solve_evidence={name: {key: value for key, value in f['solve_evidence'].items() if key != 'observed_solve_block'} for name, f in fields.items()},
        tracer_diagnostics={name: [dict(seed=p['seed'], origin=p['origin'], samples=len(p['times']), duration_s=float(p['times'][-1]), endpoint_m=p['points'][-1].tolist(),
            integrated_path_length_m=float(np.linalg.norm(np.diff(p['points'], axis=0), axis=1).sum()), stop=p['stop']) for p in values] for name, values in steady_paths.items()},
        text_QA_scope='Implemented bounding-box pair tests are supplemental; final acceptance also requires inspection of actual rendered PNGs and decoded GIF.',
        resource_floor_checks=dict(before=resources_before, after=common.resource_guard()))
    (out / 'checks.json').write_text(json.dumps(checks, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    for module, filename in ((script, 'executed-plotter.py'), (Path(elbow.__file__), 'reused-elbow-helper.py'), (Path(common.__file__), 'reused-cavity-helper.py'), (Path(common.native.__file__), 'reused-native-readers.py')):
        (out / filename).write_bytes(module.read_bytes())
    (out / 'observed-openfoam-solve-block.txt').write_text(sources['openfoam']['solve_evidence']['observed_solve_block'], encoding='utf-8')
    (out / 'input-hash-manifest.json').write_text(json.dumps(dict(command_arguments=sys.argv[1:], inputs=[dict(path=str(p), bytes=p.stat().st_size, sha256=before[str(p)]) for p in inputs]), indent=2) + '\n', encoding='utf-8')
    caption = ('# Actual steady square-obstacle wake\n\n'
        'Finite 2D confined domain x=[0,2] m, y=[-0.4,0.4] m minus the solid square x=[0.45,0.55] m, y=[-0.05,0.05] m. Body-centre distances are 5 D from inlet and 15 D to outlet. Height 8 D gives 12.5% blockage. Top and bottom are symmetry/slip; the body is stationary no-slip. Inlet U=0.1 m/s, rho=1 kg/m3, mu=0.0005 Pa s, D=0.1 m, Re=20. Native fluid fields map through complete vertex sets to 15,900 source quads; the actual body hole remains unfilled.\n\n'
        'Pressure retains the common outlet gauge reference 0 Pa. No mean subtraction or invented void values. Both codes are steady; RK4 particles only follow the fixed actual saved velocity. All four interpolation centres must be native fluid cells and every stage/segment is rejected if it touches/crosses the body. Inlet colours identify seed origin, not concentration. No forced path closure or transient vortex-shedding claim.\n\n'
        'Force histories contain actual native iterations 1..2000. Stability/symmetry windows are exactly 1801..2000. The separately declared paired Cd target is relative final-window mean difference <=1%, normalized by abs(Fluent mean Cd); a zero mean is undefined and false. Native numerical flags and the paired flag remain separate and every miss is shown. No published Cd, experimental validation, mesh independence or unbounded-flow accuracy claim.\n\n'
        + f'Tracer travel clock {TRAVEL_TIME} s, injection period {args.seed_injection_period} s, replay 48 frames at 50 ms. These parameters must be selected by an actual native-field timing probe.\n\n'
        + 'The user directs the project; Codex prepares and operates code; Fluent/OpenFOAM compute the fields. Independent user execution/reflection is not invented.\n\n' + '\n'.join(notes) + '\n')
    (out / 'CAPTIONS.md').write_text(caption, encoding='utf-8')
    outputs = sorted((p for p in out.iterdir() if p.is_file()), key=str)
    (out / 'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name, bytes=p.stat().st_size, sha256=sha(p)) for p in outputs], indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(output=str(out), animation=animation, native_numerical_failures={name: s['failed_checks'] for name, s in statuses.items()}, paired_force_comparison=paired), indent=2))


if __name__ == '__main__':
    main()
