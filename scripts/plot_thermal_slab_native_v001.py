"""Actual transient solid temperatures: 51 native states, no synthetic time.

Every one of the 101 saved temperature states is reread from native files and
compared against the independently audited NPZ before any picture is rendered.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, importlib.util, json, os, re, sys
os.environ['TEMP'] = os.environ['TMP'] = 'E:/CFD-Learning-Lab/temp'
os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-thermal-slab'
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw
import plot_cavity_native_v003 as common
import analyse_nozzle_fluent_v005 as mapping_helper
import thermal_slab_reference_v001 as reference
os.environ['MPLCONFIGDIR'] = 'E:/CFD-Learning-Lab/temp/matplotlib-thermal-slab'
require, sha, NAMES = common.require, common.sha, common.NAMES
COLOURS = {'fluent': '#087f8c', 'openfoam': '#d46725'}
WALLS = ('hot', 'cold', 'insulated_bottom', 'insulated_top')
TIMES = np.arange(101, dtype=float) * 2.
FRAME_INDICES = np.arange(0, 101, 2)
LISPY_PATH = Path('E:/CFD-Learning-Lab/environments/fluent/Lib/site-packages/ansys/fluent/core/offline/_lispy.py')
REFERENCE_SHA = '550c050d4294c08d0fb3905d06a2eeda7c82bb44397072a3dd261db1a6d9846a'


class SlabGeometry:
    def __init__(self, path):
        with np.load(path, allow_pickle=False) as source:
            self.source = {key: source[key].copy() for key in source.files}
        for key in ('nodes', 'quads', 'centres', 'areas', 'ix', 'iy', 'x_edges', 'y_edges'):
            setattr(self, key, self.source[key])
        self.quads = self.quads.astype(np.int64); self.ix = self.ix.astype(int); self.iy = self.iy.astype(int)
        self.ny, self.nx = map(int, self.source['grid_shape'])
        require((self.ny, self.nx) == (20, 100) and self.nodes.shape == (2121, 2)
                and self.quads.shape == (2000, 4) and self.centres.shape == (2000, 2), 'Expected 100x20 solid slab quads')
        require(np.allclose(self.x_edges, np.linspace(0., .1, 101), rtol=0, atol=1e-12)
                and np.allclose(self.y_edges, np.linspace(0., .02, 21), rtol=0, atol=1e-12), 'Unexpected source slab extent')
        self.polygons = self.nodes[self.quads]
        x, y = self.polygons[:, :, 0], self.polygons[:, :, 1]
        cross = x * np.roll(y, -1, axis=1) - np.roll(x, -1, axis=1) * y
        areas = cross.sum(axis=1) / 2
        centres = np.column_stack((((x + np.roll(x, -1, axis=1)) * cross).sum(axis=1),
                                   ((y + np.roll(y, -1, axis=1)) * cross).sum(axis=1))) / (6 * areas[:, None])
        require(np.isfinite(self.nodes).all() and np.all(areas > 0) and np.max(np.abs(areas - self.areas)) < 1e-14
                and abs(areas.sum() - .002) < 1e-13 and np.max(np.abs(centres - self.centres)) < 1e-10, 'Source polygon area or true centroid differs')
        slots = self.iy * self.nx + self.ix
        require(np.array_equal(np.sort(slots), np.arange(2000)), 'Source slab ordering not bijective')
        self.order = np.argsort(slots); self.xc = (self.x_edges[:-1] + self.x_edges[1:]) / 2
        require(np.max(np.abs(self.centres - np.column_stack((self.xc[self.ix], (self.y_edges[self.iy] + self.y_edges[self.iy + 1]) / 2)) )) < 1e-10, 'Source centres do not match regular slab')
        self.node_lookup = {tuple(np.round(point, 10)): i for i, point in enumerate(self.nodes)}
        self.cell_lookup = {tuple(sorted(map(int, quad))): i for i, quad in enumerate(self.quads)}

    def profile(self, values):
        array = np.asarray(values)
        grid = array[..., self.order].reshape(*array.shape[:-1], self.ny, self.nx)
        weights = self.areas[self.order].reshape(self.ny, self.nx)
        return np.sum(grid * weights, axis=-2) / weights.sum(axis=0)

    def foam_mapping(self, mesh, centres):
        node_ids = []
        for point in mesh['nodes'][:, :2]:
            key = tuple(np.round(point, 10)); require(key in self.node_lookup, 'Native solid vertex absent from source slab')
            node_ids.append(self.node_lookup[key])
        node_ids = np.asarray(node_ids, np.int64); mapping = []
        for vertices in mesh['cells']:
            key = tuple(np.unique(node_ids[vertices])); require(len(key) == 4 and key in self.cell_lookup, 'Native solid hex does not project to one source quad')
            mapping.append(self.cell_lookup[key])
        mapping = np.asarray(mapping, np.int64)
        require(np.array_equal(np.sort(mapping), np.arange(2000)), 'Native solid vertex-set map not bijective')
        error = float(np.max(np.linalg.norm(centres[:, :2] - self.centres[mapping], axis=1)))
        require(error < 1e-10 and np.max(np.abs(centres[:, 2])) < 1e-10, 'Native C differs from source slab centroids')
        return mapping, dict(vertex_set_bijection=True, source_quad_cells=2000, maximum_saved_C_error_m=error)


def fluent_time_variables(handle, parser):
    raw = handle['settings/Data Variables'][0]
    raw = raw.decode() if isinstance(raw, bytes) else str(raw)
    return {str(item[0]): item[1] for item in parser.parse(raw)[1]}


def native_temperatures(name, run, g, parser=None):
    if name == 'fluent':
        case = run / 'slab.cas.h5'
        mapping = mapping_helper.common_mapping(case, g.source)
        values, native_times, native_steps, native_dt = [], [], [], []
        for index in range(101):
            path = run / 'snapshots' / f'step-{index:06d}.dat.h5'
            with h5py.File(path, 'r') as handle:
                field = common.native.indexed(handle['results/1/phase-1/cells/SV_T'], 2000)
                variables = fluent_time_variables(handle, parser)
                time, raw_step, dt = float(variables['flow-time']), float(variables['time-step']), float(variables['physical-time-step'])
                require(time == TIMES[index] and raw_step == index and dt == 2., 'Actual Fluent time/step/dt differs from selected snapshot')
                step = int(raw_step)
                row = np.empty(2000); row[mapping] = field; values.append(row)
                native_times.append(time); native_steps.append(step); native_dt.append(dt)
        return np.array(values), dict(vertex_set_bijection=True, source_quad_cells=2000), dict(actual_times_s=native_times, actual_steps=native_steps, actual_dt_s=native_dt), g.areas.copy()
    mesh = common.native.foam_mesh(run)
    require(np.allclose(np.unique(mesh['nodes'][:, 2]), [-.005, .005], atol=1e-10, rtol=0), 'Expected native .01m solid extrusion')
    centres = common.native.foam_field(run / '0/C', 2000, 3, [0, 1, 0, 0, 0, 0, 0])
    mapping, proof = g.foam_mapping(mesh, centres)
    volumes = common.native.foam_field(run / '0/Vc', 2000, 1, [0, 3, 0, 0, 0, 0, 0])
    require(np.max(np.abs(volumes - g.areas[mapping] * .01)) < 1e-15, 'Native Vc differs from mapped source volume')
    values = []
    for time in TIMES:
        field = common.native.foam_field(run / str(int(time)) / 'T', 2000, 1, [0, 0, 0, 1, 0, 0, 0])
        row = np.empty(2000); row[mapping] = field; values.append(row)
    log = (run / 'solver.log').read_text(encoding='utf-8', errors='replace')
    actual_times = [float(value) for value in re.findall(r'^Time\s*=\s*([-+0-9.eE]+)s?\s*$', log, re.M)]
    require(actual_times == TIMES[1:].tolist() and re.search(r'Solving for (?:e|h|T),', log), 'Actual solid energy time blocks differ')
    return np.array(values), proof, dict(actual_times_s=TIMES.tolist(), actual_saved_folders=[str(int(v)) for v in TIMES], actual_solved_log_times_s=actual_times), g.areas * .01


def audited_arrays(name, path, expected_hash, T, g, volumes, refs):
    require(sha(path) == expected_hash, 'Selected audited NPZ hash differs: ' + name)
    with np.load(path, allow_pickle=False) as source:
        require(np.array_equal(source['times_s'], TIMES), 'Audited physical times differ')
        require(tuple(map(str, source['wall_names'])) == WALLS, 'Audited heat-rate column order differs')
        if name == 'fluent':
            audited = source['temperature_K']; Q = source['wall_heat_in_W'][1:].copy(); residual = source['energy_residual'].copy()
            require(np.array_equal(source['native_step'], np.arange(101)) and np.all(source['native_dt_s'] == 2.), 'Audited native Fluent step/dt differs')
            require(np.array_equal(source['wall_heat_solved'], np.arange(101) > 0), 'Initial Fluent heat rates must be explicitly unsolved')
        else:
            mapping = source['native_cell_to_source_cell'].astype(int)
            require(np.array_equal(np.sort(mapping), np.arange(2000)), 'Audited OpenFOAM map not bijective')
            audited = np.empty((101, 2000)); audited[:, mapping] = source['T']
            Q = source['native_Q_inward_W'].copy(); residual = source['final_linear_energy_residual'].copy()
            require(np.array_equal(source['heat_times_s'], TIMES[1:]) and np.array_equal(source['residual_times_s'], TIMES[1:]), 'Native solved heat/residual times differ')
        require(audited.shape == T.shape == (101, 2000) and np.isfinite(T).all()
                and np.array_equal(audited, T), 'Direct native T differs from independent audited NPZ')
        require(Q.shape == (100, 4) and residual.shape == (100,) and np.isfinite(Q).all()
                and np.isfinite(residual).all() and np.all(residual >= 0), 'Invalid audited solved-step heat rates/residuals')
        energy = 1e6 * np.sum((T - 300.) * volumes, axis=1)
        require(np.allclose(energy, source['energy_J'], rtol=1e-11, atol=2e-7), 'Audited energy differs from all native temperatures')
        power = np.diff(energy) / 2.; signed_balance = power - Q.sum(axis=1)
        depth = 1. if name == 'fluent' else .01
        scale = np.maximum.reduce([np.sum(np.abs(Q), axis=1), np.abs(power), np.full(100, 200. * depth)])
        relative = np.abs(signed_balance) / scale
        require(np.allclose(relative, source['energy_balance_relative'], rtol=2e-7, atol=5e-11), 'Audited global energy balance differs from native E and Qin')
        error = T - refs
        max_error = np.max(np.abs(error), axis=1); rms_error = np.sqrt(np.sum(error**2 * volumes, axis=1) / volumes.sum())
        require(np.allclose(max_error, source['analytic_max_error_K'], rtol=1e-8, atol=1e-9)
                and np.allclose(rms_error, source['analytic_rms_error_K' if name == 'fluent' else 'analytic_RMS_error_K'], rtol=1e-8, atol=1e-9), 'Audited analytical error metrics differ from direct native T')
    return dict(T=T, profile=g.profile(T), energy=energy, Q=Q, residual=residual, energy_balance_relative=relative,
                signed_energy_balance_W=signed_balance, depth=depth, error=error,
                max_error=max_error, rms_error=rms_error,
                audited_NPZ_sha256=expected_hash, all_native_temperature_values_match_audited_NPZ=True)


def decorate_slab(ax, g, name):
    ax.plot([0, .1, .1, 0, 0], [0, 0, .02, .02, 0], color='#304a58', lw=.8)
    ax.plot([0, 0], [0, .02], color='#d35c35', lw=2.)
    ax.plot([.1, .1], [0, .02], color='#2a86ad', lw=2.)
    ax.set(xlim=(0, .1), ylim=(0, .02), xlabel='x (m)', ylabel='y (m)')
    ax.set_aspect('equal'); ax.set_title(NAMES[name] + ' | actual solid temperature', loc='left', fontsize=11, pad=9, fontweight='bold')
    ax.set_yticks([0, .01, .02]); ax.tick_params(labelsize=9)


def render_transient(g, fields, refs, out, status_lines):
    fig, axes = plt.subplots(3, 1, figsize=(12, 10.7 + .16 * len(status_lines)), gridspec_kw={'height_ratios': [1, 1, 1.35]})
    fig.subplots_adjust(left=.10, right=.95, top=.82, bottom=.27 + .018 * len(status_lines), hspace=.65)
    fig.suptitle('Heat diffuses through a solid slab', x=.07, y=.962, ha='left', fontsize=21, fontweight='bold')
    fig.text(.07, .918, 'Actual transient Fluent + OpenFOAM temperatures | hot end 400K, cold end 300K | insulated sides', fontsize=10, color='#425e6d')
    clock = fig.text(.07, .872, '', fontsize=11, color='#304b59')
    images, artists, curves = [], {}, {}
    for ax, (name, field) in zip(axes[:2], fields.items()):
        artists[name] = common.collection(ax, g, field['T'][0], [300., 400.], 'inferno')
        decorate_slab(ax, g, name)
    profile_ax = axes[2]
    for name in fields:
        curves[name], = profile_ax.plot(g.xc, fields[name]['profile'][0], color=COLOURS[name], lw=1.9,
                                        ls='--' if name == 'openfoam' else '-', label=NAMES[name] + ': native column averages')
    exact, = profile_ax.plot(g.xc, g.profile(refs)[0], color='#303d48', ls=':', lw=1.6, label='Continuous PDE reference: exact cell averages')
    profile_ax.set(xlim=(0, .1), ylim=(298, 402), xlabel='x (m)', ylabel='Temperature (K)')
    profile_ax.grid(alpha=.18); profile_ax.spines[['top', 'right']].set_visible(False)
    profile_ax.legend(loc='upper right', frameon=False, fontsize=8.4, ncol=1)
    cax = fig.add_axes([.28, .225, .5, .016]); fig.colorbar(artists['fluent'], cax=cax, orientation='horizontal', label='Temperature (K) | fixed common 300–400K scale')
    fig.text(.07, .142, 'Matched implicit Euler, Δt = 2s. Every frame uses an actual saved physical time; no temporal interpolation.', fontsize=9, color='#304b59')
    fig.text(.07, .112, 'Continuous reference uses exact cell averages. Early coarse-step errors remain in the separate error figure.', fontsize=8.5, color='#425e6d')
    fig.text(.07, .085, 'Numerical verification of an ideal constant-property conduction model; no experimental physical-validation claim.', fontsize=8.4, color='#425e6d')
    common.add_footer(fig, status_lines)
    ref_profiles = g.profile(refs)
    for frame, state in enumerate(FRAME_INDICES):
        if frame % 8 == 0:
            common.resource_guard()
        clock.set_text(f'Actual physical time: {TIMES[state]:.0f}s | native step {state}/100 | playback 50ms per saved frame')
        for name in fields:
            artists[name].set_array(fields[name]['T'][state]); curves[name].set_ydata(fields[name]['profile'][state])
        exact.set_ydata(ref_profiles[state])
        if frame == 0:
            common.verify_text_bounds(fig, out, 'thermal-slab-native-temperature')
        fig.canvas.draw(); images.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()))
    gif = out / 'thermal-slab-actual-transient.gif'
    images[0].save(gif, save_all=True, append_images=images[1:], duration=50, loop=0, optimize=False, disposal=2)
    images[-1].save(out / 'thermal-slab-native-temperature.png'); plt.close(fig)
    with Image.open(gif) as decoded:
        require(decoded.n_frames == 51, 'Thermal GIF requires exactly 51 actual saved states')
        durations, hashes = [], []
        contact = Image.new('RGB', (images[0].width, (images[0].height + 30) * 3), '#e9eff2'); draw = ImageDraw.Draw(contact)
        for frame in range(51):
            decoded.seek(frame); rgb = decoded.convert('RGB'); durations.append(decoded.info.get('duration'))
            hashes.append(hashlib.sha256(rgb.tobytes()).hexdigest())
            if frame in (0, 25, 50):
                top = (0, 25, 50).index(frame) * (images[0].height + 30)
                draw.text((16, top + 8), f'Decoded actual native state: t={TIMES[FRAME_INDICES[frame]]:.0f}s | frame {frame + 1}/51', fill='#203c4b')
                contact.paste(rgb, (0, top + 30))
        require(durations == [50] * 51 and len(set(hashes)) == 51, 'Thermal GIF decoded timing or distinct states differ')
        contact.save(out / 'thermal-slab-decoded-GIF-contact-sheet.png')
    return dict(frames=51, frame_delay_ms=50, distinct_decoded_RGB_frames=51, playback_duration_s=2.55,
                actual_physical_times_s=TIMES[FRAME_INDICES].tolist(), native_state_indices=FRAME_INDICES.tolist(),
                actual_solver_time_evolution=True, temporal_interpolation=False, synthetic_temperature=False,
                tracer_animation=False, shared_fixed_temperature_limits_K=[300., 400.])


def render_errors(g, fields, out, status_lines):
    fig = plt.figure(figsize=(12, 9 + .16 * len(status_lines))); gs = fig.add_gridspec(2, 2)
    axes = [fig.add_subplot(gs[0, :]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    fig.subplots_adjust(left=.09, right=.96, top=.83, bottom=.18 + .018 * len(status_lines), hspace=.55, wspace=.28)
    fig.suptitle('Temperature error against the continuous slab solution', x=.09, y=.965, ha='left', fontsize=17, fontweight='bold')
    fig.text(.09, .91, 'Exact cell averages of the continuous PDE; all101 saved physical times retained, including the initial boundary jump.', fontsize=9, color='#425e6d')
    for name, field in fields.items():
        style = '--' if name == 'openfoam' else '-'
        axes[0].plot(TIMES, field['max_error'], color=COLOURS[name], ls=style, lw=1.7, label=NAMES[name] + ' max |error|')
        axes[0].plot(TIMES, field['rms_error'], color=COLOURS[name], ls='-.' if name == 'openfoam' else ':', lw=1.5, label=NAMES[name] + ' weighted RMS')
        for ax, index in zip(axes[1:], (1, 100)):
            ax.plot(g.xc, g.profile(field['error'])[index], color=COLOURS[name], ls=style, lw=1.8, label=NAMES[name])
    axes[0].scatter([200, 200], [.5, .25], color='#586b76', marker='x', s=32, label='Final-only targets: max0.5K / RMS0.25K')
    axes[0].set(xlabel='Actual physical time (s)', ylabel='Temperature error (K)', xlim=(0, 205))
    axes[0].legend(frameon=False, fontsize=8.4, ncol=3)
    for ax, title in zip(axes[1:], ('First solved state, t=2s | early Euler error', 'Final state, t=200s | accuracy targets apply here')):
        ax.set(title=title, xlabel='x (m)', ylabel='Native − exact cell average (K)', xlim=(0, .1)); ax.axhline(0, color='#80929b', lw=.7)
        ax.legend(frameon=False, fontsize=9)
    for ax in axes:
        ax.grid(alpha=.17); ax.spines[['top', 'right']].set_visible(False)
    fig.text(.09, .102 + .019 * len(status_lines), 'Final targets apply only at200s. Earlier implicit-Euler errors are reported without an all-time accuracy threshold.', fontsize=8.6, color='#425e6d')
    common.add_footer(fig, status_lines); common.save(fig, out, 'thermal-slab-temperature-error')


def render_energy(fields, refs, g, out, status_lines):
    fig, axes = plt.subplots(2, 2, figsize=(12, 9 + .16 * len(status_lines)))
    fig.subplots_adjust(left=.09, right=.96, top=.82, bottom=.18 + .018 * len(status_lines), hspace=.53, wspace=.29)
    fig.suptitle('Global thermal energy and actual boundary heat rates', x=.09, y=.965, ha='left', fontsize=17, fontweight='bold')
    fig.text(.09, .91, 'E and wall power are divided by native reference/extrusion depth for comparison; raw native W and J remain in the CSV.', fontsize=9, color='#425e6d')
    analytic_energy = 1e6 * np.sum((refs - 300.) * g.areas, axis=1)
    axes[0, 0].plot(TIMES, analytic_energy, color='#303d48', ls=':', lw=1.5, label='Continuous PDE reference')
    for name, field in fields.items():
        style = '--' if name == 'openfoam' else '-'; depth = field['depth']
        axes[0, 0].plot(TIMES, field['energy'] / depth, color=COLOURS[name], ls=style, lw=1.8, label=NAMES[name])
        for column, colour in ((0, '#cc493d'), (1, '#2c81ad')):
            axes[0, 1].plot(TIMES[1:], field['Q'][:, column] / depth, color=colour, ls=style, lw=1.5, label=NAMES[name] + ' ' + WALLS[column])
        for column, colour in ((2, '#337e70'), (3, '#9d639d')):
            axes[1, 0].plot(TIMES[1:], field['Q'][:, column] / depth, color=colour, ls=style, lw=1.5, label=NAMES[name] + ' ' + WALLS[column].replace('insulated_', ''))
        axes[1, 1].plot(TIMES[1:], field['energy_balance_relative'], color=COLOURS[name], ls=style, lw=1.5, label=NAMES[name])
    axes[0, 0].set(xlabel='Physical time (s)', ylabel='Energy above300K / depth (J/m)')
    axes[0, 1].set(xlabel='Solved physical time (s)', ylabel='Native heat rate into solid / depth (W/m)')
    axes[1, 0].set(xlabel='Solved physical time (s)', ylabel='Native insulated-wall power / depth (W/m)')
    axes[1, 1].set(xlabel='Solved physical time (s)', ylabel='Absolute global energy imbalance / scale')
    axes[1, 1].set_yscale('symlog', linthresh=1e-14); axes[1, 1].axhline(1e-5, color='#586b76', lw=.9, ls=':', label='Declared per-step target1e-5')
    for ax in axes.ravel():
        ax.grid(alpha=.17); ax.spines[['top', 'right']].set_visible(False); ax.legend(frameon=False, fontsize=8.1)
    fig.text(.09, .102 + .019 * len(status_lines), 'Backward-Euler balance: (Eₙ−Eₙ₋₁)/2s − ΣQin,n. Global only; imposed insulation alone is not conservation proof.', fontsize=8.5, color='#425e6d')
    common.add_footer(fig, status_lines); common.save(fig, out, 'thermal-slab-energy-balance')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('geometry', 'out', 'fluent-run', 'fluent-audit', 'fluent-fields', 'openfoam-run', 'openfoam-audit', 'openfoam-fields'):
        parser.add_argument('--' + flag, type=Path, required=True)
    parser.add_argument('--fluent-fields-sha256', required=True); parser.add_argument('--openfoam-fields-sha256', required=True)
    args = parser.parse_args(); out, script = args.out.resolve(), Path(__file__).resolve()
    require(out.drive.upper() == 'E:' and not out.exists(), 'Fresh E: thermal export required')
    g = SlabGeometry(args.geometry); setup_path = args.geometry.parent / 'setup.json'; setup = json.loads(setup_path.read_text(encoding='utf-8-sig'))
    expected = dict(density_kg_m3=1000., specific_heat_J_kg_K=1000., thermal_conductivity_W_m_K=10., thermal_diffusivity_m2_s=1e-5,
                    initial_temperature_K=300., hot_temperature_K=400., cold_temperature_K=300., time_step_s=2., time_steps=100,
                    final_time_s=200., reference_depth_m=1., openfoam_extrusion_depth_m=.01, solved_zone_type='solid')
    require(all(setup.get(key) == value for key, value in expected.items()) and setup.get('geometry_npz_sha256') == sha(args.geometry), 'Exact thermal setup or source geometry hash differs')
    require(sha(Path(reference.__file__)) == REFERENCE_SHA, 'Frozen analytical reference changed')
    ref_profiles = np.array([reference.cell_average_temperature(g.x_edges[:-1], g.x_edges[1:], time, terms=80) for time in TIMES])
    ref160 = np.array([reference.cell_average_temperature(g.x_edges[:-1], g.x_edges[1:], time, terms=160) for time in TIMES])
    reference_series_error = float(np.max(np.abs(ref_profiles - ref160))); require(reference_series_error < 1e-9, '80/160-term analytical cell averages differ')
    refs = ref_profiles[:, g.ix]
    runs = {'fluent': args.fluent_run.resolve(), 'openfoam': args.openfoam_run.resolve()}
    audits = {'fluent': args.fluent_audit.resolve(), 'openfoam': args.openfoam_audit.resolve()}
    data_paths = {'fluent': args.fluent_fields.resolve(), 'openfoam': args.openfoam_fields.resolve()}
    data_hashes = {'fluent': args.fluent_fields_sha256, 'openfoam': args.openfoam_fields_sha256}
    mesh = next((p for p in (runs['openfoam'] / 'constant/polyMesh', runs['openfoam'] / '0/polyMesh') if (p / 'points').is_file()), None)
    require(mesh is not None, 'Native solid polyMesh missing')
    required = {'fluent': [runs['fluent'] / 'slab.cas.h5', *(runs['fluent'] / 'snapshots' / f'step-{i:06d}.dat.h5' for i in range(101)),
                           runs['fluent'] / 'slab-final.cas.h5', runs['fluent'] / 'slab-final.dat.h5',
                           *(runs['fluent'] / name for name in ('loaded-settings.json', 'initialized-settings.json', 'first-step-settings.json', 'final-settings.json')),
                           runs['fluent'] / 'solver-transcript.txt'],
                'openfoam': [*(runs['openfoam'] / str(int(t)) / 'T' for t in TIMES), *(mesh / name for name in ('points', 'faces', 'owner', 'neighbour', 'boundary')),
                             *(runs['openfoam'] / time / field for time in ('0', '200') for field in ('C', 'Vc')),
                             *(runs['openfoam'] / 'system' / name for name in ('fvSchemes', 'fvSolution', 'controlDict')), runs['openfoam'] / 'solver.log']}
    statuses = {name: common.audit_status(name, audits[name], paths) for name, paths in required.items()}
    for name, status in statuses.items():
        if name == 'fluent':
            manifest = json.loads(Path(status['source_manifest_path']).read_text(encoding='utf-8-sig'))
            audited_hash = manifest.get('outputs', {}).get(data_paths[name].name)
        else:
            audited_hash = json.loads(audits[name].read_text(encoding='utf-8-sig')).get('audited_fields_npz_sha256')
        require(audited_hash == data_hashes[name] == sha(data_paths[name]), 'Selected NPZ is not hash-bound to the native audit output: ' + name)
    inputs = [script, args.geometry.resolve(), setup_path.resolve(), Path(common.__file__).resolve(), Path(common.native.__file__).resolve(),
              Path(mapping_helper.__file__).resolve(), Path(reference.__file__).resolve(), LISPY_PATH, *data_paths.values()]
    for status in statuses.values():
        inputs += [Path(status['path']), Path(status['source_manifest_path']), *(Path(p) for p in status['native_source_hashes'])]
    inputs = sorted(set(inputs), key=str); before = {str(path): sha(path) for path in inputs}
    resources_before = common.resource_guard(render_reserve=True)
    spec = importlib.util.spec_from_file_location('thermal_official_lispy', LISPY_PATH); lispy = importlib.util.module_from_spec(spec); spec.loader.exec_module(lispy)
    fields, mappings, time_evidence = {}, {}, {}
    for name in runs:
        T, mappings[name], time_evidence[name], volumes = native_temperatures(name, runs[name], g, parser=lispy)
        fields[name] = audited_arrays(name, data_paths[name], data_hashes[name], T, g, volumes, refs)
    status_lines = common.footer_lines(statuses)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.linewidth': .7}); out.mkdir(parents=True)
    animation = render_transient(g, fields, refs, out, status_lines)
    render_errors(g, fields, out, status_lines); render_energy(fields, refs, g, out, status_lines)
    with gzip.open(out / 'actual-native-temperatures.csv.gz', 'wt', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'actual_physical_time_s', 'source_quad_0based', 'x_m', 'y_m', 'native_temperature_K', 'continuous_PDE_exact_cell_average_K', 'native_minus_reference_K'])
        for name, field in fields.items():
            for time_index, time in enumerate(TIMES):
                writer.writerows([name, time, i, *g.centres[i], field['T'][time_index, i], refs[time_index, i], field['error'][time_index, i]] for i in range(2000))
    with (out / 'actual-time-metrics.csv').open('x', encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream); writer.writerow(['solver', 'actual_time_s', 'energy_J', 'depth_m', 'max_abs_error_K', 'weighted_RMS_error_K', 'Qhot_in_W', 'Qcold_in_W', 'Qbottom_in_W', 'Qtop_in_W', 'global_energy_balance_relative', 'own_native_thermal_residual', 'initial_unsolved_state'])
        for name, field in fields.items():
            for index, time in enumerate(TIMES):
                solved = [*field['Q'][index - 1], field['energy_balance_relative'][index - 1], field['residual'][index - 1]] if index else [''] * 6
                writer.writerow([name, time, field['energy'][index], field['depth'], field['max_error'][index], field['rms_error'][index], *solved, index == 0])
    require(before == {str(path): sha(path) for path in inputs}, 'Native thermal sources, audits or executed code changed during export')
    checks = dict(model='Transient constant-property solid slab conduction', actual_native_fields=True, synthetic_temperature=False, common_cell_count=2000,
                  actual_saved_states_per_solver=101, native_temperature_values_verified_per_solver=202000, paired_native_temperature_values_verified=404000,
                  actual_saved_physical_times_s=TIMES.tolist(), native_mapping=mappings, actual_native_time_evidence=time_evidence,
                  paired_heat_and_residual_times_s=TIMES[1:].tolist(), initial_unsolved_Fluent_heat_rates_excluded=True,
                  numerical_audits=statuses, input_hashes_stable=True, animation=animation,
                  source_bound_audited_NPZ={name: dict(path=str(data_paths[name]), sha256=data_hashes[name], all202000_native_T_values_identical=True) for name in fields},
                  actual_temperature_ranges_K={name: [float(f['T'].min()), float(f['T'].max())] for name, f in fields.items()},
                  final_errors_K={name: dict(max_abs=float(f['max_error'][-1]), weighted_RMS=float(f['rms_error'][-1])) for name, f in fields.items()},
                  all_time_errors={name: dict(max_abs_K=f['max_error'].tolist(), weighted_RMS_K=f['rms_error'].tolist(), worst_max_abs_time_s=float(TIMES[int(np.argmax(f['max_error']))])) for name, f in fields.items()},
                  final_only_accuracy_targets_K=dict(max_abs=.5, weighted_RMS=.25), all_time_accuracy_threshold_imposed=False,
                  analytical_reference=dict(kind='Continuous Dirichlet slab PDE, exact cell averages', source_sha256=REFERENCE_SHA, terms=80, max_80_vs_160_difference_K=reference_series_error,
                      initial_interior_explicit_300K=True, reference_temperatures_clipped=False),
                  global_energy_balance_only=True, native_all_cell_thermal_flux_conservation_claim=False,
                  physical_validation=False, mesh_independence=False, time_step_independence=False, fluid_motion=False,
                  resource_floor_checks=dict(before=resources_before, after=common.resource_guard()))
    (out / 'checks.json').write_text(json.dumps(checks, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    for source, name in ((script, 'executed-plotter.py'), (Path(common.__file__), 'reused-cavity-helper.py'), (Path(common.native.__file__), 'reused-native-readers.py'),
                         (Path(mapping_helper.__file__), 'reused-native-mapping-helper.py'), (Path(reference.__file__), 'analytical-reference.py')):
        (out / name).write_bytes(source.read_bytes())
    (out / 'input-hash-manifest.json').write_text(json.dumps(dict(command_arguments=sys.argv[1:], inputs=[dict(path=str(p), bytes=p.stat().st_size, sha256=before[str(p)]) for p in inputs]), indent=2) + '\n', encoding='utf-8')
    caption = ('# Actual transient solid-slab temperatures\n\n'
               'Generic constant-property model solid: rho1000kg/m3, c1000J/(kgK), k10W/(mK), alpha1e-5m2/s. Rectangle.1m by.02m, common100x20quads. Initial interior300K; hot left jumps to400K, cold right300K, insulated top/bottom. Both codes use matched first-order implicit Euler,dt2s,100steps to200s. No fluid motion, tracers or synthetic physical times.\n\n'
               'All101native temperature states are independently reread and matched to auditedNPZ,202000cell values per solver. Native cell connectivity maps by complete vertex sets. Fluent physical time, timestep index and dt are independently parsed from saved native Data Variables; OpenFOAM uses actual native time folders and solved energy-log times. GIF51frames select exactly0,4,...200s,50ms/frame,2.55s playback. Temperature colour limits stay300..400K and the rectangle retains its physical5:1aspect.\n\n'
               'Dotted reference profiles are exact cell averages of the continuous1D Dirichlet PDE, not another numerical solver result. The initial interior is explicitly300K;80vs160series terms agree within1e-9K for allsavedtimes. Reference/data values are not clipped. Errors at every actual saved time, including early coarse-implicit-Euler boundary-jump errors, remain in the CSV and separate error figure. Only final200s has accuracy targets:maxabs0.5K, volume-weightedRMS0.25K. No all-time accuracy threshold is imposed.\n\n'
               'E=sum(rho*c*V*(T-300K)); nativeQ is signed positive into the solid. Global backward-Euler balance uses(E_n-E_prev)/2s minus summed nativeQ_n. Paired heat-rate/residual figures use only solved times2..200s. Fluent initialQ is explicitly unsolved; actual OpenFOAM time0boundary powers remain in its native audit and are omitted from paired plots. No time0Euler balance is invented. Energy/power plots divide by native depth1m(Fluent)or.01m(OpenFOAM); rawJ/W remain in metrics. Insulation is reported separately and is not sole conservation proof. No all-cell thermal flux conservation is claimed.\n\n'
               'This is numerical verification of an ideal material/PDE model, not experimental physical validation or timestep/mesh independence. The user directs the project; Codex selected/prepared the model and workflow; solvers computed fields. No user reflection or independent user execution is invented.\n\n' + '\n'.join(status_lines) + '\n')
    (out / 'CAPTIONS.md').write_text(caption, encoding='utf-8')
    files = sorted((p for p in out.iterdir() if p.is_file()), key=str)
    (out / 'output-hash-manifest.json').write_text(json.dumps([dict(path=p.name, bytes=p.stat().st_size, sha256=sha(p)) for p in files], indent=2) + '\n', encoding='utf-8')
    print(json.dumps(dict(output=str(out), animation=animation, final_errors_K=checks['final_errors_K'], numerical_failures={k: v['failed_checks'] for k, v in statuses.items()}), indent=2))


if __name__ == '__main__':
    main()
