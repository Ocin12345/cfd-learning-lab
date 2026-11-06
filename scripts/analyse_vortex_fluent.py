"""Audit saved cylinder transient Fluent evidence without starting a solver.

--run is the native fluent or fluent-extended folder, --out is a fresh E: folder.
False checks are retained and return exit code 2. A startup transient need not
have constant Cd/Cl. No shedding frequency, benchmark, mesh-independence or
physical-validation acceptance is imposed on a short 2/8 second calculation.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import io
import json
import math
import re
import shutil
import zipfile

EQUATIONS = ('continuity', 'x-velocity', 'y-velocity')
REFERENCE = dict(density=1.0, viscosity=0.001, velocity=1.0, length=0.1, area=0.1, depth=1.0)
NUMBER = r'([-+0-9.eE]+)'


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, values):
    path.write_text(json.dumps(values, indent=2, allow_nan=False, default=str)+'\n', encoding='utf-8')


def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def mass_values(text):
    values = {}
    for label in ('inlet', 'outlet', 'Net'):
        matches = re.findall(r'^\s*'+label+r'\s+'+NUMBER+r'\s*$', text, re.M)
        if not matches:
            raise ValueError('Native mass report lacks '+label)
        values[label] = float(matches[-1])
    if not all(map(math.isfinite, values.values())) or values['inlet'] <= 0 or values['outlet'] >= 0:
        raise ValueError('Invalid native boundary mass values or signs')
    values['relative_imbalance'] = abs(values['Net'])/abs(values['inlet'])
    return values


def terminal_rows(text):
    """Associate the last printed residual row with each actual native time marker."""
    output, last = [], None
    residual = re.compile(r'^\s*(\d+)\s+'+NUMBER+r'\s+'+NUMBER+r'\s+'+NUMBER+r'\s+\d+:\d+:\d+\s+\d+\s*$')
    marker = re.compile(r'Flow\s+time\s*=\s*'+NUMBER+r'\s*s\s*,?\s*time\s*step\s*=\s*(\d+)', re.I)
    for line in text.splitlines():
        match = residual.match(line)
        if match:
            values = [float(match[i]) for i in (2, 3, 4)]
            if all(map(math.isfinite, values)):
                last = dict(iteration=int(match[1]), **dict(zip(EQUATIONS, values)))
        match = marker.search(line)
        if match and last is not None:
            output.append(dict(time_s=float(match[1]), step=int(match[2]), **last))
            last = None
    # A resumed transcript can repeat a reporting marker. Preserve its last
    # actual terminal report, while retaining the raw transcript untouched.
    return sorted({row['step']: row for row in output}.values(), key=lambda row: row['step'])


def scheme_flow_times(text):
    return [float(value.replace('D', 'e').replace('d', 'e')) for value in
            re.findall(r'\(\s*flow-time\s*(?:\.\s*)?([-+0-9.eEdD]+)\s*\)', text)]


def xy_rows(text):
    rows, label = {}, None
    for line in text.splitlines():
        match = re.search(r'xy/key/label\s+"([^"]+)"', line)
        if match:
            label = match[1]
        elif label in EQUATIONS:
            values = line.split()
            if len(values) == 2:
                try:
                    iteration, value = float(values[0]), float(values[1])
                    if iteration == int(iteration) and math.isfinite(value):
                        rows.setdefault(int(iteration), {})[label] = value
                except ValueError:
                    pass
    return rows


def close(value, expected, tolerance=1e-9):
    try:
        return math.isfinite(float(value)) and math.isclose(float(value), expected, rel_tol=0., abs_tol=tolerance)
    except (ValueError, TypeError):
        return False


def settings_match(state):
    solver = state['general']['solver']
    inlet = state['inlet']['momentum']
    components = inlet['velocity_components']
    methods, run = state['methods'], state['run_calculation']
    checks = dict(
        solver=solver['type'] == 'pressure-based' and solver['two_dim_space'] == 'planar'
            and solver['time'] in ('transient', 'unsteady-2nd-order'),
        laminar=state['viscous']['model'] == 'laminar',
        material=close(state['material']['density']['value'], 1.)
            and close(state['material']['viscosity']['value'], .001)
            and state['cell_zones']['fluid']['general']['material'] == 'air',
        inlet=inlet['velocity_specification_method'] == 'Components'
            and components[0]['option'] == 'profile'
            and components[0]['profile_name'] == 'parabolic-inlet'
            and components[0]['field_name'] == 'u' and close(components[1]['value'], 0.),
        outlet=close(state['outlet']['momentum']['gauge_pressure']['value'], 0.),
        walls=all(state['walls'][name]['momentum']['wall_motion'] == 'Stationary Wall'
            and state['walls'][name]['momentum']['shear_condition'] == 'No Slip'
            for name in ('walls', 'cylinder')),
        references=all(close(state['references'][key], value) for key, value in REFERENCE.items()),
        numerics=methods['p_v_coupling']['flow_scheme'] == 'PISO'
            and methods['transient_formulation'] == 'unsteady-2nd-order'
            and methods['spatial_discretization']['gradient_scheme'] == 'least-square-cell-based'
            and methods['spatial_discretization']['discretization_scheme'].get('mom') == 'second-order-upwind'
            and methods['spatial_discretization']['discretization_scheme'].get('pressure') == 'second-order',
        fixed_step=run['transient_controls']['type'] == 'Fixed'
            and close(run['parameters']['time_step_size'], .005)
            and run['parameters']['max_iter_per_time_step'] == 20,
        residual_targets=all(close(state['residual']['equations'][name]['absolute_criteria'], 1e-5, 1e-12)
            and state['residual']['equations'][name]['check_convergence'] for name in EQUATIONS))
    return checks


def hdf_residuals(path):
    import h5py
    import numpy as np
    rows, indices = {}, None
    with h5py.File(path, 'r') as handle:
        for name in EQUATIONS:
            group = handle['results/residuals/phase-1/'+name]
            raw, values = group['iterations'][:], group['data'][:]
            index = np.asarray(raw).reshape(-1)
            if (not len(index) or not np.array_equal(index, index.astype(int))
                    or np.any(np.diff(index) <= 0) or values.ndim != 2
                    or values.shape[0] != len(index) or values.shape[1] < 2
                    or not np.isfinite(values[:, :2]).all() or np.any(values[:, 1] <= 0)):
                raise ValueError('Invalid saved native residual layout/scaling: '+name)
            if indices is not None and not np.array_equal(indices, index):
                raise ValueError('Native equation iteration histories differ')
            indices = index
            scaled = values[:, 0]/values[:, 1]
            if not np.isfinite(scaled).all() or np.any(scaled < 0):
                raise ValueError('Nonfinite or negative native scaled residual')
            for iteration, value in zip(index, scaled):
                rows.setdefault(int(iteration), dict(iteration=int(iteration)))[name] = float(value)
    return [rows[index] for index in sorted(rows)]


def hdf_fields(path, expected_cells):
    import h5py
    import numpy as np
    result = {}
    with h5py.File(path, 'r') as handle:
        fields = handle['results/1/phase-1/cells']
        for name, expected in (('SV_U', None), ('SV_V', None), ('SV_P', None), ('SV_DENSITY', 1.), ('SV_MU_LAM', .001)):
            group = fields[name]
            seen = np.zeros(expected_cells, dtype=bool)
            low, high = math.inf, -math.inf
            for dataset in group.values():
                first = int(np.asarray(dataset.attrs['minId']).ravel()[0])-1
                last = int(np.asarray(dataset.attrs['maxId']).ravel()[0])
                values = np.asarray(dataset[:]).reshape(-1)
                if (not 0 <= first < last <= expected_cells or len(values) != last-first
                        or seen[first:last].any() or not np.isfinite(values).all()):
                    raise ValueError('Invalid native field IDs or nonfinite values: '+dataset.name)
                if expected is not None and not np.allclose(values, expected, rtol=1e-10, atol=1e-12):
                    raise ValueError('Saved native material field differs: '+name)
                seen[first:last] = True
                low, high = min(low, float(values.min())), max(high, float(values.max()))
            if not seen.all():
                raise ValueError('Missing native field cell IDs: '+name)
            result[name] = dict(minimum=low, maximum=high, cells=expected_cells)
    return result


def hdf_times(path):
    import h5py
    result = []
    with h5py.File(path, 'r') as handle:
        if 'settings' in handle:
            def inspect(name, item):
                if isinstance(item, h5py.Dataset) and item.size*item.dtype.itemsize <= 2*1024**2:
                    raw = item[()]
                    if isinstance(raw, bytes):
                        text = raw.decode('utf-8', errors='replace')
                    elif getattr(getattr(raw, 'dtype', None), 'kind', '') in ('S', 'O', 'U'):
                        text = ' '.join(value.decode('utf-8', errors='replace')
                            if isinstance(value, bytes) else str(value) for value in raw.flat)
                    else:
                        text = raw.tobytes().decode('utf-8', errors='replace') if hasattr(raw, 'tobytes') else str(raw)
                    for value in scheme_flow_times(text.replace('\x00', '')):
                        result.append(dict(dataset='settings/'+name, flow_time_s=value))
            handle['settings'].visititems(inspect)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    if not run.is_dir() or out.exists() or out.drive.upper() != 'E:':
        raise ValueError('Use an existing native run and fresh output directory on E:')
    out.mkdir(parents=True)
    checks, evidence, errors, before = {}, {}, {}, {}

    def remember(path):
        path = Path(path).resolve()
        if path not in before:
            before[path] = sha256(path)
        return path

    def attempt(name, operation, required=True):
        try:
            value = operation()
            if required:
                checks[name] = True
            return value
        except Exception as error:
            errors[name] = type(error).__name__+': '+str(error)
            if required:
                checks[name] = False
            return None

    def load(name):
        return attempt('read_'+name, lambda: json.loads(remember(run/name).read_text(encoding='utf-8'))) or {}

    status, summary = load('status.json'), load('summary.json')
    settings, quality = load('final-settings.json'), load('geometry-quality.json')
    geometry, provenance = load('geometry-setup.json'), load('code-provenance.json')
    setup, targets = load('mesh-setup.json'), load('declared-targets.json')
    manifest = load('snapshots.json')
    source_path = attempt('frozen_driver_present', lambda: remember(run/'executed-script.py'))
    checks['frozen_driver_matches_prelaunch_hash'] = bool(source_path and sha256(source_path) == provenance.get('source_sha256'))
    checks['driver_completed_native_final_saved'] = (status.get('stage') in ('completed', 'completed_with_export_errors')
        and summary.get('native_final_saved') is True)
    checks['declared_targets_match'] = (close(targets.get('dt_s', 0), .005)
        and close(targets.get('residual_target', 0), 1e-5, 1e-12)
        and close(targets.get('relative_mass_imbalance_target', 0), 1e-4, 1e-12)
        and targets.get('max_inner_iterations_per_time_step') == 20)
    detail = attempt('native_settings_structure', lambda: settings_match(settings)) or {}
    evidence['settings_checks'] = detail
    checks['native_settings_match_declared_physics_references_and_numerics'] = bool(detail and all(detail.values()))
    checks['native_mesh_quality_gate'] = bool(quality.get('native_numeric_gate_pass') is True
        and quality.get('checked_before_flow') is True and quality.get('minimum_orthogonal_quality', -1) >= .2
        and quality.get('minimum_cell_volume_m3', 0) > 0)
    checks['geometry_and_mesh_provenance_match'] = (geometry.get('kind') == 'cylinder-channel-Re100'
        and geometry.get('references') == REFERENCE and setup.get('mesh_sha256') == provenance.get('mesh_sha256')
        and before.get((run/'mesh-setup.json').resolve()) == provenance.get('mesh_setup_sha256'))
    evidence['native_mesh_quality'] = quality
    mass = attempt('native_mass_report_valid', lambda: mass_values(remember(run/'mass-balance.txt').read_text(encoding='utf-8'))) or {}
    checks['relative_mass_imbalance_le_1e_4'] = mass.get('relative_imbalance', math.inf) <= 1e-4
    evidence['boundary_mass_balance'] = mass
    actual_time = summary.get('end_time_s')
    evidence['physical_time'] = dict(driver_actual_rpc_end_time_s=actual_time,
        driver_status_time_s=status.get('end_time_s'), driver_summary_start_time_s=summary.get('start_time_s'))
    checks['driver_rpc_time_agrees_with_final_status'] = bool(isinstance(actual_time, (int, float))
        and math.isfinite(actual_time) and 0 < actual_time <= 8
        and close(status.get('end_time_s', -1), actual_time))
    settings_time = settings.get('run_calculation', {}).get('transient_controls', {}).get('flow_time')
    evidence['physical_time']['saved_final_settings_flow_time_s'] = settings_time
    if settings_time is not None and isinstance(actual_time, (int, float)):
        checks['available_final_settings_time_matches_rpc'] = close(settings_time, actual_time)
    final_case, final_data = run/'vortex-final.cas.h5', run/'vortex-final.dat.h5'
    attempt('native_final_case_present', lambda: remember(final_case))
    attempt('native_final_data_present', lambda: remember(final_data))
    residual = attempt('native_hdf_residual_history_valid', lambda: hdf_residuals(final_data)) or []
    if residual:
        final_residual = {name: residual[-1][name] for name in EQUATIONS}
        evidence['residuals'] = dict(final_iteration=residual[-1]['iteration'],
            final_scaled=final_residual, full_native_hdf_history_records=len(residual),
            definition='Actual native saved HDF residual data column0/column1',
            cross_solver_residual_norm_equivalence=False)
        checks['final_current_step_residuals_le_1e_5'] = all(value <= 1e-5 for value in final_residual.values())
    else:
        checks['final_current_step_residuals_le_1e_5'] = False
    write_csv(out/'fluent-residuals.csv', residual, ['iteration', *EQUATIONS])
    xy_path = run/'residuals.xy'
    xy = attempt('native_xy_read', lambda: xy_rows(remember(xy_path).read_text(encoding='utf-8', errors='replace')), required=False)
    evidence['residual_xy'] = dict(native_file_present=xy_path.is_file(), native_bytes=xy_path.stat().st_size if xy_path.is_file() else None,
        parsed_scalar_pairs=sum(len(row) for row in (xy or {}).values()),
        csv_source='Intact native final HDF history; original XY, including any empty export, is preserved')
    hdf_by_iteration = {row['iteration']: row for row in residual}
    checks['available_native_xy_matches_hdf'] = True
    if xy:
        checks['available_native_xy_matches_hdf'] = all(index in hdf_by_iteration and
            all(math.isclose(value, hdf_by_iteration[index][name], rel_tol=6e-6, abs_tol=1e-14)
                for name, value in row.items()) for index, row in xy.items())
    elif xy_path.is_file() and xy_path.stat().st_size:
        checks['available_native_xy_matches_hdf'] = False
    transcript = attempt('native_transcript_present', lambda: remember(run/'solver-transcript.txt').read_text(encoding='utf-8', errors='replace')) or ''
    terminal = terminal_rows(transcript)
    write_csv(out/'fluent-step-terminal-residuals.csv', terminal, ['time_s', 'step', 'iteration', *EQUATIONS])
    failed_steps = [row['step'] for row in terminal if any(row[name] > 1e-5 for name in EQUATIONS)]
    evidence['per_step_terminals'] = dict(actual_native_time_markers_parsed=len(terminal),
        steps_above_target=failed_steps, count_above_target=len(failed_steps),
        all_recorded_terminals_le_1e_5=bool(terminal and not failed_steps),
        expected_this_session_steps=summary.get('additional_time_steps'),
        all_this_session_steps_observed=len(terminal) == summary.get('additional_time_steps'),
        acceptance='Final current-step target is mandatory; full per-step terminal values retained as diagnostics')
    if terminal and residual:
        evidence['physical_time']['final_transcript_marker_s'] = terminal[-1]['time_s']
        checks['final_transcript_marker_matches_rpc_and_hdf_residual'] = bool(
            isinstance(actual_time, (int, float)) and close(terminal[-1]['time_s'], actual_time)
            and terminal[-1]['iteration'] == residual[-1]['iteration']
            and all(math.isclose(terminal[-1][name], residual[-1][name], rel_tol=6e-5, abs_tol=1e-14) for name in EQUATIONS))
    else:
        evidence['per_step_terminals']['limitation'] = 'Native terminal time markers not completely available; no inferred per-step residuals'
    native_times = []
    for path in (final_case, final_data):
        values = attempt('physical_time_metadata_'+path.name, lambda path=path: hdf_times(path), required=False) or []
        native_times.extend(dict(file=path.name, **value) for value in values)
    evidence['physical_time']['native_hdf_metadata'] = native_times
    evidence['physical_time']['native_hdf_metadata_available'] = bool(native_times)
    if native_times:
        checks['available_native_hdf_physical_times_match_rpc'] = bool(isinstance(actual_time, (int, float))
            and all(close(row['flow_time_s'], actual_time) for row in native_times))
    fields = attempt('final_native_fields_finite_and_materials_match', lambda: hdf_fields(final_data, int(setup['cells'])))
    evidence['final_native_field_ranges'] = fields
    snapshots, snapshot_errors = [], []
    if isinstance(manifest, list):
        for item in manifest:
            try:
                path = Path(item['data'])
                path = path if path.is_absolute() else run/path
                remember(path)
                ranges = hdf_fields(path, int(setup['cells']))
                snapshots.append(dict(time_s=float(item['time']), step=int(item['step']), data=str(path.resolve()),
                                      bytes=path.stat().st_size, field_ranges=ranges))
            except Exception as error:
                snapshot_errors.append(dict(item=item, error=type(error).__name__+': '+str(error)))
    evidence['snapshots'] = dict(count=len(snapshots), errors=snapshot_errors,
        last_time_s=snapshots[-1]['time_s'] if snapshots else None, all_checked_fields=['SV_U','SV_V','SV_P','SV_DENSITY','SV_MU_LAM'])
    checks['all_snapshot_fields_finite_and_materials_match'] = bool(snapshots and not snapshot_errors)
    expected_count = int(round(actual_time/.025))+1 if (isinstance(actual_time, (int, float))
        and math.isfinite(actual_time) and 0 < actual_time <= 8) else -1
    checks['snapshot_times_and_steps_match_full_saved_time_grid'] = bool(
        len(snapshots) == expected_count and all(close(row['time_s'], index*.025)
            and row['step'] == index*5 for index, row in enumerate(snapshots)))
    checks['final_snapshot_time_matches_driver_rpc'] = bool(snapshots and isinstance(actual_time, (int, float))
        and close(snapshots[-1]['time_s'], actual_time))
    force_path = run/'force-history.csv'
    raw_force = attempt('native_force_history_present', lambda: list(csv.DictReader(io.StringIO(remember(force_path).read_text(encoding='utf-8'))))) or []
    forces = attempt('native_force_history_numeric', lambda: [{key: int(value) if key == 'step' else float(value)
                for key, value in row.items()} for row in raw_force]) or []
    if any(not math.isfinite(value) for row in forces for value in row.values()):
        errors['native_force_history_numeric'] = 'Nonfinite native force-history value; original file retained'
        checks['native_force_history_numeric'] = False
        forces = []
    checks['force_times_finite_and_match_snapshot_grid'] = bool(len(forces) == expected_count
        and all(all(math.isfinite(value) for value in row.values()) and close(row['time_s'], index*.025)
            and row['step'] == index*5 for index, row in enumerate(forces)))
    checks['force_reference_coefficients_match'] = bool(forces and all(
        math.isclose(row['Cd'], row['drag_force_N']/.05, rel_tol=2e-4, abs_tol=1e-8)
        and math.isclose(row['Cl'], row['lift_force_N']/.05, rel_tol=2e-4, abs_tol=1e-8) for row in forces))
    write_csv(out/'fluent-forces.csv', forces, list(forces[0]) if forces else ['time_s','Cd','Cl'])
    evidence['forces'] = dict(samples=len(forces), sample_interval_s=.025,
        final=forces[-1] if forces else None, interpretation='Actual native force samples; variation is expected in an unsteady solution',
        frequency_estimate_produced=False, periodicity_established=False,
        no_separate_frequency_validation='Short startup duration does not establish a converged shedding frequency')
    remember(Path(__file__))
    checks['native_inputs_unchanged_during_analysis'] = all(sha256(path) == digest for path, digest in before.items())
    record = dict(solver='ANSYS Fluent Student v261', native_run=str(run), evidence=evidence,
        check_results=checks, errors=errors, all_declared_checks_passed=all(checks.values()),
        physical_validation=False, benchmark_acceptance=False, mesh_independence_established=False,
        time_step_independence_established=False, periodic_shedding_established=False,
        interpretation='Numerical evidence audit of an unsteady startup; no constant Cd/Cl or benchmark acceptance criterion')
    write_json(out/'fluent-checks.json', record)
    if source_path:
        shutil.copyfile(source_path, out/'executed-fluent.py')
    shutil.copyfile(Path(__file__), out/'executed-fluent-analyser.py')
    sources = dict(native_and_sources_unchanged_during_analysis=checks['native_inputs_unchanged_during_analysis'],
        sources=[dict(path=str(path), sha256=digest) for path, digest in before.items()],
        outputs=[dict(path=path.name, sha256=sha256(path)) for path in sorted(out.iterdir()) if path.is_file()],
        solver_started_or_modified=False, native_logs_archived_or_published=False,
        source_archive_matches_prelaunch_hash=checks['frozen_driver_matches_prelaunch_hash'])
    write_json(out/'fluent-source-manifest.json', sources)
    with zipfile.ZipFile(out/'fluent-evidence.zip', 'x', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(out.iterdir()):
            if path.is_file() and path.suffix != '.zip':
                archive.write(path, arcname=path.name)
    print(json.dumps(dict(all_declared_checks_passed=record['all_declared_checks_passed'],
        check_results=checks, errors=errors, checks_file=str(out/'fluent-checks.json')), indent=2, allow_nan=False))
    return 0 if record['all_declared_checks_passed'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
