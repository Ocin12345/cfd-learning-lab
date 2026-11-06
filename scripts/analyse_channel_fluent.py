"""Check saved channel Fluent results without starting or modifying a solver.

Example: --run E:/.../2026-11-01-channel-flow-v001 --out REPO/results/channel-flow-v001
Native logs stay private; only derived checks and their hashes are exported.
"""
from pathlib import Path
import argparse
import csv
import json
import os
import re

os.environ.setdefault('MPLCONFIGDIR', 'E:/CFD-Learning-Lab/temp/matplotlib-channel')
import h5py
import numpy as np
from channel_fluent import EQUATIONS, residual_rows, residual_csv
from plot_channel_flow import hdf_indexed, map_and_check, read_fluent, require, sha256


def write_json(path, values):
    path.write_text(json.dumps(values, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def hdf_residuals(path):
    with h5py.File(path, 'r') as handle:
        values = {}
        iterations = None
        for name in EQUATIONS:
            group = handle['results/residuals/phase-1/' + name]
            indices = group['iterations'][:].astype(int)
            data = group['data'][:]
            require(np.array_equal(group['iterations'][:], indices), 'Nonintegral residual iterations')
            require(data.shape == (len(indices), 4) and np.all(data[:, 1] > 0),
                    'Unexpected native residual layout or scaling')
            require(np.array_equal(indices, np.arange(1, len(indices) + 1)), 'Residual iteration gap')
            require(iterations is None or np.array_equal(iterations, indices), 'Equation iteration mismatch')
            iterations = indices
            values[name] = data[:, 0] / data[:, 1]
            require(np.isfinite(values[name]).all(), 'Nonfinite scaled residual')
    return iterations, values


def same_native_mesh(baseline, final):
    with h5py.File(baseline, 'r') as first, h5py.File(final, 'r') as last:
        names = []
        first['meshes/1'].visititems(lambda name, item: names.append(name)
                                    if isinstance(item, h5py.Dataset) else None)
        for name in names:
            require(np.array_equal(first['meshes/1/' + name][:], last['meshes/1/' + name][:]),
                    'Saved mesh dataset changed: ' + name)
    return len(names)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    native = run / 'fluent'
    names = ('fluent-residuals.csv', 'fluent-checks.json', 'fluent-source-manifest.json')
    require(not any((out / name).exists() for name in names), 'Refuse to replace existing analysis outputs')
    setup_path = run / 'mesh/setup.json'
    sources = [setup_path, run / 'mesh/channel.msh', *[native / name for name in (
        'channel-baseline.cas.h5', 'channel-baseline.dat.h5', 'channel-final.cas.h5',
        'channel-final.dat.h5', 'residuals-baseline.xy', 'residuals.xy', 'mass-balance.txt',
        'declared-targets.json', 'executed-script.py', 'code-provenance.json',
        'loaded-settings.json', 'final-settings.json', 'status.json')]]
    startup_files = list(native.glob('fluent-*.trn'))
    require(len(startup_files) == 1, 'Expected exactly one private startup log')
    private_logs = [native / 'solver-transcript.txt', startup_files[0]]
    before = {path: sha256(path) for path in sources + private_logs}
    setup = json.loads(setup_path.read_text(encoding='utf-8'))
    targets = json.loads((native / 'declared-targets.json').read_text(encoding='utf-8'))
    baseline_iterations, baseline_residuals = hdf_residuals(native / 'channel-baseline.dat.h5')
    iterations, residuals = hdf_residuals(native / 'channel-final.dat.h5')
    require(int(iterations[-1] - baseline_iterations[-1]) == 50, 'Saved endpoints are not 50 iterations apart')
    for name in EQUATIONS:
        require(np.array_equal(baseline_residuals[name], residuals[name][:len(baseline_iterations)]),
                'Native baseline/final residual histories disagree')
    xy = residual_rows(native / 'residuals-baseline.xy')
    xy_pairs = sum(len(row) - 1 for row in xy.values())
    for index, row in xy.items():
        for name in EQUATIONS:
            require(row[name] == float(format(float(np.float32(residuals[name][index - 1])), '.6g')),
                    'Native XY differs from float32-cast six-significant-figure HDF residual')

    mesh_datasets = same_native_mesh(native / 'channel-baseline.cas.h5', native / 'channel-final.cas.h5')
    final = read_fluent(run, setup)
    final_analytical = map_and_check(final, setup)
    baseline = dict(centres=final['centres'].copy(), topology='Verified identical saved native mesh and cell IDs')
    with h5py.File(native / 'channel-baseline.dat.h5', 'r') as handle:
        fields = handle['results/1/phase-1/cells']
        baseline.update({key: hdf_indexed(fields[name], setup['cells'])
                         for key, name in [('u', 'SV_U'), ('v', 'SV_V'), ('p', 'SV_P')]})
    map_and_check(baseline, setup)
    x = (np.arange(setup['nx']) + .5) * setup['length_m'] / setup['nx']
    station = int(np.argmin(np.abs(x - .8)))
    left, right = [int(np.argmin(np.abs(x - location))) for location in (.6, .9)]
    mean = setup['inlet_velocity_m_s']
    change = (final['grid']['u'][:, station] - baseline['grid']['u'][:, station]) / mean
    outlet_change = (final['grid']['u'][:, -1] - baseline['grid']['u'][:, -1]) / mean
    pressure_drop = {label: float(field['grid']['p'][:, left].mean() - field['grid']['p'][:, right].mean())
                     for label, field in [('baseline', baseline), ('final', final)]}
    pressure_change = abs(pressure_drop['final'] - pressure_drop['baseline']) / abs(pressure_drop['baseline'])
    mass_text = (native / 'mass-balance.txt').read_text(encoding='utf-8')
    mass = {label: float(re.search(r'^\s*' + label + r'\s+([-+0-9.eE]+)\s*$', mass_text, re.M)[1])
            for label in ('inlet', 'outlet', 'Net')}
    relative_mass = abs(mass['Net']) / abs(mass['inlet'])
    transcript = private_logs[0].read_text(encoding='utf-8', errors='replace')
    startup = private_logs[1].read_text(encoding='utf-8', errors='replace')
    markers = [int(value) for value in re.findall(r'!\s+(\d+)\s+solution is converged', transcript)]
    orthogonal = float(re.search(r'Minimum Orthogonal Quality\s*=\s*([-+0-9.eE]+)', transcript)[1])
    aspect = float(re.search(r'Maximum Aspect Ratio\s*=\s*([-+0-9.eE]+)', transcript)[1])
    licensed_parallel_limit = int(re.search(r'Your license enables (\d+)-way parallel', startup)[1])
    processes = len(re.findall(r'^n\d+\*?\s+', startup, re.M))
    final_printed = re.findall(r'^\s*' + str(iterations[-1])
                              + r'\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+\d+:\d+:\d+\s+\d+\s*$',
                              transcript, re.M)
    require(bool(final_printed), 'Final residual row absent from private solver evidence')
    for name, value in zip(EQUATIONS, final_printed[-1]):
        require(float(value) == float(format(float(residuals[name][-1]), '.5g')),
                'Final HDF residual differs from printed solver evidence')
    provenance = json.loads((native / 'code-provenance.json').read_text(encoding='utf-8'))
    archive = out / 'executed-fluent.py'
    require(archive.is_file() and sha256(archive) == before[native / 'executed-script.py']
            == provenance['source_sha256'], 'Original executed source archive differs')
    checks = dict(live_student_licence_and_two_process_startup=('student version' in startup.lower() and processes == 2),
                  residual_final_le_1e_6=all(values[-1] <= targets['residual_max'] for values in residuals.values()),
                  residual_last50_le_1e_6=all(np.max(values[-50:]) <= targets['residual_max'] for values in residuals.values()),
                  relative_mass_imbalance_le_1e_5=relative_mass <= targets['relative_mass_imbalance_max'],
                  downstream_profile_endpoint_change_le_0p1percent=float(np.max(np.abs(change))) <= .001,
                  outlet_profile_endpoint_change_le_0p1percent=float(np.max(np.abs(outlet_change))) <= .001,
                  pressure_drop_endpoint_change_le_0p1percent=pressure_change <= .001,
                  native_xy_pairs_equal_hdf_printed_values=True,
                  analytical_targets_passed=final_analytical['analytical_targets_passed'])
    out.mkdir(parents=True, exist_ok=True)
    recovered = residual_csv(native / 'residuals.xy', out / names[0], native / 'channel-final.dat.h5')
    require(len(recovered) == len(iterations), 'Recovered CSV iteration count differs')
    with (out / names[0]).open(newline='', encoding='utf-8') as stream:
        csv_rows = {int(row['iteration']): row for row in csv.DictReader(stream)}
    for index in iterations:
        for name in EQUATIONS:
            require(float(csv_rows[int(index)][name]) == recovered[int(index)][name]
                    == float(format(float(np.float32(residuals[name][int(index) - 1])), '.6g')), 'Recovered CSV/HDF mismatch')
    for index, row in xy.items():
        for name in EQUATIONS:
            require(float(csv_rows[index][name]) == row[name], 'Actual native XY pair changed in recovered CSV')
    unchanged = all(sha256(path) == digest for path, digest in before.items())
    require(unchanged, 'A native input changed during analysis')
    record = dict(solver='ANSYS Fluent Student 2026 R1', mesh_cells=setup['cells'],
                  startup=dict(live_licence_confirmed=checks['live_student_licence_and_two_process_startup'],
                               actual_solver_processes=processes, licensed_parallel_limit=licensed_parallel_limit),
                  mesh_quality=dict(minimum_orthogonal_quality=orthogonal, maximum_aspect_ratio=aspect,
                                    identical_baseline_final_mesh_datasets=mesh_datasets),
                  residuals=dict(initial_native_convergence_marker=min(markers), baseline_iteration=int(baseline_iterations[-1]),
                                 final_iteration=int(iterations[-1]), native_hdf_history_records_per_equation=len(iterations),
                                 final_scaled_residuals={key: float(value[-1]) for key, value in residuals.items()},
                                 maximum_scaled_residuals_last50={key: float(np.max(value[-50:])) for key, value in residuals.items()},
                                 definition='Saved HDF residual data column0 divided by column1; printed transcript final row agrees to five significant figures',
                                 cross_solver_residual_norm_equivalence=False),
                  boundary_mass_balance=dict(signed_solver_report_kg_s=mass, relative_imbalance=relative_mass,
                                             definition='abs(native reported Net)/abs(native reported inlet mass rate)'),
                  final50_endpoint_stability=dict(iterations=[int(baseline_iterations[-1]), int(iterations[-1])],
                      profile_x_m=float(x[station]), profile_rms_change_over_Umean=float(np.sqrt(np.mean(change**2))),
                      profile_maximum_change_over_Umean=float(np.max(np.abs(change))),
                      outlet_profile_last_cell_x_m=float(x[-1]), outlet_profile_rms_change_over_Umean=float(np.sqrt(np.mean(outlet_change**2))),
                      outlet_profile_maximum_change_over_Umean=float(np.max(np.abs(outlet_change))),
                      pressure_drop_interval_m=[float(x[left]), float(x[right])], pressure_drop_Pa=pressure_drop,
                      pressure_drop_relative_change=float(pressure_change),
                      limitation='Two saved endpoints only; does not establish every-iteration field stability between222 and272'),
                  final_analytical_checks=final_analytical,
                  optional_export_recovery=dict(original_native_status='completed_with_export_errors',
                      original_final_xy_bytes=(native / 'residuals.xy').stat().st_size,
                      actual_failure='Final native residual XY export was empty; labels in the intact baseline XY already matched the parser',
                      recovered_from='Complete saved final HDF scaled residual ratios, float32-cast then printed to six significant figures, matching all intact native XY pairs',
                      recovered_csv_records=len(csv_rows), recovered_csv_scalar_pairs=len(csv_rows)*len(EQUATIONS),
                      actual_baseline_xy_pairs_equal_csv=xy_pairs, all_actual_xy_pairs_checked=True,
                      original_executed_source_preserved=True, native_inputs_unchanged=True, solver_rerun_or_continuation=False),
                  check_results=checks, all_declared_checks_passed=all(checks.values()),
                  mesh_independence_established=False, experimental_validation=False)
    write_json(out / names[1], record)
    scripts = [Path(__file__).resolve(), Path(__file__).with_name('channel_fluent.py'),
               Path(__file__).with_name('plot_channel_flow.py'), archive]
    manifest = dict(inputs_unchanged_during_analysis=unchanged,
                    sources=[dict(path=str(path), sha256=before[path]) for path in sources],
                    private_log_evidence_hashes=[dict(role=role, sha256=before[path])
                        for role, path in zip(('solver numerical evidence', 'startup and licence evidence'), private_logs)],
                    scripts=[dict(path=str(path), sha256=sha256(path)) for path in scripts],
                    outputs=[dict(path=str(out / name), sha256=sha256(out / name)) for name in names[:2]],
                    executed_source_archive_matches_prelaunch_sha256=True,
                    native_logs_published=False)
    write_json(out / names[2], manifest)
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0 if all(checks.values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())