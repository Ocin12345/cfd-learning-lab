"""Check saved Fluent backward-step data; never start or modify a solver.

Run with --run the native run root and --out a results folder. Only derived
checks, residuals, hashes and executed Python sources are exported; logs stay private.
"""
from pathlib import Path
import argparse
import csv
import json
import re
import h5py
import numpy as np
from step_fluent import EQUATIONS, residual_csv, residual_rows
from plot_step_flow import (grid_coordinates, hdf_indexed, map_and_check,
                            read_fluent, require, sha256, validate_setup)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def history(path):
    values, iterations = {}, None
    with h5py.File(path, 'r') as handle:
        for name in EQUATIONS:
            group = handle['results/residuals/phase-1/' + name]
            raw, data = group['iterations'][:], group['data'][:]
            index = raw.astype(int)
            require(np.array_equal(raw, index) and np.array_equal(index, np.arange(1, len(index) + 1)),
                    'Residual history must contain every integer iteration from 1')
            require(data.shape == (len(index), 4) and np.all(data[:, 1] > 0), 'Unexpected residual scaling')
            require(iterations is None or np.array_equal(iterations, index), 'Equation history mismatch')
            iterations = index
            values[name] = data[:, 0] / data[:, 1]
            require(np.isfinite(values[name]).all() and np.all(values[name] >= 0), 'Invalid native residuals')
    return iterations, values


def same_mesh(first, last):
    with h5py.File(first, 'r') as a, h5py.File(last, 'r') as b:
        items = []
        a['meshes/1'].visititems(lambda name, obj: items.append(name))
        other = []
        b['meshes/1'].visititems(lambda name, obj: other.append(name))
        require(items == other, 'Saved native mesh object names changed')
        datasets = 0
        for name in [''] + items:
            left, right = [f['meshes/1' + ('/' + name if name else '')] for f in (a, b)]
            require(type(left) == type(right) and set(left.attrs) == set(right.attrs), 'Native mesh type/attribute mismatch')
            require(all(np.array_equal(left.attrs[k], right.attrs[k]) for k in left.attrs), 'Native mesh attributes changed')
            if isinstance(left, h5py.Dataset):
                require(np.array_equal(left[:], right[:]), 'Native mesh dataset changed: ' + name)
                datasets += 1
    return datasets


def settings_values(state):
    return dict(solver=state['general']['solver'], viscous=state['viscous']['model'],
                material=state['material'], fluid_material=state['cell_zones']['fluid']['general']['material'],
                inlet=state['inlet']['momentum'], outlet_gauge_Pa=state['outlet']['momentum']['gauge_pressure']['value'],
                walls=state['walls']['momentum'], methods=state['methods']['spatial_discretization'],
                coupling=state['methods']['p_v_coupling']['flow_scheme'],
                references=state['references'], residual_equations=state['residual']['equations'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    native = run / 'fluent'
    names = ['fluent-checks.json', 'fluent-residuals.csv', 'fluent-source-manifest.json',
             'executed-fluent.py', 'executed-fluent-analyser.py']
    require(not any((out / n).exists() for n in names), 'Refuse to replace existing Fluent analysis outputs')
    paths = [run / 'mesh/setup.json', run / 'mesh/step.msh'] + [native / name for name in (
        'step-baseline.cas.h5', 'step-baseline.dat.h5', 'step-final.cas.h5', 'step-final.dat.h5',
        'residuals-baseline.xy', 'residuals.xy', 'iteration-history.csv', 'mass-balance.txt',
        'loaded-settings.json', 'final-settings.json', 'mesh-setup.json', 'geometry-setup.json',
        'code-provenance.json', 'declared-targets.json', 'executed-script.py', 'history.json', 'status.json')]
    startup = list(native.glob('fluent-*.trn'))
    require(len(startup) == 1, 'Expected one original Fluent startup transcript')
    private_logs = [native / 'solver-transcript.txt', startup[0]]
    scripts = [Path(__file__).resolve(), Path(__file__).with_name('step_fluent.py'),
               Path(__file__).with_name('plot_step_flow.py')]
    before = {p: sha256(p) for p in paths + private_logs + scripts}
    setup = validate_setup(json.loads(paths[0].read_text(encoding='utf-8')))
    targets = json.loads((native / 'declared-targets.json').read_text(encoding='utf-8'))
    provenance = json.loads((native / 'code-provenance.json').read_text(encoding='utf-8'))
    require(before[run / 'mesh/step.msh'] == setup['mesh_sha256'] == provenance['mesh_sha256'], 'Mesh provenance mismatch')
    require(before[native / 'executed-script.py'] == before[scripts[1]] == provenance['source_sha256'], 'Executed driver differs from frozen source')
    require(before[paths[0]] == before[native / 'mesh-setup.json'] == provenance['mesh_setup_sha256'], 'Mesh setup snapshot differs')
    first_index, first_residual = history(native / 'step-baseline.dat.h5')
    index, residual = history(native / 'step-final.dat.h5')
    require(int(index[-1] - first_index[-1]) == 100, 'Saved endpoints are not 100 iterations apart')
    require(all(np.array_equal(first_residual[k], residual[k][:len(first_index)]) for k in EQUATIONS), 'Native residual histories disagree')
    xy = residual_rows(native / 'residuals-baseline.xy')
    for i, row in xy.items():
        require(all(row[k] == float(format(float(np.float32(residual[k][i - 1])), '.6g')) for k in EQUATIONS), 'Actual baseline XY differs from HDF scaled history')
    datasets = same_mesh(native / 'step-baseline.cas.h5', native / 'step-final.cas.h5')
    final = read_fluent(run, setup)
    mapping = map_and_check(final, setup)
    baseline = dict(centres=final['centres'].copy(), topology='Exact same saved native mesh and IDs as final')
    materials = {}
    for label in ('baseline', 'final'):
        with h5py.File(native / ('step-' + label + '.dat.h5'), 'r') as handle:
            fields = handle['results/1/phase-1/cells']
            materials[label] = {}
            for key, expected in (('SV_DENSITY', setup['density_kg_m3']), ('SV_MU_LAM', setup['dynamic_viscosity_Pa_s'])):
                values = hdf_indexed(fields[key], setup['cells'])
                require(np.allclose(values, expected, rtol=1e-12, atol=1e-12), 'Saved material field differs: ' + key)
                materials[label][key] = [float(values.min()), float(values.max())]
            if label == 'baseline':
                baseline.update({key: hdf_indexed(fields[name], setup['cells']) for key, name in (('u', 'SV_U'), ('v', 'SV_V'), ('p', 'SV_P'))})
    map_and_check(baseline, setup)
    speed_scale = setup['inlet_velocity_m_s']
    difference = np.column_stack([final[k] - baseline[k] for k in ('u', 'v')])
    cell, component = np.unravel_index(np.abs(difference).argmax(), difference.shape)
    velocity_change = float(np.abs(difference).max() / speed_scale)
    _, _, x, _ = grid_coordinates(setup)
    left, right = [int(np.argmin(abs(x - location))) for location in (.8, 1.8)]
    drop = {name: float(field['grid']['p'][:, left].mean() - field['grid']['p'][:, right].mean()) for name, field in (('baseline', baseline), ('final', final))}
    absolute_pressure_change = abs(drop['final'] - drop['baseline'])
    pressure_scale = setup['density_kg_m3'] * speed_scale**2
    relative_meaningful = abs(drop['baseline']) > 1e-6 * pressure_scale
    pressure_change = absolute_pressure_change / abs(drop['baseline']) if relative_meaningful else None
    pressure_pass = (pressure_change if relative_meaningful else absolute_pressure_change / pressure_scale) <= .001
    mass_text = (native / 'mass-balance.txt').read_text(encoding='utf-8')
    mass = {k: float(re.search(r'^\s*' + k + r'\s+([-+0-9.eE]+)\s*$', mass_text, re.M)[1]) for k in ('inlet', 'outlet', 'Net')}
    require(mass['inlet'] > 0 and mass['outlet'] < 0, 'Unexpected boundary mass signs')
    relative_mass = abs(mass['Net']) / abs(mass['inlet'])
    transcript, launch = [p.read_text(encoding='utf-8', errors='replace') for p in private_logs]
    def number(pattern):
        return float(re.search(pattern + r'\s*([-+0-9.eE]+)', transcript)[1])
    quality = dict(minimum_volume_m3=number(r'minimum volume \(m3\):'),
                   maximum_volume_m3=number(r'maximum volume \(m3\):'),
                   total_volume_m3=number(r'total volume \(m3\):'),
                   minimum_orthogonal_quality=number(r'Minimum Orthogonal Quality\s*='),
                   maximum_aspect_ratio=number(r'Maximum Aspect Ratio\s*='),
                   same_saved_native_mesh_datasets=datasets,
                   full_native_mesh_check_completed=bool(re.search(r'Checking mesh\.+\s*Done\.', transcript)))
    processes = len(re.findall(r'^n\d+\*?\s+', launch, re.M))
    limit = int(re.search(r'Your license enables (\d+)-way parallel', launch)[1])
    markers = [int(v) for v in re.findall(r'!\s+(\d+)\s+solution is converged', transcript)]
    printed = re.findall(r'^\s*' + str(index[-1]) + r'\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+\d+:\d+:\d+\s+\d+\s*$', transcript, re.M)
    require(printed and all(float(v) == float(format(float(residual[k][-1]), '.5g')) for k, v in zip(EQUATIONS, printed[-1])), 'Final printed residuals differ from native HDF')
    loaded, saved = [json.loads((native / (n + '-settings.json')).read_text(encoding='utf-8')) for n in ('loaded', 'final')]
    actual = settings_values(saved)
    required_settings = (saved == loaded and actual['solver']['type'] == 'pressure-based'
        and actual['solver']['two_dim_space'] == 'planar' and actual['solver']['time'] == 'steady'
        and actual['viscous'] == 'laminar' and actual['material']['density']['value'] == 1.
        and actual['material']['viscosity']['value'] == .001 and actual['fluid_material'] == 'air'
        and actual['inlet']['velocity_specification_method'] == 'Components'
        and [v['value'] for v in actual['inlet']['velocity_components']] == [1., 0.]
        and actual['outlet_gauge_Pa'] == 0. and actual['walls']['wall_motion'] == 'Stationary Wall'
        and actual['walls']['shear_condition'] == 'No Slip' and actual['coupling'] == 'SIMPLE'
        and actual['methods']['gradient_scheme'] == 'least-square-cell-based'
        and actual['methods']['discretization_scheme'] == {'mom': 'second-order-upwind', 'pressure': 'second-order'}
        and actual['references']['length'] == .1 and actual['references']['area'] == .1
        and all(actual['residual_equations'][k]['check_convergence'] and actual['residual_equations'][k]['absolute_criteria'] == 1e-6 for k in EQUATIONS))
    checks = dict(student_live_startup_two_processes='student version' in launch.lower() and processes == 2,
        full_mesh_positive_orthogonal=quality['full_native_mesh_check_completed'] and quality['minimum_volume_m3'] > 0 and quality['minimum_orthogonal_quality'] >= .99999,
        native_cells_and_materials_match=True, declared_numerical_settings_match=required_settings,
        saved_endpoints_100_iterations_apart=True, identical_baseline_final_native_mesh=True,
        final_residuals_le_1e_6=all(v[-1] <= targets['residual_max'] for v in residual.values()),
        last100_residuals_le_1e_6=all(v[-100:].max() <= targets['residual_max'] for v in residual.values()),
        relative_mass_imbalance_le_1e_5=relative_mass <= targets['relative_mass_imbalance_max'],
        all_fluid_velocity_component_endpoint_change_le_0p001=velocity_change <= .001,
        downstream_pressure_drop_endpoint_change_le_0p001=pressure_pass)
    record = dict(solver='ANSYS Fluent Student 2026 R1', mesh_cells=setup['cells'],
        startup=dict(actual_compute_processes=processes, licensed_parallel_limit=limit), mesh_quality=quality,
        governing_settings=actual, saved_material_field_ranges=materials,
        reference_viscosity_note=dict(reference_value_Pa_s=actual['references']['viscosity'], material_value_Pa_s=.001,
            reference_matches_material=actual['references']['viscosity'] == .001,
            meaning='Display/reference setting remains default; solved viscosity and Re_h use the verified material .001 Pa.s'),
        residuals=dict(native_initial_convergence_marker=min(markers) if markers else None,
            baseline_iteration=int(first_index[-1]), final_iteration=int(index[-1]), full_history_records_per_equation=len(index),
            final_scaled={k: float(v[-1]) for k, v in residual.items()}, maximum_last100_scaled={k: float(v[-100:].max()) for k, v in residual.items()},
            definition='Native HDF data[:,0]/data[:,1]; final transcript agrees at five significant figures', cross_solver_norm_equivalence=False),
        mass_balance=dict(signed_report_kg_s=mass, relative_imbalance=relative_mass, definition='abs(reported Net)/abs(reported inlet)'),
        final100_saved_endpoint_stability=dict(iterations=[int(first_index[-1]), int(index[-1])],
            maximum_all_fluid_velocity_component_change_over_Uin=velocity_change,
            maximum_change_component=('u', 'v')[component], maximum_change_cell_centre_m=final['centres'][cell].tolist(),
            maximum_velocity_vector_change_over_Uin=float(np.linalg.norm(difference, axis=1).max() / speed_scale),
            downstream_pressure_drop_x_m=[float(x[left]), float(x[right])], cross_section_cell_count=setup['ny'],
            cross_section_mean_pressure_drop_Pa=drop, absolute_pressure_drop_change_Pa=absolute_pressure_change,
            relative_pressure_drop_change=pressure_change, relative_pressure_scale_meaningful=relative_meaningful,
            pressure_change_over_rho_Uin_squared=absolute_pressure_change / pressure_scale,
            definition='Equal cell areas: cross-section mean static pressure at nearest cell-centre columns; rho*Uin^2 is the absolute fallback scale',
            limitation='Only two saved field endpoints 100 iterations apart; no every-iteration field-stability claim'),
        native_field_mapping_and_separation=mapping,
        residual_export=dict(final_native_xy_bytes=(native / 'residuals.xy').stat().st_size,
            source='Complete saved HDF scaled residuals, float32-cast then six significant figures as native XY precision',
            checked_actual_baseline_xy_scalar_pairs=sum(len(r)-1 for r in xy.values()), csv_records=len(index)),
        check_results=checks, all_declared_checks_passed=all(checks.values()),
        benchmark_matching=False, mesh_independence_established=False, physical_validation=False)
    require(all(sha256(p) == digest for p, digest in before.items()), 'Source changed during analysis')
    out.mkdir(parents=True, exist_ok=True)
    (out / 'executed-fluent.py').write_bytes((native / 'executed-script.py').read_bytes())
    (out / 'executed-fluent-analyser.py').write_bytes(Path(__file__).read_bytes())
    recovered = residual_csv(native / 'residuals.xy', out / 'fluent-residuals.csv', native / 'step-final.dat.h5')
    with (native / 'iteration-history.csv').open(newline='', encoding='utf-8') as stream:
        existing = {int(r['iteration']): r for r in csv.DictReader(stream)}
    require(len(recovered) == len(existing) == len(index), 'CSV/full HDF history length mismatch')
    for i in index:
        require(all(recovered[int(i)][k] == float(existing[int(i)][k]) == float(format(float(np.float32(residual[k][i-1])), '.6g')) for k in EQUATIONS), 'Exported/native CSV residuals differ from HDF')
    write_json(out / 'fluent-checks.json', record)
    unchanged = all(sha256(p) == digest for p, digest in before.items())
    require(unchanged, 'Input changed during analysis/output')
    manifest = dict(native_and_sources_unchanged_during_analysis=unchanged,
        sources=[dict(path=str(p), sha256=before[p]) for p in paths],
        private_log_evidence_hashes=[dict(role=role, sha256=before[p]) for role, p in zip(('native numerical evidence', 'native startup/licence evidence'), private_logs)],
        scripts=[dict(path=str(p), sha256=before[p]) for p in scripts],
        outputs=[dict(path=str(out/n), sha256=sha256(out/n)) for n in names if n != 'fluent-source-manifest.json'],
        executed_driver_archive_matches_prelaunch_sha256=sha256(out/'executed-fluent.py') == provenance['source_sha256'],
        native_logs_published=False, solver_started_or_controls_modified=False)
    write_json(out / 'fluent-source-manifest.json', manifest)
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0 if all(checks.values()) else 2


if __name__ == '__main__':
    raise SystemExit(main())
