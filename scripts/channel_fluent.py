"""Fresh, bounded Fluent calculation for the 2-D planar starter channel.

Run with --mesh PATH --out FRESH_DIRECTORY. Native baseline/final case-data
pairs support a post-run outlet-profile and pressure-drop stability comparison.
This script starts one two-process Fluent session only when explicitly executed.
"""
from pathlib import Path
import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
import re
import shutil
import tempfile
import traceback


EQUATIONS = ('continuity', 'x-velocity', 'y-velocity')
TARGETS = dict(initial_iteration_cap=800, stability_iterations=50,
               residual_max=1e-6, relative_mass_imbalance_max=1e-5,
               downstream_profile_relative_rms_max=0.02,
               developed_pressure_gradient_relative_error_max=0.03,
               final50_profile_and_pressure_drop_relative_change_max=0.001,
               reference_pressure_gradient_Pa_per_m=12.0,
               validation_scope='Educational numerical check; no physical validation')


def write_json(path, values):
    path.write_text(json.dumps(values, indent=2, default=str) + '\n', encoding='utf-8')


def residual_rows(source, native_data=None):
    """Read native XY pairs, or recover an empty export from saved HDF history."""
    rows, label = {}, None
    for line in source.read_text(encoding='utf-8', errors='replace').splitlines():
        match = re.search(r'xy/key/label\s+"([^"]+)"', line)
        if match:
            label = match.group(1)
            continue
        values = line.strip().split()
        if label and len(values) == 2:
            try:
                iteration, value = int(float(values[0])), float(values[1])
            except ValueError:
                continue
            rows.setdefault(iteration, {'iteration': iteration})[label] = value
    if not rows and native_data is not None:
        import h5py
        from struct import pack, unpack
        with h5py.File(native_data, 'r') as handle:
            for label in EQUATIONS:
                group = handle['results/residuals/phase-1/' + label]
                iterations, data = group['iterations'][:], group['data'][:]
                if len(iterations) != len(data):
                    raise ValueError('Native HDF residual lengths differ')
                for iteration, values in zip(iterations, data):
                    if values[1] <= 0 or iteration != int(iteration):
                        raise ValueError('Invalid native residual scaling or iteration')
                    # Observed native XY graph precision: float32-cast scaled
                    # ratio, then six significant figures. Full values stay in HDF.
                    plotted = unpack('f', pack('f', float(values[0] / values[1])))[0]
                    value = float(format(plotted, '.6g'))
                    index = int(iteration)
                    rows.setdefault(index, {'iteration': index})[label] = value
    if not rows:
        raise ValueError('Empty or unsupported residual XY; native file retained')
    return rows


def residual_csv(source, destination, native_data=None):
    """Retain native inputs and export their per-iteration scaled residuals."""
    rows = residual_rows(source, native_data)
    fields = ['iteration'] + sorted({key for row in rows.values() for key in row} - {'iteration'})
    with destination.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows[key] for key in sorted(rows))
    return rows

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    mesh, out = args.mesh.resolve(), args.out.resolve()
    repo = Path(__file__).resolve().parents[1]
    executable = Path(os.environ['FLUENT_EXECUTABLE'])
    if not mesh.is_file():
        raise FileNotFoundError(mesh)
    if out.exists() or out.is_relative_to(repo) or out.drive.upper() != 'E:':
        raise ValueError('Choose a fresh output directory on E: outside the repository')
    if executable.drive.upper() != 'D:' or not executable.is_file():
        raise ValueError('FLUENT_EXECUTABLE must identify the installed Fluent executable on D:')
    if shutil.disk_usage(out.anchor).free < 10 * 1024**3:
        raise RuntimeError('Less than the declared 10 GiB free-space floor on E:')
    out.mkdir(parents=True)
    temp = Path('E:/CFD-Learning-Lab/temp/channel-fluent-v001')
    temp.mkdir(parents=True, exist_ok=True)
    os.environ['TEMP'] = os.environ['TMP'] = os.environ['TMPDIR'] = str(temp)
    tempfile.tempdir = str(temp)
    source = Path(__file__).read_bytes()
    (out / 'executed-script.py').write_bytes(source)
    write_json(out / 'code-provenance.json', dict(
        source_sha256=hashlib.sha256(source).hexdigest(), mesh=str(mesh),
        mesh_sha256=hashlib.sha256(mesh.read_bytes()).hexdigest(),
        snapshot_before_solver_launch=True))
    write_json(out / 'declared-targets.json', TARGETS)
    geometry = dict(length=1.0, height=0.1, nx=240, ny=32, cells=7680,
                    rho=1.0, mu=0.01, Umean=1.0, depth=1.0,
                    inlet='uniform velocity', outlet_gauge_pressure=0.0,
                    walls='stationary no-slip', Re_height=10, Re_hydraulic_diameter=20)
    write_json(out / 'geometry-setup.json', geometry)
    events, export_errors = [], []

    def update(stage, **values):
        event = dict(stage=stage, time=datetime.now().astimezone().isoformat(), **values)
        events.append(event)
        write_json(out / 'history.json', events)
        write_json(out / 'status.json', event)
        print(json.dumps(event, default=str), flush=True)

    def export(label, operation):
        # Case/data are saved first, so an optional report/graphics error cannot
        # discard the calculation or prevent the other exports being attempted.
        try:
            operation()
        except Exception:
            error = dict(export=label, error=traceback.format_exc())
            export_errors.append(error)
            write_json(out / 'export-errors.json', export_errors)
            print(json.dumps(error), flush=True)

    session, loaded, final_saved = None, False, False
    try:
        update('importing_pyfluent')
        import ansys.fluent.core as pf
        pf.config.examples_path = str(out / 'examples')
        update('launching_single_session', processor_count=2, ui_mode='hidden_gui',
               dimension=2, precision='double', initial_iteration_cap=800)
        session = pf.launch_fluent(fluent_path=str(executable), dimension=2,
                                  precision='double', processor_count=2,
                                  ui_mode='hidden_gui', cwd=str(out),
                                  start_timeout=180, cleanup_on_exit=True)
        update('solver_connected', version=str(session.get_fluent_version()))
        session.transcript.start(file_name=str(out / 'solver-transcript.txt'))
        st = session.settings
        st.file.read_mesh(file_name=str(mesh))
        loaded = True
        st.setup.general.solver.type = 'pressure-based'
        st.setup.general.solver.two_dim_space = 'planar'
        st.setup.general.solver.time = 'steady'
        st.setup.models.viscous.model = 'laminar'
        st.setup.materials.fluid['air'].density.value = 1.0
        st.setup.materials.fluid['air'].viscosity.value = 0.01
        st.setup.cell_zone_conditions.fluid['fluid'].general.material = 'air'
        inlet = st.setup.boundary_conditions.velocity_inlet['inlet'].momentum
        inlet.velocity_specification_method = 'Components'
        inlet.velocity_components = [{'value': 1.0}, {'value': 0.0}]
        st.setup.boundary_conditions.pressure_outlet['outlet'].momentum.gauge_pressure.value = 0.0
        wall = st.setup.boundary_conditions.wall['walls'].momentum
        wall.wall_motion = 'Stationary Wall'
        wall.shear_condition = 'No Slip'
        for key, value in dict(density=1.0, velocity=1.0, length=0.1, area=0.1).items():
            setattr(st.setup.reference_values, key, value)
        st.solution.methods.p_v_coupling.flow_scheme = 'SIMPLE'
        st.solution.methods.spatial_discretization.gradient_scheme = 'least-square-cell-based'
        st.solution.methods.spatial_discretization.discretization_scheme = {
            'mom': 'second-order-upwind', 'pressure': 'second-order'}
        for name in EQUATIONS:
            equation = st.solution.monitor.residual.equations[name]
            equation.absolute_criteria = TARGETS['residual_max']
            equation.check_convergence = True
        st.mesh.check()
        st.mesh.quality()

        def settings_snapshot():
            return dict(general=st.setup.general.get_state(),
                        viscous=st.setup.models.viscous.get_state(),
                        material=st.setup.materials.fluid['air'].get_state(),
                        cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),
                        inlet=st.setup.boundary_conditions.velocity_inlet['inlet'].get_state(),
                        outlet=st.setup.boundary_conditions.pressure_outlet['outlet'].get_state(),
                        walls=st.setup.boundary_conditions.wall['walls'].get_state(),
                        references=st.setup.reference_values.get_state(),
                        methods=st.solution.methods.get_state(),
                        residual=st.solution.monitor.residual.get_state())

        write_json(out / 'loaded-settings.json', settings_snapshot())
        st.file.write_case(file_name=str(out / 'channel-setup.cas.h5'))
        update('initial_convergence', iteration_cap=800, residual_target=1e-6)
        st.solution.initialization.hybrid_initialize()
        st.solution.run_calculation.iterate(iter_count=800)
        st.file.write_case_data(file_name=str(out / 'channel-baseline.cas.h5'))
        update('baseline_saved', case='channel-baseline.cas.h5',
               data='channel-baseline.dat.h5', stability_additional_iteration=0)
        export('baseline_residuals', lambda: st.solution.monitor.residual.write(
            file_name=str(out / 'residuals-baseline.xy')))
        # The fixed 50-iteration window measures stationarity; restore all three
        # convergence checks before saving the final reusable case.
        for name in EQUATIONS:
            st.solution.monitor.residual.equations[name].check_convergence = False
        update('stability_window', additional_iterations=50,
               comparison='Saved baseline/final outlet profile and pressure drop')
        st.solution.run_calculation.iterate(iter_count=50)
        for name in EQUATIONS:
            st.solution.monitor.residual.equations[name].check_convergence = True
        st.file.write_case_data(file_name=str(out / 'channel-final.cas.h5'))
        final_saved = True
        update('final_saved', case='channel-final.cas.h5', data='channel-final.dat.h5',
               stability_additional_iteration=50)
        export('final_settings', lambda: write_json(out / 'final-settings.json', settings_snapshot()))
        export('residuals', lambda: st.solution.monitor.residual.write(
            file_name=str(out / 'residuals.xy')))
        export('iteration_history_csv', lambda: residual_csv(
            out / 'residuals.xy', out / 'iteration-history.csv', out / 'channel-final.dat.h5'))
        export('mass_balance', lambda: st.results.report.fluxes.mass_flow(
            zones=['inlet', 'outlet'], write_to_file=True,
            file_name=str(out / 'mass-balance.txt')))
        # Installed v261 API: report_of and surface_names are native arguments.
        # These reports are a convenience; saved HDF fields remain the source
        # for downstream pressure-gradient and profile checks.
        for zone in ('inlet', 'outlet'):
            export(zone + '_pressure_average', lambda zone=zone:
                   st.results.report.surface_integrals.area_weighted_avg(
                       surface_names=[zone], report_of='pressure',
                       write_to_file=True, file_name=str(out / f'{zone}-pressure-average.txt')))

        def picture(name, field):
            graphics = st.results.graphics
            graphics.picture.use_window_resolution = False
            graphics.picture.x_resolution = 1600
            graphics.picture.y_resolution = 600
            graphics.picture.invert_background = False
            graphics.views.camera.projection = 'orthographic'
            if name not in graphics.contour.get_object_names():
                graphics.contour.create(name=name)
            graphics.contour[name].field = field
            graphics.contour[name].surfaces_list = []
            graphics.contour[name].display()
            graphics.views.camera.position = [0.5, 0.05, 10.0]
            graphics.views.camera.target = [0.5, 0.05, 0.0]
            graphics.views.camera.field = [1.15, 0.43125]
            graphics.picture.save_picture(file_name=str(out / f'{name}.png'))

        export('velocity_png', lambda: picture('velocity', 'velocity-magnitude'))
        export('pressure_png', lambda: picture('pressure', 'pressure'))
        export('transcript_stop', session.transcript.stop)
        summary = dict(kind='Fresh 2-D planar laminar channel; educational numerical calculation',
                       native_case_saved=True, baseline_case='channel-baseline.cas.h5',
                       final_case='channel-final.cas.h5', initial_iteration_cap=800,
                       stability_additional_iterations=50, convergence_checks_restored=True,
                       targets=TARGETS, export_errors=export_errors,
                       numerical_checks='Pending postprocessing of residuals, mass balance and saved baseline/final fields',
                       stability_fields=['channel-baseline.dat.h5', 'channel-final.dat.h5'])
        write_json(out / 'summary.json', summary)
        update('completed_with_export_errors' if export_errors else 'completed', **summary)
    except Exception:
        failure = dict(error=traceback.format_exc(), native_final_saved=final_saved,
                       export_errors=export_errors, targets=TARGETS)
        write_json(out / 'summary.json', dict(stage='failed', **failure))
        update('failed', **failure)
        raise
    finally:
        if session is not None:
            if loaded and not final_saved:
                try:
                    session.settings.file.write_case_data(file_name=str(out / 'recovery.cas.h5'))
                    write_json(out / 'recovery.json', dict(case='recovery.cas.h5',
                                                           data='recovery.dat.h5'))
                except Exception:
                    write_json(out / 'recovery-error.json', dict(error=traceback.format_exc()))
            try:
                session.exit()
            except Exception:
                write_json(out / 'session-exit-error.json', dict(error=traceback.format_exc()))


if __name__ == '__main__':
    main()