"""Bounded Re=100 cylinder/channel transient in one hidden two-process Fluent.

Run only when authorized. --end-time is absolute physical time (default 2 s;
maximum 8 s). --resume-case requires a saved case/data pair from this driver.
All generated native, report, snapshot and temporary files stay under --out on
E:. The practical p=0 pressure outlet is not asserted equivalent to exact DFG
traction conditions. A short startup pilot is not a periodic/validated result.
"""
from pathlib import Path
import argparse
import csv
import ctypes
from datetime import datetime
import hashlib
import json
import math
import os
import re
import shutil
import tempfile
import time
import traceback

DT = 0.005
SAMPLE_STEPS = 5
MAX_INNER = 20
EQUATIONS = ('continuity', 'x-velocity', 'y-velocity')
DEFAULT_EXECUTABLE = Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
KIND = 'cylinder-channel-Re100'
TARGETS = dict(dt_s=DT, initial_pilot_end_time_s=2.0, maximum_end_time_s=8.0,
               max_inner_iterations_per_time_step=MAX_INNER,
               residual_target=1e-5, relative_mass_imbalance_target=1e-4,
               sample_interval_steps=SAMPLE_STEPS, sample_interval_s=DT*SAMPLE_STEPS,
               temporal_order='second-order implicit',
               solver_processes=2, physical_validation=False,
               mesh_independence=False, time_step_independence=False,
               acceptance='Postprocessing required; solver stopping is not a passed numerical check')
REFERENCE = dict(density=1.0, viscosity=0.001, velocity=1.0,
                 length=0.1, area=0.1, depth=1.0)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + '\n', encoding='utf-8')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def step_index(value):
    if not math.isfinite(value) or value < 0:
        raise ValueError('Physical time must be finite and nonnegative')
    step = round(value / DT)
    if not math.isclose(value, step*DT, rel_tol=0.0, abs_tol=1e-8):
        raise ValueError('Physical time must lie on the fixed 0.005 s time grid')
    return step


def remaining_steps(end_time, start_time):
    end, start = step_index(end_time), step_index(start_time)
    if end_time > 8.0 or end <= start:
        raise ValueError('End time must exceed saved physical time and be at most 8 s')
    return end-start


def force_row(text):
    rows = re.findall(r'^\s*Net\s+((?:[-+0-9.eE]+\s+){5}[-+0-9.eE]+)\s*$', text, re.M)
    if not rows:
        raise ValueError('Native force report has no six-column scalar Net row')
    values = [float(value) for value in rows[-1].split()]
    if len(values) != 6 or not all(map(math.isfinite, values)):
        raise ValueError('Native force report has nonfinite or missing values')
    return values


def mesh_quality_numbers(text):
    number = r'([-+0-9.eE]+)'
    orth = re.findall(r'Minimum Orthogonal Quality\s*=\s*'+number, text, re.I)
    volume = re.findall(r'minimum volume\s*\([^)]*\)\s*:\s*'+number, text, re.I)
    if not orth or not volume:
        raise ValueError('Native mesh report lacks minimum orthogonal quality or volume')
    orth, volume = float(orth[-1]), float(volume[-1])
    if not math.isfinite(orth) or not math.isfinite(volume) or orth < 0.2 or volume <= 0:
        raise ValueError(f'Native mesh gate failed: min orthogonal quality={orth}, min volume={volume}')
    return orth, volume


def available_memory():
    class MemoryStatus(ctypes.Structure):
        _fields_ = [('length', ctypes.c_ulong), ('load', ctypes.c_ulong)] + [
            (name, ctypes.c_ulonglong) for name in
            ('total_phys', 'available_phys', 'total_page', 'available_page',
             'total_virtual', 'available_virtual', 'available_extended')]
    memory = MemoryStatus()
    memory.length = ctypes.sizeof(memory)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(memory)):
        raise ctypes.WinError()
    return memory.available_phys


def case_data_path(case):
    if not case.name.endswith('.cas.h5'):
        raise ValueError('Resume requires a .cas.h5 file')
    return case.with_name(case.name[:-7] + '.dat.h5')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--end-time', type=float, default=2.0)
    parser.add_argument('--resume-case', type=Path)
    args = parser.parse_args()
    remaining_steps(args.end_time, 0.0)
    mesh, out = args.mesh.resolve(), args.out.resolve()
    repo = Path(__file__).resolve().parents[1]
    executable = Path(os.environ.get('FLUENT_EXECUTABLE', str(DEFAULT_EXECUTABLE))).resolve()
    resume = args.resume_case.resolve() if args.resume_case else None
    if out.exists() or out.is_relative_to(repo) or out.drive.upper() != 'E:':
        raise ValueError('Choose a fresh output directory on E: outside the repository')
    if mesh.drive.upper() != 'E:' or not mesh.is_file():
        raise ValueError('The native mesh must be an existing file on E:')
    if executable != DEFAULT_EXECUTABLE.resolve() or not executable.is_file():
        raise ValueError('Use the installed D:/ANSYS Inc/ANSYS Student/v261 Fluent executable')
    free_disk = shutil.disk_usage(out.anchor).free
    if free_disk < 10*1024**3:
        raise RuntimeError('Less than the declared 10 GiB E: free-space floor')
    free_memory = available_memory()
    if free_memory < 4*1024**3:
        raise RuntimeError('Less than the declared 4 GiB available-RAM floor')
    mesh_setup_path = mesh.parent/'setup.json'
    profile = mesh.parent/'parabolic-inlet.prof'
    inlet_points = mesh.parent/'inlet-face-profile.json'
    mesh_setup = json.loads(mesh_setup_path.read_text(encoding='utf-8'))
    mesh_hash = sha256(mesh)
    if mesh_setup.get('mesh_sha256') != mesh_hash:
        raise ValueError('Mesh differs from its setup.json SHA-256')
    if not profile.is_file() or not inlet_points.is_file():
        raise FileNotFoundError('Exact inlet-face parabolic profile and coordinates are required')
    points = json.loads(inlet_points.read_text(encoding='utf-8'))
    if not points:
        raise ValueError('Empty inlet-face profile')
    for point in points:
        y = point['y']
        expected = 6*y*(0.41-y)/0.41**2
        if (abs(point['x']) > 1e-12 or not 0 <= y <= 0.41 or
                not math.isclose(point['u'], expected, rel_tol=0.0, abs_tol=1e-10)):
            raise ValueError('Inlet profile does not match the declared parabola')
    resume_provenance = None
    if resume:
        if resume.drive.upper() != 'E:' or not resume.is_file() or not case_data_path(resume).is_file():
            raise FileNotFoundError('An existing resume case/data pair on E: is required')
        previous = json.loads((resume.parent/'geometry-setup.json').read_text(encoding='utf-8'))
        resume_provenance = json.loads((resume.parent/'code-provenance.json').read_text(encoding='utf-8'))
        if (previous.get('kind') != KIND or previous.get('references') != REFERENCE or
                resume_provenance.get('mesh_sha256') != mesh_hash):
            raise ValueError('Resume metadata does not match this mesh and physical setup')
    out.mkdir(parents=True)
    for name in ('snapshots', 'force-reports', 'temp', 'examples'):
        (out/name).mkdir()
    for variable in ('TEMP', 'TMP', 'TMPDIR'):
        os.environ[variable] = str(out/'temp')
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    tempfile.tempdir = str(out/'temp')
    os.chdir(out)
    source = Path(__file__).read_bytes()
    (out/'executed-script.py').write_bytes(source)
    shutil.copyfile(mesh_setup_path, out/'mesh-setup.json')
    shutil.copyfile(profile, out/'parabolic-inlet.prof')
    shutil.copyfile(inlet_points, out/'inlet-face-profile.json')
    provenance = dict(source_sha256=hashlib.sha256(source).hexdigest(),
                      mesh=str(mesh), mesh_sha256=mesh_hash,
                      mesh_setup_sha256=sha256(mesh_setup_path),
                      profile_sha256=sha256(profile), inlet_points_sha256=sha256(inlet_points),
                      available_ram_before_launch_bytes=free_memory,
                      free_disk_before_launch_bytes=free_disk,
                      snapshot_before_solver_launch=True,
                      resume_case=str(resume) if resume else None,
                      resume_case_sha256=sha256(resume) if resume else None,
                      resume_data_sha256=sha256(case_data_path(resume)) if resume else None,
                      parent_provenance=resume_provenance)
    write_json(out/'code-provenance.json', provenance)
    write_json(out/'declared-targets.json', TARGETS)
    geometry = dict(kind=KIND, mesh_setup=mesh_setup, references=REFERENCE,
                    rho_kg_m3=1.0, mu_Pa_s=0.001, nu_m2_s=0.001, Re_D=100,
                    inlet='u(y)=6*y*(0.41-y)/0.41^2 m/s, v=0; exact inlet-face centers',
                    inlet_discrete_flow_m3_s=sum(p['u']*p['length'] for p in points),
                    inlet_exact_flow_m3_s=0.41, cylinder_center_m=[0.2, 0.2],
                    cylinder_radius_m=0.05, channel_bounds_m=[0.0, 2.2, 0.0, 0.41],
                    walls='stationary no-slip channel walls and cylinder',
                    outlet='practical zero-gauge-pressure outlet; no exact DFG traction-equivalence claim',
                    initialization='standard zero interior velocity and gauge pressure; inlet parabola applied immediately' if not resume else 'saved transient state; no reinitialization',
                    benchmark_matching=False, physical_validation=False)
    write_json(out/'geometry-setup.json', geometry)
    write_json(out/'requested-settings.json', dict(
        targets=TARGETS, geometry=geometry, end_time_s=args.end_time,
        coupling='PISO', momentum_scheme='second-order-upwind', pressure_scheme='second-order',
        temporal_scheme='unsteady-2nd-order', ui_mode='hidden_gui', processor_count=2,
        executable=str(executable), out=str(out), temp=str(out/'temp')))
    events, export_errors, snapshots, history = [], [], [], []
    session, loaded, final_saved = None, False, False
    actual_time = None

    def update(stage, **values):
        record = dict(stage=stage, time=datetime.now().astimezone().isoformat(), **values)
        events.append(record)
        write_json(out/'history.json', events)
        write_json(out/'status.json', record)
        print(json.dumps(record, default=str), flush=True)

    def export(label, operation):
        try:
            operation()
            return True
        except Exception:
            error = dict(export=label, error=traceback.format_exc())
            export_errors.append(error)
            write_json(out/'export-errors.json', export_errors)
            print(json.dumps(error), flush=True)
            return False

    try:
        update('importing_pyfluent')
        import ansys.fluent.core as pf
        pf.config.examples_path = str(out/'examples')
        update('launching_single_session', processor_count=2, ui_mode='hidden_gui',
               dimension=2, precision='double', end_time_s=args.end_time)
        session = pf.launch_fluent(fluent_path=str(executable), dimension=2,
                                  precision='double', processor_count=2,
                                  ui_mode='hidden_gui', cwd=str(out),
                                  start_timeout=180, cleanup_on_exit=True)
        update('solver_connected', version=str(session.get_fluent_version()))
        session.transcript.start(file_name=str(out/'solver-transcript.txt'))
        st = session.settings
        if resume:
            st.file.read_case_data(file_name=str(resume))
        else:
            st.file.read_mesh(file_name=str(mesh))
        loaded = True
        st.setup.general.solver.type = 'pressure-based'
        st.setup.general.solver.two_dim_space = 'planar'
        st.setup.general.solver.time = 'transient'
        st.setup.models.viscous.model = 'laminar'
        st.setup.materials.fluid['air'].density.value = 1.0
        st.setup.materials.fluid['air'].viscosity.value = 0.001
        st.setup.cell_zone_conditions.fluid['fluid'].general.material = 'air'
        # Import on fresh mesh. Resume case embeds the profile, avoiding a duplicate
        # named profile and preserving the prior state; both are checked below.
        if not resume:
            st.file.read_profile(file_name=str(out/'parabolic-inlet.prof'))
        profile_names = st.setup.profiles.query_list_profiles()
        profile_fields = st.setup.profiles.query_list_profile_fields(profile_name='parabolic-inlet')
        write_json(out/'native-profile-list.json', dict(names=profile_names, fields=profile_fields))
        inlet = st.setup.boundary_conditions.velocity_inlet['inlet'].momentum
        inlet.velocity_specification_method = 'Components'
        inlet.velocity_components = [
            dict(option='profile', profile_name='parabolic-inlet', field_name='u'),
            dict(option='value', value=0.0)]
        components = inlet.velocity_components.get_state()
        write_json(out/'native-profile-list.json', dict(names=profile_names, fields=profile_fields,
                   assigned_velocity_components=components))
        # The typed backend readback is decisive even if a listing query prints
        # its results and returns None in a particular native build.
        if (not isinstance(components, list) or len(components) < 2 or
                components[0].get('option') != 'profile' or
                components[0].get('profile_name') != 'parabolic-inlet' or
                components[0].get('field_name') != 'u'):
            raise ValueError('Native inlet readback does not match parabolic-inlet/u')
        st.setup.boundary_conditions.pressure_outlet['outlet'].momentum.gauge_pressure.value = 0.0
        for name in ('walls', 'cylinder'):
            wall = st.setup.boundary_conditions.wall[name].momentum
            wall.wall_motion = 'Stationary Wall'
            wall.shear_condition = 'No Slip'
        for key, value in REFERENCE.items():
            setattr(st.setup.reference_values, key, value)
        methods = st.solution.methods
        methods.p_v_coupling.flow_scheme = 'PISO'
        methods.spatial_discretization.gradient_scheme = 'least-square-cell-based'
        methods.spatial_discretization.discretization_scheme = {
            'mom': 'second-order-upwind', 'pressure': 'second-order'}
        methods.transient_formulation = 'unsteady-2nd-order'
        run = st.solution.run_calculation
        run.transient_controls.type = 'Fixed'
        run.transient_controls.method = 'User-Specified'
        run.transient_controls.duration_specification_method = 'incremental-time-steps'
        run.parameters.time_step_size = DT
        run.parameters.max_iter_per_time_step = MAX_INNER
        for name in EQUATIONS:
            equation = st.solution.monitor.residual.equations[name]
            equation.absolute_criteria = TARGETS['residual_target']
            equation.check_convergence = True
        check_result = st.mesh.check()
        quality_result = st.mesh.quality()
        write_json(out/'native-mesh-reports.json', dict(check=check_result, quality=quality_result,
                    transcript='solver-transcript.txt', note='Native command output is retained in transcript'))
        quality_text = str(check_result) + '\n' + str(quality_result)
        for attempt in range(10):
            transcript = out/'solver-transcript.txt'
            combined = quality_text + '\n' + (transcript.read_text(encoding='utf-8', errors='replace')
                                               if transcript.is_file() else '')
            try:
                minimum_orthogonal_quality, minimum_volume = mesh_quality_numbers(combined)
                break
            except ValueError:
                if attempt == 9:
                    raise
                time.sleep(0.1)
        write_json(out/'geometry-quality.json', dict(
            minimum_orthogonal_quality=minimum_orthogonal_quality,
            minimum_cell_volume_m3=minimum_volume,
            orthogonal_quality_floor=0.2, positive_volume_required=True,
            native_numeric_gate_pass=True, checked_before_flow=True,
            source='native mesh.check/quality and solver-transcript.txt'))

        def settings_snapshot():
            return dict(general=st.setup.general.get_state(), viscous=st.setup.models.viscous.get_state(),
                        material=st.setup.materials.fluid['air'].get_state(),
                        cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),
                        inlet=st.setup.boundary_conditions.velocity_inlet['inlet'].get_state(),
                        outlet=st.setup.boundary_conditions.pressure_outlet['outlet'].get_state(),
                        walls={name: st.setup.boundary_conditions.wall[name].get_state()
                               for name in ('walls', 'cylinder')},
                        references=st.setup.reference_values.get_state(), methods=methods.get_state(),
                        run_calculation=run.get_state(), residual=st.solution.monitor.residual.get_state())

        write_json(out/'loaded-settings.json', settings_snapshot())
        if not resume:
            initialization = st.solution.initialization
            initialization.initialization_type = 'standard'
            for name in ('x-velocity', 'y-velocity', 'pressure'):
                initialization.defaults[name] = 0.0
            initialization.standard_initialize()
        actual_time = float(run.transient_controls.flow_time())
        remaining = remaining_steps(args.end_time, actual_time)
        start_time = actual_time
        start_step = step_index(actual_time)
        st.file.write_case(file_name=str(out/'vortex-setup.cas.h5'))
        write_json(out/'initialized-settings.json', settings_snapshot())

        def save_snapshot(step, physical_time):
            path = out/'snapshots'/f'time-{step:06d}.dat.h5'
            st.file.write_data(file_name=str(path))
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError('Native snapshot is missing or empty: ' + str(path))
            snapshots.append(dict(time=physical_time, time_s=physical_time, step=step,
                                  data=str(path.relative_to(out)), bytes=path.stat().st_size))
            write_json(out/'snapshots.json', snapshots)

        def sample_forces(step, physical_time):
            row = dict(step=step, time_s=physical_time)
            for label, direction in (('drag', [1, 0, 0]), ('lift', [0, 1, 0])):
                path = out/'force-reports'/f'{label}-{step:06d}.txt'
                st.results.report.forces(option='forces', wall_zones=['cylinder'],
                    direction_vector=direction, write_to_file=True, file_name=str(path))
                values = force_row(path.read_text(encoding='utf-8', errors='replace'))
                row[label+'_force_N'] = values[2]
                row[label+'_pressure_N'], row[label+'_viscous_N'] = values[:2]
                row['Cd' if label == 'drag' else 'Cl'] = values[5]
                # Native coefficients use rho=1,U=1,A=.1,depth=1. Report text is
                # rounded; compare to a tolerance appropriate to its precision.
                if not math.isclose(values[5], values[2]/0.05, rel_tol=2e-4, abs_tol=1e-8):
                    raise ValueError('Native force coefficient normalization differs from declared references')
            history.append(row)
            with (out/'force-history.csv').open('w', newline='', encoding='utf-8') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(history[0]))
                writer.writeheader()
                writer.writerows(history)

        if resume:
            # Carry prior pilot records into the extension manifest/history; original
            # pilot files remain untouched. Absolute ancestor paths are intentional.
            old_manifest = resume.parent/'snapshots.json'
            if old_manifest.is_file():
                for item in json.loads(old_manifest.read_text(encoding='utf-8')):
                    if float(item['time']) < actual_time-1e-8:
                        item = dict(item)
                        old_path = Path(item['data'])
                        item['data'] = str(old_path if old_path.is_absolute() else resume.parent/old_path)
                        snapshots.append(item)
            old_history = resume.parent/'force-history.csv'
            if old_history.is_file():
                with old_history.open(newline='', encoding='utf-8') as stream:
                    for row in csv.DictReader(stream):
                        if float(row['time_s']) < actual_time-1e-8:
                            history.append({key: int(value) if key == 'step' else float(value)
                                            for key, value in row.items()})
        save_snapshot(start_step, actual_time)
        export('forces_initial', lambda: sample_forces(start_step, actual_time))
        update('advancing', start_time_s=start_time, end_time_s=args.end_time,
               additional_steps=remaining, max_inner_iterations=MAX_INNER,
               sampling='native forces and full native data every five time steps (.025 s)')
        completed = 0
        while completed < remaining:
            if shutil.disk_usage(out.anchor).free < 10*1024**3:
                raise RuntimeError('E: free space fell below 10 GiB during the bounded run')
            count = min(SAMPLE_STEPS, remaining-completed)
            run.dual_time_iterate(time_step_count=count, max_iter_per_step=MAX_INNER)
            completed += count
            actual_time = float(run.transient_controls.flow_time())
            expected_time = start_time + completed*DT
            if not math.isclose(actual_time, expected_time, rel_tol=0.0, abs_tol=1e-8):
                raise RuntimeError(f'Fluent time {actual_time} differs from requested {expected_time}')
            step = start_step+completed
            save_snapshot(step, actual_time)
            export('forces_'+str(step), lambda: sample_forces(step, actual_time))
            update('advancing', time_s=actual_time, step=step, completed_steps=completed,
                   requested_additional_steps=remaining, snapshots=len(snapshots), force_samples=len(history))
        st.file.write_case_data(file_name=str(out/'vortex-final.cas.h5'))
        final_saved = True
        update('final_saved', time_s=actual_time, case='vortex-final.cas.h5', data='vortex-final.dat.h5')
        export('final_settings', lambda: write_json(out/'final-settings.json', settings_snapshot()))
        export('residuals', lambda: st.solution.monitor.residual.write(file_name=str(out/'residuals.xy')))
        mass_saved = export('mass_balance', lambda: st.results.report.fluxes.mass_flow(
            zones=['inlet', 'outlet'], write_to_file=True, file_name=str(out/'mass-balance.txt')))
        export('final_mesh_quality', lambda: st.mesh.quality())

        def picture(name, field):
            graphics = st.results.graphics
            graphics.picture.use_window_resolution = False
            graphics.picture.x_resolution = 1800
            graphics.picture.y_resolution = 600
            graphics.picture.invert_background = False
            graphics.views.camera.projection = 'orthographic'
            if name not in graphics.contour.get_object_names():
                graphics.contour.create(name=name)
            graphics.contour[name].field = field
            graphics.contour[name].surfaces_list = []
            graphics.contour[name].display()
            graphics.views.camera.position = [1.1, 0.205, 10.0]
            graphics.views.camera.target = [1.1, 0.205, 0.0]
            graphics.views.camera.field = [2.42, 0.807]
            graphics.picture.save_picture(file_name=str(out/(name+'.png')))

        for name, field in (('velocity', 'velocity-magnitude'), ('pressure', 'pressure'),
                            ('vorticity', 'z-vorticity')):
            export(name+'_png', lambda name=name, field=field: picture(name, field))
        export('transcript_stop', session.transcript.stop)
        summary = dict(kind=KIND, start_time_s=start_time, end_time_s=actual_time,
                       dt_s=DT, additional_time_steps=remaining, max_inner_iterations=MAX_INNER,
                       native_final_saved=True, snapshot_count=len(snapshots),
                       snapshots_manifest='snapshots.json', force_samples=len(history),
                       force_sampling_interval_s=SAMPLE_STEPS*DT,
                       force_history='force-history.csv' if history else None,
                       mass_report_saved=mass_saved, targets=TARGETS,
                       export_errors=export_errors, physical_validation=False,
                       numerical_checks='Pending native-field, residual, mass-balance and transient-history analysis',
                       interpretation='Startup pilot; no established periodic shedding, frequency or benchmark-accuracy claim')
        write_json(out/'summary.json', summary)
        update('completed_with_export_errors' if export_errors else 'completed', **summary)
    except Exception:
        failure = dict(error=traceback.format_exc(), native_final_saved=final_saved,
                       last_time_s=actual_time, snapshots=len(snapshots), export_errors=export_errors)
        write_json(out/'summary.json', dict(stage='failed', **failure))
        update('failed', **failure)
        raise
    finally:
        if session is not None:
            if loaded and not final_saved:
                try:
                    session.settings.file.write_case_data(file_name=str(out/'recovery.cas.h5'))
                    write_json(out/'recovery.json', dict(case='recovery.cas.h5', data='recovery.dat.h5',
                                                        last_observed_time_s=actual_time))
                except Exception:
                    write_json(out/'recovery-error.json', dict(error=traceback.format_exc()))
            try:
                session.exit()
            except Exception:
                write_json(out/'session-exit-error.json', dict(error=traceback.format_exc()))


if __name__ == '__main__':
    main()

