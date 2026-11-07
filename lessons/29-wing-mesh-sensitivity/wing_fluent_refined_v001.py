"""One bounded hidden two-process steady 3D finite-wing Fluent calculation.

--mesh points to the approved shared wing.msh; --out must be fresh on E:.
Maximum1500 fixed actual iterations; numerical targets checked afterwards. This computes steady fields. Tracers/orbits are steady-field
visualizations, never computed transient wing evolution. No automatic second run.
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

EQUATIONS=('continuity','x-velocity','y-velocity','z-velocity')
MAX_ITERATIONS=1500
REFERENCE=dict(density=1.,viscosity=.0025,velocity=1.,length=1.,area=4.)
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
TARGETS=dict(maximum_total_iterations=MAX_ITERATIONS,checkpoints=[300,600,900,1200,1500],
    residual_target=1e-5,relative_mass_imbalance_target=1e-4,minimum_orthogonal_quality=.10,
    force_stability_samples=100,force_stability_relative_range_max=.001,
    force_stability_definition='(max-min)/max(abs(mean),1e-3) for actual final100 consecutive iteration Cd and Cl',
    observation_after_convergence_maximum_iterations=100,solver_processes=2,
    mesh_independence=False,physical_validation=False,steady_visualization_only=True)


def write_json(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def available_memory():
    class MemoryStatus(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[
            (key,ctypes.c_ulonglong) for key in ('total_phys','available_phys','total_page','available_page',
                'total_virtual','available_virtual','available_extended')]
    status=MemoryStatus();status.length=ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):raise ctypes.WinError()
    return status.available_phys


def mesh_quality(text):
    patterns=dict(minimum_orthogonal_quality=r'Minimum Orthogonal Quality\s*=',
        minimum_cell_volume_m3=r'minimum volume\s*\([^)]*\)\s*:',
        maximum_cell_volume_m3=r'maximum volume\s*\([^)]*\)\s*:',
        total_cell_volume_m3=r'total volume\s*\([^)]*\)\s*:',
        maximum_aspect_ratio=r'Maximum Aspect Ratio\s*=')
    result={}
    for name,pattern in patterns.items():
        matches=re.findall(pattern+r'\s*([-+0-9.eE]+)',text,re.I)
        if matches:result[name]=float(matches[-1])
    if not all(math.isfinite(v) for v in result.values()):raise ValueError('Nonfinite native quality metric')
    if result.get('minimum_orthogonal_quality',-1)<.10 or result.get('minimum_cell_volume_m3',0)<=0:
        raise ValueError('Native positive-volume/orthogonal-quality .10 gate failed or missing')
    return result


def hdf_terminal(path):
    import h5py
    result={};last=None
    with h5py.File(path,'r') as f:
        for name in EQUATIONS:
            group=f['results/residuals/phase-1/'+name]
            iteration=float(group['iterations'][-1]);row=group['data'][-1]
            if iteration!=int(iteration) or row[1]<=0:raise ValueError('Invalid native residual terminal')
            if last is not None and int(iteration)!=last:raise ValueError('Native residual equation endpoints differ')
            value=float(row[0]/row[1])
            if not math.isfinite(value) or value<0:raise ValueError('Nonfinite native residual')
            result[name]=value;last=int(iteration)
    return last,result


def force_rows(text):
    rows={}
    for line in text.splitlines():
        fields=line.split()
        if len(fields)!=3:continue
        try:iteration,cd,cl=map(float,fields)
        except ValueError:continue
        if iteration!=int(iteration) or not all(map(math.isfinite,(iteration,cd,cl))):
            raise ValueError('Invalid actual native force-monitor row')
        rows[int(iteration)]=dict(iteration=int(iteration),Cd=cd,Cl=cl)
    if not rows:raise ValueError('No native iteration/Cd/Cl force-monitor records yet')
    return [rows[key] for key in sorted(rows)]


def force_stability(rows,actual):
    final=rows[-100:]
    complete=len(final)==100 and [row['iteration'] for row in final]==list(range(actual-99,actual+1))
    result=dict(actual_samples=len(final),consecutive_final100_available=complete,
        definition='range/max(abs(mean),1e-3), actual native per-iteration coefficients')
    for name in ('Cd','Cl'):
        values=[row[name] for row in final]
        mean=sum(values)/len(values) if values else None
        relative=(max(values)-min(values))/max(abs(mean),1e-3) if values else None
        result[name]=dict(mean=mean,relative_range=relative,
                         passed=bool(complete and relative<=.001))
    result['passed']=bool(complete and all(result[name]['passed'] for name in ('Cd','Cl')))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    args=parser.parse_args();mesh=args.mesh.resolve();out=args.out.resolve()
    repo=Path(__file__).resolve().parents[1]
    if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo):raise ValueError('Choose fresh external E: output')
    if mesh.drive.upper()!='E:' or not mesh.is_file():raise ValueError('Existing approved E: native mesh required')
    if not EXECUTABLE.is_file():raise FileNotFoundError(EXECUTABLE)
    disk=shutil.disk_usage(out.anchor).free;memory=available_memory()
    if disk<10*1024**3 or memory<4*1024**3:raise RuntimeError('Require10GiBE:disk and4GiBavailableRAM')
    setup_path=mesh.parent/'setup.json';setup=json.loads(setup_path.read_text(encoding='utf-8'))
    if setup.get('kind')!='finite-wing-steady-Re400' or not 1<=setup.get('cells',0)<=250000:
        raise ValueError('Wrong model or hard cell cap exceeded')
    if setup.get('mesh_sha256')!=sha256(mesh):raise ValueError('Native mesh SHA-256 differs from setup')
    assert setup['candidate_quality_gate']['candidate_gate_pass'], 'Candidate mesh quality not accepted'
    out.mkdir(parents=True)
    for name in ('temp','examples','force-reports'):(out/name).mkdir()
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(out/'temp')
    project_temp=Path('E:/CFD-Learning-Lab/temp/wing-refined-fluent-v001');project_temp.mkdir(exist_ok=True)
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(project_temp)
    os.environ['PYTHONDONTWRITEBYTECODE']='1';tempfile.tempdir=str(project_temp);os.chdir(out)
    source=Path(__file__).read_bytes();(out/'executed-script.py').write_bytes(source)
    shutil.copyfile(setup_path,out/'mesh-setup.json')
    write_json(out/'code-provenance.json',dict(source_sha256=hashlib.sha256(source).hexdigest(),
        mesh=str(mesh),mesh_sha256=sha256(mesh),mesh_setup_sha256=sha256(setup_path),
        geometry_npz_sha256=sha256(mesh.parent/'geometry.npz'),
        snapshot_before_solver_launch=True,available_ram_before_launch_bytes=memory,free_disk_before_launch_bytes=disk))
    write_json(out/'declared-targets.json',TARGETS)
    write_json(out/'geometry-setup.json',dict(kind=setup['kind'],mesh_setup=setup,references=REFERENCE,
        inlet_velocity_m_s=[1,0,0],outer_types=dict(inlet='velocity-inlet',farfield='velocity-inlet',outlet='pressure-outlet p_gauge0'),
        wall='stationary no-slip wing',initialization='uniform interior U=(1,0,0),gaugep=0',
        moment_origin_m=[.25,0,0],Cp_reference_gauge_pressure_Pa=0,steady_visualization_only=True,physical_validation=False))
    write_json(out/'requested-settings.json',dict(targets=TARGETS,references=REFERENCE,
        solver='3D pressure-based steady laminar Coupled, CFL20, pseudo-time off',gradient='least-square-cell-based',
        momentum='second-order-upwind',pressure='second-order',processor_count=2,ui_mode='hidden_gui',executable=str(EXECUTABLE)))
    events=[];exports=[];checkpoints=[];session=None;loaded=False;saved=False;actual=0;final_residuals={}

    def update(stage,**values):
        event=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**values)
        events.append(event);write_json(out/'history.json',events);write_json(out/'status.json',event)
        print(json.dumps(event,default=str),flush=True)

    def optional(label,operation):
        try:operation();return True
        except Exception:
            exports.append(dict(export=label,error=traceback.format_exc()))
            write_json(out/'export-errors.json',exports);return False

    try:
        update('importing_pyfluent')
        import ansys.fluent.core as pf
        pf.config.examples_path=str(out/'examples')
        update('launching_single_session',processor_count=2,ui_mode='hidden_gui',dimension=3)
        session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=3,precision='double',processor_count=2,
            ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
        update('solver_connected',version=str(session.get_fluent_version()))
        session.transcript.start(file_name=str(out/'solver-transcript.txt'))
        st=session.settings;st.file.read_mesh(file_name=str(mesh));loaded=True
        st.setup.general.solver.type='pressure-based';st.setup.general.solver.time='steady'
        st.setup.models.viscous.model='laminar'
        st.setup.materials.fluid['air'].density.value=1.;st.setup.materials.fluid['air'].viscosity.value=.0025
        st.setup.cell_zone_conditions.fluid['fluid'].general.material='air'
        for name in ('inlet','farfield'):
            inlet=st.setup.boundary_conditions.velocity_inlet[name].momentum
            inlet.velocity_specification_method='Components'
            inlet.velocity_components=[dict(value=1.),dict(value=0.),dict(value=0.)]
        st.setup.boundary_conditions.pressure_outlet['outlet'].momentum.gauge_pressure.value=0.
        wall=st.setup.boundary_conditions.wall['wing'].momentum
        wall.wall_motion='Stationary Wall';wall.shear_condition='No Slip'
        for key,value in REFERENCE.items():setattr(st.setup.reference_values,key,value)
        methods=st.solution.methods;methods.p_v_coupling.flow_scheme='Coupled'
        methods.pseudo_time_method.formulation.coupled_solver='off'
        st.solution.controls.p_v_controls.flow_courant_number=20.
        st.solution.controls.p_v_controls.explicit_pressure_under_relaxation=.5
        st.solution.controls.p_v_controls.explicit_momentum_under_relaxation=.5
        methods.spatial_discretization.gradient_scheme='least-square-cell-based'
        methods.spatial_discretization.discretization_scheme={'mom':'second-order-upwind','pressure':'second-order'}
        residual=st.solution.monitor.residual
        residual.options.criterion_type='absolute';residual.options.n_save=2000
        for name in EQUATIONS:
            residual.equations[name].absolute_criteria=1e-5;residual.equations[name].check_convergence=False
        for kind,name,direction,output in (('drag','Cd',[1,0,0],'Drag Coefficient'),('lift','Cl',[0,1,0],'Lift Coefficient')):
            collection=getattr(st.solution.report_definitions,kind);collection.create(name=name)
            report=collection[name];report.zones=['wing'];report.force_vector=direction
            report.report_output_type=output;report.average_over=1;report.create_report_file=False;report.create_report_plot=False
        reports=st.solution.monitor.report_files;reports.create(name='forces')
        monitor=reports['forces'];monitor.report_defs=['Cd','Cl'];monitor.file_name=str(out/'force-history.out')
        choices=monitor.frequency_of.allowed_values()
        iteration_choice=next((value for value in choices if 'iter' in value.lower()),None)
        if iteration_choice is None:raise ValueError('Native force monitor has no per-iteration sampling mode: '+str(choices))
        monitor.frequency_of=iteration_choice;monitor.frequency=1;monitor.active=True;monitor.print=False
        write_json(out/'native-force-monitor-settings.json',dict(report_definitions=st.solution.report_definitions.get_state(),report_files=reports.get_state(),native_column_order=['iteration','Cd','Cl']))
        check=st.mesh.check();quality=st.mesh.quality()
        for attempt in range(10):
            text=str(check)+'\n'+str(quality)+'\n'+(out/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace')
            try:metrics=mesh_quality(text);break
            except ValueError:
                if attempt==9:raise
                time.sleep(.1)
        write_json(out/'geometry-quality.json',dict(**metrics,native_numeric_gate_pass=True,checked_before_flow=True,
            orthogonal_quality_floor=.10,positive_volume_required=True,source='native mesh.check/quality transcript'))

        def state():
            return dict(general=st.setup.general.get_state(),viscous=st.setup.models.viscous.get_state(),
                material=st.setup.materials.fluid['air'].get_state(),cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),
                inlet=st.setup.boundary_conditions.velocity_inlet['inlet'].get_state(),farfield=st.setup.boundary_conditions.velocity_inlet['farfield'].get_state(),
                outlet=st.setup.boundary_conditions.pressure_outlet['outlet'].get_state(),walls=dict(wing=st.setup.boundary_conditions.wall['wing'].get_state()),
                references=st.setup.reference_values.get_state(),methods=methods.get_state(),controls=st.solution.controls.get_state(),residual=residual.get_state(),
                report_definitions=st.solution.report_definitions.get_state(),report_files=reports.get_state(),initialization=st.solution.initialization.get_state())

        write_json(out/'loaded-settings.json',state())
        initialization=st.solution.initialization;initialization.initialization_type='standard'
        for key,value in (('x-velocity',1.),('y-velocity',0.),('z-velocity',0.),('pressure',0.)):initialization.defaults[key]=value
        initialization.standard_initialize();st.file.write_case_data(file_name=str(out/'wing-setup.cas.h5'))
        write_json(out/'initialized-settings.json',state())

        def export_forces():
            rows=force_rows((out/'force-history.out').read_text(encoding='utf-8',errors='replace'))
            with (out/'force-history.csv').open('w',newline='',encoding='utf-8') as stream:
                writer=csv.DictWriter(stream,fieldnames=['iteration','Cd','Cl']);writer.writeheader();writer.writerows(rows)
            write_json(out/'force-stability.json',force_stability(rows,actual))

        def checkpoint():
            nonlocal actual,final_residuals
            # Attempt count is not actual count when native convergence stops early.
            stem=f'wing-request{len(checkpoints)+1:02d}'
            st.file.write_case_data(file_name=str(out/(stem+'.cas.h5')))
            actual,final_residuals=hdf_terminal(out/(stem+'.dat.h5'))
            if actual>MAX_ITERATIONS:raise RuntimeError('Native iteration count exceeded declared total cap')
            entry=dict(iteration=actual,case=stem+'.cas.h5',data=stem+'.dat.h5',scaled_residuals=final_residuals)
            checkpoints.append(entry);write_json(out/'checkpoints.json',checkpoints)
            optional('force_history_'+str(actual),export_forces)
            update('checkpoint_saved',**entry)

        run=st.solution.run_calculation
        while actual<MAX_ITERATIONS:
            requested=min(300,MAX_ITERATIONS-actual)
            update('iterating',current_actual_iteration=actual,requested_additional_iterations=requested,total_cap=MAX_ITERATIONS)
            run.iterate(iter_count=requested);previous=actual;checkpoint()
            if actual<=previous:raise RuntimeError('Native iterations did not advance')
        st.file.write_case_data(file_name=str(out/'wing-final.cas.h5'))
        actual,final_residuals=hdf_terminal(out/'wing-final.dat.h5');saved=True
        update('final_saved',actual_iteration=actual,case='wing-final.cas.h5',data='wing-final.dat.h5',scaled_residuals=final_residuals)
        optional('force_history_final',export_forces)
        optional('final_settings',lambda:write_json(out/'final-settings.json',state()))
        optional('residuals_xy',lambda:residual.write(file_name=str(out/'residuals.xy')))
        optional('mass_balance',lambda:st.results.report.fluxes.mass_flow(zones=['inlet','outlet','farfield'],write_to_file=True,file_name=str(out/'mass-balance.txt')))
        for name,direction in (('drag',[1,0,0]),('lift',[0,1,0])):
            optional(name+'_native_force',lambda name=name,direction=direction:st.results.report.forces(option='forces',wall_zones=['wing'],direction_vector=direction,write_to_file=True,file_name=str(out/'force-reports'/(name+'-final.txt'))))

        def picture(name,field,surfaces):
            graphics=st.results.graphics;graphics.picture.use_window_resolution=False
            graphics.picture.x_resolution=1600;graphics.picture.y_resolution=1000
            graphics.picture.invert_background=False;graphics.views.camera.projection='orthographic'
            graphics.contour.create(name=name);contour=graphics.contour[name];contour.field=field;contour.surfaces_list=surfaces;contour.display()
            graphics.views.camera.position=[4.,3.,7.];graphics.views.camera.target=[.5,0,0];graphics.views.camera.field=[6.8,4.25]
            graphics.picture.save_picture(file_name=str(out/(name+'.png')))

        optional('wing_pressure_native_png',lambda:picture('wing-pressure','pressure',['wing']))
        optional('wing_velocity_native_png',lambda:picture('wing-velocity','velocity-magnitude',['wing']))
        optional('transcript_stop',session.transcript.stop)
        summary=dict(kind='finite-wing-steady-Re400',actual_iterations=actual,total_iteration_cap=MAX_ITERATIONS,
            final_scaled_residuals=final_residuals,native_final_saved=True,native_case='wing-final.cas.h5',native_data='wing-final.dat.h5',
            native_per_iteration_force_file='force-history.out',derived_force_csv='force-history.csv',checkpoints=checkpoints,
            force_stability='force-stability.json',targets=TARGETS,export_errors=exports,physical_validation=False,
            numerical_checks='Independent saved-native audit required; failed targets must remain visible',
            visualization='Steady computed fields; animated tracers/cameras do not represent transient CFD')
        write_json(out/'summary.json',summary);update('completed_with_export_errors' if exports else 'completed',**summary)
    except Exception:
        failure=dict(error=traceback.format_exc(),native_final_saved=saved,last_observed_actual_iteration=actual,export_errors=exports)
        write_json(out/'summary.json',dict(stage='failed',**failure));update('failed',**failure);raise
    finally:
        if session is not None:
            if loaded and not saved:
                try:session.settings.file.write_case_data(file_name=str(out/'recovery.cas.h5'))
                except Exception:write_json(out/'recovery-error.json',dict(error=traceback.format_exc()))
            try:session.exit()
            except Exception:write_json(out/'session-exit-error.json',dict(error=traceback.format_exc()))


if __name__=='__main__':main()
