"""One bounded hidden two-process planar Euler nozzle calculation.

Only --authorization with the root-reviewed source and shared-mesh SHA-256 permits
launch. Steady density-based implicit, 300 first-order iterations at CFL1 followed
by 1700 first-order at CFL3 as ONE order-only control. First-order convection
has larger numerical diffusion and can broaden shocks. No automatic restart,
extension, attributed error cause, or physical claim.
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

EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
EQUATIONS=('continuity','x-velocity','y-velocity','energy')
MAX_ITERATIONS=2000
CONSTANTS=dict(gamma=1.4,R_J_kg_K=287.05,Cp_J_kg_K=1004.675,
               molecular_weight_kg_kmol=8314.46261815324/287.05,
               inlet_total_pressure_Pa=200000.,inlet_total_temperature_K=300.,
               outlet_static_pressure_Pa=140000.,operating_pressure_Pa=0.)
TARGETS=dict(maximum_total_iterations=MAX_ITERATIONS,
    stage1=dict(iterations=300,scheme='first-order-upwind',courant_number=1.),
    stage2=dict(iterations=1700,scheme='first-order-upwind',courant_number=3.),
    checkpoints=[300,600,900,1200,1500,1800,2000],
    residual_targets=dict(continuity=1e-5,**{'x-velocity':1e-5,'y-velocity':1e-5,'energy':1e-7}),
    relative_mass_imbalance_target=1e-4,minimum_orthogonal_quality=.15,
    massflow_stability_samples=100,massflow_stability_relative_range_max=.001,
    massflow_stability_definition='(max-min)/max(abs(mean),1e-9), actual final100 consecutive native monitor rows',
    solver_processes=2,mesh_independence=False,physical_validation=False,
    solver_is_steady=True,animation_is_transient_CFD=False)


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')


def sha256(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):result.update(block)
    return result.hexdigest()


def available_memory():
    class MemoryStatus(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[
            (name,ctypes.c_ulonglong) for name in ('total_phys','available_phys','total_page','available_page',
                'total_virtual','available_virtual','available_extended')]
    status=MemoryStatus();status.length=ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):raise ctypes.WinError()
    return status.available_phys


def quality_metrics(text):
    patterns=dict(minimum_orthogonal_quality=r'Minimum Orthogonal Quality\s*=',
        minimum_cell_volume_m3=r'minimum volume\s*\([^)]*\)\s*:',
        maximum_cell_volume_m3=r'maximum volume\s*\([^)]*\)\s*:',
        total_cell_volume_m3=r'total volume\s*\([^)]*\)\s*:',
        maximum_aspect_ratio=r'Maximum Aspect Ratio\s*=')
    result={}
    for name,pattern in patterns.items():
        values=re.findall(pattern+r'\s*([-+0-9.eE]+)',text,re.I)
        if values:result[name]=float(values[-1])
    if not all(map(math.isfinite,result.values())):raise ValueError('Nonfinite native mesh metrics')
    if result.get('minimum_orthogonal_quality',-1)<.15 or result.get('minimum_cell_volume_m3',0)<=0:
        raise ValueError('Native positive-volume/orthogonal-quality .15 gate failed or missing')
    return result


def hdf_terminal(path):
    import h5py
    result={};actual=None
    with h5py.File(path,'r') as handle:
        for equation in EQUATIONS:
            group=handle['results/residuals/phase-1/'+equation]
            iteration=float(group['iterations'][-1]);row=group['data'][-1]
            if iteration!=int(iteration) or row[1]<=0:raise ValueError('Invalid native residual endpoint/normalization')
            if actual is not None and actual!=int(iteration):raise ValueError('Native residual equation endpoints differ')
            value=float(row[0]/row[1])
            if not math.isfinite(value) or value<0:raise ValueError('Nonfinite native scaled residual')
            actual=int(iteration);result[equation]=value
    return actual,result


def monitor_rows(path):
    rows={}
    for line in Path(path).read_text(encoding='utf-8',errors='replace').splitlines():
        fields=line.split()
        if len(fields)!=3:continue
        try:iteration,mdot_in,mdot_out=map(float,fields)
        except ValueError:continue
        if iteration!=int(iteration) or not all(map(math.isfinite,(iteration,mdot_in,mdot_out))):
            raise ValueError('Invalid native massflow monitor row')
        rows[int(iteration)]=dict(iteration=int(iteration),mdot_in=mdot_in,mdot_out=mdot_out)
    if not rows:raise ValueError('No actual per-iteration massflow rows')
    return [rows[key] for key in sorted(rows)]


def stability(rows,actual):
    final=rows[-100:]
    complete=len(final)==100 and [row['iteration'] for row in final]==list(range(actual-99,actual+1))
    result=dict(actual_samples=len(final),consecutive_final100_available=complete,
        definition='range/max(abs(mean),1e-9), actual native per-iteration flux-massflow reports')
    for name in ('mdot_in','mdot_out'):
        values=[row[name] for row in final];mean=sum(values)/len(values) if values else None
        relative=(max(values)-min(values))/max(abs(mean),1e-9) if values else None
        result[name]=dict(mean=mean,relative_range=relative,passed=bool(complete and relative<=.001))
    result['passed']=bool(complete and all(result[name]['passed'] for name in ('mdot_in','mdot_out')))
    return result


def validate_authorization(path,mesh,out):
    auth=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    required=dict(action='run_nozzle_fluent_v004',approved=True,root='root',allowed_mode='Run')
    if any(auth.get(key)!=value for key,value in required.items()):raise ValueError('Root run authorization absent or wrong action/mode')
    if auth.get('source_sha256')!=sha256(__file__):raise ValueError('Reviewed driver SHA-256 mismatch')
    if auth.get('mesh_sha256')!=sha256(mesh):raise ValueError('Reviewed shared mesh SHA-256 mismatch')
    if Path(auth.get('out','')).resolve()!=out:raise ValueError('Authorization output path mismatch')
    return auth


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',required=True,type=Path)
    parser.add_argument('--out',required=True,type=Path)
    parser.add_argument('--authorization',required=True,type=Path)
    args=parser.parse_args();mesh=args.mesh.resolve();out=args.out.resolve()
    repo=Path(__file__).resolve().parents[1]
    if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo):raise ValueError('Fresh external E: native output required')
    if mesh.drive.upper()!='E:' or not mesh.is_file():raise ValueError('Existing approved E: mesh required')
    if not EXECUTABLE.is_file():raise FileNotFoundError(EXECUTABLE)
    authorization=validate_authorization(args.authorization,mesh,out)
    setup_path=mesh.parent/'setup.json';setup=json.loads(setup_path.read_text(encoding='utf-8-sig'))
    if setup.get('mesh_sha256')!=sha256(mesh):raise ValueError('Native mesh differs from setup SHA-256')
    if setup.get('cells')!=9600:raise ValueError('Approved fixed 240x40 shared mesh requires9600cells')
    disk=shutil.disk_usage(out.anchor).free;memory=available_memory()
    if disk<10*1024**3 or memory<4*1024**3:raise RuntimeError('Require10GiBE:disk and4GiBavailableRAM')
    out.mkdir(parents=True)
    # Keep this project's temporary files under its requested E: temp root.
    temporary=Path('E:/CFD-Learning-Lab/temp/nozzle-v001-fluent-order-control-v004');temporary.mkdir(parents=True,exist_ok=True)
    (out/'examples').mkdir();(out/'native-reports').mkdir()
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
    os.environ['PYTHONDONTWRITEBYTECODE']='1';tempfile.tempdir=str(temporary);os.chdir(out)
    source=Path(__file__).read_bytes();(out/'executed-script.py').write_bytes(source)
    shutil.copyfile(setup_path,out/'mesh-setup.json')
    write_json(out/'authorization.json',authorization)
    write_json(out/'code-provenance.json',dict(source_sha256=hashlib.sha256(source).hexdigest(),
        mesh=str(mesh),mesh_sha256=sha256(mesh),mesh_setup_sha256=sha256(setup_path),
        geometry_npz_sha256=sha256(mesh.parent/'geometry.npz'),snapshot_before_solver_launch=True,
        available_ram_before_launch_bytes=memory,free_disk_before_launch_bytes=disk,temporary_directory=str(temporary)))
    write_json(out/'declared-targets.json',TARGETS)
    write_json(out/'order-control.json',dict(control='ONE first-order convection control against v003',
        only_numerical_change='Stage2 convection remains first-order-upwind instead of second-order-upwind',
        same_total_cap=2000,same_mesh_BCs_gas_initialization_and_CFL_schedule=True,
        limitation='Larger numerical diffusion and possible shock broadening; no mesh independence or attributed cause',
        previous_v003_numeric_misses_preserved=True,physical_validation=False))
    write_json(out/'geometry-setup.json',dict(kind='planar-full-nozzle-Euler-ideal-gas',mesh_setup=setup,constants=CONSTANTS,
        inlet='pressure inlet: totalP200000Pa,totalT300K, normal inward direction',
        outlet='pressure outlet: staticP140000Pa,backflowtotalT300K',
        walls='inviscid impermeable slip,adiabatic; no turbulent or boundary-layer model',
        initialization=dict(pressure_Pa=140000.,temperature_K=300.,velocity_m_s=[50.,0.]),
        physical_validation=False,visualization='steady computed fields; animated tracers do not represent transient flow'))
    write_json(out/'requested-settings.json',dict(solver='2D planar steady density-based implicit Euler',
        constants=CONSTANTS,targets=TARGETS,flux='Roe-FDS',gradient='least-square-cell-based',
        processor_count=2,ui_mode='hidden_gui',executable=str(EXECUTABLE),automatic_convergence_stopping=False))
    events=[];errors=[];checkpoints=[];session=None;loaded=False;initialized=False;saved=False;actual=0;final_residuals={}

    def update(stage,**values):
        event=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**values)
        events.append(event);write_json(out/'history.json',events);write_json(out/'status.json',event)
        print(json.dumps(event,default=str),flush=True)

    def optional(label,operation):
        try:operation();return True
        except Exception:
            errors.append(dict(export=label,error=traceback.format_exc()))
            write_json(out/'export-errors.json',errors);return False

    try:
        update('importing_pyfluent')
        import ansys.fluent.core as pf
        pf.config.examples_path=str(out/'examples')
        update('launching_single_session',processor_count=2,ui_mode='hidden_gui',dimension=2)
        session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,
            ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
        update('solver_connected',version=str(session.get_fluent_version()),licence='Actual two-process Student session connected')
        session.transcript.start(file_name=str(out/'solver-transcript.txt'))
        st=session.settings;st.file.read_mesh(file_name=str(mesh));loaded=True
        st.setup.general.solver.type='density-based-implicit'
        st.setup.general.solver.two_dim_space='planar';st.setup.general.solver.time='steady'
        st.setup.models.energy.enabled=True;st.setup.models.viscous.model='inviscid'
        st.setup.general.operating_conditions.operating_pressure=0.
        air=st.setup.materials.fluid['air'];air.density.option='ideal-gas'
        air.specific_heat.option='constant';air.specific_heat.value=CONSTANTS['Cp_J_kg_K']
        air.molecular_weight.value=CONSTANTS['molecular_weight_kg_kmol']
        if air.thermal_conductivity.is_active():
            air.thermal_conductivity.option='constant';air.thermal_conductivity.value=0.
        st.setup.cell_zone_conditions.fluid['fluid'].general.material='air'
        st.setup.reference_values.depth=1.
        inlet=st.setup.boundary_conditions.pressure_inlet['inlet']
        inlet.momentum.gauge_total_pressure.value=200000.
        inlet.momentum.supersonic_or_initial_gauge_pressure.value=180000.
        inlet.momentum.direction_specification_method='Normal to Boundary'
        inlet.thermal.total_temperature.value=300.
        outlet=st.setup.boundary_conditions.pressure_outlet['outlet']
        outlet.momentum.gauge_pressure.value=140000.;outlet.thermal.backflow_total_temperature.value=300.
        # Euler walls are impermeable slip with no heat diffusion. The wall thermal
        # group exists but its heat-flux child is inactive in this model; do not set it.
        methods=st.solution.methods;methods.p_v_coupling.flux_type='Roe-FDS'
        methods.spatial_discretization.gradient_scheme='least-square-cell-based'
        discretization=methods.spatial_discretization.discretization_scheme
        keys=discretization.get_object_names()
        schemes={name:discretization[name].allowed_values() for name in keys}
        write_json(out/'available-spatial-schemes.json',dict(actual_keys=keys,allowed_values=schemes))
        if keys!=['amg-c']:raise ValueError('Unexpected coupled Euler discretization keys: '+str(keys))
        if not all(name in schemes['amg-c'] for name in ('first-order-upwind','second-order-upwind')):
            raise ValueError('Required native amg-c schemes unavailable')
        def set_scheme(scheme):
            discretization['amg-c']=scheme
            if discretization['amg-c']()!=scheme:raise ValueError('Native amg-c scheme readback differs')
        set_scheme('first-order-upwind');st.solution.controls.courant_number=1.
        residual=st.solution.monitor.residual;residual.options.criterion_type='absolute';residual.options.n_save=3000
        active=list(residual.equations.get_object_names())
        if not all(name in active for name in EQUATIONS):raise ValueError('Required Euler residual equations unavailable: '+str(active))
        for name in EQUATIONS:
            residual.equations[name].absolute_criteria=TARGETS['residual_targets'][name]
            residual.equations[name].check_convergence=False
        flux=st.solution.report_definitions.flux
        for report_name,zone in (('mdot_in','inlet'),('mdot_out','outlet')):
            flux.create(name=report_name);definition=flux[report_name]
            definition.report_type='flux-massflow';definition.boundaries=[zone]
            definition.average_over=1;definition.create_report_file=False;definition.create_report_plot=False
        files=st.solution.monitor.report_files;files.create(name='massflows');monitor=files['massflows']
        monitor.report_defs=['mdot_in','mdot_out'];monitor.file_name=str(out/'massflow-history.out')
        choices=monitor.frequency_of.allowed_values()
        iteration_choice=next((value for value in choices if 'iter' in value.lower()),None)
        if iteration_choice is None:raise ValueError('No native per-iteration monitor mode')
        monitor.frequency_of=iteration_choice;monitor.frequency=1;monitor.active=True;monitor.print=False
        write_json(out/'native-mass-monitor-settings.json',dict(report_definitions=st.solution.report_definitions.get_state(),
            report_files=files.get_state(),native_column_order=['iteration','mdot_in','mdot_out'],
            sign_convention='Native Fluent flux signs retained; determine inlet/outlet signs from saved rows and native report'))
        check=st.mesh.check();quality=st.mesh.quality()
        for attempt in range(10):
            text=str(check)+'\n'+str(quality)+'\n'+(out/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace')
            try:metrics=quality_metrics(text);break
            except ValueError:
                if attempt==9:raise
                time.sleep(.1)
        write_json(out/'geometry-quality.json',dict(**metrics,native_numeric_gate_pass=True,checked_before_flow=True,
            orthogonal_quality_floor=.15,positive_volume_required=True,source='native mesh.check/quality transcript'))
        def state():
            return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),
                viscous=st.setup.models.viscous.get_state(),material=air.get_state(),
                cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),inlet=inlet.get_state(),outlet=outlet.get_state(),
                walls={name:st.setup.boundary_conditions.wall[name].get_state() for name in ('upper_wall','lower_wall')},
                references=st.setup.reference_values.get_state(),
                methods=methods.get_state(),controls=st.solution.controls.get_state(),residual=residual.get_state(),
                report_definitions=st.solution.report_definitions.get_state(),report_files=files.get_state(),
                initialization=st.solution.initialization.get_state())
        write_json(out/'loaded-settings.json',state())
        material_readback=air.get_state();Cp=float(material_readback['specific_heat']['value'])
        Mw=float(material_readback['molecular_weight']['value']);R=8314.46261815324/Mw
        gas_readback=dict(density_model=material_readback['density']['option'],Cp_J_kg_K=Cp,
            molecular_weight_kg_kmol=Mw,universal_R_J_kmol_K=8314.46261815324,
            derived_R_J_kg_K=R,derived_gamma=Cp/(Cp-R),
            operating_pressure_Pa=float(st.setup.general.operating_conditions.operating_pressure()),
            reference_depth_m=float(st.setup.reference_values.depth()),viscous_model=st.setup.models.viscous.model())
        write_json(out/'gas-constant-readback.json',gas_readback)
        if (gas_readback['density_model']!='ideal-gas' or gas_readback['viscous_model']!='inviscid'
            or abs(R-287.05)>1e-8 or abs(Cp-1004.675)>1e-8
            or abs(gas_readback['derived_gamma']-1.4)>1e-12
            or gas_readback['operating_pressure_Pa']!=0. or gas_readback['reference_depth_m']!=1.):
            raise ValueError('Loaded gas/operating-pressure/depth settings differ from declared Euler model')
        initialization=st.solution.initialization;initialization.initialization_type='standard'
        for key,value in (('x-velocity',50.),('y-velocity',0.),('pressure',140000.),('temperature',300.)):
            initialization.defaults[key]=value
        initialization.standard_initialize();initialized=True
        st.file.write_case_data(file_name=str(out/'nozzle-setup.cas.h5'))
        write_json(out/'initialized-settings.json',state())
        def export_massflows():
            rows=monitor_rows(out/'massflow-history.out')
            with (out/'massflow-history.csv').open('w',newline='',encoding='utf-8') as stream:
                writer=csv.DictWriter(stream,fieldnames=['iteration','mdot_in','mdot_out']);writer.writeheader();writer.writerows(rows)
            write_json(out/'massflow-stability.json',stability(rows,actual))
        for requested_endpoint in TARGETS['checkpoints']:
            if actual==300:
                set_scheme('first-order-upwind');st.solution.controls.courant_number=3.
                write_json(out/'order-control-stage2-settings.json',state())
                update('first_order_control_stage2_started',actual_iteration=actual,courant_number=3.)
            count=requested_endpoint-actual
            if count<=0:raise RuntimeError('Nonpositive bounded iteration request')
            update('iterating',current_actual_iteration=actual,requested_additional_iterations=count,total_cap=MAX_ITERATIONS)
            st.solution.run_calculation.iterate(iter_count=count)
            stem=f'nozzle-iter{requested_endpoint:04d}'
            st.file.write_case_data(file_name=str(out/(stem+'.cas.h5')))
            actual,final_residuals=hdf_terminal(out/(stem+'.dat.h5'))
            if actual!=requested_endpoint or actual>MAX_ITERATIONS:raise RuntimeError('Actual native iteration differs from fixed bounded schedule')
            entry=dict(iteration=actual,case=stem+'.cas.h5',data=stem+'.dat.h5',scaled_residuals=final_residuals,
                spatial_scheme='first-order-upwind',courant_number=1. if actual==300 else 3.)
            checkpoints.append(entry);write_json(out/'checkpoints.json',checkpoints)
            optional('massflow_history_'+str(actual),export_massflows)
            update('checkpoint_saved',**entry)
        st.file.write_case_data(file_name=str(out/'nozzle-final.cas.h5'))
        actual,final_residuals=hdf_terminal(out/'nozzle-final.dat.h5');saved=True
        update('final_saved',actual_iteration=actual,case='nozzle-final.cas.h5',data='nozzle-final.dat.h5',scaled_residuals=final_residuals)
        optional('massflow_history_final',export_massflows)
        optional('final_settings',lambda:write_json(out/'final-settings.json',state()))
        optional('residuals_xy',lambda:residual.write(file_name=str(out/'residuals.xy')))
        optional('mass_balance',lambda:st.results.report.fluxes.mass_flow(zones=['inlet','outlet','upper_wall','lower_wall'],
            write_to_file=True,file_name=str(out/'mass-balance.txt')))
        for zone in ('inlet','outlet'):
            for field in ('pressure','total-pressure','temperature','total-temperature','mach-number'):
                optional(zone+'_'+field,lambda zone=zone,field=field:
                    st.results.report.surface_integrals.area_weighted_avg(surface_names=[zone],report_of=field,
                        write_to_file=True,file_name=str(out/'native-reports'/(zone+'-'+field+'.txt'))))
        optional('full_cell_ascii',lambda:st.file.export.ascii(file_name=str(out/'nozzle-final-cell-fields.csv'),
            surface_name_list=[],delimiter='comma',
            quantities=['pressure','density','x-velocity','y-velocity','temperature','mach-number'],location='cell-center'))
        def picture(name,field):
            graphics=st.results.graphics;graphics.picture.use_window_resolution=False
            graphics.picture.x_resolution=1800;graphics.picture.y_resolution=650;graphics.picture.invert_background=False
            graphics.views.camera.projection='orthographic';graphics.contour.create(name=name)
            contour=graphics.contour[name];contour.field=field;contour.surfaces_list=[];contour.display()
            graphics.views.camera.position=[.5,0.,10.];graphics.views.camera.target=[.5,0.,0.]
            graphics.views.camera.field=[1.15,.4152778]
            graphics.picture.save_picture(file_name=str(out/(name+'.png')))
        optional('native_mach_png',lambda:picture('nozzle-mach','mach-number'))
        optional('native_pressure_png',lambda:picture('nozzle-pressure','pressure'))
        optional('transcript_stop',session.transcript.stop)
        summary=dict(kind='planar-full-nozzle-Euler-ideal-gas',actual_iterations=actual,total_iteration_cap=MAX_ITERATIONS,
            final_scaled_residuals=final_residuals,native_final_saved=True,native_case='nozzle-final.cas.h5',native_data='nozzle-final.dat.h5',
            native_per_iteration_mass_file='massflow-history.out',checkpoints=checkpoints,targets=TARGETS,export_errors=errors,
            physical_validation=False,numerical_checks='Saved-native independent audit required; failed targets remain visible',
            visualization='Steady computed fields; animated tracers do not represent transient nozzle evolution')
        write_json(out/'summary.json',summary);update('completed_with_export_errors' if errors else 'completed',**summary)
    except Exception:
        failure=dict(error=traceback.format_exc(),native_final_saved=saved,last_observed_actual_iteration=actual,export_errors=errors)
        write_json(out/'summary.json',dict(stage='failed',**failure));update('failed',**failure);raise
    finally:
        if session is not None:
            if loaded and initialized and not saved:
                try:session.settings.file.write_case_data(file_name=str(out/'recovery.cas.h5'))
                except Exception:write_json(out/'recovery-error.json',dict(error=traceback.format_exc()))
            try:session.exit()
            except Exception:write_json(out/'session-exit-error.json',dict(error=traceback.format_exc()))


if __name__=='__main__':main()



