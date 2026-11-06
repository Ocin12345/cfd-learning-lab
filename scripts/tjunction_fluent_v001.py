"""Single authorized hidden 2-process Fluent laminar T-junction merge, 2000 iterations.

No solver starts on import. Requires a root-issued SHA-bound Run authorization.
The fixed iteration cap is deliberate: retain 1800 and2000 field snapshots.
"""
from pathlib import Path
import argparse
import ctypes
from datetime import datetime
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import tempfile
import time
import traceback

EQUATIONS=('continuity','x-velocity','y-velocity')
WALLS=('walls',)
INLETS={'inlet_left':[.1,0.],'inlet_top':[0.,-.1]}
CHECKPOINTS=(500,1000,1500,1800,2000)
MAX_ITERATIONS=2000
EXECUTABLE=Path('D:/ANSYS Inc/ANSYS Student/v261/fluent/ntbin/win64/fluent.exe')
TARGETS=dict(maximum_total_iterations=2000,checkpoints=list(CHECKPOINTS),residual_target=1e-6,
    minimum_orthogonal_quality=.15,local_dimensionless_divergence_max=1e-5,local_dimensionless_divergence_rms_max=1e-6,
    relative_mass_imbalance_max=1e-4,final200_weighted_velocity_L2_relative_change_max=.001,final200_endpoints=[1800,2000],
    rho_kg_m3=1.,U_m_s=.1,L_m=.2,reference_depth_m=1.,source_area_m2=.44,
    divergence_definition='sum signed native face mass flux/(rho*source cell area*reference depth)*L/U; L/U=2s',
    mass_imbalance_definition='abs(native Net boundary mass flux)/actual sum of two inlet mass fluxes',
    solver_processes=2,physical_validation=False,mesh_independence=False,scalar_mixing_simulated=False,
    visualization='Steady computed fields; entry seed colours are not concentration, tracer motion is not transient CFD')


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')


def sha256(path):
    value=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024**2),b''):value.update(block)
    return value.hexdigest()


def available_memory():
    class Status(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(name,ctypes.c_ulonglong) for name in
            ('total_phys','available_phys','total_page','available_page','total_virtual','available_virtual','available_extended')]
    status=Status();status.length=ctypes.sizeof(status)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status)):raise ctypes.WinError()
    return status.available_phys


def require_idle():
    query="@(Get-Process fluent*,cortex* -ErrorAction SilentlyContinue) | ForEach-Object { $_.ProcessName + ':' + $_.Id }"
    found=subprocess.run(['powershell.exe','-NoProfile','-NonInteractive','-Command',query],capture_output=True,text=True,check=True).stdout.strip()
    if found:raise RuntimeError('Existing Fluent session must be saved and closed first: '+found)
    result=subprocess.run(['wsl.exe','--list','--running','--quiet'],capture_output=True,check=True)
    raw=result.stdout
    running=raw.decode('utf-16-le' if b'\x00' in raw else 'utf-8',errors='replace').replace('\x00','').strip()
    if running:raise RuntimeError('Require no concurrently running WSL distribution: '+running)
    return dict(existing_fluent_processes=[],running_wsl_distributions=[])


def validate_authorization(path,mesh,out):
    value=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    required=dict(action='run_tjunction_fluent_v001',approved=True,root='root',allowed_mode='Run')
    if any(value.get(key)!=expected for key,expected in required.items()):raise ValueError('Root tjunction Run authorization absent or wrong mode')
    if value.get('source_sha256')!=sha256(__file__):raise ValueError('Reviewed tjunction driver SHA-256 mismatch')
    if value.get('mesh_sha256')!=sha256(mesh):raise ValueError('Reviewed shared mesh SHA-256 mismatch')
    if Path(value.get('out','')).resolve()!=out:raise ValueError('Authorized output differs')
    return value


def quality_metrics(text):
    patterns=dict(minimum_orthogonal_quality=r'Minimum Orthogonal Quality\s*=',
        minimum_cell_volume_m3=r'minimum volume\s*\([^)]*\)\s*:',
        maximum_cell_volume_m3=r'maximum volume\s*\([^)]*\)\s*:',
        total_cell_volume_m3=r'total volume\s*\([^)]*\)\s*:',maximum_aspect_ratio=r'Maximum Aspect Ratio\s*=')
    result={}
    for name,pattern in patterns.items():
        matches=re.findall(pattern+r'\s*([-+0-9.eE]+)',text,re.I)
        if matches:result[name]=float(matches[-1])
    if not all(map(math.isfinite,result.values())):raise ValueError('Nonfinite native quality metric')
    if result.get('minimum_orthogonal_quality',-1)<.15 or result.get('minimum_cell_volume_m3',0)<=0:
        raise ValueError('Native positive-volume/quality .15 gate failed or missing')
    if abs(result.get('total_cell_volume_m3',0)-.44)>1e-8:raise ValueError('Native tjunction total volume must equal0.44 m3 at reference depth1')
    return result


def hdf_terminal(path):
    import h5py
    values={};actual=None
    with h5py.File(path,'r') as handle:
        for name in EQUATIONS:
            group=handle['results/residuals/phase-1/'+name]
            iteration=float(group['iterations'][-1]);row=group['data'][-1]
            if iteration!=int(iteration) or row[1]<=0:raise ValueError('Invalid native residual endpoint/scaling')
            if actual is not None and actual!=int(iteration):raise ValueError('Native equation endpoints differ')
            value=float(row[0]/row[1])
            if not math.isfinite(value) or value<0:raise ValueError('Invalid native residual value')
            actual=int(iteration);values[name]=value
    return actual,values


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mesh',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--authorization',type=Path,required=True)
    args=parser.parse_args();mesh=args.mesh.resolve();out=args.out.resolve();repo=Path(__file__).resolve().parents[1]
    if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo):raise ValueError('Fresh native E: output outside repo required')
    if mesh.drive.upper()!='E:' or not mesh.is_file():raise ValueError('Existing approved E: native mesh required')
    if not EXECUTABLE.is_file():raise FileNotFoundError(EXECUTABLE)
    authorization=validate_authorization(args.authorization,mesh,out)
    setup_path=mesh.parent/'setup.json';geometry_path=mesh.parent/'geometry.npz'
    setup=json.loads(setup_path.read_text(encoding='utf-8-sig'))
    if (setup.get('mesh_sha256')!=sha256(mesh) or setup.get('cells')!=4400):raise ValueError('Approved4400cell shared T-junction mesh required')
    import numpy as np
    with np.load(geometry_path,allow_pickle=False) as geometry:
        if geometry['quads'].shape!=(4400,4) or geometry['nodes'].shape!=(4641,2):raise ValueError('Unexpected shared tjunction topology')
        if not np.isfinite(geometry['areas']).all() or (geometry['areas']<=0).any() or abs(float(geometry['areas'].sum())-.44)>1e-10:
            raise ValueError('Positive source cell areas summing to0.44 required')
    disk=shutil.disk_usage(out.anchor).free;memory=available_memory()
    if disk<10*1024**3 or memory<4*1024**3:raise RuntimeError('Require10GiBE: free and4GiBavailableRAM')
    idle=require_idle()
    out.mkdir(parents=True);(out/'examples').mkdir()
    temporary=Path('E:/CFD-Learning-Lab/temp/tjunction-v001-fluent');temporary.mkdir(parents=True,exist_ok=True)
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
    os.environ['PYTHONDONTWRITEBYTECODE']='1';tempfile.tempdir=str(temporary);os.chdir(out)
    source=Path(__file__).read_bytes();(out/'executed-script.py').write_bytes(source)
    shutil.copyfile(setup_path,out/'mesh-setup.json');write_json(out/'authorization.json',authorization)
    write_json(out/'code-provenance.json',dict(source_sha256=sha256(__file__),mesh=str(mesh),mesh_sha256=sha256(mesh),
        mesh_setup_sha256=sha256(setup_path),geometry_npz_sha256=sha256(geometry_path),snapshot_before_solver_launch=True,
        available_ram_before_launch_bytes=memory,free_disk_before_launch_bytes=disk,temporary_directory=str(temporary),idle_check=idle))
    write_json(out/'declared-targets.json',TARGETS)
    write_json(out/'requested-settings.json',dict(solver='2D planar steady pressure-based laminar SIMPLE',density=1.,viscosity=.001,
        operating_pressure_Pa=0.,reference_depth_m=1.,energy=False,inlet_velocities_m_s=INLETS,outlet_gauge_pressure_Pa=0.,
        gradient='least-square-cell-based',momentum='second-order-upwind',pressure='second-order',
        processor_count=2,ui_mode='hidden_gui',automatic_convergence_stopping=False,targets=TARGETS))
    events=[];errors=[];checkpoints=[];session=None;loaded=False;saved=False;actual=0;final_residuals={}
    def update(stage,**values):
        value=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**values)
        events.append(value);write_json(out/'history.json',events);write_json(out/'status.json',value)
        print(json.dumps(value,default=str),flush=True)
    def optional(label,operation):
        try:operation();return True
        except Exception:
            errors.append(dict(export=label,error=traceback.format_exc()));write_json(out/'export-errors.json',errors);return False
    try:
        update('importing_pyfluent')
        import ansys.fluent.core as pf
        pf.config.examples_path=str(out/'examples')
        update('launching_single_session',processor_count=2,dimension=2,ui_mode='hidden_gui')
        session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,
            ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
        update('solver_connected',version=str(session.get_fluent_version()))
        session.transcript.start(file_name=str(out/'solver-transcript.txt'))
        st=session.settings;st.file.read_mesh(file_name=str(mesh));loaded=True
        st.setup.general.solver.type='pressure-based';st.setup.general.solver.two_dim_space='planar';st.setup.general.solver.time='steady'
        st.setup.models.energy.enabled=False;st.setup.models.viscous.model='laminar'
        st.setup.general.operating_conditions.operating_pressure=0.
        air=st.setup.materials.fluid['air'];air.density.value=1.;air.viscosity.value=.001
        st.setup.cell_zone_conditions.fluid['fluid'].general.material='air'
        for name,components in INLETS.items():
            inlet=st.setup.boundary_conditions.velocity_inlet[name].momentum
            inlet.velocity_specification_method='Components';inlet.velocity_components=[dict(value=value) for value in components]
        st.setup.boundary_conditions.pressure_outlet['outlet'].momentum.gauge_pressure.value=0.
        wall=st.setup.boundary_conditions.wall['walls'].momentum;wall.wall_motion='Stationary Wall';wall.shear_condition='No Slip'
        for key,value in dict(density=1.,viscosity=.001,velocity=.1,length=.2,area=.2,depth=1.).items():setattr(st.setup.reference_values,key,value)
        methods=st.solution.methods;methods.p_v_coupling.flow_scheme='SIMPLE'
        methods.spatial_discretization.gradient_scheme='least-square-cell-based'
        methods.spatial_discretization.discretization_scheme={'mom':'second-order-upwind','pressure':'second-order'}
        residual=st.solution.monitor.residual;residual.options.criterion_type='absolute';residual.options.n_save=2500
        for name in EQUATIONS:
            residual.equations[name].absolute_criteria=1e-6;residual.equations[name].check_convergence=False
        definitions=st.solution.report_definitions.volume
        for name,field in (('mean_speed','velocity-magnitude'),('mean_pressure','pressure')):
            definitions.create(name=name);report=definitions[name];report.report_type='volume-average';report.field=field
            report.cell_zones=['fluid'];report.average_over=1;report.create_report_file=False;report.create_report_plot=False
        files=st.solution.monitor.report_files;files.create(name='tjunction_fields');monitor=files['tjunction_fields']
        monitor.report_defs=['mean_speed','mean_pressure'];monitor.file_name=str(out/'field-monitor-history.out')
        choices=monitor.frequency_of.allowed_values();choice=next((v for v in choices if 'iter' in v.lower()),None)
        if choice is None:raise ValueError('No native per-iteration report mode')
        monitor.frequency_of=choice;monitor.frequency=1;monitor.active=True;monitor.print=False
        flux_defs=st.solution.report_definitions.flux
        for zone in ('inlet_left','inlet_top','outlet'):
            name='mdot_'+zone;flux_defs.create(name=name);report=flux_defs[name]
            report.report_type='flux-massflow';report.boundaries=[zone];report.average_over=1;report.create_report_file=False;report.create_report_plot=False
        files.create(name='flow_balance');mass_monitor=files['flow_balance']
        mass_monitor.report_defs=['mdot_inlet_left','mdot_inlet_top','mdot_outlet'];mass_monitor.file_name=str(out/'massflow-history.out')
        mass_monitor.frequency_of=choice;mass_monitor.frequency=1;mass_monitor.active=True;mass_monitor.print=False
        check=st.mesh.check();quality=st.mesh.quality()
        for attempt in range(10):
            text=str(check)+'\n'+str(quality)+'\n'+(out/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace')
            try:metrics=quality_metrics(text);break
            except ValueError:
                if attempt==9:raise
                time.sleep(.1)
        write_json(out/'geometry-quality.json',dict(**metrics,native_numeric_gate_pass=True,checked_before_flow=True))
        def state():
            return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),
                material=air.get_state(),cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),
                walls={name:st.setup.boundary_conditions.wall[name].get_state() for name in WALLS},
                inlets={name:st.setup.boundary_conditions.velocity_inlet[name].get_state() for name in INLETS},
                outlet=st.setup.boundary_conditions.pressure_outlet['outlet'].get_state(),
                references=st.setup.reference_values.get_state(),methods=methods.get_state(),controls=st.solution.controls.get_state(),
                residual=residual.get_state(),report_definitions=st.solution.report_definitions.get_state(),report_files=files.get_state(),
                initialization=st.solution.initialization.get_state())
        readback=dict(inlets={},wall_motion=wall.wall_motion(),wall_shear_condition=wall.shear_condition(),
            outlet_gauge_pressure_Pa=st.setup.boundary_conditions.pressure_outlet['outlet'].momentum.gauge_pressure.value(),
            operating_pressure_Pa=st.setup.general.operating_conditions.operating_pressure(),reference_depth_m=st.setup.reference_values.depth(),
            density=air.density.value(),viscosity=air.viscosity.value(),viscous_model=st.setup.models.viscous.model(),energy_enabled=st.setup.models.energy.enabled())
        for name,expected in INLETS.items():
            inlet=st.setup.boundary_conditions.velocity_inlet[name].momentum
            values=[item['value'] for item in inlet.velocity_components.get_state()]
            if inlet.velocity_specification_method()!='Components' or values!=expected:raise ValueError('Native inlet components differ: '+name)
            readback['inlets'][name]=dict(velocity_specification_method=inlet.velocity_specification_method(),velocity_m_s=values)
        if (readback['wall_motion']!='Stationary Wall' or readback['wall_shear_condition']!='No Slip'
            or readback['outlet_gauge_pressure_Pa']!=0 or readback['operating_pressure_Pa']!=0 or readback['reference_depth_m']!=1
            or readback['density']!=1 or readback['viscosity']!=.001 or readback['viscous_model']!='laminar' or readback['energy_enabled']):
            raise ValueError('Native T-junction physical readback failed: '+str(readback))
        write_json(out/'native-physical-readback.json',readback);write_json(out/'loaded-settings.json',state())
        initialization=st.solution.initialization;initialization.initialization_type='standard'
        for key,value in (('x-velocity',0.),('y-velocity',0.),('pressure',0.)):initialization.defaults[key]=value
        initialization.standard_initialize();st.file.write_case_data(file_name=str(out/'tjunction-setup.cas.h5'))
        write_json(out/'initialized-settings.json',state())
        for endpoint in CHECKPOINTS:
            update('iterating',current_actual_iteration=actual,requested_additional_iterations=endpoint-actual,total_cap=MAX_ITERATIONS)
            st.solution.run_calculation.iterate(iter_count=endpoint-actual)
            stem=f'tjunction-iter{endpoint:04d}';st.file.write_case_data(file_name=str(out/(stem+'.cas.h5')))
            actual,final_residuals=hdf_terminal(out/(stem+'.dat.h5'))
            if actual!=endpoint:raise RuntimeError('Fixed checkpoint iteration mismatch: '+str(actual)+' vs '+str(endpoint))
            item=dict(iteration=actual,case=stem+'.cas.h5',data=stem+'.dat.h5',scaled_residuals=final_residuals)
            checkpoints.append(item);write_json(out/'checkpoints.json',checkpoints);update('checkpoint_saved',**item)
        st.file.write_case_data(file_name=str(out/'tjunction-final.cas.h5'));actual,final_residuals=hdf_terminal(out/'tjunction-final.dat.h5');saved=True
        optional('final_settings',lambda:write_json(out/'final-settings.json',state()))
        optional('native_residuals_xy',lambda:residual.write(file_name=str(out/'residuals.xy')))
        optional('native_total_boundary_mass_flux_report',lambda:st.results.report.fluxes.mass_flow(zones=['inlet_left','inlet_top','outlet','walls'],write_to_file=True,file_name=str(out/'mass-balance.txt')))
        def picture(name,field):
            graphics=st.results.graphics;graphics.picture.use_window_resolution=False
            graphics.picture.x_resolution=1600;graphics.picture.y_resolution=900;graphics.picture.invert_background=False
            graphics.views.camera.projection='orthographic';graphics.contour.create(name=name)
            contour=graphics.contour[name];contour.field=field;contour.surfaces_list=[];contour.display()
            graphics.views.camera.position=[.4,.3,10.];graphics.views.camera.target=[.4,.3,0.];graphics.views.camera.field=[1.84,1.035]
            graphics.picture.save_picture(file_name=str(out/(name+'.png')))
        optional('native_speed_png',lambda:picture('velocity','velocity-magnitude'))
        optional('native_gauge_pressure_png',lambda:picture('pressure','pressure'))
        optional('transcript_stop',session.transcript.stop)
        summary=dict(kind='laminar-T-junction-merge-Re20-inlets',actual_iterations=actual,total_iteration_cap=MAX_ITERATIONS,
            final_scaled_residuals=final_residuals,native_final_saved=True,native_case='tjunction-final.cas.h5',native_data='tjunction-final.dat.h5',
            checkpoints=checkpoints,targets=TARGETS,export_errors=errors,physical_validation=False,
            numerical_checks='Pending independent saved-file audit of actual inlet/outlet fluxes and all-cell native divergence',
            visualization='Steady computed fields; coloured inlet tracers identify entry only, not scalar concentration or transient CFD')
        write_json(out/'summary.json',summary);update('completed_with_export_errors' if errors else 'completed',**summary)
    except Exception:
        failure=dict(error=traceback.format_exc(),native_final_saved=saved,last_observed_actual_iteration=actual,export_errors=errors)
        write_json(out/'summary.json',dict(stage='failed',**failure));update('failed',**failure);raise
    finally:
        if session is not None:
            if loaded and not saved:
                try:session.settings.file.write_case_data(file_name=str(out/'recovery.cas.h5'))
                except Exception:write_json(out/'recovery-error.json',dict(error=traceback.format_exc()))
            try:session.exit()
            except Exception:write_json(out/'session-exit-error.json',dict(error=traceback.format_exc()))

if __name__=='__main__':main()
