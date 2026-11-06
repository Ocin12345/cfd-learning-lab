"""One root-authorized +500-iteration continuation of the saved2000 cavity.

Reads the original case/data; never initializes or overwrites a parent file.
Fresh output only, same physics/methods/targets, snapshots2300 and2500.
"""
from pathlib import Path
import argparse
from datetime import datetime
import json
import os
import shutil
import tempfile
import time
import traceback
from cavity_fluent_v001 import (EQUATIONS,WALLS,EXECUTABLE,TARGETS as INITIAL_TARGETS,
    sha256,write_json,available_memory,require_idle,quality_metrics,hdf_terminal)

HELPER=Path(__file__).with_name('cavity_fluent_v001.py')
MAX_ITERATIONS=2500
CHECKPOINTS=(2300,2500)
TARGETS=dict(INITIAL_TARGETS,maximum_total_iterations=2500,checkpoints=list(CHECKPOINTS),
    initial_budget_iterations=2000,additional_authorized_iterations=500,final200_endpoints=[2300,2500],
    original2000_continuity_residual=1.6136954255074055e-6,
    continuation_reason='Original fixed2000 run missed continuity1e-6; one separately authorized500 step continuation, unchanged criteria')


def validate_authorization(path,mesh,parent,out):
    value=json.loads(Path(path).read_text(encoding='utf-8-sig'))
    required=dict(action='continue_cavity_fluent_to2500',approved=True,root='root',allowed_mode='Run')
    if any(value.get(key)!=expected for key,expected in required.items()):raise ValueError('Root cavity continuation Run authorization absent or wrong mode')
    files=dict(source_sha256=Path(__file__),helper_source_sha256=HELPER,mesh_sha256=mesh,
        parent_case_sha256=parent/'cavity-final.cas.h5',parent_data_sha256=parent/'cavity-final.dat.h5')
    for key,file in files.items():
        if value.get(key)!=sha256(file):raise ValueError('Reviewed continuation hash mismatch: '+key)
    if Path(value.get('out','')).resolve()!=out:raise ValueError('Authorized continuation output differs')
    if value.get('additional_iteration_cap')!=500 or value.get('maximum_total_iterations')!=2500:
        raise ValueError('Require exactly500 additional iterations to total2500')
    return value


def native_state(st):
    return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),
        material=st.setup.materials.fluid['air'].get_state(),cell_zones=st.setup.cell_zone_conditions.fluid.get_state(),
        walls={name:st.setup.boundary_conditions.wall[name].get_state() for name in WALLS},references=st.setup.reference_values.get_state(),
        methods=st.solution.methods.get_state(),controls=st.solution.controls.get_state(),residual=st.solution.monitor.residual.get_state(),
        report_definitions=st.solution.report_definitions.get_state(),report_files=st.solution.monitor.report_files.get_state(),
        initialization=st.solution.initialization.get_state())


def physical_readback(st):
    air=st.setup.materials.fluid['air'];lid=st.setup.boundary_conditions.wall['lid'].momentum;values=lid.velocity_components.get_state()
    value=dict(lid_wall_motion=lid.wall_motion(),lid_relative=lid.relative(),lid_velocity_spec=lid.velocity_spec(),
        lid_velocity_components=values,lid_shear_condition=lid.shear_condition(),operating_pressure_Pa=st.setup.general.operating_conditions.operating_pressure(),
        reference_depth_m=st.setup.reference_values.depth(),density=air.density.value(),viscosity=air.viscosity.value(),
        viscous_model=st.setup.models.viscous.model(),energy_enabled=st.setup.models.energy.enabled())
    if (value['lid_wall_motion']!='Moving Wall' or value['lid_relative']!='Absolute' or value['lid_velocity_spec']!='Components'
        or len(values)!=2 or [item['value'] for item in values]!=[1,0] or value['lid_shear_condition']!='No Slip'
        or value['operating_pressure_Pa']!=0 or value['reference_depth_m']!=1 or value['density']!=1 or value['viscosity']!=.01
        or value['viscous_model']!='laminar' or value['energy_enabled']):raise ValueError('Native cavity physical settings changed')
    for name in WALLS[1:]:
        wall=st.setup.boundary_conditions.wall[name].momentum
        if wall.wall_motion()!='Stationary Wall' or wall.shear_condition()!='No Slip':raise ValueError('Stationary wall changed: '+name)
    return value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('mesh','parent','out','authorization'):parser.add_argument('--'+name,required=True,type=Path)
    args=parser.parse_args();mesh=args.mesh.resolve();parent=args.parent.resolve();out=args.out.resolve();repo=Path(__file__).resolve().parents[1]
    if out.exists() or out.drive.upper()!='E:' or out.is_relative_to(repo) or out==parent:raise ValueError('Fresh native E: continuation output outside repo/parent required')
    if mesh.drive.upper()!='E:' or not mesh.is_file() or parent.drive.upper()!='E:':raise ValueError('Existing approved E: mesh and parent required')
    if not EXECUTABLE.is_file():raise FileNotFoundError(EXECUTABLE)
    authorization=validate_authorization(args.authorization,mesh,parent,out)
    actual,parent_residuals=hdf_terminal(parent/'cavity-final.dat.h5')
    if actual!=2000:raise ValueError('Parent must contain exactly2000 actual native iterations')
    setup_path=mesh.parent/'setup.json';geometry_path=mesh.parent/'geometry.npz';setup=json.loads(setup_path.read_text(encoding='utf-8-sig'))
    if setup.get('mesh_sha256')!=sha256(mesh) or setup.get('cells')!=6400:raise ValueError('Approved6400cell shared cavity mesh required')
    parent_files=[parent/name for name in ('cavity-final.cas.h5','cavity-final.dat.h5','final-settings.json','field-monitor-history.out',
        'summary.json','declared-targets.json','code-provenance.json','executed-script.py','authorization.json','native-physical-readback.json')]
    parent_before={str(path):sha256(path) for path in parent_files}
    parent_settings=json.loads((parent/'final-settings.json').read_text(encoding='utf-8-sig'))
    disk=shutil.disk_usage(out.anchor).free;memory=available_memory()
    if disk<10*1024**3 or memory<4*1024**3:raise RuntimeError('Require10GiBE: free and4GiBavailableRAM')
    idle=require_idle();out.mkdir(parents=True);(out/'examples').mkdir()
    temporary=Path('E:/CFD-Learning-Lab/temp/cavity-v001-fluent-extended-v002');temporary.mkdir(parents=True,exist_ok=True)
    for key in ('TEMP','TMP','TMPDIR'):os.environ[key]=str(temporary)
    os.environ['PYTHONDONTWRITEBYTECODE']='1';tempfile.tempdir=str(temporary);os.chdir(out)
    shutil.copyfile(__file__,out/'executed-script.py');shutil.copyfile(HELPER,out/'driver-helper-source.py');shutil.copyfile(setup_path,out/'mesh-setup.json')
    write_json(out/'authorization.json',authorization);write_json(out/'declared-targets.json',TARGETS)
    write_json(out/'parent-input-manifest-before.json',parent_before)
    write_json(out/'code-provenance.json',dict(source_sha256=sha256(__file__),helper_source_sha256=sha256(HELPER),mesh=str(mesh),mesh_sha256=sha256(mesh),
        mesh_setup_sha256=sha256(setup_path),geometry_npz_sha256=sha256(geometry_path),parent_run=str(parent),parent_case_sha256=sha256(parent/'cavity-final.cas.h5'),
        parent_data_sha256=sha256(parent/'cavity-final.dat.h5'),snapshot_before_solver_launch=True,available_ram_before_launch_bytes=memory,
        free_disk_before_launch_bytes=disk,temporary_directory=str(temporary),idle_check=idle))
    write_json(out/'continuation-budget.json',dict(initial_cap=2000,original_final_scaled_residuals=parent_residuals,original_continuity_missed=True,
        separately_authorized_additional_cap=500,new_total_cap=2500,criteria_unchanged=True,never_reinitialize=True,parent_preserved=True))
    events=[];errors=[];checkpoints=[];session=None;loaded=False;saved=False;final_residuals={}
    def update(stage,**values):
        value=dict(stage=stage,time=datetime.now().astimezone().isoformat(),**values)
        events.append(value);write_json(out/'history.json',events);write_json(out/'status.json',value);print(json.dumps(value,default=str),flush=True)
    def optional(label,operation):
        try:operation();return True
        except Exception:
            errors.append(dict(export=label,error=traceback.format_exc()));write_json(out/'export-errors.json',errors);return False
    try:
        update('importing_pyfluent')
        import ansys.fluent.core as pf
        pf.config.examples_path=str(out/'examples');update('launching_single_session',processor_count=2,dimension=2,ui_mode='hidden_gui')
        session=pf.launch_fluent(fluent_path=str(EXECUTABLE),dimension=2,precision='double',processor_count=2,
            ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
        update('solver_connected',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(out/'solver-transcript.txt'))
        st=session.settings;st.file.read_case_data(file_name=str(parent/'cavity-final.cas.h5'));loaded=True
        # Redirect the loaded active report before any iteration, retaining the
        # original report file. Never call initialization in a continuation.
        monitor=st.solution.monitor.report_files['cavity_fields'];monitor.active=False
        monitor.file_name=str(out/'field-monitor-history.out');monitor.active=True
        loaded_state=native_state(st);write_json(out/'parent-loaded-settings.json',loaded_state)
        physics_keys=('general','energy','viscous','material','cell_zones','walls','references','methods','controls')
        changed=[key for key in physics_keys if loaded_state[key]!=parent_settings[key]]
        if changed:raise ValueError('Parent native physics/numerical settings changed on read: '+str(changed))
        readback=physical_readback(st);write_json(out/'native-physical-readback.json',readback)
        residual=st.solution.monitor.residual;residual.options.n_save=3000
        for name in EQUATIONS:
            if residual.equations[name].absolute_criteria()!=1e-6 or residual.equations[name].check_convergence():raise ValueError('Parent residual target/stopping changed')
        check=st.mesh.check();quality=st.mesh.quality()
        for attempt in range(10):
            text=str(check)+'\n'+str(quality)+'\n'+(out/'solver-transcript.txt').read_text(encoding='utf-8',errors='replace')
            try:metrics=quality_metrics(text);break
            except ValueError:
                if attempt==9:raise
                time.sleep(.1)
        write_json(out/'geometry-quality.json',dict(**metrics,native_numeric_gate_pass=True,checked_before_continuation=True))
        write_json(out/'loaded-settings.json',native_state(st))
        st.file.write_case_data(file_name=str(out/'cavity-continuation-start.cas.h5'))
        native_start,_=hdf_terminal(out/'cavity-continuation-start.dat.h5')
        if native_start!=2000:raise ValueError('Loaded parent actual iteration is not2000')
        for endpoint in CHECKPOINTS:
            update('iterating',current_actual_iteration=actual,requested_additional_iterations=endpoint-actual,initial_cap=2000,additional_cap=500,total_cap=2500)
            st.solution.run_calculation.iterate(iter_count=endpoint-actual)
            stem=f'cavity-iter{endpoint:04d}';st.file.write_case_data(file_name=str(out/(stem+'.cas.h5')))
            actual,final_residuals=hdf_terminal(out/(stem+'.dat.h5'))
            if actual!=endpoint:raise RuntimeError('Fixed continuation checkpoint mismatch')
            item=dict(iteration=actual,case=stem+'.cas.h5',data=stem+'.dat.h5',scaled_residuals=final_residuals)
            checkpoints.append(item);write_json(out/'checkpoints.json',checkpoints);update('checkpoint_saved',**item)
        st.file.write_case_data(file_name=str(out/'cavity-final.cas.h5'));actual,final_residuals=hdf_terminal(out/'cavity-final.dat.h5');saved=True
        optional('final_settings',lambda:write_json(out/'final-settings.json',native_state(st)))
        optional('native_residuals_xy',lambda:residual.write(file_name=str(out/'residuals.xy')))
        optional('native_closed_wall_mass_flux_report',lambda:st.results.report.fluxes.mass_flow(zones=list(WALLS),write_to_file=True,file_name=str(out/'wall-mass-balance.txt')))
        def picture(name,field):
            graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1200;graphics.picture.y_resolution=1200
            graphics.picture.invert_background=False;graphics.views.camera.projection='orthographic';graphics.contour.create(name=name)
            contour=graphics.contour[name];contour.field=field;contour.surfaces_list=[];contour.display()
            graphics.views.camera.position=[.5,.5,10.];graphics.views.camera.target=[.5,.5,0.];graphics.views.camera.field=[1.15,1.15]
            graphics.picture.save_picture(file_name=str(out/(name+'.png')))
        optional('native_speed_png',lambda:picture('velocity','velocity-magnitude'));optional('native_pressure_png',lambda:picture('pressure','pressure'))
        optional('transcript_stop',session.transcript.stop)
        unchanged=all(sha256(path)==digest for path,digest in parent_before.items());write_json(out/'parent-input-manifest-after.json',dict(parent_inputs_unchanged=unchanged,inputs=parent_before))
        if not unchanged:raise RuntimeError('A preserved parent input changed')
        summary=dict(kind='lid-driven-cavity-Re100-continuation',actual_iterations=actual,initial_iteration_cap=2000,additional_authorized_iterations=500,
            total_iteration_cap=2500,parent_run=str(parent),parent_inputs_unchanged=unchanged,never_reinitialized=True,
            final_scaled_residuals=final_residuals,native_final_saved=True,native_case='cavity-final.cas.h5',native_data='cavity-final.dat.h5',
            checkpoints=checkpoints,targets=TARGETS,export_errors=errors,physical_validation=False,numerical_checks='Independent saved-file audit required')
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
