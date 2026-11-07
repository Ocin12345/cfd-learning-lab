"""Run the zero-incidence wing in a new two-process hidden Fluent session."""
from pathlib import Path
import json,os,shutil,traceback
import ansys.fluent.core as pf
from wing_fluent_v001 import available_memory,sha256,hdf_terminal,force_rows,force_stability,mesh_quality
from wing_zero_common_v001 import U,LIFT,MAX_ITERATIONS
R=Path('E:/CFD-Learning-Lab');P=R/'native/runs/finite-wing-coupled-fluent-v001'
O=R/'native/runs/finite-wing-zero-fluent-v001'
EQ=('continuity','x-velocity','y-velocity','z-velocity')

def write(p,v):p.write_text(json.dumps(v,indent=2,default=str,allow_nan=False)+'\n',encoding='utf-8')
def main():
    assert not O.exists() and Path(os.environ['FLUENT_EXECUTABLE']).is_file()
    assert available_memory()>4*2**30 and shutil.disk_usage(R).free>10*2**30
    O.mkdir();(O/'source').mkdir()
    src=P/'wing-final.cas.h5';h=sha256(src);shutil.copyfile(src,O/'source/wing.cas.h5')
    assert h==sha256(O/'source/wing.cas.h5')==sha256(src)
    shutil.copyfile(__file__,O/'executed-driver.py')
    shutil.copyfile(Path(__file__).with_name('wing_zero_common_v001.py'),O/'wing_zero_common_v001.py')
    write(O/'source-hashes.json',{str(src):h})
    write(O/'declared-settings.json',dict(body_clockwise_pitch_deg=8,freestream_direction=U,lift_direction=LIFT,
        effective_incidence_deg=0,Re=400,density=1.,viscosity=.0025,area=4.,chord=1.,
        maximum_iterations=MAX_ITERATIONS,residual_target=1e-5,mass_target=1e-4,
        final100_force_relative_range_target=.001,force_reference_floor=.001,
        common_mesh_unchanged=True,physical_validation=False))
    session=None;loaded=False
    def status(stage,**kw):
        row=dict(stage=stage,**kw);write(O/'status.json',row);print(json.dumps(row),flush=True)
    try:
        status('launching')
        session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=3,precision='double',
            processor_count=2,ui_mode='hidden_gui',cwd=str(O),start_timeout=180,cleanup_on_exit=True)
        session.transcript.start(file_name=str(O/'solver-transcript.txt'))
        status('licensed_connected',version=str(session.get_fluent_version()))
        st=session.settings;st.file.read_case(file_name=str(O/'source/wing.cas.h5'));loaded=True
        sol=st.solution
        for name in sol.monitor.report_files.get_object_names():sol.monitor.report_files[name].file_name=str(O/(name+'-history.out'))
        for name in ['inlet','farfield']:
            m=st.setup.boundary_conditions.velocity_inlet[name].momentum
            m.velocity_specification_method='Components';m.velocity_components=[dict(value=v) for v in U]
        sol.report_definitions.drag['Cd'].force_vector=U
        sol.report_definitions.lift['Cl'].force_vector=LIFT
        assert sol.methods.p_v_coupling.flow_scheme()=='Coupled'
        assert sol.methods.pseudo_time_method.formulation.coupled_solver()=='off'
        sol.controls.p_v_controls.flow_courant_number=20.
        for n in EQ:
            sol.monitor.residual.equations[n].absolute_criteria=1e-5
            sol.monitor.residual.equations[n].check_convergence=False
        sol.monitor.residual.options.n_save=3000
        sol.controls.advanced.multi_grid.amg_controls.options.verbosity=0
        st.mesh.check();st.mesh.quality()
        write(O/'mesh-quality.json',mesh_quality((O/'solver-transcript.txt').read_text(errors='replace')))
        init=sol.initialization;init.initialization_type='standard'
        for n,v in zip(['x-velocity','y-velocity','z-velocity','pressure'],[*U,0.]):init.defaults[n]=v
        init.standard_initialize()
        write(O/'actual-settings.json',dict(general=st.setup.general.get_state(),viscous=st.setup.models.viscous.get_state(),
            material=st.setup.materials.fluid['air'].get_state(),references=st.setup.reference_values.get_state(),
            boundaries=st.setup.boundary_conditions.get_state(),methods=sol.methods.get_state(),controls=sol.controls.get_state(),
            reports=sol.report_definitions.get_state(),monitors=sol.monitor.report_files.get_state(),residual=sol.monitor.residual.get_state()))
        st.file.write_case_data(file_name=str(O/'wing-initial.cas.h5'))
        checkpoints=[];actual=0;observed=False
        while actual<MAX_ITERATIONS:
            n=min(100 if observed else 200,MAX_ITERATIONS-actual)
            status('iterating',current_iteration=actual,requested=n,confirmation_block=observed)
            sol.run_calculation.iterate(iter_count=n)
            stem=f'wing-{actual+n:04d}'
            st.file.write_case_data(file_name=str(O/(stem+'.cas.h5')))
            new,res=hdf_terminal(O/(stem+'.dat.h5'));assert new==actual+n
            actual=new
            rows=force_rows((O/'forces-history.out').read_text());stability=force_stability(rows,actual)
            entry=dict(iteration=actual,residuals=res,force_stability=stability,case=stem+'.cas.h5',data=stem+'.dat.h5')
            checkpoints.append(entry);write(O/'checkpoints.json',checkpoints);status('saved_checkpoint',**entry)
            meets=max(res.values())<=1e-5 and stability['passed']
            if observed and meets:break
            observed=meets
        for n in EQ:sol.monitor.residual.equations[n].check_convergence=True
        st.file.write_case_data(file_name=str(O/'wing-final.cas.h5'))
        st.results.report.fluxes.mass_flow(zones=['inlet','outlet','farfield','wing'],write_to_file=True,file_name=str(O/'mass-balance.txt'))
        for name,vector in [('drag',U),('lift',LIFT)]:
            st.results.report.forces(option='forces',wall_zones=['wing'],direction_vector=vector,write_to_file=True,file_name=str(O/(name+'-final.txt')))
        write(O/'force-stability.json',stability)
        g=st.results.graphics;g.picture.use_window_resolution=False;g.picture.x_resolution=1600;g.picture.y_resolution=1000
        g.views.camera.projection='orthographic';g.contour.create(name='zero-wing-pressure')
        q=g.contour['zero-wing-pressure'];q.field='pressure';q.surfaces_list=['wing'];q.display()
        g.views.camera.position=[4.,3.,7.];g.views.camera.target=[.5,0,0];g.views.camera.field=[6.8,4.25]
        g.picture.save_picture(file_name=str(O/'wing-zero-pressure-native.png'))
        assert sha256(src)==h
        write(O/'summary.json',dict(iteration=actual,residuals=res,force_stability=stability,numerical_audit_pending=True,physical_validation=False))
        status('saved_final',iteration=actual)
    except Exception:
        status('failed',error=traceback.format_exc())
        if loaded:
            try:session.settings.file.write_case_data(file_name=str(O/'recovery.cas.h5'))
            except Exception:write(O/'recovery-error.json',dict(error=traceback.format_exc()))
        raise
    finally:
        if session:session.exit()
    status('closed',iteration=actual)
if __name__=='__main__':main()
