"""Hidden Fluent ideal-gas cavity control and short heated qualification."""
from pathlib import Path
import json,os,shutil,traceback
import h5py,numpy as np
import ansys.fluent.core as pf
from cavity_spatial_common_v001 import *
from cavity_prelaunch_v001 import prelaunch
from cavity_step_audit_v001 import audit_steps
from fluent_gravity_work_v001 import install_gravity_work
from wing_fluent_v001 import available_memory,mesh_quality
from plot_finite_wing_native_v002 import indexed
OUT=ROOT/'native/runs/showcase-convection-spatial-v001/fluent-qualification-v001'
MESH=ROOT/'native/meshes/showcase-convection-spatial-v001/fine-v001/cavity.msh'

def saved_stats(path):
    with h5py.File(path) as h:
        a=h['results/1/phase-1/cells'];q={k:indexed(a[v],N*N) for k,v in [('u','SV_U'),('v','SV_V'),('T','SV_T'),('p','SV_P'),('rho','SV_DENSITY')]}
    assert all(np.isfinite(v).all() for v in q.values())
    return dict(max_speed_m_s=float(np.hypot(q['u'],q['v']).max()),min_T_K=float(q['T'].min()),max_T_K=float(q['T'].max()),mass_per_depth_kg_m=float(q['rho'].sum()*(L/N)**2),eos_relative_error=float(np.max(abs(q['rho']*R*q['T']/(q['p']+P0)-1))))
def main():
    assert not OUT.exists() and MESH.is_file() and Path(os.environ['FLUENT_EXECUTABLE']).is_file()
    assert available_memory()>4*2**30 and shutil.disk_usage(ROOT).free>10*2**30
    proof=json.loads((ROOT/'native/runs/vertical-energy-control-fluent-v005/energy-check.json').read_text());assert proof['slope_relative_error']<.02
    admission=prelaunch('fluent-qualification-inputs-v001.json')
    OUT.mkdir(parents=True);write(OUT/'launch-admission.json',admission);shutil.copyfile(__file__,OUT/'executed-driver.py');write(OUT/'declared-settings.json',specification())
    session=None;loaded=False
    def mark(stage,**kw):write(OUT/'status.json',dict(stage=stage,**kw));print(json.dumps(dict(stage=stage,**kw)),flush=True)
    try:
        mark('launching')
        session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True)
        mark('licensed_connected',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(OUT/'solver-transcript.txt'))
        st=session.settings;st.file.read_mesh(file_name=str(MESH));loaded=True
        st.setup.general.solver.type='pressure-based';st.setup.general.solver.time='steady';st.setup.models.viscous.model='laminar';st.setup.models.energy.enabled=True
        gas=st.setup.materials.fluid['air'];gas.density.option='ideal-gas';gas.molecular_weight.value=MW
        for attr,value in [('specific_heat',CP),('viscosity',MU),('thermal_conductivity',K)]:
            p=getattr(gas,attr);p.option='constant';p.value=value
        st.setup.cell_zone_conditions.fluid['fluid'].general.material='air'
        op=st.setup.general.operating_conditions;op.operating_pressure=P0;op.gravity.enable=True;op.gravity.components=[0.,-G]
        op.operating_density.enable=True;op.operating_density.value=0.
        for name,value in [('pressure_work',True),('kinetic_energy',True),('viscous_dissipation',False)]:
            setting=getattr(st.setup.models.energy,name)
            if setting.is_active():setting.set_state(value)
        write(OUT/'gravity-work-source.json',install_gravity_work(st))
        for name in WALLS:
            wall=st.setup.boundary_conditions.wall[name];wall.momentum.wall_motion='Stationary Wall';wall.momentum.shear_condition='No Slip'
            wall.thermal.thermal_condition='Temperature';wall.thermal.temperature.value=T0
        st.setup.reference_values.depth=1.
        sol=st.solution;sol.methods.p_v_coupling.flow_scheme='SIMPLE'
        sol.methods.spatial_discretization.gradient_scheme='least-square-cell-based'
        sol.methods.spatial_discretization.discretization_scheme={'pressure':'body-force-weighted','mom':'second-order-upwind','temperature':'second-order-upwind'}
        residual=sol.monitor.residual;residual.options.n_save=10000;residual.options.criterion_type='absolute'
        for name in residual.equations.get_object_names():
            residual.equations[name].absolute_criteria=1e-9 if name=='energy' else 1e-6
            residual.equations[name].check_convergence=False
        st.mesh.check();st.mesh.quality();write(OUT/'mesh-quality.json',mesh_quality((OUT/'solver-transcript.txt').read_text(errors='replace')))
        assert json.loads((OUT/'mesh-quality.json').read_text())['minimum_orthogonal_quality']>=.2
        def state():return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),material=gas.get_state(),boundaries=st.setup.boundary_conditions.get_state(),references=st.setup.reference_values.get_state(),methods=sol.methods.get_state(),controls=sol.controls.get_state(),residual=residual.get_state())
        write(OUT/'control-settings.json',state())
        init=sol.initialization;init.initialization_type='standard'
        for name,value in [('x-velocity',0.),('y-velocity',0.),('pressure',0.),('temperature',T0)]:init.defaults[name]=value
        init.standard_initialize();st.file.write_case_data(file_name=str(OUT/'uniform-initial.cas.h5'));initial=saved_stats(OUT/'uniform-initial.dat.h5')
        mark('isothermal_control_running',iterations=CONTROL_ITERATIONS);sol.run_calculation.iterate(iter_count=CONTROL_ITERATIONS)
        st.file.write_case_data(file_name=str(OUT/'isothermal-control.cas.h5'));control=saved_stats(OUT/'isothermal-control.dat.h5')
        control['mass_relative_change']=abs(control['mass_per_depth_kg_m']/initial['mass_per_depth_kg_m']-1)
        checks=dict(quiet=control['max_speed_m_s']<=1e-5,temperature=max(abs(control['min_T_K']-T0),abs(control['max_T_K']-T0))<=1e-4,mass=control['mass_relative_change']<=1e-5,eos=control['eos_relative_error']<=1e-6)
        write(OUT/'control-checks.json',dict(initial=initial,control=control,checks=checks))
        assert all(checks.values()),'Isothermal control not qualified; saved for diagnosis'
        for name in ['top','bottom']:
            st.setup.boundary_conditions.wall[name].thermal.thermal_condition='Heat Flux';st.setup.boundary_conditions.wall[name].thermal.heat_flux.value=0.
        st.setup.boundary_conditions.wall['left'].thermal.temperature.value=HOT;st.setup.boundary_conditions.wall['right'].thermal.temperature.value=COLD
        st.setup.general.solver.time='transient';sol.methods.transient_formulation='unsteady-1st-order';sol.methods.p_v_coupling.flow_scheme='PISO'
        for name in residual.equations.get_object_names():
            residual.equations[name].check_convergence=True
            residual.equations[name].absolute_criteria=1e-9 if name=='energy' else 1e-6
        run=sol.run_calculation;run.transient_controls.type='Fixed';run.transient_controls.method='User-Specified';run.parameters.time_step_size=DT;run.parameters.max_iter_per_time_step=INNER
        write(OUT/'heated-gravity-work-source.json',install_gravity_work(st));write(OUT/'heated-settings.json',state());
        from audit_cavity_fluent_qualification_v001 import check_settings,check_initial
        write(OUT/'matched-baseline-settings.json',check_settings(state()))
        st.file.write_case_data(file_name=str(OUT/'heated-start.cas.h5'))
        write(OUT/'initial-comparison.json',check_initial(OUT))
        frames=[];previous=0
        for step in [0,1,10,25,50,75,100]:
            if step>previous:run.dual_time_iterate(time_step_count=step-previous,max_iter_per_step=INNER)
            t=float(run.transient_controls.flow_time());assert abs(t-step*DT)<1e-9
            st.results.report.fluxes.mass_flow(zones=WALLS,write_to_file=True,file_name=str(OUT/f'wall-mass-{step:04d}.txt'))
            st.results.report.fluxes.heat_transfer(zones=WALLS,write_to_file=True,file_name=str(OUT/f'wall-heat-{step:04d}.txt'))
            p=OUT/f'heated-{step:04d}.dat.h5';st.file.write_data(file_name=str(p));row=dict(step=step,time_s=t,file=p.name,**saved_stats(p));frames.append(row);write(OUT/'frames.json',frames);previous=step;mark('heated_state_saved',**row)
            assert COLD-1e-4<=row['min_T_K']<=row['max_T_K']<=HOT+1e-4 and row['eos_relative_error']<=1e-6
            assert abs(row['mass_per_depth_kg_m']/control['mass_per_depth_kg_m']-1)<=1e-5
            if step:
                passed=audit_steps((OUT/'solver-transcript.txt').read_text(errors='replace'),step,DT)
                write(OUT/'time-step-checks.json',passed)
                assert len(passed)==step and all(x['passed'] for x in passed),'A time step missed its residual targets; result saved before stopping'

        st.file.write_case_data(file_name=str(OUT/'qualification-final.cas.h5'))
        st.results.report.fluxes.heat_transfer(zones=WALLS,write_to_file=True,file_name=str(OUT/'wall-heat-final.txt'))
        write(OUT/'final-settings.json',state());mark('saved_final',heated_time_s=HEATED_STEPS*DT)
    except Exception:
        mark('failed',error=traceback.format_exc())
        if loaded:
            try:session.settings.file.write_case_data(file_name=str(OUT/'recovery.cas.h5'))
            except Exception:pass
        raise
    finally:
        if session:session.exit()
    mark('closed',heated_time_s=HEATED_STEPS*DT)
    from audit_cavity_fluent_qualification_v001 import audit
    audit(OUT)
if __name__=='__main__':main()
