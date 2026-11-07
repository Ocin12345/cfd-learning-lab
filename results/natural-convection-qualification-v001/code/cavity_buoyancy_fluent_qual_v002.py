"""Hidden Fluent ideal-gas cavity control and short heated qualification."""
from pathlib import Path
import json,os,shutil,traceback
import h5py,numpy as np
import ansys.fluent.core as pf
from cavity_buoyancy_common_v002 import *
from wing_fluent_v001 import available_memory,mesh_quality
from plot_finite_wing_native_v002 import indexed
OUT=ROOT/'native/runs/natural-convection-qualification-fluent-v002'
PARENT=ROOT/'native/runs/natural-convection-qualification-fluent-v001'
MESH=ROOT/'native/meshes/natural-convection-qualification-v001/cavity.msh'

def saved_stats(path):
    with h5py.File(path) as h:
        a=h['results/1/phase-1/cells'];q={k:indexed(a[v],N*N) for k,v in [('u','SV_U'),('v','SV_V'),('T','SV_T'),('p','SV_P'),('rho','SV_DENSITY')]}
    assert all(np.isfinite(v).all() for v in q.values())
    return dict(max_speed_m_s=float(np.hypot(q['u'],q['v']).max()),min_T_K=float(q['T'].min()),max_T_K=float(q['T'].max()),mass_per_depth_kg_m=float(q['rho'].sum()*(L/N)**2),eos_relative_error=float(np.max(abs(q['rho']*R*q['T']/(q['p']+P0)-1))))
def main():
    assert not OUT.exists() and MESH.is_file() and Path(os.environ['FLUENT_EXECUTABLE']).is_file()
    assert available_memory()>4*2**30 and shutil.disk_usage(ROOT).free>10*2**30
    OUT.mkdir(parents=True);shutil.copyfile(__file__,OUT/'executed-driver.py');write(OUT/'declared-settings.json',specification())
    session=None;loaded=False
    def mark(stage,**kw):write(OUT/'status.json',dict(stage=stage,**kw));print(json.dumps(dict(stage=stage,**kw)),flush=True)
    try:
        mark('launching')
        session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True)
        mark('licensed_connected',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(OUT/'solver-transcript.txt'))
        st=session.settings
        import hashlib
        source_hashes={}
        (OUT/'source').mkdir()
        for n in ['isothermal-control.cas.h5','isothermal-control.dat.h5']:
            p=PARENT/n;source_hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();shutil.copyfile(p,OUT/'source'/n)
            assert hashlib.sha256((OUT/'source'/n).read_bytes()).hexdigest()==source_hashes[str(p)]
        write(OUT/'source-hashes.json',source_hashes)
        st.file.read_case_data(file_name=str(OUT/'source/isothermal-control.cas.h5'));loaded=True
        sol=st.solution;residual=sol.monitor.residual;gas=st.setup.materials.fluid['air']
        def state():return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),material=gas.get_state(),boundaries=st.setup.boundary_conditions.get_state(),references=st.setup.reference_values.get_state(),methods=sol.methods.get_state(),controls=sol.controls.get_state(),residual=residual.get_state())
        initial=saved_stats(PARENT/'uniform-initial.dat.h5');control=saved_stats(OUT/'source/isothermal-control.dat.h5')
        control['mass_relative_change']=abs(control['mass_per_depth_kg_m']/initial['mass_per_depth_kg_m']-1)
        checks=dict(quiet=control['max_speed_m_s']<=1e-5,temperature=max(abs(control['min_T_K']-T0),abs(control['max_T_K']-T0))<=1e-4,mass=control['mass_relative_change']<=1e-5,eos=control['eos_relative_error']<=1e-6)
        write(OUT/'control-checks.json',dict(initial=initial,control=control,checks=checks,reference_correction='Fluent installed UNIVERSAL_GAS_CONSTANT=8314.47247122, confirmed by saved native p/rho/T times MW; original 1ppm mismatch retained in v001'))
        assert all(checks.values()),'Isothermal control not qualified with documented constant'
        for name in ['top','bottom']:
            st.setup.boundary_conditions.wall[name].thermal.thermal_condition='Heat Flux';st.setup.boundary_conditions.wall[name].thermal.heat_flux.value=0.
        st.setup.boundary_conditions.wall['left'].thermal.temperature.value=HOT;st.setup.boundary_conditions.wall['right'].thermal.temperature.value=COLD
        st.setup.general.solver.time='transient';sol.methods.transient_formulation='unsteady-1st-order';sol.methods.p_v_coupling.flow_scheme='PISO'
        for name in residual.equations.get_object_names():residual.equations[name].check_convergence=True
        run=sol.run_calculation;run.transient_controls.type='Fixed';run.transient_controls.method='User-Specified';run.parameters.time_step_size=DT;run.parameters.max_iter_per_time_step=INNER
        write(OUT/'heated-settings.json',state());st.file.write_case_data(file_name=str(OUT/'heated-start.cas.h5'));frames=[];previous=0
        for step in [0,1,5,10,20,30,40,50]:
            if step>previous:run.dual_time_iterate(time_step_count=step-previous,max_iter_per_step=INNER)
            t=float(run.transient_controls.flow_time());assert abs(t-step*DT)<1e-9
            p=OUT/f'heated-{step:04d}.dat.h5';st.file.write_data(file_name=str(p));row=dict(step=step,time_s=t,file=p.name,**saved_stats(p));frames.append(row);write(OUT/'frames.json',frames);previous=step;mark('heated_state_saved',**row)
        st.file.write_case_data(file_name=str(OUT/'qualification-final.cas.h5'))
        st.results.report.fluxes.heat_transfer(zones=WALLS,write_to_file=True,file_name=str(OUT/'wall-heat-final.txt'))
        graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1200;graphics.picture.y_resolution=1000
        graphics.contour.create(name='cavity-temperature');c=graphics.contour['cavity-temperature'];c.field='temperature';names=session.fields.field_info.get_surfaces_info();write(OUT/'surface-info.json',names);c.surfaces_list=[n for n in names if 'interior' in n or 'fluid' in n];assert c.surfaces_list();c.display();graphics.views.auto_scale()
        graphics.picture.save_picture(file_name=str(OUT/'cavity-temperature-native.png'))
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
if __name__=='__main__':main()
