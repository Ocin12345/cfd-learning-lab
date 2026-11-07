"""Hidden Fluent ideal-gas cavity control and short heated qualification."""
from pathlib import Path
import json,os,shutil,traceback
import h5py,numpy as np
import ansys.fluent.core as pf
from plume_production_common_v001 import *
DT=.01; HEATED_STEPS=2000
from fluent_gravity_work_v001 import install_gravity_work
from wing_fluent_v001 import available_memory,mesh_quality
from plot_finite_wing_native_v002 import indexed
OUT=ROOT/'native/runs/hot-patch-plume-20s-v001/fluent'
SOURCE=ROOT/'native/runs/hot-patch-plume-qualification-v001/fluent'
import hashlib
MESH=ROOT/'native/meshes/hot-patch-plume-48x80-v001/cavity.msh'

def saved_stats(path):
    with h5py.File(path) as h:
        a=h['results/1/phase-1/cells'];q={k:indexed(a[v],NC) for k,v in [('u','SV_U'),('v','SV_V'),('T','SV_T'),('p','SV_P'),('rho','SV_DENSITY')]}
    assert all(np.isfinite(v).all() for v in q.values())
    return dict(max_speed_m_s=float(np.hypot(q['u'],q['v']).max()),min_T_K=float(q['T'].min()),max_T_K=float(q['T'].max()),mass_per_depth_kg_m=float(q['rho'].sum()*DX*DY),eos_relative_error=float(np.max(abs(q['rho']*R*q['T']/(q['p']+P0)-1))))
def main():
    assert not OUT.exists() and MESH.is_file() and Path(os.environ['FLUENT_EXECUTABLE']).is_file()
    assert available_memory()>4*2**30 and shutil.disk_usage(ROOT).free>10*2**30
    assert json.loads((ROOT/'exports/figures/2026-11-07-hot-patch-plume-qualification-v001/qualification-audit.json').read_text())['all_saved_field_checks_pass']
    proof=json.loads((ROOT/'native/runs/vertical-energy-control-fluent-v005/energy-check.json').read_text());assert proof['slope_relative_error']<.02
    OUT.mkdir(parents=True);shutil.copyfile(__file__,OUT/'executed-driver.py');write(OUT/'declared-settings.json',dict(specification(),dt_s=DT,heated_steps=HEATED_STEPS,study='20-second heater plume; qualified PISO, continuity 1e-8, cap 1600, original physical conditions and t=0 fields'))
    session=None;loaded=False
    def mark(stage,**kw):write(OUT/'status.json',dict(stage=stage,**kw));print(json.dumps(dict(stage=stage,**kw)),flush=True)
    try:
        mark('launching')
        session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True)
        mark('licensed_connected',version=str(session.get_fluent_version()));session.transcript.start(file_name=str(OUT/'solver-transcript.txt'))
        st=session.settings;st.file.read_case_data(file_name=str(SOURCE/'heated-start.cas.h5'));loaded=True
        assert abs(float(st.solution.run_calculation.transient_controls.flow_time()))<1e-12
        write(OUT/'source-hashes.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [SOURCE/'heated-start.cas.h5',SOURCE/'heated-start.dat.h5']})
        gas=st.setup.materials.fluid['air'];sol=st.solution;residual=sol.monitor.residual
        assert st.setup.general.solver.time()=='transient'
        assert sol.methods.transient_formulation()=='unsteady-1st-order'
        assert sol.methods.p_v_coupling.flow_scheme()=='PISO'
        assert st.setup.cell_zone_conditions.fluid['fluid'].sources.get_state()==json.loads((SOURCE/'heated-gravity-work-source.json').read_text())['readback']
        def state():return dict(general=st.setup.general.get_state(),energy=st.setup.models.energy.get_state(),viscous=st.setup.models.viscous.get_state(),material=gas.get_state(),boundaries=st.setup.boundary_conditions.get_state(),references=st.setup.reference_values.get_state(),methods=sol.methods.get_state(),controls=sol.controls.get_state(),residual=residual.get_state())
        original=json.loads((SOURCE/'heated-settings.json').read_text());loaded_state=state()
        differences=[]
        def equal(a,b,path):
            if isinstance(a,dict) and isinstance(b,dict):
                assert a.keys()==b.keys(),path
                for k in a:equal(a[k],b[k],path+'/'+k)
            elif isinstance(a,list) and isinstance(b,list):
                assert len(a)==len(b),path
                for i,(av,bv) in enumerate(zip(a,b)):equal(av,bv,path+'/'+str(i))
            elif isinstance(a,(float,int)) and not isinstance(a,bool) and isinstance(b,(float,int)) and not isinstance(b,bool):
                assert abs(a-b)<=1e-14*max(abs(a),abs(b),1e-20),(path,a,b)
                if a!=b:differences.append(dict(path=path,loaded=a,baseline=b))
            else:assert a==b,(path,a,b)
        for key in ['general','energy','viscous','material','boundaries','references','methods','controls']:equal(loaded_state[key],original[key],key)
        write(OUT/'restart-settings-check.json',dict(passed=True,numeric_relative_tolerance=1e-14,roundoff_differences=differences))
        write(OUT/'loaded-settings.json',loaded_state)
        assert sol.methods.p_v_coupling.neighbor_correction_itr_count()==1
        sol.controls.under_relaxation['pressure']=1.0
        sol.controls.under_relaxation['mom']=1.0
        write(OUT/'numerical-change.json',dict(pressure_relaxation=[.3,1.],momentum_relaxation=[.7,1.],continuity_limit=1e-8,prior_continuity_limit=1e-6,inner_cap=1600,reason='Retain qualified PISO settings and tighten continuity to limit accumulated mass drift in the longer heating run',source='https://ansyshelp.ansys.com/public/views/secured/corp/v252/en/flu_ug/flu_ug_sec_uns_solve_pvel_1.html',momentum_energy_criteria_unchanged=True,continuity_criterion_tightened=True,physical_parameters_unchanged=True))
        residual.options.criterion_type='absolute';residual.options.n_save=10000
        for name in residual.equations.get_object_names():
            residual.equations[name].check_convergence=True
            residual.equations[name].absolute_criteria=1e-9 if name=='energy' else (1e-8 if name=='continuity' else 1e-6)
        st.mesh.check();st.mesh.quality();write(OUT/'mesh-quality.json',mesh_quality((OUT/'solver-transcript.txt').read_text(errors='replace')))
        control=saved_stats(SOURCE/'heated-start.dat.h5')
        run=sol.run_calculation;run.transient_controls.type='Fixed';run.transient_controls.method='User-Specified';run.parameters.time_step_size=DT;run.parameters.max_iter_per_time_step=INNER
        write(OUT/'heated-gravity-work-source.json',install_gravity_work(st));write(OUT/'heated-settings.json',state());st.file.write_case_data(file_name=str(OUT/'heated-start.cas.h5'));frames=[];previous=0
        for step in [0,1]+list(range(50,HEATED_STEPS+1,50)):
            if (OUT/'STOP_AFTER_SAVE').exists():
                st.file.write_case_data(file_name=str(OUT/'requested-stop.cas.h5'))
                mark('stopped_after_save',heated_time_s=float(run.transient_controls.flow_time()))
                return
            if step>previous:run.dual_time_iterate(time_step_count=step-previous,max_iter_per_step=INNER)
            t=float(run.transient_controls.flow_time());assert abs(t-step*DT)<1e-9
            p=OUT/f'heated-{step:04d}.dat.h5';st.file.write_data(file_name=str(p));row=dict(step=step,time_s=t,file=p.name,**saved_stats(p));frames.append(row);write(OUT/'frames.json',frames);previous=step;mark('heated_state_saved',**row)
            assert COLD-1e-4<=row['min_T_K']<=row['max_T_K']<=HOT+1e-4 and row['eos_relative_error']<=1e-6
            assert abs(row['mass_per_depth_kg_m']/control['mass_per_depth_kg_m']-1)<=1e-5
            if step:
                import re
                raw=(OUT/'solver-transcript.txt').read_text(errors='replace');passed=[]
                for block in re.split(r'Updating solution at time level N',raw)[1:]:
                    event=re.search(r'Flow time = ([0-9.eE+-]+)s, time step = (\d+)',block)
                    if not event:continue
                    rows=re.findall(r'^\s*(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+\d+:',block,re.M)
                    assert rows
                    values=[float(v) for v in rows[-1][1:]];passed.append(dict(step=int(event[2]),time_s=float(event[1]),residuals=values,passed=all(v<=lim for v,lim in zip(values,[1e-8,1e-6,1e-6,1e-9]))))
                write(OUT/'time-step-checks.json',passed)
                assert len(passed)==step and all(x['passed'] for x in passed),'A time step missed its residual targets; result saved before stopping'

        st.file.write_case_data(file_name=str(OUT/'qualification-final.cas.h5'))
        st.results.report.fluxes.heat_transfer(zones=WALLS,write_to_file=True,file_name=str(OUT/'wall-heat-final.txt'))
        graphics=st.results.graphics;graphics.picture.use_window_resolution=False;graphics.picture.x_resolution=1200;graphics.picture.y_resolution=1000
        graphics.contour.create(name='cavity-temperature');c=graphics.contour['cavity-temperature'];c.field='temperature';c.surfaces_list=[];c.options.filled=True;c.options.node_values=False;c.range_options.auto_range=False;c.range_options.minimum=COLD;c.range_options.maximum=HOT;c.color_map.format='%0.2f';c.display()
        graphics.views.camera.projection='orthographic';graphics.views.camera.position=[W/2,H/2,1.];graphics.views.camera.target=[W/2,H/2,0.];graphics.views.camera.field=[.132,.11]
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
