"""Restart the original 1200-iteration wing with coupled pressure-velocity solution.
One 300-iteration diagnostic, checkpoints every 100. Physics, geometry and
second-order discretization preserved. Pseudo-time off; coupled Courant 20;
explicit pressure/momentum relaxation 0.5. Original failed runs remain intact.
"""
from pathlib import Path
import hashlib,json,os,shutil,sys,traceback
import ansys.fluent.core as pf
from wing_fluent_v001 import hdf_terminal,force_rows,force_stability,mesh_quality,available_memory

R=Path('E:/CFD-Learning-Lab')
PARENT=R/'native/runs/2026-11-01-finite-wing-v003/fluent'
OUT=R/'native/runs/finite-wing-coupled-fluent-v001'
EQ=('continuity','x-velocity','y-velocity','z-velocity')

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False,default=str)+'\n',encoding='utf-8')

def main():
    assert not OUT.exists() and Path(os.environ['FLUENT_EXECUTABLE']).is_file()
    assert shutil.disk_usage(R).free>10*2**30 and available_memory()>4*2**30
    OUT.mkdir();(OUT/'source').mkdir()
    source_hashes={}
    for n in ('wing-final.cas.h5','wing-final.dat.h5'):
        src=PARENT/n;dst=OUT/'source'/n;source_hashes[n]=sha(src)
        shutil.copyfile(src,dst);assert sha(src)==sha(dst)==source_hashes[n]
    write(OUT/'source-hashes.json',source_hashes)
    shutil.copyfile(__file__,OUT/'executed-driver.py')
    write(OUT/'declared-test.json',dict(parent_iteration=1200,additional_iterations=300,method='Coupled',pseudo_time='off',flow_courant=20,residual_target=1e-5,mass_target=1e-4,
        force_range_target=.001,force_samples=100,physical_validation=False))
    session=None;loaded=False
    def status(stage,**kw):
        entry=dict(stage=stage,**kw);write(OUT/'status.json',entry);print(json.dumps(entry),flush=True)
    try:
        status('launching')
        session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=3,precision='double',
            processor_count=2,ui_mode='hidden_gui',cwd=str(OUT),start_timeout=180,cleanup_on_exit=True)
        session.transcript.start(file_name=str(OUT/'solver-transcript.txt'))
        status('licensed_solver_connected',version=str(session.get_fluent_version()))
        st=session.settings;st.file.read_case_data(file_name=str(OUT/'source/wing-final.cas.h5'));loaded=True
        reports=st.solution.monitor.report_files
        for name in reports.get_object_names():
            reports[name].file_name=str(OUT/(name+'-history.out'))
        assert reports['forces'].active() and reports['forces'].frequency()==1
        sol=st.solution
        write(OUT/'loaded-controls.json',sol.controls.get_state())
        write(OUT/'loaded-methods.json',sol.methods.get_state())
        write(OUT/'loaded-setup.json',dict(general=st.setup.general.get_state(),
            material=st.setup.materials.fluid['air'].get_state(),references=st.setup.reference_values.get_state(),
            boundaries=st.setup.boundary_conditions.get_state(),residual=sol.monitor.residual.get_state()))
        st.mesh.check();st.mesh.quality()
        write(OUT/'mesh-quality.json',mesh_quality((OUT/'solver-transcript.txt').read_text(errors='replace')))
        assert sol.methods.p_v_coupling.flow_scheme()=='SIMPLE'
        mg=sol.controls.advanced.multi_grid
        write(OUT/'original-multigrid.json',mg.get_state())
        choices=mg.mg_controls.get_object_names();assert 'pressure' in choices,choices
        pressure=mg.mg_controls['pressure']
        old=float(pressure.termination_criteria());assert old>.001
        # Fixed-length diagnostic; record retained numerical thresholds separately.
        for name in EQ:
            assert sol.monitor.residual.equations[name].absolute_criteria()==1e-5
            sol.monitor.residual.equations[name].check_convergence=False
        sol.monitor.residual.options.n_save=4000
        mg.amg_controls.options.verbosity=1
        st.file.write_case_data(file_name=str(OUT/'restart.cas.h5'))
        checkpoints=[]
        def solve(label,n):
            status('iterating',phase=label,requested=n)
            sol.run_calculation.iterate(iter_count=n)
            stem=OUT/label
            st.file.write_case_data(file_name=str(stem.with_suffix('.cas.h5')))
            it,res=hdf_terminal(stem.with_suffix('.dat.h5'))
            row=dict(phase=label,iteration=it,residual=res,method=sol.methods.p_v_coupling.flow_scheme())
            checkpoints.append(row);write(OUT/'checkpoints.json',checkpoints)
            status('saved',**row)
        sol.methods.p_v_coupling.flow_scheme='Coupled'
        sol.methods.pseudo_time_method.formulation.coupled_solver='off'
        sol.controls.p_v_controls.flow_courant_number=20.
        sol.controls.p_v_controls.explicit_pressure_under_relaxation=.5
        sol.controls.p_v_controls.explicit_momentum_under_relaxation=.5
        write(OUT/'candidate-controls.json',sol.controls.get_state())
        write(OUT/'candidate-methods.json',sol.methods.get_state())
        for expected in (1300,1400,1500):
            solve('coupled-'+str(expected),100)
            assert checkpoints[-1]['iteration']==expected
            assert max(checkpoints[-1]['residual'].values())<1.,'Divergence guard'
        for name in EQ:sol.monitor.residual.equations[name].check_convergence=True
        st.file.write_case_data(file_name=str(OUT/'wing-final.cas.h5'))
        st.results.report.fluxes.mass_flow(zones=['inlet','outlet','farfield','wing'],write_to_file=True,file_name=str(OUT/'mass-balance.txt'))
        for label,direction in [('drag',[1,0,0]),('lift',[0,1,0])]:
            st.results.report.forces(option='forces',wall_zones=['wing'],direction_vector=direction,write_to_file=True,file_name=str(OUT/(label+'-final.txt')))
        rows=force_rows((OUT/'forces-history.out').read_text())
        write(OUT/'force-stability.json',force_stability(rows,1500))
        write(OUT/'final-controls.json',sol.controls.get_state());write(OUT/'final-methods.json',sol.methods.get_state())
        graphics=st.results.graphics;graphics.picture.use_window_resolution=False
        graphics.picture.x_resolution=1500;graphics.picture.y_resolution=1000
        graphics.views.camera.projection='orthographic'
        graphics.contour.create(name='wing-restart-pressure')
        contour=graphics.contour['wing-restart-pressure'];contour.field='pressure';contour.surfaces_list=['wing'];contour.display()
        graphics.views.camera.position=[4.,3.,7.];graphics.views.camera.target=[.5,0,0];graphics.views.camera.field=[6.8,4.25]
        graphics.picture.save_picture(file_name=str(OUT/'wing-pressure-native.png'))
        assert all(sha(PARENT/n)==h for n,h in source_hashes.items())
        status('saved_final',iteration=1500,residuals=checkpoints[-1]['residual'])
    except Exception:
        status('failed',error=traceback.format_exc())
        if loaded:
            try:session.settings.file.write_case_data(file_name=str(OUT/'recovery.cas.h5'))
            except Exception:write(OUT/'recovery-error.json',dict(error=traceback.format_exc()))
        raise
    finally:
        if session:session.exit()
    status('closed',iteration=1500)

if __name__=='__main__':main()
