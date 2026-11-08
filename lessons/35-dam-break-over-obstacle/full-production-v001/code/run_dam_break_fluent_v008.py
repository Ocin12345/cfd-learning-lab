"""Run a source-bound, two-rank hidden Fluent dam-break control or production."""
from pathlib import Path
import argparse,csv,datetime,hashlib,json,os,shutil,traceback
import h5py,numpy as np
import ansys.fluent.core as pf
from analyse_nozzle_fluent_v005 import common_mapping
from dam_break_common_v002 import *
from dam_break_runtime_guard_v001 import check_chunk
from dam_break_admission_v004 import verify,wait_resources,project_path
R=Path('E:/CFD-Learning-Lab');P=R/'local/dam-break-over-obstacle-v001';MESH=R/'native/meshes/dam-break-over-obstacle-9072-v001'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write_json(p,a):p.write_text(json.dumps(a,indent=2,default=str)+'\n',encoding='utf-8')
def verify_fluent(stage):
    proof=verify(stage);receipt=P/'fluent-inputs-v008.json';a=json.loads(receipt.read_text())
    for name,h in a['source_sha256'].items():assert sha(project_path(name))==h,'Fluent source changed: '+name
    control=P/'openfoam-flat-audit.json';a=json.loads(control.read_text());assert a['numerical_checks_pass'] and a['mode']=='flat' and a['steps']==200
    for name,h in a['native_source_sha256'].items():assert sha(project_path(name))==h,'OpenFOAM control changed: '+name
    proof[str(receipt)]=sha(receipt);proof[str(control)]=sha(control);return proof
def water_group(f):
    groups=[phase['cells/SV_VOF'] for phase in f['results/1'].values() if 'cells/SV_VOF' in phase]
    if len(groups)==1:return groups[0]
    groups=[g for g in groups if '/phase-3/' in g.name];assert len(groups)==1;return groups[0]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['flat','collapse','production']);ap.add_argument('out',type=Path);args=ap.parse_args();out=args.out;stage='fluent-'+args.mode
    assert out.resolve().is_relative_to((R/'native/runs/dam-break-over-obstacle-v001').resolve()) and not out.exists()
    assert Path(os.environ['FLUENT_EXECUTABLE']).is_file()
    # All heavy imports above precede this resource check.
    resources=wait_resources(stage);proof=verify_fluent(stage);out.mkdir(parents=True)
    shutil.copyfile(__file__,out/'executed-driver.py')
    for name in ['tank.msh','geometry.npz']:shutil.copyfile(MESH/name,out/name)
    steps=10000 if args.mode=='production' else 400;interval=100 if args.mode=='production' else 20
    g=geometry();expected=initial_fields(g,args.mode);session=None;success=False;last_time=0.;initialized=False
    def mark(state,**kw):
        write_json(out/'status.json',dict(state=state,observed_at=datetime.datetime.now().astimezone().isoformat(),mode=args.mode,**kw));print(state,kw,flush=True)
    try:
        resources=wait_resources(stage);proof=verify_fluent(stage);mark('starting',resources=resources)
        session=pf.launch_fluent(fluent_path=os.environ['FLUENT_EXECUTABLE'],dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(out),start_timeout=180,cleanup_on_exit=True)
        session.transcript.start(file_name=str(out/'solver-transcript.txt'));mark('licensed_fluent_started',version=str(session.get_fluent_version()))
        st=session.settings;st.file.read_mesh(file_name=str(out/'tank.msh'));st.mesh.check();st.mesh.quality()
        st.setup.general.solver.time='transient';st.setup.models.viscous.model='laminar'
        op=st.setup.general.operating_conditions;op.gravity.enable=True;op.gravity.components=[0.,-G]
        mats=st.setup.materials;mats.database.copy_by_name(type='fluid',name='water-liquid')
        for name,(rho,mu) in MATERIALS.items():
            fluid=mats.fluid[name];fluid.density.option='constant';fluid.density.value=rho;fluid.viscosity.option='constant';fluid.viscosity.value=mu
        mp=st.setup.models.multiphase;mp.model='vof';mp.vof_parameters.vof_formulation='explicit';mp.advanced_formulation.implicit_body_force=True
        phases=mp.phases.get_object_names();assert len(phases)==2
        for name,material in zip(phases,MATERIALS):mp.phases[name].material=material
        op.operating_density.method='user-input';op.operating_density.value=0.
        force=mp.phase_interaction.forces
        print('Surface-force state before configuration:',force.get_state(),flush=True)
        tension=force.surface_tension;assert len(tension.get_state())==1
        choices=tension[0].ppt_methods.allowed_values();constant=next(v for v in choices if v.lower()=='constant');tension[0].ppt_methods=constant;tension[0].ppt_constant=SIGMA
        force.surface_tension_model=True
        assert force.surface_tension_model.get_state() is True
        force.surface_tension_model_type='Continuum Surface Force';force.wall_adhesion=False
        def shared_section(boundary,section):
            state=boundary.get_state();print('Boundary active state:',state,flush=True)
            if section in state:return getattr(boundary,section)
            candidates=[key for key,value in state.get('phase',{}).items() if section in value]
            preferred=[key for key in candidates if key.lower()=='mixture']
            selected=preferred[0] if preferred else candidates[0] if len(candidates)==1 else None
            assert selected is not None,(section,candidates,state)
            return getattr(boundary.phase[selected],section)
        wall=st.setup.boundary_conditions.wall['walls'];wall_momentum=shared_section(wall,'momentum');wall_momentum.wall_motion='Stationary Wall';wall_momentum.shear_condition='No Slip'
        outlet=st.setup.boundary_conditions.pressure_outlet['atmosphere'];outlet_momentum=shared_section(outlet,'momentum');outlet_momentum.gauge_pressure.value=0.
        momentum=outlet_momentum.get_state();key=next(k for k in ['backflow_pressure_spec','backflow_pressure_specification'] if k in momentum);getattr(outlet_momentum,key).set_state('Static Pressure')
        outlet_momentum.backflow_dir_spec_method='Normal to Boundary'
        water_boundary=outlet.phase[phases[1]].multiphase;water_boundary.volume_frac_spec_method='Backflow Volume Fraction';water_boundary.backflow_volume_fraction.value=0.
        print('Atmosphere:',outlet.get_state(),flush=True);print('Surface tension:',force.get_state(),flush=True)
        sol=st.solution;method=sol.methods;method.p_v_coupling.flow_scheme='Coupled';
        pseudo=method.get_state().get('pseudo_time_method',{}).get('formulation',{})
        if 'coupled_solver' in pseudo:method.pseudo_time_method.formulation.coupled_solver='off'
        sol.controls.p_v_controls.flow_courant_number=20.
        sol.controls.p_v_controls.explicit_pressure_under_relaxation=.5
        sol.controls.p_v_controls.explicit_momentum_under_relaxation=.5
        scheme=method.spatial_discretization.discretization_scheme
        scheme['mom']='second-order-upwind';scheme['pressure']='body-force-weighted';scheme['mp']='geo-reconstruct'
        run=sol.run_calculation;run.transient_controls.method='User-Specified';run.parameters.time_step_size=DT;run.parameters.max_iter_per_time_step=400
        monitor=sol.monitor.residual;monitor.options.normalize=False;monitor.options.residual_values.scale_residuals=args.mode!='flat'
        targets={'continuity':1e-12,'x-velocity':1e-10,'y-velocity':1e-10} if args.mode=='flat' else {k:1e-4 for k in ['continuity','x-velocity','y-velocity']}
        for name,target in targets.items():monitor.equations[name].absolute_criteria=target;monitor.equations[name].check_convergence=True
        flux=sol.report_definitions.flux;flux.create(name='mdot_atmosphere');report=flux['mdot_atmosphere'];report.report_type='flux-massflow';report.boundaries=['atmosphere'];report.average_over=1;report.create_report_file=False;report.create_report_plot=False
        phase_choices=report.phase.allowed_values();mixture=next(v for v in phase_choices if v.lower()=='mixture');report.phase=mixture
        files=sol.monitor.report_files;files.create(name='atmosphere_mass');mass=files['atmosphere_mass'];mass.report_defs=['mdot_atmosphere'];mass.file_name=str(out/'massflow-history.out')
        choices=mass.frequency_of.allowed_values();print('Report frequency choices:',choices,flush=True)
        frequency=next(v for v in choices if 'time' in v.lower() and 'step' in v.lower());mass.frequency_of=frequency;mass.frequency=1;mass.write_instantaneous_values=True;mass.active=True;mass.print=False
        init=sol.initialization;init.initialization_type='standard'
        for name in ['x-velocity','y-velocity','pressure']:init.defaults[name]=0.
        init.standard_initialize();initialized=True;init.patch.vof_smooth_options.patch_reconstructed_interface=False;init.patch.vof_smooth_options.use_volumetric_smoothing=False
        def patch(variable,value,lower,upper,domain):
            if variable=='mp' and value<=1e-14:return
            sol.cell_registers.create(name='initial_field');region=sol.cell_registers['initial_field'].type;region.option=next(v for v in region.option.allowed_values() if 'hex' in v.lower());region.hexahedron.min_point=lower;region.hexahedron.max_point=upper;region.hexahedron.inside=True
            init.patch.calculate_patch(domain=domain,registers=['initial_field'],variable=variable,value=float(value));sol.cell_registers.delete(name_list=['initial_field'])
        x=g['x_edges'];y=g['y_edges'];height=FLAT_HEIGHT if args.mode=='flat' else WATER_HEIGHT;width=L if args.mode=='flat' else WATER_WIDTH
        ix=int(np.searchsorted(x,width,side='right')-1);iy=int(np.searchsorted(y,height,side='right')-1);ix=min(ix,NX);iy=min(iy,NY)
        patch('mp',1.,[-1e-9,-1e-9],[float(x[ix]),float(y[iy])],phases[1])
        fx=0. if ix==NX else (width-x[ix])/(x[ix+1]-x[ix]);fy=0. if iy==NY else (height-y[iy])/(y[iy+1]-y[iy])
        if fx>1e-14:patch('mp',fx,[float(x[ix]),-1e-9],[float(x[ix+1]),float(y[iy])],phases[1])
        if fy>1e-14:patch('mp',fy,[-1e-9,float(y[iy])],[float(x[ix]),float(y[iy+1])],phases[1])
        if fx>1e-14 and fy>1e-14:patch('mp',fx*fy,[float(x[ix]),float(y[iy])],[float(x[ix+1]),float(y[iy+1])],phases[1])
        # Exact shared physical gauge-pressure input, grouped by row/x fraction.
        for j in range(NY):
            yc=(y[j]+y[j+1])/2;air_p=G*(H-yc)
            value=air_p+999*G*max(height-yc,0) if args.mode=='flat' else air_p
            patch('pressure',value,[-1e-9,float(y[j])],[L+1e-9,float(y[j+1])],'mixture')
            if args.mode!='flat' and height>yc:
                patch('pressure',air_p+999*G*(height-yc),[-1e-9,float(y[j])],[float(x[ix]),float(y[j+1])],'mixture')
                if fx>1e-14:patch('pressure',air_p+999*G*fx*(height-yc),[float(x[ix]),float(y[j])],[float(x[ix+1]),float(y[j+1])],'mixture')
        mark('initial_fields_patched')
        setup=dict(solution_controls=sol.controls.get_state(),mode=args.mode,cells=NCELLS,dimensions_m=[L,H],dt_s=DT,steps=steps,save_interval=interval,surface_tension_N_m=SIGMA,materials={k:mats.fluid[k].get_state() for k in MATERIALS},model=mp.get_state(),operating_conditions=op.get_state(),outlet=outlet.get_state(),wall=wall.get_state(),patch_options=init.patch.vof_smooth_options.get_state(),residual_reporting=dict(options=monitor.options.get_state(),equations=monitor.equations.get_state(),targets=targets),numerical_methods=method.get_state(),run_calculation=run.get_state(),report_definitions=sol.report_definitions.get_state(),report_files=files.get_state(),prerequisite_sha256=proof)
        write_json(out/'setup.json',setup);st.file.write_case(file_name=str(out/'dam.cas.h5'));mapping=common_mapping(out/'dam.cas.h5',g)
        with (out/'frames.csv').open('x',newline='') as stream:
            ledger=csv.writer(stream);ledger.writerow(['step','time_s','file','water_area_m2','mixture_mass_per_m','alpha_min','alpha_max'])
            for step in range(0,steps+1,interval):
                if step:run.dual_time_iterate(time_step_count=interval,max_iter_per_step=400)
                last_time=float(run.transient_controls.flow_time());assert abs(last_time-step*DT)<1e-10
                path=out/f'dam-{step:04d}.dat.h5';st.file.write_data(file_name=str(path));fields=[]
                with h5py.File(path) as f:
                    for group in [water_group(f),f['results/1/phase-1/cells/SV_U'],f['results/1/phase-1/cells/SV_V'],f['results/1/phase-1/cells/SV_P']]:
                        raw=cell_values(group);aligned=np.empty(NCELLS);aligned[mapping]=raw;fields.append(aligned)
                alpha,u,v,pressure=fields;water=float(alpha@g['areas']);rho=1+999*alpha;inventory=float(rho@g['areas'])
                if step==0:
                    da=float(abs(alpha-expected['alpha']).max());dp=float(abs(pressure-expected['pressure']).max());passed=da<1e-10 and dp<1e-8 and abs(water-WATER_AREA)<1e-12 and np.max(np.hypot(u,v))==0
                    write_json(out/'initial-field-check.json',dict(maximum_alpha_difference=da,maximum_pressure_difference_Pa=dp,water_area_m2=water,passed=bool(passed)));assert passed,(da,dp,water)
                assert alpha.min()>=-1e-6 and alpha.max()<=1+1e-6
                ledger.writerow([step,last_time,path.name,water,inventory,float(alpha.min()),float(alpha.max())]);stream.flush();mark('saved_frame',step=step,time_s=last_time,water_area_m2=water)
                if step and args.mode=='production':
                    check=check_chunk((out/'solver-transcript.txt').read_text(errors='replace'),(out/'massflow-history.out').read_text(),completed_steps=step,start_step=step-interval+1,dt=DT,current_mass=inventory,initial_mass=expected['mixture_mass_per_m'],targets=targets)
                    with (out/'runtime-checks.jsonl').open('a',encoding='utf-8') as checks:checks.write(json.dumps(check)+'\n')

        st.file.write_case_data(file_name=str(out/'dam-final.cas.h5'));success=True;mark('saved_final',time_s=last_time,numerical_audit_pending=True)
    except Exception:
        mark('failed',error=traceback.format_exc())
        if session is not None:
            try:
                session.settings.file.write_case(file_name=str(out/'dam-partial.cas.h5'))
                if initialized:session.settings.file.write_data(file_name=str(out/'dam-partial.dat.h5'))
            except Exception as e:print('Partial save unavailable:',str(e),flush=True)
        raise
    finally:
        if session is not None:session.exit()
        if success:mark('closed',time_s=last_time,numerical_audit_pending=True)
if __name__=='__main__':main()