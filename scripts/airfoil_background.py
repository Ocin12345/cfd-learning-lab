from pathlib import Path

import sys,json,time,traceback,os

root=Path(__file__).resolve().parents[1]

import ansys.fluent.core as pf

base=Path(os.environ['CFD_OUTPUT_ROOT']).resolve()/'ansys-airfoil'

base.mkdir(parents=True,exist_ok=True)

status=base/'status.json'

def update(stage,**kw):

    status.write_text(json.dumps(dict(stage=stage,updated=time.strftime('%Y-%m-%d %H:%M:%S'),**kw),indent=2))

    print(stage,kw,flush=True)

s=None

try:

    update('starting',ui_mode='hidden_gui',solver_processes=2)

    s=pf.launch_fluent(fluent_path=os.environ.get('FLUENT_EXECUTABLE'),dimension=2,precision='double',processor_count=2,ui_mode='hidden_gui',cwd=str(base),start_timeout=180,cleanup_on_exit=True)

    for label in ['alpha0','alpha5']:

        out=base/label

        update('solving',mesh=label)

        s.transcript.start(file_name=str(out/'solver-transcript.txt'))

        st=s.settings

        st.file.read_mesh(file_name=str(out/'airfoil.msh'))

        st.mesh.check()

        st.mesh.quality()

        st.setup.models.viscous.model='laminar'

        st.setup.materials.fluid['air'].density.value=1

        st.setup.materials.fluid['air'].viscosity.value=.01

        bc=st.setup.boundary_conditions.velocity_inlet['inlet'].momentum

        bc.velocity_specification_method='Components'

        bc.velocity_components=[{'value':1},{'value':0}]

        for k,v in dict(density=1,velocity=1,length=1,area=1).items():setattr(st.setup.reference_values,k,v)

        for name in ['continuity','x-velocity','y-velocity']:

            st.solution.monitor.residual.equations[name].absolute_criteria=1e-7

        st.solution.initialization.hybrid_initialize()

        st.solution.run_calculation.iterate(iter_count=1000)

        st.results.report.fluxes.mass_flow(zones=['inlet','outlet'],write_to_file=True,file_name=str(out/'mass-balance.txt'))

        st.results.report.forces(option='forces',wall_zones=['airfoil'],direction_vector=[1,0,0],write_to_file=True,file_name=str(out/'drag.txt'))

        st.results.report.forces(option='forces',wall_zones=['airfoil'],direction_vector=[0,1,0],write_to_file=True,file_name=str(out/'lift.txt'))

        g=st.results.graphics

        g.picture.use_window_resolution=False

        g.picture.x_resolution=1800

        g.picture.y_resolution=1100

        g.picture.invert_background=False

        g.views.camera.projection='orthographic'

        for name,field in [('velocity','velocity-magnitude'),('pressure','pressure')]:

            if name not in g.contour.get_object_names():g.contour.create(name=name)

            g.contour[name].field=field

            g.contour[name].surfaces_list=[]

            g.contour[name].display()

            g.views.camera.position=[.55,0,50]

            g.views.camera.target=[.55,0,0]

            g.views.camera.field=[2.8,1.7]

            g.picture.save_picture(file_name=str(out/f'{name}.png'))

        st.file.write_case_data(file_name=str(out/'airfoil-Re100.cas.h5'))

        s.transcript.stop()

        update('mesh_saved',mesh=label)

    update('completed',meshes=['alpha0','alpha5'])

except Exception:

    update('failed',error=traceback.format_exc())

    raise

finally:

    if s is not None:s.exit()
