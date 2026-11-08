"""Run the 2-D inviscid compression corner in one hidden Fluent session."""
import json
from math import tan
import os
from pathlib import Path
import shutil
import sys
import time
import traceback
from datetime import datetime
from wedge_prelaunch_v001 import prelaunch,verify
from wing_fluent_v001 import available_memory
from wedge_endpoint_v001 import evaluate
S=Path(__file__).resolve().parent
ROOT=Path("E:/CFD-Learning-Lab")
MESH=ROOT/"native/meshes/showcase-mach2-wedge-spatial-v001/fine-v001"

import numpy as np
import ansys.fluent.core as pf

from wedge_common_fine_v001 import (CP, CORNER, H, L, LIMITS, MACH, MW, NX, NY,
                              P, P0, RESIDUALS, T, T0, THETA, U, reference)


def json_out(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=str) + "\n", encoding="utf-8")





def main(out, support):
    sys.path.insert(0, str(support))
    from nozzle_fluent_v004 import hdf_terminal, quality_metrics

    if shutil.disk_usage(out.anchor).free < 10 * 2**30:
        raise RuntimeError("Require 10 GiB free on the run drive")
    out.mkdir(parents=True)
    verify("fluent-fine-inputs-v001.json")
    shutil.copy2(Path(__file__),out/"executed-script.py")
    for name in ("wedge.msh","geometry.npz"):shutil.copy2(MESH/name,out/name)
    mesh=out/"wedge.msh"
    json_out(out / "reference.json", reference())
    session = None
    checkpoints = []
    initialized=False

    def stage(label, **kwargs):
        record=dict(stage=label,time=datetime.now().astimezone().isoformat(),pid=os.getpid(),**kwargs)
        json_out(out/"status.json",record)
        json_out(S/"fluent-fine-v001-status.json",record)
        print(label, kwargs, flush=True)

    try:
        while available_memory()<4*2**30 or shutil.disk_usage(ROOT).free<10*2**30:stage("waiting_resources");time.sleep(30)
        resources=prelaunch("fluent-fine-inputs-v001.json");json_out(out/"launch-admission.json",resources)
        pf.config.examples_path=str(out/"examples")
        stage("launching")
        session = pf.launch_fluent(fluent_path=os.environ["FLUENT_EXECUTABLE"],
                                   dimension=2, precision="double", processor_count=2,
                                   ui_mode="hidden_gui", cwd=str(out), start_timeout=180,
                                   cleanup_on_exit=True)
        stage("licensed_fluent_started", version=str(session.get_fluent_version()))
        transcript = out / "solver-transcript.txt"
        session.transcript.start(file_name=str(transcript))
        st = session.settings
        st.file.read_mesh(file_name=str(mesh))
        st.setup.general.solver.type = "density-based-implicit"
        st.setup.general.solver.two_dim_space = "planar"
        st.setup.general.solver.time = "steady"
        st.setup.models.energy.enabled = True
        st.setup.models.viscous.model = "inviscid"
        st.setup.general.operating_conditions.operating_pressure = 0.
        gas = st.setup.materials.fluid["air"]
        gas.density.option = "ideal-gas"
        gas.specific_heat.option = "constant"
        gas.specific_heat.value = CP
        gas.molecular_weight.value = MW
        if gas.thermal_conductivity.is_active():
            gas.thermal_conductivity.option = "constant"
            gas.thermal_conductivity.value = 0.
        st.setup.cell_zone_conditions.fluid["fluid"].general.material = "air"
        st.setup.reference_values.depth = 1.
        inlet = st.setup.boundary_conditions.pressure_inlet["inlet"]
        inlet.momentum.gauge_total_pressure.value = P0
        inlet.momentum.supersonic_or_initial_gauge_pressure.value = P
        inlet.momentum.direction_specification_method = "Normal to Boundary"
        inlet.thermal.total_temperature.value = T0
        outlet = st.setup.boundary_conditions.pressure_outlet["outlet"]
        outlet.momentum.gauge_pressure.value = P
        outlet.thermal.backflow_total_temperature.value = T0
        method = st.solution.methods
        method.p_v_coupling.flux_type = "Roe-FDS"
        method.spatial_discretization.gradient_scheme = "least-square-cell-based"
        disc = method.spatial_discretization.discretization_scheme
        disc["amg-c"] = "first-order-upwind"
        st.solution.controls.courant_number = 1.
        residuals = st.solution.monitor.residual
        residuals.options.n_save = 2400
        for name, target in RESIDUALS.items():
            res = residuals.equations[name]
            res.absolute_criteria = target
            res.check_convergence = False
        st.mesh.check()
        st.mesh.quality()
        for attempt in range(15):
            try:
                quality = quality_metrics(transcript.read_text(encoding="utf-8", errors="replace"))
                break
            except ValueError:
                if attempt == 14:
                    raise
                time.sleep(.2)
        assert quality["minimum_orthogonal_quality"] >= LIMITS["minimum_orthogonal_quality"]
        json_out(out / "native-mesh-quality.json", quality)
        init = st.solution.initialization
        init.initialization_type = "standard"
        for key, value in (("x-velocity", U), ("y-velocity", 0.), ("pressure", P), ("temperature", T)):
            init.defaults[key] = value
        init.standard_initialize()
        initialized=True
        actual = 0
        for end in (120, 320, 520, 720, 920, 1120, 1320, 1520):
            if actual == 120:
                disc["amg-c"] = "second-order-upwind"
                st.solution.controls.courant_number = 3.
            stage("iterating", from_iteration=actual, to_iteration=end, scheme=disc["amg-c"]())
            st.solution.run_calculation.iterate(iter_count=end - actual)
            stem = f"wedge-{end:04d}"
            st.file.write_case_data(file_name=str(out / (stem + ".cas.h5")))
            actual, residual = hdf_terminal(out / (stem + ".dat.h5"))
            assert actual == end
            st.results.report.fluxes.mass_flow(zones=["inlet", "outlet", "ramp", "top"],
                                              write_to_file=True,
                                              file_name=str(out / (stem + "-mass.txt")))
            entry = dict(iteration=actual, case=stem + ".cas.h5", data=stem + ".dat.h5",
                         mass_report=stem + "-mass.txt", residuals=residual, scheme=disc["amg-c"]())
            checkpoints.append(entry)
            json_out(out / "checkpoints.json", checkpoints)
            stage("checkpoint_saved", **entry)
            if actual>=520:
                with np.load(out/"geometry.npz") as geo:check=evaluate(out,entry,checkpoints[-2],geo)
                json_out(out/f"checkpoint-{actual:04d}-audit.json",check)
                if check["numerical_endpoint_pass"]:break
        json_out(out / "loaded-settings.json", dict(
            general=st.setup.general.get_state(), viscous=st.setup.models.viscous.get_state(),
            gas=gas.get_state(), inlet=inlet.get_state(), outlet=outlet.get_state(),
            methods=method.get_state(), reference_mach=MACH, reference=reference(),
            controls=st.solution.controls.get_state(),walls={name:st.setup.boundary_conditions.wall[name].get_state() for name in ("top","ramp")}))
        stage("saved_final", final=checkpoints[-1])
    except Exception:
        if session is not None and initialized:
            try:session.settings.file.write_case_data(file_name=str(out/"recovery.cas.h5"))
            except Exception as error:print("Recovery save error:",str(error),flush=True)
        stage("failed", error=traceback.format_exc())
        raise
    finally:
        if session is not None:
            session.exit()
    stage("closed", final=checkpoints[-1])


if __name__ == "__main__":
    out=Path(sys.argv[1]);main(out,Path(sys.argv[2]))
    from audit_wedge_fluent_fine_v001 import audit
    result=audit(out)
    json_out(S/"fluent-fine-v001-status.json",dict(stage="fluent_endpoint_pass" if result["numerical_endpoint_pass"] else "fluent_endpoint_numerical_miss",time=datetime.now().astimezone().isoformat(),pid=os.getpid(),audit=str(out/"endpoint-audit.json"),failed_checks=result["failed_checks"]))

