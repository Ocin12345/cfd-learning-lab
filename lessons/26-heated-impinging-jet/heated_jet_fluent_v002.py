"""Advance only energy from a saved Fluent impinging-jet case."""
import csv
import json
import os
from pathlib import Path
import shutil
import sys
import traceback

import h5py
import numpy as np
import ansys.fluent.core as pf

from heated_jet_common_v001 import COLD, CP, DT, ENERGY_RESIDUAL, HOT, K, LAST, STRIDE
from jet_common_v001 import ZONES


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + "\n", encoding="utf-8")


def initial_temperature(data):
    with h5py.File(data) as f:
        group = f["results/1/phase-1/cells/SV_T"]
        return [np.asarray(part) for part in group.values()]


def main(source, out):
    assert shutil.disk_usage(out.anchor).free >= 10 * 2**30
    assert source.is_file() and source.with_name(source.name.replace(".cas.", ".dat.")).is_file()
    out.mkdir(parents=True)
    session = None
    status = out / "status.json"

    def mark(stage, **extra):
        dump(status, dict(stage=stage, **extra))
        print(stage, extra, flush=True)

    try:
        mark("launching")
        session = pf.launch_fluent(fluent_path=os.environ["FLUENT_EXECUTABLE"], dimension=2,
                                   precision="double", processor_count=2, ui_mode="hidden_gui",
                                   cwd=str(out), start_timeout=180, cleanup_on_exit=True)
        session.transcript.start(file_name=str(out / "solver-transcript.txt"))
        mark("licensed_fluent_started", version=str(session.get_fluent_version()))
        st = session.settings
        st.file.read_case_data(file_name=str(source))
        energy = st.setup.models.energy
        energy.enabled = True
        for name in ("pressure_work", "kinetic_energy", "viscous_dissipation"):
            setting = getattr(energy, name)
            if setting.is_active():
                setting.set_state(False)
                assert setting() is False
        heating = st.setup.models.viscous.options.viscous_heating
        if heating.is_active():
            heating.set_state(False)
        assert energy.inlet_diffusion.is_active()
        energy.inlet_diffusion = True
        air = st.setup.materials.fluid["air"]
        air.specific_heat.option = "constant"
        air.specific_heat.value = CP
        air.thermal_conductivity.option = "constant"
        air.thermal_conductivity.value = K
        st.setup.boundary_conditions.velocity_inlet["inlet"].thermal.temperature.value = COLD
        for name in ("left", "right"):
            st.setup.boundary_conditions.pressure_outlet[name].thermal.backflow_total_temperature.value = COLD
        plate = st.setup.boundary_conditions.wall["plate"].thermal
        plate.thermal_condition = "Temperature"
        plate.temperature.value = HOT
        ceiling = st.setup.boundary_conditions.wall["ceiling"].thermal
        ceiling.thermal_condition = "Heat Flux"
        ceiling.heat_flux.value = 0.
        st.setup.general.solver.time = "transient"
        sol = st.solution
        sol.methods.transient_formulation = "unsteady-1st-order"
        sol.methods.spatial_discretization.discretization_scheme["temperature"] = "first-order-upwind"
        equations = sol.controls.get_state()["equations"]
        for name in equations:
            sol.controls.equations[name] = name in ("temperature", "energy")
        assert sol.controls.get_state()["equations"] == {"flow": False, "temperature": True}
        residual = sol.monitor.residual
        residual.options.n_save = 12000
        assert residual.equations.get_object_names() == ["energy"]
        residual.equations["energy"].check_convergence = True
        residual.equations["energy"].absolute_criteria = ENERGY_RESIDUAL
        run = sol.run_calculation
        run.transient_controls.type = "Fixed"
        run.transient_controls.method = "User-Specified"
        run.parameters.time_step_size = DT
        run.parameters.max_iter_per_time_step = 30
        sol.initialization.patch.calculate_patch(variable="temperature", value=COLD, cell_zones=["fluid"])
        dump(out / "setup.json", dict(
            energy=energy.get_state(), material=air.get_state(), equations=sol.controls.get_state()["equations"],
            methods=sol.methods.get_state(), plate=plate.get_state(), ceiling=ceiling.get_state(),
            dt=DT, source_case=source.name, last_step=LAST, initialized_temperature_K=COLD))
        st.file.write_case(file_name=str(out / "heated-jet.cas.h5"))
        checkpoints = sorted({*range(0, LAST + 1, STRIDE), LAST - 1})
        last = 0
        with (out / "frames.csv").open("x", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["step", "time_s", "file"])
            for step in checkpoints:
                if step > last:
                    run.dual_time_iterate(time_step_count=step - last, max_iter_per_step=30)
                actual = float(run.transient_controls.flow_time())
                assert abs(actual - step * DT) < 1e-9, (step, actual)
                file = out / f"heat-{step:04d}.dat.h5"
                st.file.write_data(file_name=str(file))
                if step == 0:
                    assert all(np.max(abs(T - COLD)) < 1e-10 for T in initial_temperature(file))
                writer.writerow([step, actual, file.name])
                stream.flush()
                mark("saved_thermal_frame", step=step, time_s=actual)
                last = step
        st.results.report.fluxes.heat_transfer(zones=list(ZONES), write_to_file=True,
                                              file_name=str(out / "final-heat-report.txt"))
        st.results.report.fluxes.mass_flow(zones=list(ZONES), write_to_file=True,
                                          file_name=str(out / "final-mass-report.txt"))
        st.file.write_case_data(file_name=str(out / "heated-jet-final.cas.h5"))
        mark("saved_final", step=LAST, time_s=actual)
    except Exception:
        mark("failed", error=traceback.format_exc())
        raise
    finally:
        if session is not None:
            session.exit()
    mark("closed", step=LAST, time_s=actual)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))

