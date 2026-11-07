"""Qualify/run the refined sloshing case in a hidden two-process Fluent session."""
import csv,json,os,shutil,sys,time,traceback
from pathlib import Path
import numpy as np,h5py
import ansys.fluent.core as pf
from analyse_nozzle_fluent_v005 import common_mapping
from sloshing_refinement_common_v001 import NX,NY,DX,DT,L,H,MATERIALS,WALLS,LAST_STEP,FRAME_INTERVAL,EXPECTED_AREA,RESIDUAL_TARGET,INITIAL_TILT,MODE,RAW_TARGETS,write_json,make_mesh,alpha_group,cell_values

def main(out):
    from sloshing_refinement_admission_v001 import verify_prerequisites, wait_resources
    qualification_hashes = verify_prerequisites("fluent-" + MODE)
    wait_resources("fluent-" + MODE)
    out.mkdir(parents=True)
    shutil.copyfile(__file__, out / "executed-driver.py")
    mesh = make_mesh(out)
    session = None
    status = out / "status.json"

    def mark(stage, **kwargs):
        write_json(status, dict(stage=stage, **kwargs))
        print(stage, kwargs, flush=True)

    try:
        # The memory wait and mesh preparation can outlive qualification validity.
        # Refresh native-source proof immediately before any solver launch.
        wait_resources("fluent-" + MODE)
        qualification_hashes = verify_prerequisites("fluent-" + MODE)
        mark("starting")
        session = pf.launch_fluent(fluent_path=os.environ["FLUENT_EXECUTABLE"],
                                   dimension=2, precision="double", processor_count=2,
                                   ui_mode="hidden_gui", cwd=str(out), start_timeout=180,
                                   cleanup_on_exit=True)
        session.transcript.start(file_name=str(out / "solver-transcript.txt"))
        mark("licensed_fluent_started", version=str(session.get_fluent_version()))
        st = session.settings
        st.file.read_mesh(file_name=str(mesh))
        st.mesh.check()
        st.mesh.quality()
        st.setup.general.solver.time = "transient"
        st.setup.models.viscous.model = "laminar"
        op = st.setup.general.operating_conditions
        op.gravity.enable = True
        op.gravity.components = [0., -9.81]
        mats = st.setup.materials
        mats.database.copy_by_name(type="fluid", name="water-liquid")
        for name, (rho, mu) in MATERIALS.items():
            fluid = mats.fluid[name]
            fluid.density.option = "constant"
            fluid.density.value = rho
            fluid.viscosity.option = "constant"
            fluid.viscosity.value = mu
        mp = st.setup.models.multiphase
        mp.model = "vof"
        mp.vof_parameters.vof_formulation = "explicit"
        mp.advanced_formulation.implicit_body_force = True
        phase_names = mp.phases.get_object_names()
        assert len(phase_names) == 2, phase_names
        for phase, material in zip(phase_names, MATERIALS):
            mp.phases[phase].material = material
        print("Phases:", mp.phases.get_state(), flush=True)
        print("Gravity:", op.get_state(), flush=True)
        sol = st.solution
        method = sol.methods
        method.p_v_coupling.flow_scheme = "PISO"
        scheme = method.spatial_discretization.discretization_scheme
        print("Discretization:", scheme.get_state(), flush=True)
        scheme["mom"] = "second-order-upwind"
        scheme["pressure"] = "body-force-weighted"
        run = sol.run_calculation
        run.transient_controls.type = "Fixed"
        run.transient_controls.method = "User-Specified"
        run.parameters.time_step_size = DT
        run.parameters.max_iter_per_time_step = 400
        monitor = sol.monitor.residual
        monitor.options.normalize = False
        monitor.options.residual_values.scale_residuals = MODE != "flat"
        targets = RAW_TARGETS if MODE == "flat" else {name:RESIDUAL_TARGET for name in ("continuity","x-velocity","y-velocity")}
        for name, target in targets.items():
            residual = monitor.equations[name]
            residual.absolute_criteria = target
            residual.check_convergence = True

        init = sol.initialization
        init.initialization_type = "standard"
        for name in ("x-velocity", "y-velocity", "pressure"):
            init.defaults[name] = 0.
        init.standard_initialize()
        # Fractions are already integrated from the declared interface at each cell.
        # Avoid reconstructing/smoothing the small register boxes a second time.
        init.patch.vof_smooth_options.patch_reconstructed_interface = False
        init.patch.vof_smooth_options.use_volumetric_smoothing = False
        patch_options = init.patch.vof_smooth_options.get_state()
        assert patch_options["patch_reconstructed_interface"] is False
        assert patch_options["use_volumetric_smoothing"] is False
        def patch_rectangle(lower, upper, value):
            sol.cell_registers.create(name="initial_water")
            region = sol.cell_registers["initial_water"].type
            region.option = next(option for option in region.option.allowed_values()
                                 if "hex" in option.lower())
            region.hexahedron.min_point = lower
            region.hexahedron.max_point = upper
            region.hexahedron.inside = True
            init.patch.calculate_patch(domain=phase_names[1], registers=["initial_water"],
                                       variable="mp", value=float(value))
            sol.cell_registers.delete(name_list=["initial_water"])

        # Start with the complete horizontal base, then add positive cut-cell fractions.
        base_height = 0.05 - INITIAL_TILT * L / 2
        base_rows = int(np.floor(base_height / DX + 1e-12))
        patch_rectangle([-DX, -DX], [L + DX, base_rows * DX], 1.0)
        for column in range(NX):
            x_mid = (column + 0.5) * DX
            surface_height = 0.05 + INITIAL_TILT * (x_mid - L / 2)
            full_rows = int(np.floor(surface_height / DX + 1e-12))
            fraction = surface_height / DX - full_rows
            if full_rows > base_rows:
                patch_rectangle([column * DX, base_rows * DX],
                                [(column + 1) * DX, full_rows * DX], 1.0)
            if fraction > 1e-12:
                patch_rectangle([column * DX, full_rows * DX],
                                [(column + 1) * DX, (full_rows + 1) * DX], fraction)
        mark("water_patched", initial_tilt=INITIAL_TILT, mean_height_m=0.05)
        write_json(out / "setup.json", dict(
            dimensions_m=[L,H], cells=[NX,NY], dx_m=DX, time_step_s=DT,
            case_mode=MODE,
            total_steps=LAST_STEP, frame_interval=FRAME_INTERVAL,
            qualification_audit_sha256=qualification_hashes,
            initial_mean_water_height_m=0.05, initial_tilt=INITIAL_TILT,
            expected_water_area_m2=EXPECTED_AREA, wall_patches=WALLS,
            model=mp.get_state(), operating_conditions=op.get_state(),
            vof_patch_options=patch_options,
            materials={name:mats.fluid[name].get_state() for name in MATERIALS},
            residual_reporting=dict(options=monitor.options.get_state(),
                                    equations=monitor.equations.get_state(),
                                    absolute_targets=targets),
            numerical_methods=method.get_state(), run_calculation=run.get_state()))
        case_path = out / "tank.cas.h5"
        st.file.write_case(file_name=str(case_path))
        with (out / "frames.csv").open("x", newline="", encoding="utf-8") as stream:
            ledger = csv.writer(stream)
            ledger.writerow(["step", "time_s", "file", "water_area_m2", "alpha_min", "alpha_max"])
            for step in range(0, LAST_STEP + 1, FRAME_INTERVAL):
                if step:
                    run.dual_time_iterate(time_step_count=FRAME_INTERVAL, max_iter_per_step=400)
                native_time = float(run.transient_controls.flow_time())
                assert abs(native_time - step * DT) < 1e-9, (step, native_time)
                data = out / f"tank-{step:04d}.dat.h5"
                st.file.write_data(file_name=str(data))
                with h5py.File(data) as handle:
                    group = alpha_group(handle)
                    if not step:
                        print("Water field:", group.name, flush=True)
                    alpha = cell_values(group)
                area = alpha.sum() * DX**2
                ledger.writerow([step, native_time, data.name, area, alpha.min(), alpha.max()])
                stream.flush()
                if not step:
                    with np.load(out / "geometry.npz") as geometry:
                        mapping = common_mapping(case_path, geometry)
                        centres = geometry["centres"]
                        expected = np.clip((0.05 + INITIAL_TILT * (centres[:, 0] - L / 2)
                                            - (centres[:, 1] - DX / 2)) / DX, 0, 1)
                    aligned = np.empty_like(alpha)
                    aligned[mapping] = alpha
                    max_difference = float(np.max(abs(aligned - expected)))
                    initial_pass = abs(area - EXPECTED_AREA) < 1e-10 and max_difference < 1e-10
                    write_json(out / "initial-field-check.json", dict(
                        water_area_m2=area, expected_water_area_m2=EXPECTED_AREA,
                        maximum_cell_fraction_difference=max_difference, passed=bool(initial_pass)))
                    assert initial_pass, (area, EXPECTED_AREA, max_difference)
                mark("saved_frame", step=step, time_s=native_time, water_area_m2=area)
                assert alpha.min() >= -1e-6 and alpha.max() <= 1. + 1e-6
                assert abs(area / EXPECTED_AREA - 1.) <= 0.005, area
        st.file.write_case_data(file_name=str(out / "tank-final.cas.h5"))
        mark("saved_final", final_time_s=native_time)
    except Exception:
        mark("failed", error=traceback.format_exc())
        raise
    finally:
        if session is not None:
            session.exit()
    mark("closed", final_time_s=native_time)


if __name__ == "__main__":
    main(Path(sys.argv[1]))

