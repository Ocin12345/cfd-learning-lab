"""Compute 1.2 s of small-amplitude free-surface sloshing in hidden Fluent."""
import csv
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import h5py
import numpy as np
import ansys.fluent.core as pf
from analyse_nozzle_fluent_v005 import common_mapping

NX, NY = 80, 40
DX, DT = 0.0025, 0.0005
L, H = NX * DX, NY * DX
WATER_WIDTH, WATER_HEIGHT = 0.20, 0.05
MATERIALS = {"air": (1.225, 1.7894e-5), "water-liquid": (998.2, 0.001003)}
WALLS = ["tank-walls"]
LAST_STEP = 2400
FRAME_INTERVAL = 20
EXPECTED_AREA = WATER_WIDTH * WATER_HEIGHT
RESIDUAL_TARGET = 1e-4
INITIAL_TILT = float(os.environ.get("CFD_SLOSHING_TILT", "0.025"))
assert INITIAL_TILT == 0.025, "This driver is restricted to the declared tilted surface"


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def make_mesh(out):
    x = DX * np.arange(NX + 1)
    y = DX * np.arange(NY + 1)
    nodes = np.array([(xx, yy) for yy in y for xx in x])
    quads = []
    faces = {}
    for j in range(NY):
        for i in range(NX):
            k = j * (NX + 1) + i
            q = (k, k + 1, k + NX + 2, k + NX + 1)
            quads.append(q)
            for a, b in zip(q, (*q[1:], q[0])):
                key = tuple(sorted((a, b)))
                if key in faces:
                    faces[key][3] = len(quads)
                else:
                    faces[key] = [a + 1, b + 1, len(quads), 0]
    inside = [row for row in faces.values() if row[3]]
    wall = [row for row in faces.values() if not row[3]]
    assert len(wall) == 2 * (NX + NY)
    assert len(inside) == (NX - 1) * NY + (NY - 1) * NX
    centres = nodes[np.asarray(quads)].mean(axis=1)
    np.savez_compressed(out / "geometry.npz", nodes=nodes, quads=quads,
                        centres=centres, x_edges=x, y_edges=y)
    mesh = out / "tank.msh"
    with mesh.open("x", encoding="ascii", newline="\n") as f:
        f.write('(0 "Sloshing tank, uniform 2.5 mm quads")\n(2 2)\n')
        f.write(f"(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n")
        for xx, yy in nodes:
            f.write(f"{xx:.17g} {yy:.17g}\n")
        f.write(f"))\n(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n")
        f.write(f"(13 (0 1 {len(faces):x} 0))\n")
        first = 1
        for zone, bc, rows in ((3, 2, inside), (4, 3, wall)):
            f.write(f"(13 ({zone:x} {first:x} {first + len(rows) - 1:x} {bc:x} 2)(\n")
            for row in rows:
                f.write(" ".join(f"{value:x}" for value in row) + "\n")
            f.write("))\n")
            first += len(rows)
        f.write("(45 (2 fluid fluid)())\n(45 (3 interior interior)())\n(45 (4 wall tank-walls)())\n")
    return mesh


def alpha_group(handle):
    phases = handle["results/1"]
    groups = [phase["cells/SV_VOF"] for phase in phases.values() if "cells/SV_VOF" in phase]
    if len(groups) != 1:
        matches = [group for group in groups if "/phase-3/" in group.name]
        if len(matches) != 1:
            raise ValueError(f"Expected Fluent secondary-phase VOF: {[g.name for g in groups]}")
        return matches[0]
    return groups[0]


def cell_values(group):
    result = np.full(NX * NY, np.nan)
    for ds in group.values():
        first = int(np.asarray(ds.attrs["minId"]).flat[0]) - 1
        last = int(np.asarray(ds.attrs["maxId"]).flat[0])
        assert 0 <= first < last <= len(result) and np.isnan(result[first:last]).all()
        result[first:last] = np.asarray(ds).ravel()
    assert np.isfinite(result).all()
    return result


def main(out):
    from qualification_admission_v001 import verify_qualifications
    qualification_hashes = verify_qualifications()
    from wing_fluent_v001 import available_memory
    while available_memory() < 4 * 2**30:
        print("waiting_for_launch_memory", {"available_ram_GiB": available_memory() / 2**30,
                                           "required_ram_GiB": 4}, flush=True)
        time.sleep(30)
    if shutil.disk_usage(out.anchor).free < 10 * 2**30:
        raise RuntimeError("Need 10 GiB on E:")
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
        qualification_hashes = verify_qualifications()
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
        monitor.options.residual_values.scale_residuals = True
        for name in ("continuity", "x-velocity", "y-velocity"):
            residual = monitor.equations[name]
            residual.absolute_criteria = RESIDUAL_TARGET
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
            total_steps=LAST_STEP, frame_interval=FRAME_INTERVAL,
            qualification_audit_sha256=qualification_hashes,
            initial_mean_water_height_m=0.05, initial_tilt=INITIAL_TILT,
            expected_water_area_m2=EXPECTED_AREA, wall_patches=WALLS,
            model=mp.get_state(), operating_conditions=op.get_state(),
            vof_patch_options=patch_options,
            materials={name:mats.fluid[name].get_state() for name in MATERIALS},
            residual_reporting=dict(options=monitor.options.get_state(),
                                    equations=monitor.equations.get_state(),
                                    scaled_target=RESIDUAL_TARGET),
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

