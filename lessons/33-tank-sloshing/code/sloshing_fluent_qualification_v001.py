"""Qualify small-amplitude free-surface sloshing in hidden Fluent."""
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

NX, NY = 80, 40
DX, DT = 0.0025, 0.0005
L, H = NX * DX, NY * DX
WATER_WIDTH, WATER_HEIGHT = 0.20, 0.05
MATERIALS = {"air": (1.225, 1.7894e-5), "water-liquid": (998.2, 0.001003)}
WALLS = ["tank-walls"]
LAST_STEP = 40
FRAME_INTERVAL = 1
EXPECTED_AREA = WATER_WIDTH * WATER_HEIGHT
RESIDUAL_TARGET = 1e-4
INITIAL_TILT = float(os.environ.get("CFD_SLOSHING_TILT", "0.025"))
assert INITIAL_TILT in (0.0, 0.025)


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
    from wing_fluent_v001 import available_memory
    if available_memory() < 4 * 2**30:
        raise RuntimeError("Need 4 GiB available RAM before Fluent launch")
    if shutil.disk_usage(out.anchor).free < 10 * 2**30:
        raise RuntimeError("Need 10 GiB on E:")
    out.mkdir(parents=True)
    mesh = make_mesh(out)
    session = None
    status = out / "status.json"

    def mark(stage, **kwargs):
        write_json(status, dict(stage=stage, **kwargs))
        print(stage, kwargs, flush=True)

    try:
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
        for name in ("continuity", "x-velocity", "y-velocity"):
            residual = sol.monitor.residual.equations[name]
            residual.absolute_criteria = RESIDUAL_TARGET
            residual.check_convergence = True

        init = sol.initialization
        init.initialization_type = "standard"
        for name in ("x-velocity", "y-velocity", "pressure"):
            init.defaults[name] = 0.
        init.standard_initialize()
        def patch_rectangle(name, lower, upper, value):
            sol.cell_registers.create(name=name)
            region = sol.cell_registers[name].type
            region.option = next(option for option in region.option.allowed_values()
                                 if "hex" in option.lower())
            region.hexahedron.min_point = lower
            region.hexahedron.max_point = upper
            region.hexahedron.inside = True
            init.patch.calculate_patch(domain=phase_names[1], registers=[name],
                                       variable="mp", value=float(value))

        for column in range(NX):
            x_mid = (column + 0.5) * DX
            surface_height = 0.05 + INITIAL_TILT * (x_mid - L / 2)
            row = int(np.floor(surface_height / DX))
            fraction = surface_height / DX - row
            patch_rectangle(f"water_{column:03d}",
                            [column * DX, -DX],
                            [(column + 1) * DX, row * DX], 1.0)
            patch_rectangle(f"surface_{column:03d}",
                            [column * DX, row * DX],
                            [(column + 1) * DX, (row + 1) * DX], fraction)
        mark("water_patched", initial_tilt=INITIAL_TILT, mean_height_m=0.05)
        write_json(out / "setup.json", dict(
            dimensions_m=[L,H], cells=[NX,NY], dx_m=DX, time_step_s=DT,
            initial_mean_water_height_m=0.05, initial_tilt=INITIAL_TILT,
            expected_water_area_m2=EXPECTED_AREA, wall_patches=WALLS,
            model=mp.get_state(), operating_conditions=op.get_state(),
            materials={name:mats.fluid[name].get_state() for name in MATERIALS},
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
                    assert abs(area - EXPECTED_AREA) < 1e-10, (area, EXPECTED_AREA)
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

