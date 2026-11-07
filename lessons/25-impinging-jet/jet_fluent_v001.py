"""Make the planar jet mesh and run steady Fluent in one hidden session."""
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import h5py
import numpy as np
import ansys.fluent.core as pf

from jet_common_v001 import (DX, H, INLET_MASS, MASS_BALANCE, MU, NX, NY, RE, RESIDUAL,
                             RHO, UMEAN, V_SLOT, W, X, X_SLOT, Y, ZONES)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str) + "\n", encoding="utf-8")


def make_mesh(out):
    nodes = np.array([(x, y) for y in Y for x in X])
    k = (np.arange(NY)[:, None] * (NX + 1) + np.arange(NX)).ravel()
    quads = np.column_stack((k, k + 1, k + NX + 2, k + NX + 1))
    faces = {}
    for cid, q in enumerate(quads, 1):
        for a, b in zip(q, np.roll(q, -1)):
            key = tuple(sorted((int(a), int(b))))
            if key in faces:
                assert faces[key][3] == 0
                faces[key][3] = cid
            else:
                faces[key] = [int(a) + 1, int(b) + 1, cid, 0]
    groups = {zone: [] for zone in range(3, 9)}
    for a, b, c0, c1 in faces.values():
        A, B = nodes[[a - 1, b - 1]]
        if c1:
            zone = 3
        elif A[1] == B[1] == 0:
            zone = 4
        elif A[0] == B[0] == X[0]:
            zone = 6
        elif A[0] == B[0] == X[-1]:
            zone = 7
        else:
            assert A[1] == B[1] == H
            zone = 8 if abs((A[0] + B[0]) * .5) < W / 2 else 5
        groups[zone].append((a, b, c0, c1))
    assert [len(groups[z]) for z in range(4, 9)] == [NX, NX - len(X_SLOT), NY, NY, len(X_SLOT)]
    msh = out / "jet.msh"
    with msh.open("x", encoding="ascii", newline="\n") as f:
        f.write('(0 "Planar impinging slot jet; uniform quads")\n(2 2)\n')
        f.write(f"(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n")
        for x, y in nodes:
            f.write(f"{x:.17g} {y:.17g}\n")
        f.write(f"))\n(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n")
        f.write(f"(13 (0 1 {len(faces):x} 0))\n")
        first = 1
        for zone, bc in ((3, 2), (4, 3), (5, 3), (6, 5), (7, 5), (8, 10)):
            rows = groups[zone]
            f.write(f"(13 ({zone:x} {first:x} {first+len(rows)-1:x} {bc:x} 2)(\n")
            for row in rows:
                f.write(" ".join(f"{v:x}" for v in row) + "\n")
            f.write("))\n")
            first += len(rows)
        for zone, kind, name in ((2, "fluid", "fluid"), (3, "interior", "interior"),
                                  (4, "wall", "plate"), (5, "wall", "ceiling"),
                                  (6, "pressure-outlet", "left"), (7, "pressure-outlet", "right"),
                                  (8, "velocity-inlet", "inlet")):
            f.write(f"(45 ({zone} {kind} {name})())\n")
    geo = dict(nodes=nodes, quads=quads, centres=nodes[quads].mean(axis=1),
               areas=np.full(NX * NY, DX**2), x_nodes=X, y_nodes=Y)
    np.savez_compressed(out / "geometry.npz", **geo)
    profile = out / "slot-inlet.prof"
    with profile.open("x", encoding="ascii", newline="\n") as f:
        f.write(f"((slot-inlet point {len(X_SLOT)})\n")
        for name, values in (("x", X_SLOT), ("y", np.full(len(X_SLOT), H)), ("v", V_SLOT)):
            f.write("(" + name + "\n" + "\n".join(f"{v:.17g}" for v in values) + "\n)\n")
        f.write(")\n")
    return msh, profile


def read_residual_endpoint(data):
    result = {}
    actual = None
    with h5py.File(data) as h:
        for name in ("continuity", "x-velocity", "y-velocity"):
            eq = h[f"results/residuals/phase-1/{name}"]
            iteration = int(np.asarray(eq["iterations"][-1]).flat[0])
            unscaled, scale = eq["data"][-1, :2]
            assert actual in (None, iteration) and scale > 0
            result[name] = float(unscaled / scale)
            assert math.isfinite(result[name]) and result[name] >= 0.
            actual = iteration
    return actual, result


def main(out, support):
    sys.path.insert(0, str(support))
    from nozzle_fluent_v004 import quality_metrics
    if shutil.disk_usage(out.anchor).free < 10 * 2**30:
        raise RuntimeError("Need 10 GiB free on the run drive")
    out.mkdir(parents=True)
    mesh, profile = make_mesh(out)
    session = None
    endpoints = []

    def status(stage, **extra):
        write_json(out / "status.json", dict(stage=stage, **extra))
        print(stage, extra, flush=True)

    try:
        status("launching")
        session = pf.launch_fluent(fluent_path=os.environ["FLUENT_EXECUTABLE"],
                                   dimension=2, precision="double", processor_count=2,
                                   ui_mode="hidden_gui", cwd=str(out), start_timeout=180,
                                   cleanup_on_exit=True)
        status("licensed_fluent_started", version=str(session.get_fluent_version()))
        transcript = out / "solver-transcript.txt"
        session.transcript.start(file_name=str(transcript))
        st = session.settings
        st.file.read_mesh(file_name=str(mesh))
        st.mesh.check()
        st.mesh.quality()
        for trial in range(15):
            try:
                quality = quality_metrics(transcript.read_text(encoding="utf-8", errors="replace"))
                break
            except ValueError:
                if trial == 14:
                    raise
                time.sleep(.2)
        write_json(out / "mesh-quality.json", quality)
        st.setup.general.solver.time = "steady"
        st.setup.models.viscous.model = "laminar"
        st.setup.general.operating_conditions.operating_pressure = 0.
        st.setup.general.operating_conditions.gravity.enable = False
        st.setup.materials.fluid["air"].density.option = "constant"
        st.setup.materials.fluid["air"].density.value = RHO
        st.setup.materials.fluid["air"].viscosity.option = "constant"
        st.setup.materials.fluid["air"].viscosity.value = MU
        st.setup.cell_zone_conditions.fluid["fluid"].general.material = "air"
        st.file.read_profile(file_name=str(profile))
        inlet = st.setup.boundary_conditions.velocity_inlet["inlet"].momentum
        inlet.velocity_specification_method = "Components"
        inlet.velocity_components = [dict(option="value", value=0.),
                                     dict(option="profile", profile_name="slot-inlet", field_name="v")]
        components = inlet.velocity_components.get_state()
        assert components[1]["profile_name"] == "slot-inlet" and components[1]["field_name"] == "v"
        for name in ("left", "right"):
            st.setup.boundary_conditions.pressure_outlet[name].momentum.gauge_pressure.value = 0.
        for name in ("plate", "ceiling"):
            wall = st.setup.boundary_conditions.wall[name].momentum
            wall.wall_motion = "Stationary Wall"
            wall.shear_condition = "No Slip"
        st.setup.reference_values.depth = 1.
        method = st.solution.methods
        method.p_v_coupling.flow_scheme = "SIMPLEC"
        method.spatial_discretization.gradient_scheme = "least-square-cell-based"
        method.spatial_discretization.discretization_scheme = {
            "mom": "second-order-upwind", "pressure": "second-order"}
        residual = st.solution.monitor.residual
        residual.options.n_save = 2800
        for name in ("continuity", "x-velocity", "y-velocity"):
            residual.equations[name].absolute_criteria = RESIDUAL
            residual.equations[name].check_convergence = False
        init = st.solution.initialization
        init.initialization_type = "standard"
        for name in ("pressure", "x-velocity", "y-velocity"):
            init.defaults[name] = 0.
        init.standard_initialize()
        write_json(out / "setup.json", dict(
            cells=NX * NY, nominal_slot_width_m=W, gap_m=H, Reynolds=RE,
            inlet_mean_m_s=UMEAN, discrete_expected_mass_kg_s_per_m_depth=INLET_MASS,
            material=st.setup.materials.fluid["air"].get_state(), inlet=components,
            method=method.get_state(), residual_target=RESIDUAL, mass_imbalance_target=MASS_BALANCE,
            outlet_pressure_pa=0.))
        actual = 0
        for endpoint in (250, 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500):
            status("iterating", from_iteration=actual, requested_endpoint=endpoint)
            st.solution.run_calculation.iterate(iter_count=endpoint - actual)
            stem = f"jet-{endpoint:04d}"
            st.file.write_case_data(file_name=str(out / (stem + ".cas.h5")))
            actual, residuals = read_residual_endpoint(out / (stem + ".dat.h5"))
            assert actual == endpoint
            st.results.report.fluxes.mass_flow(zones=list(ZONES), write_to_file=True,
                                              file_name=str(out / (stem + "-mass.txt")))
            endpoints.append(dict(iteration=endpoint, stem=stem, residuals=residuals))
            write_json(out / "checkpoints.json", endpoints)
            status("checkpoint_saved", **endpoints[-1])
            if actual >= 750 and max(residuals.values()) <= RESIDUAL:
                break
        status("saved_final", **endpoints[-1])
    except Exception:
        status("failed", error=traceback.format_exc())
        raise
    finally:
        if session:
            session.exit()
    status("closed", **endpoints[-1])


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))

