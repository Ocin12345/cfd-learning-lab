"""Run the 2-D inviscid compression corner in one hidden Fluent session."""
import json
from math import tan
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import numpy as np
import ansys.fluent.core as pf

from wedge_common_v001 import (CP, CORNER, H, L, LIMITS, MACH, MW, NX, NY,
                              P, P0, RESIDUALS, T, T0, THETA, U, reference)


def json_out(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=str) + "\n", encoding="utf-8")


def make_mesh(out):
    x = np.linspace(0, L, NX + 1)
    lower = np.maximum(x - CORNER, 0.) * tan(THETA)
    eta = np.linspace(0, 1, NY + 1)[:, None]
    yy = lower + eta * (H - lower)
    nodes = np.column_stack([np.broadcast_to(x, yy.shape).ravel(), yy.ravel()])
    ix = np.arange(NX)[None, :]
    iy = np.arange(NY)[:, None]
    k = (iy * (NX + 1) + ix).ravel()
    quads = np.column_stack([k, k + 1, k + NX + 2, k + NX + 1])
    polys = nodes[quads]
    nxt = np.roll(polys, -1, axis=1)
    cross = polys[..., 0] * nxt[..., 1] - nxt[..., 0] * polys[..., 1]
    areas = cross.sum(axis=1) / 2
    centres = ((polys + nxt) * cross[..., None]).sum(axis=1) / (6 * areas[:, None])
    assert np.all(areas > 0.)
    faces = {}
    for cid, q in enumerate(quads, 1):
        for a, b in zip(q, np.roll(q, -1)):
            key = tuple(sorted((int(a), int(b))))
            if key in faces:
                assert faces[key][3] == 0
                faces[key][3] = cid
            else:
                faces[key] = [int(a) + 1, int(b) + 1, cid, 0]
    zones = {zone: [] for zone in (3, 4, 5, 6, 7)}
    for a, b, c0, c1 in faces.values():
        n0, n1 = nodes[[a - 1, b - 1]]
        if c1:
            zone = 3
        elif n0[0] == n1[0] == 0:
            zone = 6
        elif n0[0] == n1[0] == L:
            zone = 7
        elif np.allclose([n0[1], n1[1]], H, atol=1e-13, rtol=0):
            zone = 4
        else:
            zone = 5
        zones[zone].append((a, b, c0, c1))
    assert [len(zones[z]) for z in (4, 5, 6, 7)] == [NX, NX, NY, NY]
    mesh = out / "wedge.msh"
    with mesh.open("x", encoding="ascii", newline="\n") as f:
        f.write('(0 "Mach 2 compression corner")\n(2 2)\n')
        f.write(f"(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n")
        for xx, y in nodes:
            f.write(f"{xx:.17g} {y:.17g}\n")
        f.write(f"))\n(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n")
        f.write(f"(13 (0 1 {len(faces):x} 0))\n")
        start = 1
        for zone, bc in [(3, 2), (4, 3), (5, 3), (6, 4), (7, 5)]:
            rows = zones[zone]
            f.write(f"(13 ({zone:x} {start:x} {start+len(rows)-1:x} {bc:x} 2)(\n")
            for row in rows:
                f.write(" ".join(f"{v:x}" for v in row) + "\n")
            f.write("))\n")
            start += len(rows)
        for zone, kind, name in ((2, "fluid", "fluid"), (3, "interior", "interior"),
                                  (4, "wall", "top"), (5, "wall", "ramp"),
                                  (6, "pressure-inlet", "inlet"), (7, "pressure-outlet", "outlet")):
            f.write(f"(45 ({zone} {kind} {name})())\n")
    np.savez_compressed(out / "geometry.npz", nodes=nodes, quads=quads, areas=areas,
                        centres=centres, x_nodes=x, y_nodes=yy)
    return mesh


def main(out, support):
    sys.path.insert(0, str(support))
    from nozzle_fluent_v004 import hdf_terminal, quality_metrics

    if shutil.disk_usage(out.anchor).free < 10 * 2**30:
        raise RuntimeError("Require 10 GiB free on the run drive")
    out.mkdir(parents=True)
    mesh = make_mesh(out)
    json_out(out / "reference.json", reference())
    session = None
    checkpoints = []

    def stage(label, **kwargs):
        json_out(out / "status.json", dict(stage=label, **kwargs))
        print(label, kwargs, flush=True)

    try:
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
            if actual >= 520 and all(residual[name] <= target for name, target in RESIDUALS.items()):
                break
        json_out(out / "loaded-settings.json", dict(
            general=st.setup.general.get_state(), viscous=st.setup.models.viscous.get_state(),
            gas=gas.get_state(), inlet=inlet.get_state(), outlet=outlet.get_state(),
            methods=method.get_state(), reference_mach=MACH, reference=reference()))
        stage("saved_final", final=checkpoints[-1])
    except Exception:
        stage("failed", error=traceback.format_exc())
        raise
    finally:
        if session is not None:
            session.exit()
    stage("closed", final=checkpoints[-1])


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))

