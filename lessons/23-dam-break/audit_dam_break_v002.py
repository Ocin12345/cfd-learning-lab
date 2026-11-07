"""Audit native logs and selected setup values for the paired dam-break run."""
import json
from pathlib import Path
import re
import sys

import numpy as np

from dam_break_fluent_v004 import DT, LAST_STEP, MATERIALS, NX, NY, RESIDUAL_TARGET

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def verify_steps(times):
    assert len(times) == LAST_STEP
    assert np.allclose(times, DT * np.arange(1, LAST_STEP + 1), atol=1e-10, rtol=0)


def fluent(run):
    log = (run / "solver-transcript.txt").read_text(encoding="utf-8", errors="replace")
    state = json.loads((run / "setup.json").read_text(encoding="utf-8"))
    assert state["model"]["model"] == "vof"
    assert state["model"]["vof_parameters"]["vof_formulation"] == "explicit"
    assert state["model"]["phases"]["phase-2"]["material"] == "water-liquid"
    assert state["model"]["phase_interaction"]["forces"]["surface_tension_model"] is False
    for name, (rho, mu) in MATERIALS.items():
        assert state["materials"][name]["density"]["value"] == rho
        assert state["materials"][name]["viscosity"]["value"] == mu
    assert state["operating_conditions"]["gravity"]["components"] == [0, -9.81]
    assert json.loads((run / "status.json").read_text())["stage"] == "closed"
    assert (run / "tank-final.cas.h5").stat().st_size > 0
    assert (run / "tank-final.dat.h5").stat().st_size > 0

    quality = float(re.search(r"Minimum Orthogonal Quality\s*=\s*(" + NUMBER + ")", log)[1])
    min_volume = float(re.search(r"minimum volume \(m3\):\s*(" + NUMBER + ")", log)[1])
    assert quality > .2 and min_volume > 0.
    step_pattern = re.compile(r"Flow time = (" + NUMBER + r")s, time step = (\d+)")
    residual_pattern = re.compile(r"^\s*(\d+)\s+(" + NUMBER + r")\s+(" + NUMBER +
                                  r")\s+(" + NUMBER + r")\s+\d+:\d+:\d+\s+\d+\s*$", re.M)
    boundaries = list(step_pattern.finditer(log))
    assert [int(x[2]) for x in boundaries] == list(range(1, LAST_STEP + 1))
    verify_steps([float(x[1]) for x in boundaries])
    terminal = []
    start = 0
    for stop in boundaries:
        rows = list(residual_pattern.finditer(log, start, stop.start()))
        assert rows
        values = [float(rows[-1][i]) for i in range(2, 5)]
        assert np.isfinite(values).all()
        terminal.append(values)
        start = stop.end()
    terminal = np.array(terminal)
    converged = np.all(terminal <= RESIDUAL_TARGET, axis=1)
    co = [float(x) for x in re.findall(r"Global Courant Number \[Explicit VOF Criteria\]\s*:\s*(" +
                                       NUMBER + r")", log)]
    assert co
    return {
        "cells": NX * NY, "steps": LAST_STEP, "last_time_s": LAST_STEP * DT,
        "minimum_orthogonal_quality": quality, "minimum_cell_area_m2": min_volume,
        "last_step_scaled_residuals": dict(zip(("continuity", "x-velocity", "y-velocity"), terminal[-1])),
        "max_terminal_scaled_residuals": dict(zip(("continuity", "x-velocity", "y-velocity"), terminal.max(0))),
        "scaled_residual_target": RESIDUAL_TARGET,
        "steps_meeting_all_scaled_residual_targets": int(converged.sum()),
        "peak_native_reported_explicit_vof_courant_number": max(co),
        "time_integration": state["numerical_methods"]["transient_formulation"],
        "interface_scheme": state["numerical_methods"]["spatial_discretization"]["discretization_scheme"]["mp"],
    }


def openfoam(run):
    log = (run / "solver.log").read_text(encoding="utf-8", errors="replace")
    mesh = (run / "mesh-check.log").read_text(encoding="utf-8", errors="replace")
    assert (run / "run-status.txt").read_text().strip() == "complete"
    assert "Mesh OK." in mesh and re.search(r"\bcells:\s+7200\b", mesh)
    assert re.search(r"^End\s*$", log, re.M) and "FOAM FATAL" not in log
    times = [float(t) for t in re.findall(r"^Time = (" + NUMBER + r")s?\s*$", log, re.M)]
    verify_steps(times)
    co = [float(x) for x in re.findall(r"^Courant Number mean: " + NUMBER + r" max: (" + NUMBER + ")", log, re.M)]
    interface_co = [float(x) for x in re.findall(r"^Interface Courant Number mean: " + NUMBER +
                                                r" max: (" + NUMBER + ")", log, re.M)]
    continuity = re.findall(r"time step continuity errors : sum local = (" + NUMBER + r"), global = (" +
                            NUMBER + r"), cumulative = (" + NUMBER + ")", log)
    assert co and interface_co and continuity
    return {
        "cells": NX * NY, "steps": LAST_STEP, "last_time_s": LAST_STEP * DT,
        "native_mesh_check": "Mesh OK.",
        "peak_native_reported_courant_number": max(co),
        "peak_native_reported_interface_courant_number": max(interface_co),
        "final_continuity_error_native_sum_local": float(continuity[-1][0]),
        "final_continuity_error_native_global": float(continuity[-1][1]),
        "time_integration": "Euler",
        "interface_scheme": "interfaceCompression vanLeer, MULES",
    }


if __name__ == "__main__":
    f, o, out = (Path(x) for x in sys.argv[1:])
    results = {"ANSYS Fluent": fluent(f), "OpenFOAM": openfoam(o)}
    with out.open("x", encoding="utf-8") as stream:
        json.dump(results, stream, allow_nan=False, indent=2)
        stream.write("\n")
    print(json.dumps(results, indent=2))

