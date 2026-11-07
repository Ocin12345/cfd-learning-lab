"""Audit the saved Fluent sloshing qualification and its flat-water control."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import sys

import h5py
import numpy as np
from sloshing_fluent_qualification_v001 import (
    DT, LAST_STEP, MATERIALS, NX, NY, DX, RESIDUAL_TARGET,
    alpha_group, cell_values,
)
from analyse_nozzle_fluent_v005 import common_mapping

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

def verify_steps(times):
    assert len(times) == LAST_STEP
    assert np.allclose(times, DT * np.arange(1, LAST_STEP + 1), atol=1e-10, rtol=0)

def fluent(run):
    log = (run / "solver-transcript.txt").read_text(encoding="utf-8", errors="replace")
    state = json.loads((run / "setup.json").read_text(encoding="utf-8"))
    assert state["vof_patch_options"]["patch_reconstructed_interface"] is False
    assert state["vof_patch_options"]["use_volumetric_smoothing"] is False
    reporting = state["residual_reporting"]
    assert state["initial_tilt"] == 0.025
    assert reporting["options"]["normalize"] is False
    assert reporting["options"]["residual_values"]["scale_residuals"] is True
    assert reporting["scaled_target"] == RESIDUAL_TARGET
    for name in ("continuity", "x-velocity", "y-velocity"):
        assert reporting["equations"][name]["check_convergence"] is True
        assert reporting["equations"][name]["absolute_criteria"] == RESIDUAL_TARGET
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--flat", action="store_true")
    args = parser.parse_args()
    assert not args.out.exists()
    assert not args.flat, "This auditor is restricted to tilted flow"
    result = fluent(args.run)
    result["auditor_source_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    residual_pass = result["steps_meeting_all_scaled_residual_targets"] == LAST_STEP
    courant_pass = result["peak_native_reported_explicit_vof_courant_number"] <= 0.25
    setup = json.loads((args.run / "setup.json").read_text())
    tilt = 0.0 if args.flat else 0.025
    assert setup["initial_tilt"] == tilt
    assert setup["cells"] == [NX, NY] and abs(setup["time_step_s"] - DT) < 1e-14
    with np.load(args.run / "geometry.npz") as geometry:
        mapping = common_mapping(args.run / "tank.cas.h5", geometry)
        centres = geometry["centres"]
        expected_alpha = np.clip((0.05 + tilt * (centres[:, 0] - 0.1)
                                  - (centres[:, 1] - DX / 2)) / DX, 0, 1)
    with (args.run / "frames.csv").open(newline="") as stream:
        frames = list(csv.DictReader(stream))
    assert [int(row["step"]) for row in frames] == list(range(LAST_STEP + 1))
    records = []
    sources = [args.run / name for name in ["tank.cas.h5", "setup.json", "solver-transcript.txt", "frames.csv"]]
    for row in frames:
        path = args.run / row["file"]
        with h5py.File(path) as native:
            from ansys.fluent.core.offline import _lispy
            clock_variables = {v[0]:v[1] for v in _lispy.parse(native["settings/Data Variables"][0].decode())[1]}
            native_time = {key:float(clock_variables[key]) for key in ("flow-time", "time-step", "physical-time-step")}
            assert np.isfinite(list(native_time.values())).all()
            assert native_time["time-step"] == int(row["step"])
            assert abs(native_time["flow-time"] - int(row["step"]) * DT) < 1e-10
            assert abs(native_time["physical-time-step"] - DT) < 1e-14
            assert abs(float(row["time_s"]) - native_time["flow-time"]) < 1e-10
            alpha = cell_values(alpha_group(native))
            velocity = native["results/1/phase-1/cells"]
            speed = np.hypot(cell_values(velocity["SV_U"]), cell_values(velocity["SV_V"]))
        aligned = np.empty_like(alpha)
        aligned[mapping] = alpha
        assert np.isfinite(speed).all()
        assert alpha.min() >= -1e-6 and alpha.max() <= 1 + 1e-6
        water_area = float(alpha.sum() * DX**2)
        assert abs(water_area / 0.01 - 1) <= 0.005
        if int(row["step"]) == 0:
            assert np.max(abs(aligned - expected_alpha)) < 1e-10
        record = dict(step=int(row["step"]), time_s=float(row["time_s"]),
                      water_area_m2=water_area, max_speed_m_s=float(speed.max()),
                      alpha_min=float(alpha.min()), alpha_max=float(alpha.max()))
        if args.flat:
            assert record["max_speed_m_s"] <= 1e-4
            assert np.max(abs(aligned - expected_alpha)) <= 0.01
        records.append(record)
        sources.append(path)
    result.update(scope="Flat-water control" if args.flat else "Tilted-surface qualification",
                  native_field_checks_pass=True, qualification_pass=residual_pass and courant_pass,
                  scaled_residual_checks_pass=residual_pass, courant_check_pass=courant_pass,
                  production_complete=False, records=records,
                  native_source_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                                        for path in sources})
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ["records", "native_source_sha256"]}, indent=2))
    print(json.dumps(dict(max_saved_speed_m_s=max(row["max_speed_m_s"] for row in records),
                          max_saved_water_area_drift=max(abs(row["water_area_m2"]/.01-1) for row in records))))
    if not result["qualification_pass"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
