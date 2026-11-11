"""Compare a saved OpenFOAM shock tube with the exact gas-dynamic solution."""
import argparse
import csv
import hashlib
import json
import re
import runpy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CELLS = 1000
LENGTH = 10.0
DX = LENGTH / CELLS
TIMES = tuple(i / 1000 for i in range(1, 8))
PRESSURE = 100_000.0
GAS_R = 8314.46261815324 / 28.96
CV = 1004.5 - GAS_R
DENSITY = PRESSURE / (GAS_R * 348.432)
LOW_DENSITY = 0.1 * PRESSURE / (GAS_R * 278.746)
VELOCITY = np.sqrt(PRESSURE / DENSITY)
EXPECTED_MASS = LENGTH * (DENSITY + LOW_DENSITY) / 2
EXPECTED_ENERGY = LENGTH * (1.1 * PRESSURE) * CV / (2 * GAS_R)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def field(path, kind):
    content = path.read_text(encoding="ascii")
    match = re.search(r"\binternalField\s+nonuniform\s+List<" + kind +
                      r">\s+(\d+)\s*\((.*?)\)\s*;", content, re.S)
    if match is None or int(match[1]) != CELLS:
        raise ValueError(f"Unexpected field count or format in {path}")
    values = np.fromstring(match[2].replace("(", " ").replace(")", " "), sep=" ")
    expected = (CELLS,) if kind == "scalar" else (CELLS, 3)
    values = values.reshape(expected)
    if not np.isfinite(values).all():
        raise ValueError(f"Nonfinite native field in {path}")
    return values


def audit(case, reference_path, output):
    input_manifest = case / "input-sha256.txt"
    for line in input_manifest.read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ", 1)
        assert sha(case / name) == expected, name
    log = (case / "solver.log").read_text(encoding="ascii")
    mesh_log = (case / "check-mesh.log").read_text(encoding="ascii")
    assert re.search(r"^\s*cells:\s*1000\s*$", mesh_log, re.M) and "Mesh OK." in mesh_log
    steps = np.array(re.findall(r"^Time = ([0-9.e+-]+)s$", log, re.M), dtype=float)
    assert np.allclose(steps, np.arange(1, 1401) * 5e-6, rtol=0, atol=5e-13)
    assert re.search(r"^End\s*$", log, re.M)
    reports = re.findall(r"^Courant Number mean:\s+\S+\s+max:\s+([0-9.e+-]+)[ \t]*$", log, re.M)
    paired = re.findall(r"^Courant Number mean:\s+\S+\s+max:\s+([0-9.e+-]+)[ \t]*\nTime = ([0-9.e+-]+)s$", log, re.M)
    assert len(paired) == 1400 and len(reports) == 1401
    assert [entry[0] for entry in paired] == reports[1:]
    assert np.array_equal(np.array([entry[1] for entry in paired], dtype=float), steps)
    courant = np.array(reports, dtype=float)
    assert np.isfinite(courant).all()

    centres = field(case / "0/C", "vector")
    index = np.argsort(centres[:, 0])
    x = centres[index, 0]
    assert np.allclose(x, -5 + (np.arange(CELLS) + 0.5)*DX, rtol=0, atol=2e-10)

    reference = runpy.run_path(str(reference_path))
    state = reference["star_state"].__globals__
    state["GAMMA"] = 1004.5 / CV
    state["RIGHT"] = (LOW_DENSITY / DENSITY, 0.0, 0.1)
    state["STAR"] = state["star_state"]()
    samples = []
    saved = []
    inputs = {name: sha(case / name) for name in ("input-sha256.txt", "solver.log", "check-mesh.log", "0/C")}
    for time in TIMES:
        name = f"{time:.3f}"
        directory = case / name
        clock = directory / "uniform/time"
        clock_text = clock.read_text(encoding="ascii")
        assert int(re.search(r"^index\s+(\d+);", clock_text, re.M)[1]) == round(time / 5e-6)
        assert abs(float(re.search(r"^value\s+([^;]+);", clock_text, re.M)[1])-time) < 1e-12
        rho = field(directory / "rho", "scalar")[index]
        p = field(directory / "p", "scalar")[index]
        temperature = field(directory / "T", "scalar")[index]
        velocity = field(directory / "U", "vector")[index]
        u = velocity[:, 0]
        assert min(rho.min(), p.min(), temperature.min()) > 0 and np.max(abs(velocity[:, 1:])) < 1e-10
        thermo_error = np.max(np.abs(p/(rho*GAS_R*temperature) - 1))
        assert thermo_error < 5e-6
        assert np.max(np.abs(u[[0,-1]])) < 1e-6
        assert np.allclose(p[[0,-1]], [PRESSURE, 0.1*PRESSURE], rtol=1e-7, atol=0)
        mass_error = abs(DX*np.sum(rho)/EXPECTED_MASS - 1)
        energy_error = abs(DX*np.sum(rho*(CV*temperature+u*u/2))/EXPECTED_ENERGY - 1)
        assert mass_error < 0.0005 and energy_error < 0.002
        physical_scale = np.array([DENSITY, VELOCITY, PRESSURE])
        exact = state["average_exact"](x/LENGTH+0.5, DX/LENGTH, time*VELOCITY/LENGTH) * physical_scale
        numerical = np.column_stack((rho, u, p))
        pressure_error = np.mean(np.abs(p-exact[:, 2])) / PRESSURE
        samples.append({"time_s": time, "pressure_l1_over_high_pressure": float(pressure_error),
                        "relative_mass_error": float(mass_error), "relative_total_energy_error": float(energy_error),
                        "maximum_relative_ideal_gas_error": float(thermo_error)})
        saved.append((numerical, exact))
        for file in (clock, *(directory/n for n in ("rho", "p", "T", "U"))):
            inputs[file.relative_to(case).as_posix()] = sha(file)

    numerical, exact = saved[-1]
    rho, _, p = numerical.T
    reference_wave_x = np.array(state["STAR"][-1])*VELOCITY*TIMES[-1]
    centres_between = (x[1:]+x[:-1])/2
    def wave_position(values, target, half_width):
        window = abs(centres_between-target) < half_width
        assert np.count_nonzero(window) > 2
        return float(centres_between[window][np.argmax(np.abs(np.diff(values))[window])])

    contact = wave_position(rho, reference_wave_x[2], .35)
    shock = wave_position(p, reference_wave_x[3], .30)
    assert samples[-1]["pressure_l1_over_high_pressure"] < 0.02
    assert abs(contact-reference_wave_x[2]) < .16 and abs(shock-reference_wave_x[3]) < .10
    checks = {"steps": len(steps), "cells": CELLS, "saved_states_checked": len(saved),
              "maximum_native_courant_number": float(courant.max()), "gamma": state["GAMMA"],
              "samples": samples, "contact_m": contact, "shock_m": shock,
              "exact_contact_m": float(reference_wave_x[2]), "exact_shock_m": float(reference_wave_x[3]),
              "native_sha256": inputs, "exact_reference_script_sha256": sha(reference_path)}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(checks, indent=2)+"\n", encoding="utf-8")
    with (output/"endpoint-profiles.csv").open("x", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(("x_m","rho_kg_m3","u_m_s","p_Pa","exact_rho_kg_m3","exact_u_m_s","exact_p_Pa"))
        writer.writerows((a,*b,*c) for a,b,c in zip(x,numerical,exact))
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7), layout="constrained")
    for j, (axis, label) in enumerate(zip(axes, ("Density (kg/m³)", "Axial velocity (m/s)", "Pressure (kPa)"))):
        scale = 1e-3 if j == 2 else 1
        axis.plot(x, exact[:,j]*scale, color="#485665", linestyle="--", label="Exact solution")
        axis.plot(x, numerical[:,j]*scale, color="#007fa3", label="OpenFOAM")
        axis.set(xlabel="Axial position (m)", ylabel=label, xlim=(-5,5))
        axis.grid(alpha=.25)
    axes[0].legend(loc="upper right")
    fig.suptitle("Shock tube at 7 ms | Native solution and exact profiles", fontweight="bold")
    for extension in ("png", "pdf"):
        fig.savefig(output/f"openfoam-shock-tube.{extension}", dpi=180)
    plt.close(fig)
    print(json.dumps({k:v for k,v in checks.items() if k!="native_sha256"}, indent=2))


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case", "reference", "output"):
        parser.add_argument(argument, type=Path)
    args=parser.parse_args()
    audit(args.case, args.reference, args.output)
