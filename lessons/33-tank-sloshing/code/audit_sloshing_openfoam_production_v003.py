"""Audit completed sloshing qualification fields and native OpenFOAM logs."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import numpy as np

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
CELLS = 3200
DX = 0.0025
DT = 0.0005
STEPS = 2400
FRAME_INTERVAL = 20

def field(path, components):
    text = path.read_text(encoding="ascii")
    uniform = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if uniform:
        value = np.fromstring(uniform[1].replace("(", " ").replace(")", " "), sep=" ")
        assert len(value) == components
        return np.tile(value, (CELLS, 1)).reshape(CELLS, components)
    match = re.search(r"internalField\s+nonuniform\s+List<\w+>\s+(\d+)\s*\(", text)
    assert match and int(match[1]) == CELLS, path
    closing = re.search(r"\n\s*\)\s*;", text[match.end():])
    assert closing, path
    body = text[match.end():match.end() + closing.start()]
    values = np.fromstring(body.replace("(", " ").replace(")", " "), sep=" ")
    assert values.size == CELLS * components
    return values.reshape(CELLS, components)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    parser.add_argument("out", type=Path)
    args = parser.parse_args()
    assert not args.out.exists()
    case = args.case
    declaration = json.loads((case / "declared-settings.json").read_text())
    tilt = declaration["initial_tilt"]
    assert tilt == .025 and declaration["momentum_predictor"]
    assert declaration["steps"] == STEPS and declaration["frame_interval"] == FRAME_INTERVAL
    assert declaration["dt_s"] == DT and declaration["duration_s"] == STEPS * DT
    assert (case / "run-status.txt").read_text().strip() == "complete"
    mesh_log = (case / "mesh-check.log").read_text()
    assert "Mesh OK." in mesh_log and re.search(r"\bcells:\s+3200\b", mesh_log)
    raw = (case / "solver.log").read_text()
    assert re.search(r"^End\s*$", raw, re.M) and "FOAM FATAL" not in raw
    chunks = re.split(r"^Time = ([0-9.eE+-]+)s?\s*$", raw, flags=re.M)
    targets = {"Ux": 1e-4, "Uy": 1e-4, "p_rgh": 1e-4}
    steps = []
    for index in range(1, len(chunks), 2):
        time = float(chunks[index])
        chunk = chunks[index + 1]
        assert "PIMPLE: Iteration 8" in chunk
        final_outer = chunk.split("PIMPLE: Iteration 8")[-1]
        residuals = {}
        for match in re.finditer(r"Solving for (\w+), Initial residual = ([-+\d.eE]+), Final residual = ([-+\d.eE]+)", final_outer):
            name = match[1]
            if name in targets:
                residuals[name] = max(residuals.get(name, 0), float(match[2]))
        assert set(residuals) == set(targets)
        alpha_rows = re.findall(r"Min\(alpha.water\)\s*=\s*(" + NUMBER + r").*?Max\(alpha.water\)\s*=\s*(" + NUMBER + r")", final_outer)
        assert alpha_rows
        amin, amax = map(float, alpha_rows[-1])
        assert np.isfinite([amin, amax]).all()
        steps.append(dict(time_s=time, residuals=residuals,
                          passed=all(value <= targets[name] for name, value in residuals.items()),
                          alpha_min_final_outer=amin,alpha_max_final_outer=amax,
                          alpha_bounds_pass=bool(amin>=-1e-6 and amax<=1+1e-6)))
    assert len(steps) == STEPS
    assert np.allclose([row["time_s"] for row in steps], DT * np.arange(1, STEPS + 1), atol=1e-10, rtol=0)
    interface_co = [float(value) for value in re.findall(r"^Interface Courant Number mean: " + NUMBER + r" max: (" + NUMBER + ")", raw, re.M)]
    assert interface_co and max(interface_co) <= .25
    centres = field(case / "0/C", 3)
    indices = np.rint(centres[:, :2] / DX - .5).astype(int)
    assert np.all((indices >= 0) & (indices < [80, 40]))
    assert np.max(abs(centres[:, :2] - (indices + .5) * DX)) < 1e-10
    assert len(np.unique(indices[:, 0] + 80 * indices[:, 1])) == CELLS
    expected = np.clip((.05 + tilt * (centres[:, 0] - .1) - indices[:, 1] * DX) / DX, 0, 1)
    records = []
    sources = [case / name for name in ["declared-settings.json", "solver.log", "mesh-check.log", "0/C"]]
    for step in range(0, STEPS + 1, FRAME_INTERVAL):
        folder = case / f"{step * DT:.8g}"
        alpha_path, velocity_path = folder / "alpha.water", folder / "U"
        alpha = field(alpha_path, 1).ravel()
        velocity = field(velocity_path, 3)
        assert np.isfinite(alpha).all() and np.isfinite(velocity).all()
        alpha_bounds_pass = bool(alpha.min() >= -1e-6 and alpha.max() <= 1 + 1e-6)
        area = float(alpha.sum() * DX**2)
        assert abs(area / .01 - 1) <= .005
        if step == 0:
            assert np.max(abs(alpha - expected)) < 1e-10
        else:
            native_time_path = folder / "uniform/time"
            native_time = re.search(r"^value\s+(" + NUMBER + r")\s*;", native_time_path.read_text(), re.M)
            assert native_time and abs(float(native_time[1]) - step * DT) < 1e-10
            native_text = native_time_path.read_text()
            native_index = re.search(r"^index\s+(\d+)\s*;", native_text, re.M)
            native_dt = re.search(r"^deltaT\s+(" + NUMBER + r")\s*;", native_text, re.M)
            assert native_index and int(native_index[1]) == step
            assert native_dt and abs(float(native_dt[1]) - DT) < 1e-14
            sources.append(native_time_path)
        speed = float(np.linalg.norm(velocity, axis=1).max())
        if tilt == 0:
            assert speed <= 1e-4 and np.max(abs(alpha - expected)) <= .01
        records.append(dict(step=step, time_s=step*DT, water_area_m2=area, max_speed_m_s=speed,
                            alpha_min=float(alpha.min()),alpha_max=float(alpha.max()),alpha_bounds_pass=alpha_bounds_pass))
        sources.extend([alpha_path, velocity_path])
    residuals_pass = all(row["passed"] for row in steps)
    bounds_pass = all(row["alpha_bounds_pass"] for row in records)
    log_bounds_pass = all(row["alpha_bounds_pass"] for row in steps)
    sources.extend([case / name for name in ["system/fvSolution","system/fvSchemes","system/controlDict","constant/g","constant/phaseProperties","constant/physicalProperties.water","constant/physicalProperties.air","constant/momentumTransport"]])
    result = dict(production_checks_pass=residuals_pass and bounds_pass and log_bounds_pass,
                  all_step_alpha_bounds_pass=log_bounds_pass,
                  steps_meeting_alpha_bounds=sum(row["alpha_bounds_pass"] for row in steps),
                  maximum_final_outer_alpha=max(row["alpha_max_final_outer"] for row in steps),
                  residual_checks_pass=residuals_pass, alpha_bounds_pass=bounds_pass,
                  alpha_bound_tolerance=1e-6, minimum_saved_alpha=min(row["alpha_min"] for row in records),
                  maximum_saved_alpha=max(row["alpha_max"] for row in records),
                  saved_states_meeting_alpha_bounds=sum(row["alpha_bounds_pass"] for row in records),
                  maximum_relative_water_area_drift=max(abs(row["water_area_m2"]/.01-1) for row in records), scope="flat control" if tilt == 0 else "tilted production",
                  declared_residual_targets=targets, steps=steps, records=records,
                  peak_interface_courant=max(interface_co), native_field_checks_pass=bounds_pass,
                  production_complete=residuals_pass and bounds_pass and log_bounds_pass,
                  auditor_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), native_source_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources})
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({key:value for key,value in result.items() if key not in ["steps", "records", "native_source_sha256"]}, indent=2))
    if not result["production_checks_pass"]:
        raise SystemExit(2)

if __name__ == "__main__":
    main()
