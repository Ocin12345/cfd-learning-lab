"""Check saved channel OpenFOAM outputs without starting solvers or rendering.

The 750-to-800 stability check compares two saved endpoints, not every iteration.
"""
from pathlib import Path
import argparse
import json
import re

import numpy as np

from plot_channel_flow import foam_field, map_and_check, require, sha256


def clean(path):
    return re.sub(r"/\*.*?\*/|//[^\n]*", "", path.read_text(encoding="utf-8"), flags=re.S)


def boundary_flux(boundary, phi):
    patches = {}
    for name, body in re.findall(r"\b(\w+)\s*\{([^{}]*)\}", clean(boundary)):
        count = re.search(r"\bnFaces\s+(\d+)\s*;", body)
        if count:
            patches[name] = dict(count=int(count[1]),
                                 kind=re.search(r"\btype\s+(\w+)\s*;", body)[1],
                                 start=int(re.search(r"\bstartFace\s+(\d+)\s*;", body)[1]))
    require(set(patches) == {"inlet", "outlet", "walls", "frontAndBackPlanes"},
            "Unexpected native boundary patch set")
    require([patches[n]["count"] for n in ("inlet", "outlet", "walls")] == [32, 32, 480],
            "Native channel boundary counts differ")
    require(patches["walls"]["kind"] == "wall"
            and patches["frontAndBackPlanes"]["kind"] == "empty",
            "Expected no-slip wall and empty depth patches")
    # Also authenticate internal phi dimensions/count; only boundary values enter net flow.
    foam_field(phi, min(p["start"] for p in patches.values()), 1, [0, 3, -1, 0, 0, 0, 0])
    raw = clean(phi).split("boundaryField", 1)[1]
    blocks = dict(re.findall(r"\b(\w+)\s*\{([^{}]*)\}", raw))
    flow = {}
    for name, patch in patches.items():
        require(name in blocks, "Boundary phi patch missing: " + name)
        if patch["kind"] == "empty":
            continue
        uniform = re.search(r"\bvalue\s+uniform\s+([^;]+);", blocks[name])
        if uniform:
            values = np.full(patch["count"], float(uniform[1]))
        else:
            value = re.search(r"\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
                              blocks[name], re.S)
            require(value is not None and int(value[1]) == patch["count"],
                    "Boundary phi count differs: " + name)
            values = np.fromstring(value[2], sep=" ")
        require(len(values) == patch["count"] and np.isfinite(values).all(),
                "Invalid boundary phi: " + name)
        flow[name] = float(values.sum())
    require(flow["inlet"] < 0 < flow["outlet"], "Expected inlet inflow and outlet outflow")
    net = sum(flow.values())
    return dict(signed_volume_flux_m3_s=flow, net_outward_volume_flux_m3_s=net,
                relative_imbalance=abs(net) / abs(flow["inlet"]),
                definition="abs(sum signed outward phi for inlet,outlet,walls)/abs(inlet phi sum)",
                empty_patches_excluded=["frontAndBackPlanes"])


def residuals(path):
    sections = re.split(r"^Time = ([\d.eE+-]+)s?\s*$", path.read_text(encoding="utf-8"), flags=re.M)
    samples = []
    for i in range(1, len(sections), 2):
        row = {}
        for name, value in re.findall(r"Solving for (Ux|Uy|p), Initial residual = ([\d.eE+-]+)",
                                      sections[i + 1]):
            row[name] = max(row.get(name, 0.), float(value))
        require(set(row) == {"Ux", "Uy", "p"}, "Missing printed native initial residual")
        samples.append((float(sections[i]), row))
    require(samples and samples[-1][0] == 800, "Expected final steady iteration800")
    require(len(samples) == 800, "Expected one bounded800-iteration execution")
    last50 = [row for time, row in samples if 750 < time <= 800]
    require(len(last50) == 50, "Incomplete final50 residual samples")
    return dict(final_iteration=800, final_initial_residuals=samples[-1][1],
                maximum_initial_residuals_last50={n: max(r[n] for r in last50) for n in ("Ux", "Uy", "p")},
                definition="Maximum printed initial residual per variable within each steady iteration",
                cross_solver_residual_norm_equivalence=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    manifest_path = out.with_name("openfoam-source-manifest.json")
    require(not out.exists() and not manifest_path.exists(), "Choose fresh check output paths")
    case = run / "openfoam"
    setup_path, declared_path = run / "mesh/setup.json", case / "declared-comparison.json"
    boundary = case / "0/polyMesh/boundary"
    script = Path(__file__).resolve()
    helper = script.with_name("plot_channel_flow.py")
    inputs = [setup_path, declared_path, run / "mesh/channel.msh", case / "input-channel.msh",
              case / "mesh-check.log", case / "solver.log", case / "run-status.json", script, helper,
              *(case / "0/polyMesh" / n for n in ("boundary", "points", "faces", "owner", "neighbour")),
              case / "800/C", *(case / str(t) / n for t in (750, 800) for n in ("U", "p", "phi"))]
    before = {str(p): sha256(p) for p in inputs}
    setup = json.loads(setup_path.read_text(encoding="utf-8"))
    declared = json.loads(declared_path.read_text(encoding="utf-8"))
    require(before[str(run / "mesh/channel.msh")] == before[str(case / "input-channel.msh")]
            == setup["mesh_sha256"] == declared["mesh_source_sha256"], "Shared mesh hash mismatch")
    for name, expected in dict(nx=240, ny=32, cells=7680, length_m=1., height_m=.1,
                               density_kg_m3=1., dynamic_viscosity_Pa_s=.01,
                               kinematic_viscosity_m2_s=.01, inlet_velocity_m_s=1.).items():
        require(setup[name] == expected, "Unsupported teaching-case input: " + name)
    status = json.loads((case / "run-status.json").read_text(encoding="utf-8"))
    require(status["stage"] == "finished_execution" and status["final_iteration"] == 800,
            "Saved run did not finish the declared cap")
    mesh_log = (case / "mesh-check.log").read_text(encoding="utf-8")
    require("checkMesh -allGeometry -allTopology" in mesh_log and "Mesh OK." in mesh_log,
            "Full native mesh check did not report Mesh OK")
    require(re.search(r"\bcells:\s+7680\b", mesh_log), "Checked native mesh count differs")
    centres = foam_field(case / "800/C", setup["cells"], 3, [0, 1, 0, 0, 0, 0, 0])
    grids, analytical = {}, {}
    for time in (750, 800):
        velocity = foam_field(case / str(time) / "U", setup["cells"], 3, [0, 1, -1, 0, 0, 0, 0])
        pressure = foam_field(case / str(time) / "p", setup["cells"], 1, [0, 2, -2, 0, 0, 0, 0])
        field = dict(centres=centres[:, :2], u=velocity[:, 0], v=velocity[:, 1],
                     p=pressure * setup["density_kg_m3"],
                     topology="Native800 C on unchanged static0/polyMesh; no imported-order assumption")
        require(np.ptp(centres[:, 2]) <= 1e-9 and np.max(np.abs(velocity[:, 2])) <= 1e-12,
                "Expected one planar cell layer")
        analytical[time] = map_and_check(field, setup)
        grids[time] = field["grid"]
    x = (np.arange(setup["nx"]) + .5) * setup["length_m"] / setup["nx"]
    column, left, right = [int(np.argmin(abs(x - v))) for v in (.8, .6, .9)]
    profile_difference = (grids[800]["u"][:, column] - grids[750]["u"][:, column]) / setup["inlet_velocity_m_s"]
    drop = {t: float(grids[t]["p"][:, left].mean() - grids[t]["p"][:, right].mean()) for t in (750, 800)}
    require(min(drop.values()) > 0, "Nonpositive downstream pressure drop")
    stability = dict(iterations=[750, 800], profile_x_m=float(x[column]),
                     profile_rms_change_over_Umean=float(np.sqrt(np.mean(profile_difference ** 2))),
                     profile_maximum_change_over_Umean=float(np.max(abs(profile_difference))),
                     pressure_drop_interval_m=[float(x[left]), float(x[right])],
                     pressure_drop_Pa=drop,
                     pressure_drop_relative_change=abs(drop[800] - drop[750]) / abs(drop[800]),
                     limitation="Two saved endpoints only; does not establish every-iteration stability between750 and800")
    flux = boundary_flux(boundary, case / "800/phi")
    residual = residuals(case / "solver.log")
    checks = dict(full_mesh_check_passed=True, shared_source_mesh_hash_matches=True,
                  residual_initial_final_le_1e_6=max(residual["final_initial_residuals"].values()) <= 1e-6,
                  residual_initial_last50_le_1e_6=max(residual["maximum_initial_residuals_last50"].values()) <= 1e-6,
                  relative_mass_imbalance_le_1e_5=flux["relative_imbalance"] <= 1e-5,
                  downstream_profile_endpoint_change_le_0p1percent=stability["profile_maximum_change_over_Umean"] <= .001,
                  pressure_drop_endpoint_change_le_0p1percent=stability["pressure_drop_relative_change"] <= .001,
                  analytical_targets_passed=analytical[800]["analytical_targets_passed"])
    require(all(sha256(p) == before[str(p)] for p in inputs), "Input changed during read-only analysis")
    record = dict(solver="OpenFOAM Foundation14 incompressibleFluid", mesh_cells=7680,
                  residuals=residual, boundary_mass_balance=flux, final50_endpoint_stability=stability,
                  final_analytical_checks=analytical[800], check_results=checks,
                  all_declared_checks_passed=all(checks.values()),
                  mesh_independence_established=False, experimental_validation=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(record, indent=2, allow_nan=False) + "\n")
    manifest = dict(inputs_unchanged_during_analysis=True, sources=[
        dict(path=p, sha256=h) for p, h in before.items()])
    with manifest_path.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0 if all(checks.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())