"""Check saved backward-step OpenFOAM fields without changing native files."""
from pathlib import Path
import argparse
import json
import re
import numpy as np
from plot_step_flow import foam_field, grid_coordinates, map_and_check, require, sha256, validate_setup


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
    require(set(patches) == {"inlet", "outlet", "walls", "frontAndBackPlanes"}, "Unexpected patch set")
    require([patches[n]["count"] for n in ("inlet", "outlet", "walls", "frontAndBackPlanes")]
            == [20, 40, 500, 17600], "Native step boundary counts differ")
    require(patches["walls"]["kind"] == "wall" and patches["frontAndBackPlanes"]["kind"] == "empty",
            "Expected wall/empty patch types")
    foam_field(phi, min(p["start"] for p in patches.values()), 1, [0, 3, -1, 0, 0, 0, 0])
    blocks = dict(re.findall(r"\b(\w+)\s*\{([^{}]*)\}", clean(phi).split("boundaryField", 1)[1]))
    flow = {}
    for name in ("inlet", "outlet", "walls"):
        require(name in blocks, "Boundary phi patch missing: " + name)
        uniform = re.search(r"\bvalue\s+uniform\s+([^;]+);", blocks[name])
        if uniform:
            values = np.full(patches[name]["count"], float(uniform[1]))
        else:
            value = re.search(r"\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", blocks[name], re.S)
            require(value is not None and int(value[1]) == patches[name]["count"], "Boundary phi count differs")
            values = np.fromstring(value[2], sep=" ")
        require(len(values) == patches[name]["count"] and np.isfinite(values).all(), "Invalid boundary phi")
        flow[name] = float(values.sum())
    require(flow["inlet"] < 0 < flow["outlet"], "Expected inlet inflow and outlet outflow")
    net = sum(flow.values())
    return dict(signed_volume_flux_m3_s=flow, net_outward_volume_flux_m3_s=net,
                relative_imbalance=abs(net) / abs(flow["inlet"]),
                definition="abs(sum outward inlet,outlet,walls phi)/abs(inlet phi sum); empty planes excluded")


def residuals(path):
    sections = re.split(r"^Time = ([\d.eE+-]+)s?\s*$", path.read_text(encoding="utf-8"), flags=re.M)
    samples = []
    for i in range(1, len(sections), 2):
        row = {}
        for name, value in re.findall(r"Solving for (Ux|Uy|p), Initial residual = ([\d.eE+-]+)", sections[i + 1]):
            row[name] = max(row.get(name, 0.), float(value))
        require(set(row) == {"Ux", "Uy", "p"}, "Missing native initial residual")
        samples.append((float(sections[i]), row))
    require([t for t, _ in samples] == list(range(1, 2001)), "Expected exactly iterations 1 through 2000")
    last100 = [row for t, row in samples if 1900 < t <= 2000]
    return dict(final_iteration=2000, final_initial_residuals=samples[-1][1],
                maximum_initial_residuals_last100={n: max(r[n] for r in last100) for n in ("Ux", "Uy", "p")},
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
    case, script = run / "openfoam", Path(__file__).resolve()
    helper = script.with_name("plot_step_flow.py")
    setup_path, declared_path = run / "mesh/setup.json", case / "declared-comparison.json"
    dictionaries = [case / p for p in ("system/controlDict", "system/fvSchemes", "system/fvSolution",
                    "system/decomposeParDict", "constant/physicalProperties", "constant/momentumTransport", "0/U", "0/p")]
    inputs = [setup_path, declared_path, run / "mesh/step.msh", case / "input-step.msh", script, helper,
              *(case / p for p in ("mesh-check.log", "solver.log", "run-status.json")), *dictionaries,
              *(case / "0/polyMesh" / n for n in ("boundary", "points", "faces", "owner", "neighbour")),
              case / "2000/C", *(case / str(t) / n for t in (1900, 2000) for n in ("U", "p", "phi"))]
    before = {p: sha256(p) for p in inputs}
    setup = validate_setup(json.loads(setup_path.read_text(encoding="utf-8")))
    declared = json.loads(declared_path.read_text(encoding="utf-8"))
    require(before[run / "mesh/step.msh"] == before[case / "input-step.msh"]
            == setup["mesh_sha256"] == declared["mesh_source_sha256"], "Shared mesh hash differs")
    require(setup["kinematic_viscosity_m2_s"] == .001, "Expected nu=.001")
    mesh = (run / "mesh/step.msh").read_text(encoding="ascii")
    require(int(re.search(r"^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)", mesh, re.M)[1], 16) == 8800,
            "Shared mesh cell count differs")
    physical, transport, control, decompose = [clean(case / p) for p in (
        "constant/physicalProperties", "constant/momentumTransport", "system/controlDict", "system/decomposeParDict")]
    require(float(re.search(r"\bnu\s+([\d.eE+-]+)\s*;", physical)[1]) == .001
            and "viscosityModel constant;" in physical and "simulationType laminar;" in transport,
            "Actual viscosity/transport dictionary differs")
    for pattern in (r"\bsolver\s+incompressibleFluid;", r"\bendTime\s+2000;", r"\bwriteInterval\s+100;", r"\bpurgeWrite\s+0;"):
        require(re.search(pattern, control), "Actual control dictionary differs: " + pattern)
    require(re.search(r"numberOfSubdomains\s+2;", decompose), "Expected two MPI ranks")
    solution, schemes = clean(case / "system/fvSolution"), clean(case / "system/fvSchemes")
    require("consistent yes;" in solution and re.search(r"\bU\s+0\.9;", solution)
            and len(re.findall(r"tolerance\s+1e-10;", solution)) == 2 and "residualControl" not in solution,
            "Actual fixed-cap SIMPLE/linear solver settings differ")
    require(all(token in schemes for token in ("default steadyState;", "default leastSquares;",
            "bounded Gauss linearUpwind grad(U);", "default Gauss linear corrected;", "default corrected;")),
            "Actual discretisation differs")
    status = json.loads((case / "run-status.json").read_text(encoding="utf-8"))
    require(status["stage"] == "finished_execution" and status["final_iteration"] == 2000, "Run cap unfinished")
    mesh_log = (case / "mesh-check.log").read_text(encoding="utf-8")
    require("checkMesh -allGeometry -allTopology" in mesh_log and "Mesh OK." in mesh_log
            and re.search(r"\bcells:\s+8800\b", mesh_log), "Full native mesh check did not pass")
    min_volume = float(re.search(r"Min volume = (\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", mesh_log)[1])
    nonorth = float(re.search(r"Mesh non-orthogonality Max:\s+([\d.eE+-]+)", mesh_log)[1])
    require(min_volume > 0 and nonorth == 0, "Expected positive orthogonal native cells")
    centres = foam_field(case / "2000/C", 8800, 3, [0, 1, 0, 0, 0, 0, 0])
    fields, mapped, velocities = {}, {}, {}
    require(np.ptp(centres[:, 2]) <= 1e-9, "Expected planar native centres")
    for time in (1900, 2000):
        require(not (case / str(time) / "polyMesh").exists(), "Unexpected endpoint mesh change")
        velocity = foam_field(case / str(time) / "U", 8800, 3, [0, 1, -1, 0, 0, 0, 0])
        pressure = foam_field(case / str(time) / "p", 8800, 1, [0, 2, -2, 0, 0, 0, 0])
        require(np.max(abs(velocity[:, 2])) <= 1e-12, "Expected zero spanwise velocity")
        fields[time] = dict(centres=centres[:, :2], u=velocity[:, 0], v=velocity[:, 1],
                            p=pressure * setup["density_kg_m3"],
                            topology="Native 2000/C mapped onto unchanged static 0/polyMesh for both endpoints")
        mapped[time] = map_and_check(fields[time], setup)
        require(mapped[time]["solid_slots_masked"] == 800, "Expected 800 masked solid slots")
        velocities[time] = velocity
    x = grid_coordinates(setup)[2]
    left, right = [int(np.argmin(abs(x - value))) for value in (.8, 1.8)]
    drops = {t: float(fields[t]["grid"]["p"][:, left].mean() - fields[t]["grid"]["p"][:, right].mean())
             for t in (1900, 2000)}
    dp_change = abs(drops[2000] - drops[1900])
    relative_dp = dp_change / abs(drops[2000]) if abs(drops[2000]) > 1e-12 else None
    component_change = np.max(abs(velocities[2000] - velocities[1900]), axis=0) / setup["inlet_velocity_m_s"]
    stability = dict(iterations=[1900, 2000], maximum_all_fluid_component_change_over_Uin=float(component_change.max()),
                     maximum_component_changes_over_Uin=dict(zip(("Ux", "Uy", "Uz"), component_change.tolist())),
                     downstream_pressure_drop_interval_m=[float(x[left]), float(x[right])],
                     pressure_drop_Pa=drops, pressure_drop_absolute_change_Pa=dp_change,
                     pressure_drop_relative_change=relative_dp,
                     pressure_drop_change_over_rho_Uin_squared=dp_change / (setup["density_kg_m3"] * setup["inlet_velocity_m_s"] ** 2),
                     limitation="Two saved endpoints only; does not establish every-iteration stability")
    flux = {t: boundary_flux(case / "0/polyMesh/boundary", case / str(t) / "phi") for t in (1900, 2000)}
    residual = residuals(case / "solver.log")
    checks = dict(full_mesh_check_passed=True, shared_source_mesh_hash_matches=True,
                  actual_dictionaries_match_declared_settings=True, native_centres_bijective_8800_and_800_solid_masked=True,
                  residual_initial_final_le_1e_6=max(residual["final_initial_residuals"].values()) <= 1e-6,
                  residual_initial_last100_le_1e_6=max(residual["maximum_initial_residuals_last100"].values()) <= 1e-6,
                  relative_mass_imbalance_le_1e_5=all(f["relative_imbalance"] <= 1e-5 for f in flux.values()),
                  maximum_all_fluid_component_endpoint_change_le_0p1percent=float(component_change.max()) <= .001,
                  pressure_drop_endpoint_change_le_0p1percent=(relative_dp <= .001 if relative_dp is not None
                    else stability["pressure_drop_change_over_rho_Uin_squared"] <= .001))
    after = {p: sha256(p) for p in inputs}
    require(before == after, "Source/input changed during read-only analysis")
    record = dict(solver="OpenFOAM Foundation14 incompressibleFluid", mesh_cells=8800,
                  mesh_check=dict(minimum_volume_m3=min_volume, maximum_nonorthogonality_degrees=nonorth),
                  residuals=residual, boundary_mass_balance=flux, final100_endpoint_stability=stability,
                  native_mapping_checks=mapped, check_results=checks, all_declared_checks_passed=all(checks.values()),
                  mesh_independence_established=False, experimental_validation=False,
                  interpretation="Bounded numerical checks on one teaching mesh; no benchmark validation claimed")
    lab = script.parents[3]
    manifest = dict(inputs_unchanged_during_analysis=True, raw_native_logs_copied_to_repository=False,
                    path_base="CFD-Learning-Lab", sources=[dict(path=p.relative_to(lab).as_posix(),
                    sha256_before=before[p], sha256_after=after[p]) for p in inputs])
    out.parent.mkdir(parents=True, exist_ok=True)
    for path, value in ((out, record), (manifest_path, manifest)):
        with path.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0 if all(checks.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
