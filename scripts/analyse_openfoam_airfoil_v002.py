"""Independently check saved matched Fluent/OpenFOAM airfoil records.

No solver, mesher, conversion, GUI or WSL process is launched. Requires numpy
and h5py to check the saved Fluent case against the imported source mesh.
Only ASCII OpenFOAM fields and the prepared single-zone triangular mesh are
supported. Inputs are read-only; outputs require a fresh directory.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import shutil

import h5py
import numpy as np

NUMBER = r"[-+0-9.eE]+"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def dictionary(path):
    return re.sub(r"/\*.*?\*/|//[^\n]*", "", path.read_text(), flags=re.S)


def block(text, name):
    match = re.search(r"(?<![\w])" + re.escape(name) + r"\s*\{", text)
    require(match is not None, "Missing dictionary block: " + name)
    start, depth = match.end(), 1
    for index in range(start, len(text)):
        depth += (text[index] == "{") - (text[index] == "}")
        if depth == 0:
            return text[start:index]
    raise ValueError("Unclosed dictionary block: " + name)


def entry(text, key):
    match = re.search(r"(?<![\w])" + re.escape(key) + r"\s+([^;{}]+);", text)
    require(match is not None, "Missing dictionary entry: " + key)
    return " ".join(match[1].split())


def scalar(text, key):
    value = float(entry(text, key))
    require(math.isfinite(value), "Nonfinite scalar: " + key)
    return value


def vector(text, key):
    raw = entry(text, key)
    require(raw.startswith("(") and raw.endswith(")"), "Expected vector: " + key)
    values = [float(v) for v in raw[1:-1].split()]
    require(len(values) == 3 and all(map(math.isfinite, values)), "Invalid vector")
    return values


def close(actual, expected, description):
    require(math.isclose(actual, expected, rel_tol=1e-10, abs_tol=1e-12),
            description + ": " + repr((actual, expected)))


def numerical_check_summary(check_results, retain_failed_checks=False):
    """Retain failed numerical criteria only by explicit opt-in; never relabel them."""
    require(check_results and all(isinstance(k, str) and type(v) is bool
                                 for k, v in check_results.items()), "Invalid check-result mapping")
    failed = [name for name, passed in check_results.items() if not passed]
    symmetry = [v for k, v in check_results.items() if k.endswith(".zero_angle_lift_diagnostic")]
    convergence = [v for k, v in check_results.items() if not k.endswith(".zero_angle_lift_diagnostic")]
    if failed and not retain_failed_checks:
        raise ValueError("Numerical acceptance checks failed: " + ", ".join(failed)
                         + "; use --retain-failed-checks only for explicitly labelled diagnostics")
    return dict(check_results=check_results, failed_checks=failed,
                numerical_convergence_passed=all(convergence),
                symmetry_diagnostic_passed=all(symmetry),
                numerical_checks_passed=not failed)


def semantic_tokens(raw):
    """Ignore comments, FoamFile metadata, whitespace and numeric spelling."""
    raw = re.sub(r"\bFoamFile\s*\{[^}]*\}", "", raw, flags=re.S)
    tokens = re.findall(r'"[^"]*"|[^\s{}();]+|[{}();]', raw)
    normalized = []
    for token in tokens:
        try:
            value = float(token)
        except ValueError:
            normalized.append(token)
        else:
            require(math.isfinite(value), "Nonfinite method token")
            normalized.append(repr(value))
    return normalized


def control_invariants(raw):
    """Permit only root restart/limit/output-interval edits, not function edits."""
    for key in ["startFrom", "startTime", "endTime", "writeInterval"]:
        raw, count = re.subn(r"(?<![\w])" + key + r"\s+[^;{}]+;",
                             key + " EXECUTION_CONTROL;", raw, count=1)
        require(count == 1, "Missing execution-control entry: " + key)
    return semantic_tokens(raw)


def execution_control_records(case, targets):
    """Verify preserved control dictionaries for every limit-only extension."""
    path = case / "execution-continuations.json"
    if not path.exists():
        return [], []
    records = json.loads(path.read_text())
    require(isinstance(records, list) and records, "Expected a nonempty continuation record list")
    declaration_sha = sha256(case / "declared-comparison.json")
    archives = []
    for record in records:
        relative = Path(record["archived_control_path"])
        require(not relative.is_absolute() and ".." not in relative.parts
                and relative.parent == Path("system")
                and relative.name.startswith("controlDict."),
                "Archived control path must stay under the native case system directory")
        archive = case / relative
        require(sha256(archive) == record["archived_control_sha256"],
                "Archived continuation control hash differs")
        require(record["declared_comparison_sha256"] == declaration_sha,
                "Continuation changed declared comparison/criteria")
        require(isinstance(record["reason"], str) and record["reason"].strip(),
                "Continuation must record its reason")
        for key in ["prior_iteration", "end_iteration"]:
            require(isinstance(record[key], int) and not isinstance(record[key], bool)
                    and record[key] > 0, "Continuation requires positive integer limits")
        require(record["end_iteration"] > record["prior_iteration"],
                "Continuation end must exceed its prior limit")
        archives.append(archive)
    method_path = case / "method-continuation-stage02.json"
    if method_path.exists():
        method = json.loads(method_path.read_text())
        require(method["stage02_end_iteration"] == records[0]["prior_iteration"]
                and method["current_controlDict_sha256"] == records[0]["archived_control_sha256"],
                "Continuation chain does not start from the recorded stage02 control")
    else:
        require(records[0]["prior_iteration"] == targets["initial_iteration_cap"],
                "Continuation chain does not start from the initial declared cap")
    for index, record in enumerate(records):
        previous = dictionary(archives[index])
        next_path = archives[index + 1] if index + 1 < len(records) else case / "system/controlDict"
        current = dictionary(next_path)
        require(sha256(next_path) == record["current_control_sha256"],
                "Continuation current-control SHA chain differs")
        require(scalar(previous, "endTime") == record["prior_iteration"]
                and scalar(current, "endTime") == record["end_iteration"],
                "Continuation control limits differ from recorded prior/end iterations")
        if index:
            require(records[index - 1]["end_iteration"] == record["prior_iteration"]
                    and records[index - 1]["current_control_sha256"]
                    == record["archived_control_sha256"],
                    "Continuation prior/current SHA or iteration chain is broken")
        require(control_invariants(previous) == control_invariants(current),
                "Continuation changed control functions/physics beyond execution controls")
        require(entry(current, "stopAt") == "endTime", "Continuation does not stop at declared end")
        start_from = entry(current, "startFrom")
        require(start_from == "latestTime"
                or (start_from == "startTime"
                    and scalar(current, "startTime") == record["prior_iteration"]),
                "Continuation is not configured to resume its saved checkpoint")
        checkpoint = [case / str(record["prior_iteration"]) / field for field in ["U", "p", "phi"]]
        require(all(p.is_file() for p in checkpoint), "Continuation checkpoint fields not preserved")
    inputs = [path, *archives]
    inputs += [case / str(record["prior_iteration"]) / field
               for record in records for field in ["U", "p", "phi"]]
    return records, inputs


def stage_method_record(case, targets):
    """Accept only the documented .7 -> .9 continuation, preserving its checkpoint."""
    path = case / "method-continuation-stage02.json"
    if not path.exists():
        return None, []
    record = json.loads(path.read_text())
    prior = int(record["prior_iteration"])
    require(prior == record["prior_iteration"] == targets["initial_iteration_cap"] == 2000,
            "Stage02 must retain the original 2000-iteration checkpoint")
    require(targets["U_equation_relaxation"] == .7
            and record["U_equation_relaxation"] == .9
            and record["n_nonorthogonal_correctors"] == targets["n_nonorthogonal_correctors"] == 2,
            "Unsupported stage02 method change")
    old_solution = case / "system/fvSolution.stage01"
    old_schemes = case / "system/fvSchemes.stage01"
    old_control = case / "system/controlDict.stage01"
    current_solution = case / "system/fvSolution"
    require(sha256(old_solution) == record["prior_fvSolution_sha256"]
            and sha256(current_solution) == record["current_fvSolution_sha256"]
            and sha256(case / "declared-comparison.json") == record["declared_comparison_sha256"],
            "Stage02 solution/declaration hash differs")
    if "prior_controlDict_sha256" in record:
        require(sha256(old_control) == record["prior_controlDict_sha256"],
                "Stage02 archived initial-control hash differs")
    if "accuracy_criteria_unchanged" in record:
        require(record["accuracy_criteria_unchanged"] is True, "Stage02 changes acceptance criteria")
    previous, current = dictionary(old_solution), dictionary(current_solution)
    for raw, expected in [(previous, .7), (current, .9)]:
        close(scalar(block(block(raw, "relaxationFactors"), "equations"), "U"),
              expected, "Stage02 relaxation")
        require(scalar(block(raw, "SIMPLE"), "nNonOrthogonalCorrectors") == 2,
                "Stage02 changed pressure corrector count")
    def without_relaxation(raw):
        result, count = re.subn(r"\bU\s+" + NUMBER + r"\s*;", "U STAGE_RELAXATION;", raw)
        require(count == 1, "Ambiguous U relaxation entry")
        return semantic_tokens(result)
    require(without_relaxation(previous) == without_relaxation(current),
            "Stage02 changed fvSolution beyond U relaxation")
    current_schemes = case / "system/fvSchemes"
    if "prior_fvSchemes_sha256" in record or "current_fvSchemes_sha256" in record:
        require(sha256(old_schemes) == record["prior_fvSchemes_sha256"]
                and sha256(current_schemes) == record["current_fvSchemes_sha256"],
                "Stage02 prior/current fvSchemes hash differs")
    require(semantic_tokens(dictionary(old_schemes)) == semantic_tokens(dictionary(current_schemes)),
            "Stage02 changed discretization schemes")
    previous_control = dictionary(old_control)
    executions, execution_inputs = execution_control_records(case, targets)
    method_control_path = (case / executions[0]["archived_control_path"] if executions
                           else case / "system/controlDict")
    current_control = dictionary(method_control_path)
    if ("current_controlDict_sha256" in record or "stage02_end_iteration" in record or executions):
        require(sha256(method_control_path) == record["current_controlDict_sha256"]
                and scalar(current_control, "endTime") == record["stage02_end_iteration"] == 2500,
                "Stage02 control hash/end limit differs from the preserved method-stage control")
    require(scalar(previous_control, "endTime") == prior
            and entry(previous_control, "startFrom") == "startTime"
            and scalar(previous_control, "startTime") == 0,
            "Archived control is not the original cold-start branch")
    require(semantic_tokens(block(block(previous_control, "functions"), "forceCoefficients"))
            == semantic_tokens(block(block(current_control, "functions"), "forceCoefficients")),
            "Stage02 changed force definitions")
    start_from = entry(current_control, "startFrom")
    require(start_from == "latestTime"
            or (start_from == "startTime" and scalar(current_control, "startTime") == prior),
            "Stage02 is not configured to resume the saved checkpoint")
    require(scalar(current_control, "endTime") > prior, "Stage02 has no continuation interval")
    checkpoint = [case / str(prior) / field for field in ["U", "p", "phi"]]
    require(all(p.is_file() for p in checkpoint), "Original cold-start checkpoint not preserved")
    return record, [path, old_solution, old_schemes, old_control, *checkpoint, *execution_inputs]


def read_source_mesh(path):
    text = path.read_text()
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("(10 (1 ")) + 1
    stop = next(i for i in range(start, len(lines)) if lines[i] == "))")
    xy = np.array([list(map(float, line.split())) for line in lines[start:stop]])
    require(xy.ndim == 2 and xy.shape[1] == 2 and np.isfinite(xy).all(),
            "Expected finite planar Fluent ASCII coordinates")
    header = re.search(r"^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)", text, re.M)
    require(header is not None, "Unsupported source cell header")
    count = int(header[1], 16)
    faces, patch_counts = [], {}
    for match in re.finditer(
            r"^\(13 \(([0-9a-fA-F]+) ([0-9a-fA-F]+) ([0-9a-fA-F]+) "
            r"([0-9a-fA-F]+) 2\)\(\s*\n(.*?)^\)\)", text, re.M | re.S):
        zone, first, last, _ = [int(v, 16) for v in match.groups()[:4]]
        rows = [[int(v, 16) for v in line.split()] for line in match[5].splitlines()
                if line.strip()]
        require(len(rows) == last - first + 1 and all(len(row) == 4 for row in rows),
                "Source face block count differs")
        faces.extend([row[:2] for row in rows])
        patch_counts[zone] = len(rows)
    require(set(patch_counts) == {3, 4, 5, 6}, "Unexpected source face zones")
    edges = np.array(faces, dtype=np.int64) - 1
    require(np.all(edges >= 0) and np.all(edges < len(xy)), "Bad source face nodes")
    return xy, edges, count, dict(airfoil=patch_counts[4], inlet=patch_counts[5],
                                  outlet=patch_counts[6])


def residual_history(logs):
    rows = {}
    for path in logs:
        text = path.read_text(errors="replace")
        require("FOAM FATAL" not in text and re.search(r"^End\s*$", text, re.M),
                "Unfinished or fatal OpenFOAM log: " + path.name)
        for time, body in re.findall(r"^Time = (" + NUMBER + r")s?\s*\n(.*?)"
                                     r"(?=^Time = |\Z)", text, re.M | re.S):
            k = float(time)
            require(k.is_integer(), "Only unit steady iteration indices are supported")
            found = {}
            for name, value in re.findall(
                    r"Solving for (Ux|Uy|p), Initial residual = (" + NUMBER + ")", body):
                found.setdefault(name, []).append(float(value))
            if not found:
                continue
            require(set(found) == {"Ux", "Uy", "p"}, "Incomplete residual iteration")
            # Use the maximum, not the final (often much smaller) pressure solve.
            row = dict(iteration=int(k), **{key: max(values) for key, values in found.items()})
            require(all(math.isfinite(v) and v >= 0 for v in row.values()),
                    "Invalid residual values")
            if int(k) in rows:
                require(rows[int(k)] == row, "Conflicting overlapping solver logs")
            rows[int(k)] = row
    return [rows[k] for k in sorted(rows)]


def force_history(directory):
    rows = {}
    files = sorted(directory.glob("*/forceCoeffs.dat"), key=lambda p: float(p.parent.name))
    require(files, "No OpenFOAM native force files")
    for path in files:
        for line in path.read_text().splitlines():
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            values = list(map(float, line.split()))
            require(len(values) >= 4 and all(map(math.isfinite, values)),
                    "Invalid native force row")
            require(values[0].is_integer(), "Only steady unit iteration indices supported")
            k = int(values[0])
            row = dict(iteration=k, Cm=values[1], Cd=values[2], Cl=values[3])
            if k in rows:
                for name in ["Cm", "Cd", "Cl"]:
                    close(row[name], rows[k][name], "Conflicting restart force sample")
            rows[k] = row
    history = [rows[k] for k in sorted(rows)]
    require(history and history[0]["iteration"] in [0, 1], "Missing force startup")
    require([r["iteration"] for r in history]
            == list(range(history[0]["iteration"], history[-1]["iteration"] + 1)),
            "Force history has missing iteration samples")
    return history, files


def boundary_flux(phi_path, patch_counts):
    text = dictionary(phi_path)
    require(entry(text, "dimensions") == "[0 3 -1 0 0 0 0]",
            "Expected volumetric phi, not a density-weighted mass flux")
    boundary = block(text, "boundaryField")
    flux = {}
    for name, count in patch_counts.items():
        body = block(boundary, name)
        uniform = re.search(r"\bvalue\s+uniform\s+(" + NUMBER + r")\s*;", body)
        if uniform:
            values = [float(uniform[1])] * count
        else:
            match = re.search(r"\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*"
                              r"\((.*?)\)\s*;", body, re.S)
            require(match is not None, "Missing saved boundary phi: " + name)
            values = list(map(float, match[2].split()))
            require(len(values) == int(match[1]) == count, "Boundary phi count differs")
        require(all(map(math.isfinite, values)), "Nonfinite boundary phi")
        flux[name] = math.fsum(values)
    require(flux["inlet"] < 0 and flux["outlet"] > 0, "Unexpected signed through-flow")
    return flux


def check_field_conditions(u, p, targets):
    require(entry(u, "dimensions") == "[0 1 -1 0 0 0 0]"
            and entry(p, "dimensions") == "[0 2 -2 0 0 0 0]", "Wrong U/p dimensions")
    for name, expected_u, expected_p in [
            ("inlet", "fixedValue", "zeroGradient"),
            ("outlet", "zeroGradient", "fixedValue"),
            ("airfoil", "fixedValue", "zeroGradient")]:
        ub, pb = block(block(u, "boundaryField"), name), block(block(p, "boundaryField"), name)
        require(entry(ub, "type") == expected_u and entry(pb, "type") == expected_p,
                "Mismatched U/p boundary type: " + name)
        if name in ["inlet", "airfoil"]:
            match = re.fullmatch(r"uniform\s+\(([^)]+)\)", entry(ub, "value"))
            require(match is not None, "Expected uniform boundary velocity")
            velocity = list(map(float, match[1].split()))
            expected = [targets["velocity_m_s"], 0, 0] if name == "inlet" else [0, 0, 0]
            require(velocity == expected, "Wrong boundary velocity: " + name)
        if name == "outlet":
            require(entry(pb, "value") == "uniform 0", "Wrong outlet kinematic pressure")
    for text in [u, p]:
        boundary = block(text, "boundaryField")
        empty_name = ("frontAndBackPlanes" if re.search(r"\bfrontAndBackPlanes\s*\{", boundary)
                      else '"(frontAndBack.*|defaultFaces)"')
        empty = block(boundary, empty_name)
        require(entry(empty, "type") == "empty", "Missing empty field condition")


def check_openfoam_settings(case, targets):
    properties = dictionary(case / "constant/physicalProperties")
    require(entry(properties, "viscosityModel") == "constant", "Nonconstant viscosity")
    close(scalar(properties, "nu"), targets["kinematic_viscosity_m2_s"], "OpenFOAM nu")
    require(entry(dictionary(case / "constant/momentumTransport"), "simulationType")
            == "laminar", "Nonlaminar OpenFOAM model")
    control = dictionary(case / "system/controlDict")
    require(entry(control, "solver") == "incompressibleFluid", "Unexpected OF solver")
    force = block(block(control, "functions"), "forceCoefficients")
    for key, expected in [("type", "forceCoeffs"),
                          ("rho", "rhoInf"), ("p", "p"), ("U", "U")]:
        require(entry(force, key) == expected, "Wrong force setting: " + key)
    require(re.fullmatch(r"\(\s*airfoil\s*\)", entry(force, "patches")),
            "Wrong force patch selection")
    for key, target in [("rhoInf", "density_kg_m3"), ("magUInf", "velocity_m_s"),
                         ("lRef", "chord_m"), ("Aref", "reference_area_m2")]:
        close(scalar(force, key), targets[target], "Force reference " + key)
    require(vector(force, "dragDir") == [1, 0, 0]
            and vector(force, "liftDir") == [0, 1, 0], "Wrong force directions")
    check_field_conditions(dictionary(case / "0/U"), dictionary(case / "0/p"), targets)
    schemes = dictionary(case / "system/fvSchemes")
    require(entry(block(schemes, "ddtSchemes"), "default") == "steadyState",
            "Expected steady discretization")
    for group, key, expected in [
            ("gradSchemes", "default", "leastSquares"),
            ("divSchemes", "div(phi,U)", "bounded Gauss linearUpwind grad(U)"),
            ("laplacianSchemes", "default", "Gauss linear corrected"),
            ("snGradSchemes", "default", "corrected")]:
        require(entry(block(schemes, group), key) == expected, "Discretization differs: " + group)
    solution = dictionary(case / "system/fvSolution")
    simple = block(solution, "SIMPLE")
    require(scalar(simple, "nNonOrthogonalCorrectors") == targets["n_nonorthogonal_correctors"]
            and entry(simple, "consistent") == "yes", "SIMPLE differs from declaration")
    stage, _ = stage_method_record(case, targets)
    expected_relaxation = stage["U_equation_relaxation"] if stage else targets["U_equation_relaxation"]
    close(scalar(block(block(solution, "relaxationFactors"), "equations"), "U"),
          expected_relaxation, "U relaxation")
    return control


def check_fluent(run, targets, source_xy, source_edges, source_cells, retain_failed_checks=False):
    settings = json.loads((run / "loaded-settings.json").read_text())
    geometry = json.loads((run / "geometry-setup.json").read_text())
    criteria = json.loads((run / "declared-targets.json").read_text())
    summary = json.loads((run / "summary.json").read_text())
    solver = settings["general"]["solver"]
    require(solver["time"] == "steady" and solver["two_dim_space"] == "planar"
            and solver["type"] == "pressure-based", "Wrong Fluent solver/model")
    require(not settings["general"]["operating_conditions"]["gravity"]["enable"],
            "Unexpected Fluent gravity")
    require(settings["viscous"]["model"] == "laminar", "Nonlaminar Fluent model")
    rho, mu = settings["material"]["density"], settings["material"]["viscosity"]
    require(rho["option"] == mu["option"] == "value", "Nonconstant Fluent material")
    close(rho["value"], targets["density_kg_m3"], "Fluent density")
    close(mu["value"], targets["dynamic_viscosity_Pa_s"], "Fluent viscosity")
    momentum = settings["inlet"]["momentum"]
    require(momentum["velocity_specification_method"] == "Components", "Wrong inlet method")
    require(all(v["option"] == "value" for v in momentum["velocity_components"]),
            "Nonconstant inlet components")
    require([v["value"] for v in momentum["velocity_components"]]
            == [targets["velocity_m_s"], 0], "Wrong Fluent inlet velocity")
    require(settings["outlet"]["momentum"]["gauge_pressure"]["value"]
            == targets["pressure_reference_Pa_gauge"], "Wrong outlet gauge pressure")
    for key, target in [("density", "density_kg_m3"), ("velocity", "velocity_m_s"),
                         ("length", "chord_m"), ("area", "reference_area_m2"),
                         ("depth", "span_m"), ("pressure", "pressure_reference_Pa_gauge")]:
        close(settings["references"][key], targets[target], "Fluent reference " + key)
    close(rho["value"] * targets["velocity_m_s"] * targets["chord_m"] / mu["value"],
          targets["reynolds_number"], "Actual Fluent Reynolds number")
    require(geometry["cells"] == source_cells and geometry["angle_degrees"] == targets["angle_degrees"],
            "Fluent geometry metadata differs")
    case_path, data_path = run / "airfoil-Re100-final.cas.h5", run / "airfoil-Re100-final.dat.h5"
    require(case_path.is_file() and data_path.is_file(), "Native Fluent saved pair missing")
    with h5py.File(case_path, "r") as file:
        mesh = file["meshes/1"]
        require(len(mesh["cells/zoneTopology/minId"]) == 1
                and int(mesh["cells/zoneTopology/minId"][0]) == 1
                and int(mesh["cells/zoneTopology/maxId"][0]) == source_cells,
                "Native Fluent cell count/zone differs")
        xy = mesh["nodes/coords/1"][:]
        require(xy.shape == source_xy.shape, "Native Fluent node count differs")
        audit_path = Path(__file__).with_name("plot_openfoam_airfoil_fields_v002.py")
        audit_spec = importlib.util.spec_from_file_location("fluent_coordinate_bijection", audit_path)
        audit_module = importlib.util.module_from_spec(audit_spec)
        audit_spec.loader.exec_module(audit_module)
        native_to_source, coordinate_metrics = audit_module.coordinate_bijection(source_xy, xy, 1e-10)
        maximum_xy = coordinate_metrics["maximum_euclidean_error_m"]
        require(int(mesh["nodes/coords/1"].attrs["minId"][0]) == 1,
                "Unsupported native Fluent node ID offset")
        require(np.all(mesh["faces/nodes/1/nnodes"][:] == 2), "Non-edge Fluent faces")
        native_edges = mesh["faces/nodes/1/nodes"][:].reshape(-1, 2).astype(np.int64) - 1
        native_edges = native_to_source[native_edges]
        def keys(edges):
            return np.min(edges, axis=1) * len(xy) + np.max(edges, axis=1)
        require(np.array_equal(np.sort(keys(native_edges)), np.sort(keys(source_edges))),
                "Saved Fluent face connectivity differs from imported mesh")
        thread = file["settings/Thread Variables"][0].decode()
        wall = re.search(r"\(39 \(\d+ wall airfoil \d+\)\((.*?)(?=\n\(39 |\Z)",
                         thread, re.S)
        require(wall is not None, "Saved airfoil zone is not a wall")
        for raw in [r"\(motion-bc\s+\.\s+0\)", r"\(shear-bc\s+\.\s+0\)",
                    r"\(moving\?\s+\.\s+#f\)"]:
            require(re.search(raw, wall[1]), "Unsupported moving/slip wall flags")
    with h5py.File(data_path, "r") as file:
        fields = file["results/1/phase-1/cells"]
        for name in ["SV_U", "SV_V", "SV_P", "SV_DENSITY", "SV_MU_LAM"]:
            dataset = fields[name + "/1"]
            require(dataset.shape == (source_cells,)
                    and int(dataset.attrs["minId"][0]) == 1
                    and int(dataset.attrs["maxId"][0]) == source_cells,
                    "Saved Fluent data cell ranges differ: " + name)
            values = dataset[:]
            require(np.isfinite(values).all(), "Nonfinite saved Fluent field: " + name)
            if name in ["SV_DENSITY", "SV_MU_LAM"]:
                expected = rho["value"] if name == "SV_DENSITY" else mu["value"]
                require(np.allclose(values, expected, rtol=1e-10, atol=1e-12),
                        "Saved Fluent data material differs: " + name)
    with (run / "force-history.csv").open(newline="") as stream:
        forces = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(stream)]
    require([r["additional_iteration"] for r in forces]
            == list(range(criteria["additional_iterations"] + 1)), "Incomplete Fluent force history")
    require(all(math.isfinite(v) for row in forces for v in row.values()), "Nonfinite Fluent forces")
    for row in forces:
        n = int(row["additional_iteration"])
        for direction in ["drag", "lift"]:
            report = (run / f"{direction}-{n:03d}.txt").read_text()
            values = re.findall(r"^\s*Net\s+((?:" + NUMBER + r"\s+){5}" + NUMBER + r")\s*$",
                                report, re.M)
            require(values, "Missing Fluent scalar force report")
            value = list(map(float, values[-1].split()))
            require(value[2] == row[direction + "_force_N"]
                    and value[5] == row[direction + "_coefficient"], "Fluent force report differs")
            if direction == "drag":
                require(value[0] == row["pressure_drag_N"] and value[1] == row["viscous_drag_N"],
                        "Fluent drag components differ")
    require(forces[-1] == summary["final"], "Stale Fluent force summary")
    window = forces[-criteria["last_window_samples"]:]
    require(len(window) == criteria["last_window_samples"], "Incomplete Fluent force window")
    cd, cl = [r["drag_coefficient"] for r in window], [r["lift_coefficient"] for r in window]
    cd_range = (max(cd) - min(cd)) / abs(math.fsum(cd) / len(cd))
    cl_range = max(cl) - min(cl)
    text = (run / "solver-transcript.txt").read_text(errors="replace")
    residuals = re.findall(r"^\s*(\d+)\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s+(" + NUMBER + ")",
                          text, re.M)
    require(residuals, "No Fluent native residual rows")
    final_residual = list(map(float, residuals[-1][1:]))
    require(all(math.isfinite(v) and v >= 0 for v in final_residual),
            "Invalid Fluent residual evidence")
    quality = float(re.search(r"Minimum Orthogonal Quality\s*=\s*(" + NUMBER + ")", text)[1])
    require(quality >= criteria["minimum_orthogonal_quality"], "Fluent mesh quality failed")
    mass = (run / "mass-balance.txt").read_text()
    net = float(re.search(r"^\s*Net\s+(" + NUMBER + ")", mass, re.M)[1])
    inlet = float(re.search(r"^\s*inlet\s+(" + NUMBER + ")", mass, re.M)[1])
    require(math.isfinite(inlet) and inlet != 0, "Invalid Fluent inlet mass flow")
    imbalance = abs(net / inlet)
    require(math.isfinite(imbalance), "Invalid Fluent mass evidence")
    checks = numerical_check_summary({
        "fluent.residuals": max(final_residual) <= criteria["residual_target"],
        "fluent.Cd_stability": cd_range <= criteria["drag_relative_range_max"],
        "fluent.Cl_stability": cl_range <= criteria["lift_absolute_range_max"],
        "fluent.zero_angle_lift_diagnostic":
            max(map(abs, cl)) <= criteria["zero_angle_abs_lift_diagnostic_target"],
        "fluent.mass_balance": imbalance <= criteria["relative_mass_imbalance_max"]
    }, retain_failed_checks)
    return dict(final_iteration=int(residuals[-1][0]), Cd=forces[-1]["drag_coefficient"],
                Cl=forces[-1]["lift_coefficient"], final_scaled_residuals=final_residual,
                force_window_samples=len(window), Cd_relative_range=cd_range,
                Cl_absolute_range=cl_range, maximum_absolute_Cl_in_window=max(map(abs, cl)),
                relative_mass_imbalance=imbalance,
                minimum_orthogonal_quality=quality, mesh_xy_max_difference_m=maximum_xy,
                native_mesh_connectivity_verified=True, saved_data_cell_ranges_and_material_verified=True,
                saved_stationary_noslip_wall_flags_verified=True,
                numerical_methods=settings["methods"], acceptance_thresholds=criteria,
                **checks), forces, residuals


def analyse(case, fluent_run, retain_failed_checks=False):
    targets = json.loads((case / "declared-comparison.json").read_text())
    status = json.loads((case / "run-status.json").read_text())
    require(status["stage"] == "finished_execution", "OpenFOAM execution is not marked finished")
    imported = case / "input-airfoil.msh"
    require(sha256(imported) == targets["mesh_source_sha256"], "Imported mesh hash differs")
    source = Path(targets["mesh_source_path"])
    require(sha256(source) == targets["mesh_source_sha256"], "Original mesh hash differs")
    require(sha256(case / "source-setup.json") == targets["mesh_setup_sha256"], "Setup hash differs")
    preparation = json.loads((case / "preparation-provenance.json").read_text())
    require(sha256(case / "executed-prepare-script.py") == preparation["script_sha256"]
            and preparation["mesh_source_sha256"] == targets["mesh_source_sha256"]
            and preparation["mesh_setup_sha256"] == targets["mesh_setup_sha256"],
            "Preparation snapshot/provenance differs")
    source_setup = json.loads((case / "source-setup.json").read_text())
    require(json.loads((fluent_run / "geometry-setup.json").read_text()) == source_setup,
            "Fluent run is not from the selected source mesh setup")
    for key, target in [("density_kg_m3", "density_kg_m3"), ("viscosity_Pa_s", "dynamic_viscosity_Pa_s"),
                         ("inlet_velocity_m_s", "velocity_m_s"), ("chord_m", "chord_m"),
                         ("reynolds_number", "reynolds_number"), ("angle_degrees", "angle_degrees")]:
        close(source_setup[key], targets[target], "Source setup " + key)
    close(targets["velocity_m_s"] * targets["chord_m"] / targets["kinematic_viscosity_m2_s"],
          targets["reynolds_number"], "Actual OpenFOAM Reynolds number")
    close(targets["dynamic_viscosity_Pa_s"] / targets["density_kg_m3"],
          targets["kinematic_viscosity_m2_s"], "Dynamic/kinematic viscosity")
    source_xy, source_edges, source_cells, patch_counts = read_source_mesh(imported)
    require(source_cells == targets["intended_cells"] == source_setup["cells"]
            and len(source_xy) == targets["intended_nodes"] == source_setup["nodes"],
            "Source mesh count differs from declaration")
    control = check_openfoam_settings(case, targets)
    method_stage, method_inputs = stage_method_record(case, targets)
    execution_stages, execution_inputs = execution_control_records(case, targets)
    mesh_log = (case / "mesh-check.log").read_text()
    require("Mesh OK." in mesh_log and "-allGeometry -allTopology" in mesh_log
            and "FOAM FATAL" not in mesh_log, "Full native OpenFOAM mesh check did not pass")
    cells = int(re.search(r"^\s*cells:\s+(\d+)", mesh_log, re.M)[1])
    require(cells == source_cells, "Converted cell count differs")
    mesh = case / "0/polyMesh"
    boundary = dictionary(mesh / "boundary")
    names = [name for name in re.findall(r"^\s*([A-Za-z]\w*)\s*\{", boundary, re.M)
             if name != "FoamFile"]
    require(set(names) == {"airfoil", "inlet", "outlet", "frontAndBackPlanes"},
            "Unexpected converted boundary patches")
    for name, count in patch_counts.items():
        patch = block(boundary, name)
        require(scalar(patch, "nFaces") == count, "Converted patch face count differs")
        require(entry(patch, "type") == ("wall" if name == "airfoil" else "patch"),
                "Converted patch type differs")
    empty = block(boundary, "frontAndBackPlanes")
    require(entry(empty, "type") == "empty" and scalar(empty, "nFaces") == 2 * cells,
            "Extrusion has missing/nonempty front/back faces")
    points = np.array([list(map(float, values)) for values in re.findall(
        r"^\(\s*(" + NUMBER + r")\s+(" + NUMBER + r")\s+(" + NUMBER + r")\s*\)$",
        (mesh / "points").read_text(), re.M)])
    require(points.shape == (2 * len(source_xy), 3), "Wrong extruded node count")
    require(set(np.round(points[:, 2], 12)) == {-targets["span_m"] / 2, targets["span_m"] / 2},
            "Wrong extrusion span or plane positions")
    maximum_xy = float(np.max(np.linalg.norm(points[:, :2] - np.tile(source_xy, (2, 1)), axis=1)))
    require(maximum_xy <= 1e-10, "Source-to-extrusion coordinate/order preservation failed")
    topology_script = Path(__file__).with_name("plot_openfoam_airfoil_fields_v002.py")
    spec = importlib.util.spec_from_file_location("openfoam_airfoil_topology_audit", topology_script)
    topology_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(topology_module)
    topology = topology_module.load_case_geometry(case)
    require(topology["metrics"]["cells"] == cells, "Topology audit cell count differs")
    forces, force_files = force_history(case / "postProcessing/forceCoefficients")
    final_k = forces[-1]["iteration"]
    require(final_k >= targets["initial_iteration_cap"], "Initial iteration cap not completed")
    require(status["final_iteration"] == final_k and scalar(control, "endTime") == final_k,
            "Final iteration differs across status/control/force files")
    logs = sorted(case.glob("solver*.log"))
    require(logs, "No native OpenFOAM solver logs")
    residuals = residual_history(logs)
    require([r["iteration"] for r in residuals] == list(range(1, final_k + 1)),
            "Incomplete OpenFOAM residual history")
    final_residual = {key: residuals[-1][key] for key in ["Ux", "Uy", "p"]}
    window = forces[-targets["final_force_window"]:]
    require(len(window) == targets["final_force_window"], "Incomplete OpenFOAM force window")
    cd, cl = [r["Cd"] for r in window], [r["Cl"] for r in window]
    cd_range = (max(cd) - min(cd)) / abs(math.fsum(cd) / len(cd))
    cl_range = max(cl) - min(cl)
    for field in ["p", "U", "phi"]:
        require((case / str(final_k) / field).is_file(), "Reconstructed final field missing")
    check_field_conditions(dictionary(case / str(final_k) / "U"),
                           dictionary(case / str(final_k) / "p"), targets)
    flux = boundary_flux(case / str(final_k) / "phi", patch_counts)
    imbalance = abs(math.fsum(flux.values())) / abs(flux["inlet"])
    require(math.isfinite(imbalance), "Invalid OpenFOAM boundary mass evidence")
    fluent, fluent_forces, fluent_residuals = check_fluent(
        fluent_run, targets, source_xy, source_edges, source_cells, retain_failed_checks)
    check_results = {
        "openfoam.residuals":
            max(final_residual.values()) <= targets["initial_p_U_residual_learning_target"],
        "openfoam.Cd_stability": cd_range <= targets["cd_relative_range_target"],
        "openfoam.Cl_stability": cl_range <= targets["cl_absolute_range_target"],
        "openfoam.zero_angle_lift_diagnostic":
            max(map(abs, cl)) <= targets["zero_angle_abs_lift_diagnostic_target"],
        "openfoam.mass_balance": imbalance <= targets["relative_boundary_mass_imbalance_target"],
        **fluent["check_results"]
    }
    checks = numerical_check_summary(check_results, retain_failed_checks)
    require(fluent["Cd"] != 0, "Cannot normalize cross-code difference by zero Fluent Cd")
    difference = 100 * abs(forces[-1]["Cd"] - fluent["Cd"]) / abs(fluent["Cd"])
    record = dict(reynolds_number=targets["reynolds_number"], angle_degrees=targets["angle_degrees"],
                  cells=cells, span_m=targets["span_m"], source_mesh_sha256=targets["mesh_source_sha256"],
                  mesh_xy_max_difference_m=maximum_xy, coordinate_preservation_tolerance_m=1e-10,
                  native_mesh_topology_checks=topology["metrics"],
                  openfoam=dict(final_iteration=final_k, Cd=forces[-1]["Cd"], Cl=forces[-1]["Cl"],
                      final_initial_residuals=final_residual,
                      pressure_residual_definition=targets["pressure_residual_definition"],
                      force_window_samples=len(window), Cd_relative_range=cd_range,
                      Cl_absolute_range=cl_range, maximum_absolute_Cl_in_window=max(map(abs, cl)),
                      relative_mass_imbalance=imbalance,
                      boundary_volume_flux_m3_s=flux,
                      boundary_mass_flux_kg_s={k: v * targets["density_kg_m3"] for k, v in flux.items()}),
                  fluent=fluent, cross_code_Cd_difference_percent=difference,
                  cross_code_learning_target_percent=targets["cross_code_cd_difference_learning_target_percent"],
                  cross_code_target_pass=difference <= targets["cross_code_cd_difference_learning_target_percent"],
                  **checks, mesh_independence_established=False,
                  domain_independence_established=False, physical_validation=False,
                  diagnostic_retention_requested=bool(retain_failed_checks),
                  diagnostic_result=bool(checks["failed_checks"]),
                  acceptance_thresholds=dict(openfoam={
                      "initial_p_U_residual_learning_target": targets["initial_p_U_residual_learning_target"],
                      "cd_relative_range_target": targets["cd_relative_range_target"],
                      "cl_absolute_range_target": targets["cl_absolute_range_target"],
                      "zero_angle_abs_lift_diagnostic_target": targets["zero_angle_abs_lift_diagnostic_target"],
                      "relative_boundary_mass_imbalance_target": targets["relative_boundary_mass_imbalance_target"]
                  }, fluent=fluent["acceptance_thresholds"]),
                  pressure_conventions=dict(openfoam="Kinematic p in m2/s2",
                                            fluent="Gauge pressure in Pa"),
                  method_continuation_stage02=method_stage,
                  execution_continuations=execution_stages,
                  final_U_equation_relaxation=(method_stage["U_equation_relaxation"] if method_stage
                                               else targets["U_equation_relaxation"]),
                  note="Independent saved-data checks; distinct residual norms and initializations. A cross-code target failure is retained as a measured result, not relabelled as physical validation.")
    inputs = [imported, case / "source-setup.json", case / "declared-comparison.json",
              case / "executed-prepare-script.py", case / "preparation-provenance.json",
              case / "run-status.json", case / "mesh-check.log", mesh / "points", mesh / "boundary",
              *method_inputs,
              *execution_inputs,
              mesh / "faces", mesh / "owner", mesh / "neighbour", topology_script,
              *(case / str(final_k) / p for p in ["p", "U", "phi"]), *logs, *force_files,
              *(case / p for p in ["system/controlDict", "system/fvSchemes", "system/fvSolution",
                                  "constant/physicalProperties", "constant/momentumTransport", "0/p", "0/U"]),
              *(fluent_run / p for p in ["loaded-settings.json", "geometry-setup.json", "declared-targets.json",
                "summary.json", "force-history.csv", "solver-transcript.txt", "mass-balance.txt",
                "airfoil-Re100-final.cas.h5", "airfoil-Re100-final.dat.h5"])]
    inputs += [fluent_run / f"{direction}-{int(row['additional_iteration']):03d}.txt"
               for row in fluent_forces for direction in ["drag", "lift"]]
    inputs = list(dict.fromkeys(inputs))
    manifest = dict(script_sha256=sha256(Path(__file__)),
                    numerical_checks_passed=checks["numerical_checks_passed"],
                    failed_checks=checks["failed_checks"],
                    diagnostic_retention_requested=bool(retain_failed_checks),
                    inputs=[dict(path=str(path.resolve()), sha256=sha256(path)) for path in inputs])
    return record, forces, residuals, fluent_forces, fluent_residuals, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--fluent-run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--retain-failed-checks", action="store_true",
                        help="Save explicitly failed numerical criteria as diagnostics; structural/evidence checks remain mandatory")
    args = parser.parse_args()
    case, fluent, out = (p.resolve() for p in (args.case, args.fluent_run, args.out))
    if out.exists():
        raise FileExistsError("Use a fresh analysis output directory")
    record, forces, residuals, _, fl_residuals, manifest = analyse(
        case, fluent, args.retain_failed_checks)
    out.mkdir(parents=True)
    for name, value in [("comparison.json", record), ("source-manifest.json", manifest)]:
        with (out / name).open("x", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
    for name, rows in [("openfoam-force-history.csv", forces),
                       ("openfoam-initial-residuals.csv", residuals)]:
        with (out / name).open("x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    with (out / "fluent-scaled-residuals.csv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["iteration", "continuity", "Ux", "Uy"])
        writer.writerows(fl_residuals)
    shutil.copy2(fluent / "force-history.csv", out / "fluent-force-history.csv")
    shutil.copy2(case / "declared-comparison.json", out / "declared-comparison.json")
    print(json.dumps(record, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
