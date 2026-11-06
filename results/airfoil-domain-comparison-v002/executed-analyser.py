"""Compare saved R20/R40 Fluent/OpenFOAM airfoil records; never run CFD.

Native input hashes, physical references, methods, relaxation schedules and exact
inner-grid preservation must agree before the four-row domain table is written.
Failed convergence, symmetry and code-comparison checks remain visible.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import shutil

for _name in ["TEMP", "TMP", "TMPDIR"]:
    os.environ[_name] = "E:/CFD-Learning-Lab/temp"

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def close(actual, expected, name):
    require(math.isclose(float(actual), float(expected), rel_tol=1e-10, abs_tol=1e-12),
            name + " differs: " + repr((actual, expected)))


def path_key(path):
    return os.path.normcase(str(Path(path).resolve()))


def verify_manifest(directory):
    manifest = json.loads((Path(directory) / "source-manifest.json").read_text())
    require(manifest.get("inputs"), "Native source manifest has no inputs")
    paths = {}
    for row in manifest["inputs"]:
        path = Path(row["path"]).resolve()
        key = path_key(path)
        require(path.is_file() and sha256(path) == row["sha256"],
                "Changed/missing native input: " + str(path))
        require(key not in paths or paths[key][1] == row["sha256"],
                "Conflicting native input hashes")
        paths[key] = (path, row["sha256"])
    for name, expected in manifest.get("outputs", {}).items():
        path = Path(directory) / name
        require(Path(name).name == name and path.is_file() and sha256(path) == expected,
                "Changed/missing derived output: " + str(path))
    return manifest, paths


def one_input(paths, name):
    found = [p for p, _ in paths.values() if p.name == name]
    require(len(found) == 1, "Expected one native input named " + name)
    return found[0]


def semantic_file(path):
    raw = re.sub(r"/\*.*?\*/|//[^\n]*", "", Path(path).read_text(), flags=re.S)
    raw = re.sub(r"\bFoamFile\s*\{[^}]*\}", "", raw, flags=re.S)
    result = []
    for token in re.findall(r'"[^"]*"|[^\s{}();]+|[{}();]', raw):
        try:
            value = float(token)
        except ValueError:
            result.append(token)
        else:
            require(math.isfinite(value), "Nonfinite dictionary token")
            result.append(repr(value))
    return result


def match_reference_values(baseline, extended):
    """Exclude only Fluent's string zone selector label; retain every other entry."""
    require(isinstance(baseline.get("zone"), str) and isinstance(extended.get("zone"), str),
            "Expected string reference-zone selector labels")
    left = {key: value for key, value in baseline.items() if key != "zone"}
    right = {key: value for key, value in extended.items() if key != "zone"}
    require(left == right, "Fluent reference values differ beyond the zone selector label")
    return {"zone": {"baseline": baseline["zone"], "extended": extended["zone"]}}


def relaxation_schedule(record, targets):
    initial = float(targets["U_equation_relaxation"])
    require(0 < initial <= 1, "Invalid initial U relaxation")
    schedule = [dict(start_iteration=0, U_equation_relaxation=initial)]
    stage = record.get("method_continuation_stage02")
    if stage:
        start, factor = int(stage["prior_iteration"]), float(stage["U_equation_relaxation"])
        require(start > 0 and 0 < factor <= 1, "Invalid relaxation continuation")
        schedule.append(dict(start_iteration=start, U_equation_relaxation=factor))
    close(record["final_U_equation_relaxation"], schedule[-1]["U_equation_relaxation"],
          "Final relaxation schedule")
    return schedule


def build_rows(label, radius, record):
    checks = record["check_results"]
    require(checks and all(type(value) is bool for value in checks.values()),
            "Check results must contain explicit booleans")
    failed = sorted(key for key, passed in checks.items() if not passed)
    require(failed == sorted(record["failed_checks"]), "Failed-check list is inconsistent")
    numerical = all(value for key, value in checks.items()
                    if not key.endswith(".zero_angle_lift_diagnostic"))
    symmetry = all(value for key, value in checks.items()
                   if key.endswith(".zero_angle_lift_diagnostic"))
    require(record["numerical_convergence_passed"] == numerical
            and record["symmetry_diagnostic_passed"] == symmetry
            and record["numerical_checks_passed"] == all(checks.values()),
            "Paired acceptance flags are inconsistent")
    if failed:
        require(record.get("diagnostic_retention_requested") and record.get("diagnostic_result"),
                "Failed record lacks explicit diagnostic retention")
    cd_fl, cd_of = record["fluent"]["Cd"], record["openfoam"]["Cd"]
    require(cd_fl != 0, "Cannot normalize a code comparison by zero Fluent drag")
    difference = 100 * abs(cd_of - cd_fl) / abs(cd_fl)
    close(record["cross_code_Cd_difference_percent"], difference, "Code-difference calculation")
    target = record["cross_code_learning_target_percent"]
    require(record["cross_code_target_pass"] == (difference <= target),
            "Code-target flag is inconsistent")
    rows = []
    for code, name in [("fluent", "Fluent"), ("openfoam", "OpenFOAM")]:
        own = {key: value for key, value in checks.items() if key.startswith(code + ".")}
        require(len(own) == 5, "Missing per-solver numerical/symmetry checks")
        values = record[code]
        require(all(math.isfinite(values[key]) for key in ["Cd", "Cl"]), "Nonfinite force")
        rows.append(dict(case=label, solver=name, outer_radius_over_chord=radius,
            cells=record.get("cells"), Cd=values["Cd"], Cl=values["Cl"],
            numerical_convergence_passed=all(v for k, v in own.items()
                if not k.endswith(".zero_angle_lift_diagnostic")),
            symmetry_diagnostic_passed=own[code + ".zero_angle_lift_diagnostic"],
            numerical_checks_passed=all(own.values()),
            failed_checks=sorted(k for k, v in own.items() if not v),
            cross_code_Cd_difference_percent=difference,
            cross_code_target_pass=record["cross_code_target_pass"],
            case_diagnostic_result=bool(failed)))
    return rows


def verify_inner_grid(source, extended):
    """Independently bind old coordinates, cell IDs, wall and interface edges."""
    n, c = len(source["xy"]), source["cells"]
    require(len(extended["xy"]) > n and extended["cells"] > c, "No added annulus")
    require(np.array_equal(extended["xy"][:n], source["xy"]),
            "Original node coordinates/IDs changed")
    require(np.array_equal(extended["triangles"][:c], source["triangles"]),
            "Original cell IDs/triangles changed")
    nn = len(extended["xy"])
    def keys(edges):
        return np.min(edges, axis=1) * nn + np.max(edges, axis=1)
    old_keys, new_keys = keys(source["edges"]), keys(extended["edges"])
    require(len(np.unique(new_keys)) == len(new_keys), "Duplicated extended edge")
    order = np.argsort(new_keys)
    positions = np.searchsorted(new_keys[order], old_keys)
    require(np.all(positions < len(order)), "Original edge missing")
    matched = order[positions]
    require(np.array_equal(new_keys[matched], old_keys), "Original edge missing")
    wall = source["zones"] == 4
    require(np.array_equal(np.sort(new_keys[extended["zones"] == 4]), np.sort(old_keys[wall])),
            "Original wall edge IDs changed")
    interior = source["zones"] == 3
    old_pairs = np.sort(np.column_stack([source["owner"], source["neighbour"]]), axis=1)
    new_pairs = np.sort(np.column_stack([extended["owner"], extended["neighbour"]]), axis=1)
    require(np.all(extended["zones"][matched[interior]] == 3)
            and np.array_equal(new_pairs[matched[interior]], old_pairs[interior]),
            "Original internal face/cell pairs changed")
    interface = np.isin(source["zones"], [5, 6])
    faces = matched[interface]
    require(interface.any() and np.all(extended["zones"][faces] == 3),
            "Former outer boundaries were not converted to interior")
    pairs = new_pairs[faces]
    require(np.all(pairs[:, 0] == source["owner"][interface])
            and np.all(pairs[:, 1] >= c),
            "Former boundary face does not join its original cell to an annulus cell")
    return dict(original_node_coordinates_exact=True, original_node_ids_preserved=True,
        original_cell_ids_and_triangles_exact=True, original_wall_edge_ids_exact=True,
        original_internal_face_cell_pairs_exact=True,
        former_boundary_edges_now_interior=int(interface.sum()),
        original_nodes=n, original_cells=c, added_nodes=len(extended["xy"]) - n,
        added_cells=extended["cells"] - c)


def source_reader():
    path = Path(__file__).with_name("plot_openfoam_airfoil_fields.py")
    spec = importlib.util.spec_from_file_location("airfoil_domain_source_reader", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.read_source_mesh, path


def load_pair(directory):
    directory = Path(directory).resolve()
    record = json.loads((directory / "comparison.json").read_text())
    targets = json.loads((directory / "declared-comparison.json").read_text())
    manifest, paths = verify_manifest(directory)
    case = one_input(paths, "source-setup.json").parent
    native_targets = case / "declared-comparison.json"
    require(path_key(native_targets) in paths
            and json.loads(native_targets.read_text()) == targets, "Public/native declaration differs")
    setup = json.loads((case / "source-setup.json").read_text())
    mesh = Path(targets["mesh_source_path"]).resolve()
    require(sha256(mesh) == targets["mesh_source_sha256"] == record["source_mesh_sha256"],
            "Selected mesh differs from paired source")
    require(sha256(mesh.parent / "setup.json") == targets["mesh_setup_sha256"],
            "Selected setup differs from paired source")
    require(json.loads((mesh.parent / "setup.json").read_text()) == setup,
            "Selected/native setup metadata differs")
    for key, target in [("outer_radius_m", "outer_radius_m"), ("chord_m", "chord_m"),
                       ("angle_degrees", "angle_degrees"), ("reynolds_number", "reynolds_number")]:
        close(setup[key], targets[target], "Source/declaration " + key)
    close(record["reynolds_number"], targets["reynolds_number"], "Recorded Reynolds number")
    close(record["angle_degrees"], targets["angle_degrees"], "Recorded incidence")
    close(record["span_m"], targets["span_m"], "Recorded span")
    settings_path = one_input(paths, "loaded-settings.json")
    settings = json.loads(settings_path.read_text())
    require(settings["methods"] == record["fluent"]["numerical_methods"],
            "Fluent native/recorded methods differ")
    for key in ["system/fvSchemes", "system/fvSolution", "0/U", "0/p"]:
        require(path_key(case / key) in paths, "Native manifest misses " + key)
    return dict(directory=directory, record=record, targets=targets, setup=setup, case=case,
        mesh=mesh, manifest=manifest, paths=paths, fluent_settings=settings,
        radius=targets["outer_radius_m"] / targets["chord_m"],
        schedule=relaxation_schedule(record, targets))


def verify_domain_provenance(baseline, extended):
    read_mesh, reader_path = source_reader()
    source, target = read_mesh(baseline["mesh"]), read_mesh(extended["mesh"])
    proof = verify_inner_grid(source, target)
    path = extended["mesh"].parent / "domain-extension-provenance.json"
    provenance = json.loads(path.read_text())
    require(path_key(provenance["source"]) == path_key(baseline["mesh"])
            and provenance["source_sha256"] == sha256(baseline["mesh"])
            and provenance["source_setup_sha256"] == sha256(baseline["mesh"].parent / "setup.json")
            and provenance["mesh_sha256"] == sha256(extended["mesh"]),
            "Extension provenance does not bind the two selected meshes")
    require(extended["setup"].get("domain_extension") == provenance,
            "Native setup/extension receipt differs")
    for key in ["preserved_source_node_coordinates", "preserved_source_triangle_cell_ids"]:
        require(extended["setup"].get(key) is True, "Missing preservation setup flag: " + key)
    close(provenance["old_radius_m"], baseline["targets"]["outer_radius_m"], "Old radius")
    close(provenance["new_radius_m"], extended["targets"]["outer_radius_m"], "New radius")
    require(provenance["centre_m"] == [.5, 0], "Unsupported outer-polygon centre")
    saved = provenance["preservation_checks"]
    for key, value in proof.items():
        if key != "original_internal_face_cell_pairs_exact":
            require(saved.get(key) == value, "Extension preservation receipt disagrees: " + key)
    for key in ["manifold_oriented_faces", "connected_cell_graph",
                "all_added_centroids_outside_source_polygon", "written_fluent_readback_exact",
                "outer_polygon_uniform_radial_scale"]:
        require(saved.get(key) is True, "Extension geometric check missing/failed: " + key)
    ratio = extended["targets"]["outer_radius_m"] / baseline["targets"]["outer_radius_m"]
    old_ids = np.unique(source["edges"][np.isin(source["zones"], [5, 6])])
    new_ids = np.unique(target["edges"][np.isin(target["zones"], [5, 6])])
    expected = (source["xy"][old_ids] - [.5, 0]) * ratio + [.5, 0]
    actual = target["xy"][new_ids]
    expected = expected[np.lexsort((expected[:, 1], expected[:, 0]))]
    actual = actual[np.lexsort((actual[:, 1], actual[:, 0]))]
    require(actual.shape == expected.shape and np.max(np.abs(actual - expected)) <= 1e-10,
            "New outer polygon is not the declared uniform scale of the original polygon")
    extra = [path, reader_path]
    for name, key in [("executed-mesher.py", "script_sha256"),
                      ("executed-source-reader.py", "source_reader_sha256")]:
        artifact = path.parent / name
        require(sha256(artifact) == provenance[key], "Extension snapshot changed: " + name)
        extra.append(artifact)
    if "annulus_mesh_sha256" in provenance:
        artifact = path.parent / "source-gmsh-annulus.msh"
        require(sha256(artifact) == provenance["annulus_mesh_sha256"], "Gmsh annulus input changed")
        extra.append(artifact)
        proof["annulus_method"] = "Gmsh"
    else:
        mesher = provenance["mesher"]
        require(mesher.get("method") == "deterministic_radial_triangles"
                and mesher.get("gmsh_used") is False, "Unrecorded non-Gmsh annulus method")
        layers = mesher["radial_layers"]
        require(type(layers) is int and layers > 0, "Invalid radial layer count")
        require(proof["added_nodes"] == len(old_ids) * layers
                and proof["added_cells"] == 2 * len(old_ids) * layers,
                "Radial annulus node/cell counts do not match recorded layers")
        name = provenance["radial_topology_file"]
        require(Path(name).name == name, "Invalid radial topology filename")
        artifact = path.parent / name
        require(sha256(artifact) == provenance["radial_topology_sha256"],
                "Radial annulus topology snapshot changed")
        with np.load(artifact, allow_pickle=False) as archive:
            require(np.array_equal(archive["added_nodes_xy"], target["xy"][proof["original_nodes"]:]),
                    "Radial snapshot/appended coordinates differ")
            require(np.array_equal(np.sort(archive["added_triangles"], axis=1),
                                   target["triangles"][proof["original_cells"]:]),
                    "Radial snapshot/appended triangles differ")
            factors = ratio ** (np.arange(1, layers + 1) / layers)
            require(archive["radius_factors"].shape == factors.shape
                    and np.allclose(archive["radius_factors"], factors, rtol=1e-13, atol=1e-14),
                    "Radial snapshot does not follow the declared geometric schedule")
        extra.append(artifact)
        proof.update(annulus_method="deterministic_radial_triangles", radial_layers=layers,
                     radial_topology_snapshot_verified=True)
    return proof, extra


def analyse(baseline_dir, extended_dir):
    baseline, extended = load_pair(baseline_dir), load_pair(extended_dir)
    require(baseline["radius"] == 20 and extended["radius"] == 40,
            "This controlled comparison requires declared R20c and R40c")
    physical = ["reynolds_number", "angle_degrees", "density_kg_m3", "dynamic_viscosity_Pa_s",
                "kinematic_viscosity_m2_s", "velocity_m_s", "chord_m", "span_m",
                "reference_area_m2", "pressure_reference_Pa_gauge", "n_nonorthogonal_correctors",
                "zero_angle_abs_lift_diagnostic_target", "cross_code_cd_difference_learning_target_percent"]
    for key in physical:
        close(baseline["targets"][key], extended["targets"][key], "Matched condition " + key)
    require(baseline["schedule"] == extended["schedule"], "OpenFOAM relaxation schedules differ")
    require(baseline["fluent_settings"]["methods"] == extended["fluent_settings"]["methods"],
            "Fluent numerical methods differ between radii")
    excluded_reference_labels = match_reference_values(
        baseline["fluent_settings"]["references"], extended["fluent_settings"]["references"])
    for name in ["fvSchemes", "fvSolution"]:
        require(semantic_file(baseline["case"] / "system" / name)
                == semantic_file(extended["case"] / "system" / name),
                "OpenFOAM method dictionary differs: " + name)
    require(baseline["record"]["acceptance_thresholds"] == extended["record"]["acceptance_thresholds"],
            "Numerical/symmetry acceptance thresholds differ between radii")
    proof, extra = verify_domain_provenance(baseline, extended)
    rows = [row for label, pair in [("baseline", baseline), ("extended", extended)]
            for row in build_rows(label, pair["radius"], pair["record"])]
    deltas = {}
    for solver in ["Fluent", "OpenFOAM"]:
        first, second = [row for row in rows if row["solver"] == solver]
        require(first["Cd"] != 0, "Zero baseline drag cannot normalize domain difference")
        change = 100 * (second["Cd"] - first["Cd"]) / abs(first["Cd"])
        deltas[solver] = change
        first["own_code_Cd_change_from_R20_percent"] = 0.
        second["own_code_Cd_change_from_R20_percent"] = change
    files = extra + [Path(__file__)]
    for pair in [baseline, extended]:
        files.extend(path for path, _ in pair["paths"].values())
        files.extend([pair["mesh"], pair["mesh"].parent / "setup.json",
                      *(pair["directory"] / name for name in
                        ["comparison.json", "declared-comparison.json", "source-manifest.json"])])
    unique = {path_key(path): Path(path).resolve() for path in files}
    result = dict(schema_version=1, reynolds_number=baseline["targets"]["reynolds_number"],
        angle_degrees=baseline["targets"]["angle_degrees"], radii_over_chord=[20., 40.],
        records=rows, own_code_Cd_change_percent=deltas,
        delta_definition="100*(Cd_R40-Cd_R20)/abs(Cd_R20); signed, within each code",
        inner_grid_preservation=proof, openfoam_relaxation_schedule=baseline["schedule"],
        zero_angle_abs_lift_diagnostic_target=baseline["targets"]["zero_angle_abs_lift_diagnostic_target"],
        native_input_hashes_verified=True, numerical_methods_and_references_matched=True,
        fluent_reference_display_labels_excluded=excluded_reference_labels,
        reference_comparison_policy="Every Fluent reference entry matches exactly except the string zone selector label; both excluded labels are recorded.",
        all_numerical_checks_passed=all(row["numerical_checks_passed"] for row in rows),
        all_cross_code_targets_passed=all(row["cross_code_target_pass"] for row in rows),
        diagnostic_result=any(not row["numerical_checks_passed"] or not row["cross_code_target_pass"]
                              for row in rows),
        mesh_independence_established=False, domain_independence_established=False,
        error_bound_established=False, convergence_order_established=False, physical_validation=False,
        note="Two measured domains with the exact original inner grid retained. Added annulus discretization also changes. Failed checks are retained; no error bound, inferred order or domain-independence claim.")
    manifest = dict(script_sha256=sha256(Path(__file__)),
        inputs=[dict(path=str(path), sha256=sha256(path)) for path in unique.values()])
    return result, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--extended", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    require(out.is_relative_to(Path("E:/CFD-Learning-Lab").resolve()),
            "Outputs must remain under the E: project")
    if out.exists():
        raise FileExistsError("Use a fresh output directory")
    result, manifest = analyse(args.baseline, args.extended)
    out.mkdir(parents=True)
    with (out / "comparison.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    with (out / "comparison.csv").open("x", newline="", encoding="utf-8") as stream:
        rows = [{**row, "failed_checks": ";".join(row["failed_checks"])} for row in result["records"]]
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    shutil.copy2(Path(__file__), out / "executed-analyser.py")
    manifest["outputs"] = {name: sha256(out / name) for name in
                           ["comparison.json", "comparison.csv", "executed-analyser.py"]}
    with (out / "source-manifest.json").open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
