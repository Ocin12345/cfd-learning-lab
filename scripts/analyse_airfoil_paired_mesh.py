"""Aggregate hash-bound medium/fine paired airfoil evidence; never run CFD.

The original polygon is checked directly on both emitted Fluent meshes. The
nonuniform grids need not retain node/cell IDs. Native pressure/viscous records
are independently replayed by the frozen force-accounting analyser.
"""
from pathlib import Path
import argparse
import csv
import importlib.util
import json
import math
import shutil

import numpy as np


def module_at(name, filename):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


domain, DOMAIN_PATH = module_at("paired_mesh_domain_helpers", "analyse_airfoil_domain.py")
force, FORCE_PATH = module_at("paired_mesh_force_helpers", "analyse_airfoil_force_components.py")
require, sha256, close, path_key = domain.require, domain.sha256, domain.close, domain.path_key
ROOT = Path("E:/CFD-Learning-Lab").resolve()
METRICS_DIR = Path(__file__).resolve().parents[1] / "results/airfoil-sensitivity-v002"
TOLERANCE_M = 1e-10


def validate_component_record(pair, component):
    """Reject laundering flags, endpoints or physical parameters in derived terms."""
    main = pair["record"]
    for key in ["check_results", "failed_checks", "numerical_checks_passed",
                "numerical_convergence_passed", "symmetry_diagnostic_passed",
                "cross_code_target_pass", "cross_code_Cd_difference_percent",
                "cross_code_learning_target_percent"]:
        require(component[key] == main[key], "Component/main " + key + " differs")
    require(component["main_acceptance_thresholds"] == main["acceptance_thresholds"],
            "Component/main acceptance thresholds differ")
    for key, expected in [("Re", main["reynolds_number"]), ("angle_degrees", main["angle_degrees"]),
                          ("cells", main["cells"]), ("outer_radius_m", pair["targets"]["outer_radius_m"]),
                          ("span_m", main["span_m"])]:
        close(component["parameters"][key], expected, "Component parameter " + key)
    require(component["evidence_consistency_checks_passed"] is True
            and component["native_HDF_evidence_hash_bound_to_main_checks"] is True,
            "Component evidence proof is missing")
    require(component["flow_iterations_executed"] == 0
            and component["numerical_cause_identified"] is False
            and component["physical_validation"] is False,
            "Unsupported component inference or new CFD")
    for code in ["fluent", "openfoam"]:
        for coefficient in ["Cd", "Cl"]:
            values = component[code]
            require(all(math.isfinite(values[term][coefficient])
                        for term in ["pressure", "viscous", "total"]), "Nonfinite component force")
            total, endpoint = values["total"][coefficient], main[code][coefficient]
            tolerance = 0 if code == "fluent" else force.COEFFICIENT_EQUIVALENCE_TOLERANCE
            require(abs(total - endpoint) <= tolerance, "Component total differs from saved " + code)
            require(force.report_close(values["pressure"][coefficient]
                                       + values["viscous"][coefficient], total),
                    "Component total closure exceeds native report rounding")


def load_components(pair, directory):
    """Recheck all source hashes and reproduce the component record read-only."""
    directory = Path(directory).resolve()
    saved = json.loads((directory / "comparison.json").read_text())
    manifest, paths = domain.verify_manifest(directory)
    for name in ["comparison.json", "source-manifest.json"]:
        expected = pair["directory"] / name
        require(path_key(expected) in paths and paths[path_key(expected)][1] == sha256(expected),
                "Component manifest is not bound to the selected paired " + name)
    raw = domain.one_input(paths, "force-components.json")
    settings = domain.one_input(pair["paths"], "loaded-settings.json")
    computed, replay_inputs = force.analyse(settings.parent, raw, pair["directory"])
    require(saved == computed, "Component record differs from independently replayed native evidence")
    validate_component_record(pair, saved)
    require(manifest["evidence_consistency_checks_passed"] is True
            and manifest["failed_checks"] == saved["failed_checks"]
            and manifest["numerical_checks_passed"] == saved["numerical_checks_passed"],
            "Component manifest flags differ")
    extra = [directory / "comparison.json", directory / "source-manifest.json", FORCE_PATH]
    extra += [p for p, _ in paths.values()] + [Path(row["path"]) for row in replay_inputs]
    return saved, extra


def original_polygons(path):
    setup_path = path.parent / "setup.json"
    setup = json.loads(setup_path.read_text())
    require(setup["chord_m"] == 1 and setup["angle_degrees"] == 0,
            "Original source must be an unrotated unit-chord polygon")
    count = int(setup["circumferential_cells"])
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("(10 (1 ")) + 1
    stop = next(i for i in range(start, len(lines)) if lines[i] == "))")
    xy = np.asarray([list(map(float, line.split())) for line in lines[start:stop]])
    require(xy.ndim == 2 and xy.shape[1] == 2 and len(xy) >= 2 * count
            and np.isfinite(xy).all(), "Invalid original polygon coordinates")
    return xy[:count], xy[-count:]


def check_polygon(edges, reference, xy):
    """Project every emitted boundary edge onto one original straight segment."""
    vector = np.roll(reference, -1, axis=0) - reference
    squared = np.sum(vector * vector, axis=1)
    require(np.all(squared > 0), "Zero original polygon segment")
    lengths = np.sqrt(squared)
    coverage = [[] for _ in reference]
    maximum_distance = 0.
    for edge in edges:
        p, q = xy[edge]
        wp, wq = p - reference, q - reference
        dp = np.abs(vector[:, 0] * wp[:, 1] - vector[:, 1] * wp[:, 0]) / lengths
        dq = np.abs(vector[:, 0] * wq[:, 1] - vector[:, 1] * wq[:, 0]) / lengths
        tp, tq = np.sum(wp * vector, axis=1) / squared, np.sum(wq * vector, axis=1) / squared
        epsilon = TOLERANCE_M / lengths
        valid = ((dp <= TOLERANCE_M) & (dq <= TOLERANCE_M)
                 & (tp >= -epsilon) & (tp <= 1 + epsilon)
                 & (tq >= -epsilon) & (tq <= 1 + epsilon))
        require(np.count_nonzero(valid) == 1, "Emitted boundary leaves/crosses original polygon")
        j = int(np.flatnonzero(valid)[0])
        coverage[j].append(tuple(sorted((float(tp[j]), float(tq[j])))))
        maximum_distance = max(maximum_distance, float(dp[j]), float(dq[j]))
    maximum_coverage_error = 0.
    for j, intervals in enumerate(coverage):
        require(intervals, "Missing original polygon segment")
        end = 0.
        for low, high in sorted(intervals):
            maximum_coverage_error = max(maximum_coverage_error, abs(low - end) * lengths[j])
            end = high
        maximum_coverage_error = max(maximum_coverage_error, abs(end - 1) * lengths[j])
    nearest = np.linalg.norm(reference[:, None, :] - xy[np.unique(edges)][None, :, :], axis=2).min(axis=1)
    require(maximum_coverage_error <= TOLERANCE_M and nearest.max() <= TOLERANCE_M,
            "Original polygon vertex/segment coverage changed")
    return dict(passed=True, original_segments=len(reference),
                maximum_edge_to_original_segment_distance_m=maximum_distance,
                maximum_segment_gap_overlap_or_endpoint_error_m=maximum_coverage_error,
                maximum_original_vertex_distance_m=float(nearest.max()))


def mesh_metric_proof(metrics_path, mesh_path):
    record = json.loads(metrics_path.read_text())
    metrics_script = Path(__file__).with_name("airfoil_mesh_metrics.py")
    require(sha256(metrics_script) == record["script_sha256"], "Changed mesh metrics script")
    require(record["checks_passed"] is True and record["patch_assignment_passed"] is True
            and record["topology_boundary_coverage_passed"] is True, "Unaccepted mesh metrics")
    require(record["definitions"]["preservation_tolerance_m"] == TOLERANCE_M,
            "Changed geometry tolerance")
    bound = []
    for source in record["inputs"].values():
        path = Path(source["path"]).resolve()
        require(path.is_file() and sha256(path) == source["sha256"], "Changed metrics source " + str(path))
        bound.append(path)
    source = Path(record["inputs"]["original_source_mesh"]["path"])
    require(path_key(record["inputs"]["original_source_setup"]["path"])
            == path_key(source.parent / "setup.json"), "Original metrics setup is not adjacent to source")
    wall_reference, outer_reference = original_polygons(source)
    reader, reader_path = domain.source_reader()
    mesh = reader(mesh_path)
    require(mesh["cells"] == record["cells"] and len(mesh["xy"]) == record["nodes"],
            "Metrics/emitted mesh cell or node count differs")
    require(mesh["patch_counts"] == record["boundary_face_counts"], "Metrics/emitted patches differ")
    xy, edges, zones = mesh["xy"], mesh["edges"], mesh["zones"]
    wall = edges[zones == 4]
    require(len(wall) == record["wall_face_count"], "Metrics/emitted wall count differs")
    wall_proof = check_polygon(wall, wall_reference, xy)
    outer_proof = check_polygon(edges[np.isin(zones, [5, 6])], outer_reference, xy)
    for name, sign in [("inlet", -1), ("outlet", 1)]:
        patch = edges[zones == {"inlet": 5, "outlet": 6}[name]]
        require(np.all(sign * (xy[patch].mean(axis=1)[:, 0] - .5) > 0), "Emitted inlet/outlet split changed")
    vertices = xy[mesh["triangles"]]
    areas = .5 * np.abs((vertices[:, 1, 0]-vertices[:, 0, 0])*(vertices[:, 2, 1]-vertices[:, 0, 1])
                        -(vertices[:, 1, 1]-vertices[:, 0, 1])*(vertices[:, 2, 0]-vertices[:, 0, 0]))
    require(np.all(areas > 0), "Nonpositive emitted triangle area")
    lengths = np.linalg.norm(xy[wall[:, 1]] - xy[wall[:, 0]], axis=1)
    heights = 2 * areas[mesh["owner"][zones == 4]] / lengths
    centroid = vertices.mean(axis=1)
    wake = (centroid[:, 0] >= 1) & (centroid[:, 0] <= 3) & (np.abs(centroid[:, 1]) < .25)
    median_length, median_height = float(np.median(lengths)), float(np.median(heights))
    median_wake = float(np.median(np.sqrt(areas[wake])))
    close(median_length, record["wall_regions"]["all"]["edge_length_m"]["median"], "Emitted wall length")
    close(median_height, record["wall_regions"]["all"]["adjacent_triangle_normal_height_m"]["median"], "Emitted wall height")
    close(median_wake, record["wake"]["triangle_sqrt_area_m"]["median"], "Emitted wake spacing")
    close(float(areas.sum()), record["fluid_area_m2"], "Emitted fluid area")
    return dict(cells=mesh["cells"], nodes=len(xy), wall_faces=len(wall),
                median_wall_edge_m=median_length, median_wall_normal_height_m=median_height,
                median_wake_sqrt_area_m=median_wake, wall=wall_proof, outer=outer_proof,
                original_source_mesh_sha256=record["inputs"]["original_source_mesh"]["sha256"],
                original_source_setup_sha256=record["inputs"]["original_source_setup"]["sha256"]), [
                    metrics_path, mesh_path, reader_path, metrics_script, *bound]


def verify_geometry_records(baseline_metrics, fine_metrics, baseline_mesh, fine_mesh):
    baseline, left = mesh_metric_proof(Path(baseline_metrics), Path(baseline_mesh))
    fine, right = mesh_metric_proof(Path(fine_metrics), Path(fine_mesh))
    for key in ["original_source_mesh_sha256", "original_source_setup_sha256"]:
        require(baseline[key] == fine[key], "Meshes do not share the same original polygon")
    require(fine["cells"] > baseline["cells"] and fine["nodes"] > baseline["nodes"]
            and fine["wall_faces"] > baseline["wall_faces"], "Fine grid is not refined")
    for key in ["median_wall_edge_m", "median_wall_normal_height_m", "median_wake_sqrt_area_m"]:
        require(fine[key] < baseline[key], "Fine regional spacing is not smaller: " + key)
    proof = dict(same_original_wall_and_outer_polygon=True,
                 exact_cell_or_node_preservation_claimed=False, uniform_refinement=False,
                 geometry_tolerance_m=TOLERANCE_M,
                 gmsh_topology_identity_claimed=False,
                 emitted_mesh_polygon_and_spacing_recomputed=True,
                 note="Direct checks use emitted Fluent triangles and wall/outer edges; hashed Gmsh metrics are corroborating records, not an assumed Gmsh-to-Fluent topology identity.")
    for label, record in [("baseline", baseline), ("fine", fine)]:
        for key, value in record.items():
            proof[label + "_" + key] = value
    return proof, left + right


def match_controls(baseline, fine):
    require(baseline["radius"] == fine["radius"] == 20,
            "This controlled comparison requires R20c in both meshes")
    physical = ["reynolds_number", "angle_degrees", "density_kg_m3", "dynamic_viscosity_Pa_s",
        "kinematic_viscosity_m2_s", "velocity_m_s", "chord_m", "span_m", "reference_area_m2",
        "pressure_reference_Pa_gauge", "outer_radius_m", "n_nonorthogonal_correctors",
        "initial_p_U_residual_learning_target", "cd_relative_range_target", "cl_absolute_range_target",
        "relative_boundary_mass_imbalance_target", "zero_angle_abs_lift_diagnostic_target",
        "cross_code_cd_difference_learning_target_percent", "initial_iteration_cap",
        "final_force_window", "mpi_processes"]
    for key in physical:
        close(baseline["targets"][key], fine["targets"][key], "Matched condition " + key)
    require(baseline["targets"]["reynolds_number"] == 100
            and baseline["targets"]["angle_degrees"] == 0, "Expected zero-incidence Re100")
    expected_schedule = [dict(start_iteration=0, U_equation_relaxation=.7),
                         dict(start_iteration=2000, U_equation_relaxation=.9)]
    require(baseline["schedule"] == fine["schedule"] == expected_schedule,
            "OpenFOAM relaxation schedules differ from the controlled schedule")
    require(baseline["fluent_settings"]["methods"] == fine["fluent_settings"]["methods"],
            "Fluent numerical methods differ")
    for key in ["general", "viscous", "material", "inlet", "outlet"]:
        require(baseline["fluent_settings"][key] == fine["fluent_settings"][key],
                "Fluent physical settings differ: " + key)
    excluded = domain.match_reference_values(baseline["fluent_settings"]["references"],
                                             fine["fluent_settings"]["references"])
    for name in ["system/fvSchemes", "system/fvSolution", "0/U", "0/p",
                 "constant/physicalProperties", "constant/momentumTransport"]:
        require(domain.semantic_file(baseline["case"] / name) == domain.semantic_file(fine["case"] / name),
                "OpenFOAM physical/method dictionary differs: " + name)
    require(baseline["record"]["acceptance_thresholds"] == fine["record"]["acceptance_thresholds"],
            "Numerical/symmetry acceptance thresholds differ")
    return excluded


def verify_fine_initialization(fluent_run, expected_geometry):
    """Enforce the saved-field continuation described by this study and caption."""
    fluent_run = Path(fluent_run).resolve()
    names = ["initial-condition-provenance.json", "reconvergence-controls.json",
             "status.json", "summary.json", "geometry-setup.json", "code-provenance.json"]
    saved = {name: json.loads((fluent_run / name).read_text()) for name in names}
    initial, controls = saved[names[0]], saved[names[1]]
    require(initial["kind"] == "Saved-case continuation"
            and initial["initialization"] == "Retain saved fields; no hybrid initialization"
            and initial["snapshot_before_solver_launch"] is True,
            "Fine Fluent must have bound Saved-case continuation provenance")
    require(controls["retained_saved_fields"] is True
            and controls["hybrid_initialization_executed"] is False
            and controls["initial_iteration_cap"] == initial["initial_iteration_cap"],
            "Fine Fluent saved-field control/provenance differs")
    require(controls["residual_equations"] == {
        name: dict(absolute_criteria=1e-8, check_convergence=True)
        for name in ["continuity", "x-velocity", "y-velocity"]},
        "Fine Fluent native 1e-8 reconvergence controls differ")
    summary, status = saved["summary.json"], saved["status.json"]
    require(status["stage"] == "completed"
            and {k: v for k, v in status.items() if k not in ["stage", "time"]} == summary
            and summary["kind"] == "Saved-case continuation; no physical validation"
            and summary["native_case_saved"] is True
            and summary["convergence_checks_restored_before_save"] is True,
            "Fine Fluent saved-field completion evidence differs")
    require(saved["geometry-setup.json"] == expected_geometry,
            "Fine Fluent continuation/source geometry differs")
    source_names = {"airfoil-Re100-final.cas.h5", "airfoil-Re100-final.dat.h5", "geometry-setup.json"}
    inputs = initial["inputs"]
    require(len(inputs) == 3 and {Path(row["path"]).name for row in inputs} == source_names,
            "Saved-field provenance must bind its source case/data/geometry")
    sources = [Path(row["path"]).resolve() for row in inputs]
    require(len({str(path.parent) for path in sources}) == 1,
            "Saved-field source case/data/geometry are not adjacent")
    for row, path in zip(inputs, sources):
        require(path.is_relative_to(ROOT) and path.is_file() and sha256(path) == row["sha256"],
                "Changed saved-field source input: " + str(path))
        if path.name == "geometry-setup.json":
            require(json.loads(path.read_text()) == expected_geometry,
                    "Initial saved-field source is not the selected fine geometry")
    code = saved["code-provenance.json"]
    executed = fluent_run / "executed-script.py"
    require(code["snapshot_before_solver_launch"] is True
            and sha256(executed) == code["source_sha256"], "Fine Fluent executed source changed")
    return dict(kind=initial["kind"], saved_source_inputs_verified=True,
        hybrid_initialization_executed=False, native_reconvergence_target=1e-8,
        initial_iteration_cap=initial["initial_iteration_cap"],
        solver_source_sha256=code["source_sha256"]), [fluent_run / name for name in names] + sources + [executed]


def assemble(baseline, fine, baseline_components, fine_components, geometry, excluded):
    rows, component_gaps = [], {}
    for label, pair, component in [("medium", baseline, baseline_components), ("fine", fine, fine_components)]:
        validate_component_record(pair, component)
        own = domain.build_rows(label, 20, pair["record"])
        for row, code in zip(own, ["fluent", "openfoam"]):
            for term in ["pressure", "viscous"]:
                for coefficient in ["Cd", "Cl"]:
                    row[term + "_" + coefficient] = component[code][term][coefficient]
        rows += own
        gap = dict(Cd_difference_percent=pair["record"]["cross_code_Cd_difference_percent"],
                   target_percent=pair["record"]["cross_code_learning_target_percent"],
                   target_pass=pair["record"]["cross_code_target_pass"])
        for term in ["pressure", "viscous", "total"]:
            for coefficient in ["Cd", "Cl"]:
                gap[term + "_" + coefficient + "_Fluent_minus_OpenFOAM"] = (
                    component["fluent"][term][coefficient] - component["openfoam"][term][coefficient])
        component_gaps[label] = gap
    changes = {}
    for solver in ["Fluent", "OpenFOAM"]:
        first, last = [next(row for row in rows if row["case"] == label and row["solver"] == solver)
                       for label in ["medium", "fine"]]
        require(first["Cd"] != 0, "Zero medium drag cannot normalize sensitivity")
        values = dict(Cd_change_absolute=last["Cd"] - first["Cd"],
                      Cd_change_percent=100 * (last["Cd"] - first["Cd"]) / abs(first["Cd"]),
                      Cl_change_absolute=last["Cl"] - first["Cl"])
        for term in ["pressure", "viscous"]:
            for coefficient in ["Cd", "Cl"]:
                key = term + "_" + coefficient
                values[key + "_change_absolute"] = last[key] - first[key]
        changes[solver] = values
    return dict(kind="Paired nonuniform airfoil mesh sensitivity; no new CFD in this analysis",
        reynolds_number=baseline["targets"]["reynolds_number"], angle_degrees=0,
        outer_radius_over_chord=20, records=rows, own_code_changes=changes,
        cross_code_gaps=component_gaps, geometry_refinement_proof=geometry,
        delta_definition="Fine minus medium; Cd percent =100*(Cd_fine-Cd_medium)/abs(Cd_medium); Cl and component changes are absolute coefficient differences.",
        cross_code_gap_definition="100*abs(Cd_OpenFOAM-Cd_Fluent)/abs(Cd_Fluent), from original paired saved endpoints; component replay totals do not redefine acceptance.",
        openfoam_relaxation_schedule=baseline["schedule"],
        zero_angle_abs_lift_diagnostic_target=baseline["targets"]["zero_angle_abs_lift_diagnostic_target"],
        acceptance_thresholds=baseline["record"]["acceptance_thresholds"],
        fluent_reference_display_labels_excluded=excluded,
        native_input_hashes_verified=True, component_native_evidence_replayed=True,
        numerical_methods_and_references_matched=True,
        all_numerical_checks_passed=all(row["numerical_checks_passed"] for row in rows),
        all_cross_code_targets_passed=all(row["cross_code_target_pass"] for row in rows),
        diagnostic_result=any(not row["numerical_checks_passed"] or not row["cross_code_target_pass"] for row in rows),
        uniform_refinement=False, GCI_established=False, error_bound_established=False,
        convergence_order_established=False, mesh_independence_established=False,
        domain_independence_established=False, physical_validation=False, numerical_cause_identified=False,
        initialization_note="New fine Fluent session continues the preserved converged fine fields; OpenFOAM starts cold. The continuation is not an independent cold Fluent calculation.",
        note="Measured sensitivity within one nonuniform mesh family. Inspect component cancellation and failed zero-angle diagnostics. This does not establish the cause of a solver difference or physical accuracy.")


def analyse(baseline_dir, fine_dir, baseline_components_dir, fine_components_dir,
            baseline_metrics=METRICS_DIR / "metrics-medium.json",
            fine_metrics=METRICS_DIR / "metrics-fine055.json"):
    baseline, fine = domain.load_pair(baseline_dir), domain.load_pair(fine_dir)
    excluded = match_controls(baseline, fine)
    geometry, extra = verify_geometry_records(baseline_metrics, fine_metrics,
                                               baseline["mesh"], fine["mesh"])
    for label, pair in [("baseline", baseline), ("fine", fine)]:
        require(pair["record"]["cells"] == geometry[label + "_cells"],
                "Paired/geometry cell count differs")
    bc, left = load_components(baseline, baseline_components_dir)
    fc, right = load_components(fine, fine_components_dir)
    initialization, initial_paths = verify_fine_initialization(
        domain.one_input(fine["paths"], "loaded-settings.json").parent, fine["setup"])
    result = assemble(baseline, fine, bc, fc, geometry, excluded)
    result["fine_fluent_initialization_proof"] = initialization
    paths = [Path(__file__), DOMAIN_PATH, FORCE_PATH, *extra, *left, *right, *initial_paths]
    for pair in [baseline, fine]:
        paths += [pair["directory"] / name for name in ["comparison.json", "source-manifest.json", "declared-comparison.json"]]
        paths += [p for p, _ in pair["paths"].values()]
        paths += [pair["mesh"], pair["mesh"].parent / "setup.json"]
    unique = {path_key(path): Path(path).resolve() for path in paths}
    manifest = dict(script_sha256=sha256(Path(__file__)),
                    inputs=[dict(path=str(path), sha256=sha256(path)) for path in unique.values()])
    return result, manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["baseline", "fine", "baseline-components", "fine-components", "out"]:
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--baseline-metrics", type=Path, default=METRICS_DIR / "metrics-medium.json")
    parser.add_argument("--fine-metrics", type=Path, default=METRICS_DIR / "metrics-fine055.json")
    args = parser.parse_args()
    out = args.out.resolve()
    require(out.is_relative_to(ROOT / "repo/cfd-learning-lab/results") or out.is_relative_to(ROOT / "local"),
            "Use a fresh project results/local output directory")
    if out.exists():
        raise FileExistsError("Use a fresh output directory")
    result, manifest = analyse(args.baseline, args.fine, args.baseline_components,
                               args.fine_components, args.baseline_metrics, args.fine_metrics)
    out.mkdir(parents=True)
    (out / "comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    rows = [{**row, "failed_checks": ";".join(row["failed_checks"])} for row in result["records"]]
    with (out / "comparison.csv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    shutil.copy2(Path(__file__), out / "executed-analyser.py")
    manifest["outputs"] = {name: sha256(out / name) for name in ["comparison.json", "comparison.csv", "executed-analyser.py"]}
    (out / "source-manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
