"""Render checked saved OpenFOAM prism-cell fields without running CFD.

The read-only mesh audit is also used by analyse_openfoam_airfoil.py. Cell
triangles are reconstructed from actual faces/owner/neighbour IDs, allowing
renumbering relative to the imported Fluent mesh. Rendering requires accepted
comparison.json and source-manifest.json, with matching native input hashes.
Failed numerical criteria require explicit --diagnostic and prominent labelling;
structural, source, geometry and evidence checks cannot be bypassed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import re
import shutil

import numpy as np


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(data)
    return digest.hexdigest()


def text(path):
    return re.sub(r"/\*.*?\*/|//[^\n]*", "", path.read_text(), flags=re.S)


def foam_list(path):
    raw = text(path)
    header = re.search(r"\bFoamFile\s*\{([^}]*)\}", raw, re.S)
    require(header is not None and re.search(r"\bformat\s+ascii\s*;", header[1]),
            "Only ASCII OpenFOAM files are supported: " + str(path))
    match = re.fullmatch(r"\s*(\d+)\s*\((.*)\)\s*;?\s*", raw[header.end():], re.S)
    require(match is not None, "Unsupported OpenFOAM list: " + str(path))
    return int(match[1]), match[2]


def read_points(path):
    count, body = foam_list(path)
    values = np.fromstring(body.replace("(", " ").replace(")", " "), sep=" ")
    require(values.size == 3 * count and np.isfinite(values).all(), "Invalid native points")
    return values.reshape(-1, 3)


def read_faces(path):
    count, body = foam_list(path)
    sizes = np.fromiter((int(m[1]) for m in re.finditer(r"^\s*(\d+)\s*\(", body, re.M)),
                        dtype=np.int64)
    require(len(sizes) == count and np.all(np.isin(sizes, [3, 4])),
            "Expected ordinary triangle/quad faceList for one-cell prisms")
    numbers = re.sub(r"^\s*\d+\s*\(", "", body, flags=re.M).replace(")", " ")
    vertices = np.fromstring(numbers, sep=" ", dtype=np.int64)
    require(len(vertices) == int(sizes.sum()), "Native face vertex count differs")
    return sizes, vertices


def read_labels(path):
    count, body = foam_list(path)
    values = np.fromstring(body, sep=" ", dtype=np.int64)
    require(len(values) == count, "Native label count differs: " + str(path))
    return values


def read_source_mesh(path):
    raw = path.read_text()
    lines = raw.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith("(10 (1 ")) + 1
    stop = next(i for i in range(start, len(lines)) if lines[i] == "))")
    xy = np.array([list(map(float, line.split())) for line in lines[start:stop]])
    require(xy.ndim == 2 and xy.shape[1] == 2 and np.isfinite(xy).all(),
            "Expected planar Fluent source nodes")
    header = re.search(r"^\(12 \(0 1 ([0-9a-fA-F]+) 0\)\)", raw, re.M)
    require(header is not None, "Unsupported source cell header")
    cells = int(header[1], 16)
    faces, zones, patch_counts = [], [], {}
    for match in re.finditer(
            r"^\(13 \(([0-9a-fA-F]+) ([0-9a-fA-F]+) ([0-9a-fA-F]+) "
            r"([0-9a-fA-F]+) 2\)\(\s*\n(.*?)^\)\)", raw, re.M | re.S):
        zone, first, last, _ = [int(v, 16) for v in match.groups()[:4]]
        rows = [[int(v, 16) for v in line.split()] for line in match[5].splitlines()
                if line.strip()]
        require(len(rows) == last - first + 1 and all(len(row) == 4 for row in rows),
                "Invalid source face block")
        faces.extend(rows)
        zones.extend([zone] * len(rows))
        patch_counts[zone] = len(rows)
    require(set(patch_counts) == {3, 4, 5, 6}, "Unexpected source face groups")
    faces = np.asarray(faces, dtype=np.int64) - 1
    edges, owner, neighbour = faces[:, :2], faces[:, 2], faces[:, 3]
    require(np.all(edges >= 0) and np.all(edges < len(xy))
            and np.all(owner >= 0) and np.all(owner < cells)
            and np.all(neighbour >= -1) and np.all(neighbour < cells), "Invalid source IDs")
    internal = neighbour >= 0
    keys = np.unique(np.concatenate([
        np.repeat(owner, 2) * len(xy) + edges.ravel(),
        np.repeat(neighbour[internal], 2) * len(xy) + edges[internal].ravel()]))
    ids, counts = np.unique(keys // len(xy), return_counts=True)
    require(np.array_equal(ids, np.arange(cells)) and np.all(counts == 3),
            "Source cells do not each have three unique vertices")
    return dict(xy=xy, edges=edges, owner=owner, neighbour=neighbour,
                zones=np.asarray(zones), cells=cells,
                triangles=(keys % len(xy)).reshape(cells, 3),
                patch_counts=dict(airfoil=patch_counts[4], inlet=patch_counts[5],
                                  outlet=patch_counts[6]))


def reconstruct_prisms(points, sizes, vertices, owner, neighbour, cells):
    """Return native-cell-ordered projected triangles and useful topology arrays."""
    require(len(owner) == len(sizes) and len(neighbour) <= len(owner), "Face/owner counts differ")
    require(np.all(vertices >= 0) and np.all(vertices < len(points)), "Invalid face vertex IDs")
    require(np.all(owner >= 0) and np.all(owner < cells)
            and np.all(neighbour >= 0) and np.all(neighbour < cells), "Invalid native cell IDs")
    require(np.all(owner[:len(neighbour)] != neighbour), "A face owns its own neighbour")
    require(np.all(sizes[:len(neighbour)] == 4), "Internal faces must be vertical quads")
    incidence = np.bincount(owner, minlength=cells) + np.bincount(neighbour, minlength=cells)
    cap_incidence = np.bincount(owner[sizes == 3], minlength=cells)
    require(np.all(incidence == 5) and np.all(cap_incidence == 2),
            "Each prism must have two triangle caps and three quad faces")
    stop = int(sizes[:len(neighbour)].sum())
    keys = np.unique(np.concatenate([
        np.repeat(owner, sizes) * len(points) + vertices,
        np.repeat(neighbour, sizes[:len(neighbour)]) * len(points) + vertices[:stop]]))
    cell_ids, counts = np.unique(keys // len(points), return_counts=True)
    require(np.array_equal(cell_ids, np.arange(cells)) and np.all(counts == 6),
            "Each native cell must have six distinct XYZ vertices")
    xyz_vertices = (keys % len(points)).reshape(cells, 6)
    xy, point_to_xy = np.unique(points[:, :2], axis=0, return_inverse=True)
    projected = np.sort(point_to_xy[xyz_vertices], axis=1)
    require(np.all(projected[:, 0] == projected[:, 1])
            and np.all(projected[:, 2] == projected[:, 3])
            and np.all(projected[:, 4] == projected[:, 5])
            and np.all(projected[:, 0] < projected[:, 2])
            and np.all(projected[:, 2] < projected[:, 4]),
            "Each prism must project to exactly three unique XY vertices")
    offsets = np.concatenate([[0], np.cumsum(sizes)])
    quad_faces = np.flatnonzero(sizes == 4)
    quads = vertices[offsets[quad_faces, None] + np.arange(4)]
    require(np.all(np.diff(np.sort(quads, axis=1), axis=1) > 0),
            "A quadrilateral face repeats an XYZ vertex")
    projected_quads = np.sort(point_to_xy[quads], axis=1)
    require(np.all(projected_quads[:, 0] == projected_quads[:, 1])
            and np.all(projected_quads[:, 2] == projected_quads[:, 3])
            and np.all(projected_quads[:, 0] < projected_quads[:, 2]),
            "A quad face is not a vertical projected edge")
    cap_faces = np.flatnonzero(sizes == 3)
    caps = vertices[offsets[cap_faces, None] + np.arange(3)]
    require(np.all(np.diff(np.sort(caps, axis=1), axis=1) > 0),
            "A triangular cap repeats an XYZ vertex")
    require(np.all(np.ptp(points[caps, 2], axis=1) == 0), "A triangle face is not a planar cap")
    planes = np.unique(points[:, 2])
    require(len(planes) == 2, "A prism layer needs two z-planes")
    for plane in planes:
        require(np.all(np.bincount(owner[cap_faces][points[caps[:, 0], 2] == plane],
                                   minlength=cells) == 1),
                "Each prism must have one cap on each z-plane")
    return dict(xy=xy, triangles=projected[:, [0, 2, 4]], quad_faces=quad_faces,
                projected_quad_edges=projected_quads[:, [0, 2]], point_to_xy=point_to_xy)


def boundary_zones(path, sizes, internal_count):
    count, body = foam_list(path)
    matches = list(re.finditer(r"\b([A-Za-z]\w*)\s*\{([^}]*)\}", body, re.S))
    require(len(matches) == count == 4, "Unexpected native boundary blocks")
    markers = np.full(len(sizes), 3, dtype=np.int64)
    assigned = np.zeros(len(sizes), dtype=bool)
    expected = dict(airfoil=(4, "wall"), inlet=(5, "patch"), outlet=(6, "patch"),
                    frontAndBackPlanes=(7, "empty"))
    require({m[1] for m in matches} == set(expected), "Unexpected patch names")
    for match in matches:
        zone, kind = expected[match[1]]
        patch = match[2]
        require(re.search(r"\btype\s+" + kind + r"\s*;", patch), "Unexpected patch type")
        first = int(re.search(r"\bstartFace\s+(\d+)\s*;", patch)[1])
        length = int(re.search(r"\bnFaces\s+(\d+)\s*;", patch)[1])
        require(first >= internal_count and first + length <= len(sizes)
                and not np.any(assigned[first:first + length]), "Invalid/overlapping patch face range")
        assigned[first:first + length] = True
        markers[first:first + length] = zone
    require(not assigned[:internal_count].any() and assigned[internal_count:].all(),
            "Boundary patches do not cover all boundary faces")
    require(np.all(markers[sizes == 3] == 7) and np.all(markers[sizes == 4] != 7),
            "Front/back patch contains a vertical face or misses a cap")
    return markers


def coordinate_bijection(source, candidate, tolerance=1e-10):
    """Unique geometric node matching; decimal serialization cannot change rank."""
    source, candidate = np.asarray(source), np.asarray(candidate)
    require(0 < tolerance <= 1e-10, "Coordinate tolerance must not exceed 1e-10 m")
    require(source.ndim == 2 and source.shape[1] == 2 and source.shape == candidate.shape
            and len(source) > 0, "Coordinate count/dimension differs")
    require(np.isfinite(source).all() and np.isfinite(candidate).all(), "Nonfinite coordinates")
    require(len(np.unique(source, axis=0)) == len(source)
            and len(np.unique(candidate, axis=0)) == len(candidate), "Duplicate coordinates")
    bins = {}
    for index, point in enumerate(source):
        key = tuple(np.floor(point / tolerance).astype(np.int64))
        bins.setdefault(key, []).append(index)
    mapping = np.empty(len(candidate), dtype=np.int64)
    for index, point in enumerate(candidate):
        bx, by = np.floor(point / tolerance).astype(np.int64)
        matches = []
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:
                for other in bins.get((bx + dx, by + dy), []):
                    delta = point - source[other]
                    if float(np.hypot(delta[0], delta[1])) <= tolerance:
                        matches.append(other)
        require(len(matches) == 1, "Missing or ambiguous coordinate match")
        mapping[index] = matches[0]
    require(len(np.unique(mapping)) == len(source), "Coordinate mapping is not bijective")
    delta = candidate - source[mapping]
    norm = np.linalg.norm(delta, axis=1)
    return mapping, dict(method="unique geometric neighbour-bucket bijection",
        source_nodes=len(source), unique_source_matches=len(np.unique(mapping)),
        tolerance_m=tolerance, maximum_component_error_m=float(np.max(np.abs(delta))),
        maximum_euclidean_error_m=float(norm.max()),
        euclidean_error_quantiles_m=dict(p50=float(np.quantile(norm, .5)),
            p90=float(np.quantile(norm, .9)), p99=float(np.quantile(norm, .99))))


def load_case_geometry(case):
    """Bind every prism, vertical face/cell pair and patch to the imported source."""
    declaration_path = case / "declared-comparison.json"
    targets = json.loads(declaration_path.read_text())
    imported = case / "input-airfoil.msh"
    require(sha256(imported) == targets["mesh_source_sha256"], "Imported source hash differs")
    source = read_source_mesh(imported)
    require(source["cells"] == targets["intended_cells"], "Imported cell count differs")
    mesh = case / "0/polyMesh"
    points = read_points(mesh / "points")
    sizes, vertices = read_faces(mesh / "faces")
    owner, neighbour = read_labels(mesh / "owner"), read_labels(mesh / "neighbour")
    require(len(points) == 2 * len(source["xy"]), "Extrusion point count differs")
    span = targets["span_m"]
    require(set(points[:, 2]) == {-span / 2, span / 2}, "Wrong extrusion planes/span")
    result = reconstruct_prisms(points, sizes, vertices, owner, neighbour, source["cells"])
    xy = result["xy"]
    require(len(xy) == len(source["xy"]), "Unexpected projected vertex count")
    xy_to_source, coordinate_metrics = coordinate_bijection(source["xy"], xy, 1e-10)
    maximum_xy = coordinate_metrics["maximum_euclidean_error_m"]
    source_order = np.lexsort((source["xy"][:, 1], source["xy"][:, 0]))
    coordinate_metrics["independent_lexsort_pair_error_m"] = float(
        np.max(np.linalg.norm(xy - source["xy"][source_order], axis=1)))
    # Each original projected vertex must have exactly one point on each z-plane.
    plane = (points[:, 2] > 0).astype(np.int64)
    projected_plane_keys = result["point_to_xy"] * 2 + plane
    require(np.array_equal(np.sort(projected_plane_keys), np.arange(2 * len(xy))),
            "Missing or duplicated extruded vertex/plane")
    source_to_xy = np.empty(len(xy), dtype=np.int64)
    source_to_xy[xy_to_source] = np.arange(len(xy))
    source_triangles = np.sort(source_to_xy[source["triangles"]], axis=1)
    native_triangles = result["triangles"]
    order_source = np.lexsort(source_triangles[:, ::-1].T)
    order_native = np.lexsort(native_triangles[:, ::-1].T)
    require(np.array_equal(native_triangles[order_native], source_triangles[order_source]),
            "Reconstructed native cell triangles differ from imported source cells")
    require(np.all(np.any(np.diff(source_triangles[order_source], axis=0) != 0, axis=1)),
            "Duplicated source triangle")
    native_to_source_cell = np.empty(source["cells"], dtype=np.int64)
    native_to_source_cell[order_native] = order_source
    quad_faces = result["quad_faces"]
    native_edges = result["projected_quad_edges"]
    source_edges = source_to_xy[source["edges"]]
    def edge_keys(edges):
        return np.min(edges, axis=1) * len(xy) + np.max(edges, axis=1)
    native_edge_keys, source_edge_keys = edge_keys(native_edges), edge_keys(source_edges)
    ni, si = np.argsort(native_edge_keys), np.argsort(source_edge_keys)
    require(len(native_edge_keys) == len(source_edge_keys)
            and len(np.unique(native_edge_keys)) == len(native_edge_keys)
            and np.array_equal(native_edge_keys[ni], source_edge_keys[si]),
            "Projected native vertical faces differ from source edges")
    full_neighbour = np.full(len(owner), -1, dtype=np.int64)
    full_neighbour[:len(neighbour)] = neighbour
    paired_neighbour = full_neighbour[quad_faces]
    valid = paired_neighbour >= 0
    paired_neighbour[valid] = native_to_source_cell[paired_neighbour[valid]]
    native_pairs = np.sort(np.column_stack(
        [native_to_source_cell[owner[quad_faces]], paired_neighbour]), axis=1)
    source_pairs = np.sort(np.column_stack([source["owner"], source["neighbour"]]), axis=1)
    require(np.array_equal(native_pairs[ni], source_pairs[si]),
            "Vertical-face owner/neighbour cell pairs differ from source")
    markers = boundary_zones(mesh / "boundary", sizes, len(neighbour))
    require(np.array_equal(markers[quad_faces][ni], source["zones"][si]),
            "Native inlet/outlet/wall face identities differ from source zones")
    result["metrics"] = dict(cells=source["cells"], points=len(points), faces=len(owner),
        internal_faces=len(neighbour), source_mesh_sha256=targets["mesh_source_sha256"],
        projected_xy_max_difference_m=maximum_xy, coordinate_tolerance_m=1e-10,
        coordinate_bijection=coordinate_metrics,
        six_xyz_and_three_xy_vertices_per_cell=True, five_faces_per_prism=True,
        one_cap_on_each_z_plane_per_cell=True,
        imported_cell_triangles_verified=True, vertical_face_cell_pairs_verified=True,
        source_patch_face_identities_verified=True)
    result["inputs"] = [imported, declaration_path,
                        *(mesh / p for p in ["points", "faces", "owner", "neighbour", "boundary"])]
    return result


def internal_field(path, kind, cells):
    raw = text(path)
    require(re.search(r"\bformat\s+ascii\s*;", raw), "Only ASCII reconstructed fields supported")
    expected = "[0 1 -1 0 0 0 0]" if kind == "vector" else "[0 2 -2 0 0 0 0]"
    dimensions = re.search(r"\bdimensions\s+(\[[^]]+\])\s*;", raw)
    require(dimensions and " ".join(dimensions[1].split()) == expected, "Wrong field dimensions")
    nonuniform = re.search(r"\binternalField\s+nonuniform\s+List<" + kind
                          + r">\s+(\d+)\s*\((.*?)\)\s*;", raw, re.S)
    width = 3 if kind == "vector" else 1
    if nonuniform:
        require(int(nonuniform[1]) == cells, "Native internalField cell count differs")
        values = np.fromstring(nonuniform[2].replace("(", " ").replace(")", " "), sep=" ")
        require(len(values) == cells * width, "Native internalField value count differs")
    else:
        uniform = re.search(r"\binternalField\s+uniform\s+([^;]+);", raw)
        require(uniform is not None, "Unsupported reconstructed internalField")
        value = np.fromstring(uniform[1].replace("(", " ").replace(")", " "), sep=" ")
        require(len(value) == width, "Invalid uniform internalField")
        values = np.tile(value, cells)
    require(np.isfinite(values).all(), "Nonfinite saved field")
    return values.reshape(cells, 3) if kind == "vector" else values


def validate_plot_record(comparison, diagnostic=False):
    failed = comparison.get("failed_checks", [])
    require(isinstance(failed, list) and all(isinstance(name, str) for name in failed),
            "Invalid recorded failed-check list")
    if diagnostic:
        results = comparison.get("check_results", {})
        require(comparison["numerical_checks_passed"] is False and failed
                and isinstance(results, dict) and results
                and all(type(v) is bool for v in results.values())
                and set(failed) == {name for name, passed in results.items() if not passed},
                "--diagnostic requires an explicitly retained failed record with named checks")
        return dict(label="DIAGNOSTIC - checks failed", diagnostic=True,
                    numerical_checks_passed=False, failed_checks=failed)
    require(comparison["numerical_checks_passed"] is True and not failed,
            "Failed numerical checks: default plotting refuses; use --diagnostic only for retained diagnostics")
    return dict(label="", diagnostic=False, numerical_checks_passed=True, failed_checks=[])


def verify_accepted_hashes(checks, required, diagnostic=False):
    comparison_path = checks / "comparison.json"
    manifest_path = checks / "source-manifest.json"
    comparison = json.loads(comparison_path.read_text())
    mode = validate_plot_record(comparison, diagnostic)
    manifest = json.loads(manifest_path.read_text())
    if diagnostic:
        require(manifest.get("numerical_checks_passed") is False
                and manifest.get("failed_checks") == mode["failed_checks"],
                "Diagnostic comparison and source-manifest failure labels differ")
    accepted = {os.path.normcase(str(Path(row["path"]).resolve())): row["sha256"]
                for row in manifest["inputs"]}
    for path in required:
        resolved = os.path.normcase(str(path.resolve()))
        require(resolved in accepted, "Accepted manifest does not cover input: " + resolved)
        require(sha256(path) == accepted[resolved], "Input changed after acceptance: " + resolved)
    return comparison, [comparison_path, manifest_path]


def check_colour_ranges(speed, cp, cp_min=-2, cp_max=2):
    require(np.isfinite(cp_min) and np.isfinite(cp_max) and cp_min < cp_max,
            "Cp colour limits must be finite and increasing")
    require(np.isfinite(speed).all() and np.isfinite(cp).all(), "Nonfinite plotted values")
    ranges = dict(speed=[float(speed.min()), float(speed.max())],
                  Cp=[float(cp.min()), float(cp.max())])
    require(speed.min() >= 0 and speed.max() <= 1.2
            and cp.min() >= cp_min and cp.max() <= cp_max,
            "Saved values need wider COMMON Fluent/OpenFOAM colour ranges; do not clip: "
            + json.dumps(ranges))
    return ranges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--checks-dir", type=Path,
                        help="Accepted or explicitly retained comparison/manifest; default: case.parent/analysis")
    parser.add_argument("--diagnostic", action="store_true",
                        help="Render retained failed numerical checks with prominent diagnostic labelling")
    parser.add_argument("--cp-min", type=float, default=-2,
                        help="Common comparison Cp minimum; all saved values must fit")
    parser.add_argument("--cp-max", type=float, default=2,
                        help="Common comparison Cp maximum; all saved values must fit")
    args = parser.parse_args()
    require(np.isfinite(args.cp_min) and np.isfinite(args.cp_max) and args.cp_min < args.cp_max,
            "Cp colour limits must be finite and increasing")
    case, out = args.case.resolve(), args.output_dir.resolve()
    checks = args.checks_dir.resolve() if args.checks_dir else case.parent / "analysis"
    if out.exists():
        raise FileExistsError("Use a fresh figure output directory")
    comparison = json.loads((checks / "comparison.json").read_text())
    mode = validate_plot_record(comparison, args.diagnostic)
    status = json.loads((case / "run-status.json").read_text())
    final = int(comparison["openfoam"]["final_iteration"])
    require(status["stage"] == "finished_execution" and status["final_iteration"] == final,
            "Live/incomplete/different OpenFOAM execution")
    targets = json.loads((case / "declared-comparison.json").read_text())
    require(all(targets[key] == 1 for key in ["density_kg_m3", "velocity_m_s", "chord_m", "span_m",
                                             "reference_area_m2"]),
            "This matched-format plot supports unit density/speed/chord/span/area only")
    require(targets["angle_degrees"] == 0 and targets["pressure_reference_Pa_gauge"] == 0,
            "This caption/crop supports zero incidence and zero gauge-pressure reference")
    require(comparison["source_mesh_sha256"] == targets["mesh_source_sha256"]
            and comparison["reynolds_number"] == targets["reynolds_number"]
            and comparison["span_m"] == targets["span_m"],
            "Accepted comparison describes a different declared case")
    geometry = load_case_geometry(case)
    cells = geometry["metrics"]["cells"]
    require(cells == comparison["cells"], "Accepted comparison cell count differs")
    field_paths = [case / str(final) / p for p in ["U", "p"]]
    required = [*geometry["inputs"], *field_paths, Path(__file__)]
    comparison, acceptance_paths = verify_accepted_hashes(checks, required, args.diagnostic)
    velocity = internal_field(field_paths[0], "vector", cells)
    pressure = internal_field(field_paths[1], "scalar", cells)
    rho, uinf = targets["density_kg_m3"], targets["velocity_m_s"]
    require(np.max(np.abs(velocity[:, 2])) <= 1e-10 * uinf, "Non-negligible out-of-plane velocity")
    speed = np.linalg.norm(velocity, axis=1) / uinf
    # OpenFOAM p is kinematic; convert the gauge-pressure reference to p/rho units.
    cp = (pressure - targets["pressure_reference_Pa_gauge"] / rho) / (.5 * uinf ** 2)
    ranges = check_colour_ranges(speed, cp, args.cp_min, args.cp_max)
    xy, triangles = geometry["xy"], geometry["triangles"]
    vertices = xy[triangles]
    ab, ac = vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]
    cross = ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]
    require(np.all(np.abs(cross) > 0), "Degenerate projected triangle")
    triangles[cross < 0] = triangles[cross < 0][:, [0, 2, 1]]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix", "font.size": 9.5,
        "axes.linewidth": .65, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True, "pdf.fonttype": 42, "svg.fonttype": "path",
        "savefig.bbox": None, "savefig.facecolor": "white"})
    triangulation = mtri.Triangulation(xy[:, 0], xy[:, 1], triangles)
    fig, axes = plt.subplots(2, 1, figsize=(7.05, 4.6), layout="constrained")
    if mode["diagnostic"]:
        fig.suptitle(mode["label"], fontsize=12, fontweight="bold", color="#A52020")
    for ax, values, limits, cmap, label, panel in zip(
            axes, [speed, cp], [(0, 1.2), (args.cp_min, args.cp_max)], ["viridis", "RdBu_r"],
            [r"$|\mathbf{u}|/U_\infty$", r"$C_p$"], ["(a)", "(b)"]):
        art = ax.tripcolor(triangulation, facecolors=values, shading="flat",
                           vmin=limits[0], vmax=limits[1], cmap=cmap, rasterized=True)
        ax.set(xlim=(-.3, 2.5), ylim=(-.55, .55), xlabel=r"$x/c$", ylabel=r"$y/c$", aspect="equal")
        ax.text(.01, 1.025, panel, transform=ax.transAxes, ha="left", va="bottom")
        fig.colorbar(art, ax=ax, pad=.02, fraction=.035).set_label(label)
    # Freeze the drawn fixed canvas. Avoid bbox_inches='tight' raster clipping.
    fig.canvas.draw()
    fig.set_layout_engine("none")
    out.mkdir(parents=True)
    base = out / "openfoam-airfoil-computed-fields"
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(base.with_suffix("." + suffix), dpi=600)
    plt.close(fig)
    shutil.copy2(Path(__file__), out / "executed-plotter.py")
    caption = (
        f"Saved OpenFOAM airfoil fields at Re={targets['reynolds_number']:g}, zero incidence, "
        f"{cells:,} triangular-prism cells, unit span and solver iteration {final}. "
        "(a) Speed normalized by U_inf=1 m/s. (b) Cp=(p_kinematic-p_ref/rho)/(U_inf^2/2), "
        "with rho=1 kg/m3 and p_ref=0 Pa gauge. The selected comparison colour ranges are "
        f"0..1.2 for normalized speed and {args.cp_min:g}..{args.cp_max:g} for Cp; "
        "every saved value lies within them. Use identical limits in the companion Fluent "
        "figure for a common-scale visual comparison. "
        "Native faces and owner/neighbour IDs determine each projected triangle and its saved "
        "cell-centred field values; source cell order is not assumed. No interpolation, smoothing "
        "or remeshing is applied; the solid interior remains blank. Only the near-body/wake crop "
        "is shown. Field layers are rasterized at 600 dpi within PDF/SVG vector axes, on a fixed "
        "7.05 x 4.6 inch canvas. OpenFOAM computed the fields; Codex reconstructed the native "
        "connectivity and prepared the figure. Solver iterations are not physical time. "
        "Numerical checks do not establish mesh/domain independence or physical validation.\n")
    if mode["diagnostic"]:
        caption = (mode["label"] + ". Failed criteria: " + ", ".join(mode["failed_checks"])
                   + ". These retained fields do not pass the declared acceptance checks. " + caption)
    (out / "CAPTION.md").write_text(caption, encoding="utf-8")
    paths = list(dict.fromkeys([*required, *acceptance_paths, case / "run-status.json"]))
    manifest = dict(script_sha256=sha256(Path(__file__)), measured_ranges=ranges,
                    diagnostic=mode["diagnostic"],
                    numerical_checks_passed=mode["numerical_checks_passed"],
                    failed_checks=mode["failed_checks"],
                    colour_limits=dict(speed=[0, 1.2], Cp=[args.cp_min, args.cp_max]),
                    software=dict(matplotlib=matplotlib.__version__, numpy=np.__version__),
                    canvas_inches=[7.05, 4.6], raster_dpi=600,
                    topology_checks=geometry["metrics"],
                    inputs=[dict(path=str(p.resolve()), sha256=sha256(p)) for p in paths])
    (out / "source-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(output=str(out), cells=cells, ranges=ranges, topology_verified=True), indent=2))


if __name__ == "__main__":
    main()
