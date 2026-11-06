"""Publish fixed-canvas paired airfoil mesh-sensitivity figures from checked data."""
from pathlib import Path
import argparse
import importlib.util
import json
import math
import os
import shutil

for _name in ["TEMP", "TMP", "TMPDIR"]:
    os.environ[_name] = "E:/CFD-Learning-Lab/temp"
os.environ["MPLCONFIGDIR"] = "E:/CFD-Learning-Lab/temp/matplotlib"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text
from matplotlib.ticker import ScalarFormatter
import numpy as np


def helpers():
    path = Path(__file__).with_name("analyse_airfoil_paired_mesh.py")
    spec = importlib.util.spec_from_file_location("paired_mesh_plot_analysis", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def validate(data):
    rows = data["records"]
    if len(rows) != 4 or {(r["case"], r["solver"]) for r in rows} != {
            (case, solver) for case in ["medium", "fine"] for solver in ["Fluent", "OpenFOAM"]}:
        raise ValueError("Expected four unique medium/fine paired results")
    if not all(math.isfinite(row[key]) for row in rows
               for key in ["Cd", "Cl", "pressure_Cd", "viscous_Cd", "pressure_Cl", "viscous_Cl"]):
        raise ValueError("Nonfinite plotted coefficient")
    for key in ["native_input_hashes_verified", "component_native_evidence_replayed",
                "numerical_methods_and_references_matched"]:
        if data[key] is not True:
            raise ValueError("Missing checked aggregation proof: " + key)
    if not data["geometry_refinement_proof"].get("same_original_wall_and_outer_polygon"):
        raise ValueError("Missing shared original polygon proof")
    initial = data.get("fine_fluent_initialization_proof", {})
    if initial.get("saved_source_inputs_verified") is not True or initial.get("hybrid_initialization_executed") is not False:
        raise ValueError("Missing fine saved-field initialization proof")
    for key in ["uniform_refinement", "GCI_established", "error_bound_established",
                "convergence_order_established", "mesh_independence_established",
                "domain_independence_established", "physical_validation", "numerical_cause_identified"]:
        if data[key] is not False:
            raise ValueError("Unsupported inference: " + key)
    for row in rows:
        if (row["numerical_checks_passed"] != (not row["failed_checks"])
                or row["numerical_checks_passed"] != (row["numerical_convergence_passed"]
                                                     and row["symmetry_diagnostic_passed"])):
            raise ValueError("Plotted row acceptance flags disagree")
    numerical = all(row["numerical_checks_passed"] for row in rows)
    cross = all(row["cross_code_target_pass"] for row in rows)
    if (data["all_numerical_checks_passed"] != numerical
            or data["all_cross_code_targets_passed"] != cross
            or data["diagnostic_result"] != (not numerical or not cross)):
        raise ValueError("Aggregate diagnostic flags disagree with plotted rows")
    if data["zero_angle_abs_lift_diagnostic_target"] != .001:
        raise ValueError("Changed zero-angle diagnostic threshold")


def make_figure(data):
    rows = data["records"]
    if len(rows) != 4 or not all(math.isfinite(r[key]) for r in rows
                               for key in ["Cd", "Cl", "pressure_Cd", "viscous_Cd"]):
        raise ValueError("Expected four finite paired result rows")
    plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": 9.5, "axes.labelsize": 10, "axes.linewidth": .65,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
        "ytick.right": True, "lines.linewidth": 1.05, "legend.fontsize": 8.5,
        "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "path",
        "savefig.bbox": None, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(2, 2, figsize=(7.05, 5.05), layout="constrained")
    keys = ["Cd", "Cl", "pressure_Cd", "viscous_Cd"]
    labels = [r"Drag coefficient, $C_D$", r"Lift coefficient, $C_L$",
              r"Pressure drag coefficient, $C_{D,p}$", r"Viscous drag coefficient, $C_{D,v}$"]
    for solver, color, marker, style in [
            ("Fluent", "#0072B2", "o", "-"), ("OpenFOAM", "#D55E00", "s", "--")]:
        own = [next(row for row in rows if row["case"] == case and row["solver"] == solver)
               for case in ["medium", "fine"]]
        for ax, key in zip(axes.flat, keys):
            ax.plot([0, 1], [row[key] for row in own], label=solver, color=color,
                    marker=marker, ls=style, ms=4.4, mfc="white", mew=.85)
    for ax, label, panel in zip(axes.flat, labels, ["(a)", "(b)", "(c)", "(d)"]):
        ax.set_xlim(-.15, 1.15)
        ax.set_xticks([0, 1], ["Medium", "Fine"])
        ax.set_ylabel(label)
        ax.text(.01, 1.035, panel, transform=ax.transAxes, ha="left", va="bottom")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_scientific(False); formatter.set_useOffset(False)
        ax.yaxis.set_major_formatter(formatter)
        ax.margins(y=.19)
    axes[0, 0].legend(loc="best")
    target = data["zero_angle_abs_lift_diagnostic_target"]
    if target != .001:
        raise ValueError("Expected the unchanged +/-0.001 lift diagnostic")
    for value in [-target, target]:
        axes[0, 1].axhline(value, color=".5", ls=":", lw=.7, zorder=0)
    axes[0, 1].axhline(0, color=".75", lw=.55, zorder=0)
    low, high = min(-target, *(r["Cl"] for r in rows)), max(target, *(r["Cl"] for r in rows))
    pad = .15 * (high - low)
    axes[0, 1].set_ylim(low-pad, high+pad)
    counts = [next(r["cells"] for r in rows if r["case"] == case) for case in ["medium", "fine"]]
    fig.supxlabel(f"Nonuniform mesh: {counts[0]:,} and {counts[1]:,} cells", fontsize=9.5)
    if data["diagnostic_result"]:
        fig.suptitle("DIAGNOSTIC: retained failed checks", fontsize=10.5, color="#9c281f")
    return fig


def freeze_and_check_layout(fig):
    """Inspect nested text, including clipping, after freezing the fixed canvas."""
    fig.canvas.draw()
    fig.set_layout_engine("none")
    renderer, canvas = fig.canvas.get_renderer(), fig.bbox
    checked = 0
    for artist in fig.findobj(match=lambda obj: isinstance(obj, Text)):
        if not artist.get_visible() or not artist.get_text().strip():
            continue
        box = artist.get_window_extent(renderer)
        if (box.x0 < canvas.x0 - 1 or box.y0 < canvas.y0 - 1
                or box.x1 > canvas.x1 + 1 or box.y1 > canvas.y1 + 1):
            raise ValueError("Figure text falls outside fixed canvas: " + artist.get_text())
        if artist.get_clip_on():
            clipping = artist.get_clip_box()
            if clipping is not None and (box.x0 < clipping.x0-1 or box.y0 < clipping.y0-1
                                        or box.x1 > clipping.x1+1 or box.y1 > clipping.y1+1):
                raise ValueError("Figure text is clipped by its clip box: " + artist.get_text())
            clip_path = artist.get_clip_path()
            if clip_path is not None:
                path, affine = clip_path.get_transformed_path_and_affine()
                corners = np.array([[box.x0, box.y0], [box.x0, box.y1],
                                    [box.x1, box.y0], [box.x1, box.y1]])
                if not path.contains_points(corners, transform=affine, radius=2).all():
                    raise ValueError("Figure text is clipped by its clip path: " + artist.get_text())
        checked += 1
    return dict(nested_text_artists_checked=checked, canvas_bounds_checked=True,
                text_clip_boxes_and_paths_checked=True, layout_frozen_after_draw=True)


def caption(data):
    counts = [next(r["cells"] for r in data["records"] if r["case"] == case) for case in ["medium", "fine"]]
    status = []
    for row in data["records"]:
        failures = "none" if not row["failed_checks"] else ", ".join(row["failed_checks"])
        status.append(f"{row['solver']} {row['case']}: failed checks {failures}; "
                      + ("passed" if row["cross_code_target_pass"] else "failed") + " code-agreement target")
    changes = data["own_code_changes"]
    gaps = data["cross_code_gaps"]
    return ("Figure 11. Paired airfoil mesh sensitivity at Re=100, zero incidence and R/c=20. "
        f"The medium and fine grids contain {counts[0]:,} and {counts[1]:,} cells. "
        "Both emitted meshes cover the same original straight wall and outer polygon; wall subdivisions, "
        "near-wall spacing, wake spacing and interior triangles change. No shared cell/node correspondence is assumed. "
        "(a) Total drag coefficient. (b) Lift coefficient; grey dotted lines show the unchanged +/-0.001 "
        "zero-angle diagnostic. (c,d) Pressure and viscous drag contributions from native force accounting. "
        "Lines connect two saved results and guide the eye; they do not establish a convergence order. "
        f"Signed fine-minus-medium total Cd changes are {changes['Fluent']['Cd_change_percent']:.6g}% "
        f"for Fluent and {changes['OpenFOAM']['Cd_change_percent']:.6g}% for OpenFOAM. "
        f"Absolute code gaps normalized by Fluent Cd are {gaps['medium']['Cd_difference_percent']:.6g}% "
        f"and {gaps['fine']['Cd_difference_percent']:.6g}% for medium and fine; the separate target remains 0.5%. "
        + "; ".join(status) + ". "
        "Component terms can cancel in total drag. Component replay does not redefine the paired saved endpoints "
        "or their acceptance flags. Every bound source hash, the per-code methods/references, relaxation schedule "
        "and unchanged thresholds were checked. The new fine Fluent session continues the preserved fine saved "
        "fields; OpenFOAM starts cold. Different algorithms, residual norms and initialization preclude a solver "
        "speed/accuracy ranking. No GCI, error bound, convergence order, mesh/domain independence, identified "
        "numerical cause or physical validation follows. Fluent and OpenFOAM computed the forces; Codex prepared "
        "the checks and figure.\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    analysis, analysis_path = helpers()
    source, out = args.data_dir.resolve(), args.output_dir.resolve()
    analysis.require(out.is_relative_to(analysis.ROOT), "Outputs must remain in the E: project")
    if out.exists():
        raise FileExistsError("Use a fresh output directory")
    manifest, paths = analysis.domain.verify_manifest(source)
    analysis.require(set(manifest.get("outputs", {})) == {
        "comparison.json", "comparison.csv", "executed-analyser.py"}, "Missing aggregate output hashes")
    data = json.loads((source / "comparison.json").read_text())
    validate(data)
    fig = make_figure(data)
    qa = freeze_and_check_layout(fig)
    out.mkdir(parents=True)
    base = out / "fig11-airfoil-paired-mesh-sensitivity"
    for extension in ["pdf", "svg", "png"]:
        fig.savefig(base.with_suffix("." + extension), dpi=600)
    plt.close(fig)
    (out / "CAPTION.md").write_text(caption(data), encoding="utf-8")
    (out / "layout-checks.json").write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")
    shutil.copy2(Path(__file__), out / "executed-plotter.py")
    inputs = [source / name for name in ["comparison.json", "comparison.csv", "source-manifest.json"]]
    inputs += [Path(__file__), analysis_path]
    record = dict(script_sha256=analysis.sha256(Path(__file__)), diagnostic_result=data["diagnostic_result"],
        all_numerical_checks_passed=data["all_numerical_checks_passed"],
        all_cross_code_targets_passed=data["all_cross_code_targets_passed"],
        verified_aggregate_inputs=len(paths), canvas_inches=[7.05, 5.05], png_dpi=600,
        software=dict(matplotlib=matplotlib.__version__), layout_checks=qa,
        inputs=[dict(path=str(path), sha256=analysis.sha256(path)) for path in inputs])
    outputs = list(out.iterdir())
    record["outputs"] = {path.name: analysis.sha256(path) for path in outputs if path.is_file()}
    (out / "source-manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(caption(data))


if __name__ == "__main__":
    main()
