"""Plot checked two-domain airfoil records, retaining every failed diagnostic."""
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
from matplotlib.ticker import ScalarFormatter


def analysis_helpers():
    path = Path(__file__).with_name("analyse_airfoil_domain.py")
    spec = importlib.util.spec_from_file_location("airfoil_domain_plot_helpers", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def make_figure(data):
    rows = data["records"]
    if len(rows) != 4 or {(r["solver"], r["outer_radius_over_chord"]) for r in rows} != {
            (solver, radius) for solver in ["Fluent", "OpenFOAM"] for radius in [20, 40]}:
        raise ValueError("Expected one saved result for each code at R20c and R40c")
    if not all(math.isfinite(r[key]) for r in rows for key in ["Cd", "Cl"]):
        raise ValueError("Nonfinite plotted force coefficient")
    plt.rcParams.update({"font.family": "STIXGeneral", "mathtext.fontset": "stix",
        "font.size": 9.5, "axes.labelsize": 10, "axes.linewidth": .65,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
        "ytick.right": True, "lines.linewidth": 1.05, "legend.fontsize": 8.5,
        "legend.frameon": False, "pdf.fonttype": 42, "svg.fonttype": "path",
        "savefig.bbox": None, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(1, 2, figsize=(7.05, 3.1), layout="constrained")
    for solver, colour, marker, style in [
            ("Fluent", "#0072B2", "o", "-"), ("OpenFOAM", "#D55E00", "s", "--")]:
        own = sorted((r for r in rows if r["solver"] == solver),
                     key=lambda r: r["outer_radius_over_chord"])
        for ax, key in zip(axes, ["Cd", "Cl"]):
            ax.plot([r["outer_radius_over_chord"] for r in own], [r[key] for r in own],
                    label=solver, color=colour, marker=marker, ls=style,
                    ms=4.2, mfc="white", mew=.85)
    target = data["zero_angle_abs_lift_diagnostic_target"]
    if target != .001:
        raise ValueError("This controlled comparison requires the unchanged +/-0.001 lift diagnostic")
    for value in [-target, target]:
        axes[1].axhline(value, color=".5", ls="--", lw=.7, zorder=0)
    axes[1].axhline(0, color=".7", lw=.55, zorder=0)
    axes[0].set_ylabel(r"Drag coefficient, $C_D$")
    axes[1].set_ylabel(r"Lift coefficient, $C_L$")
    axes[0].ticklabel_format(axis="y", style="plain", useOffset=False)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(False)
    formatter.set_useOffset(False)
    axes[1].yaxis.set_major_formatter(formatter)
    low = min(-target, *(r["Cl"] for r in rows))
    high = max(target, *(r["Cl"] for r in rows))
    pad = .15 * (high - low)
    axes[1].set_ylim(low - pad, high + pad)
    for ax, panel in zip(axes, ["(a)", "(b)"]):
        ax.set_xlim(18, 42)
        ax.set_xticks([20, 40])
        ax.set_xlabel(r"Outer boundary radius, $R/c$")
        ax.text(.01, 1.035, panel, transform=ax.transAxes, ha="left", va="bottom")
    axes[0].legend(loc="best")
    if data["diagnostic_result"]:
        fig.suptitle("DIAGNOSTIC: retained failed checks", fontsize=10.5, color="#9c281f")
    return fig


def freeze_and_check_layout(fig):
    fig.canvas.draw()
    fig.set_layout_engine("none")
    renderer = fig.canvas.get_renderer()
    canvas = fig.bbox
    artists = list(fig.texts)
    for ax in fig.axes:
        artists.extend([ax.xaxis.label, ax.yaxis.label, ax.xaxis.get_offset_text(),
                        ax.yaxis.get_offset_text(), *ax.texts,
                        *ax.get_xticklabels(), *ax.get_yticklabels()])
        legend = ax.get_legend()
        if legend:
            artists.extend(legend.get_texts())
    for artist in artists:
        if not artist.get_visible() or not artist.get_text().strip():
            continue
        box = artist.get_window_extent(renderer)
        if (box.x0 < canvas.x0 - 1 or box.y0 < canvas.y0 - 1
                or box.x1 > canvas.x1 + 1 or box.y1 > canvas.y1 + 1):
            raise ValueError("Figure text falls outside fixed canvas: " + artist.get_text())


def caption(data):
    rows = data["records"]
    status = []
    for row in rows:
        failures = row["failed_checks"]
        checks = "own numerical/symmetry checks passed" if not failures else "failed " + ", ".join(failures)
        comparison = "passed" if row["cross_code_target_pass"] else "failed"
        status.append(f"{row['solver']} at R{row['outer_radius_over_chord']:g}c: {checks}; code target {comparison}")
    deltas = data["own_code_Cd_change_percent"]
    return (f"Figure 9. Airfoil outer-domain sensitivity at Re={data['reynolds_number']:g} and zero incidence. "
        "The original R20c inner node coordinates, cell IDs, triangles and wall edges are retained exactly; "
        "the R40c model adds only the outer annulus. "
        "(a) Drag coefficient. (b) Lift coefficient; grey dashed lines mark the unchanged +/-0.001 zero-angle "
        "lift diagnostic. Each code contributes two saved results. Connecting lines guide the eye only. "
        f"Signed own-code drag changes, 100*(Cd_R40-Cd_R20)/abs(Cd_R20), are "
        f"{deltas['Fluent']:.6g}% for Fluent and {deltas['OpenFOAM']:.6g}% for OpenFOAM. "
        + " ".join(status) + ". "
        "Every native input hash and the per-code numerical methods, coefficient references and OpenFOAM "
        "relaxation schedule were checked before aggregation. Failed checks remain in the source table; "
        "this figure does not turn diagnostic results into accepted results. The added annulus discretization "
        "also changes. Two domains establish no error bound, convergence order, mesh/domain independence "
        "or physical validation. Fluent and OpenFOAM computed the forces; Codex prepared the checks and figure.\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out, source = args.output_dir.resolve(), args.data_dir.resolve()
    helper, helper_path = analysis_helpers()
    helper.require(out.is_relative_to(Path("E:/CFD-Learning-Lab").resolve()),
                   "Outputs must remain under the E: project")
    if out.exists():
        raise FileExistsError("Use a fresh output directory")
    data = json.loads((source / "comparison.json").read_text())
    manifest, paths = helper.verify_manifest(source)
    helper.require(data["native_input_hashes_verified"] and data["numerical_methods_and_references_matched"],
                   "Checked domain aggregate required")
    helper.require(data["all_numerical_checks_passed"]
                   == all(row["numerical_checks_passed"] for row in data["records"]),
                   "Aggregate acceptance flags disagree with plotted rows")
    helper.require(not any(data[key] for key in ["domain_independence_established",
                   "error_bound_established", "convergence_order_established", "physical_validation"]),
                   "Unsupported inference in domain aggregate")
    fig = make_figure(data)
    freeze_and_check_layout(fig)
    out.mkdir(parents=True)
    base = out / "fig09-airfoil-domain-sensitivity"
    for extension in ["pdf", "svg", "png"]:
        fig.savefig(base.with_suffix("." + extension), dpi=600)
    plt.close(fig)
    (out / "CAPTION.md").write_text(caption(data), encoding="utf-8")
    shutil.copy2(Path(__file__), out / "executed-plotter.py")
    inputs = [source / name for name in ["comparison.json", "comparison.csv", "source-manifest.json"]]
    inputs += [Path(__file__), helper_path]
    unique = {helper.path_key(path): path for path in inputs}
    record = dict(script_sha256=helper.sha256(Path(__file__)), diagnostic_result=data["diagnostic_result"],
        all_numerical_checks_passed=data["all_numerical_checks_passed"],
        all_cross_code_targets_passed=data["all_cross_code_targets_passed"],
        native_manifest_inputs_verified=len(paths), canvas_inches=[7.05, 3.1], png_dpi=600,
        software=dict(matplotlib=matplotlib.__version__),
        inputs=[dict(path=str(path.resolve()), sha256=helper.sha256(path)) for path in unique.values()])
    (out / "source-manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(caption(data))


if __name__ == "__main__":
    main()
