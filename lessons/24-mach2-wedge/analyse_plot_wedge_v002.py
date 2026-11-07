"""Compare native wedge fields and animate the OpenFOAM startup."""
import argparse
import csv
import json
from math import tan
from pathlib import Path
import re
import sys

import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LogNorm
import numpy as np

from plot_wedge_fluent_v001 import COLORS, GEOMETRY_LABEL, REF, REF_SLOPE, contour, fields, measure, read_fluent
from prepare_wedge_openfoam_v001 import END_TIME, FRAME_DT, MAX_CO
from wedge_common_v001 import CORNER, DEPTH, H, L, LIMITS, NX, NY, P, R, RESIDUALS, T, THETA, U

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
SL = ("gray_r", LogNorm(.3, 120), r"Density-gradient magnitude (kg/m$^4$)")
BLUE, ORANGE = "#075a74", "#a85124"


def read_openfoam(run, geo):
    from analyse_nozzle_openfoam_v002 import FIELDS, field, patch_values

    assert (run / "run-status.txt").read_text().strip() == "complete"
    mesh_log = (run / "mesh-check.log").read_text()
    assert "Mesh OK." in mesh_log and re.search(r"\bcells:\s+22400\b", mesh_log)
    log = (run / "solver.log").read_text()
    assert re.search(r"^End\s*$", log, re.M) and "Finalising parallel run" in log
    assert "FOAM FATAL" not in log
    step_times = np.array([float(t) for t in re.findall(r"^Time = (" + NUMBER + r")s?\s*$", log, re.M)])
    assert len(step_times) > 0 and np.all(np.diff(step_times) > 0)
    assert abs(step_times[-1] - END_TIME) < 1e-10
    co = np.array([float(c) for c in re.findall(
        r"^Courant Number mean: " + NUMBER + r" max: (" + NUMBER + r")", log, re.M)])
    assert len(co) == len(step_times) + 1 and co.max() <= MAX_CO * 1.02  # Includes startup report.

    count = NX * NY
    native_centres = field(run / "0.00025/C", count, *FIELDS["C"])
    x, y, z = native_centres.T
    lower = np.maximum(x - CORNER, 0) * tan(THETA)
    i = np.floor(x * NX / L).astype(int)
    j = np.floor((y - lower) * NY / (H - lower)).astype(int)
    order = j * NX + i
    assert ((i >= 0) & (i < NX) & (j >= 0) & (j < NY)).all()
    assert len(np.unique(order)) == count
    coord_error = np.linalg.norm(native_centres[:, :2] - geo["centres"][order], axis=1).max()
    assert coord_error < 2e-7 and np.max(abs(z - DEPTH / 2)) < 1e-9
    to_common = np.argsort(order)

    frame_times = np.arange(round(END_TIME / FRAME_DT) + 1) * FRAME_DT
    assert "internalField uniform 101325;" in (run / "0/p").read_text()
    assert "internalField uniform 300;" in (run / "0/T").read_text()
    assert f"internalField uniform ({U:.16g} 0 0);" in (run / "0/U").read_text()
    zeros = np.zeros((NY, NX))
    frames = [fields(dict(p=zeros + P, T=zeros + T, rho=zeros + P / (R * T),
                          u=zeros + U, v=zeros.copy()))]
    mass_rows = []
    for t in frame_times[1:]:
        directory = run / f"{t:.9g}"
        scalars = {name: field(directory / name, count, *FIELDS[name])[to_common].reshape(NY, NX)
                   for name in ("p", "rho", "T")}
        vector = field(directory / "U", count, *FIELDS["U"])[to_common].reshape(NY, NX, 3)
        assert abs(vector[..., 2]).max() < 1e-9
        frames.append(fields(dict(**scalars, u=vector[..., 0], v=vector[..., 1])))
        phi = directory / "phi"
        assert re.search(r"\bdimensions\s*\[\s*1\s+0\s+-1\s+0\s+0\s+0\s+0\s*\];", phi.read_text())
        fluxes = {name: float(patch_values(phi, name, faces).sum()) / DEPTH for name, faces
                  in (("inlet", NY), ("outlet", NY), ("top", NX), ("ramp", NX))}
        assert fluxes["inlet"] < 0 < fluxes["outlet"]
        mass_rows.append(dict(time_s=float(t), **fluxes,
                              relative_imbalance=abs(sum(fluxes.values()) / fluxes["inlet"])))
    native = dict(final_time_s=float(frame_times[-1]), comparison_time_s=float(frame_times[-5]),
                  time_steps=len(step_times), maximum_reported_courant=float(co.max()),
                  mesh_check="Mesh OK.", native_centroid_alignment_error_m=float(coord_error),
                  final_native_outward_mass_flux_kg_s_per_m_depth=mass_rows[-1],
                  mass_imbalance=mass_rows[-1]["relative_imbalance"],
                  last_five_frames_max_mass_imbalance=max(row["relative_imbalance"] for row in mass_rows[-5:]),
                  initial_frame="Prescribed uniform field")
    return frame_times, frames, mass_rows, native


def native_panels(out, geo, solutions):
    fig = plt.figure(figsize=(10.6, 8.7), facecolor="white")
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 1, .05], left=.075, right=.965,
                           top=.835, bottom=.135, wspace=.15, hspace=.29)
    for row, (name, data) in enumerate(solutions.items()):
        for col, key in enumerate(("pressure", "mach")):
            ax = fig.add_subplot(grid[row, col])
            label = "Pressure" if key == "pressure" else "Mach number"
            contour(ax, data, geo, key, f"{chr(97+2*row+col)})  {name}  ·  {label}")
            if row == 0:
                ax.set_xlabel("")
    for col, key in enumerate(("pressure", "mach")):
        cmap, norm, label = COLORS[key]
        bar = fig.colorbar(plt.cm.ScalarMappable(norm, cmap), cax=fig.add_subplot(grid[2, col]),
                           orientation="horizontal")
        bar.set_label(label)
    fig.text(.075, .956, "Mach 2 flow turns through an oblique shock", weight="bold", size=17)
    fig.text(.075, .917, GEOMETRY_LABEL + " per solver", color="#426077")
    fig.text(.075, .048, "Saved Fluent steady solution and OpenFOAM field at 8.00 ms.",
             size=9, color="#426077")
    fig.savefig(out / "mach2-wedge-comparison.png", dpi=170)
    plt.close(fig)


def pressure_profile(out, geo, solutions):
    fig, ax = plt.subplots(figsize=(7.8, 4.9), layout="constrained")
    xy = geo["centres"].reshape(NY, NX, 2)
    j = int(np.argmin(abs(xy[..., 0].mean(0) - .60)))
    x = xy[:, j, 0]
    delta = (xy[:, j, 1] - REF_SLOPE * (x - CORNER)) / np.sqrt(1 + REF_SLOPE**2)
    x_nominal = float(x.mean())
    ax.plot([-.18, 0, 0, .18], [REF["pressure_ratio"], REF["pressure_ratio"], 1, 1],
            color="#647887", lw=1.3, linestyle=":", label="Ideal weak shock")
    columns = [delta]
    for (name, data), color, style in zip(solutions.items(), (BLUE, ORANGE), ("-", "--")):
        values = data["pressure"][:, j]
        columns.append(values)
        ax.plot(delta, values, color=color, lw=2, linestyle=style, label=name)
    ax.set(xlim=(-.17, .17), ylim=(.96, 1.76), xlabel="Signed distance from ideal shock (m)",
           ylabel=r"Static pressure ratio, $p/p_1$")
    ax.set_title(f"Pressure rise across the shock at x = {x_nominal:.3f} m", loc="left", weight="bold")
    ax.grid(color="#dfe7eb", lw=.7)
    ax.legend(frameon=False, loc="center right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "mach2-wedge-pressure-profile.png", dpi=180)
    fig.savefig(out / "mach2-wedge-pressure-profile.svg")
    plt.close(fig)
    with (out / "shock-pressure-profile.csv").open("x", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["signed_distance_from_ideal_shock_m", "fluent_p_over_p1", "openfoam_p_over_p1"])
        writer.writerows(zip(*columns))


def gradient_density(frame, geo):
    xy = geo["centres"].reshape(NY, NX, 2)
    x, y = xy[..., 0], xy[..., 1]
    rho = frame["rho"]
    drdy = np.gradient(rho, axis=0) / np.gradient(y, axis=0)
    drdx = (np.gradient(rho, axis=1) - drdy * np.gradient(y, axis=1)) / np.gradient(x, axis=1)
    return np.maximum(np.hypot(drdx, drdy), .3)


def startup_movie(out, geo, times, frames):
    fig = plt.figure(figsize=(11.4, 5.65), facecolor="white")
    grid = fig.add_gridspec(2, 2, height_ratios=[1, .045], left=.075, right=.965,
                           top=.79, bottom=.15, wspace=.18, hspace=.28)
    ax0, ax1 = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    p = contour(ax0, frames[0], geo, "pressure", "a)  Static pressure")
    cmap, norm, label = SL
    x, y = geo["x_nodes"], geo["y_nodes"]
    schlieren = ax1.pcolormesh(np.broadcast_to(x, y.shape), y, gradient_density(frames[0], geo),
                              cmap=cmap, norm=norm, shading="flat")
    ax1.fill_between(x, -.03, y[0, :], color="#202f3b")
    ax1.scatter([CORNER], [0], s=12, color="#eb7251", zorder=5)
    ax1.set(xlim=(0, L), ylim=(-.025, H), aspect="equal", xlabel="$x$ (m)", ylabel="$y$ (m)")
    ax1.set_title("b)  Numerical schlieren", fontsize=10, loc="left")
    for j, (im, title) in enumerate(((p, COLORS["pressure"][2]), (schlieren, label))):
        cb = fig.colorbar(im, cax=fig.add_subplot(grid[1, j]), orientation="horizontal")
        cb.set_label(title)
    fig.text(.075, .945, "How the oblique shock forms", size=17, weight="bold")
    fig.text(.075, .892, "OpenFOAM Foundation 14  ·  transient shockFluid  ·  Mach 2 over a 10° ramp",
             color="#426077")
    clock = fig.text(.965, .944, "", ha="right", fontfamily="monospace", size=13)
    fig.text(.075, .032, "Saved fields at 0.25 ms intervals. The initial frame is the prescribed freestream.",
             size=9, color="#426077")
    gradients = [gradient_density(frame, geo) for frame in frames]

    def draw(k):
        p.set_array(frames[k]["pressure"].ravel())
        schlieren.set_array(gradients[k].ravel())
        clock.set_text(f"t = {1000 * times[k]:.2f} ms")
        return p, schlieren, clock

    # Hold the prescribed start and final snapshot so the loop is readable.
    indices = [0] * 5 + list(range(1, len(times))) + [len(times) - 1] * 7
    animation = FuncAnimation(fig, draw, frames=indices, interval=125, blit=False)
    animation.save(out / "mach2-wedge-shock-formation.gif", writer=PillowWriter(fps=8), dpi=95)
    plt.close(fig)


def main(args):
    sys.path.insert(0, str(args.support))
    with np.load(args.fluent / "geometry.npz") as geo:
        f, f_prev, f_native = read_fluent(args.fluent, args.support, geo)
        times, frames, fluxes, o_native = read_openfoam(args.openfoam, geo)
        f_checks = measure(f, f_prev, geo)
        o_checks = measure(frames[-1], frames[-5], geo)
        for checks, native in ((f_checks, f_native), (o_checks, o_native)):
            checks["mass_balance_check"] = native["mass_imbalance"] <= LIMITS["mass_imbalance"]
            checks["all_outlet_cells_supersonic"] = checks["outlet_minimum_cell_mach"] > 1.
            checks["ideal_gas_consistency_check"] = checks["maximum_ideal_gas_relative_error"] < 1e-5
        f_checks["residual_check"] = all(f_native["scaled_residuals"][name] <= target
                                         for name, target in RESIDUALS.items())
        o_checks["late_mass_balance_check"] = o_native["last_five_frames_max_mass_imbalance"] <= LIMITS["mass_imbalance"]
        payload = dict(cells_per_case=NX * NY, reference=REF, thresholds=LIMITS,
                       numerical_comparison_only=True, physical_validation=False,
                       cases={"ANSYS Fluent": dict(native=f_native, measurements=f_checks),
                              "OpenFOAM": dict(native=o_native, measurements=o_checks)})
        args.out.mkdir(parents=True)
        (args.out / "wedge-summary.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
        with (args.out / "openfoam-mass-flow.csv").open("x", newline="") as fobj:
            writer = csv.DictWriter(fobj, list(fluxes[0]))
            writer.writeheader()
            writer.writerows(fluxes)
        solutions = {"ANSYS Fluent": f, "OpenFOAM": frames[-1]}
        native_panels(args.out, geo, solutions)
        pressure_profile(args.out, geo, solutions)
        startup_movie(args.out, geo, times, frames)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for flag in ("fluent", "openfoam", "out", "support"):
        parser.add_argument("--" + flag, type=Path, required=True)
    main(parser.parse_args())

