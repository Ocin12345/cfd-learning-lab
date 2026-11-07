"""Measure and display the shock from saved Fluent fields."""
import argparse
import json
from math import atan, degrees, radians, tan
from pathlib import Path
import re
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from wedge_common_v001 import CORNER, GAMMA, H, L, LIMITS, NX, NY, P, R, THETA, reference

REF = reference()
REF_SLOPE = tan(radians(REF["shock_angle_deg"]))
P_NORM = Normalize(.98, 1.74)
M_NORM = Normalize(1.6, 2.03)
GEOMETRY_LABEL = "Mach 2 air  ·  10° compression corner  ·  22,400 cells"
COLORS = {"pressure": ("cividis", P_NORM, r"Static pressure ratio, $p/p_1$"),
          "mach": ("viridis", M_NORM, "Mach number")}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "text.color": "#142a3c", "axes.labelcolor": "#233a4b",
                     "axes.edgecolor": "#6a7d89", "xtick.color": "#233a4b", "ytick.color": "#233a4b"})


def read_fluent(run, support, geometry):
    sys.path.insert(0, str(support))
    from analyse_nozzle_fluent_v005 import chunk_array, common_mapping

    state = json.loads((run / "status.json").read_text())
    assert state["stage"] == "closed"
    checkpoints = json.loads((run / "checkpoints.json").read_text())
    last, previous = checkpoints[-1], checkpoints[-2]
    mapping = common_mapping(run / last["case"], geometry)

    def load(checkpoint):
        with h5py.File(run / checkpoint["data"], "r") as f:
            cells = f["results/1/phase-1/cells"]
            result = {}
            for field, source in (("p", "SV_P"), ("rho", "SV_DENSITY"), ("T", "SV_T"),
                                  ("u", "SV_U"), ("v", "SV_V")):
                aligned = np.empty(NX * NY)
                aligned[mapping] = chunk_array(cells[source], NX * NY)
                result[field] = aligned.reshape(NY, NX)
        return fields(result)

    text = (run / last["mass_report"]).read_text()
    native_mass = {name: float(re.search(r"^\s*" + name + r"\s+([-\d.eE+]+)\s*$", text, re.M)[1])
                   for name in ("inlet", "outlet", "top", "ramp", "Net")}
    assert native_mass["inlet"] > 0 > native_mass["outlet"]
    assert abs(native_mass["Net"] - sum(v for k, v in native_mass.items() if k != "Net")) < 1e-4
    return load(last), load(previous), dict(
        final_iteration=last["iteration"], previous_iteration=previous["iteration"],
        scaled_residuals=last["residuals"],
        mesh_quality=json.loads((run / "native-mesh-quality.json").read_text()),
        native_mass_flux_kg_s_per_m_depth=native_mass,
        mass_imbalance=abs(native_mass["Net"] / native_mass["inlet"]),
        native_case=last["case"], native_data=last["data"])


def fields(data):
    assert all(np.isfinite(value).all() for value in data.values())
    assert all((data[name] > 0).all() for name in ("p", "rho", "T"))
    data["mach"] = np.hypot(data["u"], data["v"]) / np.sqrt(GAMMA * R * data["T"])
    data["pressure"] = data["p"] / P
    return data


def measure(data, previous, geo):
    xy = geo["centres"].reshape(NY, NX, 2)
    x, y = xy[..., 0], xy[..., 1]
    pressure = data["pressure"]
    mid = .5 * (1 + REF["pressure_ratio"])
    crossings = []
    for j in np.flatnonzero((x.mean(axis=0) >= .36) & (x.mean(axis=0) <= .84)):
        ym = .5 * (y[1:, j] + y[:-1, j])
        ramp = max(x[:, j].mean() - CORNER, 0) * tan(THETA)
        p_low, p_high = pressure[:-1, j], pressure[1:, j]
        idx = np.flatnonzero((p_low >= mid) & (p_high < mid) & (ym > ramp + .03) & (ym < H - .04))
        if len(idx) == 0:
            continue
        i = idx[np.argmax((p_low - p_high)[idx])]
        frac = (p_low[i] - mid) / (p_low[i] - p_high[i])
        crossings.append((x[i, j] + frac * (x[i + 1, j] - x[i, j]),
                          y[i, j] + frac * (y[i + 1, j] - y[i, j])))
    points = np.asarray(crossings)
    assert len(points) >= 80, len(points)
    slope, intercept = np.polyfit(points[:, 0], points[:, 1], 1)
    dy = y - (slope * x + intercept)
    ramp_clearance = y - np.maximum(x - CORNER, 0) * tan(THETA)
    use = (x >= .44) & (x <= .82) & (ramp_clearance > .04) & (y < H - .035)
    upstream = use & (dy > .07) & (dy < .13)
    downstream = use & (dy < -.07) & (dy > -.13)
    assert upstream.sum() > 100 and downstream.sum() > 100
    p1 = float(np.median(data["p"][upstream]))
    p2 = float(np.median(data["p"][downstream]))
    beta = degrees(atan(slope))
    ratio = p2 / p1
    areas = geo["areas"].reshape(NY, NX)
    stability = {name: float(np.sqrt(np.sum(areas * (data[name] - previous[name])**2)
                                     / np.sum(areas * data[name]**2)))
                 for name in ("p", "mach")}
    return dict(
        shock_angle_deg=beta, shock_angle_error_deg=beta - REF["shock_angle_deg"],
        shock_fit_slope=float(slope), shock_fit_intercept_m=float(intercept),
        extrapolated_shock_foot_x_m=float(-intercept / slope), shock_fit_columns=len(points),
        shock_fit_rms_m=float(np.sqrt(np.mean((points[:, 1] - slope * points[:, 0] - intercept)**2))),
        pressure_upstream_pa=p1, pressure_downstream_pa=p2, pressure_ratio=ratio,
        pressure_ratio_relative_error=ratio / REF["pressure_ratio"] - 1,
        upstream_mach=float(np.median(data["mach"][upstream])),
        downstream_mach=float(np.median(data["mach"][downstream])),
        outlet_minimum_cell_mach=float(data["mach"][:, -1].min()),
        maximum_ideal_gas_relative_error=float(np.max(abs(data["rho"] * R * data["T"] / data["p"] - 1))),
        relative_field_change_l2=stability,
        analytical_angle_check=abs(beta - REF["shock_angle_deg"]) <= LIMITS["shock_angle_error_deg"],
        analytical_pressure_check=abs(ratio / REF["pressure_ratio"] - 1) <= LIMITS["pressure_ratio_error"],
        field_stability_check=max(stability.values()) <= LIMITS["field_stability_relative"])


def contour(ax, data, geo, key, title, theory=False):
    cmap, norm, _ = COLORS[key]
    x = geo["x_nodes"]
    y = geo["y_nodes"]
    img = ax.pcolormesh(np.broadcast_to(x, y.shape), y, data[key],
                        shading="flat", cmap=cmap, norm=norm, rasterized=True)
    ax.fill_between(x, -.03, y[0, :], color="#202f3b", zorder=3)
    ax.plot(x, y[0, :], color="#203344", lw=1.)
    ax.scatter([CORNER], [0], s=12, color="#eb7251", zorder=5)
    if theory:
        xx = np.array([CORNER, L])
        ax.plot(xx, REF_SLOPE * (xx - CORNER), "--", lw=1, color="white", zorder=5,
                dashes=(2, 3))
    ax.set(xlim=(0, L), ylim=(-.025, H), aspect="equal", xlabel="$x$ (m)", ylabel="$y$ (m)")
    ax.set_title(title, loc="left", fontsize=10)
    return img


def fluent_preview(run, support, out):
    out.mkdir(parents=True)
    with np.load(run / "geometry.npz") as geo:
        f, old, native = read_fluent(run, support, geo)
        metrics = measure(f, old, geo)
        fig = plt.figure(figsize=(11.4, 5.55), facecolor="white")
        grid = fig.add_gridspec(2, 2, height_ratios=[1, .045], left=.075, right=.965,
                               top=.79, bottom=.15, wspace=.18, hspace=.28)
        for i, (key, name) in enumerate((("pressure", "Static pressure"), ("mach", "Mach number"))):
            ax = fig.add_subplot(grid[0, i])
            img = contour(ax, f, geo, key, f"{chr(97+i)})  {name}")
            cb = fig.colorbar(img, cax=fig.add_subplot(grid[1, i]), orientation="horizontal")
            cb.set_label(COLORS[key][2])
        fig.text(.075, .94, "An attached oblique shock over a wedge", fontsize=17, weight="bold")
        fig.text(.075, .885, f"ANSYS Fluent 2026 R1  ·  {GEOMETRY_LABEL}", color="#426077")
        fig.text(.075, .038, f"Pressure-jump fit: {metrics['shock_angle_deg']:.2f}°"
                 f"    |    Ideal weak shock: {REF['shock_angle_deg']:.2f}°"
                 f"    |    Saved steady iteration {native['final_iteration']}", fontsize=9)
        fig.savefig(out / "fluent-mach2-wedge.png", dpi=160)
        plt.close(fig)
    payload = dict(native=native, measurements=metrics, weak_shock_reference=REF)
    (out / "fluent-measurements.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fluent", type=Path, required=True)
    parser.add_argument("--support", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    fluent_preview(args.fluent, args.support, args.out)

