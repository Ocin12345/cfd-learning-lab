"""Read the saved jet fields, measure impingement and render the Fluent preview."""
import argparse
import json
from pathlib import Path
import re
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
import numpy as np

from jet_common_v001 import (DX, FIELD_STABILITY, FLOW_SPLIT_ERROR, H, INLET_MASS, MASS_BALANCE,
                             NX, NY, Q_REF, RESIDUAL, UMEAN, W, X, Y, ZONES)

XC = (X[1:] + X[:-1]) / (2 * W)
YC = (Y[1:] + Y[:-1]) / (2 * W)
SPEED_NORM = Normalize(0, 1.5)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "text.color": "#172c40", "axes.labelcolor": "#293f50",
                     "axes.edgecolor": "#71828d", "xtick.color": "#293f50", "ytick.color": "#293f50"})


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def fluent_fields(run, support, geo):
    sys.path.insert(0, str(support))
    from analyse_nozzle_fluent_v005 import chunk_array, common_mapping

    assert read(run / "status.json")["stage"] == "closed"
    checkpoints = read(run / "checkpoints.json")
    a, b = checkpoints[-2:]
    mapping = common_mapping(run / (b["stem"] + ".cas.h5"), geo)

    def load(endpoint):
        with h5py.File(run / (endpoint["stem"] + ".dat.h5"), "r") as f:
            raw = f["results/1/phase-1/cells"]
            fields = {}
            for name, key in (("u", "SV_U"), ("v", "SV_V"), ("p", "SV_P")):
                arr = np.empty(NX * NY)
                arr[mapping] = chunk_array(raw[key], NX * NY)
                fields[name] = arr.reshape(NY, NX)
        return fields

    text = (run / (b["stem"] + "-mass.txt")).read_text()
    native = {name: float(re.search(r"^\s*" + name + r"\s+([-+0-9.eE]+)\s*$", text, re.M)[1])
              for name in (*ZONES, "Net")}
    assert native["inlet"] > 0 > max(native["left"], native["right"])
    mass = {name: -native[name] for name in ZONES}
    mass_imbalance = abs(native["Net"]) / native["inlet"]
    quality = read(run / "mesh-quality.json")
    assert quality["minimum_orthogonal_quality"] >= .2 and quality["minimum_cell_volume_m3"] > 0
    report = dict(final_iteration=b["iteration"], previous_iteration=a["iteration"],
                  residuals=b["residuals"], residual_check=max(b["residuals"].values()) <= RESIDUAL,
                  mesh_quality=quality, mass_imbalance=mass_imbalance, native_endpoint=b["stem"],
                  outward_mass_flux_kg_s_per_m_depth=mass)
    return load(b), load(a), report


def measure(final, previous, report):
    u, v, p = (final[key] for key in ("u", "v", "p"))
    assert all(field.shape == (NY, NX) and np.isfinite(field).all() for field in (u, v, p))
    speed_sq = u**2 + v**2
    velocity_norm = float(speed_sq.sum())
    vel_change = np.sqrt(np.sum((u - previous["u"])**2 + (v - previous["v"])**2) / velocity_norm)
    pressure_change = np.linalg.norm(p - previous["p"]) / np.linalg.norm(p)
    pressure_line = p[0, :] / Q_REF
    mid = float(pressure_line[NX // 2 - 1:NX // 2 + 1].mean())
    peak = int(np.argmax(pressure_line))
    reflection = float(np.sqrt(np.sum((u + u[:, ::-1])**2 + (v - v[:, ::-1])**2) / velocity_norm))
    flow = report["outward_mass_flux_kg_s_per_m_depth"]
    inlet_error = abs(-flow["inlet"] / INLET_MASS - 1.)
    left = flow["left"] / (flow["left"] + flow["right"])
    return dict(
        centre_pressure_coefficient_adjacent_cell=mid,
        peak_pressure_coefficient_adjacent_cell=float(pressure_line[peak]),
        peak_pressure_x_over_slot_width=float(XC[peak]),
        plate_pressure_line_integral_N_per_m_depth=float(p[0].sum() * DX),
        left_outlet_flow_percent=100. * left,
        right_outlet_flow_percent=100. * (1 - left),
        inlet_mass_relative_error=inlet_error,
        mass_imbalance_check=report["mass_imbalance"] <= MASS_BALANCE,
        inlet_profile_flow_check=inlet_error <= .001,
        outlet_split_check=abs(left - .5) <= FLOW_SPLIT_ERROR,
        mirror_velocity_relative_l2=reflection,
        mirror_symmetry_check=reflection <= .005,
        relative_endpoint_change_l2=dict(velocity=float(vel_change), pressure=float(pressure_change)),
        field_stability_check=bool(max(vel_change, pressure_change) <= FIELD_STABILITY),
        pressure_peak_near_centre_check=bool(abs(XC[peak]) <= .10),
        maximum_speed_m_s=float(np.sqrt(speed_sq.max())),
    )


def tank(ax, xlabel=True):
    ax.set(xlim=(-6, 6), ylim=(0, 3), ylabel="$y/W$", aspect="equal")
    ax.set_xticks(np.arange(-6, 7, 2))
    ax.set_yticks([0, 1, 2, 3])
    if xlabel:
        ax.set_xlabel("$x/W$")
    else:
        ax.tick_params(axis="x", labelbottom=False)
    for s in ("left", "right"):
        ax.spines[s].set_linestyle((0, (3, 3)))
    ax.spines["top"].set_visible(False)
    ax.plot([-6, -.5], [3, 3], color="#243747", lw=1.3, clip_on=False)
    ax.plot([.5, 6], [3, 3], color="#243747", lw=1.3, clip_on=False)
    ax.plot([-.5, .5], [3, 3], color="#1292a5", lw=3.0, clip_on=False)


def speed_map(ax, final, streamlines=True):
    u, v = final["u"] / UMEAN, final["v"] / UMEAN
    speed = np.hypot(u, v)
    im = ax.imshow(speed, extent=(-6, 6, 0, 3), origin="lower", cmap="viridis",
                   norm=SPEED_NORM, interpolation="nearest")
    if streamlines:
        ax.streamplot(XC, YC, u, v, color="#f4f3eb", density=(1.8, .70), arrowsize=.65,
                      linewidth=.43 + .48 * np.minimum(speed, 1.), minlength=.11)
    tank(ax, xlabel=False)
    return im


def pressure_map(ax, final, norm=None):
    cp = final["p"] / Q_REF
    if norm is None:
        norm = TwoSlopeNorm(0., vmin=min(-.05, float(cp.min())), vmax=float(cp.max()))
    im = ax.imshow(cp, extent=(-6, 6, 0, 3), origin="lower", cmap="RdBu_r",
                   norm=norm, interpolation="nearest")
    tank(ax)
    return im


def preview(run, support, out):
    with np.load(run / "geometry.npz") as geo:
        final, earlier, native = fluent_fields(run, support, geo)
    metrics = measure(final, earlier, native)
    out.mkdir(parents=True)
    fig = plt.figure(figsize=(11.8, 6.3), facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=[1, .035], left=.074, right=.90,
                           top=.815, bottom=.115, wspace=.04, hspace=.30)
    ax = fig.add_subplot(grid[0, 0])
    im = speed_map(ax, final)
    ax.set_title("a) Speed and streamlines", loc="left")
    fig.colorbar(im, cax=fig.add_subplot(grid[0, 1])).set_label(r"Speed, $|\mathbf{u}|/U_m$")
    ax = fig.add_subplot(grid[1, 0])
    im = pressure_map(ax, final)
    ax.set_title("b) Gauge pressure", loc="left")
    fig.colorbar(im, cax=fig.add_subplot(grid[1, 1])).set_label(r"Pressure coefficient, $C_p$")
    fig.text(.074, .95, "A planar jet splits along a flat plate", weight="bold", size=17)
    fig.text(.074, .895, "ANSYS Fluent 2026 R1  ·  Re = 100  ·  H/W = 3  ·  14,400 cells",
             color="#43617a")
    fig.text(.074, .038, "The cyan segment marks the slot inlet. Dashed sides are pressure outlets.",
             size=9, color="#43617a")
    fig.savefig(out / "fluent-jet-impingement.png", dpi=160)
    plt.close(fig)
    payload = dict(native=native, measurements=metrics)
    (out / "fluent-measurements.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("fluent", "support", "out"):
        parser.add_argument("--" + key, type=Path, required=True)
    args = parser.parse_args()
    preview(args.fluent, args.support, args.out)

