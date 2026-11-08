"""Read and plot native Fluent temperature snapshots without starting a solver."""
import argparse
import csv
from pathlib import Path
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np

from heated_jet_common_v001 import COLD, DT, HOT, STRIDE
from jet_common_v001 import NX, NY

CMAP = "inferno"
NORM = Normalize(COLD, HOT)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "text.color": "#182f42", "axes.labelcolor": "#203648",
                     "axes.edgecolor": "#687d8b", "xtick.color": "#203648", "ytick.color": "#203648"})


def fluent_temperatures(run, geo):
    from analyse_nozzle_fluent_v005 import chunk_array, common_mapping
    mapping = common_mapping(run / "heated-jet.cas.h5", geo)
    with (run / "frames.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    frames = {}
    for row in rows:
        step = int(row["step"])
        assert step not in frames and abs(float(row["time_s"]) - step * DT) < 1e-9
        with h5py.File(run / row["file"], "r") as h:
            native = chunk_array(h["results/1/phase-1/cells/SV_T"], NX * NY)
        T = np.empty_like(native)
        T[mapping] = native
        frames[step] = T.reshape(NY, NX)
    assert np.allclose(frames[0], COLD, atol=1e-10, rtol=0)
    assert sorted(k for k in frames if k % STRIDE == 0) == list(range(0, max(frames) + 1, STRIDE))
    return frames, mapping


def panel(ax, T, label, xlabel=False):
    im = ax.imshow(T, origin="lower", extent=(-6, 6, 0, 3), cmap=CMAP,
                   norm=NORM, interpolation="nearest")
    ax.set(xlim=(-6, 6), ylim=(0, 3), ylabel="$y/W$", aspect="equal")
    ax.set_xticks(np.arange(-6, 7, 2))
    ax.set_yticks([0, 1, 2, 3])
    ax.set_title(label, loc="left", size=10)
    ax.spines["bottom"].set_color("#db6939")
    ax.spines["bottom"].set_linewidth(2.5)
    ax.spines["top"].set_visible(False)
    ax.plot([-6, -.5], [3, 3], color="#435365", lw=1., clip_on=False)
    ax.plot([.5, 6], [3, 3], color="#435365", lw=1., clip_on=False)
    ax.plot([-.5, .5], [3, 3], color="#2ab2be", lw=2.8, clip_on=False)
    if xlabel:
        ax.set_xlabel("$x/W$")
    else:
        ax.tick_params(axis="x", labelbottom=False)
    return im


def preview(run, flow, out):
    with np.load(flow / "geometry.npz") as geo:
        frames, _ = fluent_temperatures(run, geo)
    steps = [k for k in frames if k % STRIDE == 0]
    assert len(steps) >= 4
    selected = np.linspace(0, len(steps) - 1, 4).round().astype(int)
    out.mkdir(parents=True)
    fig = plt.figure(figsize=(12., 6.4), facecolor="white")
    grid = fig.add_gridspec(2, 3, width_ratios=[1, 1, .046], left=.074, right=.914,
                           top=.82, bottom=.15, wspace=.14, hspace=.31)
    for i, index in enumerate(selected):
        k = steps[index]
        ax = fig.add_subplot(grid[i // 2, i % 2])
        im = panel(ax, frames[k], f"{chr(97+i)})  t = {k * DT:.2f} s", xlabel=i >= 2)
    fig.colorbar(im, cax=fig.add_subplot(grid[:, 2])).set_label("Air temperature (K)")
    fig.text(.074, .944, "A cold slot jet cools a hot plate", weight="bold", size=17)
    fig.text(.074, .886, "ANSYS Fluent 2026 R1  ·  320 K plate  ·  temperature evolves in a fixed airflow",
             color="#426178")
    fig.text(.074, .047, "Orange marks the heated plate. Cyan marks the 300 K jet inlet.", size=9,
             color="#426178")
    fig.savefig(out / "fluent-heated-jet.png", dpi=155)
    plt.close(fig)
    print(f"Saved actual Fluent fields through t = {steps[-1] * DT:.2f} s")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for name in ("fluent", "flow", "out", "support"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.support))
    preview(args.fluent, args.flow, args.out)

