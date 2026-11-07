"""Plot saved native dam-break volume fractions on a common cell grid."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap, Normalize
import numpy as np

from dam_break_fluent_v004 import DX, DT, FRAME_INTERVAL, H, L, LAST_STEP, NX, NY, WATER_HEIGHT, WATER_WIDTH, alpha_group, cell_values

PALETTE = LinearSegmentedColormap.from_list(
    "water", ["#f6f4e9", "#bed9d8", "#53a4ae", "#08788c", "#064562"])
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.labelcolor": "#203446", "text.color": "#142c42",
                     "xtick.color": "#203446", "ytick.color": "#203446",
                     "axes.edgecolor": "#617785"})


def load_fluent(run, scripts):
    sys.path.insert(0, str(scripts))
    from analyse_nozzle_fluent_v005 import common_mapping

    with np.load(run / "geometry.npz") as geo:
        mapping = common_mapping(run / "tank.cas.h5", geo)
        expected_initial = ((geo["centres"][:, 0] < WATER_WIDTH)
                            & (geo["centres"][:, 1] < WATER_HEIGHT))
    with (run / "frames.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    times, frames = [], []
    for row in rows:
        with h5py.File(run / row["file"], "r") as f:
            native = cell_values(alpha_group(f))
        aligned = np.empty_like(native)
        aligned[mapping] = native
        times.append(float(row["time_s"]))
        frames.append(aligned.reshape(NY, NX))
    assert len(frames) and np.allclose(frames[0].ravel(), expected_initial, atol=2e-15, rtol=0.)
    assert np.allclose(times, np.arange(len(times)) * DT * FRAME_INTERVAL, atol=1e-10, rtol=0.)
    return np.asarray(times), np.asarray(frames)


def foam_field(path, kind):
    content = path.read_text(encoding="ascii")
    match = re.search(r"internalField\s+nonuniform\s+List<" + kind + r">\s+(\d+)\s*\(", content)
    if not match:
        raise ValueError(f"Missing native internal {kind} field in {path}")
    close = re.search(r"\n\)\s*;", content[match.end():])
    if not close:
        raise ValueError(f"Unterminated field in {path}")
    end = match.end() + close.start()
    n = int(match.group(1))
    data = np.fromstring(content[match.end():end].replace("(", " ").replace(")", " "), sep=" ")
    result = data.reshape(n, 3) if kind == "vector" else data
    assert len(result) == NX * NY and np.isfinite(result).all()
    return result


def load_openfoam(run, times, initial):
    initial_dict = (run / "0/alpha.water").read_text()
    assert "type zonal" in initial_dict and "defaultValue 0" in initial_dict
    assert f"({WATER_WIDTH} {WATER_HEIGHT} 1)" in initial_dict and "value 1" in initial_dict
    centres = foam_field(run / "0.01/C", "vector")
    xy = np.rint(centres[:, :2] / DX - .5).astype(int)
    assert ((xy >= 0) & (xy < [NX, NY])).all()
    assert np.max(np.abs(centres[:, :2] - (xy + .5) * DX)) < 1e-10
    mapping = xy[:, 0] + NX * xy[:, 1]
    assert len(np.unique(mapping)) == NX * NY
    frames = [np.rint(initial)]
    for t in times[1:]:
        native = foam_field(run / f"{t:.8g}" / "alpha.water", "scalar")
        aligned = np.empty_like(native)
        aligned[mapping] = native
        frames.append(aligned.reshape(NY, NX))
    return np.asarray(frames)


def tank(ax, field, label, show_xlabel=True, show_ylabel=True):
    img = ax.imshow(field, origin="lower", extent=(0, L, 0, H),
                    interpolation="nearest", cmap=PALETTE, vmin=0, vmax=1)
    ax.set_xticks([0, .15, .30, .45, .60])
    ax.set_yticks([0, .10, .20, .30])
    ax.set_xlabel("$x$ (m)" if show_xlabel else "")
    ax.set_ylabel("$y$ (m)" if show_ylabel else "")
    ax.tick_params(axis="x", labelbottom=show_xlabel)
    ax.tick_params(axis="y", labelleft=show_ylabel)
    ax.set_title(label, loc="left", fontweight="normal", fontsize=10)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(.9)
    return img


def bar(fig, cax):
    c = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0, 1), cmap=PALETTE),
                     cax=cax, orientation="horizontal", ticks=[0, .25, .50, .75, 1.])
    c.set_label(r"Water volume fraction, $\alpha_w$")
    c.outline.set_linewidth(.6)


def overview(out, times, fields):
    selected = np.unique(np.rint(np.linspace(0, len(times) - 1, 4)).astype(int))
    paired = len(fields) == 2
    if paired:
        selected = [10, 28, 40, len(times) - 1]
    fig = plt.figure(figsize=(11, 5.55 if paired else 5.8), facecolor="white")
    if paired:
        grid = fig.add_gridspec(2, 4, left=.07, right=.98, top=.79, bottom=.24,
                               wspace=.16, hspace=.22)
        for row, (name, values) in enumerate(fields.items()):
            for col, k in enumerate(selected):
                title = f"{name}  ·  {times[k]:.2f} s"
                tank(fig.add_subplot(grid[row, col]), values[k], title,
                     show_xlabel=row == 1, show_ylabel=col == 0)
    else:
        grid = fig.add_gridspec(2, 2, left=.08, right=.96, top=.81, bottom=.23,
                               wspace=.16, hspace=.40)
        for panel, k in enumerate(selected):
            tank(fig.add_subplot(grid[panel // 2, panel % 2]), next(iter(fields.values()))[k],
                 f"{chr(97 + panel)})  t = {times[k]:.2f} s",
                 show_xlabel=panel // 2 == 1, show_ylabel=panel % 2 == 0)
    fig.text(.07, .94, "Gravity-driven collapse of a water column", weight="bold", size=17)
    subtitle = "Matched native solver snapshots" if paired else "ANSYS Fluent 2026 R1  ·  saved native volume fractions"
    fig.text(.07, .89, subtitle + "  ·  2-D closed tank  ·  5 mm cells", size=10, color="#47617a")
    bar(fig, fig.add_axes([.35, .105, .30, .028]))
    fig.savefig(out / ("dam-break-comparison.png" if paired else "dam-break-fluent.png"), dpi=160)
    plt.close(fig)


def movie(out, times, fields, ffmpeg):
    fig = plt.figure(figsize=(11, 4.35), facecolor="white")
    grid = fig.add_gridspec(1, 2, left=.075, right=.975, top=.77, bottom=.24, wspace=.22)
    artists = [tank(fig.add_subplot(grid[0, j]), values[0], f"{chr(97 + j)})  {name}")
               for j, (name, values) in enumerate(fields.items())]
    fig.text(.075, .94, "Dam break in a closed tank", weight="bold", size=17)
    fig.text(.075, .88, "Native VOF fields  ·  7,200 cells per solver  ·  0.01 s between saved frames",
             size=10, color="#47617a")
    clock = fig.text(.975, .92, "", ha="right", fontfamily="monospace", size=13)
    bar(fig, fig.add_axes([.35, .12, .30, .036]))

    def draw(k):
        for artist, values in zip(artists, fields.values()):
            artist.set_data(values[k])
        clock.set_text(f"t = {times[k]:.2f} s")
        return *artists, clock

    frames = [0] * 10 + list(range(1, len(times))) + [len(times) - 1] * 12
    anim = FuncAnimation(fig, draw, frames=frames, interval=50, blit=False)
    anim.save(out / "dam-break-paired-native.gif", writer=PillowWriter(fps=20), dpi=95)
    if ffmpeg:
        plt.rcParams["animation.ffmpeg_path"] = str(ffmpeg)
        anim.save(out / "dam-break-paired-native.mp4", dpi=120,
                  writer=FFMpegWriter(fps=20, codec="libx264", extra_args=["-pix_fmt", "yuv420p"]))
    plt.close(fig)


def diagnostics(out, times, fields):
    reference = WATER_WIDTH * WATER_HEIGHT
    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(7.6, 7.5), layout="constrained")
    fig.suptitle("Dam-break front and water-volume conservation", weight="bold", size=14)
    results = {}
    with (out / "dam-break-history.csv").open("x", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["solver", "time_s", "water_area_m2", "relative_drift_percent",
                         "floor_front_x_m", "right_wall_crest_y_m"])
        for (name, values), color, style in zip(fields.items(), ("#086878", "#ad5e25"), ("-", "--")):
            area = values.sum(axis=(1, 2)) * DX**2
            drift = 100 * (area / reference - 1)
            wet = values >= .5
            x = (np.arange(NX) + .5) * DX
            y = (np.arange(NY) + .5) * DX
            front = np.max(np.where(wet[:, 0, :], x, 0), axis=1)
            crest = np.max(np.where(wet[:, :, -1], y, 0), axis=1)
            arrival = np.flatnonzero(wet[:, 0, -1])
            results[name] = {
                "frames": len(values), "initial_water_area_m2": float(area[0]),
                "final_water_area_m2": float(area[-1]),
                "max_absolute_water_area_drift_percent": float(abs(drift).max()),
                "alpha_min": float(values.min()), "alpha_max": float(values.max()),
                "max_wall_crest_sampled_m": float(crest.max()),
                "max_wall_crest_time_s": float(times[np.argmax(crest)]),
                "first_wet_floor_cell_at_right_wall_s": float(times[arrival[0]]) if len(arrival) else None,
                "water_conservation_within_0_5_percent": bool(abs(drift).max() <= .5),
                "alpha_bounds_pass": bool(values.min() >= -1e-6 and values.max() <= 1 + 1e-6),
            }
            axes[0].plot(times, front, style, color=color, label=name, linewidth=1.9)
            axes[1].plot(times, crest, style, color=color, label=name, linewidth=1.9)
            axes[2].plot(times, drift, style, color=color, label=name, linewidth=1.9)
            writer.writerows(zip([name] * len(times), times, area, drift, front, crest))
    for ax in axes:
        ax.grid(axis="y", color="#dfe8eb")
        ax.legend(frameon=False)
    axes[0].set_ylabel("Floor front $x$ (m)")
    axes[1].set_ylabel("Right-wall crest $y$ (m)")
    axes[2].set_ylabel("Water volume drift (%)")
    axes[2].set_xlabel("Physical time (s)")
    axes[2].axhline(0., color="#8ca0a6", lw=.7, zorder=0)
    fig.savefig(out / "dam-break-diagnostics.png", dpi=180)
    plt.close(fig)
    (out / "dam-break-summary.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fluent", type=Path, required=True)
    p.add_argument("--openfoam", type=Path)
    p.add_argument("--scripts", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--ffmpeg", type=Path)
    args = p.parse_args()
    times, fluent = load_fluent(args.fluent, args.scripts)
    fields = {"ANSYS Fluent": fluent}
    if args.openfoam:
        assert len(times) == LAST_STEP // FRAME_INTERVAL + 1
        fields["OpenFOAM"] = load_openfoam(args.openfoam, times, fluent[0])
    args.out.mkdir(parents=True)
    overview(args.out, times, fields)
    if args.openfoam:
        diagnostics(args.out, times, fields)
        movie(args.out, times, fields, args.ffmpeg)


if __name__ == "__main__":
    main()

