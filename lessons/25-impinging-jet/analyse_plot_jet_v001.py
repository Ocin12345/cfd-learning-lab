"""Compare saved planar-jet runs and animate tracers through their steady fields."""
import argparse
import csv
import json
from pathlib import Path
import re
import sys

from plot_jet_fluent_v002 import XC, YC, fluent_fields, measure, speed_map

import h5py
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.colors import to_rgb
import numpy as np

from jet_common_v001 import DEPTH, DX, HALF, H, NX, NY, NSLOT, Q_REF, RESIDUAL, RHO, UMEAN, W, ZONES
from prepare_jet_openfoam_v001 import END, INTERVAL

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
N_FACE = 2 * NX * NY + NX + NY


def fluent_boundaries(run, native):
    from analyse_nozzle_fluent_v005 import chunk_array
    stem = native["native_endpoint"]
    with h5py.File(run / (stem + ".cas.h5"), "r") as f:
        root = f["meshes/1"]
        topology = root["faces/zoneTopology"]
        names = topology["name"][0].decode().split(";")
        slices = dict(zip(names, [slice(a - 1, b) for a, b in zip(
            topology["minId"][:].astype(int), topology["maxId"][:].astype(int))]))
        conn = np.empty((N_FACE, 2), dtype=int)
        for group in root["faces/nodes"].values():
            a = int(np.asarray(group.attrs["minId"]).flat[0]) - 1
            b = int(np.asarray(group.attrs["maxId"]).flat[0])
            assert np.all(group["nnodes"][:] == 2)
            conn[a:b] = group["nodes"][:].reshape(-1, 2) - 1
        nodes = chunk_array(root["nodes/coords"], (NX + 1) * (NY + 1), (2,))
        plate = slices["plate"]
        points = nodes[conn[plate]].mean(axis=1)
    order = np.argsort(points[:, 0])
    assert np.allclose(points[order, 0] / W, XC, atol=1e-10, rtol=0)
    assert abs(points[:, 1]).max() < 1e-12
    with h5py.File(run / (stem + ".dat.h5"), "r") as f:
        faces = f["results/1/phase-1/faces"]
        wall_pressure = chunk_array(faces["SV_P"], N_FACE)[plate][order]
        face_flux = chunk_array(faces["SV_FLUX"], N_FACE)
    raw_in = face_flux[slices["inlet"]].sum()
    outward = -np.sign(raw_in) * face_flux
    for name, total in native["outward_mass_flux_kg_s_per_m_depth"].items():
        assert abs(outward[slices[name]].sum() - total) < 1e-9
    native["outlet_reverse_mass_kg_s_per_m_depth"] = {
        name: float(-np.minimum(outward[slices[name]], 0.).sum()) for name in ("left", "right")}
    native["plate_pressure_source"] = "Fluent native face SV_P on the plate"
    return wall_pressure


def openfoam_fields(run, geo):
    from analyse_nozzle_openfoam_v002 import FIELDS, field, patch_values
    assert (run / "run-status.txt").read_text().strip() == "complete"
    log = (run / "solver.log").read_text()
    mesh = (run / "mesh-check.log").read_text()
    assert "Mesh OK." in mesh and re.search(r"\bcells:\s+14400\b", mesh)
    assert re.search(r"^End\s*$", log, re.M) and "Finalising parallel run" in log
    sections = re.split(r"^Time = (" + NUMBER + r")s?\s*$", log, flags=re.M)
    assert [float(t) for t in sections[1::2]] == list(range(1, END + 1))
    late = []
    for section in sections[-1:-100:-2]:
        values = {}
        for name, value in re.findall(r"Solving for (Ux|Uy|p), Initial residual = (" + NUMBER + r")", section):
            values[name] = max(values.get(name, 0.), float(value))
        assert set(values) == {"Ux", "Uy", "p"}
        late.append(values)
    assert len(late) == 50
    maxima = {name: max(row[name] for row in late) for name in late[0]}

    count = NX * NY
    centres = field(run / str(END) / "C", count, *FIELDS["C"])
    i = np.rint((centres[:, 0] + HALF) / DX - .5).astype(int)
    j = np.rint(centres[:, 1] / DX - .5).astype(int)
    assert ((i >= 0) & (i < NX) & (j >= 0) & (j < NY)).all()
    mapping = i + NX * j
    assert len(np.unique(mapping)) == count
    err = float(np.abs(centres[:, :2] - geo["centres"][mapping]).max())
    assert err < 1e-10 and abs(centres[:, 2] - DEPTH / 2).max() < 1e-10
    order = np.argsort(mapping)

    def load(it):
        path = run / str(it)
        uv = field(path / "U", count, *FIELDS["U"])[order].reshape(NY, NX, 3)
        assert abs(uv[..., 2]).max() < 1e-10
        p = field(path / "p", count, 1, [0, 2, -2, 0, 0, 0, 0])[order].reshape(NY, NX) * RHO
        return dict(u=uv[..., 0], v=uv[..., 1], p=p)

    def flux(it):
        path = run / str(it) / "phi"
        assert re.search(r"\bdimensions\s*\[\s*0\s+3\s+-1\s+0\s+0\s+0\s+0\s*\];", path.read_text())
        per_face = {name: patch_values(path, name, size) * RHO / DEPTH
                    for name, size in (("left", NY), ("right", NY), ("inlet", NSLOT),
                                       ("plate", NX), ("ceiling", NX - NSLOT))}
        mass = {name: float(arr.sum()) for name, arr in per_face.items()}
        assert mass["inlet"] < 0 < min(mass["left"], mass["right"])
        net = abs(sum(mass.values())) / -mass["inlet"]
        reverse = {name: float(-np.minimum(per_face[name], 0.).sum()) for name in ("left", "right")}
        return mass, net, reverse

    mass, net, reverse = flux(END)
    assert re.search(r"\bplate\s*\{\s*type\s+zeroGradient\s*;", (run / str(END) / "p").read_text())
    final, previous = load(END), load(END - 2 * INTERVAL)
    native = dict(final_iteration=END, previous_iteration=END - 2 * INTERVAL,
                  native_mesh_check="Mesh OK.", max_centroid_difference_m=err,
                  final_initial_residuals=late[0], max_initial_residuals_last_50=maxima,
                  residual_check=max(maxima.values()) <= RESIDUAL,
                  mass_imbalance=net, outward_mass_flux_kg_s_per_m_depth=mass,
                  outlet_reverse_mass_kg_s_per_m_depth=reverse,
                  plate_pressure_source="OpenFOAM zeroGradient pressure face equals adjacent-cell pressure")
    return final, previous, native, final["p"][0, :]


def profile_metrics(p):
    cp = p / Q_REF
    i = int(np.argmax(cp))
    return dict(centre_cp=float(cp[NX//2-1:NX//2+1].mean()), peak_cp=float(cp[i]),
                peak_x_over_W=float(XC[i]),
                centre_pressure_pa=float(p[NX//2-1:NX//2+1].mean()),
                integrated_plate_gauge_pressure_N_per_m_depth=float(p.sum() * DX))


def comparison(out, solutions, pressures):
    fig = plt.figure(figsize=(11.4, 8.2), facecolor="white")
    grid = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.3], width_ratios=[1, .029],
                           left=.082, right=.925, top=.838, bottom=.115, hspace=.49, wspace=.05)
    for i, (name, result) in enumerate(solutions.items()):
        ax = fig.add_subplot(grid[i, 0])
        im = speed_map(ax, result)
        ax.set_title(f"{chr(97+i)}) {name} speed and streamlines", loc="left", size=10)
    cbar = fig.colorbar(im, cax=fig.add_subplot(grid[:2, 1]))
    cbar.set_label(r"Speed, $|\mathbf{u}|/U_m$")
    ax = fig.add_subplot(grid[2, 0])
    for (name, p), color, style in zip(pressures.items(), ("#0b6881", "#a85329"), ("-", "--")):
        ax.plot(XC, p / Q_REF, color=color, lw=2., linestyle=style, label=name)
    ax.set(xlim=(-6, 6), xlabel="$x/W$", ylabel="$C_p$")
    ax.set_title("c) Pressure along the impingement plate", loc="left", size=10)
    ax.axhline(0, color="#99aab3", lw=.7)
    ax.grid(color="#dce7ec", linewidth=.7)
    ax.legend(loc="upper right", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(.082, .949, "A slot jet impinges on a flat plate", size=17, weight="bold")
    fig.text(.082, .901, "Re = 100  ·  H/W = 3  ·  14,400 cells per solver  ·  paired steady fields",
             color="#426077")
    fig.text(.082, .042, "The cyan line marks the inlet. Pressure is referenced to the two outlets.",
             size=9, color="#426077")
    fig.savefig(out / "impinging-jet-comparison.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.3), layout="constrained")
    for (name, p), color, style in zip(pressures.items(), ("#0b6881", "#a85329"), ("-", "--")):
        ax.plot(XC, p / Q_REF, color=color, lw=2., linestyle=style, label=name)
    ax.axvline(0, color="#a1aeb6", lw=.75)
    ax.axhline(0, color="#a1aeb6", lw=.75)
    ax.set(xlim=(-6, 6), xlabel="Distance from plate centre, $x/W$", ylabel="Pressure coefficient, $C_p$")
    ax.set_title("Plate pressure under the impinging jet", loc="left", weight="bold")
    ax.grid(color="#dce7ec", linewidth=.7)
    ax.legend(loc="upper right", frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(out / "impinging-jet-plate-pressure.png", dpi=180)
    fig.savefig(out / "impinging-jet-plate-pressure.svg")
    plt.close(fig)
    with (out / "plate-pressure.csv").open("x", newline="") as f:
        rows = csv.writer(f)
        rows.writerow(["x_over_slot_width", "fluent_plate_pressure_pa", "openfoam_plate_pressure_pa",
                       "fluent_cp", "openfoam_cp"])
        pf, po = pressures.values()
        rows.writerows(zip(XC, pf, po, pf / Q_REF, po / Q_REF))


def tracer_tracks(result, frames, substeps, dt):
    """Midpoint integration through bilinear interpolants of the steady cell fields."""
    vec = np.stack((result["u"], result["v"]), axis=-1) / UMEAN
    vec = np.pad(vec, ((1, 1), (1, 1), (0, 0)), mode="edge")
    vec[0] = -vec[1]  # Odd reflection gives zero velocity on the plate.
    # The initial tracers start one half-cell below the ceiling.
    vec[-1, abs(np.r_[XC[0] - DX/W, XC, XC[-1] + DX/W]) >= .5] *= -1

    def interp(pos):
        index = (pos + [6., 0.]) * W / DX + .5
        index = np.clip(index, 0., [NX + .999, NY + .999])
        i, j = np.floor(index).astype(int).T
        fx, fy = (index - np.column_stack((i, j))).T
        lower = (1 - fx[:, None]) * vec[j, i] + fx[:, None] * vec[j, i + 1]
        upper = (1 - fx[:, None]) * vec[j + 1, i] + fx[:, None] * vec[j + 1, i + 1]
        return (1 - fy[:, None]) * lower + fy[:, None] * upper

    jet_x = np.array([-.39, -.30, -.20, -.10, .10, .20, .30, .39])
    jet = np.column_stack((jet_x, np.full(len(jet_x), 2.973)))
    eddy = np.array([(side*x, y) for side in (-1, 1)
                     for x, y in ((2.35, 1.62), (3.05, 2.30), (3.92, 1.69))])
    pos = np.vstack((jet, eddy))
    track = np.empty((frames, len(pos), 2))
    active = np.ones(len(pos), dtype=bool)
    for k in range(frames):
        track[k] = np.where(active[:, None], pos, np.nan)
        for step in range(substeps):
            mid = pos + .5 * dt * interp(pos)
            pos = pos + dt * interp(mid)
            active &= (abs(pos[:, 0]) < 6) & (pos[:, 1] > 0) & (pos[:, 1] < 3)
            pos[~active] = [0., 0.]
    return track, len(jet)


def movie(out, solutions):
    count, substeps, dt = 77, 7, .03
    interval = substeps * dt
    fig = plt.figure(figsize=(11.6, 6.15), facecolor="white")
    grid = fig.add_gridspec(2, 2, width_ratios=[1, .030], left=.076, right=.912,
                           top=.80, bottom=.13, hspace=.33, wspace=.045)
    palette = ("#ffb1a5", "#e2fbff")
    plotted = []
    for i, (name, solution) in enumerate(solutions.items()):
        ax = fig.add_subplot(grid[i, 0])
        im = speed_map(ax, solution, streamlines=False)
        if i == 1:
            ax.set_xlabel("$x/W$")
            ax.tick_params(axis="x", labelbottom=True)
        ax.set_title(f"{chr(97+i)})  {name}", loc="left")
        track, jet_count = tracer_tracks(solution, count, substeps, dt)
        paths = []
        for color in palette:
            tail = LineCollection([], linewidths=1.45)
            ax.add_collection(tail)
            head = ax.scatter([], [], c=color, edgecolors="#ffffff", linewidths=.34, s=20, zorder=5)
            paths.append((tail, head, color))
        plotted.append((track, jet_count, paths))
    fig.colorbar(im, cax=fig.add_subplot(grid[:, 1])).set_label(r"Speed, $|\mathbf{u}|/U_m$")
    fig.text(.076, .947, "Follow the jet as it turns along the plate", weight="bold", size=16)
    fig.text(.076, .895, "Tracers advected through two saved steady velocity fields", color="#426077")
    clock = fig.text(.912, .896, "", ha="right", fontfamily="monospace")
    fig.text(.076, .041, "Peach: inlet tracers     |     White: recirculation tracers     |     Tracer time is scaled by W/Um",
             color="#426077", size=9)

    def draw(k):
        clock.set_text(f"t Um/W = {k * interval:4.1f}")
        for track, jet_count, objects in plotted:
            for points, (tail, head, color) in zip((track[:, :jet_count], track[:, jet_count:]), objects):
                segments, colors = [], []
                for old in range(max(0, k - 7), k):
                    lines = np.stack((points[old], points[old + 1]), axis=1)
                    valid = np.isfinite(lines).all(axis=(1, 2))
                    segments.extend(lines[valid])
                    colors.extend([(*to_rgb(color), .13 + .12 * (7 - k + old))] * int(valid.sum()))
                tail.set_segments(segments)
                tail.set_color(colors)
                heads = points[k]
                head.set_offsets(heads[np.isfinite(heads).all(axis=1)])

    indices = [0] * 6 + list(range(1, count)) + [count - 1] * 8
    anim = FuncAnimation(fig, draw, frames=indices, interval=85, blit=False)
    anim.save(out / "impinging-jet-tracers.gif", writer=PillowWriter(fps=12), dpi=90)
    plt.close(fig)


def main(args):
    sys.path.insert(0, str(args.support))
    with np.load(args.fluent / "geometry.npz") as geo:
        f, old_f, fn = fluent_fields(args.fluent, args.support, geo)
        fp = fluent_boundaries(args.fluent, fn)
        o, old_o, on, op = openfoam_fields(args.openfoam, geo)
    names = ("ANSYS Fluent", "OpenFOAM")
    fields = dict(zip(names, (f, o)))
    plate = dict(zip(names, (fp, op)))
    cases = {name: dict(native=native, measurements=measure(last, prev, native), plate=profile_metrics(p))
             for name, last, prev, native, p in zip(names, (f, o), (old_f, old_o), (fn, on), (fp, op))}
    payload = dict(cells_per_solver=NX * NY, Reynolds=100, gap_over_slot_width=H / W,
                   reference_speed_m_s=UMEAN, cases=cases, physical_validation=False,
                   animation_type="Passive tracers integrated through saved steady velocity fields")
    args.out.mkdir(parents=True)
    (args.out / "impinging-jet-summary.json").write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    comparison(args.out, fields, plate)
    movie(args.out, fields)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for key in ("fluent", "openfoam", "support", "out"):
        parser.add_argument("--" + key, type=Path, required=True)
    main(parser.parse_args())

