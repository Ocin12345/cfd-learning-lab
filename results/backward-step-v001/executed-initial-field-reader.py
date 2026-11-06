"""Plot saved backward-step native fields; never run a solver or invent a vortex.

Requires a fresh output folder. --solver fluent supports an explicitly labelled
early single-solver image; the default requires both final native solutions.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import re
import sys

os.environ.setdefault("MPLCONFIGDIR", "E:/CFD-Learning-Lab/temp/matplotlib-step")
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle
import numpy as np

TOL = 1e-9
NAMES = {"fluent": "ANSYS Fluent", "openfoam": "OpenFOAM"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hdf_indexed(group, count, width=1, partial=False):
    """Assemble saved ID-labelled blocks without assuming imported ordering."""
    values = np.full((count, width), np.nan)
    seen = np.zeros(count, dtype=bool)
    for dataset in group.values():
        require(isinstance(dataset, h5py.Dataset), "Unsupported nested HDF block")
        first = int(np.asarray(dataset.attrs["minId"]).ravel()[0]) - 1
        last = int(np.asarray(dataset.attrs["maxId"]).ravel()[0])
        block = np.asarray(dataset[:]).reshape(-1, width)
        require(0 <= first < last <= count and len(block) == last - first,
                "HDF block ID range or value count differs: " + dataset.name)
        require(not seen[first:last].any(), "Overlapping saved HDF IDs")
        values[first:last] = block
        seen[first:last] = True
    require(partial or seen.all(), "Saved HDF IDs do not cover the required range")
    return values[:, 0] if width == 1 else values


def read_fluent(run, setup):
    case, data = [run / ("fluent/step-final." + suffix) for suffix in ("cas.h5", "dat.h5")]
    count = setup["cells"]
    with h5py.File(case, "r") as handle:
        mesh = handle["meshes/1"]
        topology = mesh["cells/zoneTopology"]
        require(np.asarray(topology["minId"]).tolist() == [1]
                and np.asarray(topology["maxId"]).tolist() == [count],
                "Expected one contiguous saved fluid-cell zone")
        nodes = mesh["nodes/coords"]
        node_count = max(int(np.asarray(d.attrs["maxId"]).ravel()[0]) for d in nodes.values())
        width = next(iter(nodes.values())).shape[1]
        require(width in (2, 3), "Expected 2D node coordinates")
        xyz = hdf_indexed(nodes, node_count, width)
        require(np.isfinite(xyz).all() and (width == 2 or np.ptp(xyz[:, 2]) <= TOL),
                "Fluent nodes are nonfinite or not planar")
        xy = xyz[:, :2]
        blocks = mesh["faces/nodes"]
        faces = max(int(np.asarray(b.attrs["maxId"]).ravel()[0]) if "maxId" in b.attrs
                    else int(np.asarray(b["nnodes"].attrs["maxId"]).ravel()[0]) for b in blocks.values())
        edges = np.full((faces, 2), -1, dtype=np.int64)
        for block in blocks.values():
            sizes = block["nnodes"]
            first = int(np.asarray(sizes.attrs.get("minId", block.attrs.get("minId", [1]))).ravel()[0]) - 1
            values = np.asarray(block["nodes"][:], dtype=np.int64).reshape(-1, 2) - 1
            require(np.all(sizes[:] == 2) and np.all(edges[first:first + len(values)] == -1),
                    "Expected distinct two-node saved face blocks")
            edges[first:first + len(values)] = values
        require(np.all(edges >= 0) and np.all(edges < len(xy)), "Invalid saved Fluent edge vertices")
        owner = hdf_indexed(mesh["faces/c0"], faces).astype(np.int64) - 1
        other = hdf_indexed(mesh["faces/c1"], faces, partial=True)
        internal = np.isfinite(other)
        neighbour = other[internal].astype(np.int64) - 1
        require(np.all((owner >= 0) & (owner < count))
                and np.all((neighbour >= 0) & (neighbour < count)), "Invalid Fluent cell IDs")
        keys = np.unique(np.concatenate([
            np.repeat(owner, 2) * len(xy) + edges.ravel(),
            np.repeat(neighbour, 2) * len(xy) + edges[internal].ravel()]))
        cell_ids, sizes = np.unique(keys // len(xy), return_counts=True)
        require(np.array_equal(cell_ids, np.arange(count)) and np.all(sizes == 4),
                "Each Fluent cell must have exactly four distinct saved vertices")
        vertices = xy[(keys % len(xy)).reshape(count, 4)]
        centres = vertices.mean(axis=1)
        dx, dy = setup["dx_m"], setup["dy_m"]
        corners = np.sort(vertices - centres[:, None, :], axis=1)
        expected = np.array([[-dx / 2, -dy / 2], [-dx / 2, -dy / 2],
                             [dx / 2, dy / 2], [dx / 2, dy / 2]])
        require(np.max(np.abs(corners - expected)) <= TOL,
                "Saved Fluent cells are not the declared axis-aligned rectangles")
    with h5py.File(data, "r") as handle:
        fields = handle["results/1/phase-1/cells"]
        u, v, pressure = [hdf_indexed(fields[name], count) for name in ("SV_U", "SV_V", "SV_P")]
        for name, expected in (("SV_DENSITY", setup["density_kg_m3"]),
                               ("SV_MU_LAM", setup["dynamic_viscosity_Pa_s"])):
            if name in fields:
                require(np.allclose(hdf_indexed(fields[name], count), expected, rtol=1e-10, atol=1e-12),
                        "Saved Fluent material differs: " + name)
    return dict(centres=centres, u=u, v=v, p=pressure, inputs=[case, data],
                topology="Four distinct saved vertices per verified rectangular Fluent cell")


def foam_field(path, count, width, dimensions):
    raw = re.sub(r"/\*.*?\*/|//[^\n]*", "", path.read_text(encoding="utf-8"), flags=re.S)
    require(re.search(r"\bformat\s+ascii\s*;", raw), "OpenFOAM field must be ASCII: " + str(path))
    units = re.search(r"\bdimensions\s*\[([^]]+)\]\s*;", raw)
    require(units is not None and np.array_equal(np.fromstring(units[1], sep=" "), dimensions),
            "Unexpected OpenFOAM dimensions: " + str(path))
    kind = "vector" if width == 3 else "scalar"
    match = re.search(r"\binternalField\s+nonuniform\s+List<" + kind
                      + r">\s+(\d+)\s*\((.*?)\)\s*;", raw, re.S)
    if match:
        require(int(match[1]) == count, "Saved OpenFOAM cell count differs")
        values = np.fromstring(match[2].replace("(", " ").replace(")", " "), sep=" ")
    else:
        match = re.search(r"\binternalField\s+uniform\s+([^;]+);", raw)
        require(match is not None, "Unsupported OpenFOAM internalField")
        value = np.fromstring(match[1].replace("(", " ").replace(")", " "), sep=" ")
        require(len(value) == width, "Invalid uniform OpenFOAM field")
        values = np.tile(value, count)
    require(values.size == count * width and np.isfinite(values).all(),
            "Wrong-sized or nonfinite saved OpenFOAM field: " + str(path))
    return values.reshape(count, width) if width > 1 else values


def read_openfoam(run, setup):
    folder = run / "openfoam/2000"
    count = setup["cells"]
    centres = foam_field(folder / "C", count, 3, [0, 1, 0, 0, 0, 0, 0])
    velocity = foam_field(folder / "U", count, 3, [0, 1, -1, 0, 0, 0, 0])
    pressure = foam_field(folder / "p", count, 1, [0, 2, -2, 0, 0, 0, 0])
    require(np.ptp(centres[:, 2]) <= TOL and np.max(np.abs(velocity[:, 2])) <= 1e-12,
            "Expected one planar OpenFOAM cell layer and zero spanwise velocity")
    return dict(centres=centres[:, :2], u=velocity[:, 0], v=velocity[:, 1],
                p=pressure * setup["density_kg_m3"], inputs=[folder / p for p in ("C", "U", "p")],
                topology="Native writeCellCentres C field; no imported-order assumption")


def validate_setup(setup):
    expected = dict(nx=240, ny=40, cells=8800, x_min_m=-.4, x_max_m=2.,
                    height_m=.2, step_height_m=.1, density_kg_m3=1.,
                    dynamic_viscosity_Pa_s=.001, inlet_velocity_m_s=1.)
    for key, value in expected.items():
        require(key in setup and setup[key] == value, "Unsupported teaching-case setting: " + key)
    setup = dict(setup)
    for key, value in (("dx_m", .01), ("dy_m", .005)):
        require(key not in setup or abs(setup[key] - value) <= 1e-12, "Unexpected spacing: " + key)
        setup[key] = value
    return setup


def grid_coordinates(setup):
    xe = setup["x_min_m"] + np.arange(setup["nx"] + 1) * setup["dx_m"]
    ye = np.arange(setup["ny"] + 1) * setup["dy_m"]
    return xe, ye, (xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2


def map_and_check(field, setup):
    nx, ny, count = [setup[key] for key in ("nx", "ny", "cells")]
    xe, ye, x, y = grid_coordinates(setup)
    centres = field["centres"]
    require(centres.shape == (count, 2) and np.isfinite(centres).all(), "Invalid native centre array")
    require(all(np.asarray(field[k]).shape == (count,) and np.isfinite(field[k]).all()
                for k in ("u", "v", "p")), "Wrong-sized or nonfinite native field")
    ij = np.rint((centres - [setup["x_min_m"], 0]) / [setup["dx_m"], setup["dy_m"]] - .5).astype(np.int64)
    require(np.all((ij >= 0) & (ij < [nx, ny])), "Native centre outside declared grid")
    error = np.max(np.abs(centres - np.column_stack([x[ij[:, 0]], y[ij[:, 1]]])))
    solid = (x[None, :] < 0) & (y[:, None] < setup["step_height_m"])
    fluid_slots = np.flatnonzero(~solid)
    slots = ij[:, 1] * nx + ij[:, 0]
    require(len(fluid_slots) == count and error <= TOL and np.array_equal(np.sort(slots), fluid_slots),
            "Native centres must map bijectively to exactly the fluid slots within 1e-9 m")
    field["ij"] = ij
    field["grid"] = {}
    for key in ("u", "v", "p"):
        values = np.full(nx * ny, np.nan)
        values[slots] = field[key]
        field["grid"][key] = np.ma.array(values.reshape(ny, nx), mask=solid)
    field["grid"]["speed"] = np.hypot(field["grid"]["u"], field["grid"]["v"])
    u = field["grid"]["u"]
    minimum = int(np.argmin(field["u"]))
    negative = np.ma.filled(u < 0, False)
    first = u[0, x > 0].compressed()
    xf = x[x > 0]
    first_negative = np.flatnonzero(first < 0)
    proxy = dict(status="no negative first-cell axial velocity", first_cell_y_m=float(y[0]),
                 first_cell_y_over_h=float(y[0] / setup["step_height_m"]),
                 negative_first_cell_count=int(len(first_negative)), x_m=None, x_over_h=None)
    if len(first_negative):
        start = int(first_negative[0])
        crossings = np.flatnonzero((first[start:-1] < 0) & (first[start + 1:] >= 0)) + start
        proxy["first_negative_x_m"] = float(xf[start])
        if len(crossings):
            index = int(crossings[0])
            fraction = -first[index] / (first[index + 1] - first[index])
            endpoint = float(xf[index] + fraction * (xf[index + 1] - xf[index]))
            proxy.update(status="first negative-to-nonnegative crossing after first negative cell",
                         x_m=endpoint, x_over_h=endpoint / setup["step_height_m"],
                         bracket_x_m=[float(xf[index]), float(xf[index + 1])],
                         bracket_u_m_s=[float(first[index]), float(first[index + 1])])
        else:
            proxy["status"] = "negative pocket has no downstream sign-crossing before outlet"
    return dict(cells=count, native_centres_bijective=True, solid_slots_masked=int(solid.sum()),
                centre_tolerance_m=TOL, maximum_centre_mapping_error_m=float(error),
                finite_native_fields=True, topology=field["topology"],
                min_u_m_s=float(field["u"].min()), max_u_m_s=float(field["u"].max()),
                min_u_native_centre_m=field["centres"][minimum].tolist(),
                negative_u_cell_count=int(negative.sum()),
                downstream_below_step_negative_cell_count=int((negative & (x[None, :] > 0)
                                                             & (y[:, None] < setup["step_height_m"])).sum()),
                speed_range_m_s=[float(field["grid"]["speed"].min()), float(field["grid"]["speed"].max())],
                pressure_range_Pa=[float(field["p"].min()), float(field["p"].max())],
                near_bottom_first_cell_zero_u_proxy=proxy)


def reverse_seeds(field, setup, limit=9):
    """Distribute real native-centre seeds inside the saved negative-u pocket."""
    centres = field["centres"]
    select = (field["u"] < 0) & (centres[:, 0] > 0) & (centres[:, 1] < setup["step_height_m"])
    candidates = centres[select] / setup["step_height_m"]
    if not len(candidates):
        return np.empty((0, 2))
    selected = [int(np.argmin(field["u"][select]))]
    distances = np.sum((candidates - candidates[selected[0]]) ** 2, axis=1)
    for _ in range(min(limit, len(candidates)) - 1):
        index = int(np.argmax(distances))
        if distances[index] < .12 ** 2:
            break
        selected.append(index)
        distances = np.minimum(distances, np.sum((candidates - candidates[index]) ** 2, axis=1))
    return candidates[selected]


def add_flow(ax, field, setup, close_up=False):
    _, _, x, y = grid_coordinates(setup)
    h = setup["step_height_m"]
    grid = field["grid"]
    ax.streamplot(x / h, y / h, grid["u"], grid["v"], density=(2.4, .9),
                  color="#23323a99", linewidth=.55, arrowsize=.8, maxlength=8,
                  broken_streamlines=True, zorder=3)
    seeds = reverse_seeds(field, setup)
    if len(seeds):
        ax.streamplot(x / h, y / h, grid["u"], grid["v"], start_points=seeds,
                      density=(3.2, 1.4), color="#0c2437", linewidth=.9 if close_up else .7,
                      arrowsize=1.1 if close_up else .9, minlength=.015, maxlength=8,
                      broken_streamlines=False, zorder=4)
    if float(grid["u"].min()) < 0 < float(grid["u"].max()):
        ax.contour(x / h, y / h, grid["u"], levels=[0], colors=["#202328"],
                   linewidths=1.1, linestyles="--", zorder=5)
    return len(seeds)


def decorate(ax, setup, xlim, zoom=False):
    h = setup["step_height_m"]
    ax.add_patch(Rectangle((setup["x_min_m"] / h, 0), -setup["x_min_m"] / h, 1,
                           facecolor="#d4d8dc", edgecolor="none", zorder=6))
    ax.plot([setup["x_min_m"] / h, 0, 0, setup["x_max_m"] / h], [1, 1, 0, 0],
            color="#2e343b", lw=1.0, zorder=7)
    ax.axhline(setup["height_m"] / h, color="#2e343b", lw=1.0, zorder=7)
    if not zoom:
        ax.text(-2, .43, "solid step", color="#57616d", ha="center", va="center", fontsize=10, zorder=8)
    ax.set(xlim=xlim, ylim=(0, 1.35 if zoom else 2), xlabel="$x/h$", ylabel="$y/h$")
    ax.set_aspect("equal", adjustable="box")
    ax.set_yticks([0, .5, 1] if zoom else [0, 1, 2])
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=10)


def save_figure(fig, out, name):
    fig.canvas.draw()
    fig.set_layout_engine("none")
    for extension in ("png", "svg"):
        fig.savefig(out / (name + "." + extension), dpi=300, facecolor="white")
    plt.close(fig)


def plot_fields(fields, checks, setup, out):
    h = setup["step_height_m"]
    xe, ye, x, y = grid_coordinates(setup)
    lower = min(float(f["u"].min()) for f in fields.values())
    upper = max(float(f["u"].max()) for f in fields.values())
    require(upper > max(lower, 0), "Saved axial velocity has no positive range")
    negative = lower < 0
    norm = TwoSlopeNorm(vmin=lower, vcenter=0, vmax=upper) if negative else Normalize(vmin=0, vmax=upper)
    cmap = "RdBu_r" if negative else "YlOrRd"
    ticks = [lower, lower / 2, 0, upper / 3, 2 * upper / 3, upper] if negative else [0, upper / 2, upper]
    reverse_extents = [float(f["centres"][(f["u"] < 0) & (f["centres"][:, 0] > 0)
                                       & (f["centres"][:, 1] < h), 0].max() / h)
                       for f in fields.values() if np.any((f["u"] < 0) & (f["centres"][:, 0] > 0)
                                                        & (f["centres"][:, 1] < h))]
    zoom_right = min(20., max(4., max(reverse_extents, default=5.) + .9))
    seed_counts = {}
    for zoom, stem in ((False, "step-computed-flow"), (True, "step-recirculation-closeup")):
        fig, axes = plt.subplots(len(fields), 1, figsize=(14, 2.3 + len(fields) * (2.65 if zoom else 1.85)),
                                 squeeze=False, layout="constrained")
        for row, (name, field) in enumerate(fields.items()):
            ax = axes[row, 0]
            art = ax.pcolormesh(xe / h, ye / h, field["grid"]["u"], shading="flat", cmap=cmap,
                                norm=norm, rasterized=True, zorder=1)
            seeds = add_flow(ax, field, setup, close_up=zoom)
            seed_counts[name] = seeds
            decorate(ax, setup, (-.65, zoom_right) if zoom else (-4, 20), zoom=zoom)
            result = checks[name]
            ax.set_title(f"{NAMES[name]}  |  minimum $u_x$ = {result['min_u_m_s']:+.4f} m/s"
                         f"  |  {result['negative_u_cell_count']:,} cells with $u_x<0$",
                         loc="left", fontsize=11, pad=9)
            if zoom:
                proxy = result["near_bottom_first_cell_zero_u_proxy"]
                if proxy["x_over_h"] is not None:
                    ax.plot(proxy["x_over_h"], proxy["first_cell_y_over_h"], "D", ms=5,
                            color="#151b22", markeredgecolor="white", markeredgewidth=.5, zorder=9)
                    ax.annotate(f"first-cell zero-$u_x$ proxy\n$x/h$ = {proxy['x_over_h']:.2f}",
                                xy=(proxy["x_over_h"], proxy["first_cell_y_over_h"]),
                                xytext=(-10, 36), textcoords="offset points", ha="right", fontsize=8,
                                bbox=dict(boxstyle="round,pad=.3", fc="white", ec="#b9c2ca", alpha=.94),
                                arrowprops=dict(arrowstyle="-", color="#505a64", lw=.7), zorder=10)
                if not seeds:
                    ax.text(.97, .1, "No downstream negative-$u_x$ cells in saved field",
                            transform=ax.transAxes, ha="right", fontsize=9)
        cb = fig.colorbar(art, ax=list(axes[:, 0]), location="bottom", fraction=.08, pad=.08, aspect=50,
                          ticks=ticks, format="%.3g")
        cb.set_label("Axial velocity $u_x$ (m/s)   |   blue: reverse flow   |   red: forward flow" if negative
                     else "Axial velocity $u_x$ (m/s)")
        title = "Reverse flow behind the step" if zoom else "Separated laminar flow over a backward-facing step"
        qualifier = "Single-solver saved result" if len(fields) == 1 else "Two native solutions on one colour scale"
        fig.suptitle(title + "\n" + qualifier + "  |  $Re_h=100$  |  8,800 cells per solver  |  $h=0.1$ m",
                     fontsize=16, fontweight="medium")
        fig.supxlabel("Dashed line: $u_x=0$. Arrows follow the saved vector field. "
                      + ("Diamond: grid-dependent first-cell proxy; not exact wall reattachment." if zoom
                         else "Uniform 1 m/s inlet; stationary no-slip walls; no mesh-independence claim."),
                      fontsize=9, color="#4b5661")
        save_figure(fig, out, stem)
    return dict(axial_velocity_m_s=[min(lower, 0.), upper], zero_centred=negative,
                normalization="Separate linear negative/positive intervals meeting at zero" if negative
                              else "Linear from zero to saved maximum",
                recirculation_zoom_x_over_h=[-.65, zoom_right], reverse_flow_seed_counts=seed_counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True, help="Native run root containing mesh/, fluent/, openfoam/")
    parser.add_argument("--out", type=Path, required=True, help="Fresh output directory; existing paths are refused")
    parser.add_argument("--solver", choices=("both", "fluent", "openfoam"), default="both")
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    require(not out.exists(), "Choose a fresh export directory")
    setup_path = run / "mesh/setup.json"
    setup = validate_setup(json.loads(setup_path.read_text(encoding="utf-8-sig")))
    selected = ("fluent", "openfoam") if args.solver == "both" else (args.solver,)
    required = {"fluent": [run / ("fluent/step-final." + suffix) for suffix in ("cas.h5", "dat.h5")],
                "openfoam": [run / "openfoam/2000" / name for name in ("C", "U", "p")]}
    missing = [str(path) for name in selected for path in required[name] if not path.is_file()]
    require(not missing, "Native inputs are not ready: " + ", ".join(missing))
    fields = {name: (read_fluent if name == "fluent" else read_openfoam)(run, setup) for name in selected}
    checks = {name: map_and_check(field, setup) for name, field in fields.items()}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.linewidth": .65,
                         "svg.fonttype": "path", "savefig.bbox": "tight", "figure.dpi": 120})
    out.mkdir(parents=True)
    colour_limits = plot_fields(fields, checks, setup, out)
    for name, field in fields.items():
        with (out / (name + "-actual-cell-fields.csv")).open("x", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["native_cell_id_1based", "ix_0based", "iy_0based", "x_m", "y_m", "u_m_s", "v_m_s",
                             "speed_m_s", "gauge_pressure_Pa"])
            for i in range(setup["cells"]):
                writer.writerow([i + 1, *field["ij"][i], *field["centres"][i], field["u"][i], field["v"][i],
                                 np.hypot(field["u"][i], field["v"][i]), field["p"][i]])
    record = dict(solvers=checks, solver_selection=args.solver, colour_limits=colour_limits,
                  near_bottom_proxy_definition="Linear u_x zero crossing between consecutive downstream first-row cell centres, "
                  "after the first negative cell; grid dependent and not exact wall reattachment",
                  streamline_method="Matplotlib integration through saved mapped cell vectors; additional seeds are actual "
                  "negative-u native centres below the step; no invented vortex or appended wall values",
                  residual_mass_and_stability_checks="Separate native-run evidence; not assessed by this plotter",
                  analytical_poiseuille_checks_applied=False, mesh_independence_established=False,
                  physical_validation=False, scope="Bounded educational visualization of saved backward-step fields")
    if len(fields) == 2:
        diff = fields["fluent"]["grid"]["u"] - fields["openfoam"]["grid"]["u"]
        record["cross_solver_axial_difference_m_s"] = dict(max_abs=float(np.max(np.abs(diff))),
                                                          rms=float(np.sqrt(np.ma.mean(diff ** 2))))
    (out / "checks.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    script = Path(__file__).resolve()
    (out / "executed-plotter.py").write_bytes(script.read_bytes())
    inputs = [setup_path, *(p for field in fields.values() for p in field["inputs"])]
    manifest = dict(plotter_sha256=sha256(script), command_arguments=sys.argv[1:], solver_selection=args.solver,
                    inputs=[dict(path=str(path.resolve()), sha256=sha256(path)) for path in inputs])
    (out / "input-hash-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    measurements = "\n".join(f"- {NAMES[name]}: minimum axial velocity {check['min_u_m_s']:+.6g} m/s; "
                              f"{check['negative_u_cell_count']} negative-u cells; first-cell crossing proxy "
                              + (f"x/h={check['near_bottom_first_cell_zero_u_proxy']['x_over_h']:.4g}."
                                 if check['near_bottom_first_cell_zero_u_proxy']['x_over_h'] is not None
                                 else check['near_bottom_first_cell_zero_u_proxy']['status'] + ".")
                              for name, check in checks.items())
    caption = ("# Saved backward-step flow fields\n\n"
        f"Solvers shown: {', '.join(NAMES[n] for n in selected)}. "
        "Step h=0.1 m; upstream x=-0.4..0 m and y=0.1..0.2 m; downstream x=0..2 m and y=0..0.2 m. "
        "Uniform inlet 1 m/s; rho=1 kg/m3; mu=0.001 Pa.s; Re_h=100. Stationary no-slip walls; "
        "zero-gauge-pressure outlet; 8,800 quadrilateral fluid cells per code.\n\n"
        "step-computed-flow.png/svg shows flat per-cell axial velocity and integrated streamline arrows. "
        "step-recirculation-closeup.png/svg magnifies the downstream lower-channel region. Both codes and both "
        "figures share actual-data colour limits. When reverse flow exists, zero is centred in the diverging "
        "colour map; the negative and positive halves each use a separate linear interval to resolve weak "
        "negative velocities. Thus equally saturated blue and red need not have equal velocity magnitudes. "
        "A dashed contour marks u_x=0 only where that contour exists. Grey is the masked solid step. "
        "Axes x/h,y/h use equal geometric scale. Additional streamline seeds are actual cell centres where "
        "u_x<0 downstream below the step; curves integrate the native mapped vector field and introduce no "
        "synthetic vortex. No boundary velocity samples are appended.\n\n"
        "The diamond, where present, is the linear negative-to-nonnegative crossing in the first row of "
        "cell-centred axial velocities at y=0.0025 m, or y/h=0.025. This is a grid-dependent near-bottom "
        "proxy, not exact wall reattachment, and is not benchmark or physical validation.\n\n"
        + measurements + "\n\n"
        "Fluent centres are reconstructed from four unique saved face-connected vertices per verified "
        "rectangle; OpenFOAM centres come from its saved C field. Both must map bijectively to exactly the "
        "8,800 fluid slots within 1e-9 m; the 800 solid slots are masked. Imported cell ordering is never "
        "assumed. OpenFOAM kinematic pressure is multiplied by rho before CSV export. Input hashes, executed "
        "plotter, cell data and mapping diagnostics accompany the figures.\n\n"
        "Residual, mass and iteration-stability checks are separate native-run evidence. No Poiseuille "
        "analytical target applies to this geometry; no mesh-independence or physical-validation claim is "
        "made. Fluent/OpenFOAM calculated the fields; Codex wrote the workflow, reconstructed the saved data "
        "and prepared these figures. Owner reflection has not been provided.\n")
    (out / "CAPTIONS.md").write_text(caption, encoding="utf-8")
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
