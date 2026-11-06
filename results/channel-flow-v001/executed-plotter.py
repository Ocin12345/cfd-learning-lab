"""Plot saved channel-flow fields; never start or interpolate a CFD solution.

Default input is both final native solutions. --solver fluent permits an early,
explicitly single-solver export. Structural errors stop before creating output;
failed analytical targets are retained, labelled, and return exit status 2.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re
import sys

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


TOL = 1e-9
STATIONS = (0.02, 0.1, 0.3, 0.8)
COLOURS = {"fluent": "#007f86", "openfoam": "#d06b25"}
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
    """Assemble saved ID-labelled blocks; block or imported ordering is immaterial."""
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
    case = run / "fluent/channel-final.cas.h5"
    data = run / "fluent/channel-final.dat.h5"
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
        face_blocks = mesh["faces/nodes"]
        faces = max(int(np.asarray(b.attrs["maxId"]).ravel()[0])
                    if "maxId" in b.attrs else int(np.asarray(b["nnodes"].attrs["maxId"]).ravel()[0])
                    for b in face_blocks.values())
        edges = np.full((faces, 2), -1, dtype=np.int64)
        for block in face_blocks.values():
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
        # Arithmetic vertex centres are geometric centres only for these verified rectangles.
        dx, dy = setup["length_m"] / setup["nx"], setup["height_m"] / setup["ny"]
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
    folder = run / "openfoam/800"
    count = setup["cells"]
    centres = foam_field(folder / "C", count, 3, [0, 1, 0, 0, 0, 0, 0])
    velocity = foam_field(folder / "U", count, 3, [0, 1, -1, 0, 0, 0, 0])
    pressure = foam_field(folder / "p", count, 1, [0, 2, -2, 0, 0, 0, 0])
    require(np.ptp(centres[:, 2]) <= TOL and np.max(np.abs(velocity[:, 2])) <= 1e-12,
            "Expected one planar OpenFOAM cell layer and zero spanwise velocity")
    return dict(centres=centres[:, :2], u=velocity[:, 0], v=velocity[:, 1],
                p=pressure * setup["density_kg_m3"], inputs=[folder / p for p in ("C", "U", "p")],
                topology="Native writeCellCentres C field; no imported-order assumption")


def map_and_check(field, setup):
    nx, ny = setup["nx"], setup["ny"]
    length, height, mean = [setup[key] for key in ("length_m", "height_m", "inlet_velocity_m_s")]
    centres = field["centres"]
    require(centres.shape == (nx * ny, 2) and np.isfinite(centres).all(), "Invalid native centre array")
    require(all(np.asarray(field[k]).shape == (nx * ny,) and np.isfinite(field[k]).all()
                for k in ("u", "v", "p")), "Wrong-sized or nonfinite native field")
    ij = np.rint(centres / [length / nx, height / ny] - 0.5).astype(np.int64)
    require(np.all((ij >= 0) & (ij < [nx, ny])), "Native cell centre outside declared channel")
    error = np.max(np.abs(centres - (ij + 0.5) * [length / nx, height / ny]))
    slots = ij[:, 1] * nx + ij[:, 0]
    require(error <= TOL and np.array_equal(np.sort(slots), np.arange(nx * ny)),
            "Native centres do not match the structured grid bijectively within 1e-9 m")
    order = np.argsort(slots)
    field["native_order"] = order
    field["ij"] = ij
    field["grid"] = {k: field[k][order].reshape(ny, nx) for k in ("u", "v", "p")}
    field["grid"]["speed"] = np.hypot(field["grid"]["u"], field["grid"]["v"])
    speed = field["grid"]["speed"]
    require(speed.min() >= 0 and speed.max() <= 1.6, "Actual speed is outside fixed 0..1.6 m/s colour limits")
    x, y = (np.arange(nx) + .5) * length / nx, (np.arange(ny) + .5) * height / ny
    analytic = 6 * (y / height) * (1 - y / height)
    downstream = int(np.argmin(np.abs(x - .8)))
    profile = field["grid"]["u"][:, downstream] / mean
    rms = float(np.sqrt(np.mean((profile - analytic) ** 2)))
    centreline = float(np.interp(height / 2, y, field["grid"]["u"][:, downstream]))
    developed = (x >= .6) & (x <= .9)
    slope, intercept = np.polyfit(x[developed], field["grid"]["p"].mean(axis=0)[developed], 1)
    target = -12 * setup["dynamic_viscosity_Pa_s"] * mean / height ** 2
    gradient_error = float(abs((slope - target) / target))
    tests = dict(downstream_profile_rms_le_2percent_Umean=rms <= .02,
                 centreline_within_2percent_of_1p5U=abs(centreline / (1.5 * mean) - 1) <= .02,
                 developed_pressure_gradient_error_le_3percent=gradient_error <= .03)
    return dict(cells=nx * ny, native_centres_bijective=True, centre_tolerance_m=TOL,
                maximum_centre_mapping_error_m=float(error), finite_native_fields=True,
                topology=field["topology"], speed_range_m_s=[float(speed.min()), float(speed.max())],
                pressure_range_Pa=[float(field["p"].min()), float(field["p"].max())],
                profile_x_m=float(x[downstream]), profile_rms_error_over_Umean=rms,
                centreline_interpolated_u_m_s=centreline, centreline_target_m_s=1.5 * mean,
                pressure_fit_interval_m=[.6, .9], pressure_fit_intercept_Pa=float(intercept),
                developed_dp_dx_Pa_m=float(slope), analytical_dp_dx_Pa_m=float(target),
                relative_pressure_gradient_error=gradient_error, check_results=tests,
                analytical_targets_passed=all(tests.values()))


def save_figure(fig, out, name):
    fig.canvas.draw()
    fig.set_layout_engine("none")
    for extension in ("png", "svg"):
        fig.savefig(out / (name + "." + extension), dpi=300, facecolor="white")
    plt.close(fig)


def plot_fields(fields, setup, out, passed):
    names = list(fields)
    fig, axes = plt.subplots(2, len(names), figsize=(7.3 * len(names), 6.0), squeeze=False,
                             layout="constrained")
    xe = np.linspace(0, setup["length_m"], setup["nx"] + 1)
    ye = np.linspace(0, setup["height_m"], setup["ny"] + 1)
    x, y = (xe[1:] + xe[:-1]) / 2, (ye[1:] + ye[:-1]) / 2
    lower = min(0., *(float(f["p"].min()) for f in fields.values()))
    upper = max(float(f["p"].max()) for f in fields.values())
    require(upper > lower, "Pressure range has no extent")
    pad = .01 * (upper - lower)
    pressure_limits = (lower - pad if lower < 0 else 0., upper + pad)
    for row, (quantity, limits, cmap, unit) in enumerate((
            ("speed", (0., 1.6), "viridis", r"Speed $|\mathbf{u}|$ (m s$^{-1}$)"),
            ("p", pressure_limits, "cividis", "Gauge pressure (Pa)"))):
        for col, name in enumerate(names):
            ax, field = axes[row, col], fields[name]["grid"]
            art = ax.pcolormesh(xe, ye, field[quantity], shading="flat", cmap=cmap,
                                vmin=limits[0], vmax=limits[1], rasterized=True)
            if quantity == "speed":
                ax.streamplot(x, y, field["u"], field["v"], color="#ffffffb0", linewidth=.65,
                              arrowsize=.65, density=(.65, .4), broken_streamlines=False,
                              start_points=np.column_stack([np.full(5, x[0] + 1e-8),
                                                            np.linspace(.015, .085, 5)]))
            letter = chr(ord("a") + row * len(names) + col)
            ax.set_title(f"({letter})  {NAMES[name]}  |  {'Velocity and streamlines' if row == 0 else 'Pressure'}",
                         loc="left", fontsize=11, pad=9)
            ax.set(xlim=(0, setup["length_m"]), ylim=(0, setup["height_m"]), xlabel="$x$ (m)", ylabel="$y$ (m)")
            ax.set_yticks([0, .05, .1])
            ax.tick_params(labelsize=9)
            ax.spines[["top", "right"]].set_visible(False)
        cb = fig.colorbar(art, ax=list(axes[row]), fraction=.025, pad=.02, aspect=18)
        cb.set_label(unit, fontsize=10)
    fig.suptitle("Laminar flow between parallel plates\n"
                 "Uniform inlet 1 m/s  ·  height 0.1 m  ·  Re$_H$ = 10  ·  7,680 cells per solver",
                 fontsize=15, fontweight="medium")
    label = "Educational analytical checks; no mesh study or experimental validation. Vertical scale expanded for readability."
    if not passed:
        label = "ANALYTICAL TARGET(S) NOT MET — diagnostic saved fields. " + label
    fig.supxlabel(label, fontsize=8, color="#833a29" if not passed else "#555555")
    save_figure(fig, out, "channel-computed-fields")
    return {"speed_m_s": [0., 1.6], "gauge_pressure_Pa": list(pressure_limits)}


def plot_profiles(fields, setup, out, passed):
    fig, axes = plt.subplots(1, 4, figsize=(12.3, 4.4), sharey=True, layout="constrained")
    x = (np.arange(setup["nx"]) + .5) * setup["length_m"] / setup["nx"]
    eta = (np.arange(setup["ny"]) + .5) / setup["ny"]
    fine = np.linspace(0, 1, 301)
    rows = []
    for ax, station in zip(axes, STATIONS):
        index = int(np.argmin(np.abs(x - station)))
        for name, field in fields.items():
            u = field["grid"]["u"][:, index] / setup["inlet_velocity_m_s"]
            ax.plot(u, eta, color=COLOURS[name], lw=2.0, ls="-" if name == "fluent" else "--", label=NAMES[name])
            rows.extend((name, station, x[index], yy, uv, 6 * yy * (1 - yy)) for yy, uv in zip(eta, u))
        ax.plot(6 * fine * (1 - fine), fine, color="#626974", lw=1.3, ls=":", label="Fully developed planar reference")
        ax.set(xlim=(0, 1.65), ylim=(0, 1), xlabel="$u/U_{mean}$", title=f"$x$ = {x[index]:.3f} m")
        ax.set_xticks([0, .5, 1, 1.5])
        ax.grid(color="#e1e4e7", lw=.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("$y/H$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=len(handles), fontsize=10, frameon=False)
    fig.suptitle("How a uniform inlet becomes a parabola\n"
                 + ("Actual saved velocity profiles; analytical target(s) not met" if not passed
                    else "Actual saved velocity profiles compared with  $u/U_{mean}=6(y/H)(1-y/H)$"), fontsize=14)
    save_figure(fig, out, "channel-velocity-profiles")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True, help="Native run root containing mesh/, fluent/, openfoam/")
    parser.add_argument("--out", type=Path, required=True, help="Fresh output directory; existing paths are refused")
    parser.add_argument("--solver", choices=("both", "fluent", "openfoam"), default="both")
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    require(not out.exists(), "Choose a fresh export directory")
    setup_path = run / "mesh/setup.json"
    setup = json.loads(setup_path.read_text(encoding="utf-8"))
    for key, value in dict(nx=240, ny=32, cells=7680, length_m=1., height_m=.1,
                           density_kg_m3=1., dynamic_viscosity_Pa_s=.01,
                           kinematic_viscosity_m2_s=.01, inlet_velocity_m_s=1.).items():
        require(setup[key] == value, "Unsupported teaching-case setting: " + key)
    selected = ("fluent", "openfoam") if args.solver == "both" else (args.solver,)
    # Check all required files first, so a partial native run cannot produce a paired figure.
    required = {"fluent": [run / ("fluent/channel-final." + suffix) for suffix in ("cas.h5", "dat.h5")],
                "openfoam": [run / "openfoam/800" / name for name in ("C", "U", "p")]}
    missing = [str(path) for name in selected for path in required[name] if not path.is_file()]
    require(not missing, "Native inputs are not ready: " + ", ".join(missing))
    fields = {name: (read_fluent if name == "fluent" else read_openfoam)(run, setup) for name in selected}
    checks = {name: map_and_check(field, setup) for name, field in fields.items()}
    passed = all(item["analytical_targets_passed"] for item in checks.values())
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.linewidth": .65,
                         "svg.fonttype": "path", "savefig.bbox": "tight", "figure.dpi": 120})
    out.mkdir(parents=True)
    colour_limits = plot_fields(fields, setup, out, passed)
    profiles = plot_profiles(fields, setup, out, passed)
    for name, field in fields.items():
        with (out / (name + "-actual-cell-fields.csv")).open("x", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["native_cell_id_1based", "ix_0based", "iy_0based", "x_m", "y_m", "u_m_s", "v_m_s", "speed_m_s", "gauge_pressure_Pa"])
            for i in range(setup["cells"]):
                writer.writerow([i + 1, *field["ij"][i], *field["centres"][i], field["u"][i], field["v"][i],
                                 np.hypot(field["u"][i], field["v"][i]), field["p"][i]])
    with (out / "profile-comparison.csv").open("x", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["solver", "requested_x_m", "actual_column_x_m", "y_over_H", "actual_u_over_Umean", "analytical_fully_developed_u_over_Umean"])
        writer.writerows(profiles)
    record = dict(solvers=checks, analytical_targets_passed=passed, colour_limits=colour_limits,
                  pressure_conventions={"fluent": "Gauge pressure in Pa", "openfoam": "Kinematic p multiplied by density to obtain Pa"},
                  profile_rms_definition="sqrt(mean((u/Umean - 6*eta*(1-eta))**2)) at nearest native column to x=.8 m",
                  centreline_definition="Linear interpolation of the two central saved-cell values at x nearest .8 m",
                  pressure_gradient_definition="Least-squares slope of cross-channel cell-mean gauge pressure for .6 <= x <= .9 m",
                  residual_and_mass_checks="Separate native-run checks; not assessed by this plotter",
                  mesh_independence_established=False, experimental_validation=False,
                  scope="Educational comparison with the analytical fully developed planar-channel solution")
    (out / "checks.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    script = Path(__file__).resolve()
    (out / "executed-plotter.py").write_bytes(script.read_bytes())
    inputs = [setup_path, *(p for field in fields.values() for p in field["inputs"])]
    manifest = dict(plotter_sha256=sha256(script), solver_selection=args.solver,
                    inputs=[dict(path=str(path.resolve()), sha256=sha256(path)) for path in inputs])
    (out / "input-hash-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    caption = (
        "# Saved channel-flow fields\n\n"
        f"Solvers shown: {', '.join(NAMES[n] for n in selected)}. Uniform inlet 1 m/s; L = 1 m, H = 0.1 m; "
        "rho = 1 kg/m3, mu = 0.01 Pa.s; Re_H = 10; stationary no-slip plates and zero-gauge-pressure outlet. "
        "Each solver has 240 x 32 = 7,680 cells.\n\n"
        "channel-computed-fields.png/svg uses actual native cell-centred values with flat per-cell colours. "
        "Velocity uses a common 0..1.6 m/s scale; pressure limits include every shown value and are shared. "
        "The vertical display scale is expanded; axis labels retain metres. White streamline curves and arrows "
        "use Matplotlib integration through the mapped saved vector field, without changing contour values. "
        "Fluent centres are reconstructed from four unique saved face-connected vertices per rectangular cell; "
        "OpenFOAM centres come from its saved C field. Both centre sets must map bijectively within 1e-9 m. "
        "OpenFOAM kinematic pressure is multiplied by rho = 1 kg/m3.\n\n"
        "channel-velocity-profiles.png/svg shows nearest saved columns to x = 0.02, 0.1, 0.3, 0.8 m; exact column "
        "positions are printed in the panels and CSV. The dotted parabola is an analytical fully developed "
        "reference, not simulated data. Curves connect actual cell-centred samples and do not append wall values. "
        "The near-inlet profiles are developing and are not required to match that downstream reference.\n\n"
        f"Analytical targets passed: {str(passed).lower()}. See checks.json for measured profile and pressure-gradient "
        "checks. Residual, mass-balance and iteration-stability evidence is checked separately. This is a bounded "
        "educational analytical check, not mesh independence or experimental validation. Fluent/OpenFOAM calculated "
        "the fields; Codex reconstructed their saved data, checked the mappings and prepared the figures.\n")
    (out / "CAPTIONS.md").write_text(caption, encoding="utf-8")
    print(json.dumps(record, indent=2, allow_nan=False))
    return 0 if passed else 2


if __name__ == "__main__":
    sys.exit(main())
