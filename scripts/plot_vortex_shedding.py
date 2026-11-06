"""Render actual saved transient cylinder fields; never run CFD or invent motion.

Fluent uses a shared native case plus snapshots.json [{time,step,data}, ...].
OpenFOAM uses reconstructed ASCII time/U and time/p and native polyMesh files.
Every frame pairs equal physical times. Saved face/cell topology supplies every
filled triangle; no Delaunay triangulation, synthetic vortex or temporal
interpolation is used. All writes are restricted to fresh directories on E:.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import sys

os.environ["MPLCONFIGDIR"] = "E:/CFD-Learning-Lab/temp/matplotlib-vortex"
os.environ["TEMP"] = "E:/CFD-Learning-Lab/temp"
os.environ["TMP"] = "E:/CFD-Learning-Lab/temp"
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Circle
import numpy as np
from PIL import Image

NAMES = {"fluent": "ANSYS Fluent", "openfoam": "OpenFOAM"}
TIME_TOL = 1e-7
XY_TOL = 1e-8


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def e_drive(path):
    require(path.resolve().drive.upper() == "E:", "All output must remain on E: " + str(path))


def indexed(group, count, width=1, partial=False):
    """Assemble native saved IDs; dataset order and Gmsh order are irrelevant."""
    values = np.full((count, width), np.nan)
    seen = np.zeros(count, dtype=bool)
    for dataset in group.values():
        require(isinstance(dataset, h5py.Dataset), "Unsupported nested HDF block")
        first = int(np.asarray(dataset.attrs["minId"]).ravel()[0]) - 1
        last = int(np.asarray(dataset.attrs["maxId"]).ravel()[0])
        block = np.asarray(dataset[:]).reshape(-1, width)
        require(0 <= first < last <= count and len(block) == last - first,
                "Invalid saved HDF IDs: " + dataset.name)
        require(not seen[first:last].any(), "Overlapping saved HDF blocks")
        values[first:last] = block
        seen[first:last] = True
    require(partial or seen.all(), "Incomplete saved HDF IDs")
    return values[:, 0] if width == 1 else values


def fluent_mesh(case):
    with h5py.File(case, "r") as handle:
        mesh = handle["meshes/1"]
        top = mesh["cells/zoneTopology"]
        require(np.asarray(top["minId"]).tolist() == [1]
                and len(np.asarray(top["maxId"])) == 1, "Expected one native cell zone")
        count = int(top["maxId"][0])
        nodes = mesh["nodes/coords"]
        total = max(int(np.asarray(d.attrs["maxId"]).ravel()[0]) for d in nodes.values())
        width = next(iter(nodes.values())).shape[1]
        require(width in (2, 3), "Expected 2D Fluent coordinates")
        xyz = indexed(nodes, total, width)
        require(np.isfinite(xyz).all() and (width == 2 or np.ptp(xyz[:, 2]) < XY_TOL),
                "Fluent mesh is not finite and planar")
        xy = xyz[:, :2]
        blocks = mesh["faces/nodes"]
        faces = max(int(np.asarray(b["nnodes"].attrs.get("maxId", b.attrs.get("maxId"))).ravel()[0])
                    for b in blocks.values())
        edges = np.full((faces, 2), -1, dtype=np.int64)
        for block in blocks.values():
            sizes = block["nnodes"]
            first = int(np.asarray(sizes.attrs.get("minId", block.attrs.get("minId", [1]))).ravel()[0]) - 1
            require(np.all(sizes[:] == 2), "Fluent faces must have two nodes")
            values = np.asarray(block["nodes"][:], dtype=np.int64).reshape(-1, 2) - 1
            require(np.all(edges[first:first + len(values)] == -1), "Repeated saved face IDs")
            edges[first:first + len(values)] = values
        require(np.all(edges >= 0) and np.all(edges < len(xy)), "Invalid Fluent edge nodes")
        owner = indexed(mesh["faces/c0"], faces).astype(np.int64) - 1
        other = indexed(mesh["faces/c1"], faces, partial=True)
        internal = np.isfinite(other)
        neighbour = other[internal].astype(np.int64) - 1
        require(np.all((owner >= 0) & (owner < count))
                and np.all((neighbour >= 0) & (neighbour < count)), "Invalid Fluent cell IDs")
        keys = np.unique(np.concatenate([
            np.repeat(owner, 2) * len(xy) + edges.ravel(),
            np.repeat(neighbour, 2) * len(xy) + edges[internal].ravel()]))
        ids, sizes = np.unique(keys // len(xy), return_counts=True)
        require(np.array_equal(ids, np.arange(count)) and np.all(sizes == 3),
                "Every native Fluent cell must have three unique vertices")
    return validate_mesh(xy, (keys % len(xy)).reshape(count, 3), [case])


def foam_text(path):
    return re.sub(r"/\*.*?\*/|//[^\n]*", "", path.read_text(encoding="utf-8"), flags=re.S)


def foam_list(path):
    raw = foam_text(path)
    header = re.search(r"\bFoamFile\s*\{([^}]*)\}", raw, re.S)
    require(header is not None and re.search(r"\bformat\s+ascii\s*;", header[1]),
            "Only native ASCII OpenFOAM lists are supported: " + str(path))
    match = re.fullmatch(r"\s*(\d+)\s*\((.*)\)\s*;?\s*", raw[header.end():], re.S)
    require(match is not None, "Unsupported OpenFOAM list: " + str(path))
    return int(match[1]), match[2]


def foam_labels(path):
    count, body = foam_list(path)
    values = np.fromstring(body, sep=" ", dtype=np.int64)
    require(len(values) == count, "Invalid OpenFOAM label count")
    return values


def openfoam_mesh(folder):
    candidates = [folder / "constant/polyMesh", folder / "0/polyMesh"]
    mesh = next((p for p in candidates if (p / "points").is_file()), None)
    require(mesh is not None, "No reconstructed native OpenFOAM polyMesh")
    files = [mesh / name for name in ("points", "faces", "owner", "neighbour")]
    point_count, body = foam_list(files[0])
    points = np.fromstring(body.replace("(", " ").replace(")", " "), sep=" ")
    require(points.size == point_count * 3 and np.isfinite(points).all(), "Invalid native points")
    points = points.reshape(-1, 3)
    planes = np.unique(points[:, 2])
    require(len(planes) == 2, "Expected a single extruded prism layer")
    face_count, body = foam_list(files[1])
    matches = list(re.finditer(r"(\d+)\s*\(([^()]*)\)", body))
    sizes = np.array([int(m[1]) for m in matches], dtype=np.int64)
    vertices = np.concatenate([np.fromstring(m[2], sep=" ", dtype=np.int64) for m in matches])
    require(len(sizes) == face_count and np.all(np.isin(sizes, [3, 4]))
            and len(vertices) == sizes.sum(), "Expected native prism triangle/quad faces")
    owner, neighbour = foam_labels(files[2]), foam_labels(files[3])
    require(len(owner) == face_count and len(neighbour) <= len(owner)
            and np.all(vertices >= 0) and np.all(vertices < len(points)), "Invalid native face IDs")
    count = int(owner.max()) + 1
    require(np.all(owner >= 0) and np.all(neighbour >= 0) and np.all(neighbour < count)
            and np.all(sizes[:len(neighbour)] == 4), "Invalid OpenFOAM cell adjacency")
    incidence = np.bincount(owner, minlength=count) + np.bincount(neighbour, minlength=count)
    caps = np.bincount(owner[sizes == 3], minlength=count)
    require(np.all(incidence == 5) and np.all(caps == 2), "Each native cell must be a triangular prism")
    stop = int(sizes[:len(neighbour)].sum())
    keys = np.unique(np.concatenate([
        np.repeat(owner, sizes) * len(points) + vertices,
        np.repeat(neighbour, sizes[:len(neighbour)]) * len(points) + vertices[:stop]]))
    ids, counts = np.unique(keys // len(points), return_counts=True)
    require(np.array_equal(ids, np.arange(count)) and np.all(counts == 6), "Invalid native prism vertices")
    xyz_nodes = (keys % len(points)).reshape(count, 6)
    xy, projected_ids = np.unique(points[:, :2], axis=0, return_inverse=True)
    projected = np.sort(projected_ids[xyz_nodes], axis=1)
    require(np.all(projected[:, ::2] == projected[:, 1::2])
            and np.all(np.diff(projected[:, ::2], axis=1) > 0),
            "Native prisms must project to three distinct XY nodes")
    require(np.all(np.sort(points[xyz_nodes, 2], axis=1) == np.tile(np.repeat(planes, 3), (count, 1))),
            "Each prism must have three vertices on each z-plane")
    result = validate_mesh(xy, projected[:, ::2], files)
    result["span_m"] = float(planes[1] - planes[0])
    return result


def validate_mesh(xy, triangles, inputs):
    require(len(triangles) >= 100 and np.isfinite(xy).all(), "Insufficient or invalid native mesh")
    require(np.min(xy[:, 0]) >= -XY_TOL and np.max(xy[:, 0]) <= 2.2 + XY_TOL
            and np.min(xy[:, 1]) >= -XY_TOL and np.max(xy[:, 1]) <= .41 + XY_TOL,
            "Native nodes outside declared channel")
    triangles = triangles.copy()
    xyz = xy[triangles]
    ab, ac = xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]
    cross = ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]
    require(np.all(np.abs(cross) > 1e-16), "Degenerate native triangle")
    triangles[cross < 0] = triangles[cross < 0][:, [0, 2, 1]]
    area = np.abs(cross) / 2
    centres = xy[triangles].mean(axis=1)
    require(np.all(np.linalg.norm(xy - [.2, .2], axis=1) >= .05 - XY_TOL),
            "Native mesh vertex lies inside cylinder")
    require(np.all(np.linalg.norm(centres - [.2, .2], axis=1) >= .05 - 1e-4),
            "Native fluid triangle bridges cylinder")
    edges = np.sort(np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    require(np.all(np.isin(counts, [1, 2])), "Nonmanifold native mesh")
    boundary = xy[unique[counts == 1]]
    outer = np.zeros(len(boundary), dtype=bool)
    for axis, value in [(0, 0), (0, 2.2), (1, 0), (1, .41)]:
        outer |= np.all(np.abs(boundary[:, :, axis] - value) < XY_TOL, axis=1)
    cylinder = np.all(np.abs(np.linalg.norm(boundary - [.2, .2], axis=2) - .05) < XY_TOL, axis=1)
    require(np.all(outer | cylinder) and cylinder.sum() >= 16,
            "Native fluid boundaries differ from channel and polygonal cylinder")
    cylinder_edges = boundary[cylinder]
    # Boundary edges are sorted by ID, not oriented; use a circularly sorted polygon.
    cylinder_nodes = np.unique(cylinder_edges.reshape(-1, 2), axis=0)
    angles = np.arctan2(cylinder_nodes[:, 1] - .2, cylinder_nodes[:, 0] - .2)
    polygon = cylinder_nodes[np.argsort(angles)] - [.2, .2]
    polygon_area = abs(np.sum(polygon[:, 0] * np.roll(polygon[:, 1], -1)
                              - polygon[:, 1] * np.roll(polygon[:, 0], -1))) / 2
    expected = 2.2 * .41 - polygon_area
    require(abs(area.sum() - expected) < 1e-9, "Native triangles do not tile declared fluid area")
    tri = mtri.Triangulation(xy[:, 0], xy[:, 1], triangles)
    require(tri.get_trifinder()(.2, .2) == -1, "Native triangle fills cylinder centre")
    return dict(xy=xy, triangles=triangles, centres=centres, area=area,
                triangulation=tri, inputs=inputs, cells=len(triangles),
                checks=dict(cells=len(triangles), nodes=len(xy), native_topology=True,
                            cylinder_boundary_edges=int(cylinder.sum()), fluid_area_m2=float(area.sum()),
                            polygonal_cylinder_area_m2=float(polygon_area), no_solid_bridge=True))


def foam_field(path, count, width, dimensions):
    raw = foam_text(path)
    require(re.search(r"\bformat\s+ascii\s*;", raw), "OpenFOAM fields must be ASCII")
    units = re.search(r"\bdimensions\s*\[([^]]+)\]\s*;", raw)
    require(units is not None and np.array_equal(np.fromstring(units[1], sep=" "), dimensions),
            "Unexpected native dimensions: " + str(path))
    kind = "vector" if width == 3 else "scalar"
    match = re.search(r"\binternalField\s+nonuniform\s+List<" + kind + r">\s+(\d+)\s*\((.*?)\)\s*;", raw, re.S)
    if match:
        require(int(match[1]) == count, "OpenFOAM field cell count differs")
        values = np.fromstring(match[2].replace("(", " ").replace(")", " "), sep=" ")
    else:
        match = re.search(r"\binternalField\s+uniform\s+([^;]+);", raw)
        require(match is not None, "Unsupported native internalField")
        value = np.fromstring(match[1].replace("(", " ").replace(")", " "), sep=" ")
        require(value.size == width, "Invalid native uniform field")
        values = np.tile(value, count)
    require(values.size == count * width and np.isfinite(values).all(), "Invalid native field values")
    return values.reshape(count, width) if width > 1 else values


def native_field(solver, snapshot, mesh):
    count = mesh["cells"]
    if solver == "fluent":
        with h5py.File(snapshot["path"], "r") as handle:
            fields = handle["results/1/phase-1/cells"]
            u, v, p = [indexed(fields[name], count) for name in ("SV_U", "SV_V", "SV_P")]
            for name, expected in [("SV_DENSITY", 1.), ("SV_MU_LAM", .001)]:
                if name in fields:
                    require(np.allclose(indexed(fields[name], count), expected, rtol=1e-8, atol=1e-12),
                            "Native Fluent material differs: " + name)
        inputs = [snapshot["path"]]
    else:
        velocity = foam_field(snapshot["path"] / "U", count, 3, [0, 1, -1, 0, 0, 0, 0])
        require(np.max(np.abs(velocity[:, 2])) < 1e-10, "Nonzero native spanwise velocity")
        u, v = velocity[:, 0], velocity[:, 1]
        p = foam_field(snapshot["path"] / "p", count, 1, [0, 2, -2, 0, 0, 0, 0])
        inputs = [snapshot["path"] / name for name in ("U", "p")]
    require(all(np.isfinite(a).all() and a.shape == (count,) for a in (u, v, p)), "Invalid saved native fields")
    return dict(u=u, v=v, p=p, speed=np.hypot(u, v), inputs=inputs)


def reconstructed_vorticity(field, mesh):
    """Area-weighted native cell-to-node reconstruction, differentiated per cell.

    This is a postprocessed gradient, not a solver-exported vorticity or an exact
    near-wall derivative. Only incident saved fluid triangles enter each node.
    """
    xy, tri, weights = mesh["xy"], mesh["triangles"], mesh["area"]
    ids = tri.ravel()
    sums = np.bincount(ids, weights=np.repeat(weights, 3), minlength=len(xy))
    require(np.all(sums > 0), "Unconnected saved mesh node")
    nodal = [np.bincount(ids, weights=np.repeat(weights * field[k], 3), minlength=len(xy)) / sums
             for k in ("u", "v")]
    xyz = xy[tri]
    x, y = xyz[:, :, 0], xyz[:, :, 1]
    twice = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (y[:, 1] - y[:, 0]) * (x[:, 2] - x[:, 0])
    u, v = [a[tri] for a in nodal]
    dv_dx = (v[:, 0] * (y[:, 1] - y[:, 2]) + v[:, 1] * (y[:, 2] - y[:, 0])
             + v[:, 2] * (y[:, 0] - y[:, 1])) / twice
    du_dy = (u[:, 0] * (x[:, 2] - x[:, 1]) + u[:, 1] * (x[:, 0] - x[:, 2])
             + u[:, 2] * (x[:, 1] - x[:, 0])) / twice
    result = dv_dx - du_dy
    require(np.isfinite(result).all(), "Nonfinite reconstructed vorticity")
    return result


def discover(fluent, foam, manifest):
    items = json.loads(manifest.read_text(encoding="utf-8-sig"))
    if isinstance(items, dict):
        items = items.get("snapshots", items.get("frames"))
    require(isinstance(items, list) and len(items) >= 2, "Fluent snapshot manifest needs at least two entries")
    first = []
    for item in items:
        path = Path(item["data"])
        if not path.is_absolute():
            path = fluent / path
        if not path.is_file():
            alternate = manifest.parent / item["data"]
            require(alternate.is_file(), "Missing native Fluent snapshot: " + str(path))
            path = alternate
        time = float(item["time"])
        require(np.isfinite(time) and time >= 0, "Invalid physical snapshot time")
        first.append(dict(time=time, step=item.get("step"), path=path.resolve()))
    second = []
    for path in foam.iterdir():
        if path.is_dir():
            try:
                time = float(path.name)
            except ValueError:
                continue
            if np.isfinite(time) and time >= 0 and (path / "U").is_file() and (path / "p").is_file():
                second.append(dict(time=time, path=path.resolve()))
    for name, items in [("Fluent", first), ("OpenFOAM", second)]:
        items.sort(key=lambda row: row["time"])
        require(len(items) >= 2 and np.all(np.diff([r["time"] for r in items]) > TIME_TOL),
                "Insufficient or duplicate " + name + " physical times")
    paired = []
    used = set()
    for row in first:
        matches = [j for j, other in enumerate(second) if abs(row["time"] - other["time"]) <= TIME_TOL]
        require(len(matches) <= 1, "Ambiguous physical-time match")
        if matches:
            j = matches[0]
            used.add(j)
            paired.append(dict(time=row["time"], fluent=row, openfoam=second[j]))
    require(len(paired) >= 2, "Need at least two exactly paired physical times")
    checks = dict(time_tolerance_s=TIME_TOL, interpolation=False, paired_times_s=[r["time"] for r in paired],
                  fluent_only_times_s=[r["time"] for r in first if not any(abs(r["time"] - p["time"]) <= TIME_TOL for p in paired)],
                  openfoam_only_times_s=[r["time"] for j, r in enumerate(second) if j not in used])
    return paired, checks


def shared_cell_map(first, second):
    require(first["cells"] == second["cells"], "The paired native meshes have different cell counts")
    def keys(mesh):
        rows = np.rint(mesh["xy"][mesh["triangles"]] / XY_TOL).astype(np.int64)
        return [tuple(sorted(map(tuple, row))) for row in rows]
    source, target = keys(first), keys(second)
    lookup = {key: index for index, key in enumerate(target)}
    require(len(lookup) == len(source) and len(set(source)) == len(source)
            and all(key in lookup for key in source), "The codes do not share the same native triangle geometry")
    return np.array([lookup[key] for key in source], dtype=np.int64)


def decorate(ax, name):
    ax.add_patch(Circle((.2, .2), .05, facecolor="#bbc0c5", edgecolor="#495158", lw=.7, zorder=5))
    ax.set(xlim=(0, 2.2), ylim=(0, .41), xlabel="$x$ (m)", ylabel="$y$ (m)", aspect="equal")
    ax.set_xticks([0, .2, .6, 1., 1.4, 1.8, 2.2])
    ax.set_yticks([0, .2, .41])
    ax.set_title(NAMES[name], loc="left", fontsize=10, pad=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def make_panel(meshes, initial, kind, limits):
    fig, axes = plt.subplots(2, 1, figsize=(12.8, 6.1), layout="constrained")
    norm = Normalize(*limits) if kind == "speed" else TwoSlopeNorm(vmin=limits[0], vcenter=0, vmax=limits[1])
    artists = {}
    for ax, name in zip(axes, NAMES):
        artists[name] = ax.tripcolor(meshes[name]["triangulation"], facecolors=initial[name][kind],
                                    shading="flat", cmap="viridis" if kind == "speed" else "RdBu_r", norm=norm,
                                    rasterized=True)
        decorate(ax, name)
    cb = fig.colorbar(next(iter(artists.values())), ax=axes, location="bottom", fraction=.055, pad=.035, aspect=70)
    cb.set_label("Speed $|\\mathbf{u}|$ (m/s)" if kind == "speed" else "Reconstructed vorticity $\\omega_z$ (s$^{-1}$)")
    title = fig.suptitle("", fontsize=14, fontweight="medium")
    return fig, artists, title


def render_series(meshes, frames, out, kind, limits, fps, stride, no_mp4):
    fig, artists, title = make_panel(meshes, frames[0], kind, limits)
    selected = list(range(0, len(frames), stride))
    if selected[-1] != len(frames) - 1:
        selected.append(len(frames) - 1)
    require(len(selected) >= 2, "Animation needs at least two distinct physical frames")
    gifs = []
    gif_delay_ms = max(10, round(100 / fps) * 10)
    video = None
    mp4_status = "disabled by command" if no_mp4 else "imageio_ffmpeg unavailable"
    if not no_mp4:
        try:
            import imageio_ffmpeg
            require(Path(imageio_ffmpeg.get_ffmpeg_exe()).is_file(), "FFmpeg executable unavailable")
            mp4_status = "available"
        except ImportError:
            pass
    try:
        for index in selected:
            frame = frames[index]
            for name in NAMES:
                artists[name].set_array(frame[name][kind])
            title.set_text(("Cylinder wake · speed" if kind == "speed" else "Cylinder wake · signed vorticity")
                           + f"  |  $Re_D=100$  |  $t={frame['time']:.3f}$ s")
            fig.canvas.draw()
            rgba = np.asarray(fig.canvas.buffer_rgba())
            rgb = np.ascontiguousarray(rgba[:, :, :3])
            if mp4_status == "available" and video is None:
                video = imageio_ffmpeg.write_frames(str(out / ("vortex-" + kind + ".mp4")),
                            (rgb.shape[1], rgb.shape[0]), fps=fps, codec="libx264", pix_fmt_out="yuv420p",
                            macro_block_size=2, output_params=["-crf", "20", "-movflags", "+faststart"])
                video.send(None)
            if video is not None:
                video.send(rgb)
            # One palette per GIF frame preserves the actual fixed colormap values.
            gifs.append(Image.fromarray(rgb).quantize(colors=256))
            if index == len(frames) - 1:
                fig.savefig(out / ("vortex-" + kind + "-final.png"), dpi=220, facecolor="white")
                fig.savefig(out / ("vortex-" + kind + "-final.svg"), dpi=220, facecolor="white")
        gifs[0].save(out / ("vortex-" + kind + ".gif"), save_all=True, append_images=gifs[1:],
                     duration=gif_delay_ms, loop=0, optimize=False, disposal=2)
    finally:
        if video is not None:
            video.close()
        plt.close(fig)
    return dict(frames=len(selected), physical_times_s=[frames[i]["time"] for i in selected],
                requested_playback_fps=fps, gif_frame_duration_ms=gif_delay_ms,
                gif_effective_fps=1000 / gif_delay_ms, gif_duration_s=len(selected) * gif_delay_ms / 1000,
                mp4="written" if video is not None else mp4_status,
                limits=limits, no_temporal_interpolation=True)


def force_histories(args, fluent, foam, out, end):
    paths = {name: folder / "force-history.csv" for name, folder in [("fluent", fluent), ("openfoam", foam)]}
    for spec in args.force_csv:
        name, path = spec.split("=", 1)
        require(name in NAMES, "--force-csv requires fluent=path or openfoam=path")
        paths[name] = Path(path).resolve()
    found, inputs = {}, []
    for name, path in paths.items():
        if not path.is_file():
            continue
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            require(set(("time_s", "Cd", "Cl")).issubset(reader.fieldnames or []),
                    "Native force CSV needs time_s,Cd,Cl: " + str(path))
            rows = np.array([[float(r[k]) for k in ("time_s", "Cd", "Cl")] for r in reader])
        require(rows.ndim == 2 and rows.shape[1] == 3 and np.isfinite(rows).all(), "Invalid native force history")
        rows = rows[(rows[:, 0] >= 0) & (rows[:, 0] <= end + TIME_TOL)]
        require(len(rows) >= 2 and np.all(np.diff(rows[:, 0]) > 0), "Force history needs increasing physical times")
        found[name], inputs = rows, inputs + [path]
        np.savetxt(out / (name + "-actual-force-history.csv"), rows, delimiter=",", header="time_s,Cd,Cl", comments="", fmt="%.12g")
    if not found:
        return dict(status="omitted: no actual native force-history.csv provided"), inputs
    fig, axes = plt.subplots(2, 1, figsize=(9, 5), sharex=True, layout="constrained")
    for name, rows in found.items():
        for ax, column in zip(axes, (1, 2)):
            ax.plot(rows[:, 0], rows[:, column], label=NAMES[name], lw=1.2)
    for ax, ylabel in zip(axes, ("Drag coefficient $C_D$", "Lift coefficient $C_L$")):
        ax.set_ylabel(ylabel)
        ax.grid(alpha=.18)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].legend(frameon=False)
    axes[1].set_xlabel("Physical time $t$ (s)")
    fig.suptitle("Actual native force coefficient histories · $Re_D=100$")
    fig.savefig(out / "vortex-native-force-histories.png", dpi=220)
    plt.close(fig)
    return dict(status="plotted native source CSV only", solvers=list(found),
                physical_time_ranges_s={k: [float(v[0, 0]), float(v[-1, 0])] for k, v in found.items()},
                coefficient_reference="rho=1 kg/m3, mean inlet U=1 m/s, cylinder D=0.1 m; supplied native reports"), inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--fluent-dir", type=Path)
    parser.add_argument("--openfoam-dir", type=Path)
    parser.add_argument("--fluent-case", type=Path)
    parser.add_argument("--fluent-manifest", type=Path)
    parser.add_argument("--min-end-time", type=float, default=0., help="Refuse an incomplete paired pilot")
    parser.add_argument("--fps", type=float, default=20.)
    parser.add_argument("--frame-stride", type=int, default=1)
    parser.add_argument("--no-vorticity", action="store_true")
    parser.add_argument("--no-mp4", action="store_true")
    parser.add_argument("--force-csv", action="append", default=[], metavar="SOLVER=PATH")
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    e_drive(out)
    require(not out.exists(), "Use a fresh figure/export directory")
    require(np.isfinite(args.fps) and args.fps > 0 and args.frame_stride >= 1, "Invalid frame rate or stride")
    fluent = (args.fluent_dir or run / "fluent").resolve()
    foam = (args.openfoam_dir or run / "openfoam").resolve()
    case = (args.fluent_case or fluent / "vortex-setup.cas.h5").resolve()
    manifest = (args.fluent_manifest or fluent / "snapshots.json").resolve()
    paired, time_checks = discover(fluent, foam, manifest)
    require(paired[-1]["time"] >= args.min_end_time - TIME_TOL, "Paired fields have not reached --min-end-time")
    meshes = dict(fluent=fluent_mesh(case), openfoam=openfoam_mesh(foam))
    mapping = shared_cell_map(meshes["fluent"], meshes["openfoam"])
    fields, summaries, inputs = [], [], [manifest]
    for mesh in meshes.values():
        inputs += mesh["inputs"]
    for row in paired:
        frame = dict(time=row["time"])
        for name in NAMES:
            field = native_field(name, row[name], meshes[name])
            if not args.no_vorticity:
                field["vorticity"] = reconstructed_vorticity(field, meshes[name])
            inputs += field["inputs"]
            frame[name] = field
            summaries.append([row["time"], name, float(field["u"].min()), float(field["u"].max()),
                              float(field["v"].min()), float(field["v"].max()), float(field["speed"].min()),
                              float(field["speed"].max()), float(field["p"].min()), float(field["p"].max()),
                              float(field["vorticity"].min()) if "vorticity" in field else "",
                              float(field["vorticity"].max()) if "vorticity" in field else ""])
        fields.append(frame)
    speed_max = max(float(frame[name]["speed"].max()) for frame in fields for name in NAMES)
    limits = dict(speed=[0., max(speed_max, 1e-12)])
    if not args.no_vorticity:
        maximum = max(float(np.abs(frame[name]["vorticity"]).max()) for frame in fields for name in NAMES)
        maximum = max(maximum, 1e-12)
        limits["vorticity"] = [-maximum, maximum]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.linewidth": .65,
                         "svg.fonttype": "path", "figure.dpi": 100})
    out.mkdir(parents=True)
    numeric = out / "actual-cell-fields"
    numeric.mkdir()
    columns = ["time_s", "native_cell_id_1based", "x_m", "y_m", "u_m_s", "v_m_s", "speed_m_s", "gauge_pressure_Pa"]
    if not args.no_vorticity:
        columns.append("reconstructed_vorticity_z_per_s")
    for index, frame in enumerate(fields):
        for name in NAMES:
            mesh, field = meshes[name], frame[name]
            array = np.column_stack([np.full(mesh["cells"], frame["time"]),
                                     np.arange(1, mesh["cells"] + 1), mesh["centres"],
                                     *(field[k] for k in ("u", "v", "speed", "p")),
                                     *([field["vorticity"]] if not args.no_vorticity else [])])
            path = numeric / f"{name}-frame-{index:05d}-t-{frame['time']:.9f}.csv.gz"
            with gzip.open(path, "wt", compresslevel=1, encoding="utf-8", newline="") as stream:
                np.savetxt(stream, array, delimiter=",", header=",".join(columns), comments="", fmt="%.12g")
    with (out / "native-field-time-ranges.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", "solver", "u_min_m_s", "u_max_m_s", "v_min_m_s", "v_max_m_s", "speed_min_m_s",
                         "speed_max_m_s", "pressure_min_Pa", "pressure_max_Pa", "vorticity_min_per_s", "vorticity_max_per_s"])
        writer.writerows(summaries)
    renders = {kind: render_series(meshes, fields, out, kind, interval, args.fps, args.frame_stride, args.no_mp4)
               for kind, interval in limits.items()}
    forces, force_inputs = force_histories(args, fluent, foam, out, paired[-1]["time"])
    inputs += force_inputs
    differences = []
    for frame in fields:
        du = frame["fluent"]["u"] - frame["openfoam"]["u"][mapping]
        dv = frame["fluent"]["v"] - frame["openfoam"]["v"][mapping]
        differences.append([frame["time"], float(np.max(np.hypot(du, dv))), float(np.sqrt(np.mean(du ** 2 + dv ** 2)))])
    np.savetxt(out / "matched-native-cell-velocity-differences.csv", differences, delimiter=",",
               header="time_s,max_vector_difference_m_s,rms_vector_difference_m_s", comments="", fmt="%.12g")
    record = dict(geometry=dict(channel_m=[0, 2.2, 0, .41], cylinder_centre_m=[.2, .2], cylinder_radius_m=.05),
                  declared_physics=dict(Re_D=100, density_kg_m3=1, dynamic_viscosity_Pa_s=.001,
                                        mean_inlet_velocity_m_s=1, parabolic_inlet_peak_m_s=1.5),
                  topology={name: mesh["checks"] for name, mesh in meshes.items()},
                  shared_native_triangle_geometry=True, cell_matching_tolerance_m=XY_TOL, timing=time_checks,
                  animations=renders, all_colour_limits_include_all_paired_native_fields=True,
                  native_fields=True, spatial_speed_smoothing=False, synthetic_vortices=False,
                  vorticity_method="Area-weighted native cell-to-node velocity average, differentiated linearly within each native fluid triangle; no boundary values appended",
                  force_histories=forces, physical_validation=False, mesh_independence_established=False,
                  periodic_shedding_or_statistical_stationarity_established=False,
                  scope="Actual saved transient fields over the recorded interval; pilot duration alone does not establish a periodic wake")
    (out / "checks.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    script = Path(__file__).resolve()
    (out / "executed-plotter.py").write_bytes(script.read_bytes())
    unique_inputs = sorted(set(Path(p).resolve() for p in inputs), key=str)
    source = dict(command_arguments=sys.argv[1:], executed_plotter_sha256=sha256(script),
                  inputs=[dict(path=str(path), bytes=path.stat().st_size, sha256=sha256(path)) for path in unique_inputs])
    (out / "input-hash-manifest.json").write_text(json.dumps(source, indent=2) + "\n", encoding="utf-8")
    caption = ("# Actual transient cylinder fields\n\n"
               f"Paired native Fluent/OpenFOAM snapshots span {paired[0]['time']:.9g} to {paired[-1]['time']:.9g} s "
               f"at {len(paired)} matching physical times. Time interpolation is never used. The channel is "
               "2.2 by 0.41 m, with a cylinder centred at (0.2,0.2) m of diameter 0.1 m. Declared conditions are "
               "rho=1 kg/m3, mu=0.001 Pa.s, parabolic inlet of mean 1 m/s and peak 1.5 m/s, and Re_D=100. "
               "These declarations must be checked against separate native run receipts.\n\n"
               "Speed is the actual saved cell-centred hypot(u,v), shown flat on each native fluid triangle. "
               "Triangles are reconstructed from saved face/owner/neighbour IDs independently in each code; "
               "their projected geometric correspondence is verified. Native cell order is never assumed. "
               "Gray masks the polygonal cylinder with an analytic circle; no artificial flow is drawn there. "
               "The x and y axes use metres and equal geometric scale. Every speed frame and both codes use "
               f"one linear color range [{limits['speed'][0]:.9g}, {limits['speed'][1]:.9g}] m/s, containing all "
               "paired native values. No streamlines or synthetic vortex motion are added.\n\n")
    if not args.no_vorticity:
        caption += ("Signed vorticity is an explicitly postprocessed numerical reconstruction: area-weighted "
                    "cell-centred u and v are averaged to their own incident native mesh vertices, then curl "
                    "dv/dx-du/dy is calculated from the linear field within each saved fluid triangle. No "
                    "Delaunay mesh, cross-cylinder triangle, analytic vortex, imposed wall sample or synthetic "
                    "gradient enters this calculation. Boundary derivatives have reconstruction error and "
                    "are not native solver gradient exports. A single symmetric linear range "
                    f"[{limits['vorticity'][0]:.9g}, {limits['vorticity'][1]:.9g}] per second includes every "
                    "reconstructed value in both codes and every paired time. No values are clipped.\n\n")
    caption += (f"Requested playback is {args.fps:g} frames/s; GIF uses the nearest 10 ms frame delay, "
                "recorded in checks.json. Physical time is printed on every frame; playback speed is "
                "separate from physical simulation speed. The exact rendered physical "
                "times, any unpaired snapshots, native field ranges, matching-cell differences, and native "
                "source hashes are recorded in the accompanying checks and CSVs. Compressed CSVs contain "
                "the actual per-cell fields for every paired time; OpenFOAM kinematic p is expressed as Pa "
                "using rho=1 kg/m3. Forces are plotted only if actual source force CSVs were supplied.\n\n"
                "A transient picture is not evidence of a statistically stationary or periodic shedding "
                "state, a validated Strouhal number, mesh independence or physical validation. Fluent and "
                "OpenFOAM calculated the fields; Codex reconstructed native topology and prepared exports.\n")
    (out / "CAPTIONS.md").write_text(caption, encoding="utf-8")
    output_files = sorted([p for p in out.rglob("*") if p.is_file()], key=str)
    (out / "output-hash-manifest.json").write_text(json.dumps([
        dict(path=str(p.relative_to(out)), bytes=p.stat().st_size, sha256=sha256(p)) for p in output_files], indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(output=str(out), frames=len(fields), time_interval_s=[paired[0]["time"], paired[-1]["time"]],
                          animations=renders, force_histories=forces), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
