"""Audit weir overflow and draw saved water fractions on the native cell faces."""
import argparse
import csv
import hashlib
import json
import re
import runpy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np

CELLS = 5080
TIMES = range(2, 61, 2)
NUMBER = r"[-+0-9.eE]+"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compact_list(data, dtype, components=1, offset=0):
    # Native binary arrays from the 64-bit, little-endian Foundation build.
    start = re.search(rb"\n([0-9]+)\r?\n\(", data[offset:])
    assert start
    count = int(start[1])
    position = offset + start.end()
    array = np.frombuffer(data, dtype=dtype, offset=position, count=count*components)
    end = position + array.nbytes
    assert data[end:end+1] == b")"
    return array, end+1


def values(data, name, count, kind, reader):
    if re.search(rb"\bformat\s+ascii\s*;", data):
        return reader(data.decode("ascii"), name, count, kind)
    label = rb"\b"+name.encode()+rb"\s+"
    uniform = re.search(label+rb"uniform\s+([^;]+);", data)
    match = re.search(label+rb"nonuniform\s+List<"+kind.encode()+rb">\s+(\d+)\s*\(", data)
    width = 3 if kind == "vector" else 1
    if uniform and (not match or uniform.start() < match.start()):
        raw = np.fromstring(uniform[1].replace(b"(", b" ").replace(b")", b" "), sep=" ")
        assert len(raw) == width
        array = np.tile(raw, (count, 1))
    else:
        assert match and int(match[1]) == count, (name, count)
        array = np.frombuffer(data, dtype="<f8", offset=match.end(), count=count*width).reshape(count, width)
        assert data[match.end()+array.nbytes:match.end()+array.nbytes+1] == b")"
    assert np.isfinite(array).all()
    return array if width == 3 else array[:,0]


def plane_polygons(case, centres, volumes):
    mesh = case/"constant/polyMesh"
    points, _ = compact_list((mesh/"points").read_bytes(), "<f8", 3)
    points = points.reshape(-1, 3)
    owners, _ = compact_list((mesh/"owner").read_bytes(), "<i4")
    data = (mesh/"faces").read_bytes()
    offsets, tail = compact_list(data, "<i4")
    nodes, _ = compact_list(data, "<i4", offset=tail)
    assert len(offsets) == len(owners)+1 and offsets[0] == 0 and offsets[-1] == len(nodes)
    assert np.all(np.diff(offsets) == 4) and np.all((nodes >= 0) & (nodes < len(points)))
    assert np.all((owners >= 0) & (owners < CELLS)) and np.isfinite(points).all()
    faces = points[nodes.reshape(-1, 4)]
    mask = np.all(np.isclose(faces[:, :, 2], -.5, atol=1e-10, rtol=0), axis=1)
    order = np.argsort(owners[mask])
    assert np.array_equal(owners[mask][order], np.arange(CELLS))
    xy = faces[mask][order, :, :2]
    x, y = xy.transpose(2, 0, 1)
    xx, yy = np.roll(x, -1, axis=1), np.roll(y, -1, axis=1)
    cross = x*yy - y*xx
    total = cross.sum(axis=1)
    centroids = np.column_stack((((x+xx)*cross).sum(axis=1), ((y+yy)*cross).sum(axis=1)))/(3*total[:,None])
    assert np.allclose(centroids, centres[:, :2], atol=1e-4, rtol=0)
    areas = .5*np.abs(total)
    assert np.allclose(areas, volumes, rtol=1e-5, atol=1e-6)
    return xy


def audit(case, foam_reader, output):
    reader = runpy.run_path(str(foam_reader))["foam_values"]

    def field(path, kind="scalar"):
        return values((case/path).read_bytes(), "internalField", CELLS, kind, reader)

    for line in (case/"input-sha256.txt").read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ", 1)
        assert digest(case/name) == expected, name
    mesh = (case/"check-mesh.log").read_text(encoding="ascii")
    log = (case/"solver.log").read_text(encoding="ascii")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*5080\s*$", mesh, re.M)
    assert re.search(r"^End\s*$", log, re.M)
    steps = np.array(re.findall(r"^Time = ("+NUMBER+r")s$", log, re.M), dtype=float)
    bulk = np.array(re.findall(r"^Courant Number mean: "+NUMBER+r" max: ("+NUMBER+r")", log, re.M), dtype=float)
    interface = np.array(re.findall(r"^Interface Courant Number mean: "+NUMBER+r" max: ("+NUMBER+r")", log, re.M), dtype=float)
    assert len(steps) > 0 and abs(steps[-1]-60) < 1e-8 and np.all(np.diff(steps)>0)
    assert len(bulk) >= len(steps) and len(interface) >= len(steps)
    assert np.isfinite(bulk).all() and np.isfinite(interface).all()
    centres, volumes = field("0/C", "vector"), field("0/Vc")
    assert np.all(volumes > 0)
    polygons = plane_polygons(case, centres, volumes)
    x, y, _ = centres.T
    crest = (x >= 0) & (x < 15) & (y >= 30)
    downstream = (x >= 30) & (y < 30)
    boundary = (case/"constant/polyMesh/boundary").read_text(encoding="ascii")
    outlet = re.search(r"\boutlet\s*\{([^{}]*)\}", boundary, re.S)
    count = int(re.search(r"\bnFaces\s+(\d+)\s*;", outlet[1])[1])
    assert count == 64
    rows, states = [], {}
    sources = ["input-sha256.txt", "check-mesh.log", "solver.log", "0/C", "0/Vc"] + [
        f"constant/polyMesh/{name}" for name in ("points", "faces", "owner", "boundary")
    ]
    for time in TIMES:
        clock = (case/str(time)/"uniform/time").read_text(encoding="ascii")
        step = int(re.search(r"^index\s+(\d+);", clock, re.M)[1])
        value = float(re.search(r"^value\s+([^;]+);", clock, re.M)[1])
        assert abs(value-time) < 1e-8 and abs(steps[step-1]-time) < 1e-8
        alpha = field(f"{time}/alpha.water")
        field(f"{time}/U", "vector")
        field(f"{time}/p_rgh")
        assert alpha.min() >= -1e-6 and alpha.max() <= 1+1e-6
        phi = (case/str(time)/"alphaPhi.water").read_bytes()
        assert re.search(rb"\bdimensions\s+\[0\s+3\s+-1\s+0\s+0\s+0\s+0\]", phi)
        patch = re.search(rb"(?m)^\s*outlet\s*\{\s*type\s+[^;]+;\s*(?=value\s+)", phi)
        assert patch, "Missing native outlet phase flux"
        net = values(phi[patch.end():], "value", count, "scalar", reader).sum()
        rows.append({"time_s": time, "native_step": step, "alpha_min": float(alpha.min()),
                     "alpha_max": float(alpha.max()), "crest_water_m3": float(alpha[crest]@volumes[crest]),
                     "downstream_water_m3": float(alpha[downstream]@volumes[downstream]),
                     "outlet_water_m3_s": float(net)})
        states[time] = alpha
        sources.extend(f"{time}/{name}" for name in ("uniform/time", "alpha.water", "U", "p_rgh", "alphaPhi.water"))
    assert rows[-1]["native_step"] == len(steps)
    first_crest = next(row for row in rows if row["crest_water_m3"] > 1e-6)
    first_below = next(row for row in rows if row["time_s"] > first_crest["time_s"] and row["downstream_water_m3"] > 1e-6)
    first_outlet = next(row for row in rows if row["outlet_water_m3_s"] > 1)
    checks = {"cells": CELLS, "native_steps": len(steps), "saved_states": len(rows),
              "outlet_faces": count, "max_logged_bulk_courant": float(bulk.max()),
              "max_logged_interface_courant": float(interface.max()),
              "first_saved_crest_water_s": first_crest["time_s"],
              "first_later_saved_downstream_water_s": first_below["time_s"],
              "first_saved_outlet_above_1_m3_s": first_outlet["time_s"],
              "samples": rows, "native_sha256": {name: digest(case/name) for name in sources}}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(checks, indent=2)+"\n", encoding="utf-8")
    with (output/"weir-history.csv").open("x", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 6.5), layout="constrained")
    for axis, time in zip(axes[0], (first_crest["time_s"], 60)):
        collection = PolyCollection(polygons, array=states[time], cmap="viridis", edgecolors="none")
        collection.set_clim(0, 1)
        axis.add_collection(collection)
        axis.set(xlim=(-18, 90), ylim=(0, 54), aspect="equal", xlabel="x (m)", ylabel="Elevation (m)", title=f"{time} s")
    fig.colorbar(collection, ax=list(axes[0]), label="Water volume fraction", shrink=.9)
    times = [row["time_s"] for row in rows]
    for key, label in (("crest_water_m3", "Above crest"), ("downstream_water_m3", "Below crest downstream")):
        axes[1,0].plot(times, [row[key] for row in rows], "o-", markersize=3, label=label)
    axes[1,0].set(xlabel="Time (s)", ylabel="Water volume (m³)", title="Water in measured regions")
    axes[1,0].legend()
    axes[1,1].plot(times, [row["outlet_water_m3_s"] for row in rows], "o-", markersize=3)
    axes[1,1].set(xlabel="Time (s)", ylabel="Outlet water flux (m³/s)", title="Sampled net outlet flux")
    for axis in axes[1]:
        axis.grid(alpha=.25)
    fig.suptitle("Overtopping a broad-crested weir | OpenFOAM", fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(output/f"weir-overflow.{ext}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "foam_reader", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    audit(args.case, args.foam_reader, args.output)
