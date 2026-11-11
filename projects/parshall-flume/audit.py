"""Audit native Parshall flume fields, phase flux and throat monitors."""
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

CELLS = 33950
FACES = 105983
INTERNAL = 97717
SAVES = np.arange(5, 41, 5)
ZONES = ("nearInlet", "upstream", "downstream", "nearOutlet")
WIDTHS = np.array([5., 3.04, 2.4, 5.])


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, reader_path, out):
    read = runpy.run_path(str(reader_path))["foam_values"]
    sources = []

    def native(name):
        sources.append(name)
        return (case/name).read_text("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest, name = line.split("  ", 1)
        assert sha(case/name) == digest, name
    mesh, log = native("check-mesh.log"), native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*33950\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd")
    steps = np.array(re.findall(r"^Time = ([-+\d.eE]+)s$", log, re.M), float)
    assert len(steps) and steps[-1] == 40 and np.all(np.diff(steps) > 0)
    courant = {}
    for name in ("Courant Number", "Interface Courant Number"):
        reports = np.array(re.findall(r"^"+name+r" mean: [-+\d.eE]+ max: ([-+\d.eE]+)", log, re.M), float)
        assert len(reports) >= len(steps) and np.all(np.isfinite(reports)) and np.all(reports >= 0)
        courant[name] = {"reports": len(reports), "max": float(reports.max())}

    def field(name, size=CELLS, kind="scalar"):
        data = read(native(name), "internalField", size, kind)
        assert np.all(np.isfinite(data)), name
        return data

    def indices(name, size):
        contents = native(f"constant/polyMesh/{name}")
        result = np.fromstring(re.search(r"\n"+str(size)+r"\s*\(\s*([^)]*)\)", contents)[1],
                               dtype=int, sep=" ")
        assert result.shape == (size,) and np.all((result >= 0) & (result < CELLS))
        return result

    owner, neighbour = indices("owner", FACES), indices("neighbour", INTERNAL)
    points = np.array([np.fromstring(line, sep=" ") for line in re.findall(
        r"^\(([^()]+)\)$", native("constant/polyMesh/points"), re.M)])
    nodes = np.array([np.fromstring(line, dtype=int, sep=" ") for line in re.findall(
        r"^4\(([^()]+)\)$", native("constant/polyMesh/faces"), re.M)])
    assert points.shape == (38220, 3) and nodes.shape == (FACES, 4)
    assert np.all((nodes >= 0) & (nodes < len(points)))
    faces = points[nodes]
    centres, volumes = field("0/C", kind="vector"), field("0/Vc")
    assert np.all(volumes > 0)
    plane = np.flatnonzero(np.all(abs(faces[:INTERNAL, :, 1]) < 1e-11, axis=1))
    assert len(plane) == 2425
    adjacent = np.column_stack([owner[plane], neighbour[plane]])
    assert len(np.unique(adjacent)) == 4850
    assert np.all(np.prod(centres[adjacent, 1], axis=1) < 0)
    assert np.allclose(centres[adjacent[:,0]][:,[0,2]], centres[adjacent[:,1]][:,[0,2]],
                       rtol=0, atol=1e-8)
    polygons = faces[plane][:, :, [0, 2]]

    zone_dictionary = native("constant/polyMesh/faceZones")
    zone_faces = []
    zone_locations = []
    for name in ZONES:
        entry = re.search(r"\b"+name+r"\s*\{(.*?)\n\}", zone_dictionary, re.S)[1]
        ids = np.fromstring(re.search(r"faceLabels\s+List<label>\s+350\s*\(([^)]*)\)", entry)[1],
                            dtype=int, sep=" ")
        assert ids.shape == (350,) and len(np.unique(ids)) == 350
        assert np.all((ids >= 0) & (ids < INTERNAL))
        assert re.search(r"flipMap\s+List<bool>\s+350\{false\};", entry)
        coords = faces[ids]
        area_vectors = np.cross(coords[:,1]-coords[:,0], coords[:,2]-coords[:,0])
        assert np.all(area_vectors[:,0] > 0) and np.ptp(coords[:,:,0]) < 1e-10
        zone_faces.append(ids)
        zone_locations.append(float(coords[0,0,0]))

    height, froude, flow = (np.loadtxt(case/name) for name in [
        "postProcessing/0/flumeMonitor/height", "postProcessing/0/flumeMonitor/Fr",
        "postProcessing/0/flumeMonitor/flowRate"])
    sources.extend(f"postProcessing/0/flumeMonitor/{name}" for name in ("height", "Fr", "flowRate"))
    assert height.shape == froude.shape == (9,5) and flow.shape == (9,6)
    for series in (height, froude, flow):
        assert np.all(np.isfinite(series)) and np.array_equal(series[:,0], np.r_[0,SAVES])
    assert np.all(height[:,1:] > 0) and np.all(froude[:,1:] >= 0)
    assert np.allclose(6.004*height[:,2]**1.605, flow[:,5], rtol=1e-6, atol=1e-8)
    prediction = abs(flow[:,1:5])/(WIDTHS*height[:,1:]*np.sqrt(9.81*height[:,1:]))
    assert np.allclose(prediction, froude[:,1:], rtol=1e-6, atol=1e-8)
    patches = []
    first = INTERNAL
    for name, contents in re.findall(r"(\w+)\s*\{([^{}]*nFaces[^{}]*)\}",
                                     native("constant/polyMesh/boundary"), re.S):
        count = int(re.search(r"\bnFaces\s+(\d+);", contents)[1])
        start = int(re.search(r"\bstartFace\s+(\d+);", contents)[1])
        assert first == start
        first += count
        patches.append((name, start, count))
    assert first == FACES

    results, sections = [], []
    for idx, time in enumerate(SAVES,1):
        clock = native(f"{time}/uniform/time")
        step = int(re.search(r"^index\s+(\d+);", clock, re.M)[1])
        assert float(re.search(r"^value\s+([^;]+);", clock, re.M)[1]) == time
        assert steps[step-1] == time
        alpha = field(f"{time}/alpha.water")
        field(f"{time}/U", kind="vector")
        field(f"{time}/p_rgh")
        sections.append(alpha[adjacent].mean(axis=1))

        phase = field(f"{time}/alphaPhi.water", INTERNAL)
        measured = np.array([phase[ids].sum() for ids in zone_faces])
        assert np.allclose(measured, flow[idx,1:5], rtol=1e-6, atol=1e-8)
        data = native(f"{time}/phi")
        flux = np.empty(FACES)
        flux[:INTERNAL] = read(data,"internalField",INTERNAL,"scalar")
        for name, start, count in patches:
            block = re.search(r"\b"+name+r"\s*\{([^{}]+)\}", data, re.S)[1]
            flux[start:start+count] = read(block,"value",count,"scalar")
        assert np.all(np.isfinite(flux))
        imbalance = (np.bincount(owner, weights=flux, minlength=CELLS) -
                     np.bincount(neighbour, weights=flux[:INTERNAL], minlength=CELLS))
        results.append({"time_s":int(time), "step":step, "alpha_min":float(alpha.min()),
                        "alpha_max":float(alpha.max()), "above_alpha_gate":int(np.count_nonzero(alpha>1+1e-6)),
                        "upstream_depth_m":float(height[idx,2]), "downstream_depth_m":float(height[idx,3]),
                        "upstream_flux_m3_s":float(measured[1]), "downstream_flux_m3_s":float(measured[2]),
                        "empirical_flux_m3_s":float(flow[idx,5]),
                        "upstream_froude":float(froude[idx,2]), "downstream_froude":float(froude[idx,3]),
                        "cell_flux_abs_sum_m3_s":float(abs(imbalance).sum()),
                        "max_abs_flux_divergence_s_inv":float(max(abs(imbalance/volumes)))})
    assert results[-1]["step"] == len(steps)
    fraction_ok = all(r["alpha_min"] >= -1e-6 and r["alpha_max"] <= 1+1e-6 for r in results)
    flow_ok = min(results[-1]["upstream_flux_m3_s"], results[-1]["downstream_flux_m3_s"]) > .1
    checks = {"pass":fraction_ok and flow_ok, "alpha_pass":fraction_ok, "terminal_flow_pass":flow_ok,
              "cells":CELLS, "native_steps":len(steps), "saved_states":len(results),
              "zone_plane_locations_m":dict(zip(ZONES,zone_locations)), "courant":courant,
              "samples":results, "native_sha256":{name:sha(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(checks,indent=2)+"\n", encoding="utf-8")
    with (out/"flume-history.csv").open("x",newline="") as file:
        writer = csv.DictWriter(file,fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    fig, axes = plt.subplots(2, 2, figsize=(11,6.8), layout="constrained")
    for ax, section, time in ((axes[0,0],sections[0],5),(axes[1,0],sections[-1],40)):
        pc = PolyCollection(polygons,array=section,cmap="Blues",clim=(0,1),edgecolors="none")
        ax.add_collection(pc)
        for position in zone_locations[1:3]:
            ax.axvline(position,color="#c14f27",lw=1,ls="--")
        ax.set(xlim=(0,15.86),ylim=(-.18,1.05),xlabel="Distance (m)",ylabel="Elevation (m)",
               title=f"Centre-plane water fraction at {time} s")
        fig.colorbar(pc,ax=ax,label="Water fraction")
    for i,label in ((2,"Upstream"),(3,"Downstream")):
        axes[0,1].plot(height[:,0],height[:,i],"o-",ms=3,label=label)
        axes[1,1].plot(flow[:,0],flow[:,i],"o-",ms=3,label=label)
    axes[1,1].plot(flow[:,0],flow[:,5],"--",label="Head correlation")
    axes[0,1].set(xlabel="Time (s)",ylabel="Depth (m)",title="Section-equivalent water depth")
    axes[1,1].set(xlabel="Time (s)",ylabel="Water flow (m³/s)",title="Signed section discharge")
    for ax in axes[:,1]:
        ax.grid(alpha=.2)
        ax.legend(frameon=False)
    fig.suptitle("Free-surface flow through a Parshall flume",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"parshall-flume.{ext}",dpi=180,
                    metadata={"CreationDate": None, "ModDate": None} if ext == "pdf" else None)
    plt.close(fig)
    return checks["pass"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("case","reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)