"""Read the native tank discharge fields and plot early blowdown."""
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

CELLS = 3725
STEPS = 3000
DT = 5e-7


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def geometry(case, centres, volumes):
    mesh = case/"constant/polyMesh"
    points = np.array([np.fromstring(s, sep=" ") for s in
                       re.findall(r"^\(([^()]+)\)$", (mesh/"points").read_text("ascii"), re.M)])
    face_ids = np.array([np.fromstring(s, dtype=int, sep=" ") for s in
                         re.findall(r"^4\(([^()]+)\)$", (mesh/"faces").read_text("ascii"), re.M)])
    owner_list = re.search(r"\n15075\s*\(\s*([^)]*)\)", (mesh/"owner").read_text("ascii"))[1]
    owners = np.fromstring(owner_list, dtype=int, sep=" ")
    assert points.shape == (7802, 3) and face_ids.shape == (15075, 4) and owners.shape == (15075,)
    assert np.all((0 <= face_ids) & (face_ids < len(points)))
    assert np.all((0 <= owners) & (owners < CELLS))
    faces = points[face_ids]
    front = np.all(np.isclose(faces[:, :, 2], -.01, atol=1e-12, rtol=0), axis=1)
    order = np.argsort(owners[front])
    assert np.array_equal(owners[front][order], np.arange(CELLS))
    polygons = faces[front][order, :, :2]
    x, y = polygons.transpose(2, 0, 1)
    xx, yy = np.roll(x, -1, axis=1), np.roll(y, -1, axis=1)
    cross = x*yy-y*xx
    area2 = cross.sum(axis=1)
    centroid = np.column_stack((((x+xx)*cross).sum(axis=1), ((y+yy)*cross).sum(axis=1)))/(3*area2[:, None])
    assert np.allclose(centroid, centres[:, :2], rtol=0, atol=1e-9)
    assert np.allclose(.01*abs(area2), volumes, rtol=1e-8, atol=1e-13)

    patch = re.search(r"\bnozzle\s*\{([^}]+)\}", (mesh/"boundary").read_text("ascii"))[1]
    assert int(re.search(r"\bnFaces\s+(\d+);", patch)[1]) == 5
    first = int(re.search(r"\bstartFace\s+(\d+);", patch)[1])
    nozzle = faces[first:first+5]
    assert np.allclose(nozzle[:, :, 0], .15, rtol=0, atol=1e-12)
    areas = np.ptp(nozzle[:, :, 1], axis=1)*np.ptp(nozzle[:, :, 2], axis=1)
    assert np.all(areas > 0) and np.isclose(areas.sum(), .0002, rtol=0, atol=1e-12)
    return polygons, owners[first:first+5], areas


def audit(case, reader_path, out):
    read = runpy.run_path(str(reader_path))["foam_values"]
    sources = ["input-sha256.txt", "tutorial-sha256.txt", "tutorial-inputs-verified.log",
               "check-mesh.log", "solver.log", "constant/physicalProperties"]
    for line in (case/sources[0]).read_text("ascii").splitlines():
        expected, name = line.split("  ", 1)
        assert digest(case/name) == expected, name
    mesh = (case/"check-mesh.log").read_text("ascii")
    log = (case/"solver.log").read_text("ascii")
    times = np.array(re.findall(r"^Time = ([-+\d.eE]+)s$", log, re.M), dtype=float)
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*3725\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd") and len(times) == STEPS
    assert np.allclose(times, np.arange(1, STEPS+1)*DT, rtol=0, atol=1e-12)

    def field(name, kind="scalar"):
        sources.append(name)
        result = read((case/name).read_text("ascii"), "internalField", CELLS, kind)
        assert np.all(np.isfinite(result)), name
        return result

    centres, volumes = field("0/C", "vector"), field("0/Vc")
    assert np.all(volumes > 0)
    polygons, nozzle, areas = geometry(case, centres, volumes)
    sources.extend(f"constant/polyMesh/{name}" for name in ("points", "faces", "owner", "boundary"))
    thermo = (case/"constant/physicalProperties").read_text("ascii")
    assert re.search(r"equationOfState\s+perfectGas;", thermo)
    gas_r = 8314.46261815324/float(re.search(r"molWeight\s+([^;]+);", thermo)[1])
    rows = []
    for index in range(16):
        time = index/10000
        directory = f"{time:.4f}".rstrip("0").rstrip(".")
        if index:
            clock_path = f"{directory}/uniform/time"
            sources.append(clock_path)
            clock = (case/clock_path).read_text("ascii")
            assert int(re.search(r"^index\s+(\d+);", clock, re.M)[1]) == index*200
            assert abs(float(re.search(r"^value\s+([^;]+);", clock, re.M)[1])-time) < 1e-12
            assert float(re.search(r"^deltaT\s+([^;]+);", clock, re.M)[1]) == DT
        u, p, temp = field(f"{directory}/U", "vector"), field(f"{directory}/p"), field(f"{directory}/T")
        assert u.shape == (CELLS, 3) and p.shape == temp.shape == (CELLS,)
        assert np.all(p > 0) and np.all(temp > 0)
        rho = p/(gas_r*temp)
        mismatch = None
        flux = 0.
        if index:
            native_rho = field(f"{directory}/rho")
            assert np.all(native_rho > 0)
            mismatch = float(np.max(abs(native_rho/rho-1)))
            flux_file = f"{directory}/phi"
            sources.append(flux_file)
            phi = (case/flux_file).read_text("ascii")
            assert re.search(r"dimensions\s+\[1 0 -1 0 0 0 0\];", phi)
            section = re.search(r"\bnozzle\s*\{([^}]+)\}", phi)[1]
            fluxes = read(section, "value", 5, "scalar")
            assert fluxes.shape == (5,) and np.all(np.isfinite(fluxes))
            flux = float(fluxes.sum())
        else:
            assert np.all(p == 200000) and np.all(temp == 300) and np.all(u == 0)
        rows.append({"time_s": time, "bulk_pressure_pa": float(np.average(p, weights=volumes)),
                     "inventory_kg": float(rho@volumes),
                     "nozzle_adjacent_speed_m_s": float(np.average(u[nozzle, 0], weights=areas)),
                     "native_nozzle_mass_flux_kg_s": flux, "rho_crosscheck_max_relative": mismatch,
                     "min_pressure_pa": float(p.min()), "min_temperature_k": float(temp.min())})

    pressure_drop = 1-rows[-1]["bulk_pressure_pa"]/rows[0]["bulk_pressure_pa"]
    mass_drop = 1-rows[-1]["inventory_kg"]/rows[0]["inventory_kg"]
    assert pressure_drop > .01 and mass_drop > .01 and rows[-1]["nozzle_adjacent_speed_m_s"] > 0
    result = {"cells": CELLS, "native_steps": STEPS, "saved_states": 15, "time_step_s": DT,
              "pressure_reduction_fraction": pressure_drop, "inventory_reduction_fraction": mass_drop,
              "nozzle_face_area_m2": float(areas.sum()), "samples": rows,
              "native_sha256": {name: digest(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    with (out/"blowdown-history.csv").open("x", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.4), layout="constrained")
    speed = np.linalg.norm(u, axis=1)
    for ax, values, label in ((axes[0], p/1000, "Absolute pressure (kPa)"),
                               (axes[1], speed, "Speed (m/s)")):
        collection = PolyCollection(polygons*1000, array=values, cmap="viridis", edgecolors="none")
        ax.add_collection(collection)
        ax.set(xlim=(0, 150), ylim=(0, 300), aspect="equal", xlabel="x (mm)", ylabel="y (mm)",
               title=f"{label} at 1.5 ms")
        fig.colorbar(collection, ax=ax, shrink=.78)
    for key, label in (("bulk_pressure_pa", "Mean pressure"), ("inventory_kg", "Gas mass")):
        axes[2].plot([row["time_s"]*1000 for row in rows],
                     [100*row[key]/rows[0][key] for row in rows], "o-", ms=3, label=label)
    axes[2].set(xlabel="Time (ms)", ylabel="Percent of initial value", title="Tank response")
    axes[2].grid(alpha=.2)
    axes[2].legend(frameon=False)
    fig.suptitle("Early ideal-gas tank blowdown", fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(out/f"gas-tank-blowdown.{ext}", dpi=180,
                    metadata={"CreationDate": None, "ModDate": None} if ext == "pdf" else None)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for argument in ("case", "reader", "output"):
        parser.add_argument(argument, type=Path)
    args = parser.parse_args()
    audit(args.case, args.reader, args.output)