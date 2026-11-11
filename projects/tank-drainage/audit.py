"""Audit saved liquid drainage fields and integrate tank and collector inventories."""
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
import numpy as np

CELLS = 7200
TIMES = [i/20 for i in range(1, 13)]
ALPHA_TOL = 1e-6
NUMBER = r"[-+0-9.eE]+"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, foam_reader, output):
    reader = runpy.run_path(str(foam_reader))["foam_values"]

    def field(path, kind="scalar"):
        return reader((case/path).read_text(encoding="ascii"), "internalField", CELLS, kind)

    log = (case/"solver.log").read_text(encoding="ascii")
    mesh = (case/"check-mesh.log").read_text(encoding="ascii")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*7200\s*$", mesh, re.M)
    assert re.search(r"^End\s*$", log, re.M)
    for line in (case/"input-sha256.txt").read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ", 1)
        assert digest(case/name) == expected, name
    steps = np.array(re.findall(r"^Time = ("+NUMBER+r")s$", log, re.M), dtype=float)
    bulk = np.array(re.findall(r"^Courant Number mean: "+NUMBER+r" max: ("+NUMBER+r")", log, re.M), dtype=float)
    interface = np.array(re.findall(r"^Interface Courant Number mean: "+NUMBER+r" max: ("+NUMBER+r")", log, re.M), dtype=float)
    assert len(steps) > 0 and abs(steps[-1] - .6) < 1e-10 and np.all(np.diff(steps) > 0)
    assert len(bulk) >= len(steps) and len(interface) >= len(steps)
    assert np.isfinite(bulk).all() and np.isfinite(interface).all()

    centres = field("0/C", "vector")
    volumes = field("0/Vc")
    assert np.isfinite(centres).all() and np.isfinite(volumes).all() and np.all(volumes > 0)
    x, _, z = centres.T
    upper = z > 0
    initial = volumes[(z > 0) & (z < .05)].sum()
    assert abs(initial - .0001) < 1e-9
    xs, ix = np.unique(np.round(x, 8), return_inverse=True)
    zs, iz = np.unique(np.round(z, 8), return_inverse=True)
    nx, nz = len(xs), len(zs)
    assert (nx, nz) == (60, 120)
    grid = iz*nx + ix
    order = np.argsort(grid)
    assert np.array_equal(grid[order], np.arange(CELLS))
    grid_volumes = volumes[order].reshape(nz, nx)
    z_edges = np.linspace(-.1, .1, nz+1)
    assert np.allclose(zs, (z_edges[1:] + z_edges[:-1])/2, rtol=0, atol=1e-7)
    widths = grid_volumes.mean(axis=0)/(.02*(z_edges[1] - z_edges[0]))
    x_edges = np.r_[-.05, -.05 + np.cumsum(widths)]
    assert np.allclose(xs, (x_edges[1:] + x_edges[:-1])/2, rtol=0, atol=1e-7)

    rows = []
    snapshots = {}
    sources = ["input-sha256.txt", "solver.log", "check-mesh.log", "0/C", "0/Vc"]
    for time in TIMES:
        label = f"{time:.10g}"
        directory = case/label
        clock = (directory/"uniform/time").read_text(encoding="ascii")
        step = int(re.search(r"^index\s+(\d+);", clock, re.M)[1])
        assert abs(float(re.search(r"^value\s+([^;]+);", clock, re.M)[1]) - time) < 1e-8
        assert abs(steps[step-1] - time) < 1e-8
        alpha = field(f"{label}/alpha.liquid")
        assert np.isfinite(alpha).all() and alpha.min() >= -ALPHA_TOL and alpha.max() <= 1+ALPHA_TOL
        assert np.isfinite(field(f"{label}/U", "vector")).all()
        assert np.isfinite(field(f"{label}/p_rgh")).all()
        inventory = alpha*volumes
        tank, collector = inventory[upper].sum(), inventory[~upper].sum()
        rows.append({"time_s": time, "native_step": step, "tank_liquid_m3": float(tank),
                     "collector_liquid_m3": float(collector), "domain_liquid_m3": float(tank+collector),
                     "tank_net_loss_m3": float(initial-tank), "alpha_min": float(alpha.min()),
                     "alpha_max": float(alpha.max())})
        if label in {"0.15", "0.6"}:
            snapshots[label] = alpha[order].reshape(nz, nx)
        sources.extend(f"{label}/{name}" for name in ("uniform/time","alpha.liquid","U","p_rgh"))
    assert rows[-1]["native_step"] == len(steps)
    assert rows[-1]["tank_liquid_m3"] < .70*initial
    assert max(row["collector_liquid_m3"] for row in rows) >= .05*initial
    assert max(row["domain_liquid_m3"] for row in rows) <= 1.02*initial

    checks = {"cells": CELLS, "initial_tank_liquid_m3": float(initial), "native_steps": len(steps),
              "saved_states": len(rows), "max_logged_bulk_courant": float(bulk.max()),
              "max_logged_interface_courant": float(interface.max()), "samples": rows,
              "native_sha256": {name: digest(case/name) for name in sources}}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(checks, indent=2)+"\n", encoding="utf-8")
    with (output/"drainage-history.csv").open("x", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fig, (early, late, history) = plt.subplots(1, 3, figsize=(11.5, 4.5),
                                              gridspec_kw={"width_ratios": [1, 1, 2.4]},
                                              layout="constrained")
    for ax, label in ((early,"0.15"), (late,"0.6")):
        panel = ax.pcolormesh(1000*x_edges, 1000*z_edges, snapshots[label], vmin=0, vmax=1, cmap="viridis")
        ax.set(xlabel="Width coordinate (mm)", title=f"{float(label):.2f} s", aspect="equal")
        for left, right in ((-50,-10),(10,50)):
            ax.plot([left,right], [0,0], color="white", linewidth=2)
    early.set_ylabel("Elevation (mm)")
    fig.colorbar(panel, ax=[early,late], label="Liquid volume fraction", shrink=.82)
    for key, label in (("tank_liquid_m3","Upper tank"), ("collector_liquid_m3","Lower region"),
                       ("domain_liquid_m3","Whole domain")):
        first = 0 if key == "collector_liquid_m3" else initial
        history.plot([0]+TIMES, 1e6*np.r_[first,[row[key] for row in rows]], "o-", label=label)
    history.set(xlabel="Physical time (s)", ylabel="Liquid inventory (mL)", title="Liquid transfer through the slot")
    history.grid(alpha=.25)
    history.legend()
    fig.suptitle("Gravity drainage from a slotted tank | OpenFOAM", fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(output/f"tank-drainage.{ext}", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "foam_reader", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    audit(args.case, args.foam_reader, args.output)
