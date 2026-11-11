"""Audit the Euler–Euler column fields and phase slip."""
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

NX, NY = 25, 75
CELLS = NX * NY
DX, DY, DEPTH = .15 / NX, 1 / NY, .1
TIMES = np.arange(1, 11) / 2
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def audit(case, reader, out):
    read = runpy.run_path(str(reader))["foam_values"]
    sources = {}

    def native(name):
        raw = (case / name).read_bytes()
        sources[name] = hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest, name = line.split("  ", 1)
        assert hashlib.sha256((case / name).read_bytes()).hexdigest() == digest, name
        sources[name] = digest
    for name in ("pre-mesh-input-sha256.txt", "pre-mesh-inputs-verified.log",
                 "inputs-verified.log", "mesh.log", "centres.log", "volumes.log"):
        native(name)
    mesh, log = native("check-mesh.log"), native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*1875\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd") and "multiphaseEuler" in log
    times = np.array(re.findall(r"^Time = (" + NUMBER + r")s?$", log, re.M), float)
    np.testing.assert_allclose(times, np.arange(1, 1001) / 200, atol=1e-10, rtol=0)

    def field(name, kind="scalar"):
        values = read(native(name), "internalField", CELLS, kind)
        assert values.shape == ((CELLS, 3) if kind == "vector" else (CELLS,)), name
        assert np.isfinite(values).all(), name
        return values

    c, v = field("0/C", "vector"), field("0/Vc")
    np.testing.assert_allclose(v, DX*DY*DEPTH, atol=1e-14, rtol=0)
    np.testing.assert_allclose(c[:, 2], DEPTH/2, atol=1e-10, rtol=0)
    ix, iy = np.floor(c[:, :2] / [DX, DY]).astype(int).T
    assert np.all((ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY))
    assert len(np.unique(iy*NX + ix)) == CELLS
    np.testing.assert_allclose(c[:, :2], np.column_stack(((ix+.5)*DX, (iy+.5)*DY)), atol=1e-9, rtol=0)
    init = native("0/alpha.air")
    assert re.search(r"box\s+\(0 0 -0\.1\)\s+\(0\.15 0\.701 0\.1\);", init)
    assert re.search(r"defaultValue\s+1;", init)
    initial_air = (c[:, 1] > .701).astype(float)
    lower = c[:, 1] < .5
    assert not np.any(initial_air[lower])

    def grid(values):
        result = np.empty((NY, NX))
        result[iy, ix] = values
        return result

    actual = {p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?", p.name)}
    assert actual == {"0", *(f"{t:g}" for t in TIMES)}
    rows = []
    planes = {"0": grid(initial_air)}
    for t in TIMES:
        name = f"{t:g}"
        clock = native(f"{name}/uniform/time")
        step = int(re.search(r"^index\s+(\d+);", clock, re.M)[1])
        assert 1 <= step <= len(times) and abs(times[step-1]-t) < 1e-10
        assert abs(float(re.search(r"^value\s+([^;]+);", clock, re.M)[1])-t) < 1e-10
        a, w = field(f"{name}/alpha.air"), field(f"{name}/alpha.water")
        ua, uw = field(f"{name}/U.air", "vector"), field(f"{name}/U.water", "vector")
        p = field(f"{name}/p")
        ta, tw = field(f"{name}/T.air"), field(f"{name}/T.water")
        slip = ua[:, 1]-uw[:, 1]
        active = lower & (a >= .01)
        rows.append({"time_s": float(t), "step": step,
                     "phase_fraction_min": float(min(a.min(), w.min())),
                     "phase_fraction_max": float(max(a.max(), w.max())),
                     "maximum_phase_sum_error": float(np.max(abs(a+w-1))),
                     "absolute_pressure_min_pa": float(p.min()),
                     "min_temperature_k": float(min(ta.min(), tw.min())),
                     "max_temperature_k": float(max(ta.max(), tw.max())),
                     "lower_column_air_volume_fraction": float((a[lower]@v[lower])/v[lower].sum()),
                     "lower_column_cells_with_one_percent_air": int(active.sum()),
                     "mean_upward_slip_m_s": float(slip[active].mean()) if active.any() else None,
                     "max_air_speed_m_s": float(np.linalg.norm(ua, axis=1).max()),
                     "max_water_speed_m_s": float(np.linalg.norm(uw, axis=1).max())})
        if name in {"2.5", "5"}:
            planes[name] = grid(a)
    assert rows[-1]["step"] == len(times)
    terminal = rows[-1]
    gates = {"phase_fraction": all(-1e-6 <= r["phase_fraction_min"] <= r["phase_fraction_max"] <= 1+1e-6
                                  for r in rows),
             "phase_sum": max(r["maximum_phase_sum_error"] for r in rows) <= 2e-6,
             "thermodynamic_fields": all(r["absolute_pressure_min_pa"] > 0 and
                                        250 <= r["min_temperature_k"] <= r["max_temperature_k"] <= 375
                                        for r in rows),
             "lower_column_air": terminal["lower_column_cells_with_one_percent_air"] >= 10,
             "phase_slip": terminal["mean_upward_slip_m_s"] is not None and terminal["mean_upward_slip_m_s"] > .02}
    result = {"pass": all(gates.values()), "gates": gates, "cells": CELLS, "steps": len(times),
              "saved_states": len(rows), "samples": rows, "native_sha256": sources}
    out.mkdir()
    (out / "checks.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with (out / "bubble-column-history.csv").open("x", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig = plt.figure(figsize=(9.8, 6.4), layout="constrained")
    spec = fig.add_gridspec(2, 2, width_ratios=[1.35, 1])
    plane_spec = spec[:, 0].subgridspec(1, 3)
    axs = [fig.add_subplot(plane_spec[0, i]) for i in range(3)]
    for ax, (t, a) in zip(axs, planes.items()):
        im = ax.imshow(a, origin="lower", extent=(0, .15, 0, 1), cmap="viridis",
                       vmin=0, vmax=1, interpolation="nearest", aspect="equal")
        ax.axhline(.5, color="white", linestyle="--", lw=.9)
        ax.set(title=f"{t} s", xlabel="x (m)", xticks=[0,.15])
    axs[0].set_ylabel("Height (m)")
    for ax in axs[1:]:
        ax.set_yticklabels([])
    fig.colorbar(im, ax=axs, shrink=.85, label="Air volume fraction")
    ax = fig.add_subplot(spec[0, 1])
    ax.plot(np.r_[0, TIMES], [0]+[100*r["lower_column_air_volume_fraction"] for r in rows], "o-",
            color="#267466", ms=3)
    ax.set(title="Air in the lower 0.5 m", xlabel="Time (s)", ylabel="Mean air fraction (%)")
    ax.grid(alpha=.2)
    ax = fig.add_subplot(spec[1, 1])
    ax.plot(TIMES, [r["mean_upward_slip_m_s"] for r in rows], "o-", color="#245d83", ms=3)
    ax.set(title="Mean phase slip in aerated cells", xlabel="Time (s)", ylabel="Upward air minus water speed (m/s)")
    ax.grid(alpha=.2)
    fig.suptitle("Air injection into a bubble column", fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"bubble-column.{ext}", dpi=180,
                    metadata={"CreationDate": None, "ModDate": None} if ext == "pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "reader", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    raise SystemExit(0 if audit(args.case, args.reader, args.output) else 2)
