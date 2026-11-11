"""Read the native panel mesh, film thickness, fluxes and saved clocks."""
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
from matplotlib.colors import LogNorm
import numpy as np

NX, NY = 180, 240
CELLS = NX * NY
DX = DY = 1 / 240
HEIGHT = .01
WET = 1e-5
NEGATIVE_TOLERANCE_M = 1e-12
TIMES = np.arange(1, 11) / 10
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
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*43200\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd") and "isothermalFilm" in log
    step_lines = list(re.finditer(r"^Time = (" + NUMBER + r")s?$", log, re.M))
    times = np.array([float(m[1]) for m in step_lines])
    assert len(times) and np.isfinite(times).all() and (np.diff(times) > 0).all() and times[-1] == 1
    max_terminal_alpha_residual = 0.
    for i, step in enumerate(step_lines):
        next_step = step_lines[i + 1].start() if i + 1 < len(step_lines) else len(log)
        residuals = re.findall(r"Solving for alpha,\s+Initial residual = " + NUMBER
                               + r", Final residual = (" + NUMBER + ")", log[step.end():next_step])
        assert residuals and all(np.isfinite(float(r)) and float(r) >= 0 for r in residuals)
        max_terminal_alpha_residual = max(max_terminal_alpha_residual, float(residuals[-1]))

    def field(name, kind="scalar"):
        a = read(native(name), "internalField", CELLS, kind)
        assert a.shape == ((CELLS, 3) if kind == "vector" else (CELLS,)), name
        assert np.isfinite(a).all(), name
        return a

    c, v = field("0/C", "vector"), field("0/Vc")
    np.testing.assert_allclose(v, DX * DY * HEIGHT, atol=1e-14, rtol=0)
    np.testing.assert_allclose(c[:, 2], .005, atol=1e-9, rtol=0)
    ix, iy = np.floor(c[:, :2] / [DX, DY]).astype(int).T
    assert np.all((ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY))
    assert len(np.unique(iy * NX + ix)) == CELLS
    np.testing.assert_allclose(c[:, :2], np.column_stack(((ix + .5) * DX, (iy + .5) * DY)), atol=1e-9, rtol=0)
    assert np.all(field("0/delta") == 0)
    boundary = native("constant/polyMesh/boundary")
    patch_sizes = {}
    for patch, size in (("inlet", 120), ("outlet", 180), ("wall", CELLS), ("empty", CELLS)):
        body = re.search(r"\b" + patch + r"\s*\{([^}]+)\}", boundary, re.S)[1]
        patch_sizes[patch] = int(re.search(r"\bnFaces\s+(\d+);", body)[1])
        assert patch_sizes[patch] == size
        if patch == "wall":
            assert re.search(r"\btype\s+filmWall;", body)
    actual = {p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?", p.name)}
    assert actual == {"0", *(f"{t:g}" for t in TIMES)}

    rows = []
    fields = {}
    for t in TIMES:
        name = f"{t:g}"
        clock = native(f"{name}/uniform/time")
        step = int(re.search(r"^index\s+(\d+);", clock, re.M)[1])
        assert 1 <= step <= len(times) and abs(times[step - 1] - t) < 1e-9
        assert abs(float(re.search(r"^value\s+([^;]+);", clock, re.M)[1]) - t) < 1e-10
        delta, alpha = field(f"{name}/delta"), field(f"{name}/alpha")
        u, p = field(f"{name}/U", "vector"), field(f"{name}/p")
        flux = native(f"{name}/alphaRhoPhi")
        assert re.search(r"dimensions\s+\[1\s+0\s+-1\s+0\s+0\s+0\s+0\];", flux)
        sums = {}
        for patch in ("inlet", "outlet"):
            body = re.search(r"\b" + patch + r"\s*\{([^}]+)\}", flux, re.S)[1]
            values = read(body, "value", patch_sizes[patch], "scalar")
            assert values.shape == (patch_sizes[patch],) and np.isfinite(values).all()
            sums[f"{patch}_outward_mass_flux_kg_s"] = float(values.sum())
        wet = delta >= WET
        rows.append({"time_s": float(t), "step": step,
                     "min_thickness_m": float(delta.min()), "max_thickness_m": float(delta.max()),
                     "film_volume_m3": float(delta @ (v / HEIGHT)),
                     "negative_thickness_cells": int(np.count_nonzero(delta < 0)),
                     "negative_thickness_volume_m3": float(np.minimum(delta,0) @ (v / HEIGHT)),
                     "wet_area_m2": float(v[wet].sum() / HEIGHT),
                     "wet_front_distance_m": float(1 - c[wet, 1].min()) if wet.any() else 0.,
                     "alpha_to_delta_max_error_m": float(np.max(abs(HEIGHT * alpha - delta))),
                     "maximum_speed_m_s": float(np.linalg.norm(u, axis=1).max()),
                     "minimum_pressure_pa": float(p.min()), **sums})
        if name in {"0.1", "0.5", "1"}:
            grid = np.empty((NY, NX))
            grid[iy, ix] = delta
            fields[name] = grid
    assert rows[-1]["step"] == len(times)
    gates = {"thickness": all(-NEGATIVE_TOLERANCE_M <= r["min_thickness_m"] <= r["max_thickness_m"] < HEIGHT for r in rows),
             "conversion": max(r["alpha_to_delta_max_error_m"] for r in rows) <= 1e-9,
             "inlet": all(r["inlet_outward_mass_flux_kg_s"] < 0 for r in rows),
             "volume": rows[-1]["film_volume_m3"] >= 1e-6,
             "front": rows[-1]["wet_front_distance_m"] >= .1 and
                      rows[-1]["wet_front_distance_m"] > rows[0]["wet_front_distance_m"]}
    result = {"pass": all(gates.values()), "gates": gates, "cells": CELLS, "steps": len(times),
              "saved_states": len(rows), "wet_threshold_m": WET,
              "negative_thickness_tolerance_m": NEGATIVE_TOLERANCE_M,
              "strict_nonnegative_at_all_saves": all(r["min_thickness_m"] >= 0 for r in rows),
              "maximum_terminal_alpha_residual": max_terminal_alpha_residual,
              "samples": rows, "native_sha256": sources}
    out.mkdir()
    (out / "checks.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    with (out / "film-history.csv").open("x", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axs = plt.subplots(2, 3, figsize=(10.5, 6.5), layout="constrained",
                           gridspec_kw={"height_ratios": [1.65, 1]})
    max_thickness = max(r["max_thickness_m"] for r in rows) * 1000
    for ax, (time, thickness) in zip(axs[0], fields.items()):
        visible = np.ma.masked_less(thickness * 1000, WET * 1000)
        im = ax.imshow(visible, origin="lower", extent=(0, .75, 0, 1),
                       cmap="viridis", norm=LogNorm(vmin=WET * 1000, vmax=max_thickness), interpolation="nearest")
        ax.set_facecolor("#e6e9ee")
        ax.plot([.125, .625], [1, 1], color="#a95c3a", lw=2.5, clip_on=False)
        ax.set(title=f"{time} s", xlabel="x (m)", ylabel="Height (m)", aspect="equal")
    fig.colorbar(im, ax=axs[0], label="Film thickness ≥ 0.01 mm (log scale)", shrink=.84)
    for ax, key, scale, label, title in zip(axs[1],
            ("film_volume_m3", "wet_area_m2", "wet_front_distance_m"), (1e6, 1, 1000),
            ("Volume (mL)", "Area (m²)", "Distance (mm)"),
            ("Film volume", "Wetted area", "Wet front below inlet")):
        ax.plot(np.r_[0, TIMES], [0] + [scale * r[key] for r in rows], "o-", color="#245d83", ms=3)
        ax.set(xlabel="Time (s)", ylabel=label, title=title)
        ax.grid(alpha=.2)
    fig.suptitle("Water-film spreading on a vertical panel", fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"rivulet-panel.{ext}", dpi=180,
                    metadata={"CreationDate": None, "ModDate": None} if ext == "pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "reader", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    raise SystemExit(0 if audit(args.case, args.reader, args.output) else 2)
