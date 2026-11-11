"""Check the periodic cube and analyse saved velocity and enstrophy."""
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

N = 16
CELLS = N**3
DX = 1 / N
TIMES = np.arange(1, 41) / 4
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
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*4096\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd") and "Selecting OUForce" in log
    step_lines = list(re.finditer(r"^Time = (" + NUMBER + r")s?$", log, re.M))
    times = np.array([float(match[1]) for match in step_lines])
    np.testing.assert_allclose(times, np.arange(1, 401) / 40, atol=1e-10, rtol=0)
    final_residuals = []
    for i, step in enumerate(step_lines):
        next_step = step_lines[i + 1].start() if i + 1 < len(step_lines) else len(log)
        residuals = re.findall(r"Solving for p,\s+Initial residual = " + NUMBER
                               + r", Final residual = (" + NUMBER + ")", log[step.end():next_step])
        assert len(residuals) == 2
        final_residuals.append(float(residuals[-1]))
    assert np.isfinite(final_residuals).all()

    def field(name, kind="scalar"):
        a = read(native(name), "internalField", CELLS, kind)
        assert a.shape == ((CELLS, 3) if kind == "vector" else (CELLS,)), name
        assert np.isfinite(a).all(), name
        return a

    c, v = field("0/C", "vector"), field("0/Vc")
    np.testing.assert_allclose(v, DX**3, atol=1e-14, rtol=0)
    ijk = np.floor(c / DX).astype(int)
    assert np.all((ijk >= 0) & (ijk < N)) and len(np.unique(ijk, axis=0)) == CELLS
    np.testing.assert_allclose(c, (ijk + .5) * DX, atol=1e-10, rtol=0)
    assert np.all(field("0/U", "vector") == 0)

    boundary = native("constant/polyMesh/boundary")
    for pair in range(3):
        for side in range(2):
            patch = f"patch{pair}_half{side}"
            peer = f"patch{pair}_half{1-side}"
            body = re.search(r"\b" + patch + r"\s*\{([^}]+)\}", boundary, re.S)[1]
            assert re.search(r"\btype\s+cyclic;", body)
            assert re.search(r"\bneighbourPatch\s+" + peer + ";", body)
            assert int(re.search(r"\bnFaces\s+(\d+);", body)[1]) == N**2

    actual = {p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?", p.name)}
    assert actual == {"0", *(f"{t:g}" for t in TIMES)}
    rows = []
    planes = {}
    for t in TIMES:
        name = f"{t:g}"
        clock = native(f"{name}/uniform/time")
        step = int(re.search(r"^index\s+(\d+);", clock, re.M)[1])
        assert 1 <= step <= len(times) and abs(times[step - 1] - t) < 1e-10
        assert abs(float(re.search(r"^value\s+([^;]+);", clock, re.M)[1]) - t) < 1e-10
        u, p, e = field(f"{name}/U", "vector"), field(f"{name}/p"), field(f"{name}/enstrophy")
        speed = np.linalg.norm(u, axis=1)
        rows.append({"time_s": float(t), "step": step,
                     "mean_velocity_norm_m_s": float(np.linalg.norm(u.mean(axis=0))),
                     "mean_kinetic_energy_m2_s2": float(.5 * np.mean(speed**2)),
                     "mean_enstrophy_s2_inverse": float(e.mean()),
                     "minimum_enstrophy_s2_inverse": float(e.min()),
                     "maximum_speed_m_s": float(speed.max()),
                     "minimum_kinematic_pressure_m2_s2": float(p.min()),
                     "maximum_kinematic_pressure_m2_s2": float(p.max()),
                     "final_pressure_residual": final_residuals[step - 1]})
        if name in {"2.5", "10"}:
            grid = np.empty((N, N, N, 3))
            grid[tuple(ijk.T)] = u
            planes[name] = grid[:, :, N // 2, :]
    assert rows[-1]["step"] == len(times)

    k = np.fft.fftfreq(N) * N
    radius = np.sqrt(k[:, None, None]**2 + k[None, :, None]**2 + k[None, None, :]**2)
    mode_energy = .5 * np.sum(abs(np.fft.fftn(grid, axes=(0, 1, 2)) / CELLS)**2, axis=-1)
    spectrum = np.bincount(np.rint(radius).astype(int).ravel(), mode_energy.ravel())
    energy = rows[-1]["mean_kinetic_energy_m2_s2"]
    spectral_error = abs(float(spectrum.sum()) - energy) / max(energy, 1e-30)
    gates = {"enstrophy": all(r["minimum_enstrophy_s2_inverse"] >= -1e-12 for r in rows),
             "pressure": max(final_residuals) <= 1.01e-6,
             "mean_flow": max(r["mean_velocity_norm_m_s"] for r in rows) <= 1e-6,
             "forced_kinetic_energy": energy >= 1e-8,
             "parseval": spectral_error <= 1e-10}
    result = {"pass": all(gates.values()), "gates": gates, "cells": CELLS, "steps": len(times),
              "saved_states": len(rows), "parseval_relative_error": spectral_error,
              "maximum_final_pressure_residual": max(final_residuals),
              "samples": rows, "native_sha256": sources}
    out.mkdir()
    (out / "checks.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    for filename, records in (
        ("box-history.csv", rows),
        ("energy-spectrum.csv", [
            {"shell": i, "centre_wavenumber_rad_m": i * 2 * np.pi, "kinetic_energy_m2_s2": float(e)}
            for i, e in enumerate(spectrum)])):
        with (out / filename).open("x", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    fig, axs = plt.subplots(2, 2, figsize=(9.8, 7), layout="constrained")
    vmax = max(np.linalg.norm(plane, axis=-1).max() for plane in planes.values())
    x = (np.arange(N) + .5) / N
    for ax, (t, plane) in zip(axs[0], planes.items()):
        im = ax.imshow(np.linalg.norm(plane, axis=-1).T, origin="lower", extent=(0, 1, 0, 1),
                       cmap="viridis", vmin=0, vmax=vmax, interpolation="nearest")
        ax.quiver(x[::2], x[::2], plane[::2, ::2, 0].T, plane[::2, ::2, 1].T, color="white", alpha=.85,
                  angles="xy", scale_units="xy", scale=4)
        ax.set(title=f"{t} s at z = {x[N//2]:.4f} m", xlabel="x (m)", ylabel="y (m)")
    fig.colorbar(im, ax=axs[0], label="Speed (m/s)", shrink=.85)
    ax = axs[1, 0]
    other = ax.twinx()
    kplot = ax.plot(np.r_[0, TIMES], [0] + [r["mean_kinetic_energy_m2_s2"] for r in rows],
                    color="#245d83", label="Kinetic energy")
    eplot = other.plot(np.r_[0, TIMES], [0] + [r["mean_enstrophy_s2_inverse"] for r in rows],
                       color="#a95c3a", label="Enstrophy")
    ax.set(title="Volume means", xlabel="Time (s)", ylabel="Kinetic energy (m²/s²)")
    other.set_ylabel("Enstrophy (1/s²)")
    ax.legend(kplot + eplot, ["Kinetic energy", "Enstrophy"], loc="upper left", frameon=False)
    ax = axs[1, 1]
    shells = np.arange(len(spectrum))
    nonzero = (shells > 0) & (spectrum > 0)
    ax.semilogy(2*np.pi*shells[nonzero], spectrum[nonzero], "o-", ms=3, color="#245d83")
    ax.set(title="Three-dimensional spectrum at 10 s", xlabel="Shell centre (rad/m)",
           ylabel="Energy per shell (m²/s²)")
    for ax in axs[1]:
        ax.grid(alpha=.2)
    fig.suptitle("Stochastically forced flow in a periodic cube", fontweight="bold")
    for ext in ("png", "pdf"):
        fig.savefig(out / f"forced-periodic-box.{ext}", dpi=180,
                    metadata={"CreationDate": None, "ModDate": None} if ext == "pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "reader", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    raise SystemExit(0 if audit(args.case, args.reader, args.output) else 2)
