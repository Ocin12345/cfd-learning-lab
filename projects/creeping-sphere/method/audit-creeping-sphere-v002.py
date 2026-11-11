"""Check the native sphere mesh, steady flow, forces and Stokes comparison."""
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
from matplotlib.patches import Circle
import numpy as np

CELLS = 36_864
RADIUS = 1.0
U_INF = 0.05
RHO = NU = 1.0
STOKES = 6 * np.pi * RHO * NU * RADIUS * U_INF
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
EQUATIONS = ("Ux", "Uy", "Uz", "p")
PATCHES = {"inlet": 256, "outlet": 256, "farfield": 1024, "sphere": 1536}


def audit(case, reader, out):
    read = runpy.run_path(str(reader))["foam_values"]
    sources = {}

    def native(name):
        raw = (case / name).read_bytes()
        sources[name] = hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for manifest in ("pre-mesh-input-sha256.txt", "input-sha256.txt"):
        for line in native(manifest).splitlines():
            digest, name = line.split("  ", 1)
            native(name)
            assert sources[name] == digest, name
    for name in ("pre-mesh-inputs-verified.log", "inputs-verified.log", "mesh.log",
                 "centres.log", "volumes.log", "constant/polyMesh/owner",
                 "constant/polyMesh/neighbour"):
        native(name)
    mesh = native("check-mesh.log")
    assert "Mesh OK." in mesh
    assert re.search(r"^\s*cells:\s*36864\s*$", mesh, re.M)
    assert re.search(r"^\s*hexahedra:\s*36864\s*$", mesh, re.M)
    log = native("solver.log")
    assert log.rstrip().endswith("\nEnd")
    assert "OpenFOAM-14" in log or "14-" in log
    steps = list(re.finditer(r"^Time = (\d+)s?$", log, re.M))
    end = len(steps)
    assert 1 <= end <= 1000 and [int(s[1]) for s in steps] == list(range(1, end+1))
    rows = []
    for i, step in enumerate(steps):
        stop = steps[i+1].start() if i+1 < end else len(log)
        block = log[step.end():stop]
        entries = re.findall(r"Solving for (Ux|Uy|Uz|p), Initial residual = ("+NUMBER+r")", block)
        assert tuple(key for key, _ in entries) == EQUATIONS
        continuity = re.findall(r"time step continuity errors : sum local = ("+NUMBER+r")", block)
        assert len(continuity) == 1
        rows.append([float(value) for _, value in entries] + [float(continuity[0])])
    residuals = np.array(rows)
    assert np.isfinite(residuals).all() and (residuals >= 0).all()

    boundary = native("constant/polyMesh/boundary")
    start_face = {}
    for patch, count in PATCHES.items():
        body = re.search(r"\b"+patch+r"\s*\{([^}]+)\}", boundary, re.S)[1]
        assert int(re.search(r"nFaces\s+(\d+);", body)[1]) == count
        start_face[patch] = int(re.search(r"startFace\s+(\d+);", body)[1])
    points_text = native("constant/polyMesh/points")
    faces_text = native("constant/polyMesh/faces")
    for line in native("projection-sha256.txt").splitlines():
        digest, name = line.split("  ", 1)
        native(name)
        assert sources[name] == digest, name
    projection = json.loads(native("sphere-projection.json"))
    assert projection["original_points_sha256"] == sources["preprojection-points"]
    assert projection["revised_points_sha256"] == sources["constant/polyMesh/points"]
    assert projection["sphere_vertices"] == 1538 and projection["other_vertices_unchanged"] == 36912
    assert projection["maximum_radius_error_after_m"] <= 1e-14
    assert all(sources["constant/polyMesh/"+name] == digest
               for name, digest in projection["unchanged_topology_sha256"].items())
    points = np.array(re.findall(r"^\(("+NUMBER+r") ("+NUMBER+r") ("+NUMBER+r")\)$",
                                points_text, re.M), dtype=float)
    faces = np.array(re.findall(r"^4\((\d+) (\d+) (\d+) (\d+)\)$", faces_text, re.M), dtype=int)
    assert points.shape == (38450, 3) and faces.shape == (112128, 4)
    assert np.isfinite(points).all() and faces.min() == 0 and faces.max() < len(points)
    sphere_vertices = np.unique(faces[start_face["sphere"]:start_face["sphere"]+PATCHES["sphere"]])
    assert len(sphere_vertices) == 1538
    radius_error = float(np.max(abs(np.linalg.norm(points[sphere_vertices], axis=1)-RADIUS)))

    def field(name, kind="scalar"):
        values = read(native(name), "internalField", CELLS, kind)
        assert values.shape == ((CELLS, 3) if kind == "vector" else (CELLS,))
        assert np.isfinite(values).all(), name
        return values

    centres, volumes = field("0/C", "vector"), field("0/Vc")
    assert (volumes > 0).all()
    assert (np.linalg.norm(centres, axis=1) > RADIUS).all()
    assert abs(np.max(abs(points)) - 18.0) < 1e-12
    assert abs(volumes.sum() - (36**3-4*np.pi/3)) < .03
    saves = sorted(int(p.name) for p in case.iterdir()
                   if p.is_dir() and p.name.isdigit() and int(p.name) > 0)
    assert len(saves) >= 2 and saves == sorted(set(range(25, end+1, 25)) | {end})
    previous_u = None
    for step in saves:
        clock = native(f"{step}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);", clock, re.M)[1]) == step
        assert float(re.search(r"^value\s+([^;]+);", clock, re.M)[1]) == step
        u, p = field(f"{step}/U", "vector"), field(f"{step}/p")
        if step == saves[-2]:
            previous_u = u
        flux_text = native(f"{step}/phi")
        flow = {}
        sealed = 0.
        for patch, count in PATCHES.items():
            body = re.search(r"\b"+patch+r"\s*\{([^}]+)\}", flux_text, re.S)[1]
            values = read(body, "value", count, "scalar")
            assert values.shape == (count,) and np.isfinite(values).all()
            flow[patch] = float(values.sum())
            if patch in ("sphere", "farfield"):
                sealed = max(sealed, float(np.max(abs(values))))
    assert flow["inlet"] < 0 < flow["outlet"]
    imbalance = abs(sum(flow.values()))/abs(flow["inlet"])
    rms_change = float(np.sqrt(np.mean(np.sum((u-previous_u)**2, axis=1))))

    force_text = native("postProcessing/sphereForce/0/forces.dat")
    assert "forces(pressure viscous)" in force_text
    force_rows = [re.findall(NUMBER, line) for line in force_text.splitlines()
                  if line.strip() and not line.startswith("#")]
    assert all(len(row) == 13 for row in force_rows)
    forces = np.array(force_rows, dtype=float)
    assert np.isfinite(forces).all() and np.array_equal(forces[:, 0], np.arange(end+1))
    pressure, viscous = forces[:, 1:4], forces[:, 4:7]
    total = pressure + viscous
    drag = float(total[-1, 0])
    assert drag > 0
    transverse_ratio = float(np.linalg.norm(total[-1, 1:])/drag)
    drag_error = abs(drag-STOKES)/STOKES
    drag_change = abs(drag-float(total[saves[-2], 0]))/drag
    gates = {
        "native_convergence": f"PIMPLE solution converged in {end} iterations" in log,
        "residuals": bool(residuals[-1, :4].max() <= 1e-5),
        "sphere_radius": radius_error <= 5e-5,
        "boundary_balance": imbalance <= 1e-3,
        "impermeable_patches": sealed <= 1e-12,
        "transverse_force": transverse_ratio <= .01,
        "stokes_drag_agreement": bool(drag_error <= .25),
        "terminal_velocity_stability": rms_change <= .001,
        "terminal_drag_stability": drag_change <= .005,
    }
    result = {
        "pass": all(gates.values()), "gates": gates, "cells": CELLS,
        "sphere_faces": PATCHES["sphere"], "sphere_vertices": len(sphere_vertices),
        "radius_max_abs_error_m": radius_error, "cell_volume_range_m3": [float(volumes.min()), float(volumes.max())],
        "reynolds_number_diameter": 2*RADIUS*U_INF/NU, "iterations": end, "saved_states": len(saves),
        "terminal_initial_residuals": dict(zip(EQUATIONS, residuals[-1, :4])),
        "terminal_local_continuity_error": float(residuals[-1, -1]),
        "boundary_flow_m3_s": flow, "relative_boundary_imbalance": imbalance,
        "max_impermeable_flux_m3_s": sealed,
        "terminal_kinematic_pressure_range_m2_s2": [float(p.min()), float(p.max())],
        "stokes_reference_n": float(STOKES), "terminal_pressure_force_n": pressure[-1].tolist(),
        "terminal_viscous_force_n": viscous[-1].tolist(), "terminal_total_force_n": total[-1].tolist(),
        "relative_stokes_drag_difference": float(drag_error), "relative_transverse_force": transverse_ratio,
        "previous_save_iteration": saves[-2], "velocity_rms_change_m_s": rms_change,
        "relative_drag_change": drag_change, "surface_projection": projection, "native_sha256": sources,
    }
    out.mkdir()
    (out / "checks.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    with (out / "force-history.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("steady_iteration", "pressure_fx_n", "pressure_fy_n", "pressure_fz_n",
                         "viscous_fx_n", "viscous_fy_n", "viscous_fz_n", "total_fx_n",
                         "total_fy_n", "total_fz_n", "stokes_fx_n"))
        for i, (fp, fv, ft) in enumerate(zip(pressure, viscous, total)):
            writer.writerow((i, *fp, *fv, *ft, STOKES))
    with (out / "residual-history.csv").open("x", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("steady_iteration", *EQUATIONS, "local_continuity_error"))
        writer.writerows((i, *row) for i, row in enumerate(residuals, 1))

    speed = np.linalg.norm(u, axis=1)
    selected = (abs(centres[:, 2]) <= .18) & (abs(centres[:, 0]) <= 4) & (abs(centres[:, 1]) <= 3)
    assert selected.sum() > 100
    fig = plt.figure(figsize=(11, 6), layout="constrained")
    fig.suptitle("Slow flow around a fixed sphere", fontsize=17, fontweight="bold")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.4, 1])
    ax = fig.add_subplot(grid[:, 0])
    ax.set_facecolor("#f1f3f5")
    field_plot = ax.scatter(centres[selected, 0], centres[selected, 1], c=speed[selected]/U_INF,
                            s=9, cmap="viridis", vmin=0, vmax=max(1., float(speed[selected].max()/U_INF)),
                            edgecolors="none", rasterized=True)
    ax.add_patch(Circle((0, 0), RADIUS, fc="white", ec="#263347", lw=1.1))
    ax.set(xlabel="Streamwise position x (m)", ylabel="Transverse position y (m)",
           title="Native cell centres within 0.18 m of the equator", xlim=(-4.2, 4.2), ylim=(-3.2, 3.2),
           aspect="equal")
    fig.colorbar(field_plot, ax=ax, label="Speed / inlet speed", shrink=.7)
    ax = fig.add_subplot(grid[0, 1])
    late = np.arange(25, end+1)
    ax.plot(late, total[late, 0], color="#1d6996", label="Total CFD drag")
    ax.plot(late, pressure[late, 0], color="#e17c05", lw=1, label="Pressure")
    ax.plot(late, viscous[late, 0], color="#4b7e42", lw=1, label="Viscous")
    ax.axhline(STOKES, color="#263347", ls="--", lw=1, label="Stokes reference")
    ax.set(title=f"Streamwise force | final {drag:.3f} N", ylabel="Force (N)", xlim=(25, end))
    ax.legend(frameon=False, fontsize=8, ncols=2)
    ax.grid(alpha=.2)
    ax = fig.add_subplot(grid[1, 1])
    for j, name in enumerate(EQUATIONS):
        positive = residuals[:, j] > 0
        ax.semilogy(np.arange(1, end+1)[positive], residuals[positive, j], label=name, lw=1.15)
    ax.axhline(1e-5, color="#263347", ls="--", lw=1, label="Limit")
    ax.set(title="Initial equation residuals", xlabel="Steady iteration", ylabel="Initial residual")
    ax.legend(frameon=False, fontsize=8, ncols=3)
    ax.grid(alpha=.2)
    fig.savefig(out / "creeping-sphere.png", dpi=175)
    fig.savefig(out / "creeping-sphere.pdf")
    plt.close(fig)
    print(json.dumps({k: v for k, v in result.items() if k != "native_sha256"}, indent=2))
    return result["pass"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("case", "reader", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    raise SystemExit(0 if audit(args.case, args.reader, args.output) else 2)
