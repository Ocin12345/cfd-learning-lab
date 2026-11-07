"""Compute 1.2 s of small-amplitude free-surface sloshing in hidden Fluent."""
import csv
import json
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import h5py
import numpy as np
from analyse_nozzle_fluent_v005 import common_mapping

NX, NY = 160, 80
DX, DT = 0.00125, 0.0005
L, H = NX * DX, NY * DX
WATER_WIDTH, WATER_HEIGHT = 0.20, 0.05
MATERIALS = {"air": (1.225, 1.7894e-5), "water-liquid": (998.2, 0.001003)}
WALLS = ["tank-walls"]
MODE = os.environ.get("CFD_FINE_MODE", "tilted")
assert MODE in ("flat", "tilted", "production")
LAST_STEP = 2400 if MODE == "production" else 40
FRAME_INTERVAL = 20 if MODE == "production" else 1
EXPECTED_AREA = WATER_WIDTH * WATER_HEIGHT
RESIDUAL_TARGET = 1e-4
INITIAL_TILT = 0.0 if MODE == "flat" else 0.025
RAW_TARGETS = {"continuity":1e-12,"x-velocity":1e-10,"y-velocity":1e-10}


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def make_mesh(out):
    x = DX * np.arange(NX + 1)
    y = DX * np.arange(NY + 1)
    nodes = np.array([(xx, yy) for yy in y for xx in x])
    quads = []
    faces = {}
    for j in range(NY):
        for i in range(NX):
            k = j * (NX + 1) + i
            q = (k, k + 1, k + NX + 2, k + NX + 1)
            quads.append(q)
            for a, b in zip(q, (*q[1:], q[0])):
                key = tuple(sorted((a, b)))
                if key in faces:
                    faces[key][3] = len(quads)
                else:
                    faces[key] = [a + 1, b + 1, len(quads), 0]
    inside = [row for row in faces.values() if row[3]]
    wall = [row for row in faces.values() if not row[3]]
    assert len(wall) == 2 * (NX + NY)
    assert len(inside) == (NX - 1) * NY + (NY - 1) * NX
    centres = nodes[np.asarray(quads)].mean(axis=1)
    np.savez_compressed(out / "geometry.npz", nodes=nodes, quads=quads,
                        centres=centres, x_edges=x, y_edges=y)
    mesh = out / "tank.msh"
    with mesh.open("x", encoding="ascii", newline="\n") as f:
        f.write('(0 "Sloshing tank, uniform 1.25 mm quads")\n(2 2)\n')
        f.write(f"(10 (0 1 {len(nodes):x} 1 2))\n(10 (1 1 {len(nodes):x} 1 2)(\n")
        for xx, yy in nodes:
            f.write(f"{xx:.17g} {yy:.17g}\n")
        f.write(f"))\n(12 (0 1 {len(quads):x} 0))\n(12 (2 1 {len(quads):x} 1 3))\n")
        f.write(f"(13 (0 1 {len(faces):x} 0))\n")
        first = 1
        for zone, bc, rows in ((3, 2, inside), (4, 3, wall)):
            f.write(f"(13 ({zone:x} {first:x} {first + len(rows) - 1:x} {bc:x} 2)(\n")
            for row in rows:
                f.write(" ".join(f"{value:x}" for value in row) + "\n")
            f.write("))\n")
            first += len(rows)
        f.write("(45 (2 fluid fluid)())\n(45 (3 interior interior)())\n(45 (4 wall tank-walls)())\n")
    return mesh


def alpha_group(handle):
    phases = handle["results/1"]
    groups = [phase["cells/SV_VOF"] for phase in phases.values() if "cells/SV_VOF" in phase]
    if len(groups) != 1:
        matches = [group for group in groups if "/phase-3/" in group.name]
        if len(matches) != 1:
            raise ValueError(f"Expected Fluent secondary-phase VOF: {[g.name for g in groups]}")
        return matches[0]
    return groups[0]


def cell_values(group):
    result = np.full(NX * NY, np.nan)
    for ds in group.values():
        first = int(np.asarray(ds.attrs["minId"]).flat[0]) - 1
        last = int(np.asarray(ds.attrs["maxId"]).flat[0])
        assert 0 <= first < last <= len(result) and np.isnan(result[first:last]).all()
        result[first:last] = np.asarray(ds).ravel()
    assert np.isfinite(result).all()
    return result


