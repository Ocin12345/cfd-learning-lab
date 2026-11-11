"""Check the periodic duct's saved fields and plot the computed quarter section."""
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

NX, NY = 134, 70
CELLS = NX*NY
DX, DY, DEPTH = .1/NX, .0525/NY, .001


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stress_values(text):
    match = re.search(r"\binternalField\s+nonuniform\s+List<symmTensor>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if match:
        assert int(match[1]) == CELLS
        rows = re.findall(r"\(([^()]+)\)", match[2])
        assert len(rows) == CELLS
        values = np.array([np.fromstring(row, sep=" ") for row in rows])
    else:
        row = re.search(r"\binternalField\s+uniform\s*\(([^()]+)\)\s*;", text)[1]
        values = np.tile(np.fromstring(row, sep=" "), (CELLS, 1))
    assert values.shape == (CELLS, 6) and np.all(np.isfinite(values))
    return values


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
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*9380\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd") and re.search(r"Selecting RAS turbulence model LRR\b", log)
    iterations = np.array(re.findall(r"^Time = (\d+)s$", log, re.M), int)
    assert 0 < len(iterations) <= 6000
    endpoint = len(iterations)
    assert np.array_equal(iterations, np.arange(1, endpoint+1))
    assert re.search(r"solution converged", log, re.I) or endpoint == 6000
    converged = bool(re.search(r"solution converged", log, re.I))
    last = log.rsplit(f"Time = {endpoint}s",1)[1]
    local = [float(x) for x in re.findall(r"sum local = ([-+\d.eE]+)", last)]
    assert local and np.isfinite(local[-1])

    def field(name, kind="scalar"):
        values = read(native(name), "internalField", CELLS, kind)
        assert np.all(np.isfinite(values)), name
        return values

    centres, volumes = field("0/C", "vector"), field("0/Vc")
    assert np.allclose(volumes, DX*DY*DEPTH, rtol=1e-8, atol=1e-14)
    assert np.allclose(centres[:,2], DEPTH/2, rtol=0, atol=1e-12)
    ij = np.rint(centres[:,:2]/[DX,DY]-.5).astype(int)
    assert np.all((ij>=0) & (ij<[NX,NY]))
    assert np.allclose(centres[:,:2], (ij+.5)*[DX,DY], rtol=0, atol=1e-8)
    indices = ij[:,1]*NX+ij[:,0]
    assert np.array_equal(np.sort(indices), np.arange(CELLS))
    order = np.argsort(indices)

    saves = sorted(int(path.name) for path in case.iterdir()
                   if path.is_dir() and re.fullmatch(r"[1-9][0-9]*",path.name))
    assert saves == sorted(set(range(250, endpoint+1,250))|{endpoint})
    rows = []
    previous = None
    for index in saves:
        clock = native(f"{index}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);", clock, re.M)[1]) == index
        assert float(re.search(r"^value\s+([^;]+);", clock, re.M)[1]) == index
        u = field(f"{index}/U", "vector")
        field(f"{index}/p")
        epsilon = field(f"{index}/epsilon")
        stresses = stress_values(native(f"{index}/R"))
        normal_min = float(stresses[:,[0,3,5]].min())
        secondary = np.linalg.norm(u[:,:2],axis=1)
        rms_change = None if previous is None else float(np.sqrt(np.average(
            np.sum((u-previous)**2, axis=1),weights=volumes)))
        rows.append({"iteration":index,
                     "axial_mean_m_s":float(np.average(u[:,2],weights=volumes)),
                     "secondary_rms_m_s":float(np.sqrt(np.average(secondary**2,weights=volumes))),
                     "secondary_max_m_s":float(secondary.max()),
                     "velocity_rms_change_m_s":rms_change,
                     "min_normal_stress_m2_s2":normal_min,
                     "min_epsilon_m2_s3":float(epsilon.min())})
        previous = u
    axial, secondary_rms = rows[-1]["axial_mean_m_s"],rows[-1]["secondary_rms_m_s"]
    gates = {"stress_and_dissipation":all(r["min_normal_stress_m2_s2"] >= -1e-10
                                         and r["min_epsilon_m2_s3"] > 0 for r in rows),
             "local_continuity":local[-1] < 1e-4, "axial_speed":abs(axial/1.95-1) <= .005,
             "secondary_speed":.001 < secondary_rms < .15*axial,
             "convergence":converged or (len(rows)>=2 and rms_change < .01*axial)}
    summary = {"pass":all(gates.values()), "gates":gates, "cells":CELLS,
               "iterations":endpoint, "native_convergence":converged,
               "terminal_final_local_continuity":local[-1], "samples":rows,
               "native_sha256":{name:sha(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(summary, indent=2)+"\n",encoding="utf-8")
    with (out/"duct-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1,2,figsize=(9.6,4.15),layout="constrained")
    xs, ys = np.arange(NX+1)*DX*1000, np.arange(NY+1)*DY*1000
    fields = [(u[:,2],"Axial speed (m/s)","Axial velocity"),
              (secondary,"Secondary speed (m/s)","Cross-section circulation")]
    for ax,(values,label,title) in zip(axes,fields):
        im = ax.pcolormesh(xs, ys, values[order].reshape(NY,NX), cmap="viridis",shading="flat")
        ax.set(xlim=(0,100),ylim=(0,52.5),aspect="equal",xlabel="x (mm)",ylabel="y (mm)",title=title)
        fig.colorbar(im,ax=ax,label=label,shrink=.74)
    selected = (ij[:,0]%8==4) & (ij[:,1]%6==3)
    arrows = axes[1].quiver(centres[selected,0]*1000,centres[selected,1]*1000,
                            u[selected,0],u[selected,1],color="white",angles="xy",
                            scale_units="xy",scale=.012,width=.003)
    axes[1].quiverkey(arrows,.72,.91,.02,"0.02 m/s",labelpos="E",color="white",labelcolor="white")
    fig.suptitle("Reynolds-stress-driven secondary flow in a rectangular duct",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"duct-secondary-flow.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return summary["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case","reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)