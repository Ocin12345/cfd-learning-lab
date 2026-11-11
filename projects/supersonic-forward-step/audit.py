"""Audit the native Mach-3 step flow and plot the finite-time shock structure."""
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

NX, NY, STEP_X, STEP_Y = 240, 80, 48, 16
DX = .0125
CELLS = NX*NY-(NX-STEP_X)*STEP_Y
TIMES = np.arange(1, 21)/10


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def longest_true(row):
    edges = np.flatnonzero(np.diff(np.r_[False, row, False]))
    return int(np.max(edges[1::2]-edges[::2], initial=0))


def gradient(values, axis):
    # Use central differences where both adjacent cells are fluid and
    # one-sided differences at the tunnel and step boundaries.
    padded = np.pad(values, 1, constant_values=np.nan)
    before, after = ((padded[:-2,1:-1], padded[2:,1:-1]) if axis == 0 else
                     (padded[1:-1,:-2], padded[1:-1,2:]))
    return np.where(np.isfinite(before) & np.isfinite(after), (after-before)/(2*DX),
                    np.where(np.isfinite(before),(values-before)/DX,(after-values)/DX))


def audit(case, reader_path, out):
    read = runpy.run_path(str(reader_path))["foam_values"]
    sources = []

    def native(name):
        sources.append(name)
        return (case/name).read_text("ascii")

    for line in native("input-sha256.txt").splitlines():
        expected, name = line.split("  ", 1)
        assert digest(case/name) == expected, name
    mesh, log = native("check-mesh.log"), native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*16128\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd")
    steps = np.array(re.findall(r"^Time = ([-+\d.eE]+)s$", log, re.M), float)
    assert len(steps) and np.all(np.diff(steps)>0) and abs(steps[-1]-2)<1e-8
    courant = np.array(re.findall(r"^Courant Number mean:\s+\S+\s+max:\s+([-+\d.eE]+)",log,re.M),float)
    assert len(courant) == len(steps)+1 and np.all(np.isfinite(courant)) and np.all(courant>=0)

    def field(name, kind="scalar"):
        data = read(native(name), "internalField", CELLS, kind)
        assert np.all(np.isfinite(data)), name
        return data

    c, volumes = field("0/C","vector"), field("0/Vc")
    assert np.allclose(volumes, DX**2*.1, rtol=1e-8, atol=1e-12)
    assert np.allclose(c[:,2],0,rtol=0,atol=1e-12)
    ij = np.rint(c[:,:2]/DX-.5).astype(int)
    assert np.all((ij>=0)&(ij<[NX,NY]))
    assert np.allclose(c[:,:2],(ij+.5)*DX,rtol=0,atol=1e-8)
    ix, iy = ij.T
    expected = np.ones((NY,NX), dtype=bool)
    expected[:STEP_Y,STEP_X:] = False
    observed = np.zeros_like(expected)
    observed[iy,ix] = True
    assert np.array_equal(expected,observed)

    def grid(values):
        result = np.full((NY,NX), np.nan)
        result[iy,ix] = values
        return result

    thermo = native("constant/physicalProperties")
    assert re.search(r"equationOfState\s+perfectGas;",thermo)
    gas_r = 8314.46261815324/float(re.search(r"molWeight\s+([^;]+);",thermo)[1])
    cp = float(re.search(r"\bCp\s+([^;]+);",thermo)[1])
    gamma = cp/(cp-gas_r)
    initial_u = field("0/U","vector")
    initial_p = field("0/p")
    initial_t = field("0/T")
    assert np.all(initial_p==1) and np.all(initial_t==1) and np.all(initial_u==[3,0,0])
    mach = float(np.linalg.norm(initial_u[0])/np.sqrt(gamma*gas_r*initial_t[0]))
    assert abs(mach-3)<1e-3
    x, y = c[:,:2].T
    region = (x>=.6)&(x<=2)&(y>=.25)&(y<=.85)
    saved, pressure_frames = [], {}
    names = [f"{time:g}" for time in TIMES]
    actual = {p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"[1-9]\d*(\.\d+)?|0\.\d+",p.name)}
    assert actual == set(names)
    for time, name in zip(TIMES,names):
        clock = native(f"{name}/uniform/time")
        value = float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])
        step = int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert abs(value-time)<1e-8 and abs(steps[step-1]-time)<1e-8
        p, temp, u = field(f"{name}/p"),field(f"{name}/T"),field(f"{name}/U","vector")
        assert p.shape == temp.shape == (CELLS,) and u.shape == (CELLS,3)
        assert np.all(p>0) and np.all(temp>0)
        rho = p/(gas_r*temp)
        native_error = None
        if (case/f"{name}/rho").is_file():
            stored = field(f"{name}/rho")
            assert np.all(stored>0)
            native_error = float(np.max(np.abs(stored/rho-1)))
            assert native_error < 1e-6
        pressure = grid(p)
        roof = max(longest_true(row>1.25) for row in pressure[-2:])
        pressure_frames[name] = pressure
        saved.append({"time_s":float(time),"step":step,"p_min_pa":float(p.min()),"p_max_pa":float(p.max()),
                      "T_min_k":float(temp.min()),"density_min_kg_m3":float(rho.min()),
                      "rho_crosscheck_max_relative":native_error,
                      "compressed_region_cells":int(np.count_nonzero(region & (p>1.5))),
                      "longest_compressed_roof_run_cells":roof})
    assert saved[-1]["step"]==len(steps)
    first = next((row["time_s"] for row in saved if row["longest_compressed_roof_run_cells"]>=5),None)
    passed = first is not None and saved[-1]["compressed_region_cells"]>=10
    result = {"pass":passed,"cells":CELLS,"steps":len(steps),"saved_states":len(saved),
              "nominal_mach":mach,"gamma":gamma,"initial_courant_report":float(courant[0]),
              "max_logged_courant":float(courant.max()),"courant_reports":len(courant),
              "first_saved_roof_compression_s":first,"samples":saved,
              "native_sha256":{name:digest(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"shock-step-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(saved[0]))
        writer.writeheader()
        writer.writerows(saved)

    density = grid(rho)
    grad = np.hypot(gradient(density,0),gradient(density,1))
    vmax = max(np.nanmax(pressure_frames["0.5"]),np.nanmax(pressure_frames["2"]))
    fig, axes = plt.subplots(3,1,figsize=(9,8.3),layout="constrained")
    panels = ((pressure_frames["0.5"],"Pressure at 0.5 s","p / inlet pressure","magma",vmax),
              (pressure_frames["2"],"Pressure at 2 s","p / inlet pressure","magma",vmax),
              (grad,"Density gradient at 2 s","Magnitude (kg/m⁴)","inferno",np.nanpercentile(grad,99)))
    for ax,(data,title,label,cmap,upper) in zip(axes,panels):
        picture=ax.pcolormesh(np.arange(NX+1)*DX,np.arange(NY+1)*DX,data,
                              cmap=cmap,vmin=0,vmax=upper,shading="flat")
        ax.set(xlim=(0,3),ylim=(0,1),aspect="equal",ylabel="Height (m)",title=title)
        fig.colorbar(picture,ax=ax,label=label,shrink=.87,extend="max" if data is grad else "neither")
    axes[-1].set_xlabel("Distance (m)")
    fig.suptitle("Mach 3 flow over a forward-facing step",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"supersonic-forward-step.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return passed


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case","reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)