"""Audit early compressible gas-pocket expansion from native OpenFOAM output."""
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

NX, NY = 80, 160
CELLS = NX*NY
DX = DY = .0125
DEPTH = .2
TIMES = np.arange(1, 11)/1000
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def audit(case, reader, out):
    read = runpy.run_path(str(reader))["foam_values"]
    sources = {}

    def native(name):
        raw = (case/name).read_bytes()
        sources[name] = hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest, name = line.split("  ", 1)
        assert hashlib.sha256((case/name).read_bytes()).hexdigest() == digest, name
        sources[name] = digest
    for name in ("pre-mesh-input-sha256.txt", "pre-mesh-inputs-verified.log", "inputs-verified.log",
                 "mesh.log", "centres.log", "volumes.log"):
        native(name)
    mesh, log = native("check-mesh.log"), native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*12800\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd") and "compressibleVoF" in log
    times = np.array(re.findall(r"^Time = ("+NUMBER+r")s?$",log,re.M),float)
    assert len(times) and np.isfinite(times).all() and (np.diff(times)>0).all() and times[-1]==.01

    def field(name, kind="scalar"):
        data = read(native(name), "internalField", CELLS, kind)
        assert data.shape == ((CELLS,3) if kind=="vector" else (CELLS,)), name
        assert np.isfinite(data).all(), name
        return data

    c, v = field("0/C","vector"), field("0/Vc")
    assert (v>0).all()
    np.testing.assert_allclose(v, DX*DY*DEPTH, atol=1e-14, rtol=0)
    np.testing.assert_allclose(c[:,2], 0, atol=1e-10, rtol=0)
    ix, iy = np.floor(c[:,:2]/[DX,DY]).astype(int).T
    assert np.all((ix>=0)&(ix<NX)&(iy>=0)&(iy<NY)) and len(np.unique(iy*NX+ix))==CELLS
    np.testing.assert_allclose(c[:,:2], np.column_stack(((ix+.5)*DX,(iy+.5)*DY)), atol=1e-9,rtol=0)
    assert np.all(field("0/U","vector")==0)
    init = native("0/alpha.water")
    assert re.search(r"radius\s+0\.1;",init) and re.search(r"point1\s+\(0\.5 0\.5 -1\);",init)
    assert re.search(r"point\s+\(0 1 0\);",init)
    assert re.search(r"value\s+1e6;",native("0/p"))
    assert re.search(r"value\s+578;",native("0/T"))
    lower = c[:,1] < .85
    original = ((c[:,0]-.5)**2+(c[:,1]-.5)**2) <= .1**2
    ideal_volume = float(np.pi*.1**2*DEPTH)
    actual = {p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?",p.name)}
    assert actual == {"0",*(f"{t:g}" for t in TIMES)}

    def grid(values):
        out = np.empty((NY,NX))
        out[iy,ix] = values
        return out

    rows = []
    maps = {}
    for t in TIMES:
        name = f"{t:g}"
        clock = native(f"{name}/uniform/time")
        step = int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert 1<=step<=len(times) and abs(times[step-1]-t)<1e-10
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-t)<1e-10
        a, p, temp = field(f"{name}/alpha.water"), field(f"{name}/p"), field(f"{name}/T")
        u = field(f"{name}/U","vector")
        gas = 1-a
        weights = gas[lower]*v[lower]
        vg = float(weights.sum())
        assert vg>0
        rows.append({"time_s":float(t), "step":step, "alpha_water_min":float(a.min()),
                     "alpha_water_max":float(a.max()), "minimum_absolute_p_pa":float(p.min()),
                     "maximum_absolute_p_pa":float(p.max()), "minimum_temperature_k":float(temp.min()),
                     "maximum_temperature_k":float(temp.max()),
                     "gas_volume_below_y085_m3":vg,
                     "gas_volume_ratio_to_ideal_initial_cylinder":vg/ideal_volume,
                     "gas_weighted_absolute_p_pa":float(weights@p[lower]/vg),
                     "gas_volume_in_original_radius_m3":float(gas[original]@v[original]),
                     "maximum_speed_m_s":float(np.linalg.norm(u,axis=1).max())})
        if name in {"0.001","0.005","0.01"}:
            maps[name]=(grid(gas),grid(p))
    assert rows[-1]["step"]==len(times)
    gates = {"phase_fraction":all(-1e-6<=r["alpha_water_min"]<=r["alpha_water_max"]<=1+1e-6 for r in rows),
             "thermodynamic_fields":all(r["minimum_absolute_p_pa"]>0 and r["minimum_temperature_k"]>0 for r in rows),
             "gas_in_original_region":all(r["gas_volume_in_original_radius_m3"]>0 for r in rows),
             "gas_expansion":max(r["gas_volume_ratio_to_ideal_initial_cylinder"] for r in rows)>1.10,
             "gas_pressure_drop":rows[-1]["gas_weighted_absolute_p_pa"]<800000}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"steps":len(times),
            "saved_states":len(rows),"ideal_initial_cylinder_volume_m3":ideal_volume,
            "samples":rows,"native_sha256":sources}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"gas-pocket-history.csv").open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axs=plt.subplots(2,3,figsize=(11.5,7),layout="constrained")
    for ax,(name,(gas,p)) in zip(axs[0],maps.items()):
        im=ax.imshow(gas,origin="lower",extent=(0,1,0,2),vmin=0,vmax=1,
                     cmap="viridis",interpolation="nearest")
        ax.axhline(.85,color="#fafafa",ls="--",lw=1)
        ax.set(title=f"Gas at {1000*float(name):g} ms",xlabel="x (m)",ylabel="y (m)",
               xlim=(.2,.8),ylim=(.2,1.1),aspect="equal")
    fig.colorbar(im,ax=axs[0],label="Gas volume fraction",shrink=.84)
    ax=axs[1,0]
    ax.plot(1000*TIMES,[r["gas_volume_below_y085_m3"]*1000 for r in rows],"o-",ms=3,color="#286c68")
    ax.scatter(0,ideal_volume*1000,facecolors="none",edgecolors="#a75436",label="Initial cylinder geometry")
    ax.set(title="Gas below y = 0.85 m",xlabel="Time (ms)",ylabel="Gas volume (L)")
    ax.legend(frameon=False,fontsize=8)
    ax=axs[1,1]
    ax.plot(1000*TIMES,[r["gas_weighted_absolute_p_pa"]/1000 for r in rows],"o-",ms=3,color="#2a5a84")
    ax.scatter(0,1000,facecolors="none",edgecolors="#a75436",label="Prescribed pocket pressure")
    ax.set(title="Gas-weighted pressure",xlabel="Time (ms)",ylabel="Absolute pressure (kPa)")
    ax.legend(frameon=False,fontsize=8)
    ax=axs[1,2]
    for name,(gas,p) in maps.items():
        ax.plot((np.arange(NX)+.5)*DX,p[NY//4-1:NY//4+1].mean(axis=0)/1000,
                label=f"{1000*float(name):g} ms")
    ax.set(title="Pressure across the pocket",xlabel="x (m)",ylabel="Absolute pressure (kPa)")
    ax.legend(frameon=False)
    for ax in axs[1]:
        ax.grid(alpha=.2)
    fig.suptitle("Early expansion of a submerged gas pocket",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"submerged-gas-expansion.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)
