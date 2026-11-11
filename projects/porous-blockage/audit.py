"""Check the saved porous-blockage calculation and plot measured flow quantities."""
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
from matplotlib.patches import Rectangle
import numpy as np

CELLS = 2048
TIMES = tuple(range(1, 101))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, reader_path, output):
    values = runpy.run_path(str(reader_path))["foam_values"]

    def field(path, kind="scalar"):
        return values((case/path).read_text(encoding="ascii"), "internalField", CELLS, kind)

    def patch_flux(text, name):
        patch = re.search(r"\b"+name+r"\s*\{([^{}]*)\}", text, re.S)
        assert patch, name
        return float(values(patch[1], "value", 32, "scalar").sum())

    for line in (case/"input-sha256.txt").read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ", 1)
        assert sha(case/name) == expected, name
    mesh = (case/"check-mesh.log").read_text(encoding="ascii")
    log = (case/"solver.log").read_text(encoding="ascii")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*2048\s*$", mesh, re.M)
    assert re.search(r"^End\s*$", log, re.M)
    steps = np.array(re.findall(r"^Time = ([-+0-9.eE]+)s$", log, re.M), dtype=float)
    assert len(steps) == 2000 and np.allclose(steps, .05*np.arange(1, 2001), rtol=0, atol=1e-8)

    x, y, z = field("0/C", "vector").T
    assert np.allclose(z, 0, rtol=0, atol=1e-12)
    vol = field("0/Vc")
    assert np.allclose(vol, .125*.125*.2, rtol=0, atol=1e-10)
    ix = np.rint((x+2)/.125-.5).astype(int)
    iy = np.rint((y+2)/.125-.5).astype(int)
    order = np.argsort(iy*64+ix)
    assert np.array_equal((iy*64+ix)[order], np.arange(CELLS))
    assert np.allclose(x, -2+(ix+.5)*.125, rtol=0, atol=1e-8)
    assert np.allclose(y, -2+(iy+.5)*.125, rtol=0, atol=1e-8)
    zone = (np.abs(x)<.5)&(np.abs(y)<.5)
    bypass = (np.abs(x)<.5)&(np.abs(y)>1)
    up = np.isclose(x,-.5625)&(np.abs(y)<.5)
    down = np.isclose(x,.5625)&(np.abs(y)<.5)
    raw = (case/"constant/polyMesh/cellZones").read_text(encoding="ascii")
    match = re.search(r"\bporousBlockage\s*\{\s*cellLabels\s+List<label>\s+64\s*\(([^)]*)\)",raw,re.S)
    assert match
    assert np.array_equal(np.sort(np.fromstring(match[1],dtype=int,sep=" ")),np.flatnonzero(zone))
    assert zone.sum()==64 and up.sum()==down.sum()==8

    rows = []
    sources = ["input-sha256.txt","check-mesh.log","solver.log","0/C","0/Vc","constant/polyMesh/cellZones"]
    for time in TIMES:
        directory = case/str(time)
        clock = (directory/"uniform/time").read_text(encoding="ascii")
        index = int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        value = float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])
        assert index == 20*time and abs(value-time)<1e-8 and abs(steps[index-1]-time)<1e-8
        u = field(f"{time}/U","vector")
        p = field(f"{time}/p")
        phi = (directory/"phi").read_text(encoding="ascii")
        assert re.search(r"\bdimensions\s+\[0\s+3\s+-1\s+0\s+0\s+0\s+0\]",phi)
        inlet, outlet = (patch_flux(phi,name) for name in ("inlet","outlet"))
        mismatch = abs(inlet+outlet)/abs(inlet)
        assert inlet < 0 < outlet and mismatch < 1e-4, (time,mismatch)
        rows.append({"time_s":time,"native_step":index,"zone_mean_Ux_m_s":float(u[zone,0].mean()),
                     "bypass_mean_Ux_m_s":float(u[bypass,0].mean()),"centre_band_pressure_drop_m2_s2":float(p[up].mean()-p[down].mean()),
                     "inlet_m3_s":inlet,"outlet_m3_s":outlet,"relative_flow_mismatch":mismatch})
        sources.extend(f"{time}/{name}" for name in ("uniform/time","U","p","phi"))
    late = rows[80:]
    means = {key:float(np.mean([row[key] for row in late])) for key in
             ("zone_mean_Ux_m_s","bypass_mean_Ux_m_s","centre_band_pressure_drop_m2_s2")}
    ratio = means["zone_mean_Ux_m_s"]/means["bypass_mean_Ux_m_s"]
    assert 0 < ratio < .8 and means["centre_band_pressure_drop_m2_s2"] > 0
    checks = {"cells":CELLS,"zone_cells":int(zone.sum()),"native_steps":len(steps),"saved_states":len(rows),
              "late_sample_means":means,"late_zone_to_bypass_speed_ratio":ratio,
              "max_relative_flow_mismatch":max(row["relative_flow_mismatch"] for row in rows),
              "samples":rows,"native_sha256":{name:sha(case/name) for name in sources}}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    with (output/"porous-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes=plt.subplots(2,2,figsize=(11,6.5),layout="constrained")
    for ax, field_, title, label in (
        (axes[0,0],np.linalg.norm(u,axis=1),"Velocity magnitude at 100 s","Speed (m/s)"),
        (axes[0,1],p,"Kinematic pressure at 100 s","p (m²/s²)")):
        panel=ax.imshow(field_[order].reshape(32,64),origin="lower",extent=(-2,6,-2,2),
                        interpolation="nearest",cmap="viridis",aspect="equal")
        ax.add_patch(Rectangle((-.5,-.5),1,1,fill=False,color="white",linestyle="--"))
        ax.set(xlabel="x (m)",ylabel="y (m)",title=title)
        fig.colorbar(panel,ax=ax,label=label,shrink=.9)
    for key,label in (("zone_mean_Ux_m_s","Porous zone"),("bypass_mean_Ux_m_s","Bypass")):
        axes[1,0].plot(TIMES,[row[key] for row in rows],label=label)
    axes[1,0].set(xlabel="Time (s)",ylabel="Mean streamwise speed (m/s)",title="Flow through and around the blockage")
    axes[1,0].legend()
    axes[1,1].plot(TIMES,[row["centre_band_pressure_drop_m2_s2"] for row in rows])
    axes[1,1].set(xlabel="Time (s)",ylabel="Pressure difference (m²/s²)",title="Centre-band pressure drop")
    for ax in axes[1]:
        ax.axvspan(81,100,color="#8eb1cd",alpha=.18)
        ax.grid(alpha=.2)
    fig.suptitle("Flow through a Darcy porous blockage | OpenFOAM",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(output/f"porous-blockage.{ext}",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","foam_reader","output"):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    audit(args.case,args.foam_reader,args.output)
