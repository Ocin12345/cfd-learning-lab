"""Measure the saved interface and equivalent water height in the capillary."""
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

WIDTH = .001
HEIGHT = .020
NX, NY = 20, 400
CELLS = NX*NY
TIMES = [i/20 for i in range(1,11)]
ALPHA_TOL = 1e-6
NUMBER = r"[-+0-9.eE]+"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, foam_reader, output):
    reader = runpy.run_path(str(foam_reader))["foam_values"]
    def field(file, kind="scalar"):
        return reader(file.read_text(encoding="ascii"),"internalField",CELLS,kind)

    log = (case/"solver.log").read_text(encoding="ascii")
    mesh = (case/"check-mesh.log").read_text(encoding="ascii")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*8000\s*$",mesh,re.M)
    assert re.search(r"^End\s*$",log,re.M)
    for line in (case/"input-sha256.txt").read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ",1)
        assert sha(case/name)==expected,name
    logged_times = np.array(re.findall(r"^Time = ("+NUMBER+r")s$",log,re.M),dtype=float)
    bulk = np.array(re.findall(r"^Courant Number mean: "+NUMBER+r" max: ("+NUMBER+r")",log,re.M),dtype=float)
    interface = np.array(re.findall(r"^Interface Courant Number mean: "+NUMBER+r" max: ("+NUMBER+r")",log,re.M),dtype=float)
    assert len(logged_times)>0 and abs(logged_times[-1]-.5)<1e-10
    assert len(bulk)>=len(logged_times) and len(interface)>=len(logged_times)
    assert np.isfinite(bulk).all() and np.isfinite(interface).all()

    centres=field(case/"0/C","vector")
    nx=np.rint(centres[:,0]/(WIDTH/NX)-.5).astype(int)
    ny=np.rint(centres[:,1]/(HEIGHT/NY)-.5).astype(int)
    assert set(nx)==set(range(NX)) and set(ny)==set(range(NY))
    order=np.argsort(ny*NX+nx)
    assert np.array_equal((ny*NX+nx)[order],np.arange(CELLS))
    x=(np.arange(NX)+.5)*WIDTH/NX
    y=(np.arange(NY)+.5)*HEIGHT/NY
    assert np.allclose(centres[order,0].reshape(NY,NX),x,atol=1e-12,rtol=0)
    assert np.allclose(centres[order,1].reshape(NY,NX),y[:,None],atol=1e-12,rtol=0)

    gravity_text=(case/"constant/g").read_text(encoding="ascii")
    gravity_match=re.search(r"^\s*value\s+\(\s*0\s+-("+NUMBER+r")\s+0\s*\)\s*;",gravity_text,re.M)
    assert gravity_match and float(gravity_match[1])>0, "Expected downward gravity along the channel"
    gravity=float(gravity_match[1])
    ideal=2*.0707106*np.cos(np.deg2rad(45))/(999*gravity*WIDTH)
    result=[]
    fields=[]
    source_names=["input-sha256.txt","solver.log","check-mesh.log","0/C","constant/g"]
    for time in TIMES:
        label=f"{time:.10g}"
        directory=case/label
        clock=(directory/"uniform/time").read_text(encoding="ascii")
        step=int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-time)<1e-8
        alpha=field(directory/"alpha.water")[order].reshape(NY,NX)
        assert alpha.min()>=-ALPHA_TOL and alpha.max()<=1+ALPHA_TOL,(label,alpha.min(),alpha.max())
        assert np.any(alpha>.95) and np.any(alpha<.05)
        field(directory/"U","vector")
        field(directory/"p_rgh")
        centre=alpha[:,NX//2-1:NX//2+1].mean(axis=1)
        crossing=np.flatnonzero((centre[:-1]>=.5)&(centre[1:]<.5))
        assert len(crossing)>0 and centre[0]>=.5
        j=crossing[0]
        crossing_y=y[j]+(y[j+1]-y[j])*(.5-centre[j])/(centre[j+1]-centre[j])
        result.append({"time_s":time,"native_step":step,"equivalent_water_height_m":float(alpha.mean()*HEIGHT),
                       "central_alpha_half_height_m":float(crossing_y),"alpha_min":float(alpha.min()),
                       "alpha_max":float(alpha.max())})
        fields.append(alpha)
        source_names.extend(f"{label}/{name}" for name in ("uniform/time","alpha.water","U","p_rgh"))

    assert len(result)==10 and all(a["native_step"]<b["native_step"] for a,b in zip(result,result[1:]))
    assert result[-1]["native_step"]==len(logged_times)
    checks={"cells":CELLS,"native_steps":len(logged_times),"accepted_saved_states":len(result),
            "gravity_m_s2":gravity,"equilibrium_scale_m":float(ideal),"max_logged_bulk_courant":float(bulk.max()),
            "max_logged_interface_courant":float(interface.max()),"samples":result,
            "native_sha256":{name:sha(case/name) for name in source_names}}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    with (output/"rise-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(result[0]))
        writer.writeheader()
        writer.writerows(result)

    fig,(a,b)=plt.subplots(1,2,figsize=(10,4.3),gridspec_kw={"width_ratios":[1,2]},layout="constrained")
    panel=a.imshow(fields[-1],origin="lower",extent=(0,1,0,20),vmin=0,vmax=1,cmap="viridis",
                   interpolation="nearest",aspect="auto")
    a.set(xlabel="Gap coordinate (mm)",ylabel="Height (mm)",title="Water fraction at 0.50 s",ylim=(7,13))
    fig.colorbar(panel,ax=a,label="Water volume fraction")
    b.plot(TIMES,[1000*entry["equivalent_water_height_m"] for entry in result],"o-",label="Volume-equivalent height")
    b.plot(TIMES,[1000*entry["central_alpha_half_height_m"] for entry in result],"s-",label="Centreline interface")
    b.axhline(1000*ideal,linestyle="--",color="#485665",label="Static capillary scale")
    b.set(xlabel="Physical time (s)",ylabel="Height (mm)",title="Measured rise from saved fields")
    b.grid(alpha=.25)
    b.legend()
    fig.suptitle("Capillary rise between wetting walls | OpenFOAM",fontweight="bold")
    for extension in ("png","pdf"):
        fig.savefig(output/f"capillary-rise.{extension}",dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case","foam_reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    audit(args.case,args.foam_reader,args.output)
