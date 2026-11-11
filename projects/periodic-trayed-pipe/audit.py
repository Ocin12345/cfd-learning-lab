"""Audit the periodic wedge, tray faces and saved liquid redistribution."""
import argparse
from collections import Counter
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

NR, NY = 40, 120
CELLS = NR*NY
DR = DY = .000375
TIMES = np.arange(1,11)*.025
NUMBER = r"([-+\d.eE]+)"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, reader_path, out):
    read=runpy.run_path(str(reader_path))["foam_values"]
    sources=[]

    def native(name):
        sources.append(name)
        return (case/name).read_text("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest,name=line.split("  ",1)
        assert sha(case/name)==digest,name
    sources.extend(["pre-mesh-input-sha256.txt","pre-mesh-inputs-verified.log","zones.log","baffles.log"])
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*4800\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd")
    times=np.array(re.findall(r"^Time = "+NUMBER+r"s$",log,re.M),float)
    assert len(times)>0 and np.all(np.diff(times)>0) and abs(times[-1]-.25)<1e-8

    def field(name,kind="scalar"):
        values=read(native(name),"internalField",CELLS,kind)
        assert np.all(np.isfinite(values)),name
        return values

    c,v=field("0/C","vector"),field("0/Vc")
    assert c.shape==(CELLS,3) and np.all(v>0) and np.max(abs(c[:,2]))<1e-10
    ri, yi = np.floor(c[:,:2]/[DR,DY]).astype(int).T
    assert np.all((ri>=0)&(ri<NR)&(yi>=0)&(yi<NY))
    assert len(np.unique(yi*NR+ri))==CELLS
    assert np.allclose(c[:,1],(yi+.5)*DY,atol=1e-9,rtol=0)
    assert np.allclose(v,.025*((ri+1)**2-ri**2)*DR**2*DY,atol=1e-15,rtol=1e-8)

    def grid(values):
        result=np.empty((NY,NR))
        result[yi,ri]=values
        return result

    boundary=native("constant/polyMesh/boundary")
    patches={}
    for name,body in re.findall(r"(\w+)\s*\{([^{}]*nFaces[^{}]*)\}",boundary,re.S):
        count=int(re.search(r"\bnFaces\s+(\d+);",body)[1])
        start=int(re.search(r"\bstartFace\s+(\d+);",body)[1])
        patches[name]=(start,count)
    for name,other in (("top","bottom"),("bottom","top")):
        body=re.search(r"\b"+name+r"\s*\{([^}]+)\}",boundary,re.S)[1]
        assert re.search(r"\btype\s+cyclic;",body)
        assert re.search(r"\bneighbourPatch\s+"+other+r";",body)
    start,count=patches["baffles"]
    assert count>0 and count%2==0
    points=np.array([np.fromstring(line,sep=" ") for line in re.findall(
        r"^\(([^()]+)\)$",native("constant/polyMesh/points"),re.M)])
    faces=[np.fromstring(values,dtype=int,sep=" ") for values in re.findall(
        r"^\d+\(([^()]+)\)$",native("constant/polyMesh/faces"),re.M)]
    assert len(faces)==max(start+count for start,count in patches.values())
    baffle_points=[points[face] for face in faces[start:start+count]]
    keys=Counter(tuple(sorted(map(tuple,np.round(face,12)))) for face in baffle_points)
    assert len(keys)==count//2 and all(times_seen==2 for times_seen in keys.values())
    assert all(np.ptp(face[:,1])<1e-10 for face in baffle_points)
    levels=np.unique(np.round([face[0,1] for face in baffle_points],12))
    expected=(2*np.arange(6)+1)*.045/12
    assert len(levels)==6 and np.allclose(levels,expected,atol=1e-10,rtol=0)
    bounds=[]
    for level in levels:
        at_level=[face[:,0] for face in baffle_points if abs(face[0,1]-level)<1e-10]
        bounds.append((float(level),float(min(x.min() for x in at_level)),float(max(x.max() for x in at_level))))

    init=native("0/alpha.water")
    assert re.search(r"alphaMean\s+0\.25;",init)
    outside=c[:,0]>.0075
    initial=(~outside).astype(float)
    inventory=float(initial@v)
    initial_log_fraction=float(re.search(r"Phase-1 volume fraction = "+NUMBER,log)[1])
    assert abs(initial_log_fraction-inventory/v.sum())<1e-6

    saved_names=[f"{t:g}" for t in TIMES]
    actual={p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"0\.\d+",p.name)}
    assert actual==set(saved_names)
    rows=[]
    states={}
    for t,name in zip(TIMES,saved_names):
        clock=native(f"{name}/uniform/time")
        step=int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-t)<1e-8
        assert abs(times[step-1]-t)<1e-8
        a,u=field(f"{name}/alpha.water"),field(f"{name}/U","vector")
        field(f"{name}/p_rgh")
        phase_flux={}
        if (case/f"{name}/alphaPhi.water").is_file():
            phi=native(f"{name}/alphaPhi.water")
            for patch in ("top","bottom","baffles"):
                body=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",phi,re.S)[1]
                if re.search(r"\bvalue\s+",body):
                    values=read(body,"value",patches[patch][1],"scalar")
                    assert np.all(np.isfinite(values))
                    phase_flux[patch]=float(values.sum())
                    if patch=="baffles":
                        assert np.max(abs(values))<1e-12
        states[name]=grid(a)
        rows.append({"time_s":float(t),"step":step,"alpha_min":float(a.min()),"alpha_max":float(a.max()),
                     "inventory_m3":float(a@v),"relative_inventory_change":float(a@v/inventory-1),
                     "outside_fraction_of_initial_water":float(a[outside]@v[outside]/inventory),
                     "downward_water_bearing_cells":int(np.count_nonzero((a>=.5)&(u[:,1]<-.01))),
                     **{f"{patch}_phase_flux_m3_s":phase_flux.get(patch) for patch in ("top","bottom","baffles")}})
    assert rows[-1]["step"]==len(times)
    gates={"fraction":all(r["alpha_min"]>=-1e-6 and r["alpha_max"]<=1+1e-6 for r in rows),
           "inventory":max(abs(r["relative_inventory_change"]) for r in rows)<.005,
           "redistribution":any(r["outside_fraction_of_initial_water"]>=.05 and
                                r["downward_water_bearing_cells"]>=10 for r in rows)}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"steps":len(times),
            "saved_states":len(rows),"baffle_faces":count,"baffle_levels_and_radial_bounds_m":bounds,
            "initial_water_volume_m3":inventory,"samples":rows,
            "native_sha256":{name:sha(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"trayed-pipe-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig,axes=plt.subplots(1,3,figsize=(10.6,5.1),layout="constrained",gridspec_kw={"width_ratios":[1,1,1.6]})
    for ax,values,title in ((axes[0],grid(initial),"Initial water"),(axes[1],states["0.25"],"Water at 0.25 s")):
        im=ax.pcolormesh(np.arange(NR+1)*DR*1000,np.arange(NY+1)*DY*1000,values,
                         cmap="Blues",vmin=0,vmax=1,shading="flat")
        for level,lo,hi in bounds:
            ax.plot([1000*lo,1000*hi],[1000*level]*2,color="#bc562c",lw=1.5)
        ax.set(xlim=(0,15),ylim=(0,45),xlabel="Radius (mm)",ylabel="Axial position (mm)",title=title,aspect="equal")
    fig.colorbar(im,ax=axes[:2],label="Water fraction",shrink=.77)
    axes[2].plot([0]+[r["time_s"] for r in rows],[0]+[100*r["outside_fraction_of_initial_water"] for r in rows],
                 "o-",ms=3,color="#155e88")
    axes[2].set(xlabel="Time (s)",ylabel="Percent of initial water inventory",
                title="Water outside the initial radius")
    axes[2].grid(alpha=.2)
    fig.suptitle("Liquid redistribution in a periodic trayed pipe",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"periodic-trayed-pipe.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case","reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)