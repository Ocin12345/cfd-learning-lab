"""Audit and plot cellular convection in the stock horizontal layer."""
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

NX, NY, CELLS = 90, 10, 900
SAVES = np.arange(50,1001,50)


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
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*900\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd")
    assert "Boussinesq" in log and "Selecting laminar" in log
    steps=np.array(re.findall(r"^Time = (\d+)s$",log,re.M),int)
    assert np.array_equal(steps,np.arange(1,1001))
    terminal=log.rsplit("Time = 1000s",1)[1]
    continuity=float(re.findall(r"sum local = ([-+\d.eE]+)",terminal)[-1])

    def field(name,kind="scalar"):
        values=read(native(name),"internalField",CELLS,kind)
        assert np.all(np.isfinite(values)),name
        return values

    c,v=field("0/C","vector"),field("0/Vc")
    assert np.allclose(v,.02,rtol=1e-9,atol=1e-12)
    assert np.allclose(c[:,2],0,rtol=0,atol=1e-12)
    ij=np.rint(c[:,:2]/.1-.5).astype(int)
    assert np.all((ij>=0)&(ij<[NX,NY]))
    assert np.allclose(c[:,:2],(ij+.5)*.1,rtol=0,atol=1e-9)
    key=ij[:,1]*NX+ij[:,0]
    assert np.array_equal(np.sort(key),np.arange(CELLS))
    order=np.argsort(key)
    assert np.all(field("0/T")==300) and np.all(field("0/U","vector")==[1e-4,0,0])
    rows=[]
    temperatures={}
    for i in SAVES:
        clock=native(f"{i}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",clock,re.M)[1])==i
        assert float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])==i
        assert float(re.search(r"^deltaT\s+([^;]+);",clock,re.M)[1])==1
        u,temp=field(f"{i}/U","vector"),field(f"{i}/T")
        field(f"{i}/p_rgh")
        grid=temp[order].reshape(NY,NX)
        temperatures[int(i)]=grid
        vgrid=u[order,1].reshape(NY,NX)
        mid=vgrid[4:6].mean(axis=0)
        active=mid[abs(mid)>.001]
        changes=int(np.count_nonzero(np.diff(np.sign(active)))) if len(active) else 0
        rows.append({"time_s":int(i),"temperature_min_k":float(temp.min()),"temperature_max_k":float(temp.max()),
                     "vertical_speed_rms_m_s":float(np.sqrt(np.mean(u[:,1]**2))),
                     "midheight_vertical_min_m_s":float(mid.min()),"midheight_vertical_max_m_s":float(mid.max()),
                     "midheight_active_sign_changes":changes,
                     "fraction_cells_above_vertical_threshold":float(np.mean(abs(u[:,1])>.001)),
                     "out_of_plane_max_m_s":float(abs(u[:,2]).max()),
                     "floor_gradient_ratio":float((301-grid[0]).mean()/.05)})
    final=rows[-1]
    gates={"temperature":all(r["temperature_min_k"]>=299.99 and r["temperature_max_k"]<=301.01 for r in rows),
           "planarity":all(r["out_of_plane_max_m_s"]<1e-9 for r in rows),
           "continuity":np.isfinite(continuity) and continuity<1e-4,
           "roll_pattern":final["midheight_vertical_min_m_s"]<-.001 and
                          final["midheight_vertical_max_m_s"]>.001 and final["midheight_active_sign_changes"]>=2}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"steps":1000,"saved_states":len(rows),
            "terminal_local_continuity":continuity,"samples":rows,
            "native_sha256":{name:sha(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"convection-roll-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig,axes=plt.subplots(3,1,figsize=(11,6.8),layout="constrained")
    xedges,yedges=np.arange(NX+1)*.1,np.arange(NY+1)*.1
    x,y=xedges[:-1]+.05,yedges[:-1]+.05
    for ax,time in zip(axes[:2],(50,1000)):
        p=ax.pcolormesh(xedges,yedges,temperatures[time],shading="flat",cmap="coolwarm",vmin=300,vmax=301)
        ax.set(xlim=(0,9),ylim=(0,1),ylabel="Height (m)",title=f"Temperature at {time:,} s")
        fig.colorbar(p,ax=ax,label="Temperature (K)",shrink=.85)
    ux=u[order,0].reshape(NY,NX)
    speed=np.hypot(ux,vgrid)
    p=axes[2].pcolormesh(xedges,yedges,speed,cmap="viridis",shading="flat")
    axes[2].streamplot(x,y,ux,vgrid,color="white",density=(2,.6),linewidth=.65,arrowsize=.75)
    axes[2].set(xlim=(0,9),ylim=(0,1),xlabel="Distance (m)",ylabel="Height (m)",
                title="Velocity magnitude and streamlines at 1,000 s")
    fig.colorbar(p,ax=axes[2],label="Speed (m/s)",shrink=.85)
    fig.suptitle("Cellular convection in a layer heated from below",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"rayleigh-benard-rolls.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("case","reader","output"):
        parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)