"""Audit native steady flow and incident radiation in the hot-box room."""
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

CELLS=15523
NUMBER=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
LIMITS={"Ux":1e-3,"Uy":1e-3,"Uz":1e-3,"h":1e-3,"G":1e-3,"p_rgh":1e-2,"k":1e-3,"epsilon":1e-3}


def audit(case,reader,out):
    read=runpy.run_path(str(reader))["foam_values"]
    sources={}

    def native(name):
        raw=(case/name).read_bytes()
        sources[name]=hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest,name=line.split("  ",1)
        assert hashlib.sha256((case/name).read_bytes()).hexdigest()==digest,name
        sources[name]=digest
    for name in ("pre-mesh-input-sha256.txt","pre-mesh-inputs-verified.log","inputs-verified.log",
                 "mesh.log","centres.log","volumes.log","inherited-input-sha256.txt","inherited-inputs-verified.log"):
        native(name)
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*15523\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd") and "P1" in log
    prefix=native("initial-1000-solver.log")
    assert prefix.rstrip().endswith("\nEnd") and "P1" in prefix
    assert re.findall(r"^Time = (\d+)s$",prefix,re.M)==[str(i) for i in range(1,1001)]
    first=re.search(r"^Time = 1001s$",log,re.M)
    assert first
    log=prefix.rstrip()+"\n\n"+log[first.start():]
    steps=list(re.finditer(r"^Time = (\d+)s?$",log,re.M))
    end=int(steps[-1][1])
    assert 1001<=end<=1250 and [int(m[1]) for m in steps]==list(range(1,end+1))
    residuals=[]
    for i,m in enumerate(steps):
        next_step=steps[i+1].start() if i+1<len(steps) else len(log)
        block=log[m.end():next_step]
        values=re.findall(r"Solving for (Ux|Uy|Uz|h|G|p_rgh|k|epsilon), Initial residual = ("+NUMBER+r")",block)
        row={name:float(value) for name,value in values}
        assert set(row)==set(LIMITS) and all(np.isfinite(v) and v>=0 for v in row.values())
        residuals.append({"iteration":i+1,**row})

    def field(name,kind="scalar"):
        a=read(native(name),"internalField",CELLS,kind)
        assert a.shape==((CELLS,3) if kind=="vector" else (CELLS,)) and np.isfinite(a).all(),name
        return a

    c,v=field("0/C","vector"),field("0/Vc")
    assert (v>0).all() and abs(v.sum()-119.5)<1e-8
    assert np.all((c>0)&(c<[10,6,2]))
    x,y,z=c.T
    solid=(x>.5)&(x<1.5)&(y>.5)&(y<1.5)&(z<.5)
    assert not solid.any()
    saves=sorted(int(p.name) for p in case.iterdir() if p.is_dir() and p.name.isdigit() and int(p.name)>0)
    assert saves==sorted(set(range(1025,end+1,25))|{1000,end})
    for step in saves:
        clock=native(f"{step}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",clock,re.M)[1])==step
        assert float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])==step

    t,g,u,p,k,epsilon=(field(f"{end}/{name}","vector" if name=="U" else "scalar")
                     for name in ("T","G","U","p","k","epsilon"))
    low=(z>.5)&(z<1)
    above=low&(x>.5)&(x<1.5)&(y>.5)&(y<1.5)
    far=low&(x>7)
    def mean(a,mask):
        return float(np.average(a[mask],weights=v[mask]))
    tnear,tfar=mean(t,above),mean(t,far)
    gnear,gfar=mean(g,above),mean(g,far)

    phi=native(f"{end}/phi")
    boundary=native("constant/polyMesh/boundary")
    flux={}
    max_wall_flux=0.
    for patch in ("box","floor","ceiling","fixedWalls"):
        zone=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",boundary,re.S)[1]
        assert re.search(r"\btype\s+wall;",zone)
        count=int(re.search(r"\bnFaces\s+(\d+);",zone)[1])
        body=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",phi,re.S)[1]
        a=read(body,"value",count,"scalar")
        assert a.shape==(count,) and np.isfinite(a).all()
        flux[patch]=float(a.sum())
        max_wall_flux=max(max_wall_flux,float(np.max(abs(a))))
    gates={"native_convergence":f"PIMPLE solution converged in {end} iterations" in log,
           "residuals":all(residuals[-1][name]<=limit for name,limit in LIMITS.items()),
           "thermodynamic_fields":bool(p.min()>0 and 299<=t.min()<=t.max()<=501 and g.min()>=0),
           "turbulence_fields":bool(k.min()>=0 and epsilon.min()>=0),
           "radiation_contrast":gfar>0 and gnear>=1.1*gfar,
           "temperature_contrast":tnear-tfar>=1,
           "closed_walls":max_wall_flux<=1e-12}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"iterations":end,"continued_iterations":end-1000,"saved_states":39+len(saves),
            "terminal_initial_residuals":residuals[-1],"temperature_min_k":float(t.min()),
            "temperature_max_k":float(t.max()),"max_speed_m_s":float(np.linalg.norm(u,axis=1).max()),
            "above_box_temperature_k":tnear,"far_temperature_k":tfar,
            "above_box_incident_radiation_w_m2":gnear,"far_incident_radiation_w_m2":gfar,
            "closed_wall_fluxes_kg_s":flux,"max_absolute_wall_face_flux":max_wall_flux,
            "native_sha256":sources}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"radiating-room-history.csv").open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(residuals[0]))
        writer.writeheader()
        writer.writerows(residuals)

    xe=np.r_[np.linspace(0,.5,4)[:-1],np.linspace(.5,1.5,7)[:-1],np.linspace(1.5,10,41)]
    ze=np.r_[np.linspace(0,.5,4)[:-1],np.linspace(.5,2,9)]
    pick=abs(y-(1-1/12))<1e-8
    ix,iz=np.searchsorted(xe,x[pick])-1,np.searchsorted(ze,z[pick])-1
    assert pick.sum()==49*11-6*3 and len(np.unique(iz*49+ix))==pick.sum()
    np.testing.assert_allclose(x[pick],(xe[ix]+xe[ix+1])/2,atol=1e-9,rtol=0)
    np.testing.assert_allclose(z[pick],(ze[iz]+ze[iz+1])/2,atol=1e-9,rtol=0)
    fig,axs=plt.subplots(3,1,figsize=(10.8,8.6),layout="constrained",
                         gridspec_kw={"height_ratios":[1,1,1.35]})
    for ax,data,name,cmap,label in ((axs[0],t,"Air temperature and in-plane velocity","inferno","Temperature (K)"),
                                  (axs[1],g,"Incident radiation","magma","Incident radiation (W/m²)")):
        grid=np.full((11,49),np.nan)
        grid[iz,ix]=data[pick]
        im=ax.pcolormesh(xe,ze,grid,cmap=cmap,shading="flat")
        ax.add_patch(Rectangle((.5,0),1,.5,facecolor="#d9dee5",edgecolor="#374350",lw=.8))
        ax.set(title=name,xlabel="x (m)",ylabel="z (m)",aspect="equal")
        fig.colorbar(im,ax=ax,label=label)
    use=pick.copy()
    indices=np.flatnonzero(pick)
    use[indices[(ix%3!=0)|(iz%2!=0)]]=False
    axs[0].quiver(x[use],z[use],u[use,0],u[use,2],color="white",angles="xy",scale_units="xy",scale=4)
    ax=axs[2]
    for eq in ("p_rgh","h","G","Uz","k","epsilon"):
        ax.semilogy(np.arange(1,end+1),np.maximum([r[eq] for r in residuals],1e-20),label=eq)
    ax.set(title="Steady iteration history",xlabel="Iteration",ylabel="Initial scaled residual")
    ax.legend(ncol=3,frameon=False)
    ax.grid(alpha=.2)
    fig.suptitle("Buoyant flow and radiation from a hot box",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"radiating-room.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        parser.add_argument(name,type=Path)
    a=parser.parse_args()
    raise SystemExit(0 if audit(a.case,a.reader,a.output) else 2)
