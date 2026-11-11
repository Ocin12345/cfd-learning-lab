"""Check the flat-plate flow and compare native profiles to Blasius."""
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

NX,NY=180,96
CELLS=NX*NY
NU=.00025
NUMBER=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def blasius():
    eta=np.linspace(0,10,2001)
    value=np.array([0.,0.,.332057336215196])
    velocity=[0.]
    def rhs(f):
        return np.array([f[1],f[2],-.5*f[0]*f[2]])
    h=eta[1]-eta[0]
    for _ in eta[1:]:
        a=rhs(value)
        b=rhs(value+.5*h*a)
        c=rhs(value+.5*h*b)
        d=rhs(value+h*c)
        value=value+h*(a+2*b+2*c+d)/6
        velocity.append(value[1])
    assert abs(velocity[-1]-1)<1e-8
    return eta,np.array(velocity)


def audit(case,reader,out):
    read=runpy.run_path(str(reader))["foam_values"]
    sources={}

    def native(name):
        raw=(case/name).read_bytes()
        sources[name]=hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest,name=line.split("  ",1)
        native(name)
        assert sources[name]==digest,name
    for name in ("pre-mesh-input-sha256.txt","pre-mesh-inputs-verified.log","inputs-verified.log",
                 "mesh.log","centres.log","volumes.log"):
        native(name)
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*17280\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd")
    steps=list(re.finditer(r"^Time = (\d+)s?$",log,re.M))
    end=int(steps[-1][1])
    assert 1<=end<=1000 and [int(s[1]) for s in steps]==list(range(1,end+1))
    final=dict((name,float(value)) for name,value in re.findall(
        r"Solving for (Ux|Uy|Uz|p), Initial residual = ("+NUMBER+r")",log[steps[-1].end():]))
    assert {"Ux","Uy","p"}<=final.keys() and all(np.isfinite(v) and v>=0 for v in final.values())
    continuity=[float(v) for v in re.findall(r"time step continuity errors : sum local = ("+NUMBER+r")",log)]
    assert len(continuity)==end and np.isfinite(continuity).all()

    def field(name,kind="scalar"):
        a=read(native(name),"internalField",CELLS,kind)
        assert a.shape==((CELLS,3) if kind=="vector" else (CELLS,)) and np.isfinite(a).all(),name
        return a

    c,v=field("0/C","vector"),field("0/Vc")
    assert (v>0).all() and abs(v.sum()-.135)<1e-10
    assert np.max(abs(c[:,2]))<1e-10
    x,y=np.unique(c[:,0]),np.unique(c[:,1])
    assert len(x)==NX and len(y)==NY and len(np.unique(c[:,:2],axis=0))==CELLS
    ix,iy=np.searchsorted(x,c[:,0]),np.searchsorted(y,c[:,1])
    saves=sorted(int(p.name) for p in case.iterdir() if p.is_dir() and p.name.isdigit() and int(p.name)>0)
    assert saves==sorted(set(range(25,end+1,25))|{end})
    for step in saves:
        clock=native(f"{step}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",clock,re.M)[1])==step
        assert float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])==step
    u,p=field(f"{end}/U","vector"),field(f"{end}/p")
    assert np.max(abs(u[:,2]))<1e-9
    flux_text=native(f"{end}/phi")
    flux={}
    impermeable=0.
    for patch,count in (("inlet",96),("outlet",96),("farfield",180),("upstream",20),("plate",160)):
        body=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",flux_text,re.S)[1]
        values=read(body,"value",count,"scalar")
        assert values.shape==(count,) and np.isfinite(values).all()
        flux[patch]=float(values.sum())
        if patch not in ("inlet","outlet"):
            impermeable=max(impermeable,float(np.max(abs(values))))
    assert flux["inlet"]<0<flux["outlet"]
    imbalance=abs(sum(flux.values()))/abs(flux["inlet"])
    eta_ref,u_ref=blasius()
    ref99=float(np.interp(.99,u_ref,eta_ref))
    data=[]
    profiles=[]
    for desired in (.5,1,1.5):
        xvalue=float(x[np.argmin(abs(x-desired))])
        use=np.flatnonzero(c[:,0]==xvalue)
        use=use[np.argsort(c[use,1])]
        ys=c[use,1]
        assert len(use)==NY and np.array_equal(ys,y)
        ue=float(np.mean(u[use,0][ys>=.45]))
        assert ue>0
        norm=u[use,0]/ue
        eta=ys*np.sqrt(ue/(NU*xvalue))
        mask=eta<=6
        assert mask.sum()>=10 and eta[mask].max()<=eta_ref[-1]
        expected=np.interp(eta[mask],eta_ref,u_ref)
        rms=float(np.sqrt(np.mean((norm[mask]-expected)**2)))
        crossing=np.flatnonzero((norm[:-1]<.99)&(norm[1:]>=.99))
        assert len(crossing)>0
        j=int(crossing[0])
        delta=float(ys[j]+(.99-norm[j])*(ys[j+1]-ys[j])/(norm[j+1]-norm[j]))
        ref_delta=ref99*np.sqrt(NU*xvalue/ue)
        data.append({"x_m":xvalue,"outer_ux_m_s":ue,"local_reynolds_number":ue*xvalue/NU,
                     "delta99_m":delta,"blasius_delta99_m":float(ref_delta),"blasius_u_over_ue_rms":rms,
                     "bracket_y_lo_m":float(ys[j]),"bracket_y_hi_m":float(ys[j+1])})
        profiles.append((xvalue,eta,norm))
    gates={"native_convergence":f"PIMPLE solution converged in {end} iterations" in log,
           "residuals":max(final.values())<=1e-5,
           "mass_balance":imbalance<=.001,"impermeable_patches":impermeable<=1e-12,
           "profile_rms":max(row["blasius_u_over_ue_rms"] for row in data)<=.05,
           "positive_layer_growth":all(0<a["delta99_m"]<b["delta99_m"] for a,b in zip(data,data[1:]))}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"iterations":end,"saved_states":len(saves),
            "terminal_initial_residuals":final,"max_absolute_wall_or_slip_flux_m3_s":impermeable,
            "boundary_flow_m3_s":flux,"relative_mass_imbalance":imbalance,
            "terminal_p_range_m2_s2":[float(p.min()),float(p.max())],
            "stations":data,"native_sha256":sources}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"plate-profiles.csv").open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("x_m","eta","native_ux_over_ue","blasius_ux_over_ue"))
        for station,eta,norm in profiles:
            for et,val in zip(eta,norm):
                writer.writerow((station,et,val,float(np.interp(et,eta_ref,u_ref)) if et<=eta_ref[-1] else ""))
    with (out/"plate-stations.csv").open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)

    fig=plt.figure(figsize=(10.8,6.8),layout="constrained")
    spec=fig.add_gridspec(2,2,height_ratios=[.85,1.2])
    ax=fig.add_subplot(spec[0,:])
    grid=np.empty((NY,NX))
    grid[iy,ix]=u[:,0]
    yedges=np.r_[0,np.cumsum(np.bincount(iy,weights=v,minlength=NY)/(.1*2.25))]
    assert abs(yedges[-1]-.6)<1e-9
    im=ax.pcolormesh(np.linspace(-.25,2,NX+1),yedges,grid,cmap="viridis",shading="flat")
    for row in data:
        ax.axvline(row["x_m"],color="white",ls="--",lw=.7)
    ax.plot([0,2],[0,0],color="#606974",lw=3)
    ax.set(title="Native streamwise velocity near the plate",xlabel="x (m)",ylabel="y (m)",ylim=(0,.22))
    fig.colorbar(im,ax=ax,label="Streamwise velocity (m/s)")
    ax=fig.add_subplot(spec[1,0])
    ax.plot(eta_ref,u_ref,color="black",ls="--",label="Blasius")
    for station,eta,norm in profiles:
        ax.plot(eta,norm,"o-",ms=2,label=f"x = {station:.3f} m")
    ax.set(title="Velocity profiles",xlabel="Similarity coordinate η",ylabel="Streamwise velocity / outer velocity",
           xlim=(0,6),ylim=(0,1.04))
    ax.legend(frameon=False)
    ax.grid(alpha=.2)
    ax=fig.add_subplot(spec[1,1])
    positions=[r["x_m"] for r in data]
    ax.plot(positions,[1000*r["delta99_m"] for r in data],"o-",label="Native bracket")
    ax.plot(positions,[1000*r["blasius_delta99_m"] for r in data],"o--",label="Blasius")
    ax.set(title="99% boundary-layer thickness",xlabel="x (m)",ylabel="Thickness (mm)")
    ax.legend(frameon=False)
    ax.grid(alpha=.2)
    fig.suptitle("Laminar boundary layer on a flat plate",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"flat-plate-layer.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        parser.add_argument(name,type=Path)
    a=parser.parse_args()
    raise SystemExit(0 if audit(a.case,a.reader,a.output) else 2)
