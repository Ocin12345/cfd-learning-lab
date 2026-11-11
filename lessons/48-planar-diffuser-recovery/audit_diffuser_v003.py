
"""Audit OpenFOAM diffuser fields, boundary flux and pressure recovery."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NUMBER=r"[-+0-9.eE]+"
UIN=.1


def foam_values(text,name,count,kind):
    uniform=re.search(r"\b"+name+r"\s+uniform\s+([^;]+);",text)
    if uniform:
        values=np.fromstring(uniform[1].replace("("," ").replace(")"," "),sep=" ")
        values=np.tile(values,(count,1))
    else:
        field=re.search(r"\b"+name+r"\s+nonuniform\s+List<"+kind+r">\s+(\d+)\s*\((.*?)\)\s*;",text,re.S)
        assert field and int(field[1])==count,(name,count)
        values=np.fromstring(field[2].replace("("," ").replace(")"," "),sep=" ").reshape(count,-1)
    assert np.isfinite(values).all()
    return values[:,0] if kind=="scalar" else values


def audit(root,label,length,manifest):
    case=root/label
    for name,digest in manifest[label].items():
        assert hashlib.sha256((case/name).read_bytes()).hexdigest()==digest
    nx=110+round(length/.001)
    count=24*nx
    mesh=(case/"check-mesh.log").read_text()
    assert "Mesh OK." in mesh and re.search(rf"^\s*cells:\s*{count}\s*$",mesh,re.M)
    log=(case/"solver.log").read_text()
    assert re.search(r"^End\s*$",log,re.M)
    times=np.array([float(t) for t in re.findall(r"^Time = ("+NUMBER+r")s?\s*$",log,re.M)])
    assert np.array_equal(times,np.arange(1,801))
    terminal=log[log.rindex("Time = 800"):]
    continuity=np.array([float(t) for t in re.findall(r"sum local = ("+NUMBER+r")",terminal)])
    residuals=np.array([float(t) for t in re.findall(r"Final residual = ("+NUMBER+r")",terminal)])
    assert len(continuity) and continuity.max()<1e-8
    assert len(residuals)>=3 and residuals.max()<=1e-10
    c=foam_values((case/"0/C").read_text(),"internalField",count,"vector")
    bins=np.floor(c[:,0]/.001).astype(int)
    index=np.lexsort((c[:,1],bins)).reshape(nx,24)
    assert np.all(bins[index] == np.arange(nx)[:,None])
    x=c[index,0].mean(axis=1)
    u=foam_values((case/"800/U").read_text(),"internalField",count,"vector")
    p=foam_values((case/"800/p").read_text(),"internalField",count,"scalar")
    prior=foam_values((case/"750/U").read_text(),"internalField",count,"vector")
    change=float(np.sqrt(np.mean(np.sum((u-prior)**2,axis=1))))
    assert change<.001
    phi=(case/"800/phi").read_text()
    rates={}
    for name in ("inlet","outlet"):
        patch=re.search(r"\b"+name+r"\s*\{([^{}]*)\}",phi,re.S)
        assert patch,name
        rates[name]=float(foam_values(patch[1],"value",24,"scalar").sum())
    balance=abs(sum(rates.values()))/abs(rates["inlet"])
    assert rates["inlet"]<0<rates["outlet"] and balance<.001
    mean_p=p[index].mean(axis=1)
    mean_u=u[index,0].mean(axis=1)
    a=int(np.argmin(abs(x-.0445)))
    b=int(np.argmin(abs(x-(.0605+length))))
    cp=float((mean_p[b]-mean_p[a])/(.5*UIN**2))
    summary={"case":label,"cells":count,"expansion_length_m":length,"pressure_sample_x_m":[float(x[a]),float(x[b])],
             "pressure_recovery_coefficient":cp,"inlet_outlet_flux_m3_s":rates,"relative_flux_imbalance":balance,
             "last_50_iteration_velocity_change_rms_m_s":change,"terminal_local_continuity_max":float(continuity.max()),
             "terminal_linear_residual_max":float(residuals.max()),"minimum_axial_velocity_m_s":float(u[:,0].min()),
             "reverse_axial_flow_cells":int(np.count_nonzero(u[:,0]<-1e-6)),
             "source_log_sha256":hashlib.sha256((case/"solver.log").read_bytes()).hexdigest()}
    return summary,x,mean_p,mean_u,u[index,0]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("root","manifest","output"):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    manifest=json.loads(args.manifest.read_bytes())
    results=[audit(args.root,label,length,manifest) for label,length in (("short",.01),("long",.04))]
    args.output.mkdir()
    checks={d[0]["case"]:d[0] for d in results}
    (args.output/"checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    with (args.output/"pressure-sections.csv").open("x",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(("case","x_m","mean_kinematic_pressure_m2_s2","mean_axial_velocity_m_s"))
        for summary,x,p,u,_ in results:
            writer.writerows((summary["case"],*row) for row in zip(x,p,u))
    fig,ax=plt.subplots(2,2,figsize=(11,6),layout="constrained")
    fig.suptitle("Pressure recovery in two planar diffusers | OpenFOAM Foundation 14",fontweight="bold")
    minimum=min(float(d[-1].min()) for d in results)
    maximum=max(float(d[-1].max()) for d in results)
    for i,(summary,x,p,u,field) in enumerate(results):
        length=summary["expansion_length_m"]
        faces=np.arange(len(x)+1)*.001
        half=.005*(1+np.clip((faces-.05)/length,0,1))
        y=half[:,None]*np.linspace(-1,1,25)
        X=np.broadcast_to(faces[:,None],y.shape)
        m=ax[i,0].pcolormesh(X*1000,y*1000,field,cmap="viridis",vmin=minimum,vmax=maximum)
        ax[i,0].set(xlabel="x (mm)",ylabel="y (mm)",title=f"{length*1000:g} mm expansion | Axial velocity")
        ax[i,1].plot(x*1000,p/(.5*UIN**2))
        for site in summary["pressure_sample_x_m"]:
            ax[i,1].axvline(site*1000,ls=":",color=".4")
        ax[i,1].set(xlabel="x (mm)",ylabel="Normalized mean pressure",title=f"Pressure recovery = {summary['pressure_recovery_coefficient']:.3f}")
        ax[i,1].grid(alpha=.22)
    fig.colorbar(m,ax=list(ax[:,0]),label="Axial velocity (m/s)")
    fig.savefig(args.output/"diffuser-comparison.png",dpi=175)
    fig.savefig(args.output/"diffuser-comparison.pdf")
    print(json.dumps(checks,indent=2))
