"""Audit the ventilated room and its age-of-air solve on the saved flow."""
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

CELLS = 53820
NUMBER = r"([-+\d.eE]+)"
AXES = [np.r_[np.linspace(0,3.8,36),3.9,4.,4.1,4.2],
        np.r_[np.linspace(0,1.65,22),np.linspace(1.65,1.95,5)[1:],np.linspace(1.95,3.6,22)[1:]],
        np.r_[np.linspace(0,2.05,21),2.15,2.25,np.linspace(2.25,3.,9)[1:]]]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, reader_path, output):
    read = runpy.run_path(str(reader_path))["foam_values"]
    sources = []

    def native(name):
        sources.append(name)
        return (case/name).read_text("ascii")

    for manifest in ("input-sha256.txt","frozen-flow-sha256.txt"):
        for line in native(manifest).splitlines():
            expected,name=line.split("  ",1)
            assert sha(case/name)==expected,name
    mesh,log,age_log=native("check-mesh.log"),native("solver.log"),native("age.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*53820\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd") and age_log.rstrip().endswith("\nEnd")
    assert "Selecting RAS turbulence model kEpsilon" in log
    split=re.split(r"(?m)^Time = (\d+)s\s*$",log)
    actual=np.array(split[1::2],int)
    endpoint=len(actual)
    assert 0<endpoint<=2000 and np.array_equal(actual,np.arange(1,endpoint+1))
    convergence=re.search(r"SIMPLE solution converged in (\d+) iterations",log)
    assert (convergence and int(convergence[1])==endpoint) or endpoint==2000
    age_residuals=[float(value) for value in re.findall(
        r"Solving for age, Initial residual = "+NUMBER+r", Final residual = "+NUMBER,age_log)
        for value in value[1:]]
    assert age_residuals and np.all(np.isfinite(age_residuals))
    last_continuity=float(re.findall(r"sum local = "+NUMBER,split[-1])[-1])

    def field(name,kind="scalar"):
        values=read(native(name),"internalField",CELLS,kind)
        assert np.all(np.isfinite(values)),name
        return values

    c, v=field("0/C","vector"),field("0/Vc")
    shape=tuple(len(edges)-1 for edges in AXES)
    ijk=tuple(np.searchsorted(edges,c[:,i])-1 for i,edges in enumerate(AXES))
    expected_volumes=np.ones(CELLS)
    for i,(edges,index) in enumerate(zip(AXES,ijk)):
        assert np.all((index>=0)&(index<shape[i]))
        assert np.allclose(c[:,i],(edges[index]+edges[index+1])/2,rtol=0,atol=1e-5)
        expected_volumes*=edges[index+1]-edges[index]
    assert np.allclose(v,expected_volumes,rtol=5e-6,atol=1e-12)
    assert len(np.unique(np.ravel_multi_index(ijk,shape)))==CELLS

    saved=sorted(int(path.name) for path in case.iterdir() if path.is_dir() and re.fullmatch(r"[1-9]\d*",path.name))
    assert saved==sorted(set(range(250,endpoint+1,250))|{endpoint})
    previous=None
    rows=[]
    for i in saved:
        clock=native(f"{i}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",clock,re.M)[1])==i
        assert float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])==i
        u=field(f"{i}/U","vector")
        field(f"{i}/p")
        k,epsilon=field(f"{i}/k"),field(f"{i}/epsilon")
        change=None if previous is None else float(np.sqrt(np.average(np.sum((u-previous)**2,axis=1),weights=v)))
        rows.append({"iteration":i,"min_k_m2_s2":float(k.min()),"min_epsilon_m2_s3":float(epsilon.min()),
                     "velocity_rms_change_m_s":change})
        previous=u

    boundary=native("constant/polyMesh/boundary")
    phi=native(f"{endpoint}/phi")
    assert re.search(r"dimensions\s+\[0 3 -1 0 0 0 0\];",phi)
    flux={}
    for patch in ("inlet","outlet"):
        b=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",boundary)[1]
        count=int(re.search(r"nFaces\s+(\d+);",b)[1])
        assert count==8
        block=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",phi)[1]
        values=read(block,"value",count,"scalar")
        assert values.shape==(count,) and np.all(np.isfinite(values))
        flux[patch]=float(values.sum())
    assert flux["inlet"]<0<flux["outlet"]
    mismatch=abs(flux["inlet"]+flux["outlet"])/abs(flux["inlet"])
    turnover=float(v.sum()/abs(flux["inlet"]))
    age=field(f"{endpoint}/age")
    x,y,z=c.T
    near=(x<.5)&(y>=1.45)&(y<=2.15)&(z>=1.9)&(z<=2.4)
    lower=(x>=2)&(x<=3.6)&(y>=.4)&(y<=3.2)&(z>=.3)&(z<=1.7)
    assert min(near.sum(),lower.sum())>=10
    near_mean=float(np.average(age[near],weights=v[near]))
    lower_mean=float(np.average(age[lower],weights=v[lower]))
    age_grid=np.empty(shape)
    age_grid[ijk]=age

    profiles=[]
    for i in (1,2):
        files=list((case/f"postProcessing/probes{i}").glob("*/age"))
        assert len(files)==1
        body=native(files[0].relative_to(case).as_posix())
        config=native(f"system/probes{i}")
        points=np.array([np.fromstring(s,sep=" ") for s in re.findall(r"^\s*\(([^()]+)\)\s*$",config,re.M)])
        actual_points=np.array([np.fromstring(s,sep=" ") for s in re.findall(r"^# Probe \d+ \(([^)]+)\)",body,re.M)])
        assert points.shape==(30,3) and np.array_equal(points,actual_points)
        values=np.atleast_2d(np.loadtxt(files[0]))
        assert values.shape==(1,31) and values[0,0]==endpoint
        ages=values[0,1:]
        checked=0
        for point,value in zip(points,ages):
            if not np.isfinite(value) or abs(value)>1e20:
                assert point[2] in (0,3)
                continue
            adjacent=[np.flatnonzero((edges[:-1]<=coordinate+1e-10)&(edges[1:]>=coordinate-1e-10))
                      for edges,coordinate in zip(AXES,point)]
            assert all(len(a)>0 for a in adjacent)
            candidates=age_grid[np.ix_(*adjacent)]
            assert np.any(np.isclose(value,candidates,rtol=5e-6,atol=1e-4)),(point,value,candidates)
            checked+=1
        assert checked>=28
        profiles.append((points,ages))

    gates={"flow_stability":bool(convergence) or (change is not None and change<.02),
           "flow_turbulence":all(row["min_k_m2_s2"]>=0 and row["min_epsilon_m2_s3"]>0 for row in rows),
           "continuity":last_continuity<1e-4 and mismatch<1e-4,
           "age_equation":age_residuals[-1]<1e-5 and float(age.min())>=0,
           "age_contrast":lower_mean-near_mean>=.1*turnover}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"iterations":endpoint,
            "native_convergence":bool(convergence),"saved_states":len(saved),"terminal_local_continuity":last_continuity,
            "volume_m3":float(v.sum()),"volume_flux_m3_s":flux,"flux_mismatch_relative":mismatch,
            "nominal_turnover_s":turnover,"age_final_linear_residual":age_residuals[-1],
            "age_min_s":float(age.min()),"age_max_s":float(age.max()),"volume_weighted_age_s":float(np.average(age,weights=v)),
            "near_supply_cells":int(near.sum()),"lower_interior_cells":int(lower.sum()),
            "near_supply_mean_age_s":near_mean,"lower_interior_mean_age_s":lower_mean,
            "samples":rows,"native_sha256":{name:sha(case/name) for name in sources}}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (output/"age-profiles.csv").open("x",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(("x_m","y_m","height_m","age_s"))
        for points,ages in profiles:
            for point,value in zip(points,ages):
                writer.writerow((*point,float(value) if np.isfinite(value) and abs(value)<1e20 else ""))
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.6),layout="constrained")
    section=age_grid[:,[22,23],:].mean(axis=1)
    im=axes[0].pcolormesh(AXES[0],AXES[2],section.T,cmap="viridis",shading="flat")
    for points,_ in profiles:
        axes[0].axvline(points[0,0],ls="--",color="white",lw=1)
    axes[0].set(xlabel="Distance x (m)",ylabel="Height z (m)",title="Age of air at the centre section",aspect="equal")
    fig.colorbar(im,ax=axes[0],label="Age (s)",shrink=.8)
    for points,ages in profiles:
        use=np.isfinite(ages)&(abs(ages)<1e20)
        axes[1].plot(ages[use],points[use,2],"o-",ms=3,label=f"x = {points[0,0]:g} m")
    axes[1].axvline(turnover,ls="--",color="#555",label="Nominal turnover time")
    axes[1].set(xlabel="Age (s)",ylabel="Height z (m)",title="Native vertical probe profiles")
    axes[1].grid(alpha=.2)
    axes[1].legend(frameon=False)
    fig.suptitle("Age of air in a ventilated room",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(output/f"ventilation-age.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("case","reader","output"):
        parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)