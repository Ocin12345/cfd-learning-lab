"""Check native valve motion, external flux and branch selection."""
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

CELLS=30400
TIMES=np.arange(1,21)/20
NUMBER=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
PORTS=("inlet","lowerOutlet","upperOutlet")
WALLS=("pipeWalls","ballWalls","pipeNonCouple","ballNonCouple")


def audit(case,reader,out):
    read=runpy.run_path(str(reader))["foam_values"]
    sources={}

    def native(name):
        raw=(case/name).read_bytes()
        sources[name]=hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest,name=line.split("  ",1)
        data=(case/name).read_bytes()
        sources[name]=hashlib.sha256(data).hexdigest()
        assert sources[name]==digest,name
    for name in ("pre-mesh-input-sha256.txt","pre-mesh-inputs-verified.log","inputs-verified.log","mesh.log",
                 "zones.log","initial-rotation.log","couples.log","centres.log","volumes.log"):
        native(name)
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*30400\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd") and "oscillatingRotatingMotion" in log
    times=np.array(re.findall(r"^Time = ("+NUMBER+r")s?$",log,re.M),float)
    assert len(times)>0 and np.isfinite(times).all() and np.all(np.diff(times)>0) and times[-1]==1
    assert "Calculating couplings" in log

    zone_file=native("constant/polyMesh/cellZones")
    zones={}
    for name,size in (("ball",7600),("pipes",22800)):
        block=re.search(r"\b"+name+r"\s*\{\s*cellLabels\s+List<label>\s+"+str(size)+r"\s*\(([^()]+)\)",
                        zone_file,re.S)[1]
        data=np.fromstring(block,sep=" ",dtype=int)
        assert data.size==size
        zones[name]=data
    np.testing.assert_array_equal(np.sort(np.r_[zones["ball"],zones["pipes"]]),np.arange(CELLS))
    boundary=native("constant/polyMesh/boundary")
    counts={}
    for patch in (*PORTS,*WALLS):
        body=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",boundary,re.S)[1]
        counts[patch]=int(re.search(r"\bnFaces\s+(\d+);",body)[1])
        if patch in PORTS:
            assert counts[patch]==380
        else:
            assert re.search(r"\btype\s+wall;",body)
    assert "nonConformalCyclic_on_pipeNonCouple" in boundary
    assert "nonConformalCyclic_on_ballNonCouple" in boundary

    def field(name,kind="scalar"):
        a=read(native(name),"internalField",CELLS,kind)
        assert a.shape==((CELLS,3) if kind=="vector" else (CELLS,)) and np.isfinite(a).all(),name
        return a

    c0,v0=field("0/C","vector"),field("0/Vc")
    assert np.all(v0>0)
    reports={}
    for patch in PORTS:
        report=native(f"postProcessing/{patch}FlowRate/0/surfaceFieldValue.dat")
        assert re.search(r"^# Selection\s*:\s*patch\("+patch+r"\)",report,re.M)
        assert re.search(r"^# Faces\s*:\s*380",report,re.M) and "sum(phi)" in report
        table=np.array([[float(v) for v in line.split()] for line in report.splitlines()
                        if line.strip() and not line.lstrip().startswith("#")])
        assert table.shape==(len(times)+1,2) and np.isfinite(table).all()
        np.testing.assert_array_equal(table[:,0],np.r_[0,times])
        reports[patch]=table[:,1]
    peak=float(max(np.max(abs(report)) for report in reports.values()))
    assert peak>0

    actual={p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?",p.name)}
    assert actual=={"0",*(f"{t:g}" for t in TIMES)}
    rows=[]
    states={}
    for t in TIMES:
        name=f"{t:g}"
        clock=native(f"{name}/uniform/time")
        step=int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert 1<=step<=len(times) and abs(times[step-1]-t)<1e-9
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-t)<1e-10
        c,v=field(f"{name}/C","vector"),field(f"{name}/Vc")
        assert (v>0).all()
        u,p,tracer=field(f"{name}/U","vector"),field(f"{name}/p"),field(f"{name}/s")
        angle=45*np.sin(2*np.pi*t)
        theta=np.deg2rad(angle)
        cos,sin=np.cos(theta),np.sin(theta)
        matrix=np.array([[cos,-sin,0],[sin,cos,0],[0,0,1]])
        ball_error=float(np.max(abs(c[zones["ball"]]-c0[zones["ball"]]@matrix.T)))
        pipe_error=float(np.max(abs(c[zones["pipes"]]-c0[zones["pipes"]])))
        phi=native(f"{name}/phi")
        flows={}
        maximum_wall=0.
        for patch,count in counts.items():
            body=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",phi,re.S)[1]
            values=read(body,"value",count,"scalar")
            assert values.shape==(count,) and np.isfinite(values).all()
            if patch in PORTS:
                flows[patch]=float(values.sum())
                assert abs(flows[patch]-reports[patch][step])<1e-10
            else:
                maximum_wall=max(maximum_wall,float(np.max(abs(values))))
        speed=np.linalg.norm(u,axis=1)
        rows.append({"time_s":float(t),"step":step,"relative_ball_rotation_deg":float(angle),
                     "ball_pose_max_abs_error_m":ball_error,"stationary_pipe_max_abs_error_m":pipe_error,
                     "native_fluid_volume_m3":float(v.sum()),"max_abs_wall_flux_m3_s":maximum_wall,
                     "inlet_flow_m3_s":flows["inlet"],"lower_outlet_flow_m3_s":flows["lowerOutlet"],
                     "upper_outlet_flow_m3_s":flows["upperOutlet"],
                     "mass_imbalance_relative_to_peak":float(abs(sum(flows.values()))/peak),
                     "minimum_kinematic_pressure_m2_s2":float(p.min()),"max_speed_m_s":float(speed.max()),
                     "tracer_min":float(tracer.min()),"tracer_max":float(tracer.max())})
        states[name]=(c,speed)
    assert rows[-1]["step"]==len(times)
    max_lower=max(rows,key=lambda r:r["lower_outlet_flow_m3_s"])
    max_upper=max(rows,key=lambda r:r["upper_outlet_flow_m3_s"])
    def dominant(active,other):
        return any(r[active]>=1e-4 and r[active]>=2*abs(r[other]) for r in rows)
    gates={"ball_pose":max(r["ball_pose_max_abs_error_m"] for r in rows)<=1e-7,
           "static_pipe":max(r["stationary_pipe_max_abs_error_m"] for r in rows)<=1e-8,
           "sealed_walls":max(r["max_abs_wall_flux_m3_s"] for r in rows)<=1e-12,
           "saved_external_flux_balance":max(r["mass_imbalance_relative_to_peak"] for r in rows)<=.005,
           "lower_outlet_dominates":dominant("lower_outlet_flow_m3_s","upper_outlet_flow_m3_s"),
           "upper_outlet_dominates":dominant("upper_outlet_flow_m3_s","lower_outlet_flow_m3_s")}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"steps":len(times),"saved_states":len(rows),
            "peak_absolute_native_port_flow_m3_s":peak,"lower_outlet_peak_save":max_lower,
            "upper_outlet_peak_save":max_upper,"samples":rows,"native_sha256":sources}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    for filename,records in (
        ("saved-valve-states.csv",rows),
        ("native-port-flow.csv",[
            {"step":i,"time_s":float(t),**{f"{patch}_m3_s":float(report[i]) for patch,report in reports.items()}}
            for i,t in enumerate(np.r_[0,times])])):
        with (out/filename).open("x",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    fig,axs=plt.subplots(2,2,figsize=(10.8,6.6),layout="constrained",gridspec_kw={"height_ratios":[1.25,1]})
    max_speed=max(r["max_speed_m_s"] for r in rows)
    for ax,key,row in zip(axs[0],("upper","lower"),(max_upper,max_lower)):
        t=row["time_s"]
        c,speed=states[f"{t:g}"]
        pick=abs(c[:,2])<.025
        im=ax.scatter(c[pick,0],c[pick,1],c=speed[pick],s=4,linewidths=0,cmap="viridis",vmin=0,vmax=max_speed)
        ax.set(title=f"{key.title()} outlet peak save at {t:g} s",xlabel="x (m)",ylabel="y (m)",aspect="equal")
        ax.set(xlim=(c0[:,0].min()-.04,c0[:,0].max()+.04),ylim=(c0[:,1].min()-.04,c0[:,1].max()+.04))
        ax.set_facecolor("#e6e9ee")
    fig.colorbar(im,ax=axs[0],label="Speed (m/s)",shrink=.8)
    ax=axs[1,0]
    for patch,name,sign in (("inlet","Inlet", -1),("upperOutlet","Upper outlet",1),("lowerOutlet","Lower outlet",1)):
        ax.plot(np.r_[0,times],sign*1000*reports[patch],label=name)
    ax.axhline(0,color="#95a1ad",lw=.7)
    ax.set(title="Native port flow",xlabel="Time (s)",ylabel="Volumetric flow (L/s)")
    ax.legend(frameon=False)
    ax=axs[1,1]
    ax.plot(np.r_[0,TIMES],[0]+[r["relative_ball_rotation_deg"] for r in rows],"o-",color="#2a607e",ms=3)
    ax.set(title="Ball rotation from its initial position",xlabel="Time (s)",ylabel="Rotation about z (°)")
    for ax in axs[1]:
        ax.grid(alpha=.2)
    fig.suptitle("Flow switching through an oscillating ball valve",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"ball-valve.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        parser.add_argument(name,type=Path)
    a=parser.parse_args()
    raise SystemExit(0 if audit(a.case,a.reader,a.output) else 2)
