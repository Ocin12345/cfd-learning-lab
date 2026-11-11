"""Check the native high-lift mesh, solver records, forces and pressure fields."""
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
from matplotlib.collections import PolyCollection
import numpy as np

NUMBER=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
CELLS=63957
PARTS=("slat","mainElement","flap")
CAP=1600
EQS=("Ux","Uy","p","p_corrected","nuTilda")
AX=np.deg2rad(5.5)
DRAG=np.array([np.cos(AX),np.sin(AX),0])
LIFT=np.array([-np.sin(AX),np.cos(AX),0])


def evaluate(case,reader,out,stop):
    read=runpy.run_path(str(reader))["foam_values"]
    hashes={}

    def raw(name):
        data=(case/name).read_bytes()
        hashes[name]=hashlib.sha256(data).hexdigest()
        return data

    def native(name):
        return raw(name).decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        sha,name=line.split("  ",1)
        raw(name)
        assert hashes[name]==sha,name
    for name in ("source-input-sha256.txt","inputs-verified.log","centres.log","volumes.log"):
        native(name)
    mesh=native("check-mesh.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*63957\s*$",mesh,re.M)
    meta=json.loads(native("mesh-source.json"))
    assert meta["cells"]==CELLS and meta["source_points"]==64540
    assert meta["source_wall_edges"]=={"slat":266,"mainElement":254,"flap":90}
    speed,span,chord,nu=(meta[key] for key in ("velocity_m_s","extrusion_span_m","stowed_chord_m",
                                              "kinematic_viscosity_m2_s"))
    assert span==.01 and chord==.4572 and 58.58<speed<58.60
    assert abs(speed*chord/nu-1.71e6)<1e-6
    assert meta["grid_sha256"]=="dc3e70dbd11fa3e38a604705b2aa39779f99ebf3da081fd5932a75b4dcedca09"
    q=.5*speed**2
    force_ref=1.225*q*span*chord

    log=native("solver.log")
    if stop==CAP:
        assert log.rstrip().endswith("End") and "FOAM FATAL" not in log
    assert "SpalartAllmaras" in log and "14-7b05503f98a8" in log
    times=list(re.finditer(r"^Time = (\d+)s?\s*$",log,re.M))
    assert [int(t[1]) for t in times[:stop]]==list(range(1,stop+1))
    if stop==CAP:
        assert len(times)==CAP
    residuals=[]
    native_bounds=[]
    for i,t in enumerate(times[:stop]):
        end=times[i+1].start() if i+1<len(times) else len(log)
        segment=log[t.end():end]
        values=re.findall(r"Solving for (Ux|Uy|p|nuTilda), Initial residual = ("+NUMBER+r")",segment)
        assert tuple(key for key,_ in values)==("Ux","Uy","p","p","nuTilda"),i+1
        residuals.append([float(v) for _,v in values])
        if "bounding nuTilda" in segment:
            native_bounds.append(i+1)
    residuals=np.array(residuals)
    assert np.isfinite(residuals).all() and (residuals>=0).all()
    saves=sorted(int(d.name) for d in case.iterdir() if d.is_dir() and d.name.isdigit() and 0<int(d.name)<=stop)
    assert saves==list(range(100,stop+1,100))
    for step in saves:
        clock=native(f"{step}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",clock,re.M)[1])==step
        assert float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])==step

    def field(name,kind="scalar"):
        values=read(native(name),"internalField",CELLS,kind)
        assert values.shape==((CELLS,3) if kind=="vector" else (CELLS,))
        assert np.isfinite(values).all()
        return values

    def body(text,name):
        return re.search(r"\b"+name+r"\s*\{([^{}]*)\}",text,re.S)[1]

    centres,volume=field("0/C","vector"),field("0/Vc")
    assert (volume>0).all() and np.max(abs(centres[:,2]))<1e-12
    prefix="constant/polyMesh"
    pts=np.array(re.findall(r"^\(("+NUMBER+r") ("+NUMBER+r") ("+NUMBER+r")\)$",
                            native(prefix+"/points"),re.M),dtype=float)
    faces=np.array(re.findall(r"^4\((\d+) (\d+) (\d+) (\d+)\)$",native(prefix+"/faces"),re.M),dtype=int)
    own=np.array(re.findall(r"^(\d+)$",native(prefix+"/owner"),re.M),dtype=int)
    assert len(faces)==own[0]==len(own)-1 and len(pts)==129080
    assert np.array_equal(pts[:64540,:2],pts[64540:,:2])
    assert np.array_equal(np.unique(pts[:,2]),np.array([-.005,.005]))
    own=own[1:]
    assert own.min()==0 and own.max()==CELLS-1
    b=native(prefix+"/boundary")
    patches={}
    for name,count in {**meta["source_wall_edges"],"farfield":560,"frontBack":2*CELLS}.items():
        pb=body(b,name)
        start=int(re.search(r"startFace\s+(\d+);",pb)[1])
        assert int(re.search(r"nFaces\s+(\d+);",pb)[1])==count
        patches[name]=(start,count)
    geometries={}
    for name in PARTS:
        start,count=patches[name]
        polys=pts[faces[start:start+count]]
        sf=.5*np.cross(polys,np.roll(polys,-1,axis=1)).sum(axis=1)
        area=np.linalg.norm(sf,axis=1)
        xy=polys.mean(axis=1)
        owners=own[start:start+count]
        gap=abs(np.sum((centres[owners]-xy)*sf,axis=1))/area
        assert (area>0).all() and (gap>0).all()
        geometries[name]=(xy,sf,owners,gap)

    native_reports={}
    for name in ("allElements",*(p+"Force" for p in PARTS)):
        data=native(f"postProcessing/{name}/0/forces.dat")
        assert "forces(pressure viscous)" in data
        rows=[re.findall(NUMBER,line) for line in data.splitlines() if line.strip() and not line.startswith("#")]
        values=np.array(rows[:stop+1],dtype=float)
        if stop==CAP:
            assert len(rows)==CAP+1
        assert values.shape==(stop+1,13) and np.array_equal(values[:,0],np.arange(stop+1))
        assert np.isfinite(values).all()
        native_reports[name]=values
    native_all=native_reports["allElements"]
    summed=sum(native_reports[name+"Force"][:,1:7] for name in PARTS)
    assert np.max(abs(summed-native_all[:,1:7]))<1e-8

    snapshots=[]
    for step in saves:
        u,p,nt,nut=(field(f"{step}/{f}", "vector" if f=="U" else "scalar") for f in ("U","p","nuTilda","nut"))
        assert (nt>=0).all() and (nut>=0).all()
        assert np.max(abs(u[:,2]))<=1e-12
        flux=native(f"{step}/phi")
        assert re.search(r"dimensions\s+\[0 3 -1 0 0 0 0\];",flux)
        ff=read(body(flux,"farfield"),"value",560,"scalar")
        assert ff.min()<0<ff.max()
        boundary=abs(float(ff.sum()))/-float(ff[ff<0].sum())
        assert "type            empty;" in body(flux,"frontBack")
        sealed=max(float(np.max(abs(read(body(flux,n),"value",patches[n][1],"scalar")))) for n in PARTS)
        values={}
        ptxt=native(f"{step}/p")
        for name,(_,sf,owners,gap) in geometries.items():
            assert re.search(r"type\s+zeroGradient;",body(ptxt,name))
            pf=p[owners]
            derived=1.225*np.sum(pf[:,None]*sf,axis=0)
            report=native_reports[name+"Force"][step,1:4]
            discrepancy=float(np.linalg.norm(derived-report)/max(np.linalg.norm(report),.1))
            assert discrepancy<1e-6,(name,step,discrepancy)
            values[name]={"native_pressure_load_n":report.tolist(),"relative_pressure_readback_error":discrepancy}
        snapshots.append({"iteration":step,"relative_farfield_net_imbalance":boundary,"sealed_face_flux_m3_s":sealed,
                          "max_speed_m_s":float(np.linalg.norm(u,axis=1).max()),
                          "minimum_nuTilda_m2_s":float(nt.min()),"max_nut_m2_s":float(nut.max()),
                          "native_pressure_loads":values})

    total=native_all[:,1:4]+native_all[:,4:7]
    drag,lift=total@DRAG,total@LIFT
    ft=np.column_stack((drag,lift))
    average=ft[-100:].mean(axis=0)
    variation=np.ptp(ft[-100:],axis=0)
    allowed=np.maximum(.02*abs(average),.002*force_ref)
    elements={}
    for name in PARTS:
        own_force=native_reports[name+"Force"]
        vec=own_force[-100:,1:4]+own_force[-100:,4:7]
        elements[name]={"mean_lift_coefficient":float(np.mean(vec@LIFT)/force_ref),
                        "mean_drag_coefficient":float(np.mean(vec@DRAG)/force_ref),
                        "median_cell_centre_to_wall_m":float(np.median(geometries[name][3]))}
    gates={"terminal_initial_residuals":bool(residuals[-1].max()<=1e-3),
           "terminal_farfield_mass":snapshots[-1]["relative_farfield_net_imbalance"]<=1e-4,
           "impermeability":snapshots[-1]["sealed_face_flux_m3_s"]<=1e-10,
           "lift_and_drag_range":bool(average[0]>0 and 1.6<average[1]/force_ref<3.5),
           "terminal_force_stability":bool((variation<=allowed).all())}
    summary={"pass":stop==CAP and all(gates.values()),"gates":gates,"iterations":stop,
             "cells":CELLS,"saved_states":len(saves),"mesh_source":meta,
             "reference_force_n":force_ref,"mean_lift_coefficient":float(average[1]/force_ref),
             "mean_drag_coefficient":float(average[0]/force_ref),"terminal_force_range_n":variation.tolist(),
             "terminal_force_range_limits_n":allowed.tolist(),"terminal_initial_residuals":dict(zip(EQS,residuals[-1])),
             "iterations_reporting_nuTilda_bounding":len(native_bounds),
             "last_iteration_reporting_nuTilda_bounding":native_bounds[-1] if native_bounds else None,
             "element_coefficients":elements,"saved_checks":snapshots,"source_sha256":hashes}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k not in ("source_sha256","saved_checks")},indent=2))
    print("Terminal snapshot",snapshots[-1])
    if stop!=CAP:
        return True
    with (out/"coefficients.csv").open("x",newline="") as f:
        w=csv.writer(f)
        w.writerow(("iteration","drag_n","lift_n","drag_coefficient","lift_coefficient",*EQS))
        w.writerows((i,x,y,x/force_ref,y/force_ref,*residuals[i-1]) for i,(x,y) in enumerate(ft[1:],1))
    with (out/"surface-pressure.csv").open("x",newline="") as f:
        w=csv.writer(f)
        w.writerow(("element","x_over_c","y_over_c","pressure_coefficient"))
        for name,(xyz,_,ids,_) in geometries.items():
            w.writerows((name,x,y,cp) for x,y,cp in zip(xyz[:,0]/chord,xyz[:,1]/chord,p[ids]/q))

    start,count=patches["frontBack"]
    side=faces[start:start+count]
    low=side.max(axis=1)<64540
    assert low.sum()==CELLS and np.array_equal(np.sort(own[start:start+count][low]),np.arange(CELLS))
    ids=own[start:start+count][low]
    a=centres[ids,:2]/chord
    plot=(a[:,0]>-.28)&(a[:,0]<1.43)&(a[:,1]>-.46)&(a[:,1]<.51)
    display=pts[side[low][plot],:2]/chord
    us=np.linalg.norm(u[ids[plot]],axis=1)/speed
    fig=plt.figure(figsize=(11.6,7),layout="constrained")
    fig.suptitle("30P30N high-lift aircraft section",fontsize=18,fontweight="bold")
    grid=fig.add_gridspec(2,2,height_ratios=[1.32,1.])
    ax=fig.add_subplot(grid[0,:])
    map_plot=PolyCollection(display,array=us,cmap="viridis",edgecolor="none",antialiased=False,rasterized=True)
    map_plot.set_clim(0,2.5)
    ax.add_collection(map_plot)
    for name,(xyz,_,_,_) in geometries.items():
        ax.scatter(xyz[:,0]/chord,xyz[:,1]/chord,s=1.1,color="#253143")
    ax.set(xlabel="Axial position x/c",ylabel="Vertical position y/c",title="Solved velocity at 5.5° angle of attack",
           xlim=(-.27,1.4),ylim=(-.43,.48),aspect="equal")
    cax=ax.inset_axes([1.035,.10,.027,.8])
    fig.colorbar(map_plot,cax=cax,label="Speed / freestream speed")
    ax=fig.add_subplot(grid[1,0])
    for name,(xyz,_,owners,_) in geometries.items():
        ax.scatter(xyz[:,0]/chord,p[owners]/q,s=9,label={"mainElement":"Main element"}.get(name,name.title()))
    ax.invert_yaxis()
    ax.set(xlabel="Axial position x/c",ylabel="Surface pressure coefficient",title="Pressure on each element")
    ax.legend(frameon=False,fontsize=8)
    ax.grid(alpha=.18)
    ax=fig.add_subplot(grid[1,1])
    steps=np.arange(200,stop+1)
    ax.plot(steps,lift[steps]/force_ref,color="#176b88",label="Lift")
    ax.set(xlabel="Steady iteration",ylabel="Lift coefficient",title="Native force history from iteration 200")
    twin=ax.twinx()
    twin.plot(steps,drag[steps]/force_ref,color="#bc5737",label="Drag")
    twin.set_ylabel("Drag coefficient",color="#bc5737")
    ax.axvspan(stop-99,stop,color="#60778d",alpha=.13,zorder=0)
    ax.legend(ax.lines+twin.lines,["Lift","Drag"],frameon=False,fontsize=8,loc="best")
    ax.grid(alpha=.18)
    fig.savefig(out/"high-lift-section.png",dpi=180)
    fig.savefig(out/"high-lift-section.pdf")
    plt.close(fig)
    return summary["pass"]


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        cli.add_argument(name,type=Path)
    cli.add_argument("--stop",type=int,choices=(100,1600),default=1600)
    args=cli.parse_args()
    raise SystemExit(0 if evaluate(args.case,args.reader,args.output,args.stop) else 2)
