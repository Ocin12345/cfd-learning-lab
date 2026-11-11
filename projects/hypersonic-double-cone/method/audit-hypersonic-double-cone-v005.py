"""Audit the native double-cone solution, mesh, wall loads and local steps."""
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
from matplotlib.colors import LogNorm
import numpy as np

NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
CELLS = 32385
R_N2 = 8314.46261815324/28.01348
P_INF, T_INF, U_INF = 22.74, 138.9, 2713.
RHO_INF = P_INF/(R_N2*T_INF)
TURN, BASE = .09208087, .153683
PATCH_COUNTS = {"cone":255,"outlet":127,"freestream":255,"centreLeft":0,
                "wedge1":32385,"wedge2":32385}


def audit(case, reader, output, limit):
    read = runpy.run_path(str(reader))["foam_values"]
    sources = {}

    def native(name):
        raw = (case/name).read_bytes()
        sources[name] = hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        sha,name = line.split("  ",1)
        native(name)
        assert sources[name] == sha, name
    for line in native("tutorial-input-sha256.txt").splitlines():
        sha,name = line.split("  ",1)
        if name != "system/controlDict":
            assert sources[name] == sha, name
    for name in ("block-mesh.log","convert.log","collapse.log","inputs-verified.log",
                 "centres.log","volumes.log","wedge-projection-sha256.txt","native-grid-points.txt"):
        native(name)
    projection = json.loads(native("wedge-projection.json"))
    native("preprojection-wedge-points")
    assert projection["preprojection_sha256"] == sources["preprojection-wedge-points"]
    assert projection["points_sha256"] == sources["constant/polyMesh/points"]
    assert projection["off_axis_shared_nodes"] == 2 and projection["max_axis_node_change_m"] < 5e-6
    checkmesh = native("check-mesh.log")
    assert "Mesh OK." in checkmesh and re.search(r"^\s*cells:\s*32385\s*$",checkmesh,re.M)
    log = native("solver.log")
    assert "14-7b05503f98a8" in log and "Using LTS" in log
    assert "localEuler;" in native("system/fvSchemes")
    if limit == 10000:
        assert log.rstrip().endswith("End") and "FOAM FATAL" not in log
    found = list(re.finditer(r"^Time = ("+NUMBER+r")s?\s*$",log,re.M))
    assert len(found)>=limit
    if limit == 10000:
        assert len(found)==limit
    clocks = np.array([float(t[1]) for t in found[:limit]])
    assert np.max(abs(clocks-np.arange(1,limit+1)*1e-7))<1e-12
    scales = np.array(re.findall(r"^Flow time scale min/max = ("+NUMBER+r"), ("+NUMBER+r")",log,re.M)[:limit],dtype=float)
    assert scales.shape==(limit,2) and np.isfinite(scales).all() and (scales>0).all()
    assert (scales[:,0]<=scales[:,1]).all() and scales[:,1].max()<=2e-8+1e-16
    residuals=[]
    for i,match in enumerate(found[:limit]):
        end = found[i+1].start() if i+1 < len(found) else len(log)
        values = re.findall(r"Solving for (rho|Ux|Uy|Uz|e), Initial residual = ("+NUMBER+r"), Final residual = ("+NUMBER+r")",
                            log[match.end():end])
        assert tuple(a for a,_,_ in values)==("rho","Ux","Uy","Uz","e"),i+1
        residuals.append([float(x) for _,a,b in values for x in (a,b)])
    assert np.isfinite(residuals).all()

    prefix="constant/polyMesh"
    points=np.array(re.findall(r"^\(("+NUMBER+r") ("+NUMBER+r") ("+NUMBER+r")\)$",
                               native(prefix+"/points"),re.M),dtype=float)
    faces=[np.fromstring(row,sep=" ",dtype=int) for row in
           re.findall(r"^\d+\(([\d ]+)\)$",native(prefix+"/faces"),re.M)]
    owner=np.array(re.findall(r"^([0-9]+)$",native(prefix+"/owner"),re.M),dtype=int)
    assert len(points)==65408 and len(faces)==owner[0]==len(owner)-1==129795
    owner=owner[1:]
    assert owner.min()==0 and owner.max()==CELLS-1
    boundary=native(prefix+"/boundary")
    patches={}
    for name, count in PATCH_COUNTS.items():
        body=re.search(r"\b"+name+r"\s*\{([^}]+)\}",boundary,re.S)[1]
        assert int(re.search(r"nFaces\s+(\d+);",body)[1])==count
        start=int(re.search(r"startFace\s+(\d+);",body)[1])
        patches[name]=(start,count)
    cone, cone_count=patches["cone"]
    wall_faces=[points[ids] for ids in faces[cone:cone+cone_count]]
    x_wall=np.array([v[:,0].mean() for v in wall_faces])
    radius=np.array([v[:,1].mean() for v in wall_faces])
    areas=[]
    for v in wall_faces:
        # OpenFOAM face-area vector is the sum of oriented edge cross-products.
        areas.append(.5*np.linalg.norm(np.cross(v,np.roll(v,-1,axis=0)).sum(axis=0)))
    areas=np.array(areas)
    assert (areas>0).all() and np.all(np.diff(x_wall)>0)
    section = np.where(x_wall<TURN,"25-degree cone",np.where(x_wall<BASE,"55-degree cone","base"))
    vol=read(native("0/Vc"),"internalField",CELLS,"scalar")
    assert (vol>0).all()

    def field(t,name,kind="scalar"):
        return read(native(f"{t}/{name}"),"internalField",CELLS,kind)

    def patch(text,name,n,kind="scalar"):
        body=re.search(r"\b"+name+r"\s*\{([^}]+)\}",text,re.S)[1]
        return read(body,"value",n,kind)

    source_wall=native("postProcessing/wallHeatFlux/0/wallHeatFlux.dat")
    native_heat=np.array([[float(s) for s in line.replace("cone","").split()]
                          for line in source_wall.splitlines() if line.strip() and not line.startswith("#")])
    assert native_heat.shape[1]==5
    positive_saves=sorted((p for p in case.iterdir() if p.is_dir() and re.fullmatch(NUMBER,p.name)
                           and 0<float(p.name)<=limit*1e-7+1e-12),key=lambda x:float(x.name))
    assert len(positive_saves)==limit//500
    if limit==10000:
        assert len(native_heat)==21
    profiles=[]
    records=[]
    for index,directory in enumerate(positive_saves,1):
        t=directory.name
        step=index*500
        assert abs(float(t)-step*1e-7)<1e-13
        time=native(f"{t}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",time,re.M)[1])==step
        assert abs(float(re.search(r"^value\s+("+NUMBER+r");",time,re.M)[1])-step*1e-7)<1e-12
        assert float(re.search(r"^deltaT\s+("+NUMBER+r");",time,re.M)[1])==1e-7
        ptext=native(f"{t}/p")
        p=read(ptext,"internalField",CELLS,"scalar")
        T,rho,U=field(t,"T"),field(t,"rho"),field(t,"U","vector")
        local_step=1/field(t,"rDeltaT")
        speed=np.linalg.norm(U,axis=1)
        eos=float(np.max(abs(rho*R_N2*T/p-1)))
        assert U.shape==(CELLS,3) and (p>0).all() and (rho>0).all()
        assert (T>=100).all() and (T<=10000).all()
        assert (local_step>0).all() and local_step.max()<=2e-8+1e-16
        assert re.search(r"\bcone\s*\{\s*type\s+zeroGradient;\s*\}",ptext)
        pw=p[owner[cone:cone+cone_count]]
        qw=patch(native(f"{t}/wallHeatFlux"),"cone",cone_count)
        assert (pw>0).all()
        native_report=native_heat[index]
        assert abs(native_report[0]-step*1e-7)<1e-12
        assert np.allclose([qw.min(),qw.max(),qw@areas],
                           native_report[1:4],rtol=1e-8,atol=1e-10)
        phitext=native(f"{t}/phi")
        assert re.search(r"dimensions\s+\[1 0 -1 0 0 0 0\];",phitext)
        flux={name:float(patch(phitext,name,count).sum()) for name,(_,count) in patches.items() if count}
        assert flux["freestream"]<0<flux["outlet"]
        relative_flow=abs(sum(flux.values()))/abs(flux["freestream"])
        profiles.append((pw,qw))
        records.append({"step":step,"clock_label":t,"min_p_pa":float(p.min()),"max_p_pa":float(p.max()),
                        "min_t_k":float(T.min()),"max_t_k":float(T.max()),"max_velocity_m_s":float(speed.max()),
                        "relative_perfect_gas_eos_error":eos,"sector_gas_inventory_kg":float(rho@vol),
                        "boundary_mass_flow_kg_s":flux,"relative_steady_boundary_flux_mismatch":relative_flow,
                        "native_wall_heat_rate_w":float(qw@areas),"net_heat_into_wall_w":float(-qw@areas),
                        "wall_heating_peak_w_m2":float(-qw.min()),
                        "local_step_min_s":float(local_step.min()),"local_step_max_s":float(local_step.max())})

    pressure, heat=profiles[-1]
    means={name:float(np.average(pressure[section==name],weights=areas[section==name]))
           for name in ("25-degree cone","55-degree cone")}
    def profile_change(a,b):
        return float(np.sqrt(np.average((a-b)**2,weights=areas)/np.average(a**2,weights=areas)))
    p_change=profile_change(pressure,profiles[-2][0]) if len(profiles)>1 else None
    q_change=profile_change(heat,profiles[-2][1]) if len(profiles)>1 else None
    gates={"both_cones_compress":min(means.values())>P_INF,
           "positive_wall_heating":float(-heat.min())>0,
           "fast_and_slow_regions":bool(np.mean(speed>1500)>=.01 and np.mean(speed<300)>=.01),
           "terminal_pressure_profile_stable":p_change is not None and p_change<=.05}
    summary={"pass":limit==10000 and all(gates.values()),"gates":gates,"cells":CELLS,"iterations":limit,
             "positive_time_label_states":len(records),"wall_faces":cone_count,"starting_state":"installed developed solution",
             "local_time_stepping":True,"wall_section_mean_p_pa":means,
             "wall_peak_pressure_pa":float(pressure.max()),
             "wall_peak_pressure_x_m":float(x_wall[np.argmax(pressure)]),
             "wall_peak_heating_x_m":float(x_wall[np.argmin(heat)]),
             "relative_terminal_wall_pressure_profile_change":p_change,"relative_terminal_heat_flux_profile_change":q_change,
             "terminal_initial_residuals":dict(zip(("rho","Ux","Uy","Uz","e"),residuals[-1][::2])),
             "mesh_projection":projection,"snapshots":records,"native_sha256":sources}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in summary.items() if k not in ("snapshots","native_sha256")},indent=2))
    print("Terminal snapshot",records[-1])
    if limit!=10000:
        return True

    with (output/"wall-profiles.csv").open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("iteration","section","axial_x_m","radius_m","sector_face_area_m2",
                         "static_pressure_pa","native_wall_heat_flux_w_m2","wall_heating_w_m2"))
        for entry,(p_face,q_face) in zip(records,profiles):
            writer.writerows((entry["step"],*row) for row in zip(section,x_wall,radius,areas,p_face,q_face,-q_face))
    with (output/"local-step-history.csv").open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("iteration","minimum_local_step_s","maximum_local_step_s","rho_initial_residual",
                         "Ux_initial_residual","Uy_initial_residual","Uz_initial_residual","e_initial_residual"))
        writer.writerows((i,*scale,*res[::2]) for i,(scale,res) in enumerate(zip(scales,residuals),1))

    start,count=patches["wedge2"]
    display_owner=owner[start:start+count]
    assert np.array_equal(np.sort(display_owner),np.arange(CELLS))
    display_faces=[points[ids,:2]*1000 for ids in faces[start:start+count]]
    fig=plt.figure(figsize=(11,6.7),layout="constrained")
    grid=fig.add_gridspec(2,2,height_ratios=[1.4,1])
    fig.suptitle("Hypersonic flow over a 25°/55° double cone",fontweight="bold",fontsize=17)
    for col,(data,title,cmap,minimum,maximum) in enumerate((
            (rho/RHO_INF,"Density / freestream density","magma",1,18),
            (p/P_INF,"Static pressure / freestream pressure","viridis",1,400))):
        ax=fig.add_subplot(grid[0,col])
        poly=PolyCollection(display_faces,array=data[display_owner],cmap=cmap,edgecolor="none",antialiased=False,rasterized=True)
        if col==0:
            poly.set_clim(minimum,maximum)
        else:
            poly.set_norm(LogNorm(vmin=minimum,vmax=maximum))
        ax.add_collection(poly)
        ax.plot(x_wall*1000,radius*1000,color="#374151",lw=.9)
        ax.set(xlabel="Axial position (mm)",ylabel="Radial position (mm)",title=title,
               xlim=(-5,196),ylim=(0,192),aspect="equal")
        fig.colorbar(poly,ax=ax,shrink=.83,extend="both")
    for col,(prev,now,title,color) in enumerate(((profiles[-2][0]/1000,pressure/1000,"Wall static pressure (kPa)","#245a89"),
                                                (-profiles[-2][1]/1000,-heat/1000,"Heating into the wall (kW/m²)","#aa4728"))):
        ax=fig.add_subplot(grid[1,col])
        ax.plot(x_wall*1000,now,label="Iteration 10,000",color=color,lw=1.6)
        ax.plot(x_wall*1000,prev,label="Iteration 9,500",color="#64748b",ls="--",lw=1)
        for mark in (TURN,BASE):
            ax.axvline(mark*1000,color=".6",lw=.7)
        ax.set(xlabel="Axial position (mm)",ylabel=title)
        ax.legend(frameon=False,fontsize=8)
        ax.grid(alpha=.2)
    fig.savefig(output/"double-cone.png",dpi=180)
    fig.savefig(output/"double-cone.pdf")
    plt.close(fig)
    return summary["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        parser.add_argument(name,type=Path)
    parser.add_argument("--limit",type=int,choices=(500,10000),default=10000)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output,args.limit) else 2)
