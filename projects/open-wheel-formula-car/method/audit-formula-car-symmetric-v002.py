"""Audit steady SST runs and render the solved vehicle surfaces and flow."""
import argparse
import csv
import hashlib
import json
import re
import runpy
import struct
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

NUMBER=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
PARTS=("chassis","frontWheels","rearWheels","frontMain","frontFlap","rearWing")
PATCHES=("inlet","outlet","road","farfield","centrePlane",*PARTS)
REF_FORCE=.5*1.225*30**2
REF_P=.5*30**2
EQS=("Ux","Uy","Uz","p","k","omega")


def audit(case,reader,iteration_cap):
    read=runpy.run_path(str(reader))["foam_values"]
    hashes={}

    def raw(name):
        data=(case/name).read_bytes()
        hashes[name]=hashlib.sha256(data).hexdigest()
        return data

    def native(name):
        return raw(name).decode("ascii")

    for manifest in ("pre-mesh-sha256.txt","input-sha256.txt"):
        for line in native(manifest).splitlines():
            sha,name=line.split("  ",1)
            raw(name)
            assert hashes[name]==sha,name
    for name in ("block-mesh.log","snappy.log","generator-sha256.txt","original-inputs-verified.log",
                 "inputs-verified.log","centres.log","volumes.log"):
        native(name)
    for part,expected in zip(PARTS,(1,2,2,1,1,1)):
        text=native(f"surface-{part}.log")
        assert "Surface is closed." in text
        assert re.search(rf"Number of unconnected parts\s*:\s*{expected}",text)
    solid_bounds={}
    for part in PARTS:
        stl=raw(f"constant/geometry/{part}.stl")
        count=struct.unpack_from("<I",stl,80)[0]
        records=list(struct.iter_unpack("<12fH",stl[84:]))
        assert len(records)==count and count>0
        xyz=np.array([row[3:12] for row in records]).reshape((-1,3))
        solid_bounds[part]={"triangles":count,"min_xyz_m":xyz.min(axis=0).tolist(),
                            "max_xyz_m":xyz.max(axis=0).tolist()}
    assert solid_bounds["frontMain"]["max_xyz_m"][0] < solid_bounds["frontFlap"]["min_xyz_m"][0]
    assert solid_bounds["frontFlap"]["max_xyz_m"][0] < solid_bounds["chassis"]["min_xyz_m"][0]
    assert solid_bounds["chassis"]["max_xyz_m"][1] < .79
    for part in PARTS:
        assert solid_bounds[part]["min_xyz_m"][2] > 0
    mesh=native("check-mesh.log")
    assert "Mesh OK." in mesh
    cells=int(re.search(r"^\s*cells:\s*(\d+)",mesh,re.M)[1])
    assert 50000<=cells<=130000
    log=native("solver.log")
    assert log.rstrip().endswith("End") and "kOmegaSST" in log and "14-7b05503f98a8" in log
    times=list(re.finditer(r"^Time = (\d+)s?\s*$",log,re.M))
    assert [int(t[1]) for t in times]==list(range(1,iteration_cap+1))
    residual=[]
    for i,t in enumerate(times):
        stop=times[i+1].start() if i+1<len(times) else len(log)
        vals=re.findall(r"Solving for (Ux|Uy|Uz|p|k|omega), Initial residual = ("+NUMBER+r")",log[t.end():stop])
        assert len(vals)==6 and {key for key,_ in vals}==set(EQS),(i,vals)
        residual.append([float(dict(vals)[key]) for key in EQS])
    residual=np.array(residual)
    assert np.isfinite(residual).all()

    def field(name,kind="scalar"):
        a=read(native(name),"internalField",cells,kind)
        assert a.shape==((cells,3) if kind=="vector" else (cells,))
        assert np.isfinite(a).all(),name
        return a

    def body(text,name):
        start=re.search(r"\b"+name+r"\s*\{",text).end()
        depth=1
        for marker in re.finditer(r"[{}]",text[start:]):
            depth += 1 if marker[0]=="{" else -1
            if depth==0:
                return text[start:start+marker.start()]
        raise ValueError(f"Unclosed native patch {name}")

    centres,vol=field("0/C","vector"),field("0/Vc")
    assert (vol>0).all()
    faces=[np.fromstring(v,sep=" ",dtype=int) for v in re.findall(
        r"^\d+\(([0-9 ]+)\)$",native("constant/polyMesh/faces"),re.M)]
    owner=np.array(re.findall(r"^(\d+)$",native("constant/polyMesh/owner"),re.M),dtype=int)
    assert len(faces)==len(owner)-1==owner[0]
    owner=owner[1:]
    points=np.array(re.findall(r"^\(("+NUMBER+r") ("+NUMBER+r") ("+NUMBER+r")\)$",
                               native("constant/polyMesh/points"),re.M),dtype=float)
    boundary=native("constant/polyMesh/boundary")
    assert set(re.findall(r"^\s*(\w+)\s*\n\s*\{\s*type\s+",boundary,re.M))==set(PATCHES)
    patches={}
    for name in PATCHES:
        patch=body(boundary,name)
        count,start=(int(re.search(rf"\b{key}\s+(\d+);",patch)[1]) for key in ("nFaces","startFace"))
        assert count>0
        patches[name]=(start,count)
    assert points[:,1].min()>=-1e-12
    start,n=patches["centrePlane"]
    plane_vertices=np.unique(np.concatenate(faces[start:start+n]))
    assert np.max(abs(points[plane_vertices,1]))<1e-12
    assert re.search(r"type\s+symmetryPlane;",body(boundary,"centrePlane"))
    native_force={}
    for function in ("vehicleForce","chassisForce","frontWingForce","flapForce","rearWingForce"):
        reports=native(f"postProcessing/{function}/0/forces.dat")
        assert "forces(pressure viscous)" in reports
        vals=np.array([re.findall(NUMBER,line) for line in reports.splitlines()
                       if line.strip() and not line.startswith("#")],dtype=float)
        assert vals.shape==(iteration_cap+1,13) and np.array_equal(vals[:,0],np.arange(iteration_cap+1))
        assert np.isfinite(vals).all()
        native_force[function]=vals
    save_dirs=sorted(int(p.name) for p in case.iterdir() if p.is_dir() and p.name.isdigit() and int(p.name)>0)
    assert save_dirs==list(range(50,iteration_cap+1,50))
    samples=[]
    for step in save_dirs:
        t=native(f"{step}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",t,re.M)[1])==step
        assert float(re.search(r"^value\s+([^;]+);",t,re.M)[1])==step
        u,p=field(f"{step}/U","vector"),field(f"{step}/p")
        k,omega,nut=(field(f"{step}/{q}") for q in ("k","omega","nut"))
        assert (k>0).all() and (omega>0).all() and (nut>=0).all()
        flux=native(f"{step}/phi")
        assert re.search(r"dimensions\s+\[0 3 -1 0 0 0 0\];",flux)
        rates={}
        wall=0.
        for name,(_,count) in patches.items():
            f=read(body(flux,name),"value",count,"scalar")
            rates[name]=float(f.sum())
            if name not in ("inlet","outlet"):
                wall=max(wall,float(np.max(abs(f))))
        assert rates["inlet"]<0<rates["outlet"]
        samples.append({"iteration":step,"minimum_pressure_m2_s2":float(p.min()),
                        "maximum_velocity_m_s":float(np.linalg.norm(u,axis=1).max()),
                        "minimum_tke_m2_s2":float(k.min()),"minimum_omega_s_inverse":float(omega.min()),
                        "maximum_nut_m2_s":float(nut.max()),"patch_volume_flux_m3_s":rates,
                        "relative_boundary_imbalance":abs(sum(rates.values()))/abs(rates["inlet"]),
                        "maximum_impermeable_face_flux_m3_s":wall})

    geom={}
    for name in PARTS:
        start,count=patches[name]
        poly=[points[ids] for ids in faces[start:start+count]]
        sf=np.array([.5*np.cross(v,np.roll(v,-1,axis=0)).sum(axis=0) for v in poly])
        norm=np.linalg.norm(sf,axis=1)
        c=np.array([v.mean(axis=0) for v in poly])
        own=owner[start:start+count]
        gap=abs(np.sum((c-centres[own])*sf,axis=1))/norm
        assert (norm>0).all() and (gap>0).all()
        geom[name]=(poly,c,sf,gap,p[own]/REF_P)

    u_text=native(f"{iteration_cap}/U")
    road=read(body(u_text,"road"),"value",patches["road"][1],"vector")
    assert np.max(abs(road-np.array([30.,0.,0.])))<=1e-12
    wheel_motion={}
    for name in ("frontWheels","rearWheels"):
        condition=body(u_text,name)
        assert "rotatingWallVelocity;" in condition
        _,cc,sf,_,_=geom[name]
        native_u=read(condition,"value",patches[name][1],"vector")
        normals=sf/np.linalg.norm(sf,axis=1)[:,None]
        lower,upper=cc[:,2]<.12,cc[:,2]>.64
        assert lower.sum()>=10 and upper.sum()>=10
        ux_lo,ux_hi=float(native_u[lower,0].mean()),float(native_u[upper,0].mean())
        normal_speed=float(np.max(abs(np.sum(native_u*normals,axis=1))))
        assert 18<ux_lo<36 and -36<ux_hi<-18 and normal_speed<1e-8
        wheel_motion[name]={"lower_mean_ux_m_s":ux_lo,"upper_mean_ux_m_s":ux_hi,
                            "maximum_normal_velocity_m_s":normal_speed}
    _,cc,sf,_,cp=geom["chassis"]
    area=np.linalg.norm(sf,axis=1)
    pressure_zones={}
    for name,xmin,xmax in (("flat_throat",2.6,3.17),("rear_floor",3.30,4.37)):
        mask=(cc[:,0]>=xmin)&(cc[:,0]<=xmax)&(abs(cc[:,1])<=.52)&(cc[:,2]<.47)&(sf[:,2]>.75*area)
        assert mask.sum()>=20
        pressure_zones[name]={"face_count":int(mask.sum()),"area_m2":float(area[mask].sum()),
                              "area_weighted_pressure_coefficient":float(np.average(cp[mask],weights=area[mask]))}
    patch_pressure={name:1.225*np.sum(a[2]*a[4][:,None]*REF_P,axis=0) for name,a in geom.items()}
    estimated_pressure=sum(patch_pressure.values())
    actual=native_force["vehicleForce"]
    half_total=actual[:,1:4]+actual[:,4:7]
    total=half_total*np.array([2.,0.,2.])
    chassis=native_force["chassisForce"]
    chassis_total=(chassis[:,1:4]+chassis[:,4:7])*np.array([2.,0.,2.])
    pressure_error=float(np.linalg.norm(estimated_pressure-actual[-1,1:4])
                         /max(np.linalg.norm(actual[-1,1:4]),1.))
    pressure_error_body=float(np.linalg.norm(patch_pressure["chassis"]-chassis[-1,1:4])
                              /max(np.linalg.norm(chassis[-1,1:4]),1.))
    assert pressure_error<1e-5 and pressure_error_body<1e-5,(pressure_error,pressure_error_body)
    groups={"chassisForce":("chassis",),"frontWingForce":("frontMain","frontFlap"),
            "flapForce":("frontFlap",),"rearWingForce":("rearWing",)}
    component={}
    for name,parts in groups.items():
        rows=native_force[name]
        p_force=sum(patch_pressure[part] for part in parts)
        discrepancy=float(np.linalg.norm(p_force-rows[-1,1:4])/max(1.,np.linalg.norm(p_force)))
        assert discrepancy<1e-5,(name,discrepancy)
        vector=rows[-50:,1:4]+rows[-50:,4:7]
        component[name]={"mean_drag_n":float(2*vector[:,0].mean()),
                         "mean_downforce_n":float(-2*vector[:,2].mean()),
                         "pressure_force_relative_readback":discrepancy}
    component["wheelsCombined"]={
        "mean_drag_n":float(total[-50:,0].mean()-sum(component[name]["mean_drag_n"]
                           for name in ("chassisForce","frontWingForce","rearWingForce"))),
        "mean_downforce_n":float(-total[-50:,2].mean()-sum(component[name]["mean_downforce_n"]
                                for name in ("chassisForce","frontWingForce","rearWingForce")))}
    tail=total[-50:,(0,2)]
    average=tail.mean(axis=0)
    span=np.ptp(tail,axis=0)
    force_tol=np.maximum(.03*abs(average),.0015*REF_FORCE)
    gates={
        "terminal_initial_residuals":bool(residual[-1].max()<=1e-3),
        "terminal_mass":samples[-1]["relative_boundary_imbalance"]<=1e-4,
        "impermeability":samples[-1]["maximum_impermeable_face_flux_m3_s"]<=1e-10,
        "positive_drag":bool(average[0]>0),
        "terminal_force_stability":bool((span<=force_tol).all()),
    }
    summary={"pass":all(gates.values()),"gates":gates,"case":case.name,"cells":cells,
             "iterations":iteration_cap,"saved_states":len(samples),"reference_force_n":REF_FORCE,
             "last_50_mean_drag_n":float(average[0]),"last_50_mean_downforce_n":float(-average[1]),
             "last_50_mean_drag_coefficient":float(average[0]/REF_FORCE),
             "last_50_mean_downforce_coefficient":float(-average[1]/REF_FORCE),
             "terminal_full_equivalent_force_n":total[-1].tolist(),
             "terminal_native_half_force_n":half_total[-1].tolist(),
             "zero_yaw_centre_plane_symmetry":True,
             "terminal_full_equivalent_chassis_force_n":chassis_total[-1].tolist(),
             "terminal_native_half_pressure_force_n":actual[-1,1:4].tolist(),
             "terminal_native_half_viscous_force_n":actual[-1,4:7].tolist(),
             "mean_pitch_moment_n_m_about_x2p17_z0p18":float(2*np.mean(actual[-50:,8]+actual[-50:,11])),
             "last_50_force_range_n":span.tolist(),"last_50_force_range_limits_n":force_tol.tolist(),
             "terminal_initial_residuals":dict(zip(EQS,residual[-1])),
             "rotating_wall_readback":wheel_motion,"underfloor_pressure_zones":pressure_zones,
             "component_forces":component,"source_geometry_bounds":solid_bounds,
             "surface_pressure_force_relative_readback":pressure_error,
             "surface_pressure_force_chassis_relative_readback":pressure_error_body,
             "chassis_cell_centre_wall_distance_quantiles_m":np.quantile(geom["chassis"][3],[0,.5,1]).tolist(),
             "samples":samples,"source_sha256":hashes}
    return summary,geom,centres,u,p,total,residual


def output_pair(results,paths,out,second_angle):
    names=["8° front flap",f"{second_angle}° front flap"]
    out.mkdir()
    summary={"pass":all(item[0]["pass"] for item in results),"flap_angles_deg":[8,second_angle],
             "cases":[r[0] for r in results]}
    def source_manifest(result):
        text=(result/"pre-mesh-sha256.txt").read_text()
        return {path:digest for digest,path in (row.split("  ",1) for row in text.splitlines())}
    base,other=map(source_manifest,paths)
    assert base.keys()==other.keys()
    changed={key for key in base if base[key]!=other[key]}
    assert changed=={"constant/geometry/frontFlap.stl"},changed
    summary["differing_input_files"]=sorted(changed)
    q0,q1=(item[0] for item in results)
    summary["changes"]={
        "drag_n":q1["last_50_mean_drag_n"]-q0["last_50_mean_drag_n"],
        "downforce_n":q1["last_50_mean_downforce_n"]-q0["last_50_mean_downforce_n"],
        "front_wing_downforce_n":q1["component_forces"]["frontWingForce"]["mean_downforce_n"]-
                                  q0["component_forces"]["frontWingForce"]["mean_downforce_n"],
    }
    (out/"checks.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    with (out/"component-loads.csv").open("x",newline="") as f:
        w=csv.writer(f)
        w.writerow(("case","component","mean_drag_n","mean_downforce_n"))
        for name,(q,*_) in zip(names,results):
            for part,data in q["component_forces"].items():
                w.writerow((name,part,data["mean_drag_n"],data["mean_downforce_n"]))
    with (out/"forces.csv").open("x",newline="") as f:
        w=csv.writer(f)
        w.writerow(("case","iteration","full_equivalent_drag_n","full_equivalent_downforce_n","drag_coefficient","downforce_coefficient"))
        for name,(_,_,_,_,_,total,_) in zip(names,results):
            w.writerows((name,i,x,-z,x/REF_FORCE,-z/REF_FORCE) for i,(x,y,z) in enumerate(total))
    with (out/"initial-residuals.csv").open("x",newline="") as f:
        w=csv.writer(f)
        w.writerow(("case","iteration",*EQS))
        for name,(_,_,_,_,_,_,res) in zip(names,results):
            w.writerows((name,i,*row) for i,row in enumerate(res,1))
    fig=plt.figure(figsize=(12,8.2),layout="constrained")
    fig.suptitle("Front-flap angle on an open-wheel formula car",fontsize=18,fontweight="bold")
    grid=fig.add_gridspec(3,2,height_ratios=[1.5,1,1])
    cp_norm=TwoSlopeNorm(vmin=-1.8,vcenter=0.,vmax=1.2)
    panels=[]
    sections=[]
    for col,(name,(check,geom,c,u,p,_,_)) in enumerate(zip(names,results)):
        ax=fig.add_subplot(grid[0,col],projection="3d",computed_zorder=True)
        panels.append(ax)
        for polygons,_,_,_,cp in geom.values():
            reflected=[poly*np.array([1.,-1.,1.]) for poly in polygons]
            v=Poly3DCollection(polygons+reflected,cmap="coolwarm",norm=cp_norm,edgecolors="none",
                               linewidth=0,antialiased=False)
            v.set_array(np.concatenate([cp,cp]))
            ax.add_collection3d(v)
        ax.set(xlim=(-.08,4.65),ylim=(-1.12,1.12),zlim=(0,1.24),title=name)
        ax.set_box_aspect((4.73,2.24,1.24),zoom=1.22)
        ax.view_init(elev=25,azim=143)
        ax.set_axis_off()
        ax=fig.add_subplot(grid[1,col])
        sections.append(ax)
        view=(abs(c[:,1])<.052)&(c[:,0]>-.5)&(c[:,0]<4.95)&(c[:,2]<1.45)
        flow=ax.scatter(c[view,0],c[view,2],c=np.linalg.norm(u[view],axis=1)/30,
                        s=8,cmap="viridis",vmin=0,vmax=1.5,edgecolors="none",rasterized=True)
        ax.axvline(3.2,color=".4",ls=":",lw=.7)
        ax.set(xlim=(-.5,4.95),ylim=(0,1.45),xlabel="Streamwise position (m)",
               ylabel="Height (m)",title="Near-centreplane cell samples")
    fig.colorbar(v,ax=panels,label="Surface pressure coefficient (mirrored)",shrink=.8,extend="both")
    fig.colorbar(flow,ax=sections,label="Speed / inlet speed",shrink=.85,extend="max")
    for col,(label,sign,k) in enumerate((("Drag",1,0),("Vertical force, positive down",-1,2))):
        ax=fig.add_subplot(grid[2,col])
        for name,result,color in zip(names,results,("#5b6475","#146d89")):
            force=result[5]
            x=np.arange(len(force))
            ax.plot(x[20:],sign*force[20:,k],label=name,color=color,lw=1.5)
        ax.set(xlabel="Steady iteration",ylabel="Drag (N)" if col==0 else "Vertical force (N, positive down)")
        ax.legend(frameon=False,fontsize=8)
        ax.grid(alpha=.19)
    fig.savefig(out/"formula-car-comparison.png",dpi=180)
    fig.savefig(out/"formula-car-comparison.pdf")
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,7.4),layout="constrained",
                          gridspec_kw={"height_ratios":[1.,1.1]})
    fig.suptitle("Front-wing flow and vehicle load distribution",fontsize=17,fontweight="bold")
    for ax,title,(_,geom,c,u,p,_,_) in zip(axes[0],names,results):
        region=(abs(c[:,1])<.023)&(c[:,0]>-.13)&(c[:,0]<1.04)&(c[:,2]<.57)
        image=ax.scatter(c[region,0],c[region,2],c=np.linalg.norm(u[region],axis=1)/30,
                         s=8,cmap="viridis",vmin=0,vmax=1.6,edgecolors="none",rasterized=True)
        for name in ("frontMain","frontFlap"):
            _,xyz,_,_,_=geom[name]
            chosen=abs(xyz[:,1])<.024
            ax.scatter(xyz[chosen,0],xyz[chosen,2],s=5,c="#273346")
        ax.set(title=title,xlabel="Streamwise position (m)",ylabel="Height above road (m)",
               xlim=(-.13,1.04),ylim=(0,.57))
    fig.colorbar(image,ax=list(axes[0]),label="Speed / inlet speed",extend="max")
    groups=("frontWingForce","chassisForce","rearWingForce","wheelsCombined")
    labels=("Front wing","Chassis","Rear wing","Wheels")
    x=np.arange(len(groups))
    for ax,metric,label in zip(axes[1],("mean_downforce_n","mean_drag_n"),("Vertical force (N, positive down)","Drag (N)")):
        for i,(title,result,color) in enumerate(zip(names,results,("#5b6475","#146d89"))):
            ax.bar(x+(i-.5)*.36,[result[0]["component_forces"][part][metric] for part in groups],
                   .36,label=title,color=color)
        ax.set(xticks=x,xticklabels=labels,ylabel=label)
        ax.axhline(0,color=".45",lw=.65)
        ax.grid(axis="y",alpha=.18)
        ax.legend(frameon=False,fontsize=8)
    fig.savefig(out/"front-wing-detail.png",dpi=180)
    fig.savefig(out/"front-wing-detail.pdf")
    plt.close(fig)
    return summary["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline",type=Path)
    parser.add_argument("comparison",type=Path)
    parser.add_argument("--second-angle",type=int,choices=(12,18),required=True)
    parser.add_argument("reader",type=Path)
    parser.add_argument("output",type=Path)
    parser.add_argument("--iterations",type=int,default=600)
    args=parser.parse_args()
    cases=[audit(c,args.reader,args.iterations) for c in (args.baseline,args.comparison)]
    for case, *_ in cases:
        print(json.dumps({k:v for k,v in case.items() if k not in ("samples","source_sha256")},indent=2))
    raise SystemExit(0 if output_pair(cases,(args.baseline,args.comparison),args.output,args.second_angle) else 2)
