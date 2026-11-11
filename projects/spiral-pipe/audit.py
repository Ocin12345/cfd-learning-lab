"""Check the native spiral pipe, patch reports and steady-flow endpoint."""
import argparse
import csv
import gzip
import hashlib
import json
import re
import runpy
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CELLS = 112448
Q = .04
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


def patch_geometry(text):
    blocks = re.search(r"POINTS\s+(\d+)\s+\w+\s+(.+?)\s+POLYGONS\s+(\d+)\s+(\d+)\s+",
                       text,re.S)
    points = np.fromstring(blocks[2],sep=" ").reshape(int(blocks[1]),3)
    raw = np.fromstring(text[blocks.end():].split("CELL_DATA")[0],sep=" ",dtype=int)
    assert len(raw)==int(blocks[4]) and np.isfinite(points).all()
    index = 0
    vectors = []
    for _ in range(int(blocks[3])):
        size = raw[index]
        vertices = points[raw[index+1:index+1+size]]
        vectors.append(.5*np.cross(vertices,np.roll(vertices,-1,axis=0)).sum(axis=0))
        index += size+1
    assert index == len(raw)
    vectors = np.array(vectors)
    area = np.linalg.norm(vectors,axis=1)
    total = vectors.sum(axis=0)
    assert (area>0).all() and np.linalg.norm(total)/area.sum()>.999
    return {"faces":len(vectors),"area_m2":float(area.sum()),
            "outward_unit_normal":(total/np.linalg.norm(total)).tolist(),
            "point_centre_m":points.mean(axis=0).tolist()}


def audit(case,reader,out):
    read=runpy.run_path(str(reader))["foam_values"]
    sources={}

    def raw(name):
        data=(case/name).read_bytes()
        sources[name]=hashlib.sha256(data).hexdigest()
        return data

    def native(name):
        return raw(name).decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest,name=line.split("  ",1)
        assert hashlib.sha256(raw(name)).hexdigest()==digest,name
    for name in ("inputs-verified.log","extrude.log","patch-export.log","centres.log","volumes.log","mesh-sha256.txt"):
        native(name)
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and "Number of regions: 1 (OK)." in mesh
    assert re.search(r"^\s*cells:\s*112448\s*$",mesh,re.M) and log.rstrip().endswith("\nEnd")
    steps=list(re.finditer(r"^Time = (\d+)s?$",log,re.M))
    end=int(steps[-1][1])
    assert 1<=end<=1000 and [int(m[1]) for m in steps]==list(range(1,end+1))
    residuals=[]
    for i,m in enumerate(steps):
        next_step=steps[i+1].start() if i+1<len(steps) else len(log)
        block=log[m.end():next_step]
        row={name:float(value) for name,value in re.findall(
            r"Solving for (Ux|Uy|Uz|p), Initial residual = ("+NUMBER+r")",block)}
        assert set(row)=={"Ux","Uy","Uz","p"} and all(np.isfinite(v) and v>=0 for v in row.values())
        continuity=float(re.search(r"time step continuity errors : sum local = ("+NUMBER+r")",block)[1])
        assert np.isfinite(continuity) and continuity>=0
        residuals.append({"iteration":i+1,**row,"continuity_sum_local":continuity})

    path_text=gzip.decompress(raw("constant/geometry/curve.obj.gz")).decode("ascii")
    path=np.array([[float(v) for v in line.split()[1:]] for line in path_text.splitlines() if line.startswith("v ")])
    assert path.ndim==2 and path.shape[1]==3 and np.isfinite(path).all()
    length=float(np.linalg.norm(np.diff(path,axis=0),axis=1).sum())
    geometry={name:patch_geometry(native(f"VTK/{name}/{name}_0.vtk")) for name in ("inlet","outlet")}
    for g in geometry.values():
        assert g["faces"]==448 and .194<g["area_m2"]<.196
    np.testing.assert_allclose(geometry["inlet"]["outward_unit_normal"],[0,0,-1],atol=1e-8,rtol=0)
    end_tangent=path[-1]-path[-2]
    assert np.dot(geometry["outlet"]["outward_unit_normal"],end_tangent)/np.linalg.norm(end_tangent)>.999
    np.testing.assert_allclose(geometry["outlet"]["point_centre_m"],path[-1],atol=1e-5,rtol=0)

    def report(function,patch):
        folder=f"{function}(patch={patch})" if function=="patchFlowRate" else f"{function}(patch={patch},fields=(p))"
        text=native(f"postProcessing/{folder}/0/surfaceFieldValue.dat")
        assert re.search(r"^# Selection\s*:\s*patch\("+patch+r"\)",text,re.M)
        assert re.search(r"^# Faces\s*:\s*448",text,re.M)
        area=float(re.search(r"^# Area\s*:\s*("+NUMBER+")",text,re.M)[1])
        assert abs(area-geometry[patch]["area_m2"])<2e-6
        field="sum(phi)" if function=="patchFlowRate" else "areaAverage(p)"
        assert field in text
        table=np.array([[float(x) for x in line.split()] for line in text.splitlines()
                        if line.strip() and not line.lstrip().startswith("#")])
        assert table.shape==(end+1,2) and np.isfinite(table).all()
        np.testing.assert_array_equal(table[:,0],np.arange(end+1))
        return table[:,1]

    inlet,outlet=(report("patchFlowRate",patch) for patch in ("inlet","outlet"))
    p_in,p_out=(report("patchAverage",patch) for patch in ("inlet","outlet"))
    rows=[{**row,"inlet_flow_m3_s":float(inlet[i]),"outlet_flow_m3_s":float(outlet[i]),
           "relative_flux_imbalance":float(abs(inlet[i]+outlet[i])/Q),
           "inlet_kinematic_pressure_m2_s2":float(p_in[i]),"outlet_kinematic_pressure_m2_s2":float(p_out[i])}
          for i,row in enumerate(residuals,1)]

    def field(name,kind="scalar"):
        a=read(native(name),"internalField",CELLS,kind)
        assert a.shape==((CELLS,3) if kind=="vector" else (CELLS,)) and np.isfinite(a).all(),name
        return a

    centres,vol=field("0/C","vector"),field("0/Vc")
    assert (vol>0).all() and abs(vol.sum()-9.23810405016)<1e-8
    saves=sorted(int(p.name) for p in case.iterdir() if p.is_dir() and p.name.isdigit() and int(p.name)>0)
    assert saves[-1]==end and saves==sorted(set(range(25,end+1,25))|{end})
    for step in saves:
        clock=native(f"{step}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",clock,re.M)[1])==step
        assert float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])==step
    u,p=field(f"{end}/U","vector"),field(f"{end}/p")
    previous=field(f"{saves[-2]}/U","vector") if len(saves)>1 else field("0/U","vector")
    u_rms=float(np.sqrt(np.average(np.sum((u-previous)**2,axis=1),weights=vol)))
    phi=native(f"{end}/phi")
    flux={}
    for patch,count in (("inlet",448),("outlet",448),("walls",8032)):
        body=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",phi,re.S)[1]
        data=read(body,"value",count,"scalar")
        assert data.shape==(count,) and np.isfinite(data).all()
        flux[patch]=float(data.sum())
        if patch=="walls":
            assert np.max(abs(data))<=1e-12
    assert abs(flux["inlet"]-inlet[-1])<1e-10 and abs(flux["outlet"]-outlet[-1])<1e-10
    dp=float(p_in[-1]-p_out[-1])
    gates={"native_convergence":f"PIMPLE solution converged in {end} iterations" in log,
           "terminal_residuals":max(residuals[-1][eq] for eq in ("Ux","Uy","Uz","p"))<=1e-5,
           "inlet_flow":abs(inlet[-1]+Q)<=.001*Q,
           "mass_balance":abs(inlet[-1]+outlet[-1])<=1e-4*Q,
           "positive_pressure_drop":dp>0,
           "velocity_change":len(saves)>1 and u_rms<=.001}
    speed=np.linalg.norm(u,axis=1)
    gates={name:bool(value) for name,value in gates.items()}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"iterations":end,"saved_states":len(saves),
            "terminal":rows[-1],"previous_saved_iteration":saves[-2] if len(saves)>1 else 0,
            "rms_velocity_change_m_s":u_rms,"terminal_kinematic_pressure_drop_m2_s2":dp,
            "terminal_min_pressure_m2_s2":float(p.min()),"terminal_max_speed_m_s":float(speed.max()),
            "path_polyline_length_m":length,"geometry":geometry,"native_sha256":sources}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"spiral-pipe-history.csv").open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig=plt.figure(figsize=(11.8,7.1),layout="constrained")
    layout=fig.add_gridspec(2,3,width_ratios=[1.15,1.1,1.25])
    samples=np.random.default_rng(0).choice(CELLS,9500,replace=False)
    axes=[]
    for i,(coord,label) in enumerate(((0,"x"),(1,"y"))):
        ax=fig.add_subplot(layout[:,i])
        axes.append(ax)
        pts=centres[samples]
        im=ax.scatter(pts[:,coord],pts[:,2],c=speed[samples],cmap="viridis",
                      vmin=0,vmax=float(speed.max()),s=2,linewidths=0,alpha=.78,rasterized=True)
        ax.plot(path[:,coord],path[:,2],color="#263e56",alpha=.4,lw=.6)
        ax.set(xlabel=f"{label} (m)",ylabel="z (m)",title=f"{label}–z projection")
        ax.set_aspect("equal")
        ax.grid(alpha=.22)
    fig.colorbar(im,ax=axes,location="bottom",shrink=.8,label="Speed (m/s)")
    a=fig.add_subplot(layout[0,2])
    i=np.arange(1,end+1)
    a.plot(i,p_in[1:]-p_out[1:],color="#276a67")
    a.set_yscale("symlog",linthresh=4)
    a.set(title="Inlet-to-outlet pressure drop",xlabel="Steady iteration",ylabel="Kinematic pressure (m²/s²)")
    b=fig.add_subplot(layout[1,2])
    for name in ("Ux","Uy","Uz","p"):
        b.semilogy(i,np.maximum([row[name] for row in rows],1e-20),label=name)
    b.axhline(1e-5,color="#666",ls="--",lw=.8)
    b.set(title="Initial residuals",xlabel="Steady iteration",ylabel="Scaled residual")
    b.legend(frameon=False,ncol=2)
    for ax in (a,b):
        ax.grid(alpha=.2)
    fig.suptitle("Laminar flow through a three-dimensional spiral pipe",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"spiral-pipe.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        parser.add_argument(name,type=Path)
    a=parser.parse_args()
    raise SystemExit(0 if audit(a.case,a.reader,a.output) else 2)
