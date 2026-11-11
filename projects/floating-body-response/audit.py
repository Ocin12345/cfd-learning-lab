"""Check native saved water inventory, mesh and two-degree-of-freedom motion."""
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

TIMES = np.arange(1, 21) / 10
NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
INITIAL = np.array([.5, .45, .1])


def table(document, columns):
    result = []
    for line in document.splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            tokens = line.replace("(", " ").replace(")", " ").split()
            assert len(tokens) == columns and all(re.fullmatch(NUMBER, x) for x in tokens), line
            result.append(list(map(float, tokens)))
    result = np.array(result)
    assert result.ndim == 2 and result.shape[1] == columns and np.isfinite(result).all()
    return result


def audit(case, reader, out):
    read = runpy.run_path(str(reader))["foam_values"]
    sources = {}

    def native(name):
        raw = (case/name).read_bytes()
        sources[name] = hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest, name = line.split("  ", 1)
        raw = (case/name).read_bytes()
        sources[name] = hashlib.sha256(raw).hexdigest()
        assert sources[name] == digest, name
    for name in ("pre-mesh-input-sha256.txt", "pre-mesh-inputs-verified.log",
                 "inputs-verified.log", "subset.log", "centres.log", "volumes.log"):
        native(name)
    mesh, log = native("check-mesh.log"), native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*11640\s*$", mesh, re.M)
    assert log.rstrip().endswith("\nEnd")
    times = np.array(re.findall(r"^Time = ("+NUMBER+r")s$", log, re.M), float)
    assert len(times) and np.isfinite(times).all() and (np.diff(times)>0).all() and times[-1]==2
    refined = [(int(a),int(b)) for a,b in re.findall(r"^Refined from (\d+) to (\d+) cells\.",log,re.M)]
    assert refined and all(0<a<b<=200000 for a,b in refined)

    motion_text = native("postProcessing/rigidBodyState/0/floatingObject.dat")
    assert re.search(r"^# Angle Units\s*:\s*\[deg\]", motion_text, re.M)
    motion = table(motion_text, 13)
    forces = table(native("postProcessing/rigidBodyForces/0/rigidBodyForces.dat"),16)
    assert np.array_equal(motion[:,0], np.r_[0,times]) and np.array_equal(forces[:,0],motion[:,0])
    assert np.max(abs(motion[:,1:4]-forces[:,1:4])) < 1e-9
    np.testing.assert_allclose(motion[0,1:4],INITIAL,atol=1e-12,rtol=0)
    assert np.max(abs(motion[0,4:])) < 1e-12
    assert np.max(abs(motion[:,[1,3]]-INITIAL[[0,2]])) < 1e-9
    assert np.max(abs(motion[:,[4,6,7,9,10,12]])) < 1e-9

    def field(name,n,kind="scalar"):
        result=read(native(name),"internalField",n,kind)
        assert result.shape == ((n,3) if kind=="vector" else (n,)), name
        assert np.isfinite(result).all(),name
        return result

    n0=11640
    c0,v0=field("0/C",n0,"vector"),field("0/Vc",n0)
    assert np.all(v0>0) and abs(v0.sum()-.97)<=1e-6
    assert np.allclose(v0,1/12000,atol=1e-14,rtol=0)
    init=native("0/alpha.water")
    assert re.search(r"box\s+\(-100 -100 -100\)\s+\(100 100 0\.5368\);",init)
    assert re.search(r"box\s+\(0\.7 0\.8 -100\)\s+\(100 100 0\.651\);",init)
    assert re.search(r"defaultValue\s+0;",init) and len(re.findall(r"\bvalue\s+1;",init))==2
    a0=((c0[:,2]<.5368)|((c0[:,0]>.7)&(c0[:,1]>.8)&(c0[:,2]<.651))).astype(float)
    inventory=float(a0@v0)
    np.testing.assert_allclose(inventory,.5153333333333333,atol=1e-10,rtol=0)
    actual={p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?",p.name)}
    assert actual=={"0",*(f"{t:g}" for t in TIMES)}

    rows=[]
    samples={"0":(c0,a0)}
    for t in TIMES:
        name=f"{t:g}"
        clock=native(f"{name}/uniform/time")
        step=int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert 1<=step<=len(times) and abs(times[step-1]-t)<1e-9
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-t)<1e-10
        vtext=native(f"{name}/Vc")
        n=int(re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(",vtext)[1])
        assert n0<=n<=200000
        v=read(vtext,"internalField",n,"scalar")
        assert v.shape==(n,) and np.isfinite(v).all() and np.all(v>0)
        c=field(f"{name}/C",n,"vector")
        a=field(f"{name}/alpha.water",n)
        u=field(f"{name}/U",n,"vector")
        field(f"{name}/p_rgh",n)
        if name in {"1","2"}:
            samples[name]=(c,a)
        rows.append({"time_s":float(t),"step":step,"cells":n,"alpha_min":float(a.min()),
                     "alpha_max":float(a.max()),"fluid_volume_m3":float(v.sum()),
                     "water_volume_m3":float(a@v),"relative_water_inventory_change":float(a@v/inventory-1),
                     "max_speed_m_s":float(np.linalg.norm(u,axis=1).max()),
                     "rotation_reference_y_m":float(motion[step,2]),
                     "rotation_reference_y_displacement_m":float(motion[step,2]-INITIAL[1]),
                     "native_y_orientation_deg":float(motion[step,5]),
                     "force_y_n":float(forces[step,5]+forces[step,8]),
                     "moment_y_nm":float(forces[step,11]+forces[step,14])})
    lateral=float(np.max(abs(motion[:,2]-INITIAL[1])))
    pitch=float(np.max(abs(motion[:,5])))
    gates={"fraction":all(-1e-6<=r["alpha_min"]<=r["alpha_max"]<=1+1e-6 for r in rows),
           "water_inventory":max(abs(r["relative_water_inventory_change"]) for r in rows)<.005,
           "fluid_volume":max(abs(r["fluid_volume_m3"]-.97) for r in rows)<=1e-6,
           "refinement":bool(refined and any(r["cells"]>n0 for r in rows)),
           "saved_motion":any(abs(r["rotation_reference_y_displacement_m"])>=.001
                              or abs(r["native_y_orientation_deg"])>=np.degrees(.01) for r in rows)}
    result={"pass":all(gates.values()),"gates":gates,"steps":len(times),"saved_states":len(rows),
            "initial_cells":n0,"refinement_events":len(refined),"initial_water_volume_m3":inventory,
            "maximum_absolute_y_displacement_m":lateral,"maximum_absolute_native_y_orientation_deg":pitch,
            "samples":rows,"native_sha256":sources}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    for filename,records in (
        ("saved-states.csv",rows),
        ("motion-history.csv",[
            {"step":i,"time_s":m[0],"rotation_reference_y_m":m[2],"native_y_orientation_deg":m[5],
             "native_y_velocity_m_s":m[8],"native_y_angular_velocity_deg_s":m[11],
             "force_y_n":f[5]+f[8],"moment_y_nm":f[11]+f[14]}
            for i,(m,f) in enumerate(zip(motion,forces))])):
        with (out/filename).open("x",newline="") as file:
            writer=csv.DictWriter(file,fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    fig=plt.figure(figsize=(12,6.8),layout="constrained")
    spec=fig.add_gridspec(2,3,height_ratios=[1.35,1])
    panes=[fig.add_subplot(spec[0,i]) for i in range(3)]
    for ax,(label,(c,a)) in zip(panes,samples.items()):
        mask=(c[:,1]>=.4)&(c[:,1]<.5)
        points=ax.scatter(c[mask,0],c[mask,2],c=a[mask],s=3,cmap="Blues",vmin=0,vmax=1,lw=0)
        ax.set(title=f"{label} s",xlabel="x (m)",ylabel="z (m)",aspect="equal",xlim=(0,1),ylim=(0,1))
        ax.set_facecolor("#e4e7eb")
    fig.colorbar(points,ax=panes,label="Water volume fraction",shrink=.78)
    for i,ys,title,ylabel,color in (
        (0,1000*(motion[:,2]-INITIAL[1]),"Lateral displacement","Displacement (mm)","#255b8a"),
        (1,motion[:,5],"Rotation about y","Native orientation (°)","#9b4a3e"),
        (2,100*np.array([r["relative_water_inventory_change"] for r in rows]),"Water inventory","Change from initial (%)","#398376")):
        ax=fig.add_subplot(spec[1,i])
        x=motion[:,0] if i<2 else TIMES
        ax.plot(x,ys,"-" if i<2 else "o-",ms=3,color=color,lw=1.5)
        ax.axhline(0,color="#98a2ae",lw=.8)
        ax.set(title=title,xlabel="Time (s)",ylabel=ylabel,xlim=(0,2))
        ax.grid(alpha=.2)
    fig.suptitle("Floating body response to a displaced free surface",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"floating-body-response.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("case","reader","output"):
        parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)
