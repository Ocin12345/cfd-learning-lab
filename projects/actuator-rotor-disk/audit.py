"""Audit native rotor-disk fields and compare upstream and downstream speeds."""
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

NUMBER = r"([-+\d.eE]+)"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, reader_path, out):
    read = runpy.run_path(str(reader_path))["foam_values"]
    sources = []

    def native(name):
        sources.append(name)
        return (case/name).read_text("ascii")

    for line in native("input-sha256.txt").splitlines():
        expected, name = line.split("  ",1)
        assert sha(case/name)==expected,name
    sources.extend(["pre-mesh-input-sha256.txt","pre-mesh-inputs-verified.log","snappy.log"])
    mesh, log = native("check-mesh.log"), native("solver.log")
    assert "Mesh OK." in mesh and log.rstrip().endswith("\nEnd")
    cells = int(re.search(r"^\s*cells:\s*(\d+)\s*$",mesh,re.M)[1])
    assert 0 < cells <= 250000
    assert "Selecting rotorDisk with name disk" in log and "Selecting RAS turbulence model kOmegaSST" in log
    blocks = re.split(r"(?m)^Time = (\d+)s\s*$", log)
    iterations = np.array(blocks[1::2],int)
    endpoint = len(iterations)
    assert 0 < endpoint <= 1000 and np.array_equal(iterations,np.arange(1,endpoint+1))
    convergence = re.search(r"SIMPLE solution converged in (\d+) iterations",log)
    assert (convergence and int(convergence[1])==endpoint) or endpoint==1000

    def field(name, kind="scalar", key="internalField"):
        value = read(native(name), key, cells, kind)
        assert np.all(np.isfinite(value)), name
        return value

    c, v = field("0/C","vector"), field("0/Vc")
    assert np.all(v>0) and c.shape==(cells,3)
    x, y, z = c.T
    r = np.hypot(x,z)
    zones = native("constant/polyMesh/cellZones")
    record = re.search(r"rotatingZone\s*\{\s*cellLabels\s+List<label>\s+(\d+)\s*\(([^)]*)\)",zones)
    zone = np.fromstring(record[2],dtype=int,sep=" ")
    assert len(zone)==int(record[1])>0 and len(np.unique(zone))==len(zone)
    assert np.all((0<=zone)&(zone<cells))
    mask = np.zeros(cells,dtype=bool)
    mask[zone] = True
    up = (r<.25)&(y>=-.75)&(y<=-.25)
    down = (r<.25)&(y>=.25)&(y<=.75)
    assert min(up.sum(),down.sum())>=10
    model = native("constant/fvModels")
    rho_ref = float(re.search(r"\brhoRef\s+"+NUMBER+r";",model)[1])
    rpm = float(re.search(r"\bomega\s+"+NUMBER+r"\s*\[rpm\];",model)[1])
    assert rho_ref==1000 and rpm==1000

    saves = sorted(int(p.name) for p in case.iterdir() if p.is_dir() and re.fullmatch(r"[1-9]\d*",p.name))
    assert saves==sorted(set(range(50,endpoint+1,50))|{endpoint})
    previous = None
    rows = []
    for i in saves:
        clock = native(f"{i}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",clock,re.M)[1])==i
        assert float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])==i
        u = field(f"{i}/U","vector")
        field(f"{i}/p")
        k, omega = field(f"{i}/k"), field(f"{i}/omega")
        block = blocks[2*i]
        reports = {kind:float(re.findall(r"Effective "+kind+r"\s*=\s*"+NUMBER,block)[-1])
                   for kind in ("lift","drag","power")}
        assert np.all(np.isfinite(list(reports.values())))
        continuity = float(re.findall(r"sum local = "+NUMBER,block)[-1])
        change = None if previous is None else float(np.sqrt(np.average(np.sum((u-previous)**2,axis=1),weights=v)))
        force_paths = list((case/str(i)).glob("disk?rotorForce"))
        force_check = None
        if force_paths:
            assert len(force_paths)==1
            force = field(force_paths[0].relative_to(case).as_posix(),"vector","value")
            assert np.all(force[~mask]==0)
            lift = rho_ref*float(force[:,1]@v)
            torque = rho_ref*float((z*force[:,0]-x*force[:,2])@v)
            power = torque*rpm*2*np.pi/60
            assert np.isclose(lift,reports["lift"],rtol=2e-5,atol=1e-7)
            assert np.isclose(power,reports["power"],rtol=2e-5,atol=1e-7)
            force_check = {"lift_n":lift,"torque_nm":torque,"power_w":power}
        rows.append({"iteration":i,"upstream_axial_m_s":float(np.average(u[up,1],weights=v[up])),
                     "downstream_axial_m_s":float(np.average(u[down,1],weights=v[down])),
                     "terminal_local_continuity":continuity,
                     "velocity_change_rms_m_s":change,"min_k_m2_s2":float(k.min()),
                     "min_omega_s_inv":float(omega.min()),"reported_effective_lift_n":reports["lift"],
                     "reported_effective_drag_n":reports["drag"],"reported_effective_power_w":reports["power"],
                     "force_field_lift_n":None if force_check is None else force_check["lift_n"],
                     "force_field_torque_nm":None if force_check is None else force_check["torque_nm"],
                     "force_field_power_w":None if force_check is None else force_check["power_w"]})
        previous = u
    last=rows[-1]
    delta=last["downstream_axial_m_s"]-last["upstream_axial_m_s"]
    gates={"turbulence_fields":all(row["min_k_m2_s2"]>=0 and row["min_omega_s_inv"]>0 for row in rows),
           "continuity":last["terminal_local_continuity"]<1e-4,
           "stable_end":bool(convergence) or (change is not None and change<.02),
           "axial_disturbance":abs(delta)>.01}
    result={"pass":all(gates.values()),"gates":gates,"cells":cells,"rotor_zone_cells":len(zone),
            "upstream_sample_cells":int(up.sum()),"downstream_sample_cells":int(down.sum()),
            "iterations":endpoint,"native_convergence":bool(convergence),"saved_states":len(rows),
            "reference_density_kg_m3":rho_ref,"rpm":rpm,"downstream_minus_upstream_m_s":delta,
            "terminal_torque_from_reported_power_nm":last["reported_effective_power_w"]/(rpm*2*np.pi/60),
            "samples":rows,"native_sha256":{name:sha(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"rotor-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes=plt.subplots(1,2,figsize=(12,4.8),layout="constrained")
    slab=(abs(z)<.02)&(abs(x)<.38)&(abs(y)<.15)
    assert slab.sum()>100
    scatter=axes[0].scatter(y[slab],x[slab],c=u[slab,1],cmap="viridis",s=9,linewidths=0)
    axes[0].plot([0,0],[-.25,.25],color="#1d2939",lw=2,label="Nominal rotor span")
    axes[0].set(xlim=(-.15,.15),ylim=(-.38,.38),xlabel="Axial position y (m)",
                ylabel="Transverse x (m)",title="Near-disk cell centres within 20 mm of z = 0")
    axes[0].legend(loc="upper left",frameon=False)
    fig.colorbar(scatter,ax=axes[0],label="Axial velocity (m/s)",shrink=.75)
    edges=np.linspace(-1.5,1.5,31)
    profiles=[]
    for select,label in ((r<.25,"Rotor radius"),((r>=.35)&(r<.65),"Outer annulus")):
        means=[]
        for low,high in zip(edges[:-1],edges[1:]):
            use=select&(y>=low)&(y<high)
            means.append(float(np.average(u[use,1],weights=v[use])) if use.any() else np.nan)
        axes[1].plot((edges[1:]+edges[:-1])/2,means,"o-",ms=3,label=label)
        profiles.append(means)
    axes[1].axvline(0,color="#444",ls=":",lw=1)
    axes[1].set(xlabel="Axial position y (m)",ylabel="Mean axial velocity (m/s)",title="Volume averages in axial bins")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=.2)
    np.savetxt(out/"rotor-axial-profile.csv",np.column_stack([(edges[1:]+edges[:-1])/2,*profiles]),
               header="axial_position_m,rotor_radius_mean_m_s,outer_annulus_mean_m_s",delimiter=",",comments="")
    fig.suptitle("Axial flow around an actuator rotor",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"actuator-rotor-disk.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case","reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)