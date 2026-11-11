"""Audit the native OpenFOAM Venturi calculation and plot its measured fields."""
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

CELLS = 57600
NUMBER = r"[-+0-9.eE]+"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, foam_reader, output):
    read_values = runpy.run_path(str(foam_reader))["foam_values"]
    def field(path, kind):
        return read_values((case/path).read_text(encoding="ascii"), "internalField", CELLS, kind)

    for line in (case/"input-sha256.txt").read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ",1)
        assert sha(case/name)==expected, name
    mesh = (case/"check-mesh.log").read_text(encoding="ascii")
    log = (case/"solver.log").read_text(encoding="ascii")
    assert re.search(r"^\s*cells:\s*57600\s*$",mesh,re.M) and "Mesh OK." in mesh
    assert "PIMPLE solution converged in 350 iterations" in log and re.search(r"^End\s*$",log,re.M)
    steps = np.array(re.findall(r"^Time = (\d+)s$", log, re.M),dtype=int)
    assert np.array_equal(steps,np.arange(1,351))
    residuals = re.findall(r"Solving for (Ux|Uy|Uz|p), Initial residual = ("+NUMBER+")",log.split("Time = 350s",1)[1])
    assert {name for name,_ in residuals}=={"Ux","Uy","Uz","p"} and len(residuals)==4
    assert max(float(value) for _,value in residuals)<=1e-5
    continuity = np.array(re.findall(r"time step continuity errors : sum local = ("+NUMBER+")",log),dtype=float)
    assert len(continuity)==350 and np.isfinite(continuity).all() and continuity[-1]<1e-5

    c = field("0/C","vector")
    volumes = field("0/Vc","scalar")
    u = field("350/U","vector")
    p = field("350/p","scalar")
    previous = field("300/U","vector")
    assert np.all(volumes>0) and c.shape==u.shape==(CELLS,3)
    u_rms = np.sqrt(np.average(np.sum((u-previous)**2,axis=1),weights=volumes))
    assert u_rms < .001

    phi = (case/"350/phi").read_text(encoding="ascii")
    flow = {}
    for patch in ("inlet","outlet"):
        match=re.search(r"(?ms)^\s*"+patch+r"\s*\{([^}]+)\}",phi)
        assert match, patch
        flux=read_values(match[1],"value",576,"scalar")
        flow[patch]=float(flux.sum())
    imbalance=abs(flow["inlet"]+flow["outlet"])/abs(flow["inlet"])
    assert flow["inlet"]<0<flow["outlet"] and imbalance<.001

    sections = []
    for name, target in (("upstream",.05),("throat",.25),("downstream",.575)):
        nearest = np.argmin(np.abs(c[:,0]-target))
        x = c[nearest,0]
        selected = abs(c[:,0]-x)<1e-8
        assert selected.sum()==576
        weights = volumes[selected]
        sections.append({"section":name,"sample_x_m":float(x),
                         "mean_axial_speed_m_s":float(np.average(u[selected,0],weights=weights)),
                         "mean_kinematic_pressure_m2_s2":float(np.average(p[selected],weights=weights))})
    upstream,throat,downstream=sections
    dp=upstream["mean_kinematic_pressure_m2_s2"]-throat["mean_kinematic_pressure_m2_s2"]
    ideal=.5*(throat["mean_axial_speed_m_s"]**2-upstream["mean_axial_speed_m_s"]**2)
    recovery=downstream["mean_kinematic_pressure_m2_s2"]-throat["mean_kinematic_pressure_m2_s2"]
    assert throat["mean_axial_speed_m_s"]>upstream["mean_axial_speed_m_s"] and dp>0 and recovery>0

    graphs = {}
    for name in ("B","D","F"):
        file=case/f"postProcessing/graph{name}/350/line.xy"
        points=np.loadtxt(file)
        assert points.shape[1]==4 and len(points)>10 and np.isfinite(points).all()
        graphs[name]=points
    bound=("input-sha256.txt","solver.log","check-mesh.log","0/C","0/Vc","300/U","350/U","350/p","350/phi",
           *(f"postProcessing/graph{name}/350/line.xy" for name in graphs))
    checks={"mesh_cells":CELLS,"iterations":len(steps),"stopping_rule":"Native PIMPLE residual control",
            "terminal_initial_residuals":{name:float(value) for name,value in residuals},
            "terminal_local_continuity":float(continuity[-1]),"relative_flux_imbalance":imbalance,
            "boundary_flux_m3_s":flow,"velocity_rms_change_300_to_350_m_s":float(u_rms),"sections":sections,
            "upstream_to_throat_pressure_drop_m2_s2":dp,"ideal_bernoulli_drop_m2_s2":ideal,
            "throat_to_downstream_pressure_recovery_m2_s2":recovery,
            "native_sha256":{name:sha(case/name) for name in bound}}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    with (output/"velocity-sections.csv").open("x",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(("section","y_m","u_x_m_s","u_y_m_s","u_z_m_s"))
        for section, values in graphs.items():
            writer.writerows((section,*row) for row in values)

    fig,(a,b)=plt.subplots(1,2,figsize=(11,4),layout="constrained")
    labels={"B":"Upstream section","D":"Throat","F":"Downstream section"}
    for name,values in graphs.items():
        a.plot(1000*values[:,0],values[:,1],label=labels[name])
    a.set(xlabel="Transverse position (mm)",ylabel="Axial velocity (m/s)",title="Native axial-velocity profiles")
    a.axhline(0,color="#8194a4",linewidth=.7)
    a.grid(alpha=.2)
    a.legend()
    sample_x=[row["sample_x_m"]*1000 for row in sections]
    b.plot(sample_x,[row["mean_kinematic_pressure_m2_s2"] for row in sections],"o--",color="#007fa3")
    b.set(xlabel="Axial position (mm)",ylabel="Mean kinematic pressure (m²/s²)",title="Pressure at three native sections")
    b.set_xticks(sample_x, [f"{position:.1f}" for position in sample_x])
    b.grid(alpha=.2)
    fig.suptitle("Flow through a circular Venturi meter | OpenFOAM",fontweight="bold")
    for extension in ("png","pdf"):
        fig.savefig(output/f"venturi-comparison.{extension}",dpi=180)
    plt.close(fig)
    print(json.dumps({key:value for key,value in checks.items() if key!="native_sha256"},indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case","foam_reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    audit(args.case,args.foam_reader,args.output)
