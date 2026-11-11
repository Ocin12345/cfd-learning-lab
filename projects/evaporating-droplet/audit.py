"""Check one-way evaporation using saved parcel fields and native cloud reports."""
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

CELLS=125
TIMES=np.arange(1,101)/1000
NUMBER=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
M0=5.23599e-10


def audit(case,reader,out):
    read=runpy.run_path(str(reader))["foam_values"]
    sources={}

    def native(name):
        raw=(case/name).read_bytes()
        sources[name]=hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest,name=line.split("  ",1)
        native(name)
        assert sources[name]==digest,name
    for name in ("pre-mesh-input-sha256.txt","pre-mesh-inputs-verified.log","inputs-verified.log","mesh.log"):
        native(name)
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*125\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd") and "reactingCloud" in log and "liquidEvaporation" in log
    assert re.search(r"coupled\s+false;",native("constant/cloudProperties"))
    steps=list(re.finditer(r"^Time = ("+NUMBER+r")s?$",log,re.M))
    times=np.array([float(step[1]) for step in steps])
    assert np.isfinite(times).all() and np.all(np.diff(times)>0) and times[-1]==.1
    all_reports=[]
    for i,step in enumerate(steps):
        next_step=steps[i+1].start() if i+1<len(steps) else len(log)
        body=log[step.end():next_step]
        def report(key):
            matches=re.findall(re.escape(key)+r"\s*=\s*("+NUMBER+")",body)
            assert len(matches)==1,key
            return float(matches[0])
        count=report("Current number of parcels")
        mass=report("Current mass in system")
        released=report("Mass transfer phase change")
        introduced=report("mass introduced")
        assert count in (0,1) and np.isfinite([mass,released,introduced]).all()
        assert mass>=0 and released>=0 and introduced==M0
        assert re.search(r"- escape\s*=\s*0,\s*0",body) and re.search(r"- stick\s*=\s*0,\s*0",body)
        all_reports.append({"step":i+1,"time_s":float(times[i]),"parcels":int(count),
                            "mass_kg":mass,"phase_change_cumulative_kg":released,
                            "native_mass_balance_relative":abs(introduced-mass-released)/M0})
    assert all(b["parcels"]<=a["parcels"] for a,b in zip(all_reports,all_reports[1:]))

    def carrier(name,kind="scalar"):
        a=read(native(name),"internalField",CELLS,kind)
        assert a.shape==((CELLS,3) if kind=="vector" else (CELLS,)) and np.isfinite(a).all(),name
        return a

    def parcel(folder,name):
        body=re.sub(r"/\*.*?\*/|//[^\n]*","",native(f"{folder}/{name}"),flags=re.S)
        tail=body[body.index("}")+1:]
        match=re.fullmatch(r"\s*1\s*\(\s*("+NUMBER+r")\s*\)\s*;?\s*",tail)
        assert match,name
        return float(match[1])

    actual={p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?",p.name)}
    assert actual=={"0",*(f"{t:g}" for t in TIMES)}
    rows=[]
    parcel_ids=set()
    for t in TIMES:
        name=f"{t:g}"
        clock=native(f"{name}/uniform/time")
        step=int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert 1<=step<=len(times) and abs(times[step-1]-t)<1e-10
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-t)<1e-10
        u=carrier(f"{name}/U","vector")
        p,temp=carrier(f"{name}/p"),carrier(f"{name}/T")
        air,water=carrier(f"{name}/air"),carrier(f"{name}/H2O")
        native_row=all_reports[step-1]
        folder=f"{name}/lagrangian/cloud"
        diameter=rho=multiplicity=parcel_temp=mass_error=None
        if native_row["parcels"]==1:
            diameter,rho,multiplicity,parcel_temp=(parcel(folder,key) for key in ("d","rho","nParticle","T"))
            assert all(np.isfinite(v) and v>0 for v in (diameter,rho,multiplicity,parcel_temp))
            parcel_ids.add(parcel(folder,"origId"))
            mass_error=float(abs(multiplicity*rho*np.pi*diameter**3/6-native_row["mass_kg"])/M0)
        else:
            assert native_row["mass_kg"]==0 and not (case/folder/"d").is_file()
        rows.append({**native_row,"diameter_m":diameter,"droplet_temperature_k":parcel_temp,
                     "droplet_density_kg_m3":rho,"multiplicity":multiplicity,"parcel_mass_error_relative":mass_error,
                     "carrier_pressure_min_pa":float(p.min()),"carrier_temperature_min_k":float(temp.min()),
                     "carrier_temperature_max_k":float(temp.max()),"carrier_temperature_max_departure_k":float(np.max(abs(temp-473))),
                     "carrier_maximum_speed_m_s":float(np.linalg.norm(u,axis=1).max()),
                     "carrier_species_min":float(min(air.min(),water.min())),
                     "carrier_species_max":float(max(air.max(),water.max())),
                     "carrier_species_sum_error":float(np.max(abs(air+water-1)))})
    assert rows[-1]["step"]==len(times) and len(parcel_ids)==1
    live=[row for row in rows if row["parcels"]]
    gates={"carrier_fields":all(r["carrier_pressure_min_pa"]>0 and r["carrier_temperature_min_k"]>0
                                and -1e-6<=r["carrier_species_min"]<=r["carrier_species_max"]<=1+1e-6
                                and r["carrier_species_sum_error"]<=2e-6 for r in rows),
           "parcel_interval":len(live)>=5 and all(r["parcels"]==1 for r in rows[:5]),
           "diameter_and_mass_loss":len(live)>=5 and live[-1]["diameter_m"]<live[0]["diameter_m"]
                                   and live[-1]["mass_kg"]<live[0]["mass_kg"],
           "mass_crosscheck":all(r["parcel_mass_error_relative"]<=1e-6 for r in live),
           "mass_budget":max(r["native_mass_balance_relative"] for r in all_reports)<=1e-6}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"steps":len(times),"saved_states":len(rows),
            "active_saved_states":len(live),"initial_prescribed_mass_kg":M0,"samples":rows,"native_sha256":sources}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    for filename,records in (("droplet-history.csv",rows),("native-cloud-history.csv",all_reports)):
        with (out/filename).open("x",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]))
            writer.writeheader()
            writer.writerows(records)

    fig,axs=plt.subplots(1,3,figsize=(11.3,3.8),layout="constrained")
    ts=[r["time_s"]*1000 for r in live]
    for ax,key,scale,initial,title,label in (
        (axs[0],"diameter_m",1e6,100,"Droplet diameter","Diameter (μm)"),
        (axs[1],"mass_kg",1e12,M0*1e12,"Liquid mass","Mass (ng)"),
        (axs[2],"droplet_temperature_k",1,350,"Droplet temperature","Temperature (K)")):
        ax.plot(ts,[scale*r[key] for r in live],color="#2b688e",lw=1.8)
        ax.scatter(0,initial,facecolors="none",edgecolors="#a75436",label="Prescribed initial")
        ax.set(title=title,xlabel="Time (ms)",ylabel=label,xlim=(-2,102))
        ax.grid(alpha=.2)
    axs[0].legend(frameon=False,fontsize=8)
    fig.suptitle("One-way evaporation of a water droplet in hot gas",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"evaporating-droplet.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        parser.add_argument(name,type=Path)
    a=parser.parse_args()
    raise SystemExit(0 if audit(a.case,a.reader,a.output) else 2)
