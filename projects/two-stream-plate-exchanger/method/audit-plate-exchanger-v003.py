"""Read native mesh and fields and check heat transfer, pressure and mass."""
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

N = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
PORT_CP = 4180.
T_REF = 300.
DP_REF = 5.4


def audit(case, read, direction, fine):
    nx, ny, nw = (180, 24, 6) if fine else (120, 16, 4)
    sources = {}

    def native(name):
        raw = (case/name).read_bytes()
        sources[name] = hashlib.sha256(raw).hexdigest()
        return raw.decode("ascii")

    for manifest in ("pre-mesh-sha256.txt", "input-sha256.txt"):
        for row in native(manifest).splitlines():
            sha, name = row.split("  ", 1)
            native(name)
            assert sources[name] == sha, name
    for f in ("generator-sha256.txt", "mesh.log", "split.log", "inputs-verified.log",
              "original-inputs-verified.log"):
        native(f)
    log = native("solver.log")
    assert log.rstrip().endswith("End") and "14-7b05503f98a8" in log
    times = list(re.finditer(r"^\s*Time = (\d+)s?\s*$", log, re.M))
    assert [int(t[1]) for t in times] == list(range(1, 901))
    keys = [("cold","Ux"),("hot","Ux"),("cold","Uy"),("hot","Uy"),("cold","h"),
            ("hot","h"),("plate","e"),("cold","p"),("hot","p")]
    residuals = []
    for i, time in enumerate(times):
        end = times[i+1].start() if i+1 < len(times) else len(log)
        chunk = log[time.end():end]
        values = re.findall(r"^(cold|hot|plate)\s+\S+:\s+Solving for (\w+), Initial residual = ("+N+")", chunk, re.M)
        assert len(values) == len(keys) and {(r,f) for r,f,_ in values} == set(keys), i+1
        residuals.append(max(float(v) for _,_,v in values))
    assert np.isfinite(residuals).all()
    saves = sorted(int(p.name) for p in case.iterdir() if p.is_dir() and p.name.isdigit() and int(p.name)>0)
    assert saves == list(range(100, 901, 100))
    for step in saves:
        clk = native(f"{step}/uniform/time")
        assert int(re.search(r"^index\s+(\d+);",clk,re.M)[1]) == step
        assert float(re.search(r"^value\s+([^;]+);",clk,re.M)[1]) == step

    regions = {}
    for region, rows, interfaces, conductivity in (
            ("cold",ny,("plate",),.6), ("plate",nw,("hot","cold"),16.), ("hot",ny,("plate",),.6)):
        prefix = f"constant/{region}/polyMesh"
        count = nx*rows
        mesh = native(f"check-{region}.log")
        assert "Mesh OK." in mesh and re.search(rf"^\s*cells:\s*{count}\s*$",mesh,re.M)
        c = read(native(f"0/{region}/C"),"internalField",count,"vector")
        volume = read(native(f"0/{region}/Vc"),"internalField",count,"scalar")
        assert c.shape == (count,3) and (volume > 0).all()
        expected_volume = .12*.01*(.0005 if region=="plate" else .002)
        assert abs(volume.sum()-expected_volume) < 1e-14
        owner = np.array(re.findall(r"^([0-9]+)$",native(f"{prefix}/owner"),re.M),dtype=int)
        assert owner[0] == len(owner)-1
        owner = owner[1:]
        assert owner.min()==0 and owner.max()==count-1
        faces = np.array(re.findall(r"^4\((\d+) (\d+) (\d+) (\d+)\)$",
                                    native(f"{prefix}/faces"),re.M),dtype=int)
        points = np.array(re.findall(r"^\(("+N+r") ("+N+r") ("+N+r")\)$",
                                     native(f"{prefix}/points"),re.M),dtype=float)
        assert len(faces) == len(owner)
        boundary = native(f"{prefix}/boundary")
        patches = {}
        walls = ("outer_lower",) if region=="cold" else ("outer_upper",) if region=="hot" else ()
        expected_patches = {f"{region}_left": rows, f"{region}_right": rows,
                            **{f"{region}_to_{other}":nx for other in interfaces},
                            **{wall:nx for wall in walls}}
        for patch, n in expected_patches.items():
            body = re.search(r"\b"+patch+r"\s*\{([^}]+)\}",boundary,re.S)[1]
            start = int(re.search(r"startFace\s+(\d+);",body)[1])
            assert int(re.search(r"nFaces\s+(\d+);",body)[1])==n
            axis = 0 if patch.endswith(("_left","_right")) else 1
            v = points[faces[start:start+n]]
            center = v.mean(axis=1)
            area = np.linalg.norm(np.cross(v[:,1]-v[:,0],v[:,2]-v[:,1]),axis=1)
            cells = owner[start:start+n]
            gap = abs(c[cells,axis]-center[:,axis])
            assert (area>0).all() and (gap>0).all()
            if axis == 1:
                assert abs(area.sum()-.12*.01)<1e-14
            if "_to_" in patch:
                peer = patch.split("_to_")[1]
                assert "mappedWall;" in body
                assert re.search(r"neighbourRegion\s+"+peer+r";",body)
                assert re.search(r"neighbourPatch\s+"+peer+"_to_"+region+r";",body)
            patches[patch] = (center, area, gap, cells)
        regions[region] = (count,c,patches,conductivity)

    def patch_values(text,name,size):
        body = re.search(r"\b"+name+r"\s*\{([^}]+)\}",text,re.S)[1]
        return read(body,"value",size,"scalar")

    snapshots = []
    for step in saves:
        results = {}
        pairs = {}
        for region,(count,c,patches,k) in regions.items():
            t_text = native(f"{step}/{region}/T")
            tc = read(t_text,"internalField",count,"scalar")
            interfaces = {}
            for name,(centre,area,gap,cells) in patches.items():
                if "_to_" in name:
                    tp = patch_values(t_text,name,len(cells))
                    index = np.argsort(centre[:,0])
                    interfaces[name.split("_to_")[1]] = (centre[index],tp[index],(k*area*(tc[cells]-tp)/gap)[index])
            pairs[region] = interfaces
            metrics = {"min_T_k":float(tc.min()), "max_T_k":float(tc.max())}
            if region != "plate":
                p_text = native(f"{step}/{region}/p")
                p = read(p_text,"internalField",count,"scalar")
                u = read(native(f"{step}/{region}/U"),"internalField",count,"vector")
                assert p.shape == (count,) and u.shape == (count,3)
                phi = native(f"{step}/{region}/phi")
                assert re.search(r"dimensions\s+\[1 0 -1 0 0 0 0\];",phi)
                ports = {}
                wall_flux = 0.
                for name,(_,area,gap,cells) in patches.items():
                    mass = patch_values(phi,name,len(cells))
                    if not name.endswith(("_left","_right")):
                        wall_flux = max(wall_flux,float(np.max(abs(mass))))
                        continue
                    tf = patch_values(t_text,name,len(cells))
                    rate = float(mass.sum())
                    assert np.sign(mass).min() == np.sign(mass).max() == np.sign(rate) != 0
                    ports[name] = (rate,float(mass@tf/rate),float(mass@(tf-T_REF)*PORT_CP),
                                   float(np.sum(k*area*(tc[cells]-tf)/gap)),
                                   float(np.average(p[cells],weights=area)),area,gap,cells)
                inlet, outlet = sorted(ports, key=lambda name: ports[name][0])
                expected_inlet = f"{region}_right" if direction=="counter" and region=="cold" else f"{region}_left"
                assert inlet == expected_inlet
                a,b = ports[inlet],ports[outlet]
                assert a[0] < 0 < b[0]
                assert abs(-a[0]-.0003)<1e-13
                expected_temp = 340 if region=="hot" else 300
                assert abs(a[1]-expected_temp)<1e-10
                outlet_pressure = patch_values(p_text,outlet,len(b[-1]))
                assert np.max(abs(outlet_pressure)) <= 1e-14
                gauge_outlet = float(np.average(outlet_pressure,weights=b[5]))
                # Extrapolate inlet pressure by half a cell using the outlet gauge
                # and the interior-to-outlet distance in the fully developed channel.
                dp = (a[4]-gauge_outlet)*.12/(.12-float(np.mean(a[6])))
                metrics.update(mass_in_kg_s=-a[0],mass_out_kg_s=b[0],inlet_T_k=a[1],outlet_T_k=b[1],
                               relative_mass_error=abs(a[0]+b[0])/abs(a[0]),
                               maximum_sealed_mass_flux_kg_s=wall_flux,
                               sensible_advection_w=a[2]+b[2],port_conduction_w=a[3]+b[3],
                               total_port_sensible_w=a[2]+b[2]+a[3]+b[3],
                               pressure_drop_pa=dp, native_inlet_mean_p_pa=a[4],
                               native_outlet_gauge_p_pa=gauge_outlet)
            results[region] = metrics
        interfaces = {}
        for region in ("hot","cold"):
            fc,tf,qf = pairs[region]["plate"]
            sc,ts,qs = pairs["plate"][region]
            assert np.max(abs(fc-sc))<1e-12
            interfaces[region] = {"maximum_temperature_jump_k":float(np.max(abs(tf-ts))),
                                  "fluid_outward_heat_w":float(qf.sum()), "plate_outward_heat_w":float(qs.sum()),
                                  "face_l1_heat_mismatch_w":float(np.sum(abs(qf+qs)))}
        qh = -results["hot"]["total_port_sensible_w"]
        qc = results["cold"]["total_port_sensible_w"]
        qmax = min(results[x]["mass_in_kg_s"] for x in ("hot","cold"))*PORT_CP*40
        snapshots.append({"iteration":step, "regions":results, "interfaces":interfaces,
                          "hot_net_rejection_w":qh, "cold_net_gain_w":qc,
                          "hot_effectiveness":qh/qmax, "cold_effectiveness":qc/qmax})

    last, previous = snapshots[-1],snapshots[-2]
    fluid_regions = ("hot","cold")
    peak = max(last["hot_net_rejection_w"],last["cold_net_gain_w"])
    closure = abs(last["hot_net_rejection_w"]-last["cold_net_gain_w"])/peak
    duty_change = max(abs(last[k]-previous[k])/abs(last[k])
                      for k in ("hot_net_rejection_w","cold_net_gain_w"))
    gates = {
        "initial_residuals": residuals[-1] < 1e-5,
        "temperature_bounds": all(299.99<=r["min_T_k"]<=r["max_T_k"]<=340.01
                                  for row in snapshots for r in row["regions"].values()),
        "boundary_mass": all(row["regions"][r]["relative_mass_error"]<1e-4
                             and row["regions"][r]["maximum_sealed_mass_flux_kg_s"]<1e-12
                             for row in snapshots for r in fluid_regions),
        "pressure_drop": all(abs(last["regions"][r]["pressure_drop_pa"]/DP_REF-1)<.05 for r in fluid_regions),
        "positive_duties": 0<last["hot_effectiveness"]<1 and 0<last["cold_effectiveness"]<1,
        "port_energy": closure<=.01,
        "duty_stability": duty_change<=.01,
        "interface_continuity": all(p["maximum_temperature_jump_k"]<1e-5
                                    and p["face_l1_heat_mismatch_w"]/peak<.01
                                    for p in last["interfaces"].values()),
    }
    terminal_fields = {}
    for region,(count,c,_,_) in regions.items():
        t = read(native(f"900/{region}/T"),"internalField",count,"scalar")
        u = None if region=="plate" else read(native(f"900/{region}/U"),"internalField",count,"vector")
        xy = np.round(c[:,:2], 12)
        assert np.max(abs(xy-c[:,:2]))<1e-12
        x,y = np.unique(xy[:,0]),np.unique(xy[:,1])
        assert len(x)==nx and len(x)*len(y)==count
        assert len(np.unique(xy,axis=0))==count
        ix,iy = np.searchsorted(x,xy[:,0]),np.searchsorted(y,xy[:,1])
        tm = np.empty((len(y),nx))
        tm[iy,ix] = t
        bulk = None
        if u is not None:
            um = np.empty_like(tm)
            um[iy,ix] = abs(u[:,0])
            bulk = (um*tm).sum(axis=0)/um.sum(axis=0)
        terminal_fields[region] = (x,y,tm,bulk)
    summary = {"case":case.name,"direction":direction,"fine":fine,"cells":sum(r[0] for r in regions.values()),
               "iterations":900,"saves":len(saves),"gates":gates,"pass":all(gates.values()),
               "relative_sensible_port_energy_error":closure, "relative_terminal_duty_change":duty_change,
               "terminal_maximum_initial_residual":residuals[-1],"saved_metrics":snapshots,"native_sha256":sources}
    return summary,terminal_fields,residuals


def compare(root, reader, out):
    read = runpy.run_path(str(reader))["foam_values"]
    definitions = [("plate-counter-base-v001","counter",False),
                   ("plate-parallel-base-v001","parallel",False),
                   ("plate-counter-fine-v001","counter",True)]
    results = [audit(root/name,read,direction,fine) for name,direction,fine in definitions]
    for region in ("cold","hot","plate"):
        for name in ("points","faces","owner","neighbour","boundary"):
            mesh = f"constant/{region}/polyMesh/{name}"
            assert results[0][0]["native_sha256"][mesh] == results[1][0]["native_sha256"][mesh], mesh
    base,parallel,fine = [row[0]["saved_metrics"][-1] for row in results]
    duty = lambda row:(row["hot_net_rejection_w"]+row["cold_net_gain_w"])/2
    change = abs(duty(fine)/duty(base)-1)
    gates = {"counterflow_improvement": all(base[k]>parallel[k] for k in ("hot_effectiveness","cold_effectiveness")),
             "counterflow_refinement":change<=.03}
    summary = {"pass":all(row[0]["pass"] for row in results) and all(gates.values()),"comparison_gates":gates,
               "counterflow_relative_duty_mesh_change":change,"pressure_drop_analytical_pa":DP_REF,
               "cases":[r[0] for r in results]}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    labels = ("Counterflow", "Parallel flow", "Counterflow, refined")
    with (out/"port-history.csv").open("x",newline="") as stream:
        w=csv.writer(stream)
        w.writerow(("case","iteration","hot_outlet_k","cold_outlet_k","hot_net_rejection_w",
                    "cold_net_gain_w","hot_pressure_drop_pa","cold_pressure_drop_pa",
                    "hot_effectiveness","cold_effectiveness"))
        for label,(check,_,_) in zip(labels,results):
            for row in check["saved_metrics"]:
                h,c = row["regions"]["hot"],row["regions"]["cold"]
                w.writerow((label,row["iteration"],h["outlet_T_k"],c["outlet_T_k"],
                            row["hot_net_rejection_w"],row["cold_net_gain_w"],h["pressure_drop_pa"],
                            c["pressure_drop_pa"],row["hot_effectiveness"],row["cold_effectiveness"]))
    with (out/"temperature-profiles.csv").open("x",newline="") as stream:
        w=csv.writer(stream)
        w.writerow(("case","stream","x_m","mass_weighted_temperature_k"))
        for label,(_,fields,_) in zip(labels,results):
            for region in ("hot","cold"):
                x,_,_,bulk = fields[region]
                w.writerows((label,region,a,b) for a,b in zip(x,bulk))

    fig = plt.figure(figsize=(11,7),layout="constrained")
    grid = fig.add_gridspec(3,2,height_ratios=[1,1,1.3])
    fig.suptitle("Two-stream plate heat exchanger",fontweight="bold",fontsize=17)
    maps=[]
    for row,(_,fields,_) in enumerate(results[:2]):
        ax=fig.add_subplot(grid[row,:])
        maps.append(ax)
        for x,y,tm,_ in fields.values():
            image=ax.pcolormesh(x*1000,y*1000,tm,shading="nearest",cmap="inferno",vmin=300,vmax=340)
        ax.axhline(2,color="#fff",lw=.6)
        ax.axhline(2.5,color="#fff",lw=.6)
        ax.set(title=labels[row],xlabel="Axial position (mm)",ylabel="Channel height (mm)",
               xlim=(0,120),ylim=(0,4.5))
        ax.text(4,3.4,"Hot →",color="white",fontsize=9,
                bbox={"facecolor":"#202436","alpha":.72,"edgecolor":"none","pad":2})
        ax.text(4,.7,"Cold "+("←" if row==0 else "→"),color="white",fontsize=9,
                bbox={"facecolor":"#202436","alpha":.72,"edgecolor":"none","pad":2})
    fig.colorbar(image,ax=maps,label="Temperature (K)")
    ax=fig.add_subplot(grid[2,0])
    for label,(_,fields,_),style in zip(labels,results,("-", "--", ":")):
        for region,color in (("hot","#bd4b24"),("cold","#286f9f")):
            x,_,_,bulk = fields[region]
            ax.plot(x*1000,bulk,ls=style,color=color,label=f"{label}, {region}",lw=1.55)
    ax.set(xlabel="Axial position (mm)",ylabel="Bulk temperature (K)")
    ax.legend(frameon=False,ncols=2,fontsize=7)
    ax.grid(alpha=.18)
    ax=fig.add_subplot(grid[2,1])
    xa=np.arange(3)
    ax.bar(xa-.18,[a["hot_net_rejection_w"] for a in (base,parallel,fine)],.34,color="#bd4b24",label="Hot stream rejection")
    ax.bar(xa+.18,[a["cold_net_gain_w"] for a in (base,parallel,fine)],.34,color="#286f9f",label="Cold stream gain")
    ax.set(xticks=xa,xticklabels=("Counterflow","Parallel","Refined counterflow"),ylabel="Net sensible duty (W)")
    ax.legend(frameon=True,fontsize=8,loc="lower center",ncols=2,facecolor="white",framealpha=.94)
    ax.grid(axis="y",alpha=.18)
    fig.savefig(out/"plate-exchanger.png",dpi=185)
    fig.savefig(out/"plate-exchanger.pdf")
    plt.close(fig)
    for r,_,_ in results:
        print(r["case"],r["pass"],r["gates"],r["relative_sensible_port_energy_error"])
        print(r["saved_metrics"][-1])
    print("comparison",gates,change)
    return summary["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("root","reader","output"):
        parser.add_argument(arg,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if compare(args.root,args.reader,args.output) else 2)
