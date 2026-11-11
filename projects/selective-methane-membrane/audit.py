"""Audit the mapped membrane, native species fluxes and sleeve composition."""
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

SPECIES=("CH4","O2","N2")
SAVES=np.arange(1,11)*.05
PATCHES=("membranePipe","membraneSleeve")
N=r"([-+\d.eE]+)"
FLUXES=(("heatFlux","1 2 -3 0 0 0 0"),("specieFlux(O2)","1 0 -1 0 0 0 0"),
        ("specieFlux(CH4)","1 0 -1 0 0 0 0"),("specieFlux(N2)","1 0 -1 0 0 0 0"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case,reader_path,out):
    read=runpy.run_path(str(reader_path))["foam_values"]
    sources=[]

    def native(name):
        sources.append(name)
        return (case/name).read_text("ascii")

    for line in native("input-sha256.txt").splitlines():
        digest,name=line.split("  ",1)
        assert sha(case/name)==digest,name
    sources.extend(["pre-mesh-input-sha256.txt","pre-mesh-inputs-verified.log","snappy.log","baffles.log"])
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and log.rstrip().endswith("\nEnd")
    cells=int(re.search(r"^\s*cells:\s*(\d+)\s*$",mesh,re.M)[1])
    assert 0<cells<=150000
    steps=np.array(re.findall(r"^Time = "+N+r"s$",log,re.M),float)
    assert len(steps)>0 and np.all(np.diff(steps)>0) and steps[-1]==.5

    def field(name,kind="scalar"):
        value=read(native(name),"internalField",cells,kind)
        assert np.all(np.isfinite(value)),name
        return value

    c,v=field("0/C","vector"),field("0/Vc")
    assert c.shape==(cells,3) and np.all(v>0)
    zone=native("constant/polyMesh/cellZones")
    parsed=re.search(r"pipe\s*\{\s*cellLabels\s+List<label>\s+(\d+)\s*\(([^)]*)\)",zone)
    ids=np.fromstring(parsed[2],sep=" ",dtype=int)
    assert 0<len(ids)==int(parsed[1])<cells and len(np.unique(ids))==len(ids)
    assert np.all((ids>=0)&(ids<cells))
    pipe=np.zeros(cells,bool)
    pipe[ids]=True
    boundary=native("constant/polyMesh/boundary")
    points=np.array([np.fromstring(s,sep=" ") for s in re.findall(
        r"^\(([^()]+)\)$",native("constant/polyMesh/points"),re.M)])
    faces=[np.fromstring(s,sep=" ",dtype=int) for s in re.findall(
        r"^\d+\(([^()]+)\)$",native("constant/polyMesh/faces"),re.M)]
    owner=np.fromstring(re.search(r"\n"+str(len(faces))+r"\s*\(([^)]*)\)",native("constant/polyMesh/owner"))[1],
                        sep=" ",dtype=int)
    assert len(owner)==len(faces) and np.all((owner>=0)&(owner<cells))
    reports={}
    geometries={}
    for patch,partner,zone_side in zip(PATCHES,PATCHES[::-1],(True,False)):
        block=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",boundary)[1]
        assert re.search(r"\btype\s+mappedWall;",block) and re.search(r"neighbourPatch\s+"+partner+r";",block)
        count=int(re.search(r"nFaces\s+(\d+);",block)[1])
        start=int(re.search(r"startFace\s+(\d+);",block)[1])
        assert count>0 and start+count<=len(faces) and np.all(pipe[owner[start:start+count]]==zone_side)
        geometry={}
        for idx in range(start,start+count):
            coords=points[faces[idx]]
            area=.5*np.cross(coords,np.roll(coords,-1,axis=0)).sum(axis=0)
            assert np.linalg.norm(area)>0 and area@(coords.mean(axis=0)-c[owner[idx]])>0
            geometry[tuple(sorted(map(tuple,coords)))]=area
        assert len(geometry)==count
        geometries[patch]=geometry
        name=f"postProcessing/{patch}Fluxes/0/surfaceFieldValue.dat"
        body=native(name)
        assert f"patch({patch})" in body and int(re.search(r"^# Faces\s*:\s*(\d+)",body,re.M)[1])==count
        assert re.findall(r"sum\((?:[^()]|\([^()]*\))*\)",body)==[f"sum({key})" for key,_ in FLUXES]
        area_sum=sum(np.linalg.norm(area) for area in geometry.values())
        assert np.isclose(float(re.search(r"^# Area\s*:\s*"+N,body,re.M)[1]),area_sum,rtol=1e-6)
        values=np.loadtxt(case/name)
        assert values.shape==(len(steps)+1,5) and np.all(np.isfinite(values))
        assert values[0,0]==0 and np.allclose(values[1:,0],steps,rtol=0,atol=1e-8)
        reports[patch]=values
    assert geometries[PATCHES[0]].keys()==geometries[PATCHES[1]].keys()
    assert all(np.allclose(area,-geometries[PATCHES[1]][key],rtol=0,atol=1e-12)
               for key,area in geometries[PATCHES[0]].items())
    n_faces=len(geometries[PATCHES[0]])

    times=[f"{t:g}" for t in SAVES]
    actual={p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"0\.\d+",p.name)}
    assert actual==set(times)
    rows=[]
    for t,name in zip(SAVES,times):
        clock=native(f"{name}/uniform/time")
        step=int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-t)<1e-8
        assert abs(steps[step-1]-t)<1e-8
        temp,pressure=field(f"{name}/T"),field(f"{name}/p")
        field(f"{name}/U","vector")
        species=np.array([field(f"{name}/{s}") for s in SPECIES])
        assert species.shape==(3,cells)
        for col,(flux,dimensions) in enumerate(FLUXES,1):
            data=native(f"{name}/{flux}")
            assert re.search(r"dimensions\s+\["+dimensions+r"\];",data)
            for patch in PATCHES:
                section=re.search(r"\b"+patch+r"\s*\{([^}]+)\}",data,re.S)[1]
                patch_flux=read(section,"value",n_faces,"scalar")
                assert np.all(np.isfinite(patch_flux))
                assert np.isclose(patch_flux.sum(),reports[patch][step,col],rtol=1e-6,atol=1e-10)
        rows.append({"time_s":float(t),"step":step,"min_p_pa":float(pressure.min()),
                     "min_temperature_k":float(temp.min()),"max_temperature_k":float(temp.max()),
                     "min_mass_fraction":float(species.min()),"max_mass_fraction":float(species.max()),
                     "max_species_sum_error":float(abs(species.sum(axis=0)-1).max()),
                     "pipe_mean_methane_mass_fraction":float(np.average(species[0,pipe],weights=v[pipe])),
                     "sleeve_mean_methane_mass_fraction":float(np.average(species[0,~pipe],weights=v[~pipe])),
                     **{f"{patch}_{label}":float(reports[patch][step,i]) for patch in PATCHES
                        for i,label in enumerate(("heat_w","O2_kg_s","CH4_kg_s","N2_kg_s"),1)}})
    assert rows[-1]["step"]==len(steps)
    last=rows[-1]
    mf_pipe,mf_sleeve=(last[f"{patch}_CH4_kg_s"] for patch in PATCHES)
    mismatch=abs(mf_pipe+mf_sleeve)/max(abs(mf_pipe),abs(mf_sleeve))
    gates={"temperature_pressure":all(r["min_p_pa"]>0 and 270<=r["min_temperature_k"]<=r["max_temperature_k"]<=320 for r in rows),
           "species":all(r["min_mass_fraction"]>=-1e-6 and r["max_mass_fraction"]<=1+1e-6 and
                         r["max_species_sum_error"]<=2e-6 for r in rows),
           "methane_pair":mf_pipe>1e-5 and mf_sleeve<0 and mismatch<.005,
           "selectivity":all(abs(last[f"{patch}_{s}_kg_s"])<.01*abs(last[f"{patch}_CH4_kg_s"])
                             for patch in PATCHES for s in ("O2","N2")),
           "sleeve_methane":last["sleeve_mean_methane_mass_fraction"]>1e-5}
    result={"pass":all(gates.values()),"gates":gates,"cells":cells,"pipe_cells":len(ids),"steps":len(steps),
            "saved_states":len(rows),"membrane_faces_per_side":n_faces,
            "terminal_methane_flux_pair_relative_mismatch":mismatch,"samples":rows,
            "native_sha256":{name:sha(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"membrane-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    all_reports=np.column_stack([reports[PATCHES[0]],reports[PATCHES[1]][:,1:]])
    columns=["time_s"]+[f"{p}_{s}" for p in PATCHES for s in ("heat_w","O2_kg_s","CH4_kg_s","N2_kg_s")]
    np.savetxt(out/"native-membrane-flux.csv",all_reports,delimiter=",",header=",".join(columns),comments="")

    x,y,z=c.T
    slab=(z>=0)&(z<.1)
    fig,axes=plt.subplots(1,3,figsize=(12.1,4.3),layout="constrained")
    for ax,select,title,upper in ((axes[0],slab,"Pipe and sleeve",.5),
                                  (axes[1],slab&~pipe,"Sleeve detail",float(species[0,slab&~pipe].max()))):
        assert select.sum()>30
        dots=ax.scatter(x[select],y[select],c=species[0,select],s=11,linewidths=0,cmap="viridis",vmin=0,vmax=upper)
        ax.set(xlabel="x (m)",ylabel="y (m)",aspect="equal",title=title)
        fig.colorbar(dots,ax=ax,label="Methane mass fraction",shrink=.75)
    for p,label in zip(PATCHES,("Pipe side","Sleeve side")):
        axes[2].plot(reports[p][:,0],1000*reports[p][:,3],"o-",ms=2,label=label)
    axes[2].set(xlabel="Time (s)",ylabel="Outward methane flux (g/s)",title="Native membrane reports")
    axes[2].grid(alpha=.2)
    axes[2].legend(frameon=False)
    fig.suptitle("Selective methane transport through a membrane",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"selective-methane-membrane.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case","reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)