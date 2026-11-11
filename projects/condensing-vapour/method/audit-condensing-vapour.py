"""Check native phase fractions, composition and condensation source fields."""
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
from matplotlib.patches import Rectangle
import numpy as np

CELLS=8000
TIMES=np.arange(1,11)/50
NUMBER=r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"


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
    for name in ("pre-mesh-input-sha256.txt","pre-mesh-inputs-verified.log","inputs-verified.log","mesh.log",
                 "zones.log","centres.log","volumes.log"):
        native(name)
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*8000\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd")
    assert "Selecting wallCondensation with name wallCondensation" in log
    assert "Selecting phaseSurfaceCondensation with name solidCondensation" in log
    steps=list(re.finditer(r"^Time = ("+NUMBER+r")s?$",log,re.M))
    times=np.array([float(m[1]) for m in steps])
    assert len(times) and np.isfinite(times).all() and np.all(np.diff(times)>0) and times[-1]==.2
    zone=native("constant/polyMesh/cellZones")
    match=re.search(r"\bdryer\s*\{\s*cellLabels\s+List<label>\s+1000\s*\(([^()]+)\)",zone,re.S)
    indices=np.fromstring(match[1],sep=" ",dtype=int)
    assert len(indices)==len(np.unique(indices))==1000 and np.all((indices>=0)&(indices<CELLS))
    def field(name,kind="scalar",key="internalField"):
        a=read(native(name),key,CELLS,kind)
        assert a.shape==((CELLS,3) if kind=="vector" else (CELLS,)) and np.isfinite(a).all(),name
        return a

    c,v=field("0/C","vector"),field("0/Vc")
    assert np.all(v>0) and abs(v.sum()-.068)<1e-11
    bed=np.zeros(CELLS,bool)
    bed[indices]=True
    geometric_bed=(c[:,0]>.4)&(c[:,0]<.6)&(c[:,1]>0)
    np.testing.assert_array_equal(bed,geometric_bed)
    assert re.search(r"defaultValue\s+0;",native("0/alpha.liquid"))
    initial_h=field("0/H2O.gas")
    assert np.all(initial_h==.018)
    actual={p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?",p.name)}
    assert actual=={"0",*(f"{t:g}" for t in TIMES)}
    rows=[]
    for t in TIMES:
        name=f"{t:g}"
        clock=native(f"{name}/uniform/time")
        step=int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert 1<=step<=len(times) and abs(times[step-1]-t)<1e-10
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-t)<1e-10
        gas,liq,solid=(field(f"{name}/alpha.{phase}") for phase in ("gas","liquid","solid"))
        h,n,o=(field(f"{name}/{specie}.gas") for specie in ("H2O","N2","O2"))
        temperatures=[field(f"{name}/T.{phase}") for phase in ("gas","liquid","solid")]
        for phase in ("gas","liquid"):
            field(f"{name}/U.{phase}","vector")
        p=field(f"{name}/p")
        rates={}
        for source,key in (("wallCondensation","internalField"),("solidCondensation","value")):
            # WSL maps the native colon to the Windows private-use character on this volume.
            paths=list((case/name).glob(f"{source}?mDot"))
            assert len(paths)==1
            path=paths[0].relative_to(case).as_posix()
            text=native(path)
            assert re.search(r"object\s+"+source+r":mDot;",text)
            assert re.search(r"dimensions\s+\[1\s+-3\s+-1\s+0\s+0\s+0\s+0\];",text)
            mdot=read(text,key,CELLS,"scalar")
            assert mdot.shape==(CELLS,) and np.isfinite(mdot).all()
            rates[f"{source}_source_kg_s"]=float(mdot@v)
            if source=="solidCondensation":
                assert np.max(abs(mdot[~bed]))<1e-12
                stop=steps[step].start() if step<len(steps) else len(log)
                block=log[steps[step-1].end():stop]
                report=re.findall(r"^    min/average/max mDot = ("+NUMBER+r")/("+NUMBER+r")/("+NUMBER+r") ",
                                  block,re.M)
                assert report
                np.testing.assert_allclose([mdot.min(),mdot.mean(),mdot.max()],list(map(float,report[-1])),
                                           atol=1e-9,rtol=1e-7)
        rows.append({"time_s":float(t),"step":step,"phase_min":float(min(gas.min(),liq.min(),solid.min())),
                     "phase_max":float(max(gas.max(),liq.max(),solid.max())),
                     "phase_sum_max_error":float(np.max(abs(gas+liq+solid-1))),
                     "species_min":float(min(h.min(),n.min(),o.min())),
                     "species_max":float(max(h.max(),n.max(),o.max())),
                     "species_sum_max_error":float(np.max(abs(h+n+o-1))),
                     "minimum_p_pa":float(p.min()),"maximum_p_pa":float(p.max()),
                     "minimum_temperature_k":float(min(a.min() for a in temperatures)),
                     "solid_fraction_max_error":float(np.max(abs(solid-.2*bed))),
                     "liquid_volume_m3":float(liq@v),
                     "dryer_gas_water_mass_fraction":float(np.average(h[bed],weights=(v*gas)[bed])),
                     **rates})
    assert rows[-1]["step"]==len(times)
    gates={"fractions":all(-1e-6<=r["phase_min"]<=r["phase_max"]<=1+1e-6
                          and -1e-6<=r["species_min"]<=r["species_max"]<=1+1e-6
                          and max(r["phase_sum_max_error"],r["species_sum_max_error"])<=2e-6 for r in rows),
           "thermodynamics":all(r["minimum_p_pa"]>0 and r["minimum_temperature_k"]>0 for r in rows),
           "fixed_solid":max(r["solid_fraction_max_error"] for r in rows)<=1e-12,
           "liquid_production":rows[-1]["liquid_volume_m3"]>1e-8,
           "dryer_humidity_reduction":rows[-1]["dryer_gas_water_mass_fraction"]<=.0179}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"steps":len(times),"saved_states":len(rows),
            "dryer_cells":len(indices),"samples":rows,"native_sha256":sources}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"condensation-history.csv").open("x",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    xedges=np.linspace(0,1.4,71)
    yedges=np.r_[np.linspace(-.3,0,11)[:-1],np.linspace(0,.2,11)]
    sample=abs(c[:,2]-.09)<1e-9
    ix=np.searchsorted(xedges,c[sample,0])-1
    iy=np.searchsorted(yedges,c[sample,1])-1
    assert sample.sum()==800 and len(np.unique(iy*70+ix))==800
    fig,axes=plt.subplots(2,2,figsize=(11,6.2),layout="constrained")
    for ax,val,title,label,cmap in ((axes[0,0],h*100,"Gas humidity","Water vapour mass fraction (%)","viridis"),
                                   (axes[0,1],liq*1e6,"Liquid fraction","Liquid volume fraction (ppm)","magma")):
        grid=np.full((20,70),np.nan)
        grid[iy,ix]=val[sample]
        image=ax.pcolormesh(xedges,yedges,grid,cmap=cmap,shading="flat")
        ax.add_patch(Rectangle((.4,0),.2,.2,fill=False,edgecolor="#e9e9ed",ls="--",lw=1))
        ax.set(title=title+" at 0.2 s",xlabel="x (m)",ylabel="y (m)",aspect="equal")
        fig.colorbar(image,ax=ax,label=label,shrink=.82)
    ax=axes[1,0]
    ax.plot(np.r_[0,TIMES],[0]+[1e6*r["liquid_volume_m3"] for r in rows],"o-",ms=3,color="#235a82")
    ax.set(title="Liquid inventory",xlabel="Time (s)",ylabel="Volume (mL)")
    ax=axes[1,1]
    for name,label in (("wallCondensation","Wall model"),("solidCondensation","Solid-surface model")):
        ax.plot(TIMES,[1000*r[f"{name}_source_kg_s"] for r in rows],"o-",ms=3,label=label)
    ax.set(title="Saved condensation sources",xlabel="Time (s)",ylabel="Volume-integrated source (g/s)")
    ax.legend(frameon=False)
    for ax in axes[1]:
        ax.grid(alpha=.2)
    fig.suptitle("Condensation from humid air in a cooled duct",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"condensing-vapour.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","reader","output"):
        parser.add_argument(name,type=Path)
    a=parser.parse_args()
    raise SystemExit(0 if audit(a.case,a.reader,a.output) else 2)
