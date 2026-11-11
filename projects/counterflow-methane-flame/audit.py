"""Audit native counterflow combustion and plot heat release and species."""
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

NX, NY, CELLS = 100, 40, 4000
SPECIES = ("CH4","O2","N2","CO2","H2O")
SAVES = np.arange(1,11)*.05


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
    mesh,log=native("check-mesh.log"),native("solver.log")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*4000\s*$",mesh,re.M)
    assert log.rstrip().endswith("\nEnd")
    steps=np.array(re.findall(r"^Time = ([-+\d.eE]+)s$",log,re.M),float)
    assert len(steps)>0 and np.all(np.diff(steps)>0) and steps[-1]==.5
    assert 'CH4 + 2O2 = CO2 + 2H2O' in native("constant/reactions")

    def field(name,kind="scalar"):
        values=read(native(name),"internalField",CELLS,kind)
        assert np.all(np.isfinite(values)),name
        return values

    c,v=field("0/C","vector"),field("0/Vc")
    dx,dy=.02/NX,.02/NY
    assert np.allclose(v,dx*dy*.02,rtol=1e-8,atol=1e-15)
    assert np.max(abs(c[:,2]))<1e-10
    ij=np.rint((c[:,:2]-[0,-.01])/[dx,dy]-.5).astype(int)
    assert np.all((ij>=0)&(ij<[NX,NY]))
    assert np.allclose(c[:,:2],(ij+.5)*[dx,dy]+[0,-.01],rtol=0,atol=1e-10)
    key=ij[:,1]*NX+ij[:,0]
    assert np.array_equal(np.sort(key),np.arange(CELLS))
    order=np.argsort(key)

    def grid(values):
        return values[order].reshape(NY,NX)

    rows=[]
    names=[f"{time:g}" for time in SAVES]
    actual={p.name for p in case.iterdir() if p.is_dir() and re.fullmatch(r"0\.\d+",p.name)}
    assert actual==set(names)
    for time,name in zip(SAVES,names):
        clock=native(f"{name}/uniform/time")
        i=int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-time)<1e-8
        assert abs(steps[i-1]-time)<1e-8
        t,p,u=field(f"{name}/T"),field(f"{name}/p"),field(f"{name}/U","vector")
        qtext=native(f"{name}/Qdot")
        assert re.search(r"dimensions\s+\[1 -1 -3 0 0 0 0\];",qtext)
        q=read(qtext,"internalField",CELLS,"scalar")
        assert np.all(np.isfinite(q))
        species={s:field(f"{name}/{s}") for s in SPECIES}
        fractions=np.array(list(species.values()))
        peak=int(np.argmax(q))
        rows.append({"time_s":float(time),"step":i,"temperature_min_k":float(t.min()),
                     "temperature_max_k":float(t.max()),"min_pressure_pa":float(p.min()),
                     "max_speed_m_s":float(np.linalg.norm(u,axis=1).max()),
                     "minimum_species_fraction":float(fractions.min()),
                     "maximum_species_fraction":float(fractions.max()),
                     "max_species_sum_error":float(abs(fractions.sum(axis=0)-1).max()),
                     "peak_Qdot_w_m3":float(q[peak]),"integrated_Qdot_w":float(q@v),
                     "peak_Qdot_x_m":float(c[peak,0]),"peak_Qdot_y_m":float(c[peak,1]),
                     "CO2_at_peak":float(species["CO2"][peak]),"H2O_at_peak":float(species["H2O"][peak])})
    assert rows[-1]["step"]==len(steps)
    final=rows[-1]
    gates={"temperature_pressure":all(r["min_pressure_pa"]>0 and r["temperature_min_k"]>=250 and
                                       r["temperature_max_k"]<=3500 for r in rows),
           "species":all(r["minimum_species_fraction"]>=-1e-6 and r["maximum_species_fraction"]<=1+1e-6 and
                         r["max_species_sum_error"]<=2e-6 for r in rows),
           "reaction":final["peak_Qdot_w_m3"]>1e5 and final["integrated_Qdot_w"]>1e-3 and
                      .001<=final["peak_Qdot_x_m"]<=.019 and
                      min(final["CO2_at_peak"],final["H2O_at_peak"])>1e-4}
    result={"pass":all(gates.values()),"gates":gates,"cells":CELLS,"steps":len(steps),"saved_states":len(rows),
            "samples":rows,"native_sha256":{name:sha(case/name) for name in sources}}
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    with (out/"flame-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig,axes=plt.subplots(2,2,figsize=(10,7.8),layout="constrained")
    xedges=np.arange(NX+1)*dx*1000
    yedges=(np.arange(NY+1)*dy-.01)*1000
    for ax,values,title,label,cmap in ((axes[0,0],t,"Temperature at 0.5 s","Temperature (K)","inferno"),
                                       (axes[0,1],q/1e6,"Chemical heat release at 0.5 s","Qdot (MW/m³)","magma")):
        im=ax.pcolormesh(xedges,yedges,grid(values),shading="flat",cmap=cmap)
        ax.set(xlabel="Distance from fuel inlet (mm)",ylabel="Transverse position (mm)",title=title,aspect="equal")
        ax.axhline(0,ls="--",lw=.8,color="white",alpha=.6)
        fig.colorbar(im,ax=ax,label=label,shrink=.8)
    x=(xedges[:-1]+xedges[1:])/2
    centre=lambda data:grid(data)[NY//2-1:NY//2+1].mean(axis=0)
    for name in ("CH4","O2","CO2","H2O"):
        axes[1,0].plot(x,centre(species[name]),label=name)
    axes[1,0].set(xlabel="Distance from fuel inlet (mm)",ylabel="Mass fraction",title="Species across the centreline")
    axes[1,0].legend(frameon=False,ncol=2)
    axes[1,1].plot([r["time_s"] for r in rows],[r["integrated_Qdot_w"] for r in rows],"o-",ms=3,color="#bc562c")
    axes[1,1].set(xlabel="Time (s)",ylabel="Integrated release (W)",title="Sampled chemical heat release")
    for ax in axes[1]:
        ax.grid(alpha=.2)
    fig.suptitle("Opposed-jet methane diffusion flame",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(out/f"counterflow-methane-flame.{ext}",dpi=180,
                    metadata={"CreationDate":None,"ModDate":None} if ext=="pdf" else None)
    plt.close(fig)
    np.savetxt(out/"flame-centreline.csv",np.column_stack([x,centre(t),centre(q),*[centre(species[s]) for s in SPECIES]]),
               delimiter=",",comments="",header="distance_mm,T_K,Qdot_W_m3,"+",".join(SPECIES))
    return result["pass"]


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("case","reader","output"):
        parser.add_argument(argument,type=Path)
    args=parser.parse_args()
    raise SystemExit(0 if audit(args.case,args.reader,args.output) else 2)