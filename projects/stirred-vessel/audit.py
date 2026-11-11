"""Check the rotating-frame vessel and plot fields on native front-plane faces."""
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
from matplotlib.collections import PolyCollection
import numpy as np

CELLS = 3072
ITERATIONS = tuple(range(50, 501, 50))
NUMBER = r"[-+0-9.eE]+"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mesh_geometry(case, centres, volumes):
    mesh = case/"constant/polyMesh"
    points_text = (mesh/"points").read_text(encoding="ascii")
    points = np.array([np.fromstring(line, sep=" ") for line in
                       re.findall(r"^\(([0-9.eE +\-]+)\)$", points_text, re.M)])
    assert points.shape == (6528, 3)
    faces_text = (mesh/"faces").read_text(encoding="ascii")
    faces = np.array([np.fromstring(line, dtype=int, sep=" ") for line in
                      re.findall(r"^4\(([0-9 ]+)\)$", faces_text, re.M)])
    assert faces.shape == (12480, 4) and np.all((faces>=0)&(faces<len(points)))
    owner_text = (mesh/"owner").read_text(encoding="ascii")
    owners = np.fromstring(re.search(r"\n12480\s*\(\s*([^)]*)\)", owner_text)[1], dtype=int, sep=" ")
    assert len(owners) == len(faces)
    faces = points[faces]
    front = np.all(np.isclose(faces[:,:,2], -.005, rtol=0, atol=1e-12), axis=1)
    order = np.argsort(owners[front])
    assert np.array_equal(owners[front][order], np.arange(CELLS))
    polygons = faces[front][order,:,:2]
    x, y = polygons.transpose(2,0,1)
    xx, yy = np.roll(x,-1,axis=1), np.roll(y,-1,axis=1)
    cross = x*yy-y*xx
    twice_area = cross.sum(axis=1)
    centroid = np.column_stack((((x+xx)*cross).sum(axis=1), ((y+yy)*cross).sum(axis=1)))/(3*twice_area[:,None])
    assert np.allclose(centroid, centres[:,:2], rtol=0, atol=1e-6)
    assert np.allclose(.005*np.abs(twice_area),volumes,rtol=1e-5,atol=1e-12)
    zones = (mesh/"cellZones").read_text(encoding="ascii")
    groups = []
    for name in ("rotor", "stator"):
        ids = re.search(r"\b"+name+r"\s*\{\s*cellLabels\s+List<label>\s+1536\s*\(([^)]*)\)",zones,re.S)
        assert ids
        groups.append(np.fromstring(ids[1],dtype=int,sep=" "))
    assert np.array_equal(np.sort(np.concatenate(groups)),np.arange(CELLS))
    return polygons, groups


def audit(case, reader_path, output):
    read = runpy.run_path(str(reader_path))["foam_values"]

    def field(name, kind="scalar"):
        return read((case/name).read_text(encoding="ascii"), "internalField", CELLS, kind)

    for line in (case/"input-sha256.txt").read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ",1)
        assert sha(case/name)==expected,name
    mesh=(case/"check-mesh.log").read_text(encoding="ascii")
    log=(case/"solver.log").read_text(encoding="ascii")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*3072\s*$",mesh,re.M)
    assert re.search(r"^End\s*$",log,re.M)
    assert np.array_equal(np.array(re.findall(r"^Time = (\d+)s$",log,re.M),dtype=int),np.arange(1,501))
    terminal=log[log.rindex("Time = 500s"):]
    continuity=np.array(re.findall(r"sum local = ("+NUMBER+r")",terminal),dtype=float)
    assert len(continuity) and continuity[-1]<1e-4
    centres=field("0/C","vector")
    vol=field("0/Vc")
    assert np.all(vol>0)
    polygons,(rotor,stator)=mesh_geometry(case,centres,vol)
    x,y=centres[:,:2].T
    radius=np.hypot(x,y)
    assert np.all(radius>0)
    rows, states = [], {}
    sources=["input-sha256.txt","check-mesh.log","solver.log","0/C","0/Vc"]+[
        f"constant/polyMesh/{name}" for name in ("points","faces","owner","cellZones")]
    for index in ITERATIONS:
        clock=(case/str(index)/"uniform/time").read_text(encoding="ascii")
        assert int(re.search(r"^index\s+(\d+);",clock,re.M)[1])==index
        assert float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])==index
        u=field(f"{index}/U","vector")
        p=field(f"{index}/p")
        k=field(f"{index}/k")
        e=field(f"{index}/epsilon")
        assert np.all(k>=0) and np.all(e>=0)
        speed=(-y*u[:,0]+x*u[:,1])/radius
        rows.append({"iteration":index,
                     "rotor_mean_tangential_speed_m_s":float(np.average(speed[rotor],weights=vol[rotor])),
                     "stator_mean_tangential_speed_m_s":float(np.average(speed[stator],weights=vol[stator])),
                     "max_velocity_magnitude_m_s":float(np.max(np.linalg.norm(u,axis=1))),
                     "min_k_m2_s2":float(k.min()), "min_epsilon_m2_s3":float(e.min())})
        states[index]=(u,p)
        sources.extend(f"{index}/{name}" for name in ("uniform/time","U","p","k","epsilon"))
    u,p=states[500]
    rms=float(np.sqrt(np.average(np.sum((u-states[450][0])**2,axis=1),weights=vol)))
    assert rms<.02 and rows[-1]["rotor_mean_tangential_speed_m_s"]>0
    checks={"cells":CELLS,"rotor_cells":len(rotor),"stator_cells":len(stator),"iterations":500,
            "saved_states":len(rows),"terminal_local_continuity":float(continuity[-1]),
            "velocity_change_rms_450_to_500_m_s":rms,"samples":rows,
            "native_sha256":{name:sha(case/name) for name in sources}}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    with (output/"stirred-vessel-history.csv").open("x",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig,axes=plt.subplots(1,3,figsize=(12.5,4.2),gridspec_kw={"width_ratios":[1,1,1.25]},layout="constrained")
    for ax,data,title,unit in (
        (axes[0],np.linalg.norm(u,axis=1),"Velocity magnitude","Speed (m/s)"),
        (axes[1],p,"Kinematic pressure","p (m²/s²)")):
        pc=PolyCollection(1000*polygons,array=data,cmap="viridis",edgecolors="none",antialiaseds=False)
        ax.add_collection(pc)
        ax.set(xlim=(-105,105),ylim=(-105,105),aspect="equal",xlabel="x (mm)",ylabel="y (mm)",title=title)
        fig.colorbar(pc,ax=ax,label=unit,shrink=.78)
    for key,label in (("rotor_mean_tangential_speed_m_s","Rotor zone"),("stator_mean_tangential_speed_m_s","Stator zone")):
        axes[2].plot(ITERATIONS,[row[key] for row in rows],"o-",label=label)
    axes[2].set(xlabel="Iteration",ylabel="Mean tangential speed (m/s)",title="Signed zone averages")
    axes[2].grid(alpha=.2)
    axes[2].legend()
    fig.suptitle("Stirred vessel in a rotating frame | OpenFOAM",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(output/f"stirred-vessel.{ext}",dpi=180)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","foam_reader","output"):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    audit(args.case,args.foam_reader,args.output)
