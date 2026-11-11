"""Check native miscible lock-exchange results and plot sludge transport."""
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

CELLS = 8000
TIMES = tuple(range(1,9))
NUMBER = r"[-+0-9.eE]+"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(case, reader_path, output):
    read = runpy.run_path(str(reader_path))["foam_values"]

    def field(name, kind="scalar"):
        return read((case/name).read_text(encoding="ascii"), "internalField", CELLS, kind)

    for line in (case/"input-sha256.txt").read_text(encoding="ascii").splitlines():
        expected, name = line.split("  ",1)
        assert digest(case/name) == expected, name
    mesh = (case/"check-mesh.log").read_text(encoding="ascii")
    log = (case/"solver.log").read_text(encoding="ascii")
    assert "Mesh OK." in mesh and re.search(r"^\s*cells:\s*8000\s*$", mesh, re.M)
    assert re.search(r"^End\s*$", log, re.M)
    steps = np.array(re.findall(r"^Time = ("+NUMBER+r")s$",log,re.M),dtype=float)
    assert len(steps) and np.all(np.diff(steps)>0) and steps[-1] == 8
    x,y,z = field("0/C","vector").T
    vol = field("0/Vc")
    assert np.allclose(vol,.005,rtol=0,atol=1e-10)
    assert np.allclose(z,1,rtol=0,atol=1e-10)
    ix, iy = np.rint(x/.05-.5).astype(int), np.rint(y/.05-.5).astype(int)
    order = np.argsort(iy*200+ix)
    assert np.array_equal((iy*200+ix)[order], np.arange(CELLS))
    assert np.allclose(x,(ix+.5)*.05,rtol=0,atol=1e-8)
    assert np.allclose(y,(iy+.5)*.05,rtol=0,atol=1e-8)
    initial = vol[x<5].sum()*1000
    assert abs(initial-20000) < 1e-6
    rows, states = [], {}
    sources = ["input-sha256.txt","check-mesh.log","solver.log","0/C","0/Vc"]
    for t in TIMES:
        clock = (case/str(t)/"uniform/time").read_text(encoding="ascii")
        step = int(re.search(r"^index\s+(\d+);",clock,re.M)[1])
        assert abs(float(re.search(r"^value\s+([^;]+);",clock,re.M)[1])-t)<1e-8
        assert abs(steps[step-1]-t)<1e-8
        s,w,rho = (field(f"{t}/{name}") for name in ("sludge","water","rho"))
        field(f"{t}/U","vector")
        field(f"{t}/p")
        field(f"{t}/p_rgh")
        assert np.all(field(f"{t}/T")>0)
        assert min(s.min(),w.min()) >= -1e-6 and max(s.max(),w.max()) <= 1+1e-6
        closure = np.max(np.abs(s+w-1))
        assert closure <= 2e-6, (t,closure)
        mixing_error = np.max(np.abs(rho-1/(s/1000+w/990)))
        assert mixing_error <= .005, (t,mixing_error)
        mass = rho*s*vol
        mass_error = abs(mass.sum()/initial-1)
        assert mass_error < .002, (t,mass_error)
        grid = s[order].reshape(40,200)
        lower = grid[:8,100:].mean(axis=0)
        misses = np.flatnonzero(lower <= .5)
        filled = int(misses[0]) if len(misses) else len(lower)
        front = 5 + .05*(filled-.5) if filled else 5.
        right_mass = mass[x>5].sum()
        rows.append({"time_s":t,"native_step":step,"sludge_mass_kg":float(mass.sum()),
                     "sludge_mass_right_kg":float(right_mass),"fraction_sum_max_error":float(closure),
                     "density_crosscheck_max_error_kg_m3":float(mixing_error),
                     "relative_sludge_mass_drift":float(mass_error),
                     "lower_band_front_m":front,
                     "sludge_min":float(s.min()),"sludge_max":float(s.max()),
                     "water_min":float(w.min()),"water_max":float(w.max())})
        states[t] = grid
        sources.extend(f"{t}/{name}" for name in ("uniform/time","sludge","water","rho","U","T","p","p_rgh"))
    assert rows[-1]["native_step"] == len(steps)
    assert rows[-1]["sludge_mass_right_kg"] >= .05*initial
    assert rows[-1]["lower_band_front_m"]-5 >= .5
    checks = {"cells":CELLS,"native_steps":len(steps),"saved_states":len(rows),
              "initial_sludge_mass_kg":float(initial),"samples":rows,
              "native_sha256":{name:digest(case/name) for name in sources}}
    output.mkdir()
    (output/"checks.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    with (output/"lock-exchange-history.csv").open("x",newline="") as file:
        writer = csv.DictWriter(file,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(2,2,figsize=(11.5,5.8),layout="constrained")
    for ax,t in zip(axes[0],(4,8)):
        plot = ax.imshow(states[t], origin="lower", extent=(0,10,0,2), vmin=0,vmax=1,
                         cmap="viridis",interpolation="nearest",aspect="equal")
        ax.axvline(5,color="#e6e6e6",linestyle="--",linewidth=1)
        ax.set(xlabel="x (m)",ylabel="Elevation (m)",title=f"Sludge mass fraction at {t} s")
    fig.colorbar(plot,ax=list(axes[0]),label="Sludge mass fraction",shrink=.85)
    times=[0]+list(TIMES)
    axes[1,0].plot(times,[0]+[row["sludge_mass_right_kg"] for row in rows],"o-")
    axes[1,0].set(xlabel="Time (s)",ylabel="Sludge mass (kg)",title="Transport into the right half")
    axes[1,1].plot(times,[5]+[row["lower_band_front_m"] for row in rows],"o-")
    axes[1,1].set(xlabel="Time (s)",ylabel="Front x (m)",title="Front from the lowest 0.4 m")
    for ax in axes[1]:
        ax.grid(alpha=.25)
    fig.suptitle("Miscible lock exchange | OpenFOAM",fontweight="bold")
    for ext in ("png","pdf"):
        fig.savefig(output/f"lock-exchange.{ext}",dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("case","foam_reader","output"):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    audit(args.case,args.foam_reader,args.output)
