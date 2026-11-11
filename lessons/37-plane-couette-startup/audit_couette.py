
"""Compare saved OpenFOAM startup profiles with cell averages of the analytical series."""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np

WALL_SPEED = 0.1
GAP = 0.01
VISCOSITY = 0.0001
SPECS = {
    "smoke": (80, 0.00025, 0.02),
    "coarse": (40, 0.0005, 0.4),
    "fine": (80, 0.00025, 0.4),
}
NUMBER = r"[-+0-9.eE]+"


def require(condition, explanation):
    if not condition:
        raise ValueError(explanation)


def vectors(path, count):
    text = path.read_text(encoding="ascii")
    match = re.search(r"internalField\s+uniform\s+\(([^()]*)\)", text)
    if match:
        values = np.repeat(np.fromstring(match[1], sep=" ")[None, :], count, axis=0)
    else:
        match = re.search(
            r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;",
            text, re.S,
        )
        require(match and int(match[1]) == count, f"Invalid vector field {path}")
        values = np.fromstring(match[2].replace("(", " ").replace(")", " "), sep=" ").reshape(-1, 3)
    require(values.shape == (count, 3) and np.isfinite(values).all(), f"Invalid data {path}")
    return values


def analytical(y, time, width):
    """Integrate each sine over the native uniform cell centered at y."""
    if time == 0:
        return np.zeros_like(y)
    modes = np.arange(1, 201, dtype=float)
    factor = np.exp(-(modes * np.pi / GAP)**2 * VISCOSITY * time) / modes
    cell_average = np.sinc(modes * width / (2 * GAP))
    series = np.sin(np.outer(y / GAP, modes * np.pi)) @ (factor * cell_average)
    return WALL_SPEED * (1 - y / GAP - 2 / np.pi * series)


def audit(root, label, manifest):
    ny, dt, end = SPECS[label]
    case = root / label
    for name, expected in manifest[label].items():
        require(hashlib.sha256((case / name).read_bytes()).hexdigest() == expected, f"Input hash {label}/{name}")
    mesh = (case / "check-mesh.log").read_text()
    log = (case / "solver.log").read_text()
    require("Mesh OK." in mesh and re.search(rf"^\s*cells:\s*{8*ny}\s*$",mesh,re.M), f"Mesh {label}")
    require(re.search(r"^End\s*$", log, re.M), f"Missing solver end {label}")
    dt_count = round(end / dt)
    clock = np.array([float(t) for t in re.findall(r"^Time = ("+NUMBER+r")s?\s*$",log,re.M)])
    require(len(clock) == dt_count and np.allclose(clock, dt*np.arange(1,dt_count+1), rtol=0, atol=1e-10), f"Native clock {label}")
    continuity = np.array([float(x) for x in re.findall(r"sum local = ("+NUMBER+r")",log)])
    require(len(continuity) == 2*dt_count and np.isfinite(continuity).all() and continuity.max() <= 1e-9, f"Continuity {label}")
    rows = re.findall(r"Solving for (Ux|Uy|p),[^\r\n]*Final residual = ("+NUMBER+r")", log)
    require([field for field, _ in rows] == ["Ux", "Uy", "p", "p"] * dt_count, f"Linear solve sequence {label}")
    residuals = np.array([float(value) for _, value in rows])
    require(len(residuals) == 4*dt_count and np.isfinite(residuals).all() and residuals.max() <= 1e-12, f"Residual {label}")
    saved = sorted((float(p.name),p) for p in case.iterdir() if p.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?",p.name))
    require(len(saved)==round(end/.01)+1 and np.allclose([t for t,_ in saved], np.arange(len(saved))*.01, rtol=0, atol=1e-10), f"Saved times {label}")
    centres = vectors(case / "0/C",8*ny)
    by_row = np.lexsort((centres[:,0],centres[:,1].round(12))).reshape(ny,8)
    y = (np.arange(ny)+.5)*GAP/ny
    x = (np.arange(8)+.5)*.008/8
    require(np.allclose(centres[by_row,1],y[:,None],rtol=0,atol=1e-12)
            and np.allclose(centres[by_row,0],x[None,:],rtol=0,atol=1e-12)
            and np.allclose(centres[:,2],.0005,rtol=0,atol=1e-12), f"Native centres {label}")
    profiles, errors = [], []
    for t, folder in saved:
        uvw = vectors(folder / "U",8*ny)
        u = uvw[by_row,0]
        transverse = float(np.abs(uvw[:,1:]).max())
        spread = float(np.ptp(u,axis=1).max())
        require(transverse <= 1e-9 and spread <= 1e-9, f"Symmetry {label} t={t}")
        require(u.min() >= -1e-9 and u.max() <= WALL_SPEED + 1e-9, f"Velocity bounds {label} t={t}")
        values = u.mean(axis=1)
        reference = analytical(y,t,GAP/ny)
        difference = (values-reference)/WALL_SPEED
        rms = float(np.sqrt(np.mean(difference**2)))
        largest = float(np.abs(difference).max())
        require(rms < .006 and largest < .02, f"Reference errors {label} t={t}: {rms} RMS, {largest} max")
        errors.append((t,rms,largest,transverse,spread))
        profiles.extend(zip(np.full(ny,t),y,values,reference))
    return {
        "label":label, "cell_count":8*ny, "time_step_s":dt, "end_time_s":end,
        "actual_steps":dt_count, "saved_states":len(saved),
        "maximum_normalized_rms_error":max(row[1] for row in errors),
        "maximum_normalized_absolute_error":max(row[2] for row in errors),
        "maximum_transverse_velocity_m_s":max(row[3] for row in errors),
        "maximum_within_row_velocity_spread_m_s":max(row[4] for row in errors),
        "maximum_reported_local_continuity_error":float(continuity.max()),
        "maximum_reported_final_linear_residual":float(residuals.max()),
        "source_log_sha256":hashlib.sha256((case/"solver.log").read_bytes()).hexdigest(),
    },profiles,errors


def export(destination, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    destination.mkdir()
    with (destination / "checks.json").open("x",encoding="utf-8") as out:
        json.dump({r[0]["label"]:r[0] for r in results},out,indent=2)
    with (destination / "profiles.csv").open("x",newline="",encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(("case","time_s","y_m","u_openfoam_m_s","u_reference_cell_average_m_s"))
        for summary,profiles,_ in results:
            writer.writerows((summary["label"],*row) for row in profiles)
    with (destination / "errors.csv").open("x",newline="",encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(("case","time_s","rms_over_wall_speed","max_abs_over_wall_speed","max_transverse_m_s","max_within_row_spread_m_s"))
        for summary,_,errors in results:
            writer.writerows((summary["label"],*row) for row in errors)

    fig, (velocity, error) = plt.subplots(1,2,figsize=(11,4.5),layout="constrained")
    fig.suptitle("Plane Couette startup | OpenFOAM Foundation 14",fontweight="bold")
    times=(.01,.05,.15,.4)
    fine=np.array(results[-1][1])
    for t,color in zip(times,("#6b3fa0","#176c91","#007e71","#b34e24")):
        rows=fine[np.isclose(fine[:,0],t)]
        velocity.plot(rows[:,1]*1000,rows[:,3],color=color,label=f"{t:g} s")
        velocity.plot(rows[::6,1]*1000,rows[::6,2],"o",color=color,markersize=3,fillstyle="none")
    velocity.set(xlabel="Distance from moving plate (mm)",ylabel="Streamwise velocity (m/s)",title="Reference lines and fine-grid samples")
    velocity.legend(frameon=False)
    for summary,_,errors in results:
        values=np.array(errors)[1:]
        error.semilogy(values[:,0],values[:,1]*100,label=f"{summary['label'].capitalize()} | {summary['cell_count']} cells")
    error.set(xlabel="Time (s)",ylabel="Profile RMS error (% of plate speed)",title="Error at every saved time after startup")
    error.legend(frameon=False)
    for axis in (velocity,error):
        axis.grid(alpha=.22)
    fig.savefig(destination/"couette-profiles-and-errors.png",dpi=185)
    fig.savefig(destination/"couette-profiles-and-errors.pdf")
    plt.close(fig)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("manifest",type=Path)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--smoke",action="store_true")
    mode.add_argument("--export",type=Path)
    args=parser.parse_args()
    manifest=json.loads(args.manifest.read_bytes())
    if args.smoke:
        print(json.dumps(audit(args.root,"smoke",manifest)[0],indent=2))
    else:
        result=[audit(args.root,label,manifest) for label in ("coarse","fine")]
        require(result[1][0]["maximum_normalized_rms_error"] < result[0][0]["maximum_normalized_rms_error"], "Fine case must reduce maximum saved RMS error")
        export(args.export,result)
        print(json.dumps([r[0] for r in result],indent=2))
