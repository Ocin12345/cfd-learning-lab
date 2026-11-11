
"""Conservative periodic finite-volume transport of a Gaussian pulse."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SPEED = 1.0
INITIAL_CENTER = .25
INITIAL_WIDTH = .06
QUAD, WEIGHTS = np.polynomial.legendre.leggauss(8)


def wrapped_distance(x, center):
    return (x - center + .5) % 1 - .5


def reference(x, dx, t, diffusivity):
    sigma = np.sqrt(INITIAL_WIDTH**2 + 2*diffusivity*t)
    center = (INITIAL_CENTER + SPEED*t) % 1
    points = wrapped_distance(x[:,None]+dx*QUAD/2, center)
    total = sum(np.exp(-(points+shift)**2/(2*sigma**2)) for shift in (-2,-1,0,1,2))
    return total @ (WEIGHTS/2) * INITIAL_WIDTH/sigma


def slope(q):
    left = q - np.roll(q,1)
    right = np.roll(q,-1) - q
    return np.where(left*right > 0, np.sign(left)*np.minimum(np.minimum(2*abs(left),2*abs(right)),abs(left+right)/2),0)


def solve(n, method, diffusivity, end):
    dx = 1/n
    x = (np.arange(n)+.5)*dx
    q = reference(x,dx,0,diffusivity)
    mass = q.sum()*dx
    snapshots = 20
    per_snapshot = int(np.ceil(end/snapshots / (.35/(SPEED/dx+2*diffusivity/dx**2))))
    dt = end / (snapshots*per_snapshot)

    def rhs(u):
        face = u + slope(u)/2 if method=="muscl" else u
        flux = SPEED*face - diffusivity*(np.roll(u,-1)-u)/dx
        return -(flux-np.roll(flux,1))/dx

    history = []
    fields = []
    max_mass_error = 0.
    for k in range(snapshots*per_snapshot+1):
        t = k*dt
        error = abs(q.sum()*dx-mass)/mass
        max_mass_error = max(max_mass_error,float(error))
        assert np.isfinite(q).all() and q.min()>=-1e-12 and q.max()<=1+1e-12
        if k % per_snapshot == 0:
            ref = reference(x,dx,t,diffusivity)
            center = (INITIAL_CENTER+SPEED*t)%1
            r2 = wrapped_distance(x,center)**2
            variance = np.sum(q*r2)*dx/mass
            ref_variance = np.sum(ref*r2)*dx/(np.sum(ref)*dx)
            history.append((t,float(q.max()),float(ref.max()),float(variance),float(ref_variance),
                            float(np.sum(abs(q-ref))*dx/mass),float(error)))
            fields.append((t,q.copy(),ref))
        if k < snapshots*per_snapshot:
            first = q + dt*rhs(q)
            second = .75*q + .25*(first+dt*rhs(first))
            q = q/3+2/3*(second+dt*rhs(second))
    assert max_mass_error < 1e-11
    return {
        "cells":n,"method":method,"time_step":dt,"steps":snapshots*per_snapshot,
        "maximum_relative_mass_error_all_steps":max_mass_error,
        "final_peak":history[-1][1],"final_exact_peak":history[-1][2],
        "final_variance_m2":history[-1][3],"final_exact_variance_m2":history[-1][4],
        "final_normalized_l1_error":history[-1][5],
        "maximum_saved_normalized_l1_error":max(row[5] for row in history),
    },x,history,fields


def main(mode,out):
    diffusivity,end = (0.,1.) if mode==41 else (.002,.5)
    specs=[(160,"upwind"),(160,"muscl")] if mode==41 else [(160,"muscl"),(320,"muscl")]
    data=[solve(n,method,diffusivity,end) for n,method in specs]
    a,b=(item[0] for item in data)
    assert b["final_normalized_l1_error"] < a["final_normalized_l1_error"]
    if mode==41:
        assert b["final_normalized_l1_error"] < .04 and b["final_peak"] > a["final_peak"]
    else:
        assert max(case[0]["maximum_saved_normalized_l1_error"] for case in data)<.03
        assert abs(b["final_peak"]-b["final_exact_peak"])<.01
        assert abs(b["final_variance_m2"]-b["final_exact_variance_m2"])<.00015
    out.mkdir()
    checks={"study":mode,"velocity_m_s":SPEED,"diffusivity_m2_s":diffusivity,"end_time_s":end,"cases":[a,b]}
    (out/"checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    with (out/"history.csv").open("x",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(("cells","method","time_s","peak","reference_peak","variance_m2","reference_variance_m2","l1_over_initial_mass","mass_error"))
        for (info,_,rows,_) in data:
            writer.writerows((info["cells"],info["method"],*row) for row in rows)
    with (out/"profiles.csv").open("x",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(("cells","method","time_s","cell_center_m","computed","reference"))
        for info,x,_,fields in data:
            for t,q,ref in fields:
                writer.writerows((info["cells"],info["method"],t,*row) for row in zip(x,q,ref))
    fig,(profile,trace)=plt.subplots(1,2,figsize=(10.6,4.2),layout="constrained")
    fig.suptitle("One periodic crossing" if mode==41 else "Periodic advection with molecular diffusion",fontweight="bold")
    for info,x,rows,fields in data:
        label=info["method"].upper() if mode==41 else f'{info["cells"]} cells'
        profile.plot(x,fields[-1][1],label=label)
        rows=np.array(rows)
        trace.plot(rows[:,0],rows[:,1] if mode==41 else rows[:,3],label=label)
    _,x,rows,fields=data[-1]
    profile.plot(x,fields[-1][2],"k--",label="Exact cell average")
    rows=np.array(rows)
    trace.plot(rows[:,0],rows[:,2] if mode==41 else rows[:,4],"k--",label="Exact")
    profile.set(xlabel="x (m)",ylabel="Scalar concentration",title=f"Computed and exact profiles at {end:g} s")
    trace.set(xlabel="Time (s)",ylabel="Peak concentration" if mode==41 else "Variance (m²)",title="Saved scalar history")
    for ax in (profile,trace):
        ax.grid(alpha=.24)
        ax.legend(frameon=False)
    fig.savefig(out/"periodic-transport.png",dpi=180)
    fig.savefig(out/"periodic-transport.pdf")
    print(json.dumps(checks,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study",type=int,choices=(41,42))
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    main(args.study,args.out)
