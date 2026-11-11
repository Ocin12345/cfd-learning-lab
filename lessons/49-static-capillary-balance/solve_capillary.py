
"""Evolve a fixed-interface capillary balance test with a staggered pressure projection."""
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LENGTH=.04
RADIUS=.008
SIGMA=.072
DENSITY=1000.
NU=1e-5
JUMP=SIGMA/RADIUS
DT=.001
STEPS=20


def run(n,method):
    dx=LENGTH/n
    x=(np.arange(n)+.5)*dx-LENGTH/2
    xx,yy=np.meshgrid(x,x)
    distance=np.hypot(xx,yy)
    sdf=distance-RADIUS
    s=np.clip(sdf/(1.5*dx),-1,1)
    color=.5*(1-s-np.sin(np.pi*s)/np.pi)
    inside=distance<RADIUS-4*dx
    outside=distance>RADIUS+4*dx

    def gradient(f,axis):
        return (np.roll(f,-1,axis)-f)/dx

    def centered(f,axis):
        return (np.roll(f,-1,axis)-np.roll(f,1,axis))/(2*dx)

    def divergence(a,b):
        return (a-np.roll(a,1,1)+b-np.roll(b,1,0))/dx

    def laplacian(f):
        return (sum(np.roll(f,k,i) for i in (0,1) for k in (-1,1))-4*f)/dx**2

    if method=="prescribed":
        curvature=np.full_like(sdf,1/RADIUS)
    else:
        nx,ny=centered(sdf,1),centered(sdf,0)
        norm=np.hypot(nx,ny)
        nx=np.divide(nx,norm,out=np.zeros_like(nx),where=norm>1e-12)
        ny=np.divide(ny,norm,out=np.zeros_like(ny),where=norm>1e-12)
        curvature=centered(nx,1)+centered(ny,0)
    fx=SIGMA/DENSITY * (curvature+np.roll(curvature,-1,1))/2 * gradient(color,1)
    fy=SIGMA/DENSITY * (curvature+np.roll(curvature,-1,0))/2 * gradient(color,0)
    symbol=-4*np.sin(np.pi*np.fft.fftfreq(n))**2/dx**2
    modes=symbol[:,None]+symbol[None,:]
    inverse=np.divide(1,modes,out=np.zeros_like(modes),where=modes!=0)
    pkin=np.fft.ifft2(np.fft.fft2(divergence(fx,fy))*inverse).real
    pressure=(pkin-pkin[outside].mean())*DENSITY
    projected_x=fx-gradient(pkin,1)
    projected_y=fy-gradient(pkin,0)
    jump=float(pressure[inside].mean()-pressure[outside].mean())
    u=np.zeros_like(sdf)
    v=np.zeros_like(sdf)
    assert NU*DT/dx**2<.25
    rows=[]
    for k in range(STEPS+1):
        uc=(u+np.roll(u,1,1))/2
        vc=(v+np.roll(v,1,0))/2
        speed=np.hypot(uc,vc)
        rows.append((k*DT,float(speed.max()),float(abs(divergence(u,v)).max())))
        assert np.isfinite(speed).all()
        if k<STEPS:
            u+=DT*(projected_x+NU*laplacian(u))
            v+=DT*(projected_y+NU*laplacian(v))
    assert max(r[2] for r in rows)<1e-9
    if method=="prescribed":
        assert abs(jump-JUMP)<1e-9 and rows[-1][1]<1e-12
    else:
        assert abs(jump/JUMP-1)<.01 and rows[-1][1]<1e-3
    return {
        "cells_per_side":n,"curvature":method,"pressure_jump_Pa":jump,
        "relative_pressure_jump_error":float(abs(jump/JUMP-1)),
        "final_maximum_speed_m_s":rows[-1][1],
        "maximum_discrete_divergence_s_inverse":max(r[2] for r in rows)
    },rows,x,pressure,speed


if __name__=="__main__":
    out=Path(sys.argv[1])
    data=[run(n,mode) for n in (64,128) for mode in ("prescribed","reconstructed")]
    assert data[3][0]["final_maximum_speed_m_s"] < data[1][0]["final_maximum_speed_m_s"]
    assert data[3][0]["relative_pressure_jump_error"] < data[1][0]["relative_pressure_jump_error"]
    out.mkdir()
    checks={"reference_jump_Pa":JUMP,"steps":STEPS,"time_step_s":DT,"results":[row[0] for row in data]}
    (out/"checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    with (out/"speed-history.csv").open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("cells_per_side","curvature","time_s","max_speed_m_s","max_divergence_s_inverse"))
        for info,rows,*_ in data:
            writer.writerows((info["cells_per_side"],info["curvature"],*row) for row in rows)
    with (out/"fine-centerline.csv").open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("curvature","x_m","pressure_Pa","speed_m_s"))
        for info,_,x,p,speed in data[2:]:
            writer.writerows((info["curvature"],*row) for row in zip(x,p[64],speed[64]))
    fig,(a,b,c)=plt.subplots(1,3,figsize=(12.8,4),layout="constrained")
    fig.suptitle("Capillary pressure balance for a fixed circular interface",fontweight="bold")
    extent=(-20,20,-20,20)
    exact=a.imshow(data[2][3],extent=extent,origin="lower",vmin=0,vmax=JUMP,cmap="viridis")
    speed=b.imshow(data[3][4]*1e6,extent=extent,origin="lower",cmap="magma")
    fig.colorbar(exact,ax=a,label="Pressure (Pa)")
    fig.colorbar(speed,ax=b,label="Speed (µm/s)")
    a.set(xlabel="x (mm)",ylabel="y (mm)",title="128² | Prescribed curvature pressure")
    b.set(xlabel="x (mm)",ylabel="y (mm)",title="128² | Reconstructed curvature speed")
    for info,rows,*_ in data:
        t,m,_=np.array(rows).T
        c.semilogy(t[1:],np.maximum(m[1:],1e-18),label=f"{info['cells_per_side']}² {info['curvature']}")
    c.set(xlabel="Time (s)",ylabel="Peak speed (m/s)",title="Velocity after pressure projection")
    c.legend(frameon=False,fontsize=8)
    c.grid(alpha=.2)
    fig.savefig(out/"capillary-balance.png",dpi=180)
    fig.savefig(out/"capillary-balance.pdf")
    print(json.dumps(checks,indent=2))
