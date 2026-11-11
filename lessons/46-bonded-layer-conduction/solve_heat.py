
"""Three finite-volume heat-transfer benchmarks."""
import argparse
import csv
import json
from math import erfc
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

GAUSS, WEIGHTS = np.polynomial.legendre.leggauss(8)


def matrix(n, dx, coefficient):
    """Uniform cell-centered diffusion with fixed temperatures at both faces."""
    A = -2*np.eye(n)+np.eye(n,k=1)+np.eye(n,k=-1)
    A[0,0] = A[-1,-1] = -3
    return coefficient/dx**2 * A


def graetz(n, dxial=None):
    dy = 2/n
    y = -1 + (np.arange(n)+.5)*dy
    u = 1.5*(1-y**2-dy**2/12)
    A = matrix(n,dy,1.)
    scale = np.sqrt(u)
    values,V = np.linalg.eigh(A/scale[:,None]/scale[None,:])
    start = V.T@scale
    sites = np.linspace(0,.5,101)
    if dxial is None:
        modes = np.exp(sites[:,None]*values)*start
        fields = modes@V.T/scale
        balance = 0.
    else:
        steps=round(.5/dxial)
        factor=1/(1-dxial*values)
        current=start.copy()
        fields=[np.ones(n)]
        balance=0.
        old=fields[0]
        for k in range(1,steps+1):
            current*=factor
            new=V@current/scale
            wall_flux=2*(new[0]+new[-1])/dy
            error=abs(np.dot(u,new-old)*dy/dxial+wall_flux)
            balance=max(balance,error)
            assert np.isfinite(new).all() and new.min()>=-1e-10 and new.max()<=1+1e-10
            if k%(steps//100)==0:
                fields.append(new.copy())
            old=new
        fields=np.array(fields)
        assert len(fields)==101 and balance<1e-8
    bulk=fields@u/u.sum()
    flux=(fields[:,0]+fields[:,-1])/dy
    return y,sites,fields,bulk,4*flux/bulk,balance


def heat_entrance(out,profile,history):
    ref_lo=graetz(256)
    ref=graetz(512)
    use=ref[1]>=.025
    assert abs(ref_lo[3]-ref[3]).max()<.0005 and abs(ref_lo[4][use]/ref[4][use]-1).max()<.006
    checks=[]
    for n,dt in ((64,.001),(128,.0005)):
        y,s,fields,bulk,nu,balance=graetz(n,dt)
        err=float(abs(bulk-ref[3]).max())
        nu_err=float(abs(nu[use]/ref[4][use]-1).max())
        assert err<.005 and nu_err<.05 and abs(nu[-1]-7.541)<.12
        checks.append(dict(cross_stream_cells=n,axial_steps=round(.5/dt),maximum_bulk_temperature_error=err,
                           maximum_relative_nusselt_error_for_s_ge_0p025=nu_err,exit_nusselt=float(nu[-1]),
                           maximum_discrete_heat_balance_residual=float(balance)))
        history.plot(s[1:],bulk[1:],label=f"{n} cells")
        for k in (5,20,60,100):
            profile.plot(y,fields[k],label=f"s = {s[k]:.3g}" if n==128 else None,alpha=1 if n==128 else .3)
        save(out/f"axial-{n}.csv",("s","bulk_dimensionless_temperature","nusselt","reference_bulk","reference_nusselt"),
             zip(s,bulk,nu,ref[3],ref[4]))
    assert checks[1]["maximum_bulk_temperature_error"] < checks[0]["maximum_bulk_temperature_error"]
    history.plot(ref[1][1:],ref[3][1:],"k--",label="Fine reference")
    profile.set(xlabel="y / half-gap",ylabel="Dimensionless temperature",title="Transverse profiles as the fluid warms")
    history.set(xlabel="Axial coordinate s",ylabel="Bulk dimensionless temperature",title="Bulk temperature versus marching distance")
    return {"study":45,"results":checks,"reference_cells":512,"reference_exit_nusselt":float(ref[4][-1])}


def composite(out,profile,history):
    length=.04
    exact_flux=20/(.02/10+.02/1)
    interface=320-exact_flux*.02/10
    checks=[]
    for n in (20,40):
        dx=length/n
        x=(np.arange(n)+.5)*dx
        k=np.where(x<.02,10.,1.)
        reference=np.where(x<.02,320-exact_flux*x/10,300+exact_flux*(length-x))
        for scheme in ("series resistance","arithmetic conductivity"):
            face_k=2/(1/k[1:]+1/k[:-1]) if scheme=="series resistance" else (k[1:]+k[:-1])/2
            g=np.r_[2*k[0],face_k,2*k[-1]]/dx
            A=np.diag(g[:-1]+g[1:])-np.diag(g[1:-1],1)-np.diag(g[1:-1],-1)
            b=np.zeros(n)
            b[0]=320*g[0]
            b[-1]=300*g[-1]
            T=np.linalg.solve(A,b)
            flux=g*(np.r_[320,T]-np.r_[T,300])
            left=T[n//2-1]-flux[n//2]*dx/(2*k[n//2-1])
            right=T[n//2]+flux[n//2]*dx/(2*k[n//2])
            assert abs(np.ptp(flux))<1e-7 and np.isfinite(T).all() and np.min(T)>=300 and np.max(T)<=320
            if scheme=="series resistance":
                assert abs(flux.mean()-exact_flux)<1e-7 and np.max(abs(T-reference))<1e-9
                assert abs(left-interface)<1e-9 and abs(right-interface)<1e-9
            checks.append(dict(cells=n,face_model=scheme,flux_W_m2=float(flux.mean()),
                               relative_flux_error=float(abs(flux.mean()/exact_flux-1)),
                               max_temperature_error_K=float(abs(T-reference).max()),
                               extrapolated_interface_temperature_jump_K=float(left-right)))
            save(out/f"{n}-{scheme.replace(' ','-')}.csv",("x_m","temperature_K","exact_temperature_K"),zip(x,T,reference))
            if n==20:
                profile.plot(x*1000,T,"o-",label=scheme.capitalize(),markersize=3)
    assert checks[3]["relative_flux_error"]<checks[1]["relative_flux_error"]
    profile.plot([0,20,40],[320,interface,300],"k--",label="Exact")
    history.bar(["20 cells","40 cells"],[checks[i]["relative_flux_error"]*100 for i in (1,3)],color="#b34e24",label="Arithmetic")
    history.plot([0,1],[0,0],"o",color="#176c91",label="Series resistance")
    profile.axvline(20,color=".6",ls=":")
    profile.set(xlabel="x (mm)",ylabel="Temperature (K)",title="Temperature across the bonded layers")
    history.set(ylabel="Heat-flux error (%)",title="Effect of the interface face model")
    return {"study":46,"exact_heat_flux_W_m2":exact_flux,"exact_interface_temperature_K":interface,"results":checks}


def half_space(out,profile,history):
    alpha=1e-5
    L=.05
    k=10.
    rise=20.
    end=10.
    checks=[]
    for n,dt in ((80,.02),(160,.01)):
        dx=L/n
        x=(np.arange(n)+.5)*dx
        A=matrix(n,dx,alpha)
        lam,V=np.linalg.eigh(A)
        steady=1-x/L
        modes=V.T@(-steady)
        gain=1/(1-dt*lam)
        old=np.zeros(n)
        max_balance=0.
        error=[]
        rows=[]
        flux_rows=[]
        for j in range(1,round(end/dt)+1):
            t=j*dt
            modes*=gain
            T=steady+V@modes
            qL=2*(1-T[0])/dx
            qR=2*T[-1]/dx
            max_balance=max(max_balance,abs((T-old).sum()*dx/dt-alpha*(qL-qR)))
            assert np.isfinite(T).all() and T.min()>=-1e-10 and T.max()<=1+1e-10
            if j%round(.2/dt)==0:
                sample=x[:,None]+dx*GAUSS/2
                ref=np.vectorize(erfc)(sample/(2*np.sqrt(alpha*t)))@(WEIGHTS/2)
                flux_ref=k*rise/np.sqrt(np.pi*alpha*t)
                error.append((abs(T-ref).max(),abs(k*rise*qL/flux_ref-1)))
                rows.extend((t,*row) for row in zip(x,300+rise*T,300+rise*ref))
                flux_rows.append((t,k*rise*qL,flux_ref))
                if n==160 and any(np.isclose(t,s) for s in (.2,1,4,10)):
                    line,=profile.plot(x*1000,300+rise*T,label=f"{t:g} s")
                    profile.plot(x[::8]*1000,300+rise*ref[::8],"o",fillstyle="none",markersize=3,color=line.get_color())
            old=T
        errors=np.array(error)
        assert max_balance<1e-9 and errors[:,0].max()<.02 and errors[:,1].max()<.04
        checks.append(dict(cells=n,steps=round(end/dt),maximum_temperature_error_over_step=float(errors[:,0].max()),
                           maximum_relative_hot_wall_flux_error=float(errors[:,1].max()),
                           maximum_discrete_balance_residual=float(max_balance),
                           final_hot_wall_flux_W_m2=float(flux_rows[-1][1])))
        save(out/f"temperature-{n}.csv",("time_s","x_m","computed_K","halfspace_K"),rows)
        save(out/f"heat-flux-{n}.csv",("time_s","computed_W_m2","halfspace_W_m2"),flux_rows)
        history.loglog(*np.array(flux_rows).T[:2],label=f"{n} cells")
    assert checks[1]["maximum_relative_hot_wall_flux_error"]<checks[0]["maximum_relative_hot_wall_flux_error"]
    rows=np.array(flux_rows)
    history.loglog(rows[:,0],rows[:,2],"k--",label="Half-space solution")
    profile.set(xlabel="Depth (mm)",ylabel="Temperature (K)",title="Fine-grid profiles and half-space samples")
    history.set(xlabel="Time (s)",ylabel="Hot-wall flux (W/m²)",title="Heat flux after the surface step")
    return {"study":47,"results":checks,"halfspace_temperature_rise_fraction_at_far_boundary_final":erfc(L/(2*np.sqrt(alpha*end)))}


def save(path,header,rows):
    with path.open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study",type=int,choices=(45,46,47))
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    args.out.mkdir()
    fig,(profile,history)=plt.subplots(1,2,figsize=(11,4.3),layout="constrained")
    methods={
        45:(heat_entrance,"Thermal development between isothermal plates"),
        46:(composite,"Steady heat conduction through two bonded layers"),
        47:(half_space,"Transient heat penetration into a solid"),
    }
    solve,title=methods[args.study]
    checks=solve(args.out,profile,history)
    fig.suptitle(title,fontweight="bold")
    for axis in (profile,history):
        axis.grid(alpha=.22)
        axis.legend(frameon=False)
    fig.savefig(args.out/"heat-transfer.png",dpi=180)
    fig.savefig(args.out/"heat-transfer.pdf")
    (args.out/"checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    print(json.dumps(checks,indent=2))
