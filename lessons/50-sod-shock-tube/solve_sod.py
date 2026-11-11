
"""Conservative HLL finite-volume integration of the Sod shock tube."""
import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

GAMMA=1.4
LEFT=(1.,0.,1.)
RIGHT=(.125,0.,.1)
END=.2
QUAD, WEIGHTS=np.polynomial.legendre.leggauss(8)


def conserved(rho,u,p):
    return np.stack((rho,rho*u,p/(GAMMA-1)+.5*rho*u*u),axis=-1)


def primitives(Q):
    rho=Q[...,0]
    u=Q[...,1]/rho
    p=(GAMMA-1)*(Q[...,2]-.5*rho*u*u)
    if not np.isfinite(Q).all() or (rho<=0).any() or (p<=0).any():
        raise ValueError("Nonpositive or nonfinite Euler state")
    return rho,u,p


def sound(rho,p):
    return np.sqrt(GAMMA*p/rho)


def physical_flux(Q,rho,u,p):
    return np.column_stack((rho*u,rho*u*u+p,u*(Q[:,2]+p)))


def hll(Q):
    extended=np.vstack((Q[0],Q,Q[-1]))
    a,b=extended[:-1],extended[1:]
    rho_a,u_a,p_a=primitives(a)
    rho_b,u_b,p_b=primitives(b)
    left=np.minimum(np.minimum(u_a-sound(rho_a,p_a),u_b-sound(rho_b,p_b)),0)
    right=np.maximum(np.maximum(u_a+sound(rho_a,p_a),u_b+sound(rho_b,p_b)),0)
    fa=physical_flux(a,rho_a,u_a,p_a)
    fb=physical_flux(b,rho_b,u_b,p_b)
    l,r=left[:,None],right[:,None]
    return (r*fa-l*fb+l*r*(b-a))/(r-l)


def wave_function(p,state):
    rho,_,base=state
    if p>base:
        return (p-base)*np.sqrt(2/((GAMMA+1)*rho)/(p+(GAMMA-1)/(GAMMA+1)*base))
    return 2*sound(rho,base)/(GAMMA-1)*((p/base)**((GAMMA-1)/(2*GAMMA))-1)


def star_state():
    lo,hi=RIGHT[2],LEFT[2]
    for _ in range(90):
        mid=(lo+hi)/2
        if wave_function(mid,LEFT)+wave_function(mid,RIGHT)>0:
            hi=mid
        else:
            lo=mid
    p=(lo+hi)/2
    u=(wave_function(p,RIGHT)-wave_function(p,LEFT))/2
    a_left=sound(LEFT[0],LEFT[2])
    rho_left=LEFT[0]*(p/LEFT[2])**(1/GAMMA)
    a_star=sound(rho_left,p)
    ratio=(GAMMA-1)/(GAMMA+1)
    rho_right=RIGHT[0]*(p/RIGHT[2]+ratio)/(ratio*p/RIGHT[2]+1)
    shock=sound(RIGHT[0],RIGHT[2])*np.sqrt((GAMMA+1)/(2*GAMMA)*p/RIGHT[2]+(GAMMA-1)/(2*GAMMA))
    return p,u,rho_left,rho_right,(-a_left,u-a_star,u,shock)


STAR=star_state()


def exact(x,t):
    pstar,ustar,rl,rr,(head,tail,contact,shock)=STAR
    if t==0:
        return np.where((x<.5)[...,None],np.array(LEFT),np.array(RIGHT))
    speed=(x-.5)/t
    a=sound(LEFT[0],LEFT[2])
    fan_u=2/(GAMMA+1)*(a+speed)
    fan_sound=2/(GAMMA+1)*(a-(GAMMA-1)/2*speed)
    factor=np.maximum(fan_sound/a,0)
    fan=np.stack((LEFT[0]*factor**(2/(GAMMA-1)),fan_u,LEFT[2]*factor**(2*GAMMA/(GAMMA-1))),axis=-1)
    values=np.broadcast_to(np.array(RIGHT),x.shape+(3,)).copy()
    values[speed<shock]=(rr,ustar,pstar)
    values[speed<contact]=(rl,ustar,pstar)
    values[speed<tail]=fan[speed<tail]
    values[speed<head]=LEFT
    return values


def average_exact(x,dx,t):
    lower=x-dx/2
    upper=x+dx/2
    averaged=np.zeros((len(x),3))
    for wave in sorted([.5+t*v for v in STAR[-1]]+[1.]):
        right=np.maximum(lower,np.minimum(wave,upper))
        points=(lower+right)[:,None]/2+(right-lower)[:,None]*QUAD/2
        states=exact(points,t)
        Q=conserved(states[...,0],states[...,1],states[...,2])
        averaged+=(right-lower)[:,None]/dx*np.einsum("ijk,j->ik",Q,WEIGHTS/2)
        lower=right
    return np.column_stack(primitives(averaged))


def run(n):
    dx=1/n
    x=(np.arange(n)+.5)*dx
    initial=exact(x,0)
    Q=conserved(initial[:,0],initial[:,1],initial[:,2])
    totals=dx*Q.sum(axis=0)
    integrated_boundary=np.zeros(3)
    max_ledger_error=0.
    min_rho=float(Q[:,0].min())
    min_p=float(initial[:,2].min())
    t=0.
    step=0
    saved=[]
    for target in np.arange(1,11)*.02:
        while t < target-1e-14:
            rho,u,p=primitives(Q)
            dt=min(.45*dx/np.max(abs(u)+sound(rho,p)),target-t)
            flux=hll(Q)
            Q=Q-dt/dx*np.diff(flux,axis=0)
            rho,u,p=primitives(Q)
            min_rho=min(min_rho,float(rho.min()))
            min_p=min(min_p,float(p.min()))
            integrated_boundary+=dt*(flux[0]-flux[-1])
            error=abs(dx*Q.sum(axis=0)-totals-integrated_boundary).max()
            max_ledger_error=max(max_ledger_error,float(error))
            t+=dt
            step+=1
        ref=average_exact(x,dx,target)
        num=np.column_stack(primitives(Q))
        saved.append((target,*np.mean(abs(num-ref),axis=0)))
    pstar,ustar,_,_,waves=STAR
    locations=.5+END*np.array(waves)
    xf=(x[:-1]+x[1:])/2
    estimates=[]
    for col,center,width in ((0,locations[2],.055),(2,locations[3],.07)):
        allowed=np.abs(xf-center)<width
        difference=np.abs(np.diff(num[:,col]))
        estimates.append(float(xf[allowed][np.argmax(difference[allowed])]))
    plateau=(x>locations[1]+.04)&(x<locations[2]-.04)
    mean_p=float(num[plateau,2].mean())
    assert len(saved)==10 and abs(t-END)<1e-12
    assert max_ledger_error<1e-11
    assert abs((dx*Q.sum(axis=0)-totals)[1]-(LEFT[2]-RIGHT[2])*END)<1e-10
    assert saved[-1][1]<.035 and abs(mean_p-pstar)<.015
    assert abs(estimates[0]-locations[2])<.022 and abs(estimates[1]-locations[3])<.012
    summary={"cells":n,"actual_steps":step,"density_l1_at_endpoint":saved[-1][1],
             "velocity_l1_at_endpoint":saved[-1][2],"pressure_l1_at_endpoint":saved[-1][3],
             "max_conservation_ledger_error":max_ledger_error,"minimum_density":min_rho,"minimum_pressure":min_p,
             "mean_left_star_pressure":mean_p,"contact_position_m":estimates[0],"shock_position_m":estimates[1]}
    return summary,x,num,ref,saved


if __name__=="__main__":
    out=Path(sys.argv[1])
    results=[run(n) for n in (200,400)]
    assert results[1][0]["density_l1_at_endpoint"]<results[0][0]["density_l1_at_endpoint"]
    out.mkdir()
    checks={"end_time_s":END,"reference_star_pressure":STAR[0],"reference_star_velocity":STAR[1],
            "exact_wave_positions_m":list(.5+END*np.array(STAR[-1])),"cases":[d[0] for d in results]}
    (out/"checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    with (out/"profiles.csv").open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("cells","x_m","density","velocity","pressure","exact_average_density","exact_velocity_from_conserved_average","exact_pressure_from_conserved_average"))
        for summary,x,num,ref,_ in results:
            writer.writerows((summary["cells"],a,*b,*c) for a,b,c in zip(x,num,ref))
    with (out/"errors.csv").open("x",newline="") as f:
        writer=csv.writer(f)
        writer.writerow(("cells","time_s","density_l1","velocity_l1","pressure_l1"))
        for summary,*_,rows in results:
            writer.writerows((summary["cells"],*row) for row in rows)
    fig,axes=plt.subplots(1,3,figsize=(13,4.2),layout="constrained")
    fig.suptitle("Sod shock tube | Euler solution at 0.20 s",fontweight="bold")
    exact_x=np.linspace(0,1,3001)
    analytical=exact(exact_x,END)
    for j,(axis,label) in enumerate(zip(axes,("Density","Velocity","Pressure"))):
        for summary,x,num,*_ in results:
            axis.plot(x,num[:,j],label=f"{summary['cells']} cells")
        axis.plot(exact_x,analytical[:,j],"k--",label="Exact Riemann solution")
        axis.set(xlabel="x (m)",ylabel=label)
        axis.grid(alpha=.22)
        axis.legend(frameon=False,fontsize=8)
    fig.savefig(out/"shock-tube.png",dpi=180)
    fig.savefig(out/"shock-tube.pdf")
    print(json.dumps(checks,indent=2))
