
"""Steady and harmonically forced radial momentum in a circular pipe."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = .005
NU = .0001
FREQUENCY = 4.
W = 2*np.pi*FREQUENCY
FORCE = .8
GAUSS, WEIGHTS = np.polynomial.legendre.leggauss(8)


def pipe(n):
    faces = np.linspace(0,R,n+1)
    dr=R/n
    r=(faces[:-1]+faces[1:])/2
    volume=(faces[1:]**2-faces[:-1]**2)/2
    L=np.zeros((n,n))
    for i,f in enumerate(faces[1:-1]):
        coeff=f/dr
        L[i,i]-=coeff
        L[i+1,i+1]-=coeff
        L[i,i+1]+=coeff
        L[i+1,i]+=coeff
    L[-1,-1]-=2*R/dr
    return faces,r,volume,NU*L/volume[:,None]


def averages(faces, function):
    sample=(faces[:-1]+faces[1:])[:,None]/2 + np.diff(faces)[:,None]*GAUSS/2
    weight=sample*WEIGHTS
    return np.sum(weight*function(sample),axis=1)/np.sum(weight,axis=1)


def i0(z):
    term=np.ones_like(z,dtype=complex)
    result=term.copy()
    for k in range(1,42):
        term=term*(z/2)**2/k**2
        result+=term
    return result


def harmonic(r):
    z=np.sqrt(1j*W/NU)
    return FORCE/(1j*W)*(1-i0(z*r)/i0(z*R))


def steady(n):
    faces,r,volume,A=pipe(n)
    profiles=[]
    flows=[]
    for g in (.2,.4,.8,1.6):
        u=np.linalg.solve(-A,np.full(n,g))
        ref=averages(faces,lambda x:g/(4*NU)*(R**2-x**2))
        q=2*np.pi*np.dot(volume,u)
        qref=np.pi*R**4*g/(8*NU)
        wall=NU*2*u[-1]/(R/n)
        rms=np.sqrt(np.sum(volume*(u-ref)**2)/volume.sum())/(g*R**2/(4*NU))
        flows.append((g,q,qref,rms,abs(q/qref-1),abs(wall/(g*R/2)-1),abs(A@u+g).max()))
        profiles.extend((g,*row) for row in zip(r,u,ref))
        assert np.isfinite(u).all() and (u>=0).all() and (np.diff(u)<0).all()
    metrics=np.array(flows)
    assert metrics[:,3:6].max() < .002 and metrics[:,6].max() < 1e-10
    assert np.ptp(metrics[:,1]/metrics[:,0])<1e-15
    return {
        "cells":n,"maximum_profile_rms_over_exact_peak":float(metrics[:,3].max()),
        "maximum_relative_flow_error":float(metrics[:,4].max()),
        "maximum_relative_wall_shear_error":float(metrics[:,5].max()),
        "maximum_steady_equation_residual_m_s2":float(metrics[:,6].max()),
        "maximum_acceleration_m_s2":float(metrics[-1,0]),
        "flow_at_maximum_acceleration_m3_s":float(metrics[-1,1]),
        "exact_flow_at_maximum_acceleration_m3_s":float(metrics[-1,2]),
    },flows,profiles


def oscillatory(n,dt):
    faces,r,volume,A=pipe(n)
    weight=np.sqrt(volume)
    lam,basis=np.linalg.eigh(weight[:,None]*A/weight[None,:])
    source=basis.T@weight*FORCE
    decay=(1+dt*lam/2)/(1-dt*lam/2)
    amplitude=dt*source/(1-dt*lam/2)
    reference=averages(faces,harmonic)
    ref_q=2*np.pi*np.dot(volume,reference)
    transformed=np.zeros(n)
    end=1.25
    every=round(.005/dt)
    steps=round(end/dt)
    records=[]
    late=[]
    for k in range(steps+1):
        t=k*dt
        if k%every==0:
            u=basis@transformed/weight
            q=2*np.pi*np.dot(volume,u)
            qref=(ref_q*np.exp(1j*W*t)).imag
            records.append((t,q,qref))
            if t>1+1e-9:
                late.append(u)
            assert np.isfinite(u).all() and np.max(abs(u)) < 2*FORCE/W
        if k<steps:
            transformed=decay*transformed+amplitude*np.sin(W*(t+dt/2))
    final_times=np.array(records)[-50:,0]
    design=np.column_stack((np.sin(W*final_times),np.cos(W*final_times)))
    fitted=np.linalg.lstsq(design,np.array(late),rcond=None)[0]
    measured=fitted[0]+1j*fitted[1]
    measured_q=2*np.pi*np.dot(volume,measured)
    ratio=measured_q/ref_q
    err=float(abs(abs(ratio)-1))
    phase_error=float(abs(np.angle(ratio)))
    profile_error=float(np.max(abs(measured-reference))/(FORCE/W))
    assert len(late)==50 and err < .003 and phase_error < .006 and profile_error < .003
    profiles=list(zip(r,abs(measured),abs(reference),np.angle(measured),np.angle(reference)))
    return {
        "cells":n,"time_step_s":dt,"steps":steps,"flow_amplitude_m3_s":float(abs(measured_q)),
        "exact_flow_amplitude_m3_s":float(abs(ref_q)),
        "flow_phase_rad":float(np.angle(measured_q)),"exact_flow_phase_rad":float(np.angle(ref_q)),
        "relative_amplitude_error":err,"phase_error_rad":phase_error,
        "maximum_complex_profile_error_over_force_by_frequency":profile_error,
    },records,profiles


def main(study,out):
    data = [steady(n) for n in (40,80)] if study==43 else [oscillatory(n,dt) for n,dt in ((40,.001),(80,.0005))]
    metric="maximum_relative_flow_error" if study==43 else "relative_amplitude_error"
    assert data[1][0][metric]<data[0][0][metric]
    out.mkdir()
    checks={"study":study,"radius_m":R,"viscosity_m2_s":NU,"results":[row[0] for row in data]}
    if study==44:
        checks["womersley_number"]=R*np.sqrt(W/NU)
    (out/"checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    headers = (
        ("acceleration_m_s2","flow_m3_s","exact_flow_m3_s","rms_over_exact_peak","relative_flow_error","relative_wall_shear_error","max_residual_m_s2"),
        ("acceleration_m_s2","radius_m","computed_u_m_s","exact_u_m_s")
    ) if study==43 else (
        ("time_s","computed_flow_m3_s","periodic_exact_flow_m3_s"),
        ("radius_m","computed_velocity_amplitude_m_s","exact_velocity_amplitude_m_s","computed_phase_rad","exact_phase_rad")
    )
    for j,filename in enumerate(("flow.csv","profiles.csv")):
        with (out/filename).open("x",newline="") as file:
            writer=csv.writer(file)
            writer.writerow(("radial_cells",*headers[j]))
            for info,*values in data:
                writer.writerows((info["cells"],*row) for row in values[j])
    fig,(profile,flow)=plt.subplots(1,2,figsize=(10.7,4.2),layout="constrained")
    if study==43:
        fig.suptitle("Steady pressure-driven flow in a circular pipe",fontweight="bold")
        for info, rows, values in data:
            values=np.array(values)
            last=values[values[:,0]==1.6]
            profile.plot(last[:,1]*1000,last[:,2],label=f"{info['cells']} cells")
            rows=np.array(rows)
            flow.plot(rows[:,0],rows[:,1]*1e6,"o-",label=f"{info['cells']} cells")
        profile.plot(last[:,1]*1000,last[:,3],"k--",label="Exact")
        flow.plot(rows[:,0],rows[:,2]*1e6,"k--",label="Exact")
        profile.set(xlabel="Radius (mm)",ylabel="Axial velocity (m/s)",title="Profile at 1.6 m/s²")
        flow.set(xlabel="Acceleration from pressure gradient (m/s²)",ylabel="Volume flow (mL/s)",title="Flow rate versus forcing")
    else:
        fig.suptitle("Oscillating pressure-driven flow in a circular pipe",fontweight="bold")
        for info,rows,values in data:
            values=np.array(values)
            profile.plot(values[:,0]*1000,values[:,1],label=f"{info['cells']} cells")
            rows=np.array(rows)[-51:]
            flow.plot(rows[:,0],rows[:,1]*1e6,label=f"{info['cells']} cells")
        profile.plot(values[:,0]*1000,values[:,2],"k--",label="Womersley")
        flow.plot(rows[:,0],rows[:,2]*1e6,"k--",label="Womersley")
        profile.set(xlabel="Radius (mm)",ylabel="Velocity amplitude (m/s)",title="Radial amplitude from the final cycle")
        flow.set(xlabel="Time (s)",ylabel="Volume flow (mL/s)",title="Flow through one settled cycle")
    for ax in (profile,flow):
        ax.grid(alpha=.23)
        ax.legend(frameon=False)
    fig.savefig(out/"pipe-flow.png",dpi=180)
    fig.savefig(out/"pipe-flow.pdf")
    print(json.dumps(checks,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study",type=int,choices=(43,44))
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    main(args.study,args.out)
