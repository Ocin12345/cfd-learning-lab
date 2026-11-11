
"""Integrate the 2D vorticity equation and compare Taylor–Green decay."""
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NU = .01
END = 4.
WAVENUMBER = 2


def run(n, dt):
    h = 2 * np.pi / n
    x = (np.arange(n) + .5) * h
    xx, yy = np.meshgrid(x, x)
    u0 = np.sin(2 * xx) * np.cos(2 * yy)
    v0 = -np.cos(2 * xx) * np.sin(2 * yy)
    omega = 4 * np.sin(2 * xx) * np.sin(2 * yy)
    symbol = 4 * np.sin(np.pi * np.fft.fftfreq(n))**2 / h**2
    eigenvalues = symbol[:, None] + symbol[None, :]
    inverse = np.divide(1, eigenvalues, out=np.zeros_like(eigenvalues), where=eigenvalues > 0)

    def derivative(a, axis):
        return (np.roll(a, -1, axis) - np.roll(a, 1, axis)) / (2*h)

    def velocity(w):
        psi = np.fft.ifft2(np.fft.fft2(w) * inverse).real
        return derivative(psi, 0), -derivative(psi, 1)

    def rhs(w):
        u, v = velocity(w)
        lap = sum(np.roll(w, d, axis) for axis in (0, 1) for d in (-1, 1)) - 4*w
        return NU/h**2 * lap - u*derivative(w, 1) - v*derivative(w, 0)

    steps = round(END/dt)
    write_interval = round(.1/dt)
    history = []
    first_energy = None
    for step in range(steps+1):
        if step % write_interval == 0:
            u, v = velocity(omega)
            t = step * dt
            exact_amplitude = np.exp(-2*NU*WAVENUMBER**2*t)
            energy = np.mean(u*u+v*v)/2
            if first_energy is None:
                first_energy = energy
            velocity_error = np.sqrt(np.mean((u-u0*exact_amplitude)**2+(v-v0*exact_amplitude)**2))
            divergence = np.abs(derivative(u, 1)+derivative(v, 0)).max()
            history.append((t,energy/first_energy,exact_amplitude**2,velocity_error,divergence,np.mean(omega)))
        if step < steps:
            a = rhs(omega)
            b = rhs(omega+dt*a/2)
            c = rhs(omega+dt*b/2)
            d = rhs(omega+dt*c)
            omega += dt*(a+2*b+2*c+d)/6
            if not np.isfinite(omega).all():
                raise ValueError(f"Nonfinite vorticity at step {step+1}")
    return np.array(history), omega


def main(out):
    cases = [(32,.005),(64,.0025)]
    results = [run(*case) for case in cases]
    checks = {}
    for (n,dt),(history,_) in zip(cases,results):
        assert history.shape==(41,6) and np.allclose(history[:,0],np.arange(41)*.1,rtol=0,atol=1e-12)
        energy_error = float(np.abs(history[:,1]-history[:,2]).max())
        checks[str(n)] = {
            "grid_points_per_side":n,"steps":round(END/dt),
            "maximum_normalized_energy_error":energy_error,
            "maximum_velocity_rms_m_s":float(history[:,3].max()),
            "maximum_discrete_divergence_s_inverse":float(history[:,4].max()),
            "maximum_absolute_mean_vorticity_s_inverse":float(abs(history[:,5]).max()),
            "final_energy_fraction":float(history[-1,1]),
            "final_exact_energy_fraction":float(history[-1,2]),
        }
        assert energy_error < .01
        assert history[:,3].max() < .02
        assert history[:,4].max() < 1e-10 and np.abs(history[:,5]).max() < 1e-10
        assert np.max(np.diff(history[:,1])) < 1e-12
    assert checks["64"]["maximum_normalized_energy_error"] < checks["32"]["maximum_normalized_energy_error"]/2
    out.mkdir()
    (out/"checks.json").write_text(json.dumps(checks,indent=2)+"\n")
    with (out/"decay.csv").open("x",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(("grid_points_per_side","time_s","energy_fraction","exact_energy_fraction","velocity_rms_m_s","max_divergence_s_inverse","mean_vorticity_s_inverse"))
        for (n,_),(history,_) in zip(cases,results):
            writer.writerows((n,*row) for row in history)
    fig,(fields,plot)=plt.subplots(1,2,figsize=(11,4.5),layout="constrained")
    fig.suptitle("Two-dimensional Taylor–Green vortex decay",fontweight="bold")
    image=fields.imshow(results[1][1],origin="lower",extent=(0,2*np.pi,0,2*np.pi),vmin=-4,vmax=4,cmap="RdBu_r")
    fig.colorbar(image,ax=fields,label="Vorticity (1/s)")
    fields.set(xlabel="x (m)",ylabel="y (m)",title="Computed fine-grid vorticity at 4 s")
    for (n,_),(history,_) in zip(cases,results):
        plot.plot(history[:,0],history[:,1],label=f"{n} × {n}")
    history=results[-1][0]
    plot.plot(history[:,0],history[:,2],"k--",label="Exact")
    plot.set(xlabel="Time (s)",ylabel="Kinetic energy / initial energy",title="Kinetic energy from the time integration")
    plot.grid(alpha=.25)
    plot.legend(frameon=False)
    fig.savefig(out/"vortex-decay.png",dpi=180)
    fig.savefig(out/"vortex-decay.pdf")
    print(json.dumps(checks,indent=2))


if __name__=="__main__":
    main(Path(sys.argv[1]))
