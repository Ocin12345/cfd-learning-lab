# Circular pipe flow protocols

Solve the axisymmetric fully developed Newtonian axial momentum equation in annular finite volumes on 40 and 80 radial cells. The pipe radius is 5 mm and kinematic viscosity is 0.0001 m²/s. Use symmetric fluxes across the internal faces, zero centreline flux and a no-slip outer wall. Prescribed axial acceleration represents the pressure gradient divided by density.

Study 43 solves the steady momentum balance at 0.2, 0.4, 0.8 and 1.6 m/s². Compare the weighted velocity, flow rate and wall shear with Hagen–Poiseuille. Require all three normalized errors below 0.2%, algebraic residual below 10⁻¹⁰ m/s² and improved flow rate on the finer grid.

Study 44 drives the same pipe at 4 Hz with 0.8 m/s² amplitude from rest. Integrate with Crank–Nicolson using 1 ms and 0.5 ms steps through 1.25 s. Fit the last period from 50 samples and compare the complex velocity profile and flow with Womersley theory. Require flow amplitude error below 0.3%, phase error below 0.006 rad, profile error below 0.3% of forcing divided by angular frequency, and improved flow amplitude on the fine grid.
