# Saved backward-step flow fields

Solvers shown: ANSYS Fluent, OpenFOAM. Step h=0.1 m; upstream x=-0.4..0 m and y=0.1..0.2 m; downstream x=0..2 m and y=0..0.2 m. Uniform inlet 1 m/s; rho=1 kg/m3; mu=0.001 Pa.s; Re_h=100. Stationary no-slip walls; zero-gauge-pressure outlet; 8,800 quadrilateral fluid cells per code.

step-computed-flow.png/svg shows flat per-cell axial velocity and integrated streamline arrows. step-recirculation-closeup.png/svg magnifies the downstream lower-channel region. Both codes and both figures share actual-data axial-velocity limits -0.162762827883 to 1.46261412253 m/s. When reverse flow exists, zero is centred in the diverging colour map; the negative and positive halves each use a separate linear interval to resolve weak negative velocities. Thus equally saturated blue and red need not have equal velocity magnitudes. A dashed contour marks u_x=0 only where that contour exists. Grey is the masked solid step. Axes x/h,y/h use equal geometric scale. Each panel uses one streamplot call with three actual negative-u cell-centre seeds across the pocket depth and eight inlet seeds. The shared occupancy mask and broken_streamlines=True stop repeated overlapping orbits. Curves integrate the native mapped vector field and introduce no synthetic vortex. No boundary velocity samples are appended.

The first-cell crossing proxy, recorded below and in checks.json, is omitted from the image. It is the linear negative-to-nonnegative crossing in the first row of cell-centred axial velocities at y=0.0025 m, or y/h=0.025. This is a grid-dependent near-bottom proxy, not exact wall reattachment, and is not benchmark or physical validation.

- ANSYS Fluent: minimum axial velocity -0.16273 m/s; 468 negative-u cells; first-cell crossing proxy x/h=4.838.
- OpenFOAM: minimum axial velocity -0.162763 m/s; 466 negative-u cells; first-cell crossing proxy x/h=4.838.

Fluent centres are reconstructed from four unique saved face-connected vertices per verified rectangle; OpenFOAM centres come from its saved C field. Both must map bijectively to exactly the 8,800 fluid slots within 1e-9 m; the 800 solid slots are masked. Imported cell ordering is never assumed. OpenFOAM kinematic pressure is multiplied by rho before CSV export. Input hashes, executed plotter, cell data and mapping diagnostics accompany the figures.

Residual, mass and iteration-stability checks are separate native-run evidence. No Poiseuille analytical target applies to this geometry; no mesh-independence or physical-validation claim is made. Fluent/OpenFOAM calculated the fields.
