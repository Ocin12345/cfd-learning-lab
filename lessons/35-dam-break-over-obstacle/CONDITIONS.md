# Model 35: dam break over a low obstacle

Create a visually clear, actual Fluent/OpenFOAM comparison of a released water column crossing a low obstacle. This follows the completed two-grid sloshing study. The owner has opted out of prediction, reflection and design-approval pauses; continue routine technical work automatically, showing actual solver images and preserving all versions.

Use the official Foundation14 damBreakLaminar geometry as a reference: a0.584 x0.584m domain, water column0.1461 x0.292m, obstacle x=0.292..0.31599948m and y=0..0.04799896m. These obstacle coordinates retain the tutorial's actual scaled vertices. Double each block's original planar divisions to9,072 quadrilateral cells. Front/back are empty in OpenFOAM's0.0146m layer; Fluent is planar with unit depth.

Water rho1000kg/m3, mu0.001Pa s; air rho1kg/m3, mu1.48e-5Pa s; g=(0,-9.81)m/s2, surface tension0.07N/m, no wall adhesion. Both solvers are2D laminar incompressible/isothermal. This classroom idealization does not establish experimental validation, turbulent impact-load accuracy or aircraft validity.

Walls are stationary/no-slip. The top is a static0Pa-gauge atmosphere, admitting air only on backflow. Use Fluent pressure-outlet with Static Pressure backflow and user-input operating density0; use Foundation14 prghPressure with p uniform0, pressureInletOutletVelocity and alpha inletOutlet/inletValue0. This intentionally replaces the official tutorial's total-pressure atmosphere to match a static pressure condition across codes. It is a documented tutorial adaptation, not a claimed exact benchmark replay.

Initialize exact rectangular volume fractions and the same physical hydrostatic gauge pressure in both codes; convert that pressure to p_rgh in OpenFOAM. Do not use cell-centre membership as an approximate water volume. Cell areas vary, so integrate alpha and mixture inventory with actual cell area, not a uniform cell count.

Use dt0.0001s. First run a same-water-area flat control and a released-column qualification for0.02s/200steps each; save every10steps. Production may run0.5s/5,000steps with101 states every50steps only after both codes' controls, startup numerics and native open-boundary mass-flux accounting pass. Keep fixed dt and actual native times for matched visualizations.

Fluent retains explicit Geo-Reconstruct VOF/PISO/implicit body force/body-force-weighted pressure/second-upwind momentum, Euler time integration, dynamic scaled residual targets1e-4 and400 inner cap. Quiet raw continuity1e-12 and momentum1e-10; do not relax limits after seeing results. OpenFOAM starts from accepted sloshing pressure/alpha tolerances1e-10/1e-11, outer8/pressure3, MULES10/tolerance0.01, momentumPredictor yes. All final-outer initial Ux/Uy/p_rgh targets1e-4. Reported interface/VOF Courant<=0.25, finite fields, alpha within1e-6 of[0,1]. Quiet saved Umax<=1e-4m/s and alpha change<=0.01.

The top is open: area loss alone is not a conservation error. Record per-step native mixture mass flux and mixture inventory, account for sign and OpenFOAM depth, and check inventory change plus integrated outward flux within0.5% of initial mixture mass. Verify report cadence and all native clocks. Do not extend production with an unverified flux proxy.

All files stay on E:/CFD-Learning-Lab, Fluent installation on D:. Hidden/background, two ranks, one solver at a time. Require available RAM>=4GiB and E: free>=10GiB, rechecked after imports and immediately before launch. Save before closing. Preserve model34 and all originals; keep native outputs outside Git. Local record: actual figures/video, code, small inputs, numerical CSVs, conditions, limitations and honest contributions. No publication or personal reflection is invented.
