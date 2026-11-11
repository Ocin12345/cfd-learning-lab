# Periodic scalar transport protocols

Advect a periodic Gaussian of width 0.06 m with a prescribed 1 m/s velocity in a 1 m domain. Use conservative periodic face fluxes and SSPRK3. Compare numerical cell averages with eight-point quadrature of the translated or broadened periodic Gaussian. Verify finite bounded fields and relative scalar mass error below 10⁻¹¹ at every step.

Study 41 compares first-order upwind and MUSCL with a monotonized central limiter on 160 cells through one crossing, with zero physical diffusivity. Require lower final normalized L1 error and higher retained peak with MUSCL, and MUSCL final normalized L1 error below 4%.

Study 42 uses MUSCL and central physical diffusion at 0.002 m²/s through 0.5 s, using 160 and 320 cells. Require every saved L1 error below 3%, a lower final L1 error on the fine grid, final peak error below 0.01 and fine variance error below 0.00015 m².
