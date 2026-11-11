# Study 40 protocol

A periodic 2π m square starts with the two-dimensional Taylor–Green velocity field, initial peak speed 1 m/s, wavenumber 2 m⁻¹, and kinematic viscosity 0.01 m²/s. Integrate the nonlinear vorticity equation through 4 s on 32 × 32 and 64 × 64 periodic grids with time steps 0.005 and 0.0025 s. Use centered differences, an FFT inversion of the discrete Laplacian and RK4. Sample every 0.1 s.

Require finite vorticity, monotone kinetic-energy decay at saved times, velocity RMS below 0.02 m/s, maximum energy-fraction error below 0.01, discrete divergence and mean vorticity below 10⁻¹⁰ s⁻¹, and at least a factor of two reduction in maximum energy error on the finer calculation. The exact velocity amplitude decays as exp(−2νk²t). The exact energy fraction is its square.