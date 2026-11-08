# Discrete thermal energy accounting

The geometry, thermal properties, boundary conditions and first-order schemes remain those of study 26. Each code retains its own previously computed fixed airflow. The plate is prescribed at 320 K; no solid cools and temperature does not feed back on the flow.

## Storage and units

Let A_c = (0.001 m)² be a cell's planar area, rho = 1.225 kg/m³, cp = 1006 J/(kg·K) and T_ref = 300 K. All energy is reported in J per metre of out-of-plane depth.

Fluent stores specific enthalpy h with h_300 = 1861.1000000000229 J/kg. The native data must satisfy h − h_300 = cp(T − 300) within 1e-8 J/kg. Its stored energy is

`E_n = rho * A_c * sum_cells(h_n - h_300)`.

OpenFOAM uses

`E_n = rho * cp * A_c * sum_cells(T_n - 300)`.

Fluent is planar unit depth. OpenFOAM uses a 0.01 m empty extrusion; divide its integrated quantities by that depth.

## Independent boundary input

Q is positive **into the air**. Fluent's native SV_HEAT_FLUX entries are integrated face heat rates: sum them directly. Its SV_FLUX is outward mass flux. Cross-check every patch against actual native heat/mass reports using rtol = 1e-7 and atol = 1e-11, then shift the enthalpy datum:

`Q_300,patch = Q_native,patch + h_300 * mass_outward,patch`.

Plate mass flux is zero. Do not multiply already integrated native face heat by area again.

For OpenFOAM, reconstruct upwind convective input as

`-rho * cp * sum_faces(phi * (T_upwind - 300)) / depth`.

T_upwind is 300 K for inflow and the owner-cell temperature for outflow. Add orthogonal boundary diffusion: fixed 320 K at the plate, fixed 300 K at the inlet, 300 K on backflow and zero normal gradient on outflow, with an adiabatic ceiling. The source-bound uniform mesh has wall-to-owner distance DX/2 and face area DX × depth. Check native topology, phi dimensions and signs, root fields and both partition clocks. This is a scheme reconstruction, not a native OpenFOAM heat report or a continuum wall-gradient estimate.

## Every-step and cumulative checks

At every real step, `storage_rate = (E_n - E_previous) / dt`. Require

`abs(sum(Q_n) - storage_rate) / max(sum(abs(Q_n)), abs(storage_rate)) <= 0.005`.

For every prefix, define `C_n = sum_i(dt * sum_patch(Q_patch,i))` and `A_n = sum_i(dt * sum_patch(abs(Q_patch,i)))`. Require

`abs(C_n - (E_n - E_0)) / max(A_n, abs(E_n - E_0)) <= 0.005`.

An exact zero numerator and denominator gives zero error. There is no artificial denominator floor. Boundary flux is collected independently, never inferred from the storage difference. All 600 right-endpoint samples follow backward Euler; the result is discrete conservation evidence, not exact continuous integration or physical validation.

The 20-step qualification and 580-step continuation meet at the exact 0.1 s checkpoint. Duplicate restart fields are checked but state 20 is counted once. Continued step 21 must have time 0.105 s; step 600 ends at 3 s. Native current Scheme values and both OpenFOAM partition time indices are checked. Missing, duplicated, reset or nonfinite evidence fails. Residual, mass, temperature and frozen-field limits remain unchanged; see the copied spec and audits for exact definitions.

Primary sources and installed-source hashes are in [the reference record](methods/PRIMARY-SOURCES-v001.md). The online Foundation 14 branch is a contextual reference, not an asserted byte-identical copy of the installed source.
