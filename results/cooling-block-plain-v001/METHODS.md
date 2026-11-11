# Methods and numerical checks

[Return to the comparison and figures](README.md)

## Geometry, materials and solver setup

The planar channel measures 160 × 40 mm. A 20 × 1 mm block sits on the floor from x = 60 to 80 mm. The common 320 × 80 grid has 25,520 fluid cells and 80 solid cells. Air starts at 300 K and 0.1 m/s with zero gauge pressure. The inlet supplies the same temperature and speed, the outlet is held at zero gauge pressure with 300 K backflow, and the stationary channel walls are no-slip and adiabatic. The block base is adiabatic, and the exposed faces exchange heat with the fluid. Heat and energy are reported per metre of depth using Fluent's unit depth and OpenFOAM's 1 mm extrusion with empty front and back patches.

| Material | Density (kg/m³) | Specific heat (J/kg·K) | Conductivity (W/m·K) | Dynamic viscosity (Pa·s) |
| --- | ---: | ---: | ---: | ---: |
| Air | 1.225 | 1006 | 0.0262 | 1.7894 × 10⁻⁵ |
| Aluminium | 2700 | 900 | 205 | — |

Both solvers advance laminar flow and conjugate heat transfer together at fixed Δt = 0.01 s with first-order implicit Euler time discretization. The materials have constant properties. The model omits buoyancy, radiation and contact resistance. The solid starts at 320 K and cools by conduction and heat exchange with the air.

Fluent 2026 R1 uses the Coupled pressure–velocity scheme with least-squares gradients, second-order pressure and momentum discretization, first-order upwind energy, and up to 120 inner iterations per step. OpenFOAM Foundation 14 uses the fluid and solid modules of foamMultiRun, with Gauss linear gradients, linear-upwind momentum, upwind fluid enthalpy, orthogonal diffusion, 16 outer iterations and two pressure correctors.

Each solver continued from its accepted 20-step checkpoint at 0.2 s for another 1,980 steps. The combined histories contain 2,001 unique states including the initial condition. After a Windows reboot interrupted OpenFOAM post-processing, the missing output fields were reconstructed from the saved partitions with the same native utility.

## Acceptance checks

| Check | Limit | ANSYS Fluent | OpenFOAM |
| --- | ---: | ---: | ---: |
| Largest relative boundary mass imbalance | 1 × 10⁻⁶ | 5.683 × 10⁻⁸ | 4.616 × 10⁻¹² |
| Largest normalized energy mismatch in any single step | 0.005 | 4.041 × 10⁻⁶ | 4.584 × 10⁻⁶ |
| Largest normalized cumulative energy mismatch at 20 s | 0.005 | 3.494 × 10⁻⁶ | 5.112 × 10⁻⁹ |

The energy maxima in the table are taken across the fluid, solid, combined domain and interface checks. The cumulative checks also passed at every earlier time step. Fluent’s scaled continuity and velocity residuals were at most 10⁻⁵, and its scaled energy residual was at most 10⁻⁸. OpenFOAM’s initial velocity and pressure residuals in the last outer iteration were at most 10⁻⁵, and its final fluid and solid thermal residuals were at most 10⁻¹⁰. Residual values follow each solver’s native definition.

## Reading the energy and pressure plots

The energy balance compares the change in stored thermal plus kinetic energy with boundary energy transfer integrated at the right endpoint of each Euler step. Heat leaving the solid enters the fluid and cancels in the combined balance. The mismatch is normalized using independently determined energy transfer and storage terms, as appropriate to each balance.

Fluent heat comes from native face integrals checked against the solver reports, with a 300 K sensible enthalpy reference. OpenFOAM heat is reconstructed from saved fields and the upwind and orthogonal discrete operators. Fluent’s native enthalpy includes pressure departure and kinetic energy. OpenFOAM stores static enthalpy and kinetic energy separately. The storage calculation removes the incompressible pressure contribution consistently.

The pressure plot uses mean fluid-cell pressures at x = 0.25 and 159.75 mm. The difference includes developing channel flow along with the effect of the block.

## Saved outputs

The [CSV](cooling-heat-allstep.csv) gives mean and maximum solid temperature and interface heat transfer into the air for each solver. The initial state has no time-step heat report. Heat rates use W per metre of depth.

The MP4 and GIF each contain 101 frames taken from saved states at 0.2 s intervals, without interpolation. The same 300–320 K temperature scale is used throughout. All frames were decoded and checked. The [evidence JSON](comparison-evidence.json) lists the numerical thresholds and maxima, SHA-256 hashes for the source audits and exports, and the media checks.

These checks establish numerical agreement and discrete conservation for this fixed 2D case. Experimental validation and mesh, time-step and domain independence remain open.
