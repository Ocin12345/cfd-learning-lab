# CFD study index

Browse 48 numbered studies. The original study numbers identify related runs and follow-ups. Open each study for the setup, figures, and available numerical evidence. Diagnostic runs and missed targets are documented in the study records.

[Portfolio home](../README.md) · [Aerodynamics and wakes](#aerodynamics-and-wakes) · [Internal flows and transport](#internal-flows-and-transport) · [Heat transfer](#heat-transfer) · [Compressible flow](#compressible-flow) · [Free-surface flow](#free-surface-flow)

## Aerodynamics and wakes

| No. | Study | What it covers |
| --- | --- | --- |
| 01 | [Cylinder force stability](../lessons/01-cylinder-review/README.md) | Checks whether residuals and forces have settled in a saved cylinder run. |
| 02 | [Cylinder domain sensitivity](../lessons/02-cylinder-sensitivity/README.md) | Compares the effect of moving the outer boundary on drag. |
| 03 | [Airfoil mesh repair](../lessons/03-airfoil-mesh/README.md) | Repairs poor cells while preserving the airfoil geometry. |
| 04 | [Matched cylinder calculation](../lessons/04-cross-solver-cylinder/README.md) | Compares Fluent and OpenFOAM forces on a shared case. |
| 07 | [Cylinder vortex shedding](../lessons/07-vortex-shedding/README.md) | Tracks alternating vortices and force histories in a transient run. |
| 08 | [Three-dimensional finite wing](../lessons/08-finite-wing/README.md) | Compares saved steady flow fields and forces around a finite-span wing. |
| 14 | [Confined square wake](../lessons/14-square-obstacle-wake/README.md) | Compares recirculation and drag behind a square obstacle. |
| 17 | [Offset tandem-square wakes](../lessons/17-tandem-offset-wake/README.md) | Examines interacting wakes and forces on two obstacles. |
| 18 | [Rotating-cylinder startup](../lessons/18-magnus-cylinder/README.md) | Follows the development of the asymmetric flow and body forces. |
| 20 | [Cylinder wake benchmark comparison](../lessons/20-cylinder-benchmark/README.md) | Compares shedding frequency and force cycles. [Developed-cycle comparison](../lessons/20-cylinder-benchmark/convergence-v001/README.md). |
| 22 | [NACA 0012 at 0° and 5°](../lessons/22-airfoil-incidence/README.md) | Compares pressure and lift as incidence changes. [Trailing-edge sensitivity](../lessons/22-airfoil-incidence/trailing-edge-v001/README.md). |
| 27 | [Finite-wing convergence follow-up](../lessons/27-wing-convergence/README.md) | Continues the baseline wing calculation and compares lift and drag. |
| 28 | [Finite-wing incidence](../lessons/28-wing-incidence/README.md) | Compares the computed wing flow at 0° and 8° and records the unresolved zero-incidence Fluent residual target. |
| 29 | [Finite-wing mesh sensitivity](../lessons/29-wing-mesh-sensitivity/README.md) | Compares two tetrahedral meshes. [Refined-wing convergence investigation](../lessons/29-wing-mesh-sensitivity/formulation-investigation-v001/README.md). |

## Internal flows and transport

| No. | Study | What it covers |
| --- | --- | --- |
| 05 | [Laminar channel flow](../lessons/05-channel-flow/README.md) | Compares velocity profiles with the analytical channel solution. |
| 06 | [Backward-facing step](../lessons/06-backward-step/README.md) | Examines separation and reverse flow downstream of the step. |
| 10 | [Lid-driven cavity](../lessons/10-lid-driven-cavity/README.md) | Compares the recirculating flow driven by a moving wall. |
| 11 | [Two inlets at a T-junction](../lessons/11-t-junction/README.md) | Examines the velocity fields where two streams meet. |
| 12 | [Planar channel bend](../lessons/12-channel-bend/README.md) | Follows the flow through a bend in a two-dimensional channel. |
| 16 | [Annular Couette spin-up](../lessons/16-annular-couette/README.md) | Follows circular shear as motion spreads outward from a rotating wall. |
| 19 | [Tracer pulse through a T-junction](../lessons/19-tracer-tjunction/README.md) | Tracks scalar advection and diffusion through a saved flow field. |
| 21 | [Double-lid cavity](../lessons/21-double-lid-cavity/README.md) | Compares counter-rotating cells driven by two moving walls. |
| 25 | [Planar impinging jet](../lessons/25-impinging-jet/README.md) | Compares stagnation pressure and the wall jets formed after impingement. |
| 37 | [Startup of plane Couette flow](../lessons/37-plane-couette-startup/README.md) | Compares two OpenFOAM resolutions with the analytical development of flow between two plates. |
| 38 | [Startup of pressure-driven channel flow](../lessons/38-plane-poiseuille-startup/README.md) | Compares OpenFOAM profiles and flow-rate development with the transient plane Poiseuille solution. |
| 39 | [Oscillating wall in a finite gap](../lessons/39-oscillating-wall/README.md) | Compares OpenFOAM velocity attenuation and phase lag against the analytical startup solution. |
| 40 | [Taylor-Green vortex decay](../lessons/40-taylor-green-vortex/README.md) | Measures energy decay and velocity errors against the exact two-dimensional viscous solution. |
| 41 | [Numerical diffusion of a tracer pulse](../lessons/41-periodic-tracer-advection/README.md) | Compares upwind and limited MUSCL after one periodic crossing with zero physical diffusivity. |
| 42 | [Physical spreading of an advected tracer](../lessons/42-periodic-tracer-diffusion/README.md) | Compares peak decay, variance and mass against a periodic advection–diffusion solution. |
| 43 | [Steady pressure-driven flow in a circular pipe](../lessons/43-steady-pipe-friction/README.md) | Compares velocity, volume flow and wall shear with Hagen–Poiseuille theory. |
| 44 | [Pulsating pressure-driven pipe flow](../lessons/44-pulsating-pipe-flow/README.md) | Checks flow amplitude and phase lag against the Womersley solution. |
| 48 | [Pressure recovery in two planar diffusers](../lessons/48-planar-diffuser-recovery/README.md) | Compares two taper lengths at the same area ratio and documents pressure recovery and reverse flow. |

## Heat transfer

| No. | Study | What it covers |
| --- | --- | --- |
| 13 | [Conduction through a solid slab](../lessons/13-solid-conduction/README.md) | Compares transient temperatures with a conduction reference solution. |
| 15 | [Heating a flowing channel](../lessons/15-heated-channel/README.md) | Tracks temperature evolution using fixed velocity fields. |
| 26 | [Heat transfer from a plate to an impinging jet](../lessons/26-heated-impinging-jet/README.md) | Tracks warming of the air over a plate held at fixed temperature. [Discrete energy accounting](../lessons/26-heated-impinging-jet/energy-v001/README.md). |
| 30 | [Natural convection in a cavity](../lessons/30-natural-convection/README.md) | Follows buoyancy-driven circulation and heat transfer between the side walls. |
| 31 | [Natural-convection time-step sensitivity](../lessons/31-convection-time-step/README.md) | Compares the effects of a smaller time step. [Spatial sensitivity follow-up](../lessons/31-convection-time-step/spatial-v001/README.md). |
| 32 | [Plume above a floor heater](../lessons/32-heater-plume/README.md) | Compares the flow and heat transfer above a localized heater. |
| 36 | [Cooling an aluminium block in airflow](../results/cooling-block-plain-v001/README.md) | Compares 20 s of conjugate heat transfer in Fluent and OpenFOAM, with temperature fields, video and discrete energy checks. |
| 45 | [Thermal entrance region between heated plates](../lessons/45-thermal-entrance-region/README.md) | Compares bulk temperature and wall heat transfer with a fine Graetz eigenmode reference. |
| 46 | [Heat conduction across a bonded interface](../lessons/46-bonded-layer-conduction/README.md) | Tests flux conservation and temperature continuity for two materials with different conductivities. |
| 47 | [Transient heat penetration after a surface temperature step](../lessons/47-transient-heat-penetration/README.md) | Compares temperature depth profiles and wall flux with a half-space diffusion solution. |

## Compressible flow

| No. | Study | What it covers |
| --- | --- | --- |
| 09 | [Compressible nozzle flow](../lessons/09-nozzle/README.md) | Compares Mach number and pressure fields and records both solvers’ missed numerical targets. |
| 24 | [Mach 2 flow over a compression wedge](../lessons/24-mach2-wedge/README.md) | Compares the oblique shock with a weak-shock solution. [Spatial refinement](../lessons/24-mach2-wedge/spatial-v001/README.md). |

## Free-surface flow

| No. | Study | What it covers |
| --- | --- | --- |
| 23 | [Dam break in a closed tank](../lessons/23-dam-break/README.md) | Follows the released water front and its impact on the tank wall. |
| 33 | [Water-tank sloshing](../lessons/33-tank-sloshing/README.md) | Compares the motion and depth of water in the tank. |
| 34 | [Sloshing mesh sensitivity](../lessons/34-sloshing-mesh-sensitivity/README.md) | Compares two grid sizes. [Time-step sensitivity on the fine grid](../lessons/34-sloshing-mesh-sensitivity/time-step-v001/README.md). |
| 35 | [Dam break over an obstacle](../lessons/35-dam-break-over-obstacle/README.md) | Compares jet formation and breakup. [Full calculation, video and checks](../lessons/35-dam-break-over-obstacle/full-production-v001/README.md). |
