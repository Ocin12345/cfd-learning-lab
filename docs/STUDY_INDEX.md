# CFD study index

There are 35 numbered study records. Entries include geometry models, numerical sensitivity studies and comparisons; they are not 35 independently validated physical predictions. Read each record for its actual numerical status and limitations.

| Study | Record |
| --- | --- |
| 01 | [Lesson 1: when is a cylinder calculation settled?](../lessons/01-cylinder-review/README.md) |
| 02 | [Lesson 2: mesh agreement versus boundary sensitivity](../lessons/02-cylinder-sensitivity/README.md) |
| 03 | [Lesson 3: repair the mesh without changing the airfoil](../lessons/03-airfoil-mesh/README.md) |
| 04 | [Matched cylinder calculation: Fluent and OpenFOAM](../lessons/04-cross-solver-cylinder/README.md) |
| 05 | [Laminar channel flow: a visual starter model](../lessons/05-channel-flow/README.md) |
| 06 | [Backward-facing step: separation and reverse flow](../lessons/06-backward-step/README.md) |
| 07 | [Cylinder vortex shedding: an actual transient model](../lessons/07-vortex-shedding/README.md) |
| 08 | [A finite-span wing: actual 3D Fluent and OpenFOAM fields](../lessons/08-finite-wing/README.md) |
| 09 | [Compressible flow through a nozzle](../lessons/09-nozzle/README.md) |
| 10 | [Lid-driven cavity at Re = 100](../lessons/10-lid-driven-cavity/README.md) |
| 11 | [Two streams through a T-junction](../lessons/11-t-junction/README.md) |
| 12 | [Flow around a planar channel bend](../lessons/12-channel-bend/README.md) |
| 13 | [Heat spreading through a solid slab](../lessons/13-solid-conduction/README.md) |
| 14 | [A square obstacle and its steady wake](../lessons/14-square-obstacle-wake/README.md) |
| 15 | [Heating a channel with its computed flow held fixed](../lessons/15-heated-channel/README.md) |
| 16 | [Annular Couette spin-up](../lessons/16-annular-couette/README.md) |
| 17 | [Offset tandem-square wake](../lessons/17-tandem-offset-wake/README.md) |
| 18 | [Rotating-cylinder startup](../lessons/18-magnus-cylinder/README.md) |
| 19 | [A tracer pulse through a T-junction](../lessons/19-tracer-tjunction/README.md) |
| 20 | [Cylinder wake: comparison with a published periodic benchmark](../lessons/20-cylinder-benchmark/README.md) |
| 21 | [Two moving walls, two counter-rotating cells](../lessons/21-double-lid-cavity/README.md) |
| 22 | [What five degrees changes: NACA 0012 in Fluent and OpenFOAM](../lessons/22-airfoil-incidence/README.md) |
| 23 | [23. Dam break in a closed tank](../lessons/23-dam-break/README.md) |
| 24 | [24. Mach 2 flow over a compression wedge](../lessons/24-mach2-wedge/README.md) |
| 25 | [25. A planar jet impinging on a plate](../lessons/25-impinging-jet/README.md) |
| 26 | [26. A cold jet cooling a hot plate](../lessons/26-heated-impinging-jet/README.md) |
| 27 | [27. Revisiting the three-dimensional wing](../lessons/27-wing-convergence/README.md) |
| 28 | [28. How incidence changes the flow around a finite wing](../lessons/28-wing-incidence/README.md) |
| 29 | [29. Does a finer mesh change the wing results?](../lessons/29-wing-mesh-sensitivity/README.md) |
| 30 | [30. Heat drives a circulation](../lessons/30-natural-convection/README.md) |
| 31 | [31. Does the animation change with a smaller time step?](../lessons/31-convection-time-step/README.md) |
| 32 | [32. Heat above a localized floor heater](../lessons/32-heater-plume/README.md) |
| 33 | [33. Water-tank sloshing](../lessons/33-tank-sloshing/README.md) |
| 34 | [34. Sloshing mesh sensitivity](../lessons/34-sloshing-mesh-sensitivity/README.md) |
| 35 | [35. Dam break over an obstacle](../lessons/35-dam-break-over-obstacle/README.md) |


## Supplementary material

[Cooling an aluminium block in airflow](../results/cooling-block-plain-v001/README.md). Supplement with paired Fluent and OpenFOAM figures, video, methods and numerical evidence.

## Refined-wing formulation investigation — 8 November 2026

The [completed investigation record](../lessons/29-wing-mesh-sensitivity/formulation-investigation-v001/README.md) preserves all five finite trial misses, actual retained-residual maps and the source-bound SIMPLEC/OpenFOAM comparison. Refined Fluent convergence remains **unresolved**; no confirmation or additional solver run is admitted. The separate zero-incidence miss remains open.
