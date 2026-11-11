# CFD Learning Lab

Computational fluid dynamics portfolio and supplemental material. The collection contains **35 numbered studies** and a separate cooling-block comparison. It uses ANSYS Fluent and OpenFOAM to explore how numerical choices affect simulated flows and heat transfer. The study records document the setup, available results and figures, numerical checks, and limitations.

**[Browse the full study index](docs/STUDY_INDEX.md)**

## Selected studies

| Study | What to look for |
| --- | --- |
| [Cylinder wake](lessons/20-cylinder-benchmark/convergence-v001/README.md) | Compare vortex shedding and force cycles in Fluent and OpenFOAM. The record reports how the cycles settle and discusses the limits of comparison with a published benchmark. |
| [Dam break over an obstacle](lessons/35-dam-break-over-obstacle/full-production-v001/README.md) | Follow a water column as it collapses and forms a jet over an obstacle. The two solvers show different breakup details on the shared grid. |
| [Cooling an aluminium block](results/cooling-block-plain-v001/README.md) | Compare temperature fields and the cooling history of a solid block over 20 seconds, with mass and energy checks. |
| [Finite-wing convergence investigation](lessons/29-wing-mesh-sensitivity/formulation-investigation-v001/README.md) | See how five solver trials and downstream residual maps were used to investigate an unresolved Fluent convergence target. |

![Fluent and OpenFOAM water fraction and mixture speed 0.25 seconds after the dam break](lessons/35-dam-break-over-obstacle/full-production-v001/state-050.png)

*Dam break over an obstacle at 0.25 s. Both solvers are shown at the same time and on shared color scales.* [Watch the comparison (MP4)](lessons/35-dam-break-over-obstacle/full-production-v001/dam-break-development.mp4)

## Evidence and methods

The portfolio includes solver comparisons, numerical sensitivity studies and retained diagnostic runs. Each record gives the applicable checks and identifies any missed targets. Some studies vary the mesh, time step, or domain size to measure how specific outputs respond. Experimental validation remains open.

[Study updates](docs/SHOWCASE_IMPROVEMENTS_2026-11-08.md) · [Reproduction notes](docs/REPRODUCING.md) · [Corrections](docs/ERRATA.md) · [Project history](docs/PROJECT_HISTORY.md)

## Contributions

The project owner chose the topics and directed the progression of the studies. Codex prepared meshes and scripts, ran the solver workflows and numerical checks, and drafted the figures and records. ANSYS Fluent and OpenFOAM computed the simulation fields.
