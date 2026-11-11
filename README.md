# CFD Learning Lab

The portfolio contains **25 distinct CFD projects**, covering external and internal flow, heat transfer, compressible flow and free surfaces. Each project links to its simulation setup, computed figures and numerical checks.

**[Browse the CFD projects](docs/MODEL_INDEX.md)** · [Study records and benchmarks](docs/STUDY_INDEX.md)

## Selected results

| Study | What the simulations showed | Figures and video |
| --- | --- | --- |
| [Cylinder wake](lessons/20-cylinder-benchmark/convergence-v001/README.md) | The two solvers gave mean shedding periods of 0.33467 s and 0.33488 s over the last six cycles. | [Force histories](lessons/20-cylinder-benchmark/convergence-v001/forces-and-pressure.png) · [Video](lessons/20-cylinder-benchmark/convergence-v001/cylinder-developed.mp4) |
| [Dam break over an obstacle](lessons/35-dam-break-over-obstacle/full-production-v001/README.md) | Both solvers formed a jet as the water crossed the obstacle. Breakup details differed on the shared grid. The declared checks passed over the selected 0.5 s trajectory in each solver. | [Numerical checks](lessons/35-dam-break-over-obstacle/full-production-v001/numerical-checks.png) · [Video](lessons/35-dam-break-over-obstacle/full-production-v001/dam-break-development.mp4) |
| [Cooling an aluminium block](results/cooling-block-plain-v001/README.md) | The block cooled from 320 K to about 318.533 K over 20 s in both solvers. Their final mean temperatures differed by about 0.00050 K. | [Energy checks](results/cooling-block-plain-v001/figures/paired-prefix-energy-checks.png) · [Video](results/cooling-block-plain-v001/cooling-block-plain-comparison.mp4) |
| [Finite-wing convergence investigation](lessons/29-wing-mesh-sensitivity/formulation-investigation-v001/README.md) | All five Fluent trials missed the continuity target. Residual maps show where the imbalance was concentrated downstream of the wing. | [Residual maps](lessons/29-wing-mesh-sensitivity/formulation-investigation-v001/figures/three-actual-residual-maps.png) · [Convergence histories](lessons/29-wing-mesh-sensitivity/formulation-investigation-v001/figures/qualification-numerics.png) |
| [Shock tube](projects/shock-tube/README.md) | The OpenFOAM pressure profile at 7 ms had a mean absolute error of 0.110% of the high-side pressure against the exact solution. | [Profile comparison](projects/shock-tube/data/openfoam-shock-tube.png) |

| Dam break at 0.25 s | Cooling block at 20 s |
| --- | --- |
| [<img src="lessons/35-dam-break-over-obstacle/full-production-v001/state-050.png" width="290" alt="Fluent and OpenFOAM water fraction and mixture speed at 0.25 seconds">](lessons/35-dam-break-over-obstacle/full-production-v001/state-050.png) | [<img src="results/cooling-block-plain-v001/figures/paired-temperature-and-cooling.png" width="340" alt="Fluent and OpenFOAM temperature fields at 20 seconds and cooling histories">](results/cooling-block-plain-v001/figures/paired-temperature-and-cooling.png) |

## Methods and limitations

Each study describes its numerical checks and records any missed targets. Some studies vary the mesh, time step or domain size to measure how individual results respond. Experimental validation remains open.

[Reproduction notes](docs/REPRODUCING.md) · [Corrections](docs/ERRATA.md)
