# Cooling an aluminium block in airflow

**Supplementary material** · ANSYS Fluent 2026 R1 and OpenFOAM Foundation 14

A 320 K aluminium block cools in a 300 K air stream. Over 20 seconds, both solvers predict a drop of about **1.47 K** in the block’s mean temperature. The comparison uses the same two-dimensional grid, material properties, and time step in both solvers.

[Watch the 10-second comparison (MP4)](cooling-block-plain-comparison.mp4) · [View the animation (GIF)](cooling-block-plain-comparison.gif) · [Read the methods and checks](METHODS.md) · [Open the data (CSV)](cooling-heat-allstep.csv)

![Fluent and OpenFOAM temperature fields at 20 seconds, with traces of mean solid temperature and heat transfer into the air](figures/paired-temperature-and-cooling.png)

*Figure 1. Matched temperature fields at 20 s and cooling histories from both solvers. The temperature panels share a 300–320 K color scale.*

## Results at 20 seconds

| Quantity | ANSYS Fluent | OpenFOAM |
| --- | ---: | ---: |
| Mean block temperature (K) | 318.53256 | 318.53306 |
| Peak block temperature (K) | 318.54685 | 318.54732 |
| Initial sensible energy removed, relative to 300 K | 7.337% | 7.335% |
| Heat transfer into the air per metre of depth (W/m) | 3.3870 | 3.3859 |

The final mean block temperatures differ by 0.00050 K. Across the 101 matched output times, the largest cellwise temperature difference was 0.05528 K. Both solvers passed the declared residual, mass and energy checks over 2,000 time steps. The agreement describes this simulation setup. Experimental validation and sensitivity to mesh, time step and domain size remain open.

## Figures and evidence

| Figure | What it shows | Files |
| --- | --- | --- |
| 1. Temperature and cooling | Matched temperature fields, mean block temperature and interface heat transfer | [PNG](figures/paired-temperature-and-cooling.png) · [PDF](figures/paired-temperature-and-cooling.pdf) |
| 2. Velocity and pressure | Flow fields at 20 s on shared scales | [PNG](figures/paired-velocity-and-pressure.png) · [PDF](figures/paired-velocity-and-pressure.pdf) |
| 3. Cooling and pressure difference | Solid energy loss and the pressure difference across the channel | [PNG](figures/paired-cooling-and-pressure-cost.png) · [PDF](figures/paired-cooling-and-pressure-cost.pdf) |
| 4. Residuals and mass balance | Each solver’s residual history and relative mass imbalance | [PNG](figures/paired-numerical-checks.png) · [PDF](figures/paired-numerical-checks.pdf) |
| 5. Regional energy balances | Fluid and solid energy storage, external energy input and interface heat transfer | [PNG](figures/paired-regional-energy.png) · [PDF](figures/paired-regional-energy.pdf) |
| 6. Cumulative energy checks | Normalized mismatch through each time step for the fluid, solid, interface and combined domain | [PNG](figures/paired-prefix-energy-checks.png) · [PDF](figures/paired-prefix-energy-checks.pdf) |

The video uses 101 saved states at 0.2 s intervals. The [cooling history](cooling-heat-allstep.csv) contains the initial state and all 2,000 time steps. The [numerical evidence](comparison-evidence.json) records the check thresholds, observed maxima, source hashes and media verification.

Fluent and OpenFOAM computed the fields.
