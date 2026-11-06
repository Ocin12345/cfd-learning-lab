# Recorded backward-step inputs

These eight dictionaries are exact copies from the completed OpenFOAM Foundation 14 run at `native/runs/2026-11-01-backward-step-v001/openfoam`, relative to the CFD-Learning-Lab workspace. This teaching-case snapshot records inputs; native mesh, solution fields, processor directories and machine logs remain outside the repository. Preserve this directory and prepare a fresh external case for further runs.

Recorded settings: steady 2D laminar `incompressibleFluid`; step height 0.1 m and expansion 2:1; uniform inlet `(1 0 0)` m/s, stationary no-slip walls and zero outlet gauge pressure. Density is 1 kg/m3 and kinematic viscosity is 0.001 m2/s, giving step-height Reynolds number 100. OpenFOAM pressure is kinematic; multiply by density for Pa. The shared mesh has 8,800 quadrilateral fluid cells.

SIMPLE uses `consistent yes`, velocity-equation relaxation 0.9, least-squares gradients, bounded linear-upwind convection, corrected diffusion and 1e-10 absolute linear-solver tolerance. The recorded run used two MPI ranks, a fixed cap of 2,000 iterations, writes every 100 iterations and no residual early stop; checkpoints 1,900 and 2,000 were retained.

Regenerate the shared Fluent ASCII mesh with `scripts/step_mesh.py` and prepare inputs with `scripts/prepare_openfoam_step.py`. In the fresh external case, conversion uses `fluentMeshToFoam -2D 1 input-step.msh`, followed by `checkMesh -allGeometry -allTopology`. The `1` is the extrusion depth in metres. This recorded case is not a publication or benchmark validation; numerical checks are in [openfoam-checks.json](../../results/backward-step-v001/openfoam-checks.json).

The source files were unchanged during copying, and every copied dictionary matched its source SHA-256:

| File | SHA-256 |
| --- | --- |
| 0/U | b262f73625b9c6f70bdffc24c6540f2f071dcc767b7202853e62527d022e3c37 |
| 0/p | 0095f0e14797bafbd2ab440f79f6621b412ca0076c314c565ade2563915a3123 |
| constant/physicalProperties | 0c2cccd25f6ff7637351833f35ea078390e21633c99246122dc74f495fa94b0b |
| constant/momentumTransport | 7498c38a3be52e135c3f5a82c90cc81cbe88a784f678ba1fb5fa24483c164844 |
| system/controlDict | fd4b08b0eadb4b33f62e07e1bd187593a3b19f9424bccad2c71051383d665d79 |
| system/fvSchemes | 19c244039db97868db2f85b18c61d87359ae987a3bc28d448ca61cb51dcadab0 |
| system/fvSolution | d8a867e1cafbc8ca16428ee87ec6621958ccafbb879b3c81246eb6e0ba6e79bf |
| system/decomposeParDict | 5c4d636cb3070a405d3437172ffbf73a3104582ae6347170c91cc11666011797 |
