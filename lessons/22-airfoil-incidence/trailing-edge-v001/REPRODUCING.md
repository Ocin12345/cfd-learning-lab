# Reproduce or audit this record

This repository holds methods, source snapshots, small numerical records and figures. Exact native reproduction also requires the preserved E:/CFD-Learning-Lab tree, its baseline mesh/checkpoints, existing E: Python environment, D: ANSYS Fluent 2026 R1 and OpenFOAM Foundation 14 in WSL CFD-Foundation14. It is not a standalone fresh-machine bundle.

1. Restore the paths and verify selected hashes in data/selected-native-cases.json and each audit. Reuse accepted completed cases; choose fresh versioned outputs for intentional reruns.
2. The baseline is native/runs/2026-11-06-airfoil-incidence-paired-v001. The selected refiner/quality-aware mesher needs that external ASCII mesh. Current geometry proofs govern coordinates, local connectivity and outside-cell preservation; inherited incidence-morph flags are historical.
3. Fluent uses the frozen numerical template, two ranks, native quality check, finite initial solve and force observation window. OpenFOAM uses the parent study preparation script, exact mesh conversion/native checks and two ranks. Saved-state continuations must not be initialized again.
4. Read-only audits require local sources and repository scripts on PYTHONPATH. The finer audit merges solver.log, solver-continuation.log and solver-continuation-v003.log through steady iteration 32,000; the first refinement merges two logs through 9,000. These indices are not elapsed physical simulation time.
5. The selected plotter reads actual saved native fields and hashes. Export to new directories and inspect the images and numerical/source checks.

Code snapshots retain their original E: execution layout and supporting repository-script imports. Frozen manifests and native parent files remain in that external tree. Historical source-archive files are text records, not recommended commands. No installed dependencies, binaries, licence files, credentials or private application material are distributed.

Use hidden/background operation, two ranks and one solver session at a time. Recheck Windows Fluent/cortex and WSL foam/MPI processes, RAM >= 4 GiB and E: free >= 10 GiB after imports and immediately before launch, again after waits. Temp/cache paths use E: where supported; project OpenMPI FIFO metadata uses /dev/shm where E: is unsupported. Save native outputs before closing.
