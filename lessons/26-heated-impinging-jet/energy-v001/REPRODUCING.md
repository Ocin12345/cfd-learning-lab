# Inspect or reproduce this execution archive

This is a methods and execution archive, not a complete public native restart bundle. Start with the [project environment guide](../../../docs/REPRODUCING.md) and [parent thermal model](../README.md). Exact readback requires the separate E:/CFD-Learning-Lab native archive and existing D: ANSYS Fluent 261 installation. Environments, meshes, checkpoints and licence files are excluded here.

## Local source layout

The executed source folder is `E:/CFD-Learning-Lab/local/showcase-improvements-2026-11-08-v001/heated-jet-energy-v001`. Add it, `local/heated-impinging-jet-v001`, `local/impinging-jet-v001`, `repo/cfd-learning-lab/scripts` and `local/showcase-improvements-2026-11-08-v001/sloshing-time-step-v001` to PYTHONPATH. Dot-source `local/Use-Project.ps1` and use the existing `environments/fluent/Scripts/python.exe -B`.

Hash guards bind the original local source paths. Copies here document execution; they are not relocated launch entry points. Saved audits and selectors identify the exact accepted cases. Import native readers into a **new versioned read-only analysis** to avoid overwriting accepted audit filenames. The native auditors never launch a solver, but their original entry points write fixed audit paths.

For synthetic checks, select `test_jet_energy_v001`, `test_jet_energy_restart_v001` and `test_jet_energy_parent_gate_v001` (23 tests). `test_native_clock_v002` adds two tests requiring preserved native Fluent data. The intentionally failing historical clock-reader regression is under `code/history` and is not part of the current green suite.

## A new calculation

Create new versioned source/output paths and regenerate reviewed source/parent manifests. Preserve the exact native thermal time-zero case/data, original frozen-flow checkpoints and both OpenFOAM partitions/native meshes. First qualify 20 steps and inspect the paired fields and audits. Only then admit the exact 0.1 s checkpoint for 580 additional steps.

Do not rerun frozen controllers over completed directories, reset counters, reinitialize or redecompose. OpenFOAM source time 1250 denotes steady iterations, not physical seconds. Keep the original physical and numerical controls for this output-only comparison.

Run hidden, with two ranks and one solver at a time, sequential Fluent then OpenFOAM. Verify executables, licence, current source/parent hashes and Windows Fluent/cortex plus WSL foam/MPI process conflicts. Require RAM ≥ 4 GiB and E: free ≥ 10 GiB after imports, immediately before launch and after waits. Data, environments, temp and caches stay on E:; use the existing D: Fluent installation. Only the recorded isolated /dev/shm MPI FIFO metadata exception is used where E: lacks FIFO support. Save before closing; no desktop takeover or automatic publication.

The 61 display states cover 0:0.05:3 s at 10 fps without interpolation; full histories use all 600 step samples. Every MP4/GIF frame was decoded and compared with its actual source PNG. All 61 source PNGs and large native outputs remain in the separate E: export archive. Representative PNGs, complete video/GIF and verification receipts are included here.
