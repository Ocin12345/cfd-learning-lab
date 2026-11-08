# Reproduction and storage

This is an execution archive with small inputs and source-bound results, not a complete native restart bundle or a fresh-machine installation package. The exact native cases/checkpoints, processor fields, all101 PNG frames, full logs, meshes and installed environments remain outside Git on E:/CFD-Learning-Lab.

The selected cases and audit paths are recorded in data/accepted-pair-selection.json. Source receipts intentionally refer to preserved E: originals. Current scripts reject reused output directories and require those hash-bound prerequisites. Do not run them blindly against completed paths.

On the original workstation, dot-source E:/CFD-Learning-Lab/local/Use-Project.ps1, use the existing E: Python environment and include both local/dam-break-over-obstacle-v001 and repo/cfd-learning-lab/scripts on PYTHONPATH. Fluent remains on D: and OpenFOAM Foundation14 is in CFD-Foundation14 WSL. Recheck resources/licence and use a new source/output version for any future run. No new installation is required for inspecting the PNG/PDF/GIF/MP4/CSV results.

Executed full Fluent driver008/auditor006 and the source-bound runtime guard are in code/. The accepted OpenFOAM continuation and complete-trajectory audit are both in run_scheduled_openfoam_pipeline_v001.py. Its5000step prefix comes from the preserved original parent, whose later numerical miss remains in history/. The native partitioned0.25s checkpoint is required for that exact continuation; it is not included here.

Small OpenFOAM dictionaries and initial alpha/U/p_rgh are in openfoam-inputs/. Their startTime is0.25 and they require the external native checkpoint, even though the copied0/ fields document the original initialization. The relaxation schedule and continuation declaration are retained. Mesh coordinates/analytical initialization are defined by the common geometry code.

To regenerate figures from the external selected native cases, use render_dam_break_production_v002.py with the selected case/audit arguments and a fresh output folder. It explicitly checks the accepted pair selection, native hashes and nonuniform mapping. The source audit binds all101 frame hashes. The MP4/GIF exporters bind those images and decode-check every frame; export records contain the source hashes and actual physical times.

All source and result copies were SHA-256 verified. Installed dependencies, licence/connection files, credentials and large native outputs are excluded. See environment.json for actual package versions and the root project reproduction guide for baseline helpers.
