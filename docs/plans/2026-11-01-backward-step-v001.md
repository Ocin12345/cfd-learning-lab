# Backward-facing step implementation plan

> **For agentic workers:** Use superpowers:subagent-driven-development or superpowers:executing-plans for the bounded case.

**Goal:** Show actual separation and reverse flow in a small introductory Fluent/OpenFOAM model.

**Architecture:** Shared8,800-cell stair-step mesh, two sequential two-process solves, native-field reconstruction and clear common-scale images. Essential numerical checks only.

**Tech Stack:** FluentStudent2026R1/PyFluent0.43, OpenFOAMFoundation14, numpy/h5py/Matplotlib.

**Spec:** E:/CFD-Learning-Lab/local/backward-step-v001/DESIGN.md, with portable conditions copied to the final model record.

## Constraints and review focus

Use E: for project artifacts and temp, existing D: ANSYS install; preserve prior projects. Uniform upstream inlet is explicit; no benchmark matching or aircraft validation. Correct solid mask/native cell mapping, reverse-flow sign, residual/mass checks and100-iteration saved-endpoint stability. Two processes, one solver session at a time, current resource/licence checks. No publication/commit.

- [x] Generate and geometrically check the shared mesh and record conditions.
- [x] Run/save/close Fluent with native images and basic checks.
- [x] Convert/check/run/reconstruct/save OpenFOAM with basic checks.
- [x] Plot and visually inspect actual separation/recirculation figures; show the images.
- [x] Save small methods/data/code/contribution record and verify source hashes.


Delivered: saved calculations, checked actual fields and paired images. Final
residual/mass/saved-endpoint limits pass; an additional Fluent last-100
continuity diagnostic narrowly fails and remains false in the result record.
Completing these implementation tasks does not mean all diagnostics passed or
physical validation was achieved.
