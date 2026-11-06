# Basic laminar channel-flow implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Run one small familiar internal-flow model in Fluent and OpenFOAM and show actual, readable images.

**Architecture:** Shared uniform 2D quadrilateral mesh; sequential hidden two-rank solves; reconstruct fields by saved native cell connectivity and plot common scales. Essential numerical checks and short reproducible records only.

**Tech Stack:** PyFluent0.43, FluentStudent2026R1, OpenFOAMFoundation14, Python/numpy/h5py/Matplotlib, PowerShell, WSL.

**Spec:** Local channel-flow-v001/DESIGN.md; public run README will carry the portable conditions.

## Global constraints

All project outputs/env/temp on E:/CFD-Learning-Lab; ANSYS installation on D:. Preserve originals and fresh output directories. No continuation of paused airfoil fine study. User wants images first and concise explanation. No GitHub publication/commit.

## Review focus

Correct planar rather than circular profile; mesh cell identities reconstructed from native connectivity; OpenFOAM pressure is kinematic with rho=1; actual residual/mass/profile checks; provisional results never labelled converged.

- [x] Prepare shared7,680-cell mesh, check cell areas/boundaries and save conditions/code.
- [x] Check resources/idle solvers, launch hidden two-process Fluent, check quality/convergence/mass, save case/data/native images and close.
- [x] Convert same mesh in OpenFOAM, run full checkMesh, execute two ranks after fresh resource check, reconstruct/save.
- [x] Plot actual native fields/profiles with common scales, visually inspect and show the images.
- [x] Save concise methods/results/source hashes/contribution record in Git-ready paths; keep native data separate.

Executed outcome: both native runs saved and closed. Full mesh, residual, mass, saved-endpoint and analytical-profile/gradient checks passed. Actual PNG/SVG images were visually inspected. The optional empty final Fluent XY export was recovered from saved native HDF history without another solve; executed source and original native files are preserved. Per-code available-RAM admission was4GiB for Fluent and2.5GiB for the small OpenFOAM case, chosen/documented before its first solve; paused fine-airfoil6GiB gate is unchanged. Resource holds and receipts remain private on E:. No publication/commit occurred.
