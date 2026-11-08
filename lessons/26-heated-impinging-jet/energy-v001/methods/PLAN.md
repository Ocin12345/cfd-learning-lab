# Heated-jet discrete energy implementation plan

> For agentic workers: use superpowers:executing-plans inline. One fresh whole-improvement review at the end, not per task.

Goal: collect and audit a complete600-step discrete thermal energy ledger without changing the accepted frozen-flow model.
Architecture: immutable baseline snapshots; tested energy arithmetic and native readers; finite sequential qualifications; only after actual acceptance admit exact-checkpoint finite continuations; read-only final analysis and record.
Tech stack: existing E: Python/PyFluent/h5py/numpy/matplotlib, D: Fluent261, Foundation14 WSL.
Spec: SPEC.md. All constraints and numerical limits there apply unchanged.

Review focus: heat sign and enthalpy datum; per-depth normalization; native integrated versus flux-density units; missing/duplicate restart steps; output-only controls changing neither initial fields nor numerical settings.

- [x] Task1 evidence inventory and primary-source check. Baseline/public persistent hashes and32-state counts saved. No all-step energy series in originals.
- [ ] Task2 tested ledger helpers. Write synthetic conservation/datum-shift/depth/restart/coverage tests, observe RED, implement minimal helpers, observe GREEN. Explicit right-endpoint sums; zero denominator and nonfinite rejection. New files test_jet_energy_v001.py, jet_energy_ledger_v001.py.
- [ ] Task3 finite qualifications. Prepare byte-verified exact initial copies and unchanged settings; new native dirs fluent-qualification-v001/openfoam-qualification-v001. Freeze sources, native inputs and parent hashes. Tests for step parsers/restart/output-only changes. Guarded two-rank20-step F then O, closed full native audits with actual reports; no full job from mere solver End. Real paired qualification images and numerical plots inspected before admission.
- [ ] Task4 bounded full admission. New production drivers/auditors and declared exact20->600 global steps. Source/prelaunch/restart gates,580 new steps sequential F then O. Allstep thermal states/reports, complete merged audits, no duplicated boundary rows. Preserve any miss, diagnose native evidence, new version only.
- [ ] Task5 presentation and existing-study26 record. Actual source-bound fixed scales, wall Q and cumulative/per-step ledger plots; source checks and every encoded/decoded frame if media. Compare all32 old states at original times, disclose separate original frozen flows. Update navigation only with exact backups; baseline, manifest and private-file checks.
- [ ] Task6 one fresh whole-improvement review. Apply one Critical/Important fix pass with tests, ledger Minor. Final source/preservation/exclusion/process checks and completion proof before claiming complete. Show actual checked images/media. No GitHub upload without authorization/remote verification.

Interfaces: Task2 consumes explicit step/time/stored-energy/patch-input arrays fromTask3/4 native readers; arrays MUST include initial E and each real step, no interpolation. Task3 acceptance is the admission gate forTask4. Task5 consumes only accepted full selection. Task6 reviews the entire source-to-visual-to-record chain.
