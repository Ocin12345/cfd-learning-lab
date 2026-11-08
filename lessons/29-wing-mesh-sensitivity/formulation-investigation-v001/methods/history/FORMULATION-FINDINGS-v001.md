# Wing solver-formulation investigation, v001

The owner explicitly chose to investigate the solver formulation after the three failed repair hypotheses. This is a read-only analysis of existing native results; no new Fluent session, solver iteration, changed numerical control or OpenFOAM rerun occurred.

## What the native evidence establishes

A tested strict parser read all400new coupled-solve blocks from the four100iteration pilots. In every block, the printed initial coupled-continuity value matches the unscaled HDF residual numerator, and the reported outer residual matches that numerator divided by its stored scale (within the precision of printed tables). The continuity scale remains exactly2.2344914086503613 throughout; no parser/scaling reset explains the miss. In all20saved checkpoints, the final printed linear continuity residual matches the independently reconstructed sum of absolute native cell mass imbalances from signed face fluxes (within printed precision).

The tighter AMG solve therefore works on the assembled linear system: original-control cases use2cycles while the tight cases use4..5. In the tight-AMG/Courant20case, the next iteration starts with continuity imbalance a median2595.55times the preceding inner terminal value (range1048.58..3886.70). The other tight cases also rebound by more than1000. The original-control median is57.07. This is an observed distinction between assembled/initial and corrected/final stages in THESE logs, not a universal claim that every Fluent residual is sampled at the same stage.

Across20checkpoints, the largest volume-weighted velocity L2 change relative to the original1500field is0.00010457282092648344%; the largest pressure RMS change is2.156773568751815e-7Pa. Small field/force changes and good corrected mass closure do not replace the original1e-5 residual test: every pilot remains MISS. These steady iteration oscillations do not establish physical unsteadiness or round-off as the cause.

## What remains unknown

No per-cell residual fields were retained in the saved data, and save_cell_residuals isfalse in all recorded settings. We cannot spatially localize the PRE-correction residual from corrected flux fields, pressure differences, or a visually smooth wing. Coupled operator reassembly/pressure-velocity-flux consistency is the remaining investigation focus. Mesh/skewness, interpolation, the flux reconstruction, or an implementation detail have NOT been proved as the unique cause. No direct transition to another algorithm is yet justified by a localized cause.

## Proposed bounded next diagnostic

Use the unchanged original1500checkpoint and original physical/spatial/reference/numerical controls. Enable residual retention and output verbosity only. First close and inspect a ONE-iteration1501probe, proving initial-field/settings/source-prefix identity and obtaining actual per-cell residual fields. If that retention/identity audit passes, admit a NEW exact-checkpoint NINE-iteration diagnostic to1510, saving every iteration. Map the retained continuity and momentum residuals, compare their actual sums/stages with the matrix transcript and corrected native flux, and report which mesh/boundary regions dominate. No new convergence qualification, acceptance reset, full run or physics claim follows this diagnostic. The finite specification is RETAINED-RESIDUAL-PROPOSAL-v001.md; implementation/launch have not occurred.

## Primary sources verified

- [Fluent v261 Monitoring Solution Convergence](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_sec_solve_monitor.html): defines raw/scaled residual storage and the startup-dependent continuity scale. Its residual-postprocessing section requires enabling retained cell residuals and advancing at least one iteration; the postprocessed values are unscaled. This supports the diagnostic, not changing the acceptance norm.
- [Fluent v261 Pressure-Velocity Coupling theory](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_uns_solve_pv.html): the Coupled formulation assembles pressure-based continuity and momentum together with implicit pressure/face-flux coupling; segregated correction methods use another solution architecture. This is context for the observed stages, not proof of a defect.
- [Fluent v261 Algebraic Multigrid theory](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_th/flu_th_sec_amg_multigrid.html): the AMG system is linearized; nonlinear effects enter when the fine-level operator is updated. This distinguishes inner-system reduction from outer convergence.
- [Fluent v261 Judging Convergence](https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/flu_ug/flu_ug_sec_solve_conv_stab.html): residuals and relevant integrated quantities should both be considered. We retain the original stricter1e-5 requirement; no relaxation or renormalization is proposed.

## Reproduction and attribution

Run investigate_formulation_v001.py with the existing E:environment and explicit wing_import_paths_v001 module. Its strict parser has7tests watched RED then GREEN; all19helper/control/deployment/transcript tests pass together. CSVs retain400actual matrix/outer rows and20actual saved-field cross-checks. investigation.json and persistent-source-bindings.json record the evidence and hashes; figures are under exports/figures/2026-11-08-wing-formulation-investigation-v001. Actual PNG inspected. Native inputs are rehashed before/after analysis. Allpriorfailures and repository content preserved; no upload.

Owner selected this investigation and authorized the portfolio. Codex performed the read-only analysis and plotting. ANSYS computed the existing fields/transcripts; no independent owner execution or reflection is invented. This is not wing convergence completion, mesh/time/domain independence or physical validation. The separate zero-incidence MISS remains open. No whole-improvement review was spawned for this read-only subtask.
