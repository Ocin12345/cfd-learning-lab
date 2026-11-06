# Paired fine-airfoil comparison implementation plan

> For agentic workers: execute one case at a time with independent source and final-result review. User authorization permits automatic routine work without lesson checkpoints.

**Goal:** Compute the existing fine Re=100, zero-incidence airfoil in both solvers and measure medium-to-fine sensitivity of the remaining drag/lift discrepancy.

**Architecture:** Reuse the verified radius-20 fine mesh, continue the existing fine Fluent saved fields in a fresh session/output, and run OpenFOAM cold on the exact same imported mesh. Compare both with their original radius-20 medium results; keep all prior cases and failed diagnostics.

**Tech stack:** Existing E: Python/PyFluent, D: Fluent Student v261, Foundation 14 under WSL, NumPy/Matplotlib/HDF5, Git.

**Spec:** Owner's continuing paired-model request and docs/AIRFOIL_DOMAIN_STUDY.md. The controlled case choice is Codex's decision, not an independent owner contribution.

## Design and alternatives

The radius extension changed drag by about -0.30% in both codes but retained a 1.755% code gap. Saved-field force accounting locates about 87.7% of that gap in the viscous contribution without identifying its cause. Refining the existing mesh family tests discretization sensitivity directly. A mirrored mesh would also change construction and symmetry; a nonzero-angle comparison would add a physical parameter before the zero-angle diagnostic is understood. Those alternatives remain subsequent work.

Use fine055: 944,026 triangles, 473,166 nodes, Re=100, alpha=0, R/c=20. Source mesh SHA256 is 6bf10e6d695f80e27f64f510e21de7d628b0e1f70dec560fd85e4151853458d9. Setup source_mesh_sha256 describes upstream geometry lineage, not the emitted fine mesh. The old fine Fluent endpoint at iteration172 is a reused initial state, not a new independent cold calculation. Report the new continuation and this provenance explicitly.

## Global constraints

- All project outputs, native fields, environments and supported temp/cache paths remain on E:; originals remain unchanged.
- Use a fresh native/runs/2026-11-01-airfoil-fine-paired-v001 and local/airfoil-fine-paired-v001. Public large native solver files remain excluded.
- Use one hidden Fluent session with 2 processes. Run two-process solves sequentially. Require current available RAM >=6 GiB before either fine solve and E: free >=10 GiB; do not lower gates to launch. Preserve saved states before exit.
- Student mesh admission: cells and nodes <=950,000; current fine source qualifies. Query licence before launch and record unsupported-status responses honestly; actual successful Student startup establishes token availability before flow iterations.
- Fluent continuation restores native 1e-8 stopping controls, then records 100 additional observations. Acceptance remains p/U residuals <=1e-7, last20 Cd relative range <=1e-4, Cl absolute range <=1e-5, |Cl|<=0.001, mass imbalance <=1e-6, orthogonal quality >=0.2.
- OpenFOAM retains cold uniform U=(1,0,0), p=0; identical methods/references and U relaxation0.7 to2000,0.9 to2500 and5000. Use full checkMesh -allGeometry -allTopology. Acceptance remains max initial p/U <=1e-7, last100 Cd range <=1e-4, Cl absolute range <=1e-5, |Cl|<=0.001, mass imbalance <=1e-6.
- Preserve the separate0.5% code-agreement target and all failed diagnostic flags. No physical validation, GCI, convergence-order or mesh/domain-independence claim.
- Frozen prior analyzers/helpers remain unchanged. Use the versioned coordinate-bijection analyzer/plotter and unchanged1e-10 geometry tolerance.

## Review focus

Check exact emitted mesh versus upstream lineage hash; saved-field initialization and new iteration IDs; pressure maximum across all correctors; finite fields/topology/mass and every native force sample; unchanged methods/references/thresholds; component cancellation, false diagnostic flags and honest owner/Codex/solver roles.

## Execution

- [ ] Bind source mesh/setup and Fluent case/data/settings hashes. Recheck current licence, disk, RAM and single-session admission immediately before launch.
- [ ] Run fresh hidden Fluent saved-field continuation; save native case/data/images and checks, close, show actual image.
- [ ] Prepare exact fine OpenFOAM mesh/case, full native mesh checks, then run recorded two-rank cold-start/continuation stages. Save before closure; extend in separately recorded stages only if convergence/stability remains incomplete.
- [ ] Analyze actual saved geometry/fields/mass/forces using v002 helper. Re-evaluate pressure/viscous forces in isolated saved-field scratch with zero flow iterations and unchanged1e-8 equivalence bound.
- [ ] Aggregate original medium and new fine paired results with matching physics/methods and source hashes. Test any new aggregation logic against real evidence and refusal guards; retain failures.
- [ ] Produce common-scale fields and academic sensitivity/history charts with PDF/SVG/600dpiPNG, source snapshots/hashes and rendered PDF inspection. Show actual outputs.
- [ ] Save methods, code, results, limits and contribution roles; run final source/hash/privacy/link checks and byte-verify staged contents. No publication or university submission.
