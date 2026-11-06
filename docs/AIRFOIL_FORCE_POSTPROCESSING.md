# Saved-field airfoil force decomposition

`scripts/airfoil_force_components.py` prepares an isolated scratch case and verifies pressure/viscous force outputs. It never starts postprocessing during preparation. Native cases remain unchanged.

Foundation 14's [foamPostProcess](https://cpp.openfoam.org/v14/foamPostProcess_8C.html) can instantiate the incompressible solver model and evaluate function objects at selected saved times without a flow-solving loop. The [force-function guide](https://doc.cfd.direct/openfoam/user-guide-v14/post-processing-functionality) distinguishes kinematic-pressure incompressible forces. The scratch uses the actual original density, velocity, area, chord and pressure references, with the same momentum-transport model and schemes.

After R40 execution and reconstruction are finished at 5000, prepare a **fresh** scratch:

```powershell
E:/CFD-Learning-Lab/environments/fluent/Scripts/python.exe -B E:/CFD-Learning-Lab/repo/cfd-learning-lab/scripts/airfoil_force_components.py prepare --case E:/CFD-Learning-Lab/native/runs/2026-11-01-airfoil-domain40-v002/openfoam --time 5000 --out E:/CFD-Learning-Lab/native/postprocessing/2026-11-01-airfoil-forces-r40-v001
```

Record host available RAM >= 4 GiB and E free space >= 10 GiB. Run only one serial postprocessor at a time, assuming up to 2 GiB extra RAM:

```powershell
wsl -d CFD-Foundation14 -- bash /mnt/e/CFD-Learning-Lab/native/postprocessing/2026-11-01-airfoil-forces-r40-v001/run-postprocess.sh
```

Set the launcher working directory to scratch as well; the generated shell explicitly enters it before OpenFOAM utilities. Afterwards record available RAM >= 2 GiB and the 10 GiB disk floor.

The helper copies only explicit saved fields, topology, dictionaries and evidence, checks original/copy SHA256 values before and after, and requires one selected processing time with no equation-solving log evidence. Both pressure+viscous totals and a coefficient replay must reproduce saved Cd/Cl/Cm within an absolute 1e-8 equivalence tolerance. A mismatch is retained and verification fails. Outputs include `force-components.json`, CSV, raw function output, runtime source hashes and the complete log.

This diagnostic locates differences between force contributions. It does not establish their cause, accept failed convergence/symmetry checks, or provide physical validation.

## Aggregate records

`scripts/analyse_airfoil_force_components.py` reads the Fluent native scalar pressure/viscous/total force reports, the verified OpenFOAM decomposition and the matched-pair analysis. It checks source hashes, all Fluent observations, reference normalization, the OpenFOAM raw vectors/replay and selected times. Saved HDF fields are bound to the prior matched-pair native-field checks by their exact hashes; this standard-library aggregate does not parse HDF itself.

Its required flags are `--fluent-run`, `--openfoam-components`, `--comparison-dir` and `--out`. The output directory must be fresh under project `local` or repository `results`.

Outputs are `comparison.json`, `force-components.csv`, `force-differences.csv`, `source-manifest.json` and `executed-analyser.py`. The component CSV has one row per solver and pressure/viscous/total term. The difference CSV uses Fluent minus OpenFOAM and signed shares of the total gap. Shares can be negative or exceed 100% when terms cancel; a zero total gap has no defined share. Fluent report rounding limits closure and share precision, so publication should round shares appropriately.

The aggregate retains the main `check_results`, `failed_checks`, `numerical_checks_passed`, `numerical_convergence_passed`, `symmetry_diagnostic_passed` and `cross_code_target_pass` fields. `evidence_consistency_checks_passed` describes data consistency, not acceptance of the underlying solutions. Numerical cause and physical validation remain explicitly unestablished.
