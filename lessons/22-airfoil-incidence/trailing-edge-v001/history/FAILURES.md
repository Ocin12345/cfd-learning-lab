# Preserved unsuccessful stages

- An initial driver attempt failed on Python local scope before a native solver/output case was created.
- The first local split mesh failed native quality at 0.115009, below the unchanged 0.20 floor. No flow ran on it. A new quality-aware transition closure repaired the mesh.
- An initial quality predictor failed its native comparison. Its log was prematurely called GREEN before checking the exit code; the ledger explicitly corrects that statement. Failed source/logs remain preserved. The corrected predictor does not replace native quality checks.
- Both OpenFOAM 3,000-iteration endpoints missed residual and force stability. Native geometry, fields, settings and mass passed. Exact checkpoint continuations preserved these misses.
- The finer 24,000 endpoint passed mass and force stability but missed Uy residual: 2.79768610138e-7 against 1e-7. Its exact checkpoint was copied into a new directory and continued to a finite 32,000 endpoint without changing methods or limits.

Native unsuccessful cases and full logs stay in the external archive. JSON here retains failed-check labels. Historical source snapshots are archival text, not selected executable versions.
